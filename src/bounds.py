"""Upper and lower bounds on the Dubins TSP.

Both bounds use the same pipeline:

    clustered graph  ->  integer costs  ->  Noon-Bean ATSP  ->  solve  ->  decode

They differ only in what a node means and what an edge costs:

* Upper bound: node = one discrete heading; edge = point-to-point Dubins path.
  The decoded tour arrives at and leaves every target with the same heading,
  so it is a real flyable tour  ->  length >= optimum, whatever solver is used.

* Lower bound: node = a heading interval; edge = interval-problem minimum.
  Every leg may use a different heading at each end (the tour can "kink" at
  targets), so it is not flyable, and its *optimal* value is <= the optimum.
  The GTSP is solved with LKH, a heuristic: it may miss the optimal GTSP
  tour, and then the value it reports is too high. So the lower bound is an
  ESTIMATE, not a certified bound (the exact MILP solver was too slow).
"""

from __future__ import annotations

import sys
import time
from concurrent.futures import Future
from dataclasses import dataclass, field

import numpy as np

from . import parallel
from .config import Config
from .tsp import noon_bean, solvers
from .tsp.gtsp import (
    ClusteredGraph,
    IntervalEdgeCosts,
    build_heading_graph,
    build_interval_graph,
    uniform_headings,
    uniform_intervals,
)
from .utils.angles import HeadingInterval


@dataclass
class Leg:
    from_target: int
    to_target: int
    depart_heading: float
    arrive_heading: float
    length: float
    from_option: str
    to_option: str


@dataclass
class BoundResult:
    kind: str  # "upper" or "lower"
    value: float
    is_valid_bound: bool  # True for the upper bound; the LKH lower bound is an estimate
    status: str  # "valid upper bound" or "estimate (LKH)"
    solver: str
    options_per_target: list[int]
    order: list[int]  # targets in visiting order
    legs: list[Leg]
    selected: list  # chosen heading / interval of every target (index = target)
    build_seconds: float
    solve_seconds: float
    notes: list[str] = field(default_factory=list)
    interval_solves: int = 0  # interval problems actually solved for this bound
    # seconds per stage: dubins | interval | transform | lkh | other
    timings: dict[str, float] = field(default_factory=dict)


def _progress(done: int, total: int):
    if sys.stderr.isatty() or done == total:
        sys.stderr.write(f"\r  interval costs: {done}/{total} target pairs")
        if done == total:
            sys.stderr.write("\n")
        sys.stderr.flush()


def _integer_costs(graph: ClusteredGraph, scale: float, round_down: bool) -> np.ndarray:
    """Scale to integers for the ATSP solver.

    For the lower bound we round *down*, so the integer problem can only be
    cheaper than the real one and the bound stays valid.
    """
    scaled = np.where(np.isfinite(graph.cost), graph.cost * scale, 0.0)
    return (np.floor(scaled) if round_down else np.rint(scaled)).astype(np.int64)


def _timed_lkh(atsp: noon_bean.NoonBeanATSP, runs: int):
    """LKH solve and its duration (runs in a worker process with use_all_cores)."""
    t = time.perf_counter()
    tour = solvers.solve_lkh(atsp, runs)
    return tour, time.perf_counter() - t


@dataclass
class _Started:
    """A bound whose LKH solve has been started (possibly in a worker process)."""

    graph: ClusteredGraph
    kind: str
    atsp: noon_bean.NoonBeanATSP
    build_seconds: float
    transform_seconds: float
    lkh: Future  # -> (tour, LKH seconds)


def _start_solve(graph: ClusteredGraph, kind: str, cfg: Config, build_seconds: float, pool):
    """Integer costs + Noon-Bean, then start LKH (cfg.lkh_runs restarts)."""
    t0 = time.perf_counter()
    int_cost = _integer_costs(graph, cfg.cost_scale, round_down=kind == "lower")
    atsp = noon_bean.transform(int_cost, graph.clusters)
    transform_seconds = time.perf_counter() - t0
    lkh = parallel.submit(pool, _timed_lkh, atsp, cfg.lkh_runs)
    return _Started(graph, kind, atsp, build_seconds, transform_seconds, lkh)


def _finish_solve(started: _Started, cfg: Config) -> BoundResult:
    """Wait for LKH, decode the tour and build the result."""
    graph, kind, atsp = started.graph, started.kind, started.atsp
    is_lower = kind == "lower"
    tour, lkh_seconds = started.lkh.result()
    t = time.perf_counter()
    chosen = noon_bean.decode(tour, atsp)
    timings = {
        "transform": started.transform_seconds + time.perf_counter() - t,
        "lkh": lkh_seconds,
    }
    # Stage times only: waiting for a solve that ran in parallel is not counted.
    solve_seconds = sum(timings.values())

    n = len(chosen)
    legs = []
    for idx in range(n):
        u, v = chosen[idx], chosen[(idx + 1) % n]
        legs.append(
            Leg(
                from_target=graph.node_target[u],
                to_target=graph.node_target[v],
                depart_heading=float(graph.depart_heading[u, v]),
                arrive_heading=float(graph.arrive_heading[u, v]),
                length=float(graph.cost[u, v]),
                from_option=graph.label(u),
                to_option=graph.label(v),
            )
        )
    selected = [None] * graph.num_targets
    for node in chosen:
        selected[graph.node_target[node]] = graph.node_option[node]

    notes = []
    if is_lower:
        # Integer optimum / scale <= real GTSP optimum <= true Dubins TSP optimum,
        # but LKH's tour may cost more than the integer optimum: an estimate only.
        value = atsp.gtsp_cost(tour) / cfg.cost_scale
        valid = False
        notes.append(
            "solved with LKH (heuristic): this is an ESTIMATE, not a certified "
            "lower bound. If LKH misses the optimal GTSP tour it is too high."
        )
        notes.append(
            "legs may use different headings when arriving at and leaving a "
            "target, so this tour is not flyable (that is why it is a lower bound)."
        )
    else:
        value = sum(leg.length for leg in legs)  # real length of a flyable tour
        valid = True

    return BoundResult(
        kind=kind,
        value=value,
        is_valid_bound=valid,
        status="valid upper bound" if valid else "estimate (LKH)",
        solver="lkh",
        options_per_target=[len(c) for c in graph.clusters],
        order=[leg.from_target for leg in legs],
        legs=legs,
        selected=selected,
        build_seconds=started.build_seconds,
        solve_seconds=solve_seconds,
        notes=notes,
        timings=timings,
    )


def _finish_timings(result: BoundResult) -> None:
    """Put whatever the named stages don't cover under 'other'."""
    total = result.build_seconds + result.solve_seconds
    result.timings["other"] = max(0.0, total - sum(result.timings.values()))


def initial_headings(cfg: Config) -> list[list[float]]:
    """Uniform heading grid at every target, from config."""
    grid = uniform_headings(cfg.upper.resolution, cfg.heading_offset)
    return [list(grid) for _ in cfg.points]


def initial_intervals(cfg: Config) -> list[list[HeadingInterval]]:
    """Uniform interval partition at every target, from config."""
    grid = uniform_intervals(cfg.lower.resolution, cfg.heading_offset)
    return [list(grid) for _ in cfg.points]


def _build_upper(cfg: Config, headings, cache, pool):
    t0 = time.perf_counter()
    headings = headings or initial_headings(cfg)
    graph = build_heading_graph(cfg.points, headings, cfg.turning_radius, cache, pool)
    return graph, time.perf_counter() - t0


def _finish_upper(started: _Started, cfg: Config) -> BoundResult:
    result = _finish_solve(started, cfg)
    result.timings["dubins"] = result.build_seconds  # building = Dubins path lengths
    _finish_timings(result)
    return result


def _build_lower(cfg: Config, intervals, costs: IntervalEdgeCosts, pool):
    t0 = time.perf_counter()
    intervals = intervals or initial_intervals(cfg)
    graph = build_interval_graph(
        cfg.points, intervals, cfg.turning_radius, costs, progress=_progress, pool=pool
    )
    return graph, time.perf_counter() - t0


def _finish_lower(started: _Started, cfg: Config, costs: IntervalEdgeCosts, before) -> BoundResult:
    result = _finish_solve(started, cfg)
    result.interval_solves = costs.solves - before[0]
    result.timings["interval"] = costs.solve_seconds - before[1]
    _finish_timings(result)
    return result


def compute_upper_bound(
    cfg: Config,
    headings: list[list[float]] | None = None,
    cache: dict | None = None,
    pool=None,
) -> BoundResult:
    """Upper bound with the given headings per target (default: uniform grid)."""
    graph, build_seconds = _build_upper(cfg, headings, cache, pool)
    return _finish_upper(_start_solve(graph, "upper", cfg, build_seconds, pool), cfg)


def compute_lower_bound(
    cfg: Config,
    intervals: list[list[HeadingInterval]] | None = None,
    costs: IntervalEdgeCosts | None = None,
    pool=None,
) -> BoundResult:
    """Lower bound with the given intervals per target (default: uniform partition).

    Pass the same `costs` object between refinement iterations to reuse every
    interval problem solved so far.
    """
    costs = costs or IntervalEdgeCosts(cfg.points, cfg.turning_radius)
    before = (costs.solves, costs.solve_seconds)
    graph, build_seconds = _build_lower(cfg, intervals, costs, pool)
    return _finish_lower(_start_solve(graph, "lower", cfg, build_seconds, pool), cfg, costs, before)


def compute_bounds(
    cfg: Config,
    headings: list[list[float]],
    intervals: list[list[HeadingInterval]],
    upper_cache: dict,
    lower_costs: IntervalEdgeCosts,
    pool=None,
) -> tuple[BoundResult, BoundResult]:
    """Both bounds. With a pool, the upper-bound LKH solve runs in a worker
    while this process builds the lower-bound graph and starts its solve."""
    graph, build_seconds = _build_upper(cfg, headings, upper_cache, pool)
    upper = _start_solve(graph, "upper", cfg, build_seconds, pool)
    before = (lower_costs.solves, lower_costs.solve_seconds)
    graph, build_seconds = _build_lower(cfg, intervals, lower_costs, pool)
    lower = _start_solve(graph, "lower", cfg, build_seconds, pool)
    return _finish_upper(upper, cfg), _finish_lower(lower, cfg, lower_costs, before)
