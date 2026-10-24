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
from dataclasses import dataclass, field

import numpy as np

from .config import Config
from .tsp import noon_bean, solvers
from .tsp.gtsp import (ClusteredGraph, IntervalEdgeCosts, build_heading_graph,
                       build_interval_graph, refresh_edge, uniform_headings, uniform_intervals)
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
    kind: str                 # "upper" or "lower"
    value: float
    is_valid_bound: bool      # False only for a lower bound from a heuristic solver
    solver: str
    options_per_target: list[int]
    order: list[int]          # targets in visiting order
    legs: list[Leg]
    selected: list            # chosen heading / interval of every target (index = target)
    build_seconds: float
    solve_seconds: float
    notes: list[str] = field(default_factory=list)
    solve_rounds: int = 1        # > 1 when lazy edge evaluation had to re-solve
    interval_solves: int = 0     # interval problems actually solved for this bound
    # seconds per stage: dubins | interval | transform | lkh | exact | other
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


def _refresh_tour_neighbourhood(graph: ClusteredGraph, chosen: list[int],
                                edge_costs: IntervalEdgeCosts) -> bool:
    """If the tour uses an estimated edge, make a batch of edges exact and
    return True (solve again). Return False when the tour is all exact."""
    tour_edges = [(chosen[k], chosen[(k + 1) % len(chosen)]) for k in range(len(chosen))]
    if all(graph.exact[u, v] for u, v in tour_edges):
        return False
    for u, v in _lazy_batch(graph, tour_edges):
        refresh_edge(graph, edge_costs, u, v)
    return True


def _lazy_batch(graph: ClusteredGraph, tour_edges) -> set[tuple[int, int]]:
    """Estimated edges to make exact after one lazy round.

    Only fixing the tour's own edges would let the solver swap to a
    neighbouring estimated edge next round, and the next, one at a time. So
    for every tour edge u -> v we also fix the edges u -> (any option of v's
    target) and (any option of u's target) -> v: the alternatives the solver
    would most likely try next.
    """
    batch = set()
    for u, v in tour_edges:
        tu, tv = graph.node_target[u], graph.node_target[v]
        batch.update((u, w) for w in graph.clusters[tv])
        batch.update((w, v) for w in graph.clusters[tu])
    return {(a, b) for a, b in batch if not graph.exact[a, b]}


def _solve_graph(graph: ClusteredGraph, kind: str, cfg: Config,
                 build_seconds: float, edge_costs: IntervalEdgeCosts | None = None):
    """Noon-Bean + LKH + decode.

    Lazy edge evaluation (lower bound only): some graph costs may be lower
    estimates (graph.exact == False). Solve with them; if the tour uses any
    estimated edge, compute its exact cost and solve again. Stop when the tour
    uses exact edges only. With an optimal solver that tour would be optimal
    for the exact costs too:
        estimated optimum  <=  exact optimum  <=  exact cost of this tour
                                              =   its estimated cost = estimated optimum.
    LKH is not optimal, so here it is only LKH's best tour on the exact costs.

    Screening: the lazy rounds only decide which edges to make exact, so they
    use a cheap LKH (cfg.lkh_runs_screening restarts). Once a screening tour
    uses exact edges only, solve once more with the full cfg.lkh_runs. If that
    tour is all exact too, stop; otherwise make its edges exact and go back
    to screening.
    """
    is_lower = kind == "lower"
    lazy = graph.exact is not None and not graph.exact.all()
    rounds = 0
    timings = {"transform": 0.0, "lkh_screen": 0.0, "lkh": 0.0}
    t0 = time.perf_counter()

    def solve_once(runs, stage):
        t = time.perf_counter()
        int_cost = _integer_costs(graph, cfg.cost_scale, round_down=is_lower)
        atsp = noon_bean.transform(int_cost, graph.clusters)
        timings["transform"] += time.perf_counter() - t
        t = time.perf_counter()
        tour = solvers.solve_lkh(atsp, runs)
        timings[stage] += time.perf_counter() - t
        return atsp, tour, noon_bean.decode(tour, atsp)

    no_screening = not lazy or cfg.lkh_runs_screening >= cfg.lkh_runs
    final = no_screening
    while True:
        rounds += 1
        if final:
            atsp, tour, chosen = solve_once(cfg.lkh_runs, "lkh")
        else:
            atsp, tour, chosen = solve_once(cfg.lkh_runs_screening, "lkh_screen")
        if lazy and _refresh_tour_neighbourhood(graph, chosen, edge_costs):
            final = no_screening        # tour used estimated edges: screen again
        elif final:
            break
        else:
            final = True                # screening tour is all exact: full solve
    solve_seconds = time.perf_counter() - t0

    n = len(chosen)
    legs = []
    for idx in range(n):
        u, v = chosen[idx], chosen[(idx + 1) % n]
        legs.append(Leg(
            from_target=graph.node_target[u], to_target=graph.node_target[v],
            depart_heading=float(graph.depart_heading[u, v]),
            arrive_heading=float(graph.arrive_heading[u, v]),
            length=float(graph.cost[u, v]),
            from_option=graph.label(u), to_option=graph.label(v),
        ))
    selected = [None] * graph.num_targets
    for node in chosen:
        selected[graph.node_target[node]] = graph.node_option[node]

    notes = []
    if is_lower:
        # Integer optimum / scale <= real GTSP optimum <= true Dubins TSP optimum,
        # but LKH's tour may cost more than the integer optimum: an estimate only.
        value = atsp.gtsp_cost(tour) / cfg.cost_scale
        valid = False
        notes.append("solved with LKH (heuristic): this is an ESTIMATE, not a certified "
                     "lower bound. If LKH misses the optimal GTSP tour it is too high.")
        notes.append("legs may use different headings when arriving at and leaving a "
                     "target, so this tour is not flyable (that is why it is a lower bound).")
    else:
        value = sum(leg.length for leg in legs)  # real length of a flyable tour
        valid = True

    return BoundResult(
        kind=kind, value=value, is_valid_bound=valid, solver="lkh",
        options_per_target=[len(c) for c in graph.clusters],
        order=[leg.from_target for leg in legs], legs=legs, selected=selected,
        build_seconds=build_seconds, solve_seconds=solve_seconds, notes=notes,
        solve_rounds=rounds, timings=timings,
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


def compute_upper_bound(cfg: Config, headings: list[list[float]] | None = None,
                        cache: dict | None = None) -> BoundResult:
    """Upper bound with the given headings per target (default: uniform grid)."""
    t0 = time.perf_counter()
    headings = headings or initial_headings(cfg)
    graph = build_heading_graph(cfg.points, headings, cfg.turning_radius, cache)
    result = _solve_graph(graph, "upper", cfg, time.perf_counter() - t0)
    result.timings["dubins"] = result.build_seconds   # building = Dubins path lengths
    _finish_timings(result)
    return result


def compute_lower_bound(cfg: Config, intervals: list[list[HeadingInterval]] | None = None,
                        costs: IntervalEdgeCosts | None = None) -> BoundResult:
    """Lower bound with the given intervals per target (default: uniform partition).

    Pass the same `costs` object between refinement iterations to reuse every
    interval problem solved so far.
    """
    t0 = time.perf_counter()
    intervals = intervals or initial_intervals(cfg)
    if costs is None:
        costs = IntervalEdgeCosts(cfg.points, cfg.turning_radius, cfg.sweep_samples)
        costs.register_initial(intervals)
    solves_before, seconds_before = costs.solves, costs.solve_seconds
    graph = build_interval_graph(cfg.points, intervals, cfg.turning_radius, costs,
                                 lazy=cfg.lower.lazy, progress=_progress)
    result = _solve_graph(graph, "lower", cfg, time.perf_counter() - t0, costs)
    result.interval_solves = costs.solves - solves_before
    # interval solves happen both while building and in the lazy re-solve rounds
    result.timings["interval"] = costs.solve_seconds - seconds_before
    _finish_timings(result)
    return result
