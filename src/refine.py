"""Iterative refinement driven by the lower bound.

Idea: the lower-bound tour tells us, for every target, which heading interval
the best tours seem to need. So we spend extra resolution exactly there.

Iteration 0: solve both bounds on the uniform grids from config.yml.

Iteration k = 1, 2, ... for every target t:
  1. Take the interval I_t chosen at t by the previous lower-bound tour.
  2. Place lambda evenly spaced headings strictly inside I_t:
         start + j * width / (lambda + 1),   j = 1 .. lambda
  3. Upper bound: add these headings to t's heading set.
  4. Lower bound: split I_t at the same headings into lambda + 1 pieces.
  Then solve both bounds again.

Stop at the first iteration whose gap is <= the gap target (percent).

Why both are refined: if only headings were added, the next lower bound
would pick the same intervals and every later iteration would add the same
headings again. Splitting the intervals tightens the lower bound and lets the
next iteration "zoom in" further (or move to another interval if the
lower-bound tour changes).

With optimal GTSP solutions both bounds would be monotone:
  * the heading sets only grow  -> the upper bound can only go down;
  * each new interval is inside an old one, so its edge costs can only go up
                                -> the lower bound can only go up.
Both are solved with LKH, which can break this. Every upper bound is a real
tour, so the smallest one seen is kept. The lower bound is only an estimate
(LKH may be above the GTSP optimum), so the gap always uses the current
iteration's estimate, never the maximum over iterations.
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass, field

from . import parallel
from .bounds import BoundResult, compute_bounds, initial_headings, initial_intervals
from .config import Config
from .tsp.gtsp import IntervalEdgeCosts
from .utils.angles import HeadingInterval, mod2pi
from .utils.report import LB_WARNING


@dataclass
class Iteration:
    index: int
    upper: BoundResult
    lower: BoundResult
    best_upper: float  # smallest upper bound so far (every UB is a real tour)
    max_lb_estimate: float  # largest LB estimate so far; NOT used for the gap, since
    # an LKH estimate can be too high
    headings: list[list[float]]  # heading set of every target
    intervals: list[list[HeadingInterval]]  # interval partition of every target
    added: list[list[float]] = field(default_factory=list)  # headings added per target
    wall_seconds: float = 0.0  # elapsed time for both bounds (stages may overlap)

    @property
    def gap_percent(self) -> float:
        """(best upper - this iteration's LB estimate) / this iteration's LB estimate."""
        return (self.best_upper - self.lower.value) / self.lower.value * 100

    @property
    def lb_exceeds_ub(self) -> bool:
        """An LB estimate above a real tour's length means LKH missed the GTSP optimum."""
        return self.lower.value > self.best_upper + 1e-9


def _same_heading(a: float, b: float, tol: float = 1e-9) -> bool:
    diff = abs(float(mod2pi(a - b)))
    return diff < tol or diff > 2 * math.pi - tol


def refine(headings, intervals, selected, lam: int):
    """One refinement step. Returns new heading sets, new intervals, added headings."""
    new_headings = [list(h) for h in headings]
    new_intervals = [list(iv) for iv in intervals]
    added = []
    for t, chosen in enumerate(selected):
        step = chosen.width / (lam + 1)
        cuts = [float(mod2pi(chosen.start + j * step)) for j in range(1, lam + 1)]

        fresh = [c for c in cuts if not any(_same_heading(c, h) for h in new_headings[t])]
        new_headings[t] = sorted(new_headings[t] + fresh)
        added.append(fresh)

        pieces = [
            HeadingInterval(float(mod2pi(chosen.start + j * step)), step) for j in range(lam + 1)
        ]
        k = new_intervals[t].index(chosen)
        new_intervals[t][k : k + 1] = pieces
    return new_headings, new_intervals, added


def run(
    cfg: Config,
    lam: int,
    gap_target: float | None,
    on_iteration=None,
    max_iterations: int | None = None,
) -> list[Iteration]:
    """Refine until the gap (%) is <= gap_target; calls on_iteration(it) after each iteration.

    gap_target None: iteration 0 only. max_iterations caps the loop (used by
    the tests; the pipeline has no cap).
    """
    headings, intervals = initial_headings(cfg), initial_intervals(cfg)
    upper_cache = {}
    lower_costs = IntervalEdgeCosts(cfg.points, cfg.turning_radius)
    iterations: list[Iteration] = []
    added: list[list[float]] = [[] for _ in cfg.points]

    with parallel.pool(cfg.use_all_cores) as pool:
        k = 0
        while True:
            print(f"iteration {k}: bounds ...", file=sys.stderr)
            t0 = time.perf_counter()
            upper, lower = compute_bounds(cfg, headings, intervals, upper_cache, lower_costs, pool)
            wall = time.perf_counter() - t0

            best_upper = min([upper.value] + [it.best_upper for it in iterations])
            max_lb = max([lower.value] + [it.max_lb_estimate for it in iterations])
            it = Iteration(
                k,
                upper,
                lower,
                best_upper,
                max_lb,
                [list(h) for h in headings],
                [list(iv) for iv in intervals],
                added,
                wall,
            )
            iterations.append(it)
            if it.lb_exceeds_ub:
                print(LB_WARNING.format(k=k, lb=lower.value, ub=best_upper), file=sys.stderr)
            if on_iteration:
                on_iteration(it)

            print(f"iteration {k}: gap {it.gap_percent:.3f}%", file=sys.stderr)

            if gap_target is None or it.gap_percent <= gap_target:
                return iterations
            if max_iterations is not None and len(iterations) >= max_iterations:
                return iterations
            headings, intervals, added = refine(headings, intervals, lower.selected, lam)
            k += 1
