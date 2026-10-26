"""Iterative refinement driven by the lower bound.

Idea: the lower-bound tour tells us, for every target, which heading interval
the best tours seem to need. So we spend extra resolution exactly there.

Iteration 0: solve both bounds on the uniform grids from config.yml.

Iteration k = 1 .. mu, for every target t:
  1. Take the interval I_t chosen at t by the previous lower-bound tour.
  2. Place lambda evenly spaced headings strictly inside I_t:
         start + j * width / (lambda + 1),   j = 1 .. lambda
  3. Upper bound: add these headings to t's heading set.
  4. Lower bound: split I_t at the same headings into lambda + 1 pieces.
  Then solve both bounds again.

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
from dataclasses import dataclass, field

from .bounds import (BoundResult, compute_lower_bound, compute_upper_bound,
                     initial_headings, initial_intervals)
from .config import Config
from .tsp.gtsp import IntervalEdgeCosts
from .utils.report import LB_WARNING
from .utils.angles import HeadingInterval, mod2pi


@dataclass
class Iteration:
    index: int
    upper: BoundResult
    lower: BoundResult
    best_upper: float          # smallest upper bound so far (every UB is a real tour)
    max_lb_estimate: float     # largest LB estimate so far; NOT used for the gap, since
                               # an LKH estimate can be too high
    headings: list[list[float]]                  # heading set of every target
    intervals: list[list[HeadingInterval]]       # interval partition of every target
    added: list[list[float]] = field(default_factory=list)  # headings added per target

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

        pieces = [HeadingInterval(float(mod2pi(chosen.start + j * step)), step)
                  for j in range(lam + 1)]
        k = new_intervals[t].index(chosen)
        new_intervals[t][k:k + 1] = pieces
    return new_headings, new_intervals, added


def run(cfg: Config, lam: int, mu: int, on_iteration=None) -> list[Iteration]:
    """Iteration 0 plus `mu` refinement iterations; calls on_iteration(it) after each."""
    headings, intervals = initial_headings(cfg), initial_intervals(cfg)
    upper_cache = {}
    lower_costs = IntervalEdgeCosts(cfg.points, cfg.turning_radius)
    iterations: list[Iteration] = []
    added: list[list[float]] = [[] for _ in cfg.points]

    for k in range(mu + 1):
        print(f"iteration {k}/{mu}: upper bound ...", file=sys.stderr)
        upper = compute_upper_bound(cfg, headings, upper_cache)
        print(f"iteration {k}/{mu}: lower bound ...", file=sys.stderr)
        lower = compute_lower_bound(cfg, intervals, lower_costs)

        best_upper = min([upper.value] + [it.best_upper for it in iterations])
        max_lb = max([lower.value] + [it.max_lb_estimate for it in iterations])
        it = Iteration(k, upper, lower, best_upper, max_lb,
                       [list(h) for h in headings], [list(iv) for iv in intervals], added)
        iterations.append(it)
        if it.lb_exceeds_ub:
            print(LB_WARNING.format(k=k, lb=lower.value, ub=best_upper), file=sys.stderr)
        if on_iteration:
            on_iteration(it)

        if k < mu:
            headings, intervals, added = refine(headings, intervals, lower.selected, lam)
    return iterations
