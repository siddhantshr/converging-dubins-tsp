"""Two-point Dubins interval problem.

Problem: the vehicle starts at point A with any heading in interval I0 and must
reach point B with any heading in interval I1. Find the shortest such path.

    L*(A, I0, B, I1) = min over th0 in I0, th1 in I1 of  Dubins(A, th0 -> B, th1)

Why this gives a lower bound for the Dubins TSP: any real tour passes each
target with a heading that lies in *some* interval, and on every leg its length
is >= this minimum (a minimum over a set is <= every member).

How we solve it. The best (th0, th1) lies in one of three places of the
rectangle I0 x I1:

  1. a corner: both headings on a boundary -> ordinary Dubins path;
  2. an edge: one heading on a boundary, the other free and strictly inside
     its interval. Pontryagin's transversality condition then says the free
     end behaves like a "free heading" end, so the path is CS or CC
     (Dubins to a point), see free_heading.paths_to_point. The case with the
     *end* heading fixed is handled by reversing the path;
  3. the interior: both headings free. Candidates are S, single arcs C, and
     CC paths that switch at the midpoint of AB, all in closed form, see
     free_heading.paths_between_points.

This is the candidate list of Manyam et al. for the two-point interval problem
(normal extremals, lambda_0 = 1). We keep every C and CC path, not only those
with arcs longer than pi: the extra ones are real paths, so they cannot push
the minimum too low, and the pi-length arcs are where the abnormal extremals
of this problem live.

Every candidate is a real, flyable path, so the result can never be *below*
the true minimum by more than floating-point error. `brute_force_interval`
gives an independent check used by the tests.

Speed: `solve_interval_problem` runs the C++ port in _interval.cpp (~3 us per
solve; build it with `python tools/build_interval.py`). This Python version,
`solve_interval_problem_py`, is the reference: the tests check that the two
agree. If the C++ module is not built, the Python version is used (~0.5 ms
per solve) and a warning is printed.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

from ..utils.angles import HeadingInterval, mod2pi
from .free_heading import paths_between_points, paths_to_point
from .path import shortest_length

try:
    from . import _interval as _ext
except ImportError:
    _ext = None

Point = tuple[float, float]

# Order of the case codes returned by the C++ solver.
CASES = ("corner", "edge(start fixed)", "edge(end fixed)", "interior(S/C/CC)")


@dataclass(frozen=True)
class IntervalSolution:
    length: float
    depart_heading: float  # heading at A, inside I0
    arrive_heading: float  # heading at B, inside I1
    case: str  # which candidate won, for debugging


def _better(best: IntervalSolution | None, length, th0, th1, case) -> IntervalSolution | None:
    if not np.isfinite(length):
        return best
    if best is None or length < best.length:
        return IntervalSolution(float(length), float(mod2pi(th0)), float(mod2pi(th1)), case)
    return best


def _corners(a, i0, b, i1, r, best):
    for th0 in i0.boundaries:
        for th1 in i1.boundaries:
            length = shortest_length((a[0], a[1], th0), (b[0], b[1], th1), r)
            best = _better(best, length, th0, th1, "corner")
    return best


def _start_fixed_edges(a, i0, b, i1, r, best):
    """th0 on a boundary of I0, th1 free inside I1."""
    for th0 in i0.boundaries:
        for length, th1 in paths_to_point(a, th0, b, r):
            if np.isfinite(length[0]) and i1.contains(th1[0]):
                best = _better(best, length[0], th0, th1[0], "edge(start fixed)")
    return best


def _end_fixed_edges(a, i0, b, i1, r, best):
    """th1 on a boundary of I1, th0 free inside I0.

    Reversing a path swaps its ends and flips headings by pi but keeps its
    length, so we fly backwards from (B, th1 + pi) to the point A.
    """
    for th1 in i1.boundaries:
        for length, rev_heading in paths_to_point(b, th1 + math.pi, a, r):
            th0 = rev_heading[0] + math.pi
            if np.isfinite(length[0]) and i0.contains(th0):
                best = _better(best, length[0], th0, th1, "edge(end fixed)")
    return best


def _interior_closed_form(a, i0, b, i1, r, best):
    for length, th0, th1 in paths_between_points(a, b, r):
        if i0.contains(th0) and i1.contains(th1):
            best = _better(best, length, th0, th1, "interior(S/C/CC)")
    return best


def solve_interval_problem(
    a: Point, i0: HeadingInterval, b: Point, i1: HeadingInterval, r: float
) -> IntervalSolution:
    """Shortest Dubins path from (A, any heading in I0) to (B, any heading in I1)."""
    if _ext is None:
        return solve_interval_problem_py(a, i0, b, i1, r)
    length, th0, th1, case = _ext.solve(
        a[0], a[1], i0.start, i0.width, b[0], b[1], i1.start, i1.width, r
    )
    return IntervalSolution(length, th0, th1, CASES[case])


def solve_interval_problem_py(
    a: Point, i0: HeadingInterval, b: Point, i1: HeadingInterval, r: float
) -> IntervalSolution:
    """Pure-Python reference for solve_interval_problem (same candidates, same order)."""
    best = None
    best = _corners(a, i0, b, i1, r, best)
    best = _start_fixed_edges(a, i0, b, i1, r, best)
    best = _end_fixed_edges(a, i0, b, i1, r, best)
    best = _interior_closed_form(a, i0, b, i1, r, best)
    if best is None:  # cannot happen: corners always give a path
        raise RuntimeError("no Dubins path found")
    return best


if _ext is None:
    warnings.warn(
        "C++ interval solver not built, using the pure-Python one (~150x slower). "
        "Build it with: python tools/build_interval.py",
        stacklevel=2,
    )


def brute_force_interval(
    a: Point, i0: HeadingInterval, b: Point, i1: HeadingInterval, r: float, grid: int = 120
) -> float:
    """Reference answer: minimum over a grid of heading pairs (slow, for tests)."""
    best = math.inf
    for th0 in i0.start + i0.width * np.linspace(0, 1, grid):
        for th1 in i1.start + i1.width * np.linspace(0, 1, grid):
            best = min(best, shortest_length((a[0], a[1], th0), (b[0], b[1], th1), r))
    return best
