"""The interval solver must never be above brute force (needed for a valid
lower bound) and must always return a path that really exists."""

import math
import random

import pytest

from src.dubins import interval
from src.dubins.interval import (
    brute_force_interval,
    solve_interval_problem,
    solve_interval_problem_py,
)
from src.dubins.path import shortest_length
from src.utils.angles import HeadingInterval


def _random_case(rng):
    r = rng.uniform(0.5, 2.0)
    dist = r * rng.choice(
        [rng.uniform(0.05, 1), rng.uniform(1, 2), rng.uniform(2, 4), rng.uniform(4, 8)]
    )
    ang = rng.uniform(0, 2 * math.pi)
    a = (rng.uniform(-3, 3), rng.uniform(-3, 3))
    b = (a[0] + dist * math.cos(ang), a[1] + dist * math.sin(ang))
    i0 = HeadingInterval(rng.uniform(0, 2 * math.pi), rng.uniform(0.05, 2 * math.pi))
    i1 = HeadingInterval(rng.uniform(0, 2 * math.pi), rng.uniform(0.05, 2 * math.pi))
    return a, i0, b, i1, r


@pytest.mark.parametrize("seed", range(4))
def test_against_brute_force(seed):
    rng = random.Random(seed)
    for _ in range(15):
        a, i0, b, i1, r = _random_case(rng)
        sol = solve_interval_problem(a, i0, b, i1, r)
        assert sol.length <= brute_force_interval(a, i0, b, i1, r, grid=50) + 1e-6
        # The returned headings are inside their intervals ...
        assert i0.contains(sol.depart_heading) and i1.contains(sol.arrive_heading)
        # ... and a real path of that length exists between them.
        real = shortest_length((*a, sol.depart_heading), (*b, sol.arrive_heading), r)
        assert real <= sol.length + 1e-6


def test_interval_wrapping_past_2pi():
    # [315 deg, 360 deg] given as start=5.5 with the wrap; the Julia version
    # mishandled this case (heading interval [-45 deg, 0]).
    a, b, r = (0.0, 0.0), (0.5, -1.2), 1.0
    i0 = HeadingInterval(7 * math.pi / 4, math.pi / 4)
    i0_neg = HeadingInterval(-math.pi / 4, math.pi / 4)  # same set of headings
    i1 = HeadingInterval(math.pi, math.pi / 2)
    s1 = solve_interval_problem(a, i0, b, i1, r)
    s2 = solve_interval_problem(a, i0_neg, b, i1, r)
    assert math.isclose(s1.length, s2.length, abs_tol=1e-9)
    assert s1.length <= brute_force_interval(a, i0, b, i1, r, grid=80) + 1e-6


def test_full_circle_interval_is_straight_line():
    full = HeadingInterval(0.0, 2 * math.pi)
    sol = solve_interval_problem((0, 0), full, (3, 4), full, 1.0)
    assert math.isclose(sol.length, 5.0, abs_tol=1e-9)


@pytest.mark.skipif(interval._ext is None, reason="C++ solver not built (tools/build_interval.py)")
@pytest.mark.parametrize("seed", range(4))
def test_cpp_matches_python(seed):
    """The C++ port gives the Python reference's length. Headings may differ
    only on ties, so check those give a real path of that length instead."""
    rng = random.Random(100 + seed)
    for _ in range(500):
        a, i0, b, i1, r = _random_case(rng)
        cpp = solve_interval_problem(a, i0, b, i1, r)
        ref = solve_interval_problem_py(a, i0, b, i1, r)
        assert math.isclose(cpp.length, ref.length, rel_tol=1e-12, abs_tol=1e-12)
        assert i0.contains(cpp.depart_heading) and i1.contains(cpp.arrive_heading)
        real = shortest_length((*a, cpp.depart_heading), (*b, cpp.arrive_heading), r)
        assert real <= cpp.length + 1e-9


@pytest.mark.skipif(interval._ext is None, reason="C++ solver not built (tools/build_interval.py)")
def test_cpp_shortest_length_matches_python():
    """Upper-bound edges use the C++ port of path.shortest_length."""
    rng = random.Random(7)
    for _ in range(2000):
        r = rng.uniform(0.5, 2.0)
        q0 = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-7, 14))
        q1 = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-7, 14))
        assert math.isclose(
            interval.dubins_length(q0, q1, r),
            shortest_length(q0, q1, r),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
