"""Checks of the lower-bound-driven refinement."""
import math

from src import refine
from src.utils.angles import TWO_PI, HeadingInterval
from tests.test_bounds import POINTS, _config
from src.config import BoundConfig


def test_refine_step_adds_headings_inside_selected_interval():
    headings = [[0.0, math.pi / 2, math.pi, 3 * math.pi / 2]]
    intervals = [[HeadingInterval(k * math.pi / 2, math.pi / 2) for k in range(4)]]
    chosen = intervals[0][1]                                  # [90°, 180°]
    new_h, new_iv, added = refine.refine(headings, intervals, [chosen], lam=2)

    assert [round(math.degrees(a), 6) for a in added[0]] == [120.0, 150.0]
    assert len(new_h[0]) == 6
    # The selected interval is split into 3 equal pieces; the partition still
    # covers the full circle exactly once.
    assert len(new_iv[0]) == 6
    assert math.isclose(sum(iv.width for iv in new_iv[0]), TWO_PI)
    pieces = [iv for iv in new_iv[0] if math.isclose(iv.width, math.pi / 6)]
    assert [round(math.degrees(iv.start), 6) for iv in pieces] == [90.0, 120.0, 150.0]


def test_bounds_move_the_right_way():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    its = refine.run(cfg, lam=1, mu=3)
    lowers = [it.lower.value for it in its]
    uppers = [it.upper.value for it in its]
    tol = len(POINTS) / cfg.cost_scale  # integer rounding of the costs
    # Exact solvers: refined intervals only raise the LB, extra headings only lower the UB.
    assert all(b >= a - tol for a, b in zip(lowers, lowers[1:]))
    assert all(b <= a + tol for a, b in zip(uppers, uppers[1:]))
    assert all(it.best_lower <= it.best_upper for it in its)
    assert its[-1].gap_percent <= its[0].gap_percent


def test_lazy_lower_bound_equals_eager():
    # Lazy evaluation reuses parent intervals' costs as lower estimates and
    # only computes edges the solver wants; the bound must not change.
    def lower_bounds(lazy):
        cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh", lazy))
        return [it.lower.value for it in refine.run(cfg, lam=2, mu=3)]

    assert lower_bounds(lazy=True) == lower_bounds(lazy=False)


def test_child_interval_inherits_exact_cost():
    from src.tsp.gtsp import IntervalEdgeCosts

    costs = IntervalEdgeCosts(POINTS, 1.0)
    parent = HeadingInterval(0.0, math.pi / 2)
    other = HeadingInterval(math.pi, math.pi / 2)
    cost, dep, _ = costs.exact(0, parent, 1, other)
    pieces = [HeadingInterval(k * math.pi / 6, math.pi / 6) for k in range(3)]
    costs.register_split(parent, pieces)
    solves = costs.solves
    holder = next(p for p in pieces if p.contains(dep))
    found = costs.estimate(0, holder, 1, other)
    assert found[3] and found[0] == cost and costs.solves == solves   # exact, no new solve
    for p in pieces:
        if p != holder and not p.contains(dep):
            est = costs.estimate(0, p, 1, other)
            assert not est[3] and est[0] <= costs.exact(0, p, 1, other)[0] + 1e-12
