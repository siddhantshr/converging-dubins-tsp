"""Checks of the lower-bound-driven refinement."""
import dataclasses
import math

from src import refine
from src.utils.angles import TWO_PI, HeadingInterval
from tests.test_bounds import POINTS, _config
from src.bounds import compute_lower_bound
from src.config import BoundConfig
from src.utils.report import format_convergence, format_iteration


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


def test_reported_values_on_lkh_run():
    # LKH gives no monotonicity guarantee, so only check what always holds.
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    its = refine.run(cfg, lam=1, mu=3)
    for k, it in enumerate(its):
        assert it.lower.status == "estimate (LKH)" and not it.lower.is_valid_bound
        assert it.best_upper == min(i.upper.value for i in its[:k + 1])
        assert it.max_lb_estimate == max(i.lower.value for i in its[:k + 1])
        # gap uses this iteration's LB estimate, not the max over iterations
        assert it.gap_percent == (it.best_upper - it.lower.value) / it.lower.value * 100
        assert it.lower.value <= it.best_upper and not it.lb_exceeds_ub


def test_warns_when_lb_estimate_exceeds_best_upper():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    it = refine.run(cfg, lam=1, mu=0)[0]
    too_high = dataclasses.replace(it, lower=dataclasses.replace(it.lower, value=it.best_upper + 1))
    assert too_high.lb_exceeds_ub
    assert "WARNING" in format_iteration(too_high) and "WARNING" in format_convergence([too_high])


def test_lazy_lower_bound_equals_eager():
    # Lazy evaluation reuses parent intervals' costs as lower estimates and
    # only computes edges the solver wants. Both runs use LKH; elkai has no
    # seed argument but LKH's default SEED = 1 makes it deterministic. On this
    # small instance LKH finds the GTSP optimum, so the values must agree.
    def lower_bounds(lazy):
        cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh", lazy))
        return [it.lower.value for it in refine.run(cfg, lam=2, mu=3)]

    assert lower_bounds(lazy=True) == lower_bounds(lazy=False)


def test_lkh_is_deterministic():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    runs = [[(it.upper.value, it.lower.value) for it in refine.run(cfg, lam=1, mu=2)]
            for _ in range(2)]
    assert runs[0] == runs[1]


def test_lazy_rounds_end_with_a_full_lkh_solve():
    cfg = _config(lower=BoundConfig(True, 4, "lkh", True), lkh_runs=5, lkh_runs_screening=1)
    lower = compute_lower_bound(cfg)
    assert lower.solve_rounds >= 2                     # >= 1 screening round + the full solve
    assert lower.timings["lkh_screen"] > 0 and lower.timings["lkh"] > 0


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
