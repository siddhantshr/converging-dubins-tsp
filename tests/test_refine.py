"""Checks of the lower-bound-driven refinement."""

import dataclasses
import math

from src import refine
from src.config import BoundConfig
from src.utils.angles import TWO_PI, HeadingInterval
from src.utils.report import format_convergence, format_iteration
from tests.test_bounds import POINTS, _config


def test_refine_step_adds_headings_inside_selected_interval():
    headings = [[0.0, math.pi / 2, math.pi, 3 * math.pi / 2]]
    intervals = [[HeadingInterval(k * math.pi / 2, math.pi / 2) for k in range(4)]]
    chosen = intervals[0][1]  # [90°, 180°]
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
    its = refine.run(cfg, lam=1, gap_target=1e-9, max_iterations=4)
    for k, it in enumerate(its):
        assert it.lower.status == "estimate (LKH)" and not it.lower.is_valid_bound
        assert it.best_upper == min(i.upper.value for i in its[: k + 1])
        assert it.max_lb_estimate == max(i.lower.value for i in its[: k + 1])
        # gap uses this iteration's LB estimate, not the max over iterations
        assert it.gap_percent == (it.best_upper - it.lower.value) / it.lower.value * 100
        assert it.lower.value <= it.best_upper and not it.lb_exceeds_ub


def test_warns_when_lb_estimate_exceeds_best_upper():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    it = refine.run(cfg, lam=1, gap_target=None)[0]
    too_high = dataclasses.replace(it, lower=dataclasses.replace(it.lower, value=it.best_upper + 1))
    assert too_high.lb_exceeds_ub
    assert "WARNING" in format_iteration(too_high) and "WARNING" in format_convergence([too_high])


def test_lkh_is_deterministic():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    runs = [
        [
            (it.upper.value, it.lower.value)
            for it in refine.run(cfg, lam=1, gap_target=1e-9, max_iterations=3)
        ]
        for _ in range(2)
    ]
    assert runs[0] == runs[1]


def test_interval_costs_are_cached_between_iterations():
    from src.tsp.gtsp import IntervalEdgeCosts

    costs = IntervalEdgeCosts(POINTS, 1.0)
    a, b = HeadingInterval(0.0, math.pi / 2), HeadingInterval(math.pi, math.pi / 2)
    first = costs.edge(0, a, 1, b)
    assert costs.edge(0, HeadingInterval(0.0, math.pi / 2), 1, b) == first and costs.solves == 1


def test_stops_at_first_iteration_within_gap_target():
    cfg = _config(upper=BoundConfig(True, 4, "lkh"), lower=BoundConfig(True, 4, "lkh"))
    gaps = [it.gap_percent for it in refine.run(cfg, lam=1, gap_target=1e-9, max_iterations=6)]
    target = sorted(gaps)[len(gaps) // 2]  # reached somewhere in the middle
    its = refine.run(cfg, lam=1, gap_target=target)
    assert its[-1].gap_percent <= target
    assert all(it.gap_percent > target for it in its[:-1])
    assert len(refine.run(cfg, lam=1, gap_target=1000.0)) == 1  # met at iteration 0


def test_config_rejects_removed_mu(tmp_path):
    from src.config import ConfigError, load_config

    path = tmp_path / "c.yml"
    path.write_text("points: [[0, 0], [3, 0], [0, 3]]\nrefinement: {mu: 3}\n")
    try:
        load_config(path)
    except ConfigError as err:
        assert "gap_target" in str(err)
    else:
        raise AssertionError("mu should be rejected")
