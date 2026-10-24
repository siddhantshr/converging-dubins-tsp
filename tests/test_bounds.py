"""End-to-end checks of the two bounds."""
import dataclasses
import itertools
import math

from src.bounds import compute_lower_bound, compute_upper_bound
from src.config import BoundConfig, Config

POINTS = [(0.0, 0.0), (4.0, 1.0), (5.0, 5.0), (1.0, 4.0), (2.0, 2.0)]


def _config(**kw):
    base = Config(points=POINTS, turning_radius=1.0, heading_offset=0.0,
                  upper=BoundConfig(True, 6, "lkh"), lower=BoundConfig(True, 6, "lkh"),
                  sweep_samples=64, lkh_runs=5, cost_scale=1000,
                  refine_lambda=1, refine_mu=0, plot=False, results_dir="")
    return dataclasses.replace(base, **kw)


def test_lower_below_upper():
    cfg = _config()
    assert compute_lower_bound(cfg).value <= compute_upper_bound(cfg).value


def test_one_interval_gives_euclidean_tsp():
    # With a single 360-degree interval every heading is free, each leg is a
    # straight line, so the lower bound is the Euclidean TSP length.
    cfg = _config(lower=BoundConfig(True, 1, "lkh"))
    best = min(
        sum(math.dist(POINTS[o[i]], POINTS[o[(i + 1) % len(o)]]) for i in range(len(o)))
        for o in ((0,) + p for p in itertools.permutations(range(1, len(POINTS)))))
    lb = compute_lower_bound(cfg).value
    assert lb <= best + 1e-9 and lb >= best - len(POINTS) / cfg.cost_scale


def test_upper_bound_tour_is_flyable():
    res = compute_upper_bound(_config())
    for leg, nxt in zip(res.legs, res.legs[1:] + res.legs[:1]):
        assert math.isclose(leg.arrive_heading, nxt.depart_heading)  # no kink
