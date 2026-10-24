"""Noon-Bean + ATSP solver must reproduce the brute-force one-in-a-set optimum."""
import itertools

import numpy as np
import pytest

from src.tsp import noon_bean, solvers


def brute_force_gtsp(cost, clusters):
    n, best = len(clusters), None
    for perm in itertools.permutations(range(1, n)):
        order = (0,) + perm
        for choice in itertools.product(*[clusters[c] for c in order]):
            total = sum(cost[choice[i], choice[(i + 1) % n]] for i in range(n))
            best = total if best is None or total < best else best
    return best


@pytest.mark.parametrize("method", ["lkh", "exact"])
def test_matches_brute_force(method):
    rng = np.random.default_rng(3)
    for _ in range(10):
        n, k = int(rng.integers(3, 6)), int(rng.integers(1, 5))
        cost = rng.integers(1, 100, (n * k, n * k))
        clusters = [list(range(c * k, (c + 1) * k)) for c in range(n)]
        atsp = noon_bean.transform(cost, clusters)
        tour = solvers.solve(atsp, method, lkh_runs=10, time_limit=60)
        chosen = noon_bean.decode(tour, atsp)
        real = sum(cost[chosen[i], chosen[(i + 1) % n]] for i in range(n))
        assert real == atsp.gtsp_cost(tour)          # successor shift is exact
        assert real == brute_force_gtsp(cost, clusters)


def test_decode_rejects_broken_ring():
    clusters = [[0, 1], [2, 3], [4, 5]]
    atsp = noon_bean.transform(np.ones((6, 6), dtype=int), clusters)
    with pytest.raises(noon_bean.DecodeError):
        noon_bean.decode([0, 2, 1, 3, 4, 5], atsp)  # leaves cluster 0 twice
