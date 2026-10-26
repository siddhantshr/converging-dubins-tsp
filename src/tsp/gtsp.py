"""Clustered graphs for the one-in-a-set (generalised) TSP.

Each target becomes a cluster of nodes, one per heading option. Targets may
have different numbers of options (needed by the iterative refinement).

* Upper bound: an option is a single heading.
  Edge cost = ordinary point-to-point Dubins length.
* Lower bound: an option is a heading interval.
  Edge cost = two-point Dubins interval problem.

Edges between two nodes of the same target are not used (cost = inf).

Edge costs are expensive (especially interval ones), so both builders accept a
cache dict that survives between calls: when the refinement adds a few options,
only the new edges are computed.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..dubins.interval import solve_interval_problem
from ..dubins.path import shortest_length
from ..utils.angles import TWO_PI, HeadingInterval, mod2pi

Point = tuple[float, float]
ProgressFn = Callable[[int, int], None]


@dataclass
class ClusteredGraph:
    cost: np.ndarray            # (N, N) float, inf inside a cluster
    depart_heading: np.ndarray  # heading used when leaving u towards v
    arrive_heading: np.ndarray  # heading used when arriving at v from u
    node_target: list[int]      # target of every node
    node_option: list           # heading (float) or HeadingInterval of every node
    clusters: list[list[int]]   # node ids of every target

    @property
    def num_targets(self) -> int:
        return len(self.clusters)

    def label(self, node: int) -> str:
        option = self.node_option[node]
        if isinstance(option, HeadingInterval):
            return str(option)
        return f"{math.degrees(option) % 360:.1f}°"


# --- option grids -------------------------------------------------------------

def uniform_headings(count: int, offset: float = 0.0) -> list[float]:
    return [float(mod2pi(offset + k * TWO_PI / count)) for k in range(count)]


def uniform_intervals(count: int, offset: float = 0.0) -> list[HeadingInterval]:
    width = TWO_PI / count
    return [HeadingInterval(float(mod2pi(offset + k * width)), width) for k in range(count)]


# --- graph construction ---------------------------------------------------------

def _layout(options_per_target):
    node_target, node_option, clusters = [], [], []
    for target, options in enumerate(options_per_target):
        start = len(node_target)
        clusters.append(list(range(start, start + len(options))))
        node_target.extend([target] * len(options))
        node_option.extend(options)
    return node_target, node_option, clusters


def _build(options_per_target, edge_fn, progress=None) -> ClusteredGraph:
    """edge_fn(i, option_a, j, option_b) -> (cost, depart, arrive)."""
    node_target, node_option, clusters = _layout(options_per_target)
    size = len(node_target)
    cost = np.full((size, size), np.inf)
    dep, arr = np.zeros((size, size)), np.zeros((size, size))
    n = len(clusters)
    total, done = n * (n - 1), 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for u in clusters[i]:
                for v in clusters[j]:
                    cost[u, v], dep[u, v], arr[u, v] = edge_fn(
                        i, node_option[u], j, node_option[v])
            done += 1
            if progress:
                progress(done, total)
    return ClusteredGraph(cost, dep, arr, node_target, node_option, clusters)


def build_heading_graph(points: list[Point], headings_per_target: list[list[float]], r: float,
                        cache: dict | None = None) -> ClusteredGraph:
    """Graph for the upper bound."""
    cache = {} if cache is None else cache

    def edge(i, th_a, j, th_b):
        key = (i, round(th_a, 12), j, round(th_b, 12))
        if key not in cache:
            cache[key] = shortest_length((*points[i], th_a), (*points[j], th_b), r)
        return cache[key], th_a, th_b

    return _build(headings_per_target, edge)


class IntervalEdgeCosts:
    """Interval-problem edge costs, cached between refinement iterations.

    An edge depends only on its two targets and two intervals, so when the
    refinement splits a few intervals, every other edge is a cache hit.
    """

    def __init__(self, points: list[Point], r: float):
        self.points, self.r = points, r
        self.cache: dict = {}
        self.solves = 0  # number of interval problems actually solved
        self.solve_seconds = 0.0  # time spent inside the interval solver

    @staticmethod
    def _key(i, a, j, b):
        return (i, round(a.start, 12), round(a.width, 12), j, round(b.start, 12), round(b.width, 12))

    def edge(self, i, a, j, b):
        """(cost, depart heading, arrive heading), solving the interval problem if needed."""
        key = self._key(i, a, j, b)
        if key not in self.cache:
            if a.is_full_circle and b.is_full_circle:   # both headings free: straight line
                (x0, y0), (x1, y1) = self.points[i], self.points[j]
                heading = float(mod2pi(math.atan2(y1 - y0, x1 - x0)))
                self.cache[key] = (math.hypot(x1 - x0, y1 - y0), heading, heading)
            else:
                t0 = time.perf_counter()
                sol = solve_interval_problem(self.points[i], a, self.points[j], b, self.r)
                self.solve_seconds += time.perf_counter() - t0
                self.solves += 1
                self.cache[key] = (sol.length, sol.depart_heading, sol.arrive_heading)
        return self.cache[key]


def build_interval_graph(points: list[Point], intervals_per_target: list[list[HeadingInterval]],
                         r: float, costs: IntervalEdgeCosts,
                         progress: ProgressFn | None = None) -> ClusteredGraph:
    """Graph for the lower bound: every edge is an interval-problem minimum."""
    return _build(intervals_per_target, costs.edge, progress)
