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
    exact: np.ndarray | None = None  # False where `cost` is only a lower estimate

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
    """edge_fn(i, option_a, j, option_b) -> (cost, depart, arrive, is_exact)."""
    node_target, node_option, clusters = _layout(options_per_target)
    size = len(node_target)
    cost = np.full((size, size), np.inf)
    dep, arr = np.zeros((size, size)), np.zeros((size, size))
    exact = np.ones((size, size), dtype=bool)
    n = len(clusters)
    total, done = n * (n - 1), 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for u in clusters[i]:
                for v in clusters[j]:
                    cost[u, v], dep[u, v], arr[u, v], exact[u, v] = edge_fn(
                        i, node_option[u], j, node_option[v])
            done += 1
            if progress:
                progress(done, total)
    return ClusteredGraph(cost, dep, arr, node_target, node_option, clusters, exact)


def build_heading_graph(points: list[Point], headings_per_target: list[list[float]], r: float,
                        cache: dict | None = None) -> ClusteredGraph:
    """Graph for the upper bound."""
    cache = {} if cache is None else cache

    def edge(i, th_a, j, th_b):
        key = (i, round(th_a, 12), j, round(th_b, 12))
        if key not in cache:
            cache[key] = shortest_length((*points[i], th_a), (*points[j], th_b), r)
        return cache[key], th_a, th_b, True

    return _build(headings_per_target, edge)


FULL_CIRCLE = HeadingInterval(0.0, TWO_PI)


class IntervalEdgeCosts:
    """Interval-problem edge costs, cached, with lower estimates from ancestors.

    When the refinement splits an interval, the pieces remember their parent
    (`parent[piece] = interval`). For an edge between a child interval a and
    an interval b, every computed edge between an ancestor-or-self of a and an
    ancestor-or-self of b gives information, because a child is a subset of
    its parent:

      * cost(a, b) >= cost(ancestor pair)          -> a valid LOWER ESTIMATE;
      * if the ancestor's optimal path departs inside a and arrives inside b,
        that path is feasible for (a, b), so cost(a, b) equals it EXACTLY.

    Intervals of the first iteration have the full circle as parent; the full
    circle edge is just the straight line (Euclidean distance), computed
    without calling the interval solver.
    """

    def __init__(self, points: list[Point], r: float, sweep_samples: int = 128):
        self.points, self.r, self.sweep_samples = points, r, sweep_samples
        self.exact_cache: dict = {}
        self.parent: dict[HeadingInterval, HeadingInterval] = {}
        self.solves = 0  # number of interval problems actually solved
        self.solve_seconds = 0.0  # time spent inside the interval solver

    @staticmethod
    def _key(i, a, j, b):
        return (i, round(a.start, 12), round(a.width, 12), j, round(b.start, 12), round(b.width, 12))

    def register_split(self, parent: HeadingInterval, pieces: list[HeadingInterval]) -> None:
        for piece in pieces:
            if piece != parent:
                self.parent[piece] = parent

    def register_initial(self, intervals_per_target: list[list[HeadingInterval]]) -> None:
        """First-iteration intervals descend from the full circle."""
        for intervals in intervals_per_target:
            self.register_split(FULL_CIRCLE, [iv for iv in intervals if not iv.is_full_circle])

    def _ancestors(self, interval):
        chain = [interval]
        while chain[-1] in self.parent:
            chain.append(self.parent[chain[-1]])
        return chain

    def exact(self, i, a, j, b):
        """(cost, depart heading, arrive heading), solving the interval problem if needed."""
        key = self._key(i, a, j, b)
        if key not in self.exact_cache:
            if a.is_full_circle and b.is_full_circle:   # both headings free: straight line
                (x0, y0), (x1, y1) = self.points[i], self.points[j]
                heading = float(mod2pi(math.atan2(y1 - y0, x1 - x0)))
                self.exact_cache[key] = (math.hypot(x1 - x0, y1 - y0), heading, heading)
            else:
                t0 = time.perf_counter()
                sol = solve_interval_problem(self.points[i], a, self.points[j], b, self.r,
                                             self.sweep_samples)
                self.solve_seconds += time.perf_counter() - t0
                self.solves += 1
                self.exact_cache[key] = (sol.length, sol.depart_heading, sol.arrive_heading)
        return self.exact_cache[key]

    def estimate(self, i, a, j, b):
        """(cost, depart, arrive, is_exact) from known ancestor edges, or None."""
        key = self._key(i, a, j, b)
        if key in self.exact_cache:
            return (*self.exact_cache[key], True)
        best = None
        for anc_a in self._ancestors(a):
            for anc_b in self._ancestors(b):
                if anc_a.is_full_circle and anc_b.is_full_circle:
                    known = self.exact(i, anc_a, j, anc_b)   # straight line, no solver
                else:
                    known = self.exact_cache.get(self._key(i, anc_a, j, anc_b))
                if known is None:
                    continue
                cost, dep, arr = known
                if a.contains(dep) and b.contains(arr):     # ancestor's path fits: exact
                    self.exact_cache[key] = known
                    return (cost, dep, arr, True)
                if best is None or cost > best[0]:          # tightest lower estimate
                    best = (cost, dep, arr, False)
        return best


def build_interval_graph(points: list[Point], intervals_per_target: list[list[HeadingInterval]],
                         r: float, costs: IntervalEdgeCosts, lazy: bool = False,
                         progress: ProgressFn | None = None) -> ClusteredGraph:
    """Graph for the lower bound.

    lazy=False: every edge cost is exact.
    lazy=True:  edges are filled from `costs.estimate` where possible; graph.exact
                marks which entries are exact. The caller must make the edges of
                its final tour exact (see bounds.py).
    """
    def edge(i, iv_a, j, iv_b):
        if lazy:
            found = costs.estimate(i, iv_a, j, iv_b)
            if found is not None:
                return found
        return (*costs.exact(i, iv_a, j, iv_b), True)

    return _build(intervals_per_target, edge, progress)


def refresh_edge(graph: ClusteredGraph, costs: IntervalEdgeCosts, u: int, v: int) -> None:
    """Replace a lower estimate in the graph by the exact interval cost."""
    i, j = graph.node_target[u], graph.node_target[v]
    cost, dep, arr = costs.exact(i, graph.node_option[u], j, graph.node_option[v])
    graph.cost[u, v], graph.depart_heading[u, v], graph.arrive_heading[u, v] = cost, dep, arr
    graph.exact[u, v] = True
