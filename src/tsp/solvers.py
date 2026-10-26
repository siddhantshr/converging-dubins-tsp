"""Solvers for the (transformed) ATSP.

* "lkh":   Lin-Kernighan-Helsgaun heuristic, via the `elkai` package (LKH is
           bundled, no binary to install). Fast and usually optimal, but gives
           no guarantee. Fine for the upper bound: any tour is a valid upper
           bound.
* "exact": integer program solved with CBC (via `pulp`). Proves optimality.
           No longer used by the bounds (too slow); kept only as the
           reference solver in tests/test_noon_bean.py.

Both return a tour as a list of node ids (each node once, no repeat at the end).
"""

from __future__ import annotations

import time

import numpy as np

from .noon_bean import NoonBeanATSP


def solve_lkh(atsp: NoonBeanATSP, runs: int = 10) -> list[int]:
    import elkai

    tour = elkai.DistanceMatrix(atsp.matrix.tolist()).solve_tsp(runs=runs)
    return [int(v) for v in tour[:-1]]  # elkai repeats the start node at the end


def _cycles(successor: dict[int, int]) -> list[list[int]]:
    seen, cycles = set(), []
    for start in successor:
        if start in seen:
            continue
        cycle, node = [], start
        while node not in seen:
            seen.add(node)
            cycle.append(node)
            node = successor[node]
        cycles.append(cycle)
    return cycles


class ExactATSP:
    """Exact ATSP by integer programming with lazily added subtour cuts.

    Start from the assignment problem (every node has one edge in and one
    out). Whenever the solution splits into several cycles, add the
    Dantzig-Fulkerson-Johnson cut  sum_{u,v in S} x_uv <= |S| - 1  for each
    cycle S and solve again. The ring of every cluster is cut from the start,
    because those zero-cost rings are the first subtours the solver would find.

    The model can be solved again after the costs change: the constraints
    only depend on which edges exist,
    so all cuts found so far stay valid and are kept, and the previous tour
    is given to CBC as a starting solution.
    """

    def __init__(self, atsp: NoonBeanATSP):
        import pulp

        self._pulp = pulp
        size = atsp.matrix.shape[0]
        self.arcs = [
            (u, v)
            for u in range(size)
            for v in range(size)
            if u != v and atsp.matrix[u, v] < atsp.forbidden
        ]
        self.prob = pulp.LpProblem("noon_bean_atsp", pulp.LpMinimize)
        self.x = {(u, v): pulp.LpVariable(f"x_{u}_{v}", cat="Binary") for u, v in self.arcs}
        self.last_tour: list[int] | None = None

        out_arcs = {u: [] for u in range(size)}
        in_arcs = {v: [] for v in range(size)}
        for (u, v), var in self.x.items():
            out_arcs[u].append(var)
            in_arcs[v].append(var)
        for node in range(size):
            self.prob += pulp.lpSum(out_arcs[node]) == 1
            self.prob += pulp.lpSum(in_arcs[node]) == 1

        for c in range(atsp.num_clusters):
            members = list(np.flatnonzero(atsp.cluster_of == c))
            if 1 < len(members) < size:
                self._add_subtour_cut(members)

    def _add_subtour_cut(self, nodes):
        members = set(nodes)
        self.prob.addConstraint(
            self._pulp.lpSum(var for (u, v), var in self.x.items() if u in members and v in members)
            <= len(members) - 1
        )

    def solve(self, atsp: NoonBeanATSP, time_limit: float = 120.0) -> list[int]:
        pulp = self._pulp
        self.prob.setObjective(
            pulp.lpSum(int(atsp.matrix[u, v]) * var for (u, v), var in self.x.items())
        )
        deadline = time.monotonic() + time_limit
        while True:
            warm = self.last_tour is not None
            if warm:
                used = {
                    (self.last_tour[i], self.last_tour[(i + 1) % len(self.last_tour)])
                    for i in range(len(self.last_tour))
                }
                for arc, var in self.x.items():
                    var.setInitialValue(1 if arc in used else 0)
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("exact ATSP solve hit the time limit")
            status = self.prob.solve(
                pulp.PULP_CBC_CMD(msg=False, timeLimit=max(1, int(remaining)), warmStart=warm)
            )
            # sol_status distinguishes "proved optimal" from "stopped with some
            # integer solution" (e.g. on a time limit); only the first is a proof.
            if self.prob.sol_status != pulp.LpSolutionOptimal:
                raise TimeoutError(
                    f"exact ATSP solve did not prove optimality ({pulp.LpStatus[status]})"
                )

            successor = {u: v for (u, v), var in self.x.items() if var.value() > 0.5}
            cycles = _cycles(successor)
            if len(cycles) == 1:
                self.last_tour = _cycle_from(successor, start=0)
                return self.last_tour
            for cycle in cycles:
                self._add_subtour_cut(cycle)


def solve_exact(atsp: NoonBeanATSP, time_limit: float = 120.0) -> list[int]:
    return ExactATSP(atsp).solve(atsp, time_limit)


def _cycle_from(successor: dict[int, int], start: int) -> list[int]:
    tour, node = [start], successor[start]
    while node != start:
        tour.append(node)
        node = successor[node]
    return tour


def solve(atsp: NoonBeanATSP, method: str, lkh_runs: int, time_limit: float) -> list[int]:
    if method == "lkh":
        return solve_lkh(atsp, lkh_runs)
    if method == "exact":
        return solve_exact(atsp, time_limit)
    raise ValueError(f"unknown solver {method!r} (use 'lkh' or 'exact')")
