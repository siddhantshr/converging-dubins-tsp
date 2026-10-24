"""Noon-Bean transformation: one-in-a-set ATSP  ->  ordinary ATSP.

(Noon & Bean 1993; explained in Oberlin, Rathinam & Darbha, "Today's
Traveling Salesman Problem", IEEE RAM 2010.)

An ATSP solver must visit *every* node. We want *one* node per cluster.
Three tricks make the two problems equivalent:

1. Zero-cost ring: inside each cluster add edges
       u -> successor(u)  with cost 0
   so the solver can sweep through the unused options for free.
   (These edges are bookkeeping only; the vehicle never flies them.)

2. Successor shift: a tour that enters a cluster at node e goes round the
   ring and leaves from the node just before e. So every edge leaving the
   cluster is charged as if it left from the successor:
       C_hat(u -> w) = C(successor(u) -> w)
   With this, the cluster costs exactly as if only e had been visited.

3. Big M: add M to every edge between different clusters. An honest tour
   uses exactly n such edges; a tour that leaves a cluster's ring early needs
   at least n + 1. With M larger than any real tour cost, honest tours are
   always cheaper.

Decoding: the first node visited in each cluster is the chosen option, and
    GTSP cost = ATSP cost - n * M.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class DecodeError(RuntimeError):
    """The ATSP tour does not visit each cluster as one consecutive block."""


@dataclass
class NoonBeanATSP:
    matrix: np.ndarray      # int64 ATSP costs
    big_m: int
    forbidden: int          # cost of edges that must never be used
    successor: np.ndarray   # ring successor of every node
    cluster_of: np.ndarray  # cluster index of every node
    num_clusters: int

    def tour_cost(self, tour: list[int]) -> int:
        return int(sum(self.matrix[tour[i], tour[(i + 1) % len(tour)]] for i in range(len(tour))))

    def gtsp_cost(self, tour: list[int]) -> int:
        """Cost of the corresponding one-in-a-set tour (removes the n*M offset)."""
        return self.tour_cost(tour) - self.num_clusters * self.big_m


def transform(cost: np.ndarray, clusters: list[list[int]]) -> NoonBeanATSP:
    """Build the ATSP matrix from integer GTSP costs.

    `cost[u, v]` is only read for u and v in different clusters.
    """
    size = cost.shape[0]
    n = len(clusters)

    cluster_of = np.empty(size, dtype=int)
    successor = np.empty(size, dtype=int)
    for c, members in enumerate(clusters):
        for idx, u in enumerate(members):
            cluster_of[u] = c
            successor[u] = members[(idx + 1) % len(members)]

    different = cluster_of[:, None] != cluster_of[None, :]
    max_cost = int(cost[different].max())
    big_m = n * max_cost + 1          # strictly larger than any GTSP tour cost
    forbidden = (n + 1) * big_m       # larger than any honest ATSP tour

    if forbidden >= 2**31 - 1:
        raise OverflowError(
            "ATSP costs too large for the solver; lower `solver.cost_scale` in config.yml")

    # Trick 2 + 3: rows of the successor, plus M, for edges between clusters.
    matrix = np.where(different, cost[successor, :] + big_m, forbidden).astype(np.int64)
    # Trick 1: free ring inside each cluster.
    for u in range(size):
        if successor[u] != u:
            matrix[u, successor[u]] = 0
    np.fill_diagonal(matrix, forbidden)
    return NoonBeanATSP(matrix, big_m, forbidden, successor, cluster_of, n)


def decode(tour: list[int], atsp: NoonBeanATSP) -> list[int]:
    """Chosen node of each cluster, in visiting order.

    The chosen node is the one where the tour *enters* the cluster.
    """
    entries = [tour[i] for i in range(len(tour))
               if atsp.cluster_of[tour[i]] != atsp.cluster_of[tour[i - 1]]]
    if len(entries) != atsp.num_clusters:
        raise DecodeError(
            f"tour enters clusters {len(entries)} times, expected {atsp.num_clusters}; "
            "the solver broke a zero-cost ring")
    return entries
