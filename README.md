# Dubins TSP bounds (Python)

Computes an **upper bound** and a **lower bound** on the shortest closed tour of
a Dubins vehicle (forward-only, minimum turning radius `r`) through a set of points.
Both bounds turn the problem into a one-in-a-set (generalised) TSP and solve it
with the **Noon–Bean transformation** to an ordinary asymmetric TSP.

Python rewrite of the Julia `Dubins-TSP` code, with the issues found in that
version fixed (see the end of this file).

## Quick start

```bash
pip install -r requirements.txt
python -m src                    # uses ./config.yml
python -m src --config other.yml --no-plot
python -m pytest                 # run the checks
```

Edit `config.yml` to set the points, the turning radius, the number of
discrete headings (upper bound), the number of heading intervals (lower bound)
and the refinement settings `lambda` and `mu`.

## Output

```
results/
  iter_00/   summary.txt  results.json  tour.png    uniform grids from config.yml
  iter_01/   ...                                    after the 1st refinement
  ...
  convergence.txt    table of every iteration (bounds, best-so-far, gap, nodes, time)
  convergence.json   everything, machine-readable
  convergence.png    bounds and gap against iteration
```

## Iterative refinement (lambda, mu)

The lower-bound tour picks one heading interval at every target. Each
refinement iteration, for every target:

1. place `lambda` evenly spaced headings strictly inside that interval;
2. add them to the target's upper-bound heading set;
3. split the interval at the same headings (lower bound gets finer there);

then both bounds are solved again. Repeated `mu` times.

Both parts are needed: if only headings were added, the next lower bound would
pick the same intervals and every later iteration would repeat the first.

With exact solvers the upper bound can only go down (heading sets only grow)
and the lower bound can only go up (each new interval is inside an old one, so
its edge costs can only rise). LKH can wobble slightly, so the gap is computed
from the best upper and best lower bound seen so far. Edge costs are cached
between iterations; only edges touching new options are computed.

## What the two bounds are

| | Upper bound | Lower bound |
|---|---|---|
| node | target + one heading `k·360°/N` | target + heading interval `[k, k+1]·360°/N` |
| edge cost | point-to-point Dubins path | two-point Dubins **interval problem** |
| tour | flyable (same heading in and out) | may "kink" at targets, not flyable |
| guarantee | always ≥ optimum, any solver | ≤ optimum, **only with `solver: exact`** |

Increasing `N` tightens both bounds.

## Code layout

```
config.yml               points + settings
src/
  __main__.py            entry point: python -m src
  config.py              load / validate config.yml
  bounds.py              pipeline: graph -> Noon-Bean -> solve -> decode
  refine.py              lower-bound-driven refinement (lambda, mu)
  dubins/
    path.py              classical Dubins shortest path (6 words)
    free_heading.py      paths with free end headings (CS, CC, S, C)
    interval.py          two-point Dubins interval problem (+ brute-force reference)
  tsp/
    gtsp.py              clustered graphs for both bounds
    noon_bean.py         Noon-Bean transformation and decoding
    solvers.py           LKH (elkai) and exact MILP (pulp/CBC)
  utils/
    angles.py            angle wrapping, HeadingInterval
    plotting.py          plots of the tours
    report.py            text report, JSON export
tests/                   correctness checks (see below)
```

## Interval solver (how the lower-bound edge costs are computed)

The best heading pair `(θ0, θ1)` in the rectangle `I0 × I1` is at

1. a **corner**: ordinary Dubins path between boundary headings;
2. an **edge**: one heading on a boundary, the other free inside its interval.
   Transversality makes the free end a "free heading" end, so the path is CS
   or CC. The end-fixed case is done by flying the path backwards;
3. the **interior**: both headings free. S and single arcs C in closed form;
   the CC family is swept over `θ0` and every local minimum is refined.

Every candidate is a real path, so the solver can't undershoot the true
minimum. The tests check it never exceeds a brute-force grid search.

## Reusing previous work (`lower_bound.lazy`)

Most of the runtime is solving interval problems for the lower-bound edges.
Between refinement iterations they are reused in three ways:

1. **Cache.** Every interval problem solved is kept; unchanged edges are free.
2. **Inheritance (exact).** A split piece is a subset of its parent, so its
   edge cost is >= the parent's. If the parent's optimal path departs and
   arrives with headings inside the pieces, the piece's cost is *exactly* the
   parent's: no solve needed.
3. **Lazy evaluation.** Every other new edge starts with its parent's cost as
   a lower estimate (for the first iteration: the straight-line distance).
   Solve, make the edges around the tour exact, and repeat until the tour uses
   exact edges only. That tour is then optimal for the exact costs as well:

       estimated optimum <= exact optimum <= exact cost of the tour = estimated optimum

   so the lower bound is identical to computing every edge. Cheap LKH rounds
   pick which edges to compute; the final answer is always the exact solver's.
   The exact model keeps its subtour cuts between rounds.

On the example (6 targets, lambda = 2, 10 iterations) this cuts interval
solves from 34,320 to about 6,900 and runtime from ~121 s to ~96 s, with
identical bounds at every iteration. In later iterations the exact MILP solve
becomes the bottleneck.

## Choosing the solver

* `lkh`: fast heuristic (LKH via `elkai`). Use for the upper bound.
* `exact`: integer program with subtour cuts (CBC via `pulp`). Proves
  optimality, so it is needed for a *certified* lower bound. It becomes slow
  for large `targets × intervals`; if the time limit is hit the program
  stops instead of reporting an uncertified value.

For the lower bound, costs are rounded **down** to integers, and the reported
value is the integer optimum divided by `cost_scale`. This keeps the bound valid.

## Tests

* every Dubins word really ends at the goal configuration;
* the interval solver is never above brute force, and its path really exists;
* Noon–Bean + each solver reproduces the brute-force one-in-a-set optimum;
* with one 360° interval the lower bound equals the Euclidean TSP;
* the upper-bound tour has no kinks;
* a refinement step adds headings inside the chosen interval and keeps the
  interval partition covering the circle exactly once;
* with exact solvers the lower bound never drops and the upper bound never rises;
* lazy evaluation gives exactly the same lower bounds as computing every edge;
* a split piece containing the parent's optimal headings inherits its cost exactly.

## Changes from the Julia version

* Successor shift `(V1 + 1) % N` produced heading interval `[-360°/N, 0]` for
  one option; the interval check then rejected most candidates, making some
  lower-bound edge costs too high. Intervals now wrap properly.
* The 3-segment search never updated its running minimum and used a global
  variable that leaked between calls.
* The turning radius was hard-coded to 1 in several places.
* A lower bound computed with LKH was reported as a bound; it is now flagged
  as an estimate unless `solver: exact` is used.
* Big M is strictly larger than any tour cost (`n · max + 1`).
* Decoding checks the tour really visits each cluster as one block.
* Dubins formulas are guarded against rounding noise (a zero-length turn
  coming out as a full circle, a zero-length straight being rejected).
* Runs on Linux / macOS / Windows (no bundled `.exe`).
