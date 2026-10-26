# Dubins TSP bounds (Python)

Computes an **upper bound** and a **lower bound** on the shortest closed tour of
a Dubins vehicle (forward-only, minimum turning radius `r`) through a set of points.
Both bounds turn the problem into a one-in-a-set (generalised) TSP and solve it
with the **Noon–Bean transformation** to an ordinary asymmetric TSP, solved by LKH.

**The lower bound is an estimate.** LKH is a heuristic: when it misses the
optimal GTSP tour, the reported lower bound is too high. It is labelled
"estimate (LKH)" everywhere. (An exact MILP solver used to certify it, but
it dominated the runtime; it was removed.)

Python rewrite of the Julia `Dubins-TSP` code, with the issues found in that
version fixed (see the end of this file).

## Quick start

```bash
pip install -r requirements.txt
python tools/build_interval.py   # compile the C++ interval solver (needs a C++17 compiler)
python -m src                    # uses ./config.yml
python -m src --config other.yml --no-plot
python -m pytest                 # run the checks
python tools/lb_stability.py 12  # LB estimate of iteration 12 with 1 / 10 / 50 LKH restarts
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
  convergence.txt    table of every iteration (bounds, best UB, max LB estimate, gap, nodes, time)
  summary.txt        time per iteration, split by stage (interval problems, LKH, ...)
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

With optimal GTSP solutions the upper bound could only go down (heading sets
only grow) and the lower bound only go up (each new interval is inside an old
one, so its edge costs can only rise). LKH can break both. Every upper bound
is a real tour, so the best one so far is kept. The lower bound is only an
estimate, so the gap is

    gap = (best upper bound - this iteration's LB estimate) / this iteration's LB estimate

The maximum LB estimate over iterations is shown ("max LB est.") but not used:
an estimate that is too high would then stick. If an LB estimate ever exceeds
the best upper bound, a warning is printed (LKH missed the optimum).
Edge costs are cached between iterations; only edges touching new options are computed.

## What the two bounds are

| | Upper bound | Lower bound |
|---|---|---|
| node | target + one heading `k·360°/N` | target + heading interval `[k, k+1]·360°/N` |
| edge cost | point-to-point Dubins path | two-point Dubins **interval problem** |
| tour | flyable (same heading in and out) | may "kink" at targets, not flyable |
| guarantee | always ≥ optimum | ≤ optimum only if the GTSP is solved optimally: **estimate** with LKH |

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
    interval.py          two-point Dubins interval problem (Python reference + brute force)
    _interval.cpp        C++ port of interval.py, used by the pipeline
  tsp/
    gtsp.py              clustered graphs for both bounds
    noon_bean.py         Noon-Bean transformation and decoding
    solvers.py           LKH (elkai); exact MILP (pulp/CBC) kept only as a test reference
  utils/
    angles.py            angle wrapping, HeadingInterval
    plotting.py          plots of the tours
    report.py            text report, JSON export
tools/
  lb_stability.py        LB estimate vs number of LKH restarts for one iteration
tests/                   correctness checks (see below)
```

## Interval solver (how the lower-bound edge costs are computed)

The best heading pair `(θ0, θ1)` in the rectangle `I0 × I1` is at

1. a **corner**: ordinary Dubins path between boundary headings;
2. an **edge**: one heading on a boundary, the other free inside its interval.
   Transversality makes the free end a "free heading" end, so the path is CS
   or CC. The end-fixed case is done by flying the path backwards;
3. the **interior**: both headings free. S, single arcs C, and CC paths that
   switch at the midpoint of AB (`λθ` vanishes on a line through A, the
   switch point and B), all in closed form.

These are the candidates of Manyam et al. for the two-point interval problem;
no numerical search is needed.

Every candidate is a real path, so the solver can't undershoot the true
minimum. The tests check it never exceeds a brute-force grid search.

The pipeline runs the C++ port (`src/dubins/_interval.cpp`, about 2–3 µs per
solve against about 0.5 ms in Python). Build it with
`python tools/build_interval.py`; without it, the Python solver is used and a
warning is printed. The tests check that the two give the same lengths.

## Reusing previous work

Every interval problem solved is cached. An edge depends only on its two
targets and two intervals, so when the refinement splits one interval per
target, all other edges are cache hits and only edges touching the new
pieces are solved. Every edge of the lower-bound graph is exact: no
estimates, one LKH solve per bound.

## Solver

Both bounds use LKH (via `elkai`, deterministic: LKH's default seed is 1).
`solver: lkh` is the only accepted value. For the lower bound, costs are
rounded **down** to integers and the value is the tour's integer cost divided
by `cost_scale`, so rounding never pushes the estimate up; LKH missing the
optimum still can. `tools/lb_stability.py` shows how much the estimate moves
with more LKH restarts.

## Tests

* every Dubins word really ends at the goal configuration;
* the interval solver is never above brute force, and its path really exists;
* Noon–Bean + each solver (LKH, and the exact MILP as reference) reproduces the
  brute-force one-in-a-set optimum;
* with one 360° interval the lower bound equals the Euclidean TSP;
* the upper-bound tour has no kinks;
* a refinement step adds headings inside the chosen interval and keeps the
  interval partition covering the circle exactly once;
* the LB estimate stays below the upper bound on the small test instance;
* best UB / max LB estimate / gap are reported as defined above, with a
  warning when an LB estimate exceeds the best UB;
* LKH is deterministic;
* interval costs are cached between iterations;
* the C++ interval solver gives the Python reference's lengths.

## Changes from the Julia version

* Successor shift `(V1 + 1) % N` produced heading interval `[-360°/N, 0]` for
  one option; the interval check then rejected most candidates, making some
  lower-bound edge costs too high. Intervals now wrap properly.
* The 3-segment search never updated its running minimum and used a global
  variable that leaked between calls.
* The turning radius was hard-coded to 1 in several places.
* A lower bound computed with LKH was reported as a bound; it is now
  labelled as an estimate.
* Big M is strictly larger than any tour cost (`n · max + 1`).
* Decoding checks the tour really visits each cluster as one block.
* Dubins formulas are guarded against rounding noise (a zero-length turn
  coming out as a full circle, a zero-length straight being rejected).
* Runs on Linux / macOS / Windows (no bundled `.exe`).
