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
version fixed (see the end of this file). The Dubins edge costs of both bounds
are computed by a small C++ extension; everything else is Python.

## Quick start

```bash
pip install -r requirements.txt
python tools/build_interval.py   # compile the C++ extension (needs a C++17 compiler)
python -m src                    # uses ./config.yml
python -m src --config other.yml --no-plot
python -m src --config examples/big_75.yml   # 75 targets
python -m pytest                 # run the checks
ruff check . && ruff format .    # lint + format (settings in ruff.toml)
python tools/lb_stability.py 8   # LB estimate of iteration 8 with 1 / 10 / 50 LKH restarts
```

Edit `config.yml` to set the points, the turning radius, the number of
discrete headings (upper bound), the number of heading intervals (lower bound)
and the refinement settings `lambda` and `gap_target`. `parallelise: true`
(the default) runs the two LKH solves of an iteration at the same time on big
problems; the results are the same.

The two shipped configs:

| | `config.yml` | `examples/big_75.yml` |
|---|---|---|
| targets | 6 | 75 (random, seed 75, 30 × 30 square) |
| headings / intervals per target | 8 / 8 | 6 / 6 |
| `lambda` | 2 | 1 |
| `gap_target` | 2 % | 3 % |
| result | 9 iterations, gap 1.90 %, ~1 s | 5 iterations, gap 1.98 %, ~23 s |

## Output

Every run writes to the repository's `results/` folder, whichever config it
uses, and replaces the previous run's results there. The folder is not
tracked by git.

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

## Iterative refinement (lambda, gap_target)

The lower-bound tour picks one heading interval at every target. Each
refinement iteration, for every target:

1. place `lambda` evenly spaced headings strictly inside that interval;
2. add them to the target's upper-bound heading set;
3. split the interval at the same headings (lower bound gets finer there);

then both bounds are solved again. This repeats until the gap is at most
`gap_target` percent (5 if the key is missing). There is no iteration limit: a run ends
only when the target is met, so a very small target can run for a long time.
Results of each iteration are written as soon as it finishes.

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

## Architecture

### The big picture

One run of `python -m src` is a loop. Each pass of the loop is one
**iteration**, and every iteration does the same three things:

1. compute the **upper bound** (a real, flyable tour);
2. compute the **lower bound estimate** (a tour that may "kink" at targets);
3. **refine**: use the lower-bound tour to decide where to add detail for the
   next iteration.

```
config.yml
    │
    ▼
 start: num_headings headings, num_intervals intervals per target
    │
    ▼
┌──────────────────────── one iteration ────────────────────────┐
│                                                               │
│   upper bound  ── same 5-step pipeline ──►  real tour length  │
│   lower bound  ── same 5-step pipeline ──►  estimate          │
│                                                               │
│   refine: at every target, split the interval the lower-bound │
│   tour used, and add the same headings to the upper bound     │
│                                                               │
└───────────────────────────────┬───────────────────────────────┘
                                │  repeat until gap <= gap_target
                                ▼
                 results/ (text, JSON, plots)
```

Iteration 0 uses the uniform grids. Each later iteration adds `lambda`
options per target, so the problems grow steadily. Default config (6
targets, 8 options each, `lambda = 2`): 48 options in iteration 0, 144 in
iteration 8. `examples/big_75.yml` (75 targets, 6 options each,
`lambda = 1`): 450 in iteration 0, then 525, 600, 675, 750.

### The pipeline (shared by both bounds)

Both bounds are the same kind of problem: **at every target, pick one option
(a heading, or a heading interval), and pick the order of the targets, so
that the closed tour is as short as possible.** This is a "one-in-a-set" (or
generalised) TSP. Both bounds go through the same five steps, in
`src/bounds.py`:

```
 1. build graph  ─►  2. integer costs  ─►  3. Noon–Bean  ─►  4. LKH  ─►  5. decode
   (gtsp.py)          (bounds.py)          (noon_bean.py)    (solvers.py)  (noon_bean.py)
```

**1. Build the graph** (`src/tsp/gtsp.py`). One node per (target, option).
Every target's nodes form a **cluster**. There is an edge between every two
nodes of different targets; its cost is the shortest way to fly between
them. The two bounds differ only here:

| | Upper bound | Lower bound |
|---|---|---|
| option at a target | one heading | a heading interval |
| edge cost | Dubins path between two fixed headings (`dubins/path.py`, run in C++) | shortest Dubins path with *any* heading inside each interval: the **interval problem** (`dubins/interval.py`, run in C++) |
| cache | dict keyed by (targets, headings) | `IntervalEdgeCosts`, keyed by (targets, intervals) |

Both caches live for the whole run. When the refinement splits one interval
per target, only edges touching the new pieces are new; everything else is
a cache hit. (`examples/big_75.yml`: 199,800 interval problems in
iteration 0; in iteration 4, 199,800 of the 555,000 edges.)

**2. Integer costs** (`bounds._integer_costs`). LKH works with integers, so
every cost is multiplied by `solver.cost_scale` (default 1000). The upper
bound rounds to the nearest integer. The lower bound rounds **down**, so
rounding can only make the estimate smaller, never too high.

**3. Noon–Bean transformation** (`src/tsp/noon_bean.py`). LKH solves the
ordinary TSP, which visits *every* node. We want *one* node per cluster.
The Noon–Bean transformation rewrites the costs so that the two problems
have the same answer:

- inside each cluster, the nodes are joined in a free (zero-cost) ring, so
  the tour can pass through the unused options for nothing;
- each edge leaving a cluster is charged as if it left from the ring
  neighbour, so a cluster costs exactly what its first-visited option costs;
- a large constant M is added to every edge between clusters, so a tour that
  enters a cluster more than once is always worse. M only has to beat any
  real tour, so it is the sum over targets of each target's largest
  outgoing cost, plus 1;
- every other edge inside a cluster is forbidden, with cost 2M.

The resulting matrix has the same number of nodes as the graph. It is
asymmetric (cost A→B ≠ B→A). LKH multiplies every cost by 100 internally
in 32-bit integers and crashes when an entry reaches about 21 million
(measured), so entries are kept at most about 10.7 million. If a problem
needs more, the run stops and asks for a smaller `cost_scale`. (75 targets
in a 30 × 30 square: largest entry about 5 million.)

**4. LKH** (`src/tsp/solvers.py::solve_lkh`). The matrix goes to LKH, a
fast heuristic TSP solver, through the `elkai` package (which hands it to
LKH as a text problem). `solver.lkh_runs` sets the number of restarts
(default 1). LKH is deterministic (fixed seed), so runs are reproducible.
LKH is the one heuristic step: it usually finds the optimal tour, but not
always.

**5. Decode** (`noon_bean.decode`). The first node the tour visits in each
cluster is the chosen option; the order of the clusters is the visiting
order. Then:

- **upper bound** = sum of the real Dubins lengths of the legs. This tour
  can be flown, so it is always a valid upper bound, whatever LKH did;
- **lower bound estimate** = the tour's integer cost ÷ `cost_scale`. It is a
  true lower bound only if LKH found the optimal tour. When LKH misses it,
  the value is too high, which is why it is always labelled
  "estimate (LKH)".

### The refinement loop (`src/refine.py`)

`refine.run` holds the state that carries over between iterations: each
target's heading list, each target's interval list, and the two edge-cost
caches. After both bounds are solved, it:

1. records the iteration: both values, the best upper bound so far, and
   the gap `(best UB − this LB estimate) / this LB estimate`. If the
   estimate is above the best upper bound, LKH must have missed the optimum,
   and a warning is printed;
2. stops if the gap is at most `gap_target`. Because the gap uses the
   lower-bound *estimate*, an estimate that is too high (LKH missed the
   optimum) makes the gap look smaller and can stop the loop early;
3. otherwise calls `refine.refine`: at every target, take the interval the lower-bound
   tour chose, cut it into `lambda + 1` equal pieces, and add the `lambda`
   cut points as new upper-bound headings.

If only one bound is enabled in `config.yml`, there is no loop: that bound
is computed once (`__main__._run_single`).

### The interval solver (lower-bound edges)

Each lower-bound edge asks: what is the shortest Dubins path from point A,
leaving with any heading in interval I0, to point B, arriving with any
heading in interval I1? The answer is always one of a short list of
candidate paths (see "Interval solver" below for why): the solver builds all
of them, throws away those whose headings fall outside the intervals, and
keeps the shortest. There is no search or iteration, so a solve is fast and
exact.

- `src/dubins/_interval.cpp` is the C++ version used by the pipeline: about
  2–3 µs per solve. It also computes the upper bound's point-to-point Dubins
  lengths (`interval.dubins_length`): about 1.7 µs per edge including the
  Python bookkeeping, against about 13 µs in pure Python. Build it once with
  `python tools/build_interval.py`.
- `src/dubins/interval.py` is the same algorithm in Python (about 0.5 ms
  per solve). It is the reference the tests compare against, and the
  fallback (with a warning) if the C++ module has not been built.

### Parallel LKH (`parallelise: true`)

LKH is the slowest step, and each iteration calls it twice (upper bound,
then lower bound). The two calls don't depend on each other, so
`bounds.compute_bounds` sends the upper-bound solve to a worker process and
builds and solves the lower bound in the meantime.

- **Only for big graphs.** The worker only starts once an upper-bound graph
  has at least 300 nodes (`MIN_PARALLEL_NODES` in `src/bounds.py`), so small
  problems, like the default 6-target config, never start it. Starting it
  took about 0.2 s here (it can be slower on Windows or on a cold first
  run), and it overlaps the lower-bound work.
- **One worker, a process.** elkai keeps Python locked while LKH runs, and
  LKH keeps global state in C, so a thread would not run in parallel.
- **Same results.** Every computation is the same as with
  `parallelise: false`; only the elapsed time changes. On
  `examples/big_75.yml` it cut the run from 44.9 s to about 30 s (measured
  before the upper-bound edges moved to C++).
- **Timings.** `summary.txt` still shows each stage's own time; the `iter`
  column and `time s` in `convergence.txt` are elapsed time.

### Where the time goes

`examples/big_75.yml` (75 targets, 6 headings / 6 intervals, `lambda = 1`,
3 % target; 5 iterations: gap 7.17 % -> 5.36 % -> 4.76 % -> 4.62 % ->
1.98 %), from `results/summary.txt`, with `parallelise: true`:

| Stage | Seconds |
|---|---|
| Upper bound: Dubins edge costs (C++) | 1.51 |
| Upper bound: LKH (in the worker, overlaps the lower bound) | 16.98 |
| Lower bound: interval problems (C++) | 1.99 |
| Lower bound: LKH | 16.83 |
| Lower bound: everything else (graph building, bookkeeping) | 2.53 |
| **Elapsed** | **23.0** |

LKH is most of the work: one call takes about 2.3 s on 450 nodes and 5–6 s
on 750. Before the upper bound's Dubins edges moved to C++ they took 8.2 s
(about 13 µs per edge in Python); the whole run took 29.9 s. Most of the
remaining graph-building time is Python bookkeeping per edge, not the C++
calls.

The default 6-target config needs about 1 s of computing (9 iterations to
its 2 % target). Its per-iteration plots are not in these numbers and add a
little more.

### Output (`src/utils/report.py`, `src/utils/plotting.py`)

After every iteration, `iter_XX/` gets a text summary, a JSON record (it
includes the interval list, which `tools/lb_stability.py` reads) and a plot
of both tours. At the end, `convergence.txt|json|png` show the bounds and
gap across iterations and `summary.txt` shows the time per stage.

### Files

```
config.yml                 points and settings
src/
  __main__.py              entry point (python -m src): runs the loop, writes results/
  config.py                reads and checks config.yml
  refine.py                the iteration loop and the refinement step
  bounds.py                the 5-step pipeline, for either bound; the LKH worker
  dubins/
    path.py                Dubins shortest path between two fixed headings (Python reference)
    free_heading.py        Dubins paths when one or both end headings are free
    interval.py            interval problem: Python reference, brute-force check, C++ calls
    _interval.cpp          C++ version used by the pipeline: interval problem + Dubins length
  tsp/
    gtsp.py                builds the clustered graphs, caches interval costs
    noon_bean.py           Noon–Bean transformation and decoding
    solvers.py             LKH via elkai (+ an exact MILP solver, used only in tests)
  utils/
    angles.py              angle wrapping, HeadingInterval
    report.py              text tables and JSON
    plotting.py            tour and convergence plots
tools/
  build_interval.py        compiles _interval.cpp
  lb_stability.py          re-solves one iteration with more LKH restarts
tests/                     correctness checks (see "Tests")
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
solve against about 0.5 ms in Python). The same module also computes the
upper bound's point-to-point Dubins lengths (`interval.dubins_length`). Build it with
`python tools/build_interval.py`; without it, the Python solver is used and a
warning is printed. The tests check that C++ and Python give the same lengths.

## Tests

`python -m pytest` (about 5 s) checks that:

* every Dubins word really ends at the goal configuration, with no phantom
  full circles;
* the interval solver is never above brute force, its path really exists,
  and intervals that wrap past 360° work;
* the C++ interval solver and the C++ Dubins length match the Python
  reference;
* Noon–Bean + each solver (LKH, and the exact MILP as reference) reproduces the
  brute-force one-in-a-set optimum, and decoding rejects broken tours
  (a cluster left twice, a forbidden edge used);
* with one 360° interval the lower bound equals the Euclidean TSP;
* the upper-bound tour has no kinks;
* the LB estimate stays below the upper bound on the small test instance;
* a refinement step adds headings inside the chosen interval and keeps the
  interval partition covering the circle exactly once;
* the loop stops at the first iteration within `gap_target`, and a config
  that still sets the removed `mu` is rejected;
* best UB / max LB estimate / gap are reported as defined above, with a
  warning when an LB estimate exceeds the best UB;
* LKH is deterministic, interval costs are cached between iterations, and
  `parallelise: true` gives exactly the same results as `false`.

## Changes from the Julia version

* Successor shift `(V1 + 1) % N` produced heading interval `[-360°/N, 0]` for
  one option; the interval check then rejected most candidates, making some
  lower-bound edge costs too high. Intervals now wrap properly.
* The 3-segment search never updated its running minimum and used a global
  variable that leaked between calls.
* The turning radius was hard-coded to 1 in several places.
* A lower bound computed with LKH was reported as a bound; it is now
  labelled as an estimate.
* Big M is strictly larger than any tour cost (sum of each cluster's
  largest outgoing cost, + 1), and small enough for LKH's integer limits.
* Decoding checks the tour really visits each cluster as one block.
* Dubins formulas are guarded against rounding noise (a zero-length turn
  coming out as a full circle, a zero-length straight being rejected).
* Runs on Linux / macOS / Windows (no bundled `.exe`).
