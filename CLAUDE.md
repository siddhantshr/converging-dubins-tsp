# CLAUDE.md

Upper and lower bounds for the Dubins TSP via one-in-a-set TSP + Noon–Bean +
LKH, with lower-bound-driven refinement (lambda, mu). See README.md for the method.

## Commands

```bash
source venv/bin/activate
python tools/build_interval.py                     # compile the C++ interval solver
python -m src [--config config.yml] [--no-plot]   # writes results/
python -m pytest -q                                # ~20 s
python tools/lb_stability.py <iteration>           # needs results/ from a run
```

## The lower bound is an LKH estimate

- Both bounds are solved with LKH (elkai) only. `solver:` accepts only `lkh`.
- Why: the exact MILP (CBC via pulp) made the lower bound certified, but it
  took over half the runtime and grew fast with the refinement
  (14–19 s per iteration at 192 nodes). It was removed from the pipeline.
- Consequence: when LKH misses the optimal GTSP tour, the "lower bound" is too
  high. Always label it "estimate (LKH)" (`BoundResult.status`), never as a bound.
- Gap = (best UB − **this iteration's** LB estimate) / this iteration's LB
  estimate. Do not use the max of LB estimates across iterations for the gap
  (it is kept as `max_lb_estimate`, for reference only). If an LB estimate
  exceeds the best UB, `refine.run` and the reports print a warning.
- `ExactATSP` in `src/tsp/solvers.py` (and the pulp dependency) stays only as
  the reference solver in `tests/test_noon_bean.py`. Don't reconnect it to `bounds.py`.
- elkai has no seed argument. LKH's default `SEED = 1` makes runs
  deterministic (`test_lkh_is_deterministic`).

## Lazy loop (`src/bounds.py::_solve_graph`)

With `lower_bound.lazy: true`, edges start as parent-interval lower estimates.
Rounds that only pick which edges to make exact use LKH with
`lower_bound.lkh_runs_screening` restarts. Once a screening tour is all-exact,
the graph is re-solved with `solver.lkh_runs`. Stop if that tour is all-exact,
else refresh its edges and go back to screening.

## Config keys

- `solver.lkh_runs` (default 1): upper bound and the final lower-bound solve.
  When it equals `lkh_runs_screening` (the default: both 1) there is no separate
  final solve: every lazy round is one 1-restart LKH call, stop when all-exact.
- `lower_bound.lkh_runs_screening` (default 1): lazy screening rounds.
- `lower_bound.lazy` (default true).
- Removed: `solver.exact_time_limit`, `solver: exact`.

## Output

- `results/summary.txt`: per-iteration time by stage (Dubins, interval
  problems, Noon–Bean, LKH screening / final, other), #interval solves, LKH rounds.
- `results/convergence.txt|json|png`: bounds, gap, warnings.
- `results/iter_XX/results.json`: includes `intervals_rad_per_target`, which
  `tools/lb_stability.py` reads.

## Interval solver: C++ + Python reference

- `src/dubins/_interval.cpp` is a line-by-line port of `interval.py`,
  `path.py::shortest_length` and `free_heading.py`. `solve_interval_problem`
  uses it; `solve_interval_problem_py` is the reference
  (`test_cpp_matches_python`). Change both together: same candidates, order,
  tolerances and strict-`<` tie-break. Rebuild after editing the .cpp.
- The .so is not committed (`.gitignore`). Without it, the Python solver runs
  (~150x slower) with a warning.
- Never build with `-ffast-math` (the solver relies on NaN checks).

## Don't change without asking

The interval solver (`src/dubins/interval.py`, `src/dubins/_interval.cpp`) and
the refinement rule (`src/refine.py::refine`).
