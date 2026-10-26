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

## Lower-bound graph: eager, cached

Every lower-bound edge is solved exactly (`IntervalEdgeCosts.edge`, cached
across refinement iterations), then one LKH solve per bound.
- Why: the lazy loop (parent-interval lower estimates + LKH screening rounds)
  was removed. With the C++ interval solver (~2–3 µs per solve) its extra
  LKH rounds cost more than the solves it saved (default config: 9.0 s lazy
  vs 2.7 s eager, same final bounds). Don't reintroduce it.

## Config keys

- `solver.lkh_runs` (default 1): LKH restarts, both bounds.
- Removed: `solver.exact_time_limit`, `solver: exact`, `lower_bound.lazy`,
  `lower_bound.lkh_runs_screening` (old configs that set them still load;
  the keys are ignored).

## Output

- `results/summary.txt`: per-iteration time by stage (Dubins, interval
  problems, Noon–Bean, LKH, other) and #interval solves (cache misses).
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
