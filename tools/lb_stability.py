"""How stable is the LKH lower-bound estimate? Solve one iteration's GTSP with
more and more LKH restarts and compare.

    python tools/lb_stability.py 6                 # iteration 6 of the last run
    python tools/lb_stability.py 12 --runs 1 10 50 --config config.yml

Reads the interval set of that iteration from results/iter_XX/results.json
(written by `python -m src`), computes the interval graph once (all LKH
calls see the same matrix) and solves it with each number of LKH restarts.
If the values differ, LKH is not finding the optimal GTSP tour reliably, and
the lower-bound estimate is too high by at least the difference.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.bounds import compute_lower_bound  # noqa: E402
from src.config import RESULTS_DIR, ConfigError, load_config  # noqa: E402
from src.tsp.gtsp import IntervalEdgeCosts  # noqa: E402
from src.utils.angles import HeadingInterval  # noqa: E402


def load_intervals(results_dir: Path, iteration: int):
    path = results_dir / f"iter_{iteration:02d}" / "results.json"
    if not path.exists():
        raise SystemExit(f"{path} not found: run `python -m src` first")
    record = json.loads(path.read_text())
    if "intervals_rad_per_target" not in record:
        raise SystemExit(f"{path} has no intervals_rad_per_target: rerun `python -m src`")
    intervals = [
        [HeadingInterval(start, width) for start, width in ivs]
        for ivs in record["intervals_rad_per_target"]
    ]
    return intervals, record["lower_bound"]["value"], record["best_upper"]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("iteration", type=int)
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--runs", type=int, nargs="+", default=[1, 10, 50])
    args = parser.parse_args(argv)

    try:
        cfg = load_config(args.config)
    except ConfigError as err:
        print(f"config error: {err}", file=sys.stderr)
        return 2
    intervals, run_value, best_upper = load_intervals(RESULTS_DIR, args.iteration)

    # One cost cache for all solves: the interval problems are solved once, and
    # every LKH call below gets exactly the same matrix.
    costs = IntervalEdgeCosts(cfg.points, cfg.turning_radius)
    rows = []
    for runs in args.runs:
        t0 = time.perf_counter()
        res = compute_lower_bound(dataclasses.replace(cfg, lkh_runs=runs), intervals, costs)
        rows.append((runs, res.value, res.timings.get("lkh", 0.0), time.perf_counter() - t0))

    print(
        f"\nLB stability, iteration {args.iteration}: "
        f"{sum(len(ivs) for ivs in intervals)} intervals "
        f"({costs.solves} interval problems)"
    )
    print(
        f"value from the run ({cfg.lkh_runs} restarts): {run_value:.6f}"
        f"    best UB: {best_upper:.6f}\n"
    )
    print(f"{'lkh_runs':>8}  {'LB estimate':>12}  {'LKH s':>7}  {'total s':>7}")
    for runs, value, lkh_s, total_s in rows:
        print(f"{runs:>8}  {value:>12.6f}  {lkh_s:>7.2f}  {total_s:>7.2f}")
    values = [v for _, v, _, _ in rows]
    spread = max(values) - min(values)
    print(
        f"\nspread (max - min) = {spread:.6f}"
        + ("  -> LKH result depends on the restarts" if spread > 0 else "  -> identical")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
