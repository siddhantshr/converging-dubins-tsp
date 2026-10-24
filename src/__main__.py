"""Run with:  python -m src  [--config config.yml] [--no-plot]

Output goes to the results folder (output.results_dir in config.yml):

    results/
      iter_00/  summary.txt  results.json  tour.png     <- uniform grids
      iter_01/  ...                                     <- after 1st refinement
      ...
      convergence.txt   convergence.json   convergence.png
      summary.txt                                       <- time per iteration, by stage
"""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import asdict
from pathlib import Path

from . import refine
from .bounds import compute_lower_bound, compute_upper_bound
from .config import ConfigError, load_config
from .utils.plotting import plot_convergence, plot_results
from .utils.report import (format_convergence, format_iteration, format_result,
                           format_timings, iteration_record, save_records, save_text)


def _clean(results_dir: Path) -> None:
    """Remove iteration folders from an earlier run so old results don't linger."""
    results_dir.mkdir(parents=True, exist_ok=True)
    for old in results_dir.glob("iter_*"):
        if old.is_dir():
            shutil.rmtree(old)


def _run_refinement(cfg, results_dir: Path, plot: bool) -> None:
    records = []

    def save_iteration(it):
        folder = results_dir / f"iter_{it.index:02d}"
        folder.mkdir()
        text = format_iteration(it)
        print(text + "\n")
        save_text(text, folder / "summary.txt")
        records.append(iteration_record(it))
        save_records(records[-1], folder / "results.json")
        if plot:
            plot_results(cfg.points, [it.upper, it.lower], cfg.turning_radius,
                         str(folder / "tour.png"),
                         suptitle=f"Iteration {it.index}  |  gap {it.gap_percent:.2f}%  "
                                  f"(lower bound is an LKH estimate)")

    iterations = refine.run(cfg, cfg.refine_lambda, cfg.refine_mu, on_iteration=save_iteration)

    table = format_convergence(iterations)
    print(table)
    save_text(f"lambda = {cfg.refine_lambda}, mu = {cfg.refine_mu}\n\n{table}",
              results_dir / "convergence.txt")
    timings = format_timings(iterations)
    print("\n" + timings)
    save_text(timings, results_dir / "summary.txt")
    save_records(records, results_dir / "convergence.json")
    if plot:
        plot_convergence(iterations, str(results_dir / "convergence.png"))


def _run_single(cfg, results_dir: Path, plot: bool) -> None:
    """Only one bound enabled: compute it once, no refinement."""
    results = [compute_upper_bound(cfg)] if cfg.upper.enabled else [compute_lower_bound(cfg)]
    folder = results_dir / "iter_00"
    folder.mkdir()
    print(format_result(results[0]))
    save_text(format_result(results[0]), folder / "summary.txt")
    save_records([asdict(r) for r in results], folder / "results.json")
    if plot:
        plot_results(cfg.points, results, cfg.turning_radius, str(folder / "tour.png"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Upper and lower bounds for the Dubins TSP")
    parser.add_argument("--config", default="config.yml", help="path to config file")
    parser.add_argument("--no-plot", action="store_true", help="skip the plots")
    args = parser.parse_args(argv)

    try:
        cfg = load_config(args.config)
    except ConfigError as err:
        print(f"config error: {err}", file=sys.stderr)
        return 2

    results_dir = Path(cfg.results_dir)
    _clean(results_dir)
    plot = cfg.plot and not args.no_plot
    print(f"{len(cfg.points)} targets, turning radius {cfg.turning_radius}, "
          f"lambda = {cfg.refine_lambda}, mu = {cfg.refine_mu}\n")

    if cfg.upper.enabled and cfg.lower.enabled:
        _run_refinement(cfg, results_dir, plot)
    else:
        _run_single(cfg, results_dir, plot)
    print(f"\nresults written to {results_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
