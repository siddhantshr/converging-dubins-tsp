"""Load and validate config.yml."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

SOLVERS = ("lkh",)  # the exact MILP solver was removed: too slow for refinement

# Every run writes here (the repository's results/ folder), whatever the config.
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


class ConfigError(ValueError):
    pass


@dataclass
class BoundConfig:
    enabled: bool
    resolution: int  # number of headings (upper) or intervals (lower)
    solver: str


@dataclass
class Config:
    points: list[tuple[float, float]]
    turning_radius: float
    heading_offset: float  # radians; rotates every heading / interval grid
    upper: BoundConfig
    lower: BoundConfig
    lkh_runs: int
    cost_scale: float
    refine_lambda: int  # headings added per target per iteration
    gap_target: float | None  # stop once gap (%) <= this; None = single iteration
    plot: bool
    use_all_cores: bool = False  # worker processes for edge costs and the two LKH solves


def _require(cond: bool, message: str):
    if not cond:
        raise ConfigError(message)


def _bound(raw: dict, name: str, count_key: str, default_solver: str) -> BoundConfig:
    section = raw.get(name, {}) or {}
    cfg = BoundConfig(
        enabled=bool(section.get("enabled", True)),
        resolution=int(section.get(count_key, 8)),
        solver=str(section.get("solver", default_solver)).lower(),
    )
    _require(cfg.resolution >= 1, f"{name}.{count_key} must be >= 1")
    _require(
        cfg.solver in SOLVERS,
        f"{name}.solver must be 'lkh' (the exact MILP solver was removed; "
        f"the lower bound is now an LKH estimate)",
    )
    return cfg


def load_config(path: str | Path) -> Config:
    path = Path(path)
    _require(path.exists(), f"config file not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}

    points = [tuple(float(c) for c in p) for p in raw.get("points", [])]
    _require(len(points) >= 2, "need at least 2 points")
    _require(all(len(p) == 2 for p in points), "every point must be [x, y]")
    _require(len(set(points)) == len(points), "points must be distinct")

    vehicle = raw.get("vehicle", {}) or {}
    solver = raw.get("solver", {}) or {}
    output = raw.get("output", {}) or {}
    refinement = raw.get("refinement", {}) or {}

    cfg = Config(
        points=points,
        turning_radius=float(vehicle.get("turning_radius", 1.0)),
        heading_offset=math.radians(float(raw.get("heading_offset_deg", 0.0))),
        upper=_bound(raw, "upper_bound", "num_headings", "lkh"),
        lower=_bound(raw, "lower_bound", "num_intervals", "lkh"),
        lkh_runs=int(solver.get("lkh_runs", 1)),
        cost_scale=float(solver.get("cost_scale", 1000)),
        refine_lambda=int(refinement.get("lambda", 2)),
        gap_target=(
            float(refinement.get("gap_target", 5.0)) if refinement.get("enabled", True) else None
        ),
        plot=bool(output.get("plot", True)),
        use_all_cores=bool(raw.get("use_all_cores", False)),
    )
    _require(cfg.turning_radius > 0, "vehicle.turning_radius must be > 0")
    _require(cfg.lkh_runs >= 1, "solver.lkh_runs must be >= 1")
    _require(cfg.cost_scale > 0, "solver.cost_scale must be > 0")
    _require(cfg.upper.enabled or cfg.lower.enabled, "enable at least one bound")
    _require(cfg.refine_lambda >= 1, "refinement.lambda must be >= 1")
    _require(
        "results_dir" not in output, "output.results_dir was removed: every run writes to results/"
    )
    _require(
        "mu" not in refinement,
        "refinement.mu was removed: the loop now runs until the gap is <= "
        "refinement.gap_target (percent)",
    )
    _require(cfg.gap_target is None or cfg.gap_target > 0, "refinement.gap_target must be > 0")
    _require(
        cfg.gap_target is None or (cfg.upper.enabled and cfg.lower.enabled),
        "refinement needs both upper_bound and lower_bound enabled",
    )
    return cfg
