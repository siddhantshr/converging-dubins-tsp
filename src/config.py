"""Load and validate config.yml."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

SOLVERS = ("lkh",)  # the exact MILP solver was removed: too slow for refinement


class ConfigError(ValueError):
    pass


@dataclass
class BoundConfig:
    enabled: bool
    resolution: int        # number of headings (upper) or intervals (lower)
    solver: str
    lazy: bool = False     # lower bound only: lazy edge evaluation


@dataclass
class Config:
    points: list[tuple[float, float]]
    turning_radius: float
    heading_offset: float           # radians; rotates every heading / interval grid
    upper: BoundConfig
    lower: BoundConfig
    sweep_samples: int
    lkh_runs: int
    cost_scale: float
    refine_lambda: int              # headings added per target per iteration
    refine_mu: int                  # number of refinement iterations (0 = single run)
    plot: bool
    results_dir: str
    lkh_runs_screening: int = 1     # LKH restarts in the lazy screening rounds (lower bound)


def _require(cond: bool, message: str):
    if not cond:
        raise ConfigError(message)


def _bound(raw: dict, name: str, count_key: str, default_solver: str) -> BoundConfig:
    section = raw.get(name, {}) or {}
    cfg = BoundConfig(
        enabled=bool(section.get("enabled", True)),
        resolution=int(section.get(count_key, 8)),
        solver=str(section.get("solver", default_solver)).lower(),
        lazy=bool(section.get("lazy", True)),
    )
    _require(cfg.resolution >= 1, f"{name}.{count_key} must be >= 1")
    _require(cfg.solver in SOLVERS,
             f"{name}.solver must be 'lkh' (the exact MILP solver was removed; "
             f"the lower bound is now an LKH estimate)")
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
        sweep_samples=int((raw.get("lower_bound", {}) or {}).get("sweep_samples", 128)),
        lkh_runs=int(solver.get("lkh_runs", 1)),
        cost_scale=float(solver.get("cost_scale", 1000)),
        refine_lambda=int(refinement.get("lambda", 2)),
        refine_mu=int(refinement.get("mu", 0)) if refinement.get("enabled", True) else 0,
        plot=bool(output.get("plot", True)),
        results_dir=str(output.get("results_dir", "results")),
        lkh_runs_screening=int((raw.get("lower_bound", {}) or {}).get("lkh_runs_screening", 1)),
    )
    _require(cfg.turning_radius > 0, "vehicle.turning_radius must be > 0")
    _require(cfg.sweep_samples >= 8, "lower_bound.sweep_samples must be >= 8")
    _require(cfg.lkh_runs >= 1, "solver.lkh_runs must be >= 1")
    _require(cfg.lkh_runs_screening >= 1, "lower_bound.lkh_runs_screening must be >= 1")
    _require(cfg.cost_scale > 0, "solver.cost_scale must be > 0")
    _require(cfg.upper.enabled or cfg.lower.enabled, "enable at least one bound")
    _require(cfg.refine_lambda >= 1, "refinement.lambda must be >= 1")
    _require(cfg.refine_mu >= 0, "refinement.mu must be >= 0")
    _require(cfg.refine_mu == 0 or (cfg.upper.enabled and cfg.lower.enabled),
             "refinement needs both upper_bound and lower_bound enabled")
    return cfg
