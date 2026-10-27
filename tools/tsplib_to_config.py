"""Turn a TSPLIB instance into a config for this project.

    python tools/tsplib_to_config.py examples/tsplib/berlin52.tsp      # -> examples/berlin52.yml

TSPLIB: https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/tspindex.html

The Dubins TSP needs planar points, a turning radius and a cost scale; TSPLIB
gives only the first, and not always:

* EUC_2D: the coordinates are used as they are.
* GEO: the coordinates (latitude, longitude in TSPLIB's DDD.MM form) are used
  as plain x, y. The instance's GEO distances do not apply.
* EXPLICIT (FULL_MATRIX): there are no coordinates. Classical multidimensional
  scaling (MDS) places the cities in the plane so that straight-line distances
  match the matrix as closely as 2-D allows. The result is an approximation of
  the instance, and the config header says how close it is.

Turning radius: r = RADIUS_FACTOR x the mean nearest-neighbour distance
(examples/big_75.yml has the same ratio: r = 1, mean nearest neighbour 2.48),
rounded to 2 significant figures.

Cost scale: the largest power of 10 (at most 1000) that keeps the Noon-Bean
matrix at most half of LKH's integer limit (see src/tsp/noon_bean.py). The
lower bound rounds every edge down to a multiple of 1 / cost_scale, so the
estimate can be at most n / cost_scale too low; the header gives that number.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.tsp.noon_bean import LKH_MAX_ENTRY  # noqa: E402

RADIUS_FACTOR = 0.4
SETTINGS = """\
upper_bound:
  enabled: true
  num_headings: 6
  solver: lkh

lower_bound:
  enabled: true
  num_intervals: 6
  solver: lkh

refinement:
  enabled: true
  lambda: 1
  gap_target: 3.0

solver:
  lkh_runs: 1
  cost_scale: {cost_scale:g}

output:
  plot: true
"""


def read_tsplib(path: Path) -> tuple[dict, list[str]]:
    """Header fields (KEY: value) and the lines of the data section."""
    header, data, in_data = {}, [], False
    for line in path.read_text().splitlines():
        text = line.strip()
        if not text or text == "EOF":
            continue
        if text.endswith("SECTION"):
            header["SECTION"] = text
            in_data = True
        elif in_data:
            data.append(text)
        else:
            key, _, value = text.partition(":")
            header[key.strip()] = value.strip()
    return header, data


def classical_mds(dist: np.ndarray) -> tuple[np.ndarray, str]:
    """2-D points whose distances approximate `dist`, and a note on the fit."""
    n = len(dist)
    centre = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centre @ (dist**2) @ centre
    values, vectors = np.linalg.eigh(gram)
    top = np.argsort(values)[::-1][:2]
    coords = vectors[:, top] * np.sqrt(np.maximum(values[top], 0.0))
    coords *= np.sign(coords[np.argmax(np.abs(coords), axis=0), [0, 1]])  # fixed orientation
    embedded = np.hypot(*(coords[:, None, :] - coords[None, :, :]).transpose(2, 0, 1))
    off = ~np.eye(n, dtype=bool)
    rel = np.abs(embedded[off] - dist[off]) / dist[off]
    positive = values[values > 0]
    note = (
        f"MDS fit: the 2 kept dimensions carry {values[top].sum() / positive.sum():.0%} of the "
        f"spread; distance error median {np.median(rel):.0%}, 90th percentile "
        f"{np.percentile(rel, 90):.0%}, max {rel.max():.0%}"
    )
    return coords, note


def points_from(header: dict, data: list[str]) -> tuple[np.ndarray, list[str]]:
    kind = header.get("EDGE_WEIGHT_TYPE")
    n = int(header["DIMENSION"])
    if kind in ("EUC_2D", "GEO"):
        coords = np.array([[float(v) for v in line.split()[1:3]] for line in data[:n]])
        if kind == "EUC_2D":
            return coords, ["EUC_2D: TSPLIB coordinates, used as they are."]
        return coords, [
            "GEO: TSPLIB's latitude / longitude (DDD.MM) used as plain x, y;",
            "the instance's GEO distances do not apply here.",
        ]
    if kind == "EXPLICIT" and header.get("EDGE_WEIGHT_FORMAT") == "FULL_MATRIX":
        values = np.array([float(v) for line in data for v in line.split()])
        coords, note = classical_mds(values.reshape(n, n))
        return coords, [
            "EXPLICIT: TSPLIB gives only a distance matrix, no coordinates. These",
            "points come from classical MDS, so this is an APPROXIMATION of the instance.",
            note + ".",
        ]
    raise SystemExit(f"unsupported instance: EDGE_WEIGHT_TYPE {kind}")


def nice(x: float, digits: int = 2) -> float:
    """x rounded to `digits` significant figures."""
    return float(f"{x:.{digits}g}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("tsp", type=Path)
    parser.add_argument("--out", type=Path, help="default: examples/<name>.yml")
    args = parser.parse_args(argv)

    header, data = read_tsplib(args.tsp)
    name = header.get("NAME", args.tsp.stem)
    coords, notes = points_from(header, data)
    n = len(coords)
    decimals = max(0, 2 - int(math.floor(math.log10(np.ptp(coords)))) + 2)
    coords = np.round(coords, decimals)
    if len({tuple(p) for p in coords}) != n:
        raise SystemExit("two points coincide after rounding")

    dist = np.hypot(*(coords[:, None, :] - coords[None, :, :]).transpose(2, 0, 1))
    np.fill_diagonal(dist, np.inf)
    mean_nn = dist.min(axis=1).mean()
    r = nice(RADIUS_FACTOR * mean_nn)

    # Noon-Bean's M is about the sum over targets of the longest edge leaving it;
    # Dubins edges are at most the distance plus a few turns. Keep 2M <= half the limit.
    np.fill_diagonal(dist, 0.0)
    m_estimate = float(dist.max(axis=1).sum() + n * (2 * math.pi + 4) * r)
    cost_scale = 1000.0
    while 2 * m_estimate * cost_scale > LKH_MAX_ENTRY / 2:
        cost_scale /= 10
    if cost_scale < 1e-3:
        raise SystemExit("instance too large for LKH's integer costs")

    out = args.out or Path(__file__).resolve().parent.parent / "examples" / f"{name}.yml"
    lines = [
        f"# TSPLIB {name}: {header.get('COMMENT', '')} ({n} targets).",
        "# Source: https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/tspindex.html",
        f"# Generated by tools/tsplib_to_config.py from examples/tsplib/{args.tsp.name}.",
        *[f"# {note}" for note in notes],
        f"# Turning radius = {RADIUS_FACTOR} x mean nearest-neighbour distance "
        f"({mean_nn:.4g}) = {r:g}; TSPLIB has none.",
        f"# cost_scale {cost_scale:g} keeps LKH's integers in range; the LB estimate can be at",
        f"# most {n / cost_scale:g} too low from rounding.",
        "# Run with:",
        f"#     python -m src --config examples/{name}.yml",
        "",
        "vehicle:",
        f"  turning_radius: {r:g}",
        "",
        "points:",
        *[f"  - [{x:.{decimals}f}, {y:.{decimals}f}]" for x, y in coords],
        "",
        "heading_offset_deg: 0.0",
        "parallelise: true",
        "",
        SETTINGS.format(cost_scale=cost_scale),
    ]
    out.write_text("\n".join(lines))
    print(f"{out}: {n} targets, r = {r:g}, cost_scale = {cost_scale:g}")
    for note in notes:
        print(f"  {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
