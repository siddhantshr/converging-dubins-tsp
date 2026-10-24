"""Classical Dubins shortest path between two configurations.

A configuration is (x, y, heading). The vehicle moves forward at unit speed and
cannot turn tighter than radius r. Dubins (1957) showed the shortest path is one
of six "words" made of left turns (L), right turns (R) and straight lines (S):
LSL, RSR, LSR, RSL, RLR, LRL. We compute all six and keep the shortest.

The formulas follow the standard normalised form (as in A. Walker's dubins
library): rotate and scale so the start is at the origin, the goal is on the
+x axis at distance d = D / r, and the headings become alpha and beta.
Each word returns (t, p, q):
    * for a turn segment: the turned angle in radians,
    * for a straight segment: its length divided by r.
So the path length is r * (t + p + q).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..utils.angles import mod2pi, turn_angle as _turn

Config = tuple[float, float, float]  # (x, y, heading)

WORDS = ("LSL", "RSR", "LSR", "RSL", "RLR", "LRL")


@dataclass(frozen=True)
class DubinsPath:
    start: Config
    word: str                          # e.g. "LSR"
    params: tuple[float, float, float]  # see module docstring
    radius: float

    @property
    def length(self) -> float:
        return self.radius * sum(self.params)


_ROUND = 1e-10  # p^2 values this close below 0 are rounding noise, not "no path"


def _sqrt(p_sq: float):
    """sqrt that treats tiny negative rounding noise as 0; None if truly negative."""
    if p_sq < -_ROUND:
        return None
    return math.sqrt(max(p_sq, 0.0))


def _word_params(word: str, alpha: float, beta: float, d: float):
    """Normalised segment parameters (t, p, q) for one word, or None if it doesn't exist."""
    sa, sb = math.sin(alpha), math.sin(beta)
    ca, cb = math.cos(alpha), math.cos(beta)
    c_ab = math.cos(alpha - beta)

    if word == "LSL":
        p = _sqrt(2 + d * d - 2 * c_ab + 2 * d * (sa - sb))
        if p is None:
            return None
        tmp = math.atan2(cb - ca, d + sa - sb)
        return _turn(tmp - alpha), p, _turn(beta - tmp)

    if word == "RSR":
        p = _sqrt(2 + d * d - 2 * c_ab + 2 * d * (sb - sa))
        if p is None:
            return None
        tmp = math.atan2(ca - cb, d - sa + sb)
        return _turn(alpha - tmp), p, _turn(tmp - beta)

    if word == "LSR":
        p = _sqrt(-2 + d * d + 2 * c_ab + 2 * d * (sa + sb))
        if p is None:
            return None
        tmp = math.atan2(-ca - cb, d + sa + sb) - math.atan2(-2.0, p)
        return _turn(tmp - alpha), p, _turn(tmp - mod2pi(beta))

    if word == "RSL":
        p = _sqrt(-2 + d * d + 2 * c_ab - 2 * d * (sa + sb))
        if p is None:
            return None
        tmp = math.atan2(ca + cb, d - sa - sb) - math.atan2(2.0, p)
        return _turn(alpha - tmp), p, _turn(beta - tmp)

    if word == "RLR":
        tmp = (6.0 - d * d + 2 * c_ab + 2 * d * (sa - sb)) / 8.0
        if abs(tmp) > 1 + _ROUND:
            return None
        tmp = max(-1.0, min(1.0, tmp))
        phi = math.atan2(ca - cb, d - sa + sb)
        p = 2 * math.pi - math.acos(tmp)  # middle arc, in [pi, 2*pi]
        t = _turn(alpha - phi + mod2pi(p / 2.0))
        return t, p, _turn(alpha - beta - t + mod2pi(p))

    if word == "LRL":
        tmp = (6.0 - d * d + 2 * c_ab + 2 * d * (sb - sa)) / 8.0
        if abs(tmp) > 1 + _ROUND:
            return None
        tmp = max(-1.0, min(1.0, tmp))
        phi = math.atan2(ca - cb, d + sa - sb)
        p = 2 * math.pi - math.acos(tmp)  # middle arc, in [pi, 2*pi]
        t = _turn(-alpha - phi + p / 2.0)
        return t, p, _turn(mod2pi(beta) - alpha - t + mod2pi(p))

    raise ValueError(f"unknown Dubins word {word!r}")


def all_paths(q0: Config, q1: Config, r: float) -> list[DubinsPath]:
    """Every Dubins word that connects q0 to q1 (between 4 and 6 of them)."""
    dx, dy = q1[0] - q0[0], q1[1] - q0[1]
    d = math.hypot(dx, dy) / r
    theta = mod2pi(math.atan2(dy, dx))
    alpha = mod2pi(q0[2] - theta)
    beta = mod2pi(q1[2] - theta)

    paths = []
    for word in WORDS:
        params = _word_params(word, alpha, beta, d)
        if params is not None:
            paths.append(DubinsPath(tuple(q0), word, tuple(float(v) for v in params), r))
    return paths


def shortest_path(q0: Config, q1: Config, r: float) -> DubinsPath:
    return min(all_paths(q0, q1, r), key=lambda path: path.length)


def shortest_length(q0: Config, q1: Config, r: float) -> float:
    return shortest_path(q0, q1, r).length


# --- Following a path -------------------------------------------------------

def _advance(x: float, y: float, th: float, kind: str, amount: float, r: float):
    """Move along one segment. `amount` is an angle for L/R and a length for S."""
    if kind == "L":
        return (x + r * (math.sin(th + amount) - math.sin(th)),
                y - r * (math.cos(th + amount) - math.cos(th)),
                th + amount)
    if kind == "R":
        return (x - r * (math.sin(th - amount) - math.sin(th)),
                y + r * (math.cos(th - amount) - math.cos(th)),
                th - amount)
    return x + amount * math.cos(th), y + amount * math.sin(th), th


def _segment_amount(kind: str, param: float, r: float) -> float:
    return param if kind in "LR" else param * r


def endpoint(path: DubinsPath) -> Config:
    """Configuration reached at the end of the path."""
    x, y, th = path.start
    for kind, param in zip(path.word, path.params):
        x, y, th = _advance(x, y, th, kind, _segment_amount(kind, param, path.radius), path.radius)
    return x, y, float(mod2pi(th))


def sample(path: DubinsPath, step: float) -> np.ndarray:
    """Points (N x 2) along the path, roughly `step` apart, for plotting."""
    points = [path.start[:2]]
    x, y, th = path.start
    r = path.radius
    for kind, param in zip(path.word, path.params):
        seg_len = param * r  # arc length or straight length
        n = max(1, int(math.ceil(seg_len / step)))
        for i in range(1, n + 1):
            frac = i / n
            px, py, _ = _advance(x, y, th, kind, _segment_amount(kind, param * frac, r), r)
            points.append((px, py))
        x, y, th = _advance(x, y, th, kind, _segment_amount(kind, param, r), r)
    return np.asarray(points)
