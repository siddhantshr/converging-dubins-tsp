"""Dubins paths where one or both end headings are free.

These are the building blocks of the interval problem (see interval.py).

* `paths_to_point`: start configuration fixed, end heading free.
  Candidates are CS (LS, RS) and CC (LR, RL) paths. This is the classical
  "Dubins to a point" result (Bui, Soueres, Boissonnat, Laumond, 1994).
  The function is vectorised over the start heading so it can be swept cheaply.

* `paths_between_points`: both headings free.
  Candidates are the straight line S, the four single arcs C through both
  points, and the CC paths whose turns switch at the midpoint of AB
  (Manyam et al.: S, C_psi, C_psi C_psi).

Notation: a turn of side s (+1 = left, -1 = right) from a point p with
heading th goes round the circle centred at  c = p + s * r * (-sin th, cos th).
"""

from __future__ import annotations

import math

import numpy as np

from ..utils.angles import mod2pi, turn_angle

HALF_PI = 0.5 * math.pi
LEFT, RIGHT = +1, -1

# Order of the families returned by paths_to_point (kept fixed so that a
# family can be re-evaluated later by index).
FAMILIES = ("LS", "RS", "LR+", "LR-", "RL+", "RL-")


def _turn_centre(px, py, theta, side, r):
    return px - side * r * np.sin(theta), py + side * r * np.cos(theta)


def _curve_straight(ax, ay, theta0, bx, by, side, r):
    """Turn (side) then go straight to point B. Returns (length, final heading)."""
    cx, cy = _turn_centre(ax, ay, theta0, side, r)
    wx, wy = bx - cx, by - cy
    dist = np.hypot(wx, wy)
    gap = dist * dist - r * r
    # NaN when B is strictly inside the circle; tiny negatives are rounding noise.
    straight = np.where(gap >= -1e-10, np.sqrt(np.clip(gap, 0.0, None)), np.nan)
    # The straight segment leaves the circle tangentially and hits B.
    # In the frame of the straight heading h, (B - c) has components
    # (straight, -side * r), hence  h = angle(B - c) + side * atan2(r, straight).
    heading = np.arctan2(wy, wx) + side * np.arctan2(r, straight)
    arc = turn_angle(side * (heading - theta0))
    return r * arc + straight, mod2pi(heading)


def _curve_curve(ax, ay, theta0, bx, by, side, branch, r):
    """Turn (side), then turn the other way, ending at point B.

    The second circle has radius r, passes through B, and touches the first
    circle, so its centre c2 satisfies |c2 - c1| = 2r and |c2 - B| = r.
    The two intersection points of those circles are the two branches (+1/-1).
    """
    c1x, c1y = _turn_centre(ax, ay, theta0, side, r)
    vx, vy = bx - c1x, by - c1y
    dist = np.hypot(vx, vy)
    feasible = (dist >= r - 1e-12) & (dist <= 3 * r + 1e-12)
    dist_safe = np.where(dist > 0, dist, 1.0)
    ux, uy = vx / dist_safe, vy / dist_safe
    along = (3 * r * r + dist * dist) / (2 * dist_safe)  # from circle-circle intersection
    across = np.sqrt(np.clip(4 * r * r - along * along, 0.0, None))
    c2x = c1x + along * ux - branch * across * uy
    c2y = c1y + along * uy + branch * across * ux

    jx, jy = 0.5 * (c1x + c2x), 0.5 * (c1y + c2y)  # junction point
    heading_j = np.arctan2(jy - c1y, jx - c1x) + side * HALF_PI
    heading_b = np.arctan2(by - c2y, bx - c2x) - side * HALF_PI
    arc1 = turn_angle(side * (heading_j - theta0))
    arc2 = turn_angle(-side * (heading_b - heading_j))
    length = np.where(feasible, r * (arc1 + arc2), np.nan)
    return length, mod2pi(heading_b)


def paths_to_point(start_xy, theta0, target_xy, r):
    """All CS and CC paths from (start_xy, theta0) to target_xy, end heading free.

    `theta0` may be a float or an array. Returns a list (one entry per family,
    in the order of FAMILIES) of (length, final_heading) arrays; length is NaN
    where that family doesn't exist.
    """
    ax, ay = start_xy
    bx, by = target_xy
    theta0 = np.atleast_1d(np.asarray(theta0, dtype=float))
    return [
        _curve_straight(ax, ay, theta0, bx, by, LEFT, r),
        _curve_straight(ax, ay, theta0, bx, by, RIGHT, r),
        _curve_curve(ax, ay, theta0, bx, by, LEFT, +1, r),
        _curve_curve(ax, ay, theta0, bx, by, LEFT, -1, r),
        _curve_curve(ax, ay, theta0, bx, by, RIGHT, +1, r),
        _curve_curve(ax, ay, theta0, bx, by, RIGHT, -1, r),
    ]


def paths_between_points(a_xy, b_xy, r) -> list[tuple[float, float, float]]:
    """S, C and CC paths from A to B, both headings free.

    Returns a list of (length, heading at A, heading at B).
    """
    candidates = _line_and_arcs(a_xy, b_xy, r)
    # CC: with both headings free, lambda_theta vanishes at A and at B, and it
    # is zero along a straight line, so A, the switch point J and B are
    # collinear. The two tangent circles are point reflections of each other
    # through J, hence J is the midpoint of AB, the second arc is the first
    # one reflected through J, and the heading at B equals the heading at A.
    mid = (0.5 * (a_xy[0] + b_xy[0]), 0.5 * (a_xy[1] + b_xy[1]))
    for length, head_a, _ in _line_and_arcs(a_xy, mid, r)[1:]:
        candidates.append((2.0 * length, head_a, head_a))
    return candidates


def _line_and_arcs(a_xy, b_xy, r) -> list[tuple[float, float, float]]:
    """Straight line S and the four single arcs C from A to B."""
    ax, ay = a_xy
    bx, by = b_xy
    dx, dy = bx - ax, by - ay
    chord = math.hypot(dx, dy)
    line_heading = float(mod2pi(math.atan2(dy, dx)))
    candidates = [(chord, line_heading, line_heading)]  # S

    half = 0.5 * chord
    if 0 < half <= r:
        # Two circle centres of radius r through A and B.
        offset = math.sqrt(max(r * r - half * half, 0.0))
        mx, my = ax + 0.5 * dx, ay + 0.5 * dy
        nx, ny = -dy / chord, dx / chord
        for sign in (+1, -1):
            cx, cy = mx + sign * offset * nx, my + sign * offset * ny
            for side in (LEFT, RIGHT):
                head_a = math.atan2(ay - cy, ax - cx) + side * HALF_PI
                head_b = math.atan2(by - cy, bx - cx) + side * HALF_PI
                arc = float(turn_angle(side * (head_b - head_a)))
                candidates.append((r * arc, float(mod2pi(head_a)), float(mod2pi(head_b))))
    return candidates
