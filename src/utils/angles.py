"""Angle helpers.

Convention used everywhere in this package: a heading is an angle in radians,
measured counter-clockwise from the +x axis.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

TWO_PI = 2.0 * math.pi


def mod2pi(angle):
    """Wrap an angle (float or numpy array) into [0, 2*pi)."""
    return np.mod(angle, TWO_PI)


def turn_angle(angle, snap: float = 1e-7):
    """Wrap a turn angle into [0, 2*pi), snapping values within `snap` of
    2*pi to 0.

    Without the snap, a turn that should be exactly 0 (e.g. when a heading is
    exactly tangent to a turning circle) can come out as a phantom full circle
    because of rounding.
    """
    value = mod2pi(angle)
    return (
        np.where(value > TWO_PI - snap, 0.0, value)
        if np.ndim(value)
        else (0.0 if value > TWO_PI - snap else float(value))
    )


@dataclass(frozen=True)
class HeadingInterval:
    """The closed set of headings [start, start + width], swept counter-clockwise.

    The interval may wrap past 2*pi, e.g. start=5.5, width=1.0 covers
    [5.5, 2*pi) together with [0, 0.217].
    """

    start: float
    width: float

    @property
    def end(self) -> float:
        return self.start + self.width

    @property
    def is_full_circle(self) -> bool:
        return self.width >= TWO_PI - 1e-12

    @property
    def boundaries(self) -> tuple[float, ...]:
        """Boundary headings. A full circle has a single (arbitrary) boundary."""
        if self.is_full_circle:
            return (self.start,)
        return (self.start, self.end)

    def contains(self, theta, tol: float = 1e-9):
        """True where theta lies in the interval (works on floats and arrays)."""
        if self.is_full_circle:
            return np.ones_like(np.asarray(theta, dtype=float), dtype=bool)
        offset = mod2pi(np.asarray(theta, dtype=float) - self.start)
        # An offset just below 2*pi means theta sits just below `start`,
        # i.e. on the boundary up to rounding.
        return (offset <= self.width + tol) | (offset >= TWO_PI - tol)

    def __str__(self) -> str:
        return f"[{math.degrees(self.start):.1f}°, {math.degrees(self.end):.1f}°]"
