"""Every Dubins word must actually end at the goal configuration."""
import math
import random

from src.dubins.path import all_paths, endpoint, shortest_path


def _angle_diff(a, b):
    return abs((a - b + math.pi) % (2 * math.pi) - math.pi)


def test_every_word_reaches_goal():
    rng = random.Random(0)
    for _ in range(3000):
        q0 = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(0, 2 * math.pi))
        q1 = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(0, 2 * math.pi))
        r = rng.uniform(0.3, 2.0)
        for path in all_paths(q0, q1, r):
            x, y, th = endpoint(path)
            assert math.hypot(x - q1[0], y - q1[1]) < 1e-7, path.word
            assert _angle_diff(th, q1[2]) < 1e-7, path.word


def test_straight_ahead_is_a_line():
    path = shortest_path((0, 0, 0), (5, 0, 0), 1.0)
    assert math.isclose(path.length, 5.0, abs_tol=1e-9)


def test_no_phantom_full_circle():
    # Goal heading exactly tangent: the last turn must be 0, not 2*pi.
    path = shortest_path((0, 0, math.pi / 2), (-1, 1, math.pi), 1.0)  # quarter left turn
    assert math.isclose(path.length, math.pi / 2, abs_tol=1e-9)
