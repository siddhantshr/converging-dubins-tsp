// Two-point Dubins interval problem in C++: a line-by-line port of
// interval.py, path.py (shortest_length) and free_heading.py.
//
// The Python version stays the reference (tests/test_interval.py compares
// the two). Keep both in sync: same candidates, same order, same tolerances,
// same tie-break (a candidate wins only if strictly shorter).
//
// Build: python tools/build_interval.py
// Never compile with -ffast-math: the NaN checks below rely on IEEE semantics.

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cmath>
#include <initializer_list>
#include <limits>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double HALF_PI = 0.5 * PI;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr int LEFT = +1, RIGHT = -1;

// --- angles.py ---------------------------------------------------------------

// np.mod(angle, 2*pi): result has the sign of 2*pi (can round up to 2*pi).
inline double mod2pi(double a) {
    double m = std::fmod(a, TWO_PI);
    if (m != 0.0) {
        if (m < 0.0) m += TWO_PI;
    } else {
        m = 0.0;
    }
    return m;
}

inline double turn_angle(double a, double snap = 1e-7) {
    double v = mod2pi(a);
    return v > TWO_PI - snap ? 0.0 : v;
}

struct Interval {
    double start, width;
    bool full() const { return width >= TWO_PI - 1e-12; }
    int boundaries(double out[2]) const {
        out[0] = start;
        if (full()) return 1;
        out[1] = start + width;
        return 2;
    }
    bool contains(double theta, double tol = 1e-9) const {
        if (full()) return true;
        double offset = mod2pi(theta - start);
        return (offset <= width + tol) || (offset >= TWO_PI - tol);
    }
};

// --- path.py: classical Dubins shortest path --------------------------------

constexpr double ROUND = 1e-10;

// sqrt that treats tiny negative rounding noise as 0; false if truly negative.
inline bool safe_sqrt(double p_sq, double& p) {
    if (p_sq < -ROUND) return false;
    p = std::sqrt(p_sq > 0.0 ? p_sq : 0.0);
    return true;
}

enum Word { LSL, RSR, LSR, RSL, RLR, LRL };

// Normalised (t, p, q) for one word; false if the word doesn't exist.
bool word_params(Word word, double alpha, double beta, double d, double& t, double& p, double& q) {
    double sa = std::sin(alpha), sb = std::sin(beta);
    double ca = std::cos(alpha), cb = std::cos(beta);
    double c_ab = std::cos(alpha - beta);
    double tmp, phi;
    switch (word) {
    case LSL:
        if (!safe_sqrt(2 + d * d - 2 * c_ab + 2 * d * (sa - sb), p)) return false;
        tmp = std::atan2(cb - ca, d + sa - sb);
        t = turn_angle(tmp - alpha);
        q = turn_angle(beta - tmp);
        return true;
    case RSR:
        if (!safe_sqrt(2 + d * d - 2 * c_ab + 2 * d * (sb - sa), p)) return false;
        tmp = std::atan2(ca - cb, d - sa + sb);
        t = turn_angle(alpha - tmp);
        q = turn_angle(tmp - beta);
        return true;
    case LSR:
        if (!safe_sqrt(-2 + d * d + 2 * c_ab + 2 * d * (sa + sb), p)) return false;
        tmp = std::atan2(-ca - cb, d + sa + sb) - std::atan2(-2.0, p);
        t = turn_angle(tmp - alpha);
        q = turn_angle(tmp - mod2pi(beta));
        return true;
    case RSL:
        if (!safe_sqrt(-2 + d * d + 2 * c_ab - 2 * d * (sa + sb), p)) return false;
        tmp = std::atan2(ca + cb, d - sa - sb) - std::atan2(2.0, p);
        t = turn_angle(alpha - tmp);
        q = turn_angle(beta - tmp);
        return true;
    case RLR:
        tmp = (6.0 - d * d + 2 * c_ab + 2 * d * (sa - sb)) / 8.0;
        if (std::fabs(tmp) > 1 + ROUND) return false;
        tmp = std::fmax(-1.0, std::fmin(1.0, tmp));
        phi = std::atan2(ca - cb, d - sa + sb);
        p = 2 * PI - std::acos(tmp);  // middle arc, in [pi, 2*pi]
        t = turn_angle(alpha - phi + mod2pi(p / 2.0));
        q = turn_angle(alpha - beta - t + mod2pi(p));
        return true;
    case LRL:
        tmp = (6.0 - d * d + 2 * c_ab + 2 * d * (sb - sa)) / 8.0;
        if (std::fabs(tmp) > 1 + ROUND) return false;
        tmp = std::fmax(-1.0, std::fmin(1.0, tmp));
        phi = std::atan2(ca - cb, d + sa - sb);
        p = 2 * PI - std::acos(tmp);  // middle arc, in [pi, 2*pi]
        t = turn_angle(-alpha - phi + p / 2.0);
        q = turn_angle(mod2pi(beta) - alpha - t + mod2pi(p));
        return true;
    }
    return false;
}

double shortest_length(double x0, double y0, double th0, double x1, double y1, double th1, double r) {
    double dx = x1 - x0, dy = y1 - y0;
    double d = std::hypot(dx, dy) / r;
    double theta = mod2pi(std::atan2(dy, dx));
    double alpha = mod2pi(th0 - theta);
    double beta = mod2pi(th1 - theta);
    double best = std::numeric_limits<double>::infinity();
    for (Word w : {LSL, RSR, LSR, RSL, RLR, LRL}) {
        double t, p, q;
        if (word_params(w, alpha, beta, d, t, p, q)) {
            double length = r * (t + p + q);
            if (length < best) best = length;
        }
    }
    return best;
}

// --- free_heading.py ---------------------------------------------------------

struct Candidate { double length, th0, th1; };

// Turn (side) then straight to B. Returns (length, final heading); NaN if none.
inline void curve_straight(double ax, double ay, double theta0, double bx, double by,
                           int side, double r, double& length, double& heading_out) {
    double cx = ax - side * r * std::sin(theta0), cy = ay + side * r * std::cos(theta0);
    double wx = bx - cx, wy = by - cy;
    double dist = std::hypot(wx, wy);
    double gap = dist * dist - r * r;
    double straight = gap >= -1e-10 ? std::sqrt(gap > 0.0 ? gap : 0.0) : NaN;
    double heading = std::atan2(wy, wx) + side * std::atan2(r, straight);
    double arc = turn_angle(side * (heading - theta0));
    length = r * arc + straight;
    heading_out = mod2pi(heading);
}

// Turn (side), then turn the other way, ending at B (branch +1 / -1).
inline void curve_curve(double ax, double ay, double theta0, double bx, double by,
                        int side, int branch, double r, double& length, double& heading_out) {
    double c1x = ax - side * r * std::sin(theta0), c1y = ay + side * r * std::cos(theta0);
    double vx = bx - c1x, vy = by - c1y;
    double dist = std::hypot(vx, vy);
    bool feasible = (dist >= r - 1e-12) && (dist <= 3 * r + 1e-12);
    double dist_safe = dist > 0 ? dist : 1.0;
    double ux = vx / dist_safe, uy = vy / dist_safe;
    double along = (3 * r * r + dist * dist) / (2 * dist_safe);
    double across_sq = 4 * r * r - along * along;
    double across = std::sqrt(across_sq > 0.0 ? across_sq : 0.0);
    double c2x = c1x + along * ux - branch * across * uy;
    double c2y = c1y + along * uy + branch * across * ux;

    double jx = 0.5 * (c1x + c2x), jy = 0.5 * (c1y + c2y);  // junction point
    double heading_j = std::atan2(jy - c1y, jx - c1x) + side * HALF_PI;
    double heading_b = std::atan2(by - c2y, bx - c2x) - side * HALF_PI;
    double arc1 = turn_angle(side * (heading_j - theta0));
    double arc2 = turn_angle(-side * (heading_b - heading_j));
    length = feasible ? r * (arc1 + arc2) : NaN;
    heading_out = mod2pi(heading_b);
}

// paths_to_point: the 6 families LS, RS, LR+, LR-, RL+, RL-.
inline void paths_to_point(double ax, double ay, double theta0, double bx, double by, double r,
                           double length[6], double heading[6]) {
    curve_straight(ax, ay, theta0, bx, by, LEFT, r, length[0], heading[0]);
    curve_straight(ax, ay, theta0, bx, by, RIGHT, r, length[1], heading[1]);
    curve_curve(ax, ay, theta0, bx, by, LEFT, +1, r, length[2], heading[2]);
    curve_curve(ax, ay, theta0, bx, by, LEFT, -1, r, length[3], heading[3]);
    curve_curve(ax, ay, theta0, bx, by, RIGHT, +1, r, length[4], heading[4]);
    curve_curve(ax, ay, theta0, bx, by, RIGHT, -1, r, length[5], heading[5]);
}

// _line_and_arcs: S then up to 4 single arcs. Returns the number written.
int line_and_arcs(double ax, double ay, double bx, double by, double r, Candidate out[5]) {
    double dx = bx - ax, dy = by - ay;
    double chord = std::hypot(dx, dy);
    double line_heading = mod2pi(std::atan2(dy, dx));
    out[0] = {chord, line_heading, line_heading};  // S
    int n = 1;
    double half = 0.5 * chord;
    if (0 < half && half <= r) {
        double off_sq = r * r - half * half;
        double offset = std::sqrt(off_sq > 0.0 ? off_sq : 0.0);
        double mx = ax + 0.5 * dx, my = ay + 0.5 * dy;
        double nx = -dy / chord, ny = dx / chord;
        for (int sign : {+1, -1}) {
            double cx = mx + sign * offset * nx, cy = my + sign * offset * ny;
            for (int side : {LEFT, RIGHT}) {
                double head_a = std::atan2(ay - cy, ax - cx) + side * HALF_PI;
                double head_b = std::atan2(by - cy, bx - cx) + side * HALF_PI;
                double arc = turn_angle(side * (head_b - head_a));
                out[n++] = {r * arc, mod2pi(head_a), mod2pi(head_b)};
            }
        }
    }
    return n;
}

// paths_between_points: S, C, and CC switching at the midpoint of AB.
int paths_between_points(double ax, double ay, double bx, double by, double r, Candidate out[9]) {
    int n = line_and_arcs(ax, ay, bx, by, r, out);
    Candidate half[5];
    int m = line_and_arcs(ax, ay, 0.5 * (ax + bx), 0.5 * (ay + by), r, half);
    for (int k = 1; k < m; ++k)
        out[n++] = {2.0 * half[k].length, half[k].th0, half[k].th0};
    return n;
}

// --- interval.py -------------------------------------------------------------

enum Case { NONE = -1, CORNER = 0, START_FIXED = 1, END_FIXED = 2, INTERIOR = 3 };

struct Best {
    double length = 0, th0 = 0, th1 = 0;
    int which = NONE;
    void offer(double len, double a0, double a1, int c) {
        if (!std::isfinite(len)) return;
        if (which == NONE || len < length) {
            length = len; th0 = mod2pi(a0); th1 = mod2pi(a1); which = c;
        }
    }
};

Best solve(double ax, double ay, Interval i0, double bx, double by, Interval i1, double r) {
    Best best;
    double bd0[2], bd1[2];
    int n0 = i0.boundaries(bd0), n1 = i1.boundaries(bd1);
    double length[6], heading[6];

    for (int u = 0; u < n0; ++u)            // corners
        for (int v = 0; v < n1; ++v)
            best.offer(shortest_length(ax, ay, bd0[u], bx, by, bd1[v], r), bd0[u], bd1[v], CORNER);

    for (int u = 0; u < n0; ++u) {          // th0 on a boundary, th1 free
        paths_to_point(ax, ay, bd0[u], bx, by, r, length, heading);
        for (int f = 0; f < 6; ++f)
            if (std::isfinite(length[f]) && i1.contains(heading[f]))
                best.offer(length[f], bd0[u], heading[f], START_FIXED);
    }

    for (int v = 0; v < n1; ++v) {          // th1 on a boundary, th0 free: fly backwards
        paths_to_point(bx, by, bd1[v] + PI, ax, ay, r, length, heading);
        for (int f = 0; f < 6; ++f) {
            double th0 = heading[f] + PI;
            if (std::isfinite(length[f]) && i0.contains(th0))
                best.offer(length[f], th0, bd1[v], END_FIXED);
        }
    }

    Candidate cands[9];                     // both free
    int n = paths_between_points(ax, ay, bx, by, r, cands);
    for (int k = 0; k < n; ++k)
        if (i0.contains(cands[k].th0) && i1.contains(cands[k].th1))
            best.offer(cands[k].length, cands[k].th0, cands[k].th1, INTERIOR);
    return best;
}

// --- Python binding ----------------------------------------------------------

// solve(ax, ay, start0, width0, bx, by, start1, width1, r)
//   -> (length, depart_heading, arrive_heading, case)
PyObject* py_solve(PyObject*, PyObject* const* args, Py_ssize_t nargs) {
    if (nargs != 9) {
        PyErr_SetString(PyExc_TypeError,
                        "solve(ax, ay, start0, width0, bx, by, start1, width1, r) takes 9 arguments");
        return nullptr;
    }
    double v[9];
    for (int k = 0; k < 9; ++k) {
        v[k] = PyFloat_AsDouble(args[k]);
        if (v[k] == -1.0 && PyErr_Occurred()) return nullptr;
    }
    Best best = solve(v[0], v[1], {v[2], v[3]}, v[4], v[5], {v[6], v[7]}, v[8]);
    if (best.which == NONE) {  // cannot happen: corners always give a path
        PyErr_SetString(PyExc_RuntimeError, "no Dubins path found");
        return nullptr;
    }
    return Py_BuildValue("(dddi)", best.length, best.th0, best.th1, best.which);
}

PyMethodDef methods[] = {
    {"solve", reinterpret_cast<PyCFunction>(reinterpret_cast<void (*)()>(py_solve)), METH_FASTCALL,
     "solve(ax, ay, start0, width0, bx, by, start1, width1, r) -> "
     "(length, depart_heading, arrive_heading, case)"},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {PyModuleDef_HEAD_INIT, "_interval",
                      "Two-point Dubins interval problem (C++ port of interval.py).", -1, methods};

}  // namespace

PyMODINIT_FUNC PyInit__interval() { return PyModule_Create(&module); }
