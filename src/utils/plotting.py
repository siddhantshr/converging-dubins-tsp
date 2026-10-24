"""Plot the tours behind each bound."""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")  # write to file, no window needed
import matplotlib.pyplot as plt

from ..dubins.path import sample, shortest_path
from .report import describe_options


def _arrow(ax, x, y, heading, length, **style):
    ax.annotate("", xy=(x + length * math.cos(heading), y + length * math.sin(heading)),
                xytext=(x, y), arrowprops=dict(arrowstyle="->", lw=1.4, **style))


def plot_results(points, results, r: float, filename: str, suptitle: str = ""):
    fig, axes = plt.subplots(1, len(results), figsize=(6.5 * len(results), 6), squeeze=False)
    xs, ys = zip(*points)
    span = max(max(xs) - min(xs), max(ys) - min(ys), 4 * r)
    arrow_len = 0.08 * span

    for ax, res in zip(axes[0], results):
        for leg in res.legs:
            a, b = points[leg.from_target], points[leg.to_target]
            path = shortest_path((*a, leg.depart_heading), (*b, leg.arrive_heading), r)
            pts = sample(path, step=r / 20)
            ax.plot(pts[:, 0], pts[:, 1], color="tab:blue", lw=1.6)
            _arrow(ax, *a, leg.depart_heading, arrow_len, color="tab:green")
            if res.kind == "lower":
                # Arrival heading can differ from the next departure: show both.
                _arrow(ax, *b, leg.arrive_heading, arrow_len, color="tab:red", linestyle="--")

        ax.scatter(xs, ys, color="black", zorder=3)
        for idx, (x, y) in enumerate(points):
            ax.annotate(str(idx), (x, y), textcoords="offset points", xytext=(6, 6))

        label = "Upper bound" if res.kind == "upper" else "Lower bound"
        flag = "" if res.is_valid_bound else " (estimate)"
        ax.set_title(f"{label}{flag}: {res.value:.4f}\n"
                     f"{describe_options(res)}, solver = {res.solver}")
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(alpha=0.3)

    if any(res.kind == "lower" for res in results):
        fig.text(0.5, 0.01, "green = leaving heading, red dashed = arriving heading "
                 "(they may differ in the lower bound)", ha="center", fontsize=9)
    if suptitle:
        fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def plot_convergence(iterations, filename: str):
    """Bounds and gap against the iteration number."""
    idx = [it.index for it in iterations]
    lower = [it.lower.value for it in iterations]
    upper = [it.upper.value for it in iterations]
    best_lower = [max(lower[: k + 1]) for k in range(len(lower))]
    best_upper = [min(upper[: k + 1]) for k in range(len(upper))]
    gap = [(u - l) / l * 100 for u, l in zip(best_upper, best_lower)]

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    top.plot(idx, upper, "o-", color="tab:blue", label="upper bound")
    top.plot(idx, lower, "s-", color="tab:orange", label="lower bound")
    top.fill_between(idx, best_lower, best_upper, color="grey", alpha=0.15,
                     label="optimum lies here")
    top.set_ylabel("tour length")
    top.legend()
    top.grid(alpha=0.3)

    bottom.plot(idx, gap, "o-", color="tab:red")
    for x, y in zip(idx, gap):
        bottom.annotate(f"{y:.2f}%", (x, y), textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=8)
    bottom.margins(y=0.15)  # room for the labels
    bottom.set_ylabel("gap (best UB - best LB) / best LB  [%]")
    bottom.set_xlabel("iteration")
    bottom.set_xticks(idx)
    bottom.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)
