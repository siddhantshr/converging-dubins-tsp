"""Print and save results."""

from __future__ import annotations

import json
import math
from dataclasses import asdict


def _deg(rad: float) -> str:
    return f"{math.degrees(rad) % 360:6.1f}°"


def describe_options(res) -> str:
    """E.g. '8 headings per target' or '8-12 headings per target (56 total)'."""
    what = "headings" if res.kind == "upper" else "intervals"
    counts = res.options_per_target
    if min(counts) == max(counts):
        return f"{counts[0]} {what} per target"
    return f"{min(counts)}-{max(counts)} {what} per target ({sum(counts)} total)"


def format_result(res) -> str:
    title = "UPPER BOUND" if res.kind == "upper" else "LOWER BOUND"
    status = "valid bound" if res.is_valid_bound else "ESTIMATE ONLY"
    lines = [
        f"{title}: {res.value:.6f}   ({status})",
        f"  {describe_options(res)}, solver = {res.solver}, "
        f"build {res.build_seconds:.2f}s, solve {res.solve_seconds:.2f}s",
        f"  order: {' -> '.join(map(str, res.order + res.order[:1]))}",
        "  legs:   from -> to   leave     arrive    length",
    ]
    for leg in res.legs:
        lines.append(f"          {leg.from_target:>4} -> {leg.to_target:<4} "
                     f"{_deg(leg.depart_heading)}  {_deg(leg.arrive_heading)}  {leg.length:9.4f}")
    lines += [f"  note: {note}" for note in res.notes]
    return "\n".join(lines)


def format_gap(upper, lower) -> str:
    gap = (upper.value - lower.value) / lower.value * 100 if lower.value > 0 else float("inf")
    caveat = "" if lower.is_valid_bound else "  (lower bound is only an estimate)"
    return (f"{lower.value:.6f} <= optimal Dubins tour <= {upper.value:.6f}\n"
            f"gap = (UB - LB) / LB = {gap:.2f}%{caveat}")


def save_json(results, filename: str):
    with open(filename, "w") as fh:
        json.dump([asdict(res) for res in results], fh, indent=2)


# --- refinement output ------------------------------------------------------------

def _deg_list(headings):
    return [round(math.degrees(h) % 360, 4) for h in headings]


def format_iteration(it) -> str:
    lines = [f"ITERATION {it.index}", ""]
    if it.index > 0:
        lines.append("headings added (deg) per target, inside the interval the previous "
                     "lower bound picked:")
        for t, fresh in enumerate(it.added):
            lines.append(f"  target {t}: {', '.join(f'{d:.2f}' for d in _deg_list(fresh)) or '-'}")
        lines.append("")
    lines += [format_result(it.upper), "", format_result(it.lower), "",
              f"this iteration:  {it.lower.value:.6f} <= optimum <= {it.upper.value:.6f}",
              f"best so far:     {it.best_lower:.6f} <= optimum <= {it.best_upper:.6f}",
              f"gap (best)     = {it.gap_percent:.3f}%"]
    return "\n".join(lines)


def iteration_record(it) -> dict:
    return {
        "iteration": it.index,
        "upper_bound": asdict(it.upper),
        "lower_bound": asdict(it.lower),
        "best_upper": it.best_upper,
        "best_lower": it.best_lower,
        "gap_percent": it.gap_percent,
        "headings_deg_per_target": [_deg_list(h) for h in it.headings],
        "intervals_deg_per_target": [[str(iv) for iv in ivs] for ivs in it.intervals],
        "added_headings_deg_per_target": [_deg_list(a) for a in it.added],
    }


def format_convergence(iterations) -> str:
    header = (f"{'iter':>4}  {'nodes UB':>8}  {'nodes LB':>8}  {'upper':>11}  {'lower':>11}  "
              f"{'best UB':>11}  {'best LB':>11}  {'gap %':>7}  {'time s':>7}")
    lines = [header, "-" * len(header)]
    for it in iterations:
        secs = sum(r.build_seconds + r.solve_seconds for r in (it.upper, it.lower))
        lines.append(
            f"{it.index:>4}  {sum(it.upper.options_per_target):>8}  "
            f"{sum(it.lower.options_per_target):>8}  {it.upper.value:>11.6f}  "
            f"{it.lower.value:>11.6f}  {it.best_upper:>11.6f}  {it.best_lower:>11.6f}  "
            f"{it.gap_percent:>7.3f}  {secs:>7.1f}")
    if not all(it.lower.is_valid_bound for it in iterations):
        lines.append("\nnote: lower bounds were solved with LKH, so they are estimates only.")
    return "\n".join(lines)


def format_timings(iterations) -> str:
    """Seconds per iteration, split by stage, for both bounds."""
    ub_cols = [("dubins", "dubins"), ("transform", "NB"), ("lkh", "lkh"), ("other", "other")]
    lb_cols = [("interval", "interval"), ("transform", "NB"), ("lkh", "lkh"),
               ("exact", "exact"), ("other", "other")]

    def total(res):
        return res.build_seconds + res.solve_seconds

    def row(label, ub, lb, solves, rounds):
        cells = [f"{label:>5}"]
        cells += [f"{ub[k]:>8.2f}" for k, _ in ub_cols] + [f"{ub['total']:>8.2f}", " |"]
        cells += [f"{lb[k]:>8.2f}" for k, _ in lb_cols] + [f"{lb['total']:>8.2f}"]
        cells += [f"{solves:>7}", f"{rounds:>6}", " |", f"{ub['total'] + lb['total']:>8.2f}"]
        return " ".join(cells)

    head = " ".join([f"{'iter':>5}"] + [f"{h:>8}" for _, h in ub_cols] + [f"{'total':>8}", " |"]
                    + [f"{h:>8}" for _, h in lb_cols] + [f"{'total':>8}", f"{'#solves':>7}",
                    f"{'rounds':>6}", " |", f"{'iter':>8}"])
    ub_width = 6 + 9 * (len(ub_cols) + 1)
    group = (f"{'':6}{' UPPER BOUND ':-^{ub_width - 6}}  |"
             f"{' LOWER BOUND ':-^{9 * (len(lb_cols) + 1) + 15}}  |")
    lines = ["TIME PER ITERATION (seconds)", "", group, head, "-" * len(head)]

    sums_ub, sums_lb, sum_solves, sum_rounds = {}, {}, 0, 0
    for it in iterations:
        ub = {k: it.upper.timings.get(k, 0.0) for k, _ in ub_cols} | {"total": total(it.upper)}
        lb = {k: it.lower.timings.get(k, 0.0) for k, _ in lb_cols} | {"total": total(it.lower)}
        lines.append(row(it.index, ub, lb, it.lower.interval_solves, it.lower.solve_rounds))
        for sums, vals in ((sums_ub, ub), (sums_lb, lb)):
            for k, v in vals.items():
                sums[k] = sums.get(k, 0.0) + v
        sum_solves += it.lower.interval_solves
        sum_rounds += it.lower.solve_rounds
    lines += ["-" * len(head), row("total", sums_ub, sums_lb, sum_solves, sum_rounds)]

    grand = sums_ub["total"] + sums_lb["total"]
    if grand > 0:
        share = [f"UB {h} {sums_ub[k] / grand:.1%}" for k, h in ub_cols] + \
                [f"LB {h} {sums_lb[k] / grand:.1%}" for k, h in lb_cols]
        lines += ["", "share of total time: " + ", ".join(share)]

    lines += [
        "",
        "columns:",
        "  dubins    building the upper-bound graph (point-to-point Dubins lengths, cached)",
        "  interval  time inside the interval-problem solver (graph build + lazy refreshes)",
        "  NB        integer cost scaling + Noon-Bean transform (+ decode), all solve rounds",
        "  lkh       LKH solve (lower bound: all rounds if solver = lkh, screening rounds if exact)",
        "  exact     exact ATSP solve (CBC), incl. building the MILP and subtour cuts",
        "  other     everything else (lazy cost estimates from parents, bookkeeping)",
        "  #solves   interval problems actually solved in that iteration",
        "  rounds    lower-bound solve rounds (one per lazy re-solve)",
        "  iter      upper total + lower total (plotting and file output not included)",
    ]
    return "\n".join(lines)


def save_text(text: str, filename) -> None:
    with open(filename, "w") as fh:
        fh.write(text + "\n")


def save_records(records, filename) -> None:
    with open(filename, "w") as fh:
        json.dump(records, fh, indent=2)
