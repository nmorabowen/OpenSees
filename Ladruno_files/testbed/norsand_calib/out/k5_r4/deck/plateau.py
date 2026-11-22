"""WP-144 P3: the TIMs plateau rule (WP1), computed from a leg's steps.csv.

SOURCE OF THE RULE.  TIMs Workbench, Tries/2d-model/04 -- How the Strip Is Built.md, section 11 ("the plateau
criterion, WP1", read-only): over the last window of settlement ds = 0.01 B the tail tangent k_tail is measured and
normalised by the tangent of the FIRST step, k_0:  k_norm = k_tail / k_0 < eps = 0.02.  The verdict reads three facts:
the tail tangent, its trend over the previous window, and why the run ended.  A run that ends on the step floor or the
wall clock with a tail still rising is a WALL, not a limit point.  (A refused commit is a verdict of its own and is
never put through this rule.  This deck does not tell a refused commit from a failed step: its ladder halves ds on any
failed analyze, so read engine.log / log.log for refusal lines before trusting a verdict.)

WHAT THIS FILE DECIDES THAT 04 DOES NOT SAY (stated, not hidden):
  * k_tail is the ordinary-least-squares slope of q against s/B over the rows whose s/B lies in the last window
    [s_end - 0.01, s_end] (04 says "measured", not how); the secant over the same window is reported next to it.
  * k_0 is (q_1 - q_start) / (s_1 - 0): the first accepted step against the push-start state.  q_start is
    summary.json["seating"]["q0_kPa"] (the gross dead-load pressure at the start of the push) when summary.json is
    given.  Without it k_0 falls back to the secant of steps 1 -> 2 and says so ("k0_source").
  * the tangents are in kPa per unit s/B; k_norm is a ratio, so the unit cancels.
  * trend = k_tail against the slope of the previous window of the same width: falling = the tail slope is smaller than
    the previous window's (the curve is flattening), rising = larger.
  * k_norm is SIGNED.  A post-peak softening tail (k_tail < 0) satisfies k_norm < eps; the verdict then names it
    "plateau_softening" so a reader does not mistake it for a flat-topped curve.

VERDICTS: plateau | plateau_softening | not_reached (run reached its s/B target with the tail still rising: the target
was too small) | wall (run stopped on FLOOR/WALL/PROBETEST with the tail not flat) | too_short (fewer than 3 steps
in the window, or the leg is shorter than the window).

Usage:  python plateau.py <steps.csv> [<summary.json>] [--eps 0.02] [--window 0.01]      (prints one JSON object)
"""
from __future__ import annotations

import csv
import json
import sys


def _ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 0.0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def _window(sB, q, hi, w):
    idx = [k for k in range(len(sB)) if hi - w - 1e-15 <= sB[k] <= hi + 1e-15]
    return idx


def plateau_from_rows(sB, q, q_start=None, mode=None, eps=0.02, window=0.01):
    """sB, q: lists (accepted steps, ascending s/B).  Returns the verdict dict."""
    out = dict(rule="WP1 (04 sec. 11)", eps=eps, window_sB=window, steps=len(sB), mode=mode)
    if len(sB) < 3:
        out.update(verdict="too_short", reason="fewer than 3 accepted steps")
        return out
    # k_0: first step against the push start (s = 0, q = q_start), else steps 1 -> 2
    if q_start is not None and sB[0] > 0.0:
        k0 = (q[0] - q_start) / sB[0]
        out["k0_source"] = "step 1 vs push start (q_start from summary.json seating.q0_kPa)"
    else:
        k0 = (q[1] - q[0]) / (sB[1] - sB[0])
        out["k0_source"] = "secant of steps 1 -> 2 (no q_start given)"
    out["k0_kPa_per_sB"] = k0
    hi = sB[-1]
    wi = _window(sB, q, hi, window)
    out["s_end_over_B"] = hi
    out["q_end_kPa"] = q[-1]
    qmax = max(q)
    out["q_max_kPa"] = qmax
    out["s_at_q_max_over_B"] = sB[q.index(qmax)]
    if len(wi) < 3 or hi < window:
        out.update(verdict="too_short",
                   reason=f"{len(wi)} step(s) in the last window of {window} s/B (need >= 3)")
        return out
    kt = _ols([sB[k] for k in wi], [q[k] for k in wi])
    kt_sec = (q[wi[-1]] - q[wi[0]]) / (sB[wi[-1]] - sB[wi[0]])
    out.update(k_tail_kPa_per_sB=kt, k_tail_secant_kPa_per_sB=kt_sec, n_in_window=len(wi),
               k_norm=kt / k0, k_norm_secant=kt_sec / k0)
    # previous window
    wp = _window(sB, q, hi - window, window)
    wp = [k for k in wp if sB[k] < sB[wi[0]] - 1e-15]
    if len(wp) >= 3:
        kp = _ols([sB[k] for k in wp], [q[k] for k in wp])
        out.update(k_prev_kPa_per_sB=kp, k_prev_norm=kp / k0,
                   trend=("falling" if kt < kp else "rising" if kt > kp else "flat"))
    else:
        out.update(trend="n/a (previous window has < 3 steps)")
    kn = kt / k0
    ended = (mode or "TARGET")
    if kn < eps:
        out["verdict"] = "plateau_softening" if kn < -eps else "plateau"
    elif ended in ("FLOOR", "WALL", "PROBETEST"):
        out["verdict"] = "wall"
    else:
        out["verdict"] = "not_reached"
    return out


def plateau_from_files(steps_csv, summary_json=None, eps=0.02, window=0.01):
    sB, q = [], []
    with open(steps_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            sB.append(float(r["s_over_B"]))
            q.append(float(r["q_kPa"]))
    q_start, mode = None, None
    if summary_json:
        try:
            with open(summary_json) as fh:
                sm = json.load(fh)
            q_start = (sm.get("seating") or {}).get("q0_kPa")
            mode = sm.get("mode")
        except (OSError, ValueError):
            pass
    return plateau_from_rows(sB, q, q_start, mode, eps, window)


def main(argv):
    eps, window, pos = 0.02, 0.01, []
    it = iter(argv[1:])
    for a in it:
        if a == "--eps":
            eps = float(next(it))
        elif a == "--window":
            window = float(next(it))
        else:
            pos.append(a)
    if not pos:
        print(__doc__)
        return 2
    res = plateau_from_files(pos[0], pos[1] if len(pos) > 1 else None, eps, window)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
