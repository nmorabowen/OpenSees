"""WP-138: the wall table and the cost per unit s/B, from the committed run records.

    python wall_table.py [LEG ...]      (default: the Esmeralda arms + the DP control)

Reads only runs/<leg>/steps.csv, runs/<leg>/summary.json and runs/<leg>/logs/log.log.

* Refusal counts by code are the sum of the per-step census lines
  'refusals step N: k by code {...}' in the log. Those lines are written for every
  CONVERGED step (they cover the failed attempts before it too); the final ladder at
  the floor never converges, so its refusals are not in the sum (see the
  floor-ladder census in esmeralda/analysis/tables/walls_summary.json).
* ModifiedEuler (IntScheme 1) has no refusal path: it force-accepts at dt_min and
  caps at -maxSubsteps. For it the columns are capHits (steps.csv cap_step) and
  forcedAtDTmin (steps.csv dtmin_step).
* 'first NonPosH' = s/B of the first converged step whose census contains
  loadingNonPosH.
* Wall clock is the leg's own wall_total_s; each Esmeralda arm had a node to itself,
  MKL/OMP = 1 thread.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEGS = ["E_A", "E_B", "E_D", "E_C2", "E_B16", "ctrl_dp38"]
REF = re.compile(r"refusals step (\d+): (\d+) by code (\{[^}]*\})")
BINS = [(0.0, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.051)]


def load(leg):
    d = os.path.join(HERE, "runs", leg)
    S = list(csv.DictReader(open(os.path.join(d, "steps.csv"), newline="")))
    sm = json.load(open(os.path.join(d, "summary.json")))
    codes, first = {}, None
    lp = os.path.join(d, "logs", "log.log")
    if os.path.exists(lp):
        for m in REF.finditer(open(lp, errors="replace").read()):
            for k, v in eval(m.group(3)).items():      # noqa: S307 (our own log)
                codes[k] = codes.get(k, 0) + v
            if first is None and "loadingNonPosH" in m.group(3):
                first = int(m.group(1))
    return S, sm, codes, first


def at(S, key, lo, hi):
    """sum of `key` over steps with lo < s/B <= hi, and the s/B actually covered."""
    tot, s0, s1 = 0.0, None, None
    prev = 0.0
    for r in S:
        s = float(r["s_over_B"])
        if lo < s <= hi:
            tot += float(r[key])
            s0 = prev if s0 is None else s0
            s1 = s
        prev = s
    return tot, (s1 - s0) if s0 is not None else 0.0


def main():
    legs = sys.argv[1:] or LEGS
    print("| leg | scheme / TanType / TolR | mode | s/B end | q end (kPa) | steps | failed attempts "
          "| first NonPosH s/B (step) | refusals by code (converged-step census) "
          "| capHits / forcedAtDTmin | push wall (h) | h per 0.01 s/B |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    cost = {}
    for leg in legs:
        S, sm, codes, first = load(leg)
        a = sm.get("args", {})
        s_end, q_end = float(S[-1]["s_over_B"]), float(S[-1]["q_kPa"])
        wall_h = float(S[-1]["wall_total_s"]) / 3600.0
        cap = sum(int(r.get("cap_step", 0) or 0) for r in S)
        dtm = sum(int(r.get("dtmin_step", 0) or 0) for r in S)
        fs = f"{float(S[first - 1]['s_over_B']):.5f} ({first})" if first else "none"
        sch = (f"{a.get('scheme')} / {a.get('tantype')} / {a.get('tolr')}"
               + (f", mesh {a.get('mesh')}" if a.get("mesh", "b8") != "b8" else "")
               + (f", maxSub {a.get('maxsub')}" if a.get("maxsub", 2000) != 2000 else "")
               + (f", Krylov tol x{a.get('krylov_tol_mult'):g}"
                  if a.get("krylov_tol_mult", 10) != 10 else "")
               ) if a.get("mat") == "sanisand" else "DruckerPrager 38 deg"
        ref = ", ".join(f"{k} {v}" for k, v in sorted(codes.items(), key=lambda kv: -kv[1])) or "-"
        capcol = f"{cap} / {dtm}" if a.get("scheme") == 1 else "-"
        print(f"| {leg} | {sch} | {sm.get('mode', '?')} | {s_end:.5f} | {q_end:.1f} | {len(S)} "
              f"| {sm.get('failed_attempts', '?')} | {fs} | {ref} | {capcol} | {wall_h:.2f} "
              f"| {wall_h / (s_end / 0.01):.2f} |")
        cost[leg] = S
    print()
    print("Wall clock (h) / substeps (1e9) per 0.01 s/B, by s/B interval "
          "(blank = under 0.0005 of the interval reached; a partial interval is prorated):")
    print()
    print("| leg | " + " | ".join(f"{lo:.2f}-{min(hi, 0.05):.2f}" for lo, hi in BINS) + " |")
    print("|---|" + "---|" * len(BINS))
    for leg, S in cost.items():
        prev_wall = 0.0
        for r in S:
            w = float(r["wall_total_s"])
            r["_wstep"] = w - prev_wall
            prev_wall = w
        cells = []
        for lo, hi in BINS:
            w, ds = at(S, "_wstep", lo, hi)
            sub, _ = at(S, "sub_step_total", lo, hi) if "sub_step_total" in S[0] else (0.0, 0)
            if ds < 5e-4:
                cells.append("")
                continue
            cells.append(f"{w / 3600 / (ds / 0.01):.2f} / {sub / 1e9 / (ds / 0.01):.2f}")
        print(f"| {leg} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main()
