"""ADR-92 CP1 -- read the surcharge legs against ADR-86b T8's zero-surcharge ones.

The owner took ADR-86b handoff sec.5b option **B (a small free-surface surcharge)** at
CP1. This rolls the resulting legs up next to T8's, which are the same deck, the same
engine family and the same cap with `--surcharge 0`.

    python3.12 adr92_cp1_rollup.py [--cp1 adr92_cp1] [--t8 <path to adr86b_t8>]

READING RULE, applied by this script and not left to the reader: a surcharge `Q` adds
`Q*Nq` to the COLLAPSE load and nothing like it before the mechanism forms. At this
soil's `M_c = 1.3309` (`sin phi = 3Mc/(6+Mc)` => `phi = 33.0 deg`, `Nq = 26.1`) that is
+52 kPa at `Q = 2` and +261 kPa at `Q = 10`. Every capacity column therefore prints BOTH
the raw `q` and `q - Q*Nq`; the pre-peak matched-settlement columns print raw only, with
the correction marked `n/a`, because `Nq` is a collapse quantity.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
M_C = 1.3309
SIN_PHI = 3.0 * M_C / (6.0 + M_C)
PHI = math.degrees(math.asin(SIN_PHI))
NQ = math.exp(math.pi * math.tan(math.radians(PHI))) * \
    math.tan(math.radians(45.0 + PHI / 2.0)) ** 2


def load(pattern):
    out = []
    for f in sorted(glob.glob(pattern)):
        try:
            with open(f, encoding="utf-8") as fh:
                d = json.load(fh)
        except (OSError, ValueError):
            continue
        if "tag" in d:
            d["_file"] = f
            out.append(d)
    return out


def q_at(rec, target):
    for c in rec.get("checkpoints", []):
        if abs(float(c.get("cp_target", -1)) - target) < 1e-12:
            for k in ("q_foot", "q", "q_base"):
                if k in c:
                    return abs(float(c[k]))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cp1", default=os.path.join(HERE, "adr92_cp1"))
    ap.add_argument("--t8", default=os.path.abspath(os.path.join(
        HERE, "..", "..", "..", "..", "agent-aeef311decf9449d1",
        "Ladruno_files", "testbed", "hypo_bearing", "adr86b_t8")))
    a = ap.parse_args()

    cp1 = load(os.path.join(a.cp1, "*", "a2_*.json"))
    t8 = load(os.path.join(a.t8, "*", "a2_*.json"))
    if not cp1:
        raise SystemExit(f"no legs under {a.cp1}")

    print("=" * 108)
    print("ADR-92 CP1 -- the surcharge legs (owner decision B) against ADR-86b T8's "
          "zero-surcharge legs")
    print("=" * 108)
    print(f"  phi = {PHI:.2f} deg from M_c = {M_C}; Nq = {NQ:.2f}; "
          f"Q*Nq = {2 * NQ:.0f} kPa at Q=2, {10 * NQ:.0f} kPa at Q=10")
    builds = sorted({r.get("build", "?")[:9] for r in cp1})
    print(f"  CP1 engine   : {', '.join(builds)}")
    if len(builds) > 1:
        print("  *** REFUSING to compare: CP1 legs span more than one build hash")
        return
    if t8:
        print(f"  T8 engine    : {', '.join(sorted({r.get('build','?')[:9] for r in t8}))}")
        print("  (a build difference between the two sets is a caveat on every "
              "delta below, not a defect)")

    print("\n--- TERMINATION AND DEPTH ---")
    hdr = (f"{'leg':>22} {'Q':>4} {'cap':>6} {'mode':>8} {'s/B end':>9} "
           f"{'plateau':>8} {'peaked':>7} {'tail %':>8} {'subdiv':>8} "
           f"{'clamp':>6} {'wall s':>8}")
    print(hdr); print("-" * len(hdr))
    rows = []
    for r in sorted(cp1 + t8, key=lambda x: (x.get("e_init", 0),
                                             x.get("surcharge_kpa", 0.0))):
        q0 = float(r.get("surcharge_kpa", 0.0) or 0.0)
        tail = r.get("tail_pct")
        rows.append((r, q0))
        print(f"{r['tag']:>22} {q0:>4.0f} {int(r.get('max_substeps', 0)):>6} "
              f"{str(r.get('mode', '?')):>8} {float(r.get('s_end_over_B', 0)):>9.5f} "
              f"{str(bool(r.get('plateau'))):>8} {str(bool(r.get('peaked'))):>7} "
              f"{('nan' if tail is None or tail != tail else f'{float(tail):.0f}'):>8} "
              f"{int(r.get('nsub', 0)):>3}/{int(r.get('subdiv_budget', 0)):<4} "
              f"{int(r.get('n_clamping', 0)):>6} {float(r.get('wall_s', 0)):>8.0f}")

    print("\n--- CAPACITY (only rows whose leg is admissible mean anything) ---")
    hdr = (f"{'leg':>22} {'Q':>4} {'admissible':>11} {'q_u raw':>10} "
           f"{'- Q*Nq':>10} {'s_peak/B':>9}")
    print(hdr); print("-" * len(hdr))
    for r, q0 in rows:
        qu = r.get("q_u")
        qu = abs(float(qu)) if qu is not None else float("nan")
        print(f"{r['tag']:>22} {q0:>4.0f} "
              f"{str(bool(r.get('admissible_as_capacity'))):>11} "
              f"{qu:>10.2f} {qu - q0 * NQ:>10.2f} "
              f"{float(r.get('s_peak_over_B') or float('nan')):>9.5f}")
    print("  `q_u raw` on an INADMISSIBLE leg is where the RUN stopped, not where the "
          "soil failed.\n  The `- Q*Nq` column is only meaningful on an admissible one.")

    print("\n--- MATCHED SETTLEMENT (pre-peak; Nq does NOT apply, raw q only) ---")
    cps = [0.002, 0.005, 0.010, 0.020, 0.040]
    hdr = f"{'leg':>22} {'Q':>4}" + "".join(f"{f's/B {c:g}':>12}" for c in cps)
    print(hdr); print("-" * len(hdr))
    for r, q0 in rows:
        line = f"{r['tag']:>22} {q0:>4.0f}"
        for c in cps:
            v = q_at(r, c)
            line += f"{(f'{v:.2f}' if v is not None else '-'):>12}"
        print(line)

    print("\n--- LOCALIZATION WIDTH w2 (m) at the two probe depths ---")
    hdr = f"{'leg':>22} {'Q':>4} {'w2 z=-0.625':>13} {'w2 z=-1.375':>13} {'yield ele':>10}"
    print(hdr); print("-" * len(hdr))
    for r, q0 in rows:
        print(f"{r['tag']:>22} {q0:>4.0f} {float(r.get('w2_z-0.625', float('nan'))):>13.4f} "
              f"{float(r.get('w2_z-1.375', float('nan'))):>13.4f} "
              f"{int(r.get('n_yield_ele', 0)):>10}")


if __name__ == "__main__":
    main()
