"""ADR-92 P1 gate 1 -- the BVP gate, stated so it can FAIL.

CP1 measured where the wall clock goes on the real footing deck
(`_adr92_cp1_surcharge_results` §6): only ~30-35 steps per leg converge on plain
Newton, 61-83 % need a fallback rung, and **89-93 % of all Newton iterations are
spent in rungs that fail and are discarded**. IMPL-EX's structural claim is that a
frozen SPD operator makes the step LINEAR, so there is no ladder to fail.

    PREDICTION: with `-implex`, `past rung 1` falls from 61 % to near zero and the
    failed-rung share of iterations falls from 89 % to single digits.

    IF THE LADDER STILL FIRES AT CP1's RATE, the cost case is REFUTED and P1 closes
    as a correctness-only feature. This script prints the table either way and says
    which happened. It is not allowed to be silent about a refutation.

Usage -- it does not run anything itself, it READS leg records:

    python3.12 adr92_bvp_gate.py --implex <dir> [--baseline <dir> ...]

`<dir>` is any directory holding `a2_*.json` leg records (recursively). The CP1
baselines live in the 92b worktree at
`Ladruno_files/testbed/hypo_bearing/adr92_cp1/{g0.6944_q10,...}`.

THE DECOMPOSITION, and why it is recoverable. The push ladder is
`Newton`(25) -> `NewtonLineSearch`(40) -> `KrylovNewton`(60, relaxed). `nfail`
counts failed RUNGS, so: a step converging on rung 1 adds 0, on rung 2 adds 1, on
rung 3 adds 2, and a step that fails all three adds 3 and subdivides. Hence

    rung3 = nrelax                       (the relaxed rung is rung 3 by construction)
    rung2 = nfail - 3*nsub - 2*nrelax
    rung1 = steps - rung2 - rung3

The iteration split costs each FAILED rung at its cap -- an UPPER bound, since a
rung that diverges early costs less -- and a converged step at `CONVERGED_ITERS`.
The ratio is insensitive to that last number; it is stated, not hidden.
"""
from __future__ import annotations

import argparse
import glob
import json
import os

RUNG_CAPS = (25, 40, 60)          # Newton, NewtonLineSearch, KrylovNewton
CONVERGED_ITERS = 5               # assumed cost of a step that converges; see docstring
REFUTE_PAST_RUNG1 = 40.0          # % -- at or above CP1's band is a refutation
PASS_PAST_RUNG1 = 10.0            # % -- the prediction is "near zero"
PASS_FAILED_SHARE = 10.0          # % -- "single digits"


def load(root):
    out = []
    for f in sorted(glob.glob(os.path.join(root, "**", "a2_*.json"), recursive=True)):
        try:
            with open(f, encoding="utf-8") as fh:
                d = json.load(fh)
        except (OSError, ValueError):
            continue
        if "tag" in d and d.get("steps"):
            d["_dir"] = os.path.basename(os.path.dirname(f))
            out.append(d)
    return out


def decompose(d):
    st, nf, ns, nr = d["steps"], d["nfail"], d["nsub"], d["nrelax"]
    rung3 = nr
    rung2 = nf - 3 * ns - 2 * nr
    rung1 = st - rung2 - rung3
    fail_it = ns * sum(RUNG_CAPS) + rung2 * RUNG_CAPS[0] + rung3 * (RUNG_CAPS[0] + RUNG_CAPS[1])
    ok_it = st * CONVERGED_ITERS
    return dict(steps=st, rung1=rung1, rung2=rung2, rung3=rung3, nsub=ns,
                past1=(st - rung1) / st * 100.0 if st else float("nan"),
                failshare=fail_it / (fail_it + ok_it) * 100.0 if (fail_it + ok_it) else 0.0)


def table(title, legs):
    print(f"\n--- {title} ---")
    hdr = (f"{'leg':>26}{'implex':>8}{'steps':>7}{'rung1':>7}{'rung2':>7}{'rung3':>7}"
           f"{'nsub':>6}{'past rung1':>12}{'failed-rung iters':>19}")
    print(hdr); print("-" * len(hdr))
    rows = []
    for d in legs:
        k = decompose(d)
        rows.append(k)
        flag = "yes" if d.get("implex") else ("?" if "implex" not in d else "no")
        print(f"{d['_dir']:>26}{flag:>8}{k['steps']:>7}{k['rung1']:>7}{k['rung2']:>7}"
              f"{k['rung3']:>7}{k['nsub']:>6}{k['past1']:>11.1f}%{k['failshare']:>18.1f}%")
    if rows:
        p = sum(r["past1"] for r in rows) / len(rows)
        f = sum(r["failshare"] for r in rows) / len(rows)
        print(f"{'MEAN':>26}{'':>8}{'':>7}{'':>7}{'':>7}{'':>7}{'':>6}{p:>11.1f}%{f:>18.1f}%")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--implex", required=True, help="dir of leg records run WITH -implex")
    ap.add_argument("--baseline", action="append", default=[],
                    help="dir(s) of CP1 leg records run WITHOUT it; repeatable")
    a = ap.parse_args()

    impl = load(a.implex)
    if not impl:
        raise SystemExit(f"no leg records under {a.implex}")
    base = [d for r in a.baseline for d in load(r)]

    print("=" * 104)
    print("ADR-92 P1 GATE 1 -- the BVP gate: does a linear step remove the ladder?")
    print("=" * 104)
    print(f"  rung caps {RUNG_CAPS}, converged step costed at {CONVERGED_ITERS} iterations")
    print(f"  builds: implex {sorted({d.get('build','?')[:9] for d in impl})}"
          + (f", baseline {sorted({d.get('build','?')[:9] for d in base})}" if base else ""))

    brows = table("BASELINE (CP1, no -implex)", base) if base else []
    irows = table("WITH -implex", impl)

    ip = sum(r["past1"] for r in irows) / len(irows)
    ifs = sum(r["failshare"] for r in irows) / len(irows)

    print("\n" + "=" * 104)
    if brows:
        bp = sum(r["past1"] for r in brows) / len(brows)
        bfs = sum(r["failshare"] for r in brows) / len(brows)
        print(f"  past rung 1        {bp:6.1f} %  ->  {ip:6.1f} %")
        print(f"  failed-rung iters  {bfs:6.1f} %  ->  {ifs:6.1f} %")
    if ip <= PASS_PAST_RUNG1 and ifs <= PASS_FAILED_SHARE:
        print("\n  VERDICT: PREDICTION MET. The ladder is gone -- the step is effectively")
        print("  linear, which is IMPL-EX's structural claim measured on the real deck.")
    elif ip >= REFUTE_PAST_RUNG1:
        print("\n  VERDICT: **REFUTED**. The ladder still fires at or near CP1's rate with")
        print("  -implex on. IMPL-EX's COST case does not hold on this deck: P1 closes as a")
        print("  correctness-only feature and the ADR §2 item 1 must be rewritten to say so.")
        print("  Do NOT reach for a tuning knob before reporting this number.")
    else:
        print("\n  VERDICT: PARTIAL -- between the prediction and a refutation. Report the")
        print("  table as measured; do not round it toward either.")
    print("=" * 104)


if __name__ == "__main__":
    main()
