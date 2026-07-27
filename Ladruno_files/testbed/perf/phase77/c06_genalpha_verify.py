"""ADR-77 C0-6 — acceptance test for the GeneralizedAlpha::update() fix.

The fix passes the alphaM-weighted acceleration (Ualphadotdot) to setResponse
instead of the unweighted Udotdot. Two things must be true afterwards, and they
pull in opposite directions, which is what makes this a real test:

  REGRESSION GUARD: at alphaM == 1.0 the result must be BIT-IDENTICAL, because
  there Ualphadotdot == Udotdot identically and the fix is a no-op by
  construction. If this moves, the change did something it should not have.

  CORRECTION: at every other alphaM the result must CHANGE, and change toward
  convergence -- a positive observed order instead of the pre-fix negative one
  (error growing as dt shrinks).

Also checks LadrunoGeneralizedAlpha, which inherits update() and so must track
GeneralizedAlpha exactly, before and after.

Usage:
  ... python3.12 -S c06_genalpha_verify.py before   # writes c06_baseline.json
  ... python3.12 -S c06_genalpha_verify.py after    # compares against it
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "1")   # 1 thread: bit-exact oracle (P1f)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees"
    r"\.claude\worktrees\pardisio-profiling-0a03b1\dist\bin",
)
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import t3_scheme_shootout as t3  # noqa: E402  (reuses the T3 deck verbatim)

MODE = sys.argv[1] if len(sys.argv) > 1 else "before"
BASELINE = HERE / "c06_baseline.json"

t3.N = 5
DTS = [2.0e-5, 1.0e-5]          # two dt so an observed ORDER can be formed
CASES = [
    ("GeneralizedAlpha", 1.0, 2.0 / 3.0, "REGRESSION GUARD (alphaM=1.0, rho=0.5)"),
    ("GeneralizedAlpha", 0.8, 0.6, "corrected (rho=0.667)"),
    ("GeneralizedAlpha", 0.65, 0.55, "corrected (rho=0.818)"),
    ("GeneralizedAlpha", 0.5, 0.5, "corrected (rho=1.0, no damping)"),
    ("LadrunoGeneralizedAlpha", 1.0, 2.0 / 3.0, "fork class, guard"),
    ("LadrunoGeneralizedAlpha", 0.65, 0.55, "fork class, corrected"),
]


def main():
    print(f"ADR-77 C0-6 verify [{MODE}] — 1 thread, n={t3.N}, T3 ramp deck")
    _, ref, _, _ = t3.run(1.25e-6, ["Newmark", 0.5, 0.25])

    out = {}
    for cls, aM, aF, note in CASES:
        rec = {"note": note, "dts": {}}
        for dt in DTS:
            _, h, it, ns = t3.run(dt, [cls, aM, aF])
            rec["dts"][f"{dt:g}"] = {
                "uz_end": h[-1], "l2": t3.l2_err(h, ref), "it_per_step": it / ns,
            }
        e = [rec["dts"][f"{d:g}"]["l2"] for d in DTS]
        rec["order"] = (math.log(e[0] / e[1]) / math.log(DTS[0] / DTS[1])
                        if e[0] > 0 and e[1] > 0 else None)
        out[f"{cls}({aM:g},{aF:g})"] = rec
        o = rec["order"]
        print(f"  {cls}({aM:g},{aF:g})  {note}")
        for dt in DTS:
            d = rec["dts"][f"{dt:g}"]
            print(f"      dt={dt:g}  uz_end={d['uz_end']:.17e}  L2={d['l2']:.4e}"
                  f"  it/step={d['it_per_step']:.2f}")
        print(f"      observed order = {o:.2f}" if o else "      order n/a")

    if MODE == "before":
        BASELINE.write_text(json.dumps(out, indent=2))
        print(f"\nwrote {BASELINE.name} (pre-fix baseline)")
        return

    if not BASELINE.exists():
        print("\nno baseline to compare against")
        return
    base = json.loads(BASELINE.read_text())
    print("\n=== VERDICT ===")
    fail = 0
    for key, rec in out.items():
        b = base.get(key)
        if not b:
            continue
        guard = "alphaM=1" in rec["note"] or rec["note"].endswith("guard")
        same = all(rec["dts"][d]["uz_end"] == b["dts"][d]["uz_end"] for d in rec["dts"])
        if guard:
            ok = same
            print(f"  {'PASS' if ok else 'FAIL'}  {key:<42} bit-identical "
                  f"(required) -> {same}")
        else:
            improved = (rec["order"] or -9) > 0 and not same
            ok = improved
            print(f"  {'PASS' if ok else 'FAIL'}  {key:<42} changed={not same} "
                  f"order {b['order']:.2f} -> {rec['order']:.2f} (must be > 0)")
        fail += (0 if ok else 1)
    # the fork class must track the base class exactly
    for a, bkey in (("GeneralizedAlpha(1,0.666667)", "LadrunoGeneralizedAlpha(1,0.666667)"),
                    ("GeneralizedAlpha(0.65,0.55)", "LadrunoGeneralizedAlpha(0.65,0.55)")):
        if a in out and bkey in out:
            same = all(out[a]["dts"][d]["uz_end"] == out[bkey]["dts"][d]["uz_end"]
                       for d in out[a]["dts"])
            print(f"  {'PASS' if same else 'FAIL'}  fork tracks base: {a} == {bkey} -> {same}")
            fail += (0 if same else 1)
    print(f"\n{'ALL CHECKS PASSED' if not fail else f'{fail} CHECK(S) FAILED'}")


if __name__ == "__main__":
    main()
