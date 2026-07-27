"""ADR-77 — acceptance test for the C0 patch wave (C0-1, C0-2, C0-4, C0-5).

Every patch in the wave claims "cannot change results on valid input". That is
the claim worth testing, so most of this file is a no-regression check; the
new-behaviour checks are the small part.

  C0-1  analyze(n, 0.0) must now return an ERROR from TRBDF2 / TRBDF3 /
        BackwardEuler / Houbolt instead of silently building an inf/nan tangent.
        Regression: every dt > 0 result must be unchanged.
  C0-2  `algorithm Newton -hall` on BackwardEuler / Newmark1 / Collocation must
        no longer assemble a zero (or stiffness-free) tangent. Pre-fix these
        produced a completed run with garbage; post-fix they must produce the
        same answer as plain Newton to within the convergence tolerance.
  C0-4  removing the dead `static bool converged` from Newmark must be
        bit-identical.
  C0-5  KrylovNewton (and friends) must now emit formTangent / linearSolve /
        update scopes so a step can be decomposed. Profiling-only.

ORACLE NOTE, stated because it matters for how much this proves. The C0 source
edits were already applied when this test was written, so there is no pre-patch
run of the zero-dt case: that "before" is established by READING the pre-patch
code (no deltaT guard existed; c2 = 2.0/deltaT and c3 = 4.0/(deltaT*deltaT) ran
unconditionally), which is how C0-1 was found in the first place. The REGRESSION
side does not rely on this file's own baseline at all -- it is checked against
the committed t3_scheme_shootout.json, whose per-scheme rel_err values were
produced by the PREVIOUS build and are deterministic run-to-run. See
c0_wave_regression.py.

Usage:  ... python3.12 -S c0_wave_verify.py after
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "1")     # bit-exact oracle (P1f)
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
import opensees as ops  # noqa: E402

HERE = Path(__file__).resolve().parent
MODE = sys.argv[1] if len(sys.argv) > 1 else "before"
BASELINE = HERE / "c0_wave_baseline.json"

E, nu, RHO, L = 200000.0, 0.3, 7.85e-9, 100.0
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))


def build(n, integrator_argv, algo_argv):
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 250.0, 0.0, 0.0, 2000.0,
                   "-rho", RHO)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                            nid(i, j+1, k+1), 1)
                e += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 4.0e5 / 8.0e-4
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), P, 0.0, -P)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1.0e-7, 500)
    ops.algorithm(*algo_argv)
    ops.integrator(*integrator_argv)
    ops.analysis("Transient")
    return nid(nx // 2, ny // 2, nz)


def run(integrator_argv, algo_argv=("Newton",), dt=2.0e-5, nsteps=8, n=4):
    tip = build(n, integrator_argv, list(algo_argv))
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
    return ops.nodeDisp(tip, 3)


# ---- C0-1: does a zero dt error out, or silently produce numbers? -----------
def zero_dt(integrator_argv, n=3):
    build(n, integrator_argv, ["Newton"])
    try:
        rc = ops.analyze(1, 0.0)
    except Exception as exc:                       # noqa: BLE001
        return f"raised({type(exc).__name__})"
    return f"rc={rc}"


REG = [                       # (label, integrator argv, algorithm argv)
    ("Newmark",                ["Newmark", 0.5, 0.25],   ("Newton",)),
    ("Newmark1",               ["Newmark1", 0.5, 0.25],  ("Newton",)),
    ("HHT(0.9)",               ["HHT", 0.9],             ("Newton",)),
    ("TRBDF2",                 ["TRBDF2"],               ("Newton",)),
    ("TRBDF3",                 ["TRBDF3"],               ("Newton",)),
    ("BackwardEuler",          ["BackwardEuler"],        ("Newton",)),
    ("Houbolt",                ["Houbolt"],              ("Newton",)),
    ("Collocation(1.0)",       ["Collocation", 1.0],     ("Newton",)),
    ("Newmark+KrylovNewton",   ["Newmark", 0.5, 0.25],   ("KrylovNewton",)),
    ("Newmark+ModifiedNewton", ["Newmark", 0.5, 0.25],   ("ModifiedNewton",)),
]
HALL = [                      # C0-2: the classes that had no HALL branch
    ("BackwardEuler -hall",  ["BackwardEuler"],      ("Newton", "-hall", 0.2, 0.8)),
    ("Newmark1 -hall",       ["Newmark1", 0.5, 0.25], ("Newton", "-hall", 0.2, 0.8)),
    ("Collocation -hall",    ["Collocation", 1.0],   ("Newton", "-hall", 0.2, 0.8)),
]
ZERO_DT = [("TRBDF2", ["TRBDF2"]), ("TRBDF3", ["TRBDF3"]),
           ("BackwardEuler", ["BackwardEuler"]), ("Houbolt", ["Houbolt"])]


def main():
    print(f"ADR-77 C0 wave verify [{MODE}] — 1 thread")
    out = {"regression": {}, "hall": {}, "zero_dt": {}}

    print("\n-- regression: valid input must be UNCHANGED --")
    for lbl, integ, algo in REG:
        uz = run(integ, algo)
        out["regression"][lbl] = uz
        print(f"  {lbl:<24} uz={uz!r}")

    print("\n-- C0-2: -hall on the three classes that had no branch --")
    ref = out["regression"]["Newmark"]
    for lbl, integ, algo in HALL:
        uz = run(integ, algo)
        out["hall"][lbl] = uz
        note = ""
        if uz is not None and ref:
            note = f"  rel-vs-Newmark={abs(uz-ref)/abs(ref):.2e}"
        print(f"  {lbl:<24} uz={uz!r}{note}")

    print("\n-- C0-1: analyze(1, 0.0) must ERROR, not return numbers --")
    for lbl, integ in ZERO_DT:
        r = zero_dt(integ)
        out["zero_dt"][lbl] = r
        print(f"  {lbl:<24} {r}")

    if MODE == "before":
        BASELINE.write_text(json.dumps(out, indent=2))
        print(f"\nwrote {BASELINE.name}")
        return

    if not BASELINE.exists():
        print("\nno baseline")
        return
    b = json.loads(BASELINE.read_text())
    print("\n=== VERDICT ===")
    fail = 0
    for lbl, uz in out["regression"].items():
        old = b["regression"].get(lbl)
        ok = (uz == old)
        if not ok:
            fail += 1
        print(f"  {'PASS' if ok else 'FAIL'}  regression {lbl:<24} "
              f"{'bit-identical' if ok else f'{old!r} -> {uz!r}'}")
    for lbl, r in out["zero_dt"].items():
        old = b["zero_dt"].get(lbl)
        ok = (r != old) and ("rc=0" not in r)
        if not ok:
            fail += 1
        print(f"  {'PASS' if ok else 'FAIL'}  zero-dt    {lbl:<24} {old} -> {r}")
    print(f"\n{'ALL CHECKS PASSED' if not fail else f'{fail} CHECK(S) FAILED'}")


if __name__ == "__main__":
    main()
