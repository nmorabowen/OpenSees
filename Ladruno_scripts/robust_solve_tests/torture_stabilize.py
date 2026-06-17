"""Torture fixture #3 - the Layer-1.5 stabilization-energy observability seam.

Exercises the rung-4 gate getters the robust-solve driver reads (ADR-31
"gate vs audit"; LadrunoArcLength 33004):

    ladrunoArcLength dissipationRatio | dissipatedEnergy | referenceEnergy
    ladrunoArcLength scaleCVisc <factor>          (the R-RAMPDOWN actuator)

The model is the softening Concrete02 truss (same physics as torture_softening),
driven under `-stabilize` viscous load control. We do NOT assert that stabilize
traces the softening branch - it cannot (pure monotone softening has no
equilibrium above the peak, so stabilized LOAD control stalls at the peak exactly
like plain load control; that is rung-3's job). What we assert is that the
observability seam reports physically-sane dissipation and that the actuator
ramps the viscous coefficient down (verified against the live build, not guessed).

Run standalone:
    python Ladruno_scripts/robust_solve_tests/torture_stabilize.py
"""
import os
import sys

DIST = os.environ.get("LADRUNO_DIST",
                      r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

FPC, EPSC0, FPCU, EPSU = -30.0, -0.002, -6.0, -0.012
PEAK = 1.0 * abs(FPC)              # peak compressive load = 30
FTARGET = 2.0e-4                   # Abaqus-default dissipated-energy fraction


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Concrete02", 1, FPC, EPSC0, FPCU, EPSU, 0.1, 3.0, 300.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _arm(adaptStab):
    dlam = PEAK / 200.0
    args = ["LadrunoArcLength", dlam, 1.0, "-stabilize", FTARGET]
    if adaptStab:
        args.append("-adaptStab")
    ops.integrator(*args)


def run_stabilized(adaptStab=False, n_steps=200):
    """Drive the softening truss under -stabilize; return the seam's readout."""
    build()
    _arm(adaptStab)
    committed = 0
    for _ in range(n_steps):
        if ops.analyze(1) != 0:
            break
        committed += 1
    return {
        "committed": committed,
        "ratio": ops.ladrunoArcLength("dissipationRatio"),
        "dissip": ops.ladrunoArcLength("dissipatedEnergy"),
        "ref": ops.ladrunoArcLength("referenceEnergy"),
    }


def peak_load_stabilized(f, n_steps=300):
    """Drive the softening truss under -stabilize at fraction f; return peak |load|.

    The metric for the c-reduction diffusion bound (R-DIFFUSION): on pure softening
    the converged load tops out at the physical strength peak (~PEAK), which is
    essentially independent of the artificial dissipation fraction f — so halving f
    barely moves it, and `diffusion_drift` is small (the answer is trustworthy).
    """
    build()
    dlam = PEAK / 200.0
    ops.integrator("LadrunoArcLength", dlam, 1.0, "-stabilize", f)
    peak = 0.0
    for _ in range(n_steps):
        if ops.analyze(1) != 0:
            break
        peak = max(peak, abs(ops.getLoadFactor(1)))
    return peak


def measure_rampdown(factor=0.1, warmup=40, window=5):
    """Verify scaleCVisc ramps the per-window artificial dissipation down."""
    build()
    _arm(adaptStab=False)
    for _ in range(warmup):
        if ops.analyze(1) != 0:
            break
    d0 = ops.ladrunoArcLength("dissipatedEnergy")
    for _ in range(window):
        ops.analyze(1)
    inc_full = ops.ladrunoArcLength("dissipatedEnergy") - d0
    ops.ladrunoArcLength("scaleCVisc", factor)          # R-RAMPDOWN
    d1 = ops.ladrunoArcLength("dissipatedEnergy")
    for _ in range(window):
        ops.analyze(1)
    inc_ramped = ops.ladrunoArcLength("dissipatedEnergy") - d1
    return inc_full, inc_ramped


if __name__ == "__main__":
    a = run_stabilized(adaptStab=True)
    b = run_stabilized(adaptStab=False)
    inc_full, inc_ramped = measure_rampdown()

    print(f"{'config':<28}{'committed':>10}{'ratio':>14}{'ref(Estrain0)':>16}")
    print("-" * 68)
    print(f"{'-stabilize -adaptStab':<28}{a['committed']:>10}{a['ratio']:>14.4e}"
          f"{a['ref']:>16.4e}")
    print(f"{'-stabilize (cumulative)':<28}{b['committed']:>10}{b['ratio']:>14.4e}"
          f"{b['ref']:>16.4e}")
    print(f"\nscaleCVisc(0.1) ramp-down: {inc_full:.4e} -> {inc_ramped:.4e} "
          f"({'OK' if inc_ramped < inc_full else 'FAIL'})")

    ident = abs(a["ratio"] - a["dissip"] / a["ref"]) <= 1e-9 * a["ratio"]
    bounded = a["ratio"] < 100 * FTARGET          # adaptStab holds it near fTarget
    drifts = b["ratio"] > a["ratio"]              # cumulative (no adaptStab) is larger
    ramps = inc_ramped < inc_full
    ok = ident and bounded and drifts and ramps and a["ref"] > 0.0
    print("\nSELF-TEST:", "PASS - Layer-1.5 dissipation seam reports sane values "
          "and scaleCVisc ramps c down" if ok else "FAIL - inspect above")
    raise SystemExit(0 if ok else 1)
