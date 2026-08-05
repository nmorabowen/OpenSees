"""Calibrate an elastic-perfectly-plastic Drucker-Prager surrogate onto the
cone `cone_probe.py` MEASURED inside PressureDependMultiYield, and verify it.

Why a surrogate at all. `cone_probe.py` established that PDMY's failure surface
is a Lode-independent cone of slope alpha = sqrt(J2)/(3p) = 0.2436, and that
Vesic's N_gamma — a Mohr-Coulomb formula — evaluated at the cone's PLANE-STRAIN
equivalent (53.7 deg) is 52x the anchor the project quotes. That factor rests on
Chen & Han matching, which assumes ASSOCIATED flow (this PDMY set is
non-dilatant: dil1 = dil2 = 0) and treats a SQUARE footing as plane strain.
Neither assumption is safe, so the number is a hypothesis, not a measurement.

The measurement is a numerical collapse load on the actual 3D footing with the
actual cone. PDMY itself cannot deliver it: its 20 nested yield surfaces
approach the envelope asymptotically, so it hardens for as long as you push
(11 h to s/B = 0.15, still climbing at 31 % of the envelope). An
elastic-PERFECTLY-plastic material sharing the same failure surface reaches its
limit load in a few percent of settlement, and by the limit theorems the
collapse load depends on the FAILURE SURFACE and the flow rule — not on the
hardening path that leads there. So the surrogate is the right instrument, and
this script is where it earns the right to be called "the same cone".

Conventions. OpenSees `DruckerPrager` yields at

    ||s|| + rho*I1 - sqrt(2/3)*sigma_y = 0,     I1 = tr(sigma) (compression < 0)

With p = -I1/3 and ||s|| = sqrt(2*J2) that is sqrt(J2) = (3*rho/sqrt(2))*p +
sigma_y/sqrt(3), i.e.

    alpha = rho / sqrt(2)          <=>   rho = sqrt(2) * alpha
    cohesion intercept in sqrt(J2) k = sigma_y / sqrt(3)

sigma_y is NOT physical here — PDMY is cohesionless. It is a small numerical
regularizer that keeps the cone apex off the origin where the surface DOFs of a
free surface live. It is swept, and the collapse run reports the sweep.

Instrument. Two independent probes of the SAME surface:

  * `strain` — lateral stress held at p0 by constant face loads, axial
    displacement DRIVEN. Perfect plasticity gives a genuine PLATEAU in sigma_zz,
    which locates the surface exactly rather than from just inside it.
  * `stress` — the identical stress-controlled probe `cone_probe.py` ran on
    PDMY: load until non-convergence, read the last converged state. Biased LOW
    by construction; running it here MEASURES that bias (PDMY's own alpha is
    biased the same way, so the two materials are compared like with like).

Paths txc (axial compression, Lode -30 deg) and txe (axial unloading, +30 deg)
are the two Lode extremes — the pair that makes a cone test a cone test.

Run:  py -3.12 dp_cone_check.py
"""
import os
import sys

import numpy as np

_WT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
_BIN = os.path.join(_WT, "dist", "bin")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
import opensees as ops  # noqa: E402

if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
        os.path.normcase(_BIN):
    raise SystemExit(f"WRONG ENGINE: {ops.__file__}, expected {_BIN}")

# --- the cone, as measured on PDMY by cone_probe.py -------------------------
ALPHA = 0.24362                  # sqrt(J2)/(3p), txc + txe fit, 3.9 % spread
RHO_DP = np.sqrt(2.0) * ALPHA    # OpenSees DruckerPrager `rho`
# PDMY's elastic reference moduli. DruckerPrager cannot combine pressure-
# dependent moduli with plasticity (mElastFlag == 1 is elastic-only), so these
# are constants. Collapse loads do not depend on elastic constants.
K_EL, G_EL = 1.5e5, 5.5e4        # kPa
SY_LIST = [0.2, 2.0]             # sqrt(3)*k: numerical cohesion sweep, kPa
P0 = [50.0, 100.0, 200.0]
XF, YF, ZF = (2, 3, 6, 7), (3, 4, 7, 8), (5, 6, 7, 8)


def invariants(s):
    xx, yy, zz, xy, yz, zx = s
    T = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
    w = np.sort(np.linalg.eigvalsh(T))[::-1]
    I1 = xx + yy + zz
    p = -I1 / 3.0
    dev = T - (I1 / 3.0) * np.eye(3)
    J2 = 0.5 * np.tensordot(dev, dev)
    den = -(w[0] + w[2])
    sphi = (w[0] - w[2]) / den if den > 1e-9 else np.nan
    return p, np.sqrt(J2), sphi


def build(sy, rho_bar):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
                                  start=1):
        ops.node(i, float(x), float(y), float(z))
    for i in (1, 4, 5, 8):
        ops.fix(i, 1, 0, 0)
    for i in (1, 2, 5, 6):
        ops.fix(i, 0, 1, 0)
    for i in (1, 2, 3, 4):
        ops.fix(i, 0, 0, 1)
    # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta <rho_mass> <atm>
    ops.nDMaterial("DruckerPrager", 1, K_EL, G_EL, sy, RHO_DP, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1, "-geom", "linear")


def consolidate(p0):
    """Isotropic consolidation to p0, elastic by construction (K0 = 1 is deep
    inside any frictional cone)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in XF:
        ops.load(i, -p0 / 4.0, 0.0, 0.0)
    for i in YF:
        ops.load(i, 0.0, -p0 / 4.0, 0.0)
    for i in ZF:
        ops.load(i, 0.0, 0.0, -p0 / 4.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-6 * p0, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    assert ops.analyze(10) == 0, f"consolidation failed at p0={p0}"


def stress_of(): return np.array(ops.eleResponse(1, "stress")).reshape(-1, 6).mean(axis=0)


def run_strain(path, p0, sy, rho_bar):
    """Lateral stress held; axial displacement driven -> a real plateau."""
    build(sy, rho_bar)
    consolidate(p0)
    ops.loadConst("-time", 0.0)
    sgn = -1.0 if path == "txc" else +1.0          # compress / unload axially
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for i in ZF:
        ops.sp(i, 3, sgn * 1.0)                    # 1.0 m per unit load factor
    ops.test("NormUnbalance", 1e-6 * p0, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 2.0e-5)          # 20 microstrain per step
    # yield arrives at ~8e-4 axial strain (dq/3G with q_yield ~ 1.27 p0), so
    # 300 steps = 0.6 % strain is several times past it — the tail IS plateau.
    hist = []
    for _ in range(300):
        if ops.analyze(1) != 0:
            break
        hist.append(stress_of())
    # the plateau: last 10 % of the history, and its spread is the evidence
    tail = np.array(hist[-max(3, len(hist) // 10):])
    inv = np.array([invariants(s)[:2] for s in tail])
    a = inv[:, 1] / (3.0 * inv[:, 0])
    p, sq, sphi = invariants(tail[-1])
    return dict(alpha=float(a.mean()), drift=float(a.max() - a.min()),
                p=p, sqrtJ2=sq, sphi=sphi, n=len(hist))


def run_stress(path, p0, sy, rho_bar, nstep=400):
    """cone_probe.py's instrument, unchanged, applied to the surrogate."""
    build(sy, rho_bar)
    consolidate(p0)
    # HOLD the consolidation pattern. Without this its Linear series keeps
    # scaling with the same load factor as the deviatoric pattern below, and
    # the very first "increment" applies the deviator at full amplitude —
    # the probe then reports failure at step 0 with zero deviator. (Measured.)
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    amp = 4.0 * p0
    if path == "txc":
        for i in ZF:
            ops.load(i, 0.0, 0.0, -amp / 4.0)
    else:
        for i in ZF:
            ops.load(i, 0.0, 0.0, +0.98 * p0 / 4.0)
    ops.test("NormUnbalance", 1e-6 * p0, 60, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nstep)
    last, reached = stress_of(), False
    for _ in range(nstep):
        if ops.analyze(1) != 0:
            reached = True
            break
        last = stress_of()
    p, sq, sphi = invariants(last)
    return dict(alpha=sq / (3.0 * p), p=p, sqrtJ2=sq, sphi=sphi,
                reached=reached)


def main():
    print(f"target cone (measured on PDMY): alpha = {ALPHA:.5f}")
    print(f"  => DruckerPrager rho = sqrt(2)*alpha = {RHO_DP:.6f}\n")
    for sy in SY_LIST:
        for rho_bar, lbl in ((RHO_DP, "associated"), (0.0, "non-dilatant")):
            print(f"--- sigma_y = {sy} kPa   rho_bar = {rho_bar:.4f} ({lbl})")
            print(f"{'probe':>7} {'path':>5} {'p0':>6} {'p_f':>9} "
                  f"{'sqrtJ2':>9} {'alpha':>9} {'err%':>7} {'phi_mob':>8} "
                  f"{'drift':>9}")
            acc = {}
            for path in ("txc", "txe"):
                for p0 in P0:
                    r = run_strain(path, p0, sy, rho_bar)
                    acc.setdefault(("strain", path), []).append(r["alpha"])
                    print(f"{'strain':>7} {path:>5} {p0:6.0f} {r['p']:9.2f} "
                          f"{r['sqrtJ2']:9.2f} {r['alpha']:9.5f} "
                          f"{100 * (r['alpha'] / ALPHA - 1):7.2f} "
                          f"{np.degrees(np.arcsin(min(r['sphi'], 1.0))):8.2f} "
                          f"{r['drift']:9.2e}")
                for p0 in P0:
                    r = run_stress(path, p0, sy, rho_bar)
                    acc.setdefault(("stress", path), []).append(r["alpha"])
                    print(f"{'stress':>7} {path:>5} {p0:6.0f} {r['p']:9.2f} "
                          f"{r['sqrtJ2']:9.2f} {r['alpha']:9.5f} "
                          f"{100 * (r['alpha'] / ALPHA - 1):7.2f} "
                          f"{np.degrees(np.arcsin(min(r['sphi'], 1.0))):8.2f} "
                          + ("" if r["reached"] else "  NOT AT FAILURE"))
            for probe in ("strain", "stress"):
                v = np.concatenate([acc[(probe, p)] for p in ("txc", "txe")])
                print(f"    {probe}: alpha = {v.mean():.5f} "
                      f"({100 * (v.mean() / ALPHA - 1):+.2f} % vs target), "
                      f"Lode spread = "
                      f"{100 * (v.max() - v.min()) / v.mean():.2f} %")
            print()


if __name__ == "__main__":
    main()
