"""ADR-41 C3.2 — consistent (symmetric) mortar friction TANGENT: implicit frictional Newton.

C3.1 shipped the friction FORCE but no friction tangent, so a free tangential DOF under an
implicit solve is SINGULAR (the ADR-39 P3.5 lesson). C3.2 assembles the symmetric friction
tangent K[A][B] += Σ_I (b_IA b_IB / a_I) K_ss_I (b = the B̃=[D,−M] operator row, K_ss the 3×3
LadrunoFrictionKernel::frictionTangentBlock at the per-node slip) into addMortarTang, so static /
Newmark frictional Newton CONVERGES. Default SYMMETRIC (drop the non-symmetric Coulomb Csl —
solver-safe on any SOE, the P3.5 Q2 rule; the full mortar Csl is deferred with the geometric terms).

Gates:
  - static frictional mortar Newton CONVERGES (THE gate — singular without the tangent) for both
    STICK (Q < cap ⇒ tiny tangential disp ≈ Q/(epsT·A)) and SLIP (Q > cap ⇒ slides, reaction = cap);
  - the SAME model converges on a SYMMETRIC solver (ProfileSPD) — the default tangent is symmetric;
  - μ=0 byte-identical to frictionless C2 (the tangent short-circuits too).

Assembled-tangent FD validated oracle-first in proto_c3_mortar_friction.py T6 (rel ~1e-10).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e4
DELTA = 1.0e-4    # seeded normal penetration (contact active from step 1)
COH = 50.0        # Tresca cohesion ⇒ cap = COH·A, pressure-independent (Uzawa-immune)


def _master(z=0.0):
    for t, x, y in [(1, 0, 0), (2, 1, 0), (3, 1, 1), (4, 0, 1)]:
        ops.node(t, float(x), float(y), z)
        ops.fix(t, 1, 1, 1)


def _slave(z):
    for t, x, y in [(5, 0, 0), (6, 1, 0), (7, 1, 1), (8, 0, 1)]:
        ops.node(t, float(x), float(y), z)
        ops.fix(t, 0, 0, 0)            # all free: z held by the normal penalty, x,y by friction


def _solve(system="FullGeneral"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1.0e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    return ops.analyze(1)


def _frictional_static(Qx, Pz=200.0, mu=0.0, cohesion=COH, system="FullGeneral"):
    """A slave facet pressed onto a fixed master (−Pz) + a tangential load Qx, Tresca friction.
    Returns (rc, x-disp of node 5). Newton MUST converge (the friction tangent is what makes the
    free tangential DOFs non-singular)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master(0.0)
    _slave(-DELTA)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    opts = ["-mortar", "-epsN", EPS, "-epsT", EPS, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts = ["-mortar", "-epsN", EPS, "-mu", mu, "-epsT", EPS, "-outward", 0.0, 0.0, 1.0]
    if cohesion > 0.0:
        opts += ["-cohesion", cohesion]
    ops.contact(1, 1, 2, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, Qx / 4.0, 0.0, -Pz / 4.0)
    rc = _solve(system)
    return rc, ops.nodeDisp(5, 1)


def test_c3_2_static_stick_converges():
    """STICK: a tangential load below the friction cap (cap = COH·A = 50). Static Newton CONVERGES
    (the free x,y DOFs are non-singular ONLY because the friction stick tangent kt·P_t is assembled
    — without C3.2 this solve is singular), and the slave barely moves (sticks)."""
    Qx = 20.0                                       # < cap = 50 ⇒ stick
    rc, ux = _frictional_static(Qx)
    assert rc == 0, "static frictional Newton did not converge (missing/wrong friction tangent?)"
    # stick: the elastic tangential slip is ≈ (Q/A)/epsT — small. NOT free sliding.
    assert abs(ux) < 1.0e-2, f"slave slid {ux} — expected stick (tiny disp)"
    assert abs(ux) == pytest.approx((Qx / 1.0) / EPS, rel=0.10), f"stick disp {ux} != Q/(epsT·A)"


def test_c3_2_static_slip_converges():
    """SLIP: a tangential load above the cap. A free block under Q>cap has NO static equilibrium
    (it slides indefinitely — the slip-direction friction tangent is zero by construction), so the
    slave is restrained by a small tangential SPRING (kx) that supplies the slip-direction stiffness.
    Static Newton then CONVERGES on the slip tangent, the friction saturates at the cap, and the
    spring carries the excess: kx·x = Q − cap·A. This exercises the slip branch in a real solve."""
    Qx, Pz, kx = 120.0, 200.0, 1.0e3               # Q=120 > cap=COH·A=50
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master(0.0)
    _slave(-DELTA)
    ops.uniaxialMaterial("Elastic", 1, kx)
    for i, t in enumerate((5, 6, 7, 8)):           # a tangential (x) spring to ground per slave node
        gt = 200 + t
        ops.node(gt, *[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)][i])
        ops.fix(gt, 1, 1, 1)
        ops.element("zeroLength", 300 + t, gt, t, "-mat", 1, "-dir", 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-epsT", EPS, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, Qx / 4.0, 0.0, -Pz / 4.0)
    assert _solve() == 0, "slip-regime frictional Newton did not converge"
    ux = ops.nodeDisp(5, 1)
    # at slip the friction saturates at cap=COH·A=50; the spring (4·kx total) carries Q−cap.
    cap = COH * 1.0
    x_expected = (Qx - cap) / (4.0 * kx)           # 4 springs of kx, total stiffness 4·kx
    assert abs(ux) > 1.0e-2, f"slave did not slide ({ux}) — expected slip past the cap"
    assert ux == pytest.approx(x_expected, rel=0.05), f"slip disp {ux} != (Q−cap)/(4kx)={x_expected}"


def test_c3_2_symmetric_solver_safe():
    """The DEFAULT friction tangent is SYMMETRIC, so it is correct on a symmetric SOE (ProfileSPD
    reads only the upper triangle). The stick solve converges identically on ProfileSPD and
    FullGeneral — a non-symmetric default would silently corrupt the symmetric solver (P3.5 Q2)."""
    rc_sym, ux_sym = _frictional_static(20.0, system="ProfileSPD")
    rc_full, ux_full = _frictional_static(20.0, system="FullGeneral")
    assert rc_sym == 0, "frictional Newton failed on the symmetric ProfileSPD solver"
    assert rc_full == 0
    assert abs(ux_sym - ux_full) < 1.0e-9, f"symmetric vs full disagree: {ux_sym} vs {ux_full}"


def test_c3_2_coulomb_static_converges():
    """Coulomb (μ>0) static Newton converges too — the cap = μN tracks the (penalty + one-commit
    Uzawa) normal pressure, but a single LoadControl step holds λ_N=0 during the Newton, so the
    cap is μ·epsN·δ-stable within the solve. Confirms the Coulomb friction tangent (not just Tresca)."""
    rc, ux = _frictional_static(5.0, Pz=200.0, mu=0.5, cohesion=0.0)
    assert rc == 0, "Coulomb frictional Newton did not converge"
    assert abs(ux) < 1.0, "unphysical tangential motion under Coulomb stick"
