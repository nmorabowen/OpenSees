"""ADR-41 D2 — viscous contact stabilization (p_visc = μ_c·v_rel).

A velocity-proportional NORMAL contact pressure that damps contact chatter / snap-through in the
pounding/rocking/uplift regime. It reuses the normal-penalty operator driven by velocity:

    gap rate ġ = B·v ;  viscous force f_visc = −μ_c·ġ·B ;  damping matrix C_visc = μ_c·B Bᵀ

active ONLY while in contact (gap<0); velocity ≡ 0 in statics ⇒ the term is identically inert
(static byte-identical). Force-only under the mass-only CentralDifferenceLadruno (explicit); force +
the consistent C_visc tangent (addCtoTang) under implicit Newmark.

It is a Kelvin–Voigt (spring+dashpot) contact, so a normal impact has the exact restitution
e = exp(−ζπ/√(1−ζ²)), ζ = μ_c/(2√(kₙ·m)) — pinned to machine tol in the numpy oracle
(contact_prototypes/proto_d2_viscous.py); here we confirm the C++ adapter reproduces it.

Gates (D2.1, NTS lane = RIGID_PLANE + SEGMENT):
  - explicit impact: μ_c=0 ⇒ e≈1 (elastic); μ_c>0 ⇒ e<1, monotone-decreasing — chatter energy bleed;
  - restitution matches the Kelvin–Voigt analytic e=exp(−ζπ/√(1−ζ²));
  - static solve with μ_c>0 == without (velocity ≡ 0 ⇒ viscous identically zero — byte-identical);
  - implicit Newmark with -visc CONVERGES (the consistent C_visc tangent) and dissipates energy;
  - viscous also acts on a SEGMENT (NTS faceted) contact (the second D2.1 mode);
  - `-visc` on a mortar contact ⇒ clean refusal (D2.2 scope).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
M = 1.0


def _anchor_truss(node, X0):
    """A weak truss to ground so a separated node is well-posed (no contact-only singularity)."""
    ops.node(900 + node, X0[0], X0[1], X0[2])
    ops.fix(900 + node, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, 1.0e-6)        # negligible vs KN
    ops.element("truss", 900 + node, 900 + node, node, 1.0, 900)


def _zeta(muc):
    return muc / (2.0 * math.sqrt(KN * M))


def _e_kelvin_voigt(muc):
    z = _zeta(muc)
    return math.exp(-z * math.pi / math.sqrt(1.0 - z * z)) if z < 1.0 else 0.0


def _explicit_impact(muc, v_mag=2.0, gap0=0.01, fine=True):
    """Mass node gap0 above a rigid plane (normal +z), flung down at v_mag. CentralDifferenceLadruno
    (the viscous term is FORCE-only here — mass-only LHS). Returns (blew_up, restitution, max_pen)."""
    n = [0.0, 0.0, 1.0]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    X0 = [gap0 * c for c in n]
    ops.node(1, X0[0], X0[1], X0[2])
    ops.mass(1, M, M, M)
    _anchor_truss(1, X0)
    ops.fix(1, 1, 1, 0)                                  # uniaxial along z (x,y fixed)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, n[0], n[1], n[2], 0.0, 0.0, 0.0, KN, "-visc", muc)
    ops.setNodeVel(1, 3, -v_mag, "-commit")
    T = 2.0 * math.pi * math.sqrt(M / KN)                # penalty period
    dt = (T / 60.0) if fine else (math.sqrt(M / KN))     # fine ⇒ ~30 steps in the contact half-period
    nsteps = int(40.0 * T / dt) + 200
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    max_pen, v_out, blew = 0.0, 0.0, False
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        z = X0[2] + ops.nodeDisp(1, 3)
        vz = ops.nodeVel(1, 3)
        if z < 0.0:
            max_pen = max(max_pen, -z)
        if z > 0.5 * gap0 and vz > 0.0:                  # separated and moving away
            v_out = vz
    return blew, abs(v_out / v_mag), max_pen


def test_d2_explicit_impact_dissipates_monotone():
    """μ_c=0 ⇒ near-elastic (e≈1); μ_c>0 ⇒ e<1; more μ_c ⇒ less rebound (monotone). The headline
    chatter/restitution control, through the real CentralDifferenceLadruno + the contact adapter."""
    b0, e0, _ = _explicit_impact(0.0)
    b1, e1, _ = _explicit_impact(200.0)                  # ζ≈0.10
    b2, e2, _ = _explicit_impact(500.0)                  # ζ≈0.25
    assert not (b0 or b1 or b2), "an impact blew up"
    assert e0 > 0.95, f"μ_c=0 should be ~elastic (e={e0:.3f})"
    assert e1 < e0 and e2 < e1, f"restitution not monotone in μ_c: e=({e0:.3f},{e1:.3f},{e2:.3f})"
    assert e2 < 0.7, f"μ_c=500 should clearly dissipate (e={e2:.3f})"


def test_d2_restitution_matches_kelvin_voigt():
    """The explicit impact restitution matches the analytic Kelvin–Voigt e=exp(−ζπ/√(1−ζ²))
    (the spring+dashpot contact) — the quantitative gate the numpy oracle pins to machine tol."""
    for muc in (200.0, 500.0):
        _, e, _ = _explicit_impact(muc, fine=True)
        e_kv = _e_kelvin_voigt(muc)
        assert abs(e - e_kv) < 0.06, (
            f"μ_c={muc}: e={e:.4f} vs Kelvin–Voigt {e_kv:.4f} (ζ={_zeta(muc):.3f})")


def _static_penetration(muc):
    """A node pressed into the rigid plane by a static downward load; returns the converged z-disp.
    Velocity ≡ 0 in statics, so the viscous term must be identically inert (μ_c cannot matter)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN, "-visc", muc)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)                        # push down into the plane
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "static contact did not converge"
    return ops.nodeDisp(1, 3)


def test_d2_static_byte_identical():
    """THE REGRESSION CONTRACT: in statics the nodal velocity is identically zero, so the viscous
    term contributes nothing — a static solve with -visc μ_c>0 is BIT-IDENTICAL to -visc 0
    (penetration = −P/kₙ either way)."""
    z_off = _static_penetration(0.0)
    z_on = _static_penetration(5.0e3)
    assert z_off == z_on, f"viscous perturbed the static solution: {z_off!r} vs {z_on!r}"
    assert abs(z_off - (-1.0e3 / KN)) < 1.0e-9, f"static penetration {z_off:.3e} != −P/kₙ"


def test_d2_implicit_newmark_converges_and_damps():
    """Implicit Newmark dynamic impact: -visc must (a) CONVERGE every step — exercising the
    consistent C_visc=μ_c·B Bᵀ damping tangent assembled in addCtoTang(c2) — and (b) dissipate
    energy vs -visc 0 (lower rebound velocity). Validates the implicit path (CDL is force-only)."""
    def newmark_impact(muc):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        gap0 = 0.01
        ops.node(1, 0.0, 0.0, gap0)
        ops.mass(1, M, M, M)
        ops.fix(1, 1, 1, 0)
        ops.contactSurface(1, "-slave", 1)
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN, "-visc", muc)
        ops.setNodeVel(1, 3, -2.0, "-commit")
        T = 2.0 * math.pi * math.sqrt(M / KN)
        dt = T / 40.0
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-10, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        v_out = 0.0
        for _ in range(int(60.0 * T / dt)):
            if ops.analyze(1, dt) != 0:
                return None                              # non-convergence ⇒ fail
            z = gap0 + ops.nodeDisp(1, 3)
            vz = ops.nodeVel(1, 3)
            if z > 0.5 * gap0 and vz > 0.0:
                v_out = vz
        return abs(v_out / 2.0)

    e_off = newmark_impact(0.0)
    e_on = newmark_impact(500.0)
    assert e_off is not None, "Newmark -visc 0 impact did not converge"
    assert e_on is not None, "Newmark -visc>0 did NOT converge (the C_visc tangent is the gate)"
    assert e_on < e_off, f"implicit viscous did not dissipate: e_on={e_on:.3f} !< e_off={e_off:.3f}"


def test_d2_segment_viscous_dissipates():
    """The SECOND D2.1 mode: viscous on a faceted NTS (SEGMENT) contact. A mass dropped onto a fixed
    master quad with -visc dissipates energy (lower rebound) vs -visc 0 — the B=[n|−Nᵢ n] operator
    path, distinct from the rigid-plane n-only path."""
    def seg_impact(muc):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        # fixed master quad in the z=0 plane
        for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
            ops.node(t, float(x), float(y), 0.0)
            ops.fix(t, 1, 1, 1)
        gap0 = 0.01
        ops.node(5, 0.0, 0.0, gap0)                      # slave above the quad centre
        ops.mass(5, M, M, M)
        ops.fix(5, 1, 1, 0)
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave", 5)
        # NTS contact: -kn KN, -outward +z, -visc muc
        ops.contact(1, 1, 2, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-visc", muc)
        ops.setNodeVel(5, 3, -2.0, "-commit")
        T = 2.0 * math.pi * math.sqrt(M / KN)
        dt = T / 60.0
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
        ops.analysis("Transient")
        v_out, blew = 0.0, False
        for _ in range(int(40.0 * T / dt) + 200):
            if ops.analyze(1, dt) != 0:
                blew = True
                break
            z = gap0 + ops.nodeDisp(5, 3)
            vz = ops.nodeVel(5, 3)
            if z > 0.5 * gap0 and vz > 0.0:
                v_out = vz
        return blew, abs(v_out / 2.0)

    b0, e0 = seg_impact(0.0)
    b1, e1 = seg_impact(500.0)
    assert not (b0 or b1), "segment impact blew up"
    assert e0 > 0.9, f"segment μ_c=0 should be ~elastic (e={e0:.3f})"
    assert e1 < 0.7, f"segment viscous should dissipate (e={e1:.3f} vs e0={e0:.3f})"


def test_d2_visc_refused_on_mortar():
    """-visc is NTS-only in D2.1; combined with -mortar it must be REFUSED (not silently ignored) —
    mortar viscous is the D2.2 scope."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(4 + t, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-visc", 100.0)
