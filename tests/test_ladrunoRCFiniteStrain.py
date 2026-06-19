"""LadrunoRCFiniteStrain — finite-strain (Hencky) view of the RC plastic-damage
material (classTag 33018, ADR 19 Phase 4b) — Zone-A, driven end-to-end through
`element LadrunoBrick ... -geom finite` (which calls setTrialF(F)).

The material is a FiniteStrainNDMaterial subclass, so unlike the LogStrain wrapper
it is passed DIRECTLY to the element (no `nDMaterial LogStrain` indirection):

    nDMaterial LadrunoRCFiniteStrain 1 E nu -Ce {..} -Cs {..} -Te {..} -Ts {..} ...
    element LadrunoBrick 1 ...nodes... 1 -formulation std -geom finite

Seam:  B = F Fᵀ  →  εᵉ = ½ln B  →  the SHARED LadrunoRCKernel returnMap3D  →
Cauchy σ = τ/J + spatial tangent c = (1/2J)[D:L:B]. Because the RC spine carries
no tensorial plastic strain, the elastic left Cauchy–Green is the TOTAL B (no bᵉ
to track) — the whole correctness of the stress seam therefore reduces to: the
finite Cauchy stress must equal the SMALL-STRAIN material evaluated at the same
Hencky strain, divided by J. That is the patch oracle below — an exact, same-binary
cross-check needing no separate numpy reference.

Gates:
  * tiny strain  → finite reduces to the small-strain LadrunoBrick -geom linear;
  * ELASTIC finite deformation → consistent spatial tangent K == FD of the force
    (validates the LogStrainKernel spatial-tangent seam; the damage-range tangent
    is a deliberate secant, so the FD gate is restricted to the elastic regime);
  * finite stretch INTO DAMAGE → GP Cauchy σ == (small-strain σ at ½ln B)/J;
  * multi-step committed path → same final Cauchy as a single jump (the committed Hencky
    strain in RCHist carrying correctly across commitState → next setTrialF);
  * rigid rotation of the unloaded element → stress-free;
  * ISOTROPIC-spine objectivity (no interlock): σ(QF) == Q σ(F) Qᵀ at two rotation
    magnitudes (the OBJECTIVE half of §14.11; the directional fixed-crack/interlock state
    is NOT objective under a committed-then-rotated path — documented in the class header
    + ADR, with the corotational ELEMENT as the objective route, tested for the small-strain
    material in test_ladrunoRCConcrete_objectivity.py);
  * IMPL-EX completes and tracks the implicit run;
  * det F ≤ 0 is rejected (the element cuts the step);
  * a committed finite state survives the FE-database round-trip.
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

# distorted 8-node brick (good tangent test) — same geometry family as the J2 finite battery
_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010, 0.0040]
TS = [0.0, 3.0,    0.5,    0.0]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0, 0.999]


def _mat_fin(tag, interlock=False, cyclic=False, implex=False):
    a = ["LadrunoRCFiniteStrain", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
         "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if interlock:
        a += ["-interlock", "-agg", 16.0, "-crackSpacing", 50.0]
        if cyclic:
            a += ["-cyclic"]
    if implex:
        a += ["-implex"]
    ops.nDMaterial(*a)


def _mat_small(tag):
    ops.nDMaterial("LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC)


def _build(geom, mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "std", "-geom", geom)


def _impose_and_solve(u, geom, mat_fn, nsteps=1):
    _build(geom, mat_fn)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    rc = 0
    for _ in range(nsteps):
        rc = ops.analyze(1)
        if rc != 0:
            break
    return rc


def _resisting_force(u, mat_fn, geom="finite"):
    assert _impose_and_solve(u, geom, mat_fn) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _affine_disp(Fbar):
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        d = (np.asarray(Fbar) - I) @ np.array([x, y, z])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


def _hencky(F):
    """Hencky strain tensor εᵉ = ½ ln(F Fᵀ)."""
    B = np.asarray(F) @ np.asarray(F).T
    w, V = np.linalg.eigh(B)
    return V @ np.diag(0.5 * np.log(w)) @ V.T


def _gp_cauchy(geom, mat_fn, u, nsteps=1):
    assert _impose_and_solve(u, geom, mat_fn, nsteps) == 0, "solve failed"
    return np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)


def _to_tensor(v):
    return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]])


# --------------------------------------------------------------------------
def test_finite_rc_reduces_to_small_strain():
    """At tiny strain the finite (Hencky) view reduces to the small-strain RC on
    LadrunoBrick -geom linear: the resisting forces coincide."""
    eps = 1.0e-6
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    u = _affine_disp(Fbar)
    f_fin = _resisting_force(u, _mat_fin, "finite")
    f_lin = _resisting_force(u, lambda t: _mat_small(t), "linear")
    scale = max(np.abs(f_lin).max(), 1.0e-30)
    assert np.abs(f_fin - f_lin).max() <= 1.0e-3 * scale, "finite did not reduce to small-strain"


def test_finite_rc_elastic_tangent_matches_fd():
    """In the ELASTIC finite regime the spatial tangent (LogStrainKernel channel-A)
    must equal the FD of the resisting force — validates the finite-strain seam. The
    damage-range tangent is a deliberate secant, so the FD gate stays elastic."""
    # a modest deviatoric stretch that stays below the tensile crack + compressive
    # damage thresholds (max principal Hencky strain well under ft/E ~ 1e-4 in tension)
    Fbar = [[1.00003, 0.00002, 0.0],
            [0.00002, 0.99997, 0.0],
            [0.0, 0.0, 1.00001]]
    u = _affine_disp(Fbar)
    assert _impose_and_solve(u, "finite", _mat_fin) == 0
    dmg = list(ops.eleResponse(1, "material", 1, "damage"))
    assert max(abs(v) for v in dmg) < 1.0e-9, f"tangent-FD gate not elastic (damage {dmg})"
    K = np.array(ops.eleResponse(1, "stiff"), dtype=float).reshape(24, 24)
    assert np.isfinite(K).all()
    h = 1.0e-8
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, _mat_fin) - _resisting_force(um, _mat_fin)) / (2.0 * h)
    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    assert err <= 1.0e-4 * scale, f"finite elastic tangent != FD: err {err:.3e} vs {scale:.3e}"


def test_finite_rc_patch_matches_smallstrain_oracle():
    """The stress seam: a homogeneous finite stretch INTO tensile damage. Every GP
    Cauchy σ must equal (the small-strain RC stress at the same Hencky strain ½ln B)
    divided by J — an exact same-binary cross-check (same kernel, two strain measures)."""
    lam = 1.01                                   # ~1% axial Hencky => well past cracking
    lat = 1.0 / math.sqrt(lam)
    Fm = np.diag([lam, lat, lat])
    J = float(np.linalg.det(Fm))
    sig = _gp_cauchy("finite", _mat_fin, _affine_disp(Fm.tolist()))
    s0 = sig[0]
    assert np.abs(sig - s0).max() <= 1.0e-6 * max(np.abs(s0).max(), 1.0), "Cauchy not constant across GPs"
    # guard the cross-check is meaningful: the state is genuinely damaged (read on the finite run)
    dmg = list(ops.eleResponse(1, "material", 1, "damage"))
    assert max(abs(v) for v in dmg) > 0.0, "patch oracle is vacuous — no damage formed"

    # oracle: small-strain RC at the Hencky strain (impose eps_h affinely on -geom linear)
    eps_h = _hencky(Fm)
    tau = _gp_cauchy("linear", lambda t: _mat_small(t), _affine_disp((np.eye(3) + eps_h).tolist()))[0]
    sig_oracle = tau / J
    tol = 1.0e-5 * max(np.abs(sig_oracle).max(), 1.0)
    assert np.allclose(s0, sig_oracle, rtol=1.0e-5, atol=tol), (
        f"finite Cauchy {s0} != small-strain(½lnB)/J {sig_oracle}")


def test_finite_rc_multistep_committed_matches_single_step():
    """Multi-step committed finite path: ramping to Fm in N committed LoadControl steps must
    reach the same final Cauchy as a single jump to Fm (monotone proportional stretch ⇒ the
    damage threshold = max equivalent strain is path-independent). Exercises the compiled
    commitState → next setTrialF round-trip (the committed Hencky strain in RCHist carrying
    across steps) — a step-2 accumulation slip would pass every single-step gate but fail here."""
    lam = 1.012
    Fm = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)])
    u = _affine_disp(Fm.tolist())
    s_single = _gp_cauchy("finite", _mat_fin, u, nsteps=1)[0]
    s_multi = _gp_cauchy("finite", _mat_fin, u, nsteps=10)[0]
    # guard the path actually damaged (else this is a trivial elastic check)
    dmg = list(ops.eleResponse(1, "material", 1, "damage"))
    assert max(abs(v) for v in dmg) > 0.0, "multistep path is vacuous — no damage"
    scale = max(np.abs(s_single).max(), 1.0)
    assert np.abs(np.array(s_multi) - np.array(s_single)).max() <= 1.0e-4 * scale, (
        f"multi-step committed finite path diverged from single-step: {s_multi} vs {s_single}")


def test_finite_rc_rigid_rotation_stress_free():
    """Pure rigid rotation of the unloaded element induces no stress / force."""
    cz, sz = math.cos(0.4), math.sin(0.4)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(0.3), math.sin(0.3)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    R = Ry @ Rz
    assert _impose_and_solve(_affine_disp(R), "finite", _mat_fin) == 0
    s = ops.eleResponse(1, "stresses")
    assert max(abs(v) for v in s) <= 1.0e-6 * max(TS), "rigid rotation induced stress"


def test_finite_rc_isotropic_objectivity_holds():
    """ISOTROPIC-damage-spine objectivity (no interlock ⇒ no stored crack frame):
    σ(Q F) == Q σ(F) Qᵀ for a damaging stretch. The spectral split of B and the
    scalar damage thresholds are frame-indifferent, so the wrapper IS objective here."""
    lam = 1.008
    Fm = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)])
    s_ref = _gp_cauchy("finite", _mat_fin, _affine_disp(Fm.tolist()))[0]
    th = 0.5
    c, s = math.cos(th), math.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    s_rot = _gp_cauchy("finite", _mat_fin, _affine_disp((Q @ Fm).tolist()))[0]

    def to_tensor(v):
        return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]])
    expected = Q @ to_tensor(s_ref) @ Q.T
    got = to_tensor(s_rot)
    tol = 1.0e-5 * max(np.abs(got).max(), 1.0)
    assert np.allclose(got, expected, rtol=1.0e-5, atol=tol), (
        "isotropic spine not objective under finite rotation")


def test_finite_rc_isotropic_objectivity_under_two_rotations():
    """Reinforce the isotropic-spine objectivity across a second rotation magnitude: σ(QF)
    == Q σ(F) Qᵀ for a damaging stretch, Q a 90° rotation about z (swaps x↔y axes). This is
    the OBJECTIVE half of the §14.11 split — the scalar-damage spine is frame-indifferent.
    (The directional fixed-crack/interlock state is NOT objective under a committed-then-
    rotated path; that boundary is documented in the class header + ADR and the objective
    route is the corotational ELEMENT, validated for the small-strain material in
    test_ladrunoRCConcrete_objectivity.py. It is not re-tested here at the element level
    because a committed-crack-then-large-rotation static solve fights the softening interlock
    tangent — the same reason the cyclic wall validation is quasi-static explicit.)"""
    lam = 1.007
    Fm = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)])
    s_ref = _gp_cauchy("finite", _mat_fin, _affine_disp(Fm.tolist()))[0]
    th = math.pi / 2.0
    c, s = math.cos(th), math.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    s_rot = _gp_cauchy("finite", _mat_fin, _affine_disp((Q @ Fm).tolist()))[0]
    expected = Q @ _to_tensor(s_ref) @ Q.T
    tol = 1.0e-5 * max(np.abs(expected).max(), 1.0)
    assert np.allclose(_to_tensor(s_rot), expected, rtol=1.0e-5, atol=tol), (
        "isotropic spine not objective under a 90° finite rotation")


def test_finite_rc_implex_runs_and_tracks():
    """-implex completes a finite stretch into damage and reports a finite implexError;
    the committed stress stays within a loose band of the implicit run (O(dt) lag)."""
    lam = 1.006
    Fm = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)])
    u = _affine_disp(Fm.tolist())
    s_imp = _gp_cauchy("finite", lambda t: _mat_fin(t), u, nsteps=8)[0]
    s_iex = _gp_cauchy("finite", lambda t: _mat_fin(t, implex=True), u, nsteps=8)[0]
    err = list(ops.eleResponse(1, "material", 1, "implexError"))
    assert err and math.isfinite(err[0])
    scale = max(np.abs(s_imp).max(), 1.0)
    assert np.abs(np.array(s_iex) - np.array(s_imp)).max() <= 0.30 * scale, (
        f"implex strayed too far from implicit: {s_iex} vs {s_imp}")


def test_finite_rc_detF_guard():
    """A strongly inverting deformation (det F → 0/negative) must be rejected by
    setTrialF (the element then fails the step rather than returning sign-flipped σ)."""
    # squash z to near-zero thickness => det F tiny/negative under the affine map
    Fbad = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -0.5]]
    rc = _impose_and_solve(_affine_disp(Fbad), "finite", _mat_fin)
    assert rc != 0, "inverted det F was not rejected"


def test_finite_rc_database_roundtrip():
    """A committed finite damaged state survives sendSelf/recvSelf through the FE database."""
    lam = 1.008
    Fm = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)])
    u = _affine_disp(Fm.tolist())

    def build():
        assert _impose_and_solve(u, "finite", _mat_fin) == 0

    database_roundtrip(build, probe_nodes=[7], ndf=3)
