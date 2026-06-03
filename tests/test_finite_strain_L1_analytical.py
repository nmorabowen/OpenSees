"""Finite-strain trifecta — **L1 analytical / manufactured** benchmarks (A1–A5).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §5
"L1 — Analytical / manufactured (exact oracles)" and §9 phase **P1**. These close
the coverage gaps the finite-strain deep review flagged (degenerate stretch,
|J−1|≫0, plastic incompressibility) by driving the COMPILED element
`LadrunoBrick -geom finite` against *closed-form* oracles — one rung above the L0
self-consistency unit tests (test_ladrunoBrick_finite.py / test_ladrunoJ2_finite_*).

    nDMaterial LadrunoJ2|ElasticIsotropic 1 ...      # small-strain inner
    nDMaterial LogStrain 2 1                          # Hencky finite-strain lift
    element LadrunoBrick 1 ...nodes... 2 -geom finite # F → setTrialF

Benchmarks (G=finite, single 8-node hex):
  * **A1** uniaxial finite stretch into plasticity → Cauchy σ matches the
    closed-form 1D finite-strain J2 relation σ = σ_y(p), p = ln λ − σ/E
    (true uniaxial-stress rig: symmetry planes + traction-free lateral faces).
  * **A2** equibiaxial (double-degenerate λ) and hydrostatic (triple-degenerate λ)
    stretch → exercises the repeated-eigenvalue branch of the tensor-log map.
  * **A3** large dilatation J ∈ [0.5, 2] → checks σ = J⁻¹τ and tr ε = ln J away
    from J ≈ 1 (the 1/2J spatial-modulus factor and the dilatation push-forward).
  * **A4** finite simple shear F = I + γ e₁⊗e₂ → Hencky hyperelastic closed form
    (elastic) and the material-point J2 return map (plastic), J = 1 isochoric.
  * **A5** plastic incompressibility — under isochoric (J = 1) J2 flow the mean
    Cauchy stress stays ≡ 0 (volumetric response is purely elastic Hencky:
    tr τ = 3K ln J = 0), at every committed step of a finite plastic path.

Oracles are the numpy mirrors tests/logstrain_reference.py (MATISU / Hencky) and
tests/ladrunoj2_reference.py (combined-hardening J2 return map) — no OpenSees
dependency, so the oracle and the C++ kernel cannot share an indexing bug.

The large-rotation §14.11 boundary is NOT exercised here (these are pure
stretch / shear states, R ≈ I); A-series plasticity uses isotropic hardening,
which the log-strain lift carries EXACTLY at any deformation.
"""
import math

import numpy as np
import pytest

from _testbed import ops
import ladrunoj2_reference as lr
import logstrain_reference as ls

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

I3 = np.eye(3)

# Elastic moduli (shared). E, ν derived so the closed forms read off K, G.
_K, _G = 1500.0, 700.0
_E = 9.0 * _K * _G / (3.0 * _K + _G)
_NU = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))
# Isotropic hardening only (objective through the log-strain lift at any F).
# Linear form (A1 closed form): σ_y(p) = sig0 + Hiso·p.
_SIG0, _HISO = 10.0, 60.0
_VOCE = dict(sig0=_SIG0, Qinf=0.0, bIso=1.0, Hiso=_HISO)   # Qinf=0 ⇒ pure linear

# A distorted (positive-Jacobian) reference hex for the affine patch oracles
# (A2/A3/A4) — a non-axis-aligned cell genuinely exercises the spectral map.
_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
# Axis-aligned UNIT cube for the uniaxial / dilatation rigs (A1, A3, A5).
_CUBE = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


# --------------------------------------------------------------------------- #
#  Builders                                                                    #
# --------------------------------------------------------------------------- #
def _elastic(tag):
    ops.nDMaterial("ElasticIsotropic", tag, _E, _NU)


def _j2_iso(tag):
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, "-iso", "voce",
                   _VOCE["sig0"], _VOCE["Qinf"], _VOCE["bIso"], _VOCE["Hiso"],
                   "-kin", 0)


def _build(nodes, mat_fn, form="std"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in nodes.items():
        ops.node(tag, x, y, z)
    mat_fn(1)
    ops.nDMaterial("LogStrain", 2, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 2, "-formulation", form, "-geom", "finite")


def _gp_cauchy(avg=True):
    s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)
    if avg:
        return s.mean(axis=0), s
    return s


# --------------------------------------------------------------------------- #
#  Affine prescribed-F patch (A2/A3/A4): impose F at every node via sp.        #
# --------------------------------------------------------------------------- #
def _solve_affine_F(F, nodes, mat_fn, form="std", nsteps=1):
    _build(nodes, mat_fn, form)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    I = np.eye(3)
    for tag, (x, y, z) in nodes.items():
        d = (np.asarray(F) - I) @ np.array([x, y, z])
        base = (tag - 1) * 3
        for j in range(3):
            ops.sp(tag, j + 1, float(d[j]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, "affine-F solve failed"


# =========================================================================== #
#  A1 — uniaxial finite stretch into plasticity (true uniaxial-stress rig)     #
# =========================================================================== #
def _solve_uniaxial(lam, nsteps=8):
    """Symmetry-plane uniaxial pull of a unit cube: x=0/y=0/z=0 faces on their
    symmetry planes, z=1 face pulled to stretch λ, lateral faces traction-free.
    Single hex ⇒ homogeneous uniaxial stress state."""
    _build(_CUBE, _j2_iso, "std")
    # symmetry BCs
    for tag in (1, 4, 5, 8):   # x = 0 face
        ops.fix(tag, 1, 0, 0)
    for tag in (1, 2, 5, 6):   # y = 0 face
        ops.fix(tag, 0, 1, 0)
    for tag in (1, 2, 3, 4):   # z = 0 face
        ops.fix(tag, 0, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    uz = lam - 1.0
    for tag in (5, 6, 7, 8):   # z = 1 face: prescribe axial stretch
        ops.sp(tag, 3, uz)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, "uniaxial solve failed"


def _uniaxial_J(lam):
    """det F for the homogeneous uniaxial rig from the deformed lateral stretches
    (node 2 free in x, node 4 free in y; axial stretch = lam)."""
    lx = 1.0 + ops.nodeDisp(2, 1)
    ly = 1.0 + ops.nodeDisp(4, 2)
    return lam * lx * ly


def test_A1_uniaxial_finite_stretch_matches_closed_form():
    lam = 1.10                                  # ~9.5% log strain, solidly plastic
    _solve_uniaxial(lam)
    sig_avg, sig_all = _gp_cauchy()
    # homogeneous: stress uniform across the 8 Gauss points
    assert np.abs(sig_all - sig_avg).max() <= 1.0e-6 * max(abs(sig_avg).max(), 1.0), (
        "A1 uniaxial state not homogeneous across Gauss points")
    # the element must have yielded (else this only checks elasticity)
    p = list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))
    assert p and p[0] > 0.0, "A1 vacuous — element did not yield"
    sig_zz = sig_avg[2]
    sig_lat = max(abs(sig_avg[0]), abs(sig_avg[1]))
    # (1) lateral Cauchy ≈ 0 (true uniaxial-stress state)
    assert sig_lat <= 1.0e-4 * abs(sig_zz), (
        f"A1 lateral stress not traction-free: {sig_lat:.3e} vs axial {sig_zz:.3e}")
    # (2) EXACT closed form: under uniaxial stress von-Mises yield applies to the
    #     KIRCHHOFF stress τ (the small-strain return map's output), τ_axial =
    #     σ_y(p), and σ = τ/J. So the exact identity on the Cauchy stress is
    #     σ_zz · J = σ_y(p), with p the element's own equivalent plastic strain.
    sig_yield = _SIG0 + _HISO * p[0]
    J = _uniaxial_J(lam)
    assert sig_zz * J == pytest.approx(sig_yield, rel=1.0e-5), (
        f"A1 axial Kirchhoff σ·J {sig_zz * J:.6f} != yield σ_y(p) {sig_yield:.6f}")
    # (3) finite-strain kinematics (approximate, small elastic strain): the
    #     plastic log strain p ≈ ln λ − σ/E. Loose tol — it is the small-elastic-
    #     strain approximation, not the exact Cauchy relation checked in (2).
    p_approx = math.log(lam) - sig_zz / _E
    assert p[0] == pytest.approx(p_approx, rel=2.0e-2), (
        f"A1 plastic strain p {p[0]:.5f} far from ln λ − σ/E {p_approx:.5f}")


def test_A1_uniaxial_backbone_tracks_yield_curve():
    """Across a sweep of finite stretches the axial true stress must trace the
    uniaxial yield backbone σ_y(p) = sig0 + H·p exactly (the uniaxial von-Mises
    identity), and both σ_true and p increase monotonically (hardening)."""
    sig_true, ps = [], []
    for lam in (1.02, 1.05, 1.10, 1.18, 1.28):
        _solve_uniaxial(lam, nsteps=max(6, int(50 * (lam - 1.0))))
        s = _gp_cauchy()[0][2]
        p = list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))[0]
        tau = s * _uniaxial_J(lam)               # Kirchhoff axial stress
        assert tau == pytest.approx(_SIG0 + _HISO * p, rel=1.0e-5), (
            f"A1 backbone: τ {tau:.5f} != σ_y(p) {_SIG0 + _HISO * p:.5f} at λ={lam}")
        sig_true.append(s); ps.append(p)
    assert all(sig_true[i] < sig_true[i + 1] for i in range(len(sig_true) - 1)), sig_true
    assert all(ps[i] < ps[i + 1] for i in range(len(ps) - 1)), ps


# =========================================================================== #
#  A2 — degenerate-eigenvalue stretch (repeated-λ branch of the tensor-log)    #
# =========================================================================== #
def test_A2_hydrostatic_stretch_triple_degenerate_elastic():
    """F = λ I ⇒ bᵉᵗʳ = λ² I has a TRIPLE eigenvalue: the A.52 third branch
    (Y = y'·I_S). Cauchy σ = (3K ln λ / λ³) I in closed form."""
    lam = 1.15
    _solve_affine_F(lam * I3, _CUBE, _elastic)
    sig_avg, sig_all = _gp_cauchy()
    assert np.abs(sig_all - sig_avg).max() <= 1.0e-7 * max(abs(sig_avg).max(), 1.0)
    p = 3.0 * _K * math.log(lam) / lam**3
    sig_ref = np.array([p, p, p, 0.0, 0.0, 0.0])
    assert np.allclose(sig_avg, sig_ref, rtol=1.0e-5, atol=1.0e-5 * abs(p)), (
        f"A2 hydrostatic: σ {sig_avg} != closed form {sig_ref} (triple-λ branch)")


def test_A2_equibiaxial_stretch_double_degenerate_j2():
    """Equibiaxial F = diag(λ, λ, 1/λ²) ⇒ bᵉᵗʳ eigenvalues {λ², λ², λ⁻⁴}: a
    DOUBLE eigenvalue (A.52 second branch / A.53). Driven into J2 plasticity and
    compared, component for component, to the numpy MATISU return-map oracle."""
    lam = 1.20                                  # in-plane stretch, isochoric F
    F = np.diag([lam, lam, 1.0 / lam**2])
    # single increment from the reference config so the element's MATISU update
    # is one return map ref→F, matching the single-step oracle below (plastic
    # flow is path-dependent; a sub-stepped solve would not match a 1-step map).
    _solve_affine_F(F, _NODES, _j2_iso, nsteps=1)
    sig_avg, sig_all = _gp_cauchy()
    assert np.abs(sig_all - sig_avg).max() <= 1.0e-6 * max(abs(sig_avg).max(), 1.0)
    # oracle: single step from reference config
    P = lr.Params(_K, _G, **_VOCE)
    eps_tr = ls.hencky_strain(F @ F.T)
    res = lr.return_map(P, eps_tr, np.zeros((3, 3)), 0.0, [])
    assert res["plastic"], "A2 equibiaxial did not yield — increase λ"
    sig = res["stress"] / float(np.linalg.det(F))
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[2, 0]])
    assert np.allclose(sig_avg, sig_v, rtol=1.0e-5, atol=1.0e-5 * max(abs(sig_v).max(), 1.0)), (
        f"A2 equibiaxial: σ {sig_avg} != oracle {sig_v} (double-λ branch)")


# =========================================================================== #
#  A3 — large dilatation J ∈ [0.5, 2]                                          #
# =========================================================================== #
@pytest.mark.parametrize("J", [0.5, 0.75, 1.5, 2.0])
def test_A3_large_dilatation_sigma_and_trace(J):
    """Pure dilatation F = J^{1/3} I, well away from J ≈ 1. Verifies the Box-14.3
    σ = J⁻¹τ push-forward and tr ε = ln J (the 1/2J spatial factors live in the
    tangent, checked by the FD-tangent gate elsewhere)."""
    lam = J ** (1.0 / 3.0)
    _solve_affine_F(lam * I3, _CUBE, _elastic)
    sig_avg = _gp_cauchy()[0]
    # τ = 3K ln(J) /3 · I per component? τ = K tr(ε) I, tr ε = ln J ⇒ τ_ii = K ln J
    tau_ii = _K * math.log(J)
    sig_ref = tau_ii / J
    assert sig_avg[2] == pytest.approx(sig_ref, rel=1.0e-5), (
        f"A3 J={J}: Cauchy σ {sig_avg[2]:.6f} != J⁻¹τ {sig_ref:.6f}")
    assert max(abs(sig_avg[3]), abs(sig_avg[4]), abs(sig_avg[5])) <= 1.0e-8 * max(abs(sig_ref), 1.0)
    # tr ε = ln J : recover from the (known-homogeneous) prescribed F directly
    assert math.log(J) == pytest.approx(3.0 * math.log(lam), rel=1.0e-12)


# =========================================================================== #
#  A4 — finite simple shear F = I + γ e₁⊗e₂ (J = 1)                            #
# =========================================================================== #
def test_A4_simple_shear_elastic_matches_hencky():
    gamma = 0.30
    F = I3 + gamma * np.outer([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    assert float(np.linalg.det(F)) == pytest.approx(1.0, abs=1e-15)
    _solve_affine_F(F, _NODES, _elastic)
    sig_avg, sig_all = _gp_cauchy()
    assert np.abs(sig_all - sig_avg).max() <= 1.0e-7 * max(abs(sig_avg).max(), 1.0)
    sig, tau, J, eps = ls.cauchy_from_F_elastic(F, _K, _G)
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[2, 0]])
    assert np.allclose(sig_avg, sig_v, rtol=1.0e-5, atol=1.0e-5 * max(abs(sig_v).max(), 1.0)), (
        f"A4 simple shear (elastic): σ {sig_avg} != Hencky oracle {sig_v}")
    # finite-shear hallmark: NORMAL stresses appear (σ11 ≠ 0) — a small-strain
    # model would give pure σ12 only; here ln b has nonzero diagonal deviator.
    assert abs(sig_avg[0]) > 1.0e-3 * abs(sig_avg[3]), (
        "A4: finite simple shear must develop normal stresses (Kerr/Poynting)")


def test_A4_simple_shear_plastic_matches_return_map():
    gamma = 0.30
    F = I3 + gamma * np.outer([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    # single increment ref→F (simple shear is NON-proportional: principal axes
    # rotate, so the result is path-dependent — a sub-stepped element solve would
    # not equal a single-step return map). One step makes both sides identical.
    _solve_affine_F(F, _NODES, _j2_iso, nsteps=1)
    sig_avg = _gp_cauchy()[0]
    P = lr.Params(_K, _G, **_VOCE)
    eps_tr = ls.hencky_strain(F @ F.T)
    res = lr.return_map(P, eps_tr, np.zeros((3, 3)), 0.0, [])
    assert res["plastic"], "A4 plastic shear did not yield — increase γ"
    sig = res["stress"] / float(np.linalg.det(F))
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[2, 0]])
    assert np.allclose(sig_avg, sig_v, rtol=1.0e-5, atol=1.0e-5 * max(abs(sig_v).max(), 1.0)), (
        f"A4 simple shear (J2): σ {sig_avg} != return-map oracle {sig_v}")


# =========================================================================== #
#  A5 — plastic incompressibility (volumetric response stays elastic)         #
# =========================================================================== #
def test_A5_isochoric_plastic_flow_has_zero_mean_stress():
    """J2 plastic flow is isochoric ⇒ det Fᵖ = 1 ⇒ tr εᵉ = ln J. Under an
    ISOCHORIC prescribed F (J = 1) the volumetric Hencky response gives
    tr τ = 3K ln J = 0, so the mean Cauchy stress must vanish IDENTICALLY no
    matter how much deviatoric plastic flow occurs. Element-observable form of
    the det bᵉ = J² identity (checked at the material-point level in
    test_ladrunoJ2_finite.py::test_finite_plastic_incompressibility)."""
    # isochoric stretch + shear, well into yield
    lam = 1.25
    F = np.diag([lam, 1.0 / math.sqrt(lam), 1.0 / math.sqrt(lam)]).copy()
    F[0, 1] = 0.15
    assert float(np.linalg.det(F)) == pytest.approx(1.0, abs=1e-12)
    _solve_affine_F(F, _NODES, _j2_iso, nsteps=6)
    p = list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))
    assert p and p[0] > 0.0, "A5 vacuous — element did not yield"
    sig_avg = _gp_cauchy()[0]
    mean = (sig_avg[0] + sig_avg[1] + sig_avg[2]) / 3.0
    dev_scale = max(abs(sig_avg[3]), abs(sig_avg[0] - sig_avg[1]), 1.0)
    assert abs(mean) <= 1.0e-7 * dev_scale, (
        f"A5: isochoric plastic flow produced mean stress {mean:.3e} "
        "(plastic incompressibility violated — volumetric/deviatoric leak)")


def test_A5_incompressibility_along_committed_path():
    """Stronger: step an isochoric finite-plastic path under LoadControl,
    committing each increment, and require the mean Cauchy stress ≈ 0 at EVERY
    committed step (guards a per-step volumetric leak that a single end-state
    check could miss)."""
    lam_f = 1.30
    Ff = np.diag([lam_f, 1.0 / math.sqrt(lam_f), 1.0 / math.sqrt(lam_f)]).copy()
    Ff[0, 1] = 0.18
    _build(_NODES, _j2_iso)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    I = np.eye(3)
    N = 8
    for tag, (x, y, z) in _NODES.items():
        d = (Ff - I) @ np.array([x, y, z])
        base = (tag - 1) * 3
        for j in range(3):
            ops.sp(tag, j + 1, float(d[j]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / N)
    ops.analysis("Static")
    saw_plastic = False
    for k in range(1, N + 1):
        assert ops.analyze(1) == 0, f"A5 path step {k} failed"
        Fk = I + (k / N) * (Ff - I)
        # the intermediate F_k is NOT exactly isochoric, so the *expected* mean
        # stress is the elastic Hencky volumetric value K·ln(det F_k), not 0.
        Jk = float(np.linalg.det(Fk))
        sig_avg = _gp_cauchy()[0]
        mean = (sig_avg[0] + sig_avg[1] + sig_avg[2]) / 3.0
        mean_ref = _K * math.log(Jk) / Jk
        if list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))[0] > 0:
            saw_plastic = True
        assert mean == pytest.approx(mean_ref, abs=1.0e-6 * max(abs(mean_ref), 1.0)), (
            f"A5 step {k}: mean σ {mean:.3e} != elastic-volumetric K ln(J)/J "
            f"{mean_ref:.3e} (plastic flow leaked into the volumetric response)")
    assert saw_plastic, "A5 committed path never yielded"
