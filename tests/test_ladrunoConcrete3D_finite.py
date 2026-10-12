"""LadrunoConcrete3D under finite strain via `nDMaterial LogStrain` (P4).

    nDMaterial LadrunoConcrete3D 1 E nu fc ft Gf Gc      # small-strain inner
    nDMaterial LogStrain 2 1                              # Hencky finite-strain lift
    element LadrunoBrick 1 ...nodes... 2 -geom finite     # F -> setTrialF

The 3D material lifts to finite strain "for free" through the generic Hencky wrapper
(LogStrainNDMaterial 33010): the element passes F, LogStrain forms b=F Fᵀ, the
logarithmic strain εᵉ=½ln b, feeds it to the UNCHANGED small-strain LadrunoConcrete3D
as the trial strain, treats the returned stress as the Kirchhoff stress τ, and reports
the Cauchy stress σ=τ/J. The kernel exposes `sigEffImplicit` (the UNDAMAGED effective
stress) precisely so LogStrain recovers bᵉ from the non-degraded stress — sidestepping
the (1−d) bᵉ-drift that forces a damage inner like LadrunoRCConcrete to use a native
FiniteStrainNDMaterial subclass; for THIS material the generic wrapper is exact.

The headline: LadrunoConcrete3D is an ISOTROPIC plastic-damage law (Menétrey–Willam +
two SCALAR damage variables ωt/ωc, NO tensorial/kinematic internal variable), so it is
**fully objective under arbitrary rotation** — σ(QF)=Qσ(F)Qᵀ holds EXACTLY. This is the
clean contrast with LadrunoJ2Finite, whose Chaboche backstress does NOT co-rotate in the
LogStrain wrapper (the de Souza Neto §14.11 boundary, pinned xfail there).

Gates (self-referential — the small-strain material IS the reference; no numpy oracle):
  * reduce-to-small-strain: a tiny stretch reproduces the small-strain LadrunoConcrete3D
    stress (Hencky ε≈engineering ε, J≈1);
  * rigid rotation of the unloaded element is stress-free;
  * OBJECTIVITY: a damaged state rotated by Q gives σ(QF)=Qσ(F)Qᵀ (exact, isotropic);
  * uniaxial finite stretch into DAMAGE: the Gauss-point Cauchy stress == the small-strain
    LadrunoConcrete3D stress evaluated at the Hencky strain ½ln b, pushed back by /J (the
    LogStrain seam through the damage law).

Requires a local/CI build (`from _testbed import ops`).
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# a deliberately irregular hex (no GP degeneracy) — same connectivity family as the J2 finite test
_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

# oracle-gate concrete (compression NEGATIVE; fc/ft positive magnitudes)
_E, _NU, _FC, _FT, _GF, _GC = 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0


def _mat(tag):
    ops.nDMaterial("LadrunoConcrete3D", tag, _E, _NU, _FC, _FT, _GF, _GC)


def _build(geom):
    """Hex with the small-strain material (geom='linear') or LogStrain-wrapped (geom='finite')."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _mat(1)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    ops.element("LadrunoBrick", 1, *_CONN, mtag, "-formulation", "std", "-geom", geom)
    return mtag


def _affine_disp(Fbar):
    """Nodal displacements imposing a homogeneous deformation gradient Fbar (d = (F−I)·X)."""
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        X = np.array([x, y, z])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = (np.asarray(Fbar) - I) @ X
    return u


def _impose_and_solve(u, geom):
    _build(geom)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")                       # CDPM2 tangent is NON-SYMMETRIC
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    # Solver ladder (WP concrete3d-hang-diagnosis #877 follow-up): the whole prescribed deformation is imposed
    # in ONE load step, so plain Newton's basin of attraction decides whether the (physically fine) state is
    # reached -- test_finite_objectivity_through_damage documents this ("either the F or the Q@F state,
    # unpredictably"), and the return-map changes of #877 moved it across the basin boundary. Retrying the SAME
    # step with a more globalized algorithm changes the path, not the physics or the tolerance.
    rc = -1
    for alg in (("Newton",), ("NewtonLineSearch", "-type", "Bisection"), ("KrylovNewton",),
                ("Newton", "-initial")):
        ops.algorithm(*alg)
        rc = ops.analyze(1)
        if rc == 0:
            break
    return rc


def _gp_cauchy(geom="finite"):
    """The 8 Gauss points' Cauchy stress (6-Voigt each)."""
    return np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)


def _rot(ax, ang):
    """Rotation matrix: angle `ang` about unit axis `ax` (Rodrigues)."""
    ax = np.asarray(ax, float); ax = ax / np.linalg.norm(ax)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1.0 - math.cos(ang)) * (K @ K)


def _voigt_to_mat(s):
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def _small_strain_cauchy(F):
    """Reference: the LogStrain seam = small-strain LadrunoConcrete3D evaluated at the Hencky strain
    εᵉ=½ln(F Fᵀ) (treated as Kirchhoff τ), pushed to Cauchy σ=τ/J. Driven via NDTest on a fresh
    small-strain material — the same binary, so it certifies the wrapper's Hencky seam end-to-end."""
    F = np.asarray(F, float)
    B = F @ F.T
    w, V = np.linalg.eigh(B)
    H = V @ np.diag(0.5 * np.log(w)) @ V.T                 # ½ ln b (symmetric)
    eng = [H[0, 0], H[1, 1], H[2, 2], 2.0 * H[0, 1], 2.0 * H[1, 2], 2.0 * H[0, 2]]   # engineering shear ×2
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    _mat(1)
    ops.NDTest("SetStrain", 1, *eng); ops.NDTest("CommitState", 1)
    tau = np.array(ops.NDTest("GetStress", 1), dtype=float)   # treated as Kirchhoff
    return tau / float(np.linalg.det(F))


# --------------------------------------------------------------------------- #
# 1. reduce-to-small-strain: a tiny stretch reproduces the small-strain material
# --------------------------------------------------------------------------- #
def test_finite_reduces_to_small_strain():
    e = 4.0e-5                                       # < onset ε0=ft/E≈1e-4 ⇒ elastic, ½lnB≈ε, J≈1
    Fm = np.diag([1.0 + e, 1.0, 1.0])
    assert _impose_and_solve(_affine_disp(Fm.tolist()), "finite") == 0
    s_fin = _gp_cauchy()[0]
    # small-strain reference at the SAME engineering strain (no log/J correction needed at this size)
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3); _mat(1)
    ops.NDTest("SetStrain", 1, e, 0, 0, 0, 0, 0); ops.NDTest("CommitState", 1)
    s_ss = np.array(ops.NDTest("GetStress", 1), dtype=float)
    assert np.allclose(s_fin, s_ss, rtol=2.0e-3, atol=1.0e-3), (
        f"finite at tiny strain {s_fin} != small-strain {s_ss}")


# --------------------------------------------------------------------------- #
# 2. rigid rotation of the unloaded element is stress-free
# --------------------------------------------------------------------------- #
def test_finite_rigid_rotation_stress_free():
    Q = _rot([0.3, -0.7, 0.65], 0.9)                 # ~51° rotation, no stretch
    assert _impose_and_solve(_affine_disp(Q.tolist()), "finite") == 0
    s = _gp_cauchy()
    assert np.abs(s).max() < 1.0e-6 * _FC, f"rigid rotation not stress-free: max|σ|={np.abs(s).max():.3e}"


# --------------------------------------------------------------------------- #
# 3. OBJECTIVITY (the headline): a DAMAGED state rotated by Q gives σ(QF)=Qσ(F)Qᵀ
#    EXACTLY — LadrunoConcrete3D is isotropic (scalar damage, no co-rotating internal var),
#    so there is NO de Souza Neto §14.11 boundary (contrast LadrunoJ2Finite's backstress).
# --------------------------------------------------------------------------- #
def _objectivity_err(scale):
    """Solve the damaged stretch F = I + scale*(F_full - I) and its rotation Q@F; return (err, omega_t) with
    err = ||sigma(QF) - Q sigma(F) Q^T||_max / max(||sigma(QF)||_max, 1)."""
    F_full = np.array([[1.10, 0.03, 0.02],
                        [0.0, 0.97, 0.015],
                        [0.0, 0.0, 0.98]])
    F = np.eye(3) + scale * (F_full - np.eye(3))
    Q = _rot([0.2, 0.5, -0.84], 1.1)                 # ~63 deg
    assert _impose_and_solve(_affine_disp(F.tolist()), "finite") == 0
    sF = _voigt_to_mat(_gp_cauchy()[0])
    wt = list(ops.eleResponse(1, "material", 1, "damage"))[0]
    assert _impose_and_solve(_affine_disp((Q @ F).tolist()), "finite") == 0
    sQF = _voigt_to_mat(_gp_cauchy()[0])
    pushed = Q @ sF @ Q.T
    return np.abs(sQF - pushed).max() / max(np.abs(sQF).max(), 1.0), wt


def test_finite_objectivity_through_damage():
    """sigma(QF) = Q sigma(F) Q^T through damage, EXACTLY (isotropic scalar damage, no co-rotating variable).

    SCALE, derived from measurement (WP concrete3d-hang-diagnosis review M1). The Hencky strains of F and Q@F are
    rotations of each other, equal to ~1e-16 (and the global solve reproduces them to its 1e-10 tolerance), and the
    kernel is exactly objective for one and the same branch of the deterministic sub-incrementation map (kernel probe,
    virgin state, scales 0.005..0.03: ||sigma(QF) - Q sigma(F) Q^T|| = 1e-15). The two solves can only differ if a 1e-12 change of
    the strain changes the BRANCH: the piece count n = ceil(f_tr/0.3) (the trial sits >= 0.13 from an integer at these
    scales, so no), or the failure ladder n -> 2n -> 4n. The ladder flickers because ONE piece of the chain -- the first
    plastic tensile piece, sigma_xx ~ 2.8 = 0.92 ft, kappa_p ~ 0.05, f_tr = 0.088 -- fails its direct return in ~30 % of
    1e-12 perturbations (Newton at the edge of convergence in the locally indefinite tension regime), so a chain of n pieces
    escalates with probability ~ 1 - (1 - p)^n. Measured by perturbing the Hencky strain by 1e-12, 40 draws each:
        scale 0.01 (n = 9):  ladder level differs in  0/40, stress changes <= 4e-9   -> branch-safe
        scale 0.02 (n = 42): differs in 12/40, stress changes up to 4.6e-3
        scale 0.03 (n = 64, saturated, subInfo 128): differs in 29/40, up to 5.0e-3
        scale 0.10 (saturated, subInfo 256): F and QF land on kappa_p 58.7 vs 54.2, error 3.2e-2 (below)
    Scale 0.01 is the largest of these with no observed flip and still damages (omega_t = 0.29); the earlier 0.03 passed
    only by luck of the ladder level (0.04, 0.06 likewise; 0.05, 0.10 did not)."""
    err, wt = _objectivity_err(0.01)
    assert wt > 0.01, f"objectivity test is vacuous — no damage developed (ωt={wt})"
    assert err < 1.0e-8, f"NOT objective through damage: ‖σ(QF)−Qσ(F)Qᵀ‖={err:.3e} (isotropic ⇒ should be exact)"


def test_finite_objectivity_saturated_branch_jump_bound():
    """Where the trial saturates the piece count (scale 0.10: f_tr/c ~ 5e3 >> nmax = 64) F and Q@F may land on different
    ladder levels of the sub-incremented return -- two consistent integrations of the same increment, not a violation of
    isotropy. The difference is BOUNDED by the branch jump of the deterministic map: measured 3.2e-2 in the normalized
    stress (kappa_p 58.7 vs 54.2, omega_t 0.394 vs 0.375); the review measured 0.4-1.5 % of sigma_eff at n boundaries
    (0.06-0.40 MPa on 15-41 MPa) and 1-2 % of sigma_eff with kappa_p jumps of 1.1-2.5 at ladder switches near first
    cracking. Bound 6e-2 = 2x the measured stress difference."""
    err, wt = _objectivity_err(0.10)
    assert wt > 0.01
    assert err < 6.0e-2, f"branch jump exceeds the documented bound: {err:.3e}"


# --------------------------------------------------------------------------- #
# 4. the LogStrain SEAM through the damage law: uniaxial finite stretch into damage →
#    Gauss-point Cauchy == small-strain stress at the Hencky strain ½lnB, pushed by /J.
# --------------------------------------------------------------------------- #
def test_finite_uniaxial_stretch_matches_hencky_seam():
    # Under the current CDPM2-default kernel (B1 flow potential, B2 damage drive), an ISOCHORIC
    # uniaxial stretch (lat = 1/sqrt(lam), the pre-B1/B2 choice here) has zero trace Hencky strain,
    # so the trial state is purely deviatoric (sigma_V = 0) and the M-W return map for this triaxiality
    # lands entirely in the COMPRESSIVE regime (ωc grows, ωt stays exactly 0 — verified numerically:
    # at lam=1.012 the Gauss-point Cauchy is [-1.8, -44.7, -44.7] MPa, all compressive). "exact value
    # irrelevant" no longer holds: pin lat = 1.0 (no lateral contraction) so the Hencky strain keeps a
    # net tensile trace and the return map actually damages the TENSILE side (ωt), matching this test's
    # intent (the seam through TENSILE damage; ωc growth is already covered elsewhere).
    lam = 1.0003                                      # log axial ≈3e-4 ≫ onset ε0≈1e-4 ⇒ tensile damage
    lat = 1.0
    Fm = np.diag([lam, lat, lat])
    assert _impose_and_solve(_affine_disp(Fm.tolist()), "finite") == 0
    s = _gp_cauchy()
    s0 = s[0]
    assert np.abs(s - s0).max() <= 1.0e-6 * max(np.abs(s0).max(), 1.0), (
        "uniaxial finite patch: Cauchy not constant across Gauss points")
    # damage must actually be active, else this is just an elastic seam check
    assert list(ops.eleResponse(1, "material", 1, "damage"))[0] > 0.01, "no damage — increase lam"

    s_ref = _small_strain_cauchy(Fm)
    tol = 1.0e-5 * max(np.abs(s_ref).max(), 1.0)
    assert np.allclose(s0, s_ref, rtol=1.0e-5, atol=tol), (
        f"finite GP Cauchy {s0} != Hencky-seam reference {s_ref}")
