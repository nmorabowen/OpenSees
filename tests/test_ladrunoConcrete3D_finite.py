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
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


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
def test_finite_objectivity_through_damage():
    # a non-symmetric stretch driven WELL past tensile onset so damage develops (ωt>0), then rotate it
    F = np.array([[1.10, 0.03, 0.02],
                  [0.0, 0.97, 0.015],
                  [0.0, 0.0, 0.98]])
    Q = _rot([0.2, 0.5, -0.84], 1.1)                 # ~63°

    assert _impose_and_solve(_affine_disp(F.tolist()), "finite") == 0
    sF = _voigt_to_mat(_gp_cauchy()[0])
    wt = list(ops.eleResponse(1, "material", 1, "damage"))[0]
    assert wt > 0.01, f"objectivity test is vacuous — no damage developed (ωt={wt})"

    assert _impose_and_solve(_affine_disp((Q @ F).tolist()), "finite") == 0
    sQF = _voigt_to_mat(_gp_cauchy()[0])

    pushed = Q @ sF @ Q.T
    err = np.abs(sQF - pushed).max() / max(np.abs(sQF).max(), 1.0)
    assert err < 1.0e-8, f"NOT objective through damage: ‖σ(QF)−Qσ(F)Qᵀ‖={err:.3e} (isotropic ⇒ should be exact)"


# --------------------------------------------------------------------------- #
# 4. the LogStrain SEAM through the damage law: uniaxial finite stretch into damage →
#    Gauss-point Cauchy == small-strain stress at the Hencky strain ½lnB, pushed by /J.
# --------------------------------------------------------------------------- #
def test_finite_uniaxial_stretch_matches_hencky_seam():
    lam = 1.012                                       # log axial ≈1.2e-2 ≫ onset ε0≈1e-4 ⇒ tensile damage
    lat = 1.0 / math.sqrt(lam)                        # a prescribed (near-isochoric) F — exact value irrelevant
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
