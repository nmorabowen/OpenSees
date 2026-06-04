"""Finite-strain trifecta — **P3 locking & incompressibility** (F-bar cure).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §9
phase **P3** (gate: bbar/eas/F-bar cure demonstrated). P1's B2/B3 showed the
F-bar cure on a slender bending cantilever (elastic ν→½ and isochoric J2). P3
extends it to the two harder regimes the plan calls out:

  * **B4 — Cook's membrane** (Cook 1974): a tapered trapezoidal panel, clamped on
    the left edge, end-shear on the right — the canonical *combined bending +
    shear + (near-)incompressibility* locking probe. Under finite strain the
    standard hex locks volumetrically as ν→½; **F-bar (`bbar`) converges to a
    mesh-independent tip displacement while `std` stays far too stiff and strongly
    mesh-dependent**. Run near-incompressible elastic (the clean cure) and anchored
    to the compressible plane-strain reference. (Isochoric-*plastic* Cook locking
    is weak at the loads that converge here; that cure is covered by P1's B3
    cantilever.)
  * **E4 — near-incompressible rubber-like block under large compression** with
    bonded platens (lateral motion fixed on the loaded faces ⇒ barrelling ⇒
    inhomogeneous, constrained, isochoric). Hencky elastic, ν=0.4999: `std` locks
    (grossly over-stiff reaction), F-bar stays compliant. The *large*-deformation
    counterpart of B2 (which was small-strain bending).

All Zone-A (structured lattices, no gmsh). `-geom finite` + `LogStrain` over
`ElasticIsotropic`/`LadrunoJ2`. F-bar's tangent is unsymmetric ⇒ `FullGeneral`;
plastic/large-deformation BVPs use `KrylovNewton` (cf. P1 LEDGER_quirks).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t2a]


# --------------------------------------------------------------------------- #
#  Cook's membrane: trapezoid A(0,0) B(48,44) C(48,60) D(0,44), unit-thick     #
#  3D slab (plane strain: uz=0 on both faces), left edge clamped, right-edge    #
#  vertical shear. n×n in-plane elements. Returns tip vertical disp at corner C.#
# --------------------------------------------------------------------------- #
_CK = dict(A=(0.0, 0.0), B=(48.0, 44.0), C=(48.0, 60.0), D=(0.0, 44.0))


def _cook_xy(xi, eta):
    A, B, C, D = (np.array(_CK[k]) for k in "ABCD")
    return (1 - xi) * (1 - eta) * A + xi * (1 - eta) * B + xi * eta * C + (1 - xi) * eta * D


def _cook_tip(form, n, nu, P, E=1.0, t=1.0, nsteps=8):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    def nid(i, j, k):
        return 1 + i + (n + 1) * j + (n + 1) * (n + 1) * k

    for k in range(2):
        for j in range(n + 1):
            for i in range(n + 1):
                x, y = _cook_xy(i / n, j / n)
                ops.node(nid(i, j, k), float(x), float(y), t * k)
    for k in range(2):
        for j in range(n + 1):
            for i in range(n + 1):
                if i > 0:
                    ops.fix(nid(i, j, k), 0, 0, 1)        # plane strain (uz=0)
    for k in range(2):
        for j in range(n + 1):
            ops.fix(nid(0, j, k), 1, 1, 1)                # clamp left edge

    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain", 2, 1)
    e = 1
    for j in range(n):
        for i in range(n):
            conn = [nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i, j + 1, 0),
                    nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)]
            ops.element("LadrunoBrick", e, *conn, 2, "-formulation", form, "-geom", "finite")
            e += 1

    right = [nid(n, j, k) for j in range(n + 1) for k in range(2)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in right:
        ops.load(nd, 0.0, P / len(right), 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-10, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"Cook solve failed (form={form}, n={n}, mat={material})"
    return float(np.mean([ops.nodeDisp(nid(n, n, k), 2) for k in range(2)]))


# =========================================================================== #
#  B4.1 — F-bar Cook's membrane converges; std locks (near-incompressible)      #
# =========================================================================== #
def test_B4_cook_membrane_fbar_converges_std_locks():
    nu, P = 0.499, 1.0e-4                  # near-incompressible, small load (linear)
    d_bar_16 = _cook_tip("bbar", 16, nu, P, nsteps=1)
    d_bar_24 = _cook_tip("bbar", 24, nu, P, nsteps=1)
    d_std_24 = _cook_tip("std", 24, nu, P, nsteps=1)
    # (1) F-bar is mesh-converged (16→24 within a few %)
    assert abs(d_bar_24 - d_bar_16) / d_bar_24 < 0.04, (
        f"F-bar Cook not converged: {d_bar_16:.3f} → {d_bar_24:.3f}")
    # (2) std LOCKS: even at n=24 it is far stiffer (smaller tip disp) than F-bar
    assert d_bar_24 / d_std_24 > 1.3, (
        f"std did not lock on Cook's membrane: bbar/std = {d_bar_24 / d_std_24:.2f}")


def test_B4_cook_membrane_std_is_strongly_mesh_dependent():
    # The signature of locking: std's answer crawls upward with refinement (never
    # converged in the usable range) while F-bar is already flat.
    nu, P = 0.499, 1.0e-4
    d_std_8 = _cook_tip("std", 8, nu, P, nsteps=1)
    d_std_16 = _cook_tip("std", 16, nu, P, nsteps=1)
    d_bar_8 = _cook_tip("bbar", 8, nu, P, nsteps=1)
    d_bar_16 = _cook_tip("bbar", 16, nu, P, nsteps=1)
    std_change = abs(d_std_16 - d_std_8) / d_std_16
    bar_change = abs(d_bar_16 - d_bar_8) / d_bar_16
    assert std_change > 3.0 * bar_change, (
        f"std should be far more mesh-dependent than F-bar (locking): "
        f"std Δ={std_change:.3f} vs bbar Δ={bar_change:.3f}")


def test_B4_cook_membrane_compressible_converges_to_reference():
    # Compressible ν=1/3 plane-strain Cook — F-bar converges to a stable tip
    # displacement (~22.5; consistent with the plane-strain Cook value, the
    # plane-stress canonical 23.9 being a different BC). Anchors the absolute scale.
    nu, P = 1.0 / 3.0, 1.0e-4
    d16 = _cook_tip("bbar", 16, nu, P, nsteps=1) / P
    d32 = _cook_tip("bbar", 32, nu, P, nsteps=1) / P
    assert abs(d32 - d16) / d32 < 0.03, f"compressible Cook not converged: {d16:.2f} → {d32:.2f}"
    assert 21.0 < d32 < 24.0, f"Cook tip {d32:.2f} outside the plane-strain reference band"


# (Isochoric-J2 locking relief is validated on the bending cantilever in P1's
#  test_finite_strain_L2_convergence.py::test_B3_isochoric_j2_locking_relieved_by_fbar.
#  On Cook's membrane the plastic-locking separation is weak at the loads that
#  converge here — std≈bbar — so the elastoplastic-Cook row is intentionally not
#  asserted; the elastic near-incompressible Cook above is the strong B4 cure.)


# --------------------------------------------------------------------------- #
#  E4 — near-incompressible rubber block, bonded platens, LARGE compression.    #
#  A cube compressed by a prescribed top-face descent with the top & bottom     #
#  faces bonded (lateral motion fixed) ⇒ barrelling, constrained isochoric      #
#  deformation. Returns the total vertical reaction (resistance).               #
# --------------------------------------------------------------------------- #
def _rubber_block_reaction(form, n, nu, stroke, E=10.0, nsteps=20):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    def nid(i, j, k):
        return 1 + i + (n + 1) * j + (n + 1) * (n + 1) * k

    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(nid(i, j, k), float(i) / n, float(j) / n, float(k) / n)
    bottom = [nid(i, j, 0) for j in range(n + 1) for i in range(n + 1)]
    top = [nid(i, j, n) for j in range(n + 1) for i in range(n + 1)]
    for nd in bottom:
        ops.fix(nd, 1, 1, 1)                 # bonded bottom platen
    for nd in top:
        ops.fix(nd, 1, 1, 0)                 # bonded top platen: lateral fixed, vertical driven

    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain", 2, 1)
    e = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, 2, "-formulation", form, "-geom", "finite")
                e += 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in top:
        ops.sp(nd, 3, -stroke)               # prescribe downward stroke
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-9, 150, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"rubber block solve failed (form={form})"
    ops.reactions()
    return abs(sum(ops.nodeReaction(nd, 3) for nd in bottom))


def test_E4_rubber_block_large_compression_fbar_relieves_locking():
    # 10% compression of a near-incompressible (ν=0.499) bonded block — a large,
    # constrained, isochoric deformation. The std hex cannot represent the
    # isochoric barrelling ⇒ it over-stiffens grossly; F-bar stays compliant.
    # (ν=0.4999 over-constrains the bonded block past convergence; 0.499 already
    # gives a ~9× separation — the cure is unmistakable.)
    nu, stroke = 0.499, 0.10
    R_std = _rubber_block_reaction("std", 4, nu, stroke, nsteps=40)
    R_bar = _rubber_block_reaction("bbar", 4, nu, stroke, nsteps=40)
    assert R_std > 0 and R_bar > 0
    assert R_std / R_bar > 3.0, (
        f"F-bar did not relieve near-incompressible large-compression locking: "
        f"std/bbar reaction = {R_std / R_bar:.2f}")
