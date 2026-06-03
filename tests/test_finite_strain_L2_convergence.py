"""Finite-strain trifecta — **L2 convergence & locking** benchmarks (B1–B3).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §5
"L2 — Convergence & locking" and §9 phase **P1**. Where L1 (the A-series) pins
the *constitutive* answer at a point, L2 pins the *discretization*: that the
8-node hex converges at its theoretical rate and that the volumetric-locking
cures (F-bar) actually relieve locking under finite strain — both elastic
(ν→½) and plastic (isochoric J2).

  * **B1** h-convergence — a smooth finite-strain elastic bending cantilever,
    refined through the thickness, converges monotonically toward a
    Richardson-extrapolated tip deflection (solution verification).
  * **B2** near-incompressible elastic block (ν = 0.4999) — `std` (plain F) locks
    (tip deflection collapses), `bbar` (F-bar) stays compliant. The locking
    SIGNATURE is reference-free: the cross-ν trend and the bbar/std ratio.
  * **B3** isochoric J2 block — plastic incompressibility (det Fᵖ = 1) locks the
    standard hex in bending exactly as ν→½ does elastically; F-bar relieves it.
    This is the *plastic* counterpart of B2 and the review-flagged gap that no
    existing test covers.

All Zone-A (structured lattices, no gmsh). Material: `LogStrain` over
`ElasticIsotropic` (B1/B2) or `LadrunoJ2` isotropic hardening (B3). F-bar's
tangent is unsymmetric ⇒ a full (`FullGeneral`) solver throughout.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t2a]

_E = 200.0
_SIG0, _HISO = 0.40, 10.0          # J2 (B3): yield below bending stress, mild hardening


# --------------------------------------------------------------------------- #
#  A 1×1×nz column of LadrunoBricks, base clamped, transverse tip shear.       #
#  Returns the mean tip transverse (x) displacement.                           #
# --------------------------------------------------------------------------- #
def _cantilever(form, nu, nz, tip_force, material="elastic", nsteps=4):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    def nid(i, j, k):
        return 1 + i + 2 * j + 4 * k

    for k in range(nz + 1):
        for j in range(2):
            for i in range(2):
                ops.node(nid(i, j, k), float(i), float(j), float(k) * (4.0 / nz))
    for j in range(2):
        for i in range(2):
            ops.fix(nid(i, j, 0), 1, 1, 1)

    if material == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    else:
        ops.nDMaterial("LadrunoJ2", 1, _E / 3.0 / (1 - 2 * nu), _E / 2.0 / (1 + nu),
                       "-iso", "voce", _SIG0, 0.0, 1.0, _HISO, "-kin", 0)
    ops.nDMaterial("LogStrain", 2, 1)
    for k in range(nz):
        conn = [nid(0, 0, k), nid(1, 0, k), nid(1, 1, k), nid(0, 1, k),
                nid(0, 0, k + 1), nid(1, 0, k + 1), nid(1, 1, k + 1), nid(0, 1, k + 1)]
        ops.element("LadrunoBrick", k + 1, *conn, 2, "-formulation", form, "-geom", "finite")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(2):
        for i in range(2):
            ops.load(nid(i, j, nz), tip_force / 4.0, 0.0, 0.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"cantilever solve failed (form={form}, nu={nu}, nz={nz}, mat={material})"
    return float(np.mean([ops.nodeDisp(nid(i, j, nz), 1) for j in range(2) for i in range(2)]))


def _cantilever_disp_controlled(form, nu, tip_disp, nx=2, ny=2, nz=6, L=4.0, nsteps=60):
    """A refined a×a×L cantilever (nx×ny elements across the section, nz up the
    length) with the transverse tip displacement PRESCRIBED — stable into the
    fully-plastic range (no collapse runaway). Returns (support reaction_x,
    max ε̄ᵖ). A volumetrically-locked element over-resists at the same imposed
    displacement ⇒ higher reaction; F-bar relieves it. The 2×2 section resolves
    bending well enough for stable plastic Newton iterations."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    a = 1.0
    nxp, nyp = nx + 1, ny + 1

    def nid(i, j, k):
        return 1 + i + nxp * j + nxp * nyp * k

    for k in range(nz + 1):
        for j in range(nyp):
            for i in range(nxp):
                ops.node(nid(i, j, k), a * i / nx, a * j / ny, L * k / nz)
    base_nodes = [nid(i, j, 0) for j in range(nyp) for i in range(nxp)]
    for n in base_nodes:
        ops.fix(n, 1, 1, 1)

    K = _E / 3.0 / (1 - 2 * nu)
    G = _E / 2.0 / (1 + nu)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", _SIG0, 0.0, 1.0, _HISO, "-kin", 0)
    ops.nDMaterial("LogStrain", 2, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, 2, "-formulation", form, "-geom", "finite")
                e += 1
    nel = e - 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tip_nodes = [nid(i, j, nz) for j in range(nyp) for i in range(nxp)]
    for n in tip_nodes:
        ops.sp(n, 1, tip_disp)                 # prescribe transverse tip displacement

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-6, 200, 0)
    # KrylovNewton: plain Newton (even with line search) diverges on this
    # bending-into-plasticity BVP; the accelerated subspace Newton is robust.
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"disp-controlled cantilever failed (form={form})"

    ops.reactions()
    R = sum(ops.nodeReaction(n, 1) for n in base_nodes)
    return abs(R), _max_eqplastic(nel)


# =========================================================================== #
#  B1 — h-convergence of a smooth finite-strain elastic cantilever            #
# =========================================================================== #
def test_B1_h_convergence_monotone_and_richardson():
    """Refine the cantilever through-thickness nz ∈ {2,4,8} (compressible
    ν = 0.3, modest tip load so the response is smooth/near-linear). The tip
    deflection must (a) increase monotonically toward the converged value as
    the bending-locking of a coarse single-layer mesh relaxes, and (b) show a
    DIMINISHING increment (Richardson: successive differences shrink by a
    roughly geometric factor), evidencing convergence rather than drift."""
    nzs = [2, 4, 8]
    d = [_cantilever("std", 0.30, nz, tip_force=0.05) for nz in nzs]
    assert all(v > 0 for v in d)
    # monotone increase (coarse mesh is too stiff in bending → under-predicts)
    assert d[0] < d[1] < d[2], f"tip deflection not monotone under refinement: {d}"
    # diminishing returns: the 4→8 change is markedly smaller than 2→4
    inc1, inc2 = d[1] - d[0], d[2] - d[1]
    assert inc2 < 0.6 * inc1, (
        f"refinement increments not shrinking (no convergence): {inc1:.3e} → {inc2:.3e}")
    # Richardson extrapolation overshoots the finest mesh by a bounded amount
    # (the limit is finite, not running away)
    rich = d[2] + inc2 * inc2 / (inc1 - inc2)        # geometric-series limit
    assert d[2] < rich < d[2] * 1.5, f"Richardson limit unreasonable: {rich:.4f} vs {d[2]:.4f}"


# =========================================================================== #
#  B2 — near-incompressible ELASTIC block: F-bar relieves volumetric locking  #
# =========================================================================== #
def test_B2_near_incompressible_elastic_locking_relieved_by_fbar():
    """At ν = 0.4999 the standard hex bending cantilever LOCKS volumetrically;
    F-bar does not. Reference-free discriminators:
      (1) std stiffens sharply going ν 0.30 → 0.4999 (deflection collapses);
      (2) F-bar deflection is ~ν-insensitive;
      (3) at ν = 0.4999, F-bar is dramatically softer than the locked std hex.
    Pushes one decade past the existing 0.499 locking test (a harder regime)."""
    nz, P = 6, 0.02
    d_std_c = _cantilever("std", 0.30, nz, P)
    d_bar_c = _cantilever("bbar", 0.30, nz, P)
    d_std_i = _cantilever("std", 0.4999, nz, P)
    d_bar_i = _cantilever("bbar", 0.4999, nz, P)
    assert min(d_std_c, d_bar_c, d_std_i, d_bar_i) > 0.0

    std_lock = d_std_i / d_std_c
    assert std_lock < 0.4, f"std hex did not lock at ν=0.4999 (ratio {std_lock:.3f})"
    bar_lock = d_bar_i / d_bar_c
    assert 0.7 < bar_lock < 1.5, f"F-bar should be ν-insensitive (ratio {bar_lock:.3f})"
    unlock = d_bar_i / d_std_i
    assert unlock > 4.0, (
        f"F-bar did not relieve near-incompressible locking: bbar/std = {unlock:.2f}")


# =========================================================================== #
#  B3 — isochoric J2 block: plastic incompressibility locking, cured by F-bar #
# =========================================================================== #
def test_B3_isochoric_j2_locking_relieved_by_fbar():
    """J2 plastic flow is isochoric (det Fᵖ = 1), so once the cantilever yields
    through the thickness the standard hex locks volumetrically just as an
    incompressible elastic solid does. DISPLACEMENT-controlled (stable past the
    plastic-collapse load, no runaway): impose the same transverse tip
    displacement for std vs bbar and compare the support REACTION. A locked
    element over-resists ⇒ higher reaction; F-bar relieves it ⇒ lower reaction.
    Both must be solidly plastic, else the test is vacuous."""
    delta = 0.3
    R_std, eqp_std = _cantilever_disp_controlled("std", 0.30, delta)
    R_bar, eqp_bar = _cantilever_disp_controlled("bbar", 0.30, delta)

    assert eqp_std > 1.0e-3 and eqp_bar > 1.0e-3, (
        f"B3 vacuous — not plastic (ε̄ᵖ std={eqp_std:.2e}, bbar={eqp_bar:.2e})")
    assert R_std > 0 and R_bar > 0
    # isochoric-plastic locking: the std hex is artificially stiff ⇒ higher
    # reaction than the F-bar element at the same imposed displacement.
    ratio = R_std / R_bar
    assert ratio > 1.15, (
        f"F-bar did not relieve isochoric-J2 locking: std/bbar reaction "
        f"= {ratio:.3f} (std should over-stiffen once plastic flow is isochoric)")


def _max_eqplastic(nel):
    """Largest equivalent plastic strain over the column's elements (post-solve)."""
    vals = []
    for e in range(1, nel + 1):
        r = ops.eleResponse(e, "material", 1, "equivalentPlasticStrain")
        if r:
            vals.append(max(abs(v) for v in r))
    return max(vals) if vals else 0.0
