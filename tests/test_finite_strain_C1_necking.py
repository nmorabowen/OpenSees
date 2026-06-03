"""Finite-strain trifecta — **C1 Simo–Armero necking bar** (L3 literature).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §5
"L3 — Literature benchmarks" row C1, and §9 phase **P1** (the flagship
finite-strain-J2 validation). A round bar pulled axially past the Considère
point develops a NECK: a localized reduction of the cross-section at a seeded
imperfection, driven by geometric softening even though the material hardens.
This is THE canonical finite-strain elastoplastic benchmark (Simo 1992; de Souza
Neto, Perić & Owen §14.10).

Model — `element LadrunoBrick -geom finite -formulation bbar` (F-bar; necking is
isochoric-plastic, so the standard hex would lock — see the L2 B3 test) driving
`nDMaterial LogStrain` over `nDMaterial LadrunoJ2` with **isotropic** saturation
hardening (Simo's law σ_y = σ0 + (σ∞−σ0)(1−e^{−δp}) + Hp; objective through the
log-strain lift at any deformation — combined hardening would hit the §14.11
boundary). Half-bar symmetry (z=0 mid-plane), axisymmetric axis line pinned, end
pulled by prescribed displacement. The round cross-section is a structured
"squircle" hex lattice (FG-squircle map of a square grid onto the disk — no
gmsh, fully deterministic, all-hex, no axis degeneracy), so the whole study is
**Zone-A** and travels with the PR (plan §6).

> [!scope] What this test validates vs defers
> It establishes the necking PHYSICS to engineering rigor:
>   * the neck LOCALIZES at the imperfection plane (r_neck < r_end, min on-axis);
>   * plastic strain localizes there (ε̄ᵖ_neck ≫ ε̄ᵖ_end);
>   * the reaction–elongation curve SOFTENS (incremental stiffness collapses
>     toward the Considère peak) — geometric softening under a hardening material;
>   * a PERFECT bar (no imperfection) deforms uniformly — the neck is physical,
>     not a mesh artifact.
> The tight QUANTITATIVE neck-ratio match to Simo (r/r0 ≈ 0.5 at 7 mm, tol 3 %)
> needs a mesh-converged FINE model and is the **Zone-B C1 contract** (plan P6 /
> §6); a coarse structured lattice necks only mildly before element distortion
> stalls it. See Ladruno_implementation/18_finite_strain_validation_report.md §4.

Solver recipe (from B3 / the report): KrylovNewton + EnergyIncr + an unsymmetric
sparse solver (F-bar tangent is unsymmetric); plain Newton diverges.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t2a]

# Simo (1992) / dSNPO §14.10 round bar, consistent N–mm–MPa units.
_R0 = 6.413          # mm, nominal radius
_LHALF = 26.667      # mm, half length (z=0 is the symmetry mid-plane)
_E, _NU = 206900.0, 0.29
_SIG0, _SINF, _DELTA, _HLIN = 450.0, 715.0, 16.93, 129.24   # MPa saturation law


def _squircle(u, v, R):
    """FG-squircle map of the square [-1,1]² onto the radius-R disk (orientation
    preserving, smooth, no center degeneracy)."""
    return (R * u * math.sqrt(max(0.0, 1.0 - 0.5 * v * v)),
            R * v * math.sqrt(max(0.0, 1.0 - 0.5 * u * u)))


def _build(n, nz, imp):
    """Structured squircle-hex half-bar. imp = fractional radius dip at z=0
    (localized Gaussian imperfection); imp=0 ⇒ perfect prismatic bar.
    Returns (idx, nel, end_nodes, z0_nodes, boundary_ij, neck_elems, end_elems)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    c = n // 2
    zs = [_LHALF * (k / nz) ** 2.0 for k in range(nz + 1)]      # bias toward neck
    Rz = lambda z: _R0 * (1.0 - imp * math.exp(-(z / (0.22 * _LHALF)) ** 2))
    idx = {}
    tag = 1
    for k, z in enumerate(zs):
        R = Rz(z)
        for j in range(n + 1):
            for i in range(n + 1):
                x, y = _squircle(-1 + 2 * i / n, -1 + 2 * j / n, R)
                ops.node(tag, x, y, z)
                idx[(i, j, k)] = tag
                tag += 1
    K = _E / 3.0 / (1.0 - 2.0 * _NU)
    G = _E / 2.0 / (1.0 + _NU)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce",
                   _SIG0, _SINF - _SIG0, _DELTA, _HLIN, "-kin", 0)
    ops.nDMaterial("LogStrain", 2, 1)
    e = 1
    neck_elems, end_elems = [], []
    for k in range(nz):
        for j in range(n):
            for i in range(n):
                conn = [idx[(i, j, k)], idx[(i + 1, j, k)], idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                        idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)], idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                ops.element("LadrunoBrick", e, *conn, 2, "-formulation", "bbar", "-geom", "finite")
                if k == 0:
                    neck_elems.append(e)
                if k == nz - 1:
                    end_elems.append(e)
                e += 1
    nel = e - 1
    for j in range(n + 1):
        for i in range(n + 1):
            ops.fix(idx[(i, j, 0)], 0, 0, 1)        # z=0 symmetry plane
    for k in range(nz + 1):
        ops.fix(idx[(c, c, k)], 1, 1, 0)            # axisymmetric axis line
    end_nodes = [idx[(i, j, nz)] for j in range(n + 1) for i in range(n + 1)]
    z0_nodes = [idx[(i, j, 0)] for j in range(n + 1) for i in range(n + 1)]
    bnd = [(i, j) for j in range(n + 1) for i in range(n + 1) if i in (0, n) or j in (0, n)]
    return idx, nel, end_nodes, z0_nodes, bnd, neck_elems, end_elems


def _mean_radius(idx, bnd, k):
    rs = []
    for (i, j) in bnd:
        t = idx[(i, j, k)]
        rs.append(math.hypot(ops.nodeCoord(t, 1) + ops.nodeDisp(t, 1),
                              ops.nodeCoord(t, 2) + ops.nodeDisp(t, 2)))
    return sum(rs) / len(rs)


def _max_eqp(elems):
    m = 0.0
    for e in elems:
        r = ops.eleResponse(e, "material", 1, "equivalentPlasticStrain")
        if r:
            m = max(m, max(abs(v) for v in r))
    return m


def _run(imp, end_disp, n=4, nz=8, nsteps=40):
    """Pull the half-bar end by end_disp, returning per-step
    (elongation, reaction, r_neck/R0, r_end/R0) plus the final state handle."""
    idx, nel, end_nodes, z0_nodes, bnd, neck_elems, end_elems = _build(n, nz, imp)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nN in end_nodes:
        ops.sp(nN, 3, end_disp)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("UmfPack")              # unsymmetric sparse (F-bar tangent)
    except Exception:
        ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-6, 150, 0)
    ops.algorithm("KrylovNewton")
    step = 1.0 / nsteps
    ops.integrator("LoadControl", step)
    ops.analysis("Static")
    curve = []
    done, cuts = 0.0, 0
    while done < 1.0 - 1e-9 and cuts < 8:
        if ops.analyze(1) == 0:
            done += step
            ops.reactions()
            R = abs(sum(ops.nodeReaction(nN, 3) for nN in z0_nodes))
            curve.append((done * end_disp, R,
                          _mean_radius(idx, bnd, 0) / _R0,
                          _mean_radius(idx, bnd, nz) / _R0))
        else:
            step *= 0.5
            cuts += 1
            ops.integrator("LoadControl", step)
    return curve, (idx, bnd, nz, neck_elems, end_elems)


# Shared run: a 2.8 mm half-elongation (5.6 mm total) neck — reliably reached
# before coarse-mesh element distortion would stall it, yet solidly necking
# (r_neck/R0 ≈ 0.91, well past the Considère stiffness collapse).
@pytest.fixture(scope="module")
def necking_run():
    curve, state = _run(imp=0.02, end_disp=2.8, n=4, nz=8, nsteps=44)
    assert curve and curve[-1][0] >= 2.5, (
        f"necking run did not reach a meaningful elongation: {curve[-1][0] if curve else 0:.2f} mm")
    return curve, state


# =========================================================================== #
#  C1.1 — the neck LOCALIZES at the imperfection plane                         #
# =========================================================================== #
def test_C1_neck_localizes_at_imperfection_plane(necking_run):
    curve, (idx, bnd, nz, _, _) = necking_run
    r_neck = curve[-1][2]
    r_end = curve[-1][3]
    # (1) the section shrank at the neck
    assert r_neck < 0.97, f"no necking developed (r_neck/R0={r_neck:.3f})"
    # (2) and it localized: the neck radius is clearly below the end radius
    assert r_neck < r_end - 0.01, (
        f"neck did not localize: r_neck/R0={r_neck:.4f} not < r_end/R0={r_end:.4f}")
    # (3) the minimum cross-section is AT the imperfection plane (z=0), scanning z
    profile = [_mean_radius(idx, bnd, k) for k in range(nz + 1)]
    assert profile[0] == min(profile), (
        f"minimum radius not at the symmetry/imperfection plane: {profile}")
    # (4) radius reduction is monotone in elongation (the neck only deepens)
    rn = [row[2] for row in curve]
    assert all(rn[i] >= rn[i + 1] - 1e-6 for i in range(len(rn) - 1)), (
        f"neck radius not monotonically decreasing: {rn}")


# =========================================================================== #
#  C1.2 — geometric softening: the reaction–elongation curve flattens          #
#         (Considère necking onset) despite a hardening material               #
# =========================================================================== #
def test_C1_reaction_softens_toward_considere_peak(necking_run):
    curve, _ = necking_run
    el = np.array([row[0] for row in curve])
    R = np.array([row[1] for row in curve])
    assert R[-1] > 0
    # secant incremental stiffness over the first vs the last fifth of the path
    nseg = max(2, len(curve) // 5)
    k_init = (R[nseg] - R[0]) / (el[nseg] - el[0])
    k_late = (R[-1] - R[-1 - nseg]) / (el[-1] - el[-1 - nseg])
    assert k_init > 0, "reaction did not rise initially (no elastic/hardening branch)"
    # geometric softening: late stiffness collapses to a small fraction of the
    # initial — the hallmark of approaching the Considère (necking) instability.
    # Measured k_late/k_init ≈ 0.07 at 2.8 mm; 0.12 leaves comfortable margin.
    assert k_late < 0.12 * k_init, (
        f"no geometric softening: late dR/dl {k_late:.1f} not ≪ initial {k_init:.1f}")


# =========================================================================== #
#  C1.3 — plastic strain localizes at the neck                                 #
# =========================================================================== #
def test_C1_plastic_strain_localizes_at_neck(necking_run):
    _, (idx, bnd, nz, neck_elems, end_elems) = necking_run
    eqp_neck = _max_eqp(neck_elems)
    eqp_end = _max_eqp(end_elems)
    assert eqp_neck > 1.0e-2, f"neck did not yield substantially (ε̄ᵖ={eqp_neck:.3e})"
    # neck plastic strain clearly exceeds the end (measured ratio ≈ 1.9 at
    # 2.8 mm; the whole bar yields, so the localization is real but modest at
    # this coarse-mesh elongation — it sharpens in the fine Zone-B run).
    assert eqp_neck > 1.5 * max(eqp_end, 1.0e-12), (
        f"plastic strain did not localize: neck ε̄ᵖ {eqp_neck:.3e} vs end {eqp_end:.3e}")


# =========================================================================== #
#  C1.4 — control: a PERFECT bar (no imperfection) deforms UNIFORMLY            #
#         (proves the neck is imperfection-seeded, not a mesh artifact)         #
# =========================================================================== #
def test_C1_perfect_bar_deforms_uniformly():
    curve, (idx, bnd, nz, _, _) = _run(imp=0.0, end_disp=2.0, n=4, nz=8, nsteps=30)
    assert curve and curve[-1][0] >= 1.5
    profile = [_mean_radius(idx, bnd, k) for k in range(nz + 1)]
    spread = (max(profile) - min(profile)) / np.mean(profile)
    # a perfect prismatic bar contracts uniformly (Poisson + plastic): the radius
    # is nearly z-independent — no localized neck. (A few % residual spread comes
    # from the free end vs the symmetry plane, not from a neck.)
    assert spread < 0.03, (
        f"perfect bar developed a spurious neck (radius spread {spread:.3f} along z): {profile}")
    assert curve[-1][2] == pytest.approx(curve[-1][3], abs=0.02), (
        "perfect bar: neck-plane and end-plane radii should match (uniform contraction)")
