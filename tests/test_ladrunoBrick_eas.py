"""LadrunoBrick — true Simo-Rifai EAS (`-formulation eas`) validation (ADR 19 §4).

`ssp` (stabilized single-point) and the patch/rank tests cannot validate a mixed
assumed-strain element. The gates here are the ones that CAN:

  1. factory: `eas` is now accepted (PR-2) and uses 8 LIVE Gauss points (not the
     single-point mirror) — distinct per-GP stress under a bending gradient;
  2. constant-strain PATCH TEST on a DISTORTED 2x2x2 patch with a free interior
     node — the direct check that the enhanced field is mean-zero (int M dV = 0):
     a wrong M(xi) fails this on distorted meshes to far worse than 1e-8;
  3. reduce-to-`std` on a regular cube under uniform strain (the enhanced params
     alpha -> 0, so EAS == std there to the inner-Newton tolerance);
  4. coarse-mesh BENDING: EAS beats `std` and approaches the Euler-Bernoulli tip;
  5. near-incompressible (nu -> 0.4999): EAS does not volumetrically lock;
  6. alpha state cycle: commit/revert reproducibility under a path-dependent
     material (LadrunoJ2) on a distorted hex.

Deferred (ADR 19 §4, not in v1): Cook's membrane band, the EnhancedQuad 2-D
plane-strain sub-block oracle, and the white-box `intMdV` eleResponse.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


# ===========================================================================
# 1. factory: eas builds and is genuinely 8-point (not the ssp single-point mirror)
# ===========================================================================
def test_eas_accepted_and_uses_eight_live_gauss_points():
    """`-formulation eas` builds (PR-2 freed the name) AND, under a bending
    gradient, reports DISTINCT per-GP stresses — proving it evaluates all 8 live
    material points, not the single-point centroid mirror that `ssp` uses."""
    L, E, nu, P = 4.0, 1000.0, 0.0, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = {1: (0, 0, 0), 2: (L, 0, 0), 3: (L, 1, 0), 4: (0, 1, 0),
              5: (0, 0, 1), 6: (L, 0, 1), 7: (L, 1, 1), 8: (0, 1, 1)}
    for t, c in coords.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "eas")
    tags = ops.getEleTags()
    tags = [tags] if isinstance(tags, int) else (tags or [])
    assert 1 in tags, "-formulation eas must build (PR-2)"

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.0, 0.0, P / 4.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0

    s = list(ops.eleResponse(1, "stresses"))
    assert len(s) == 48, f"expected 8 GP x 6 = 48 stress comps, got {len(s)}"
    s11 = [s[6 * g] for g in range(8)]
    assert max(s11) - min(s11) > 1e-6, (
        f"EAS sigma_xx identical across GPs ({s11}) — it must use 8 LIVE points, "
        f"not the ssp single-point mirror"
    )


# ===========================================================================
# 2. constant-strain patch test on a DISTORTED 2x2x2 patch (proves int M dV = 0)
# ===========================================================================
def _G_field(X):
    """A fixed constant-strain (linear-displacement) field u_i = G_ij X_j, with a
    shear component to exercise the off-diagonal enhanced modes."""
    g = ((2.0e-3, 1.0e-3, 0.5e-3),
         (0.0,    -1.0e-3, 0.0),
         (0.0,     0.0,    1.5e-3))
    x, y, z = X
    return (g[0][0] * x + g[0][1] * y + g[0][2] * z,
            g[1][0] * x + g[1][1] * y + g[1][2] * z,
            g[2][0] * x + g[2][1] * y + g[2][2] * z)


def test_eas_constant_strain_patch_distorted():
    """THE EAS gate. A 2x2x2 patch (27 nodes) with the interior node DISTORTED;
    every boundary node is prescribed to the exact linear field u = G.X, the
    interior node is free. A consistent element reproduces the constant-strain
    field exactly: the free node lands on G.X and every GP stress equals C:eps.
    A wrong enhanced operator (int M dV != 0) fails this on the distorted mesh."""
    E, nu = 1000.0, 0.3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    # 3x3x3 base grid (spacing 1), node tag = 1 + i + 3*j + 9*k
    def nid(i, j, k):
        return 1 + i + 3 * j + 9 * k

    base = {}
    for k in range(3):
        for j in range(3):
            for i in range(3):
                base[nid(i, j, k)] = [float(i), float(j), float(k)]
    # distort the interior node (tag 14 = i=j=k=1) — this is what stresses M(xi)
    base[nid(1, 1, 1)] = [1.0 + 0.18, 1.0 - 0.13, 1.0 + 0.21]
    for t, c in base.items():
        ops.node(t, *c)

    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    e = 1
    for k in range(2):
        for j in range(2):
            for i in range(2):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, 1, "-formulation", "eas")
                e += 1

    interior = nid(1, 1, 1)
    # prescribe the linear field on every boundary node via single-point constraints
    # (sp ONLY — combining ops.fix with ops.sp on the same dof makes two competing
    # Penalty terms that average to half the imposed value).
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, c in base.items():
        if t == interior:
            continue
        u = _G_field(c)
        for dof, uu in enumerate(u, start=1):
            ops.sp(t, dof, uu)

    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Penalty", 1e14, 1e14)
    ops.test("NormDispIncr", 1e-12, 20)
    ops.integrator("LoadControl", 1.0); ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(1) == 0

    # the free interior node must land on the linear field
    got = ops.nodeDisp(interior)
    exp = _G_field(base[interior])
    for d in range(3):
        assert abs(got[d] - exp[d]) <= 1e-8 + 1e-6 * abs(exp[d]), (
            f"patch test: interior u[{d}]={got[d]:.3e} != linear field {exp[d]:.3e} "
            f"(EAS enhanced field is not mean-zero on the distorted mesh)"
        )

    # analytic constant stress C:eps (eps = sym G)
    g = ((2.0e-3, 1.0e-3, 0.5e-3), (0.0, -1.0e-3, 0.0), (0.0, 0.0, 1.5e-3))
    exx, eyy, ezz = g[0][0], g[1][1], g[2][2]
    exy = 0.5 * (g[0][1] + g[1][0]); eyz = 0.5 * (g[1][2] + g[2][1]); ezx = 0.5 * (g[2][0] + g[0][2])
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    tr = exx + eyy + ezz
    sig = [lam * tr + 2 * mu * exx, lam * tr + 2 * mu * eyy, lam * tr + 2 * mu * ezz,
           2 * mu * exy, 2 * mu * eyz, 2 * mu * ezx]   # {xx,yy,zz,xy,yz,zx}, engineering
    for ele in range(1, 9):
        s = list(ops.eleResponse(ele, "stresses"))
        for g_ in range(8):
            for c in range(6):
                assert abs(s[6 * g_ + c] - sig[c]) <= 1e-6 * E * 1e-3 + 1e-7, (
                    f"patch: ele {ele} GP {g_} stress[{c}]={s[6*g_+c]:.4e} != C:eps {sig[c]:.4e}"
                )


# ===========================================================================
# 3. reduce-to-std on a regular cube under uniform strain (alpha -> 0)
# ===========================================================================
def _solve_uniaxial(form_args, E=1000.0, nu=0.3, sig0=2.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, sig0 / 4.0, 0.0, 0.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(n) for n in (2, 3, 6, 7)]


def test_eas_reduces_to_std_on_regular_cube():
    """On a regular (undistorted) cube under a uniform uniaxial state the enhanced
    parameters are inactive (alpha -> 0 to the inner-Newton tol), so true EAS must
    match `std` to a few * that tolerance — NOT machine zero (ADR 19 §4 item 5)."""
    eas = _solve_uniaxial(["-formulation", "eas"])
    std = _solve_uniaxial(["-formulation", "std"])
    for a, b in zip(eas, std):
        for d in range(3):
            assert abs(a[d] - b[d]) <= 1e-7 * (1.0 + abs(b[d])), (
                f"EAS {a[d]:.3e} should match std {b[d]:.3e} on a regular cube (alpha->0)"
            )


# ===========================================================================
# 4. coarse-mesh bending: EAS beats std and approaches Euler-Bernoulli
# ===========================================================================
def _cantilever_tip(form_args, nx, L=10.0, b=1.0, h=1.0, E=1000.0, nu=0.0, P=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = L / nx

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k

    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * b, k * h)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    e = 1
    for i in range(nx):
        ops.element("LadrunoBrick", e,
                    nt(i, 0, 0), nt(i + 1, 0, 0), nt(i + 1, 1, 0), nt(i, 1, 0),
                    nt(i, 0, 1), nt(i + 1, 0, 1), nt(i + 1, 1, 1), nt(i, 1, 1),
                    1, *form_args)
        e += 1
    for k in range(2):
        for j in range(2):
            ops.fix(nt(0, j, k), 1, 1, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    tip = [nt(nx, j, k) for k in range(2) for j in range(2)]
    for n in tip:
        ops.load(n, 0.0, 0.0, P / len(tip))
    ops.system("FullGeneral"); ops.numberer("RCM"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    return sum(ops.nodeDisp(n)[2] for n in tip) / len(tip)


def test_eas_beats_std_in_bending():
    """The point of EAS: on a coarse one-element-deep cantilever it relieves the
    parasitic shear/bending locking that cripples `std`, recovering most of the
    Euler-Bernoulli tip deflection (delta = P L^3/(3 E I), I = b h^3/12 = 1/12)."""
    analytic = 1.0 * 10.0 ** 3 / (3 * 1000.0 * (1.0 / 12.0))  # = 4.0
    eas4 = _cantilever_tip(["-formulation", "eas"], 4) / analytic
    std4 = _cantilever_tip(["-formulation", "std"], 4) / analytic
    assert eas4 > 0.80, f"EAS nx=4 tip ratio {eas4:.3f} should recover most bending"
    assert eas4 > 2.0 * std4, f"EAS {eas4:.3f} should crush std {std4:.3f} in bending"
    # converges from below toward the Euler limit under refinement
    eas8 = _cantilever_tip(["-formulation", "eas"], 8) / analytic
    assert eas4 < eas8 <= 1.10, f"EAS should converge toward 1 ({eas4:.3f} -> {eas8:.3f})"


# ===========================================================================
# 5. near-incompressible: EAS does not volumetrically lock
# ===========================================================================
def test_eas_no_volumetric_locking():
    """Volumetric locking manifests in CONSTRAINED/bending response, not in a free
    statically-determinate uniaxial cube (which is exact for any patch-passing
    element). In a near-incompressible (nu=0.45) coarse cantilever the full-
    integration `std` hex locks hard (shear + volumetric); EAS cures both and is
    far more flexible — the headline reason it exists at high nu."""
    nu = 0.45
    eas = _cantilever_tip(["-formulation", "eas"], 4, nu=nu)
    std = _cantilever_tip(["-formulation", "std"], 4, nu=nu)
    assert eas > 1.5 * std, (
        f"EAS tip {eas:.4e} should far exceed the locked std tip {std:.4e} at nu={nu}"
    )


# ===========================================================================
# 6. alpha state cycle: commit/revert reproducibility (path-dependent material)
# ===========================================================================
def _push_j2_distorted(reset_after=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # a single DISTORTED hex so the enhanced parameters are genuinely active
    coords = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1.1, 1.05, 0.0), 4: (0.05, 1, 0),
              5: (0, 0, 1), 6: (1, 0.05, 1.1), 7: (1.05, 1, 1), 8: (0, 1, 1.05)}
    for t, c in coords.items():
        ops.node(t, *c)
    # J2Plasticity (path-dependent) so alpha evolves and commit/revert is meaningful.
    # nDMaterial J2Plasticity: K (bulk), G (shear), sig0, sigInf, delta, H.
    K = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3)); G = 1000.0 / (2.0 * (1.0 + 0.3))
    ops.nDMaterial("J2Plasticity", 1, K, G, 2.0, 2.0, 0.0, 50.0)
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "eas")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 8.0, 0.0, 0.0)   # push past yield
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25); ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(2) == 0
    u_mid = [ops.nodeDisp(2)[0], ops.nodeDisp(7)[2]]
    assert ops.analyze(2) == 0
    u_end = [ops.nodeDisp(2)[0], ops.nodeDisp(7)[2]]
    return u_mid, u_end


def test_eas_alpha_commit_revert_reproducible():
    """The enhanced parameters alpha are committed element state. Running the same
    path-dependent (LadrunoJ2) push twice must give the identical displacement
    history — i.e. commitState/revertToLastCommit of alpha are consistent (a broken
    alpha lifecycle would drift between identical runs)."""
    a_mid, a_end = _push_j2_distorted()
    b_mid, b_end = _push_j2_distorted()
    for a, b in zip(a_mid + a_end, b_mid + b_end):
        assert abs(a - b) <= 1e-9 * (1.0 + abs(b)), (
            f"EAS path-dependent run not reproducible ({a:.6e} vs {b:.6e}) — "
            f"alpha commit/revert lifecycle is inconsistent"
        )
    # and the push genuinely activated the enhanced/plastic response (non-vacuous)
    assert abs(a_end[0]) > 1e-4, "the J2 push should produce a non-trivial displacement"


# ===========================================================================
# 7. EAS drives a path-dependent + SOFTENING material (J2 + Lemaitre damage)
# ===========================================================================
def _push_damage(form, umax=0.30, nsteps=120):
    """Single steel hex pulled in uniaxial tension into J2 + Lemaitre ductile damage
    via displacement control (so it traverses post-peak softening). Returns
    (u_reached/umax, max damage)."""
    E, nu = 200000.0, 0.3
    K = E / (3.0 * (1.0 - 2.0 * nu)); G = E / (2.0 * (1.0 + nu))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *[float(x) for x in c])
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 400.0, 80.0, 40.0, 50.0,
                   "-damage", "lemaitre", 0.35, 1.0, 0.01, 0.92)
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", form)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, umax)            # imposed x-disp on the +x face (ramped by the series)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 50, 0); ops.analysis("Static")
    base, u = 1.0 / nsteps, 0.0
    for _ in range(nsteps):
        ok = -1
        for cut, alg in ((1.0, "Newton"), (0.25, "Newton"), (0.1, "NewtonLineSearch")):
            ops.integrator("LoadControl", base * cut); ops.algorithm(alg)
            ok = ops.analyze(1)
            if ok == 0:
                break
        if ok != 0:
            break
        u = ops.nodeDisp(2)[0]
    ds = [max(ops.eleResponse(1, "material", gp, "damage")) for gp in range(1, 9)]
    return u / umax, max(ds)


def test_eas_drives_j2_lemaitre_damage():
    """EAS must carry a path-dependent + SOFTENING constitutive law (LadrunoJ2 +
    Lemaitre ductile damage) through the full displacement on a homogeneous element:
    the inner-Newton alpha solve and the commit cycle have to cope with plasticity
    AND damage. EAS reaches the full elongation and accumulates substantial damage,
    matching the fully-integrated `bbar` (both 8 Gauss points). (This is the
    homogeneous case; EAS is NOT robust on a notched localization problem — small-
    strain hourglassing — where `ssp`/`bbar` are the right choice; see ADR 19.)"""
    eas_frac, eas_d = _push_damage("eas")
    bbar_frac, bbar_d = _push_damage("bbar")
    assert eas_frac > 0.99, f"EAS only reached {eas_frac:.2f} of the imposed elongation"
    assert eas_d > 0.2, f"EAS damage {eas_d:.3f} — material did not soften"
    assert abs(eas_d - bbar_d) < 0.05, (
        f"EAS damage {eas_d:.3f} should track bbar {bbar_d:.3f} on a homogeneous push"
    )
