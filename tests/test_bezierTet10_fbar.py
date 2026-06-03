"""BezierTet10 ``-geom finite -bbar`` — F-bar (bbar + finite) battery (3D tet).

F-bar (dSNPO §15.1) is the large-strain volumetric-locking cure: every GP is
driven by F̄ = (J₀/J)^(1/3) F, sharing the element-centroid dilatation J₀ (tet
centroid, barycentric L=(¼,¼,¼)). The residual is UNCHANGED (eq 15.9 — only σ̄
changes); the tangent gains the eq 15.10 coupling

    K_{(a,i)(b,k)} += ∫ (Σ_j g_aj q_ij)(G₀_kb − g_kb) dv ,
    q_ij = (1/3) a_ijpp − (2/3) σ̄_ij     (eq 15.11, GENERALLY UNSYMMETRIC)

using the SAME a4 = c̄ − σ̄δ as the std term. The −(2/3)σ̄ is the spatial
initial-stress part (NOT a (1/3)c shortcut) — the COEFFICIENT TRAP.

Gates:
  * FD-consistent-tangent headline (UNSYMMETRIC-aware: assert the FD match AND
    that the tangent is genuinely unsymmetric — proves the eq 15.10 coupling is
    live, G₀ ≠ g);
  * FD-tangent at ν=0.499 — the regression guard for the −(2/3)σ̄ coefficient
    (it blows up as the bulk modulus diverges);
  * residual reduces to std finite under HOMOGENEOUS F (F̄=F), but the tangent
    DIFFERS (centroid G₀ ≠ per-GP g — F-bar is a distinct volumetric treatment);
  * F-bar relieves the near-incompressible PRESSURE spike (tet displacement
    locking is mild; the payoff is the stress/pressure field — Kadapa §6);
  * parse: -geom finite -bbar accepted; -pressure + finite still rejected.

All F-bar tests run under system('FullGeneral') (the tangent is unsymmetric).
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

_E = 200.0
_NU = 0.3

# mid-edge (vertexA, vertexB) order for nodes 5..10 (TenNodeTetrahedron / element)
_EDGE_V = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

# A distorted (positive-volume) reference tet, so the centroid F₀/G₀ and the
# per-GP g genuinely differ (G₀ ≠ g ⇒ the F-bar coupling is exercised).
_VERTS = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.12, 0.05),
    3: (0.18, 1.00, 0.08),
    4: (0.07, 0.15, 1.00),
}

_WIGGLE = [
    (0.000, 0.000, 0.000), (0.011, -0.008, 0.006), (-0.007, 0.010, -0.009),
    (0.008, 0.005, 0.011), (-0.010, 0.007, -0.006), (0.006, -0.011, 0.008),
    (-0.009, 0.006, 0.010), (0.012, -0.005, -0.007), (-0.006, 0.009, 0.011),
    (0.010, -0.007, -0.005),
]


def _tet_nodes():
    coord = dict(_VERTS)
    for e, (a, b) in enumerate(_EDGE_V):
        va, vb = _VERTS[a + 1], _VERTS[b + 1]
        coord[5 + e] = tuple(0.5 * (va[i] + vb[i]) for i in range(3))
    return coord


_COORD = _tet_nodes()


# the two F-bar variants (the -fbar switch); centroid is the default.
_FBAR_MODES = ["centroid", "mean_dilatation"]


def _build(bbar, nu=_NU, fbar="centroid"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    ops.nDMaterial("LogStrain", 2, 1)
    args = ["BezierTet10", 1, *range(1, 11), 2]
    if bbar:
        args.append("-bbar")
        args += ["-fbar", fbar]
    args += ["-geom", "finite"]
    ops.element(*args)


def _impose_and_solve(u, bbar, nu=_NU, fbar="centroid"):
    _build(bbar, nu, fbar)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _COORD:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")               # F-bar tangent is unsymmetric
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _resisting_force(u, bbar, nu=_NU, fbar="centroid"):
    assert _impose_and_solve(u, bbar, nu, fbar) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, bbar, nu=_NU, fbar="centroid"):
    assert _impose_and_solve(u, bbar, nu, fbar) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 30 * 30, "BezierTet10 must answer eleResponse('stiff') with 30x30"
    return np.array(K, dtype=float).reshape(30, 30)


def _affine_disp(Fbar, wiggle=False):
    u = np.zeros(30)
    I = np.eye(3)
    F = np.asarray(Fbar)
    for tag, (x, y, z) in _COORD.items():
        X = np.array([x, y, z])
        d = (F - I) @ X
        if wiggle:
            d = d + np.array(_WIGGLE[tag - 1])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


def _rot(thz, thy):
    cz, sz = math.cos(thz), math.sin(thz)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(thy), math.sin(thy)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return Ry @ Rz


# 4-point tet rule barycentric (L1,L2,L3) — must match BezierTet10::GP4_L exactly
# (so GP index 0 of eleResponse('stresses') is the F̄ we reconstruct here).
_A = 0.585410196624968
_B = 0.138196601125011
_GP4_L = [(_A, _B, _B), (_B, _A, _B), (_B, _B, _A), (_B, _B, _B)]


def _shape_derivs_bary(L1, L2, L3):
    """∂Nₐ/∂L_i for the 10-node quadratic Bernstein tet — byte-mirror of
    BezierTet10::shapeDerivatives (node order N1=L1²…N10=2L2L4, L4=1−ΣL)."""
    L4 = 1.0 - L1 - L2 - L3
    dN = np.zeros((3, 10))
    dN[0] = [2*L1, 0.0, 0.0, -2*L4, 2*L2, 0.0, 2*L3, 2*(L4 - L1), -2*L3, -2*L2]
    dN[1] = [0.0, 2*L2, 0.0, -2*L4, 2*L1, 2*L3, 0.0, -2*L1, -2*L3, 2*(L4 - L2)]
    dN[2] = [0.0, 0.0, 2*L3, -2*L4, 0.0, 2*L2, 2*L1, -2*L1, 2*(L4 - L3), -2*L2]
    return dN


def _F_at(u, L):
    """Deformation gradient at barycentric point L, computed INDEPENDENTLY of the
    element (reference control points = _COORD; u are the Bernstein control
    coefficients). F_iJ = δ_iJ + Σ_a u_a[i] ∂N_a/∂X_J."""
    X = np.array([_COORD[a] for a in range(1, 11)])    # 10×3 reference coords
    dN = _shape_derivs_bary(*L)                          # 3×10 (∂N/∂L)
    Jm = dN @ X                                          # 3×3 J[i][j] = Σ_a ∂N_i/∂L · X_aj
    dN_dX = np.linalg.inv(Jm) @ dN                       # 3×10 ∂N/∂X
    F = np.eye(3)
    U = u.reshape(10, 3)
    for a in range(10):
        F += np.outer(U[a], dN_dX[:, a])                 # F_iJ += u_a[i] ∂N_a/∂X_J
    return F


# --------------------------------------------------------------------------- #
#  Headline: consistent tangent == FD of the resisting force, UNSYMMETRIC.     #
# --------------------------------------------------------------------------- #
def _fd_tangent_gate(nu, fbar):
    # Non-affine (wiggled) state so the bar gradient Ĝ differs from the per-GP g
    # and the eq 15.10 coupling is genuinely exercised (for BOTH variants).
    Fbar = [[1.15, 0.08, -0.05],
            [0.04, 0.92, 0.06],
            [-0.03, 0.05, 1.10]]
    u = _affine_disp(Fbar, wiggle=True)

    K = _element_tangent(u, True, nu, fbar)
    f0 = _resisting_force(u, True, nu, fbar)
    assert np.isfinite(K).all() and np.isfinite(f0).all()

    h = 1.0e-6
    Kfd = np.zeros((30, 30))
    for d in range(30):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, True, nu, fbar)
                     - _resisting_force(um, True, nu, fbar)) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    return K, scale, err


@pytest.mark.parametrize("fbar", _FBAR_MODES)
def test_fbar_consistent_tangent_matches_finite_difference(fbar):
    # THE arbiter — for centroid AND mean_dilatation. The mean-dilatation tangent
    # (Ĵ→J̄, Ĝ→Ḡ volume averages) has no prior validated baseline, so this FD match
    # is the proof that the (J₀,G₀)→(J̄,Ḡ) derivation is the complete consistent
    # tangent. The FD match is the SOLE arbiter of the eq 15.11 q-coupling.
    K, scale, err = _fd_tangent_gate(_NU, fbar)
    assert err <= 1.0e-4 * scale, (
        f"F-bar[{fbar}] tangent != FD of resisting force: max abs err {err:.3e} vs "
        f"scale {scale:.3e} (eq 15.11 coupling q_ij = (1/3)a_ijpp − (2/3)σ_ij?)"
    )
    # F-bar carries a genuinely unsymmetric coupling on this non-affine state:
    # confirm the test is actually exercising it (else it proves nothing).
    asym = np.abs(K - K.T).max()
    assert asym > 1.0e-6 * scale, (
        f"F-bar[{fbar}] tangent came out symmetric on a non-affine state — the "
        "eq 15.10 coupling term is not being exercised (Ĝ == g?)"
    )


@pytest.mark.parametrize("fbar", _FBAR_MODES)
def test_fbar_consistent_tangent_matches_fd_near_incompressible(fbar):
    # The −(2/3)σ̄ term feeds the volumetric block whose magnitude blows up as
    # ν→1/2 (bulk modulus ~ E/(3(1−2ν))). Re-running the FD gate at ν=0.499 guards
    # that term where it matters most — for both variants.
    _, scale, err = _fd_tangent_gate(0.499, fbar)
    assert err <= 1.0e-4 * scale, (
        f"near-incompressible F-bar[{fbar}] tangent != FD: max abs err {err:.3e} vs "
        f"scale {scale:.3e} (the −(2/3)σ̄ term in eq 15.11 is most consequential at ν→1/2)"
    )


# --------------------------------------------------------------------------- #
#  Residual reduces to std finite under HOMOGENEOUS F (F̄ = F when J₀ = J),     #
#  but the TANGENT does NOT (centroid G₀ ≠ per-GP g).                           #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fbar", _FBAR_MODES)
def test_fbar_residual_matches_std_on_homogeneous_deformation(fbar):
    Fm = [[1.10, 0.05, -0.02],
          [0.03, 0.96, 0.04],
          [-0.01, 0.02, 1.07]]
    u = _affine_disp(Fm, wiggle=False)            # affine => homogeneous F (Ĵ = J)

    f_bar = _resisting_force(u, True, fbar=fbar)
    f_std = _resisting_force(u, False)
    scale = max(np.abs(f_std).max(), 1.0e-30)
    assert np.abs(f_bar - f_std).max() <= 1.0e-9 * scale, (
        f"F-bar[{fbar}] residual did not reduce to std finite on a homogeneous "
        f"deformation (F̄ should equal F): max diff {np.abs(f_bar - f_std).max():.3e}"
    )

    # The F-bar tangent legitimately DIFFERS even here (bar gradient Ĝ ≠ per-GP g):
    # confirm the coupling is present (not a silent no-op).
    K_bar = _element_tangent(u, True, fbar=fbar)
    K_std = _element_tangent(u, False)
    ks = np.abs(K_std).max()
    assert np.abs(K_bar - K_std).max() > 1.0e-6 * ks, (
        f"F-bar[{fbar}] tangent is identical to std on a distorted tet — the "
        "eq 15.10 coupling term appears to be missing"
    )


# --------------------------------------------------------------------------- #
#  Parse: -geom finite -bbar accepted; -pressure + finite still rejected.       #
# --------------------------------------------------------------------------- #
def _ele_absent_after(make):
    try:
        make()
    except Exception:
        pass
    tags = ops.getEleTags()
    if tags is None:
        return
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags


def test_fbar_parse_accept_and_pressure_still_rejected():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.nDMaterial("LogStrain", 2, 1)
    ops.element("BezierTet10", 1, *range(1, 11), 2, "-bbar", "-geom", "finite")
    tags = ops.getEleTags()
    tags = [tags] if isinstance(tags, int) else (tags or [])
    assert 1 in tags, "-geom finite -bbar (F-bar) should be accepted"

    # -fbar centroid|mean_dilatation both accepted with -bbar -geom finite
    for fb in _FBAR_MODES:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for tag, (x, y, z) in _COORD.items():
            ops.node(tag, x, y, z)
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
        ops.nDMaterial("LogStrain", 2, 1)
        ops.element("BezierTet10", 1, *range(1, 11), 2, "-bbar", "-fbar", fb, "-geom", "finite")
        tags = ops.getEleTags()
        tags = [tags] if isinstance(tags, int) else (tags or [])
        assert 1 in tags, f"-fbar {fb} should be accepted"

    # an unknown -fbar value is rejected
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.nDMaterial("LogStrain", 2, 1)
    _ele_absent_after(
        lambda: ops.element("BezierTet10", 1, *range(1, 11), 2,
                            "-bbar", "-fbar", "bogus", "-geom", "finite")
    )

    # -pressure + finite (with or without -bbar) is still rejected
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.nDMaterial("LogStrain", 2, 1)
    _ele_absent_after(
        lambda: ops.element("BezierTet10", 1, *range(1, 11), 2,
                            "-bbar", "-pressure", 1.0, "-geom", "finite")
    )


# --------------------------------------------------------------------------- #
#  F-bar relieves near-incompressible PRESSURE locking — the mechanism, made    #
#  directly observable. F-bar pins det F̄ = J₀ at EVERY GP, so for a near-       #
#  incompressible material (pressure p = −tr(σ)/3 is a function of J only) the   #
#  GP pressure is UNIFORM across the 4 Gauss points under F-bar, while std lets  #
#  J vary GP-to-GP and so its pressure SPREADS — the spurious volumetric         #
#  oscillation that is the very thing F-bar cures. Imposed (deterministic) non-  #
#  homogeneous field at ν=0.499 — no solver, so std cannot just diverge.         #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fbar", _FBAR_MODES)
def test_fbar_pins_uniform_pressure_under_inhomogeneous_near_incompressible(fbar):
    nu = 0.499
    # affine stretch + per-node wiggle ⇒ a NON-homogeneous F (det J varies across
    # the 4 GPs). Moderate so det F stays well positive at every GP for std too.
    Fbar = [[1.06, 0.03, -0.02],
            [0.02, 0.97, 0.03],
            [-0.01, 0.02, 1.04]]
    u = _affine_disp(Fbar, wiggle=True)

    def gp_pressures(bbar, fb="centroid"):
        assert _impose_and_solve(u, bbar, nu, fb) == 0, "imposed near-incompressible solve failed"
        s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(4, 6)
        return -(s[:, 0] + s[:, 1] + s[:, 2]) / 3.0       # hydrostatic pressure per GP

    p_std = gp_pressures(False)
    p_bar = gp_pressures(True, fbar)                       # both variants pin det F̄ uniform
    spread_std = float(p_std.max() - p_std.min())
    spread_bar = float(p_bar.max() - p_bar.min())

    # sanity: std genuinely spreads its pressure across GPs (the locking artefact)
    assert spread_std > 1.0e-3 * abs(float(p_std.mean())), (
        f"std GP pressures unexpectedly uniform (spread {spread_std:.3e}); the test "
        "is not exercising the inhomogeneity"
    )
    # F-bar pins det F̄ = J₀ at every GP ⇒ pressure (a function of J only for the
    # Hencky-elastic material) is uniform across GPs to MACHINE precision — a far
    # stronger, spread_std-independent bound than the 5%-of-std ratio.
    assert spread_bar < 1.0e-9 * max(abs(float(p_bar.mean())), 1.0), (
        f"F-bar pressure not uniform across GPs to machine precision: "
        f"spread_bar={spread_bar:.3e}, mean p={p_bar.mean():.3e} — F-bar should "
        "pin det F̄=J₀ ⇒ a single dilatation ⇒ a single pressure at every GP"
    )
    # and (redundantly) far tighter than the std spread it cures.
    assert spread_bar < 0.05 * spread_std


# --------------------------------------------------------------------------- #
#  Load-driven, multi-step Newton through the UNSYMMETRIC F-bar tangent. The     #
#  point is to drive K (the eq 15.10 unsymmetric coupling) through solver        #
#  convergence with FullGeneral, not just match it by FD at a frozen state — so  #
#  this runs in a WELL-CONDITIONED regime (ν=0.3, the proven finite-std setup).  #
#  Near-incompressible load-driven solves are ill-conditioned by nature (the     #
#  Newton predictor overshoots regardless of F-bar); that regime's BEHAVIOUR is  #
#  covered deterministically by the uniform-pressure test above.                 #
# --------------------------------------------------------------------------- #
def _kuhn_tets(n0, n1, n2, n3, n4, n5, n6, n7):
    return [
        (n0, n1, n2, n6), (n0, n2, n3, n6),
        (n0, n3, n7, n6), (n0, n7, n4, n6),
        (n0, n4, n5, n6), (n0, n5, n1, n6),
    ]


def test_fbar_load_driven_converges():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nz, side, nu = 4, 0.5, _NU              # ν=0.3 — well-conditioned (proven setup)

    def vid(i, j, k):
        return 1 + i + 2 * j + 4 * k

    for k in range(nz + 1):
        for j in range(2):
            for i in range(2):
                ops.node(vid(i, j, k), side * i, side * j, float(k))
    for j in range(2):
        for i in range(2):
            ops.fix(vid(i, j, 0), 1, 1, 1)

    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    ops.nDMaterial("LogStrain", 2, 1)

    etag = 1
    nc = vid(1, 1, nz) + 1
    for k in range(nz):
        hexc = [vid(0, 0, k), vid(1, 0, k), vid(1, 1, k), vid(0, 1, k),
                vid(0, 0, k + 1), vid(1, 0, k + 1), vid(1, 1, k + 1), vid(0, 1, k + 1)]
        for tet in _kuhn_tets(*hexc):
            v = list(tet)
            mids = []
            for (a, b) in _EDGE_V:
                xm = 0.5 * (np.array(ops.nodeCoord(v[a])) + np.array(ops.nodeCoord(v[b])))
                ops.node(nc, *[float(c) for c in xm]); mids.append(nc); nc += 1
            ops.element("BezierTet10", etag, *v, *mids, 2, "-bbar", "-geom", "finite")
            etag += 1

    top = [vid(i, j, nz) for j in range(2) for i in range(2)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.01 / 4.0, 0.0, 0.0)        # tip shear (the proven finite-std magnitude)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")                    # F-bar tangent is unsymmetric
    ops.test("NormDispIncr", 1.0e-9, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25)
    ops.analysis("Static")

    tip_prev = 0.0
    for step in range(4):
        assert ops.analyze(1) == 0, f"F-bar cantilever step {step} failed to converge"
        tip = float(np.mean([ops.nodeDisp(t, 1) for t in top]))
        assert np.isfinite(tip) and tip > tip_prev, (
            f"F-bar tip deflection should grow monotonically (step {step}: "
            f"{tip:.4e} <= prev {tip_prev:.4e})"
        )
        tip_prev = tip
    assert tip_prev > 0.0


# --------------------------------------------------------------------------- #
#  Oracle leg: pin the bar dilatation Ĵ (and the (Ĵ/J)^(1/3) scaling) against    #
#  an INDEPENDENT closed-form, for BOTH variants. The FD-of-self gates and the   #
#  spread-only / homogeneous-F gates structurally cannot catch a consistent      #
#  wrong-Ĵ (it shifts σ̄ and the tangent together, and keeps det F̄ uniform). We   #
#  reconstruct F̄ at GP 0 from scratch — Ĵ = det F at the tet centroid (centroid  #
#  variant) or the GP-mean of det F (mean_dilatation; the simple mean is exact   #
#  for a straight-sided tet where detJ_ref is constant) — and check the GP-0      #
#  Cauchy against the LogStrain oracle. A wrong centroid/J̄ value fails this.      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fbar", _FBAR_MODES)
def test_fbar_dilatation_matches_oracle(fbar):
    import logstrain_reference as o
    nu = _NU
    Fimp = [[1.08, 0.05, -0.03],
            [0.02, 0.95, 0.04],
            [-0.02, 0.03, 1.06]]
    u = _affine_disp(Fimp, wiggle=True)          # non-affine ⇒ Ĵ ≠ J(GP)
    assert _impose_and_solve(u, True, nu, fbar) == 0, "imposed F-bar solve failed"

    J_gp = np.array([float(np.linalg.det(_F_at(u, L))) for L in _GP4_L])
    if fbar == "mean_dilatation":
        Jhat = float(J_gp.mean())                # J̄ = mean det F (straight-sided tet)
    else:
        Jhat = float(np.linalg.det(_F_at(u, (0.25, 0.25, 0.25))))   # centroid J₀
    F_gp0 = _F_at(u, _GP4_L[0])
    Jg = J_gp[0]
    assert Jhat > 0.0 and Jg > 0.0
    Fbar0 = (Jhat / Jg) ** (1.0 / 3.0) * F_gp0   # F̄ at GP 0 (det = Ĵ)

    Kbulk = _E / (3.0 * (1.0 - 2.0 * nu))
    G = _E / (2.0 * (1.0 + nu))
    sig, tau, J, eps = o.cauchy_from_F_elastic(Fbar0, Kbulk, G)
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2],
                      sig[0, 1], sig[1, 2], sig[2, 0]])

    s0 = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(4, 6)[0]
    tol = 1.0e-6 * max(np.abs(sig_v).max(), 1.0)
    assert np.allclose(s0, sig_v, rtol=1.0e-6, atol=tol), (
        f"F-bar[{fbar}] GP-0 Cauchy {s0} != oracle from hand-built F̄ {sig_v} — "
        "the bar dilatation Ĵ or the (Ĵ/J)^(1/3) scaling is wrong"
    )


# --------------------------------------------------------------------------- #
#  Objectivity through the F-bar path: a pure rigid rotation is stress-free.    #
#  (F=R ⇒ J=J₀=1 ⇒ F̄=F=R ⇒ σ̄=0 and the eq 15.10 coupling vanishes; this also    #
#  confirms J₀=1 uniformly, i.e. the centroid evaluation is sane on a rotation.)  #
# --------------------------------------------------------------------------- #
def test_fbar_rigid_rotation_is_stress_free():
    R = _rot(0.9, 0.5)
    u = _affine_disp(R)
    assert _impose_and_solve(u, True) == 0
    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 24
    assert max(abs(v) for v in s) <= 1.0e-7 * _E, "F-bar rigid rotation induced stress"
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, "F-bar rigid rotation induced force"


# --------------------------------------------------------------------------- #
#  The two variants are DISTINCT: under a non-homogeneous state the centroid     #
#  gradient G₀ and the volume-averaged Ḡ differ, so their tangents differ. (If    #
#  they were identical the -fbar switch would be a no-op.)                       #
# --------------------------------------------------------------------------- #
def test_fbar_centroid_and_mean_dilatation_differ():
    Fimp = [[1.12, 0.06, -0.04],
            [0.03, 0.93, 0.05],
            [-0.02, 0.04, 1.08]]
    u = _affine_disp(Fimp, wiggle=True)           # non-affine ⇒ G₀ ≠ Ḡ
    K_cen = _element_tangent(u, True, fbar="centroid")
    K_mean = _element_tangent(u, True, fbar="mean_dilatation")
    scale = np.abs(K_cen).max()
    assert np.abs(K_cen - K_mean).max() > 1.0e-6 * scale, (
        "centroid and mean_dilatation F-bar tangents are identical on a "
        "non-homogeneous state — the -fbar switch is not taking effect"
    )


# --------------------------------------------------------------------------- #
#  Consistency payoff of mean_dilatation: at SMALL strain it reduces to the      #
#  element's volume-averaged small-strain -bbar (-geom linear -bbar), which the   #
#  CENTROID variant does NOT (centroid samples one point, not the volume mean).   #
#  This is the reason mean_dilatation exists.                                     #
# --------------------------------------------------------------------------- #
def _linear_bbar_force(u, nu):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)   # small-strain mean-dilatation bbar
    ops.element("BezierTet10", 1, *range(1, 11), 1, "-bbar", "-geom", "linear")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _COORD:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "linear -bbar solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def test_fbar_mean_dilatation_reduces_to_linear_bbar_at_small_strain():
    nu = _NU
    eps = 1.0e-5
    Fm = [[1.0 + eps, 0.5 * eps, 0.0],
          [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
          [0.0, 0.0, 1.0 + 0.2 * eps]]
    u = _affine_disp(Fm, wiggle=False)            # tiny pure strain (small-strain regime)

    f_lin_bbar = _linear_bbar_force(u, nu)        # small-strain volume-averaged B-bar
    f_mean = _resisting_force(u, True, nu, "mean_dilatation")
    scale = max(np.abs(f_lin_bbar).max(), 1.0e-30)
    # mean_dilatation finite ⇒ small-strain mean-dilatation B-bar as F→I (O(strain) gap)
    assert np.abs(f_mean - f_lin_bbar).max() <= 1.0e-3 * scale, (
        "mean_dilatation F-bar did not reduce to the small-strain -bbar at small "
        f"strain: max diff {np.abs(f_mean - f_lin_bbar).max():.3e} vs scale {scale:.3e}"
    )


# --------------------------------------------------------------------------- #
#  Serialization: the -fbar variant (iData slot 16) must survive recvSelf, or a  #
#  parallel/DB worker silently reverts mean_dilatation → centroid. nodeDisp is   #
#  sp-pinned and vacuous here, so we probe eleResponse('stiff') — the mean and   #
#  centroid tangents differ, so a reverted fbarMode changes it.                  #
# --------------------------------------------------------------------------- #
def test_fbar_mode_serialization_roundtrip():
    Fimp = [[1.10, 0.05, -0.03],
            [0.02, 0.95, 0.04],
            [-0.02, 0.03, 1.07]]
    u = _affine_disp(Fimp, wiggle=True)

    def build():
        assert _impose_and_solve(u, True, _NU, "mean_dilatation") == 0, "roundtrip build failed"

    build()
    probe = list(_COORD)
    assert max(abs(ops.nodeDisp(t, 1)) for t in probe) > 1.0e-6
    database_roundtrip(build, probe, ndf=3,
                       probe_fn=lambda: ops.eleResponse(1, "stiff"))
