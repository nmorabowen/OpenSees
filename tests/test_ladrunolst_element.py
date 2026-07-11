"""LadrunoLST (T6, classTag 33016, ADR 70 P3) — linear-geometry Zone-A battery.

The headline gate is reduce-to-upstream: LadrunoLST ``-formulation std -geom
linear`` reproduces ``SixNodeTri`` (~1e-9) — same quadratic shape functions,
same 3-point interior rule (barycentric (2/3,1/6),(1/6,2/3),(1/6,1/6), w=1/6).

On top of that:
  * rank / zero-energy: a free std T6 has exactly 3 rigid-body modes;
  * the B-BAR REFUTATION PIN (the P3 finding): constant element-mean
    dilatation B-bar/F-bar on the T6 is RANK-DEFICIENT — the two quadratic
    conformal modes (Re/Im z²) carry zero deviatoric strain and zero MEAN
    dilatation, so a free element would get 5 zero-energy modes (stacked-B̄
    rank 7 of the required 9; verified numerically AND on a compiled build
    during P3). LadrunoLST therefore ships std-only and '-formulation bbar'
    must stay parser-refused; this battery pins the refusal;
  * constant-strain patch: every integration point reports the closed-form
    stress under a linear displacement field;
  * LINEAR-strain-field reproduction (the reason the LST exists, ADR 70 §3.3):
    under an exact quadratic displacement field the T6 interpolates the exact
    LINEAR strain at every integration point — the field a constant-strain
    CST cannot carry (this is the mechanism behind off-mesh-line shear bands);
  * charLength = sqrt(2*area) (triangle crack-band convention);
  * parse guards (bbar refused with the rank message, quad-only formulations
    refused).

Plan: Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress

pytestmark = [pytest.mark.zone_a]

# A distorted (positive-Jacobian) T6: straight edges, midsides at midpoints
# (the ADR's distortion caveat — curved/off-center midsides degrade the T6;
# the gates run on the well-posed configuration).
_C = {1: (0.00, 0.00), 2: (2.00, 0.20), 3: (0.30, 1.60)}
_T6N = {
    1: _C[1], 2: _C[2], 3: _C[3],
    4: ((_C[1][0] + _C[2][0]) / 2, (_C[1][1] + _C[2][1]) / 2),   # mid 1-2
    5: ((_C[2][0] + _C[3][0]) / 2, (_C[2][1] + _C[3][1]) / 2),   # mid 2-3
    6: ((_C[3][0] + _C[1][0]) / 2, (_C[3][1] + _C[1][1]) / 2),   # mid 3-1
}
_CONN = [1, 2, 3, 4, 5, 6]
_BASE = [1, 2, 4]              # bottom edge (corner-mid-corner) fully fixed
_LOADS = {                     # mixed so the full 12x12 K is exercised
    3: (6.0, -5.0),
    5: (2.0, 3.0),
    6: (-1.5, 4.0),
}
_THK = 0.7


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x - y):.3e} > {tol:.3e})"
        )


def _static_solve(place_fn, E=1000.0, nu=0.3):
    """Build a single 6-node element via place_fn(matTag), static-solve, return
    (disps, stresses, forces)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in _BASE:
        ops.fix(n, 1, 1)
    place_fn(1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy) in _LOADS.items():
        ops.load(n, fx, fy)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    forces = list(ops.eleResponse(1, "forces"))
    stresses = list(ops.eleResponse(1, "stresses"))
    return disps, stresses, forces


# --------------------------------------------------------------------------
# headline reduce-to-upstream gates
# --------------------------------------------------------------------------
@pytest.mark.parametrize("ptype", ["PlaneStrain", "PlaneStress"])
def test_lst_std_matches_SixNodeTri(ptype):
    """LadrunoLST std/linear == upstream SixNodeTri to ~1e-9 (disps, per-GP
    stresses, resisting forces) on a distorted mesh with mixed loads."""
    # upstream SixNodeTri registers as "tri6n" in the shared functionMap
    ref = _static_solve(
        lambda m: ops.element("tri6n", 1, *_CONN, _THK, ptype, m))
    got = _static_solve(
        lambda m: ops.element("LadrunoLST", 1, *_CONN, m, "-thick", _THK,
                              "-type", ptype))
    _assert_close(got[0], ref[0], f"disps({ptype})")
    _assert_close(got[1], ref[1], f"stresses({ptype})")
    _assert_close(got[2], ref[2], f"forces({ptype})")


def test_lst_default_formulation_is_std():
    """No -formulation flag == std (SixNodeTri overlap)."""
    ref = _static_solve(
        lambda m: ops.element("LadrunoLST", 1, *_CONN, m, "-thick", _THK,
                              "-type", "PlaneStrain", "-formulation", "std"))
    got = _static_solve(
        lambda m: ops.element("LadrunoLST", 1, *_CONN, m, "-thick", _THK,
                              "-type", "PlaneStrain"))
    _assert_close(got[0], ref[0], "disps(default)")


# --------------------------------------------------------------------------
# rank / zero-energy (no spurious mechanisms)
# --------------------------------------------------------------------------
def test_lst_no_spurious_modes():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", "std")
    assert_zero_energy(_CONN, ndf=2, ndm=2)


# --------------------------------------------------------------------------
# patch tests: constant strain (exact) + LINEAR strain (the LST raison d'être)
# --------------------------------------------------------------------------
def _prescribe_field(u_fn, form="std"):
    """Single T6 with EVERY node displacement prescribed to u_fn(X, Y) via the
    Transformation handler — sp constraints condense EXACTLY (no penalty
    error, which sits at ~k/α ≈ 1e-7 and fails the locked 1e-9 patch rtol).
    A fully-prescribed model is a FullGenLinSOE fatal (LEDGER quirk), so a
    dummy grounded spring supplies the one free equation."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", form)
    # dummy free DOF: grounded x-spring far from the element
    ops.node(100, 10.0, 10.0)
    ops.node(101, 10.0, 10.0)
    ops.fix(101, 1, 1)
    ops.fix(100, 0, 1)
    ops.uniaxialMaterial("Elastic", 900, 1.0)
    ops.element("zeroLength", 900, 101, 100, "-mat", 900, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y) in _T6N.items():
        ux, uy = u_fn(x, y)
        ops.sp(tag, 1, ux)
        ops.sp(tag, 2, uy)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def test_lst_constant_stress_patch():
    """Linear displacement field -> every integration point reports the
    closed-form constant plane-strain stress."""
    E, nu = 1000.0, 0.3
    exx, eyy, gxy = 1e-3, -4e-4, 5e-4
    _prescribe_field(lambda x, y: (exx * x + 0.5 * gxy * y,
                                   eyy * y + 0.5 * gxy * x))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    sxx = (lam + 2 * mu) * exx + lam * eyy
    syy = lam * exx + (lam + 2 * mu) * eyy
    sxy = mu * gxy
    check_constant_stress(1, 3, [sxx, syy, sxy], E, exx)


def test_lst_carries_linear_strain_field():
    """Quadratic displacement patch: u = (a x^2, b y^2 + c x y). The T6's
    quadratic interpolation reproduces the exact LINEAR strain field at every
    integration point (eps_xx = 2ax, eps_yy = 2by + cx, gamma = cy). This is
    the capability gap vs the CST (constant strain) that lets the LST carry an
    inclined localization band (ADR 70 §3.3); a wrong shape-function derivative
    or B assembly fails it immediately."""
    a, bq, c = 2e-4, -1.5e-4, 1e-4
    _prescribe_field(lambda x, y: (a * x * x, bq * y * y + c * x * y))

    # integration-point physical coords: midsides at midpoints => affine map
    # from the corners (barycentric L1=s at node1, L2=t at node2, L3 at node3)
    gps = [(2 / 3, 1 / 6), (1 / 6, 2 / 3), (1 / 6, 1 / 6)]
    strains = list(ops.eleResponse(1, "strains"))       # 3 GP x [exx, eyy, gxy]
    for i, (s, t) in enumerate(gps):
        X = s * _C[1][0] + t * _C[2][0] + (1 - s - t) * _C[3][0]
        Y = s * _C[1][1] + t * _C[2][1] + (1 - s - t) * _C[3][1]
        exx_ref = 2 * a * X
        eyy_ref = 2 * bq * Y + c * X
        gxy_ref = c * Y
        got = strains[3 * i:3 * i + 3]
        for comp, (g, r) in enumerate(zip(got, (exx_ref, eyy_ref, gxy_ref))):
            tol = 1e-9 + 1e-8 * max(abs(g), abs(r))
            assert abs(g - r) <= tol, (
                f"GP{i} strain[{comp}]: {g!r} != exact linear field {r!r}"
            )


# --------------------------------------------------------------------------
# HRZ lumped mass: strictly positive (deliberate divergence from SixNodeTri)
# --------------------------------------------------------------------------
def test_lst_hrz_mass_positive_dynamics():
    """The LST uses HRZ lumping (∫ρN² rescaled) instead of upstream
    SixNodeTri's plain N-lumping, which yields EXACTLY ZERO corner masses
    (negative on distorted T6s) — unusable for explicit dynamics. Gate: a
    free-free element with element mass only has a well-posed generalized
    eigenproblem — 12 FINITE eigenvalues, exactly 3 rigid-body zeros. With a
    singular (zero-corner-mass) M this blows up instead."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-rho", 2.5)
    eigs = ops.eigen("-fullGenLapack", 12)
    assert all(math.isfinite(e) for e in eigs), (
        f"non-finite eigenvalues — singular lumped mass? {eigs}"
    )
    lam_max = max(abs(e) for e in eigs)
    nzero = sum(1 for e in eigs if abs(e) < 1e-8 * lam_max)
    assert nzero == 3, f"expected 3 rigid-body modes, got {nzero}: {eigs}"


# --------------------------------------------------------------------------
# consistent quadratic-edge pressure: closed-contour resultant is ZERO
# --------------------------------------------------------------------------
def test_lst_pressure_closed_contour():
    """Uniform -pressure on the full (closed) element boundary must produce a
    ZERO net force resultant (∮ p·n dΓ = 0). This discriminates the fork's
    side-61 fix from upstream SixNodeTri's 'dx61 = x4-x6' copy-paste typo,
    which breaks contour closure and leaks a spurious net force."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-pressure", 3.7)
    f = list(ops.eleResponse(1, "forces"))       # −pressureLoad at zero disp
    fx = sum(f[0::2])
    fy = sum(f[1::2])
    fmag = max(abs(v) for v in f)
    assert fmag > 1e-3, "pressure load should be nonzero at the nodes"
    assert abs(fx) <= 1e-12 + 1e-12 * fmag, f"net Fx {fx} != 0 (contour not closed)"
    assert abs(fy) <= 1e-12 + 1e-12 * fmag, f"net Fy {fy} != 0 (contour not closed)"


# --------------------------------------------------------------------------
# characteristic length
# --------------------------------------------------------------------------
def test_lst_charlen_sqrt_2area():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK, "-type", "PlaneStrain")
    x1, y1 = _C[1]; x2, y2 = _C[2]; x3, y3 = _C[3]
    A = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    lch = ops.eleResponse(1, "charLength")[0]
    assert abs(lch - math.sqrt(2 * A)) <= 1e-12 + 1e-9 * lch


# --------------------------------------------------------------------------
# parse guards
# --------------------------------------------------------------------------
def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)


def test_lst_bbar_refused_rank_deficient():
    """THE refutation pin (ADR 70 P3): '-formulation bbar' must be refused —
    constant element-mean dilatation on the T6 admits the two conformal
    zero-energy modes (5 zero modes on a free element, stacked-B̄ rank 7 < 9).
    If someone re-enables it without a rank-sufficient projection (P4), this
    trips."""
    _fresh()
    try:
        ops.element("LadrunoLST", 1, *_CONN, 1, "-type", "PlaneStrain",
                    "-formulation", "bbar")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), (
        "bbar must stay parser-refused on the LST (rank-deficient, ADR 70 P3)"
    )
    ops.element("LadrunoLST", 2, *_CONN, 1, "-type", "PlaneStrain",
                "-formulation", "std")
    assert 2 in ops.getEleTags()


@pytest.mark.parametrize("form", ["ssp", "eas"])
def test_lst_quad_only_formulations_refused(form):
    """ssp/eas are quad-only (single-point stabilization / Q1E4 modes have no
    T6 analog here) — the factory must refuse, and the interpreter survives."""
    _fresh()
    try:
        ops.element("LadrunoLST", 1, *_CONN, 1, "-type", "PlaneStrain",
                    "-formulation", form)
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), f"{form} must be refused for the LST"
    ops.element("LadrunoLST", 2, *_CONN, 1, "-type", "PlaneStrain")
    assert 2 in ops.getEleTags()


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
