"""StagedStrainNDMaterial — small-strain, dimension-general staged-activation battery.

StagedStrain is the small-strain sibling of InitDefGrad: it makes a small-strain
continuum element born STRESS-FREE at the deformed geometry when appended mid-stage.
At the first setTrialStrain it captures the birth strain ε0 and thereafter feeds the
inner material the relative strain ε − ε0 (additive; exact for small strain). It fixes
the three ways the upstream InitStrain bites staged use: 3D-only (it is order-general,
2D+3D), fixed ε0 (it auto-captures), and additive-prestrain intent (it subtracts the
birth strain). It also composes with InitStrain — nested OUTERMOST — for staged-birth-
plus-prestrain.

Driven through real small-strain elements (the material has no Python setTrialStrain):
FourNodeQuad (2D, PlaneStrain/PlaneStress) and stdBrick (3D).

Design note: Ladruno_implementation/staged_deformation_gradiend.md.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_E = 200.0
_NU = 0.3

# 3D: a distorted positive-Jacobian hex.
_N3 = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.05, 0.00), 3: (1.05, 1.00, 0.00),
    4: (0.00, 0.95, 0.05), 5: (0.00, 0.00, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.00, 1.05, 1.00), 8: (0.05, 1.00, 1.00),
}
_C3 = [1, 2, 3, 4, 5, 6, 7, 8]

# 2D: a single slightly-distorted quad.
_N2 = {1: (0.0, 0.0), 2: (1.0, 0.05), 3: (1.05, 1.0), 4: (0.0, 0.95)}
_C2 = [1, 2, 3, 4]


# --------------------------------------------------------------------------- #
#  Model builders (single element), parametrized by dimension.                #
# --------------------------------------------------------------------------- #
def _mk_inner(itag, kind, planar):
    """Build the inner small-strain nDMaterial; return its tag."""
    if kind == "elastic":
        ops.nDMaterial("ElasticIsotropic", itag, _E, _NU)
    elif kind == "j2":
        # small kinematic/iso hardening J2 so yielding is exercised
        ops.nDMaterial("J2Plasticity", itag, _E / (3 * (1 - 2 * _NU)),
                       _E / (2 * (1 + _NU)), 0.20, 0.30, 0.0, 1.0)
    return itag


def _wrap(wtag, itag, spec):
    kind = spec[0]
    if kind == "bare":
        return itag
    if kind == "noinit":
        ops.nDMaterial("StagedStrain", wtag, itag, "-noInit")
    elif kind == "auto":
        ops.nDMaterial("StagedStrain", wtag, itag)
    elif kind == "eps0":
        ops.nDMaterial("StagedStrain", wtag, itag, "-eps0", *[float(v) for v in spec[1]])
    elif kind == "compose":  # StagedStrain( InitStrain(inner, eps_pre) )  -- Staged OUTER
        pre = spec[1]
        ops.nDMaterial("InitStrain", wtag + 1, itag, *[float(v) for v in pre])
        ops.nDMaterial("StagedStrain", wtag, wtag + 1)
    else:
        raise ValueError(kind)
    return wtag


def _build(dim, spec, inner="elastic", eltag=1, mbase=100):
    nodes = _N3 if dim == 3 else _N2
    conn = _C3 if dim == 3 else _C2
    itag = _mk_inner(mbase + 1, inner, dim == 2)
    mtag = _wrap(mbase + 3, itag, spec)
    if dim == 3:
        ops.element("stdBrick", eltag, *conn, mtag)
    else:
        ops.element("quad", eltag, *conn, 1.0, "PlaneStrain", mtag)
    return eltag


def _setup(dim):
    ops.wipe()
    ops.model("basic", "-ndm", dim, "-ndf", dim)
    nodes = _N3 if dim == 3 else _N2
    for t, xyz in nodes.items():
        ops.node(t, *xyz)


def _affine(dim, strain):
    """Uniform infinitesimal strain field -> nodal displacements u = ε·X (symmetric ε
    given as a small matrix). Returns the flat dof vector."""
    nodes = _N3 if dim == 3 else _N2
    eps = np.asarray(strain, float)
    u = np.zeros(len(nodes) * dim)
    for t, xyz in nodes.items():
        X = np.array(xyz)
        u[(t - 1) * dim:(t - 1) * dim + dim] = eps @ X
    return u


def _impose_solve(dim, u):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    nodes = _N3 if dim == 3 else _N2
    for t in nodes:
        b = (t - 1) * dim
        for d in range(dim):
            ops.sp(t, d + 1, float(u[b + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _stresses(dim, spec, u, inner="elastic"):
    _setup(dim)
    _build(dim, spec, inner)
    assert _impose_solve(dim, u) == 0, "solve failed"
    s = ops.eleResponse(1, "stresses")
    ncomp = 6 if dim == 3 else 3
    ngp = 8 if dim == 3 else 4
    assert s and len(s) == ncomp * ngp, "stresses response missing"
    return np.array(s, float).reshape(ngp, ncomp)


def _force(dim, spec, u, inner="elastic"):
    _setup(dim)
    _build(dim, spec, inner)
    assert _impose_solve(dim, u) == 0, "solve failed"
    return np.array(ops.eleForce(1), float)


# A small (small-strain) symmetric strain state per dimension.
_EPS3 = [[1.0e-3, 2.0e-4, -1.0e-4],
         [2.0e-4, -5.0e-4, 1.5e-4],
         [-1.0e-4, 1.5e-4, 8.0e-4]]
_EPS2 = [[1.0e-3, 3.0e-4], [3.0e-4, -6.0e-4]]


def _eps(dim):
    return _EPS3 if dim == 3 else _EPS2


# --------------------------------------------------------------------------- #
#  -noInit is a pure pass-through, in 2D and 3D.                              #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dim", [2, 3])
def test_noinit_is_passthrough(dim):
    u = _affine(dim, _eps(dim))
    s_bare = _stresses(dim, ("bare",), u)
    s_wrap = _stresses(dim, ("noinit",), u)
    f_bare = _force(dim, ("bare",), u)
    f_wrap = _force(dim, ("noinit",), u)
    sc = max(np.abs(s_bare).max(), 1.0)
    fc = max(np.abs(f_bare).max(), 1.0)
    assert np.abs(s_wrap - s_bare).max() <= 1.0e-10 * sc, "-noInit changed stress"
    assert np.abs(f_wrap - f_bare).max() <= 1.0e-10 * fc, "-noInit changed force"


# --------------------------------------------------------------------------- #
#  Explicit -eps0 reduces to the relative strain: StagedStrain(-eps0=g) at ε    #
#  equals the bare inner at ε − g.                                            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dim", [2, 3])
def test_explicit_eps0_reduces_to_relative(dim):
    eps = np.asarray(_eps(dim), float)
    g_mat = 0.5 * eps                       # a birth strain to subtract
    ncomp = 6 if dim == 3 else 3
    # voigt order matches the element 'stresses'/strain: 3D {11,22,33,12,23,31},
    # 2D {11,22,12}. Build g in that engineering-shear voigt form.
    if dim == 3:
        g = [g_mat[0, 0], g_mat[1, 1], g_mat[2, 2],
             2 * g_mat[0, 1], 2 * g_mat[1, 2], 2 * g_mat[2, 0]]
    else:
        g = [g_mat[0, 0], g_mat[1, 1], 2 * g_mat[0, 1]]
    u = _affine(dim, eps)
    u_rel = _affine(dim, eps - g_mat)
    s_wrap = _stresses(dim, ("eps0", g), u)
    s_bare = _stresses(dim, ("bare",), u_rel)
    sc = max(np.abs(s_bare).max(), 1.0)
    assert np.abs(s_wrap - s_bare).max() <= 1.0e-8 * sc, (
        f"explicit-eps0 wrapper != bare inner at relative strain (dim {dim})"
    )


# --------------------------------------------------------------------------- #
#  Genuine two-phase staged construction, in 2D AND 3D: deform+commit a holder, #
#  append an auto-capture StagedStrain element on the same nodes -> born ~0.    #
#  (This also IS the dimensional-generality proof: order 3 capture for the quad, #
#   order 6 for the brick, both through one class.)                            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dim", [2, 3])
def test_two_phase_auto_capture_births_stress_free(dim):
    u = _affine(dim, _eps(dim))
    _setup(dim)
    _build(dim, ("bare",), eltag=1, mbase=110)        # holder
    assert _impose_solve(dim, u) == 0, "phase-1 holder solve failed"
    f_holder = np.array(ops.eleForce(1), float)
    assert np.abs(f_holder).max() > 1.0e-6 * _E, "holder should be stressed"

    _build(dim, ("auto",), eltag=2, mbase=120)        # appended, same nodes
    ops.integrator("LoadControl", 0.0)                # no extra load
    assert ops.analyze(1) == 0, "phase-2 staged-append solve failed"

    f_new = np.array(ops.eleForce(2), float)
    assert np.abs(f_new).max() <= 1.0e-7 * _E, (
        f"appended element not born stress-free (dim {dim}): {np.abs(f_new).max():.3e}"
    )


# --------------------------------------------------------------------------- #
#  PlaneStress also captures correctly (order-3, different inner view).        #
# --------------------------------------------------------------------------- #
def test_plane_stress_birth_stress_free():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, xyz in _N2.items():
        ops.node(t, *xyz)
    ops.nDMaterial("ElasticIsotropic", 111, _E, _NU)
    ops.element("quad", 1, *_C2, 1.0, "PlaneStress", 111)
    u = _affine(2, _EPS2)
    assert _impose_solve(2, u) == 0
    f_holder = np.array(ops.eleForce(1), float)
    assert np.abs(f_holder).max() > 1.0e-6 * _E

    ops.nDMaterial("ElasticIsotropic", 121, _E, _NU)
    ops.nDMaterial("StagedStrain", 123, 121)
    ops.element("quad", 2, *_C2, 1.0, "PlaneStress", 123)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    f_new = np.array(ops.eleForce(2), float)
    assert np.abs(f_new).max() <= 1.0e-7 * _E, "PlaneStress staged birth not stress-free"


# --------------------------------------------------------------------------- #
#  J2 inner: the appended element is born GENUINELY VIRGIN (zero stress and no   #
#  inherited plastic strain), then yields relative to birth.                   #
# --------------------------------------------------------------------------- #
def test_j2_inner_born_virgin():
    # Phase 1: drive a holder well past yield so a naively-born element would       #
    # inherit large plastic strain.                                               #
    big = [[5.0e-3, 0, 0], [0, -2.0e-3, 0], [0, 0, -2.0e-3]]
    u = _affine(3, big)
    _setup(3)
    _mk_inner(110 + 1, "j2", False)
    ops.element("stdBrick", 1, *_C3, 111)             # bare J2 holder
    assert _impose_solve(3, u) == 0
    s_holder = np.array(ops.eleResponse(1, "stresses"), float)
    # E=200, strain ~5e-3, yield 0.2 -> yielded holder carries stress ~O(0.2-0.4)
    assert np.abs(s_holder).max() > 0.1, "holder should carry real stress"

    _mk_inner(120 + 1, "j2", False)
    ops.nDMaterial("StagedStrain", 123, 121)
    ops.element("stdBrick", 2, *_C3, 123)             # appended staged J2
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    s_new = np.array(ops.eleResponse(2, "stresses"), float)
    assert np.abs(s_new).max() <= 1.0e-6 * _E, (
        f"staged J2 element not born virgin: stress {np.abs(s_new).max():.3e}"
    )


# --------------------------------------------------------------------------- #
#  Composition with the upstream InitStrain (StagedStrain OUTER): born carrying #
#  exactly the prestress σ(ε_pre), INDEPENDENT of the birth deformation.       #
# --------------------------------------------------------------------------- #
def test_composition_with_initstrain_keeps_only_prestrain():
    pre = [1.0e-3, 1.0e-3, 1.0e-3, 0.0, 0.0, 0.0]     # volumetric prestrain (3D voigt)

    def born_stress(birth_eps):
        u = _affine(3, birth_eps)
        _setup(3)
        _build(3, ("bare",), eltag=1, mbase=110)      # holder deforms the nodes
        assert _impose_solve(3, u) == 0
        _build(3, ("compose", pre), eltag=2, mbase=120)
        ops.integrator("LoadControl", 0.0)
        assert ops.analyze(1) == 0
        return np.array(ops.eleResponse(2, "stresses"), float).reshape(8, 6)[0]

    s_small = born_stress([[1.0e-3, 0, 0], [0, 0, 0], [0, 0, 0]])
    s_large = born_stress([[4.0e-3, 1.0e-3, 0], [1.0e-3, -2.0e-3, 0], [0, 0, 1.0e-3]])

    # The born stress must be the SAME for both very different birth states (only the
    # prestrain survives the staged re-reference) -> proves StagedStrain captured the
    # geometric birth strain and InitStrain's prestrain rode through.
    sc = max(np.abs(s_small).max(), 1.0)
    assert np.abs(s_small - s_large).max() <= 1.0e-6 * sc, (
        "composed born stress depends on birth deformation -> staged re-reference "
        "failed (or nesting order wrong)"
    )
    # And it is genuinely nonzero (the prestress σ(ε_pre) ~ K·tr(ε_pre) ~ 0.5), not zero.
    assert np.abs(s_small).max() > 0.05, "prestress vanished -> InitStrain cancelled"


# --------------------------------------------------------------------------- #
#  Graceful failure where InitStrain would exit(-1): a missing inner tag makes  #
#  the command fail (return null), not kill the interpreter.                   #
# --------------------------------------------------------------------------- #
def test_missing_inner_fails_gracefully():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    try:
        ops.nDMaterial("StagedStrain", 5, 999)        # inner 999 does not exist
    except Exception:
        pass
    # The point is simply that we got here — the interpreter is alive.
    ops.node(1, 0.0, 0.0, 0.0)
    assert ops.nodeCoord(1) == [0.0, 0.0, 0.0]


# --------------------------------------------------------------------------- #
#  An explicit -eps0 whose size does not match the element's order is gracefully #
#  dropped (with a warning) and the element falls back to AUTO-capture — so it   #
#  is still born stress-free in a two-phase staged build, never crashing.       #
#  (Review fix: getCopy warns instead of silently ignoring the prescribed eps0.) #
# --------------------------------------------------------------------------- #
def test_mismatched_explicit_eps0_falls_back_to_autocapture_2d():
    Fdef = [[8.0e-4, 2.0e-4], [2.0e-4, -5.0e-4]]
    u = _affine(2, Fdef)
    _setup(2)
    # Phase 1: bare PlaneStrain holder deforms + commits.
    ops.nDMaterial("ElasticIsotropic", 111, _E, _NU)
    ops.element("quad", 1, *_C2, 1.0, "PlaneStrain", 111)
    assert _impose_solve(2, u) == 0
    assert np.abs(np.array(ops.eleForce(1), float)).max() > 1.0e-6 * _E

    # Phase 2: append a quad whose StagedStrain has a WRONG-size explicit -eps0
    # (6 components for an order-3 PlaneStrain element). getCopy("PlaneStrain")
    # drops it (warns) -> auto-capture at birth -> born stress-free.
    ops.nDMaterial("ElasticIsotropic", 121, _E, _NU)
    ops.nDMaterial("StagedStrain", 123, 121, "-eps0", 1e-3, 1e-3, 1e-3, 0.0, 0.0, 0.0)
    ops.element("quad", 2, *_C2, 1.0, "PlaneStrain", 123)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0, "mismatched-eps0 staged append failed to solve"
    f_new = np.array(ops.eleForce(2), float)
    assert np.abs(f_new).max() <= 1.0e-7 * _E, (
        f"mismatched-eps0 element not born stress-free: {np.abs(f_new).max():.3e}"
    )
