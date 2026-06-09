"""InitStrainNDMaterial — dimension-general fixed-prestrain battery (Ladruno).

Upstream InitStrain (Petracca, 2024) is a fine FIXED additive prestrain (thermal /
lack-of-fit), but it was 3D-ONLY: getOrder() was hardcoded to 6 and getCopy(type)
only answered "ThreeDimensional", so the base NDMaterial::getCopy() returned 0 for
PlaneStrain / AxiSymmetric — the DEFAULT view of LadrunoQuad / LadrunoCST (and the
upstream FourNodeQuad). A 2D element wrapping an InitStrain material therefore got a
null material and failed to construct.

The Ladruno change makes it dimension-general: PlaneStrain / AxiSymmetric views are
built natively (inner = that view of the stored 3D material; eps0 reduced via the
same vmap LadrunoJ2 uses), while PlaneStress / PlateFiber / 3D keep their existing
paths and 3D behaviour stays byte-identical.

These tests pin the headline contract: a PlaneStrain element now BUILDS and IMPOSES
the prestrain, and InitStrain stays exactly ADDITIVE — InitStrain(eps0=g) at element
strain ε equals the bare inner at ε + g — in 2D and 3D. Driven through real elements
(the material has no Python setTrialStrain): FourNodeQuad + LadrunoQuad (2D), stdBrick
(3D).
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

# A small (small-strain) symmetric tensor strain state per dimension.
_EPS3 = [[1.0e-3, 2.0e-4, -1.0e-4],
         [2.0e-4, -5.0e-4, 1.5e-4],
         [-1.0e-4, 1.5e-4, 8.0e-4]]
_EPS2 = [[1.0e-3, 3.0e-4], [3.0e-4, -6.0e-4]]


def _nodes(dim):
    return _N3 if dim == 3 else _N2


def _conn(dim):
    return _C3 if dim == 3 else _C2


def _setup(dim):
    ops.wipe()
    ops.model("basic", "-ndm", dim, "-ndf", dim)
    for t, xyz in _nodes(dim).items():
        ops.node(t, *xyz)


def _affine(dim, tensor_strain):
    """Uniform infinitesimal tensor strain field -> nodal displacements u = ε·X."""
    eps = np.asarray(tensor_strain, float)
    u = np.zeros(len(_nodes(dim)) * dim)
    for t, xyz in _nodes(dim).items():
        u[(t - 1) * dim:(t - 1) * dim + dim] = eps @ np.array(xyz)
    return u


def _impose_solve(dim, u):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _nodes(dim):
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


def _mk_element(dim, eltag, conn, mtag, eltype):
    if dim == 3:
        ops.element("stdBrick", eltag, *conn, mtag)
    elif eltype == "ladruno":
        ops.element("LadrunoQuad", eltag, *conn, mtag, "-type", "PlaneStrain")
    else:
        ops.element("quad", eltag, *conn, 1.0, "PlaneStrain", mtag)


def _voigt_g(dim, g_mat):
    """3D-Voigt engineering-shear initial-strain vector {11,22,33,12,23,31} from a
    tensor matrix. (InitStrain's eps0 input is always full 3D Voigt.)"""
    g = np.asarray(g_mat, float)
    if dim == 3:
        return [g[0, 0], g[1, 1], g[2, 2], 2 * g[0, 1], 2 * g[1, 2], 2 * g[2, 0]]
    return [g[0, 0], g[1, 1], 0.0, 2 * g[0, 1], 0.0, 0.0]


def _stresses(dim, u, wrap_g=None, eltype="upstream", mbase=100):
    """Build one element (optionally InitStrain-wrapped), solve, return GP stresses."""
    _setup(dim)
    itag = mbase + 1
    ops.nDMaterial("ElasticIsotropic", itag, _E, _NU)
    if wrap_g is None:
        mtag = itag
    else:
        mtag = mbase + 3
        ops.nDMaterial("InitStrain", mtag, itag, *[float(v) for v in _voigt_g(dim, wrap_g)])
    _mk_element(dim, 1, _conn(dim), mtag, eltype)
    assert _impose_solve(dim, u) == 0, "solve failed"
    s = ops.eleResponse(1, "stresses")
    ncomp = 6 if dim == 3 else 3
    ngp = 8 if dim == 3 else 4
    assert s and len(s) == ncomp * ngp, "stresses response missing"
    return np.array(s, float).reshape(ngp, ncomp)


def _eps(dim):
    return _EPS3 if dim == 3 else _EPS2


# --------------------------------------------------------------------------- #
#  Headline: a PlaneStrain element wrapping InitStrain now BUILDS and solves.  #
#  (Before the fix this returned a null material and the element failed.)      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("eltype", ["upstream", "ladruno"])
def test_planestrain_initstrain_builds_and_solves(eltype):
    u = _affine(2, _eps(2))
    g = 0.5 * np.asarray(_eps(2), float)
    s = _stresses(2, u, wrap_g=g, eltype=eltype)   # asserts solve == 0 internally
    assert np.isfinite(s).all() and np.abs(s).max() > 0.0, "no stress developed"


# --------------------------------------------------------------------------- #
#  InitStrain stays exactly ADDITIVE in every dimension: InitStrain(g) at ε     #
#  equals the bare inner at ε + g. This is the correctness pin (2D was broken). #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dim,eltype", [(2, "upstream"), (2, "ladruno"), (3, "upstream")])
def test_initstrain_is_additive(dim, eltype):
    eps = np.asarray(_eps(dim), float)
    g_mat = 0.5 * eps
    u = _affine(dim, eps)
    u_plus = _affine(dim, eps + g_mat)             # bare inner at ε + g
    s_wrap = _stresses(dim, u, wrap_g=g_mat, eltype=eltype, mbase=100)
    s_bare = _stresses(dim, u_plus, wrap_g=None, eltype=eltype, mbase=200)
    sc = max(np.abs(s_bare).max(), 1.0)
    assert np.abs(s_wrap - s_bare).max() <= 1.0e-8 * sc, (
        f"InitStrain not additive (dim {dim}, {eltype}): "
        f"{np.abs(s_wrap - s_bare).max():.3e}"
    )


# --------------------------------------------------------------------------- #
#  PlaneStress (the legacy base-wrapper path) still builds and solves -> the    #
#  fix did not regress the views the base could already construct.             #
# --------------------------------------------------------------------------- #
def test_planestress_initstrain_still_builds():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, xyz in _N2.items():
        ops.node(t, *xyz)
    ops.nDMaterial("ElasticIsotropic", 11, _E, _NU)
    ops.nDMaterial("InitStrain", 13, 11, 1e-3, 1e-3, 0.0, 0.0, 0.0, 0.0)
    ops.element("quad", 1, *_C2, 1.0, "PlaneStress", 13)
    assert _impose_solve(2, _affine(2, _EPS2)) == 0, "PlaneStress InitStrain solve failed"
    s = np.array(ops.eleResponse(1, "stresses"), float)
    assert np.isfinite(s).all()


# --------------------------------------------------------------------------- #
#  setParameter eps0 reaches a live PlaneStrain element through the canonical   #
#  3D store: an in-plane component (eps0_11) takes effect; an out-of-plane one  #
#  (eps0_23, absent from the PlaneStrain view) is a graceful no-op, not a crash.#
# --------------------------------------------------------------------------- #
def test_setparameter_eps0_on_planestrain():
    def stress_with_param(name, value):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for t, xyz in _N2.items():
            ops.node(t, *xyz)
        ops.nDMaterial("ElasticIsotropic", 11, _E, _NU)
        ops.nDMaterial("InitStrain", 13, 11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.element("quad", 1, *_C2, 1.0, "PlaneStrain", 13)
        # set the prestrain component on the element's material, then solve at zero disp
        ops.setParameter("-val", float(value), "-ele", 1, name)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in _N2:
            ops.sp(t, 1, 0.0)
            ops.sp(t, 2, 0.0)
        ops.constraints("Lagrange")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-11, 30, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return np.array(ops.eleResponse(1, "stresses"), float)

    s_in = stress_with_param("eps0_11", 1.0e-3)    # in-plane -> develops stress
    s_out = stress_with_param("eps0_23", 1.0e-3)   # out-of-plane in PlaneStrain -> no-op
    assert np.abs(s_in).max() > 1.0e-3 * _E, "eps0_11 did not impose a PlaneStrain prestress"
    assert np.abs(s_out).max() <= 1.0e-9 * _E, (
        f"eps0_23 (absent from PlaneStrain) should be a no-op, got {np.abs(s_out).max():.3e}"
    )
