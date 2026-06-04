"""LadrunoQuad (classTag 33006) & LadrunoCST (classTag 33007) — Zone-A battery.

ADR 25, Phase 1: small-strain (geometrically linear) ``-geom linear``. The
headline gate is bit-for-bit overlap with the upstream elements:

  * LadrunoQuad ``-formulation std`` reproduces ``FourNodeQuad``  (~1e-9)
  * LadrunoCST                       reproduces ``Tri31``         (~1e-9)

On top of that:
  * rank / zero-energy: a free quad/CST has exactly 3 rigid-body modes (2D) for
    every shipped formulation — no spurious mechanisms;
  * volumetric-locking relief: near-incompressible (nu=0.4999) PlaneStrain, the
    quad's ``bbar`` is strictly more flexible than ``std`` (std locks);
  * constant-strain patch: every Gauss point reports the closed-form stress;
  * reserved formulations ``ssp``/``eas`` are refused at the factory (ADR 25
    P2/P3) and ``bbar`` is PlaneStrain-only.

Plan: Ladruno_implementation/25_ladruno_plane_elements_adr.md.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress

pytestmark = [pytest.mark.zone_a]

# A deliberately distorted (positive-Jacobian) quad. CCW node order.
_QNODES = {
    1: (0.00, 0.00),
    2: (1.00, 0.10),
    3: (1.10, 1.00),
    4: (0.05, 0.95),
}
_QCONN = [1, 2, 3, 4]
_QBASE = [1, 2]            # fully fixed
_QLOADS = {                # mixed so the full 8x8 K is exercised
    3: (10.0, -4.0),
    4: (3.0, 7.0),
}

# A distorted triangle.
_TNODES = {
    1: (0.00, 0.00),
    2: (1.00, 0.05),
    3: (0.10, 1.00),
}
_TCONN = [1, 2, 3]
_TBASE = [1, 2]
_TLOADS = {3: (6.0, -5.0)}

_THK = 0.7


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _static_solve(nodes, conn, base, loads, place_fn, E=1000.0, nu=0.3):
    """Build a single element via place_fn(matTag), static-solve, return
    (disps, stresses, forces)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in nodes.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in base:
        ops.fix(n, 1, 1)

    place_fn(1)  # matTag = 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy) in loads.items():
        ops.load(n, fx, fy)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    disps = []
    for n in conn:
        disps.extend(ops.nodeDisp(n))
    forces = list(ops.eleResponse(1, "forces"))
    stresses = list(ops.eleResponse(1, "stresses"))
    return disps, stresses, forces


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x - y):.3e} > {tol:.3e})"
        )


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


# --------------------------------------------------------------------------
# headline regression gates
# --------------------------------------------------------------------------
def test_quad_std_matches_FourNodeQuad():
    """LadrunoQuad -formulation std reproduces upstream FourNodeQuad to ~1e-9.
    (Upstream FourNodeQuad is registered under the command name 'quad'.)"""
    ref = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("quad", 1, *_QCONN, _THK, "PlaneStrain", m))
    ours = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain",
                                              "-formulation", "std"))
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


def test_quad_std_matches_FourNodeQuad_planestress():
    ref = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("quad", 1, *_QCONN, _THK, "PlaneStress", m))
    ours = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStress"))
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")


def test_quad_default_formulation_is_std():
    ref = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-formulation", "std"))
    bare = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m, "-thick", _THK))
    _assert_close(bare[0], ref[0], "disp")


def test_cst_matches_Tri31():
    """LadrunoCST reproduces upstream Tri31 to ~1e-9."""
    ref = _static_solve(_TNODES, _TCONN, _TBASE, _TLOADS,
                        lambda m: ops.element("Tri31", 1, *_TCONN, _THK, "PlaneStrain", m))
    ours = _static_solve(_TNODES, _TCONN, _TBASE, _TLOADS,
                        lambda m: ops.element("LadrunoCST", 1, *_TCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain"))
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


# --------------------------------------------------------------------------
# rank / zero-energy (no spurious mechanisms)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar"])
def test_quad_no_spurious_modes(formulation):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QNODES.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoQuad", 1, *_QCONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", formulation)
    assert_zero_energy(list(_QNODES.keys()), ndf=2, ndm=2)


def test_cst_no_spurious_modes():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _TNODES.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoCST", 1, *_TCONN, 1, "-thick", _THK, "-type", "PlaneStrain")
    assert_zero_energy(list(_TNODES.keys()), ndf=2, ndm=2)


# --------------------------------------------------------------------------
# volumetric-locking relief: bbar more flexible than std at nu -> 0.5
# --------------------------------------------------------------------------
def test_bbar_relieves_volumetric_locking():
    nu = 0.4999
    std = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain",
                                              "-formulation", "std"), nu=nu)
    bbar = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain",
                                              "-formulation", "bbar"), nu=nu)
    # bbar must be strictly more flexible (larger displacement) — std locks.
    assert _norm(bbar[0]) > _norm(std[0]) * (1.0 + 1e-6), (
        f"bbar |u|={_norm(bbar[0]):.6e} not > std |u|={_norm(std[0]):.6e} "
        "(bbar should relieve volumetric locking)"
    )


# --------------------------------------------------------------------------
# constant-strain patch: closed-form stress at every Gauss point
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar"])
def test_quad_constant_stress_patch(formulation):
    """Constant-stress (force-driven) patch test on a unit square: x-fix the left
    edge, y-fix the bottom edge, pull the right edge with a uniform x-load. The
    exact field u = exx*x, v = eyy*y (linear) gives uniform uniaxial stress
    sxx = P/(h*t), syy = sxy = 0. A complete element reproduces it identically at
    every Gauss point. Both std and bbar must pass (B-bar preserves completeness)."""
    E, nu, t, P = 1000.0, 0.25, 0.7, 5.0
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-thick", t,
                "-type", "PlaneStrain", "-formulation", formulation)
    ops.fix(1, 1, 1)   # left-bottom: u=v=0
    ops.fix(4, 1, 0)   # left-top:    u=0
    ops.fix(2, 0, 1)   # right-bottom: v=0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.5 * P, 0.0)   # right edge (x=1), consistent split
    ops.load(3, 0.5 * P, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    sxx = P / (1.0 * t)   # total load over height*thickness
    check_constant_stress(1, 4, [sxx, 0.0, 0.0], E, sxx / E)


# --------------------------------------------------------------------------
# reserved formulations / option guards refused in this build
# --------------------------------------------------------------------------
def _quad_refused(*extra):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QNODES.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    try:
        ops.element("LadrunoQuad", 1, *_QCONN, 1, "-thick", _THK, *extra)
    except Exception:
        pass  # raising is an acceptable refusal
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    return 1 not in tags


def test_formulation_ssp_reserved():
    assert _quad_refused("-formulation", "ssp"), "ssp must be refused (ADR 25 P2)"


def test_formulation_eas_reserved():
    assert _quad_refused("-formulation", "eas"), "eas must be refused (ADR 25 P3)"


def test_bbar_planestress_refused():
    assert _quad_refused("-type", "PlaneStress", "-formulation", "bbar"), (
        "bbar is PlaneStrain-only"
    )
