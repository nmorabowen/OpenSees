"""LadrunoDispBeamColumn3d (Ladruno regularized displacement-based fiber frame,
3D, classTag 33014) — Zone-A element battery.

3D sibling of LadrunoDispBeamColumn2d: same Tier-1 per-IP characteristic-length
(`lch`) channel for crack-band regularization, plus the `-nl` ½(θy²+θz²) biaxial
bowing strain (DispBeamColumnNL3d). Driven through the same `LadrunoDispBeamColumn`
command (ndm-dispatch picks the 3D class for ndm 3 / ndf 6).

Anchors: reduce-to-stock (elastic == dispBeamColumn), -lch accept/reject, large
displacement via Corotational == stock, and the -nl biaxial bowing (reduces to
linear at small deformation; stiffens under axial restraint at finite rotation).

Plan: Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# elastic 3D section / member properties
_E = 29000.0
_G = 11000.0
_A = 10.0
_Iz = 100.0
_Iy = 80.0
_J = 50.0
_L = 100.0


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def _section_transf():
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)


def _cantilever_tip(ele_name, transf="Linear", n=2):
    _model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    for i in range(1, n + 1):
        ops.element(ele_name, i, i, i + 1, 1, 1)
    tip = n + 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # mixed 6-dof tip load exercising both bending axes + torsion
    ops.load(tip, 1.0e3, -2.0e3, 1.5e3, 3.0e2, 4.0e2, -2.5e2)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(tip, d) for d in range(1, 7)]


def test_reduce_to_stock_elastic_3d():
    d_lad = _cantilever_tip("LadrunoDispBeamColumn")
    d_ref = _cantilever_tip("dispBeamColumn")
    for a, b in zip(d_lad, d_ref):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def _try_build_lch(lch_args):
    _model()
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    _section_transf()
    try:
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, *lch_args)
    except Exception:
        return False
    return 1 in ops.getEleTags()


@pytest.mark.parametrize("args", [(), ("-lch", "ip"), ("-lch", "element"), ("-lch", "12.5"), ("-nl",)])
def test_lch_nl_accepts_valid_3d(args):
    assert _try_build_lch(args) is True


@pytest.mark.parametrize("bad", ["inf", "nan", "-5", "0"])
def test_lch_rejects_invalid_3d(bad):
    assert _try_build_lch(("-lch", bad)) is False


# ---------------------------------------------------------------------------
# Large displacement via Corotational == stock dispBeamColumn (3D)
# ---------------------------------------------------------------------------
def _corot_largedef_tip(ele_name, n=4, nsteps=30):
    _model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Corotational", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    for i in range(1, n + 1):
        ops.element(ele_name, i, i, i + 1, 1, 1)
    tip = n + 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 0.0, -300.0, 0.0, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0
    return [ops.nodeDisp(tip, d) for d in range(1, 7)]


def test_corotational_large_displacement_matches_stock_3d():
    d_lad = _corot_largedef_tip("LadrunoDispBeamColumn")
    d_ref = _corot_largedef_tip("dispBeamColumn")
    assert abs(d_lad[1]) > 5.0  # genuinely large deflection
    for a, b in zip(d_lad, d_ref):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# -nl ½(θy²+θz²) bowing (Linear transform, axially restrained)
# ---------------------------------------------------------------------------
def _lin_transf_tip(nl, P, n=4, nsteps=20):
    _model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    tip = n + 1
    ops.fix(tip, 1, 0, 0, 0, 0, 0)  # restrain axial DOF at the free end
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    extra = ("-nl",) if nl else ()
    for i in range(1, n + 1):
        ops.element("LadrunoDispBeamColumn", i, i, i + 1, 1, 1, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 0.0, -float(P), 0.0, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0
    return ops.nodeDisp(tip, 2)


def test_nl_reduces_to_linear_at_small_deformation_3d():
    assert math.isclose(_lin_transf_tip(True, 1.0), _lin_transf_tip(False, 1.0), rel_tol=1e-5)


def test_nl_stiffens_under_finite_rotation_3d():
    d_nl = _lin_transf_tip(True, 300.0)
    d_lin = _lin_transf_tip(False, 300.0)
    assert abs(d_lin) > 5.0
    assert abs(d_nl) < 0.98 * abs(d_lin)
