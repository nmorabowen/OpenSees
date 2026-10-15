"""Prescribed displacements under constraints Auto, sp -subtractInit, and Path -useLast.

  auto_sp        two-bar chain, node 1 fixed, node 3 driven to u3 = 3.0 by an sp,
                 Newton with NormDispIncr. Equal bars give the closed form
                 u2 = u3 / 2 = 1.5 for every constraint handler.

Each case builds its own model after ops.wipe().
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E = 100.0
A = 1.0
L = 1.0


def bar_chain():
    # nodes 1-2-3 on the x axis, two equal trusses, y DOFs fixed
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 2.0 * L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.element("Truss", 2, 2, 3, A, 1)


def static_analysis(handler, test_type="NormDispIncr"):
    if handler == "Penalty":
        ops.constraints("Penalty", 1.0e12, 1.0e12)
    else:
        ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test(test_type, 1.0e-10, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


@pytest.mark.parametrize("handler", ["Auto", "Transformation", "Penalty"])
@pytest.mark.parametrize("test_type", ["NormDispIncr", "EnergyIncr"])
def test_auto_sp_prescribed_displacement(handler, test_type):
    """u3 = 3.0 prescribed; u2 = 1.5 (equal bars). Auto returned u2 = 0.0."""
    bar_chain()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(3, 1, 3.0)
    static_analysis(handler, test_type)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 1) == pytest.approx(3.0, rel=1.0e-6)
    assert ops.nodeDisp(2, 1) == pytest.approx(1.5, rel=1.0e-6)
