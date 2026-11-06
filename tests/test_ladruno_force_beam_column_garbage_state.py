"""Zone-A regression: `LadrunoForceBeamColumn2d/3d` never accept or commit a state that is not
a solution of the element equations.

Found by frame-ladruno-validation (recipe 02, blind cyclic RC columns): three runs committed
tip axial displacements of 1e11 - 1e13 mm. Two element defects let such states through.

1. Blind energy test. `evalBulkState` accepted on `|dW| < tol`, `dW = dv . kv . dv`. Once a
   section sits on a flat plateau (crushed concrete, yielded bars, tangent ~ 1e-10), `fs` is huge
   and `kv -> 0`: the linearized correction `fs (Ss - Ssr)` inside `vr` absorbs any target `v`,
   so dW passes while the sections carry the wrong force at deformations incompatible with `v`.
   The element now also requires the section deformations actually set to integrate to `v`.

2. Ignored update failure. `TransformationConstraintHandler::enforceSPs` updates every element
   when it applies the SP values and drops the return code. When that update failed, the sections
   kept the failed iterate, the Newton correction was zero, and `update()`'s
   `|dv| <= DBL_EPSILON` early return reported success: the step converged with the element force
   unchanged and the failed section states were committed. A failed update now disables the early
   returns, puts the sections back on the last solved state and refuses the commit.

Every case drives a cantilever with prescribed tip axial displacement and rotation(s), lateral DOFs
free, from rest, and checks one property: either the analyze fails and the element is back on its
committed state, or the sections are in equilibrium with the element force.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "dist" / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import opensees as ops  # noqa: E402

pytestmark = pytest.mark.zone_a

L = 1000.0


def _materials():
    ops.uniaxialMaterial("Concrete02", 1, -30.0, -0.002, -6.0, -0.006, 0.1, 3.0, 1500.0)
    ops.uniaxialMaterial("Steel01", 2, 400.0, 200000.0, 1.0e-9)   # flat yield plateau


def _analysis():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-8, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _cantilever_2d(elem, iter_args=("-iter", 10, 1.0e-4), elastic=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, L)
    ops.fix(1, 1, 1, 1)
    if elastic:
        ops.section("Elastic", 1, 30000.0, 9.0e4, 6.75e8)
    else:
        _materials()
        ops.section("Fiber", 1)
        ops.patch("rect", 1, 10, 1, -150.0, -150.0, 150.0, 150.0)
        for y in (-120.0, 120.0):
            ops.fiber(y, 0.0, 600.0, 2)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(elem, 1, 1, 2, 1, 1, *iter_args)
    ops.timeSeries("Linear", 1)
    _analysis()
    return (2, 3)          # prescribed tip DOFs: axial, rotation


def _cantilever_3d(elem, iter_args=("-iter", 10, 1.0e-4), elastic=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    if elastic:
        ops.section("Elastic", 1, 30000.0, 9.0e4, 6.75e8, 6.75e8, 12500.0, 1.1e9)
    else:
        _materials()
        ops.section("Fiber", 1, "-GJ", 1.0e12)
        ops.patch("rect", 1, 10, 10, -150.0, -150.0, 150.0, 150.0)
        for y in (-120.0, 120.0):
            for z in (-120.0, 120.0):
                ops.fiber(y, z, 300.0, 2)
    ops.geomTransf("Linear", 1, 1.0, 0.0, 0.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(elem, 1, 1, 2, 1, 1, *iter_args)
    ops.timeSeries("Linear", 1)
    _analysis()
    ops.pattern("Plain", 9, 1)
    ops.sp(2, 6, 0.0)       # no twist
    return (3, 4, 5)       # prescribed tip DOFs: axial, rotations


def _impose(dofs, values):
    """Prescribe absolute tip values `values` on `dofs` in one LoadControl step."""
    try:
        ops.remove("loadPattern", 1)
    except Exception:  # noqa: BLE001 - first step has no pattern yet
        pass
    ops.pattern("Plain", 1, 1)
    for dof, val in zip(dofs, values):
        ops.sp(2, dof, val)
    ops.setTime(0.0)
    return ops.analyze(1)


def _section_axial(i):
    return ops.eleResponse(1, "section", i, "force")[0]


def _section_def(i):
    return list(ops.eleResponse(1, "section", i, "deformation"))


def _assert_solution_or_reverted(rc, committed_defs):
    """rc == 0: sections in axial equilibrium with the basic force (no span load).
    rc != 0: the element is back on its committed section state."""
    if rc == 0:
        N = ops.eleResponse(1, "basicForce")[0]
        secN = [_section_axial(i) for i in range(1, 6)]
        scale = max([abs(N)] + [abs(s) for s in secN] + [1.0])
        worst = max(abs(s - N) for s in secN)
        assert worst <= 1.0e-6 * scale, (
            f"converged with sections out of equilibrium: N = {N}, section N = {secN}")
    else:
        assert [_section_def(i) for i in range(1, 6)] == committed_defs


# Single jumps from rest found by an element-level random search on the pre-fix build
# (e7d489c8b): every one converged with rc 0 and section forces off by up to 30 %.
CASES_2D = [(-1.1232, 0.0024), (1.5411, 0.004027), (-4.4667, -0.0076)]
CASES_3D = [(-0.11642, -0.0020910, -0.0016523), (0.86219, 0.00051486, 0.0043134),
            (1.05391, -0.00094002, 0.0052974)]


@pytest.mark.parametrize("target", CASES_2D)
def test_lfbc2d_plateau_jump_is_solution_or_failure(target):
    dofs = _cantilever_2d("LadrunoForceBeamColumn")
    committed = [_section_def(i) for i in range(1, 6)]
    rc = _impose(dofs, target)
    _assert_solution_or_reverted(rc, committed)


@pytest.mark.parametrize("target", CASES_3D)
def test_lfbc3d_plateau_jump_is_solution_or_failure(target):
    dofs = _cantilever_3d("LadrunoForceBeamColumn")
    committed = [_section_def(i) for i in range(1, 6)]
    rc = _impose(dofs, target)
    _assert_solution_or_reverted(rc, committed)


@pytest.mark.parametrize("builder", [_cantilever_2d, _cantilever_3d], ids=["2d", "3d"])
def test_lfbc_failed_update_is_not_committed(builder):
    """An element that cannot converge (`-iter 1 1e-300`) on an elastic jump.

    Pre-fix: the update inside enforceSPs failed silently, Newton's zero correction took the
    early return, and the step 'converged' with the tip displaced and zero element force."""
    dofs = builder("LadrunoForceBeamColumn", iter_args=("-iter", 1, 1.0e-300), elastic=True)
    committed = [_section_def(i) for i in range(1, 6)]
    rc = _impose(dofs, (-0.5,) + (0.002,) * (len(dofs) - 1))
    assert rc != 0
    assert ops.nodeDisp(2, dofs[0]) == 0.0
    assert [_section_def(i) for i in range(1, 6)] == committed


@pytest.mark.parametrize("elem", ["LadrunoForceBeamColumn", "forceBeamColumn"])
@pytest.mark.parametrize("builder", [_cantilever_2d, _cantilever_3d], ids=["2d", "3d"])
def test_failed_analyze_restores_committed_section_deformation(builder, elem):
    """Issue B check: revertToLastCommit puts the sections back on the committed state that
    eleResponse(..., 'section', i, 'deformation') reads, as vanilla forceBeamColumn does."""
    dofs = builder(elem)
    assert _impose(dofs, (-0.3,) + (0.001,) * (len(dofs) - 1)) == 0
    committed = [_section_def(i) for i in range(1, 6)]
    assert any(abs(d[0]) > 1.0e-5 for d in committed)
    ops.test("NormDispIncr", 1.0e-8, 1, 0)       # one iteration: this jump cannot converge
    rc = _impose(dofs, (-6.0,) + (0.02,) * (len(dofs) - 1))
    assert rc != 0
    assert [_section_def(i) for i in range(1, 6)] == committed
