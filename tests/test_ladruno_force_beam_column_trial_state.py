"""Zone-A regression: `LadrunoForceBeamColumn2d/3d` iterate from the last TRIAL state.

Two defects, both from `evalBulkState` restarting every update from the COMMITTED
state (vanilla `ForceBeamColumn` iterates from the last trial state):

1. Newton floor. Each update re-solved the whole step to the element tolerance, so
   the resisting force carried a tolerance-sized error that no global iteration
   removed. A plain RC cantilever stalled at |R| ~ 1e-2 N and failed the 4th
   displacement step under `NormDispIncr 1e-8`; vanilla converges in 3 iterations.

2. Committed-state tangent. A load re-application at the committed state
   (`eleLoadDirty`: `DisplacementControl::domainChanged` re-applies the loads and
   calls `formUnbalance` at the start of every `analyze`) re-solved span-loaded
   members at zero increment, pushed the sub-tolerance section residual into the
   fibres and flipped them onto their unloading tangent. The next step's predictor
   was then built on a stiffer tangent (up to 46 % off on the portal below, 29 % on
   the 04a frame of frame-ladruno-validation); vanilla's predictor tangent is the
   last-iteration tangent.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "dist" / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import opensees as ops  # noqa: E402

pytestmark = pytest.mark.zone_a

FC, FY, ES = 35.0, 420.0, 200000.0   # N, mm
EC = 4700.0 * FC ** 0.5


def _materials():
    eps0 = -2.0 * FC / EC
    ops.uniaxialMaterial("Concrete02", 1, -1.2 * FC, 1.2 * eps0, -0.5 * FC, -0.025, 0.10, 3.0, 1500.0)
    ops.uniaxialMaterial("Concrete02", 2, -FC, eps0, -0.35 * FC, -0.018, 0.10, 2.5, 1500.0)
    ops.uniaxialMaterial("Steel02", 3, FY, ES, 0.01, 18.0, 0.925, 0.15)


def _rc_section(tag, b, h, As, nb, three_d=False):
    cov = 40.0
    ys = 0.5 * h - cov
    if three_d:
        ops.section("Fiber", tag, "-GJ", 1.0e13)
    else:
        ops.section("Fiber", tag)
    ops.patch("rect", 1, 10, 2 if three_d else 1, -ys, -0.5 * b + cov, ys, 0.5 * b - cov)
    ops.patch("rect", 2, 2, 1, -0.5 * h, -0.5 * b, -ys, 0.5 * b)
    ops.patch("rect", 2, 2, 1, ys, -0.5 * b, 0.5 * h, 0.5 * b)
    ops.layer("straight", 3, nb, As / nb, -ys, -0.5 * b + cov, -ys, 0.5 * b - cov)
    ops.layer("straight", 3, nb, As / nb, ys, -0.5 * b + cov, ys, 0.5 * b - cov)


def _gravity_then_lateral_pattern(lateral):
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.2)
    ops.analysis("Static")
    assert ops.analyze(5) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for node, values in lateral:
        ops.load(node, *values)


def _cantilever_2d(elem):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 5000.0, 0.0)
    ops.fix(1, 1, 1, 1)
    _materials()
    _rc_section(1, 350.0, 550.0, 1080.0, 3)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(elem, 1, 1, 2, 1, 1, "-iter", 30, 1.0e-10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    _gravity_then_lateral_pattern([(2, (0.0, 1.0, 0.0))])


def _portal_2d(elem):
    """RC portal: P-Delta columns, Linear beam with a -beamUniform gravity load."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 6000.0, 0.0)
    ops.node(3, 0.0, 3000.0)
    ops.node(4, 6000.0, 3000.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    _materials()
    _rc_section(1, 450.0, 450.0, 1600.0, 4)
    _rc_section(2, 300.0, 500.0, 1080.0, 3)
    ops.geomTransf("PDelta", 1)
    ops.geomTransf("Linear", 2)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.beamIntegration("Lobatto", 2, 2, 5)
    ops.element(elem, 1, 1, 3, 1, 1, "-iter", 30, 1.0e-10)
    ops.element(elem, 2, 2, 4, 1, 1, "-iter", 30, 1.0e-10)
    ops.element(elem, 3, 3, 4, 2, 2, "-iter", 30, 1.0e-10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 3, "-type", "-beamUniform", -30.0, 0.0)
    ops.load(3, 0.0, -400.0e3, 0.0)
    ops.load(4, 0.0, -400.0e3, 0.0)
    _gravity_then_lateral_pattern([(3, (1.0, 0.0, 0.0)), (4, (1.0, 0.0, 0.0))])


def _cantilever_3d(elem):
    """Biaxially pushed RC cantilever with an axial load and a biaxial -beamUniform load."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 5000.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    _materials()
    _rc_section(1, 350.0, 550.0, 1080.0, 3, three_d=True)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(elem, 1, 1, 2, 1, 1, "-iter", 30, 1.0e-10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", -20.0, -6.0, 0.0)
    ops.load(2, -200.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
    _gravity_then_lateral_pattern([(2, (0.0, 1.0, 0.3, 0.0, 0.0, 0.0))])


def _push(node, dof, du, nsteps, on_converged=None):
    """Plain Newton, NormDispIncr 1e-8; returns (iterations per step, load factors)."""
    iters, lams = [], []
    for k in range(nsteps):
        ops.test("NormDispIncr", 1.0e-8, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", node, dof, du)
        if ops.analyze(1) != 0:
            pytest.fail(f"Newton failed at displacement step {k + 1} (iterations so far {iters})")
        iters.append(ops.testIter())
        lams.append(ops.getLoadFactor(2))
        if on_converged is not None:
            on_converged(k)
    return iters, np.array(lams)


def _tangent():
    return np.array(ops.printA("-ret"))


def _reapplied_load_tangent_gap():
    """Relative change of the committed-state tangent across what the next analyze does
    before its predictor: re-apply the loads, update, formUnbalance."""
    k_last = _tangent()                               # last-iteration tangent (no update)
    ops.ladrunoTrialResidualNorm(ops.getTime() + 1.0)  # applyLoad + Domain::update + formUnbalance
    k_commit = _tangent()
    return np.linalg.norm(k_commit - k_last) / np.linalg.norm(k_last)


def test_lfbc2d_newton_converges_quadratically_like_forceBeamColumn():
    """No element-tolerance floor: every step converges and tracks vanilla."""
    _cantilever_2d("forceBeamColumn")
    it_van, lam_van = _push(2, 2, 2.0, 15)
    _cantilever_2d("LadrunoForceBeamColumn")
    it_lad, lam_lad = _push(2, 2, 2.0, 15)
    assert max(it_lad) <= max(it_van) + 1, (it_lad, it_van)
    np.testing.assert_allclose(lam_lad, lam_van, rtol=1.0e-9)


@pytest.mark.parametrize("model,node,dof,nsteps", [
    ("portal_2d", 3, 1, 40),
    ("cantilever_3d", 2, 2, 20),
])
def test_lfbc_committed_tangent_equals_last_iteration_tangent(model, node, dof, nsteps):
    """On every converged step, the tangent at the committed state after the loads are
    re-applied is the last-iteration tangent, bit for bit, and the pushover tracks
    vanilla forceBeamColumn."""
    build = {"portal_2d": _portal_2d, "cantilever_3d": _cantilever_3d}[model]
    build("forceBeamColumn")
    _, lam_van = _push(node, dof, 2.0, nsteps)

    gaps = []
    build("LadrunoForceBeamColumn")
    _, lam_lad = _push(node, dof, 2.0, nsteps, on_converged=lambda k: gaps.append(_reapplied_load_tangent_gap()))
    assert max(gaps) == 0.0, ["%.2e" % g for g in gaps]
    np.testing.assert_allclose(lam_lad, lam_van, rtol=1.0e-9)


def _failed_solve_then_revert():
    """A solve that fails by construction; analyze then reverts the domain to the last commit."""
    ops.test("NormUnbalance", 1.0e-99, 1, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) != 0


def _push_history(node, dof, du, nsteps, revert_first):
    iters, lams = [], []
    for k in range(nsteps):
        if revert_first:
            _failed_solve_then_revert()
        ops.test("NormDispIncr", 1.0e-8, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", node, dof, du)
        if ops.analyze(1) != 0:
            pytest.fail(f"Newton failed at displacement step {k + 1} (revert_first={revert_first})")
        iters.append(ops.testIter())
        lams.append(ops.getLoadFactor(2))
    return iters, np.array(lams)


@pytest.mark.parametrize("model,node,dof,nsteps", [
    ("portal_2d", 3, 1, 40),
    ("cantilever_3d", 2, 2, 20),
])
def test_lfbc_revert_to_last_commit_leaves_newton_history_unchanged(model, node, dof, nsteps):
    """revertToLastCommit restores the committed state exactly: a failed solve + revert
    before every step must give the same Newton history as no revert at all."""
    build = {"portal_2d": _portal_2d, "cantilever_3d": _cantilever_3d}[model]
    build("LadrunoForceBeamColumn")
    it_ref, lam_ref = _push_history(node, dof, 2.0, nsteps, revert_first=False)
    build("LadrunoForceBeamColumn")
    it_rev, lam_rev = _push_history(node, dof, 2.0, nsteps, revert_first=True)
    assert it_rev == it_ref
    np.testing.assert_allclose(lam_rev, lam_ref, rtol=1.0e-12)


def _elastic_cantilever(ndm, hinge):
    """Elastic cantilever (L = 2 m), optionally with a cohesive -hingeI at node i (01b G-3 model)."""
    L, E, A, Iz, Mc, Gf, Kpen = 2000.0, 30000.0, 60000.0, 4.5e8, 60.0e6, 1.2e6, 1.0e14
    ops.wipe()
    if ndm == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iz, 12500.0, 2.0 * Iz)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, Mc, Gf, "-linear", "-penalty", Kpen)
    ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, *(("-hingeI", 10) if hinge else ()))
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, *([0.0, 1.0, 0.0] if ndm == 2 else [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]))
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.analysis("Static")


def _reaction_trajectory(ndm, hinge, excursion):
    _elastic_cantilever(ndm, hinge)
    ndf = 3 if ndm == 2 else 6
    out = []
    for i in range(20):
        if excursion:
            u = list(ops.nodeDisp(2))
            u[1] += (-1.0) ** i * 5.0
            ops.ladrunoSetNodeTrial(2, *u, *([0.0] * ndf), *([0.0] * ndf))
            ops.ladrunoTrialResidualNorm()
            _failed_solve_then_revert()
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 2, 2, 0.5)
        assert ops.analyze(1) == 0
        ops.reactions()
        out.append(-ops.nodeReaction(1, 2))
    return out


@pytest.mark.parametrize("ndm", [2, 3])
@pytest.mark.parametrize("hinge", [False, True])
def test_lfbc_trial_excursion_and_revert_is_bit_exact(ndm, hinge):
    """A trial excursion followed by revertToLastCommit before every step leaves the
    trajectory bit-identical (01b G-3). The committed state must restore the deformation
    the trial state was actually solved at (vTrialBulkCommit): when the last Newton
    correction is below DBL_EPSILON, update() returns early and that deformation differs
    from the committed basic deformation by round-off."""
    clean = _reaction_trajectory(ndm, hinge, excursion=False)
    probed = _reaction_trajectory(ndm, hinge, excursion=True)
    assert probed == clean
