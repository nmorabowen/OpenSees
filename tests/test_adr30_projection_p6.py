"""ADR-30 Phase 6 — two robustness guards for LadrunoProjection (closes ADR O1).

3a  Near-singular L^T M L condition-number gate (projector). buildMass() used to probe
    singularity with a single exact-pivot test solve, which a near-singular (ill-conditioned)
    L^T M L passes — then project() amplifies round-off into a garbage acceleration. P6 estimates
    cond = lambda_max/lambda_min via a self-contained symmetric Jacobi eigensolve and REFUSES
    above 1e12 (also catching the exact-singular case), WARNS above 1e8.

3b  Frozen-Ccr staleness guard (handler). A rigidLink-beam / rigidDiaphragm bakes the master's
    lever arm into L from Ccr once (small-rotation). For accumulated master rotation >~0.1 rad the
    tie silently drifts. P6 flags the master-rotation -> slave-translation couplings at handle()
    and warns ONCE in applyLoad() when the master's rotation drift crosses 0.1 rad.

  TP6-1  well-conditioned 2-master EQ group (equal masses) runs (gate does not false-refuse).
  TP6-2  ill-conditioned 2-master EQ group (mass ratio 1e13, cond ~1e13) -> REFUSED at analyze.
  TP6-3  moderately ill-conditioned (ratio 1e10, cond in (1e8,1e12)) -> WARNS but RUNS.
  TP6-4  staleness: a rigidLink-beam master rotated past 0.1 rad emits the one-time NOTE.
  TP6-5  staleness: the same model kept at small rotation does NOT emit the NOTE.
"""
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

PROJ = "LadrunoProjection"
CDL = "CentralDifferenceLadruno"


# --------------------------------------------------------------------------
# 3a — condition-number gate. A 2-master equationConstraint puts BOTH retained
# DOFs in one group (nRet=2), so L^T M L is 2x2 and its conditioning is driven by
# the retained-mass ratio: cond ~ m1/m2. equal -> ~1, 1e13 -> refused, 1e10 -> warn.
# --------------------------------------------------------------------------
def _run_eq_group(m1, m2, m3, nsteps=3, dt=1.0e-3):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.node(3, 2.0); ops.mass(3, m3)
    # u3 = 0.5 u1 + 0.5 u2   (nodes 1,2 retained -> nRet=2 in one group, node 3 slave)
    ops.equationConstraint(3, 1, 1.0, 1, 1, -0.5, 2, 1, -0.5)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    rc = 0
    for _ in range(nsteps):
        rc = ops.analyze(1, dt)
        if rc != 0:
            break
    return rc


def test_TP6_1_well_conditioned_group_runs():
    # equal masses -> cond ~ 1.5 -> the gate must NOT refuse a healthy group
    assert _run_eq_group(1.0, 1.0, 1.0) == 0


def test_TP6_2_near_singular_group_refused():
    # mass ratio 1e13 -> cond ~ 1e13 > 1e12 -> buildMass refuses -> analyze fails.
    # The ONLY difference from TP6-1 is the retained-mass ratio, so a failure here is the gate.
    assert _run_eq_group(1.0e13, 1.0, 1.0) != 0


def test_TP6_3_ill_conditioned_group_warns_but_runs():
    # ratio 1e10 -> cond ~ 1e10 in (1e8, 1e12) -> WARN but still RUN (no refusal)
    assert _run_eq_group(1.0e10, 1.0, 1.0) == 0


# --------------------------------------------------------------------------
# 3b — frozen-Ccr staleness warning. 2D rigidLink -beam: node 2 (slave) is tied to
# node 1 (master) with a lever arm, so the master rz (dof 2, rotational) couples into
# the slave ux,uy (dofs 0,1, translational) -> flagged. Drive the master rotation with
# an applied moment; the NOTE fires once the drift crosses 0.1 rad.
# --------------------------------------------------------------------------
def _run_rigidlink_beam(moment, nsteps, dt, logpath):
    ops.wipe()
    ops.logFile(logpath, "-noEcho")          # capture opserr (the NOTE) to a file
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)                     # offset slave -> nonzero lever arm
    ops.mass(1, 1.0, 1.0, 1.0)               # translational + rotational mass on the master
    ops.mass(2, 1.0, 1.0, 1.0)               # tied DOFs need lumped mass (projection rule)
    ops.rigidLink("beam", 1, 2)              # ties node2 (ux,uy,rz) to node1 with the lever arm
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, moment)            # moment on the master rz
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
    rz = abs(ops.nodeDisp(1, 3))
    # opserr's FileStream buffers; redirecting to a throwaway file CLOSES (flushes) logpath.
    ops.logFile(logpath + ".flush", "-noEcho")
    ops.wipe()
    txt = ""
    if os.path.exists(logpath):
        with open(logpath, "r", errors="ignore") as fh:
            txt = fh.read()
    return rz, txt


def test_TP6_4_staleness_note_fires_on_large_rotation(tmp_path):
    log = str(tmp_path / "stale_big.log")
    # large moment + many steps -> master rz blows past 0.1 rad
    rz, txt = _run_rigidlink_beam(50.0, 400, 5.0e-3, log)
    assert rz > 0.1, f"master did not rotate enough to trigger (rz={rz:.4f})"
    assert ("stale" in txt.lower()) or ("rotated" in txt.lower()), (
        f"frozen-Ccr staleness NOTE not emitted; log:\n{txt[-800:]}")


def test_TP6_5_no_staleness_note_on_small_rotation(tmp_path):
    log = str(tmp_path / "stale_small.log")
    # tiny moment + few steps -> rz stays well under 0.1 rad -> no NOTE
    rz, txt = _run_rigidlink_beam(1.0e-4, 20, 1.0e-3, log)
    assert rz < 0.1, f"control rotated too much (rz={rz:.4f})"
    assert "stale" not in txt.lower(), f"staleness NOTE fired on a small rotation; log:\n{txt}"
