"""ADR-74 review follow-up (F1) — positive coverage for the PlainHandler EQ-sweep.

#598 rewrote PlainHandler's per-node `getEQs()` sweep into an `allEQs` multimap
(the same while->equal_range transform proven byte-identical for the MP sweep by
the tie gate). But the adversarial review (2026-07-23) found that gate has ZERO
EQ coverage: `equalDOF` makes MP_Constraints, `EQ_Constraint` is produced only by
LadrunoTie, and no gate deck used it — so the rewritten EQ loop NEVER EXECUTED
under test. A regression in its `equal_range` bounds or body would have sailed
through 18/18.

This gate makes the loop execute. The witness is the per-node "constraint matrix
not identity" warning: it is emitted from INSIDE the rewritten loop, for the
specific constrained node, so its presence proves `equal_range(nodeID)` reached
that node and ran the body. A loop that skipped the node (wrong bounds) would go
silent -> test fails. Runs under both `numberer Plain` and `numberer RCM` (the
two -4-fixup variants), and checks the free-DOF numbering is a valid bijection.

`equationConstraint cNode cDOF cCoef rNode rDOF rCoef` (1-based dof): the handler
computes Ccr = -rCoef/cCoef; Ccr==[1.0] is the trivial-identity branch (marks
the dof -4), anything else is warn-and-ignore.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _build(numberer):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, 6):
        ops.node(i, float(i - 1), 0.0)
        ops.mass(i, 1.0, 1.0)
    ops.fix(1, 1, 1)
    for e in range(2, 6):
        ops.element("truss", e, e - 1, e, 1.0, 1)
    # NON-trivial EQ on node 3 and node 4 (Ccr = 0.5 != 1.0 -> warn+ignore each).
    # Two distinct constrained nodes so equal_range must resolve more than one key.
    ops.equationConstraint(3, 1, 1.0, 2, 1, -0.5)
    ops.equationConstraint(4, 1, 1.0, 2, 1, -0.5)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 1.0, 0.0)
    ops.constraints("Plain")
    ops.numberer(numberer)
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-8, 20, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


@pytest.mark.parametrize("numberer", ["Plain", "RCM"])
def test_plainhandler_eq_sweep_executes_and_warns(numberer, capfd):
    _build(numberer)
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    assert rc == 0, f"analyze failed (rc={rc}) under numberer {numberer}"
    # the rewritten EQ loop must reach BOTH constrained nodes and run the
    # non-identity warn branch for each
    for node in (3, 4):
        assert "not identity" in err and f"node {node}" in err, (
            f"PlainHandler EQ-sweep did not warn for the non-trivial EQ on node "
            f"{node} under numberer {numberer} — the rewritten equal_range loop "
            f"may be skipping nodes. stderr tail={err[-500:]!r}")
    # free-DOF numbering is a valid dense bijection 0..K-1 (no leaked -4/-3/-2)
    eqns = []
    for n in range(1, 6):
        eqns += [d for d in ops.nodeDOFs(n) if d >= 0]
    assert len(eqns) == len(set(eqns)), f"duplicate equation numbers: {sorted(eqns)}"
    assert sorted(eqns) == list(range(len(eqns))), (
        f"free-DOF numbering is not a dense bijection 0..K-1: {sorted(eqns)}")


def test_trivial_identity_eq_marks_and_stays_healthy(capfd):
    """Exercise the OTHER branch of the loop: a trivial-identity EQ (Ccr==[1.0])
    takes the setID(-4) path instead of the warn path. Assert the run stays
    healthy (this documents the stock PlainHandler+EQ behavior the rewrite
    preserves; the -4 resolution itself is a pre-existing upstream limitation)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, 5):
        ops.node(i, float(i - 1), 0.0)
        ops.mass(i, 1.0, 1.0)
    ops.fix(1, 1, 1)
    for e in range(2, 5):
        ops.element("truss", e, e - 1, e, 1.0, 1)
    ops.equationConstraint(3, 1, 1.0, 2, 1, -1.0)   # Ccr = 1.0 -> trivial-identity branch
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 1.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-8, 20, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    # the trivial branch must NOT emit the non-identity warning
    assert "not identity" not in err, (
        f"trivial-identity EQ wrongly took the non-identity warn branch: {err[-400:]!r}")
    assert rc == 0, f"trivial-identity EQ under constraints Plain broke the analysis (rc={rc})"
