"""ADR-78 removal lane — a contact survives the analysis removing its nodes.

Before this, `handle()` re-ran on every `domainChanged()` and any node missing
from a contact surface was FATAL (ADR-78 P1), because a deck typo and a node the
analysis had removed were indistinguishable at that point. That made
`recorder Collapse` + a declared contact abort mid-run, which blocked the fork's
own progressive-collapse roadmap (ADR-51/54).

The fix splits the two cases by WHEN the node went missing:

  * absent at `contactSurface` declaration  -> deck typo   -> refuse, at the deck line
  * present then, absent at handle()        -> the analysis removed it -> prune, continue

so the abort is preserved exactly where it was right, and dropped exactly where
it was wrong. The `-kn auto` / `-soft` / `kn <= 0` aborts are untouched — apeGmsh
ADR 0092 INV-1 leans on those, not on this path.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _master_block():
    """The fixed bottom brick. Its top face (5,6,7,8) is the contact master."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)
    for t, (x, y, z) in zip(
            (1, 2, 3, 4, 5, 6, 7, 8),
            [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
             (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]):
        ops.node(t, float(x), float(y), float(z))
    ops.element("stdBrick", 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    for n in (5, 6, 7, 8):
        ops.fix(n, 1, 1, 0)
        ops.mass(n, 1.0, 1.0, 1.0)   # explicit lane needs mass on every free DOF


def _slaves(grounded, spare):
    """Slave node-set just above the master face, slightly penetrating.

    `grounded` get a soft z spring to an anchor so the static system is never
    singular once contact separates (the ADR-76 idiom the P4 soft tests use).
    `spare` get NO element and NO load, which is what makes them removable:
    removing a node that still carries a NodalLoad leaves that load pointing at
    freed memory and OpenSees ACCESS-VIOLATES on the next analyze -- measured,
    and it does so with `constraints Plain` and no contact at all, so it is a
    stock lifecycle hazard rather than anything contact does.
    """
    xy = {11: (0, 0), 12: (1, 0), 13: (1, 1), 14: (0, 1)}
    anchor = 900
    for t in sorted(xy):
        ops.node(t, float(xy[t][0]), float(xy[t][1]), 0.995)
        # A `spare` node is fixed in z as well. Its z DOF would otherwise be held
        # up by the contact penalty ALONE, so the tangent goes singular the moment
        # contact separates -- measured as `U(i,i) = 0`. Fixing it keeps the node a
        # full member of the contact surface (which is all this test needs of it)
        # while leaving it removable: no element, no load, no free DOF.
        ops.fix(t, 1, 1, 0 if t in grounded else 1)
        ops.mass(t, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 90, 1.0e4)
    for t in grounded:
        ops.node(anchor + t, float(xy[t][0]), float(xy[t][1]), 0.995)
        ops.fix(anchor + t, 1, 1, 1)
        ops.element("zeroLength", anchor + t, anchor + t, t, "-mat", 90, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave", *sorted(xy))
    ops.contact(1, 1, 2, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    return list(grounded), list(spare)


def _analysis():
    """Explicit transient, deliberately.

    The property under test is that handle() no longer ABORTS when the analysis
    removes a surface node -- not how a toy penalty deck converges. An implicit
    Newton run here chatters across the contact activation threshold and fails on
    the convergence test, which would make this suite report a solver artefact as
    a removal-lane defect. Central difference takes the iteration out entirely,
    and it is the same idiom the P3 missing-node test already uses.
    """
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _step():
    return ops.analyze(1, 1.0e-5)


# The declaration-time refusal that makes this whole lane sound lives in the
# validation battery where the other refusals are:
#   test_contact_review_p3_validation.py::test_missing_node_refused_at_declaration
# Kept there rather than duplicated here.


def test_removed_slave_node_is_pruned_and_the_run_continues(capfd):
    """The analysis removes ONE slave node mid-run. Pre-fix this was a mid-run
    FATAL (`node 14 ... does not exist; ABORTING`); now it is pruned and the next
    step still solves."""
    _master_block()
    grounded, spare = _slaves(grounded=(11, 12, 13), spare=(14,))
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in grounded:
        ops.load(n, 0.0, 0.0, -10.0)
    _analysis()
    assert _step() == 0, "baseline (nothing removed yet) must solve"

    ops.remove("node", 14)                  # the analysis retires a node
    assert _step() == 0, "removed slave node must prune and continue"
    # Assert the lane actually RAN. Without this the test would also pass if the
    # contact had quietly stopped doing anything, which is the failure mode this
    # whole ADR exists to forbid -- and the one an earlier ADR-41 deck fell into.
    out = capfd.readouterr()
    blob = out.out + out.err
    assert "removal lane" in blob and "contactSurface 2" in blob, blob[-2000:]


def test_contact_retires_when_its_surface_is_emptied(capfd):
    """Remove EVERY slave node. The interaction has nothing left to act on, so it
    is retired with a named notice and the run continues -- rather than aborting,
    and rather than silently pretending a contact is still active."""
    _master_block()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):                  # load the MASTER, so no slave is loaded
        ops.load(n, 0.0, 0.0, -10.0)
    _slaves(grounded=(), spare=(11, 12, 13, 14))
    _analysis()
    assert _step() == 0

    for n in (11, 12, 13, 14):
        ops.remove("node", n)
    assert _step() == 0, "emptied surface must retire, not abort"
    out = capfd.readouterr()
    blob = out.out + out.err
    assert "RETIRING contact 1" in blob, blob[-2000:]


def test_kn_auto_failure_is_still_fatal():
    """Guard on the guard. ADR 0092 INV-1 picks the owner rank by a master-node
    majority PROXY, and ADR-78 states that proxy is safe ONLY because `-kn auto`
    failure is fatal. The removal lane must not have softened it: a master face
    with no owning solid still ABORTS."""
    L, P, KSPRING = 1.0, 1.0e2, 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    top = [1, 2, 3, 4]
    for t, (x, y, z) in zip(top, [(0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]):
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 1, 1, 1)
    ops.node(99, L / 2, L / 2, L - 1.0e-8)
    ops.fix(99, 1, 1, 0)
    ops.node(98, L / 2, L / 2, L - 1.0e-8)
    ops.fix(98, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KSPRING)
    ops.element("zeroLength", 1, 98, 99, "-mat", 1, "-dir", 3)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, "auto", "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    _analysis()
    assert _step() == -1, "-kn auto with no owning solid must still ABORT"
