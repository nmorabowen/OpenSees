"""Contact review fix P4 (2026-07) — EQ_Constraint parity with upstream PlainHandler.

The contact handler REPLICATES PlainHandler's DOF/FE loop (so the contact-FE start tag
is knowable), but upstream PlainHandler gained an EQ_Constraint block AFTER the replica
was written: enforce a trivial-identity EQ (DOF mark -4) and loudly warn-and-ignore a
non-trivial one. The replica had NEITHER — every `equationConstraint`, including the
fork's own LadrunoTie (ADR-62) rows, was SILENTLY dropped under `constraints
LadrunoContact`: parts of the structure ran unconnected with no diagnostic (review
HIGH-3, the exact silent-zero class the handler's non-homogeneous-SP warning guards).

Post-fix: the non-trivial EQ is still not enforced (actual EQ enforcement is the
projection handler's job — contact + non-trivial ties in one analysis stays
unsupported) but the handler now SAYS SO loudly. The warning assert below FAILS on the
pre-fix build (silence).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _model_with_contact_and_eq():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # a minimal live contact: slave node over a rigid plane
    ops.node(1, 0.0, 0.0, 1.0e-4)
    ops.mass(1, 1.0, 1.0, 1.0)
    # two extra free nodes tied by a NON-TRIVIAL equationConstraint u2x = 0.5*u3x
    ops.node(2, 5.0, 0.0, 0.0); ops.mass(2, 1.0, 1.0, 1.0)
    ops.node(3, 6.0, 0.0, 0.0); ops.mass(3, 1.0, 1.0, 1.0)
    ops.equationConstraint(2, 1, 1.0, 3, 1, -0.5)
    ops.contactSurface(20, "-slave", 1)
    ops.contactPlane(1, 20, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e6)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -100.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")


def test_nontrivial_eq_warns_loudly_under_contact_handler(capfd):
    """A non-trivial equationConstraint under `constraints LadrunoContact` must produce
    the loud NOT-ENFORCED warning (pre-fix: dead silence) and the analysis still runs."""
    _model_with_contact_and_eq()
    assert ops.analyze(1, 1.0e-4) == 0
    err = capfd.readouterr().err
    assert "EQ_Constraint at node 2" in err and "NOT ENFORCED" in err, (
        "the contact handler silently dropped a non-trivial equationConstraint "
        f"(no warning in stderr; got: {err[-300:]!r})")


def test_eq_is_not_enforced_but_analysis_is_healthy():
    """Documented behavior: the tie is warned-and-IGNORED (node 2 does not follow
    0.5*u3 — both drift freely), and the contact itself still works."""
    _model_with_contact_and_eq()
    ops.load(3, 50.0, 0.0, 0.0)              # drive node 3 laterally
    for _ in range(200):
        assert ops.analyze(1, 1.0e-3) == 0
    u2, u3 = ops.nodeDisp(2, 1), ops.nodeDisp(3, 1)
    assert abs(u3) > 1.0e-6, "driver node never moved (rig broken)"
    assert abs(u2 - 0.5 * u3) > 1.0e-9, (
        "u2 == 0.5*u3: the EQ appears ENFORCED — the warn-and-ignore contract changed")
    # the contact slave still settles at ~F/kn (the handler kept doing its real job)
    z = 1.0e-4 + ops.nodeDisp(1, 3)
    assert abs(z + 100.0 / 1.0e6) < 5.0e-4
