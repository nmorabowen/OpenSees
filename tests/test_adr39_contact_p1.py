"""ADR-39 P1a — LadrunoContactHandler injection is graph-neutral (bitwise gate).

P1a injects ONE empty-connectivity zero LadrunoContactFE adapter via the custom
LadrunoContactHandler. The gate verdict: with EMPTY connectivity the adapter adds
no DOF-graph edges, so the numberer permutation is untouched and the result is
**bitwise** identical to no-contact — under BOTH an explicit (CentralDifference-
Ladruno) and an implicit (Newmark) integrator. This is the clean true/false
plumbing gate (P2 switches to declared connectivity + a 1e-12 tolerance).

MP constraints (equalDOF) are NOT enforced by the P1a handler (Plain-style) — that
delegation lands in P1b — so this battery uses SP (fix) constraints only.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _run(handler, integrator, nsteps=40):
    """A small 3D truss chain, tip given an initial velocity, run dynamically."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n = 5
    for i in range(n):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, n + 1):
        ops.fix(i, 0, 1, 1)            # free along x only (SP constraints only)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, n):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    ops.setNodeVel(n, 1, 0.5, "-commit")

    ops.constraints(handler)
    ops.numberer("Plain")
    if integrator == "Newmark":
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1e-12, 25, 0)
        ops.algorithm("Newton")
    else:
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
    ops.analysis("Transient")

    out = []
    dt = 1.0e-3
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        out.append(tuple(ops.nodeDisp(i + 1, 1) for i in range(n)))
    return out


@pytest.mark.parametrize("integrator", ["CentralDifferenceLadruno", "Newmark"])
def test_p1a_contact_handler_bitwise_identical(integrator):
    """constraints LadrunoContact (empty zero adapter) == constraints Plain,
    BITWISE, under explicit and implicit."""
    ref = _run("Plain", integrator)
    con = _run("LadrunoContact", integrator)
    assert con == ref, f"{integrator}: LadrunoContact deviated from Plain (should be graph-neutral)"


def test_p1a_handler_runs_and_rebuilds():
    """Smoke: the handler runs, and a forced domainChanged (remove an element
    mid-run) leaves the analysis healthy — exercises adapter destroy/rebuild."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(4):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, 5):
        ops.fix(i, 0, 1, 1)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, 4):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(10, 1.0e-3) == 0
    ops.remove("element", 3)          # forced domainChanged -> handler re-runs
    assert ops.analyze(10, 1.0e-3) == 0
