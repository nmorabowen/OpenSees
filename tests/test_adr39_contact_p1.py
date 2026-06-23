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


# ---------------------------------------------------------------------------
# P1b — LadrunoContactDomain attachment + lifecycle hooks (still ZERO force)
# ---------------------------------------------------------------------------
def _run_with_contact(integrator, nsteps=40):
    """Same truss model, but with a LadrunoContactDomain attached (one contact
    definition over zero-force empty-connectivity adapters)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n = 5
    for i in range(n):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, n + 1):
        ops.fix(i, 0, 1, 1)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, n):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    ops.setNodeVel(n, 1, 0.5, "-commit")

    # --- define contact (topology only; P1b is zero-force) ---
    ops.contactSurface(1, "-master", 3, 1, 2, 3)   # master segment surface
    ops.contactSurface(2, "-slave", 4, 5)           # slave node-set
    ops.contact(1, 1, 2, 1.0e6, 5.0e5, 0.3)         # one contact -> one adapter

    ops.constraints("LadrunoContact")
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
    for _ in range(nsteps):
        assert ops.analyze(1, 1.0e-3) == 0
        out.append(tuple(ops.nodeDisp(i + 1, 1) for i in range(n)))
    return out


@pytest.mark.parametrize("integrator", ["CentralDifferenceLadruno", "Newmark"])
def test_p1b_contact_defined_still_bitwise_and_commits(integrator):
    """A defined (zero-force) contact pair: still bitwise-identical to no-contact
    (empty-connectivity adapters), and Domain::commit() drives ContactDomain.commit()
    every step (the design-gate BLOCKER-1 lifecycle hook actually fires)."""
    ref = _run("Plain", integrator)
    con = _run_with_contact(integrator)
    assert con == ref, f"{integrator}: defined zero-force contact changed the result"
    info = ops.ladrunoContactInfo()              # [numContacts, numCommits, numReverts]
    assert info[0] == 1, f"expected 1 contact, got {info}"
    assert info[1] >= 40, f"commit hook should fire each step, got {info}"


def test_p1b_wipe_clears_contact_engine():
    """ops.wipe() (Domain::clearAll) must delete the contact engine — a fresh
    model must not inherit the old one (the ADR-30 theEQs-leak class of bug)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, 2.0, 0.0, 0.0)
    ops.contactSurface(1, "-master", 3, 1, 2, 3)
    ops.contactSurface(2, "-slave", 3)
    ops.contact(1, 1, 2, 1.0e6, 0.0, 0.0)
    assert ops.ladrunoContactInfo()[0] == 1
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)      # fresh model, no contact
    # [numContacts, numCommits, numReverts, numMortarContacts] (4th added ADR-41 C2)
    assert ops.ladrunoContactInfo() == [0, 0, 0, 0], "contact engine leaked across wipe"


def test_p1b_null_contact_is_pure_plain():
    """constraints LadrunoContact with NO contact defined injects ZERO adapters
    -> identical to constraints Plain (the nullable-engine guard)."""
    ref = _run("Plain", "Newmark")
    con = _run("LadrunoContact", "Newmark")        # no contactSurface/contact calls
    assert con == ref
    assert ops.ladrunoContactInfo() == [0, 0, 0, 0]


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
