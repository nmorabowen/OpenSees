"""ADR-57 E6 — edge-edge one-scalar commit-cycle ALM (augmented Lagrange).

The EDGE_EDGE mode gains the OPTIONAL per-pair normal multiplier λ_N — the point-like analogue of the
mortar C2.2 Uzawa-on-commit. The augmented pressure is p = min(0, λ_N + εN·gN), traction tN = −p; one
Uzawa step per Domain::commit() drives λ_N ← min(0, λ_N + εN·gN_committed), so a HELD load augmented
to convergence reaches an εN-INDEPENDENT penetration at FINITE εN (penalty alone leaves O(P/εN)).
Point-like ⇒ ONE scalar per edge pair, NO shared-node accumulator (the whole simplification vs mortar).
The held-load `analyze_augmented` proc (ADR-41 D1) works UNCHANGED with query=ops.ladrunoEdgePenetration.
Implicit-only (the mass-only explicit tangent degenerates Uzawa to single-pass penalty — the disclosed
ADR-47 limitation, inherited). OFF by default (λ_N≡0 ⇒ the E2 penalty path) ⇒ byte-identical. Mechanics
validated oracle-first in contact_prototypes/proto_e6_alm.py.

Two crossed bars (the E2 geometry): a slave bottom edge ∥x crossing perpendicular over a fixed master
top edge ∥y, n=ẑ. Pressed down, held by the edge-edge contact; soft z-springs regularize the point-pair
differential mode (as in E2/E3).

Gates:
  - WITHIN-STEP CONVERGENCE — held-load augmentation drives ladrunoEdgePenetration < augTol, εN-INDEPENDENT
    at finite εN (penalty-only leaves O(P/εN), εN-DEPENDENT);
  - RELEASE — reversing the load lifts the bar off: the contact OPENS (penetration → 0, no locked-in pull);
  - OPT-IN BYTE-IDENTITY — the first solve (λ_N still 0) is BITWISE identical with -edgeAlm ON vs OFF (the
    ALM adds nothing until augmentation — the penalty path).
"""
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

KSPRING = 500.0      # soft regularizing spring (anchored SEPARATED so release fully opens the contact)
PTOT = 400.0         # total downward load on the slave bottom edge
ZC = -1.0e-4         # slave bottom edge seeded slightly penetrating (contact active from step 1)
ZSPRING = 0.01       # spring anchor z (ABOVE the master ⇒ removing the load opens the contact)


def _build(epsN, alm, mu=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.0, ZC + 0.3), (8, -1.0, 0.0, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)               # x,y fixed; z free
    ops.uniaxialMaterial("Elastic", 1, KSPRING)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, ZSPRING)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    args = [1, 1, 2, "-mortar", "-epsN", epsN, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if alm:
        args += ["-edgeAlm", "-edgeAugTol", 1.0e-9]
    ops.contact(*args)


def _load_down():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -PTOT / 2.0)


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def test_e6_augment_converges_eps_independent():
    """Held-load augmentation drives the edge penetration → augTol at FINITE εN, εN-INDEPENDENTLY —
    while the penalty-only penetration is O(P/εN) (spans orders of magnitude across εN)."""
    AUGTOL = 1.0e-9
    pen_penalty, pen_alm = {}, {}
    for epsN in (1.0e5, 1.0e7):
        # -edgeAlm on (so ladrunoEdgePenetration is live). On the FIRST solve λ_N is still 0, so the
        # reported penetration is the PURE PENALTY penetration O(P/εN) — εN-DEPENDENT. analyze_augmented
        # then drives it to an εN-INDEPENDENT < augTol via the held-load Uzawa.
        _build(epsN, alm=True)
        _load_down()
        _static()
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, "ALM edge-edge push did not converge"
        pen_penalty[epsN] = ops.ladrunoEdgePenetration()    # λ_N=0 on this first solve ⇒ penalty pen
        status, n_aug, pen = analyze_augmented(ops, maxAug=40, augTol=AUGTOL,
                                               query=ops.ladrunoEdgePenetration)
        assert status == 0, f"εN={epsN:.0e}: augmentation did not converge (status {status}, pen {pen})"
        pen_alm[epsN] = pen
        assert pen < AUGTOL, f"εN={epsN:.0e}: ALM penetration {pen} not < augTol {AUGTOL}"
    # the penalty penetration is εN-DEPENDENT (shrinks ~1/εN ⇒ spans ≥10× across the two εN)
    spread = max(pen_penalty.values()) / max(min(pen_penalty.values()), 1e-300)
    assert spread > 10.0, f"penalty penetration not εN-dependent (spread {spread}): {pen_penalty}"
    # the ALM penetration is εN-INDEPENDENT (both < augTol)
    assert max(pen_alm.values()) < AUGTOL, f"ALM penetration not εN-independent: {pen_alm}"


def test_e6_release_opens():
    """RELEASE: augment under the down-load (λ_N carries the reaction), then REVERSE the load to lift
    the bar off. The contact OPENS — penetration → 0, no locked-in λ_N pull keeping it engaged."""
    epsN = 1.0e6
    _build(epsN, alm=True)
    _load_down()
    _static()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    status, _, pen = analyze_augmented(ops, maxAug=40, augTol=1e-9, query=ops.ladrunoEdgePenetration)
    assert status == 0 and pen < 1e-9, f"pre-release augmentation failed (status {status}, pen {pen})"
    # REVERSE: a second pattern applies an equal-and-opposite UP load + extra ⇒ net upward ⇒ the bar
    # lifts to the SEPARATED spring rest (z≈ZSPRING>0). If λ_N were locked negative, p=min(0,λ_N+εN·gN)
    # would keep the contact active and pull the bar back; a clean open proves the KKT release.
    ops.pattern("Plain", 2, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, +PTOT)              # cancel the −P/2 each + push up
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "release solve did not converge"
    # drive the released state to convergence (the Uzawa clamps λ_N → 0 as gN goes positive)
    analyze_augmented(ops, maxAug=20, augTol=1e-9, query=ops.ladrunoEdgePenetration)
    z5 = ZC + ops.nodeDisp(5, 3)
    assert z5 > 0.0, f"bar did not lift off the master (z5={z5}); contact stayed locked?"
    assert ops.ladrunoEdgePenetration() < 1e-9, "edge-edge still reports active penetration after release"


def test_e6_optin_byte_identical():
    """OPT-IN BYTE-IDENTITY: on the FIRST solve λ_N is still 0, so -edgeAlm ON is BITWISE identical to
    OFF (the ALM contributes nothing until the commit-cycle augmentation runs — it IS the penalty path)."""
    def _first_step(alm):
        _build(1.0e6, alm=alm)
        _load_down()
        _static()
        assert ops.analyze(1) == 0
        return tuple(ops.nodeDisp(t, 3) for t in (5, 6))
    assert _first_step(True) == _first_step(False)
