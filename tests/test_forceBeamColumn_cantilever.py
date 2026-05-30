"""T1 worked example (feature class: derivative) — end-loaded 2D cantilever.

Demonstrates the canonical Zone-A pattern:
  * physics pulled from the SHARED params table (D1) the Tcl twin also reads,
  * closed-form asserts at isclose rel_tol 1e-9,
  * the universal equilibrium guard (G7),
  * the Axis-1 Python serialization proxies (D4): setResponse smoke + database
    round-trip.

Mirror: tests/tcl/forceBeamColumn_cantilever.tcl (same params file).
"""
from math import isclose

import pytest

from _testbed import ops
from _testbed.params import load_params
from _testbed.equilibrium import assert_equilibrium
from _testbed.roundtrip import setresponse_smoke, database_roundtrip

P = load_params("forceBeamColumn_cantilever")
L, E, A, I, Pload = P["L"], P["E"], P["A"], P["I"], P["P"]

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]


def _build():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.node(2, L, 0.0)
    ops.section("Elastic", 1, E, A, I)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.geomTransf("Linear", 1)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, Pload, 0.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)
    ops.reactions()


def test_tip_deflection():
    _build()
    assert isclose(ops.nodeDisp(2, 2), Pload * L**3 / (3 * E * I), rel_tol=1e-9)


def test_tip_rotation():
    _build()
    assert isclose(ops.nodeDisp(2, 3), Pload * L**2 / (2 * E * I), rel_tol=1e-9)


def test_base_reactions():
    _build()
    assert isclose(ops.nodeReaction(1, 2), -Pload, rel_tol=1e-9)
    assert isclose(ops.nodeReaction(1, 3), -Pload * L, rel_tol=1e-9)


def test_global_equilibrium():
    _build()
    # 2D beam -> 2 leading translational DOFs (Fx, Fy); the moment DOF is not
    # summed. Total applied force = P in global Y at the tip.
    assert_equilibrium([1, 2], n_force_dofs=2, applied_total=[0.0, Pload], atol=1e-6)


def test_setresponse_smoke():
    _build()
    setresponse_smoke(1, [("force",), ("section", 1, "force")])


def test_database_roundtrip():
    database_roundtrip(_build, probe_nodes=[2], ndf=3)
