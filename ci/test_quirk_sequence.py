"""Self-test for quirk-lint rule `sequence` (alias L7, WP-140) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_sequence.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _run(tmp_path, body, stamped=False, where="SRC/element/Element.cpp"):
    src = (kit.STAMP if stamped else "") + body
    root = kit.tree(tmp_path, {where: src})
    used = set()
    out = cq.check_sequence(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def _resp(body):
    return "int Element::getResponse(int responseID, Information &eleInfo)\n{\n" + body + "\n}\n"


def test_flags_the_c15_incident_in_a_vanilla_file(tmp_path):
    """The exact pre-a2004e0a7 line of vanilla Element.cpp (WP-124 C15)."""
    out = _run(tmp_path, _resp(
        "  switch (responseID) {\n  case 444444:\n"
        "    return eleInfo.setVector(this->getResistingForceIncInertia()-this->getRayleighDampingForces()"
        "-this->getResistingForce());\n  default:\n    return -1;\n  }"))
    assert len(out) == 1 and out[0].startswith("L7 SRC/element/Element.cpp:4:")
    assert "getResistingForceIncInertia(), getRayleighDampingForces(), getResistingForce()" in out[0]


def test_passes_the_c15_fix(tmp_path):
    assert _run(tmp_path, _resp(
        "  switch (responseID) {\n  case 444444: {\n"
        "    Vector inertial(this->getResistingForceIncInertia());\n"
        "    inertial -= this->getRayleighDampingForces();\n"
        "    inertial -= this->getResistingForce();\n"
        "    return eleInfo.setVector(inertial);\n  }\n  default:\n    return -1;\n  }")) == []


def test_flags_vector_plus_matrix_times_vector(tmp_path):
    """Cross-type too: a getMass() that forms the mass can refill the residual storage as a
    side effect (LadrunoBrick's formInertiaTerms(1) writes resid)."""
    out = _run(tmp_path, "const Vector &E::getResistingForceIncInertia(void)\n{\n"
                        "  res = this->getResistingForce() + this->getMass() * accel;\n  return res;\n}\n")
    assert len(out) == 1 and "getResistingForce(), getMass()" in out[0]


@pytest.mark.parametrize("body", [
    "  theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);",                 # one call, args
    "  foo(this->getMass(), this->getTangentStiff());",                               # separate ARGUMENTS
    "  theVector->addMatrixVector(0.0, this->getMass(), vel, alphaM);",
    "  res = this->getMass() * accel;",                                             # one call, arithmetic
    "  res = this->getResistingForce();\n  res += this->getRayleighDampingForces();", # separate statements
    "  // res = this->getResistingForce() - this->getRayleighDampingForces();\n  x = 1;",
    "  opserr << \"getResistingForce() - getMass()\" << endln;",
    "  K = theEle->getTangentStiff();\n  M = theEle->getMass();\n  A = K - M;",      # owned copies
])
def test_passes_non_arithmetic_or_sequenced_uses(tmp_path, body):
    assert _run(tmp_path, "void E::f(void)\n{\n" + body + "\n}\n") == []


def test_flags_a_pointer_receiver_and_an_other_file(tmp_path):
    out = _run(tmp_path, "void FE::g(void)\n{\n  r = myEle->getResistingForceIncInertia() - "
                        "myEle->getResistingForce();\n}\n", where="SRC/analysis/fe_ele/FE.cpp")
    assert len(out) == 1 and out[0].startswith("L7 SRC/analysis/fe_ele/FE.cpp:3:")


def test_waiver(tmp_path):
    ok = ("void E::f(void)\n{\n  // ladruno-lint: sequence-ok getMass never writes resid in this element\n"
          "  res = this->getResistingForce() + this->getMass() * a;\n}\n")
    assert _run(tmp_path, ok, stamped=True) == []
    short = ("void E::f(void)\n{\n  res = this->getResistingForce() + this->getMass() * a;"
             "   // ladruno-lint: sequence-ok ok\n}\n")
    out = _run(tmp_path, short, stamped=True)
    assert len(out) == 1 and "reason too short" in out[0]


def test_stale_waiver(tmp_path):
    stale = ("void E::f(void)\n{\n  // ladruno-lint: sequence-ok this used to combine two accessors\n"
             "  res = this->getResistingForce();\n}\n")
    out = _run(tmp_path, stale, stamped=True)
    assert len(out) == 1 and out[0].startswith("W ") and "stale sequence-ok" in out[0]
