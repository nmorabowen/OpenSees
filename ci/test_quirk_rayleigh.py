"""Self-test for quirk-lint rule `rayleigh` (alias L1, WP-115) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_rayleigh.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _run(tmp_path, body, stamped=True):
    src = (kit.STAMP if stamped else "") + body
    root = kit.tree(tmp_path, {"SRC/element/Elem.cpp": src})
    used = set()
    out = cq.check_rayleigh(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def _inc(body):
    return "const Vector &Elem::getResistingForceIncInertia(void)\n{\n" + body + "\n}\n"



GOOD = {
    "snapshot_static_local": _inc(
        "  static Vector res(6);\n  res = this->getResistingForce();\n"
        "  if (betaK != 0.0)\n    res += this->getRayleighDampingForces();\n  return res;"),
    "cst_copy_back": _inc(
        "  static Vector res(6);\n  if (r == 0.0) {\n    this->getResistingForce();\n    res = P;\n"
        "    if (betaK != 0.0)\n      res += this->getRayleighDampingForces();\n    P = res;\n    return P;\n  }\n"
        "  return P;"),
    "seeded_then_bound_ref": _inc(
        "  static Vector res(6);\n  res = P;\n  const Vector &v = this->getRayleighDampingForces();\n"
        "  res += v;\n  return res;"),
    "addVector_into_local": _inc(
        "  static Vector res(12);\n  res = this->getResistingForce();\n"
        "  res.addVector(1.0, this->getRayleighDampingForces(), 1.0);\n  return res;"),
    "no_rayleigh_at_all": _inc("  P.addVector(1.0, Q, -1.0);\n  return P;"),
}


@pytest.mark.parametrize("case", sorted(GOOD))
def test_passes(tmp_path, case):
    assert _run(tmp_path, GOOD[case]) == []



BAD = {
    # the #562 pattern itself
    "static_P_plus_equals": _inc(
        "  this->getResistingForce();\n  if (betaK != 0.0)\n    P += this->getRayleighDampingForces();\n  return P;"),
    "addVector_into_static": _inc(
        "  P = this->getResistingForce();\n  P.addVector(1.0, this->getRayleighDampingForces(), 1.0);\n  return P;"),
    # review holes
    "snapshot_after_bind (L3-1)": _inc(
        "  this->getResistingForce();\n  const Vector &v = this->getRayleighDampingForces();\n"
        "  static Vector res(6);\n  res = P;\n  res += v;\n  P = res;\n  return P;"),
    "local_reference_alias": _inc(
        "  Vector &buf = P;\n  buf = this->getResistingForce();\n  buf += this->getRayleighDampingForces();\n  return buf;"),
    "self_assign_plus": _inc(
        "  this->getResistingForce();\n  P = P + this->getRayleighDampingForces();\n  return P;"),
    "bound_far_above": _inc(
        "  const Vector &v = this->getRayleighDampingForces();\n  int a = 0;\n  int b = 1;\n  int c = 2;\n"
        "  int d = 3;\n  P += v;\n  return P;"),
    "plain_assign_bind (L3-2)": _inc(
        "  static Vector damp(6);\n  damp = this->getRayleighDampingForces();\n  P += damp;\n  return P;"),
    "copy_ctor_bind (L3-2)": _inc(
        "  Vector damp(this->getRayleighDampingForces());\n  P.addVector(1.0, damp, 1.0);\n  return P;"),
    "pointer_arrow_target (L3-3)": _inc(
        "  theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);\n  return *theVector;"),
    "pointer_deref_target (L3-3)": _inc(
        "  *theVector += this->getRayleighDampingForces();\n  return *theVector;"),
    "this_member_target (L3-3)": _inc(
        "  this->P += this->getRayleighDampingForces();\n  return P;"),
    "one_line_if (L3-4)": _inc(
        "  if (betaK != 0.0) P += this->getRayleighDampingForces();\n  return P;"),
    "braced_one_line_if (L3-4)": _inc(
        "  if (betaK != 0.0) { P += this->getRayleighDampingForces(); }\n  return P;"),
    "multi_line_statement (L3-6)": _inc(
        "  P.addVector(1.0,\n              this->getRayleighDampingForces(),\n              1.0);\n  return P;"),
    "unseeded_local": _inc(
        "  static Vector res(6);\n  res += this->getRayleighDampingForces();\n  return res;"),
}


@pytest.mark.parametrize("case", sorted(BAD))
def test_flags(tmp_path, case):
    out = _run(tmp_path, BAD[case])
    assert len(out) == 1 and out[0].startswith("L1 "), out


def test_indented_header_does_not_borrow_previous_locals(tmp_path):
    # L3-7: an inline (indented) method must not inherit locals of the function above it.
    body = ("void Elem::other(void)\n{\n  static Vector res(6);\n  res = P;\n}\n"
            "class Inner {\n  public:\n    const Vector &f(void)\n    {\n"
            "      res += this->getRayleighDampingForces();\n      return res;\n    }\n};\n")
    out = _run(tmp_path, body)
    assert len(out) == 1 and "Inner::f" in out[0], out


def test_ignores_comments_and_vanilla(tmp_path):
    commented = _inc("  // P += this->getRayleighDampingForces();\n  /* P += this->getRayleighDampingForces(); */\n  return P;")
    assert _run(tmp_path / "a", commented) == []
    assert _run(tmp_path / "b", BAD["static_P_plus_equals"], stamped=False) == []


def test_waiver(tmp_path):
    base = BAD["static_P_plus_equals"]
    stmt = "    P += this->getRayleighDampingForces();"
    good = base.replace(stmt, "    // ladruno-lint: rayleigh-ok P is never refilled by getTangentStiff here\n" + stmt)
    short = base.replace(stmt, "    // ladruno-lint: rayleigh-ok fine\n" + stmt)
    assert _run(tmp_path / "a", good) == []
    out = _run(tmp_path / "b", short)
    assert len(out) == 1 and "too short" in out[0]


def test_stale_waiver_is_a_finding(tmp_path):
    # L3-9: a waiver left above code that no longer needs it
    body = GOOD["snapshot_static_local"].replace(
        "    res += this", "    // ladruno-lint: rayleigh-ok left over from the old code path\n    res += this")
    out = _run(tmp_path, body)
    assert len(out) == 1 and "stale rayleigh-ok" in out[0], out
