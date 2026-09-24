"""Self-test for ci/check_quirk_patterns.py (WP-115).

Each rule must flag the pattern it exists for and pass the sanctioned fix, on
synthetic fork-stamped sources. Every hole found by the WP-115 adversarial
review has a case here, so it cannot silently reopen.
Run: pytest -q ci/test_check_quirk_patterns.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_quirk_patterns as cq  # noqa: E402

STAMP = "// LADRUNO-HEADER-START\n// LADRUNO-HEADER-END\n"


def _tree(tmp_path, files):
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return tmp_path


def _rel(root):
    return lambda p: p.resolve().relative_to(root.resolve()).as_posix()


def _l1(tmp_path, body, stamped=True):
    src = (STAMP if stamped else "") + body
    root = _tree(tmp_path, {"SRC/element/Elem.cpp": src})
    used = set()
    out = cq.check_rayleigh(root, _rel(root), used)
    return out + cq.check_stale_waivers(root, _rel(root), used)


def _inc(body):
    return "const Vector &Elem::getResistingForceIncInertia(void)\n{\n" + body + "\n}\n"


# ---------------------------------------------------------------- L1: must pass
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
def test_l1_passes(tmp_path, case):
    assert _l1(tmp_path, GOOD[case]) == []


# ---------------------------------------------------------------- L1: must flag
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
def test_l1_flags(tmp_path, case):
    out = _l1(tmp_path, BAD[case])
    assert len(out) == 1 and out[0].startswith("L1 "), out


def test_l1_indented_header_does_not_borrow_previous_locals(tmp_path):
    # L3-7: an inline (indented) method must not inherit locals of the function above it.
    body = ("void Elem::other(void)\n{\n  static Vector res(6);\n  res = P;\n}\n"
            "class Inner {\n  public:\n    const Vector &f(void)\n    {\n"
            "      res += this->getRayleighDampingForces();\n      return res;\n    }\n};\n")
    out = _l1(tmp_path, body)
    assert len(out) == 1 and "Inner::f" in out[0], out


def test_l1_ignores_comments_and_vanilla(tmp_path):
    commented = _inc("  // P += this->getRayleighDampingForces();\n  /* P += this->getRayleighDampingForces(); */\n  return P;")
    assert _l1(tmp_path / "a", commented) == []
    assert _l1(tmp_path / "b", BAD["static_P_plus_equals"], stamped=False) == []


def test_l1_waiver(tmp_path):
    base = BAD["static_P_plus_equals"]
    stmt = "    P += this->getRayleighDampingForces();"
    good = base.replace(stmt, "    // ladruno-lint: rayleigh-ok P is never refilled by getTangentStiff here\n" + stmt)
    short = base.replace(stmt, "    // ladruno-lint: rayleigh-ok fine\n" + stmt)
    assert _l1(tmp_path / "a", good) == []
    out = _l1(tmp_path / "b", short)
    assert len(out) == 1 and "too short" in out[0]


def test_stale_waiver_is_a_finding(tmp_path):
    # L3-9: a waiver left above code that no longer needs it
    body = GOOD["snapshot_static_local"].replace(
        "    res += this", "    // ladruno-lint: rayleigh-ok left over from the old code path\n    res += this")
    out = _l1(tmp_path, body)
    assert len(out) == 1 and "stale rayleigh-ok" in out[0], out


# ---------------------------------------------------------------- L2
SINGLETON = STAMP + """
namespace {
class Globals
{
  public:
    static Globals &instance(void)
    {
        static Globals theInstance;
        return theInstance;
    }
    void resetOnWipe(void) { n = 0; }
    void clear(void) { Globals::instance().resetOnWipe(); }
    int n = 0;
};
}
"""
FREE_RESET = "\nvoid\nresetGlobalsOnWipe(void)\n{\n    Globals::instance().resetOnWipe();\n}\n"


def _l2(tmp_path, files):
    root = _tree(tmp_path, files)
    used = set()
    return cq.check_wipe(root, _rel(root), used) + cq.check_stale_waivers(root, _rel(root), used)


def test_l2_flags_a_singleton_never_reset(tmp_path):
    out = _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON})
    assert len(out) == 1 and "'Globals' is not reset on wipe" in out[0]


def test_l2_reset_through_a_free_helper_called_by_an_OPS_clearAll_hook(tmp_path):
    # the #841 shape
    hook = "void OPS_clearAllNDMaterial(void)\n{\n  resetGlobalsOnWipe();\n}\n"
    out = _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + FREE_RESET, "SRC/material/NDMaterial.cpp": hook})
    assert out == []


def test_l2_reset_directly_in_domain_clearAll(tmp_path):
    hook = "int\nDomain::clearAll(void)\n{\n  Globals::instance().resetOnWipe();\n  return 0;\n}\n"
    assert _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/domain/Domain.cpp": hook}) == []


def test_l2_reset_outside_any_wipe_path_does_not_count(tmp_path):
    assert len(_l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + FREE_RESET})) == 1


def test_l2_method_sharing_a_name_with_a_hook_call_does_not_count(tmp_path):
    # L2-1: Globals::clear() resets, and a hook calls some *other* clear(); not wired.
    hook = "void OPS_clearAllNDMaterial(void)\n{\n  theMap.clear();\n}\n"
    out = _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/material/NDMaterial.cpp": hook})
    assert len(out) == 1


def test_l2_unrelated_clearAll_is_not_a_wipe_hook(tmp_path):
    # L2-2: only Domain/PartitionedDomain::clearAll and OPS_clearAll* are hooks.
    other = "void\nSomeCache::clearAll(void)\n{\n  Globals::instance().resetOnWipe();\n}\n"
    out = _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/misc/SomeCache.cpp": other})
    assert len(out) == 1


def test_l2_reset_in_comment_or_if0_does_not_count(tmp_path):
    # L2-3
    hook = ("int\nDomain::clearAll(void)\n{\n  /* Globals::instance().resetOnWipe(); */\n"
            "#if 0\n  Globals::instance().resetOnWipe();\n#endif\n  return 0;\n}\n")
    out = _l2(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/domain/Domain.cpp": hook})
    assert len(out) == 1


def test_l2_waiver(tmp_path):
    waived = SINGLETON.replace("    static Globals &instance(void)",
                               "    // ladruno-lint: wipe-ok owner-scoped, cleared in the owner dtor\n"
                               "    static Globals &instance(void)")
    assert _l2(tmp_path, {"SRC/material/Mat.cpp": waived}) == []


# ---------------------------------------------------------------- L3
def test_l3_flags_an_orphaned_pointer(tmp_path):
    root = _tree(tmp_path, {
        "Ladruno_implementation/LEDGER_quirks.md": "### `wipe()` does NOT recreate the Domain\n",
        ".claude/skills/g/SKILL.md": 'Quirks: "`wipe()` does NOT recreate the Domain", "renamed heading".\n',
    })
    out = cq.check_pointers(root, _rel(root))
    assert len(out) == 1 and "renamed heading" in out[0]
