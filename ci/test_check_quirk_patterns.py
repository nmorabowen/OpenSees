"""Self-test for ci/check_quirk_patterns.py (WP-115).

Each rule must flag the pattern it exists for and pass the sanctioned fix, on
synthetic fork-stamped sources. A lint nobody has seen fail is not a gate.
Run: pytest -q ci/test_check_quirk_patterns.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_quirk_patterns as cq  # noqa: E402

STAMP = "// LADRUNO-HEADER-START\n// LADRUNO-HEADER-END\n"

BUGGY_562 = STAMP + """
const Vector &Elem::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  if (betaK != 0.0)
    P += this->getRayleighDampingForces();
  return P;
}
"""

SNAPSHOT_FIX = STAMP + """
const Vector &Elem::getResistingForceIncInertia(void)
{
  static Vector res(6);
  res = P;
  if (betaK != 0.0)
    res += this->getRayleighDampingForces();
  P = res;
  return P;
}
"""

VIA_REF = STAMP + """
const Vector &Elem::getResistingForceIncInertia()
{
    if (betaK != 0.0) {
        const Vector &v = this->getRayleighDampingForces();
        P_return += v;
    }
    return P_return;
}
"""

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
    int n = 0;
};
}
"""

RESET_HELPER = """
void
resetGlobalsOnWipe(void)
{
    Globals::instance().resetOnWipe();
}
"""

WIPE_HOOK = """
void OPS_clearAllNDMaterial(void)
{
  resetGlobalsOnWipe();
}
"""


def _tree(tmp_path, files):
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return tmp_path


def _rel(root):
    return lambda p: p.resolve().relative_to(root.resolve()).as_posix()


def test_l1_flags_the_562_pattern(tmp_path):
    root = _tree(tmp_path, {"SRC/element/Elem.cpp": BUGGY_562})
    out = cq.check_rayleigh(root, _rel(root))
    assert len(out) == 1 and "'P' is not declared" in out[0]


def test_l1_passes_the_snapshot_fix(tmp_path):
    root = _tree(tmp_path, {"SRC/element/Elem.cpp": SNAPSHOT_FIX})
    assert cq.check_rayleigh(root, _rel(root)) == []


def test_l1_follows_a_bound_reference(tmp_path):
    root = _tree(tmp_path, {"SRC/element/Elem.cpp": VIA_REF})
    assert len(cq.check_rayleigh(root, _rel(root))) == 1


def test_l1_waiver_needs_a_reason(tmp_path):
    short = BUGGY_562.replace("    P += this->getRayleighDampingForces();",
                              "    // ladruno-lint: rayleigh-ok fine\n    P += this->getRayleighDampingForces();")
    good = BUGGY_562.replace("    P += this->getRayleighDampingForces();",
                             "    // ladruno-lint: rayleigh-ok P is never refilled by getTangentStiff\n"
                             "    P += this->getRayleighDampingForces();")
    root = _tree(tmp_path / "a", {"SRC/element/Elem.cpp": short})
    assert "too short" in cq.check_rayleigh(root, _rel(root))[0]
    root = _tree(tmp_path / "b", {"SRC/element/Elem.cpp": good})
    assert cq.check_rayleigh(root, _rel(root)) == []


def test_l1_ignores_vanilla_files(tmp_path):
    root = _tree(tmp_path, {"SRC/element/Elem.cpp": BUGGY_562.replace(STAMP, "")})
    assert cq.check_rayleigh(root, _rel(root)) == []


def test_l2_flags_a_singleton_never_reset(tmp_path):
    root = _tree(tmp_path, {"SRC/material/Mat.cpp": SINGLETON})
    out = cq.check_wipe(root, _rel(root))
    assert len(out) == 1 and "'Globals' is not reset on wipe" in out[0]


def test_l2_reset_reached_from_a_wipe_hook_through_one_helper(tmp_path):
    # The #841 shape: the reset lives in a helper that OPS_clearAllNDMaterial calls.
    root = _tree(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + RESET_HELPER,
                            "SRC/material/NDMaterial.cpp": WIPE_HOOK})
    assert cq.check_wipe(root, _rel(root)) == []


def test_l2_reset_outside_any_wipe_path_does_not_count(tmp_path):
    root = _tree(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + RESET_HELPER})
    assert len(cq.check_wipe(root, _rel(root))) == 1


def test_l2_waiver(tmp_path):
    waived = SINGLETON.replace("    static Globals &instance(void)",
                               "    // ladruno-lint: wipe-ok owner-scoped, cleared in the owner dtor\n"
                               "    static Globals &instance(void)")
    root = _tree(tmp_path, {"SRC/material/Mat.cpp": waived})
    assert cq.check_wipe(root, _rel(root)) == []


def test_l3_flags_an_orphaned_pointer(tmp_path):
    root = _tree(tmp_path, {
        "Ladruno_implementation/LEDGER_quirks.md": "### `wipe()` does NOT recreate the Domain\n",
        ".claude/skills/g/SKILL.md": 'Quirks: "`wipe()` does NOT recreate the Domain", "renamed heading".\n',
    })
    out = cq.check_pointers(root, _rel(root))
    assert len(out) == 1 and "renamed heading" in out[0]
