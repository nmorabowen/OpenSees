"""Self-test for quirk-lint rule `wipe` (alias L2, WP-115) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_wipe.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


SINGLETON = kit.STAMP + """
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


def _run(tmp_path, files):
    root = kit.tree(tmp_path, files)
    used = set()
    return cq.check_wipe(root, kit.rel(root), used) + cq.check_stale_waivers(root, kit.rel(root), used)


def test_flags_a_singleton_never_reset(tmp_path):
    out = _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON})
    assert len(out) == 1 and "'Globals' is not reset on wipe" in out[0]


def test_reset_through_a_free_helper_called_by_an_OPS_clearAll_hook(tmp_path):
    # the #841 shape
    hook = "void OPS_clearAllNDMaterial(void)\n{\n  resetGlobalsOnWipe();\n}\n"
    out = _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + FREE_RESET, "SRC/material/NDMaterial.cpp": hook})
    assert out == []


def test_reset_directly_in_domain_clearAll(tmp_path):
    hook = "int\nDomain::clearAll(void)\n{\n  Globals::instance().resetOnWipe();\n  return 0;\n}\n"
    assert _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/domain/Domain.cpp": hook}) == []


def test_reset_outside_any_wipe_path_does_not_count(tmp_path):
    assert len(_run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON + FREE_RESET})) == 1


def test_method_sharing_a_name_with_a_hook_call_does_not_count(tmp_path):
    # L2-1: Globals::clear() resets, and a hook calls some *other* clear(); not wired.
    hook = "void OPS_clearAllNDMaterial(void)\n{\n  theMap.clear();\n}\n"
    out = _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/material/NDMaterial.cpp": hook})
    assert len(out) == 1


def test_unrelated_clearAll_is_not_a_wipe_hook(tmp_path):
    # L2-2: only Domain/PartitionedDomain::clearAll and OPS_clearAll* are hooks.
    other = "void\nSomeCache::clearAll(void)\n{\n  Globals::instance().resetOnWipe();\n}\n"
    out = _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/misc/SomeCache.cpp": other})
    assert len(out) == 1


def test_reset_in_comment_or_if0_does_not_count(tmp_path):
    # L2-3
    hook = ("int\nDomain::clearAll(void)\n{\n  /* Globals::instance().resetOnWipe(); */\n"
            "#if 0\n  Globals::instance().resetOnWipe();\n#endif\n  return 0;\n}\n")
    out = _run(tmp_path, {"SRC/material/Mat.cpp": SINGLETON, "SRC/domain/Domain.cpp": hook})
    assert len(out) == 1


def test_waiver(tmp_path):
    waived = SINGLETON.replace("    static Globals &instance(void)",
                               "    // ladruno-lint: wipe-ok owner-scoped, cleared in the owner dtor\n"
                               "    static Globals &instance(void)")
    assert _run(tmp_path, {"SRC/material/Mat.cpp": waived}) == []
