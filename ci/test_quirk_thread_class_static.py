"""Self-test for quirk-lint rule `thread-class-static` (WP-258) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_thread_class_static.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _header(members, optin="true"):
    return ("class SoilPS : public NDMaterial {\npublic:\n"
            "  const Vector& getStress(void);\n"
            "  bool ladrunoThreadSafeUpdate(void) const { return " + optin + "; }\n"
            "private:\n" + members + "\n};\n")


def _run(tmp_path, header, where="SRC/material/nD/SoilPS"):
    root = kit.tree(tmp_path, {where + ".h": header, where + ".cpp": "int x;\n"})
    used = set()
    out = cq.check_thread_class_static(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def test_flags_the_wp258_druckerprager_planestrain_buffers(tmp_path):
    """The pre-WP-258 DruckerPragerPlaneStrain return buffers."""
    out = _run(tmp_path, _header("  static Vector strain;\n  static Vector stress;\n  static Matrix tangent;"))
    assert len(out) == 3, out
    assert out[1].startswith("thread-class-static SRC/material/nD/SoilPS.h:7:") and "'stress'" in out[1], out


def test_passes_the_wp258_fix(tmp_path):
    assert _run(tmp_path, _header("  static thread_local Vector stress;\n  thread_local static Matrix t;")) == []


@pytest.mark.parametrize("decl", [
    "  static const double one3;",
    "  static constexpr int n = 6;",
    "  static std::atomic<int> warnCount;",
    "  static std::mutex m;",
    "  static int numInstances(void);",
    "  static const Vector& helper(const Vector &v);",
    "  Vector stress;",
])
def test_exempt_declarations(tmp_path, decl):
    assert _run(tmp_path, _header(decl)) == []


def test_function_scope_static_is_thread_static_not_this_rule(tmp_path):
    hdr = _header("  int y;").replace("  const Vector& getStress(void);",
                                      "  int f(void) {\n    static int k;\n    return k;\n  }")
    assert _run(tmp_path, hdr) == []


def test_a_refusing_material_is_not_scanned(tmp_path):
    assert _run(tmp_path, _header("  static Vector stress;", optin="false")) == []


def test_waiver_with_reason(tmp_path):
    hdr = _header("  // ladruno-lint: thread-class-ok stage flag, written only on the main thread\n"
                  "  static char unsigned mElastFlag;")
    assert _run(tmp_path, hdr) == []


def test_waiver_reason_too_short(tmp_path):
    out = _run(tmp_path, _header("  static Vector s; // ladruno-lint: thread-class-ok ok"))
    assert len(out) == 1 and "reason too short" in out[0], out


def test_stale_waiver_is_reported(tmp_path):
    hdr = _header("  // ladruno-lint: thread-class-ok stage flag, written only on the main thread\n"
                  "  static thread_local Vector s;")
    out = _run(tmp_path, kit.STAMP + hdr)
    assert len(out) == 1 and "stale thread-class-ok waiver" in out[0], out
