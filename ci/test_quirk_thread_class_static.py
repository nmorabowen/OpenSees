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


# --- the opt-in is inherited (critic of record, WP-258): wrappers and bases are scanned too ------
BASE_H = ("class SoilBase : public NDMaterial {\npublic:\n"
          "  bool ladrunoThreadSafeUpdate(void) const { return true; }\n};\n")


def _tree(tmp_path, files):
    root = kit.tree(tmp_path, files)
    used = set()
    out = cq.check_thread_class_static(root, kit.rel(root), used) + \
        cq.check_thread_static(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def test_flags_an_inheriting_wrapper(tmp_path):
    """The ManzariDafalias3D shape: the wrapper never spells the opt-in, it inherits it."""
    out = _tree(tmp_path, {
        "SRC/material/nD/SoilBase.h": BASE_H,
        "SRC/material/nD/Soil3D.h": "class Soil3D : public SoilBase {\nprivate:\n  static Vector mSigma_M;\n};\n",
    })
    assert len(out) == 1 and out[0].startswith("thread-class-static SRC/material/nD/Soil3D.h:3:"), out


def test_a_refusing_descendant_is_not_scanned(tmp_path):
    """The ManzariDafaliasRO shape: a descendant that answers false keeps its statics."""
    out = _tree(tmp_path, {
        "SRC/material/nD/SoilBase.h": BASE_H,
        "SRC/material/nD/SoilRO.h": ("class SoilRO : public SoilBase {\npublic:\n"
                                     "  bool ladrunoThreadSafeUpdate(void) const { return false; }\n"
                                     "private:\n  static Vector s;\n};\n"),
    })
    assert out == [], out


def test_a_base_whose_leaves_opt_in_is_scanned(tmp_path):
    """The DruckerPrager shape: the leaves opt in; the base's .cpp holds the integrator."""
    out = _tree(tmp_path, {
        "SRC/material/nD/DPBase.h": "class DPBase : public NDMaterial {\n};\n",
        "SRC/material/nD/DPBase.cpp": ("void\nDPBase::plastic_integrator(void)\n{\n"
                                       "  static int count = 0;\n  count++;\n}\n"),
        "SRC/material/nD/DPLeaf.h": ("class DPLeaf : public DPBase {\npublic:\n"
                                     "  bool ladrunoThreadSafeUpdate(void) const { return true; }\n};\n"),
    })
    assert len(out) == 1 and out[0].startswith("thread-static SRC/material/nD/DPBase.cpp:4:"), out


def test_a_qualified_call_does_not_pull_in_another_class_file(tmp_path):
    out = _tree(tmp_path, {
        "SRC/material/nD/SoilBase.h": BASE_H,
        "SRC/material/nD/Other.cpp": ("int\nOther::f(void)\n{\n  static int k = 0;\n"
                                      "  opserr << \"SoilBase::f\";\n  SoilBase::g();\n  return k;\n}\n"),
    })
    assert out == [], out
