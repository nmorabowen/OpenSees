"""Self-test for quirk-lint rule `thread-static` (WP-239) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_thread_static.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402

OPT_IN = ("class Steel : public UniaxialMaterial {\npublic:\n"
          "  bool ladrunoThreadSafeUpdate(void) const { return !params.slip; }\n"
          "  bool ladrunoThreadSafeForm(void) const { return !params.slip; }\n};\n")
REFUSED = ("class Soil : public NDMaterial {\npublic:\n"
           "  bool ladrunoThreadSafeUpdate(void) const { return false; }\n};\n")


def _run(tmp_path, cpp, header=OPT_IN, stamped=False, where="SRC/material/uniaxial/Steel"):
    files = {where + ".h": header, where + ".cpp": (kit.STAMP if stamped else "") + cpp}
    root = kit.tree(tmp_path, files)
    used = set()
    out = cq.check_thread_static(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def _fn(body, name="namespace {\n  int rve_process_element(int iter)"):
    return name + "\n  {\n" + body + "\n    return 0;\n  }\n" + ("}\n" if name.startswith("namespace") else "")


def test_flags_the_wp239_asdsteel1d_incident(tmp_path):
    """The pre-WP-239 ASDSteel1D implicit-correction history, in an anonymous-namespace helper."""
    out = _run(tmp_path, _fn("    static std::array<Vector, 10> U_iterative = { Vector(7), Vector(7) };\n"
                             "    U_iterative[iter] = x;"))
    assert len(out) == 1 and out[0].startswith("thread-static SRC/material/uniaxial/Steel.cpp:4:"), out
    assert "rve_process_element" in out[0] and "U_iterative" in out[0]


def test_flags_the_wp239_shell_gettype_latch(tmp_path):
    out = _run(tmp_path, _fn("    static bool initialized = false;\n    if (!initialized) { array(0) = 1; initialized = true; }",
                             name="const ID& Steel::getType(void)"))
    assert len(out) == 1 and "static bool initialized" in out[0], out


def test_passes_the_wp239_fix(tmp_path):
    assert _run(tmp_path, _fn("    static thread_local std::array<Vector, 10> U_iterative = { Vector(7) };\n"
                              "    U_iterative[iter] = x;")) == []


@pytest.mark.parametrize("decl", [
    "    static const double tol = 1.0e-8;",
    "    static constexpr int n = 6;",
    "    static thread_local Vector a(3);",
    "    thread_local static Vector a(3);",
    "    static std::atomic<int> warnCount(0);",
    "    static std::mutex m;",
    "    Vector a(3);",
    "    // static Vector a(3);",
])
def test_passes_safe_declarations(tmp_path, decl):
    assert _run(tmp_path, _fn(decl)) == []


@pytest.mark.parametrize("name", [
    "void *OPS_Steel(void)",
    "int Steel::sendSelf(int tag, Channel &ch)",
    "int Steel::recvSelf(int tag, Channel &ch, FEM_ObjectBroker &b)",
    "Response *Steel::setResponse(const char **argv, int argc, OPS_Stream &s)",
    "int Steel::getResponse(int id, Information &info)",
    "void Steel::Print(OPS_Stream &s, int flag)",
    "int Steel::setParameter(const char **argv, int argc, Parameter &p)",
    "const Vector &Steel::getStressSensitivity(int g, bool c)",
    "Steel::Steel(int tag)",
])
def test_skips_functions_that_never_run_on_a_worker(tmp_path, name):
    assert _run(tmp_path, _fn("    static ID data(3);", name=name)) == []


def test_skips_a_class_that_refuses(tmp_path):
    assert _run(tmp_path, _fn("    static Vector n(6);"), header=REFUSED) == []


def test_scans_a_class_that_delegates(tmp_path):
    hdr = ("class Wrap : public NDMaterial {\npublic:\n  bool ladrunoThreadSafeUpdate(void) const "
           "{ return theMat != 0 && theMat->ladrunoThreadSafeUpdate(); }\n};\n")
    out = _run(tmp_path, _fn("    static Vector n(6);"), header=hdr)
    assert len(out) == 1 and out[0].startswith("thread-static "), out


def test_ignores_files_outside_src_material(tmp_path):
    assert _run(tmp_path, _fn("    static Matrix B(3, 8);"), where="SRC/element/quad/Quad") == []


def test_waiver(tmp_path):
    ok = _fn("    // ladruno-lint: thread-ok guarded by its own mutex on every access\n    static std::vector<int> v;")
    assert _run(tmp_path, ok, stamped=True) == []
    short = _fn("    static std::vector<int> v;   // ladruno-lint: thread-ok fine")
    out = _run(tmp_path, short, stamped=True)
    assert len(out) == 1 and "reason too short" in out[0], out


def test_stale_waiver(tmp_path):
    stale = _fn("    // ladruno-lint: thread-ok this used to be a shared scratch vector\n    Vector v(3);")
    out = _run(tmp_path, stale, stamped=True)
    assert len(out) == 1 and out[0].startswith("W ") and "stale thread-ok" in out[0], out
