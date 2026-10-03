"""Self-test for quirk-lint rule `dead-decl` (alias L9, WP-158) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_dead_decl.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _run(tmp_path, body, stamped=False, where="SRC/material/nD/Mat.cpp"):
    src = (kit.STAMP if stamped else "") + body
    root = kit.tree(tmp_path, {where: src})
    used = set()
    out = cq.check_dead_decl(root, kit.rel(root), used)
    return out + cq.check_stale_waivers(root, kit.rel(root), used)


def _fe(body):
    return ("void ManzariDafalias::ForwardEuler(const Vector& CurStress, double& G)\n{\n"
            "    double p = one3 * GetTrace(CurStress);\n" + body + "\n}\n")


def test_flags_the_wp158_incident_in_a_vanilla_file(tmp_path):
    """The exact pre-WP-158 lines of vanilla ManzariDafalias::ForwardEuler."""
    out = _run(tmp_path, _fe("    Vector r(6);\n    if (p > small)\n        Vector r = GetDevPart(CurStress) / p;\n"
                            "    double x = DoubleDot2_2_Contr(n, r);"))
    assert len(out) == 1 and out[0].startswith("L9 SRC/material/nD/Mat.cpp:5:"), out
    assert "re-declares 'r'" in out[0] and "outer 'r' is never assigned" in out[0]


def test_passes_the_wp158_fix(tmp_path):
    assert _run(tmp_path, _fe("    Vector r(6);\n    if (p > small)\n        r = GetDevPart(CurStress) / p;")) == []


@pytest.mark.parametrize("body", [
    "    Vector r(6);\n    if (p > small) Vector r = GetDevPart(CurStress) / p;",         # one line
    "    Vector r(6);\n    if (a) { }\n    else\n        Vector r(GetDevPart(CurStress));",  # else, ctor form
    "    for (int i = 0; i < 3; i++)\n        double G = 2.0 * i;",                     # shadows a PARAMETER
    "    Matrix M(3,3);\n    while (k-- > 0)\n        const Matrix &M = foo();",              # reference form
])
def test_flags_other_shadowing_forms(tmp_path, body):
    out = _run(tmp_path, _fe(body))
    assert len(out) == 1 and out[0].startswith("L9 "), out


@pytest.mark.parametrize("body", [
    # vanilla Domain::initialize: a dead copy kept ON PURPOSE (it forms Ki), nothing shadowed
    "    while ((e = it()) != 0)\n        Matrix initM(e->getInitialStiff());",
    "    double x = 1.0;\n    if (p > small)\n        return x;",                   # a keyword, not a type
    "    Vector r(6);\n    if (p > small)\n        delete r;",
    "    Vector r(6);\n    if (p > small) {\n        Vector r = GetDevPart(CurStress) / p;\n        use(r);\n    }",
    "    if (p > small)\n        q = p * 2.0;",
    "    if (p > small)\n        delete thePtr;",
    "    if (p > small)\n        opserr << \"Vector r = x\" << endln;",
    "    // if (p > small) Vector r = x;\n    Vector r(6);",
    # vanilla elementAPI_TCL.cpp / PythonWrapper.cpp: strip_prefix leaves `if constexpr (...) {`
    "    double v;\n    if constexpr (A) {\n        v = 1;\n    } else if constexpr (B) {\n        v = 2;\n    }",
    "    Vector r(6);\n    {\n        Vector r = x;\n        use(r);\n    }",       # nested block, no control prefix
])
def test_passes_non_shadowing_or_braced_forms(tmp_path, body):
    assert _run(tmp_path, _fe(body)) == []


def test_waiver(tmp_path):
    ok = _fe("    Vector r(6);\n    // ladruno-lint: decl-ok a scoped copy is wanted here for the destructor\n"
             "    if (p > small) Vector r = x;")
    assert _run(tmp_path, ok, stamped=True) == []
    short = _fe("    Vector r(6);\n    if (p > small) Vector r = x;   // ladruno-lint: decl-ok ok")
    out = _run(tmp_path, short, stamped=True)
    assert len(out) == 1 and "reason too short" in out[0]


def test_stale_waiver(tmp_path):
    stale = _fe("    Vector r(6);\n    // ladruno-lint: decl-ok this used to shadow the outer r\n"
                "    if (p > small) r = x;")
    out = _run(tmp_path, stale, stamped=True)
    assert len(out) == 1 and out[0].startswith("W ") and "stale decl-ok" in out[0]
