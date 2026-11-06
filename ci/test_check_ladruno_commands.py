"""Self-test for ci/check_ladruno_commands.py (WP-168, the command-hook gate).
Run: pytest -q ci/test_check_ladruno_commands.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_ladruno_commands as cc  # noqa: E402

TABLE = (
    "// LADRUNO_COMMAND(\"name\", dl, classic)  <- a comment, not a row\n"
    'LADRUNO_COMMAND("contactSurface", OPS_LadrunoContactSurface, ladrunoContactSurface)\n'
    'LADRUNO_COMMAND("ladrunoNumbering", LADRUNO_NONE, TclCommand_ladrunoNumbering)\n'
)
CLASSIC = ('#include "LadrunoCommandsClassicTcl.h"   // Ladruno WP-168\n'
           "int OpenSeesAppInit(Tcl_Interp *interp) {\n"
           "    Ladruno_registerCommands(interp);   // Ladruno WP-168\n"
           '    Tcl_CreateCommand(interp, "wipe", &wipeModel,\n'
           "        (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);\n"
           "}\n")
TCLWRAP = ('#include "LadrunoCommandsTclWrapper.h"\n'
           "void TclWrapper::addOpenSeesCommands(Tcl_Interp* interp) {\n"
           '    addCommand(interp,"wipe", &Tcl_ops_wipe);\n'
           "    Ladruno_registerCommands(this, interp);\n}\n")
PYWRAP = ('#include "LadrunoCommandsPython.h"\n'
          "void PythonWrapper::addOpenSeesCommands() {\n"
          '    addCommand("wipe", &Py_ops_wipe);\n'
          "    Ladruno_registerCommands(this);\n}\n")


def _tree(tmp_path, override=None):
    files = {
        "SRC/interpreter/LadrunoCommandTable.h": TABLE,
        "SRC/tcl/commands.cpp": CLASSIC,
        "SRC/interpreter/TclWrapper.cpp": TCLWRAP,
        "SRC/interpreter/PythonWrapper.cpp": PYWRAP,
    }
    files.update(override or {})
    return kit.tree(tmp_path, files)


def _run(tmp_path, override=None):
    rows, findings, _ = cc.run(_tree(tmp_path, override))
    return rows, findings


def test_clean_tree_passes_and_comment_rows_are_not_rows(tmp_path):
    rows, findings = _run(tmp_path)
    assert findings == []
    assert [r[0] for r in rows] == ["contactSurface", "ladrunoNumbering"]


def test_hand_registration_of_a_table_row_is_flagged(tmp_path):
    # the pre-WP-168 shape: the same verb registered by hand in the engine file
    bad = PYWRAP.replace("    Ladruno_registerCommands(this);",
                         '    addCommand("contactSurface", &Py_ops_LadrunoContactSurface);\n'
                         "    Ladruno_registerCommands(this);")
    _, f = _run(tmp_path, {"SRC/interpreter/PythonWrapper.cpp": bad})
    assert len(f) == 1 and f[0].startswith("H1 SRC/interpreter/PythonWrapper.cpp:4:") and "'contactSurface'" in f[0]


def test_a_new_ladruno_named_command_outside_the_table_is_flagged(tmp_path):
    bad = CLASSIC.replace("    Tcl_CreateCommand(interp, \"wipe\"",
                          "    Tcl_CreateCommand(interp, \"ladrunoFoo\", &ladrunoFoo,\n"
                          "        (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);\n"
                          "    Tcl_CreateCommand(interp, \"wipe\"")
    _, f = _run(tmp_path, {"SRC/tcl/commands.cpp": bad})
    assert len(f) == 1 and "'ladrunoFoo' is a Ladruno command name" in f[0]


def test_a_ladruno_marked_registration_on_the_continuation_line_is_flagged(tmp_path):
    # `profiler` carried its `// Ladruno` mark on the second line of the statement
    src = ('void f(Tcl_Interp *interp) {\n'
           '    Tcl_CreateCommand(interp, "profiler", &TclCommand_profiler,\n'
           '        (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); // Ladruno\n}\n')
    _, f = _run(tmp_path, {"SRC/tcl/Other.cpp": src})
    assert len(f) == 1 and "carries a `// Ladruno` mark" in f[0]


def test_other_registration_surfaces_are_scanned(tmp_path):
    src = 'void g(Tcl_Interp *interp) { Tcl_CreateObjCommand(interp, "contactSurface", &x, NULL, NULL); }\n'
    _, f = _run(tmp_path, {"SRC/element/Foo.cpp": src})
    assert len(f) == 1 and f[0].startswith("H1 SRC/element/Foo.cpp:1:")


def test_commented_out_and_waived_registrations_are_not_flagged(tmp_path):
    src = ('void f(Tcl_Interp *interp) {\n'
           '    // addCommand(interp,"contact", &Tcl_ops_LadrunoContact);\n'
           '    /* Tcl_CreateCommand(interp, "ladrunoX", &x, NULL, NULL); */\n'
           '    // ladruno-command-hook-ok vendored test harness\n'
           '    Tcl_CreateCommand(interp, "ladrunoY", &y, NULL, NULL);\n}\n')
    _, f = _run(tmp_path, {"SRC/tcl/Other.cpp": src})
    assert f == []


def test_upstream_registration_without_a_mark_is_not_flagged(tmp_path):
    src = 'void f(Tcl_Interp *interp) { Tcl_CreateCommand(interp, "eigen", &eigenAnalysis, NULL, NULL); }\n'
    _, f = _run(tmp_path, {"SRC/tcl/Other.cpp": src})
    assert f == []


def test_an_engine_that_drops_the_hook_is_flagged(tmp_path):
    _, f = _run(tmp_path, {"SRC/interpreter/TclWrapper.cpp": TCLWRAP.replace(
        "    Ladruno_registerCommands(this, interp);\n", "")})
    assert len(f) == 1 and f[0].startswith("H2 SRC/interpreter/TclWrapper.cpp: Ladruno_registerCommands( called 0")


def test_an_engine_that_calls_the_hook_twice_is_flagged(tmp_path):
    twice = CLASSIC.replace("    Ladruno_registerCommands(interp);   // Ladruno WP-168\n",
                            "    Ladruno_registerCommands(interp);\n    Ladruno_registerCommands(interp);\n")
    _, f = _run(tmp_path, {"SRC/tcl/commands.cpp": twice})
    assert len(f) == 1 and "called 2 times" in f[0]


def test_an_engine_without_the_hook_header_is_flagged(tmp_path):
    _, f = _run(tmp_path, {"SRC/interpreter/PythonWrapper.cpp": PYWRAP.replace(
        '#include "LadrunoCommandsPython.h"\n', "")})
    assert len(f) == 1 and 'does not #include "LadrunoCommandsPython.h"' in f[0]


def test_malformed_table_rows_are_flagged(tmp_path):
    table = TABLE + (
        'LADRUNO_COMMAND("contactSurface", OPS_X, LADRUNO_NONE)\n'
        'LADRUNO_COMMAND("orphan", LADRUNO_NONE, LADRUNO_NONE)\n'
        "LADRUNO_COMMAND(oops, OPS_Y)\n")
    _, f = _run(tmp_path, {"SRC/interpreter/LadrunoCommandTable.h": table})
    text = "\n".join(f)
    assert len(f) == 3, text
    assert "duplicate command 'contactSurface'" in text
    assert "'orphan' is registered in no engine" in text
    assert "does not parse" in text


def test_name_on_a_continuation_line_is_flagged(tmp_path):
    src = ('void f(Tcl_Interp *interp) {\n'
           '    Tcl_CreateCommand(interp,\n'
           '                      "ladrunoSneaky", &x, NULL, NULL);\n}\n')
    _, f = _run(tmp_path, {"SRC/tcl/Other.cpp": src})
    assert len(f) == 1 and f[0].startswith("H1 SRC/tcl/Other.cpp:2:") and "'ladrunoSneaky'" in f[0]


def test_a_call_expression_as_the_first_argument_is_flagged(tmp_path):
    src = ('void f() {\n'
           '    Tcl_CreateCommand(theInterp(), "ladrunoSneaky2", &x, NULL, NULL);\n'
           '    addCommand(getInterp(0), "contactSurface", &y);\n}\n')
    _, f = _run(tmp_path, {"SRC/tcl/Other.cpp": src})
    assert len(f) == 2, f
    assert "'ladrunoSneaky2'" in f[0] and "'contactSurface'" in f[1]


def test_registrations_are_counted(tmp_path):
    _, _, scanned = cc.run(_tree(tmp_path))
    assert scanned == 3   # the three `wipe` registrations of the fixture engines


def test_the_real_tree_is_clean():
    root = Path(__file__).resolve().parent.parent
    rows, findings, scanned = cc.run(root)
    assert findings == [], "\n".join(findings)
    assert len(rows) >= 30
    # non-vacuity: a REGISTER regex that matched nothing would also report no findings
    assert scanned >= 600, f"only {scanned} registrations scanned"
