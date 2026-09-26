#!/usr/bin/env python3
"""WP-124 mutation rows: prove tests/test_ladruno_element_shell_helpers.py fails on one-line
breaks of each shared helper in SRC/element/LadrunoElementShell.h.

Reuses the WP-123 driver (../wp123_undamped/mutation_rows.py) unchanged: exact count-checked
edits, `build.bat OpenSeesPy` per row, each test file in its own process, every edited source
restored in a `finally`. Only the ROWS and the TESTS differ. After the rows, run a full
5-target build.bat (a named target leaves the other dist binaries stale -- AGENTS.md rule 2).

    <py3.12> mutation_rows.py [ROW ...]        (default: every row)
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("wp123_rows", HERE.parent / "wp123_undamped" / "mutation_rows.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

ROOT = M.ROOT
SHELL = ROOT / "SRC/element/LadrunoElementShell.h"
MUT_H = ROOT / "SRC/Ladruno_mutation.h"
BRICK = ROOT / "SRC/element/ladrunoBrick/LadrunoBrick.cpp"

ROWS = {
    # --- stage 1: Ki cache
    "K1": ("cacheKi stores a ZERO matrix instead of the formed K0",
           [(SHELL, "lit", "  Ki = new Matrix(formed);\n",
             "  Ki = new Matrix(formed.noRows(), formed.noCols());   // MUTATION K1\n")]),
    "K3": ("dropKi keeps the stale Ki",
           [(SHELL, "lit", "  if (Ki != 0) {\n    delete Ki;\n    Ki = 0;\n  }\n",
             "  (void)Ki;   // MUTATION K3\n")]),
    "K4": ("cacheKi caches but returns the SCRATCH (the C1 shape) -- expected EQUIVALENT in a "
           "default build: every consumer copies the reference at once",
           [(SHELL, "lit", "  Ki = new Matrix(formed);\n  return *Ki;\n",
             "  Ki = new Matrix(formed);\n  return formed;   // MUTATION K4\n")]),
    # --- C1 evidence: the ADR-87 CONTINUUM tangent gate in IDENT mode, with and without the fix
    "C1a": ("CONTINUUM=IDENT mutant + the C1 fix REVERTED (std/bbar returns the scratch) -- "
            "the probe must FAIL: the mutation never reaches the caller",
            [(MUT_H, "lit", "#  define LADRUNO_MUTATE_CONTINUUM LADRUNO_MUT_NONE\n",
              "#  define LADRUNO_MUTATE_CONTINUUM LADRUNO_MUT_IDENT   // MUTATION C1a\n"),
             (BRICK, "lit", "  return *Ki;\n}\n\n//----------------------------------------------------------------------\n"
                            "void  LadrunoBrick::zeroLoad",
              "  return stiff;   // MUTATION C1a: pre-fix\n}\n\n"
              "//----------------------------------------------------------------------\n"
              "void  LadrunoBrick::zeroLoad")]),
    "C1b": ("CONTINUUM=IDENT mutant with the C1 fix -- the probe must PASS",
            [(MUT_H, "lit", "#  define LADRUNO_MUTATE_CONTINUUM LADRUNO_MUT_NONE\n",
              "#  define LADRUNO_MUTATE_CONTINUUM LADRUNO_MUT_IDENT   // MUTATION C1b\n")]),
}

TESTS = ["tests/test_ladruno_element_shell_helpers.py"]
ROW_TESTS = {"C1a": ["Ladruno_implementation/wp124_shells/probe_c1_ident.py"],
             "C1b": ["Ladruno_implementation/wp124_shells/probe_c1_ident.py"]}


def _crlf_aware(edits):
    """The WP-123 apply() matches literals verbatim; our sources are CRLF on disk."""
    out = []
    for path, kind, old, new in edits:
        if kind == "lit" and b"\r\n" in path.read_bytes():
            old, new = old.replace("\n", "\r\n"), new.replace("\n", "\r\n")
        out.append((path, kind, old, new))
    return out


if __name__ == "__main__":
    M.ROWS = {k: (d, _crlf_aware(e)) for k, (d, e) in ROWS.items()}
    rows = sys.argv[1:] or list(ROWS)
    # the WP-123 driver runs ONE test list for all rows: group the rows by their list
    groups = {}
    for r in rows:
        groups.setdefault(tuple(ROW_TESTS.get(r, TESTS)), []).append(r)
    rc = 0
    for tests, grp in groups.items():
        M.TESTS = list(tests)
        sys.argv[1:] = grp
        rc |= M.main()
    raise SystemExit(rc)
