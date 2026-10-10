#!/usr/bin/env python3
"""WP-139 mutation rows: prove tests/test_ladrunoBrick_lumped_inertia.py fails on one-line breaks
of the lumped inertia residual. Reuses the WP-123 driver (count-checked edits, `build.bat
OpenSeesPy` per row, sources restored in a `finally`); only ROWS and TESTS differ.

    <py3.12> mutation_rows.py [ROW ...]        (default: every row)
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("wp123_rows", HERE.parent / "wp123_undamped" / "mutation_rows.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
BRICK = M.ROOT / "SRC/element/ladrunoBrick/LadrunoBrick.cpp"

ROWS = {
    "L1": ("the pre-WP-139 hybrid: the -lumped residual is the CONSISTENT mass again",
           [(BRICK, "lit", "      if (massType == 1)                  // Ladruno (WP-139): lumped -- see after the GP loop\n",
             "      if (false)   // MUTATION L1\n"),
            (BRICK, "lit", "  if (massType == 1) {\n    int jj = 0;\n",
             "  if (false) {   // MUTATION L1\n    int jj = 0;\n")]),
    "L2": ("the lumped residual reads the COMMITTED acceleration instead of the trial one",
           [(BRICK, "lit", "      const Vector &aj = nodePointers[j]->getTrialAccel();\n      for (int p = 0; p < ndf; p++)\n        resid(jj + p) += mL[j] * aj(p);\n",
             "      const Vector &aj = nodePointers[j]->getAccel();   // MUTATION L2\n      for (int p = 0; p < ndf; p++)\n        resid(jj + p) += mL[j] * aj(p);\n")]),
}
TESTS = ["tests/test_ladrunoBrick_lumped_inertia.py"]


def _crlf_aware(edits):
    out = []
    for path, kind, old, new in edits:
        if kind == "lit" and b"\r\n" in path.read_bytes():
            old, new = old.replace("\n", "\r\n"), new.replace("\n", "\r\n")
        out.append((path, kind, old, new))
    return out


if __name__ == "__main__":
    M.ROWS = {k: (d, _crlf_aware(e)) for k, (d, e) in ROWS.items()}
    M.TESTS = TESTS
    if len(sys.argv) == 1:
        sys.argv += list(ROWS)
    raise SystemExit(M.main())
