#!/usr/bin/env python3
"""Gate the Tcl verification suite (G1).

The upstream suite only `puts` PASSED/FAILED into results.out and never sets an
exit code, so a failing Tcl check passes CI silently. Run this AFTER
`tclsh runVerificationSuite.tcl` (and any ladruno tcl tests) to turn a FAILED
line into a nonzero exit.

    python ci/check_tcl_results.py [path/to/results.out]
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "EXAMPLES" / "verification" / "results.out"


def main(argv):
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT
    if not path.exists():
        print(f"check_tcl_results: ERROR results file not found: {path}")
        return 1
    failed = [ln.rstrip() for ln in path.read_text().splitlines() if "FAILED" in ln]
    total = sum(1 for ln in path.read_text().splitlines() if ("PASSED" in ln or "FAILED" in ln))
    if failed:
        print(f"check_tcl_results: {len(failed)} FAILED of {total} checks:")
        for ln in failed:
            print("  " + ln)
        return 1
    if total == 0:
        # A gate that counts failures must refuse an empty result set. The
        # suite driver clears results.out with core-Tcl `open ... w` BEFORE any
        # OpenSees command runs, so a binary that dies on its first model
        # command (TCL_LIBRARY unset -> "Can't find a usable init.tcl" -> no
        # commands registered) leaves a fresh empty file — which this checker
        # blessed as "OK (0 checks passed)" while zero decks executed. See
        # Ladruno_internal/WORKFLOW_GOTCHAS.md §7.
        print("check_tcl_results: ERROR 0 checks ran — results.out is empty. "
              "The suite died before its first check (TCL_LIBRARY / init.tcl?). "
              "An empty result set is a FAILURE, not a pass.")
        return 1
    print(f"check_tcl_results: OK ({total} checks passed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
