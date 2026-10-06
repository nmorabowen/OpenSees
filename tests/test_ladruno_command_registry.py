"""Parity gate for the fork command table (WP-168, apeGmsh#1490, panel P9 R10).

Every fork-only interpreter command is one row of
SRC/interpreter/LadrunoCommandTable.h and is registered by the one
`Ladruno_registerCommands` hook each engine calls. The row's two columns say
where it must exist: `dl` (TclWrapper + Python) and `classic` (classic Tcl:
OpenSees.exe, OpenSeesSP, OpenSeesMP); LADRUNO_NONE means "not registered here".

This file proves the binaries agree with the table, in BOTH directions:

* every row the table puts in classic Tcl answers `info commands` in the classic
  exe, and every row it leaves out does not (a gap that closes without the
  table being updated is a finding too);
* the same for the Python module, via `hasattr`;
* no fork command is registered outside the hook (ci/check_ladruno_commands.py
  on the real tree), so the table cannot drift from the code.

Before WP-168 the three registrations were hand-maintained, and a verb could
silently miss an engine: the contact family was absent from classic Tcl until
ADR-78 P0.5, `ladrunoDR`/`ladrunoArcLength` until #729.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "SRC" / "interpreter" / "LadrunoCommandTable.h"
NONE = "LADRUNO_NONE"
ROW = re.compile(r'^\s*LADRUNO_COMMAND\s*\(\s*"([^"]+)"\s*,\s*(\w+)\s*,\s*(\w+)\s*\)', re.M)
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

sys.path.insert(0, str(REPO / "ci"))
import check_ladruno_commands  # noqa: E402


def _rows():
    rows = ROW.findall(TABLE.read_text(encoding="utf-8"))
    assert rows, f"no LADRUNO_COMMAND rows parsed from {TABLE}"
    return rows


ROWS = _rows()
CLASSIC_IN = sorted(n for n, _, c in ROWS if c != NONE)
CLASSIC_OUT = sorted(n for n, _, c in ROWS if c == NONE)
DL_IN = sorted(n for n, d, _ in ROWS if d != NONE)
DL_OUT = sorted(n for n, d, _ in ROWS if d == NONE)


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    name = "OpenSees.exe" if os.name == "nt" else "OpenSees"
    for cand in (REPO / "dist" / "bin" / name,           # build.bat output
                 REPO / "build" / "Release" / name):     # Zone-A CI build tree
        if cand.exists():
            return str(cand)
    return None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(
    TCL_EXE is None,
    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)",
)


def _require_fork_module():
    # Decided by WHERE the module came from, never by a table row: a broken Python
    # hook that dropped `ladrunoBuild` must FAIL here, not skip (WP-168 review F1).
    if ops.__name__.startswith("openseespy"):
        pytest.skip("openseespy wheel fallback in use — fork build required")


def test_table_shape():
    names = [n for n, _, _ in ROWS]
    assert len(names) == len(set(names)), "duplicate command in the table"
    assert len(ROWS) >= 30, f"only {len(ROWS)} rows: the table lost commands"
    # the anchors every engine must carry (the card's named examples)
    for anchor in ("ladrunoBuild", "ladrunoThreads", "contactSurface", "contact", "ladrunoArcLength"):
        assert anchor in CLASSIC_IN and anchor in DL_IN, anchor


def test_no_ladruno_command_is_registered_outside_the_hook():
    rows, findings, scanned = check_ladruno_commands.run(REPO)
    assert len(rows) == len(ROWS)
    assert findings == [], "\n".join(findings)
    # non-vacuity (WP-168 review F2): upstream commands.cpp + TclWrapper.cpp +
    # PythonWrapper.cpp alone carry ~800 registrations; a REGISTER regex that
    # matched nothing would pass the assertion above silently.
    assert scanned >= 600, f"only {scanned} registrations scanned"


@needs_exe
def test_classic_tcl_registers_exactly_the_table_rows(tmp_path):
    deck = tmp_path / "commands.tcl"
    deck.write_text('puts "LADRUNO_CMDS:[join [lsort [info commands]] ,]"\n'
                    'puts "LADRUNO_BUILD:[ladrunoBuild]"\n', encoding="ascii")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: under pytest capture, an inherited stdin handle is invalid on Windows.
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(tmp_path), stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    line = next((ln for ln in out.splitlines() if ln.startswith("LADRUNO_CMDS:")), None)
    assert line is not None, f"deck did not print the command list (rc={proc.returncode}):\n{out}"
    have = set(line.split(":", 1)[1].split(","))
    assert "wipe" in have and "analyze" in have, "non-vacuity: the core commands must be listed"

    missing = [n for n in CLASSIC_IN if n not in have]
    assert not missing, f"table rows NOT registered in classic Tcl ({TCL_EXE}): {missing}"
    unexpected = [n for n in CLASSIC_OUT if n in have]
    assert not unexpected, (
        f"registered in classic Tcl but the table says LADRUNO_NONE: {unexpected} — "
        "update the row's classic column")

    stamp = next((ln for ln in out.splitlines() if ln.startswith("LADRUNO_BUILD:")), "")
    assert _HEX40.match(stamp.split(":", 1)[-1].strip()), f"ladrunoBuild smoke failed:\n{out}"


def test_python_registers_exactly_the_table_rows():
    _require_fork_module()
    missing = [n for n in DL_IN if not hasattr(ops, n)]
    assert not missing, f"table rows NOT registered in {ops.__file__}: {missing}"
    unexpected = [n for n in DL_OUT if hasattr(ops, n)]
    assert not unexpected, (
        f"registered in Python but the table says LADRUNO_NONE: {unexpected} — "
        "update the row's dl column")
    assert _HEX40.match(ops.ladrunoBuild())


def test_python_generated_bridge_propagates_errors_and_results():
    """The generated bridge keeps the hand-written one's contract: a command that
    fails (OPS_* < 0) raises, and a query returns its value."""
    _require_fork_module()
    n = ops.ladrunoThreads()
    assert isinstance(n, int) and n >= 1
    ops.wipe()
    with pytest.raises(Exception):
        ops.ladrunoArcLength("not-a-subcommand")   # no active LadrunoArcLength -> error
