"""`ladrunoBuild` -> machine-readable engine build stamp (git commit hash).

Born from the TIMs T1 incident (2026-08-10): a tangent-symmetry probe silently
ran against the WRONG engine build because the only provenance signal was the
splash banner's "Ladruno OpenSees build: <hash>" line, which (a) disappears
when the banner is suppressed, (b) breaks cp1252 text-mode subprocess capture
via its UTF-8 box glyphs, and (c) is not queryable in-process.

`ops.ladrunoBuild()` returns the SAME compile-time constant the banner prints
(CMake `git log -1 --format=%H` -> the `OPENSEES_VERSION` define), so equality
with the banner is structural — the gates here pin the contract that made the
incident possible: a 40-char hex hash, retrievable in-process, that survives
`LADRUNO_OPENSEES_QUIET=1` plus a cp1252 text-mode capture.
"""

import os
import re
import subprocess
import sys

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

_HEX40 = re.compile(r"^[0-9a-f]{40}$")


def _require_fork():
    if not hasattr(ops, "ladrunoBuild"):
        pytest.skip("openseespy wheel fallback in use — fork build required")


def test_ladrunoBuild_returns_40_char_hex_hash():
    _require_fork()
    stamp = ops.ladrunoBuild()
    assert isinstance(stamp, str)
    assert _HEX40.match(stamp), (
        f"ladrunoBuild() returned {stamp!r}, not a 40-char lowercase hex git hash "
        "(a build without git stamps 'unknown' — CI/dev builds must not)"
    )


def test_ladrunoBuild_survives_quiet_banner_and_cp1252_capture():
    """The incident scenario, end to end: banner suppressed, subprocess output
    captured in text mode with a cp1252 console encoding. The stamp must come
    through (pure ASCII) and match this process's engine."""
    _require_fork()
    mod_dir = os.path.dirname(os.path.abspath(ops.__file__))
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    env["PYTHONPATH"] = mod_dir + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONIOENCODING"] = "cp1252"
    code = (
        "import opensees as ops, os\n"
        "print('MODFILE=' + os.path.abspath(ops.__file__))\n"
        "print('STAMP=' + ops.ladrunoBuild())\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        stdin=subprocess.DEVNULL,  # parent stdin may be closed (headless/CI)
        capture_output=True, text=True, encoding="cp1252", errors="strict",
        env=env, timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    lines = dict(
        line.split("=", 1) for line in proc.stdout.splitlines() if "=" in line
    )
    # Same module file => same OPENSEES_VERSION constant => same banner line.
    assert os.path.dirname(lines["MODFILE"]) == mod_dir, (
        "subprocess imported a different opensees module — PYTHONPATH hijacked "
        "(see the installed-.pth quirk in LEDGER_quirks)"
    )
    assert lines["STAMP"] == ops.ladrunoBuild()
