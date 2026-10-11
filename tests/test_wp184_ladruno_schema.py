"""`ladrunoSchema` v0 (WP-184, apeGmsh#1667, program link F3).

`ladrunoSchema` returns one JSON string describing what this build offers, so
apeGmsh's BackendInfo can detect capabilities instead of keying on the build sha
alone. The contract this file pins:

* valid JSON; compact and key-sorted at every level (byte-identical to
  ``json.dumps(..., sort_keys=True, separators=(",", ":"))``), the same bytes on
  every call, and no filesystem paths;
* ``schema_version == 0`` and ``build == ladrunoBuild()``;
* ``commands`` is exactly this engine's rows of SRC/interpreter/LadrunoCommandTable.h,
  sorted (``engine`` says which: ``dl`` for Python/TclWrapper, ``classic`` for
  OpenSees.exe);
* ``options`` is GENERATED from SRC/material/LadrunoOptSpec.h: it equals the
  declarations parsed from the header, so a token removed from (or added to) a
  table changes the schema, and the pinned anchors below fail if one goes;
* every parser that refuses through ``ladruno_opt::reportUnknown`` has its
  (name, table) pair in ``kOptTables``, so a new fail-closed parser cannot be
  missing from the schema;
* any argument to ``ladrunoSchema`` fails closed (raises / Tcl error).
"""
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "SRC" / "interpreter" / "LadrunoCommandTable.h"
SPEC = REPO / "SRC" / "material" / "LadrunoOptSpec.h"
NONE = "LADRUNO_NONE"
ROW = re.compile(r'^\s*LADRUNO_COMMAND\s*\(\s*"([^"]+)"\s*,\s*(\w+)\s*,\s*(\w+)\s*\)', re.M)
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_PATHY = re.compile(r"[A-Za-z]:[\\/]|\\\\|/(?:Users|home|tmp|usr)/")

# Anchors: a token removed from kLadrunoRCOptions must fail this file even after a
# rebuild (where the header cross-check alone would agree with the new binary).
PINNED_RC = {
    "-Ce": {"nargs": "list", "since": "PR-155"},
    "-rho": {"nargs": 1, "since": "PR-155"},
    "-secant": {"nargs": 0, "since": "PR-192"},
    "-implexControl": {"nargs": 2, "since": "PR-263"},
    "-shearRetention": {"nargs": 1, "since": "ADR19-2b.2c.1"},
    "-autoRegularization": {"nargs": 1, "since": "PR-277"},
}
RC_FAMILY = ("LadrunoRCConcrete", "LadrunoRCFiniteStrain")


# ---- what the source declares ------------------------------------------------

def _table_rows():
    rows = ROW.findall(TABLE.read_text(encoding="utf-8"))
    assert rows, f"no LADRUNO_COMMAND rows parsed from {TABLE}"
    return rows


def _spec_tables():
    """{table symbol: {token: {"nargs", "since"}}} parsed from LadrunoOptSpec.h."""
    text = re.sub(r"//[^\n]*", "", SPEC.read_text(encoding="utf-8"))
    out = {}
    for sym, body in re.findall(r"static const OptSpec (k\w+)\[\]\s*=\s*\{(.*?)\};", text, re.S):
        ents = re.findall(r'\{"(-[A-Za-z]\w*)",\s*(LIST|\d+),\s*"([^"]+)"\}', body)
        assert ents, f"no entries parsed from {sym}"
        names = [n for n, _, _ in ents]
        assert len(names) == len(set(names)), f"duplicate token in {sym}: {names}"
        out[sym] = {n: {"nargs": "list" if k == "LIST" else int(k), "since": s}
                    for n, k, s in ents}
    assert out, f"no OptSpec tables parsed from {SPEC}"
    return out


def _spec_registry():
    """[(owner, table symbol)] parsed from kOptTables."""
    text = re.sub(r"//[^\n]*", "", SPEC.read_text(encoding="utf-8"))
    body = text.split("kOptTables[] = {", 1)[1].split("};", 1)[0]
    rows = re.findall(r'\{"(\w+)",\s*(k\w+),\s*count\((k\w+)\)\}', body)
    assert rows, "no kOptTables rows parsed"
    for owner, sym, counted in rows:
        assert sym == counted, f"kOptTables row {owner}: counts {counted} but points at {sym}"
    return [(o, s) for o, s, _ in rows]


def _expected_options():
    tables = _spec_tables()
    return {owner: tables[sym] for owner, sym in _spec_registry()}


def _check_shape(text, engine, commands):
    doc = json.loads(text)
    assert text == json.dumps(doc, sort_keys=True, separators=(",", ":")), \
        "schema is not compact + key-sorted (not deterministic)"
    assert not _PATHY.search(text), f"schema leaks a path: {_PATHY.search(text).group(0)!r}"
    assert set(doc) == {"build", "commands", "engine", "options", "schema_version"}
    assert doc["schema_version"] == 0
    assert doc["engine"] == engine
    assert _HEX40.match(doc["build"]), doc["build"]
    assert doc["commands"] == sorted(doc["commands"])
    assert doc["commands"] == commands, (
        f"commands differ from the table: missing {sorted(set(commands) - set(doc['commands']))}, "
        f"extra {sorted(set(doc['commands']) - set(commands))}")
    assert "ladrunoBuild" in doc["commands"] and "ladrunoSchema" in doc["commands"]

    opts = doc["options"]
    assert opts == _expected_options(), "options are not the LadrunoOptSpec.h declarations"
    for mat in RC_FAMILY:
        for tok, want in PINNED_RC.items():
            assert opts[mat].get(tok) == want, f"{mat} {tok}: {opts[mat].get(tok)} != {want}"
        # refused until fork PR #877 declares them (WP-167)
        assert "-crackedNu" not in opts[mat] and "-betaC" not in opts[mat]
    return doc


# ---- source-only checks (no binary) -------------------------------------------

def test_every_failclosed_parser_is_in_the_registry():
    registry = set(_spec_registry())
    call = re.compile(r'ladruno_opt::reportUnknown\(\s*"(\w+)"\s*,[^;]*?ladruno_opt::(k\w+)\s*\)')
    found = set()
    for path in (REPO / "SRC").rglob("*.cpp"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "ladruno_opt::reportUnknown" in text:
            found |= set(call.findall(text))
    assert found, "non-vacuity: no reportUnknown call sites found"
    missing = sorted(found - registry)
    assert not missing, (
        f"parsers refuse through these tables but kOptTables has no row, so "
        f"ladrunoSchema does not advertise them: {missing}")
    assert set(RC_FAMILY) <= {o for o, _ in registry}


# ---- Python (DL engine) ---------------------------------------------------------

def _require_fork_module():
    if ops.__name__.startswith("openseespy"):
        pytest.skip("openseespy wheel fallback in use — fork build required")


def test_python_schema():
    _require_fork_module()
    text = ops.ladrunoSchema()
    assert isinstance(text, str)
    dl_in = sorted(n for n, d, _ in _table_rows() if d != NONE)
    doc = _check_shape(text, "dl", dl_in)
    assert doc["build"] == ops.ladrunoBuild()
    assert ops.ladrunoSchema() == text, "two calls differ"


@pytest.mark.parametrize("bad", [("v1",), (1,), ("-commands",)])
def test_python_schema_refuses_arguments(bad):
    _require_fork_module()
    with pytest.raises(Exception):
        ops.ladrunoSchema(*bad)


# ---- classic Tcl ----------------------------------------------------------------

def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    name = "OpenSees.exe" if os.name == "nt" else "OpenSees"
    for cand in (REPO / "dist" / "bin" / name, REPO / "build" / "Release" / name):
        if cand.exists():
            return str(cand)
    return None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_schema(tmp_path):
    deck = tmp_path / "schema.tcl"
    deck.write_text(
        'puts "SCHEMA:[ladrunoSchema]"\n'
        'puts "BUILD:[ladrunoBuild]"\n'
        'puts "BADARG:[catch {ladrunoSchema v1}]"\n'
        'puts "DONE"\n', encoding="ascii")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True, timeout=300,
                          env=env, cwd=str(tmp_path), stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    lines = {ln.split(":", 1)[0]: ln.split(":", 1)[1] for ln in out.splitlines()
             if ln.split(":", 1)[0] in ("SCHEMA", "BUILD", "BADARG")}
    assert "DONE" in out and "SCHEMA" in lines, f"deck did not run to the end:\n{out}"
    classic_in = sorted(n for n, _, c in _table_rows() if c != NONE)
    doc = _check_shape(lines["SCHEMA"], "classic", classic_in)
    assert doc["build"] == lines["BUILD"].strip()
    assert lines["BADARG"].strip() == "1", "ladrunoSchema accepted an argument in classic Tcl"
