"""`ladrunoSchema` v1 (WP-185, apeGmsh#1667, program link F3; v1: apeGmsh#2517).

`ladrunoSchema` returns one JSON string describing what this build offers, so
apeGmsh's BackendInfo can detect capabilities instead of keying on the build sha
alone. The contract this file pins:

* valid JSON; compact and key-sorted at every level (byte-identical to
  ``json.dumps(..., sort_keys=True, separators=(",", ":"))``), the same bytes on
  every call, and no filesystem paths;
* ``schema_version == 1`` and ``build == ladrunoBuild()``; the top-level keys
  are exactly build, capabilities, commands, engine, options, schema_version;
* ``capabilities`` is GENERATED from ``kCapabilities`` in
  SRC/interpreter/LadrunoSchema.h, sorted, and keeps the pinned anchors below;
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
SCHEMA_H = REPO / "SRC" / "interpreter" / "LadrunoSchema.h"
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

# Anchors: a capability is a promise apeGmsh gates on; dropping it must fail here
# even after a rebuild. Each one's behaviour test is named in kCapabilities.
PINNED_CAPABILITIES = ("mdStageFlagPreserved",)   # apeGmsh#2517 / #2444


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


def _header_capabilities():
    """The string rows of kCapabilities, parsed from LadrunoSchema.h."""
    text = re.sub(r"//[^
]*", "", SCHEMA_H.read_text(encoding="utf-8"))
    body = text.split("kCapabilities[] = {", 1)[1].split("};", 1)[0]
    caps = re.findall(r'"([^"]+)"', body)
    assert caps, "no kCapabilities rows parsed"
    assert len(caps) == len(set(caps)), f"duplicate capability in kCapabilities: {caps}"
    return caps


def _expected_options():
    tables = _spec_tables()
    return {owner: tables[sym] for owner, sym in _spec_registry()}


def _check_shape(text, engine, commands):
    doc = json.loads(text)
    assert text == json.dumps(doc, sort_keys=True, separators=(",", ":")), \
        "schema is not compact + key-sorted (not deterministic)"
    assert not _PATHY.search(text), f"schema leaks a path: {_PATHY.search(text).group(0)!r}"
    assert set(doc) == {"build", "capabilities", "commands", "engine", "options",
                        "schema_version"}
    assert doc["schema_version"] == 1
    assert doc["capabilities"] == sorted(_header_capabilities()), (
        "capabilities are not the kCapabilities declarations")
    for cap in PINNED_CAPABILITIES:
        assert cap in doc["capabilities"], f"pinned capability {cap!r} missing"
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

_SCAN_SUFFIXES = (".cpp", ".cc", ".cxx", ".c", ".h", ".hpp", ".hh")
_COMMENT = re.compile(r'//[^\n]*|/\*.*?\*/|("(?:\\.|[^"\\])*")', re.S)
_CALL = re.compile(r"\breportUnknown\s*\(")
_DEFN = re.compile(r"\bvoid\s+$")    # `void reportUnknown(` is the definition, not a call


def _strip_comments(text):
    # keep string literals (group 1), drop comments; keep newlines for line numbers
    return _COMMENT.sub(lambda m: m.group(1) or re.sub(r"[^\n]", " ", m.group(0)), text)


def _split_args(text, start):
    """Top-level comma-split of the argument list opening at text[start] == '('."""
    depth, args, cur, i = 0, [], [], start
    while i < len(text):
        c = text[i]
        if c == '"':                                   # skip a string literal whole
            j = i + 1
            while j < len(text) and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            cur.append(text[i:j + 1]); i = j + 1; continue
        if c in "([{":
            depth += 1
            if depth > 1: cur.append(c)
        elif c in ")]}":
            depth -= 1
            if depth == 0:
                args.append("".join(cur).strip()); return args
            cur.append(c)
        elif c == "," and depth == 1:
            args.append("".join(cur).strip()); cur = []
        else:
            cur.append(c)
        i += 1
    return None


def registry_findings(sources, registry, declared):
    """Every reportUnknown( call in `sources` ({name: text}) must parse into a
    (string-literal owner, table declared in LadrunoOptSpec.h) pair that has a
    kOptTables row. Returns (findings, number of calls seen)."""
    findings, seen = [], 0
    for name, raw in sources.items():
        text = _strip_comments(raw)
        for m in _CALL.finditer(text):
            if _DEFN.search(text[max(0, m.start() - 40):m.start()]):
                continue                                # the definition itself
            seen += 1
            where = f"{name}:{text.count(chr(10), 0, m.start()) + 1}"
            args = _split_args(text, m.end() - 1)
            if args is None or len(args) != 4:
                findings.append(f"{where}: reportUnknown call does not parse into 4 arguments: {args}")
                continue
            owner = re.fullmatch(r'"(\w+)"', args[0])
            table = re.fullmatch(r"(?:::)?(?:ladruno_opt\s*::\s*)?(k\w+)", args[3])
            if not owner:
                findings.append(f"{where}: owner {args[0]!r} is not a string literal")
            elif not table or table.group(1) not in declared:
                findings.append(f"{where}: table {args[3]!r} is not an OptSpec table declared in "
                                "SRC/material/LadrunoOptSpec.h (tables must live there)")
            elif (owner.group(1), table.group(1)) not in registry:
                findings.append(f"{where}: ({owner.group(1)}, {table.group(1)}) has no kOptTables row, "
                                "so ladrunoSchema does not advertise it")
    return findings, seen


def _src_sources():
    out = {}
    for path in (REPO / "SRC").rglob("*"):
        if path.suffix.lower() in _SCAN_SUFFIXES and path.is_file():
            text = path.read_text(encoding="utf-8", errors="replace")
            if "reportUnknown" in text:
                out[str(path.relative_to(REPO))] = text
    return out


def test_every_failclosed_parser_is_in_the_registry():
    registry, declared = set(_spec_registry()), set(_spec_tables())
    findings, seen = registry_findings(_src_sources(), registry, declared)
    assert seen >= 2, f"non-vacuity: only {seen} reportUnknown call sites found"
    assert not findings, "\n".join(findings)
    assert set(RC_FAMILY) <= {o for o, _ in registry}


@pytest.mark.parametrize("snippet", [
    # the header's documented usage: unprefixed, .cpp-local table
    'if (s == 0) { ladruno_opt::reportUnknown("MyMat", tag, opt, kOpts); return 0; }',
    # using-namespace form, unregistered owner on a declared table
    'using namespace ladruno_opt;\nif (!s) { reportUnknown("MyMat", tag, opt, kLadrunoRCOptions); }',
    # non-literal owner
    'ladruno_opt::reportUnknown(name, tag, opt, ladruno_opt::kLadrunoRCOptions);',
    # call split across lines, in a header
    'ladruno_opt::reportUnknown(\n  "Other",\n  tag, opt,\n  ladruno_opt::kLadrunoRCOptions);',
])
def test_registry_guard_flags_unregistered_calls(snippet):
    registry, declared = set(_spec_registry()), set(_spec_tables())
    findings, seen = registry_findings({"synthetic.h": snippet}, registry, declared)
    assert seen == 1 and len(findings) == 1, (seen, findings)


def test_registry_guard_ignores_comments_and_the_definition():
    registry, declared = set(_spec_registry()), set(_spec_tables())
    text = ('// ladruno_opt::reportUnknown("MyMat", tag, opt, kOpts);\n'
            '/* reportUnknown("X", t, o, kY); */\n'
            'inline void reportUnknown(const char* material, int tag, const char* tok,\n'
            '                          const OptSpec (&table)[N]) {}\n'
            'ladruno_opt::reportUnknown("LadrunoRCConcrete", tag, opt, ladruno_opt::kLadrunoRCOptions);\n')
    findings, seen = registry_findings({"ok.h": text}, registry, declared)
    assert (seen, findings) == (1, [])

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
