#!/usr/bin/env python3
"""Command-hook gate (WP-168, apeGmsh#1490). Dependency-free.

Every fork-only interpreter command is a row of SRC/interpreter/LadrunoCommandTable.h
and is registered by the one Ladruno_registerCommands hook each engine calls:
classic Tcl (SRC/tcl/commands.cpp), the DL Tcl engine (SRC/interpreter/TclWrapper.cpp)
and Python (SRC/interpreter/PythonWrapper.cpp). Before WP-168 each command was
registered by hand in all three files, and a command could silently miss an engine
(contact verbs missing from classic Tcl until ADR-78 P0.5; ladrunoDR/ladrunoArcLength
until #729; the ADR-44 modal family until its classic bridge).

  H1  a command registered outside the hook: a `Tcl_CreateCommand`,
      `Tcl_CreateObjCommand` or `addCommand` call with a string-literal name that
      is a table row, contains "ladruno" (any case), or carries a `Ladruno` comment
      on its statement. Waive one registration (same statement, or the line above):
          // ladruno-command-hook-ok <reason>
  H2  an engine file that does not call `Ladruno_registerCommands(` exactly once.
  H3  a malformed table: a row that does not parse, a duplicate name, or a row
      registered in no engine (both columns LADRUNO_NONE).

Run: python ci/check_ladruno_commands.py   (exit 1 on any finding)
Self-test: pytest -q ci/test_check_ladruno_commands.py; tests/test_ladruno_command_registry.py
runs it on the real tree in Zone-A.
"""
import argparse
import re
import sys
from pathlib import Path

TABLE = "SRC/interpreter/LadrunoCommandTable.h"
HOOK_FILES = {
    TABLE,
    "SRC/interpreter/LadrunoCommandsTclWrapper.h",
    "SRC/interpreter/LadrunoCommandsPython.h",
    "SRC/tcl/LadrunoCommandsClassicTcl.h",
}
ENGINES = {
    "SRC/tcl/commands.cpp": "LadrunoCommandsClassicTcl.h",
    "SRC/interpreter/TclWrapper.cpp": "LadrunoCommandsTclWrapper.h",
    "SRC/interpreter/PythonWrapper.cpp": "LadrunoCommandsPython.h",
}
NONE = "LADRUNO_NONE"

ROW = re.compile(r'^\s*LADRUNO_COMMAND\s*\(\s*"([^"]+)"\s*,\s*(\w+)\s*,\s*(\w+)\s*\)')
ROW_START = re.compile(r"^\s*LADRUNO_COMMAND\s*\(")
REGISTER = re.compile(
    r'\b(?:Tcl_CreateCommand|Tcl_CreateObjCommand)\s*\(\s*[\w>.\-]+\s*,\s*"([^"]+)"'
    r'|\baddCommand\s*\(\s*(?:[\w>.\-]+\s*,\s*)?"([^"]+)"')
HOOK_CALL = re.compile(r"\bLadruno_registerCommands\s*\(")
WAIVER = re.compile(r"//\s*ladruno-command-hook-ok\b\s*(.*)$")


def parse_table(text):
    """[(name, dl, classic, line_no)] plus H3 findings (as strings without a path)."""
    rows, problems = [], []
    for i, line in enumerate(text.splitlines(), 1):
        if not ROW_START.match(line):
            continue
        m = ROW.match(line)
        if not m:
            problems.append(f"{i}: row does not parse as LADRUNO_COMMAND(\"name\", dl, classic)")
            continue
        rows.append((m.group(1), m.group(2), m.group(3), i))
    seen = {}
    for name, dl, classic, i in rows:
        if name in seen:
            problems.append(f"{i}: duplicate command '{name}' (first at line {seen[name]})")
        seen.setdefault(name, i)
        if dl == NONE and classic == NONE:
            problems.append(f"{i}: '{name}' is registered in no engine (both columns {NONE})")
    if not rows:
        problems.append("no LADRUNO_COMMAND rows found")
    return rows, problems


def load_table(root):
    path = root / TABLE
    if not path.exists():
        return [], [f"H3 {TABLE}: command table not found"]
    rows, problems = parse_table(path.read_text(encoding="utf-8", errors="replace"))
    return rows, [f"H3 {TABLE}:{p}" for p in problems]


def _code_and_comment(line):
    """Split a line at the first // that is not inside a string literal."""
    in_str = False
    i = 0
    while i < len(line):
        c = line[i]
        if c == "\\" and in_str:
            i += 2
            continue
        if c == '"':
            in_str = not in_str
        elif not in_str and line.startswith("//", i):
            return line[:i], line[i:]
        i += 1
    return line, ""


def _blank_block_comments(text):
    """Replace /* ... */ comments with spaces (newlines kept), so commented-out code is not scanned."""
    return re.sub(r"/\*.*?\*/", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)


def check_registrations(root, rel, names):
    """H1: a fork command registered outside the hook."""
    findings = []
    src = root / "SRC"
    if not src.exists():
        return findings
    for path in sorted(src.rglob("*")):
        if path.suffix not in (".cpp", ".cc", ".c", ".h", ".hpp") or not path.is_file():
            continue
        r = rel(path)
        if r in HOOK_FILES:
            continue
        try:
            raw = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if "Command" not in raw:
            continue
        lines = _blank_block_comments(raw).splitlines()
        for i, line in enumerate(lines):
            code, comment = _code_and_comment(line)
            m = REGISTER.search(code)
            if not m:
                continue
            name = m.group(1) or m.group(2)
            # the statement's comments: this line until the `;`, at most 3 lines
            span = [comment]
            j = i
            while ";" not in _code_and_comment(lines[j])[0] and j + 1 < len(lines) and j < i + 3:
                j += 1
                span.append(_code_and_comment(lines[j])[1])
            above = _code_and_comment(lines[i - 1])[1] if i > 0 else ""
            if any(WAIVER.search(c) for c in span + [above]):
                continue
            why = None
            if name in names:
                why = "is a LadrunoCommandTable.h row"
            elif "ladruno" in name.lower():
                why = "is a Ladruno command name"
            elif any("Ladruno" in c for c in span):
                why = "carries a `// Ladruno` mark"
            if why:
                findings.append(
                    f"H1 {r}:{i + 1}: command '{name}' {why} but is registered by hand; add a row to "
                    f"{TABLE} instead (or waive: // ladruno-command-hook-ok <reason>)")
    return findings


def check_engines(root):
    """H2: each engine calls the hook exactly once and includes its header."""
    findings = []
    for r, header in ENGINES.items():
        path = root / r
        if not path.exists():
            findings.append(f"H2 {r}: engine file not found")
            continue
        code = [_code_and_comment(line)[0]
                for line in _blank_block_comments(path.read_text(encoding="utf-8", errors="replace")).splitlines()]
        calls = [i + 1 for i, c in enumerate(code) if HOOK_CALL.search(c)]
        if len(calls) != 1:
            where = ", ".join(map(str, calls)) or "none"
            findings.append(f"H2 {r}: Ladruno_registerCommands( called {len(calls)} times (lines: {where}); "
                            "expected exactly once")
        if not any(re.search(r'#\s*include\s*[<"]' + re.escape(header) + r'[>"]', c) for c in code):
            findings.append(f"H2 {r}: does not #include \"{header}\"")
    return findings


def run(root):
    root = root.resolve()

    def rel(p):
        return p.resolve().relative_to(root).as_posix()

    rows, findings = load_table(root)
    findings += check_engines(root)
    findings += check_registrations(root, rel, {n for n, _, _, _ in rows})
    return rows, findings


def main(argv=None):
    ap = argparse.ArgumentParser(description="Command-hook gate (WP-168).")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    args = ap.parse_args(argv)
    rows, findings = run(args.root)
    for f in findings:
        print(f)
    print(f"check_ladruno_commands: {len(rows)} table row(s), {len(findings)} finding(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
