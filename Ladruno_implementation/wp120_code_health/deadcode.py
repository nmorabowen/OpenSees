#!/usr/bin/env python3
"""WP-120 step 3: dead-code survey, fork-authored vs vanilla (read-only).

 (a) unbuilt files: a .cpp named in no CMakeLists.txt / *.cmake; a .h #included
     by no source file and named in no CMake file.
 (b) unreferenced functions (fork files only): a function defined in a fork file
     whose bare name occurs nowhere in vanilla SRC, even in comments (so it is not a
     virtual override or a vanilla-called hook), and which has ZERO uses in fork code
     other than its own definition header(s) and prototype declarations.
     Candidates are hand-checked before the report calls anything dead.
 (c) always-false guards: `const|constexpr bool [Q::]X = false;` flags with the
     `if (X` / `X &&` uses they guard (counted over every file of the same group that
     names X); `#define X 0` flags used in `#if X`/`if (X)`; `if (0)`/`if (false)`;
     `#if 0` block lines; `#ifdef X` / `#if defined(X)` where X is #defined nowhere
     in the repo and is not named in any CMake file (dead unless a user passes -DX).
 (d) commented-out statements: `// call(...);` or `// x = ...;` lines, per KLOC.

    python deadcode.py [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, clean, functions, read, rel, split_fork_vanilla  # noqa: E402

IDENT = re.compile(r"[A-Za-z_]\w*")
SKIP_DIRS = ("build/", "mumps-", "dist/", ".claude/", "Win32/", "Win64/", "OTHER/")
INC = re.compile(r'#\s*include\s*[<"](?:[^>"]*/)?([^/>"]+)[>"]')
FALSE_CONST = re.compile(r"\b(?:const|constexpr)\s+bool\s+(?:\w+\s*::\s*)*(\w+)\s*=\s*false\s*;")
DEFINE0 = re.compile(r"^\s*#\s*define\s+(\w+)\s+0\s*$", re.M)
IF_LIT = re.compile(r"\bif\s*\(\s*(?:0|false)\s*\)")
IFDEF = re.compile(r"^\s*#\s*(?:ifdef\s+(\w+)|if\s+defined\s*\(?\s*(\w+))", re.M)
IF_ZERO = re.compile(r"^\s*#\s*if\s+0\b")
COMMENTED_STMT = re.compile(r"^\s*//\s*[\w.\->\[\]*]+\s*(?:\(.*\)|[-+*/]?=\s*[^=].*)\s*;\s*$")
PROTO_TAIL = re.compile(r"^\s*\([^;{}]*\)\s*(?:const\s*)?(?:noexcept\s*)?(?:override\s*)?(?:final\s*)?"
                        r"(?:=\s*(?:0|default|delete)\s*)?;")
TYPED_PREV = re.compile(r"(?:[A-Za-z_]\w*|[*&]|[^-]>)$")
NOT_TYPE_PREV = re.compile(r"(?:\b(?:return|else|new|delete|throw|case|do)|->)$")


def cmake_text():
    parts = []
    for p in ROOT.rglob("*"):
        r = p.relative_to(ROOT).as_posix()
        if r.startswith(SKIP_DIRS) or not p.is_file():
            continue
        if p.name == "CMakeLists.txt" or p.suffix == ".cmake":
            parts.append(read(p))
    return "\n".join(parts)


def unbuilt(fork, van, cm, includes):
    names = Counter(re.findall(r"[\w.+-]+\.(?:cpp|cc|cxx|c|h|hpp)\b", cm))
    out = {"fork": {"cpp": [], "h": []}, "vanilla": {"cpp": [], "h": []}}
    for group, files in (("fork", fork), ("vanilla", van)):
        for p in files:
            if p.suffix in (".cpp", ".cc", ".cxx"):
                if names[p.name] == 0:
                    out[group]["cpp"].append(rel(p))
            elif names[p.name] == 0 and includes[p.name] == 0:
                out[group]["h"].append(rel(p))
    return out


def if0_lines(raw_lines):
    n, depth = 0, 0
    for ln in raw_lines:
        s = ln.strip()
        if depth:
            n += 1
            if re.match(r"#\s*if", s):
                depth += 1
            elif re.match(r"#\s*endif", s):
                depth -= 1
        elif IF_ZERO.match(ln):
            depth, n = 1, n + 1
    return n


def guards(files, cleaned, defined_names):
    res = {"false_consts": [], "define0": [], "if_literal": 0, "if0_blocks_lines": 0,
           "undefined_ifdef": [], "commented_stmts": 0, "raw_lines": 0,
           "per_file_commented": []}
    # a flag is scoped to its class (`const bool Q::X = false;` -> files naming Q) or,
    # unqualified, to its own file and same-stem companion; type traits are not guards
    for p in files:
        for m in FALSE_CONST.finditer(cleaned[p]):
            name = m.group(1)
            if name == "value":
                continue
            q = re.search(r"(\w+)\s*::\s*" + re.escape(name) + r"\s*=", m.group(0))
            if q:
                scope = [x for x in files if q.group(1) in cleaned[x] and name in cleaned[x]]
            else:
                scope = [x for x in files if x.stem == p.stem and x.parent == p.parent]
            pat = re.compile(r"\bif\s*\(\s*" + re.escape(name) + r"\b|\b" + re.escape(name) + r"\s*&&")
            uses = sum(len(pat.findall(cleaned[x])) for x in scope)
            res["false_consts"].append((rel(p), (q.group(1) + "::" if q else "") + name, uses))
    undef = Counter()
    for p in files:
        text = read(p)
        raw = text.splitlines()
        res["raw_lines"] += len(raw)
        code = cleaned[p]
        for m in DEFINE0.finditer(text):
            name = m.group(1)
            uses = len(re.findall(r"#\s*if\s+" + re.escape(name) + r"\b", text)) + \
                len(re.findall(r"\bif\s*\(\s*" + re.escape(name) + r"\s*\)", code))
            if uses:
                res["define0"].append((rel(p), name, uses))
        res["if_literal"] += len(IF_LIT.findall(code))
        res["if0_blocks_lines"] += if0_lines(raw)
        for m in IFDEF.finditer(text):
            name = m.group(1) or m.group(2)
            if name and name not in defined_names and not name.startswith("_") and \
                    not re.match(r"(?:WIN|APPLE|linux|unix|MSC|GNUC|clang)", name):
                undef[(rel(p), name)] += 1
        n = sum(1 for ln in raw if COMMENTED_STMT.match(ln))
        if n:
            res["commented_stmts"] += n
            res["per_file_commented"].append((rel(p), n))
    res["per_file_commented"].sort(key=lambda r: -r[1])
    res["undefined_ifdef"] = sorted(((f, n, c) for (f, n), c in undef.items()), key=lambda r: -r[2])
    return res


def _uses(cl, name, def_lines):
    """Line numbers of occurrences of `name` in cleaned lines that are neither in a
    definition header nor a prototype declaration `T name(args) [const] [override];`."""
    uses, off = [], 0
    text = "\n".join(cl)
    line_of = []
    for li, ln in enumerate(cl):
        line_of.extend([li] * (len(ln) + 1))
    for m in re.finditer(r"\b" + re.escape(name) + r"\b", text):
        li = line_of[m.start()]
        if li in def_lines:
            continue
        prev = text[max(0, m.start() - 200):m.start()].rstrip()
        if TYPED_PREV.search(prev) and not NOT_TYPE_PREV.search(prev) and \
                PROTO_TAIL.match(text[m.end():m.end() + 400]):
            continue
        uses.append(li + 1)
    return uses


def unreferenced(fork, cleaned, van_ident):
    defs = defaultdict(list)
    for p in fork:
        cl = cleaned[p].split("\n")
        for f in functions(cl):
            parts = f.name.split("::")
            base = parts[-1]
            if base.startswith(("operator", "~", "OPS_")) or base == "main":
                continue
            if len(parts) >= 2 and parts[-1] == parts[-2]:
                continue                               # constructor
            defs[base].append((p, f.start, f.open, f.name))
    out = []
    for base, sites in defs.items():
        if van_ident[base]:
            continue                                   # override / vanilla-called / shared name
        uses = []
        for p in fork:
            if base not in cleaned[p]:
                continue
            dl = {li for sp, st, op, _ in sites if sp == p for li in range(st, op + 1)}
            uses += [(rel(p), u) for u in _uses(cleaned[p].split("\n"), base, dl)]
        if not uses:
            out.append({"name": base, "sites": [(rel(sp), st + 1, nm) for sp, st, _, nm in sites]})
    out.sort(key=lambda d: (d["sites"][0][0], d["sites"][0][1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json")
    a = ap.parse_args()
    fork, van = split_fork_vanilla()
    includes, van_ident, defined = Counter(), Counter(), set()
    cleaned = {}
    for p in fork + van:
        t = read(p)
        includes.update(INC.findall(t))
        defined.update(re.findall(r"#\s*define\s+(\w+)", t))
        cleaned[p] = "\n".join(clean(t))
    for p in van:
        van_ident.update(IDENT.findall(read(p)))     # raw text: comments count (conservative)
    cm = cmake_text()
    defined.update(re.findall(r"\b([A-Za-z_]\w{2,})\b", cm))   # any name CMake mentions may be -D'd
    res = {
        "n_fork": len(fork), "n_vanilla": len(van),
        "unbuilt": unbuilt(fork, van, cm, includes),
        "guards_fork": guards(fork, cleaned, defined),
        "guards_vanilla": guards(van, cleaned, defined),
        "unreferenced_fork": unreferenced(fork, cleaned, van_ident),
    }
    u = res["unbuilt"]
    print(f"unbuilt: fork .cpp={len(u['fork']['cpp'])} .h={len(u['fork']['h'])} | "
          f"vanilla .cpp={len(u['vanilla']['cpp'])} .h={len(u['vanilla']['h'])}")
    for x in u["fork"]["cpp"] + u["fork"]["h"]:
        print("   ", x)
    for g in ("guards_fork", "guards_vanilla"):
        r = res[g]
        kloc = r["raw_lines"] / 1000.0
        print(f"{g}: false-const flags={len(r['false_consts'])} guarding {sum(x[2] for x in r['false_consts'])} uses; "
              f"#define-0 flags used={len(r['define0'])}; if(0/false)={r['if_literal']}; "
              f"#if0 lines={r['if0_blocks_lines']}; undefined #ifdef={len(r['undefined_ifdef'])}; "
              f"commented-out stmts={r['commented_stmts']} ({r['commented_stmts'] / kloc:.2f}/KLOC)")
        for x in r["false_consts"]:
            print("   F", x)
        for x in r["define0"][:20]:
            print("   0", x)
        for x in r["undefined_ifdef"][:40 if g == "guards_fork" else 10]:
            print("   D", x)
        for x in r["per_file_commented"][:10]:
            print("   C", x)
    print(f"unreferenced fork functions: {len(res['unreferenced_fork'])}")
    for c in res["unreferenced_fork"]:
        print(f"   {c['name']}  {c['sites'][:3]}")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, default=list), encoding="utf-8")


if __name__ == "__main__":
    main()
