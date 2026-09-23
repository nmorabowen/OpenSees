#!/usr/bin/env python3
"""Quirk-pattern gate (WP-115). Dependency-free.

Turns LEDGER_quirks entries that name a greppable pattern into checks, so a
known trap fails CI instead of relying on someone re-reading the ledger.
Scans only fork-authored sources (files carrying the LADRUNO-HEADER-START
stamp); vanilla code is out of scope (vanilla-footprint rule).

  L1 rayleigh   In an element, adding getRayleighDampingForces() into a buffer
                that is NOT declared in the same function. betaK Rayleigh calls
                getTangentStiff(), which may refill that shared buffer and
                silently drop inertia (and -Q). LEDGER_quirks: "MUST snapshot
                the shared static `resid`" (#562 recurrence).
                Waive one site, on the statement line or the line above:
                    // ladruno-lint: rayleigh-ok <reason>

  L2 wipe       A process-wide singleton (a `static X &instance(...)`) whose
                state is not reset on `wipe`. Reset = a call
                `X::instance().reset*(...)` inside a wipe hook
                (Domain::clearAll, PartitionedDomain::clearAll, OPS_clearAll*)
                or inside a function a wipe hook calls. LEDGER_quirks:
                "`wipe()` does NOT recreate the Domain".
                Waive at the declaration (same line or up to two lines above):
                    // ladruno-lint: wipe-ok <reason>

  L3 pointers   Every `Quirks: "..."` pointer in .claude/skills/*/SKILL.md
                must still match text in LEDGER_quirks.md, so a renamed
                heading cannot silently orphan a guide item.

A waiver needs a reason of at least 12 characters.

Usage:
    python ci/check_quirk_patterns.py              # all checks, exit 1 on any finding
    python ci/check_quirk_patterns.py --only L1,L2
    python ci/check_quirk_patterns.py --root DIR   # scan another tree (e.g. a git archive)
    python ci/check_quirk_patterns.py --list-waivers
"""
import argparse
import re
import sys
from pathlib import Path

STAMP = "LADRUNO-HEADER-START"
MIN_REASON = 12

WAIVER = re.compile(r"//\s*ladruno-lint:\s*(rayleigh-ok|wipe-ok)\b(.*)$")
RAYLEIGH_CALL = re.compile(r"getRayleighDampingForces\s*\(\s*\)")
# const Vector &v = this->getRayleighDampingForces();
RAYLEIGH_BIND = re.compile(
    r"\bconst\s+Vector\s*&\s*(\w+)\s*=\s*(?:this\s*->\s*)?getRayleighDampingForces\s*\(\s*\)")
ACCUM_PLUS = re.compile(r"^\s*([A-Za-z_]\w*)\s*\+=\s*(.+)$")
ACCUM_ADDV = re.compile(r"^\s*([A-Za-z_]\w*)\s*\.\s*addVector\s*\((.+)$")
SINGLETON = re.compile(r"\bstatic\s+([A-Za-z_]\w*)\s*&\s*instance\s*\(")
FUNC_HEAD = re.compile(r"^(?:[A-Za-z_][\w:<>\*&,\s]*?[\s\*&])?(~?[A-Za-z_][\w:]*)\s*\(")
CONTROL = {"if", "for", "while", "switch", "return", "sizeof", "catch", "else", "do"}
WIPE_HOOK = re.compile(r"(?:^|::)(clearAll|OPS_clearAll\w*)$")


def strip_comment(line):
    """Drop a // comment (good enough for these patterns; no string literals involved)."""
    i = line.find("//")
    return line if i < 0 else line[:i]


def _sources(root, stamped_only):
    for p in sorted((root / "SRC").rglob("*")):
        if p.suffix not in (".cpp", ".h", ".hpp") or not p.is_file():
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if stamped_only and STAMP not in text:
            continue
        yield p, text.splitlines()


def fork_sources(root):
    return _sources(root, stamped_only=True)


def all_sources(root):
    return _sources(root, stamped_only=False)


def enclosing_function(lines, idx):
    """(name, header_index) of the column-0 function definition enclosing lines[idx]."""
    for j in range(idx, -1, -1):
        raw = lines[j]
        if not raw or raw[0] in " \t#/{}*":
            continue
        code = strip_comment(raw).rstrip()
        if code.endswith(";"):
            continue
        m = FUNC_HEAD.match(code)
        if m and m.group(1).split("::")[-1] not in CONTROL:
            return m.group(1), j
    return None, 0


def waiver_near(lines, idx, kind, above):
    """(True, reason) valid waiver, (False, reason) too-short reason, (None, "") none."""
    for j in range(idx, max(-1, idx - above - 1), -1):
        m = WAIVER.search(lines[j])
        if m and m.group(1) == kind:
            reason = m.group(2).strip()
            return len(reason) >= MIN_REASON, reason
    return None, ""


def declared_in(lines, start, end, name):
    decl = re.compile(r"\b(?:static\s+)?(?:const\s+)?Vector\s*&?\s*" + re.escape(name) + r"\s*[\(;=\{]")
    return any(decl.search(strip_comment(lines[k])) for k in range(start, end + 1))


def check_rayleigh(root, rel):
    findings = []
    for path, lines in fork_sources(root):
        bound = {}  # local ref name -> line index
        for i, raw in enumerate(lines):
            code = strip_comment(raw)
            m = RAYLEIGH_BIND.search(code)
            if m:
                bound[m.group(1)] = i
            for pat in (ACCUM_PLUS, ACCUM_ADDV):
                a = pat.match(code)
                if not a:
                    continue
                target, rhs = a.group(1), a.group(2)
                via_ref = any(re.search(r"\b" + re.escape(v) + r"\b", rhs) and 0 < i - k <= 3
                              for v, k in bound.items())
                if not (RAYLEIGH_CALL.search(rhs) or via_ref):
                    continue
                fname, head = enclosing_function(lines, i)
                if declared_in(lines, head, i, target):
                    continue
                ok, _ = waiver_near(lines, i, "rayleigh-ok", above=1)
                if ok:
                    continue
                why = ("rayleigh-ok waiver reason too short" if ok is False else
                       f"'{target}' is not declared in {fname or 'the enclosing function'}; "
                       "snapshot into a function-local vector before adding Rayleigh forces, "
                       "or waive with '// ladruno-lint: rayleigh-ok <reason>'")
                findings.append(f"L1 {rel(path)}:{i + 1}: {why}")
    return findings


def reset_calls(root):
    """{class: [function name]} for X::instance().reset*(...) calls anywhere in SRC."""
    call = re.compile(r"\b([A-Za-z_]\w*)::instance\s*\(\s*\)\s*(?:\.|->)\s*reset\w*\s*\(")
    out = {}
    for _, lines in all_sources(root):
        for i, raw in enumerate(lines):
            for m in call.finditer(strip_comment(raw)):
                fname, _ = enclosing_function(lines, i)
                out.setdefault(m.group(1), []).append(fname or "")
    return out


def wipe_hook_bodies(root):
    """Concatenated source text of every wipe-hook function body."""
    bodies = []
    for _, lines in all_sources(root):
        for i, raw in enumerate(lines):
            if not raw or raw[0] in " \t#/{}*":
                continue
            code = strip_comment(raw).rstrip()
            if code.endswith(";"):
                continue
            m = FUNC_HEAD.match(code)
            if not (m and WIPE_HOOK.search(m.group(1))):
                continue
            depth, started, body = 0, False, []
            for raw2 in lines[i:]:
                c = strip_comment(raw2)
                body.append(c)
                depth += c.count("{") - c.count("}")
                started = started or "{" in c
                if started and depth <= 0:
                    break
            bodies.append("\n".join(body))
    return "\n".join(bodies)


def check_wipe(root, rel):
    findings = []
    resets = reset_calls(root)
    hooks = wipe_hook_bodies(root)
    for path, lines in fork_sources(root):
        for i, raw in enumerate(lines):
            m = SINGLETON.search(strip_comment(raw))
            if not m:
                continue
            cls = m.group(1)
            ok, _ = waiver_near(lines, i, "wipe-ok", above=2)
            if ok:
                continue
            if ok is False:
                findings.append(f"L2 {rel(path)}:{i + 1}: wipe-ok waiver reason too short for '{cls}'")
                continue
            wired = False
            for fname in resets.get(cls, []):
                short = fname.split("::")[-1]
                if WIPE_HOOK.search(fname) or (short and re.search(r"\b" + re.escape(short) + r"\s*\(", hooks)):
                    wired = True
                    break
            if not wired:
                findings.append(
                    f"L2 {rel(path)}:{i + 1}: singleton '{cls}' is not reset on wipe; call "
                    f"{cls}::instance().reset...() from Domain::clearAll or an OPS_clearAll* hook, "
                    "or waive with '// ladruno-lint: wipe-ok <reason>'")
    return findings


def check_pointers(root, rel):
    findings = []
    ledger = root / "Ladruno_implementation" / "LEDGER_quirks.md"
    guides = sorted((root / ".claude" / "skills").glob("*/SKILL.md"))
    if not guides:
        return findings
    if not ledger.exists():
        return [f"L3 {rel(ledger)}: ledger not found"]
    flat = re.sub(r"\s+", " ", ledger.read_text(encoding="utf-8", errors="replace"))
    quote = re.compile(r'\s*(?:,|and)?\s*"([^"]+)"')
    for guide in guides:
        text = re.sub(r"\s+", " ", guide.read_text(encoding="utf-8", errors="replace"))
        for m in re.finditer(r"Quirks:", text):
            pos = m.end()
            while True:
                q = quote.match(text, pos)
                if not q:
                    break
                if q.group(1) not in flat:
                    findings.append(f'L3 {rel(guide)}: pointer not found in LEDGER_quirks.md: "{q.group(1)}"')
                pos = q.end()
    return findings


def list_waivers(root, rel):
    for path, lines in fork_sources(root):
        for i, raw in enumerate(lines):
            m = WAIVER.search(raw)
            if m:
                print(f"{rel(path)}:{i + 1}: {m.group(1)} {m.group(2).strip()}")


def main():
    ap = argparse.ArgumentParser(description="Quirk-pattern gate (WP-115).")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--only", default="L1,L2,L3", help="comma list of L1,L2,L3")
    ap.add_argument("--list-waivers", action="store_true")
    args = ap.parse_args()
    root = args.root.resolve()

    def rel(p):
        try:
            return p.resolve().relative_to(root).as_posix()
        except ValueError:
            return str(p)

    if args.list_waivers:
        list_waivers(root, rel)
        return 0
    wanted = {s.strip().upper() for s in args.only.split(",") if s.strip()}
    checks = [("L1", check_rayleigh), ("L2", check_wipe), ("L3", check_pointers)]
    findings = []
    for key, fn in checks:
        if key in wanted:
            findings += fn(root, rel)
    for f in findings:
        print(f)
    print(f"check_quirk_patterns: {len(findings)} finding(s) [{','.join(sorted(wanted))}]")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
