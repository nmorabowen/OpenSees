#!/usr/bin/env python3
"""Quirk-pattern gate (WP-115). Dependency-free.

Turns LEDGER_quirks entries that name a greppable pattern into checks, so a
known trap fails CI instead of relying on someone re-reading the ledger.
L1/L2 scan only fork-authored sources (files carrying the LADRUNO-HEADER-START
stamp); vanilla code is out of scope (vanilla-footprint rule).

The scanner is a small C++ tokenizer, not a line grep: comments, string
literals and `#if 0` blocks are blanked first; functions are found by brace
matching (so inline class methods and indented headers resolve correctly);
statements are split on `;` (so multi-line statements and one-line `if (...)`
bodies are seen whole).

  L1 rayleigh   The Rayleigh force (getRayleighDampingForces(), or a variable
                bound to it) is accumulated into a buffer that is not a
                function-local vector SEEDED BEFORE the first Rayleigh call in
                that function. betaK Rayleigh re-enters getTangentStiff(),
                which may refill a shared buffer and silently drop inertia and
                -Q; a snapshot taken after the call is too late.
                LEDGER_quirks: "MUST snapshot the shared static `resid`" (#562).
                Non-local = class/file static, member, `this->X`, a pointer
                target (`*p`, `p->`), or a local REFERENCE alias (`Vector &r`).
                Also catches `P = P + R`.
                Waive one statement (comment anywhere in its line span, or on
                the line above it):
                    // ladruno-lint: rayleigh-ok <reason>

  L2 wipe       A process-wide singleton (`static X &instance(...)`) whose state
                is not reset on `wipe`. Reset = `X::instance().reset*(...)`
                called directly in a wipe hook -- Domain::clearAll,
                PartitionedDomain::clearAll, or a free OPS_clearAll* function --
                or in a FREE function that a wipe hook calls.
                LEDGER_quirks: "`wipe()` does NOT recreate the Domain".
                Waive at the declaration (same line or up to two lines above):
                    // ladruno-lint: wipe-ok <reason>

  L3 pointers   Every `Quirks: "..."` pointer in .claude/skills/*/SKILL.md must
                still match text in LEDGER_quirks.md.

A waiver needs a reason of at least 12 characters, and a waiver that no longer
suppresses anything is itself a finding (stale).

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
SUFFIXES = (".cpp", ".h", ".hpp", ".cc", ".cxx")

WAIVER = re.compile(r"//\s*ladruno-lint:\s*(rayleigh-ok|wipe-ok)\b(.*)$")
RAYLEIGH = re.compile(r"(?:\bthis\s*->\s*)?\bgetRayleighDampingForces\s*\(\s*\)")
SINGLETON = re.compile(r"\bstatic\s+([A-Za-z_]\w*)\s*&\s*instance\s*\(")
RESET_CALL = re.compile(r"\b([A-Za-z_]\w*)::instance\s*\(\s*\)\s*(?:\.|->)\s*reset\w*\s*\(")
CONTROL = {"if", "for", "while", "switch", "catch", "return", "sizeof", "else", "do", "new", "delete"}
SCOPE_KW = re.compile(r"\b(namespace|class|struct|union|enum)\b")
HOOKS_QUALIFIED = {"Domain::clearAll", "PartitionedDomain::clearAll"}
HOOK_FREE = re.compile(r"^OPS_clearAll\w*$")


# --------------------------------------------------------------------------
# source cleaning: blank comments, string/char literals and `#if 0` blocks,
# keeping every newline so line numbers survive.
# --------------------------------------------------------------------------
def clean(text):
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        nxt = text[i + 1] if i + 1 < n else ""
        if c == "/" and nxt == "/":
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c == "/" and nxt == "*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append(c + " " * max(0, j - i - 2) + (c if j - i >= 2 else ""))
            i = j
        else:
            out.append(c)
            i += 1
    lines = "".join(out).split("\n")
    # blank `#if 0 ... #endif` (nesting-aware), and every other preprocessor line
    depth, res = 0, []
    for ln in lines:
        s = ln.strip()
        if depth:
            if re.match(r"#\s*if", s):
                depth += 1
            elif re.match(r"#\s*endif", s):
                depth -= 1
            res.append("")
            continue
        if re.match(r"#\s*if\s+0\b", s):
            depth = 1
            res.append("")
            continue
        res.append("" if s.startswith("#") else ln)
    return res


def _balanced_end(s, open_idx):
    depth = 0
    for k in range(open_idx, len(s)):
        if s[k] == "(":
            depth += 1
        elif s[k] == ")":
            depth -= 1
            if depth == 0:
                return k
    return -1


FUNC_NAME = re.compile(r"([~A-Za-z_][\w:~]*)\s*\(")


def _function_name(header):
    """Qualified name if `header` (text before a `{`) is a function definition."""
    h = " ".join(header.split())
    if not h or "=" in h.split("(")[0] or SCOPE_KW.search(h.split("(")[0]):
        return None
    for m in FUNC_NAME.finditer(h):
        name = m.group(1)
        if name.split("::")[-1] in CONTROL:
            return None
        end = _balanced_end(h, m.end() - 1)
        if end < 0:
            return None
        rest = h[end + 1:].strip()
        # allow qualifiers, and a constructor initializer list after ':'
        if re.fullmatch(r"(?:const|override|final|noexcept|volatile|&|\s)*", rest) or \
                re.match(r"(?:const|noexcept|\s)*:(?!:)", rest):
            return name
        return None
    return None


class Func:
    __slots__ = ("name", "start", "open", "end")

    def __init__(self, name, start, open_, end):
        self.name, self.start, self.open, self.end = name, start, open_, end


def functions(cl):
    """All function definitions in cleaned lines, by brace matching.
    Returns Func(name=qualified name, start=header line, open=line of '{', end)."""
    found, stack, header, header_line = [], [], [], None
    for li, line in enumerate(cl):
        for ch in line:
            if ch == "{":
                htxt = "".join(header)
                name = _function_name(htxt)
                if name is not None and not any(kind == "func" for kind, *_ in stack):
                    scopes = [nm for kind, nm, *_ in stack if kind == "type" and nm]
                    if scopes and "::" not in name:
                        name = "::".join(scopes + [name])
                    stack.append(("func", name, header_line if header_line is not None else li, li))
                else:
                    m = SCOPE_KW.search(htxt)
                    tname = None
                    if m and m.group(1) in ("class", "struct", "union"):
                        mm = re.search(r"\b(?:class|struct|union)\s+(?:\w+\s+)*?([A-Za-z_]\w*)\s*(?::[^{]*)?$",
                                       " ".join(htxt.split()))
                        tname = mm.group(1) if mm else None
                    stack.append(("type" if m else "block", tname, li, li))
                header, header_line = [], None
            elif ch == "}":
                if stack:
                    kind, nm, hl, ol = stack.pop()
                    if kind == "func":
                        found.append(Func(nm, hl, ol, li))
                header, header_line = [], None
            elif ch == ";":
                header, header_line = [], None
            else:
                if not ch.isspace() and header_line is None:
                    header_line = li
                header.append(ch)
        header.append(" ")
    return found


def statements(cl, f):
    """(text, first_line, last_line) for each `;`-terminated statement in function
    f's body, with block braces as boundaries."""
    out, buf, first = [], [], None
    for li in range(f.open, f.end + 1):
        line = cl[li]
        start = line.index("{") + 1 if li == f.open else 0
        for ch in line[start:]:
            if ch in ";{}":
                if ch == ";" and buf:
                    out.append(("".join(buf).strip(), first, li))
                elif ch == "{" and buf:
                    out.append(("".join(buf).strip() + " {", first, li))
                buf, first = [], None
            else:
                if first is None and not ch.isspace():
                    first = li
                buf.append(ch)
        buf.append(" ")
    return [(t, a, b) for t, a, b in out if t]


def strip_prefix(stmt):
    """Drop leading `if (...)`, `else`, `for (...)`, `while (...)` so a one-line
    body is analysed as the statement it is."""
    s = stmt.strip()
    while True:
        m = re.match(r"(?:else\b\s*)?(?:if|for|while)\s*\(", s)
        if m:
            end = _balanced_end(s, m.end() - 1)
            if end < 0:
                return s
            s = s[end + 1:].strip()
            continue
        if re.match(r"else\b", s):
            s = s[4:].strip()
            continue
        return s


# --------------------------------------------------------------------------
# file iteration
# --------------------------------------------------------------------------
def _sources(root, stamped_only, needles=None):
    """(path, raw_lines, cleaned_lines). `needles`: skip files containing none of them
    (cheap prefilter before the character-level clean)."""
    for p in sorted((root / "SRC").rglob("*")):
        if p.suffix not in SUFFIXES or not p.is_file():
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if stamped_only and STAMP not in text:
            continue
        if needles and not any(n in text for n in needles):
            continue
        yield p, text.splitlines(), clean(text)


def waiver_at(raw, idx, kind, above):
    """(line_index, reason) of a `kind` waiver on raw[idx] or up to `above` lines before."""
    for j in range(idx, max(-1, idx - above - 1), -1):
        m = WAIVER.search(raw[j])
        if m and m.group(1) == kind:
            return j, m.group(2).strip()
    return None, ""


# --------------------------------------------------------------------------
# L1
# --------------------------------------------------------------------------
TARGET = r"(?P<target>\(\s*\*\s*\w+\s*\)|\*\s*\w+|this\s*->\s*\w+|\w+\s*->\s*\w+|\w+)"
ACC_PLUS = re.compile(r"^" + TARGET + r"\s*\+=\s*(?P<rhs>.+)$", re.S)
ACC_ADDV = re.compile(r"^" + TARGET + r"\s*(?:\.|->)\s*addVector\s*\((?P<rhs>.+)$", re.S)
ACC_SELF = re.compile(r"^(?P<target>\w+)\s*=\s*(?P<rhs>.+)$", re.S)
BIND_ASSIGN = re.compile(r"(?:\b(?:const\s+)?(?:Vector|auto)\s*&?\s*)?\b(\w+)\s*(?<![+\-*/<>!=])=(?!=)\s*" + RAYLEIGH.pattern)
BIND_CTOR = re.compile(r"\b(?:const\s+)?(?:Vector|auto)\s+(\w+)\s*\(\s*" + RAYLEIGH.pattern + r"\s*\)")
LOCAL_DECL = r"\b(?:static\s+)?(?:const\s+)?Vector\s+{name}\s*(?P<init>[\(=\{{]|$)"
REF_DECL = r"\b(?:static\s+)?(?:const\s+)?Vector\s*&\s*{name}\b"


def _analyze_rayleigh(cl, f):
    """Yield (stmt_first_line, message) for L1 findings in function f."""
    stmts = statements(cl, f)
    ray_idx = [k for k, (s, *_) in enumerate(stmts) if RAYLEIGH.search(s)]
    if not ray_idx:
        return
    first_ray = ray_idx[0]
    bound = set()
    for k, (stmt, line, last) in enumerate(stmts):
        span = (line, last)
        body = strip_prefix(stmt)
        for rx in (BIND_ASSIGN, BIND_CTOR):
            m = rx.search(body)
            if m:
                bound.add(m.group(1))

        def carries_rayleigh(rhs):
            return bool(RAYLEIGH.search(rhs)) or any(re.search(r"\b" + re.escape(b) + r"\b", rhs) for b in bound)

        m = ACC_PLUS.match(body) or ACC_ADDV.match(body)
        target = None
        if m and carries_rayleigh(m.group("rhs")):
            target = m.group("target")
        else:
            s = ACC_SELF.match(body)
            if s and s.group("target") not in bound and carries_rayleigh(s.group("rhs")) and \
                    re.search(r"\b" + re.escape(s.group("target")) + r"\b", s.group("rhs")):
                target = s.group("target")             # P = P + R
        if target is None:
            continue
        t = " ".join(target.split())
        if "*" in t or "->" in t:
            yield span, (f"accumulates Rayleigh forces into '{t}' (member/pointer target) in {f.name}; "
                         "snapshot into a function-local vector first")
            continue
        name = t
        prior = "\n".join(s for s, *_ in stmts[:k])
        if re.search(REF_DECL.format(name=re.escape(name)), prior):
            yield span, (f"'{name}' is a local REFERENCE (an alias, not a snapshot) in {f.name}; "
                         "declare a function-local Vector instead")
            continue
        decl_k = next((j for j, (s, *_) in enumerate(stmts[:k])
                       if re.search(LOCAL_DECL.format(name=re.escape(name)), s)), None)
        if decl_k is None:
            yield span, (f"'{name}' is not declared in {f.name}; snapshot into a function-local "
                         "vector before the first getRayleighDampingForces() call")
            continue
        dm = re.search(LOCAL_DECL.format(name=re.escape(name)), stmts[decl_k][0])
        seeded_at = decl_k if dm.group("init") == "=" else None
        if seeded_at is None:
            seeded_at = next((j for j in range(decl_k + 1, k + 1)
                              if re.search(r"\b" + re.escape(name) + r"\b", stmts[j][0])), None)
        if seeded_at is None or seeded_at >= first_ray:
            yield span, (f"local '{name}' in {f.name} is seeded AFTER the first getRayleighDampingForces() "
                         "call; the re-entry has already happened -- seed the snapshot before it")


def check_rayleigh(root, rel, used_waivers=None):
    findings = []
    used = set() if used_waivers is None else used_waivers
    for path, raw, cl in _sources(root, stamped_only=True, needles=("getRayleighDampingForces",)):
        for f in functions(cl):
            for (line, last), msg in _analyze_rayleigh(cl, f):
                # waiver: anywhere in the statement's line span, or on the line above it
                wl, reason = waiver_at(raw, last, "rayleigh-ok", above=last - line + 1)
                if wl is not None:
                    used.add((str(path), wl))
                    if len(reason) >= MIN_REASON:
                        continue
                    msg = "rayleigh-ok waiver reason too short"
                findings.append(f"L1 {rel(path)}:{line + 1}: {msg}")
    return findings


# --------------------------------------------------------------------------
# L2
# --------------------------------------------------------------------------
def _enclosing(funcs, li):
    best = None
    for f in funcs:
        if f.open <= li <= f.end and (best is None or f.open >= best.open):
            best = f
    return best


def check_wipe(root, rel, used_waivers=None):
    findings = []
    used = set() if used_waivers is None else used_waivers
    resets, hook_text = {}, []
    for _, _, cl in _sources(root, stamped_only=False, needles=("::instance", "clearAll")):
        funcs = functions(cl)
        for f in funcs:
            short = f.name.split("::")[-1]
            if f.name in HOOKS_QUALIFIED or ("::" not in f.name and HOOK_FREE.match(short)):
                hook_text.append("\n".join(cl[f.open:f.end + 1]))
        for li, line in enumerate(cl):
            for m in RESET_CALL.finditer(line):
                f = _enclosing(funcs, li)
                resets.setdefault(m.group(1), []).append(f.name if f else "")
    hooks = "\n".join(hook_text)

    def wired(cls):
        for fname in resets.get(cls, []):
            short = fname.split("::")[-1]
            if fname in HOOKS_QUALIFIED or ("::" not in fname and HOOK_FREE.match(short)):
                return True
            if fname and "::" not in fname and re.search(r"(?<![\w:.>])" + re.escape(fname) + r"\s*\(", hooks):
                return True
        return False

    for path, raw, cl in _sources(root, stamped_only=True, needles=("instance",)):
        for li, line in enumerate(cl):
            m = SINGLETON.search(line)
            if not m:
                continue
            cls = m.group(1)
            wl, reason = waiver_at(raw, li, "wipe-ok", above=2)
            if wl is not None:
                used.add((str(path), wl))
                if len(reason) >= MIN_REASON:
                    continue
                findings.append(f"L2 {rel(path)}:{li + 1}: wipe-ok waiver reason too short for '{cls}'")
                continue
            if not wired(cls):
                findings.append(
                    f"L2 {rel(path)}:{li + 1}: singleton '{cls}' is not reset on wipe; call "
                    f"{cls}::instance().reset...() from Domain::clearAll or an OPS_clearAll* hook "
                    "(directly, or via a free function the hook calls), or waive with "
                    "'// ladruno-lint: wipe-ok <reason>'")
    return findings


def check_stale_waivers(root, rel, used):
    findings = []
    for path, raw, _ in _sources(root, stamped_only=True, needles=("ladruno-lint",)):
        for li, line in enumerate(raw):
            m = WAIVER.search(line)
            if m and (str(path), li) not in used:
                findings.append(f"W  {rel(path)}:{li + 1}: stale {m.group(1)} waiver -- it no longer "
                                "suppresses anything; delete it")
    return findings


# --------------------------------------------------------------------------
# L3
# --------------------------------------------------------------------------
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
    for path, raw, _ in _sources(root, stamped_only=True):
        for li, line in enumerate(raw):
            m = WAIVER.search(line)
            if m:
                print(f"{rel(path)}:{li + 1}: {m.group(1)} {m.group(2).strip()}")


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
    used = set()
    findings = []
    if "L1" in wanted:
        findings += check_rayleigh(root, rel, used)
    if "L2" in wanted:
        findings += check_wipe(root, rel, used)
    if {"L1", "L2"} <= wanted:
        findings += check_stale_waivers(root, rel, used)
    if "L3" in wanted:
        findings += check_pointers(root, rel)
    for f in findings:
        print(f)
    print(f"check_quirk_patterns: {len(findings)} finding(s) [{','.join(sorted(wanted))}]")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
