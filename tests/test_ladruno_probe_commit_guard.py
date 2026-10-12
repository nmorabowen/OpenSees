"""Guard: throwaway material copies must never trip a REAL commit (WP concrete3d-hang-diagnosis #877 follow-up).

LadrunoConcrete3D::commitState() declares a refused trial through the process-wide WP-99 counter
(ladrunoNoteCommitRefusal, SRC/material/LadrunoMaterialStatus.h), and Domain::commit() fails the whole step
when the counter is non-zero. LadrunoBrick and LadrunoQuad own hourglass-SHADOW material copies that they commit
through the same commitState() path; a refusal in such a probe says nothing about the element's real state.
Measured: K&R coarse (LadrunoQuad ssp) aborted at 0.6 mm on a probe refusal while every real Gauss point was
healthy -- the shadow commit was unguarded in LadrunoQuad.

The seam is LadrunoProbeCommitScope (an RAII guard that turns ladrunoNoteCommitRefusal() into a no-op). It is an
opt-in per owning element, so this test turns "a third element with shadows forgets it" into a red test instead of
a silent hang/abort.

STRUCTURAL (review nit): the check parses the source instead of grepping for a substring. Comments and string
literals are stripped, braces are matched, and for every CALL of commitHgShadows() the test requires a
`LadrunoProbeCommitScope <name>;` declaration that is still IN SCOPE at the call (declared earlier in an enclosing
block), so a scope declared in an inner block that has already closed, or mentioned only in a comment, a string or a
different function, does not satisfy the guard. test_the_checker_itself_is_structural proves the checker rejects those.
"""
import os
import re

import pytest

pytestmark = [pytest.mark.zone_a]

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SRC_ELEMENT = os.path.join(ROOT, "SRC", "element")
STATUS_H = os.path.join(ROOT, "SRC", "material", "LadrunoMaterialStatus.h")

_DECL = re.compile(r"\bLadrunoProbeCommitScope\s+[A-Za-z_]\w*\s*;")
# a CALL of the shadow-commit helper: `->commitHgShadows(`, `.commitHgShadows(` or a bare statement `commitHgShadows(`
# (a definition `X::commitHgShadows(` is preceded by `::` and excluded)
_CALL = re.compile(r"(?:->|\.)\s*commitHgShadows\s*\(|(?<![:\w])commitHgShadows\s*\(")


def strip_comments_and_strings(text):
    """Blank out // and /* */ comments and string / char literals (newlines kept, so offsets and lines survive)."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i)); i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j])); i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j])); i = j
        else:
            out.append(c); i += 1
    return "".join(out)


def unguarded_calls(text):
    """Line numbers of every commitHgShadows() CALL that has no LadrunoProbeCommitScope declared in an enclosing,
    still-open block before it."""
    code = strip_comments_and_strings(text)
    events = [(m.start(), "decl") for m in _DECL.finditer(code)] + [(m.start(), "call") for m in _CALL.finditer(code)]
    events += [(i, ch) for i, ch in enumerate(code) if ch in "{}"]
    events.sort()
    stack, bad = [False], []                      # stack[k] = "a scope guard is declared in open block k"
    for pos, kind in events:
        if kind == "{":
            stack.append(False)
        elif kind == "}":
            if len(stack) > 1:
                stack.pop()
        elif kind == "decl":
            stack[-1] = True
        elif kind == "call" and not any(stack):
            bad.append(code.count("\n", 0, pos) + 1)
    return bad


def calls(text):
    return len(_CALL.findall(strip_comments_and_strings(text)))


def _cpp_files(root):
    for d, _dirs, files in os.walk(root):
        for f in files:
            if f.endswith((".cpp", ".cc")):
                yield os.path.join(d, f)


def test_every_element_committing_shadow_probes_opens_the_probe_scope():
    offenders, seen = [], []
    for path in _cpp_files(SRC_ELEMENT):
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        if calls(text):
            seen.append(os.path.relpath(path, ROOT))
            bad = unguarded_calls(text)
            if bad:
                offenders.append((os.path.relpath(path, ROOT), bad))
    assert seen, "the scan found no element committing shadow probes -- the guard pattern is stale"
    assert not offenders, (
        "these elements call commitHgShadows() with no LadrunoProbeCommitScope in scope at the call "
        f"(a probe refusal would fail their real commit): {offenders}"
    )


def test_the_checker_itself_is_structural():
    ok = "void E::commitState() { LadrunoProbeCommitScope probe; this->commitHgShadows(); }"
    assert unguarded_calls(ok) == []
    in_comment = "void E::commitState() { // LadrunoProbeCommitScope probe;\n this->commitHgShadows(); }"
    assert unguarded_calls(in_comment) == [2]
    in_string = 'void E::commitState() { const char *s = "LadrunoProbeCommitScope probe;"; this->commitHgShadows(); }'
    assert unguarded_calls(in_string) == [1]
    closed_block = "void E::commitState() { { LadrunoProbeCommitScope probe; } this->commitHgShadows(); }"
    assert unguarded_calls(closed_block) == [1]
    other_function = ("void A::f() { LadrunoProbeCommitScope probe; }\n"
                      "void B::commitState() { this->commitHgShadows(); }")
    assert unguarded_calls(other_function) == [2]
    enclosing = "void E::g() { LadrunoProbeCommitScope probe; if (x) { this->commitHgShadows(); } }"
    assert unguarded_calls(enclosing) == []
    definition_only = "void E::commitHgShadows() { for (;;) {} }"
    assert calls(definition_only) == 0


def test_probe_scope_actually_suppresses_the_refusal_note():
    with open(STATUS_H, encoding="utf-8", errors="replace") as fh:
        code = strip_comments_and_strings(fh.read())
    assert re.search(r"struct\s+LadrunoProbeCommitScope\b", code)
    start = code.index("inline void ladrunoNoteCommitRefusal")
    body_start = code.index("{", start)
    depth, i = 0, body_start
    while True:                                   # match the braces of the function body
        depth += code[i] == "{"
        depth -= code[i] == "}"
        if depth == 0:
            break
        i += 1
    body = code[body_start:i + 1]
    assert re.search(r"if\s*\(\s*ladrunoProbeScopeDepth\(\)\s*==\s*0\s*\)\s*\+\+\s*ladrunoCommitRefusalCounter\(\)", body), (
        "ladrunoNoteCommitRefusal must count a refusal only when no probe scope is open"
    )
