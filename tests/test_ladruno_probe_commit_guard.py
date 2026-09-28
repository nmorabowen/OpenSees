"""Guard: throwaway material copies must never trip a REAL commit (WP concrete3d-hang-diagnosis #877 follow-up).

LadrunoConcrete3D::commitState() declares a refused trial through the process-wide WP-99 counter
(ladrunoNoteCommitRefusal, SRC/material/LadrunoMaterialStatus.h), and Domain::commit() fails the whole step
when the counter is non-zero. LadrunoBrick and LadrunoQuad own hourglass-SHADOW material copies that they commit
through the same commitState() path; a refusal in such a probe says nothing about the element's real state.
Measured: K&R coarse (LadrunoQuad ssp) aborted at 0.6 mm on a probe refusal while every real Gauss point was
healthy -- the shadow commit was unguarded in LadrunoQuad.

The seam is LadrunoProbeCommitScope (an RAII guard that turns ladrunoNoteCommitRefusal() into a no-op). It is an
opt-in per owning element, so this test turns "a third element with shadows forgets it" into a red test instead of
a silent hang/abort: every source file that calls commitHgShadows() must open the scope, and the counter seam
must honour it.
"""
import os
import re

import pytest

pytestmark = [pytest.mark.zone_a]

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SRC_ELEMENT = os.path.join(ROOT, "SRC", "element")
STATUS_H = os.path.join(ROOT, "SRC", "material", "LadrunoMaterialStatus.h")

# a CALL of a shadow-commit helper (not its definition `X::commitHgShadows(`): `->commitHgShadows(` or
# ` commitHgShadows();` at statement start / after `this->`
_CALL = re.compile(r"(->|\.)\s*commitHgShadows\s*\(|^\s*commitHgShadows\s*\(", re.M)


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
        if _CALL.search(text):
            seen.append(os.path.relpath(path, ROOT))
            if "LadrunoProbeCommitScope" not in text:
                offenders.append(os.path.relpath(path, ROOT))
    assert seen, "the scan found no element committing shadow probes -- the guard pattern is stale"
    assert not offenders, (
        "these elements commit throwaway shadow material copies without a LadrunoProbeCommitScope "
        f"(a probe refusal would fail their real commit): {offenders}"
    )


def test_probe_scope_actually_suppresses_the_refusal_note():
    with open(STATUS_H, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    assert "struct LadrunoProbeCommitScope" in text
    m = re.search(r"inline void ladrunoNoteCommitRefusal\(void\)\s*\{(.*?)\n\}", text, re.S)
    assert m, "ladrunoNoteCommitRefusal definition not found"
    body = m.group(1)
    assert "ladrunoProbeScopeDepth()" in body and "++ladrunoCommitRefusalCounter()" in body, (
        "ladrunoNoteCommitRefusal must be gated on the probe-scope depth"
    )
