"""`ladrunoSchema` capability `mdStageFlagPreserved` keeps its promise (apeGmsh #2517).

The capability says: ManzariDafalias' class-wide stage flag (the static
`mElastFlag`) survives the construction of any further ManzariDafalias or
LadrunoSANISAND instance (fix 43b235bbb, fork PR #1005, apeGmsh #2444). apeGmsh
lifts its #2443 refusal on builds that list it, so a build that lists the name
while a constructor resets the flag would give a silently elastic answer there.

tests/test_manzari_stage_static_activation.py covers element activation, whose
`getCopy(const char*)` builds the copy through a FULL constructor. This file adds
two more checks:
  * behaviour: a NEW `nDMaterial` command after the flip (full constructor, no
    element, no copy) keeps the flag -- Python and classic Tcl;
  * source: no constructor of the ManzariDafalias / LadrunoSANISAND family
    assigns `mElastFlag` at all. This is the only tripwire for the classTag
    constructor `ManzariDafalias(int classTag)`: no deck here reaches it (only
    `*3D::getCopy()` and the object broker do), so a reset put back there alone
    passes every behaviour test while the build still lists the capability
    (review of fork PR #1017, measured).
Before 43b235bbb both full constructors and the classTag one wrote
`mElastFlag = 0`.

Deck: one laterally confined stdBrick cube, settled under a top load in stage 0
(elastic), `loadConst`, then optionally flipped to stage 1, then optionally a
second material of the same class is constructed (tag 2, unused), then loaded
further. uz at the top (node 5) is read.

Oracles:
  * flip + construct == flip alone (rel 1e-9): the new material is attached to
    nothing, so a kept flag gives the same answer bit for bit;
  * flip alone, for ManzariDafalias, matches stock openseespy 3.8.0 on this exact
    deck (-3.4067686516e-3, rel 1e-4; the fork's MD integrator changes move it
    by ~3e-6), and no flip matches stock's elastic answer (-2.8970896097e-3,
    rel 1e-9) -- both values from tests/test_manzari_stage_static_activation.py;
  * the elastic answer is 17 % away from the plastic one, so a reset (which
    returns the elastic answer) cannot pass the 1e-9 check.

The schema check and the behaviour check run in one test per class: the name
is listed AND the behaviour holds. Classic Tcl runs the same deck through
tests/tcl/md_stage_flag_capability.tcl.
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
DECK = REPO / "tests" / "tcl" / "md_stage_flag_capability.tcl"
CAPABILITY = "mdStageFlagPreserved"
MARKER = "MD_STAGE_FLAG_CAP_2517"

# Toyoura-sand set (OpenSees docs example), as tests/test_manzari_safety_pack.py.
MD = [125, 0.05, 0.8, 1.25, 0.712, 0.019, 0.934, 0.7, 100, 0.01, 7.05,
      0.968, 1.1, 0.704, 3.5, 4, 600, 1.42]
XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]

STOCK_380_PLASTIC = -3.4067686516e-3    # stock openseespy 3.8.0, flip, no activation
STOCK_380_ELASTIC = -2.8970896097e-3    # stock openseespy 3.8.0, no flip

# ---- source: no family constructor writes the class-wide flag -----------------

_ND = REPO / "SRC" / "material" / "nD"
# Every class that inherits ManzariDafalias' static mElastFlag, with the file that
# defines its constructors and the minimum number of constructors it has (so a
# renamed file or a parser miss cannot pass vacuously).
FAMILY = {
    "ManzariDafalias": (_ND / "UWmaterials" / "ManzariDafalias.cpp", 4),
    "ManzariDafalias3D": (_ND / "UWmaterials" / "ManzariDafalias3D.cpp", 2),
    "ManzariDafaliasPlaneStrain": (_ND / "UWmaterials" / "ManzariDafaliasPlaneStrain.cpp", 2),
    "ManzariDafaliasRO": (_ND / "UWmaterials" / "ManzariDafaliasRO.cpp", 2),
    "ManzariDafalias3DRO": (_ND / "UWmaterials" / "ManzariDafalias3DRO.cpp", 2),
    "ManzariDafaliasPlaneStrainRO": (_ND / "UWmaterials" / "ManzariDafaliasPlaneStrainRO.cpp", 2),
    "LadrunoSANISAND": (_ND / "LadrunoSANISAND.cpp", 4),
    "LadrunoSANISAND3D": (_ND / "LadrunoSANISAND3D.cpp", 2),
    "LadrunoSANISANDPlaneStrain": (_ND / "LadrunoSANISANDPlaneStrain.cpp", 2),
}
_COMMENT_OR_STRING = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'',
                                re.S)
# any write: =, op=, ++/-- on either side, with or without a Class:: / this-> prefix
_WRITE = re.compile(r"\bmElastFlag\s*(?:[-+*/%|&^]|<<|>>)?=(?!=)|\bmElastFlag\s*(?:\+\+|--)"
                    r"|(?:\+\+|--)\s*(?:this\s*->\s*|\w+\s*::\s*)?mElastFlag\b")


def _blank(text):
    """Comments and string/char literals replaced by spaces (newlines kept)."""
    return _COMMENT_OR_STRING.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


def _match_close(text, i, open_ch, close_ch):
    depth = 0
    for j in range(i, len(text)):
        if text[j] == open_ch:
            depth += 1
        elif text[j] == close_ch:
            depth -= 1
            if depth == 0:
                return j
    raise AssertionError(f"unbalanced {open_ch}{close_ch} from offset {i}")


def constructor_spans(text, cls):
    """[(line, span text)] for each `cls::cls(...) [: init-list] { body }`
    definition in `text` (already blanked): signature, initialiser list and body."""
    out = []
    for m in re.finditer(r"\b%s\s*::\s*%s\s*\(" % (cls, cls), text):
        close = _match_close(text, m.end() - 1, "(", ")")
        k = close + 1
        while k < len(text) and text[k] != "{":
            if text[k] == ";":            # a declaration / call, not a definition
                break
            if text[k] == "(":            # an initialiser's argument list
                k = _match_close(text, k, "(", ")")
            k += 1
        if k >= len(text) or text[k] != "{":
            continue
        end = _match_close(text, k, "{", "}")
        out.append((text.count("\n", 0, m.start()) + 1, text[m.start():end + 1]))
    return out


def stage_flag_writes_in_constructors(sources):
    """{cls: (n constructors, [findings])} for {cls: source text}."""
    result = {}
    for cls, raw in sources.items():
        spans = constructor_spans(_blank(raw), cls)
        bad = [f"{cls} constructor at line {line}: writes mElastFlag "
               f"({_WRITE.search(span).group(0).strip()!r})"
               for line, span in spans if _WRITE.search(span)]
        result[cls] = (len(spans), bad)
    return result


def test_no_family_constructor_writes_the_stage_flag():
    sources = {cls: path.read_text(encoding="utf-8", errors="replace")
               for cls, (path, _) in FAMILY.items()}
    result = stage_flag_writes_in_constructors(sources)
    findings = []
    for cls, (n, bad) in result.items():
        assert n >= FAMILY[cls][1], (
            f"non-vacuity: found {n} {cls} constructor(s) in {FAMILY[cls][0].name}, "
            f"expected >= {FAMILY[cls][1]}")
        findings += bad
    assert not findings, (
        f"{CAPABILITY!r} promises no constructor resets the class-wide stage flag "
        "(apeGmsh #2517 / #2444):\n" + "\n".join(findings))


@pytest.mark.parametrize("body, writes", [
    ("ManzariDafalias::ManzariDafalias(int classTag)\n  : NDMaterial(0, classTag), mA(6)\n"
     "{\n  mTolF = 1e-7;\n  mElastFlag = 0;\n}\n", True),
    ("ManzariDafalias ::ManzariDafalias(int classTag) : NDMaterial(0, classTag)\n"
     "{ ManzariDafalias::mElastFlag=0; }\n", True),
    ("ManzariDafalias::ManzariDafalias() : NDMaterial(0, 1) { ++mElastFlag; }\n", True),
    ("ManzariDafalias::ManzariDafalias(int t) : NDMaterial(t, 1) { this->mElastFlag |= 1; }\n", True),
    # reads, comments and strings are not writes
    ("ManzariDafalias::ManzariDafalias(int t) : NDMaterial(t, 1)\n{\n"
     "  // vanilla reset mElastFlag = 0 here\n  /* mElastFlag = 0; */\n"
     "  opserr << \"mElastFlag = 0\";\n  if (mElastFlag == 0) mTolF = 1;\n}\n", False),
    # a write outside any constructor is not this rule's business
    ("void ManzariDafalias::setFlag() { mElastFlag = 1; }\n"
     "ManzariDafalias::ManzariDafalias() : NDMaterial(0, 1) { mTolF = 0; }\n", False),
])
def test_constructor_scan_self_test(body, writes):
    (n, bad), = stage_flag_writes_in_constructors({"ManzariDafalias": body}).values()
    assert n >= 1 and bool(bad) is writes, (n, bad)


def _require_fork_module():
    if ops.__name__.startswith("openseespy"):
        pytest.skip("openseespy wheel fallback in use -- fork build required")


def _run(mat, flip, construct):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(mat, 1, *MD)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in range(1, 9):
        ops.fix(n, 1, 1, 1 if n <= 4 else 0)
    ops.updateMaterialStage('-material', 1, '-stage', 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n in range(5, 9):
        ops.load(n, 0., 0., -25.)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25)
    ops.analysis('Static')
    assert ops.analyze(4) == 0, 'stage-0 (elastic) gravity failed'
    ops.loadConst('-time', 0.0)
    if flip:
        ops.updateMaterialStage('-material', 1, '-stage', 1)
    if construct:
        # a second instance through the FULL constructor -- the reset site
        # 43b235bbb removed; attached to no element
        ops.nDMaterial(mat, 2, *MD)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for n in range(5, 9):
        ops.load(n, 0., 0., -25.)
    ops.integrator('LoadControl', 0.05)
    assert ops.analyze(20) == 0, 'loading stage failed'
    return ops.nodeDisp(5, 3)


@pytest.mark.parametrize('mat', ['ManzariDafalias', 'LadrunoSANISAND'])
def test_capability_listed_and_kept(mat):
    _require_fork_module()
    caps = json.loads(ops.ladrunoSchema())["capabilities"]
    assert CAPABILITY in caps, f"{CAPABILITY!r} not in ladrunoSchema capabilities {caps}"

    plastic = _run(mat, flip=True, construct=False)
    plastic_new = _run(mat, flip=True, construct=True)
    elastic = _run(mat, flip=False, construct=False)

    assert abs(plastic_new - plastic) <= 1e-9 * abs(plastic), (
        f'{mat}: constructing another instance after the stage-1 flip moved uz5 '
        f'from {plastic:.10e} to {plastic_new:.10e} (elastic answer {elastic:.10e}) '
        f'-- a constructor reset the class-wide stage flag, yet the build lists '
        f'{CAPABILITY!r} (apeGmsh #2517 / #2444)')
    assert abs(plastic - elastic) > 0.1 * abs(elastic), (
        f'{mat}: the plastic and elastic legs are not distinguishable '
        f'({plastic:.10e} vs {elastic:.10e}); the deck lost its power')


def test_manzari_legs_match_stock_380():
    plastic = _run('ManzariDafalias', flip=True, construct=False)
    elastic = _run('ManzariDafalias', flip=False, construct=False)
    assert plastic == pytest.approx(STOCK_380_PLASTIC, rel=1e-4)
    assert elastic == pytest.approx(STOCK_380_ELASTIC, rel=1e-9)


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


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_capability_listed_and_kept():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO), stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert f"{MARKER}: PASS" in out, out[-3000:]
