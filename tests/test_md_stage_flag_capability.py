"""`ladrunoSchema` capability `mdStageFlagPreserved` keeps its promise (apeGmsh #2517).

The capability says: ManzariDafalias' class-wide stage flag (the static
`mElastFlag`) survives the construction of any further ManzariDafalias or
LadrunoSANISAND instance (fix 43b235bbb, fork PR #1005, apeGmsh #2444). apeGmsh
lifts its #2443 refusal on builds that list it, so a build that lists the name
while a constructor resets the flag would give a silently elastic answer there.

tests/test_manzari_stage_static_activation.py covers element activation, which
copies the material through the classTag constructor (`getCopy()`). This file
covers the other path the capability names: a NEW `nDMaterial` command after
the flip, which runs the full constructor and attaches the material to no
element. Before 43b235bbb both full constructors wrote `mElastFlag = 0`.

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
