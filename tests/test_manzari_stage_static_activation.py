"""A later-constructed ManzariDafalias must not reset the class-wide stage (apeGmsh #2444).

`ManzariDafalias::mElastFlag` is a STATIC: `updateMaterialStage -stage 1` flips
every MD instance at once. Before the fix, the fork's two full constructors and
vanilla's classTag constructor (the one `getCopy()` runs) all wrote
`mElastFlag = 0`, so a stage that ACTIVATED more MD elements after the flip
silently sent every MD material in the domain back to elastic. apeGmsh saw it as
uy -1.6597e-3 -> -1.4923e-3, the elastic answer. Stock 3.8.0 never reset it.

Deck: a laterally confined column of stdBrick cubes. Stage 0 settles one cube
under a top load (elastic). The material is flipped to stage 1, then a second
cube is activated on top of it (its new top nodes ride on the settled interface,
so it starts stress-free) and the load is increased. The new cube only carries
the extra load down, so uz at the interface (node 5) is the same with or without
the activation.

Oracles:
  * flip + activate == flip without activation (rel 1e-4) -- independent of
    any reference binary. Stock openseespy 3.8.0 on this exact deck gives
    -3.4067686468e-3 / -3.4067686516e-3 for the two; the fork's MD legs agree
    to ~1e-9, LadrunoSANISAND's to ~2e-6 (its newly activated material runs
    its own one-shot flip priming on its first step);
  * the fork's MD matches stock 3.8.0 to rel 1e-4 (measured 3e-6: the fork carries
    integrator changes to MD, which move the plastic path slightly) and the
    elastic answer to 1e-9;
  * the elastic answer (no flip) is -2.8970896097e-3 on stock 3.8.0 and is
    17 % away, so the 1e-4 check cannot pass on an elastic run. Pre-fix the flip +
    activate leg returned exactly that elastic answer;
  * the default stays elastic: no stage call at all reproduces the stage-0
    answer, and so does the next model after `wipe` following a flip.

Same deck for LadrunoSANISAND (it inherits the static from its base), and once
through the classic Tcl engine (tests/tcl/manzari_stage_static_activation.tcl).
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "manzari_stage_static_activation.tcl"

# Toyoura-sand set (OpenSees docs example), as tests/test_manzari_safety_pack.py.
MD = [125, 0.05, 0.8, 1.25, 0.712, 0.019, 0.934, 0.7, 100, 0.01, 7.05,
      0.968, 1.1, 0.704, 3.5, 4, 600, 1.42]
XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]

STOCK_380_PLASTIC = -3.4067686468e-3    # stock openseespy 3.8.0, flip + activate
STOCK_380_ELASTIC = -2.8970896097e-3    # stock openseespy 3.8.0, no flip


def _run(mat, flip, activate, stage0_call=True):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(mat, 1, *MD)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in range(1, 9):
        ops.fix(n, 1, 1, 1 if n <= 4 else 0)
    if stage0_call:
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
    if activate:
        for j, (x, y) in enumerate(XY):
            ops.node(9 + j, x, y, 2.0)
            ops.fix(9 + j, 1, 1, 0)
            ops.setNodeDisp(9 + j, 3, ops.nodeDisp(5 + j, 3), '-commit')
        # constructs a new MD material from the prototype -- the reset site
        ops.element('stdBrick', 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for n in (range(9, 13) if activate else range(5, 9)):
        ops.load(n, 0., 0., -25.)
    ops.integrator('LoadControl', 0.05)
    assert ops.analyze(20) == 0, 'loading stage failed'
    return ops.nodeDisp(5, 3)


@pytest.mark.parametrize('mat', ['ManzariDafalias', 'LadrunoSANISAND'])
def test_activation_after_flip_stays_plastic(mat):
    plastic_no_add = _run(mat, flip=True, activate=False)
    plastic_add = _run(mat, flip=True, activate=True)
    elastic = _run(mat, flip=False, activate=False)

    assert abs(plastic_add - plastic_no_add) <= 1e-4 * abs(plastic_no_add), (
        f'{mat}: activating an element after the stage-1 flip moved uz5 from '
        f'{plastic_no_add:.10e} to {plastic_add:.10e} (elastic answer '
        f'{elastic:.10e}) -- a constructor reset the class-wide stage flag (#2444)')
    assert abs(plastic_no_add - elastic) > 0.1 * abs(elastic), (
        f'{mat}: the plastic and elastic legs are not distinguishable '
        f'({plastic_no_add:.10e} vs {elastic:.10e}); the deck lost its power')


def test_manzari_matches_stock_380():
    plastic_add = _run('ManzariDafalias', flip=True, activate=True)
    elastic = _run('ManzariDafalias', flip=False, activate=False)
    assert plastic_add == pytest.approx(STOCK_380_PLASTIC, rel=1e-4)
    assert elastic == pytest.approx(STOCK_380_ELASTIC, rel=1e-9)


@pytest.mark.parametrize('mat', ['ManzariDafalias', 'LadrunoSANISAND'])
def test_default_is_elastic_fresh_and_after_wipe(mat):
    elastic = _run(mat, flip=False, activate=False)
    _run(mat, flip=True, activate=True)          # leaves the static at 1 ...
    after_wipe = _run(mat, flip=False, activate=False, stage0_call=False)
    assert after_wipe == pytest.approx(elastic, rel=1e-9), (
        f'{mat}: the next model after `wipe` did not start elastic '
        f'({after_wipe:.10e} vs {elastic:.10e})')


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
def test_classic_tcl_activation_after_flip():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO), stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert "MD_STAGE_2444: PASS" in out, out[-3000:]
