"""WP-112 (TIMs F14): `-flipAlphaIn init` is the LadrunoSANISAND default, and
the first plastic steps after `updateMaterialStage` no longer depend on the MKL
thread count.

THE DEFECT (F14, TIMs fork request 2026-09-18). `ManzariDafalias::integrate()`
(`ManzariDafalias.cpp:1022-1026`) resets `alpha_in := alpha_n` when
`(alpha_n - alpha_in_n) : (Ce : d_eps) < 0`, with no magnitude guard on either
factor, and it runs that test in the ELASTIC stage too (it sits before the
`mElastFlag` branch). A `LoadControl(0)` hold's `d_eps` is solver noise, so a
hold sets `alpha_in := alpha_n` at a coin-flip of points, and on the next step
`alpha_n - alpha_in_n` is only the round-off by which `alpha` moved since. On
the first PLASTIC increment the SIGN of that round-off vector against a real
`d_eps` picks the branch, and the branch picks the plastic modulus. The TIMs
strip read the first push step at 1.511 / 1.824 / 1.824 / 1.489 kPa at
1 / 2 / 4 / 8 MKL threads under the then-default `vanilla`, and 1.824 / 14.339 /
36.586 kPa at rows 1 / 8 / 15 on every thread count under `init`.

THE DECK. A small SELF-WEIGHT PLANE-STRAIN strip: 12 x 6 `LadrunoQuad
-formulation bbar -type PlaneStrain` over 6 m x 3 m, Gorini's `_PARAMS`,
`IntScheme 1` / `TanType 2`, gravity (gamma = 9.81 as nodal loads) plus a
20 kPa surface surcharge ramped at stage 0 in 10 steps, `loadConst`, TWO
`LoadControl(0)` holds (the F14 trigger -- measured below), `updateMaterialStage
1`, then a rigid 1 m footing (`equalDOF` on its surface nodes) pushed by
`DisplacementControl` at 5e-5 m per step, `system Pardiso` (MKL), `KrylovNewton`,
`NormDispIncr 1e-8`. The output per step is the footing load (the load factor
of the unit reference load on the footing master node, kN/m), compared by
`float.hex`, i.e. bit for bit.

MEASURED BEFORE THE BUILD (2026-09-18, release build 48c0e99bc, which still
defaults to `vanilla`; the new default was simulated by passing `-flipAlphaIn
init` explicitly -- WP-112 changes nothing else):

  * at the flip, after the two holds, 229 of 288 Gauss points carry
    0 < ||alpha - alpha_in|| / max(||alpha||, m) < 1e-12 (round-off) and the
    other 59 read exactly 1.0 (alpha_in = 0): the deck carries the F14 state.
    With no holds all 288 read 1.0 -- the holds are the trigger.
  * `init`: the ten push steps read 9.659111 13.158149 16.291804 19.268990
    22.171316 25.042673 27.860072 30.632409 33.366440 36.067358 kN/m, and are
    BIT-IDENTICAL at MKL_NUM_THREADS = 1 / 2 / 4 / 8.
  * the hold-count lottery, at one thread: the first push step under
    `vanilla` reads 4.107 / FAIL / 8.332 / 9.483 / FAIL kN/m after 0 / 1 / 2 /
    3 / 4 holds, i.e. the branch is decided by round-off the holds leave
    behind; under `init` it reads 9.659111 after every one of them (to
    3e-14 relative: the holds are real solves and move the state by ULPs,
    which is not a branch).
    (Pardiso on this 170-DOF system does not change its arithmetic with the
    thread count, so the lottery shows here through the hold count instead;
    the TIMs strip, at 9 720 Gauss points, showed it through the threads.
    The thread-count leg below is therefore a regression guard, and the
    hold-count leg is what makes this file non-vacuous.)
  * STATED, not asserted: under `init` the curve from step 4 on still moves
    with the hold count (step 10: 35.72 / 36.06 / 36.07 / 35.54 / 35.53 kN/m
    after 0..4 holds). The holds perturb the committed state by round-off
    and later branch decisions amplify it; `init` removes the step-1 sign
    lottery at the flip, not every sensitivity of the model to its state.

Each thread count runs in its OWN subprocess because MKL reads
MKL_NUM_THREADS once, at its first initialisation. `system Pardiso -stats`
prints `threads=<mkl_get_max_threads()>`, which the test reads back to prove
the setting reached MKL.
"""
import json
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="system Pardiso requires MKL (Windows/oneAPI build)"),
]

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_THREADS = (1, 2, 4, 8)
_N_PUSH = 10
_HOLDS = 2
# Hold-count leg: the default's first step may move by round-off (measured
# 3e-14 relative across 0..3 holds); 1e-9 is five decades above that and seven
# below vanilla's branch spread (measured 0.57 = (9.483 - 4.107)/9.483, plus two
# outright failures). A spread above 1e-2 is a branch, not round-off.
_ROUNDOFF_SPREAD = 1.0e-9
_BRANCH_SPREAD = 1.0e-2

# The child: builds the deck once per (mode, holds) config, in order, and
# prints ONE `RESULT <json>` line. `mode` is 'default' (no -flipAlphaIn token),
# 'init' or 'vanilla'.
_CHILD = r'''
import json, math, sys
sys.path.insert(0, sys.argv[1])
from _testbed import ops

PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005,
          1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
NX, NZ, W, H, B = 12, 6, 6.0, 3.0, 1.0
GAMMA, Q_SUR, NG, DS = 9.81, 20.0, 10, 5.0e-5
configs = json.loads(sys.argv[2])
n_push = int(sys.argv[3])
stats = sys.argv[4] == '1'


def run(mode, holds):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    dx, dz = W / NX, H / NZ
    nid = lambda i, k: 1 + i + (NX + 1) * k
    for k in range(NZ + 1):
        for i in range(NX + 1):
            ops.node(nid(i, k), -0.5 * W + i * dx, -H + k * dz)
    flip = () if mode == 'default' else ('-flipAlphaIn', mode)
    ops.nDMaterial('LadrunoSANISAND', 1, *PARAMS, 1, 2, 1, 1e-7, 1e-7,
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * 101.0, *flip)
    e = 0
    for k in range(NZ):
        for i in range(NX):
            e += 1
            ops.element('LadrunoQuad', e, nid(i, k), nid(i + 1, k), nid(i + 1, k + 1),
                        nid(i, k + 1), 1, '-thick', 1.0, '-type', 'PlaneStrain',
                        '-formulation', 'bbar')
    for k in range(1, NZ + 1):
        ops.fix(nid(0, k), 1, 0)
        ops.fix(nid(NX, k), 1, 0)
    for i in range(NX + 1):
        ops.fix(nid(i, 0), 1, 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    fz = GAMMA * dx * dz / 4.0
    for k in range(NZ):
        for i in range(NX):
            for (a, c) in ((i, k), (i + 1, k), (i + 1, k + 1), (i, k + 1)):
                ops.load(nid(a, c), 0.0, -fz)
    for i in range(NX + 1):
        ops.load(nid(i, NZ), 0.0, -Q_SUR * (dx if 0 < i < NX else 0.5 * dx))
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('Pardiso', *(('-stats',) if stats else ()))
    ops.test('NormDispIncr', 1.0e-8, 100, 0)
    ops.algorithm('KrylovNewton')
    ops.integrator('LoadControl', 1.0 / NG)
    ops.analysis('Static')
    ops.updateMaterialStage('-material', 1, '-stage', 0)
    for s in range(NG):
        assert ops.analyze(1) == 0, 'gravity step %d failed' % (s + 1)
    ops.loadConst('-time', 0.0)
    ops.integrator('LoadControl', 0.0)
    for s in range(holds):
        assert ops.analyze(1) == 0, 'hold %d failed' % (s + 1)
    ops.updateMaterialStage('-material', 1, '-stage', 1)
    # census of the committed (alpha, alpha_in) pair, read before any plastic
    # trial (elements 2..: the stage dispatch reaches element 1 eagerly).
    n_round = 0
    n_tot = 0
    for ee in range(2, NX * NZ + 1):
        for gp in range(1, 5):
            a = ops.eleResponse(ee, 'material', gp, 'alpha')
            ai = ops.eleResponse(ee, 'material', gp, 'alpha_in')
            d = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, ai)))
            na = math.sqrt(sum(x * x for x in a))
            n_tot += 1
            if 0.0 < d <= 1.0e-8 * max(na, PARAMS[9]):
                n_round += 1
    foot = [nid(i, NZ) for i in range(NX + 1)
            if abs(-0.5 * W + i * dx) <= 0.5 * B + 1e-9]
    master = nid(NX // 2, NZ)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for n in foot:
        if n != master:
            ops.equalDOF(master, n, 2)
    ops.load(master, 0.0, -1.0)
    ops.integrator('DisplacementControl', master, 2, -DS)
    steps = []
    for s in range(n_push):
        r = ops.analyze(1)
        q = ops.getLoadFactor(2)
        steps.append([r, q.hex(), q])
        if r != 0:
            break
    return {'mode': mode, 'holds': holds, 'roundoff': n_round, 'census': n_tot,
            'steps': steps}


out = [run(m, h) for m, h in configs]
print('RESULT ' + json.dumps({'build': ops.ladrunoBuild(), 'runs': out}))
'''


def _run_child(configs, threads=1, n_push=_N_PUSH, stats=False):
    """One fresh interpreter with MKL_NUM_THREADS / OMP_NUM_THREADS = threads.
    Returns (parsed RESULT dict, combined output). Same Windows-safe flags as
    `_testbed.subprocess_run.run_python_script` (stdin=DEVNULL, utf-8)."""
    env = dict(os.environ)
    env['MKL_NUM_THREADS'] = str(threads)
    env['OMP_NUM_THREADS'] = str(threads)
    env.setdefault('LADRUNO_OPENSEES_QUIET', '1')
    p = subprocess.run(
        [sys.executable, '-u', '-c', _CHILD, _TESTS_DIR, json.dumps(configs),
         str(n_push), '1' if stats else '0'],
        env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace',
        timeout=600)
    text = p.stdout + (p.stderr or '')
    lines = [l for l in p.stdout.splitlines() if l.startswith('RESULT ')]
    assert p.returncode == 0 and lines, (
        'child failed (threads=%d, rc=%d)' % (threads, p.returncode), text[-3000:])
    return json.loads(lines[-1][len('RESULT '):]), text


def _pardiso_threads(text):
    vals = {int(tok.split('=')[1]) for tok in text.split()
            if tok.startswith('threads=') and tok.split('=')[1].isdigit()}
    return vals


def test_first_ten_push_steps_bit_identical_across_mkl_threads():
    """F14's own test: under the DEFAULT (no -flipAlphaIn token), the first ten
    push steps after `updateMaterialStage` are bit-identical at MKL_NUM_THREADS
    1 / 2 / 4 / 8."""
    parent_build = ops.ladrunoBuild()
    ref = None
    for t in _THREADS:
        res, text = _run_child([['default', _HOLDS]], threads=t, stats=True)
        assert res['build'] == parent_build, (
            'the child imported a different engine from the parent -- a stale '
            'or shadowing opensees.pyd; nothing below would mean anything',
            res['build'], parent_build)
        assert _pardiso_threads(text) == {t}, (
            'system Pardiso -stats did not report threads=%d -- MKL_NUM_THREADS '
            'did not reach MKL, so this leg would not test what it claims' % t,
            sorted(_pardiso_threads(text)))
        run = res['runs'][0]
        assert run['roundoff'] > 0, (
            'no Gauss point carries a round-off alpha - alpha_in after the holds '
            '-- the deck no longer reaches the F14 state (measured 229/288 on '
            '48c0e99bc)', run)
        steps = run['steps']
        assert len(steps) == _N_PUSH and all(s[0] == 0 for s in steps), (
            'a push step failed to converge at threads=%d under the default' % t,
            steps)
        hexes = [s[1] for s in steps]
        if ref is None:
            ref = hexes
            continue
        assert hexes == ref, (
            'the first %d push steps under the DEFAULT -flipAlphaIn differ '
            'between MKL_NUM_THREADS=1 and %d -- the F14 thread-count dependence '
            'is back' % (_N_PUSH, t), ref, hexes)


def test_default_first_step_is_immune_to_the_hold_lottery():
    """Non-vacuity for the thread-count leg: on THIS deck the round-off the
    holds leave behind decides vanilla's first push step (measured 4.107 /
    FAIL / 8.332 / 9.483 kN/m after 0..3 holds), while the default reads
    the same value to round-off after every hold count. One subprocess, one
    thread."""
    configs = [['default', h] for h in (0, 1, 2, 3)] + \
              [['vanilla', h] for h in (0, 1, 2, 3)]
    res, _ = _run_child(configs, threads=1, n_push=1)
    runs = {(r['mode'], r['holds']): r for r in res['runs']}

    # The holds are real solves, so they move the committed state by round-off
    # and the default's first step moves with it -- by ULPs (measured spread
    # 3e-14 relative: 0x1.351770cc47098p+3 .. 0x1.351770cc47124p+3), never by
    # a branch. Vanilla's moves by the branch: 4.107 / FAIL / 8.332 / 9.483.
    default_first = {h: runs[('default', h)]['steps'][0] for h in (0, 1, 2, 3)}
    assert all(v[0] == 0 for v in default_first.values()), (
        'the default leg failed its first push step', default_first)
    q = [v[2] for v in default_first.values()]
    spread = (max(q) - min(q)) / max(abs(x) for x in q)
    assert spread <= _ROUNDOFF_SPREAD, (
        'under the DEFAULT the first push step moved by %.3e (relative) with '
        'the number of zero-load holds before the flip -- more than round-off, '
        'so a branch was picked by the holds' % spread, default_first)

    vanilla_first = {h: runs[('vanilla', h)]['steps'][0] for h in (0, 1, 2, 3)}
    failed = [h for h, v in vanilla_first.items() if v[0] != 0]
    qv = [v[2] for v in vanilla_first.values() if v[0] == 0]
    v_spread = ((max(qv) - min(qv)) / max(abs(x) for x in qv)) if len(qv) > 1 else 0.0
    assert failed or v_spread >= _BRANCH_SPREAD, (
        'explicit -flipAlphaIn vanilla gave the same first push step (spread '
        '%.3e, no failure) after 0..3 holds -- the deck no longer exhibits the '
        'sign-at-round-off lottery, so the default leg above proves nothing'
        % v_spread, vanilla_first)
    assert runs[('vanilla', 0)]['roundoff'] == 0 and runs[('vanilla', 2)]['roundoff'] > 0, (
        'the holds are supposed to be the trigger: zero round-off points '
        'without them, many with two',
        runs[('vanilla', 0)]['roundoff'], runs[('vanilla', 2)]['roundoff'])


_WARN = 'is at round-off'


def test_vanilla_roundoff_warning_and_default_silence():
    """The WP-112 warning: under explicit `vanilla` on the F14 state it fires,
    once per Gauss point, within the 10-per-process print budget (then one
    suppression line); it is silent under the default and under vanilla
    without the holds (alpha_in = 0, a genuine difference). Fresh processes,
    so the process-wide budget starts at zero in each."""
    _, text_v = _run_child([['vanilla', _HOLDS]], threads=1, n_push=1)
    n_warn = text_v.count(_WARN)
    assert 1 <= n_warn <= 10, (
        'expected the round-off alpha_in warning between 1 and 10 times '
        '(budget 10 per process) under -flipAlphaIn vanilla after two holds',
        n_warn, text_v[-3000:])
    if n_warn == 10:
        assert text_v.count('further round-off alpha_in warnings suppressed') == 1, (
            'the budget was spent but the suppression line did not print once',
            text_v[-3000:])

    _, text_v0 = _run_child([['vanilla', 0]], threads=1, n_push=1)
    assert _WARN not in text_v0, (
        'the warning fired on a deck whose alpha - alpha_in is a genuine '
        'difference (alpha_in = 0, no holds) -- the threshold is wrong',
        text_v0[-3000:])

    _, text_d = _run_child([['default', _HOLDS]], threads=1, n_push=1)
    assert _WARN not in text_d, (
        'the warning fired under the DEFAULT (init), where alpha_in := alpha '
        'at the flip makes the difference exactly zero', text_d[-3000:])
