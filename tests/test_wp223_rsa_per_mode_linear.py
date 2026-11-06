"""WP-223: `responseSpectrumPerMode -linear`, one element pass per mode.

ADR: Ladruno_implementation/214_ladruno_rsa_per_mode_output_adr.md, D10.

Without -linear the command runs one revert / update / set / update pass per
(direction, mode); with it, one pass per mode at the direction of largest
|Gamma|, scaled to the others by Gamma_d / Gamma_ref, and a self-check pass
at -1/2 the amplitude that refuses the call when the touched elements are not
linear about the committed state.

  (a) flag off: byte-identical to the build before WP-223. The baseline
      tests/data/wp223_rsa_flagoff_golden.json was recorded on ladruno
      d56fbfc48 (no WP-223); it also stores a fingerprint of the modal input
      (eigenvalues, eigenvectors, modalProperties), and a case whose input
      differs on this platform is skipped, not compared;
  (b) flag on agrees with flag off on linear models: the shell / Timoshenko
      tower (ASDShellQ4 with and without -pDelta, over a gravity preload and
      loadConst) and elastic frames with Linear and PDelta geomTransf; the
      scaled directions are exactly Gamma_d / Gamma_ref x the reference one;
  (c) the self-check refuses nonlinear models (Corotational frames under a
      large preload, a Steel01 base spring past yield, symmetric, which a
      +/- check would miss) and leaves the domain at the committed state;
  (d) options: -modes subsets, one direction, -linearCheck list / all,
      -linearTol, refusals; Tcl parity.
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops
from _testbed.rsa_tower import build_tower, spectrum

pytestmark = [pytest.mark.zone_a]

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = Path(__file__).resolve().parent / 'data' / 'wp223_rsa_flagoff_golden.json'
TN, SA = spectrum()
SPEC = ['-Tn', *TN, '-Sa', *SA]
ABOUT = (15.0, 15.0, 0.0)
NMODES = 6


@pytest.fixture(autouse=True)
def _quiet(monkeypatch):
    monkeypatch.setenv('LADRUNO_OPENSEES_QUIET', '1')


# --------------------------------------------------------------------------- models
def _gravity(loads):
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n, f in loads:
        ops.load(n, *f)
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)


def _tower(transf='Linear', shell_args=(), preload=True, h=1.5, eig='-genBandArpack'):
    t = build_tower(2, h, transf=transf, shell_args=shell_args)
    if preload:
        loads = []
        for s in t['storey_nodes']:
            for n in t['storey_nodes'][s]:
                if ops.nodeCoord(n)[2] == pytest.approx(s * 3.0):
                    loads.append((n, (0.0, 0.0, -2.0e4, 0.0, 0.0, 0.0)))
        _gravity(loads)
    ops.eigen(eig, NMODES)
    ops.modalProperties('-unorm')
    return t


def _tower_args(t):
    core = [n for n in t['supports'] if 11.0 <= ops.nodeCoord(n)[0] <= 19.0
            and 11.0 <= ops.nodeCoord(n)[1] <= 19.0]
    e1, c1 = t['piers'][1]
    e2, c2 = t['piers'][2]
    every = sorted(set().union(*t['storey_nodes'].values()))
    return ['-reactions', '-about', *ABOUT, '-reactionSet', 'core', *core,
            '-eleGroup', 'pier1', '-ele', *e1, '-cut', *c1, '-about', *ABOUT,
            '-eleGroup', 'pier2', '-ele', *e2, '-cut', *c2,
            '-inertia', 'all', *every, '-inertia', 'st2', *t['storey_nodes'][2], '-about', *ABOUT]


def _frame3d(transf='PDelta', preload=True, nst=4, P=4.0e5):
    """2x2-bay elastic frame, lumped lateral mass, gravity on the columns"""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    xs, h = (0.0, 6.0, 12.0), 3.0
    tag = lambda i, j, s: 1 + i + 3 * j + 9 * s
    for s in range(nst + 1):
        for j in range(3):
            for i in range(3):
                n = tag(i, j, s)
                ops.node(n, xs[i], xs[j] + 0.4 * i, s * h)     # a little plan skew
                if s == 0:
                    ops.fix(n, 1, 1, 1, 1, 1, 1)
                else:
                    ops.mass(n, 2.0e4, 2.0e4, 0.0, 0.0, 0.0, 0.0)
    ops.geomTransf(transf, 1, 1.0, 0.0, 0.0)
    ops.geomTransf('Linear', 2, 0.0, 0.0, 1.0)
    e = 0
    for s in range(nst):
        for j in range(3):
            for i in range(3):
                e += 1
                ops.element('elasticBeamColumn', e, tag(i, j, s), tag(i, j, s + 1),
                            0.25, 30e9, 12.5e9, 0.009, 0.0052, 0.0052, 1)
    for s in range(1, nst + 1):
        for j in range(3):
            for i in range(3):
                for a, b in ((tag(i, j, s), tag(i + 1, j, s)) if i < 2 else (None, None),
                             (tag(i, j, s), tag(i, j + 1, s)) if j < 2 else (None, None)):
                    if a is not None:
                        e += 1
                        ops.element('elasticBeamColumn', e, a, b,
                                    0.28, 30e9, 12.5e9, 0.006, 0.0114, 0.0037, 2)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
    ops.test('NormDispIncr', 1e-12, 20); ops.algorithm('Newton')
    if preload:
        _gravity([(tag(i, j, s), (0.0, 0.0, -P, 0.0, 0.0, 0.0))
                  for s in range(1, nst + 1) for j in range(3) for i in range(3)])
    else:
        ops.integrator('LoadControl', 0.0)
        ops.analysis('Static')
    ops.eigen('-fullGenLapack', NMODES)
    ops.modalProperties('-unorm')
    base = [tag(i, j, 0) for j in range(3) for i in range(3)]
    cols1 = list(range(1, 10))
    every = [tag(i, j, s) for s in range(1, nst + 1) for j in range(3) for i in range(3)]
    return ['-reactions', '-about', 6.0, 6.0, 0.0, '-reactionSet', 'edge', *base[:3],
            '-eleGroup', 'cols1', '-ele', *cols1, '-cut', *base, '-about', 6.0, 6.0, 0.0,
            '-eleGroup', 'cols1top', '-ele', *cols1,
            '-inertia', 'all', *every]


def _frame2d(kind, P=0.0, Fy=None, scale=1.0):
    """two-storey 2D column line: 'corot' Corotational columns under P;
    'spring' a zeroLength Steel01 base spring (yield Fy) under elastic columns"""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    E, A, I, L, m = 200e9, 0.01, 1.0e-4, 3.0, 1.0e4
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.node(3, 0.0, L); ops.mass(3, m, 0.0, 0.0)
    ops.node(4, 0.0, 2 * L); ops.mass(4, m, 0.0, 0.0)
    ops.geomTransf('Corotational' if kind == 'corot' else 'Linear', 1)
    if kind == 'corot':
        ops.equalDOF(1, 2, 1, 2, 3)
    else:
        ops.uniaxialMaterial('Steel01', 1, Fy, 5.0e7, 0.01)
        ops.uniaxialMaterial('Elastic', 2, 1.0e12)
        ops.element('zeroLength', 100, 1, 2, '-mat', 1, 2, 2, '-dir', 1, 2, 3)
    ops.element('elasticBeamColumn', 1, 2, 3, A, E, I, 1)
    ops.element('elasticBeamColumn', 2, 3, 4, A, E, I, 1)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
    ops.test('NormDispIncr', 1e-12, 50); ops.algorithm('Newton')
    if P:
        _gravity([(3, (0.0, -P, 0.0)), (4, (0.0, -P, 0.0))])
    else:
        ops.integrator('LoadControl', 0.0)
        ops.analysis('Static')
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties('-unorm')
    return ['-Tn', 0.0, 1.0, 10.0, '-Sa', 10.0 * scale, 10.0 * scale, 10.0 * scale,
            '-dir', 1, '-reactions', '-eleGroup', 'c1', '-ele', 1, '-cut', 2, '-inertia', 'top', 3, 4]


# --------------------------------------------------------------------------- helpers
def _run(*args):
    out = ops.responseSpectrumPerMode(*args)
    nd, nm, nc = (int(v) for v in out['shape'])
    return {k: (np.array(v).reshape(nd, nm, nc) if ':' in k or k == 'reactions'
                else np.array(v)) for k, v in out.items()}


def _kind_err(on, off, ncomp):
    """per key: max over component kind (forces, moments) of max |on - off| /
    max |off| of that kind over every (dir, mode). Per (dir, mode) would be
    round-off over round-off where a resultant cancels by equilibrium (the
    net base shear of a torsional mode)"""
    nF = 3 if ncomp == 6 else 2
    worst = {}
    for key, ref in off.items():
        if not (':' in key or key == 'reactions'):
            continue
        err = 0.0
        for lo, hi in ((0, nF), (nF, ncomp)):
            s = np.abs(ref[..., lo:hi]).max()
            d = np.abs(on[key][..., lo:hi] - ref[..., lo:hi]).max()
            err = max(err, d / s if s > 0.0 else (np.inf if d > 0.0 else 0.0))
        worst[key] = err
    return worst


def _assert_scaled_exactly(on):
    """the non-reference directions are Gamma_d / Gamma_ref x the reference
    direction, bit for bit (IEEE division and product, as in the command);
    a direction with Gamma = 0 is exact zeros"""
    g = on['gamma'].reshape(len(on['dirs']), len(on['modes']))
    for a in range(g.shape[1]):
        if not np.any(g[:, a]):
            continue
        r = int(np.argmax(np.abs(g[:, a])))
        for d in range(g.shape[0]):
            if d == r:
                continue
            ratio = g[d, a] / g[r, a]
            for key, v in on.items():
                if ':' in key and key.startswith('inertia:'):
                    continue
                if ':' in key or key == 'reactions':
                    want = np.zeros(v.shape[2]) if ratio == 0.0 else ratio * v[r, a]
                    assert np.array_equal(v[d, a], want), (key, d, a)


def _fingerprint():
    """sha256 of the input the command reads: modal properties, eigenvectors
    and the committed displacements"""
    h = hashlib.sha256()
    mp = ops.modalProperties('-unorm', '-return')
    for k in sorted(mp):
        v = mp[k]
        for x in (v if isinstance(v, (list, tuple)) else [v]):
            for y in (x if isinstance(x, (list, tuple)) else [x]):
                if isinstance(y, float):
                    h.update(y.hex().encode())
    nmodes = len(mp['eigenLambda'])
    for n in ops.getNodeTags():
        # the committed (gravity) state: its solve can move in the last bit
        # under load, and the outputs are increments over it
        for y in ops.nodeDisp(n):
            h.update(float(y).hex().encode())
        for k in range(1, nmodes + 1):
            for y in ops.nodeEigenvector(n, k):
                h.update(float(y).hex().encode())
    return h.hexdigest()


def _hex(out):
    return {k: [float(x).hex() for x in np.ravel(v)] for k, v in out.items()}


# the flag-off cases of the golden baseline: (name, builder) -> call args
def _golden_cases():
    # coarse towers with the dense eigensolver: ARPACK + UmfPack eigenvectors
    # moved in the last bits between runs of the same build on this machine
    def tower_pdelta():
        t = _tower('PDelta', ('-pDelta',), h=3.0, eig='-fullGenLapack')
        return [*SPEC, '-dir', 1, 2, *_tower_args(t)]

    def tower_subset():
        t = _tower('Linear', (), preload=False, h=3.0, eig='-fullGenLapack')
        return [*SPEC, '-dir', 2, '-modes', 4, 1, 3, *_tower_args(t)]

    def frame_pdelta():
        args = _frame3d('PDelta')
        return [*SPEC, '-dir', 1, 2, '-scale', 0.8, *args]

    return [('tower_pdelta', tower_pdelta), ('tower_subset', tower_subset),
            ('frame_pdelta', frame_pdelta)]


def write_golden():
    """record the baseline on a build WITHOUT WP-223 (see the module doc)"""
    rows = {}
    for name, build in _golden_cases():
        args = build()
        fp = _fingerprint()
        rows[name] = {'fingerprint': fp, 'out': _hex(ops.responseSpectrumPerMode(*args))}
    GOLDEN.write_text(json.dumps({
        'note': 'responseSpectrumPerMode flag-off outputs, float.hex, recorded on a build '
                'without WP-223; tests/test_wp223_rsa_per_mode_linear.py write_golden()',
        'build': str(ops.ladrunoBuild()), 'platform': sys.platform, 'rows': rows}, indent=1))


# --------------------------------------------------------------------------- (a)
@pytest.mark.parametrize('name', [c[0] for c in _golden_cases()])
def test_a_flag_off_is_byte_identical_to_the_build_before(name):
    gold = json.loads(GOLDEN.read_text())
    row = gold['rows'][name]
    args = dict(_golden_cases())[name]()
    if _fingerprint() != row['fingerprint']:
        pytest.skip(f'{name}: the modal input differs from the {gold["platform"]} '
                    f'baseline on this platform; byte identity is checked on the same input only')
    got = _hex(ops.responseSpectrumPerMode(*args))
    assert set(got) == set(row['out'])
    for k in got:
        assert got[k] == row['out'][k], k


# --------------------------------------------------------------------------- (b)
@pytest.mark.parametrize('transf,shell', [('Linear', ()), ('Linear', ('-pDelta',)),
                                          ('PDelta', ('-pDelta',))])
def test_b_shell_tower_flag_on_equals_flag_off(transf, shell):
    t = _tower(transf, shell)
    args = [*SPEC, '-dir', 1, 2, *_tower_args(t)]
    off = _run(*args)
    on = _run(*args, '-linear')
    err = _kind_err(on, off, 6)
    tol = 1e-9 if transf == 'Linear' else PDELTA_TOL
    assert max(err.values()) <= tol, err
    for k in ('inertia:all', 'inertia:st2'):
        assert np.array_equal(on[k], off[k]), k               # inertia: same code
    _assert_scaled_exactly(on)
    assert list(on['linearCheckModes']) == [1]
    assert on['linearCheckError'][0] <= on['linearTol'][0] == 1e-6
    assert 'linearCheckError' not in off


# The PDelta transformation is not linear about the committed state: it
# multiplies the trial axial force (committed + increment) by the trial chord
# drift, so the increment carries dP x dDelta (ledger quirk WP-223). In the
# tower the walls carry the overturning and the tolerance holds; in a pure
# frame it does not, see test_c_pdelta_frame...
PDELTA_TOL = 1e-6


def test_b_elastic_frame_flag_on_equals_flag_off():
    args = [*SPEC, '-dir', 1, 2, *_frame3d('Linear')]
    off = _run(*args)
    on = _run(*args, '-linear', '-linearCheck', 'all')
    err = _kind_err(on, off, 6)
    assert max(err.values()) <= 1e-12, err
    _assert_scaled_exactly(on)
    assert list(on['linearCheckModes']) == list(range(1, NMODES + 1))
    assert on['linearCheckError'][0] < 1e-8


def test_b_mode_subset_one_direction_and_time_series():
    args = _frame3d('Linear')
    ops.timeSeries('Path', 77, '-time', *TN, '-values', *SA)
    call = ['-ts', 77, '-scale', 0.7, '-dir', 2, '-modes', 5, 2, *args]
    off = _run(*call)
    on = _run(*call, '-linear', '-linearCheck', 2)
    assert max(_kind_err(on, off, 6).values()) <= 1e-12       # one direction: same passes
    assert list(on['modes']) == [5, 2] and list(on['linearCheckModes']) == [2]
    on2 = _run('-ts', 77, '-scale', 0.7, '-dir', 2, 1, '-modes', 5, 2, *args, '-linear')
    assert list(on2['linearCheckModes']) == [5]               # default: first listed
    _assert_scaled_exactly(on2)


def test_b_zero_participation_direction_gives_exact_zeros():
    """a mode with Gamma_d = 0 exactly: zeros there, not round-off"""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2, 0.0, 3.0); ops.mass(2, 1.0e4, 0.0, 0.0)
    ops.node(3, 5.0, 0.0); ops.fix(3, 1, 1, 1)
    ops.node(4, 5.0, 0.0); ops.fix(4, 1, 0, 1); ops.mass(4, 0.0, 2.0e3, 0.0)
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', 1, 1, 2, 0.01, 200e9, 1e-4, 1)
    ops.uniaxialMaterial('Elastic', 1, 1.0e6)
    ops.element('zeroLength', 2, 3, 4, '-mat', 1, '-dir', 2)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 0.0); ops.algorithm('Linear'); ops.analysis('Static')
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties('-unorm')
    on = _run('-Tn', 0.0, 1.0, 10.0, '-Sa', 2.0, 3.0, 0.5, '-dir', 1, 2, '-reactions', '-linear')
    g = on['gamma'].reshape(2, 2)
    zero = [(d, a) for d in range(2) for a in range(2) if g[d, a] == 0.0]
    if not zero:
        pytest.skip('the eigensolver left round-off in the decoupled block')
    for d, a in zero:
        assert np.all(on['reactions'][d, a] == 0.0)
    _assert_scaled_exactly(on)


# --------------------------------------------------------------------------- (c)
def _refused(call):
    with pytest.raises(Exception):
        ops.responseSpectrumPerMode(*call)


def test_c_corotational_frame_under_a_large_preload_is_refused():
    call = _frame2d('corot', P=6.0e5)       # ~0.27 of the cantilever Euler load
    off = _run(*call)
    _refused([*call, '-linear'])
    again = _run(*call)                      # the refusal left the domain committed
    for k in off:
        assert np.array_equal(off[k], again[k]), k
    # ...and the same model at a tiny amplitude passes: the check measures
    # the nonlinearity actually reached
    ok = _run(*_frame2d('corot', P=6.0e5, scale=1e-7), '-linear', '-linearTol', 1e-5)
    assert ok['linearCheckError'][0] <= 1e-5


def test_c_symmetric_yielding_spring_is_refused():
    """no preload, Steel01 base spring yielding under the mode: the response
    is odd in the amplitude, so only the -1/2 (not a -1) check sees it"""
    elastic = _run(*_frame2d('spring', Fy=1.0e8), '-linear')
    assert elastic['linearCheckError'][0] < 1e-9
    call = _frame2d('spring', Fy=2.0e4)
    off = _run(*call)
    assert abs(off['reactions'][0, 0, 0]) > 2.0e4 * 0.99    # the spring did yield
    _refused([*call, '-linear'])
    again = _run(*call)
    for k in off:
        assert np.array_equal(off[k], again[k]), k


def test_c_pdelta_frame_is_second_order_and_refused():
    """the exact path itself is not linear in the amplitude (dP x dDelta,
    the same for any preload P), so -linear refuses it at the default
    tolerance; with the tolerance raised past that measured gap it agrees
    with the exact path to about that gap"""
    for P in (4.0e5, 1.0e5):
        args = _frame3d('PDelta', P=P)
        one = _run(*SPEC, '-dir', 1, 2, *args)
        half = _run(*SPEC, '-scale', 0.5, '-dir', 1, 2, *args)
        x, y = one['reaction:edge'][..., :3], 2.0 * half['reaction:edge'][..., :3]
        gap = np.abs(x - y).max() / np.abs(x).max()
        assert 1e-4 < gap < 1e-2, gap                 # second order, not round-off
        _refused([*SPEC, '-dir', 1, 2, *args, '-linear'])
        on = _run(*SPEC, '-dir', 1, 2, *args, '-linear', '-linearTol', 0.05)
        assert 1e-4 < on['linearCheckError'][0] < 0.05
        assert max(_kind_err(on, one, 6).values()) < 0.02


def test_c_check_only_named_modes_and_tolerance():
    call = _frame2d('corot', P=6.0e5)
    _refused([*call, '-linear', '-linearCheck', 2])
    loose = _run(*call, '-linear', '-linearTol', 10.0)        # an explicit opt-out
    assert loose['linearCheckError'][0] > 1e-6


# --------------------------------------------------------------------------- (d)
def test_d_option_refusals():
    call = _frame2d('spring', Fy=1.0e8)
    for bad in (['-linearTol', 1e-6], ['-linearCheck', 1],
                ['-linear', '-linearTol', 0.0], ['-linear', '-linearTol'],
                ['-linear', '-linearCheck', 7], ['-linear', '-linearCheck']):
        _refused([*call, *bad])
    out = _run(*call, '-linear', '-linearCheck', 'all', '-linearTol', 1e-8)
    assert list(out['linearCheckModes']) == [1, 2]


def test_d_inertia_only_needs_no_element_pass():
    _frame2d('spring', Fy=2.0e4)              # would be refused with element groups
    call = ['-Tn', 0.0, 1.0, 10.0, '-Sa', 10.0, 10.0, 10.0, '-dir', 1, '-inertia', 'top', 3, 4]
    on = _run(*call, '-linear')
    off = _run(*call)
    assert np.array_equal(on['inertia:top'], off['inertia:top'])
    assert len(on['linearCheckModes']) == 0


def test_d_classic_tcl_parity(tmp_path):
    exe = ROOT / 'dist' / 'bin' / 'OpenSees.exe'
    if not exe.exists():
        pytest.skip('dist/bin/OpenSees.exe not built')
    deck = ['model basic -ndm 2 -ndf 3',
            'node 1 0 0', 'fix 1 1 1 1', 'node 2 0 3', 'mass 2 2.0e4 0 0',
            'node 3 0 6', 'mass 3 1.0e4 0 0', 'node 4 4 0', 'fix 4 1 1 1',
            'node 5 4 3', 'node 6 4 6', 'mass 5 2.0e4 0 0', 'mass 6 1.0e4 0 0',
            'geomTransf PDelta 1',
            'element elasticBeamColumn 1 1 2 0.01 200e9 1e-4 1',
            'element elasticBeamColumn 2 2 3 0.01 200e9 1e-4 1',
            'element elasticBeamColumn 3 4 5 0.01 200e9 1e-4 1',
            'element elasticBeamColumn 4 5 6 0.01 200e9 1e-4 1',
            'element elasticBeamColumn 5 2 5 0.01 200e9 2e-4 1',
            'element elasticBeamColumn 6 3 6 0.01 200e9 2e-4 1',
            'system UmfPack', 'numberer RCM', 'constraints Plain',
            'integrator LoadControl 0.0', 'algorithm Linear', 'analysis Static',
            'eigen -fullGenLapack 4', 'modalProperties -unorm',
            'set Tn {0.0 1.0 10.0}', 'set Sa {2.0 3.0 0.5}',
            'set r [responseSpectrumPerMode -Tn $Tn -Sa $Sa -dir 1 -reactions '
            '-eleGroup c -ele 1 3 -cut 1 4 -linear -linearCheck all]',
            'set bad [catch {responseSpectrumPerMode -Tn $Tn -Sa $Sa -dir 1 -reactions -linearTol 1e-6}]',
            'set out [open res.txt w]',
            'puts $out [dict get $r reactions]', 'puts $out [dict get $r ele:c]',
            'puts $out [dict get $r linearCheckModes]', 'puts $out $bad', 'close $out']
    (tmp_path / 'deck.tcl').write_text('\n'.join(deck) + '\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET='1')
    subprocess.run([str(exe), 'deck.tcl'], cwd=tmp_path, env=env, check=True,
                   stdin=subprocess.DEVNULL, capture_output=True, timeout=120)
    lines = (tmp_path / 'res.txt').read_text().split('\n')
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    for n, x, y in ((1, 0, 0), (2, 0, 3), (3, 0, 6), (4, 4, 0), (5, 4, 3), (6, 4, 6)):
        ops.node(n, float(x), float(y))
    ops.fix(1, 1, 1, 1); ops.fix(4, 1, 1, 1)
    for n, m in ((2, 2.0e4), (3, 1.0e4), (5, 2.0e4), (6, 1.0e4)):
        ops.mass(n, m, 0.0, 0.0)
    ops.geomTransf('PDelta', 1)
    for e, i, j, I in ((1, 1, 2, 1e-4), (2, 2, 3, 1e-4), (3, 4, 5, 1e-4), (4, 5, 6, 1e-4),
                       (5, 2, 5, 2e-4), (6, 3, 6, 2e-4)):
        ops.element('elasticBeamColumn', e, i, j, 0.01, 200e9, I, 1)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 0.0); ops.algorithm('Linear'); ops.analysis('Static')
    ops.eigen('-fullGenLapack', 4); ops.modalProperties('-unorm')
    py = ops.responseSpectrumPerMode('-Tn', 0.0, 1.0, 10.0, '-Sa', 2.0, 3.0, 0.5, '-dir', 1,
                                     '-reactions', '-eleGroup', 'c', '-ele', 1, 3, '-cut', 1, 4,
                                     '-linear', '-linearCheck', 'all')
    for line, key in zip(lines[:2], ('reactions', 'ele:c')):
        got = np.array(line.split(), float)
        np.testing.assert_allclose(got, py[key], rtol=1e-12, atol=1e-12 * np.abs(py[key]).max())
    assert [int(float(v)) for v in lines[2].split()] == [1, 2, 3, 4]
    assert lines[3].strip() == '1'


if __name__ == '__main__':
    if '--write-golden' in sys.argv:
        write_golden()
        print('wrote', GOLDEN)
