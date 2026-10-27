"""WP-216: `responseSpectrumPerMode` — per-mode RSA resultants in one call.

ADR: Ladruno_implementation/214_ladruno_rsa_per_mode_output_adr.md.

The oracle is today's per-mode loop, run through the stock commands on the
same model: for each direction d and mode k
`responseSpectrumAnalysis d .. -mode k -increment`, then `reactions` +
`nodeReaction` (base reactions), `eleForce` (element-group resultants at a
cut) and lambda_k * M (u - u_commit) from the IncrDisp node response and
`nodeMass` (modal inertia; `nodeDisp` is the committed displacement). Every quantity is the increment over the committed state.

  (a) equal to the loop to 1e-12 per mode, relative to the scale of the
      summed terms (sum |term| per component; the loop itself is not
      bit-reproducible, see test_the_loop_itself_moves...), under Transformation,
      Penalty and Lagrange, with -about moments, -ts and -Tn/-Sa spectra and
      a gravity preload;
  (b) the same on a residualModes set (with the -period override) and on a
      ritzVectors set, and the term / saPeriod / kind keys follow ADR 213 D3b;
  (c) bad groups, modes, directions and options are refused, and a refused
      call leaves the domain as it was.

Independent of the loop: for eigenmodes the total base reaction balances the
total modal inertia (iota^T K phi = 0 => R = -lambda iota^T M phi).
Model: tests/_testbed/rsa_tower.py (ASDShellQ4 walls, ElasticTimoshenkoBeam
frames, LadrunoMembrane floors, rigidDiaphragm, lateral lumped mass).
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops
from _testbed.rsa_tower import build_tower, spectrum

pytestmark = [pytest.mark.zone_a]

ROOT = Path(__file__).resolve().parent.parent
NST, H, NMODES = 3, 1.0, 6
TN, SA = spectrum()
SPEC = ['-Tn', *TN, '-Sa', *SA]
ABOUT = (15.0, 15.0, 0.0)
RTOL = 1e-12


@pytest.fixture(autouse=True)
def _quiet(monkeypatch):
    monkeypatch.setenv('LADRUNO_OPENSEES_QUIET', '1')


def _tower(handler='Transformation', preload=False):
    t = build_tower(NST, H, handler=handler)
    if preload:
        ops.timeSeries('Constant', 1)
        ops.pattern('Plain', 1, 1)
        for s, m in enumerate(t['masters'], start=1):
            for n in t['storey_nodes'][s]:
                if ops.nodeCoord(n)[2] == pytest.approx(s * 3.0):
                    ops.load(n, 0.0, 0.0, -1.0e4, 0.0, 0.0, 0.0)
        ops.integrator('LoadControl', 1.0)
        assert ops.analyze(1) == 0
        ops.loadConst('-time', 0.0)
    ops.eigen('-genBandArpack', NMODES)
    ops.modalProperties('-unorm')
    return t


def _groups(t):
    """(new-command args, oracle description) for a representative request."""
    core_base = [n for n in t['supports'] if 11.0 <= ops.nodeCoord(n)[0] <= 19.0
                 and 11.0 <= ops.nodeCoord(n)[1] <= 19.0]
    e1, c1 = t['piers'][1]
    e2, c2 = t['piers'][2]
    req = dict(sets={'core': (core_base, None)},
               eles={'pier1': (e1, c1, ABOUT), 'pier2': (e2, c2, None),
                     'pier2all': (e2, None, None)},
               inertia={'all': (sorted(set().union(*t['storey_nodes'].values())), None),
                        'st2': (t['storey_nodes'][2], ABOUT)})
    args = ['-reactions', '-about', *ABOUT]
    for name, (nodes, p) in req['sets'].items():
        args += ['-reactionSet', name, *nodes] + (['-about', *p] if p else [])
    for name, (eles, cut, p) in req['eles'].items():
        args += ['-eleGroup', name, '-ele', *eles] + (['-cut', *cut] if cut else []) \
            + (['-about', *p] if p else [])
    for name, (nodes, p) in req['inertia'].items():
        args += ['-inertia', name, *nodes] + (['-about', *p] if p else [])
    return args, req


def _moment(f, x, p):
    if p is None:
        return f
    f = f.copy()
    r = np.asarray(x) - np.asarray(p)
    f[3:6] += np.cross(r, f[0:3])
    return f


def _snapshot():
    """committed-state reactions, element forces and displacements"""
    ops.reactions()
    R = {n: np.array(ops.nodeReaction(n)) for n in ops.getFixedNodes()}
    F = {e: np.array(ops.eleForce(e)) for e in ops.getEleTags()}
    return R, F


def _loop(req, dirs, modes, spec_rsa, extra_calls=True):
    """today's per-mode loop through the stock commands.

    Returns ({key: [d][k][c]}, {key: [d][k][c]}): the values and the scale of
    the summed terms, sum |term| per component (the round-off floor of a sum
    that cancels, e.g. the vertical reactions of a core under overturning).
    For reactions the terms are the element force blocks at the support nodes
    (what calculateNodalReactions adds up: several shells meet at a base node
    and largely cancel), for element groups the blocks at the cut, for
    inertia the nodal forces. extra_calls=False skips every eleForce read
    (values only, no reaction scale), only to measure how much the call
    history alone moves the loop (ASDShellQ4, see test below)."""
    R0, F0 = _snapshot()
    # same normalization as _tower: -return re-stores the modal properties
    lam = np.array(ops.modalProperties('-unorm', '-return')['eigenLambda'])
    fixed = ops.getFixedNodes()
    # a reaction exists on the SP-constrained DOFs only
    mask = {n: np.isin(np.arange(1, 7), ops.getFixedDOFs(n)) for n in fixed}
    at = {n: [] for n in fixed}                       # support node -> (ele, block)
    for e in ops.getEleTags():
        for i, n in enumerate(ops.eleNodes(e)):
            if n in at:
                at[n].append((e, i))
    out, scl = {}, {}

    def put(key, d, a, terms, scale=None):
        terms = np.array(terms).reshape(-1, 6)
        out.setdefault(key, np.zeros((len(dirs), len(modes), 6)))[d, a] = terms.sum(axis=0)
        scl.setdefault(key, np.zeros((len(dirs), len(modes), 6)))[d, a] = (
            np.abs(terms).sum(axis=0) if scale is None else scale)

    def block_scale(F, nodes, p):
        """sum over the support nodes of |r x f| / |f| of every element block"""
        sc = np.zeros(6)
        for n in nodes:
            for e, i in at[n]:
                f = np.where(mask[n], F[e][6 * i:6 * i + 6], 0.0)
                sc += np.abs(_moment(f, ops.nodeCoord(n), p))
        return sc

    for di, d in enumerate(dirs):
        for a, k in enumerate(modes):
            ops.responseSpectrumAnalysis(*spec_rsa(d), '-mode', k, '-increment')
            ops.reactions()
            dR = {n: np.where(mask[n], np.array(ops.nodeReaction(n)) - R0[n], 0.0) for n in fixed}
            F = ({e: np.array(ops.eleForce(e)) - F0[e] for e in {e for n in fixed for e, _ in at[n]}}
                 if extra_calls else None)
            put('reactions', di, a, [_moment(dR[n], ops.nodeCoord(n), ABOUT) for n in fixed],
                block_scale(F, fixed, ABOUT) if extra_calls else None)
            for name, (nodes, p) in req['sets'].items():
                put('reaction:' + name, di, a, [_moment(dR[n], ops.nodeCoord(n), p) for n in nodes],
                    block_scale(F, nodes, p) if extra_calls else None)
            for name, (eles, cut, p) in (req['eles'].items() if extra_calls else ()):
                terms = []
                for e in eles:
                    f = np.array(ops.eleForce(e)) - F0[e]
                    for i, n in enumerate(ops.eleNodes(e)):
                        if cut is None or n in cut:
                            terms.append(_moment(f[6 * i:6 * i + 6], ops.nodeCoord(n), p))
                put('ele:' + name, di, a, terms)
            for name, (nodes, p) in req['inertia'].items():
                terms = []
                for n in nodes:
                    m = np.array(ops.nodeMass(n))
                    # nodeDisp is the COMMITTED displacement; -increment's trial
                    # minus committed is the IncrDisp response (4)
                    du = np.array([ops.nodeResponse(n, i, 4) for i in range(1, 7)])
                    terms.append(_moment(lam[k - 1] * m * du, ops.nodeCoord(n), p))
                put('inertia:' + name, di, a, terms)
    ops.responseSpectrumAnalysis('-revert')
    return out, scl


def _new(args, dirs, modes=None, spec=SPEC):
    extra = ['-modes', *modes] if modes else []
    out = ops.responseSpectrumPerMode(*spec, '-dir', *dirs, *extra, *args)
    nd, nm, nc = (int(v) for v in out['shape'])
    return {k: (np.array(v).reshape(nd, nm, nc) if ':' in k or k == 'reactions'
                else np.array(v)) for k, v in out.items()}


def _rel_err(got, ref, scale):
    """per (dir, mode): max over the components of |got - ref| / the term
    scale OF THAT COMPONENT (forces against forces, moments against moments,
    review of #1032). A component whose terms are all exactly zero must match
    exactly."""
    diff = np.abs(got - ref)
    with np.errstate(divide='ignore', invalid='ignore'):
        err = np.where(scale > 0.0, diff / scale, np.where(diff > 0.0, np.inf, 0.0))
    return err.max(axis=2)


def _assert_same(new, oracle, rtol=RTOL, over_modes=False):
    """new == loop to rtol x term scale, per (dir, mode) and component.
    over_modes=True takes each component's term scale as its largest over the
    direction's modes: the round-off a mode can inherit from the larger modes
    updated before it (ASDShellQ4, test_the_loop_itself_moves...)"""
    old, scl = oracle
    for key, ref in old.items():
        got = new[key]
        assert got.shape == ref.shape, key
        sc = scl[key]
        if over_modes:
            sc = np.broadcast_to(sc.max(axis=1, keepdims=True), sc.shape)
        err = _rel_err(got, ref, sc)
        d, a = np.unravel_index(np.argmax(err), err.shape)
        assert err.max() <= rtol, f'{key} dir {d} mode {a}: rel err {err[d, a]:.2e}'


def _rsa_lists(d):
    return [d, *SPEC]


# --------------------------------------------------------------------------- (a)
@pytest.mark.parametrize('handler', ['Transformation', 'Penalty', 'Lagrange'])
def test_a_equals_the_per_mode_loop(handler):
    t = _tower(handler)
    args, req = _groups(t)
    modes = list(range(1, NMODES + 1))
    new = _new(args, (1, 2))
    old = _loop(req, (1, 2), modes, _rsa_lists)
    _assert_same(new, old)
    assert list(new['modes']) == modes and list(new['dirs']) == [1, 2]
    # supports + pier shells only, not the whole domain
    assert new['numUpdated'][0] < 0.25 * t['n_ele']


def test_a_mode_subset_time_series_spectrum_and_scale():
    t = _tower()
    args, req = _groups(t)
    ops.timeSeries('Path', 77, '-time', *TN, '-values', *SA)
    new = _new(args, (2, 1), modes=[5, 2], spec=['-ts', 77, '-scale', 0.7])
    old = _loop(req, (2, 1), [5, 2], lambda d: [77, d, '-scale', 0.7])
    _assert_same(new, old)


def test_a_increment_over_a_gravity_preload():
    """with a committed gravity state the output is the loop minus that state"""
    t = _tower(preload=True)
    args, req = _groups(t)
    R0, _ = _snapshot()
    assert abs(sum(R0[n][2] for n in t['supports'])) > 1e5    # the preload is there
    new = _new(args, (1,))
    old = _loop(req, (1,), list(range(1, NMODES + 1)), _rsa_lists)
    # reactions/forces are differences of preloaded totals: cancellation
    _assert_same(new, old, rtol=1e-9)


def test_the_loop_itself_moves_with_its_call_history():
    """Why the gate measures against the term scale, not bit equality.
    ASDShellQ4's AGQI update starts from the KQU / Q_residual of the element's
    LAST calculateAll call, which revertToLastCommit does not restore, so the
    stock loop read with and without the eleForce calls differs at round-off.
    On a sum that cancels (the vertical reactions of a core under
    overturning) round-off is relative to the terms, not to the result."""
    t = _tower()
    _, req = _groups(t)
    modes = list(range(1, NMODES + 1))
    a = _loop(req, (1, 2), modes, _rsa_lists, extra_calls=True)
    b = _loop(req, (1, 2), modes, _rsa_lists, extra_calls=False)
    key = 'reaction:core'
    assert not np.array_equal(a[0][key], b[0][key])          # not bit-reproducible
    assert _rel_err(a[0][key], b[0][key], a[1][key]).max() < RTOL


def test_eigen_mode_reactions_balance_the_modal_inertia():
    """independent of the loop: iota^T K phi_k = 0, so per eigenmode the base
    reaction is minus the total modal inertia (force and moment)"""
    t = _tower()
    args, _ = _groups(t)
    new = _new(args, (1, 2))
    R = new['reactions']                       # about ABOUT
    I = _new(['-inertia', 'all', *sorted(set().union(*t['storey_nodes'].values())),
              '-about', *ABOUT], (1, 2))['inertia:all']
    scale = np.abs(I[..., :3]).max(axis=2, keepdims=True)
    assert np.abs(R[..., :3] + I[..., :3]).max() / scale.min() < 1e-6
    mscale = np.abs(I[..., 3:]).max()
    assert np.abs(R[..., 3:] + I[..., 3:]).max() / mscale < 1e-6


def test_domain_left_at_the_committed_state():
    t = _tower()
    args, _ = _groups(t)
    ops.responseSpectrumAnalysis(*_rsa_lists(1), '-mode', 2, '-increment')  # a trial left behind
    first = _new(args, (1,))
    n = t['storey_nodes'][NST][-1]
    assert all(ops.nodeResponse(n, i, 4) == 0.0 for i in range(1, 7))   # trial == committed
    again = _new(args, (1,))
    for k in first:
        assert np.array_equal(first[k], again[k]), k


# --------------------------------------------------------------------------- (b)
def test_b_residual_mode_set_with_period_override():
    t = _tower()
    ops.residualModes('-dir', 1, 2, '-period', 0.05)
    ops.modalProperties('-unorm')
    n = len(ops.modalProperties('-unorm', '-return')['eigenLambda'])
    assert n > NMODES
    args, req = _groups(t)
    modes = list(range(1, n + 1))
    new = _new(args, (1, 2))
    old = _loop(req, (1, 2), modes, _rsa_lists)
    # the stiff pseudo-modes (T < 1 ms here) follow modes of much larger
    # amplitude, whose round-off ASDShellQ4 carries into them (see
    # test_the_loop_itself_moves...): the loop run in reverse order differs
    # from itself by up to ~7e-12 of a pseudo-mode's OWN term scale, on Windows
    # and Linux alike. So on this set each component is measured against its
    # largest term scale over the direction's modes. (A bound of a few times
    # the loop's measured spread was tried first; one round-off sample is too
    # noisy an estimate of its size, and it failed on Linux.)
    _assert_same(new, old, over_modes=True)
    kind = new['kind']
    assert np.all(kind[:NMODES] == 0) and np.all(kind[NMODES:] == 1)
    # D3b: the overridden modes read Sa(0.05) and share one combination term
    assert np.all(new['saPeriod'][NMODES:] == 0.05)
    assert np.all(new['period'][NMODES:] != 0.05)
    assert len(set(new['term'][NMODES:])) == 1
    assert len(set(new['term'])) == NMODES + 1
    assert np.allclose(new['omega'][NMODES:], 2 * np.pi / 0.05)


def test_b_ritz_vector_set():
    t = build_tower(NST, H)
    ops.ritzVectors(NMODES, '-dir', 1, 2)
    ops.modalProperties('-unorm')
    args, req = _groups(t)
    new = _new(args, (1, 2))
    old = _loop(req, (1, 2), list(range(1, NMODES + 1)), _rsa_lists)
    _assert_same(new, old)
    assert np.all(new['kind'] == 2)
    assert np.array_equal(new['saPeriod'], new['period'])


# --------------------------------------------------------------------------- (c)
def test_c_refusals_leave_the_domain_alone():
    t = _tower()
    good = ['-dir', 1, '-reactions']
    ref = ops.responseSpectrumPerMode(*SPEC, *good)
    free = t['storey_nodes'][NST][-1]          # a top-floor node: no SP, on no pier
    e1, c1 = t['piers'][1]
    bad = [
        [*SPEC, '-reactions'],                                    # no -dir
        [*SPEC, '-dir', 7, '-reactions'],                         # dir out of range
        [*SPEC, '-dir', 1, 1, '-reactions'],                      # dir twice
        [*SPEC, '-dir', 1],                                       # nothing requested
        ['-dir', 1, '-reactions'],                                # no spectrum
        ['-Tn', 0.0, 1.0, '-Sa', 1.0, '-dir', 1, '-reactions'],   # length mismatch
        ['-Tn', 1.0, 0.5, '-Sa', 1.0, 1.0, '-dir', 1, '-reactions'],  # Tn not increasing
        ['-ts', 999, '-dir', 1, '-reactions'],                    # unknown series
        [*SPEC, '-dir', 1, '-modes', 0, '-reactions'],            # mode 0
        [*SPEC, '-dir', 1, '-modes', NMODES + 1, '-reactions'],   # mode out of range
        [*SPEC, '-dir', 1, '-modes', 2, 2, '-reactions'],         # mode twice
        [*SPEC, '-dir', 1, '-reactionSet', 'x', 99999],           # unknown node
        [*SPEC, '-dir', 1, '-reactionSet', 'x', free],            # not a supported node
        [*SPEC, '-dir', 1, '-reactionSet', 'x'],                  # empty set
        [*SPEC, '-dir', 1, '-eleGroup', 'p', '-ele', 99999],      # unknown element
        [*SPEC, '-dir', 1, '-eleGroup', 'p', '-ele', *e1, '-cut', free],  # cut not on group
        [*SPEC, '-dir', 1, '-eleGroup', 'p', 1, 2],               # missing -ele
        [*SPEC, '-dir', 1, '-eleGroup', 'p', '-ele', 1, 1],       # element twice
        [*SPEC, '-dir', 1, '-inertia', 'a', 1, '-inertia', 'a', 2],   # duplicate name
        [*SPEC, '-dir', 1, '-inertia', 'a', 99999],               # unknown node
        [*SPEC, '-dir', 1, '-about', 0.0, 0.0, 0.0, '-reactions'],    # -about first
        [*SPEC, '-dir', 1, '-reactions', '-about', 0.0, 0.0],     # 2 coords in 3D
        [*SPEC, '-dir', 1, '-reactions', '-bogus'],               # unknown option
    ]
    for args in bad:
        with pytest.raises(Exception):
            ops.responseSpectrumPerMode(*args)
    assert ops.responseSpectrumPerMode(*SPEC, *good) == ref
    # a new eigen without modalProperties is refused (stale modal properties)
    ops.eigen('-genBandArpack', NMODES - 1)
    with pytest.raises(Exception):
        ops.responseSpectrumPerMode(*SPEC, *good)


def test_c_non_positive_eigenvalue_is_refused():
    """a mechanism (unrestrained floor rotation) gives lambda ~ 0: refused"""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2, 0.0, 3.0); ops.mass(2, 1.0, 1.0, 0.0)
    ops.node(3, 1.0, 3.0); ops.mass(3, 1.0, 1.0, 0.0)
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', 1, 1, 2, 1.0, 1.0e3, 1.0, 1)
    ops.uniaxialMaterial('Elastic', 1, 1.0e3)
    ops.element('zeroLength', 2, 2, 3, '-mat', 1, '-dir', 1)
    ops.fix(3, 0, 1, 1)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    lam = ops.modalProperties('-return')['eigenLambda']
    bad = [k + 1 for k, v in enumerate(lam) if not v > 0.0]
    good = [k + 1 for k, v in enumerate(lam) if v > 0.0]
    if bad:
        with pytest.raises(Exception):
            ops.responseSpectrumPerMode(*SPEC, '-dir', 1, '-reactions')
    out = ops.responseSpectrumPerMode(*SPEC, '-dir', 1, '-modes', *good, '-reactions')
    assert [int(v) for v in out['shape']] == [1, len(good), 3]


# --------------------------------------------------------------------------- Tcl
def test_classic_tcl_bridge(tmp_path):
    exe = ROOT / 'dist' / 'bin' / 'OpenSees.exe'
    if not exe.exists():
        pytest.skip('dist/bin/OpenSees.exe not built')
    lines = ['model basic -ndm 2 -ndf 3',
             'node 1 0 0', 'fix 1 1 1 1', 'node 2 0 3', 'mass 2 2.0 0 0',
             'node 3 0 6', 'mass 3 1.0 0 0', 'geomTransf Linear 1',
             'element elasticBeamColumn 1 1 2 1.0 1000.0 1.0 1',
             'element elasticBeamColumn 2 2 3 1.0 1000.0 1.0 1',
             'fix 2 0 1 0', 'fix 3 0 1 0',
             'system UmfPack', 'analysis Static', 'eigen -fullGenLapack 2', 'modalProperties',
             # unexpanded Tcl lists, as responseSpectrumAnalysis takes them
             'set Tn {0.0 1.0 10.0}', 'set Sa {2.0 3.0 0.5}', 'set top {2 3}',
             'set r [responseSpectrumPerMode -Tn $Tn -Sa $Sa -dir 1 -reactions '
             '-eleGroup c1 -ele 1 -cut 1 -inertia top $top]',
             # no free mass in direction 2: refused, as -mode k -increment is
             'set bad [catch {responseSpectrumPerMode -Tn $Tn -Sa $Sa -dir 2 -reactions}]',
             'set out [open res.txt w]',
             'puts $out [dict get $r reactions]', 'puts $out [dict get $r ele:c1]',
             'puts $out [dict get $r inertia:top]', 'puts $out $bad', 'close $out']
    (tmp_path / 'deck.tcl').write_text('\n'.join(lines) + '\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET='1')
    subprocess.run([str(exe), 'deck.tcl'], cwd=tmp_path, env=env, check=True,
                   stdin=subprocess.DEVNULL, capture_output=True, timeout=120)
    lines = (tmp_path / 'res.txt').read_text().split('\n')
    got = [np.array(l.split(), float) for l in lines[:3]]
    assert lines[3].strip() == '1'                       # -dir 2 refused
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0, 0); ops.fix(1, 1, 1, 1)
    ops.node(2, 0, 3); ops.mass(2, 2.0, 0, 0)
    ops.node(3, 0, 6); ops.mass(3, 1.0, 0, 0)
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', 1, 1, 2, 1.0, 1000.0, 1.0, 1)
    ops.element('elasticBeamColumn', 2, 2, 3, 1.0, 1000.0, 1.0, 1)
    ops.fix(2, 0, 1, 0); ops.fix(3, 0, 1, 0)
    ops.system('UmfPack'); ops.analysis('Static')
    ops.eigen('-fullGenLapack', 2); ops.modalProperties()
    py = ops.responseSpectrumPerMode('-Tn', 0.0, 1.0, 10.0, '-Sa', 2.0, 3.0, 0.5, '-dir', 1,
                                     '-reactions', '-eleGroup', 'c1', '-ele', 1, '-cut', 1,
                                     '-inertia', 'top', 2, 3)
    # WP-212 guard per direction (review of #1032): no mass in Y -> refused,
    # alone or next to a good direction, like stock -mode k -increment
    spec = ['-Tn', 0.0, 1.0, 10.0, '-Sa', 2.0, 3.0, 0.5]
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(2, *spec, '-mode', 1, '-increment')
    for dirs in ((2,), (1, 2)):
        with pytest.raises(Exception):
            ops.responseSpectrumPerMode(*spec, '-dir', *dirs, '-reactions')
    for g, key in zip(got, ('reactions', 'ele:c1', 'inertia:top')):
        np.testing.assert_allclose(g, py[key], rtol=1e-12, atol=1e-12 * np.abs(py[key]).max())
    # 2D: 3 components, base shear = -inertia for the eigenmodes
    R = np.array(py['reactions']).reshape(2, 3)
    I = np.array(py['inertia:top']).reshape(2, 3)
    np.testing.assert_allclose(R[:, 0], -I[:, 0], rtol=1e-8)


# ------------------------------------------- TransformationFE > 64 nodes (found here)
_MANY_NODES = r'''
import numpy as np
from _testbed import ops
out = []
for handler in ('Transformation', 'Plain') * 4:
    ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
    pts = [(0, 0), (10, 0), (10, 10), (0, 10)]
    pts += [(1 + 8 * (i % 10) / 9.0, 1 + 8 * (i // 10) / 9.0) for i in range(100)]
    for k, (x, y) in enumerate(pts, start=1):
        ops.node(k, x, y, 0.0)
        ops.fix(k, int(k <= 4), int(k <= 4), 1, 1, 1, 1)
        ops.mass(k, 1.0, 1.0, 0, 0, 0, 0)
    ops.element('LadrunoMembrane', 1, '-nodes', 1, 2, 3, 4,
                '-diaphragmNodes', *range(5, len(pts) + 1), '-stiff', 30e9, 0.2, 0.2)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints(handler)
    ops.integrator('LoadControl', 0.0); ops.algorithm('Linear'); ops.analysis('Static')
    out.append(ops.eigen('-fullGenLapack', 3))
ops.wipe()
lam = np.array(out)
assert np.allclose(lam, lam[1], rtol=1e-10), lam
print('MANY-NODES-OK')
'''


def test_transformation_fe_element_with_more_than_64_nodes():
    """WP-216 found it: TransformationFE kept one count per NODE in a static
    buffer of MAX_NUM_DOF = 64 ints, so a LadrunoMembrane floor with more nodes
    corrupted the heap on every assembly under Transformation (0xc0000374 at
    the next wipe). 104-node membrane: Transformation = Plain, no crash.
    NOT a deterministic regression: heap corruption need not crash, so this
    can pass on an unfixed build (review of #1032); the eigenvalue match and
    the gate tower (which crashed every time before the fix) are the evidence."""
    from _testbed.subprocess_run import run_python_script
    rc, text = run_python_script(_MANY_NODES, cwd=str(Path(__file__).resolve().parent))
    assert rc == 0 and 'MANY-NODES-OK' in text, text[-2000:]
