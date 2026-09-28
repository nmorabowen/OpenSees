"""PlateRebarMaterial response forwarding to its wrapped uniaxial (WP-142 P1b, oracle O6).

PlateRebarMaterial (vanilla) had no setResponse, so the wrapped bar's own state (plastic
strain, back stress, damage, buckling state) was unreachable inside a LayeredShell. The fork
adds (``// Ladruno`` in PlateRebarMaterial.cpp):

  * ``material <args...>`` -> the wrapped uniaxial's ``setResponse(<args...>)``, for ANY key,
    so the keys that collide with the NDMaterial base reach the bar:
    ``section k fiber j material stress`` is the scalar bar stress;
  * keys the NDMaterial base answers (``stress``, ``strain``, ``tangent``, ...) keep their
    5-component plate response, unchanged;
  * any other key is forwarded unchanged (``plasticStrain``, ``backStress``, ``PLE``, ...).

Oracle (T3, identity, bit-exact). One ASDShellQ4 on a LayeredShell (elastic concrete + bars
at 0 and 30 degrees wrapping a plastic uniaxial) is driven through a non-proportional cyclic
membrane history by prescribed displacements (Penalty). For every Gauss point and bar, the
bar strain eps11 c^2 + eps22 s^2 + gamma12 c s of the layer's committed plate strain drives a
standalone copy of the same uniaxial (a unit Truss whose end displacement is set exactly).
Every forwarded response must equal the standalone one bit for bit, and the base
``stress`` must still be the 5-component plate stress built from the bar stress.

Two drivers. ``Linear`` (one trial strain per step): the standalone is driven by the formula
strain itself, and the forwarded bar strain must equal the formula bit for bit. ``Newton``
(several trials per step): many uniaxials ignore a trial change below an ABSOLUTE dead-band
(LadrunoUniaxialJ2: ``fabs(Tstrain - strain) < DBL_EPSILON`` -> return), so the committed bar
state can belong to an earlier iterate that differs from the final plate strain by < 2.2e-16.
There the standalone replays the bar's own forwarded strain, and the formula is held to that
dead-band.

The forwarded responses are read twice: by ``eleResponse`` (a fresh Response each step) and
by Element recorders created BEFORE the history. Only the recorders keep one Response object
alive across steps, so only they catch a Response bound to a stale or foreign uniaxial copy
(named mutant "forward to a stale copy": ``theMat->getCopy()->setResponse(...)``).
"""
import itertools
import math
import os
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E_C, NU_C = 30000.0, 0.2
E_S = 200000.0
T_C = 0.06
T_S = 0.002
BAR_FIBERS = {2: 0.0, 4: 30.0}               # fiber index -> PlateRebar angle
LAYERS = [(1, T_C), (20, T_S), (1, T_C), (30, T_S), (1, T_C)]
GPS = (1, 2, 3, 4)
NODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}

# non-proportional cyclic membrane history (exx, eyy, gxy); both bars yield both ways
STATES = [(0.0, 0.0, 0.0), (4.0e-3, 1.0e-3, 3.0e-3), (-3.0e-3, 2.0e-3, -4.0e-3),
          (2.0e-3, -3.0e-3, 5.0e-3), (-1.0e-3, 0.0, -2.0e-3)]
NSUB = 5

# (key path, standalone key path) for the scalar responses forwarded to the bar
FWD_COMMON = [
    (('material', 'stress'), ('stress',)),
    (('material', 'strain'), ('strain',)),
    (('material', 'tangent'), ('tangent',)),
]


def _j2(tag):
    # Voce + linear isotropic, one Armstrong-Frederick back stress; yield strain 2e-3
    ops.uniaxialMaterial('LadrunoUniaxialJ2', tag, E_S, '-iso', 'voce', 400.0, 100.0, 20.0,
                         1000.0, '-kin', 1, 20000.0, 150.0)


def _asdsteel(tag):
    ops.uniaxialMaterial('ASDSteel1D', tag, E_S, 400.0, 540.0, 0.1)


MATS = {
    'LadrunoUniaxialJ2': (_j2, FWD_COMMON + [
        (('plasticStrain',), ('plasticStrain',)),
        (('backStress',), ('backStress',)),
        (('equivalentPlasticStrain',), ('equivalentPlasticStrain',)),
        (('material', 'plasticStrain'), ('plasticStrain',)),
    ]),
    'ASDSteel1D': (_asdsteel, FWD_COMMON + [
        (('PLE',), ('PLE',)),
        (('Damage',), ('Damage',)),
        (('material', 'PLE'), ('PLE',)),
    ]),
}


def _cs(angle):
    # the same expression, in the same order, as the PlateRebarMaterial constructor
    rang = angle * 4.0 * math.asin(1.0) / 360.0
    return math.cos(rang), math.sin(rang)


def _bar_strain(angle, e):
    # the same branches and operation order as PlateRebarMaterial::setTrialStrain
    if angle == 0:
        return e[0]
    if angle == 90:
        return e[1]
    c, s = _cs(angle)
    return e[0] * c * c + e[1] * s * s + e[2] * c * s


def _plate_stress(angle, sig):
    # PlateRebarMaterial::getStress
    out = [0.0] * 5
    if angle == 0:
        out[0] = sig
    elif angle == 90:
        out[1] = sig
    else:
        c, s = _cs(angle)
        out[0], out[1], out[2] = sig * c * c, sig * s * s, sig * c * s
    return out


def _path():
    out = []
    for a, b in zip(STATES[:-1], STATES[1:]):
        for k in range(1, NSUB + 1):
            out.append(tuple(a[i] + (b[i] - a[i]) * k / NSUB for i in range(3)))
    return out


def _impose(exx, eyy, gxy, first):
    if not first:
        ops.remove('loadPattern', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y) in NODES.items():
        ops.sp(n, 1, exx * x + 0.5 * gxy * y)       # symmetric gradient: no drilling
        ops.sp(n, 2, 0.5 * gxy * x + eyy * y)


def _resp(*args):
    return [float(v) for v in ops.eleResponse(*args)]


def run_shell(mat_fn, fwd, rec_dir, algo):
    """Drive the shell; return per-step readings + the recorder file map."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for n, (x, y) in NODES.items():
        ops.node(n, x, y, 0.0)
        ops.fix(n, 0, 0, 1, 1, 1, 1)
    mat_fn(11)
    ops.nDMaterial('ElasticIsotropic', 1, E_C, NU_C)
    ops.nDMaterial('PlateRebar', 20, 11, BAR_FIBERS[2])
    ops.nDMaterial('PlateRebar', 30, 11, BAR_FIBERS[4])
    ops.section('LayeredShell', 10, len(LAYERS), *[v for lay in LAYERS for v in lay])
    ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 10)

    files = {}
    for gp in GPS:
        for fib in BAR_FIBERS:
            keys = [k for k, _ in fwd] + [('stress',), ('strain',)]
            for key in keys:
                fn = os.path.join(rec_dir, f'g{gp}_f{fib}_{"_".join(key)}.out')
                files[(gp, fib, key)] = fn
                ops.recorder('Element', '-file', fn, '-precision', 17, '-ele', 1,
                             'section', str(gp), 'fiber', str(fib), *key)

    ops.timeSeries('Constant', 1)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-14, 50, 0)
    ops.algorithm(algo)          # Linear: exactly one trial strain per committed step
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')

    steps = []
    for i, (exx, eyy, gxy) in enumerate(_path()):
        _impose(exx, eyy, gxy, first=(i == 0))
        assert ops.analyze(1) == 0, f'shell step {i} failed'
        rd = {}
        for gp in GPS:
            for fib in BAR_FIBERS:
                pre = (1, 'section', str(gp), 'fiber', str(fib))
                rd[(gp, fib, 'plate_strain')] = _resp(*pre, 'strain')
                rd[(gp, fib, 'plate_stress')] = _resp(*pre, 'stress')
                for key, _ in fwd:
                    rd[(gp, fib, key)] = _resp(*pre, *key)
        steps.append(rd)
    ops.remove('recorders')
    return steps, files


def run_standalone(mat_fn, fwd, bar_strains):
    """bar_strains[step][(gp, fib)] -> exact strain; returns per-step standalone readings."""
    pairs = sorted(bar_strains[0])
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    mat_fn(11)
    ops.uniaxialMaterial('Elastic', 99, 1.0)
    ops.node(1, 0.0)
    ops.fix(1, 1)
    ops.node(2, 1.0)                       # free dummy equation: keeps the SOE non-empty
    ops.element('Truss', 1000, 1, 2, 1.0, 99)
    ele = {}
    for k, pair in enumerate(pairs):
        a, b = 10 + 2 * k, 11 + 2 * k
        ops.node(a, 0.0)
        ops.node(b, 1.0)
        ops.fix(a, 1)
        ops.fix(b, 1)                      # displacement set exactly by setNodeDisp
        ops.element('Truss', k + 1, a, b, 1.0, 11)   # A = 1, L = 1 -> strain = u_b
        ele[pair] = (k + 1, b)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-14, 10, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')

    out = []
    for step in bar_strains:
        for pair, (e, b) in ele.items():
            ops.setNodeDisp(b, 1, step[pair], '-commit')
        assert ops.analyze(1) == 0
        rd = {}
        for pair, (e, b) in ele.items():
            for key, skey in fwd:
                rd[(pair[0], pair[1], key)] = _resp(e, 'material', *skey)
        out.append(rd)
    return out


def _read(fn):
    with open(fn) as f:
        return [[float(v) for v in line.split()] for line in f if line.strip()]


@pytest.fixture(scope='module', params=list(itertools.product(MATS, ('Linear', 'Newton'))),
                ids=lambda p: f'{p[0]}-{p[1]}')
def histories(request, tmp_path_factory):
    name, algo = request.param
    mat_fn, fwd = MATS[name]
    rec_dir = str(tmp_path_factory.mktemp(f'rec_{name}_{algo}'))
    shell, files = run_shell(mat_fn, fwd, rec_dir, algo)
    formula = [{(gp, fib): _bar_strain(BAR_FIBERS[fib], rd[(gp, fib, 'plate_strain')])
                for gp in GPS for fib in BAR_FIBERS} for rd in shell]
    if algo == 'Linear':
        drive = formula
    else:   # the bar's own committed strain (see the dead-band note in the module docstring)
        drive = [{pair: rd[(*pair, ('material', 'strain'))][0] for pair in f}
                 for rd, f in zip(shell, formula)]
    alone = run_standalone(mat_fn, fwd, drive)
    recs = {key: _read(fn) for key, fn in files.items()}
    return dict(name=name, algo=algo, fwd=fwd, shell=shell, formula=formula, alone=alone,
                recs=recs)


def test_history_is_plastic_both_ways(histories):
    # the gate is only meaningful if both bars yield in tension and compression
    for fib in BAR_FIBERS:
        sig = [rd[(1, fib, ('material', 'stress'))][0] for rd in histories['shell']]
        assert max(sig) > 400.0 and min(sig) < -400.0, (histories['name'], fib, max(sig), min(sig))


def test_forwarded_bar_strain_is_the_rotated_plate_strain(histories):
    tol = 0.0 if histories['algo'] == 'Linear' else sys.float_info.epsilon
    for rd, formula in zip(histories['shell'], histories['formula']):
        for (gp, fib), e in formula.items():
            got = rd[(gp, fib, ('material', 'strain'))]
            assert len(got) == 1 and abs(got[0] - e) <= tol, (histories['algo'], gp, fib, got, e)


def test_forwarded_responses_equal_standalone_bit_exact(histories):
    n = 0
    for i, (rd, al) in enumerate(zip(histories['shell'], histories['alone'])):
        for key, got in al.items():
            assert rd[key] == got, (histories['name'], i, key, rd[key], got)
            n += 1
    assert n == len(histories['shell']) * len(GPS) * len(BAR_FIBERS) * len(histories['fwd'])


def test_recorders_bind_to_the_live_copy(histories):
    # one Response object per recorder, alive across the whole history
    nstep = len(histories['alone'])
    for (gp, fib, key), rows in histories['recs'].items():
        if key in (('stress',), ('strain',)):
            continue
        assert len(rows) == nstep, (key, len(rows))
        for i, row in enumerate(rows):
            assert row == histories['alone'][i][(gp, fib, key)], (histories['name'], i, gp, fib, key)


def test_base_plate_stress_is_unchanged(histories):
    nstep = len(histories['shell'])
    for i, rd in enumerate(histories['shell']):
        for gp in GPS:
            for fib, ang in BAR_FIBERS.items():
                sig = rd[(gp, fib, ('material', 'stress'))][0]
                want = _plate_stress(ang, sig)
                assert rd[(gp, fib, 'plate_stress')] == want
                assert histories['recs'][(gp, fib, ('stress',))][i] == want
                assert histories['recs'][(gp, fib, ('strain',))][i] == rd[(gp, fib, 'plate_strain')]
                assert len(rd[(gp, fib, 'plate_strain')]) == 5
    assert nstep == len(STATES[1:]) * NSUB


def test_reachability_and_null_keys():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for n, (x, y) in NODES.items():
        ops.node(n, x, y, 0.0)
    _j2(11)
    ops.uniaxialMaterial('LadrunoRebarBuckling', 12, 11, '-lsr', 8.0, '-fy', 400.0, '-E', E_S)
    ops.nDMaterial('ElasticIsotropic', 1, E_C, NU_C)
    ops.nDMaterial('PlateRebar', 20, 11, 0.0)
    ops.nDMaterial('PlateRebar', 40, 12, 90.0)
    ops.section('LayeredShell', 10, 3, 1, T_C, 20, T_S, 40, T_S)
    ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 10)
    for head in ('section', 'material'):          # ASDShellQ4 takes either keyword
        pre = (1, head, '1', 'fiber', '2')
        assert len(ops.eleResponse(*pre, 'stress')) == 5
        assert len(ops.eleResponse(*pre, 'strain')) == 5
        assert len(ops.eleResponse(*pre, 'material', 'stress')) == 1
        assert len(ops.eleResponse(*pre, 'plasticStrain')) == 1
        assert list(ops.eleResponse(*pre, 'noSuchResponseKey')) == []
        assert list(ops.eleResponse(*pre, 'material')) == []
        assert list(ops.eleResponse(*pre, 'material', 'noSuchResponseKey')) == []
    # a wrapper as the bar: LadrunoRebarBuckling's own state, and a key it forwards to its J2
    pre = (1, 'section', '1', 'fiber', '3')
    assert len(ops.eleResponse(*pre, 'stress')) == 5
    assert len(ops.eleResponse(*pre, 'buckling')) == 1
    assert len(ops.eleResponse(*pre, 'material', 'reduction')) == 1
    assert len(ops.eleResponse(*pre, 'plasticStrain')) == 1


# ---------------------------------------------------------------------------
# recorder metadata, MPCO end-to-end, database round trip
# ---------------------------------------------------------------------------
def _small_shell(angles=(0.0,)):
    """ASDShellQ4 + [concrete, bar(angles[0]), concrete, bar(angles[1]), concrete ...]."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for n, (x, y) in NODES.items():
        ops.node(n, x, y, 0.0)
        ops.fix(n, 0, 0, 1, 1, 1, 1)
    _j2(11)
    ops.nDMaterial('ElasticIsotropic', 1, E_C, NU_C)
    layers = [(1, T_C)]
    for k, ang in enumerate(angles):
        ops.nDMaterial('PlateRebar', 20 + k, 11, ang)
        layers += [(20 + k, T_S), (1, T_C)]
    ops.section('LayeredShell', 10, len(layers), *[v for lay in layers for v in lay])
    ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 10)
    ops.timeSeries('Constant', 1)
    return layers


def _analysis(algo='Newton'):
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-14, 50, 0)
    ops.algorithm(algo)
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def test_recorder_metadata_is_one_material_block_per_fiber(tmp_path):
    """-xml element recorder metadata: each fiber carries exactly ONE material block.

    NDMaterial::setResponse writes an NdMaterialOutput tag even when it returns null, so a
    plain "base first, then forward" left an EMPTY NdMaterialOutput beside the bar's own
    block. PlateRebar probes the base on a silent stream and wraps a forwarded request in
    one NdMaterialOutput enclosing the bar's block. The wrapper stays even when the bar
    emits no tags (LadrunoUniaxialJ2 plasticStrain): MPCO needs a material node per fiber,
    see test_mpco_section_fiber_forwarded_key.
    """
    import xml.etree.ElementTree as ET
    _small_shell()
    pre = ('-ele', 1, 'section', '1', 'fiber', '2')
    cases = {'fwd': ('stressStrain',), 'prefixed': ('material', 'stress'),
             'tagless': ('plasticStrain',), 'base': ('stress',)}
    files = {}
    for name, key in cases.items():
        files[name] = str(tmp_path / f'{name}.xml')
        ops.recorder('Element', '-xml', files[name], *pre, *key)
    _impose(1.0e-3, 0.0, 0.0, first=True)
    _analysis()
    assert ops.analyze(1) == 0
    ops.remove('recorders')
    ops.wipe()
    kids = {}
    for name, fn in files.items():
        (fiber,) = list(ET.parse(fn).getroot().iter('FiberOutput'))
        kids[name] = [(c.tag, [g.tag for g in c]) for c in fiber]
    assert kids['fwd'] == [('NdMaterialOutput', ['UniaxialMaterialOutput'])], kids['fwd']
    assert kids['prefixed'] == [('NdMaterialOutput', ['UniaxialMaterialOutput'])], kids['prefixed']
    assert kids['tagless'] == [('NdMaterialOutput', [])], kids['tagless']
    (base,) = kids['base']
    assert base[0] == 'NdMaterialOutput' and base[1] and set(base[1]) == {'ResponseType'}, base


def test_mpco_section_fiber_forwarded_key(tmp_path):
    h5py = pytest.importorskip('h5py')
    os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')
    layers = _small_shell(angles=(0.0, 30.0))
    nfib = len(layers)
    bars = [i + 1 for i, (m, _) in enumerate(layers) if m != 1]
    fn = str(tmp_path / 'o6.mpco')
    ops.recorder('mpco', fn, '-E', 'section.fiber.stress', 'section.fiber.plasticStrain')
    _impose(4.0e-3, 1.0e-3, 3.0e-3, first=True)
    _analysis()
    assert ops.analyze(1) == 0
    want_pl, want_st = [], []
    for gp in GPS:
        for fib in bars:
            want_pl += _resp(1, 'section', str(gp), 'fiber', str(fib), 'plasticStrain')
        for fib in range(1, nfib + 1):
            want_st += _resp(1, 'section', str(gp), 'fiber', str(fib), 'stress')
    assert max(abs(v) for v in want_pl) > 1.0e-4          # the bars did yield
    ops.remove('recorders')
    ops.wipe()

    def bucket(f, result):
        grp = f[f'MODEL_STAGE[1]/RESULTS/ON_ELEMENTS/{result}']
        (name,) = list(grp.keys())                          # one element class -> one bucket
        return grp[name]

    def step(b):
        # one analysis step -> one dataset, but its STEP_<n> name is not always STEP_0
        # (in the full module run, after other analyses, STEP_0 was absent): read the one
        (ds,) = list(b['DATA'].values())
        return [float(v) for v in ds[0]]

    with h5py.File(fn, 'r') as f:
        pl = bucket(f, 'section.fiber.plasticStrain')
        # the rebar-layer bucket: per Gauss point only the bar fibers (the concrete layers
        # answer nothing), one scalar each. LadrunoUniaxialJ2 emits no ResponseType tag for
        # this key, so MPCO groups a point's bars as ONE entry of len(bars) components
        # ("C1,C2"); the product is what is fixed.
        assert int(pl.attrs['NUM_COLUMNS'][0]) == len(GPS) * len(bars)
        mult = [int(v) for v in pl['META/MULTIPLICITY'][:].ravel()]
        ncmp = [int(v) for v in pl['META/NUM_COMPONENTS'][:].ravel()]
        assert [m * c for m, c in zip(mult, ncmp)] == [len(bars)] * len(GPS), (mult, ncmp)
        assert step(pl) == want_pl
        st = bucket(f, 'section.fiber.stress')
        # the plate-stress bucket: every layer, the 5-component plate stress
        assert int(st.attrs['NUM_COLUMNS'][0]) == len(GPS) * nfib * 5
        assert list(st['META/MULTIPLICITY'][:].ravel()) == [nfib] * len(GPS)
        assert list(st['META/NUM_COMPONENTS'][:].ravel()) == [5] * len(GPS)
        assert step(st) == want_st


RT_ANGLES = (30.0, 17.0)                   # off-axis: exercises c, s (0/90 take shortcuts)
RT_PATH = [(2.0e-3, 5.0e-4, 1.5e-3), (4.0e-3, 1.0e-3, 3.0e-3), (1.0e-3, 2.0e-3, -2.0e-3)]
RT_NEXT = (-2.0e-3, 1.0e-3, -4.0e-3)


def _rt_shell(angles):
    """Every shell DOF fixed and set exactly by setNodeDisp: the bar strains are then pure
    functions of the imposed displacements, with no global solve in between. (A Penalty
    solve has condition number ~alpha/K, so it amplifies any ulp-level difference in the
    restored state -- e.g. LadrunoUniaxialJ2::recvSelf re-derives Tstress, 1 ulp off -- to
    ~1e-9 and would hide the angle bug this test is for.) One free dummy equation keeps the
    system non-empty. EAS off: this gates the material's send/recv, not the element's.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for n, (x, y) in NODES.items():
        ops.node(n, x, y, 0.0)
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    ops.node(98, 0.0, 0.0, 1.0)
    ops.node(99, 1.0, 0.0, 1.0)
    ops.fix(98, 1, 1, 1, 1, 1, 1)
    ops.fix(99, 0, 1, 1, 1, 1, 1)
    _j2(11)
    ops.uniaxialMaterial('Elastic', 99, 1.0)
    ops.element('Truss', 99, 98, 99, 1.0, 99)
    ops.nDMaterial('ElasticIsotropic', 1, E_C, NU_C)
    layers = [(1, T_C)]
    for k, ang in enumerate(angles):
        ops.nDMaterial('PlateRebar', 20 + k, 11, ang)
        layers += [(20 + k, T_S), (1, T_C)]
    ops.section('LayeredShell', 10, len(layers), *[v for lay in layers for v in lay])
    ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 10, '-noeas')


def _rt_analysis():
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _rt_step(exx, eyy, gxy):
    for n, (x, y) in NODES.items():
        ops.setNodeDisp(n, 1, exx * x + 0.5 * gxy * y, '-commit')
        ops.setNodeDisp(n, 2, 0.5 * gxy * x + eyy * y, '-commit')
    assert ops.analyze(1) == 0


def _rt_read(state_only=False):
    out = []
    for gp in GPS:
        for fib in (2, 4):
            pre = (1, 'section', str(gp), 'fiber', str(fib))
            out += _resp(*pre, 'material', 'strain') + _resp(*pre, 'plasticStrain')
            if not state_only:
                out += _resp(*pre, 'material', 'stress') + _resp(*pre, 'stress')
    return out


def test_database_roundtrip_off_axis_bar_is_bit_identical(tmp_path):
    """save/wipe/restore (FE_Datastore File), then one more step == the uninterrupted run.

    PlateRebarMaterial::recvSelf used to rebuild (c, s) with angle * 0.0174532925 while the
    constructor uses angle * 4 asin(1)/360, so a restored off-axis bar saw a different bar
    strain (~1e-9 relative) from the first step after the restore.
    """
    # A: uninterrupted
    _rt_shell(RT_ANGLES)
    _rt_analysis()
    for e in RT_PATH + [RT_NEXT]:
        _rt_step(*e)
    ref = _rt_read()
    assert max(abs(v) for v in ref[1::8]) > 1.0e-4        # the bars carry plastic strain
    ops.wipe()

    # B: same history, save after RT_PATH, wipe, restore, then RT_NEXT
    _rt_shell(RT_ANGLES)
    _rt_analysis()
    for e in RT_PATH:
        _rt_step(*e)
    before = _rt_read(state_only=True)
    db = str(tmp_path / 'o6rt')
    try:
        ops.database('File', db)
    except Exception as exc:  # noqa: BLE001 - build without FE_Datastore
        pytest.skip(f'database() unsupported in this build: {exc}')
    ops.save(1)
    ops.wipe()
    _rt_shell(RT_ANGLES)                   # skeleton for the restore to land on
    ops.database('File', db)
    ops.restore(1)
    assert _rt_read(state_only=True) == before     # committed bar state came back exactly
    _rt_analysis()
    _rt_step(*RT_NEXT)
    got = _rt_read()
    ops.wipe()
    assert got == ref
