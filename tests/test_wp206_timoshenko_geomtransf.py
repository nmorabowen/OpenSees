"""WP-206: ElasticTimoshenkoBeam3d honours its geomTransf (P-Delta, Corotational,
-jntOffset rigid arms).

The stock element built its own rotation-only Tgl from the transformation's axes
and took its geometric flag only from -geomNonlinear, so a PDelta transf gave no
P-Delta and -jntOffset rigid arms were silently dropped (YOO SLE port: frames
with joint offsets must keep them). It now works in the 6-dof basic system and
lets the CrdTransf map to global, like ElasticBeam3d.

Oracles:
  * shear-rigid limit (Avy, Avz -> huge): must equal elasticBeamColumn with the
    SAME transformation (offsets / PDelta / Corotational), to round-off;
  * finite shear area: closed-form Timoshenko cantilever (bending + shear);
  * rigid arms: the same frame modelled with explicit rigidLink beams;
  * -geomNonlinear with a PDelta transf equals the legacy Linear+-geomNonlinear
    path (no offsets), which carries the full consistent geometric matrix.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, G = 30.0e6, 12.0e6
A, Jx, Iy, Iz = 0.12, 2.0e-3, 9.0e-4, 1.6e-3
RIGID = 1.0e6          # shear area making Timoshenko -> Euler-Bernoulli
L = 4.0


def _model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)


def _static(load_fn, nsteps=1, algo='Newton'):
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    load_fn()
    ops.system('FullGeneral')
    ops.numberer('Plain')
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1.0e-12, 50)
    ops.algorithm(algo)
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    assert ops.analyze(nsteps) == 0


def _cantilever(ele, transf, offsets=None, av=RIGID, extra=()):
    """vertical cantilever along global Z, base fixed, tip node 2."""
    _model()
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    args = [transf, 1, 1.0, 0.0, 0.0]
    if offsets is not None:
        args += ['-jntOffset', *offsets]
    ops.geomTransf(*args)
    if ele == 'timo':
        ops.element('ElasticTimoshenkoBeam', 1, 1, 2, E, G, A, Jx, Iy, Iz, av, av, 1, *extra)
    else:
        ops.element('elasticBeamColumn', 1, 1, 2, A, E, G, Jx, Iy, Iz, 1, *extra)


def _tip(load, ele, transf, nsteps=1, **kw):
    _cantilever(ele, transf, **kw)
    _static(lambda: ops.load(2, *load), nsteps=nsteps)
    return [ops.nodeDisp(2, i) for i in range(1, 7)]


def _close(a, b, rel=1e-9, abs_=1e-14):
    for x, y in zip(a, b):
        assert x == pytest.approx(y, rel=rel, abs=abs_), (a, b)


# ---------------------------------------------------------------- linear path
def test_timoshenko_closed_form_cantilever_shear():
    """finite shear area: tip = PL^3/(3EI) + PL/(G Av), both bending planes."""
    Av = 0.05
    P = 10.0
    d = _tip((P, P, 0, 0, 0, 0), 'timo', 'Linear', av=Av)
    # which global direction bends about Iy vs Iz depends on vecxz; compare
    # the pair of tip deflections as a set
    ex = abs(d[0]); ey = abs(d[1])
    cand = sorted([P * L**3 / (3 * E * Iy) + P * L / (G * Av),
                   P * L**3 / (3 * E * Iz) + P * L / (G * Av)])
    assert sorted([ex, ey]) == pytest.approx(cand, rel=1e-10)


@pytest.mark.parametrize("transf", ['Linear', 'PDelta'])
def test_linear_loading_matches_elastic_beam_shear_rigid(transf):
    load = (3.0, -2.0, 5.0, 0.7, -0.4, 0.3)
    _close(_tip(load, 'timo', transf), _tip(load, 'ebc', transf), rel=1e-8)


# ---------------------------------------------------------------- joint offsets
@pytest.mark.parametrize("transf", ['Linear', 'PDelta'])
def test_joint_offsets_match_elastic_beam(transf):
    off = (0.2, -0.1, 0.3, -0.15, 0.25, -0.4)
    load = (3.0, -2.0, 5.0, 0.7, -0.4, 0.3)
    t = _tip(load, 'timo', transf, offsets=off)
    e = _tip(load, 'ebc', transf, offsets=off)
    _close(t, e, rel=1e-8)
    # and the offsets matter (the stock element ignored them)
    t0 = _tip(load, 'timo', transf)
    assert abs(t[0] - t0[0]) > 1e-3 * abs(t0[0])


def test_joint_offsets_match_explicit_rigid_links():
    """rigid arms == rigidLink beams to two extra nodes, finite shear area."""
    Av = 0.04
    offI = (0.0, 0.0, 0.5)
    offJ = (0.3, 0.0, -0.6)
    load = (4.0, 3.0, -2.0, 0.2, 0.1, -0.3)
    # offsets
    _cantilever('timo', 'Linear', offsets=(*offI, *offJ), av=Av)
    _static(lambda: ops.load(2, *load))
    d_off = [ops.nodeDisp(2, i) for i in range(1, 7)]
    # explicit: flexible element between the arm ends, rigid links to the nodes
    _model()
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.node(11, *offI)
    ops.node(12, offJ[0], offJ[1], L + offJ[2])
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.rigidLink('beam', 1, 11)
    ops.rigidLink('beam', 2, 12)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    ops.element('ElasticTimoshenkoBeam', 1, 11, 12, E, G, A, Jx, Iy, Iz, Av, Av, 1)
    _static(lambda: ops.load(2, *load))
    d_exp = [ops.nodeDisp(2, i) for i in range(1, 7)]
    _close(d_off, d_exp, rel=1e-7)


# ---------------------------------------------------------------- P-Delta
def _pdelta_tip(ele, transf, P, H=1.0, extra=(), offsets=None):
    _cantilever(ele, transf, extra=extra, offsets=offsets)
    _static(lambda: ops.load(2, H, 0.0, -P, 0, 0, 0), nsteps=1)
    return ops.nodeDisp(2, 1)


def test_pdelta_transf_now_softens():
    P = 0.4 * math.pi**2 * E * min(Iy, Iz) / (4 * L * L)
    u0 = _pdelta_tip('timo', 'Linear', P)
    up = _pdelta_tip('timo', 'PDelta', P)
    ue = _pdelta_tip('ebc', 'PDelta', P)
    assert up > 1.2 * u0                      # stock: up == u0 (ignored)
    assert up == pytest.approx(ue, rel=1e-8)  # same chord P-Delta as ElasticBeam3d


def test_pdelta_with_offsets_matches_elastic_beam():
    P = 0.3 * math.pi**2 * E * min(Iy, Iz) / (4 * L * L)
    off = (0.0, 0.0, 0.3, 0.0, 0.0, -0.4)
    up = _pdelta_tip('timo', 'PDelta', P, offsets=off)
    ue = _pdelta_tip('ebc', 'PDelta', P, offsets=off)
    assert up == pytest.approx(ue, rel=1e-8)


def test_geomnonlinear_pdelta_equals_legacy_linear_geomnonlinear():
    """chord (PDelta transf) + element P-delta (kgb) == legacy full klgeo."""
    P = 0.3 * math.pi**2 * E * min(Iy, Iz) / (4 * L * L)
    legacy = _pdelta_tip('timo', 'Linear', P, extra=('-geomNonlinear',))
    new = _pdelta_tip('timo', 'PDelta', P, extra=('-geomNonlinear',))
    assert new == pytest.approx(legacy, rel=1e-9)


def test_pdelta_eigen_period_ratio():
    """cantilever + tip mass + axial load: T(P)/T(0) ~ 1/sqrt(1 - P/Pcr_chord)."""
    m = 50.0
    def period(P):
        _cantilever('timo', 'PDelta')
        ops.mass(2, m, m, 0.0, 0.0, 0.0, 0.0)
        _static(lambda: ops.load(2, 0, 0, -P, 0, 0, 0))
        ops.loadConst('-time', 0.0)
        lam = ops.eigen('-fullGenLapack', 1)[0]
        return 2 * math.pi / math.sqrt(lam)
    k = 3 * E * min(Iy, Iz) / L**3   # the softest direction is mode 1
    # one element, PDelta chord term on the translations only: k(P) = k - P/L
    P = 0.25 * k * L                 # chord stiffness P/L = 0.25 k
    ratio = period(P) / period(0.0)
    assert ratio == pytest.approx(1.0 / math.sqrt(1.0 - 0.25), rel=1e-6)


# ---------------------------------------------------------------- corotational
def test_corotational_large_rotation_matches_elastic_beam():
    load = (0.0, 0.0, 0.0, 0.0, 0.25 * E * Iy / L, 0.0)   # ~quarter-radian tip rotation
    t = _tip(load, 'timo', 'Corotational', nsteps=10)
    e = _tip(load, 'ebc', 'Corotational', nsteps=10)
    _close(t, e, rel=1e-7, abs_=1e-12)
    assert max(abs(x) for x in t[3:]) > 0.1


# ---------------------------------------------------------------- element loads
@pytest.mark.parametrize("transf,offsets", [('Linear', None), ('PDelta', None),
                                            ('Linear', (0.0, 0.0, 0.2, 0.0, 0.0, -0.3))])
def test_uniform_and_partial_loads_match_elastic_beam(transf, offsets):
    def run(ele):
        _cantilever(ele, transf, offsets=offsets)
        def ld():
            ops.eleLoad('-ele', 1, '-type', '-beamUniform', 2.0, -1.5, 0.8)
            ops.eleLoad('-ele', 1, '-type', '-beamUniform', 1.0, 0.5, 0.0, 0.2, 0.7)
        _static(ld)
        ops.reactions()
        return [ops.nodeDisp(2, i) for i in range(1, 7)] + \
               [ops.nodeReaction(1, i) for i in range(1, 7)]
    _close(run('timo'), run('ebc'), rel=1e-8, abs_=1e-10)


# ---------------------------------------------------------------- responses
def test_local_and_basic_forces():
    off = (0.0, 0.0, 0.2, 0.0, 0.0, -0.3)
    load = (3.0, -2.0, 5.0, 0.7, -0.4, 0.3)
    _cantilever('timo', 'Linear', offsets=off)
    _static(lambda: ops.load(2, *load))
    lt = ops.eleResponse(1, 'localForce')
    bt = ops.eleResponse(1, 'basicForce')
    _cantilever('ebc', 'Linear', offsets=off)
    _static(lambda: ops.load(2, *load))
    le = ops.eleResponse(1, 'localForce')
    be = ops.eleResponse(1, 'basicForce')
    _close(lt, le, rel=1e-8, abs_=1e-9)
    assert len(bt) == 6
    _close(bt, be, rel=1e-8, abs_=1e-9)


def test_save_restore_keeps_flag_and_path(tmp_path):
    P = 0.3 * math.pi**2 * E * min(Iy, Iz) / (4 * L * L)
    _cantilever('timo', 'PDelta', extra=('-geomNonlinear',))
    ops.database('File', str(tmp_path / 'timo'))
    ops.save(1)
    ref = _pdelta_tip('timo', 'PDelta', P, extra=('-geomNonlinear',))
    _cantilever('timo', 'PDelta', extra=('-geomNonlinear',))
    ops.database('File', str(tmp_path / 'timo'))
    ops.restore(1)
    _static(lambda: ops.load(2, 1.0, 0.0, -P, 0, 0, 0))
    assert ops.nodeDisp(2, 1) == pytest.approx(ref, rel=1e-9)
