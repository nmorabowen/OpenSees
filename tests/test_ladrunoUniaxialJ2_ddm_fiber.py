"""WP-262: DDM sensitivities of LadrunoUniaxialJ2 inside a fiber-section beam.

Before WP-262 LadrunoUniaxialJ2 inherited UniaxialMaterial's sensitivity stubs
(getStressSensitivity = 0, commitSensitivity = -1), so a fork fiber in a
DDM-capable beam contributed d(sigma)/d(theta) = 0 without a warning. WP-262
implements them with dual numbers on the templated kernel
(LadrunoUniaxialJ2Kernel.h): one Dual<1> pass seeded on the active parameter and
on the committed history sensitivities.

The gate: a 2D dispBeamColumn cantilever with a fiber section of LadrunoUniaxialJ2
(Voce + linear isotropic + one Armstrong-Frederick term) under a CYCLIC lateral
load (load control, into plasticity both ways, so the history sensitivities
matter). The DDM tip-displacement sensitivity sensNodeDisp(tip, 1, param) must
equal a central finite difference of the whole run, at every step, for sig0, E
and Hiso. A Steel01 run of the same harness is the control: if it fails, the
harness is wrong, not the material.
"""
import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

L = 3000.0                 # mm
B, H = 300.0, 300.0        # mm
E0 = 200000.0              # MPa
S0, QINF, BISO, HISO = 250.0, 100.0, 20.0, 1500.0
CKIN, GKIN = 8000.0, 80.0
P_MAX = 640.0e3            # N: past the plastic moment Mp/L = 562 kN (hardening carries it);
                           # measured tip response 1.21-1.27x the elastic one, both ways
N_UP, N_DOWN = 20, 40      # steps up to +P_MAX, then down to -P_MAX


def _protocol():
    up = np.linspace(P_MAX / N_UP, P_MAX, N_UP)
    down = np.linspace(P_MAX - 2.0 * P_MAX / N_DOWN, -P_MAX, N_DOWN)
    return np.concatenate([up, down])


def _run(mat, vals, sens_param=None):
    """Run the cyclic protocol. vals: material constants. Returns tip
    displacement per step and, if sens_param, the DDM d(u_tip)/d(param)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, L)
    ops.fix(1, 1, 1, 1)
    if mat == 'LadrunoUniaxialJ2':
        ops.uniaxialMaterial('LadrunoUniaxialJ2', 1, vals['E'],
                             '-iso', 'voce', vals['sig0'], vals['Qinf'], vals['b'], vals['Hiso'],
                             '-kin', 1, CKIN, GKIN)
    else:
        ops.uniaxialMaterial('Steel01', 1, vals['sig0'], vals['E'], 0.01)
    ops.section('Fiber', 1)
    ops.patch('rect', 1, 24, 1, -H / 2, -B / 2, H / 2, B / 2)
    ops.geomTransf('Linear', 1)
    ops.beamIntegration('Lobatto', 1, 1, 5)
    ops.element('dispBeamColumn', 1, 1, 2, 1, 1)
    prot = _protocol()
    ops.timeSeries('Path', 1, '-time', *range(len(prot) + 1),
                   '-values', 0.0, *[float(v) for v in prot])
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    if sens_param is not None:
        key = {'sig0': 'sig0' if mat == 'LadrunoUniaxialJ2' else 'Fy',
               'E': 'E', 'Hiso': 'Hiso'}[sens_param]
        ops.parameter(1, 'element', 1, 'material', 1, key)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 100)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)       # pseudo-time; the Path series is the load
    ops.analysis('Static')
    if sens_param is not None:
        ops.sensitivityAlgorithm('-computeAtEachStep')
    us, ss = [], []
    for k in range(len(prot)):
        assert ops.analyze(1) == 0, f"{mat}: step {k} failed"
        us.append(ops.nodeDisp(2, 1))
        ss.append(ops.sensNodeDisp(2, 1, 1) if sens_param is not None else 0.0)
    ops.wipe()
    return np.array(us), np.array(ss)


def _check(mat, base, par, rel_h=1.0e-6, tol=1.0e-5):
    u0, ddm = _run(mat, base, par)
    h = rel_h * abs(base[par])
    vp = dict(base); vp[par] += h
    vm = dict(base); vm[par] -= h
    up, _ = _run(mat, vp)
    um, _ = _run(mat, vm)
    fd = (up - um) / (2.0 * h)
    scale = np.abs(fd).max()
    assert scale > 0.0, "the response does not depend on the parameter"
    err = np.abs(ddm - fd).max() / scale
    return err, ddm, fd, u0


BASE_LJ2 = dict(E=E0, sig0=S0, Qinf=QINF, b=BISO, Hiso=HISO)
BASE_S01 = dict(E=E0, sig0=S0)


def test_protocol_is_plastic_both_ways():
    u, _ = _run('LadrunoUniaxialJ2', BASE_LJ2)
    # elastic tip flexibility under the peak load vs the reached displacement
    I = B * H**3 / 12.0
    u_el = P_MAX * L**3 / (3.0 * E0 * I)
    assert u.max() > 1.15 * u_el and u.min() < -1.15 * u_el, "the protocol never yields"


@pytest.mark.parametrize("par", ["sig0", "E"])
def test_control_steel01_harness(par):
    """Steel01 has hand-coded DDM: the harness must reproduce its FD."""
    err, *_ = _check('Steel01', BASE_S01, par)
    assert err <= 1e-5, f"Steel01 control: DDM vs FD {err:.3e} -- the harness is wrong"


@pytest.mark.parametrize("par", ["sig0", "E", "Hiso"])
def test_ladruno_uniaxial_j2_ddm_matches_fd(par):
    err, ddm, fd, _ = _check('LadrunoUniaxialJ2', BASE_LJ2, par)
    assert np.abs(ddm).max() > 0.0, "DDM sensitivity is identically zero (the pre-WP-262 stub)"
    assert err <= 1e-5, f"LadrunoUniaxialJ2 d(u_tip)/d({par}): DDM vs FD {err:.3e}"
