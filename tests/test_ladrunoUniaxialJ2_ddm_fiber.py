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


# ---------------------------------------------------------------------------
# Truss-chain harness (WP-262 review findings 1 and 3): Lemaitre damage, the
# crack-band regularization -lchRef on elements of DIFFERENT lengths, IMPL-EX,
# and two parameters in one run (the SHV block per gradIndex). Finding 1 was
# that the DDM pass re-read ops_TheActiveElement -- stale outside
# Domain::update -- so every material replayed its map with the LAST element's
# lch: 1.7e-2 .. 1.2e-1 off FD on L = 1, 2 with -lchRef 1.
# ---------------------------------------------------------------------------
T_E, T_S0, T_HISO = 200000.0, 250.0, 2000.0


def _truss(vals, xs, lchref, sens=None, implex=False, pmax=330.0, nup=20):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    for i, x in enumerate(xs):
        ops.node(i + 1, x)
    ops.fix(1, 1)
    args = ['LadrunoUniaxialJ2', 1, vals['E'], '-iso', 'voce', vals['sig0'], 0.0, 0.0, vals['Hiso'],
            '-damage', 'lemaitre', 0.5, 1.0, 0.0, 0.99]
    if implex:
        args += ['-implex']
    if lchref > 0.0:
        args += ['-lchRef', lchref]
    ops.uniaxialMaterial(*args)
    ne = len(xs) - 1
    for e in range(ne):
        ops.element('Truss', e + 1, e + 1, e + 2, 1.0, 1)
    prot = ([pmax * (k + 1) / nup for k in range(nup)]
            + [pmax * (1.0 - 0.1 * (k + 1) / 40.0) for k in range(40)])
    ops.timeSeries('Path', 1, '-time', *range(len(prot) + 1), '-values', 0.0, *prot)
    ops.pattern('Plain', 1, 1)
    ops.load(len(xs), 1.0)
    names = list(sens) if sens else []
    for j, nm in enumerate(names):
        ops.parameter(j + 1, 'element', 1, nm)
        for e in range(1, ne):
            ops.addToParameter(j + 1, 'element', e + 1, nm)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 100)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    if names:
        ops.sensitivityAlgorithm('-computeAtEachStep')
    us, ss = [], []
    for k in range(len(prot)):
        assert ops.analyze(1) == 0, f"truss: step {k} failed"
        us.append(ops.nodeDisp(len(xs), 1))
        ss.append([ops.sensNodeDisp(len(xs), 1, j + 1) for j in range(len(names))])
    ops.wipe()
    return np.array(us), np.array(ss)


def _truss_check(xs, lchref, names, implex=False):
    base = dict(E=T_E, sig0=T_S0, Hiso=T_HISO)
    _, ddm = _truss(base, xs, lchref, sens=names, implex=implex)
    errs = {}
    for j, nm in enumerate(names):
        h = 1.0e-6 * abs(base[nm])
        vp = dict(base); vp[nm] += h
        vm = dict(base); vm[nm] -= h
        up, _ = _truss(vp, xs, lchref, implex=implex)
        um, _ = _truss(vm, xs, lchref, implex=implex)
        fd = (up - um) / (2.0 * h)
        scale = np.abs(fd).max()
        assert scale > 0.0
        errs[nm] = np.abs(ddm[:, j] - fd).max() / scale
    return errs


@pytest.mark.parametrize("xs", [[0.0, 1.0, 3.0], [0.0, 2.0, 3.0]], ids=["L1-2", "L2-1"])
def test_ddm_crack_band_unequal_lengths(xs):
    """Regression for review finding 1: -lchRef on elements of different length."""
    errs = _truss_check(xs, 1.0, ["sig0"])
    assert errs["sig0"] <= 1e-6, f"lchRef DDM vs FD {errs['sig0']:.3e} (stale active element?)"


def test_ddm_two_parameters_one_run_damage():
    """Two gradients in one run: each SHV block is indexed by its gradIndex."""
    errs = _truss_check([0.0, 1.0, 3.0], 1.0, ["E", "Hiso"])
    for nm, e in errs.items():
        assert e <= 1e-6, f"two-parameter run, {nm}: DDM vs FD {e:.3e}"


def test_ddm_damage_implex():
    errs = _truss_check([0.0, 1.0, 3.0], 1.0, ["sig0"], implex=True)
    assert errs["sig0"] <= 1e-6, f"IMPL-EX DDM vs FD {errs['sig0']:.3e}"
