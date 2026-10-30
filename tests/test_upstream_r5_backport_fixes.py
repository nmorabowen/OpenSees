"""Fixes found while porting to jaabell/OpenSees (round 5, jaabell#72-#74).

1. LadrunoHHT re-forms the current tangent before the DDM solves. With
   ModifiedNewton on a yielding Steel01 truss the stale factorization gave
   d(u)/dfy about 1e-3 off a central finite difference.
2. An optional trailing number that is absent no longer eats the next
   token under openseespy (a failed OPS_GetDoubleInput advances there):
   LadrunoArcLength -stabilize / -mass lumped, LadrunoDynamicRelaxation
   -mass lumped / -damping viscous, LadrunoIndirectControl -iter.
3. InitDefGrad refuses a singular -F0 at the command instead of exiting
   the process from the constructor (run in a child interpreter).

LogStrain's revert reset is covered by the existing LogStrain protocol
tests only as a regression.
"""
import subprocess
import sys
import textwrap

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---- 1. LadrunoHHT DDM with ModifiedNewton ---------------------------------

E0, A, L, M = 2.0e3, 1.5, 1.0, 1.0
DT, NSTEPS = 5.0e-4, 80
ALPHA = 0.8
FY0, B0 = 20.0, 0.02


def _hht_run(fy, sens, algo):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.mass(2, M, M)
    ops.uniaxialMaterial('Steel01', 1, fy, E0, B0)
    ops.element('Truss', 1, 1, 2, A, 1)
    ops.timeSeries('Trig', 1, 0.0, DT*NSTEPS, DT*NSTEPS, '-factor', 2.0*FY0*A)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0)
    if sens:
        ops.parameter(1, 'element', 1, 'fy')
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 100)
    ops.algorithm(algo)
    ops.integrator('LadrunoHHT', ALPHA)
    ops.analysis('Transient')
    if sens:
        ops.sensitivityAlgorithm('-computeAtEachStep')
    umax = 0.0
    for _ in range(NSTEPS):
        assert ops.analyze(1, DT) == 0
        umax = max(umax, abs(ops.nodeDisp(2, 1)))
    u = ops.nodeDisp(2, 1)
    ddm = ops.sensNodeDisp(2, 1, 1) if sens else None
    ops.wipe()
    return u, umax, ddm


@pytest.mark.parametrize('algo', ['Newton', 'ModifiedNewton'])
def test_hht_ddm_fy_matches_finite_difference(algo):
    _, umax, ddm = _hht_run(FY0, True, algo)
    assert umax > FY0/E0*L                     # the truss yields
    dfy = 1.0e-6*FY0
    up, _, _ = _hht_run(FY0 + dfy, False, algo)
    um, _, _ = _hht_run(FY0 - dfy, False, algo)
    fd = (up - um)/(2.0*dfy)
    assert abs(fd) > 0.0
    assert abs(ddm - fd) <= 1.0e-6*abs(fd)


# ---- 2. optional trailing numbers ------------------------------------------

APEX = 2


def _truss():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, -1.0, 0.0)
    ops.node(APEX, 0.0, 0.1)
    ops.node(3, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(APEX, 1, 0)
    ops.uniaxialMaterial('Elastic', 1, 1.0e4)
    ops.element('corotTruss', 1, 1, APEX, 1.0, 1)
    ops.element('corotTruss', 2, APEX, 3, 1.0, 1)
    ops.mass(APEX, 1.0, 1.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(APEX, 0.0, -1.0)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-10, 50, 0)
    ops.algorithm('Newton')


@pytest.mark.parametrize('args', [
    ('LadrunoArcLength', 0.1, 1.0, '-stabilize', '-mass', 'unity'),
    ('LadrunoArcLength', 0.1, 1.0, '-stabilize', '-mass', 'lumped',
     '-cVisc', 0.5),
    ('LadrunoDynamicRelaxation', '-mass', 'lumped', '-damping', 'viscous'),
    ('LadrunoDynamicRelaxation', '-damping', 'viscous', '-mass', 'unity'),
], ids=['al-stabilize', 'al-lumped', 'dr-lumped', 'dr-viscous'])
def test_absent_optional_number_keeps_next_option(args, capfd):
    _truss()
    capfd.readouterr()
    ops.integrator(*args)
    err = capfd.readouterr().err
    assert 'unknown option' not in err, err


def test_indirect_control_iter_without_bounds_keeps_dof():
    _truss()
    ops.integrator('LadrunoIndirectControl', -0.004, '-iter', 4,
                   '-dof', APEX, 2, 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(APEX, 2) < 0.0


# ---- 3. InitDefGrad singular -F0 -------------------------------------------

def test_initdefgrad_singular_f0_is_refused():
    script = textwrap.dedent("""
        try:
            import opensees as ops
        except ModuleNotFoundError:
            import openseespy.opensees as ops
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
        ops.nDMaterial("LogStrain", 2, 1)
        try:
            ops.nDMaterial("InitDefGrad", 3, 2, "-F0", 1, 0, 0, 0, 1, 0, 0, 0, 0)
            print("ACCEPTED")
        except Exception:
            print("REFUSED")
        ops.nDMaterial("InitDefGrad", 4, 2)
        print("ALIVE")
    """)
    r = subprocess.run([sys.executable, "-c", script], capture_output=True,
                       text=True, timeout=120)
    assert "REFUSED" in r.stdout, (r.stdout, r.stderr)
    assert "ALIVE" in r.stdout, (r.stdout, r.stderr)
