"""`FluidSolidPorous` over the UW/Manzari family must BUILD, in 3D and in plane strain.

Gate for the one-line fix in `FluidSolidPorousMaterial`'s constructor (ADR 86
follow-up; TIMs `Work Plan/02` OQ29, `Objectives/11 -- The Solution Path` sec.9).

The defect: the wrapper cloned its soil with `soilMat.getCopy()`, i.e. the
dimension-free overload, on the BASE object the parser registered.  For the UW
family that overload is "subclass responsibility" and calls `exit(-1)`, so
`nDMaterial FluidSolidPorous $tag 3 <ManzariDafalias tag> ...` killed the process
at model-build time -- and `LadrunoSANISAND` inherited the same line.  The fix
asks for `getCopy("ThreeDimensional")` or `getCopy("PlaneStrain")` from the
wrapper's own `nd`, which every soil nDMaterial implements.

Why every case runs in a SUBPROCESS: the failure mode is a process exit, which
no in-process assertion can catch.  A subprocess that dies is the "fails before"
half of the gate; a returncode of 0 plus the identity below is the "passes after"
half.

The physics check is the wrapper's own plumbing identity, exact by construction
(`getStress` accumulates precisely this): after the undrained switch,

    p_exc == K_comb * (tr(eps) - tr(eps)|switch)

to 1e-6 relative, on a short imposed-strain leg.  Any departure means the strain
reaching the wrapper is not the strain reaching the soil -- which is exactly
what a wrong-dimension copy would produce -- or the cavitation clamp bound.

Wall time: three subprocesses, a few seconds each (zone_a, not slow).
"""
import os
import subprocess
import sys

import pytest

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

_HERE = os.path.dirname(os.path.abspath(__file__))

# Gorini's calibrated set, kPa (test_ladruno_sanisand.py carries the same list).
_PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
           0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
_TAIL = [1, 0, 1, 1.0e-10, 1.0e-10]          # ModifiedEuler, elastic tangent, analytic Jacobian
_E_INIT = 0.6944
_KW = 2.2e6
_POR = _E_INIT / (1.0 + _E_INIT)
_KCOMB = _KW / _POR                            # the wrapper takes K_f/n RAW (no porosity division inside)


_SCRIPT_3D = r'''
import sys, os
sys.path.insert(0, {here!r})
from _testbed import ops
P = {params!r}; T = {tail!r}
ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)
for n, (x, y, z) in enumerate([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)], start=1):
    ops.node(n, x, y, z)
for n in (1,4,5,8): ops.fix(n, 1,0,0)
for n in (1,2,5,6): ops.fix(n, 0,1,0)
for n in (1,2,3,4): ops.fix(n, 0,0,1)
if {ladruno!r}:
    ops.nDMaterial('LadrunoSANISAND', 1, *P, *T, '-Presidual', 0.0, '-Pmin', 0.0101, '-honorTolR', 0)
else:
    ops.nDMaterial('ManzariDafalias', 1, *P, *T)
ops.nDMaterial('FluidSolidPorous', 2, 3, 1, {kcomb!r}, 101.0)
ops.element('stdBrick', 1, 1,2,3,4,5,6,7,8, 2)
print('BUILT')
# isotropic consolidation to 100 kPa (static, load control)
ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
q4 = -100.0/4.0
for n in (2,3,6,7): ops.load(n, q4, 0, 0)
for n in (3,4,7,8): ops.load(n, 0, q4, 0)
for n in (5,6,7,8): ops.load(n, 0, 0, q4)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('UmfPack')
ops.test('NormUnbalance', 1e-6, 100, 0); ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1); ops.analysis('Static')
assert ops.analyze(10) == 0, 'consolidation failed'
ops.updateMaterialStage('-material', 1, '-stage', 1)   # the SOIL's plasticity
ops.updateMaterialStage('-material', 2, '-stage', 1)   # the WRAPPER's undrained switch
ops.integrator('LoadControl', 0.0)
assert ops.analyze(2) == 0, 'stage switch failed'
e0 = ops.eleResponse(1, 'material', 1, 'strain'); trE0 = e0[0] + e0[1] + e0[2]
ops.loadConst('-time', 0.0)
ops.timeSeries('Linear', 2); ops.pattern('Plain', 2, 2)
for n in (5,6,7,8): ops.sp(n, 3, -1.0)
ops.wipeAnalysis()
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormUnbalance', 1e-2, 200, 0); ops.algorithm('KrylovNewton')
ops.integrator('LoadControl', 0.0005); ops.analysis('Static')
worst = 0.0
for i in range(5):
    assert ops.analyze(1) == 0, 'push step %d failed' % (i+1)
    e = ops.eleResponse(1, 'material', 1, 'strain'); trE = e[0] + e[1] + e[2]
    pw = ops.eleResponse(1, 'material', 1, 'pressure')[0]
    pred = {kcomb!r} * (trE - trE0)
    rel = abs(pw - pred) / (abs(pred) + 1e-9)
    worst = max(worst, rel)
print('IDENTITY_WORST_REL', repr(worst))
print('PEXC_LAST', repr(pw))
'''

_SCRIPT_2D = r'''
import sys, os
sys.path.insert(0, {here!r})
from _testbed import ops
P = {params!r}; T = {tail!r}
ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
for n, (x, y) in enumerate([(0,0),(1,0),(1,1),(0,1)], start=1):
    ops.node(n, x, y)
ops.fix(1, 1, 1); ops.fix(2, 0, 1); ops.fix(4, 1, 0)
ops.nDMaterial('ManzariDafalias', 1, *P, *T)
ops.nDMaterial('FluidSolidPorous', 2, 2, 1, {kcomb!r}, 101.0)
ops.element('quad', 1, 1,2,3,4, 1.0, 'PlaneStrain', 2)
print('BUILT')
s = ops.eleResponse(1, 'material', 1, 'stress')
print('NSTRESS', len(s))
'''


def _run(script):
    p = subprocess.run([sys.executable, '-c', script], cwd=_HERE,
                       capture_output=True, text=True, timeout=300)
    return p.returncode, p.stdout + p.stderr


def _identity(out):
    for line in out.splitlines():
        if line.startswith('IDENTITY_WORST_REL'):
            return float(line.split()[1])
    raise AssertionError('identity line missing from subprocess output:\n' + out)


@pytest.mark.parametrize('ladruno', [False, True], ids=['ManzariDafalias', 'LadrunoSANISAND'])
def test_fspm_3d_builds_and_carries_the_identity(ladruno):
    """Before the fix this subprocess dies at `nDMaterial FluidSolidPorous` with
    "ManzariDafalias::getCopy -- subclass responsibility" (or the LadrunoSANISAND
    twin) and a non-zero returncode; after it, the model builds and the wrapper's
    plumbing identity holds to 1e-6 on five undrained steps."""
    rc, out = _run(_SCRIPT_3D.format(here=_HERE, params=_PARAMS, tail=_TAIL,
                                     ladruno=ladruno, kcomb=_KCOMB))
    assert 'subclass responsibility' not in out, out
    assert rc == 0, out
    assert 'BUILT' in out, out
    assert _identity(out) < 1.0e-6, out


def test_fspm_planestrain_builds_with_a_2d_copy():
    """The nd = 2 lane: the wrapper must ask the soil for a PlaneStrain copy, and
    the element must then see a 3-component stress from the wrapped material
    (a ThreeDimensional copy inside a plane-strain wrapper would report 6)."""
    rc, out = _run(_SCRIPT_2D.format(here=_HERE, params=_PARAMS, tail=_TAIL, kcomb=_KCOMB))
    assert 'subclass responsibility' not in out, out
    assert rc == 0, out
    assert 'BUILT' in out, out
    n = [int(l.split()[1]) for l in out.splitlines() if l.startswith('NSTRESS')]
    assert n and n[0] == 3, out
