"""ADR-52 W2-E1 — explicit bulk viscosity on LadrunoQuad / LadrunoCST (2D runtime).

The 2D companions to test_adr52_w2e1_bulk_viscosity.py (Brick). A single plane
element (PlaneStrain, rho>0) constrained for uniaxial motion -> nonzero volumetric
strain rate, kicked with an initial velocity. Integrated with EXPLICIT central
difference (bulk viscosity's intended regime; the lumped element mass + `system
Diagonal` give the diagonal mass it needs). Implicit Newmark would treat the
velocity-dependent viscous force explicitly within the implicit step and blow up at
large b1, so it is deliberately not used here.

Per element, two decisive checks:
  1. OFF-PATH IDENTITY: `-bulkViscosity 0 0` is bit-identical to no option.
  2. DISSIPATION + SIGN: b1>0 decays the oscillation and never amplifies it.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, RHO = 1000.0, 0.0, 1.0           # nu=0 -> D00=E, c_d=sqrt(E/rho)=31.6
DT, NSTEPS, V0 = 1.0e-3, 1600, -1.0     # dt << dt_cr (~0.03); ~11 axial periods
BV = (1.0, 1.2)                          # strong b1 for an unambiguous decay signal


def _analyse(rec_node, rec_dof):
    # explicit central difference -- bulk viscosity's intended regime; the lumped
    # element mass (getMass) + `system Diagonal` give the diagonal mass it needs.
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    hist = np.empty(NSTEPS)
    for k in range(NSTEPS):
        assert ops.analyze(1, DT) == 0, f"step {k} failed"
        hist[k] = ops.nodeDisp(rec_node, rec_dof)
    return hist


def _run_quad(bulk):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y))
    ops.fix(1, 1, 1); ops.fix(2, 1, 1)      # bottom fixed
    ops.fix(3, 1, 0); ops.fix(4, 1, 0)      # top: x fixed, y free (uniaxial)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    args = ["LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0]
    if bulk is not None:
        args += ["-bulkViscosity", float(bulk[0]), float(bulk[1])]
    ops.element(*args)
    ops.setNodeVel(3, 2, V0, "-commit")
    ops.setNodeVel(4, 2, V0, "-commit")
    return _analyse(3, 2)


def _run_cst(bulk):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in zip((1, 2, 3), [(0, 0), (1, 0), (0, 1)]):
        ops.node(t, float(x), float(y))
    ops.fix(1, 1, 1)                         # corner fixed
    ops.fix(2, 0, 1)                         # node 2: x free, y fixed
    ops.fix(3, 1, 0)                         # node 3: x fixed, y free
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    args = ["LadrunoCST", 1, 1, 2, 3, 1, "-type", "PlaneStrain"]
    if bulk is not None:
        args += ["-bulkViscosity", float(bulk[0]), float(bulk[1])]
    ops.element(*args)
    ops.setNodeVel(2, 1, V0, "-commit")
    ops.setNodeVel(3, 2, V0, "-commit")
    return _analyse(2, 1)


def _assert_dissipates(undamped, damped):
    q = len(undamped) // 4
    amp_early = np.max(np.abs(damped[:q]))
    amp_late = np.max(np.abs(damped[-q:]))
    amp_undamped_late = np.max(np.abs(undamped[-q:]))
    assert np.max(np.abs(undamped[-q:])) > 0.5 * np.max(np.abs(undamped[:q])), \
        "undamped baseline unexpectedly decayed -- check setup"
    assert amp_late < 0.9 * amp_early, \
        f"no decay: early={amp_early:.4e} late={amp_late:.4e}"
    assert amp_late < amp_undamped_late, \
        f"damped late {amp_late:.4e} not below undamped late {amp_undamped_late:.4e}"
    assert np.max(np.abs(damped)) <= 1.05 * np.max(np.abs(undamped)), \
        "bulk viscosity AMPLIFIED the response -- wrong sign"


def test_quad_bulk_viscosity_off_path_bit_identical():
    assert np.array_equal(_run_quad(None), _run_quad((0.0, 0.0))), \
        "LadrunoQuad: b1=b2=0 not bit-identical to no -bulkViscosity"


def test_quad_bulk_viscosity_dissipates():
    _assert_dissipates(_run_quad(None), _run_quad(BV))


def test_cst_bulk_viscosity_off_path_bit_identical():
    assert np.array_equal(_run_cst(None), _run_cst((0.0, 0.0))), \
        "LadrunoCST: b1=b2=0 not bit-identical to no -bulkViscosity"


def test_cst_bulk_viscosity_dissipates():
    _assert_dissipates(_run_cst(None), _run_cst(BV))
