"""Ground-motion inertia load counted ONCE on the vanilla elastic beams -- WP-119.

Vanilla ElasticBeam2d, ElasticTimoshenkoBeam2d and ElasticTimoshenkoBeam3d each
subtract their load vector (which holds only the UniformExcitation inertia
load, -M*R*a_g) in getResistingForce() when rho != 0, and then AGAIN in
getResistingForceIncInertia(), which calls getResistingForce() first. With
element mass, a ground-motion run drove the beam with 2*a_g. ElasticBeam3d was
correct. Still present in upstream OpenSees master (checked 2026-09-24).

Gates (ElasticBeam3d rides along as a control that was always right):
  t1  rigid-body probe -- a beam free only in y (x and rotations fixed) has
      zero stiffness in rigid y-translation, so under a constant ground
      acceleration a_g every node's RELATIVE acceleration must be exactly
      -a_g (the bug gave -2*a_g). Lumped and consistent (-cMass) mass.
  t1  deformable differential -- a cantilever with element -mass (lumped) under
      UniformExcitation must reproduce, step for step, the same cantilever
      carrying the identical lumped masses as nodal masses.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

AG, E, G, A, I, J, RHO = 2.0, 2.0e5, 8.0e4, 0.02, 1.5e-4, 3.0e-4, 0.8
AV = 0.8 * A

ELEMENTS = [("elasticBeamColumn", 2), ("elasticBeamColumn", 3),
            ("ElasticTimoshenkoBeam", 2), ("ElasticTimoshenkoBeam", 3)]
IDS = ["ElasticBeam2d", "ElasticBeam3d(control)", "ElasticTimoshenkoBeam2d", "ElasticTimoshenkoBeam3d"]


def _model(dim):
    ops.wipe()
    if dim == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf("Linear", 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)


def _section(name, dim):
    if name == "elasticBeamColumn":
        return (A, E, I) if dim == 2 else (A, E, G, J, I, I)
    # ElasticTimoshenkoBeam: 2D E G A Iz Avy ; 3D E G A J Iz Iy Avy Avz
    return (E, G, A, I, AV) if dim == 2 else (E, G, A, J, I, I, AV, AV)


def _element(name, dim, tag, i, j, element_mass, cmass=False):
    args = [name, tag, i, j, *_section(name, dim), 1]
    if element_mass:
        args += ["-mass", RHO]
    if cmass:
        args.append("-cMass")
    ops.element(*args)


def _transient():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.analysis("Transient")


# --------------------------------------------------------------------------
# gate 1: rigid-body probe
# --------------------------------------------------------------------------
@pytest.mark.parametrize("cmass", [False, True], ids=["lumped", "consistent"])
@pytest.mark.parametrize("name,dim", ELEMENTS, ids=IDS)
def test_rigid_body_relative_accel_is_minus_ag(name, dim, cmass):
    _model(dim)
    for k in range(3):
        ops.node(k + 1, *((1.5 * k, 0.0) if dim == 2 else (1.5 * k, 0.0, 0.0)))
        ops.fix(k + 1, *((1, 0, 1) if dim == 2 else (1, 0, 1, 1, 1, 1)))    # free in y only
    for e in (1, 2):
        _element(name, dim, e, e, e + 1, element_mass=True, cmass=cmass)
    ops.timeSeries("Constant", 1, "-factor", AG)
    ops.pattern("UniformExcitation", 1, 2, "-accel", 1)
    _transient()
    assert ops.analyze(1, 0.01) == 0
    acc = [ops.nodeAccel(k + 1, 2) for k in range(3)]
    worst = max(abs(a + AG) for a in acc)
    assert worst <= 1e-9 * AG, (
        f"[{name} {dim}D] relative y-acceleration {min(acc):+.6f}..{max(acc):+.6f}, expected {-AG} "
        "(-2*a_g means the ground-motion load is subtracted twice)")


# --------------------------------------------------------------------------
# gate 2: element -mass under UniformExcitation == the same masses as nodal masses
# --------------------------------------------------------------------------
NEL, L_BEAM = 4, 3.0


def _cantilever(name, dim, element_mass):
    _model(dim)
    h = L_BEAM / NEL
    for k in range(NEL + 1):
        ops.node(k + 1, *((k * h, 0.0) if dim == 2 else (k * h, 0.0, 0.0)))
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))
    for e in range(1, NEL + 1):
        _element(name, dim, e, e, e + 1, element_mass=element_mass)
    if not element_mass:
        for k in range(2, NEL + 2):
            m = RHO * h * (0.5 if k == NEL + 1 else 1.0)
            ops.mass(k, *([m, m, 0.0] if dim == 2 else [m, m, m, 0.0, 0.0, 0.0]))
    ops.timeSeries("Constant", 1, "-factor", AG)
    ops.pattern("UniformExcitation", 1, 2, "-accel", 1)
    _transient()
    w1 = math.sqrt(ops.eigen("-fullGenLapack", 1)[0])
    nstep = 100
    dt = 1.5 * (2.0 * math.pi / w1) / nstep
    out = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
        out.append(ops.nodeDisp(NEL + 1, 2))
    return out


@pytest.mark.parametrize("name,dim", ELEMENTS, ids=IDS)
def test_element_mass_equals_nodal_mass_under_ground_motion(name, dim):
    elem = _cantilever(name, dim, element_mass=True)
    nodal = _cantilever(name, dim, element_mass=False)
    scale = max(abs(u) for u in nodal)
    assert scale > 0.0
    worst = max(abs(a - b) for a, b in zip(elem, nodal))
    assert worst <= 1e-9 * scale, (
        f"[{name} {dim}D] element -mass departs from the same masses as nodal masses by "
        f"{worst / scale:.3e} of peak under UniformExcitation")
