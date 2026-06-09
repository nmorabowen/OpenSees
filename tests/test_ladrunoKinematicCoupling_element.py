"""LadrunoKinematicCoupling (RBE2 / kinematic coupling) — Zone-A battery.

A reference node R (6-DOF) that RIGIDLY drives a set of slave nodes: each slave follows
R's rigid-body motion u_i = u_R + theta_R x d_i (and theta_i = theta_R where the slave
rotation is tied), d_i = x_i - x_R. The rigid sibling of RBE3 — R is the MASTER (adds
stiffness to the slaves), not the dependent. See SRC/element/ladrunoKinematicCoupling/
and Ladruno_implementation/29_ladruno_kinematic_coupling_rbe2_adr.md.

Kinematic tests EXPLOIT a clean fact: with the slaves otherwise free, their only stiffness
is the penalty tie, so equilibrium drives each slave onto R's rigid prediction (gap -> 0,
independent of K). Force/moment tests fix the slaves, load R, and read the slave reactions
(global equilibrium: Sum reaction = -applied, Sum x_i x reaction = -applied moment).
"""
import math

import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------- helpers
def _face(ref_xyz=(0.0, 0.0, 0.0), a=1.0, slave_ndf=3, dof=None, kt=1.0e7, kr=None,
          enforce=None, bip=None, ref_mass=None, slave_mass=None):
    """Reference node (tag 1, ndf=6) rigidly driving 4 corner slaves (tags 2..5) at
    (+-a, +-a, 0). slave_ndf selects 3-DOF solid slaves or 6-DOF frame/shell slaves."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, *ref_xyz, "-ndf", 6)                       # reference node (master, 6-DOF)
    if ref_mass is not None:
        ops.mass(1, *ref_mass)
    corners = [(a, a, 0.0), (-a, a, 0.0), (-a, -a, 0.0), (a, -a, 0.0)]
    for i, c in enumerate(corners):
        ops.node(2 + i, *c, "-ndf", slave_ndf)             # slaves
        if slave_mass is not None:
            ops.mass(2 + i, *slave_mass)
    args = ["LadrunoKinematicCoupling", 1, 1, 4, 2, 3, 4, 5]
    if dof is not None:
        args += ["-dof", *dof]
    args += ["-k", kt]
    if kr is not None:
        args += ["-kr", kr]
    if enforce is not None:
        args += ["-enforce", enforce]
    if bip is not None:
        args += ["-bipenalty", "-dtcr", bip]
    ops.element(*args)
    return corners


def _solve_static():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _prescribe_ref(u6):
    """Prescribe all 6 DOFs of the reference master node via sp (incl. zeros)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for d in range(6):
        ops.sp(1, d + 1, u6[d])


# ----------------------------------------- 1. rigid translation -> slaves follow, zero gap
def test_rigid_translation_zero_gap():
    corners = _face()
    c = (0.003, -0.002, 0.001)
    _prescribe_ref((c[0], c[1], c[2], 0.0, 0.0, 0.0))
    _solve_static()
    for i in range(4):
        d = ops.nodeDisp(2 + i)
        for k in range(3):
            assert d[k] == pytest.approx(c[k], rel=1e-7, abs=1e-12)
    g = ops.eleResponse(1, "gap")
    assert max(abs(v) for v in g) < 1e-9                   # gap driven to ~0


# ----------------------------------------- 2. rigid rotation + transport (offset R) + SIGN
def test_rotation_transport_with_sign():
    """Offset reference, prescribed z-rotation phi: each slave must TRANSLATE by
    theta x d_i = phi*(-d_iy, d_ix, 0). The SIGN of the transport couple is asserted
    (a flipped transOp would reverse it). d_i = x_i - x_R."""
    xoff = 2.0
    corners = _face(ref_xyz=(xoff, 0.0, 0.0))
    phi = 1.0e-3
    _prescribe_ref((0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    for i, c in enumerate(corners):
        dx, dy = c[0] - xoff, c[1] - 0.0                   # d_i
        u = ops.nodeDisp(2 + i)
        assert u[0] == pytest.approx(-phi * dy, rel=1e-6, abs=1e-12)   # (theta x d)_x = -phi*d_y
        assert u[1] == pytest.approx(phi * dx, rel=1e-6, abs=1e-12)    # (theta x d)_y = +phi*d_x
        assert abs(u[2]) < 1e-9


# ----------------------------------------- 3a. -dof full-rigid drives slave rotation
def test_dof_full_rigid_drives_slave_rotation():
    """6-DOF slaves, default -dof (all 6): a prescribed reference rotation drives each
    slave's rotation theta_i = theta_R."""
    _face(slave_ndf=6)
    phi = 1.0e-3
    _prescribe_ref((0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    for i in range(4):
        assert ops.nodeDisp(2 + i)[5] == pytest.approx(phi, rel=1e-6)   # slave theta_z = phi
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(24.0)     # 6 DOF x 4 slaves


# ----------------------------------------- 3b. -dof translation-only leaves slave rotation free
def test_dof_translation_only_leaves_rotation_free():
    """6-DOF slaves, -dof 1 2 3: only translations tied. Slave rotations are NOT driven
    (fixed to 0 here to keep the system well-posed); a prescribed reference rotation still
    drives the slave TRANSLATIONS via transport, but the slave spin stays 0."""
    corners = _face(slave_ndf=6, dof=[1, 2, 3])
    for i in range(4):
        ops.fix(2 + i, 0, 0, 0, 1, 1, 1)                   # rotations fixed (untied -> would float)
    phi = 1.0e-3
    _prescribe_ref((0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    for i, c in enumerate(corners):
        u = ops.nodeDisp(2 + i)
        assert u[0] == pytest.approx(-phi * c[1], rel=1e-6, abs=1e-12)  # translation follows transport
        assert u[1] == pytest.approx(phi * c[0], rel=1e-6, abs=1e-12)
        assert abs(u[5]) < 1e-12                            # rotation NOT driven (stays fixed 0)
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(12.0)     # 3 DOF x 4 slaves


# ----------------------------------------- 4. force transfer (3-DOF slaves), global balance
def test_force_balance():
    """A force at R, slaves fixed: the slave reactions sum to -F (global equilibrium)."""
    _face(kt=1.0e7)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Fx = 1000.0
    ops.load(1, Fx, 0.0, 0.0, 0.0, 0.0, 0.0)
    _solve_static()
    ops.reactions()
    sx = sum(ops.nodeReaction(nd)[0] for nd in (2, 3, 4, 5))
    assert sx == pytest.approx(-Fx, rel=1e-5)


# ----------------------------------------- 5. moment transfer into a 3-DOF face (force couple)
def test_moment_into_solid():
    """A moment at R enters the 3-DOF slave face as a self-equilibrated force couple:
    zero net force, correct moment about R. The ndf-mismatch moment-transfer driver,
    rigid variant (the patch is held rigid — contrast RBE3)."""
    corners = _face(kt=1.0e7)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Mz = 1000.0
    ops.load(1, 0.0, 0.0, 0.0, 0.0, 0.0, Mz)
    _solve_static()
    ops.reactions()
    R = {nd: ops.nodeReaction(nd) for nd in (2, 3, 4, 5)}
    sx = sum(R[nd][0] for nd in (2, 3, 4, 5))
    sy = sum(R[nd][1] for nd in (2, 3, 4, 5))
    assert abs(sx) < 1e-3 and abs(sy) < 1e-3              # pure moment -> zero net force
    mz = sum(corners[i][0] * R[2 + i][1] - corners[i][1] * R[2 + i][0] for i in range(4))
    assert mz == pytest.approx(-Mz, rel=1e-5)            # moment balance about R (at origin)


# ----------------------------------------- 6. reduce to a single rigid link (N=1)
def test_single_rigid_link():
    """N=1 + a prescribed reference: the lone slave follows u_R + theta_R x d exactly —
    the rigidLink kinematics, generalized with a moment arm."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.5, 0.5, 0.0)
    ops.element("LadrunoKinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7)
    ux, phi = 0.002, 1.0e-3
    _prescribe_ref((ux, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    u = ops.nodeDisp(2)
    assert u[0] == pytest.approx(ux - phi * 0.5, rel=1e-6, abs=1e-12)   # u_R + theta x d, d=(1.5,0.5,0)
    assert u[1] == pytest.approx(phi * 1.5, rel=1e-6, abs=1e-12)


# ----------------------------------------- 7. augmented-Lagrangian solves (free slaves)
def test_al_runs():
    _face(kt=1.0e6, enforce="al")
    c = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3)
    _prescribe_ref(c)
    _solve_static()
    for i, cc in enumerate([(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)]):
        u = ops.nodeDisp(2 + i)
        assert u[0] == pytest.approx(-1.0e-3 * cc[1], rel=1e-5, abs=1e-9)   # gap ~ 0 at finite K


# ----------------------------------------- 8. bipenalty: massless reference self-report
def test_bipenalty_massless_ref():
    dt = 1.0e-3
    _face(kt=1.0e7, bip=dt)
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(dt, rel=1e-6)


# ----------------------------------------- 9. bipenalty: massless SLAVE is detected (E2)
def test_bipenalty_massless_slave_scanned():
    """The R-centric bipenalty of RBE3 would miss a massless dependent slave. RBE2 scans
    every tied DOF of R AND slaves: with R prescribed and one free massless slave, the
    slave's tied DOFs get the penalty mass and the self-report bounds them."""
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)               # -ndf 6 so R is massable (mass sizes by model ndf)
    ops.node(1, 0.0, 0.0, 0.0)                             # R (6-DOF)
    ops.mass(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)             # R massed -> only the slave is the hazard
    ops.node(2, 1.0, 0.0, 0.0, "-ndf", 3)                  # free, massless 3-DOF slave
    ops.element("LadrunoKinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7, "-bipenalty", "-dtcr", dt)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(dt, rel=1e-6)


# ----------------------------------------- 10. bipenalty: massed R is NOT double-counted (E1)
def test_bipenalty_massed_ref_not_doublecounted():
    """When every tied DOF already carries mass (R AND slaves massed), the element lumps
    NO penalty mass — R's own mass is not double-counted and the massed slaves are skipped
    — so the self-report has 'no opinion' (0). (A massless DOF WOULD be lumped, test 9.)"""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)             # -ndf 6 so every node is massable
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0)             # R fully massed
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, -1.0, 0.0, 0.0)
    ops.node(4, 0.0, 1.0, 0.0)
    for nd in (2, 3, 4):
        ops.mass(nd, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)       # slaves massed too -> nothing massless
    ops.element("LadrunoKinematicCoupling", 1, 1, 3, 2, 3, 4, "-k", 1.0e7, "-bipenalty", "-dtcr", 1.0e-3)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(0.0, abs=1e-12)   # nothing lumped


# ----------------------------------------- 11. TRANSIENT smoke (getDamp regression guard)
def test_transient_newmark_runs():
    """A pure penalty coupling overrides setRayleighDampingFactors to a no-op; the base
    Element::getDamp then dereferences an unallocated damping slot the first time a
    transient integrator forms the C-tangent -> hard crash. Guards the getDamp /
    getRayleighDampingForces overrides. Reference mass comes from bipenalty."""
    _face(kt=1.0e7, bip=1.0e-3)
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, 0.0, 0.0, 0.0, 1000.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-8, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(3, 1.0e-4) == 0                     # would CRASH without getDamp override
    assert math.isfinite(ops.nodeDisp(1)[5])


# ----------------------------------------- 12. degeneracy: self-tie refused (inert)
def test_self_tie_refused():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.0, 0.0, 0.0)
    # slave list includes the reference node -> refused inertly (valid=false, nGap stays 0)
    ops.element("LadrunoKinematicCoupling", 1, 1, 2, 1, 2, "-k", 1.0e7)
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(0.0)   # inert


# ----------------------------------------- 13. degeneracy: all-coincident, no Inf / floored Kr
def test_all_coincident_no_inf_and_floored_kr():
    """All slaves coincident with R. Translation-only + bipenalty must NOT report +Inf for
    the reference rotation; and with a rotation tied the default K_r must be floored > 0."""
    # (a) translation-only, all coincident, bipenalty -> finite dtcr (no 2*sqrt(I_p/0))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.node(3, 0.0, 0.0, 0.0)
    ops.element("LadrunoKinematicCoupling", 1, 1, 2, 2, 3, "-dof", 1, 2, 3,
                "-k", 1.0e7, "-bipenalty", "-dtcr", 1.0e-3)
    dt = ops.eleResponse(1, "dtcr")[0]
    assert math.isfinite(dt) and dt >= 0.0

    # (b) all coincident, rotation tied -> default K_r floored strictly positive
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.element("LadrunoKinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7)   # default -dof = all 6
    assert ops.eleResponse(1, "kr")[0] > 0.0


# ----------------------------------------- 14. derived rotational penalty K_r = K_t*ell^2
def test_derived_rotational_penalty():
    """Default K_r = K_t * ell^2, ell^2 floored to the largest |d_i|^2. Flat face a=1:
    |d_i|^2 = 2 each, mean = 2, so ell^2 = 2."""
    _face(kt=1.0e7, slave_ndf=6)
    assert ops.eleResponse(1, "kt")[0] == pytest.approx(1.0e7, rel=1e-9)
    assert ops.eleResponse(1, "kr")[0] == pytest.approx(1.0e7 * 2.0, rel=1e-6)


# ----------------------------------------- 15. serialization round-trip (sendSelf/recvSelf)
def test_database_roundtrip():
    """FE_Datastore round-trip exercises sendSelf/recvSelf + the broker. The probe reads
    element-owned, geometry-derived state (K_r and tiedDOFs), so it fails if recvSelf did
    not reconstruct the element (geometry/layout recomputed from coords on recv)."""
    def build():
        _face(kt=1.0e7, slave_ndf=6)
        _prescribe_ref((0.001, 0.0, 0.0, 0.0, 0.0, 0.0))
        _solve_static()

    database_roundtrip(
        build, probe_nodes=[1], ndf=6,
        probe_fn=lambda: [ops.eleResponse(1, "kr")[0], ops.eleResponse(1, "tiedDOFs")[0]],
    )
