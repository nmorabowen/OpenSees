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


# ----------------------------------------- 14. u-p slaves: ambiguous ndf refused without -dof
# TIMs 2026-09-07 no-ask finding 1. With the DEFAULT component list the element ties
# components 1..ndm+nrot on every slave and only checks that the node HAS that DOF
# index, so an ndf-4 u-p slave had its pressure DOF (index 3) tied to theta_x of the
# master, silently. The parser now refuses a slave whose ndf is neither ndm nor
# ndm+nrot when -dof is omitted; an explicit -dof keeps working.
def _up_slaves(with_dof):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)                   # master
    for i, c in enumerate([(1.0, 1.0, 0.0), (-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0), (1.0, -1.0, 0.0)]):
        ops.node(2 + i, *c, "-ndf", 4)                       # u-p slaves
    args = ["LadrunoKinematicCoupling", 1, 1, 4, 2, 3, 4, 5]
    if with_dof:
        args += ["-dof", 1, 2, 3]
    args += ["-k", 1.0e7]
    try:
        ops.element(*args)
        return True
    except Exception:
        return False


def test_up_slave_without_dof_is_refused(capfd):
    created = _up_slaves(with_dof=False)
    text = "".join(capfd.readouterr())
    assert "REFUSED" in text and "ndf = 4" in text and "-dof" in text, text
    # whether the binding raises or returns, the element must NOT exist
    assert (not created) or (1 not in ops.getEleTags())


def test_up_slave_with_explicit_dof_is_accepted():
    assert _up_slaves(with_dof=True)
    assert 1 in ops.getEleTags()
    # 4 slaves x 3 translational components tied, nothing on the pressure DOF
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(12.0)
    # and the tie is mechanically live: a master translation moves every slave
    _prescribe_ref([0.01, 0.0, 0.0, 0.0, 0.0, 0.0])
    for n in range(2, 6):
        ops.fix(n, 0, 0, 0, 1)                                # pin p; translations follow R
    _solve_static()
    for n in range(2, 6):
        assert ops.nodeDisp(n, 1) == pytest.approx(0.01, abs=1e-9)
        assert ops.nodeDisp(n, 4) == 0.0


def test_ndf3_and_ndf6_slaves_default_still_accepted():
    for ndf in (3, 6):
        _face(slave_ndf=ndf)                                  # default component list
        assert 1 in ops.getEleTags()


# ===================================================================== WP-101
# The RIGIDITY GATE: a real elastic block with a rigid footing skin driven by a
# master node.  This is the configuration the penalty formulation is actually
# asked to hold rigid (a footing on soil), and the one on which the pre-WP-101
# `-enforce al` was measured to be no better than plain penalty: its Uzawa
# recursion advanced ONCE PER COMMITTED STEP, so within a single push the tie
# was penalty-only and the constraint never converged.
#
# Metric: err = max |gap| over the footprint / |push|.  With the master's 6 DOFs
# prescribed the only thing resisting the tie is the soil, so the residual gap is
# exactly the tie force divided by K_t => err = c / K_t, c a property of the
# fixture (measured below: c ~ 1.66e4 for this block).
#
# Model: 2x2x2 stdBrick cube, B = 1.5 m, E = 45 000 kPa / nu = 0.3 (the TIMs
# strip's soil), base fixed, the 9 top-face nodes are the footing skin.
_E, _NU, _B, _N = 45000.0, 0.3, 1.5, 2
_PUSH = -0.005
_MASTER = 90001
_COUPLE = 90002


def _bn(i, j, k):
    return 1000 + (k * (_N + 1) + j) * (_N + 1) + i


def _soil_block():
    """Elastic block, base fixed. Returns (skin node tags, representative host ele tag)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    h = _B / _N
    for k in range(_N + 1):
        for j in range(_N + 1):
            for i in range(_N + 1):
                ops.node(_bn(i, j, k), i * h, j * h, k * h)
    ele = 1
    for k in range(_N):
        for j in range(_N):
            for i in range(_N):
                ops.element("stdBrick", ele,
                            _bn(i, j, k), _bn(i + 1, j, k), _bn(i + 1, j + 1, k), _bn(i, j + 1, k),
                            _bn(i, j, k + 1), _bn(i + 1, j, k + 1), _bn(i + 1, j + 1, k + 1),
                            _bn(i, j + 1, k + 1), 1)
                ele += 1
    for j in range(_N + 1):
        for i in range(_N + 1):
            ops.fix(_bn(i, j, 0), 1, 1, 1)
    skin = [_bn(i, j, _N) for j in range(_N + 1) for i in range(_N + 1)]
    return skin, 1


_TOL = {"NormUnbalance": 1.0e-6, "NormDispIncr": 1.0e-12}


def _build_gate(kt="auto", enforce=None, alupdate=None, host=False, host_first=True):
    """Elastic block + footing skin tied to a 6-DOF master. Returns the skin node tags."""
    skin, hostele = _soil_block()
    ops.node(_MASTER, _B / 2.0, _B / 2.0, _B, "-ndf", 6)
    args = ["LadrunoKinematicCoupling", _COUPLE, _MASTER, len(skin)] + skin
    args += ["-dof", 1, 2, 3, "-k", kt]
    if host or kt == "auto":
        args += ["-host", hostele]
    if enforce is not None:
        args += ["-enforce", enforce]
    if alupdate is not None:
        args += ["-alUpdate", alupdate]
    ops.element(*args)
    return skin


def _drive(integ="LoadControl", push=_PUSH):
    """Prescribe the master (LoadControl) or load it for DisplacementControl."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if integ == "LoadControl":
        for d, v in enumerate((0.0, 0.0, push, 0.0, 0.0, 0.0)):
            ops.sp(_MASTER, d + 1, v)
    else:
        for d in (1, 2, 4, 5, 6):
            ops.sp(_MASTER, d, 0.0)
        ops.load(_MASTER, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0)


def _gap_err(push=_PUSH):
    g = ops.eleResponse(_COUPLE, "gap")
    return max(abs(v) for v in g) / abs(push)


def _reaction_z():
    ops.reactions()
    return sum(ops.nodeReaction(_bn(i, j, 0))[2] for j in range(_N + 1) for i in range(_N + 1))


def _gate(kt="auto", enforce=None, alupdate=None, host=False, push=_PUSH):
    """Build the gate, push the MASTER by `push`, return (analyze status, err, sum Rz)."""
    _build_gate(kt, enforce, alupdate, host)
    _drive("LoadControl", push)
    ok = _solve_static_gate()
    return ok, _gap_err(push), _reaction_z()


def _solve_static_gate(tol=1e-14, maxiter=60, algo="Newton", integ="LoadControl",
                       testname="NormDispIncr", nsteps=1, push=_PUSH):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test(testname, tol, maxiter, 0)
    ops.algorithm(algo)
    if integ == "LoadControl":
        ops.integrator("LoadControl", 1.0)
    else:
        ops.integrator("DisplacementControl", _MASTER, 3, push / nsteps)
    ops.analysis("Static")
    return sum(1 for _ in range(nsteps) if ops.analyze(1) != 0)     # number of FAILED steps


def _augment_sweep(algo="Newton", integ="LoadControl", testname="NormDispIncr",
                   kt=1.0e6, n_aug=12, augtol=1.0e-9):
    """The ADR-41 D1 held-load within-step augmentation sweep = an OUTER Uzawa loop.

    Each held-load analyze() is an inner solve at FIXED lambda (so the residual is a
    genuine function of u and every algorithm works); the Domain::commit() that ends it
    performs the outer update lambda += D g via the element's commitState.
    """
    _build_gate(kt, "al", "commit")
    _drive(integ)
    if _solve_static_gate(_TOL[testname], 60, algo, integ, testname, 1) != 0:
        return -1, _gap_err(), 0
    ops.ladrunoBeginAugment()
    ops.integrator("LoadControl", 0.0)          # hold EVERYTHING, whatever drove the step
    npass = 0
    for _ in range(n_aug):
        if ops.analyze(1) != 0:
            ops.ladrunoEndAugment()
            return -2, _gap_err(), npass
        npass += 1
        if _gap_err() <= augtol:
            break
    ops.ladrunoEndAugment()
    return 0, _gap_err(), npass


def _direct_push_reaction():
    """Leg B: push the SAME footprint nodes directly (no coupling element)."""
    skin, _ = _soil_block()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in skin:
        ops.sp(n, 1, 0.0)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, _PUSH)
    assert _solve_static_gate() == 0
    ops.reactions()
    return sum(ops.nodeReaction(_bn(i, j, 0))[2] for j in range(_N + 1) for i in range(_N + 1))


# ------------------------- 20. penalty rigidity error is exactly c/K_t
def test_penalty_rigidity_error_scales_as_inverse_k():
    """err = c/K_t over three decades — the reason a "rigid" penalty tie needs a
    K_t nobody can condition. c is a fixture property, not a tolerance."""
    cs = []
    for kt in (1.0e6, 1.0e7, 1.0e8):
        ok, err, _ = _gate(kt=kt)
        assert ok == 0
        cs.append(err * kt)
    for c in cs[1:]:
        assert c == pytest.approx(cs[0], rel=0.05)             # err*K_t constant
    assert 1.0e4 < cs[-1] < 1.0e5                              # measured c ~ 1.66e4


# ------------------------- 21. -enforce al closes the constraint WITHIN one step (WP-101)
def test_al_closes_constraint_within_one_step():
    """THE WP-101 gate. At a MODERATE, host-order K_t = 1e6 (~1.3e2 x the host
    element's diagonal stiffness) the penalty alone leaves ~1.6e-2 of the push
    un-transmitted. With the per-iteration Uzawa update `-enforce al` drives the
    same gap below 1e-8 in ONE step, and the reaction matches a direct push of
    the same footprint exactly. `iter` is OPT-IN and only valid here because this
    fixture uses full Newton under LoadControl -- see test_iter_cadence_refused."""
    ok_p, err_p, _ = _gate(kt=1.0e6)
    assert ok_p == 0 and err_p > 1.0e-3                        # penalty: ~1.6e-2

    ok_a, err_a, rz_a = _gate(kt=1.0e6, enforce="al", alupdate="iter")
    assert ok_a == 0
    assert err_a <= 1.0e-8, f"AL left gap/push = {err_a:.3e} after one step"
    assert err_a < err_p / 1.0e5                               # >5 decades better

    rz_direct = _direct_push_reaction()
    assert rz_a == pytest.approx(rz_direct, rel=1e-9)


# ------------------------- 22. -alUpdate commit reproduces the legacy (pre-WP-101) behaviour
def test_al_update_commit_reproduces_penalty_within_a_step():
    """`-alUpdate commit` is the pre-WP-101 cadence: one Uzawa step per COMMITTED
    step, so within a single push the tie is penalty-only and the gap is exactly
    the penalty gap. This pins the escape hatch AND documents the defect."""
    _, err_pen, _ = _gate(kt=1.0e7)
    ok, err_commit, _ = _gate(kt=1.0e7, enforce="al", alupdate="commit")
    assert ok == 0
    assert err_commit == pytest.approx(err_pen, rel=1e-9)
    # ...and `commit` is the DEFAULT, so a bare -enforce al matches it
    ok_d, err_default, _ = _gate(kt=1.0e7, enforce="al")
    assert ok_d == 0
    assert err_default == pytest.approx(err_commit, rel=1e-12)


# ------------------------- 23. -k auto -host cannot hold a rigid footing
def test_k_auto_host_cannot_hold_a_rigid_footing():
    """`-k auto` scales K_t to the HOST's stiffness (kAlpha=1e3 x max|K_host(i,i)|)
    — that is its job: a tie that does not wreck conditioning. It is NOT a rigidity
    setting: on this gate it resolves to ~8e6 and leaves ~2e-3 of the push
    un-transmitted, 5 decades short of rigid. Use -enforce al (or a larger -k)
    when the footing must actually be rigid."""
    ok, err, _ = _gate(kt="auto")
    assert ok == 0
    kt = ops.eleResponse(_COUPLE, "kt")[0]
    assert 1.0e6 < kt < 1.0e8                                  # host-order, not 1e12
    assert err > 1.0e-4, f"-k auto unexpectedly rigid ({err:.3e})"


# ------------------------- 24. a failed step must not keep the trial-state multipliers
def test_al_lambda_reverts_on_failed_step():
    """With the per-iteration recursion lambda moves INSIDE the step, so
    revertToLastCommit (previously a bare `return 0`) must roll it back — otherwise
    a retried step inherits the multipliers of a discarded trial state."""
    assert _gate(kt=1.0e6, enforce="al", alupdate="iter")[0] == 0
    lam0 = list(ops.eleResponse(_COUPLE, "lambda"))
    assert max(abs(v) for v in lam0) > 0.0                     # AL is actually live

    # a second step that CANNOT converge (1 iteration, impossible tolerance)
    ops.sp(_MASTER, 3, 2.0 * _PUSH)
    ops.test("NormDispIncr", 1.0e-30, 1, 0)
    assert ops.analyze(1) != 0                                 # fails -> domain reverts
    lam1 = list(ops.eleResponse(_COUPLE, "lambda"))
    assert lam1 == pytest.approx(lam0, rel=1e-12, abs=1e-12)


# ------------------------- 25. conditioning guard: a huge K_t against a named -host warns
def test_high_kt_against_host_warns(capfd):
    """K_t far above the host's stiffness scale buys no rigidity (the error floors
    on round-off) and costs conditioning. With a -host named, say so."""
    ok, _, _ = _gate(kt=1.0e12, host=True)
    assert ok == 0
    text = "".join(capfd.readouterr())
    assert "x the -host element's stiffness scale" in text, text

    # ...and stay quiet inside the recommended band
    ok, _, _ = _gate(kt=1.0e6, host=True)
    assert ok == 0
    assert "stiffness scale" not in "".join(capfd.readouterr())


# ------------------------- 26. -alUpdate without -enforce al is called out
def test_al_update_without_al_warns(capfd):
    _face(kt=1.0e7)                                            # fresh model
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.element("LadrunoKinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7, "-alUpdate", "commit")
    text = "".join(capfd.readouterr())
    assert "-alUpdate has no effect without -enforce al" in text, text
    assert 1 in ops.getEleTags()                               # warned, not refused

# ------------------------- 27. algorithm x integrator x test sweep (WP-101 review r1)
@pytest.mark.parametrize("algo", ["Newton", "ModifiedNewton", "KrylovNewton"])
@pytest.mark.parametrize("integ", ["LoadControl", "DisplacementControl"])
@pytest.mark.parametrize("testname", ["NormUnbalance", "NormDispIncr"])
def test_enforcement_sweep(algo, integ, testname):
    """The 12-cell sweep that the first cut of WP-101 did not have, and which caught
    two hard regressions. Per cell:
      penalty  -- converges (0 failed steps), error at the c/K_t floor;
      commit   -- converges, and under LoadControl matches penalty exactly (one step,
                  one Uzawa update at the END of it, so the step itself is penalty);
      iter     -- EITHER converges to <=1e-8 (full Newton + LoadControl only) OR is
                  refused (the analysis fails loudly) -- never a silent wrong answer.
    """
    nsteps = 1 if integ == "LoadControl" else 5
    tol = _TOL[testname]

    _build_gate(1.0e6)
    _drive(integ)
    nfail_pen = _solve_static_gate(tol, 60, algo, integ, testname, nsteps)
    err_pen = _gap_err()
    assert nfail_pen == 0, f"penalty failed {nfail_pen}/{nsteps} steps"
    assert err_pen > 1.0e-3                                  # the c/K_t floor, ~1.6e-2

    _build_gate(1.0e6, "al", "commit")
    _drive(integ)
    nfail_com = _solve_static_gate(tol, 60, algo, integ, testname, nsteps)
    err_com = _gap_err()
    assert nfail_com == 0, f"-alUpdate commit failed {nfail_com}/{nsteps} steps"
    if integ == "LoadControl":
        assert err_com == pytest.approx(err_pen, rel=1e-9)   # one step => penalty exactly
    else:
        assert err_com <= err_pen                            # across-step Uzawa, 5 steps

    _build_gate(1.0e6, "al", "iter")
    _drive(integ)
    nfail_it = _solve_static_gate(tol, 60, algo, integ, testname, nsteps)
    supported = (algo == "Newton" and integ == "LoadControl")
    if supported:
        assert nfail_it == 0
        assert _gap_err() <= 1.0e-8
    else:
        assert nfail_it == nsteps, (
            f"-alUpdate iter must be REFUSED under {algo}/{integ}; it 'succeeded' "
            f"({nfail_it}/{nsteps} failed) with gap {_gap_err():.3e}"
        )


# ------------------------- 28. the SUPPORTED within-step route: the ADR-41 augment sweep
@pytest.mark.parametrize("algo", ["Newton", "ModifiedNewton", "KrylovNewton", "BFGS", "Broyden"])
@pytest.mark.parametrize("integ", ["LoadControl", "DisplacementControl"])
def test_augment_sweep_closes_constraint(algo, integ):
    """The consistent way to close the constraint within a step: an OUTER Uzawa loop
    with a FIXED residual per inner solve. The fork already ships the mechanism
    (ADR-41 D1: ladrunoBeginAugment / held-load analyze / ladrunoEndAugment), and
    LadrunoKinematicCoupling participates with NO new hook -- Domain::commit() still
    runs every element's commitState() during the sweep (contactAugmenting only
    suppresses the recorder loop and the commitTag bump), and with -alUpdate commit
    that commitState IS the outer update lambda += D g.

    Unlike -alUpdate iter this works with EVERY algorithm and integrator, because
    lambda is frozen inside each inner solve."""
    status, err, npass = _augment_sweep(algo=algo, integ=integ)
    assert status == 0, f"augment sweep failed (status {status}) after {npass} passes"
    assert err <= 1.0e-8, f"augment sweep left gap/push = {err:.3e} after {npass} passes"
    assert 1 <= npass <= 10


# ------------------------- 29. iter is refused, loudly, outside its validity window
def test_iter_cadence_refused(capfd):
    """-alUpdate iter makes the residual path-dependent, so it is refused anywhere but
    full Newton + LoadControl. The refusal must be LOUD (the analysis fails) rather than
    a silently wrong answer."""
    _build_gate(1.0e6, "al", "iter")
    _drive("DisplacementControl")
    nfail = _solve_static_gate(_TOL["NormDispIncr"], 60, "KrylovNewton",
                               "DisplacementControl", "NormDispIncr", 2)
    assert nfail == 2
    text = "".join(capfd.readouterr())
    assert "-alUpdate iter is REFUSED" in text, text
    assert "ladrunoBeginAugment" in text                     # points at the supported route


# ------------------------- 30. -alUpdate iter echoes its validity window at parse time
def test_iter_echo_at_parse(capfd):
    _build_gate(1.0e6, "al", "iter")
    text = "".join(capfd.readouterr())
    assert "EXPERT/opt-in" in text and "LoadControl" in text, text


# ------------------------- 31. declaring the coupling BEFORE its -host must still resolve
def test_host_declared_after_coupling_still_resolves(capfd):
    """Domain::addElement calls update() -> resolveAutoKt() at DECLARATION time. Latching
    ktResolved on that failed host lookup left `-k auto` silently at the 1e12 default and
    made the conditioning warning unreachable. The lookup must be retried, not latched."""
    # coupling FIRST, host element only afterwards
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    h = _B / _N
    for k in range(_N + 1):
        for j in range(_N + 1):
            for i in range(_N + 1):
                ops.node(_bn(i, j, k), i * h, j * h, k * h)
    ops.node(_MASTER, _B / 2.0, _B / 2.0, _B, "-ndf", 6)
    skin = [_bn(i, j, _N) for j in range(_N + 1) for i in range(_N + 1)]
    ops.element("LadrunoKinematicCoupling", _COUPLE, _MASTER, len(skin), *skin,
                "-dof", 1, 2, 3, "-k", "auto", "-host", 1)          # host ele 1 does not exist yet
    e = 1
    for k in range(_N):
        for j in range(_N):
            for i in range(_N):
                ops.element("stdBrick", e, _bn(i, j, k), _bn(i + 1, j, k), _bn(i + 1, j + 1, k),
                            _bn(i, j + 1, k), _bn(i, j, k + 1), _bn(i + 1, j, k + 1),
                            _bn(i + 1, j + 1, k + 1), _bn(i, j + 1, k + 1), 1)
                e += 1
    for j in range(_N + 1):
        for i in range(_N + 1):
            ops.fix(_bn(i, j, 0), 1, 1, 1)
    _drive("LoadControl")
    assert _solve_static_gate() == 0
    kt = ops.eleResponse(_COUPLE, "kt")[0]
    assert 1.0e6 < kt < 1.0e8, f"-k auto did not resolve off the late -host (K_t = {kt:.3e})"


# ------------------------- 32. round-trip must preserve the AL multiplier
def test_database_roundtrip_preserves_lambda():
    """The v1 round-trip probed kr + tiedDOFs only. lambda is now ALSO committed state
    (and lambdaCommitted rides the payload), and Domain::recv() calls update() right
    after recvSelf -- which under `iter` used to advance lambda by ~6.6e-9 on the
    just-restored state. Probe lambda so both are pinned."""
    def build():
        _build_gate(1.0e6, "al", "commit")
        _drive("LoadControl")
        assert _solve_static_gate() == 0

    database_roundtrip(
        build, probe_nodes=[_MASTER], ndf=6,
        probe_fn=lambda: list(ops.eleResponse(_COUPLE, "lambda")),
    )
