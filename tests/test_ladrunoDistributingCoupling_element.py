"""LadrunoDistributingCoupling (RBE3 / distributing coupling) — Zone-A battery.

A reference node R (6-DOF in 3D) whose motion is the weighted least-squares rigid-body
fit of a set of independent nodes (3-DOF) — equivalently, a load at R distributes to the
set without stiffening it. See SRC/element/ladrunoDistributingCoupling/ and
Ladruno_implementation/28_ladruno_distributing_coupling_rbe3_adr.md.

The kinematic tests EXPLOIT a clean fact: with the reference node otherwise free, its
only force is the penalty traction, so equilibrium forces the gap to ZERO exactly
(independent of K) and R lands exactly on the LSQ fit. The force tests fix the
independents and read their reactions (the distributed force set).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------- helpers
def _flat_face(weights=None, ref_xyz=(0.0, 0.0, 0.0), a=1.0, kt=1.0e7, kr=None,
               enforce=None, bip=None):
    """Reference node (tag 1, ndf=6) coupled to 4 independent corner nodes of a square
    face (tags 2..5) at (±a, ±a, 0). Equal weights -> centroid = origin."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, *ref_xyz, "-ndf", 6)                      # reference node (6-DOF)
    corners = [(a, a, 0.0), (-a, a, 0.0), (-a, -a, 0.0), (a, -a, 0.0)]
    for i, c in enumerate(corners):
        ops.node(2 + i, *c)                              # independents (3-DOF)
    args = ["LadrunoDistributingCoupling", 1, 1, 4, 2, 3, 4, 5]
    if weights is not None:
        args += ["-w", *weights]
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


def _prescribe_indep(field):
    """field: dict {nodeTag: (ux,uy,uz)} prescribed via sp under a Plain/Linear pattern.

    EVERY translational DOF of each independent node is prescribed (including zeros),
    so the only free DOFs are the reference node's — otherwise an unconstrained
    independent DOF (the RBE3 element gives it stiffness only through the gap) leaves a
    spurious mechanism and a singular system."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd, u in field.items():
        for d in range(3):
            ops.sp(nd, d + 1, u[d])


# ----------------------------------------- 1. weighted-average translation
def test_translation_weighted_average():
    """2 collinear independents, w=[0.25,0.75]; reference translation = weighted mean
    (rotations fixed so the dropped collinear axis is irrelevant)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, -1.0, 0.0, 0.0)
    ops.node(3, 1.0, 0.0, 0.0)
    ops.fix(1, 0, 0, 0, 1, 1, 1)                # free translations; fix ref rotations
    ops.element("LadrunoDistributingCoupling", 1, 1, 2, 2, 3, "-w", 0.25, 0.75, "-k", 1.0e7)
    _prescribe_indep({2: (0.004, 0.0, 0.0), 3: (0.0, 0.0, 0.0)})
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[0] == pytest.approx(0.25 * 0.004, rel=1e-8)   # 0.25*u2x + 0.75*u3x
    assert abs(d[1]) < 1e-12 and abs(d[2]) < 1e-12


# ----------------------------------------- 2. LSQ rotation fit (flat face)
def test_rotation_fit_flat_face():
    """Independents given a pure z-rotation field u_i = phi_z (z) x r_i; the reference
    (at the centroid, free) recovers theta_z = phi_z and zero translation."""
    corners = _flat_face()
    phi = 1.0e-3
    field = {2 + i: (-phi * c[1], phi * c[0], 0.0) for i, c in enumerate(corners)}
    _prescribe_indep(field)
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[5] == pytest.approx(phi, rel=1e-7)            # theta_z
    for k in (0, 1, 2, 3, 4):
        assert abs(d[k]) < 1e-9                            # no translation / x,y rotation


# ----------------------------------------- 3. offset-reference transport term
def test_offset_reference_transport():
    """Reference offset in x by xoff; under the same z-rotation field the reference must
    TRANSLATE by u_R = -theta x d = (0, phi*xoff, 0) (the transport term). A formulation
    that ties u_R directly to the centroid would give zero translation — the discriminator."""
    xoff = 2.0
    corners = _flat_face(ref_xyz=(xoff, 0.0, 0.0))
    phi = 1.0e-3
    field = {2 + i: (-phi * c[1], phi * c[0], 0.0) for i, c in enumerate(corners)}
    _prescribe_indep(field)
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[5] == pytest.approx(phi, rel=1e-7)            # theta_z = phi
    assert d[1] == pytest.approx(phi * xoff, rel=1e-6)     # u_R_y = phi*xoff (transport!)
    assert abs(d[0]) < 1e-9 and abs(d[2]) < 1e-9


# ----------------------------------------- 4. rigid translation -> zero gap
def test_rigid_translation_zero_gap():
    corners = _flat_face()
    c = (0.003, -0.002, 0.001)
    field = {2 + i: c for i in range(4)}
    _prescribe_indep(field)
    _solve_static()
    d = ops.nodeDisp(1)
    for k in range(3):
        assert d[k] == pytest.approx(c[k], rel=1e-7, abs=1e-12)
    for k in range(3, 6):
        assert abs(d[k]) < 1e-9
    g = ops.eleResponse(1, "gap")
    assert max(abs(v) for v in g) < 1e-9                   # gap driven to ~0


# ----------------------------------------- 5. force balance (distributed load)
def test_force_balance():
    """A force at the (centroid) reference distributes to the independents; their
    reactions sum to -F (global equilibrium of a statically-equivalent set)."""
    _flat_face(kt=1.0e7)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Fx = 1000.0
    ops.load(1, Fx, 0.0, 0.0, 0.0, 0.0, 0.0)
    _solve_static()
    ops.reactions()
    sx = sum(ops.nodeReaction(nd)[0] for nd in (2, 3, 4, 5))
    assert sx == pytest.approx(-Fx, rel=1e-5)              # Sum f_i = F
    # equal weights -> equal split
    for nd in (2, 3, 4, 5):
        assert ops.nodeReaction(nd)[0] == pytest.approx(-Fx / 4.0, rel=1e-4)


# ----------------------------------------- 6. moment -> distributed force couple
def test_moment_into_solid():
    """A moment at the reference enters the 3-DOF independent set as a self-equilibrated
    force couple: zero net force, correct moment about the centroid. This is the
    ndf-mismatch moment-transfer driver (6-DOF ref -> 3-DOF solid face)."""
    corners = _flat_face(kt=1.0e7)
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
    assert mz == pytest.approx(-Mz, rel=1e-5)            # moment balance about centroid


# ----------------------------------------- 7. flat face is full rank (debunked premise)
def test_flat_face_full_rank():
    _flat_face()
    assert ops.eleResponse(1, "nKept")[0] == pytest.approx(3.0)   # all 3 rotations resolved


# ----------------------------------------- 8. collinear set drops one rotation axis
def test_collinear_drops_one_axis():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, -1.0, 0.0, 0.0)
    ops.node(3, 0.0, 0.0, 0.0)
    ops.node(4, 1.0, 0.0, 0.0)
    ops.element("LadrunoDistributingCoupling", 1, 1, 3, 2, 3, 4, "-k", 1.0e7)
    nk = ops.eleResponse(1, "nKept")[0]
    assert nk == pytest.approx(2.0)                        # the line axis is dropped, no NaN


# ----------------------------------------- 9. derived rotational penalty K_r = K_t*ell^2
def test_derived_rotational_penalty():
    """Default K_r = K_t * (Sum w_i|r_i|^2 / W). Flat face a=1: |r_i|^2 = 2 each, ell^2 = 2."""
    _flat_face(kt=1.0e7)
    assert ops.eleResponse(1, "kt")[0] == pytest.approx(1.0e7, rel=1e-9)
    assert ops.eleResponse(1, "kr")[0] == pytest.approx(1.0e7 * 2.0, rel=1e-6)


# ----------------------------------------- 10. bipenalty dt_cr self-report
def test_bipenalty_dtcr():
    """-dtcr dt gives m_p = K_t(dt/2)^2 and I_p = K_r(dt/2)^2, so the self-reported
    critical step is exactly dt for both classes."""
    dt = 1.0e-3
    _flat_face(kt=1.0e7, bip=dt)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(dt, rel=1e-6)


# ----------------------------------------- 11. augmented-Lagrangian solves (free ref)
def test_al_runs():
    corners = _flat_face(kt=1.0e6, enforce="al")
    field = {2 + i: (0.002, 0.0, 0.0) for i in range(4)}
    _prescribe_indep(field)
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[0] == pytest.approx(0.002, rel=1e-6)          # free ref -> gap ~ 0 even at finite K
