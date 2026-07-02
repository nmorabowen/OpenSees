"""ADR-64 (LadrunoTie P4) — shell-to-solid plane-section tie (-shellSolid), Zone-A.

Ties an ndf-6 shell EDGE (MASTER polyline) to an ndf-3 solid FACE (SLAVES) with the
linearized rigid plane-section map (direction b-B, Abaqus shell-to-solid coupling):

    u_{s,i} = Σ_{j=A,B} N_j(ξ) [ u_{j,i} + (θ_j × d)_i ],   d = x_s − x̄(ξ)  FROZEN

Three mixed EQ rows per solid face node; drilling (θ·n) drops out via the 1e-12
near-zero filter (free, not clamped); tied DOFs fixed to the solid translations
{1,2,3} — BLOCKER-2 satisfied by -rho, no rotary-mass precondition. Cross-ndf rows
are Transformation-INCOMPATIBLE: static tests run Lagrange, explicit runs the
shipped LadrunoProjection. Oracle: proto_p4_shell_solid_tie.py (T1–T8, 22/22).

Mesh (shared): solid block x∈[0,1] (SSPbrick, 4×4×2, thickness t=0.1 in z,
y-grid NESTED 2:1 in the shell's) butting a shell strip x∈[1,2] (2 ShellMITC4,
y-split at 0.5, mid-surface z=0). The 15 solid face nodes at x=1
(y∈{0,.25,.5,.75,1} × z∈{−0.05,0,0.05}) are the slaves — the z=0 row exercises
the |d|=0 identity degeneracy; the 0.25/0.75 rows exercise genuine 0.5/0.5
interpolation on the 2-segment edge polyline. (Nested 2:1 keeps the collocation
force transfer CONSISTENT so the force-balance patch is exact — the same nesting
P1's patch test uses; a non-nested grid is a documented collocation limitation,
mortar's territory. Kinematic exactness holds for ANY grid — oracle T1–T5.)

  MEMBRANE  a uniform pull (ν=0, matching EA) transmits the constant-stress state
            exactly: interface at Δ/2, constant σ_xx in the solid.
  RIGID     all 6 master DOFs prescribed to a rigid rotation+translation: every
            solid node follows the rigid field through the tie (arm engaged).
  SECTION   a prescribed constant-moment section state (w0, θy0) lands the face
            nodes exactly on u = N(u_j + θ_j×d) — the plane-section kinematics,
            linear-in-thickness u_x, between DISSIMILAR element types.
  HEADLINE  an imposed far-edge rotation on matching-EI members produces the pure
            constant-moment state: the moment CROSSES the tie (interface rotation
            ≈ θ/2), and the face reproduces the plane-section rows exactly.
  CONTRAST  the same load through hand-emitted TRANSLATIONS-ONLY rows is a hinge:
            the shell rides freely (interface rotation ≈ θ), the solid stays
            unbent — the θ columns are the load-bearing delta of this rung.
  REFUSE    -mortar/-hermite/-dof combos, ndf<6 master edge, out-of-polyline
            projection, the -thickness arm guard, massless slave, non-disjoint,
            -masterEdge without -shellSolid.
  EXPLICIT  CentralDifferenceLadruno + LadrunoProjection + system Diagonal:
            dt_cr-neutral; the tie holds to machine precision every step under a
            transverse (moment-carrying) load. MEASURED CORRECTION to ADR-64
            Decision 1 win 4: the projector keeps every group DOF in the explicit
            equation set, so the tied masters' θx/θy DO need a small nodal rotary
            mass (generic to every rotational tie under handler 33001 — rigidLink
            -beam precedent; ShellMITC4 mass is CONSISTENT and refused outright,
            ASDShellQ4 lumps). The GENERATOR precondition is still gone: BLOCKER-2
            lives on the solid slaves (-rho), and static Lagrange needs no rotary
            mass at all. Total linear momentum is conserved through the tie
            (kick–coast, probe-calibrated nodal masses).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"

E, NU, T, RHO = 1000.0, 0.0, 0.1, 1.0
C = T / 2.0
MR = 1.0e-4                     # nodal rotary mass (only where a θ DOF is FREE+OUTSIDE the tie)

XS = [0.0, 0.25, 0.5, 0.75, 1.0]
YS = [0.0, 0.25, 0.5, 0.75, 1.0]   # nested 2:1 in the shell's 0.5 split (see docstring)
ZS = [-C, 0.0, C]

# shell nodes: (x, y) on z=0; edge nodes at x=1 are the MASTERS
_SHELL = {201: (1.0, 0.0), 202: (1.0, 0.5), 203: (1.0, 1.0),
          211: (2.0, 0.0), 212: (2.0, 0.5), 213: (2.0, 1.0)}
EDGE_MASTERS = [201, 202, 203]
FAR_EDGE = [211, 212, 213]
_SHELL_QUADS = [[201, 211, 212, 202], [202, 212, 213, 203]]
MASTER_EDGE = [2, 201, 202, 202, 203]      # -masterEdge <nseg=2> a1 b1 a2 b2


def _solid_tag(i, j, k):
    return 1000 + i * 100 + j * 10 + k      # i over XS, j over YS, k over ZS


def _mesh(tie=True, shell_rotary="none", shell_ele="ShellMITC4"):
    """Solid 4×4×2 SSPbrick + shell 2×ShellMITC4, tied with -shellSolid.
    shell_rotary: 'none' | 'far' (full rotary on far shell nodes + DRILLING-only
    on the masters — their θx/θy stay massless: inside the projection group) |
    'all' (full rotary everywhere — untied explicit reference).
    MIXED-NDF GOTCHA: the `mass` command sizes its matrix by the BUILDER's current
    ndf (OPS_GetNDF), not the node's, so the model is re-issued at ndf 6 before the
    shell side — the per-node '-ndf' override alone breaks `mass` on shell nodes."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for i, x in enumerate(XS):
        for j, y in enumerate(YS):
            for k, z in enumerate(ZS):
                ops.node(_solid_tag(i, j, k), x, y, z)
    eid = 1
    for i in range(len(XS) - 1):
        for j in range(len(YS) - 1):
            for k in range(len(ZS) - 1):
                n = [_solid_tag(i, j, k), _solid_tag(i + 1, j, k),
                     _solid_tag(i + 1, j + 1, k), _solid_tag(i, j + 1, k),
                     _solid_tag(i, j, k + 1), _solid_tag(i + 1, j, k + 1),
                     _solid_tag(i + 1, j + 1, k + 1), _solid_tag(i, j + 1, k + 1)]
                ops.element("SSPbrick", eid, *n, 1)
                eid += 1
    ops.model("basic", "-ndm", 3, "-ndf", 6)          # switch builder ndf (see gotcha)
    for t, (x, y) in _SHELL.items():
        ops.node(t, x, y, 0.0)
    ops.section("ElasticMembranePlateSection", 1, E, NU, T, RHO)
    for e, conn in enumerate(_SHELL_QUADS, start=100):
        ops.element(shell_ele, e, *conn, 1)
    if shell_rotary == "all":
        # EXPLICIT rule (measured, ADR-64 correction): shell getMass() has NO rotary
        # inertia, and the projection handler keeps every group DOF in the explicit
        # equation set — so the tied masters' θx/θy need a small nodal rotary mass
        # exactly like every rotational tie under handler 33001 (rigidLink -beam
        # precedent). The GENERATOR imposes no rotary-mass precondition (BLOCKER-2
        # lives on the solid slaves, auto-satisfied by -rho) and static Lagrange
        # needs none — the requirement is explicit-projector-only.
        for t in _SHELL:
            ops.mass(t, 0.0, 0.0, 0.0, MR, MR, MR)
    if tie:
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def _face_nodes():
    return [_solid_tag(len(XS) - 1, j, k) for j in range(len(YS)) for k in range(3)]


def _face_coords():
    return {_solid_tag(len(XS) - 1, j, k): (YS[j], ZS[k])
            for j in range(len(YS)) for k in range(3)}


def _plane_section_pred(y, z):
    """The b-B row RHS for a face slave at (1, y, z) from the SOLVED master DOFs:
    u_pred = Σ N_j [u_j + θ_j×d], d = (0,0,z), θ×d = (θy·z, −θx·z, 0)."""
    if y <= 0.5:
        (mA, mB), xi = (201, 202), y / 0.5
    else:
        (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
    pred = [0.0, 0.0, 0.0]
    for m, N in ((mA, 1.0 - xi), (mB, xi)):
        u = [ops.nodeDisp(m, d) for d in (1, 2, 3)]
        th = [ops.nodeDisp(m, d) for d in (4, 5, 6)]
        pred[0] += N * (u[0] + th[1] * z)
        pred[1] += N * (u[1] - th[0] * z)
        pred[2] += N * u[2]
    return pred


def _assert_face_satisfies_rows(tol):
    for s, (y, z) in _face_coords().items():
        pred = _plane_section_pred(y, z)
        for d in (1, 2, 3):
            err = abs(ops.nodeDisp(s, d) - pred[d - 1])
            assert err <= tol, f"face node {s} dof {d} violates the tie by {err}"


def _static_solve():
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# --------------------------------------------------------------------------
# MEMBRANE — constant stress crosses the tie (ν=0, matching EA)
# --------------------------------------------------------------------------
def test_membrane_patch_crosses_tie():
    DELTA = 1e-3
    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(len(YS)):
        for k in range(3):
            ops.fix(_solid_tag(0, j, k), 1, 1, 1)      # clamp the solid root (ν=0: clean)
    for t in FAR_EDGE:
        ops.sp(t, 1, DELTA)                            # pull the shell far edge in x
    _static_solve()
    # uniform strain over matching EA members: interface displaces Δ/2
    for s in _face_nodes():
        assert abs(ops.nodeDisp(s, 1) - 0.5 * DELTA) <= 1e-8 * DELTA, \
            f"face node {s} ux != Δ/2 (constant-stress state did not cross)"
    for m in EDGE_MASTERS:
        assert abs(ops.nodeDisp(m, 1) - 0.5 * DELTA) <= 1e-8 * DELTA
    # constant σ_xx in every solid element (one GP each — SSPbrick)
    sig = E * 0.5 * DELTA
    for e in range(1, 33):
        ops.eleResponse(e, "forces")
        s = list(ops.eleResponse(e, "stresses"))
        assert abs(s[0] - sig) <= 1e-6 * sig, f"element {e} σ_xx {s[0]} != {sig}"
    _assert_face_satisfies_rows(1e-12 * DELTA + 1e-15)


# --------------------------------------------------------------------------
# RIGID — a fully prescribed rigid master motion carries the whole solid
# --------------------------------------------------------------------------
def test_rigid_rotation_crosses_tie():
    U0 = (0.003, -0.002, 0.001)
    OM = (0.010, -0.020, 0.015)                        # incl. drilling-carrying θz
    P0 = (0.3, -0.2, 0.05)

    def cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for m, (x, y) in _SHELL.items():                   # prescribe EVERY shell DOF
        r = (x - P0[0], y - P0[1], 0.0 - P0[2])
        u = cross(OM, r)
        vals = [U0[0]+u[0], U0[1]+u[1], U0[2]+u[2], OM[0], OM[1], OM[2]]
        for d in range(6):
            ops.sp(m, d + 1, vals[d])
    _static_solve()
    err = 0.0
    for i, x in enumerate(XS):
        for j, y in enumerate(YS):
            for k, z in enumerate(ZS):
                r = (x - P0[0], y - P0[1], z - P0[2])
                u = cross(OM, r)
                for d in range(3):
                    err = max(err, abs(ops.nodeDisp(_solid_tag(i, j, k), d + 1)
                                       - (U0[d] + u[d])))
    assert err < 1e-10, f"solid did not follow the rigid field (err={err})"


# --------------------------------------------------------------------------
# SECTION — prescribed constant-moment section state: plane-section kinematics
# exact between dissimilar element types (incl. the |d|=0 z=0 identity row)
# --------------------------------------------------------------------------
def test_plane_section_state_exact_on_face():
    W0, TH0 = 2.0e-3, -4.0e-3                          # u=(θy·z, 0, w0) on the face
    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for m in _SHELL.keys():
        vals = [0.0, 0.0, W0, 0.0, TH0, 0.0]
        for d in range(6):
            ops.sp(m, d + 1, vals[d])
    _static_solve()
    for s, (y, z) in _face_coords().items():
        exp = (TH0 * z, 0.0, W0)
        for d in (1, 2, 3):
            err = abs(ops.nodeDisp(s, d) - exp[d - 1])
            assert err <= 1e-10, f"face node {s} dof {d} off the plane section by {err}"


# --------------------------------------------------------------------------
# HEADLINE — the moment CROSSES between dissimilar element types
# --------------------------------------------------------------------------
THETA = 0.01


def _bend_bcs_and_solve():
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(len(YS)):
        for k in range(3):
            ops.fix(_solid_tag(0, j, k), 1, 1, 1)      # solid root clamp
    for t in FAR_EDGE:
        ops.sp(t, 5, THETA)                            # imposed far-edge rotation θy
    _static_solve()


def test_bending_crosses_tie():
    """Uniform-EI members (same E, t) under an imposed end rotation with a free tip
    = a PURE constant-moment state: θ(x) ≈ THETA·x/2, so the shell edge must rotate
    ≈ THETA/2 — the moment crossed into the solid — and the face nodes must sit
    exactly on the plane-section rows (linear-in-thickness u_x)."""
    _mesh()
    _bend_bcs_and_solve()
    th_edge = sum(ops.nodeDisp(m, 5) for m in EDGE_MASTERS) / 3.0
    assert 0.3 * THETA < th_edge < 0.7 * THETA, \
        f"interface rotation {th_edge} not ~THETA/2 — the moment did not cross"
    _assert_face_satisfies_rows(1e-10)
    # the face carries the through-thickness couple: u_x odd in z, gradient = θy_edge
    fc = _face_coords()
    g = {}
    for s, (y, z) in fc.items():
        g.setdefault(y, {})[z] = ops.nodeDisp(s, 1)
    for y, col in g.items():
        grad = (col[C] - col[-C]) / (2.0 * C)
        assert abs(grad) > 0.2 * THETA, f"no through-thickness couple at y={y}"


def test_translations_only_contrast_is_hinge():
    """CONTRAST (oracle T8's FE twin): the same load through hand-emitted
    TRANSLATIONS-ONLY rows (u_s = Σ N_j u_j — the tie available before this rung)
    is a moment HINGE: the shell rides freely (θ_edge ≈ THETA, rigid rotation about
    the interface), and the solid face stays unbent."""
    _mesh(tie=False)
    for s, (y, z) in _face_coords().items():
        if y <= 0.5:
            (mA, mB), xi = (201, 202), y / 0.5
        else:
            (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
        for d in (1, 2, 3):
            args = [s, d, 1.0]
            for m, N in ((mA, 1.0 - xi), (mB, xi)):
                if abs(N) > 1e-12:
                    args += [m, d, -N]
            ops.equationConstraint(*args)
    _bend_bcs_and_solve()
    th_edge = sum(ops.nodeDisp(m, 5) for m in EDGE_MASTERS) / 3.0
    assert abs(th_edge - THETA) < 0.05 * THETA, \
        f"hinge should ride freely (θ_edge={th_edge} vs THETA={THETA})"
    fc = _face_coords()
    for y in YS:
        top = ops.nodeDisp([s for s, (yy, zz) in fc.items() if yy == y and zz == C][0], 1)
        bot = ops.nodeDisp([s for s, (yy, zz) in fc.items() if yy == y and zz == -C][0], 1)
        grad = (top - bot) / (2.0 * C)
        assert abs(grad) < 0.05 * THETA, \
            f"hinge leaked a through-thickness couple at y={y} (grad={grad})"


# --------------------------------------------------------------------------
# REFUSE — the named refusals
# --------------------------------------------------------------------------
def test_refuse_mortar_combo():
    _mesh(tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-mortar", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE)


def test_refuse_hermite_combo():
    _mesh(tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-hermite", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE)


def test_refuse_dof_combo():
    """Tied DOFs are FIXED to {1,2,3} by the plane-section operator — -dof refused."""
    _mesh(tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE, "-dof", 3, 1, 2, 3)


def test_refuse_low_ndf_master_edge():
    """The master edge must be ndf-6 shell nodes (the rows put moments there)."""
    _mesh(tie=False)
    a = _solid_tag(0, 0, 0)                            # ndf-3 solid nodes as "edge"
    b = _solid_tag(0, 2, 0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", 1, a, b)


def test_refuse_out_of_polyline_projection():
    """A slave beyond the edge polyline span (ξ out of [0,1] on every segment) is
    refused — the shell edge must span the solid face."""
    _mesh(tie=False)
    ops.node(999, 1.0, 1.5, 0.0)                       # beyond y=1 — off the polyline
    ops.mass(999, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999,
                       "-masterEdge", *MASTER_EDGE)


def test_refuse_arm_guard():
    """-thickness arms the |d| ≤ (0.5+tol)·t guard: a slave 2t off the mid-surface
    is not this shell's through-thickness material."""
    _mesh(tie=False)
    ops.node(999, 1.0, 0.5, 2.0 * T)
    ops.mass(999, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999,
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def test_arm_guard_skipped_without_thickness():
    """OQ-3 (signed off): -thickness omitted ⇒ warn and SKIP the arm guard — the
    same far slave is ACCEPTED (footgun-catcher, not a correctness precondition)."""
    _mesh(tie=False)
    ops.node(999, 1.0, 0.5, 2.0 * T)
    ops.mass(999, 1.0, 1.0, 1.0)
    ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999, "-masterEdge", *MASTER_EDGE)


def test_refuse_massless_slave():
    """BLOCKER-2 per-DOF: a tied solid DOF with no lumped mass is refused."""
    _mesh(tie=False)
    ops.node(999, 1.0, 0.5, 0.02)                      # on the face, NO mass, NO element
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999,
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def test_refuse_non_disjoint():
    """BLOCKER-1: a master edge node listed as a slave is refused (MP-chain)."""
    _mesh(tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 201,
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def test_refuse_masterEdge_without_shellSolid():
    """-masterEdge/-thickness are -shellSolid vocabulary."""
    _mesh(tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 15, *_face_nodes(), "-masterEdge", *MASTER_EDGE)


# --------------------------------------------------------------------------
# EXPLICIT — LadrunoProjection: dt_cr-neutral, machine-zero tie violation under a
# moment-carrying load with NO rotary mass on the tied masters, momentum conserved
# --------------------------------------------------------------------------
def _explicit_analysis():
    ops.constraints(PROJ)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")


def test_explicit_dt_cr_neutral():
    """Adding the plane-section tie does not change dt_cr (no penalty stiffness)."""
    dt0 = 1e-5
    _mesh(tie=True, shell_rotary="all", shell_ele="ASDShellQ4")
    _explicit_analysis()
    assert ops.analyze(1, dt0) == 0
    dt_tied = ops.criticalTimeStep()

    _mesh(tie=False, shell_rotary="all", shell_ele="ASDShellQ4")               # untied reference (free shell)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    assert ops.analyze(1, dt0) == 0
    dt_free = ops.criticalTimeStep()
    assert dt_tied > 0 and dt_free > 0
    assert abs(dt_tied / dt_free - 1.0) < 1e-9, f"dt_cr not neutral: {dt_tied} vs {dt_free}"


def test_explicit_zero_violation_moment_load():
    """A transverse far-edge load makes the moment cross the tie DYNAMICALLY: the
    plane-section rows hold to machine precision every step. (Masters carry the
    small nodal rotary mass every rotational tie needs under the explicit
    projector — see _mesh; the generator itself imposed no rotary precondition.)"""
    _mesh(tie=True, shell_rotary="all", shell_ele="ASDShellQ4")
    for j in range(len(YS)):
        for k in range(3):
            ops.fix(_solid_tag(0, j, k), 1, 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in FAR_EDGE:
        ops.load(t, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0)    # transverse tip load
    _explicit_analysis()
    assert ops.analyze(1, 1e-5) == 0
    dt = 0.4 * ops.criticalTimeStep()
    viol = 0.0
    for _ in range(300):
        assert ops.analyze(1, dt) == 0
        for s, (y, z) in _face_coords().items():
            pred = _plane_section_pred(y, z)
            for d in (1, 2, 3):
                viol = max(viol, abs(ops.nodeDisp(s, d) - pred[d - 1]))
    assert viol < 1e-10, f"tie violation {viol} (kinematic tie must be machine-zero)"
    assert abs(ops.nodeDisp(_solid_tag(4, 1, 2), 3)) > 1e-12, \
        "transverse load did not cross the tie into the solid"


def test_explicit_momentum_conserved():
    """Kick–coast momentum gate: probe-calibrate every node's lumped x-mass (untied,
    force at every node, one CDL step ⇒ v_i ∝ F/m_i with one universal factor that
    cancels), then kick the TIED free-free model on the solid root face and coast:
    Σ m̂_i·v_i,x is constant through the coast — the projection transfers momentum
    across the tie without leaking any."""
    all_tags = ([_solid_tag(i, j, k) for i in range(5) for j in range(len(YS)) for k in range(3)]
                + list(_SHELL.keys()))

    # probe: untied, full rotary mass everywhere so Diagonal is regular, Fx=1 on all
    _mesh(tie=False, shell_rotary="all", shell_ele="ASDShellQ4")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in all_tags:
        if t in _SHELL:
            ops.load(t, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        else:
            ops.load(t, 1.0, 0.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    dtp = 1e-5
    assert ops.analyze(1, dtp) == 0
    inv_m = {t: ops.nodeVel(t, 1) for t in all_tags}   # v_i = c/m_i, one universal c
    assert all(v > 0 for v in inv_m.values()), "probe velocities must be positive"

    # kick-coast on the TIED free-free model (no fixes at all)
    _mesh(tie=True, shell_rotary="all", shell_ele="ASDShellQ4")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    kick = [_solid_tag(0, j, k) for j in range(len(YS)) for k in range(3)]
    for t in kick:
        ops.load(t, 1.0, 0.0, 0.0)
    _explicit_analysis()
    assert ops.analyze(1, 1e-5) == 0
    dt = 0.4 * ops.criticalTimeStep()

    def momentum():
        return sum(ops.nodeVel(t, 1) / inv_m[t] for t in all_tags)

    for _ in range(30):
        assert ops.analyze(1, dt) == 0
    ops.remove("loadPattern", 1)                       # end of kick — coast
    # settle the leapfrog bookkeeping: the last half-step of forcing lands in the
    # first post-removal steps (verified identical tied vs UNTIED — not a tie effect)
    for _ in range(5):
        assert ops.analyze(1, dt) == 0
    p0 = momentum()
    assert p0 > 0.0, "kick produced no momentum"
    for _ in range(200):
        assert ops.analyze(1, dt) == 0
    p1 = momentum()
    assert abs(p1 - p0) <= 1e-5 * abs(p0), \
        f"momentum leaked through the tie: {p0} -> {p1}"
