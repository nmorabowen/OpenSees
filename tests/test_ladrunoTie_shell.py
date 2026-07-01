"""ADR-62 P3 — LadrunoTie shell / ROTATIONAL mesh-tie (ndf-6), Zone-A.

P1/P2 tie the translational DOFs of a non-conforming interface. P3 extends the SAME
per-slave transfer P (collocation shape weights, or the mortar D⁻¹M row) to the three
ROTATIONAL DOFs of an ndf-6 shell node:  θ_s = Σ_k P_sk θ_{m,k}. The generator now
defaults `-dof` to 1..6 for a 3D ndf-6 node (translations + rotations), emits the extra
rotational EQ_Constraint rows with the identical P coefficients, and the shipped
LadrunoProjectionHandler enforces them unchanged (every (node,dof) is an independent
union-find vertex — a rotational row is just more master-only vertices).

Scope v1 = shell-to-shell (same ndf, same P on rotations). Shell-to-solid (θ=½∇×u) is
a later rung. The EXACT bending case is CYLINDRICAL bending with the curvature axis
PARALLEL to the (straight) interface — then w and θ are constant along the interface and
the linear-complete tie reproduces the constant-moment state exactly (oracle T4). A
curvature varying ALONG the interface leaves an O(h²) residual on the quadratic w (never
on θ) — deferred as the Hermite w–θ transfer P3.1.

  BENDING   a constant-moment cantilever whose moment must cross a NON-CONFORMING shell
            tie: the tied split strip reproduces the monolithic conforming strip
            (interface rotation + transverse displacement, tip displacement) exactly.
  CONTRAST  tying translations ONLY (-dof 1 2 3) leaves the interface a hinge: the moment
            does NOT cross (the master stays flat) — proving the rotational tie carries it.
  MORTAR    a rigid-body rotation (w linear, θ constant — both linearly complete) crosses
            a genuinely non-matching mortar shell tie exactly, on ALL 6 DOFs.
  REFUSE    ShellMITC4 / ASDShellQ4 getMass() neglect rotational inertia, so tying a
            rotational DOF without explicit nodal rotary mass is refused (per-DOF
            BLOCKER-2); -dof 1 2 3 (translations only) is still accepted (escape hatch).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, H, RHO = 1000.0, 0.0, 0.1, 1.0
THETA = 0.01           # imposed tip rotation (rad) — constant-curvature cantilever
UZ, THY = 3, 5         # 1-based DOF indices: uz = 3, θy (rot about y) = 5
MR = 1.0e-4            # nodal rotary mass to satisfy the per-DOF BLOCKER-2 on tied rots


def _shell_section():
    ops.section("ElasticMembranePlateSection", 1, E, NU, H, RHO)


# ==========================================================================
# collocation butt-joint strip: master shell x∈[0,1], slave shells x∈[1,2]
# split in y at 0.5 (the y=0.5 slave interface node is NON-CONFORMING).
# ==========================================================================
_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0)}
MASTER_FACET = [1, 2, 3, 4]
_SLAVE = {11: (1, 0, 0), 12: (2, 0, 0), 13: (2, 0.5, 0), 14: (1, 0.5, 0),
          15: (2, 1, 0), 16: (1, 1, 0)}
_SLAVE_QUADS = [[11, 12, 13, 14], [14, 13, 15, 16]]
SLAVE_IFACE = [11, 14, 16]                     # x=1 slave nodes (14 is non-conforming)
_ALL = {**_MASTER, **_SLAVE}


def _mesh_butt(rotary_mass=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    _shell_section()
    ops.element("ShellMITC4", 1, *MASTER_FACET, 1)
    for e, conn in enumerate(_SLAVE_QUADS, start=2):
        ops.element("ShellMITC4", e, *conn, 1)
    if rotary_mass:
        for n in SLAVE_IFACE:                  # shells give NO rotary mass — add it
            ops.mass(n, 0.0, 0.0, 0.0, MR, MR, MR)


def _static_solve(handler="Lagrange"):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _apply_bending_bcs(nodes):
    """Cantilever constant-moment state: clamp the x=0 root fully, impose θy=THETA on the
    x=2 tip. Everything else (incl. the interface) free."""
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in nodes.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1, 1, 1, 1)               # root clamp
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, THY, THETA)                      # tip: θy prescribed (rest free)


def _monolithic_bending():
    """Conforming reference strip x∈[0,2] (2 quads, single y-span). Returns (θy, uz) at
    the x=1 interface and uz at the tip — the constant-curvature field the tie must match."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nd = {1: (0, 0, 0), 2: (0, 1, 0), 3: (1, 0, 0), 4: (1, 1, 0), 5: (2, 0, 0), 6: (2, 1, 0)}
    for t, (x, y, z) in nd.items():
        ops.node(t, float(x), float(y), float(z))
    _shell_section()
    ops.element("ShellMITC4", 1, 1, 3, 4, 2, 1)        # x∈[0,1]
    ops.element("ShellMITC4", 2, 3, 5, 6, 4, 1)        # x∈[1,2]
    _apply_bending_bcs(nd)
    _static_solve("Transformation")
    return ops.nodeDisp(3, THY), ops.nodeDisp(3, UZ), ops.nodeDisp(5, UZ)


# --------------------------------------------------------------------------
# BENDING — constant moment crosses the non-conforming shell tie exactly
# --------------------------------------------------------------------------
def test_collocation_bending_crosses_tie():
    th_ref, uz_ref, tip_ref = _monolithic_bending()
    assert abs(th_ref) > 1e-6, "reference did not bend — check the BCs"

    _mesh_butt(rotary_mass=True)
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET)
    _apply_bending_bcs(_ALL)
    _static_solve("Lagrange")

    tol = 1e-6 * max(abs(th_ref), abs(uz_ref), abs(tip_ref)) + 1e-9
    # the FREE master interface node (2) is pulled to the reference bending state by the
    # tie — the moment crossed from the slave (loaded) side into the master (clamped) side
    assert abs(ops.nodeDisp(2, THY) - th_ref) <= tol, "interface rotation did not cross"
    assert abs(ops.nodeDisp(2, UZ) - uz_ref) <= tol, "interface transverse disp did not cross"
    # the non-conforming slave node (14, y=0.5) gets the SAME rotation (constant along the
    # interface) via the 0.5/0.5 collocation weights
    assert abs(ops.nodeDisp(14, THY) - th_ref) <= tol, "non-conforming node rotation wrong"
    # the slave free end reproduces the reference tip displacement
    assert abs(ops.nodeDisp(12, UZ) - tip_ref) <= tol, "slave tip disp wrong"


def test_collocation_translations_only_blocks_moment():
    """CONTRAST: tie translations ONLY (-dof 1 2 3) ⇒ the interface is a moment hinge ⇒
    the tip rotation does NOT bend the (clamped) master. Proves the ROTATIONAL tie is what
    carries the moment across the interface."""
    th_ref, _, _ = _monolithic_bending()

    _mesh_butt(rotary_mass=True)
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET,
                   "-dof", 3, 1, 2, 3)                 # translations only
    _apply_bending_bcs(_ALL)
    _static_solve("Lagrange")

    # the master interface node barely rotates: the moment did not cross the hinge
    assert abs(ops.nodeDisp(2, THY)) < 1e-3 * abs(th_ref), \
        "moment crossed a translations-only tie — the interface should be a hinge"


# --------------------------------------------------------------------------
# MORTAR — a rigid rotation crosses a genuinely non-matching mortar shell tie
# --------------------------------------------------------------------------
def _shell_grid(n, tag0):
    """(n+1)² nodes / n² quads over the unit square [0,1]², z=0. Returns (nodes, quads)."""
    idx = lambda i, j: tag0 + j * (n + 1) + i
    nodes = {idx(i, j): (i / n, j / n, 0.0) for j in range(n + 1) for i in range(n + 1)}
    quads = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
             for j in range(n) for i in range(n)]
    return nodes, quads


def test_mortar_rigid_rotation_crosses_tie():
    """Overlapping coplanar shells, GENUINELY non-matching (master 2×2 vs slave 3×3),
    tied with -mortar. A rigid-body rotation about y (uz=-φx linear, θy=φ constant — both
    linearly complete) prescribed on the master perimeter must reproduce EXACTLY at every
    slave node on all 6 DOFs (the mortar rotational transfer)."""
    phi = 0.02
    m_nodes, m_quads = _shell_grid(2, 1)               # nodes 1..9, center = 5
    s_nodes, s_quads = _shell_grid(3, 101)             # nodes 101..116 (non-matching)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in {**m_nodes, **s_nodes}.items():
        ops.node(t, x, y, z)
    _shell_section()
    eid = 1
    for conn in m_quads + s_quads:
        ops.element("ShellMITC4", eid, *conn, 1)
        eid += 1
    for t in s_nodes:                                  # rotary mass on every tied slave node
        ops.mass(t, 0.0, 0.0, 0.0, MR, MR, MR)

    sf = [n for q in s_quads for n in q]
    mf = [n for q in m_quads for n in q]
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, len(s_quads), *sf,
                   "-masterFacets", 4, len(m_quads), *mf, "-outward", 0.0, 0.0, 1.0)

    # prescribe the rigid rotation on the master PERIMETER; leave the master center (5) free
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in m_nodes.items():
        if t == 5:
            continue                                   # free interior — solves to rigid mode
        ops.fix(t, 1, 1, 0, 1, 0, 1)                   # ux,uy,θx,θz = 0; uz(3),θy(5) via sp
        ops.sp(t, UZ, -phi * x)                        # uz = -φ x  (linear)
        ops.sp(t, THY, phi)                            # θy = φ     (constant)
    _static_solve("Lagrange")

    tol = 1e-6
    for t, (x, y, z) in s_nodes.items():
        assert abs(ops.nodeDisp(t, UZ) - (-phi * x)) <= tol, f"slave {t} uz wrong"
        assert abs(ops.nodeDisp(t, THY) - phi) <= tol, f"slave {t} θy wrong"
        for d in (1, 2, 4, 6):                         # ux, uy, θx, θz stay ~0
            assert abs(ops.nodeDisp(t, d)) <= tol, f"slave {t} spurious DOF {d}"


# --------------------------------------------------------------------------
# REFUSE — the per-DOF rotational-mass precondition
# --------------------------------------------------------------------------
def test_collocation_refuse_missing_rotary_mass():
    """A default shell tie (DOFs 1..6) with NO nodal rotary mass is refused: ShellMITC4
    getMass() neglects rotational inertia, so the tied rotational DOFs are massless."""
    _mesh_butt(rotary_mass=False)                      # translational mass only (from -rho)
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET)


def test_collocation_translations_only_accepted_without_rotary_mass():
    """Escape hatch: -dof 3 1 2 3 (translations only) is ACCEPTED without rotary mass —
    the per-DOF guard only requires mass on the DOFs actually tied (1,2,3 have -rho mass)."""
    _mesh_butt(rotary_mass=False)
    # accepted = returns without raising (mirrors test_mortar_nonaffine_slave_accepted)
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET,
                   "-dof", 3, 1, 2, 3)


def test_mortar_refuse_missing_rotary_mass():
    """Same per-DOF rotational-mass refusal on the -mortar path."""
    m_nodes, m_quads = _shell_grid(2, 1)
    s_nodes, s_quads = _shell_grid(3, 101)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in {**m_nodes, **s_nodes}.items():
        ops.node(t, x, y, z)
    _shell_section()
    eid = 1
    for conn in m_quads + s_quads:
        ops.element("ShellMITC4", eid, *conn, 1)
        eid += 1
    # NO rotary mass on the slave nodes
    sf = [n for q in s_quads for n in q]
    mf = [n for q in m_quads for n in q]
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, len(s_quads), *sf,
                       "-masterFacets", 4, len(m_quads), *mf, "-outward", 0.0, 0.0, 1.0)
