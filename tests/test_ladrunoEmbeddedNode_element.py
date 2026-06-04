"""LadrunoEmbeddedNode (general node-to-host coupling element) — Zone-A battery.

Phase 1 = the ISOTROPIC translational tie (U). The element ties a constrained node's
translations to a host element's nodes via shape-function weights N_i (penalty/AL,
`-k auto`, bipenalty), sharing LadrunoEmbeddedKernel with LadrunoEmbeddedRebar. See
SRC/element/ladrunoEmbeddedNode/ and Ladruno_implementation/23_ladruno_embedded_node_adr.md.

Like the rebar battery, these checks pin the coupling MECHANICS without a host continuum:
we play the "host" with explicit nodes whose displacements we prescribe (fix), so the gap
g = u_c - sum N_i u_host is fully controlled. The -host/-xi legs use a real LadrunoBrick.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------- helpers
def _single_host_node(ku=1.0e5, ndm=3, enforce=None):
    """cNode (1) tied to a single coincident host node (2, fixed) with N=[1.0], so the
    gap is g = u_1 — an isotropic 3D/2D spring of stiffness ku in every direction."""
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    ops.node(1, *([0.0] * ndm))       # constrained node
    ops.node(2, *([0.0] * ndm))       # host node (coincident)
    ops.fix(2, *([1] * ndm))          # prescribe host -> u_host = 0
    args = ["LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku]
    if enforce is not None:
        args += ["-enforce", enforce]
    ops.element(*args)


def _push1(dof, u, nstep=1):
    """Prescribe ONE dof of the constrained node and return its reactions. The other
    cNode dofs are left FREE (they spring to ~0 through the isotropic penalty), so the
    system always keeps >=1 free DOF (don't fix-host AND sp-all-cNode -> 0 free DOFs)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, dof, u)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 20); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nstep); ops.analysis("Static")
    assert ops.analyze(nstep) == 0
    ops.reactions()
    return ops.nodeReaction(1)


# --------------------------------------------------- 1. isotropic spring (3D + 2D)
def test_isotropic_x_3d():
    ku = 1.0e5
    _single_host_node(ku=ku)
    R = _push1(1, 2.0e-3)                      # push x; y,z free -> 0
    assert abs(R[0]) == pytest.approx(ku * 2.0e-3, rel=1e-6)
    assert abs(R[1]) < 1e-3 and abs(R[2]) < 1e-3   # no cross-coupling (isotropic)


def test_isotropic_y_3d():
    """Same stiffness in y as in x (isotropic — unlike the anisotropic rebar)."""
    ku = 1.0e5
    _single_host_node(ku=ku)
    R = _push1(2, 1.5e-3)                      # push y; x,z free -> 0
    assert abs(R[1]) == pytest.approx(ku * 1.5e-3, rel=1e-6)
    assert abs(R[0]) < 1e-3 and abs(R[2]) < 1e-3


def test_isotropic_2d():
    ku = 1.0e5
    _single_host_node(ku=ku, ndm=2)
    R = _push1(1, 2.0e-3)                      # push x; y free -> 0
    assert abs(R[0]) == pytest.approx(ku * 2.0e-3, rel=1e-6)
    assert abs(R[1]) < 1e-3


# --------------------------------------------------- 2. partition of unity on hosts
def test_force_split_by_N():
    """Two host nodes with N=[0.25,0.75]; the reaction the element pushes onto the
    hosts splits as N_i*(cNode force) and balances the cNode (P = B^T t)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.node(3, 0.0, 0.0, 0.0)
    for nd in (2, 3):
        ops.fix(nd, 1, 1, 1)
    ops.element("LadrunoEmbeddedNode", 1, 1, 2, 2, 3, "-shape", 0.25, 0.75, "-k", 1.0e5)
    R = _push1(1, 1.0e-3)
    ops.reactions()
    fr = R[0]
    f2 = ops.nodeReaction(2)[0]
    f3 = ops.nodeReaction(3)[0]
    assert (fr + f2 + f3) == pytest.approx(0.0, abs=1e-6 * abs(fr))   # equilibrium
    assert f2 == pytest.approx(0.25 * (-fr), rel=1e-6)
    assert f3 == pytest.approx(0.75 * (-fr), rel=1e-6)


# --------------------------------------------- 3. -host / -xi with a real LadrunoBrick
_CUBE = {  # unit-cube nodes in LadrunoBrick natCoord order ((c+1)/2 mapping)
    11: (0.0, 0.0, 0.0), 12: (1.0, 0.0, 0.0), 13: (1.0, 1.0, 0.0), 14: (0.0, 1.0, 0.0),
    15: (0.0, 0.0, 1.0), 16: (1.0, 0.0, 1.0), 17: (1.0, 1.0, 1.0), 18: (0.0, 1.0, 1.0),
}
_CUBE_NAT = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
             (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def _trilinear(xi):
    return [0.125 * (1 + cx * xi[0]) * (1 + cy * xi[1]) * (1 + cz * xi[2])
            for (cx, cy, cz) in _CUBE_NAT]


def _cube_host_node(trailing):
    """Single LadrunoBrick host (tag 100, nodes 11..18 all fixed) + a constrained node
    (1) embedded via `trailing` (e.g. ['-host', 100, '-xi', 0, 0, 0, '-k', 1e5])."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
        ops.fix(tag, 1, 1, 1)
    ops.node(1, 0.5, 0.5, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element("LadrunoEmbeddedNode", 1, 1, *trailing)


def _assert_split(fr, weights):
    ops.reactions()
    for tag, Ni in zip(_CUBE, weights):
        assert ops.nodeReaction(tag)[0] == pytest.approx(Ni * (-fr), rel=1e-6)


def test_host_form_resolves_nodes():
    """`-host eleTag` pulls the 8 host nodes off the LadrunoBrick automatically."""
    ku = 1.0e5
    _cube_host_node(["-host", 100, "-shape", *([0.125] * 8), "-k", ku])
    R = _push1(1, 2.0e-3)
    assert abs(R[0]) == pytest.approx(ku * 2.0e-3, rel=1e-6)
    _assert_split(R[0], [0.125] * 8)


def test_xi_centroid_weights():
    """-xi 0 0 0 queries the host -> all weights 0.125 (centroid of the cube)."""
    ku = 1.0e5
    _cube_host_node(["-host", 100, "-xi", 0.0, 0.0, 0.0, "-k", ku])
    R = _push1(1, 2.0e-3)
    _assert_split(R[0], [0.125] * 8)


def test_xi_offcentroid_matches_trilinear():
    xi = (0.3, -0.2, 0.5)
    _cube_host_node(["-host", 100, "-xi", *xi, "-k", 1.0e5])
    R = _push1(1, 2.0e-3)
    _assert_split(R[0], _trilinear(xi))


# --------------------------------------------------- 4. -k auto (host-stiffness scaled)
def test_k_auto_scales_linearly_with_alpha():
    """K_u = kAlpha * max|K_host(i,i)| (ADR 23 D3/M2): doubling kAlpha doubles K_u."""
    _cube_host_node(["-host", 100, "-xi", 0.0, 0.0, 0.0, "-k", "auto", "-kAlpha", 1.0e3])
    _push1(1, 1.0e-3)
    ku1 = ops.eleResponse(1, "k")[0]
    _cube_host_node(["-host", 100, "-xi", 0.0, 0.0, 0.0, "-k", "auto", "-kAlpha", 2.0e3])
    _push1(1, 1.0e-3)
    ku2 = ops.eleResponse(1, "k")[0]
    assert ku1 > 0.0
    assert ku2 == pytest.approx(2.0 * ku1, rel=1e-9)


# --------------------------------------------------- 5. augmented Lagrangian
def test_enforce_al_drives_gap_to_zero():
    """cNode FREE in x, loaded by F, held by the penalty tie to a fixed host. Penalty
    leaves a residual gap F/K_u; AL's per-step Uzawa drives it ~0 at MODERATE K_u."""
    F, Ku = 50.0, 1.0e4

    def _build(enf):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 0.0, 0.0, 0.0)
        ops.fix(2, 1, 1, 1)
        ops.fix(1, 0, 1, 1)            # cNode free in x only
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", Ku,
                    "-enforce", enf)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, F, 0.0, 0.0)
        ops.constraints("Transformation"); ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1e-12, 30); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0); ops.analysis("Static")

    # penalty: ramp to F in one step -> residual gap F/K_u
    _build("penalty")
    assert ops.analyze(1) == 0
    gap_pen = abs(ops.nodeDisp(1, 1))
    assert gap_pen == pytest.approx(F / Ku, rel=1e-3)

    # AL: ramp to F (step 1), freeze the load, then re-solve at CONSTANT F so each
    # commit fires a Uzawa update that drives the gap toward zero.
    _build("al")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.integrator("LoadControl", 0.0)          # no further load increment
    assert ops.analyze(8) == 0                   # 8 constant-load Uzawa updates
    gap_al = abs(ops.nodeDisp(1, 1))
    assert gap_al < 0.05 * gap_pen               # AL near-exact at moderate K_u


# --------------------------------------------------- 6. responses
def test_responses():
    ku = 1.0e5
    _single_host_node(ku=ku)
    _push1(1, 1.0e-3)                                    # u_c = (1e-3, 0, 0)
    assert ops.eleResponse(1, "gap")[0] == pytest.approx(1.0e-3, rel=1e-6)
    assert ops.eleResponse(1, "force")[0] == pytest.approx(ku * 1.0e-3, rel=1e-6)
    assert ops.eleResponse(1, "k")[0] == pytest.approx(ku, rel=1e-12)
    assert ops.eleResponse(1, "constraintViolation")[0] == pytest.approx(1.0e-3, rel=1e-6)
    assert ops.eleResponse(1, "penaltyEnergy")[0] == pytest.approx(
        0.5 * ku * (1.0e-3) ** 2, rel=1e-6)


# --------------------------------------------------- 7. bipenalty (explicit dt_cr)
def test_bipenalty_dtcr_self_report():
    """-bipenalty -dtcr dt: m_p = k_eff*(dt/2)^2 and the self-reported dt_cr == dt."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1)
    Ku, dt = 1.0e6, 1.0e-3
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", Ku,
                "-bipenalty", "-dtcr", dt)
    assert ops.eleResponse(1, "mpenalty")[0] == pytest.approx(Ku * (dt / 2.0) ** 2, rel=1e-9)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(dt, rel=1e-9)
