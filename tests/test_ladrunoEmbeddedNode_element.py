"""LadrunoEmbeddedNode (general node-to-host coupling element) — Zone-A battery.

Phase 1 = the ISOTROPIC translational tie (U). The element ties a constrained node's
translations to a host element's nodes via shape-function weights N_i (penalty/AL,
`-k auto`, bipenalty), sharing LadrunoEmbeddedKernel with LadrunoEmbeddedRebar. See
SRC/element/ladrunoEmbeddedNode/ and Ladruno_implementation/23_ladruno_embedded_node_adr.md.

Like the rebar battery, these checks pin the coupling MECHANICS without a host continuum:
we play the "host" with explicit nodes whose displacements we prescribe (fix), so the gap
g = u_c - sum N_i u_host is fully controlled. The -host/-xi legs use a real LadrunoBrick.
"""
import math

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
def test_bipenalty_dtcr_host_reduced_mass():
    """ADR review #2 — with a QUERYABLE host that carries mass, the self-reported dt_cr
    tightens to the REDUCED-mass value 2√(μ/k_eff), μ = m_p·M_h/(m_p+M_h), BELOW the
    slave-only 2√(m_p/k_eff): the penalty stiffness loads the (penalty-massless) host DOFs,
    so the true relative-mode step is governed by the reduced mass, not m_p alone."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.node(3, 1.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1); ops.fix(3, 1, 1, 1)
    rho, A = 2.0, 0.5
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    ops.element("Truss", 100, 2, 3, A, 1, "-rho", rho)        # host element carrying mass
    Ku, dt = 1.0e6, 1.0e-3
    ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-shape", 1.0, 0.0,
                "-k", Ku, "-bipenalty", "-dtcr", dt)
    m_p = Ku * (dt / 2.0) ** 2
    slave_only = 2.0 * math.sqrt(m_p / Ku)                    # == dt (the old, unconservative report)
    reported = ops.eleResponse(1, "dtcr")[0]
    assert 0.0 < reported < slave_only                        # host mass tightens the bound
    # Truss -rho is mass PER UNIT LENGTH, lumped ρ·L/2 per node ⇒ M_h = min positive host
    # mass diagonal = ρ·L/2 (L = 1.0 here); A does not enter the Truss mass.
    Mh = rho * 1.0 / 2.0
    mu = m_p * Mh / (m_p + Mh)
    assert reported == pytest.approx(2.0 * math.sqrt(mu / Ku), rel=1e-6)


def test_k_nonnumeric_rejected():
    """ADR review #4 — `-k` followed by a non-numeric token (a forgotten value) is rejected
    at parse time, instead of atof silently setting Ku=0 (a disabled tie)."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    try:
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", "-enforce", "penalty")
    except Exception:
        pass
    assert 1 not in (ops.getEleTags() or [])


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


# =================================================== Phase 1b — pressure tie (UP)
def _up_host_node(ku=1.0e5, kp=2.0e5, ndm=3):
    """cNode (1) tied to a single coincident u-p host node (2, all DOFs fixed). u-p
    nodes carry ndf = ndm+1 with the pressure DOF at index ndm (3D: dof 4, 2D: dof 3)."""
    ndf = ndm + 1
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)
    ops.node(1, *([0.0] * ndm))
    ops.node(2, *([0.0] * ndm))
    ops.fix(2, *([1] * ndf))          # host fully fixed (u AND p)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-pressure", "-kp", kp)


def _push_pressure(p_dof, p_val):
    """Prescribe the cNode pressure DOF; translations stay free (spring to 0)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, p_dof, p_val)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 20); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()


def test_pressure_tie_3d():
    """3D u-p (ndf=4): the pressure DOF (4) couples as t_p = K_p*(p_c - sum N_i p_host)."""
    ku, kp = 1.0e5, 2.0e5
    _up_host_node(ku=ku, kp=kp, ndm=3)
    p_val = 3.0e-3
    _push_pressure(4, p_val)                     # prescribe cNode pressure (dof 4)
    assert abs(ops.nodeReaction(1)[3]) == pytest.approx(kp * p_val, rel=1e-6)
    assert ops.eleResponse(1, "pgap")[0] == pytest.approx(p_val, rel=1e-6)
    assert ops.eleResponse(1, "pforce")[0] == pytest.approx(kp * p_val, rel=1e-6)


def test_pressure_tie_2d():
    """2D u-p (ndf=3): -pressure disambiguates dof 3 as PRESSURE (not drilling Rz)."""
    ku, kp = 1.0e5, 2.0e5
    _up_host_node(ku=ku, kp=kp, ndm=2)
    p_val = 3.0e-3
    _push_pressure(3, p_val)                     # pressure DOF = index ndm = dof 3
    assert abs(ops.nodeReaction(1)[2]) == pytest.approx(kp * p_val, rel=1e-6)
    assert ops.eleResponse(1, "pgap")[0] == pytest.approx(p_val, rel=1e-6)


def test_pressure_does_not_disturb_translations():
    """The pressure tie is on a SEPARATE DOF (index ndm); the isotropic U coupling is
    unchanged (the ndf>ndm scatter writes pressure/translation to disjoint slots)."""
    ku, kp = 1.0e5, 2.0e5
    _up_host_node(ku=ku, kp=kp, ndm=3)
    R = _push1(1, 2.0e-3)                         # push x; y,z,p free -> 0
    assert abs(R[0]) == pytest.approx(ku * 2.0e-3, rel=1e-6)
    assert abs(R[1]) < 1e-3 and abs(R[2]) < 1e-3
    assert abs(R[3]) < 1e-3                       # pressure DOF untouched by a U push


def test_pressure_force_split_by_N():
    """Partition of unity on the PRESSURE DOF: two u-p hosts, N=[0.4,0.6]."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for nd in (1, 2, 3):
        ops.node(nd, 0.0, 0.0, 0.0)
    for nd in (2, 3):
        ops.fix(nd, 1, 1, 1, 1)
    ops.element("LadrunoEmbeddedNode", 1, 1, 2, 2, 3, "-shape", 0.4, 0.6,
                "-k", 1.0e5, "-pressure", "-kp", 2.0e5)
    _push_pressure(4, 1.0e-3)
    fr = ops.nodeReaction(1)[3]
    f2 = ops.nodeReaction(2)[3]
    f3 = ops.nodeReaction(3)[3]
    assert (fr + f2 + f3) == pytest.approx(0.0, abs=1e-6 * abs(fr))
    assert f2 == pytest.approx(0.4 * (-fr), rel=1e-6)
    assert f3 == pytest.approx(0.6 * (-fr), rel=1e-6)


def test_pressure_degrades_when_not_up():
    """-pressure requested but nodes are NOT u-p (ndf=3 in 3D) -> degrade to U-only."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e5,
                "-pressure", "-kp", 1.0e5)
    R = _push1(1, 2.0e-3)
    assert abs(R[0]) == pytest.approx(1.0e5 * 2.0e-3, rel=1e-6)   # still works as U
    assert ops.eleResponse(1, "pgap")[0] == 0.0                   # UP inactive


# =================================================== Phase 2 — rotation tie (UR)
# UR ties the constrained node's ROTATIONS to the host CONTINUUM rotation
# θ = ½ curl(u) = skew(∇u) at the embedded point, read from the host cartesian
# gradients ∂N_i/∂x. As in Phase 1 we pin the MECHANICS without a host continuum by
# supplying ∂N/∂x explicitly via -dNdx (the gradient analog of -shape) and prescribing
# the host translations; the real-host legs use a LadrunoBrick / BezierTet10 + -xi.

def _solve_static():
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-11, 30); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0


def test_rotation_2d_drilling_via_dNdx():
    """2D drilling: θ_z = ½(∂N/∂x·u_y − ∂N/∂y·u_x). Host translation (a,b) and ∂N/∂x
    are prescribed; the cNode drilling R_z is prescribed; check the moment reaction =
    K_r·(R_z − θ_host) and the rgap response (signed)."""
    ku, kr = 1.0e5, 1.0e6
    gx, gy = 1.5, 2.5
    a, b, rz = 0.02, 0.03, 0.05
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)        # cNode (u_x,u_y,R_z)
    ops.node(2, 0.0, 0.0)        # host
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-rot", "-kr", kr, "-dNdx", gx, gy)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, a); ops.sp(2, 2, b); ops.sp(2, 3, 0.0)   # host fully prescribed
    ops.sp(1, 3, rz)                                       # cNode drilling; u_x,u_y free
    _solve_static()
    ops.reactions()
    thz = 0.5 * (gx * b - gy * a)
    assert abs(ops.nodeReaction(1)[2]) == pytest.approx(kr * abs(rz - thz), rel=1e-6)
    assert ops.eleResponse(1, "rgap")[0] == pytest.approx(rz - thz, rel=1e-6)


def test_rotation_3d_continuum_via_dNdx():
    """3D: θ_host = ½ Σ_i (∇N_i × u_i). One host node (N=1, ∇N=(gx,gy,gz)), host
    translation (a,b,c); prescribe the three cNode rotations; check each moment
    reaction = K_r·(θ_c − θ_host) and the 3-vector rgap (signed)."""
    ku, kr = 1.0e5, 1.0e6
    gx, gy, gz = 1.0, 2.0, 3.0
    a, b, c = 0.01, 0.0, 0.0
    rx, ry, rz = 0.1, 0.0, 0.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)   # cNode (u + r)
    ops.node(2, 0.0, 0.0, 0.0)   # host
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-rot", "-kr", kr, "-dNdx", gx, gy, gz)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for d, val in ((1, a), (2, b), (3, c), (4, 0.0), (5, 0.0), (6, 0.0)):
        ops.sp(2, d, val)                      # host fully prescribed
    ops.sp(1, 4, rx); ops.sp(1, 5, ry); ops.sp(1, 6, rz)   # cNode rotations; u free
    _solve_static()
    ops.reactions()
    thx = 0.5 * (gy * c - gz * b)
    thy = 0.5 * (gz * a - gx * c)
    thz = 0.5 * (gx * b - gy * a)
    R = ops.nodeReaction(1)
    assert abs(R[3]) == pytest.approx(kr * abs(rx - thx), rel=1e-6, abs=1e-9)
    assert abs(R[4]) == pytest.approx(kr * abs(ry - thy), rel=1e-6, abs=1e-9)
    assert abs(R[5]) == pytest.approx(kr * abs(rz - thz), rel=1e-6, abs=1e-9)
    rg = ops.eleResponse(1, "rgap")
    assert rg[0] == pytest.approx(rx - thx, rel=1e-6, abs=1e-9)
    assert rg[1] == pytest.approx(ry - thy, rel=1e-6, abs=1e-9)
    assert rg[2] == pytest.approx(rz - thz, rel=1e-6, abs=1e-9)


def test_rotation_3d_brick_shear_field():
    """Real LadrunoBrick host (validates getInterpolationGradients + the -xi auto-
    gradient path + mixed ndf): prescribe a shear field u_x = α·y on the 8 host nodes
    ⇒ host continuum rotation θ = (0, 0, −α/2). The embedded node's rotations (free,
    no applied moment) follow the host continuum rotation through the penalty."""
    alpha, kr = 0.01, 1.0e6
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)            # 3-dof solid nodes
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.model("basic", "-ndm", 3, "-ndf", 6)            # switch ndf (no wipe): 6-dof cNode
    ops.node(1, 0.5, 0.5, 0.5)
    ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
                "-k", 1.0e5, "-rot", "-kr", kr)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _CUBE.items():
        ops.sp(tag, 1, alpha * y); ops.sp(tag, 2, 0.0); ops.sp(tag, 3, 0.0)
    _solve_static()                                      # cNode 6 DOFs free
    assert ops.nodeDisp(1, 4) == pytest.approx(0.0, abs=1e-6)        # θ_x
    assert ops.nodeDisp(1, 5) == pytest.approx(0.0, abs=1e-6)        # θ_y
    assert ops.nodeDisp(1, 6) == pytest.approx(-alpha / 2, rel=1e-4)  # θ_z = −α/2


# straight-sided reference tet (vertices + edge midpoints, BezierTet10 node order)
_TET_V = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]


def _mid(a, b):
    return tuple(0.5 * (a[k] + b[k]) for k in range(3))


# N1..N4 vertices, N5 e(1-2), N6 e(2-3), N7 e(1-3), N8 e(1-4), N9 e(3-4), N10 e(2-4)
_TET = {
    21: _TET_V[0], 22: _TET_V[1], 23: _TET_V[2], 24: _TET_V[3],
    25: _mid(_TET_V[0], _TET_V[1]), 26: _mid(_TET_V[1], _TET_V[2]),
    27: _mid(_TET_V[0], _TET_V[2]), 28: _mid(_TET_V[0], _TET_V[3]),
    29: _mid(_TET_V[2], _TET_V[3]), 30: _mid(_TET_V[1], _TET_V[3]),
}


def test_rotation_3d_bezierTet10_shear_field():
    """Real BezierTet10 host (validates its getInterpolationGradients override): same
    shear field u_x = α·y ⇒ θ_z = −α/2 at the centroid (L1=L2=L3=0.25). The embedded
    node's drilling rotation follows the host continuum rotation."""
    alpha, kr = 0.01, 1.0e6
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _TET.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("BezierTet10", 100, *_TET.keys(), 1)
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    cen = tuple(sum(v[k] for v in _TET_V) / 4.0 for k in range(3))
    ops.node(1, *cen)
    ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.25, 0.25, 0.25,
                "-k", 1.0e5, "-rot", "-kr", kr)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _TET.items():
        ops.sp(tag, 1, alpha * y); ops.sp(tag, 2, 0.0); ops.sp(tag, 3, 0.0)
    _solve_static()
    assert ops.nodeDisp(1, 6) == pytest.approx(-alpha / 2, rel=1e-4)
    assert ops.nodeDisp(1, 4) == pytest.approx(0.0, abs=1e-6)
    assert ops.nodeDisp(1, 5) == pytest.approx(0.0, abs=1e-6)


def test_rot_and_pressure_mutually_exclusive():
    """ADR 23 M5/D7-1 — parse-time guard: -rot and -pressure together are rejected
    (the extra DOF is rotation OR pressure, ambiguous in 2D ndf=3); the element is
    NOT created."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    try:
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e5,
                    "-rot", "-dNdx", 1.0, 0.0, 0.0, "-pressure")
    except Exception:
        pass
    assert 1 not in (ops.getEleTags() or [])


def test_rot_degrades_when_no_rotation_dofs():
    """-rot but the cNode has no rotation DOFs (3D ndf=3) ⇒ degrade to U-only (warn);
    the translational tie still works and rgap is inactive (zero)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e5,
                "-rot", "-dNdx", 1.0, 0.0, 0.0)
    R = _push1(1, 2.0e-3)
    assert abs(R[0]) == pytest.approx(1.0e5 * 2.0e-3, rel=1e-6)   # U still works
    assert ops.eleResponse(1, "rgap")[0] == 0.0                   # UR inactive


def test_rotation_bipenalty_per_class_dtcr():
    """ADR 23 D5 / M1·ES-1 — per-DOF-class (stiffness, inertia) pairs: m_p = K_u·(dt/2)²
    (translation), I_p = K_r·(dt/2)² (rotation), so the self-reported dt_cr = min over
    classes = dt (the rotation mode is bounded too, not left free)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    ku, kr, dt = 1.0e6, 4.0e6, 1.0e-3
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-rot", "-kr", kr, "-dNdx", 1.0, 0.0, 0.0, "-bipenalty", "-dtcr", dt)
    assert ops.eleResponse(1, "mpenalty")[0] == pytest.approx(ku * (dt / 2.0) ** 2, rel=1e-9)
    assert ops.eleResponse(1, "dtcr")[0] == pytest.approx(dt, rel=1e-9)


def test_enforce_al_drives_rotation_gap_to_zero():
    """AL on the rotation tie: cNode R_z FREE under an applied moment M, held by K_r to
    a fixed host (θ_host=0). Penalty leaves a residual M/K_r; AL's per-step Uzawa drives
    it ~0 at moderate K_r."""
    M, Kr = 50.0, 1.0e4

    def _build(enf):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
        ops.fix(2, 1, 1, 1, 1, 1, 1)            # host fully fixed ⇒ θ_host = 0
        ops.fix(1, 1, 1, 1, 1, 1, 0)            # cNode free in R_z only
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e6,
                    "-rot", "-kr", Kr, "-dNdx", 1.0, 0.0, 0.0, "-enforce", enf)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, 0.0, 0.0, 0.0, 0.0, M)
        ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
        ops.test("NormDispIncr", 1e-12, 30); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0); ops.analysis("Static")

    _build("penalty")
    assert ops.analyze(1) == 0
    gap_pen = abs(ops.nodeDisp(1, 6))
    assert gap_pen == pytest.approx(M / Kr, rel=1e-3)

    _build("al")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(8) == 0
    gap_al = abs(ops.nodeDisp(1, 6))
    assert gap_al < 0.05 * gap_pen


# =============================================== Phase 2b — D9 material interface
# Any local direction (normal e_0, tangents e_1/e_2) can carry a UNIAXIAL material
# instead of the bare penalty K_u — turning the translational tie into a node↔host
# INTERFACE: cohesive, unilateral gap, elastic bedding, bond. We pin the mechanics with
# a single coincident fixed host (N=[1], host fixed ⇒ g = u_c) and an axis-aligned
# -normal so the local frame coincides with the global axes (e_0=normal, e_1/e_2 the
# tangents). v1 uses the REFERENCE frame (-corot deferred).

def test_matN_elastic_reproduces_penalty():
    """-matN Elastic(K) on the normal reproduces the bare penalty (force = K·g_n); a
    direction with no -mat* keeps the K_u penalty (the tangents here)."""
    ku, ke = 1.0e5, 2.0e5
    def _build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", 10, ke)
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                    "-normal", 1.0, 0.0, 0.0, "-matN", 10)
    _build(); R = _push1(1, 2.0e-3)                       # push x = normal (Elastic ke)
    assert abs(R[0]) == pytest.approx(ke * 2.0e-3, rel=1e-6)
    _build(); R = _push1(2, 1.5e-3)                       # push y = tangent (penalty ku)
    assert abs(R[1]) == pytest.approx(ku * 1.5e-3, rel=1e-6)


def test_matN_ent_unilateral_gap():
    """ENT (elastic-no-tension) on the normal: opens (≈0 force) in tension/separation,
    holds stiff in compression — the unilateral-contact/gap behaviour."""
    ku, E = 1.0e5, 3.0e5
    def _build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
        ops.uniaxialMaterial("ENT", 11, E)
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                    "-normal", 1.0, 0.0, 0.0, "-matN", 11)
    _build(); R = _push1(1, 2.0e-3)                       # +x tension -> gap opens
    assert abs(R[0]) < 1.0
    _build(); R = _push1(1, -2.0e-3)                      # -x compression -> holds
    assert abs(R[0]) == pytest.approx(E * 2.0e-3, rel=1e-6)


def test_matT_elasticpp_fixed_slip():
    """ElasticPP on a tangent gives a FIXED slip force F_y past yield (the documented
    UNCOUPLED approximate-friction behaviour — not μ·N; rigorous → LadrunoContact)."""
    ku, E, epsy = 1.0e5, 1.0e6, 1.0e-3
    Fy = E * epsy
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("ElasticPP", 12, E, epsy)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 1.0, 0.0, 0.0, "-matT1", 12)
    R = _push1(2, 5.0e-3)                                 # push y (tangent e_1) past yield
    assert abs(R[1]) == pytest.approx(Fy, rel=1e-4)


def test_material_interface_local_responses():
    """localGap / localForce / normal responses report the frame components."""
    ku, ke = 1.0e5, 2.0e5
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 13, ke)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 1.0, 0.0, 0.0, "-matN", 13)
    _push1(1, 2.0e-3)
    assert ops.eleResponse(1, "localGap")[0] == pytest.approx(2.0e-3, rel=1e-6)
    assert ops.eleResponse(1, "localForce")[0] == pytest.approx(ke * 2.0e-3, rel=1e-6)
    assert ops.eleResponse(1, "normal")[0] == pytest.approx(1.0, rel=1e-9)


def test_al_with_material_direction_reprojection():
    """AL + one material direction (normal Elastic) + one penalty direction (a tangent
    under load): AL drives the PENALTY gap → 0; the material direction carries physical
    force and the multiplier is re-projected off it (M4/D9-3)."""
    F, ku, ke = 50.0, 1.0e4, 2.0e5

    def _build(enf):
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
        ops.fix(1, 0, 0, 1)                               # cNode free in x,y; z fixed
        ops.uniaxialMaterial("Elastic", 14, ke)
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                    "-normal", 0.0, 1.0, 0.0, "-matN", 14, "-enforce", enf)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        ops.load(1, F, 0.0, 0.0)                          # load along x = penalty tangent
        ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
        ops.test("NormDispIncr", 1e-12, 30); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0); ops.analysis("Static")

    _build("penalty")
    assert ops.analyze(1) == 0
    gap_pen = abs(ops.nodeDisp(1, 1))
    assert gap_pen == pytest.approx(F / ku, rel=1e-3)
    _build("al")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(8) == 0
    gap_al = abs(ops.nodeDisp(1, 1))
    assert gap_al < 0.05 * gap_pen
    assert abs(ops.nodeDisp(1, 2)) < 1e-6                 # material (y) dir stays put


def test_material_bipenalty_keff_from_initial_tangents():
    """D9 bipenalty k_eff = max(K_u, material initial tangents) (M6); m_p = k_eff·(dt/2)²."""
    ku, ke, dt = 1.0e6, 4.0e6, 1.0e-3
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 15, ke)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 1.0, 0.0, 0.0, "-matN", 15, "-bipenalty", "-dtcr", dt)
    assert ops.eleResponse(1, "mpenalty")[0] == pytest.approx(ke * (dt / 2.0) ** 2, rel=1e-9)


def test_material_interface_2d():
    """2D material interface: Elastic on the in-plane normal."""
    ku, ke = 1.0e5, 2.0e5
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 0.0); ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 16, ke)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 1.0, 0.0, "-matN", 16)
    R = _push1(1, 2.0e-3)
    assert abs(R[0]) == pytest.approx(ke * 2.0e-3, rel=1e-6)


def test_matN_requires_normal():
    """-matN without -normal is rejected at parse time (the frame is undefined)."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    ops.uniaxialMaterial("Elastic", 17, 1.0e5)
    try:
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e5,
                    "-matN", 17)
    except Exception:
        pass
    assert 1 not in (ops.getEleTags() or [])


# ====================================== Phase 2b v2 — D9 material-frame co-rotation (-corot)
# -corot makes the material frame track the host CONTINUUM rotation θ = skew(∇u) at the
# embedded point: frameCur = R(θ)·frame (3D Rodrigues / 2D drilling). A directional contact
# normal then follows the deformed host. It needs the host gradients ∂N/∂x (here via -dNdx,
# the gradient analog of -shape) and applies to the material interface only (the isotropic /
# penalty tie is already frame-objective). The dropped ∂e_d/∂u tangent term keeps the
# residual exact (tangent inexact — R7). As elsewhere we pin the mechanics by prescribing the
# host translations (⇒ θ_host is known and held), with the cNode free.

def test_corot_normal_rotates_2d():
    """2D: a prescribed host translation + ∂N/∂x induce a continuum drilling rotation
    θ_z = ½(∂N/∂x·u_y − ∂N/∂y·u_x); under -corot the material normal rotates with it
    (normal response = R(θ_z)·n = (cos θ_z, sin θ_z)). WITHOUT -corot it stays the
    reference normal (1,0) — the negative control proving -corot does the rotating."""
    ku, ke = 1.0e5, 2.0e5
    gx, gy = 0.0, 2.0
    a, b = 0.1, 0.0
    thz = 0.5 * (gx * b - gy * a)                 # = −gy·a/2 = −0.1

    def _build(corot):
        ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 0.0)
        ops.uniaxialMaterial("Elastic", 20, ke)
        args = ["LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 1.0, 0.0, "-matN", 20, "-dNdx", gx, gy]
        if corot:
            args.append("-corot")
        ops.element(*args)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, a); ops.sp(2, 2, b)          # host prescribed ⇒ θ_z = −0.1
        _solve_static()

    _build(corot=True)
    nrm = ops.eleResponse(1, "normal")
    assert nrm[0] == pytest.approx(math.cos(thz), rel=1e-4)
    assert nrm[1] == pytest.approx(math.sin(thz), rel=1e-4)

    _build(corot=False)                            # negative control: reference frame
    nrm = ops.eleResponse(1, "normal")
    assert nrm[0] == pytest.approx(1.0, abs=1e-9)
    assert nrm[1] == pytest.approx(0.0, abs=1e-9)


def test_corot_requires_material_mode():
    """-corot without any -mat* is rejected at parse time (the isotropic/penalty tie is
    already frame-objective; there is no frame to co-rotate)."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    try:
        ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e5,
                    "-normal", 1.0, 0.0, 0.0, "-dNdx", 1.0, 0.0, 0.0, "-corot")
    except Exception:
        pass
    assert 1 not in (ops.getEleTags() or [])


def test_corot_al_material_reprojection_3d():
    """ADR 23 item 11 — AL × material × corot. A host translation + ∂N/∂x induce a continuum
    rotation about z (θ = ½ g×u_host), so the material normal (reference y) co-rotates toward
    x. Under -enforce al the multiplier λ augments only the PENALTY directions and is
    re-projected off the (rotated) material direction each commit (M4/D9-3): λ·e_0^cur stays
    ~0 even though frame rotation alone could leak it there."""
    F, ku, ke = 50.0, 1.0e4, 2.0e5
    gx, gy, gz = 2.0, 0.0, 0.0
    a, b, c = 0.0, 0.05, 0.0                       # θ = ½ g×u = (0,0,0.05) ⇒ θ_z = 0.05
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0)
    ops.uniaxialMaterial("Elastic", 21, ke)
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", ku,
                "-normal", 0.0, 1.0, 0.0, "-matN", 21, "-dNdx", gx, gy, gz,
                "-corot", "-enforce", "al")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, a); ops.sp(2, 2, b); ops.sp(2, 3, c)   # host prescribed ⇒ θ_z = 0.05
    ops.fix(1, 0, 0, 1)                            # cNode free in x,y; z fixed
    ops.load(1, F, 0.0, 0.0)                       # load along x (≈ penalty tangent e_1)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-11, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(8) == 0                      # constant-load Uzawa updates

    nrm = ops.eleResponse(1, "normal")             # the frame really rotated (corot active)
    assert abs(nrm[0]) == pytest.approx(math.sin(0.05), rel=1e-2)
    lam = ops.eleResponse(1, "lambda")             # λ purged off the rotated material normal
    leak = sum(lam[k] * nrm[k] for k in range(3))
    lam_mag = math.sqrt(sum(l * l for l in lam))
    assert abs(leak) < 1e-6 * max(lam_mag, 1.0)


# =================================== Initial-gap (offset) capture — stress-free staged activation
# By default the element captures the gap at activation (setDomain) and drives ALL traction from
# the RELATIVE gap (g − g0), so an element added to an ALREADY-DEFORMED host (staged construction)
# is born stress-free instead of yanking the slave to chase the host's accumulated displacement
# (the parent ASDEmbeddedNodeElement m_U0 behavior). `-absolute` opts out (the legacy absolute
# tie / a deliberate snap-to-host). We deform a real LadrunoBrick (stage 1), then tie a fresh
# slave node into the deformed host (stage 2) and inspect the element at activation.

def _deformed_brick_host(load=10.0):
    """LadrunoBrick (tag 100): fix the z=0 face, shear-load the z=1 face in +x, solve + freeze.
    Leaves the domain holding the DEFORMED, committed host ready for a stage-2 element add."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
    for tag in (11, 12, 13, 14):                  # z=0 face fully fixed
        ops.fix(tag, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag in (15, 16, 17, 18):                  # shear the top face in +x
        ops.load(tag, load, 0.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 30); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)                    # freeze stage-1 load; host stays deformed


def test_staged_activation_born_stress_free():
    """Default capture: a slave tied into the deformed host at -xi centroid is born with ZERO
    force and ZERO relative gap; the captured offset equals −(host centroid disp)."""
    Ku = 1.0e6
    _deformed_brick_host()
    uc = sum(ops.nodeDisp(t, 1) for t in _CUBE) / 8.0     # host centroid x-disp (weights 0.125)
    assert abs(uc) > 1e-4                                  # the host really moved
    ops.node(1, 0.5, 0.5, 0.5)                             # fresh slave (u = 0)
    ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0, "-k", Ku)
    assert abs(ops.eleResponse(1, "force")[0]) < 1e-3      # born force-free (no jolt)
    assert abs(ops.eleResponse(1, "gap")[0]) < 1e-9        # relative gap zero at activation
    assert ops.eleResponse(1, "initGap")[0] == pytest.approx(-uc, rel=1e-6)   # offset captured
    ops.integrator("LoadControl", 0.0)                     # one constant-load step
    assert ops.analyze(1) == 0
    assert abs(ops.nodeDisp(1, 1)) < 1e-4                  # slave NOT yanked (tracks increments)


def test_staged_activation_absolute_yanks():
    """-absolute (negative control): the same staged add WITHOUT capture is born with the full
    accumulated gap → a large jolt force, and the slave is yanked to the deformed host point."""
    Ku = 1.0e6
    _deformed_brick_host()
    uc = sum(ops.nodeDisp(t, 1) for t in _CUBE) / 8.0
    ops.node(1, 0.5, 0.5, 0.5)
    ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
                "-k", Ku, "-absolute")
    assert abs(ops.eleResponse(1, "gap")[0]) == pytest.approx(abs(uc), rel=1e-6)     # full gap
    assert abs(ops.eleResponse(1, "force")[0]) == pytest.approx(Ku * abs(uc), rel=1e-3)  # jolt
    assert ops.eleResponse(1, "initGap")[0] == 0.0          # no capture
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(1, 1) == pytest.approx(uc, rel=2e-2)   # slave yanked to the host point


def test_staged_capture_noop_when_undeformed():
    """Capture is a no-op when the element is added at the UNDEFORMED state (g0 = 0): the
    isotropic tie behaves exactly as before — so the whole v1 battery stays byte-identical."""
    ku = 1.0e5
    _single_host_node(ku=ku)                                # added at rest ⇒ g0 = 0
    assert ops.eleResponse(1, "initGap")[0] == 0.0
    R = _push1(1, 2.0e-3)
    assert abs(R[0]) == pytest.approx(ku * 2.0e-3, rel=1e-6)   # identical to the absolute tie
    assert ops.eleResponse(1, "gap")[0] == pytest.approx(2.0e-3, rel=1e-6)


def test_staged_serialization_roundtrip():
    """sendSelf/recvSelf round-trip (FE_Datastore) of the element carrying a captured g0 — proves
    the bumped header/payload buffers are self-consistent and the broker reconstructs the element.
    Skips if this build lacks database support."""
    from _testbed.roundtrip import database_roundtrip

    def _build():
        Ku = 1.0e6
        ops.wipe()                       # the roundtrip helper calls _build before any wipe
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for tag, (x, y, z) in _CUBE.items():
            ops.node(tag, x, y, z)
        for tag in (11, 12, 13, 14):
            ops.fix(tag, 1, 1, 1)
        ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
        ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        for tag in (15, 16, 17, 18):
            ops.load(tag, 10.0, 0.0, 0.0)
        ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
        ops.test("NormDispIncr", 1e-10, 30); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0); ops.analysis("Static")
        assert ops.analyze(1) == 0
        ops.node(1, 0.5, 0.5, 0.5)                          # stage-2 slave ⇒ nonzero captured g0
        ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0, "-k", Ku)

    database_roundtrip(_build, probe_nodes=[15, 17], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, "initGap")))


# =================================================== implicit-transient damping (getDamp crash)
def test_transient_newmark_smoke():
    """REGRESSION — the pure penalty coupling no-ops setRayleighDampingFactors, which suppresses
    the base Element::getDamp LAZY allocation of its Rayleigh scratch. Without an explicit
    getDamp/getRayleighDampingForces override the base indexes theMatrices[-1] (index stays -1)
    and HARD-CRASHES on the first IMPLICIT-transient step: Newmark's velocity coefficient
    c2 = γ/(βΔt) is always nonzero, so FE_Element::addCtoTang → Element::getDamp fires every step
    (the residual addD_Force path hits it too). The rest of this battery is quasi-static (Static)
    or — for the bipenalty legs — explicit, whose effective tangent dodges getDamp when there is
    no Rayleigh damping, so they all miss it. Build with -bipenalty so the cNode carries m_p, run
    a few Newmark steps, and assert they return 0 (would segfault/exit before the getDamp fix)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)              # constrained node (free) — carries m_p
    ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)   # host fixed ⇒ g = u_c
    Ku, dt = 1.0e6, 1.0e-3
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", Ku,
                "-bipenalty", "-dtcr", dt)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, 0.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-8, 30); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    assert ops.analyze(5, dt) == 0          # crashes in getDamp without the override
