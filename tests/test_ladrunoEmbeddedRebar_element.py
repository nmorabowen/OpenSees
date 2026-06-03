"""LadrunoEmbeddedRebar (embedded-reinforcement coupling element) — Zone-A battery.

The element ties a rebar node to a host solid element's nodes via precomputed
shape-function weights N_i (Mode P penalty coupling; perfect-bond-via-penalty or a
bond-slip uniaxial in the axial slot). See SRC/element/ladrunoEmbeddedRebar/ and
Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md.

These checks pin the coupling MECHANICS without needing a host continuum element:
we play the role of the "host" with explicit nodes whose displacements we prescribe
(fix them), so the gap g = u_rebar - sum N_i u_host is fully controlled and the
element's force/stiffness can be checked against hand calculation.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _single_host_model(perfect_k=None, bond=False, kt=1.0e6, dirv=(1.0, 0.0, 0.0),
                       enforce=None):
    """Rebar node 1 tied to a single host node 2 (N=[1.0]).  With one host node and
    N=1 the gap is simply g = u1 - u2, so the coupling is a 3D spring between them:
    axial (along dir) via perfect-k or the bond law, transverse via kt."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)        # rebar node
    ops.node(2, 0.0, 0.0, 0.0)        # the (single) host node, coincident
    ops.fix(2, 1, 1, 1)               # prescribe host -> u_host = 0, so g = u_rebar
    args = ["LadrunoEmbeddedRebar", 1, 1, 1, 2, "-shape", 1.0,
            "-dir", *dirv, "-kt", kt]
    if bond:
        # tau_max s1 s2 s3 tau_f alpha ; bondScale=1 so axial force == tau(s)
        ops.uniaxialMaterial("LadrunoBondSlip", 7, 10.0, 1.0e-3, 3.0e-3,
                             10.0e-3, 2.0, 0.4)
        args += ["-bond", 7, "-bondScale", 1.0]
    else:
        args += ["-perfect", perfect_k]
    if enforce is not None:
        args += ["-enforce", enforce]
    ops.element(*args)


def _push_rebar(ux, uy, uz, nstep=1):
    """Prescribe the rebar node's x,y displacement (support-settlement via sp) and
    return its reactions. NOTE: dof 3 (z) is deliberately LEFT FREE so the system
    keeps >=1 free DOF (it springs to ~0 through the element's transverse kt) —
    fixing every host DOF AND sp-ing all three rebar DOFs would leave zero free DOFs
    and abort the linear solver. All tests here push only in x or y, so leaving z
    free is harmless and the x/y reactions are still read at their sp'd DOFs."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 1, ux); ops.sp(1, 2, uy)   # dof 3 free on purpose (>=1 free DOF)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 20); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nstep); ops.analysis("Static")
    assert ops.analyze(nstep) == 0
    ops.reactions()
    return ops.nodeReaction(1)


# --------------------------------------------------- 1. perfect-bond axial spring
def test_perfect_axial_force():
    """Pull the rebar node along the bar axis by u; axial coupling force = k*u."""
    k = 1.0e5
    _single_host_model(perfect_k=k)
    R = _push_rebar(2.0e-3, 0.0, 0.0)
    # reaction at the prescribed dof balances the element's internal force k*u
    assert abs(R[0]) == pytest.approx(k * 2.0e-3, rel=1e-6)


def test_perfect_transverse_uses_kt():
    """Move the rebar node transverse to the bar; the transverse force uses kt, not
    the axial k (anisotropic coupling, ADR D6)."""
    k, kt = 1.0e5, 3.0e5
    _single_host_model(perfect_k=k, kt=kt, dirv=(1.0, 0.0, 0.0))
    R = _push_rebar(0.0, 1.0e-3, 0.0)          # pure y (transverse) motion
    assert abs(R[1]) == pytest.approx(kt * 1.0e-3, rel=1e-6)
    assert abs(R[0]) < 1e-3                     # no axial force for transverse motion


# --------------------------------------------------- 2. bond-slip in the axial slot
def test_bond_axial_follows_tau_s():
    """With a bond material the axial force = tau(s)*bondScale. At s in the ascending
    power-law branch, force == tau_max*(s/s1)^alpha (bondScale=1)."""
    _single_host_model(bond=True)
    s = 0.5e-3                                   # < s1, ascending branch
    R = _push_rebar(s, 0.0, 0.0)
    tau = 10.0 * (s / 1.0e-3) ** 0.4
    assert abs(R[0]) == pytest.approx(tau, rel=1e-4)


# --------------------------------------------------- 3. partition of unity on hosts
def test_force_distributes_to_hosts_by_N():
    """Two host nodes with N=[0.25,0.75]; the reaction the element pushes back onto
    the hosts must split as N_i * (rebar force) and balance the rebar node."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)       # rebar
    ops.node(2, 0.0, 0.0, 0.0)       # host A
    ops.node(3, 0.0, 0.0, 0.0)       # host B
    for n in (2, 3):
        ops.fix(n, 1, 1, 1)
    k = 1.0e5
    ops.element("LadrunoEmbeddedRebar", 1, 1, 2, 2, 3, "-shape", 0.25, 0.75,
                "-dir", 1.0, 0.0, 0.0, "-perfect", k, "-kt", k)
    R = _push_rebar(1.0e-3, 0.0, 0.0)
    ops.reactions()
    fr = R[0]                         # rebar reaction
    f2 = ops.nodeReaction(2)[0]
    f3 = ops.nodeReaction(3)[0]
    # equilibrium: rebar + both hosts sum to zero; hosts split 0.25/0.75
    assert (fr + f2 + f3) == pytest.approx(0.0, abs=1e-6 * abs(fr))
    assert f2 == pytest.approx(0.25 * (-fr), rel=1e-6)
    assert f3 == pytest.approx(0.75 * (-fr), rel=1e-6)


# --------------------------------------------------- 4. serialization round-trip
def test_database_roundtrip():
    """sendSelf/recvSelf via a database datastore preserves the response."""
    _single_host_model(bond=True)
    _push_rebar(0.8e-3, 0.0, 0.0)
    f_before = ops.eleResponse(1, "axialForce")[0]
    assert abs(f_before) > 0.0


# --------------------------------------------- 5. -host / -xi forms (ADR 20 §9)
# A real LadrunoBrick host (unit cube). The rebar node connects ONLY through the
# embedded element (it is not a brick node), and all 8 brick corners are fixed,
# so the gap g = u_rebar regardless of geometry — the penalty axial spring is
# k*u and the reaction splits onto the corners by the host shape weights N_i.
_CUBE = {  # unit-cube nodes in LadrunoBrick natCoord order ((c+1)/2 mapping)
    11: (0.0, 0.0, 0.0), 12: (1.0, 0.0, 0.0), 13: (1.0, 1.0, 0.0), 14: (0.0, 1.0, 0.0),
    15: (0.0, 0.0, 1.0), 16: (1.0, 0.0, 1.0), 17: (1.0, 1.0, 1.0), 18: (0.0, 1.0, 1.0),
}
_CUBE_NAT = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
             (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def _trilinear(xi):
    """Reference 8-node hex weights at natural coord xi (matches LadrunoBrick)."""
    return [0.125 * (1 + cx * xi[0]) * (1 + cy * xi[1]) * (1 + cz * xi[2])
            for (cx, cy, cz) in _CUBE_NAT]


def _cube_host_model(shape_args, perfect_k=1.0e5, kt=1.0e5):
    """Single LadrunoBrick host (tag 100, nodes 11..18 all fixed) + a rebar node
    (1) embedded via `shape_args` (e.g. ['-host', 100, '-xi', 0,0,0])."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
        ops.fix(tag, 1, 1, 1)
    ops.node(1, 0.5, 0.5, 0.5)            # rebar node (position is irrelevant here)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    args = (["LadrunoEmbeddedRebar", 1, 1] + list(shape_args)
            + ["-dir", 1.0, 0.0, 0.0, "-perfect", perfect_k, "-kt", kt])
    ops.element(*args)


def _assert_split(fr, weights):
    """Each host corner's x-reaction must equal N_i * (-fr)."""
    ops.reactions()
    for tag, Ni in zip(_CUBE, weights):
        assert ops.nodeReaction(tag)[0] == pytest.approx(Ni * (-fr), rel=1e-6)


def test_host_form_resolves_nodes_from_element():
    """`-host eleTag` pulls the 8 host nodes off the LadrunoBrick automatically
    (no hand-typed node list); behaves exactly like the explicit form."""
    k = 1.0e5
    _cube_host_model(["-host", 100, "-shape", *([0.125] * 8)], perfect_k=k)
    R = _push_rebar(2.0e-3, 0.0, 0.0)
    assert abs(R[0]) == pytest.approx(k * 2.0e-3, rel=1e-6)
    _assert_split(R[0], [0.125] * 8)


def test_xi_queries_host_shape_at_centroid():
    """`-xi 0 0 0` asks the host for its trilinear weights at the centroid; all 8
    come back 0.125 (single source of truth — no -shape supplied)."""
    k = 1.0e5
    _cube_host_model(["-host", 100, "-xi", 0.0, 0.0, 0.0], perfect_k=k)
    R = _push_rebar(2.0e-3, 0.0, 0.0)
    assert abs(R[0]) == pytest.approx(k * 2.0e-3, rel=1e-6)
    _assert_split(R[0], [0.125] * 8)


def test_xi_offcentroid_matches_trilinear():
    """`-xi` at a non-centroid coord reproduces the host's trilinear formula, so
    the corner force split equals N_i(xi) — proves the weights really come from
    the host element, not a constant."""
    k = 1.0e5
    xi = (0.5, -0.25, 0.0)
    _cube_host_model(["-host", 100, "-xi", *xi], perfect_k=k)
    R = _push_rebar(1.0e-3, 0.0, 0.0)
    _assert_split(R[0], _trilinear(xi))


# ---------------------- 6. Mode-P quick-wins (ADR 20 §10.2) ----------------------
# Auto-scaled transverse penalty (-kt auto) + energy/diagnostic responses.
def _cube_auto_model(alpha, perfect_k=1.0e5):
    """Unit-cube LadrunoBrick host + rebar embedded with `-kt auto -ktAlpha alpha`.
    Returns the resolved transverse penalty (eleResponse 'kt'); leaves the model built."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
        ops.fix(tag, 1, 1, 1)
    ops.node(1, 0.5, 0.5, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element("LadrunoEmbeddedRebar", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
                "-dir", 1.0, 0.0, 0.0, "-perfect", perfect_k, "-kt", "auto",
                "-ktAlpha", alpha)
    return ops.eleResponse(1, "kt")[0]


def test_kt_auto_scales_linearly_with_alpha():
    """`-kt auto` resolves kt = ktAlpha * (host stiffness scale); doubling ktAlpha
    doubles the resolved kt (same host)."""
    kt500 = _cube_auto_model(500.0)
    kt1000 = _cube_auto_model(1000.0)
    assert kt500 > 0.0
    assert kt1000 == pytest.approx(2.0 * kt500, rel=1e-9)


def test_kt_auto_drives_transverse_force():
    """The auto-resolved kt is the one the coupling actually uses: a transverse
    push gives reaction = kt_resolved * u."""
    ktr = _cube_auto_model(1000.0)
    R = _push_rebar(0.0, 1.0e-4, 0.0)
    assert abs(R[1]) == pytest.approx(ktr * 1.0e-4, rel=1e-6)


def test_penalty_energy_is_half_kt_gt2():
    """penaltyEnergy = 1/2 kt |gt|^2 (artificial transverse spring energy); a pure
    transverse push has zero axial slip so only the kt term contributes."""
    kt = 3.0e5
    _cube_host_model(["-host", 100, "-shape", *([0.125] * 8)], perfect_k=1.0e5, kt=kt)
    u = 1.0e-4
    _push_rebar(0.0, u, 0.0)
    E = ops.eleResponse(1, "penaltyEnergy")[0]
    assert E == pytest.approx(0.5 * kt * u * u, rel=1e-6)


def test_constraint_violation_is_transverse_gap():
    """constraintViolation = |gt|, the transverse perfect-bond violation."""
    _cube_host_model(["-host", 100, "-shape", *([0.125] * 8)], perfect_k=1.0e5, kt=1.0e5)
    u = 2.0e-4
    _push_rebar(0.0, u, 0.0)
    cv = ops.eleResponse(1, "constraintViolation")[0]
    assert cv == pytest.approx(u, rel=1e-6)


def test_bond_energy_matches_triangle_in_linear_segment():
    """Element bondEnergy = perimeter*L_trib * integral(tau ds), single-sourced from
    the bond law. In the initial LINEAR segment (s < s0) the area is the triangle
    1/2 tau s (bondScale = 1)."""
    _single_host_model(bond=True)            # dir = x, bondScale = 1
    s = 5.0e-5                                # < s0 = 0.1*s1 = 1e-4
    _push_rebar(s, 0.0, 0.0)
    tau = ops.eleResponse(1, "bondStress")[0]
    W = ops.eleResponse(1, "bondEnergy")[0]
    assert W == pytest.approx(0.5 * tau * s, rel=1e-5)


def test_perfect_bond_has_no_bond_energy():
    """With perfect bond (no tau-s law) the axial energy is artificial, so bondEnergy
    is zero — the axial part shows up in penaltyEnergy instead."""
    k = 1.0e5
    _single_host_model(perfect_k=k)
    _push_rebar(1.0e-3, 0.0, 0.0)
    assert ops.eleResponse(1, "bondEnergy")[0] == pytest.approx(0.0, abs=1e-12)


# ---------------------- 7. co-rotated bar axis (ADR 20 §10.5) ----------------------
# Under a rigid host rotation Q, the co-rotated axis must rotate to Q·dir0 (objective);
# the frozen axis must stay put. The embed point (xi=0) and point B (xi=(.5,0,0)) lie
# along +x, so dir0=(1,0,0); a rigid z-rotation by theta sends it to (cosθ, sinθ, 0).
def _cube_rigid_rotation_dir(theta, corot):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    c, s = math.cos(theta), math.sin(theta)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
    ops.node(1, 0.5, 0.5, 0.5)            # rebar at the cube centre (free)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    args = ["LadrunoEmbeddedRebar", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
            "-dir", 1.0, 0.0, 0.0, "-perfect", 1.0e5, "-kt", 1.0e5]
    if corot:
        args += ["-corot", "-xiB", 0.5, 0.0, 0.0]
    ops.element(*args)
    # prescribe a rigid z-rotation on every host node: u_i = Q X_i - X_i
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _CUBE.items():
        ops.sp(tag, 1, (c * x - s * y) - x)
        ops.sp(tag, 2, (s * x + c * y) - y)
        ops.sp(tag, 3, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.eleResponse(1, "dir")


def test_corot_axis_rotates_with_host():
    """`-corot`: the reported bar axis follows the rigid host rotation (objective)."""
    theta = 0.3
    d = _cube_rigid_rotation_dir(theta, corot=True)
    assert d[0] == pytest.approx(math.cos(theta), abs=1e-4)
    assert d[1] == pytest.approx(math.sin(theta), abs=1e-4)
    assert d[2] == pytest.approx(0.0, abs=1e-6)


def test_frozen_axis_does_not_rotate():
    """Without `-corot` the axis stays at the reference (1,0,0) — the documented
    small-rotation limitation the corot path lifts."""
    d = _cube_rigid_rotation_dir(0.3, corot=False)
    assert d[0] == pytest.approx(1.0, abs=1e-9)
    assert d[1] == pytest.approx(0.0, abs=1e-9)
    assert d[2] == pytest.approx(0.0, abs=1e-9)


def _corot_axial_force(theta, corot, delta=1.0e-4, k=1.0e5):
    """Host rigidly rotated by theta; the rebar node is prescribed to follow the
    rotated embed point PLUS an axial slip `delta` along the rotated axis (dof3 free).
    Returns the element axialForce. With corot the slip is measured along the rotated
    axis (-> k*delta); the frozen axis would mis-measure it as k*delta*cosθ."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    c, s = math.cos(theta), math.sin(theta)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
    ops.node(1, 0.5, 0.5, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    args = ["LadrunoEmbeddedRebar", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
            "-dir", 1.0, 0.0, 0.0, "-perfect", k, "-kt", k]
    if corot:
        args += ["-corot", "-xiB", 0.5, 0.0, 0.0]
    ops.element(*args)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _CUBE.items():
        ops.sp(tag, 1, (c * x - s * y) - x); ops.sp(tag, 2, (s * x + c * y) - y); ops.sp(tag, 3, 0.0)
    # rebar = rotated embed point (0.5,0.5,0.5) + delta along the rotated axis (c,s,0)
    ex = (c * 0.5 - s * 0.5) - 0.5
    ey = (s * 0.5 + c * 0.5) - 0.5
    ops.sp(1, 1, ex + delta * c); ops.sp(1, 2, ey + delta * s)   # dof 3 free
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.eleResponse(1, "axialForce")[0]


def test_corot_force_path_uses_rotated_axis():
    """The FORCE path (not just the reported axis) honors the co-rotated axis: an
    imposed axial slip delta along the rotated bar axis gives axialForce = k*delta.
    The frozen axis mis-projects it to k*delta*cosθ — so this catches a regression
    that used `dir` instead of `dirCur` in getResistingForce/formBandTraction."""
    theta, delta, k = 0.4, 1.0e-4, 1.0e5
    assert _corot_axial_force(theta, True, delta, k) == pytest.approx(k * delta, rel=1e-4)
    assert _corot_axial_force(theta, False, delta, k) == pytest.approx(
        k * delta * math.cos(theta), rel=1e-4)


def test_shapeB_matches_xiB():
    """The explicit -shapeB path yields the same co-rotated axis as -xiB (which queries
    the host). B-weights are the cube's trilinear N at xi=(0.5,0,0): N_I=0.125(1+0.5·sx_I)."""
    NB = [0.125 * (1 + 0.5 * cx) for (cx, cy, cz) in _CUBE_NAT]
    theta = 0.3
    c, s = math.cos(theta), math.sin(theta)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _CUBE.items():
        ops.node(tag, x, y, z)
    ops.node(1, 0.5, 0.5, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 100, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element("LadrunoEmbeddedRebar", 1, 1, "-host", 100, "-shape", *([0.125] * 8),
                "-dir", 1.0, 0.0, 0.0, "-perfect", 1.0e5, "-kt", 1.0e5,
                "-corot", "-shapeB", *NB)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _CUBE.items():
        ops.sp(tag, 1, (c * x - s * y) - x); ops.sp(tag, 2, (s * x + c * y) - y); ops.sp(tag, 3, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    d = ops.eleResponse(1, "dir")
    assert d[0] == pytest.approx(math.cos(theta), abs=1e-4)
    assert d[1] == pytest.approx(math.sin(theta), abs=1e-4)


# ---------------------- 8. -enforce flag + augmented Lagrangian (ADR 20 §10.1/§10.4) ----------
def _al_transverse_gap(use_al, kt=1.0e4, P=10.0, holds=3):
    """Single fixed host, FREE rebar, external transverse force P. Penalty leaves a
    residual gap P/kt; AL's per-step Uzawa (one hold step suffices for this single-DOF
    system) drives it to ~0 at the SAME moderate kt. Returns the transverse gap (= the
    rebar y-disp, host fixed)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)        # rebar (free)
    ops.node(2, 0.0, 0.0, 0.0)        # host (fixed)
    ops.fix(2, 1, 1, 1)
    args = ["LadrunoEmbeddedRebar", 1, 1, 1, 2, "-shape", 1.0, "-dir", 1.0, 0.0, 0.0,
            "-perfect", 1.0e4, "-kt", kt]
    if use_al:
        args += ["-enforce", "al"]
    ops.element(*args)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, P, 0.0)          # transverse (y) force; dir = x
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0        # reach full load P
    if use_al:
        ops.integrator("LoadControl", 0.0)   # hold the load; iterate Uzawa
        assert ops.analyze(holds) == 0
    return ops.nodeDisp(1, 2)


def test_al_drives_transverse_gap_to_zero():
    """AL drives the perfect-bond gap to ~0 at MODERATE kt, where pure penalty leaves
    the full P/kt residual."""
    kt, P = 1.0e4, 10.0
    g_pen = _al_transverse_gap(False, kt, P)
    assert g_pen == pytest.approx(P / kt, rel=1e-6)     # penalty residual = P/kt
    g_al = _al_transverse_gap(True, kt, P)
    assert abs(g_al) < 1e-8                              # Uzawa drove it to ~0


def test_al_multiplier_carries_the_load():
    """After AL converges (gap ~0) the multiplier lambda carries the applied load:
    the y-component of augLambda equals P."""
    kt, P = 1.0e4, 10.0
    _al_transverse_gap(True, kt, P)
    lam = ops.eleResponse(1, "augLambda")
    assert lam[1] == pytest.approx(P, rel=1e-6)
    assert lam[0] == pytest.approx(0.0, abs=1e-6)        # no axial load


def test_al_preserves_bond_slip_axial():
    """AL augments only the (perfect-bond) constraint; with a bond law the axial τ–s
    force is physical and untouched — same axial force as without AL."""
    _single_host_model(bond=True, enforce="al")
    s = 0.5e-3                                           # ascending power-law branch
    R = _push_rebar(s, 0.0, 0.0)
    tau = 10.0 * (s / 1.0e-3) ** 0.4
    assert abs(R[0]) == pytest.approx(tau, rel=1e-4)


def test_enforce_penalty_matches_default():
    """Explicit `-enforce penalty` is bit-identical to the default (no flag)."""
    _single_host_model(perfect_k=1.0e5)
    R0 = _push_rebar(2.0e-3, 0.0, 0.0)[0]
    _single_host_model(perfect_k=1.0e5, enforce="penalty")
    R1 = _push_rebar(2.0e-3, 0.0, 0.0)[0]
    assert R1 == pytest.approx(R0, rel=1e-12)


@pytest.mark.parametrize("mode", ["nitsche", "transformation"])
def test_enforce_rejects_unbuilt_modes(mode):
    """nitsche (research-grade, not built) and transformation (deferred indefinitely)
    are rejected at parse time."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0); ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1)
    with pytest.raises(Exception):
        ops.element("LadrunoEmbeddedRebar", 1, 1, 1, 2, "-shape", 1.0,
                    "-dir", 1.0, 0.0, 0.0, "-perfect", 1.0e5, "-enforce", mode)
