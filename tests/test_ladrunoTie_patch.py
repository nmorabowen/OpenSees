"""ADR-62 — LadrunoTie kinematic mesh-tie generator (Zone-A).

LadrunoTie pairs each slave node to its nearest master facet (bucket-sort broad
phase + closest-point projection), evaluates the master shape weights N_i, and
emits one EQ_Constraint per slave node per translational DOF (u_s = Σ N_i u_{m,i}).
Any EQ-aware handler enforces it: Lagrange/Penalty (implicit/static, exact) or the
shipped LadrunoProjection (explicit). These tests cover both paths.

Mesh (shared): a NON-CONFORMING solid–solid interface at x=1.
  * master block: one stdBrick hex, x∈[0,1] (nodes 1-8).
  * slave  block: two stdBrick hexes stacked in z, x∈[1,2], z∈{0,0.5,1}
                  (nodes 11-22). Its x=1 face has 6 nodes; the z=0.5 nodes
                  (15,18) are NON-CONFORMING — they land on the MIDDLE of the
                  master facet z-edges (shape weights 0.5/0.5), so the tests
                  genuinely exercise the interpolation, not just corner ties.
The master interface facet is the quad-4 [2,3,7,6] = (1,0,0)(1,1,0)(1,1,1)(1,0,1).

We tie all 3 translational DOFs. The slave interface nodes are left otherwise
free (no SP on a tied/slave DOF — the projection handler refuses SP-on-slave);
their y,z are carried to 0 by the tie from the y,z-fixed master facet nodes.

  PATCH    a linear field u_x=γz (constant ε_xz, ν=0) must transmit EXACTLY across
           the tie: every Gauss point of both blocks reports the same constant
           stress = the analytic value. Fails unless the 0.5/0.5 weights are right.
  SPLIT    the non-conforming tied column reproduces the monolithic (conforming)
           column's stress and free cut-plane displacement under a uniform stretch.
  EXPLICIT on the shipped CentralDifferenceLadruno + LadrunoProjection: the tie
           holds to machine precision every step (ZERO penetration, vs a penalty
           SOFT tie's residual δ/h>0), and adding it does NOT change dt_cr.
  REFUSE   the named BLOCKER refusals (non-disjoint, duplicate, off-manifold,
           massless).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"

MASTER_FACET = [2, 3, 7, 6]
SLAVE_NODES = [11, 14, 15, 18, 19, 22]
# slave : list of (master_node, weight) the generator should emit
TIE_RELATIONS = {
    11: [(2, 1.0)], 14: [(3, 1.0)], 19: [(6, 1.0)], 22: [(7, 1.0)],
    15: [(2, 0.5), (6, 0.5)], 18: [(3, 0.5), (7, 0.5)],
}

_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_SLAVE = {11: (1, 0, 0.0), 12: (2, 0, 0.0), 13: (2, 1, 0.0), 14: (1, 1, 0.0),
          15: (1, 0, 0.5), 16: (2, 0, 0.5), 17: (2, 1, 0.5), 18: (1, 1, 0.5),
          19: (1, 0, 1.0), 20: (2, 0, 1.0), 21: (2, 1, 1.0), 22: (1, 1, 1.0)}
_ALL = {**_MASTER, **_SLAVE}
_MASTER_CUT = {2, 3, 6, 7}                 # master facet nodes (retained; x free)


def _mesh(E=1000.0, nu=0.0, rho=1.0, tie=True, ele="stdBrick"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.element(ele, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element(ele, 2, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element(ele, 3, 15, 16, 17, 18, 19, 20, 21, 22, 1)
    if tie:
        ops.LadrunoTie("-slaveNodes", len(SLAVE_NODES), *SLAVE_NODES,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def _gp_stresses(ele):
    ops.eleResponse(ele, "forces")           # force lazy-strain elements to evaluate
    s = list(ops.eleResponse(ele, "stresses"))
    return [s[i * 6:(i + 1) * 6] for i in range(len(s) // 6)]


def _static_solve(handler):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# --------------------------------------------------------------------------
# PATCH — constant-stress transmission across the non-conforming tie (static)
# --------------------------------------------------------------------------
def test_constant_stress_patch_nonconforming():
    """u_x = γ z (constant ε_xz, ν=0) gives the SAME constant stress at every Gauss
    point of both blocks. The z=0.5 slave nodes need the 0.5/0.5 facet weights."""
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, tie=True)

    # Prescribe the exact linear field u_x=γz on EVERY non-tied node (including the
    # master facet nodes — a shear field needs the full boundary driven, else a free
    # lateral edge cannot carry the shear). The 6 slave interface nodes stay free and
    # tied: the tie must deliver u_x=γz_s to them — the z=0.5 nodes via the 0.5/0.5
    # weights — or the slave-block stress comes out non-constant.
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_NODES:
            continue                          # tied: governed by the tie
        ops.fix(t, 0, 1, 1)                   # fix y,z; x free
        ops.sp(t, 1, gamma * float(z))        # impose u_x = γ z

    _static_solve("Lagrange")

    allgp = _gp_stresses(1) + _gp_stresses(2) + _gp_stresses(3)
    ref = allgp[0]
    atol = 1e-8 * E * gamma
    for v in allgp:                            # (a) constant across both blocks
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant stress: {v} vs {ref}"
    target = E * gamma / 2.0                    # (b) = analytic σ_xz (ν=0), uniaxial shear
    comps = [abs(c) for c in ref]
    assert abs(max(comps) - target) <= 1e-6 * target + atol, f"stress {ref} != {target}"
    assert sum(1 for c in comps if c > 1e-6 * target) == 1, f"not uniaxial shear: {ref}"


# --------------------------------------------------------------------------
# SPLIT — non-conforming tied column == monolithic conforming column
# --------------------------------------------------------------------------
def _monolithic_axial(E, U):
    """Conforming reference: a 2-hex column x∈[0,2], confined y/z, u_x=0 at x=0 and
    u_x=U at x=2. Returns (gp stress, u_x at the x=1 mid-plane node)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nodes, tag = {}, 1
    for xi in (0.0, 1.0, 2.0):
        for (y, z) in [(0, 0), (1, 0), (1, 1), (0, 1)]:
            nodes[tag] = (xi, y, z)
            ops.node(tag, float(xi), float(y), float(z))
            tag += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 1.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in nodes.items():
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Transformation")
    return _gp_stresses(1)[0], ops.nodeDisp(5, 1)   # node 5 = (1,0,0) mid-plane


def test_split_bar_equivalence():
    """A non-conforming split column tied with LadrunoTie reproduces the monolithic
    conforming column: GP stress + the FREE mid-plane (x=1) displacement."""
    E, U = 1000.0, 0.02
    ref_stress, ref_mid = _monolithic_axial(E, U)

    _mesh(E=E, nu=0.0, tie=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_NODES:
            continue
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)                 # master x=0 face
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)                   # slave x=2 face
    _static_solve("Lagrange")

    # master cut node 2 (x=1) is FREE — pulled to the linear-field value U/2 via the tie
    mid = ops.nodeDisp(2, 1)
    assert abs(mid - ref_mid) <= 1e-7 + 1e-6 * abs(ref_mid), f"mid-plane {mid} != {ref_mid}"
    got = _gp_stresses(1)[0]
    for gi, ri in zip(got, ref_stress):
        assert abs(gi - ri) <= 1e-6 * E * U + 1e-9, f"stress {got} != {ref_stress}"


# --------------------------------------------------------------------------
# EXPLICIT — zero penetration + dt_cr-neutral (the shipped projection path)
# --------------------------------------------------------------------------
def _build_explicit(tie, with_load):
    # SSPbrick: lumped (diagonal) mass — required by the projection handler on tied
    # DOFs (stdBrick's consistent mass is refused) and by the explicit `system Diagonal`.
    _mesh(E=1000.0, nu=0.0, rho=1.0, tie=tie, ele="SSPbrick")
    for t, (x, y, z) in _ALL.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)               # clamp the master base
        elif tie and t in SLAVE_NODES:
            continue                          # tied slave: all DOFs free (tie governs)
        else:
            ops.fix(t, 0, 1, 1)               # confine y,z; x free
    if with_load:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in (20, 21):                    # pull the slave free end (+x)
            ops.load(t, 5.0, 0.0, 0.0)
    ops.constraints(PROJ if tie else "Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")


def test_explicit_dt_cr_neutral():
    """Adding the kinematic tie does NOT change dt_cr (no penalty stiffness)."""
    dt0 = 1e-4
    _build_explicit(tie=True, with_load=False)
    assert ops.analyze(1, dt0) == 0
    dt_tied = ops.criticalTimeStep()
    _build_explicit(tie=False, with_load=False)
    assert ops.analyze(1, dt0) == 0
    dt_free = ops.criticalTimeStep()
    assert dt_tied > 0 and dt_free > 0, f"dt_cr not populated ({dt_tied}, {dt_free})"
    assert abs(dt_tied / dt_free - 1.0) < 1e-9, f"dt_cr not neutral: {dt_tied} vs {dt_free}"


def test_explicit_zero_penetration():
    """The tie holds to machine precision every step (gap = 0), and the load crosses
    the non-conforming interface (the slave end moves)."""
    dt0 = 1e-4
    _build_explicit(tie=True, with_load=True)
    assert ops.analyze(1, dt0) == 0
    dt = 0.4 * ops.criticalTimeStep()

    tie_err = 0.0
    for _ in range(200):
        assert ops.analyze(1, dt) == 0
        for s, rel in TIE_RELATIONS.items():
            us = ops.nodeDisp(s, 1)
            um = sum(w * ops.nodeDisp(m, 1) for (m, w) in rel)
            tie_err = max(tie_err, abs(us - um))
    assert tie_err < 1e-10, f"tie penetration {tie_err} (kinematic tie must be 0)"
    assert abs(ops.nodeDisp(20, 1)) > 1e-9, "slave end did not move — tie not load-carrying"

    # contrast: a penalty SOFT tie leaves a residual penetration
    #   δ/h = ε_iface·CFL²/(4·α_m·SOFSCL) > 0 (ADR-61 P0 closed form), whereas the
    # kinematic tie's gap is identically zero (asserted above).
    soft_pen = 1.0 * 0.9 ** 2 / (4.0 * 1.0 * 1.0)
    assert soft_pen > tie_err, "kinematic tie should beat the SOFT penalty residual"


# --------------------------------------------------------------------------
# REFUSE — the named BLOCKER refusals
# --------------------------------------------------------------------------
def _bare_master(E=1000.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 1.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def test_refuse_non_disjoint_surfaces():
    """BLOCKER-1: a node that is both slave and master is refused (MP-chain)."""
    _bare_master()
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 1, 2, "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_duplicate_slave():
    """BLOCKER-1: the same slave node listed twice is refused (double constraint)."""
    _bare_master()
    ops.node(99, 1.0, 0.5, 0.5)
    ops.mass(99, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 2, 99, 99, "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_off_manifold_slave():
    """OQ-3: a slave node off the master surface is refused (no IC snapping)."""
    _bare_master()
    ops.node(99, 1.5, 0.5, 0.5)               # 0.5 off the x=1 facet plane
    ops.mass(99, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_massless_slave():
    """BLOCKER-2: a tied slave DOF with no lumped mass is refused."""
    _bare_master()
    ops.node(99, 1.0, 0.5, 0.5)               # on the facet, but NO mass, NO element
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)
