"""ADR-78 — quadratic (serendipity) facets for the LadrunoTie integral mortar, Zone-A.

The driving case apeGmsh ADR 0086 measured: a hex20 face (quad8) tied onto a hex8
plate. The kernel keeps ALL geometry on the corner polygon and evaluates the full
serendipity basis at the mapped Gauss points (exact for flat facets with midsides
at edge midpoints — guarded, status -2).

Mesh — the quad8-on-quad4 patch (mirrors test_ladrunoTie_mortar.py's rig):
  * master block: ONE stdBrick hex, x∈[0,1] (nodes 1-8); its x=1 face is ONE quad4.
  * slave  block: ONE LadrunoBrick20 (-lumped = HRZ mass, the ADR-78 BLOCKER-1
    recipe — row-sum lumping would give zero/negative serendipity corner masses),
    x∈[1,2]; its x=1 face is ONE quad8 facet (4 corners + 4 midsides). The corner
    coords coincide with the master's but the MIDSIDES have no master counterpart —
    their tie rows are the genuinely quadratic transfer.

  PATCH   a linear field u_x=γz imposed on every non-tied node must transfer
          EXACTLY to the tied quad8 face (midsides included) with -dual — the
          linear-completeness gate the numpy oracle (proto_p2_2) proves at 2e-15.
  SPLIT   the quad8-mortar-tied split column matches the series-bar solution
          (mid-plane u = U/2, master stress = E·U/2).
  REFUSE  tri6 slave (named: the dual zero-rowsum degeneracy); curved quad8
          (midside off its edge midpoint, status -2); quadratic facets in the
          collocation mode; a quad8 slave off the master plane (the L1 gap guard).
  ACCEPT  tri6 MASTER facets (refused only on the slave side).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# master hex x∈[0,1]; x=1 face = quad facet [2,3,7,6]
_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
MASTER_FACET = [2, 3, 7, 6]

# slave hex20 x∈[1,2], standard 20-node serendipity ordering (corners 1-8, bottom
# edge mids 9-12, top edge mids 13-16, vertical mids 17-20), tags 101..120.
_S = {101: (1, 0, 0), 102: (2, 0, 0), 103: (2, 1, 0), 104: (1, 1, 0),
      105: (1, 0, 1), 106: (2, 0, 1), 107: (2, 1, 1), 108: (1, 1, 1),
      109: (1.5, 0, 0), 110: (2, 0.5, 0), 111: (1.5, 1, 0), 112: (1, 0.5, 0),
      113: (1.5, 0, 1), 114: (2, 0.5, 1), 115: (1.5, 1, 1), 116: (1, 0.5, 1),
      117: (1, 0, 0.5), 118: (2, 0, 0.5), 119: (2, 1, 0.5), 120: (1, 1, 0.5)}

# the slave x=1 face as a quad8 facet: corners-first, midsides in edge order
QUAD8_FACET = [101, 104, 108, 105, 112, 120, 116, 117]
SLAVE_IFACE = set(QUAD8_FACET)
_ALL = {**_MASTER, **_S}


def _mesh(E=1000.0, nu=0.0, rho=1.0, dual=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("LadrunoBrick20", 2, *range(101, 121), 1, "-lumped")
    args = ["-mortar", "-slaveFacets", 8, 1, *QUAD8_FACET,
            "-masterFacets", 4, 1, *MASTER_FACET]
    if dual:
        args.append("-dual")
    ops.LadrunoTie(*args)


def _static_solve(handler="Lagrange", numberer="Plain"):
    ops.constraints(handler)
    ops.numberer(numberer)
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _gp_stresses(ele):
    s = list(ops.eleResponse(ele, "stresses"))
    return [s[i * 6:(i + 1) * 6] for i in range(len(s) // 6)]


# --------------------------------------------------------------------------
# PATCH — linear field transfers exactly through the quad8 tie (midsides too)
# --------------------------------------------------------------------------
def test_quad8_patch_linear_exact_dual():
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, dual=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue                              # tied: governed by the mortar tie
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))            # impose u_x = γ z everywhere else
    _static_solve("Lagrange")

    # the headline gate: the tied quad8 face — CORNERS AND MIDSIDES — carries the
    # exact linear field (linear completeness of the quadratic dual transfer)
    for t in QUAD8_FACET:
        z = _ALL[t][2]
        got = ops.nodeDisp(t, 1)
        assert abs(got - gamma * z) <= 1e-9 * max(gamma, 1.0), \
            f"tied node {t} (z={z}): u_x={got} != {gamma * z}"
    # and the master block sees the constant stress state
    gps = _gp_stresses(1)
    ref = gps[0]
    target = E * gamma / 2.0
    atol = 1e-7 * E * gamma
    for v in gps:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant master stress: {v} vs {ref}"
    assert abs(max(abs(c) for c in ref) - target) <= 1e-5 * target + atol


def test_quad8_patch_linear_exact_standard():
    """Same patch with the STANDARD (dense) basis — both bases are linearly complete."""
    gamma = 1e-3
    _mesh(E=1000.0, nu=0.0, dual=False)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))
    _static_solve("Lagrange")
    for t in QUAD8_FACET:
        z = _ALL[t][2]
        assert abs(ops.nodeDisp(t, 1) - gamma * z) <= 1e-9


# --------------------------------------------------------------------------
# SPLIT — quad8-tied series column == the series-bar solution
# --------------------------------------------------------------------------
@pytest.mark.parametrize("numberer", ["Plain", "RCM"])
def test_quad8_split_bar_equivalence_dual(numberer):
    """Series-bar equivalence through the quad8 tie.

    The two numberers are the well-posedness gate, not decoration: this rig once
    fixed y/z on the TIED nodes as well, which the tie already ties to masters
    whose y/z are fixed — 16 linearly dependent Lagrange rows (EQ = SP_slave −
    Σ P_m·SP_master). u was still exact, only the multipliers were indeterminate,
    so LAPACK failed or not purely on pivot order: green on Plain/MKL, `U(i,i)=0`
    on CI. A second DOF ordering makes that class of over-constraint fail loudly
    instead of intermittently. See LEDGER_quirks (redundant SP on a tied DOF).
    """
    E, U = 1000.0, 0.02
    _mesh(E=E, nu=0.0, dual=True)
    x2_face = [t for t, (x, y, z) in _S.items() if abs(x - 2.0) < 1e-12]
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t not in SLAVE_IFACE:      # tied: y/z arrive through the mortar tie,
            ops.fix(t, 0, 1, 1)       # fixing them too is a redundant constraint
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif t in x2_face:
            ops.sp(t, 1, U)
    _static_solve("Lagrange", numberer)

    # equal E, A, L on both sides ⇒ the interface sits at U/2; the tied quad8 face
    # (master facet node 2 is free, pulled by the tie) must land there too
    assert abs(ops.nodeDisp(2, 1) - 0.5 * U) <= 1e-6 * U + 1e-9
    for t in QUAD8_FACET:
        assert abs(ops.nodeDisp(t, 1) - 0.5 * U) <= 1e-6 * U + 1e-9, \
            f"tied node {t}: {ops.nodeDisp(t, 1)} != {0.5 * U}"
    got = _gp_stresses(1)[0]
    assert abs(got[0] - E * U / 2.0) <= 1e-5 * E * U, f"master sigma_xx {got[0]}"


# --------------------------------------------------------------------------
# REFUSE / ACCEPT
# --------------------------------------------------------------------------
def _bare_master(E=1000.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 1.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def test_refuse_tri6_slave():
    """tri6 SLAVE facets: the corner shape functions integrate to zero ⇒ the dual
    scaling divides by zero — refused by name in BOTH bases (ADR-78 D2)."""
    _bare_master()
    tri6 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1),
            204: (1, 0.5, 0), 205: (1, 1, 0.5), 206: (1, 0.5, 0.5)}
    for t, (x, y, z) in tri6.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 6, 1, *sorted(tri6),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_curved_quad8():
    """A quad8 with a midside OFF its edge midpoint (curved face) is an INVALID
    quadratic facet (kernel status -2) — hard error, never a silent skip/mis-integration."""
    _bare_master()
    q8 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1), 204: (1, 0, 1),
          205: (1, 0.5, 0.1),   # lifted off the (201,202) midpoint => curved
          206: (1, 1, 0.5), 207: (1, 0.5, 1), 208: (1, 0, 0.5)}
    for t, (x, y, z) in q8.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 8, 1, *sorted(q8),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_quad8_off_manifold():
    """A (valid) quad8 slave off the master plane trips the L1 gap guard."""
    _bare_master()
    q8 = {201: (1.4, 0, 0), 202: (1.4, 1, 0), 203: (1.4, 1, 1), 204: (1.4, 0, 1),
          205: (1.4, 0.5, 0), 206: (1.4, 1, 0.5), 207: (1.4, 0.5, 1), 208: (1.4, 0, 0.5)}
    for t, (x, y, z) in q8.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 8, 1, *sorted(q8),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_collocation_quadratic_master():
    """Collocation stays linear-only: a quad8 -masterFacets without -mortar is refused
    (needs the quadratic closest-point projection — the message points at -mortar)."""
    _bare_master()
    q8 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1), 204: (1, 0, 1),
          205: (1, 0.5, 0), 206: (1, 1, 0.5), 207: (1, 0.5, 1), 208: (1, 0, 0.5)}
    for t, (x, y, z) in q8.items():
        ops.node(t, float(x), float(y), float(z))
    ops.node(300, 1.0, 0.4, 0.4)
    ops.mass(300, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-slaveNodes", 1, 300, "-masterFacets", 8, 1, *sorted(q8))


def test_refuse_accordion_gap_L1():
    """The ADR-78 unified L1 gap guard (OQ-2 behaviour change, quirks-ledger entry):
    an accordion slave — two facets tilted oppositely with the ridge ON the master
    plane — cancels the pre-ADR-78 SIGNED per-node gap signal at the shared ridge
    nodes, but the per-facet ∫|g_N|/area measure sees the full tilt ⇒ refuse.
    (Oracle T11 quantifies: signed ~1e-20, L1 = d/2.)"""
    _bare_master()
    d = 0.2                                # >> tolFrac·facetScale ⇒ must refuse
    acc = {201: (1 - d, 0, 0), 202: (1, 0.5, 0), 203: (1, 0.5, 1), 204: (1 - d, 0, 1),
           205: (1 + d, 1, 0), 206: (1 + d, 1, 1)}
    for t, (x, y, z) in acc.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar",
                       "-slaveFacets", 4, 2, 201, 202, 203, 204, 202, 205, 206, 203,
                       "-masterFacets", 4, 1, *MASTER_FACET,
                       "-outward", 1.0, 0.0, 0.0)


def test_refuse_duplicate_master_facets():
    """ADR-78 review fix (MAJOR-1): a doubly-listed master facet double-counts the
    per-facet coverage sum and could mask an uncovered slave region — the master
    self-overlap guard refuses coincident master-master overlap by name."""
    _bare_master()
    sn = {201: (1, 0.2, 0.2), 202: (1, 0.8, 0.2), 203: (1, 0.8, 0.8), 204: (1, 0.2, 0.8)}
    for t, (x, y, z) in sn.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 2, *MASTER_FACET, *MASTER_FACET)


# --------------------------------------------------------------------------
# MASTER-basis correctness (ADR-78 review MAJOR-3): a spring-rig transfer test.
# Master nodes are sp'd to a linear field; each slave node gets nodal mass and a
# soft grounding spring (so K is nonsingular); the tie then SLAVES the node to
# u_s = P u_m exactly (Lagrange), so u_s must reproduce the linear field — this
# fails loudly for any wrong quadratic master basis (e.g. permuted tri6 midsides,
# which preserve partition of unity and pass every internal generator gate).
# --------------------------------------------------------------------------
def _spring_rig_solve(master_nodes, master_facets, npsM, field):
    """master_nodes: {tag:(x,y,z)}; field(y,z) -> imposed u_x. Returns 4 slave tags."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in master_nodes.items():
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 0, 1, 1)
    sv = {501: (1, 0.2, 0.2), 502: (1, 0.8, 0.2), 503: (1, 0.8, 0.8), 504: (1, 0.2, 0.8)}
    ops.uniaxialMaterial("Elastic", 91, 10.0)
    for t, (x, y, z) in sv.items():
        ops.node(t, x, y, z)
        ops.mass(t, 1.0, 1.0, 1.0)
        g = t + 100                               # coincident fixed ground node
        ops.node(g, x, y, z)
        ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", t, g, t, "-mat", 91, 91, 91, "-dir", 1, 2, 3)
    ops.LadrunoTie("-mortar", "-dual", "-slaveFacets", 4, 1, 501, 502, 503, 504,
                   "-masterFacets", npsM, len(master_facets), *sum(master_facets, []),
                   "-outward", 1.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in master_nodes.items():
        ops.sp(t, 1, field(y, z))
    _static_solve("Lagrange")
    return sv


def test_tri6_master_linear_field_exact():
    mn = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 4: (1, 0, 1),
          5: (1, 0.5, 0), 6: (1, 1, 0.5), 7: (1, 0.5, 0.5),
          8: (1, 0.5, 1), 9: (1, 0, 0.5)}
    tri_a = [1, 2, 3, 5, 6, 7]
    tri_b = [1, 3, 4, 7, 8, 9]
    field = lambda y, z: 0.3 + 1.7 * y - 0.9 * z
    sv = _spring_rig_solve(mn, [tri_a, tri_b], 6, field)
    for t, (x, y, z) in sv.items():
        got = ops.nodeDisp(t, 1)
        assert abs(got - field(y, z)) <= 1e-9, \
            f"tri6-master transfer wrong at slave {t}: {got} != {field(y, z)}"


def test_quad8_master_linear_field_exact():
    mn = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 4: (1, 0, 1),
          5: (1, 0.5, 0), 6: (1, 1, 0.5), 7: (1, 0.5, 1), 8: (1, 0, 0.5)}
    q8 = [1, 2, 3, 4, 5, 6, 7, 8]
    field = lambda y, z: -0.2 + 0.8 * y + 1.3 * z
    sv = _spring_rig_solve(mn, [q8], 8, field)
    for t, (x, y, z) in sv.items():
        got = ops.nodeDisp(t, 1)
        assert abs(got - field(y, z)) <= 1e-9, \
            f"quad8-master transfer wrong at slave {t}: {got} != {field(y, z)}"


# --------------------------------------------------------------------------
# EXPLICIT — the quad8 tie on the shipped LadrunoProjection path (the handler the
# emitted message tells users to run; review MAJOR-3 noted no test exercised it).
# --------------------------------------------------------------------------
def test_quad8_explicit_projection():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 1.0)
    ops.element("SSPbrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("LadrunoBrick20", 2, *range(101, 121), 1, "-lumped")
    ops.LadrunoTie("-mortar", "-dual", "-slaveFacets", 8, 1, *QUAD8_FACET,
                   "-masterFacets", 4, 1, *MASTER_FACET)
    x2_face = [t for t, (x, y, z) in _S.items() if abs(x - 2.0) < 1e-12]
    for t, (x, y, z) in _ALL.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
        elif t in SLAVE_IFACE:
            continue
        else:
            ops.fix(t, 0, 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (102, 103):
        ops.load(t, 5.0, 0.0, 0.0)
    ops.constraints("LadrunoProjection")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    assert ops.analyze(1, 1e-4) == 0
    dt = 0.4 * ops.criticalTimeStep()
    assert dt > 0
    for _ in range(150):
        assert ops.analyze(1, dt) == 0
    assert abs(ops.nodeDisp(102, 1)) > 1e-9, "slave end did not move — quad8 tie inert"


def test_accept_tri6_master():
    """tri6 is refused only on the SLAVE side — a quad4 slave tied onto tri6 MASTER
    facets (two triangles covering the unit square) must be ACCEPTED."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # two tri6 masters on x=1 covering (y,z) ∈ [0,1]^2, corners-first ordering
    mn = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 4: (1, 0, 1),
          5: (1, 0.5, 0), 6: (1, 1, 0.5), 7: (1, 0.5, 0.5),
          8: (1, 0.5, 1), 9: (1, 0, 0.5)}
    for t, (x, y, z) in mn.items():
        ops.node(t, float(x), float(y), float(z))
    tri_a = [1, 2, 3, 5, 6, 7]      # (0,0)-(1,0)-(1,1) with edge mids
    tri_b = [1, 3, 4, 7, 8, 9]      # (0,0)-(1,1)-(0,1) with edge mids
    # slave: one quad4 strictly inside, nodal mass (no element needed)
    sv = {501: (1, 0.2, 0.2), 502: (1, 0.8, 0.2), 503: (1, 0.8, 0.8), 504: (1, 0.2, 0.8)}
    for t, (x, y, z) in sv.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    # must return without raising
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 501, 502, 503, 504,
                   "-masterFacets", 6, 2, *tri_a, *tri_b, "-outward", 1.0, 0.0, 0.0)
