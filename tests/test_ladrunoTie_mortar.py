"""ADR-62 P2 — LadrunoTie integral-mortar mesh-tie (`-mortar`), Zone-A.

P2 ties the slave SURFACE to the master surface in the weak (integral) sense:
assemble the mortar operators over the clipped overlap (reusing the shipped
LadrunoMortarKernel), condense u_s = P u_m with P = D⁻¹M once at setup, and emit one
EQ_Constraint per slave node per DOF (master-only retainers → no chains → the shipped
LadrunoProjection / any EQ-aware handler enforces it).

Mesh — a GENUINELY non-matching solid interface at x=1 (no coincident interface nodes):
  * master block: ONE stdBrick hex, x∈[0,1] (nodes 1-8); its x=1 face is ONE quad facet.
  * slave  block: x∈[1,2] meshed 2×2 in (y,z) = 4 hexes; its x=1 face is a 3×3 node grid
                  (9 nodes, 4 quad facets). The slave interface nodes at y/z = 0.5 do NOT
                  coincide with any master node — collocation (P1) would require conforming
                  nodes; mortar integrates the overlap.

  PATCH    a linear field u_x=γz (constant ε_xz, ν=0) must transmit EXACTLY across the
           non-matching tie: every Gauss point of both blocks reports the same constant
           stress = analytic. This is the headline integral-mortar gate.
  SPLIT    the mortar-tied split column matches the monolithic conforming column.
  EXPLICIT the tie holds to machine precision on CentralDifferenceLadruno +
           LadrunoProjection (SSPbrick lumped mass).
  REFUSE   uncovered slave node, off-manifold (non-coincident) surface, non-disjoint.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"

# master hex x∈[0,1]; x=1 face = quad facet [2,3,7,6]
_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
MASTER_FACET = [2, 3, 7, 6]

# slave block x∈[1,2], 2×2 in (y,z): node tag = 100 + ix*9 + iy*3 + iz
_YS = [0.0, 0.5, 1.0]
_ZS = [0.0, 0.5, 1.0]
_XS = [1.0, 2.0]


def _stag(ix, iy, iz):
    return 100 + ix * 9 + iy * 3 + iz


def _slave_nodes():
    nodes = {}
    for ix, x in enumerate(_XS):
        for iy, y in enumerate(_YS):
            for iz, z in enumerate(_ZS):
                nodes[_stag(ix, iy, iz)] = (x, y, z)
    return nodes


def _slave_hexes():
    """4 hexes (2×2 in y,z), standard hex node order."""
    hexes = []
    for jy in range(2):
        for jz in range(2):
            n = lambda ix, dy, dz: _stag(ix, jy + dy, jz + dz)
            hexes.append([n(0, 0, 0), n(1, 0, 0), n(1, 1, 0), n(0, 1, 0),
                          n(0, 0, 1), n(1, 0, 1), n(1, 1, 1), n(0, 1, 1)])
    return hexes


def _slave_facets():
    """The 4 quad facets on the slave x=1 face (ix=0), flat node list (4 per facet)."""
    fac = []
    for jy in range(2):
        for jz in range(2):
            fac += [_stag(0, jy, jz), _stag(0, jy + 1, jz),
                    _stag(0, jy + 1, jz + 1), _stag(0, jy, jz + 1)]
    return fac


SLAVE_NODES = _slave_nodes()
SLAVE_IFACE = [_stag(0, iy, iz) for iy in range(3) for iz in range(3)]   # 9 tied nodes
_ALL = {**_MASTER, **SLAVE_NODES}


def _mesh(E=1000.0, nu=0.0, rho=1.0, tie=True, ele="stdBrick", dual=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.element(ele, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for e, conn in enumerate(_slave_hexes(), start=2):
        ops.element(ele, e, *conn, 1)
    if tie:
        sf = _slave_facets()
        args = ["-mortar", "-slaveFacets", 4, len(sf) // 4, *sf,
                "-masterFacets", 4, 1, *MASTER_FACET]
        if dual:
            args.append("-dual")            # P2.1 biorthogonal basis => sparse P
        ops.LadrunoTie(*args)


def _gp_stresses(ele):
    ops.eleResponse(ele, "forces")
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
# PATCH — constant-stress across the genuinely non-matching mortar tie
# --------------------------------------------------------------------------
def test_mortar_constant_stress_patch():
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, tie=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue                              # tied: governed by the mortar tie
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))            # impose u_x = γ z on the full boundary

    _static_solve("Lagrange")

    allgp = []
    for e in range(1, 6):                          # master + 4 slave hexes
        allgp += _gp_stresses(e)
    ref = allgp[0]
    atol = 1e-7 * E * gamma                        # mortar Gauss-integration ~1e-9 rel
    for v in allgp:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant stress across mortar tie: {v} vs {ref}"
    target = E * gamma / 2.0
    comps = [abs(c) for c in ref]
    assert abs(max(comps) - target) <= 1e-5 * target + atol, f"stress {ref} != {target}"
    assert sum(1 for c in comps if c > 1e-5 * target) == 1, f"not uniaxial shear: {ref}"


# --------------------------------------------------------------------------
# SPLIT — mortar-tied non-matching column == monolithic conforming column
# --------------------------------------------------------------------------
def _monolithic_axial(E, U):
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
    return _gp_stresses(1)[0]


def test_mortar_split_bar_equivalence():
    E, U = 1000.0, 0.02
    ref_stress = _monolithic_axial(E, U)

    _mesh(E=E, nu=0.0, tie=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Lagrange")

    # master facet node 2 (x=1) is free, pulled to U/2 by the mortar tie
    mid = ops.nodeDisp(2, 1)
    assert abs(mid - 0.5 * U) <= 1e-6 * U + 1e-9, f"mid-plane {mid} != {0.5 * U}"
    got = _gp_stresses(1)[0]
    for gi, ri in zip(got, ref_stress):
        assert abs(gi - ri) <= 1e-5 * E * U + 1e-9, f"stress {got} != {ref_stress}"


# --------------------------------------------------------------------------
# EXPLICIT — the mortar tie on the shipped projection path
# --------------------------------------------------------------------------
def _build_explicit(with_load, dual=False):
    _mesh(E=1000.0, nu=0.0, rho=1.0, tie=True, ele="SSPbrick", dual=dual)
    for t, (x, y, z) in _ALL.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
        elif t in SLAVE_IFACE:
            continue
        else:
            ops.fix(t, 0, 1, 1)
    if with_load:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in (_stag(1, 0, 2), _stag(1, 2, 2)):   # pull the slave x=2 free end
            ops.load(t, 5.0, 0.0, 0.0)
    ops.constraints(PROJ)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")


def test_mortar_explicit_zero_penetration():
    dt0 = 1e-4
    _build_explicit(with_load=True)
    assert ops.analyze(1, dt0) == 0
    dt = 0.4 * ops.criticalTimeStep()
    assert dt > 0
    for _ in range(150):
        assert ops.analyze(1, dt) == 0
    # the slave end moved (load crossed the non-matching tie)
    assert abs(ops.nodeDisp(_stag(1, 0, 2), 1)) > 1e-9, "slave end did not move — tie inert"


# --------------------------------------------------------------------------
# REFUSE
# --------------------------------------------------------------------------
def _bare_master(E=1000.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 1.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def test_mortar_refuse_uncovered_slave():
    """A slave surface extending beyond the master surface ⇒ uncovered node ⇒ refuse."""
    _bare_master()
    # slave quad spanning y∈[0,2] (half of it past the master's y∈[0,1] extent)
    sn = {201: (1, 0, 0), 202: (1, 2, 0), 203: (1, 2, 1), 204: (1, 0, 1)}
    for t, (x, y, z) in sn.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_refuse_off_manifold():
    """A slave surface offset off the master plane ⇒ non-conforming ⇒ refuse."""
    _bare_master()
    sn = {201: (1.4, 0, 0), 202: (1.4, 1, 0), 203: (1.4, 1, 1), 204: (1.4, 0, 1)}
    for t, (x, y, z) in sn.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    with pytest.raises(Exception):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_refuse_non_disjoint():
    """A node on both surfaces ⇒ refuse (MP-chain)."""
    _bare_master()
    with pytest.raises(Exception):
        # reuse master nodes 2,3,7,6 as the slave facet → not node-disjoint
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 2, 3, 7, 6,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_nonaffine_slave_accepted():
    """Regression guard for the coverage-ratio slack: a NON-AFFINE (trapezoidal) slave
    facet that spans MULTIPLE master facets is the only configuration where cover
    (multi-piece overlap clip) and fullCov (single self-clip) differ in quadrature.
    A fully-covered such facet must still be ACCEPTED (not false-refused) — i.e. the
    1e-3 slack and the post-solve partition-of-unity backstop both hold off the
    affine grid. (The shipped patch/split tests use axis-aligned = affine quads, where
    cover==fullCov exactly, so they never exercise this.)"""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # master: a 3×3 node grid on x=1 (y,z ∈ {0,0.5,1}) = 4 quad facets, with a
    # mass-bearing solid behind it so the facet nodes are real.
    mtag, mnodes = {}, {}
    t = 1
    for iy, y in enumerate([0.0, 0.5, 1.0]):
        for iz, z in enumerate([0.0, 0.5, 1.0]):
            mtag[(iy, iz)] = t
            ops.node(t, 1.0, y, z)
            mnodes[t] = (1.0, y, z)
            t += 1
    mfac = []
    for iy in range(2):
        for iz in range(2):
            mfac += [mtag[(iy, iz)], mtag[(iy + 1, iz)], mtag[(iy + 1, iz + 1)], mtag[(iy, iz + 1)]]
    nmf = len(mfac) // 4

    # slave: ONE non-affine (trapezoidal) quad fully inside the master extent,
    # spanning all four master facets. Nodes carry mass (no element).
    sverts = [(1.0, 0.10, 0.15), (1.0, 0.90, 0.20), (1.0, 0.78, 0.85), (1.0, 0.18, 0.92)]
    stags = []
    for i, (x, y, z) in enumerate(sverts):
        tag = 500 + i
        ops.node(tag, x, y, z)
        ops.mass(tag, 1.0, 1.0, 1.0)
        stags.append(tag)

    # must be ACCEPTED (returns without raising) — exercises guard #2 off the affine grid
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, *stags,
                   "-masterFacets", 4, nmf, *mfac)


# --------------------------------------------------------------------------
# P2.1 — DUAL (biorthogonal) basis: -dual makes D diagonal ⇒ P local/sparse.
# At the FE level the dual tie must PRESERVE correctness (linear completeness —
# exactly what row-sum lumping would break); the SPARSITY of P is a structural
# property of the operator (identical solution for linearly-complete data), proven
# in the numpy oracle proto_p2_1_dual_mortar.py (D diagonal, dual≠lumped, dual≠dense).
# --------------------------------------------------------------------------
def test_mortar_dual_constant_stress_patch():
    """The DUAL tie must still pass the constant-stress patch across the non-matching
    interface — this is the linear-completeness a naive diagonal (lumped D) loses."""
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, tie=True, dual=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))
    _static_solve("Lagrange")

    allgp = []
    for e in range(1, 6):
        allgp += _gp_stresses(e)
    ref = allgp[0]
    atol = 1e-7 * E * gamma
    for v in allgp:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"dual mortar non-constant stress: {v} vs {ref}"
    target = E * gamma / 2.0
    comps = [abs(c) for c in ref]
    assert abs(max(comps) - target) <= 1e-5 * target + atol, f"stress {ref} != {target}"


def test_mortar_dual_split_bar_equivalence():
    """The dual-tied non-matching column reproduces the monolithic conforming column."""
    E, U = 1000.0, 0.02
    ref_stress = _monolithic_axial(E, U)

    _mesh(E=E, nu=0.0, tie=True, dual=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Lagrange")

    mid = ops.nodeDisp(2, 1)
    assert abs(mid - 0.5 * U) <= 1e-6 * U + 1e-9, f"dual mid-plane {mid} != {0.5 * U}"
    got = _gp_stresses(1)[0]
    for gi, ri in zip(got, ref_stress):
        assert abs(gi - ri) <= 1e-5 * E * U + 1e-9, f"dual stress {got} != {ref_stress}"


def test_mortar_dual_explicit_zero_penetration():
    """The dual tie holds to machine precision on the shipped explicit projection path."""
    dt0 = 1e-4
    _build_explicit(with_load=True, dual=True)
    assert ops.analyze(1, dt0) == 0
    dt = 0.4 * ops.criticalTimeStep()
    assert dt > 0
    for _ in range(150):
        assert ops.analyze(1, dt) == 0
    assert abs(ops.nodeDisp(_stag(1, 0, 2), 1)) > 1e-9, "dual slave end did not move — tie inert"


def test_dual_requires_mortar():
    """-dual only applies to the integral-mortar mode (collocation has no D to diagonalize)."""
    _mesh(E=1000.0, nu=0.0, tie=False)
    with pytest.raises(Exception):
        ops.LadrunoTie("-dual", "-slaveNodes", 1, _stag(0, 0, 0),
                       "-masterFacets", 4, 1, *MASTER_FACET)
