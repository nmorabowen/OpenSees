"""LadrunoBrick20 ``-formulation uri`` (2x2x2 reduced integration, the C3D20R
analog) — Zone-A P2 battery. ADR 72 §6 P2 row; signed specs:
Ladruno_implementation/_adr72_p2_test_specs.md (S0-S9).

THE P2 ANCHOR: upstream has no reduced-integration 20-node element, so the
reduce-to anchor does not exist for uri. Correctness is anchored by the P0
sympy oracle (tests/hex20_reference.py): numpy-assembled K / f / sigma at the
GP8 table (z-fastest lexicographic, the Brick order) vs the element, on the
unit cube AND the P1 distorted hex (S0). The contract tests (S1-S5) then pin
the ADR §3.2 spurious-mode semantics:

  S1  census: a lone free uri element has EXACTLY 12 zero-energy modes
      (6 RBM + 6 spurious), the spurious subspace matching the oracle's;
  S2  non-communicability: a determinately-restrained 2x2x2 block has ZERO
      spurious global modes — while a lone element under the same 3-2-1
      restraint is outright unstable (Bathe Fig. 5.46);
  S3  the single-stack pathology DEMONSTRATED (1x1x4 cantilever stack carries
      a near-mechanism orders softer than std's softest physical mode) —
      the honest counter-example, kept as living documentation;
  S4  coarse bending (Onate Fig. 8.23 class): 1-layer uri cantilever >= 0.98
      of beam theory where the H8 LadrunoBrick std locks;
  S5  nu=0.4999: uri relieves the volumetric lock, std locks (escalation
      trigger if it does not — ADR §3.4);
  S6  Barlow: uri GP stresses beat std GP stresses vs the analytic bending
      field (each set at its OWN stations);
  S7  recorder seam (debt a): a uri element self-describes via basisInfo and
      the LadrunoRecorder writes 8 GP stations (GP_PARAM == GP8, GLOBAL
      coords == the isoparametric map);
  S8  surface: uri/reduced accepted, JSON says uri, DB round-trip preserves
      ordinal 1, -hourglass still refused, advisory probe alive under uri;
  S9  cost: uri/std tangent-assembly time ratio reported (loose CI-safe gate).

S3 and S5 carry ESCALATION semantics (plan task 2.4): if the stack shows no
soft mode or std does not lock, that contradicts ADR §3.2/§3.4 — the test
failing IS the alarm; do not weaken the assertion to make it pass.
"""
import glob
import math
import os
import time

import pytest

from _testbed import ops
from _testbed.fem_checks import zero_energy_mode_count, n_rigid_body_modes
from _testbed.roundtrip import database_roundtrip

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")           # the oracle derives symbolically
import hex20_reference as ref          # noqa: E402  (P0 oracle — debt c: no re-transcription)

pytestmark = [pytest.mark.zone_a]

E, NU = 1000.0, 0.3
RHO = 2.0e-3

# --------------------------------------------------------------------------
# fixtures — shared with the P1 battery (same distorted straight-edged hex)
# --------------------------------------------------------------------------
_CORNERS = [
    (0.00, 0.00, 0.00), (1.00, 0.10, 0.00), (1.10, 1.00, 0.10), (0.05, 0.95, 0.00),
    (0.00, 0.05, 1.00), (1.00, 0.00, 1.05), (1.05, 1.00, 1.10), (0.00, 1.00, 0.95),
]
_CUBE_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_CONN = list(range(1, 21))
_BOTTOM = [1, 2, 3, 4, 9, 10, 11, 12]
_LOADED = {
    5: (10.0, 0.0, 0.0), 6: (0.0, 8.0, 0.0), 7: (0.0, 0.0, -12.0),
    8: (5.0, 5.0, 5.0), 13: (-3.0, 2.0, 0.0), 15: (0.0, -4.0, 6.0),
    18: (1.0, 1.0, 1.0), 20: (0.0, 0.0, -5.0),
}

_URI = ["-formulation", "uri"]
_STD = ["-formulation", "std"]


def _hex20_nodes(corners):
    """{tag: coords} straight-edged hex20 from 8 corners (brcshl order) —
    edge midpoints from the oracle's EDGES table (debt c)."""
    nodes = {i + 1: tuple(corners[i]) for i in range(8)}
    for k, (a, b) in enumerate(ref.EDGES):        # 0-based corner pairs
        ca, cb = corners[a], corners[b]
        nodes[9 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    return nodes


def _build_single(corners, extra, rho=0.0, e_mod=E, nu=NU):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nodes = _hex20_nodes(corners)
    for tag, c in nodes.items():
        ops.node(tag, *c)
    if rho:
        ops.nDMaterial("ElasticIsotropic", 1, e_mod, nu, rho)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, e_mod, nu)
    ops.element("LadrunoBrick20", 1, *_CONN, 1, *extra)
    return nodes


def _static_rig():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


def _dof_map(node_tags):
    """{(node, dof): equation number} AFTER an analysis is set up (numberer
    has run) — lets the printA matrix be compared entry-for-entry against an
    oracle matrix in element DOF order, with no numbering assumptions."""
    m = {}
    for n in node_tags:
        eqns = ops.nodeDOFs(n)
        for d, e in enumerate(eqns):
            if e >= 0:
                m[(n, d)] = e
    return m


def _oracle_reduced(K_or, dof_map, node_tags):
    """Reduce the 60x60 oracle K to the free-DOF printA layout."""
    nfree = len(dof_map)
    K = np.zeros((nfree, nfree))
    idx = {}
    for a, n in enumerate(node_tags):
        for d in range(3):
            if (n, d) in dof_map:
                idx[dof_map[(n, d)]] = 3 * a + d
    for i in range(nfree):
        for j in range(nfree):
            K[i, j] = K_or[idx[i], idx[j]]
    return K


def _mesh_eigs(all_nodes, n_free):
    """Eigenvalues of the assembled free-DOF K (unit masses => K's spectrum),
    ascending absolute values."""
    for n in all_nodes:
        ops.mass(n, 1.0, 1.0, 1.0)
    eigs = ops.eigen("-fullGenLapack", n_free)
    return sorted(abs(e) for e in eigs)


# --------------------------------------------------------------------------
# hex20 block mesher (dedup by coordinate) — used by S2/S3/S4/S5/S6/S9
# --------------------------------------------------------------------------
def _hex20_block(ele_name, mat_tag, extra, nx, ny, nz, sx, sy, sz, tag0=1):
    """Build an nx x ny x nz block of straight-edged hex20 (or hex8 for
    LadrunoBrick) elements. Returns (conns, coords_of) — conns[e] is the
    element's node-tag list, coords_of maps tag -> coords."""
    tag_of, coords_of = {}, {}
    next_tag = [1]

    def add_node(c):
        key = tuple(round(v, 10) for v in c)
        if key not in tag_of:
            ops.node(next_tag[0], *c)
            tag_of[key] = next_tag[0]
            coords_of[next_tag[0]] = c
            next_tag[0] += 1
        return tag_of[key]

    hex8 = "Brick" in ele_name and "20" not in ele_name
    conns = []
    etag = tag0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x0, x1 = i * sx, (i + 1) * sx
                y0, y1 = j * sy, (j + 1) * sy
                z0, z1 = k * sz, (k + 1) * sz
                corners = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                           (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
                if hex8:
                    conn = [add_node(c) for c in corners]
                else:
                    n20 = _hex20_nodes(corners)
                    conn = [add_node(n20[t]) for t in range(1, 21)]
                ops.element(ele_name, etag, *conn, mat_tag, *extra)
                conns.append(conn)
                etag += 1
    return conns, coords_of


# ==========================================================================
# S0 — oracle anchor: uri K / f / per-GP sigma vs the P0 sympy oracle
# ==========================================================================
@pytest.mark.parametrize("corners", [_CUBE_CORNERS, _CORNERS],
                         ids=["cube", "distorted"])
def test_S0_uri_stiffness_matches_oracle(corners):
    """Assembled free-DOF K (printA) == oracle K at GP8, rel-max <= 1e-12."""
    nodes = _build_single(corners, _URI)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    _static_rig()
    assert ops.analyze(1) == 0
    K_ele = np.asarray(ops.printA("-ret"), dtype=float)

    X = ref.straight_edge_nodes(corners)
    K_or = ref.stiffness(X, ref.iso_D(E, NU), ref.GP8_F)
    dof_map = _dof_map(_CONN)
    K_red = _oracle_reduced(K_or, dof_map, _CONN)
    assert K_ele.size == K_red.size
    diff = np.max(np.abs(K_ele.reshape(K_red.shape) - K_red))
    assert diff <= 1e-12 * np.max(np.abs(K_red)), (
        f"uri K vs oracle: rel-max {diff / np.max(np.abs(K_red)):.3e} > 1e-12"
    )


def test_S0_uri_force_and_stress_match_oracle():
    """Resisting force == K_oracle @ u and the 8 per-GP stresses ==
    D B(GP8_L) u, station-by-station (pins materialPointers[L] <-> GP8[L])."""
    nodes = _build_single(_CORNERS, _URI)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADED.items():
        ops.load(n, fx, fy, fz)
    _static_rig()
    assert ops.analyze(1) == 0

    u = np.zeros(60)
    for a, n in enumerate(_CONN):
        d = ops.nodeDisp(n)
        u[3 * a: 3 * a + 3] = d

    X = ref.straight_edge_nodes(_CORNERS)
    D = ref.iso_D(E, NU)
    K_or = ref.stiffness(X, D, ref.GP8_F)
    f_or = K_or @ u
    forces = np.asarray(ops.eleResponse(1, "forces"), dtype=float)
    scale = np.max(np.abs(f_or))
    assert np.max(np.abs(forces - f_or)) <= 1e-12 * scale, (
        f"uri resisting force vs oracle: "
        f"{np.max(np.abs(forces - f_or)) / scale:.3e} > 1e-12"
    )

    sig = np.asarray(ops.eleResponse(1, "stresses"), dtype=float)
    assert sig.size == 8 * 6, f"uri must expose 8 GP stations, got {sig.size // 6}"
    smax = np.max(np.abs(sig))
    for L in range(8):
        _, dN = ref.shape_numeric(tuple(ref.GP8_F[L, :3]))
        dNdx, _ = ref.cart_grad(X, dN)
        s_or = D @ (ref.fill_B(dNdx) @ u)
        assert np.max(np.abs(sig[6 * L: 6 * L + 6] - s_or)) <= 1e-12 * smax, (
            f"uri GP {L}: stress != oracle (pairing/order broken?)"
        )


def _build_any(ele_name, extra, rho=0.0):
    """Single distorted hex20 of any 20-node element type."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, c in _hex20_nodes(_CORNERS).items():
        ops.node(tag, *c)
    if rho:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU, rho)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element(ele_name, 1, *_CONN, 1, *extra)


def _transient_tangent(ele_name, extra, dt=1.0e-2):
    _build_any(ele_name, extra, rho=RHO)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, dt) == 0
    return np.asarray(ops.printA("-ret"), dtype=float)


def test_S0_uri_mass_equals_std_and_upstream():
    """Newmark tangent K+c3*M: (uri - K_uri) == (std - K_std) == upstream —
    i.e. the consistent mass is 27-pt and formulation-INDEPENDENT by design.
    Compared via (K+c3M) - K printA pairs so only M survives."""
    def K_static(ele, extra):
        _build_any(ele, extra)
        for n in _BOTTOM:
            ops.fix(n, 1, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(7, 0.0, 0.0, -1.0)
        _static_rig()
        assert ops.analyze(1) == 0
        return np.asarray(ops.printA("-ret"), dtype=float)

    M = {}
    for name, (ele, extra) in {
        "uri": ("LadrunoBrick20", _URI),
        "std": ("LadrunoBrick20", _STD),
        "upstream": ("20NodeBrick", []),
    }.items():
        M[name] = _transient_tangent(ele, extra) - K_static(ele, extra)

    # M is extracted as (K+c3M) - K from two printA assemblies, so the floor
    # is the CANCELLATION noise |K|*eps, not machine eps on M itself — the
    # in-element M is bitwise formulation-independent (same 27-pt code path);
    # 1e-12 on the c3M scale is the honest printA-level pin.
    scale = np.max(np.abs(M["std"]))
    assert np.max(np.abs(M["uri"] - M["std"])) <= 1e-12 * scale, (
        "uri consistent M must equal std M (both 27-pt)"
    )
    assert np.max(np.abs(M["uri"] - M["upstream"])) <= 1e-12 * scale, (
        "uri consistent M must equal upstream Twenty_Node_Brick M"
    )


# ==========================================================================
# S1 — mode census: exactly 12 zero modes (6 RBM + 6 spurious), catalogued
# ==========================================================================
def _rbm_basis(coords):
    """60x6 rigid-body basis (3 translations + 3 centroid rotations)."""
    X = np.asarray(coords, dtype=float)
    c = X.mean(axis=0)
    R = np.zeros((60, 6))
    for a in range(20):
        x, y, z = X[a] - c
        R[3 * a + 0, 0] = 1.0
        R[3 * a + 1, 1] = 1.0
        R[3 * a + 2, 2] = 1.0
        R[3 * a + 1, 3], R[3 * a + 2, 3] = -z, y      # rot x
        R[3 * a + 0, 4], R[3 * a + 2, 4] = z, -x      # rot y
        R[3 * a + 0, 5], R[3 * a + 1, 5] = -y, x      # rot z
    return R


def _spurious_complement(V, R):
    """Orthonormal basis of the part of span(V) orthogonal to span(R) —
    the 6 spurious directions of a 12-dim zero-energy space containing the
    6 rigid-body modes."""
    QR, _ = np.linalg.qr(R)
    P = V - QR @ (QR.T @ V)
    U, sv, _ = np.linalg.svd(P, full_matrices=False)
    assert sv[5] > 1e-3, f"spurious complement is rank-deficient: sv={sv[:8]}"
    return U[:, :6]


@pytest.mark.parametrize("corners", [_CUBE_CORNERS, _CORNERS],
                         ids=["cube", "distorted"])
def test_S1_uri_census_exactly_12_zero_modes(corners):
    nodes = _build_single(corners, _URI)
    count, eigs = zero_energy_mode_count(_CONN, ndf=3, tol_rel=1e-8)
    assert count == 12, (
        f"uri census: {count} zero modes (want 12 = 6 RBM + 6 spurious); "
        f"eigs[:14]={sorted(abs(e) for e in eigs)[:14]}"
    )
    es = sorted(abs(e) for e in eigs)
    sep = es[12] / max(es[11], 1e-300)
    assert sep >= 1e10, f"zero/physical separation {sep:.3e} < 1e10"
    print(f"[S1:{'cube' if corners is _CUBE_CORNERS else 'distorted'}] "
          f"12 zero modes, separation {sep:.3e}, eigs[:13]={es[:13]}")


def test_S1_std_control_exactly_6():
    _build_single(_CUBE_CORNERS, _STD)
    count, _ = zero_energy_mode_count(_CONN, ndf=3, tol_rel=1e-8)
    assert count == n_rigid_body_modes(3)


def test_S1_spurious_subspace_matches_oracle():
    """The 6 non-RBM zero modes span the SAME subspace as the P0 oracle's
    (recomputed via hex20_reference, not a checked-in dump): the projection
    residual sin(theta_max) < 1e-7 for every principal direction."""
    nodes = _build_single(_CUBE_CORNERS, _URI)
    count, eigs = zero_energy_mode_count(_CONN, ndf=3, tol_rel=1e-8)
    assert count == 12
    V = np.zeros((60, 12))
    for m in range(12):
        for a, n in enumerate(_CONN):
            V[3 * a: 3 * a + 3, m] = ops.nodeEigenvector(n, m + 1)

    coords = [nodes[t] for t in _CONN]
    R = _rbm_basis(coords)
    Q1 = _spurious_complement(V, R)

    X = ref.straight_edge_nodes(_CUBE_CORNERS)
    K_or = ref.stiffness(X, ref.iso_D(E, NU), ref.GP8_F)
    w, vec = np.linalg.eigh(K_or)
    V_or = vec[:, :12]                     # 12 smallest of the oracle K8
    assert w[12] / max(abs(w[11]), 1e-300) > 1e8, "oracle census sanity"
    Q2 = _spurious_complement(V_or, R)

    # sin of the largest principal angle between the two 6-dim subspaces
    s = np.linalg.svd(Q1.T @ Q2, compute_uv=False)
    sin_max = math.sqrt(max(0.0, 1.0 - min(s) ** 2))
    assert sin_max < 1e-7, (
        f"element spurious subspace != oracle's: sin(theta_max) = {sin_max:.3e}"
    )
    print(f"[S1] spurious-subspace match: sin(theta_max) = {sin_max:.3e}")


# ==========================================================================
# S2 — non-communicability: determinate 2x2x2 block clean; lone element dies
# ==========================================================================
def _fix_321(tag_A, tag_B, tag_C):
    """Statically determinate 3-2-1 restraint (exactly 6 DOFs)."""
    ops.fix(tag_A, 1, 1, 1)
    ops.fix(tag_B, 0, 1, 1)
    ops.fix(tag_C, 0, 0, 1)


def _tag_at(coords_of, xyz, tol=1e-9):
    for t, c in coords_of.items():
        if all(abs(c[d] - xyz[d]) < tol for d in range(3)):
            return t
    raise AssertionError(f"no node at {xyz}")


@pytest.mark.parametrize("extra", [_URI, _STD], ids=["uri", "std"])
def test_S2_block_2x2x2_zero_spurious_modes(extra, request):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    conns, coords = _hex20_block("LadrunoBrick20", 1, extra, 2, 2, 2, 1, 1, 1)
    _fix_321(_tag_at(coords, (0, 0, 0)), _tag_at(coords, (2, 0, 0)),
             _tag_at(coords, (0, 2, 0)))
    n_free = 3 * len(coords) - 6
    es = _mesh_eigs(list(coords), n_free)
    lam_min, lam_max = es[0], es[-1]
    n_zero = sum(1 for e in es if e < 1e-10 * lam_max)
    assert n_zero == 0, (
        f"{request.node.callspec.id}: {n_zero} near-zero global modes in the "
        f"determinate 2x2x2 block (non-communicability violated); "
        f"lam_min/lam_max = {lam_min / lam_max:.3e}"
    )
    # stash for the cross-formulation ratio check
    _S2_LAMBDA_MIN[request.node.callspec.id] = lam_min
    print(f"[S2:{request.node.callspec.id}] lam_min={lam_min:.6e} "
          f"lam_min/lam_max={lam_min / lam_max:.3e}")


_S2_LAMBDA_MIN = {}


def test_S2_uri_softest_mode_is_physical():
    """Runs after the two parametrized cases: the softest uri block mode is a
    physical mode, not a near-mechanism."""
    assert "uri" in _S2_LAMBDA_MIN and "std" in _S2_LAMBDA_MIN, (
        "parametrized S2 cases must run first (file order)"
    )
    ratio = _S2_LAMBDA_MIN["uri"] / _S2_LAMBDA_MIN["std"]
    assert ratio >= 0.05, f"lam_min(uri)/lam_min(std) = {ratio:.3e} < 0.05"
    print(f"[S2] lam_min(uri)/lam_min(std) = {ratio:.4f}")


def test_S2_lone_element_determinate_restraint_contrast():
    """Bathe Fig. 5.46: one uri element under the SAME 3-2-1 determinate
    restraint is outright unstable (>=1 zero mode; expect 6); std is clean."""
    for extra, expect_unstable in ((_URI, True), (_STD, False)):
        nodes = _build_single(_CUBE_CORNERS, extra)
        _fix_321(1, 2, 4)          # corners (0,0,0), (1,0,0), (0,1,0)
        es = _mesh_eigs(_CONN, 3 * 20 - 6)
        n_zero = sum(1 for e in es if e < 1e-8 * es[-1])
        if expect_unstable:
            assert n_zero >= 1, (
                "lone uri element under determinate restraint must be "
                f"unstable; census = {n_zero}"
            )
            print(f"[S2] lone uri element: {n_zero} zero modes (unstable, as "
                  "documented)")
        else:
            assert n_zero == 0, f"std lone element must be stable, got {n_zero}"


# ==========================================================================
# S3 — single-stack pathology (ESCALATION trigger; do not weaken)
# ==========================================================================
def _stack_ratio(extra):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    conns, coords = _hex20_block("LadrunoBrick20", 1, extra, 1, 1, 4, 1, 1, 1)
    fixed = [t for t, c in coords.items() if abs(c[2]) < 1e-9]
    assert len(fixed) == 8            # the z=0 face of a 1x1 hex20 column
    for t in fixed:
        ops.fix(t, 1, 1, 1)
    n_free = 3 * (len(coords) - len(fixed))
    es = _mesh_eigs(list(coords), n_free)
    return es[0] / es[-1], es, coords


def test_S3_single_stack_pathology_demonstrated():
    r_uri, es_uri, coords = _stack_ratio(_URI)
    # capture the offending mode BEFORE the next model wipes the domain
    tip = [t for t, c in coords.items() if abs(c[2] - 4.0) < 1e-9]
    tip_pat = {t: ops.nodeEigenvector(t, 1) for t in tip}
    r_std, es_std, _ = _stack_ratio(_STD)
    assert r_uri < 1e-3 * r_std, (
        f"ESCALATION (plan 2.4): the 1x1x4 uri stack shows NO near-mechanism "
        f"(lam_min/lam_max: uri {r_uri:.3e} vs std {r_std:.3e}) — this "
        f"contradicts ADR 72 §3.2's propagation claim; STOP and adjudicate."
    )
    # living documentation: the tip-face pattern of the offending mode
    print(f"[S3] pathology ratio (uri/std normalized lam_min) = "
          f"{r_uri / r_std:.3e}; uri lam_min/lam_max = {r_uri:.3e}")
    for t, v in sorted(tip_pat.items()):
        print(f"[S3]   tip node {t}: u = ({v[0]:+.3e}, {v[1]:+.3e}, {v[2]:+.3e})")


# ==========================================================================
# S4 / S5 / S6 — cantilever family (coarse bending, nu-relief, Barlow)
# ==========================================================================
_L, _B, _H = 10.0, 1.0, 1.0


def _tip_face_loads(conns, coords, load_z):
    """Consistent nodal forces of a uniform tip-face traction (total load_z)
    on the LAST element's x=L face, from the oracle basis at xi=+1."""
    conn = conns[-1]
    X = np.asarray([coords[t] for t in conn], dtype=float)
    face = [a for a in range(20) if abs(X[a, 0] - _L) < 1e-9]
    assert len(face) == 8
    tau = load_z / (_B * _H)
    # 3x3 GL on (eta, zeta); face jacobian for the straight rectangular face
    g = math.sqrt(3.0 / 5.0)
    pts = [-g, 0.0, g]
    wts = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    dy_deta = 0.5 * _B / 1.0   # element spans the full cross-section
    dz_dzeta = 0.5 * _H / 1.0
    F = {}
    for a in face:
        s = 0.0
        for ge, we in zip(pts, wts):
            for gz, wz in zip(pts, wts):
                N, _ = ref.shape_numeric((1.0, ge, gz))
                s += we * wz * N[a] * dy_deta * dz_dzeta
        F[conn[a]] = tau * s
    total = sum(F.values())
    assert abs(total - load_z) < 1e-9 * abs(load_z)
    return F


def _cantilever(ele_name, extra, n_ele, e_mod, nu, P=-1.0):
    """Meshed 1-layer cantilever, tip shear P (z); returns (mean tip u_z,
    conns, coords)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, e_mod, nu)
    conns, coords = _hex20_block(ele_name, 1, extra, n_ele, 1, 1,
                                 _L / n_ele, _B, _H)
    for t, c in coords.items():
        if abs(c[0]) < 1e-9:
            ops.fix(t, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if "20" in ele_name:
        for t, f in _tip_face_loads(conns, coords, P).items():
            ops.load(t, 0.0, 0.0, f)
    else:                                    # hex8: 4 tip nodes, equal shares
        tips = [t for t, c in coords.items() if abs(c[0] - _L) < 1e-9]
        assert len(tips) == 4
        for t in tips:
            ops.load(t, 0.0, 0.0, P / 4.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    tips = [t for t, c in coords.items() if abs(c[0] - _L) < 1e-9]
    uz = np.mean([ops.nodeDisp(t)[2] for t in tips])
    return uz, conns, coords


def _delta_timoshenko(e_mod, nu, P=1.0):
    I = _B * _H ** 3 / 12.0
    A = _B * _H
    G = e_mod / (2.0 * (1.0 + nu))
    kappa = 5.0 / 6.0
    return P * _L ** 3 / (3.0 * e_mod * I) + P * _L / (kappa * G * A)


def test_S4_coarse_bending_onate():
    """1-layer 5-element cantilever, nu=0: Brick20 uri AND std >= 0.98 of
    Timoshenko; the H8 LadrunoBrick std locks (< 0.90); H8 eas reported."""
    delta = _delta_timoshenko(E, 0.0)
    got = {}
    got["uri"] = abs(_cantilever("LadrunoBrick20", _URI, 5, E, 0.0)[0]) / delta
    got["std"] = abs(_cantilever("LadrunoBrick20", _STD, 5, E, 0.0)[0]) / delta
    got["h8-std"] = abs(_cantilever("LadrunoBrick", [], 5, E, 0.0)[0]) / delta
    try:
        got["h8-eas"] = abs(_cantilever("LadrunoBrick",
                                        ["-formulation", "eas"], 5, E, 0.0)[0]) / delta
    except Exception:
        got["h8-eas"] = float("nan")
    print(f"[S4] tip/analytic: {got}")
    assert got["uri"] >= 0.98, f"uri coarse bending {got['uri']:.4f} < 0.98"
    assert got["uri"] <= 1.05
    assert got["std"] >= 0.98, f"std coarse bending {got['std']:.4f} < 0.98"
    assert got["h8-std"] < 0.90, (
        f"H8 std should shear-lock here ({got['h8-std']:.4f}) — the S4 contrast"
    )


def test_S5_near_incompressible_relief():
    """nu=0.4999 cantilever (4x1x1): std LOCKS (<= 0.90 — ESCALATION trigger
    per ADR §3.4 if it does not); uri relieves (>= 0.85 and uri-std >= 0.10 —
    "mostly relieves ... not a mixed element", §3.4). nu=0.3 control: both
    >= 0.95 with NO formulation gap (|uri-std| <= 0.03) — isolates nu as the
    moved variable; the shared ~3.5% shortfall at this coarse mesh is
    discretization + beam-reference error, both -> ~0.99 at n=16 (spec
    amendment 5)."""
    for nu_ctrl in (0.3,):
        d = _delta_timoshenko(E, nu_ctrl)
        r_uri = abs(_cantilever("LadrunoBrick20", _URI, 4, E, nu_ctrl)[0]) / d
        r_std = abs(_cantilever("LadrunoBrick20", _STD, 4, E, nu_ctrl)[0]) / d
        print(f"[S5] nu={nu_ctrl}: uri {r_uri:.4f} std {r_std:.4f}")
        assert r_uri >= 0.95 and r_std >= 0.95, "nu=0.3 control failed"
        assert abs(r_uri - r_std) <= 0.03, (
            f"formulation gap at nu=0.3 ({abs(r_uri - r_std):.4f}) — locking "
            f"should not discriminate at moderate nu"
        )

    nu = 0.4999
    d = _delta_timoshenko(E, nu)
    r_uri = abs(_cantilever("LadrunoBrick20", _URI, 4, E, nu)[0]) / d
    r_std = abs(_cantilever("LadrunoBrick20", _STD, 4, E, nu)[0]) / d
    print(f"[S5] nu={nu}: uri {r_uri:.4f} std {r_std:.4f}")
    assert r_std <= 0.90, (
        f"ESCALATION (plan 2.4): std does NOT lock at nu=0.4999 "
        f"({r_std:.4f} > 0.90) — contradicts ADR 72 §3.4; adjudicate before "
        f"changing this test."
    )
    assert r_uri >= 0.85, f"uri at nu=0.4999: {r_uri:.4f} < 0.85 (relief floor)"
    assert r_uri - r_std >= 0.10, (
        f"relief contrast uri-std = {r_uri - r_std:.4f} < 0.10 at nu=0.4999"
    )


def test_S6_barlow_superconvergence():
    """uri GP stresses beat std GP stresses against the analytic bending
    field, each set at its OWN stations (support element excluded)."""
    results = {}
    for name, extra, tbl in (("uri", _URI, ref.GP8_F), ("std", _STD, ref.GP27_F)):
        uz, conns, coords = _cantilever("LadrunoBrick20", extra, 5, E, 0.0)
        num, ana = [], []
        for e, conn in enumerate(conns):
            if e == 0:
                continue                      # St-Venant support zone
            X = np.asarray([coords[t] for t in conn], dtype=float)
            sig = np.asarray(ops.eleResponse(e + 1, "stresses"), dtype=float)
            ngp = tbl.shape[0]
            assert sig.size == 6 * ngp
            for L in range(ngp):
                N, _ = ref.shape_numeric(tuple(tbl[L, :3]))
                xg = N @ X
                # |sigma_xx| = P (L - x) |z - c| / I, sign fitted below
                a = 1.0 * (_L - xg[0]) * (xg[2] - _H / 2.0) / (_B * _H ** 3 / 12.0)
                num.append(sig[6 * L])
                ana.append(a)
        num, ana = np.asarray(num), np.asarray(ana)
        # one global sign (bending convention) fitted per formulation-shared field
        s = 1.0 if np.linalg.norm(num - ana) < np.linalg.norm(num + ana) else -1.0
        rel = np.linalg.norm(num - s * ana) / np.linalg.norm(ana)
        results[name] = rel
    print(f"[S6] RMS rel sigma_xx error: uri {results['uri']:.4e} "
          f"vs std {results['std']:.4e}")
    assert results["uri"] < results["std"], (
        f"Barlow: uri ({results['uri']:.3e}) must beat std ({results['std']:.3e})"
    )


# ==========================================================================
# S7 — recorder seam (debt a): uri self-describes, 8 stations round-trip
# ==========================================================================
def test_S7_uri_recorder_roundtrip(tmp_path):
    h5py = pytest.importorskip("h5py")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    out = str(tmp_path / "lb20_uri.ladruno")

    nodes = _build_single(_CORNERS, _URI)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADED.items():
        ops.load(n, fx, fy, fz)
    ops.recorder("ladruno", out, "-E", "stresses")
    _static_rig()
    assert ops.analyze(1) == 0

    ops.eleResponse(1, "forces")
    live = list(ops.eleResponse(1, "stresses"))
    assert len(live) == 48
    ops.wipe()

    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no .ladruno output produced"
    with h5py.File(files[0], "r") as h:
        buckets, model_groups = {}, {}

        def visit(name, obj):
            if not isinstance(obj, h5py.Group):
                return
            if "ON_ELEMENTS/stresses/" in name and name.count("/") == 4:
                buckets[name.rsplit("/", 1)[-1]] = obj
            if "/MODEL/ELEMENTS/" in name and name.rsplit("/", 2)[-2] == "ELEMENTS":
                model_groups[name.rsplit("/", 1)[-1]] = obj

        h.visititems(visit)
        keys = [k for k in buckets if "LadrunoBrick20" in k]
        assert keys, f"no LadrunoBrick20 bucket in {sorted(buckets)}"
        data = np.asarray(buckets[keys[0]]["DATA"][...], dtype=float)
        assert data.shape[-1] == 6 * 8, (
            f"uri must record 8 GP stations, got shape {data.shape} — the "
            "recorder fell through to the hardwired 27-pt class-tag arm?"
        )
        assert data[-1].reshape(-1) == pytest.approx(live, rel=1e-12, abs=1e-14)

        mkeys = [k for k in model_groups if "LadrunoBrick20" in k]
        assert mkeys
        mg = model_groups[mkeys[0]]
        assert int(np.atleast_1d(mg.attrs["NUM_GP"])[0]) == 8

        # GP_PARAM == the GP8 table (z-fastest lexicographic — the Brick /
        # Hexahedron_GaussLegendre_2 order, == Ladruno::hex20::GP8)
        gp = np.asarray(mg["QUADRATURE"]["GP_PARAM"][...], dtype=float).reshape(8, 3)
        assert np.allclose(gp, ref.GP8_F[:, :3], atol=1e-14), (
            "GP_PARAM != GP8 (ordering seam broken)"
        )

        ggp = np.asarray(mg["GLOBAL_GP_COORDS"][...], dtype=float).reshape(8, 3)
        X = np.asarray([nodes[t] for t in _CONN], dtype=float)
        for L in range(8):
            N, _ = ref.shape_numeric(tuple(ref.GP8_F[L, :3]))
            assert np.allclose(ggp[L], N @ X, atol=1e-12)


def test_S7_battery_has_no_bare_27_pin_for_uri():
    """Hygiene: every uri response in THIS battery already asserts 8-sized
    trees; this test pins the element-side count directly."""
    _build_single(_CUBE_CORNERS, _URI)
    assert len(ops.eleResponse(1, "stresses")) == 48
    assert len(ops.eleResponse(1, "strains")) == 48
    _build_single(_CUBE_CORNERS, _STD)
    assert len(ops.eleResponse(1, "stresses")) == 162


# ==========================================================================
# S8 — surface & serialization
# ==========================================================================
def test_S8_uri_accepted_and_reported(tmp_path):
    _build_single(_CUBE_CORNERS, _URI)
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags, "-formulation uri must now be ACCEPTED (P2)"
    # printModel('-JSON') with no -file emits NOTHING (upstream quirk: the
    # interpreter's -JSON branch only sets the flag; Domain::Print only runs
    # in the filename branch) — assert on the file form instead.
    out = tmp_path / "model.json"
    ops.printModel("-JSON", "-file", str(out))
    text = out.read_text(errors="replace")
    assert '"formulation": "uri"' in text, (
        f"Print JSON must report uri; got {text[:400]!r}"
    )


def test_S8_reduced_alias_is_uri():
    nodes = _build_single(_CUBE_CORNERS, _URI)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    _static_rig()
    assert ops.analyze(1) == 0
    K_uri = list(ops.printA("-ret"))

    nodes = _build_single(_CUBE_CORNERS, ["-formulation", "reduced"])
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    _static_rig()
    assert ops.analyze(1) == 0
    K_red = list(ops.printA("-ret"))
    assert K_uri == K_red, "'reduced' must be bit-identical to 'uri'"


def test_S8_unknown_formulation_refused():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, c in _hex20_nodes(_CUBE_CORNERS).items():
        ops.node(tag, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    try:
        ops.element("LadrunoBrick20", 1, *_CONN, 1, "-formulation", "i14")
    except Exception:
        pass
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags


def test_S8_hourglass_still_hard_error_under_uri():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, c in _hex20_nodes(_CUBE_CORNERS).items():
        ops.node(tag, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    try:
        ops.element("LadrunoBrick20", 1, *_CONN, 1,
                    "-formulation", "uri", "-hourglass", "stiffness", 0.05)
    except Exception:
        pass
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags, "-hourglass stays a hard error under uri (ADR §2.2)"


def test_S8_database_roundtrip_preserves_uri():
    """sendSelf/recvSelf keeps ordinal 1: the restored element must produce
    the same 8-station stresses and the same solution (probe_fn compares)."""
    def build():
        nodes = _build_single(_CORNERS, _URI)
        for n in _BOTTOM:
            ops.fix(n, 1, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for n, (fx, fy, fz) in _LOADED.items():
            ops.load(n, fx, fy, fz)
        _static_rig()
        assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[5, 6, 7, 8, 13, 18], ndf=3,
                       dbname="lb20_uri_rt",
                       probe_fn=lambda: ops.eleResponse(1, "stresses"))


def test_S8_advisory_probe_alive_under_uri(capfd):
    """The U1 damage-advisory probe is formulation-independent; under uri the
    non-firing branch must stay silent (the fires-once branch is pinned by the
    P1 battery — the flag is PROCESS-once, so it cannot be re-asserted here
    without cross-file ordering flakiness; code path identical by design)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, c in _hex20_nodes(_CUBE_CORNERS).items():
        ops.node(tag, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, 1, *_URI)
    out = capfd.readouterr()
    assert "theory-shaky" not in (out.out + out.err)


# ==========================================================================
# S9 — assembly cost (REPORT + loose CI-safe gate)
# ==========================================================================
def test_S9_assembly_cost_ratio():
    """Tangent-formation time uri/std on a 4x4x4 block. The pure-assembly
    expectation is ~0.3-0.4x (8 vs 27 material points); the eleResponse
    round-trip adds equal absolute overhead per call to both sides, so the
    gate is a LOOSE < 0.8 (CI-safe) and the measured number is the report."""
    times = {}
    for name, extra in (("std", _STD), ("uri", _URI)):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
        conns, coords = _hex20_block("LadrunoBrick20", 1, extra, 4, 4, 4, 1, 1, 1)
        # warm-up (geometry caches, material clones touched once)
        for e in range(1, len(conns) + 1):
            ops.eleResponse(e, "stiffness")
        t0 = time.perf_counter()
        for _ in range(3):
            for e in range(1, len(conns) + 1):
                ops.eleResponse(e, "stiffness")
        times[name] = time.perf_counter() - t0
    ratio = times["uri"] / times["std"]
    print(f"[S9] tangent-assembly time: uri {times['uri']:.3f}s "
          f"std {times['std']:.3f}s ratio {ratio:.3f} (expect ~0.3-0.4 pure)")
    assert ratio < 0.8, f"uri/std assembly ratio {ratio:.3f} >= 0.8"
