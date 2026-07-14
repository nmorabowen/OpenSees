"""LadrunoBrick20 (Ladruno 20-node serendipity quadratic hex, classTag 33018) —
Zone-A element battery. ADR 72 P1.

The headline gate is the **reduce-to anchor**: ``-formulation std`` (27-pt full
Gauss) must reproduce upstream ``Twenty_Node_Brick`` (``20NodeBrick``) on the
SAME distorted 20-node hex + ElasticIsotropic to ~1e-12:

  * solved displacements, resisting forces, and per-GP stresses (162 values,
    GP L <-> GP L — the brcshl pairing pinned at P0),
  * the assembled free-DOF stiffness matrix (printA, static tangent),
  * the assembled transient tangent K + c3*M (printA under Newmark), which —
    given K equality — pins the consistent mass to the same tolerance.

On top of that (ADR 72 §6 P1 row):

  * constant-strain patch on a 2-element distorted straight-edged mesh
    (all 6 strain components at every one of the 54 GPs),
  * rank/6-RBM: a free element has exactly 6 zero-energy modes (rank 54),
  * SelfWeight/-b vs closed form, ASSERTING the corner-share sign pathology:
    the consistent H20 body load pulls the 8 corners BACKWARDS (-V/8 each,
    mid-edges +V/6; Hughes §7.3.2 / Cook Fig. 13.3-4 3D analog) — users
    hand-lumping gravity onto quadratic-mesh nodes get the wrong answer,
  * lch = cbrt(V),
  * sendSelf/recvSelf database round-trip,
  * parser surface: 'uri' ACCEPTED (P2 — contract battery in
    test_ladrunoBrick20_uri.py), '-hourglass' rejected (never, by design),
    '-lumped' (accepted at parse with a process-once NOTICE; getMass errors
    "lands P3" once per process and falls back), trailing '-damp' (refused),
  * U1 advisory: a PROCESS-once opserr advisory when a "damage"-channel material
    (ASDConcrete3D-class) is attached; NEVER for ElasticIsotropic/LadrunoJ2,
  * R3 hardening: non-positive detJ (mid-edge node past the quarter point)
    fails LOUDLY at setDomain + analyze; a non-3D material is refused at the
    factory; 'material' eleResponse with no GP number is guarded (no crash),
  * LadrunoRecorder round-trip: GP_PARAM == the brcshl 27-pt table and
    GLOBAL_GP_COORDS == sum_i N_i(xi_L) X_i from an independent Python
    serendipity evaluation.

Plan: Ladruno_implementation/72_ladruno_second_order_brick_adr.md (+ the
_adr72_implementation_plan companion).
"""
import glob
import os

import pytest

from _testbed import ops
from _testbed.fem_checks import zero_energy_mode_count, n_rigid_body_modes
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

E, NU = 1000.0, 0.3
RHO = 2.0e-3

# --------------------------------------------------------------------------
# fixtures: a deliberately distorted (positive-J, straight-edged) 20-node hex.
# Corners in the upstream Twenty_Node_Brick / brcshl order (lower face CCW,
# upper face CCW), mid-edge nodes at EXACT edge midpoints (straight edges —
# quadratic completeness in physical coords is only guaranteed there; ADR §3.1).
# --------------------------------------------------------------------------
_CORNERS = [
    (0.00, 0.00, 0.00),   # 1
    (1.00, 0.10, 0.00),   # 2
    (1.10, 1.00, 0.10),   # 3
    (0.05, 0.95, 0.00),   # 4
    (0.00, 0.05, 1.00),   # 5
    (1.00, 0.00, 1.05),   # 6
    (1.05, 1.00, 1.10),   # 7
    (0.00, 1.00, 0.95),   # 8
]
# edge list in the brcshl mid-node order: lower ring 1-2,2-3,3-4,4-1;
# upper ring 5-6,6-7,7-8,8-5; vertical 1-5,2-6,3-7,4-8  (1-based corners)
_EDGES = [(1, 2), (2, 3), (3, 4), (4, 1),
          (5, 6), (6, 7), (7, 8), (8, 5),
          (1, 5), (2, 6), (3, 7), (4, 8)]


def _hex20_nodes(corners):
    """{tag: coords} for a straight-edged 20-node hex from its 8 corners
    (tags 1..20 in brcshl order; mid-edges at exact midpoints)."""
    nodes = {i + 1: corners[i] for i in range(8)}
    for k, (a, b) in enumerate(_EDGES):
        ca, cb = corners[a - 1], corners[b - 1]
        nodes[9 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    return nodes


_NODES = _hex20_nodes(_CORNERS)
_CONN = list(range(1, 21))
_BOTTOM = [1, 2, 3, 4, 9, 10, 11, 12]          # the zeta=-1 face (fixed)
_LOADED = {                                     # mixed loads: corners + mids
    5: (10.0, 0.0, 0.0),
    6: (0.0, 8.0, 0.0),
    7: (0.0, 0.0, -12.0),
    8: (5.0, 5.0, 5.0),
    13: (-3.0, 2.0, 0.0),
    15: (0.0, -4.0, 6.0),
    18: (1.0, 1.0, 1.0),
    20: (0.0, 0.0, -5.0),
}


def _build_common(rho=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    if rho:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU, rho)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)


def _static_rig():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


def _solve_hex20(ele_name, extra_args, rho=0.0):
    """Single distorted hex20, bottom face fixed, mixed nodal loads; returns
    (disps[60], stresses[162], forces[60], K_free flattened via printA)."""
    _build_common(rho=rho)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)

    ops.element(ele_name, 1, *_CONN, 1, *extra_args)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADED.items():
        ops.load(n, fx, fy, fz)

    _static_rig()
    assert ops.analyze(1) == 0

    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    forces = list(ops.eleResponse(1, "forces"))       # before stresses (lazy-strain safe)
    stresses = list(ops.eleResponse(1, "stresses"))
    K = list(ops.printA("-ret"))
    return disps, stresses, forces, K


def _transient_tangent(ele_name, extra_args, dt=1.0e-2):
    """Assembled Newmark tangent K + c3*M (free DOFs, flattened) — with K
    pinned separately this pins the consistent M."""
    _build_common(rho=RHO)
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.element(ele_name, 1, *_CONN, 1, *extra_args)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, dt) == 0
    return list(ops.printA("-ret"))


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x-y):.3e} > {tol:.3e})"
        )


def _assert_matrix_close(a, b, label, rtol=5e-13):
    """Matrix gate: max |a-b| <= rtol * max|b| (conditioning-aware absolute
    floor pinned to the matrix scale — entrywise relative is meaningless for
    near-zero off-diagonal entries)."""
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    ref = max(abs(x) for x in b)
    worst = max(abs(x - y) for x, y in zip(a, b))
    assert worst <= rtol * ref, (
        f"{label}: max|diff|={worst:.3e} > {rtol:.0e}*max|ref|={rtol * ref:.3e} "
        f"(relative {worst / ref:.3e})"
    )


# ==========================================================================
# 1. HEADLINE: reduce-to Twenty_Node_Brick (~1e-12)
# ==========================================================================
def test_reduce_to_twenty_node_brick_solution():
    """std reproduces upstream 20NodeBrick: displacements, resisting force,
    per-GP stresses (GP L <-> GP L, 27x6) to ~1e-12 on a distorted hex."""
    ref = _solve_hex20("20NodeBrick", [])
    ours = _solve_hex20("LadrunoBrick20", ["-formulation", "std"])
    _assert_close(ours[0], ref[0], "disp", rtol=1e-12, atol=1e-15)
    _assert_close(ours[2], ref[2], "force", rtol=1e-12, atol=1e-10)
    # per-GP pairing: 27 GPs x 6 components, same brcshl ordering both sides
    assert len(ours[1]) == 27 * 6 and len(ref[1]) == 27 * 6
    _assert_close(ours[1], ref[1], "stressGP", rtol=1e-12, atol=1e-12)


def test_reduce_to_twenty_node_brick_stiffness():
    """Assembled free-DOF K (printA static tangent) matches to ~1e-12."""
    ref = _solve_hex20("20NodeBrick", [])
    ours = _solve_hex20("LadrunoBrick20", [])
    assert len(ref[3]) == 36 * 36
    _assert_matrix_close(ours[3], ref[3], "K")


def test_reduce_to_twenty_node_brick_mass():
    """Assembled Newmark tangent K + c3*M matches to ~1e-12; with K pinned by
    the static gate this pins the 27-pt consistent mass."""
    ref = _transient_tangent("20NodeBrick", [])
    ours = _transient_tangent("LadrunoBrick20", [])
    _assert_matrix_close(ours, ref, "K+c3*M")


def test_default_formulation_is_std():
    ref = _solve_hex20("LadrunoBrick20", ["-formulation", "std"])
    bare = _solve_hex20("LadrunoBrick20", [])
    _assert_close(bare[0], ref[0], "disp", rtol=0.0, atol=0.0)


# ==========================================================================
# 2. constant-strain patch — 2-element distorted straight-edged mesh
# ==========================================================================
# imposed constant strain (all six components exercised)
_EPS = [1.0e-4, -2.0e-4, 1.5e-4, 8.0e-5, -6.0e-5, 4.0e-5]  # xx yy zz gxy gyz gzx


def _u_exact(x, y, z):
    exx, eyy, ezz, gxy, gyz, gzx = _EPS
    return (exx * x + 0.5 * gxy * y + 0.5 * gzx * z,
            0.5 * gxy * x + eyy * y + 0.5 * gyz * z,
            0.5 * gzx * x + 0.5 * gyz * y + ezz * z)


def test_constant_strain_patch_two_elements():
    """Prescribe the exact linear field on every node of a 2-element distorted
    straight-edged mesh; every one of the 2x27 GPs must carry the imposed
    constant strain in all 6 components (and the matching constant stress)."""
    # element 1: the distorted hex; element 2: stacked on its top face
    lower = _CORNERS
    upper = [
        (0.00, 0.05, 1.00), (1.00, 0.00, 1.05),   # = lower top corners
        (1.05, 1.00, 1.10), (0.00, 1.00, 0.95),
        (0.10, 0.00, 2.00), (1.05, 0.05, 2.10),   # new distorted top
        (1.00, 1.05, 2.05), (0.05, 1.00, 1.95),
    ]
    n1 = _hex20_nodes(lower)
    n2raw = _hex20_nodes(upper)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords_of = {}
    tag_of = {}
    next_tag = 1

    def add_node(c):
        nonlocal next_tag
        key = tuple(round(v, 12) for v in c)
        if key not in tag_of:
            ops.node(next_tag, *c)
            tag_of[key] = next_tag
            coords_of[next_tag] = c
            next_tag += 1
        return tag_of[key]

    conn1 = [add_node(n1[k]) for k in range(1, 21)]
    conn2 = [add_node(n2raw[k]) for k in range(1, 21)]

    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoBrick20", 1, *conn1, 1)
    ops.element("LadrunoBrick20", 2, *conn2, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in coords_of.items():
        ux, uy, uz = _u_exact(x, y, z)
        ops.sp(t, 1, ux)
        ops.sp(t, 2, uy)
        ops.sp(t, 3, uz)

    # Penalty: fully-prescribed rig — Transformation/Plain would condense to a
    # zero-size SOE and fatally exit (LEDGER_quirks)
    ops.constraints("Penalty", 1.0e15, 1.0e15)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # closed-form constant stress for the imposed strain
    lam = E * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E / (2 * (1 + NU))
    tr = _EPS[0] + _EPS[1] + _EPS[2]
    sig_exact = [lam * tr + 2 * mu * _EPS[0],
                 lam * tr + 2 * mu * _EPS[1],
                 lam * tr + 2 * mu * _EPS[2],
                 mu * _EPS[3], mu * _EPS[4], mu * _EPS[5]]

    for ele in (1, 2):
        ops.eleResponse(ele, "forces")
        eps = list(ops.eleResponse(ele, "strains"))
        sig = list(ops.eleResponse(ele, "stresses"))
        assert len(eps) == 162 and len(sig) == 162
        for gp in range(27):
            _assert_close(eps[6 * gp: 6 * gp + 6], _EPS,
                          f"patch strain ele{ele} gp{gp}", rtol=1e-8, atol=1e-13)
            _assert_close(sig[6 * gp: 6 * gp + 6], sig_exact,
                          f"patch stress ele{ele} gp{gp}", rtol=1e-8, atol=1e-9)


# ==========================================================================
# 3. rank / 6 rigid-body modes (rank 54 at the 27-pt rule)
# ==========================================================================
def test_free_element_has_six_rigid_body_modes():
    _build_common()
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    count, eigs = zero_energy_mode_count(_CONN, ndf=3)
    assert count == n_rigid_body_modes(3), (
        f"{count} zero-energy modes (want 6 => rank 54); eigs={eigs[:10]}..."
    )


# ==========================================================================
# 4. SelfWeight / -b vs closed form — the corner-share sign pathology
# ==========================================================================
# On the UNIT CUBE the consistent H20 body-load shares are exact closed form:
#   corners:   b * integral(N_corner) dV = b * (-V/8)   <- SIGN-FLIPPED
#   mid-edges: b * integral(N_mid)    dV = b * (+V/6)
# (row-sum-mass analog pinned at P0: -M/8 / +M/6; Hughes §7.3.2.) A downward
# body force therefore pushes the corners UP — this is why hand-lumped gravity
# on quadratic meshes is wrong by construction (ADR 72 §3.8). We assert the
# element's resisting force at u ~= 0 equals -F_consistent, sign included.
def _selfweight_forces(use_eleload):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    cube = _hex20_nodes([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    for t, c in cube.items():
        ops.node(t, *c)
        ops.fix(t, 1, 1, 1)                      # fully fixed => u ~= 0
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoBrick20", 1, *_CONN, 1, "-b", 0.0, 0.0, -1.0)

    if use_eleload:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-selfWeight", 1.0, 1.0, 1.0)

    # Penalty: fully-prescribed rig (LEDGER_quirks)
    ops.constraints("Penalty", 1.0e15, 1.0e15)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return list(ops.eleResponse(1, "forces"))


@pytest.mark.parametrize("use_eleload", [False, True])
def test_selfweight_consistent_corner_share_sign(use_eleload):
    forces = _selfweight_forces(use_eleload)
    bz, V = -1.0, 1.0
    # resisting force at u=0 is MINUS the consistent load
    f_corner = -bz * (-V / 8.0)     # = -0.125  (load +0.125 UP under gravity DOWN)
    f_mid = -bz * (V / 6.0)         # = +1/6
    for a in range(8):              # corners
        z = forces[3 * a + 2]
        assert abs(z - f_corner) < 1e-9, f"corner {a}: {z} != {f_corner}"
        assert z < 0.0, "corner share must be SIGN-FLIPPED (load pulls corners up)"
    for a in range(8, 20):          # mid-edges
        z = forces[3 * a + 2]
        assert abs(z - f_mid) < 1e-9, f"mid {a}: {z} != {f_mid}"
        assert z > 0.0
    # x/y untouched; total balances b*V exactly
    assert all(abs(forces[3 * a + d]) < 1e-9 for a in range(20) for d in (0, 1))
    total_z = sum(forces[3 * a + 2] for a in range(20))
    assert abs(total_z - (-bz * V)) < 1e-9      # sum(-F) = -b*V = +1


# ==========================================================================
# 5. lch = cbrt(V)
# ==========================================================================
@pytest.mark.parametrize("dims", [(1.0, 1.0, 1.0), (2.0, 0.5, 1.25)])
def test_char_length_is_cbrt_volume(dims):
    a, b, c = dims
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    box = _hex20_nodes([(0, 0, 0), (a, 0, 0), (a, b, 0), (0, b, 0),
                        (0, 0, c), (a, 0, c), (a, b, c), (0, b, c)])
    for t, crd in box.items():
        ops.node(t, *crd)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    lch = ops.eleResponse(1, "charLength")[0]
    assert lch == pytest.approx((a * b * c) ** (1.0 / 3.0), rel=1e-12)


# ==========================================================================
# 6. sendSelf/recvSelf database round-trip
# ==========================================================================
def test_database_roundtrip():
    def build():
        _build_common()
        for n in _BOTTOM:
            ops.fix(n, 1, 1, 1)
        ops.element("LadrunoBrick20", 1, *_CONN, 1, "-b", 0.5, -0.25, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for n, (fx, fy, fz) in _LOADED.items():
            ops.load(n, fx, fy, fz)
        _static_rig()
        assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[5, 6, 7, 8, 13, 18], ndf=3,
                       dbname="lb20_rt",
                       probe_fn=lambda: ops.eleResponse(1, "stresses"))


# ==========================================================================
# 7. parser rejections / API honesty
# ==========================================================================
def _refused(*extra):
    _build_common()
    try:
        ops.element("LadrunoBrick20", 1, *_CONN, 1, *extra)
    except Exception:
        pass  # raising is an acceptable refusal
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    return 1 not in tags


def test_formulation_uri_accepted_since_p2():
    """P2: '-formulation uri' is LIVE — the element builds and exposes the
    8-station response tree (the full uri contract lives in
    test_ladrunoBrick20_uri.py)."""
    _build_common()
    ops.element("LadrunoBrick20", 1, *_CONN, 1, "-formulation", "uri")
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags, "'uri' must be ACCEPTED (P2 shipped)"
    assert len(ops.eleResponse(1, "stresses")) == 8 * 6


def test_hourglass_hard_error_by_design():
    assert _refused("-hourglass", "stiffness", 0.05), (
        "-hourglass must be a hard error (ADR 72 §2.2: no control exists on "
        "H20 by design)"
    )
    assert _refused("-hourglass", "viscous"), "-hourglass (any flavour) refused"


def test_lumped_accepted_at_parse_but_getmass_errors(capfd):
    """'-lumped' parses (API stable) but getMass errors 'lands P3' and falls
    back to consistent — asserted via the opserr channel (capfd, not capsys)."""
    _build_common(rho=RHO)
    ops.element("LadrunoBrick20", 1, *_CONN, 1, "-lumped")
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags, "'-lumped' must be ACCEPTED at parse time"

    capfd.readouterr()                      # drop construction chatter
    for n in _CONN:
        ops.mass(n, 1e-12, 1e-12, 1e-12)    # tiny floor so eigen is well posed
    ops.eigen("-fullGenLapack", 6)          # triggers getMass
    out = capfd.readouterr()
    text = out.out + out.err
    assert "lands in P3" in text and "-lumped" in text, (
        f"getMass with massType=1 must error clearly; got: {text[:400]!r}"
    )


# ==========================================================================
# 8. U1 advisory — softening (damage-channel) material on a quadratic element
# ==========================================================================
_FT = 3.0
_TE = [_FT / E, 1.0e-3, 8.0e-3]
_TS = [_FT, 1.0, 0.05]
_TD = [0.0, 1.0 - 1.0 / (E * 1.0e-3), 1.0 - 0.05 / (E * 8.0e-3)]


# UNIQUE process-wide tag: ASDConcrete3D keeps a STATIC HardeningLawStorage
# keyed by material tag with store-if-absent semantics
# (ASDConcrete3DMaterial.cpp::HardeningLawStorage::store) — it survives
# ops.wipe(), so the FIRST material built with a tag pins that tag's hardening
# law for the whole pytest process. Reusing a common tag (e.g. 1) here would
# poison test_ladrunoBrick_asdconcrete.py's tag-1 material (different Gf) and
# break its mesh-objectivity gate. See LEDGER_quirks.md.
_ASD_TAG = 337218


def _asdconcrete(tag):
    ops.nDMaterial("ASDConcrete3D", tag, E, NU, "-rho", 0.0,
                   "-Te", *_TE, "-Ts", *_TS, "-Td", *_TD,
                   "-Ce", -_FT / E, -1.0e-3,
                   "-Cs", -_FT, -1.0,
                   "-Cd", 0.0, 0.0,
                   "-autoRegularization", 1.0)


_ADVISORY_KEY = "ADR 72"     # distinctive advisory fragment


def _advisory_count(text):
    return text.count("theory-shaky")


# ORDERING CONSTRAINT (F3): the advisory flag is PROCESS-once (a class-level
# static — a 100k-element mesh prints once, not 100k times). Non-firing
# materials never trip the flag, but the never-fires test runs FIRST in file
# order anyway so it exercises the probe with the flag guaranteed untripped;
# the ASD fires-once test below is then the first (and only) trip in the
# process. Do not reorder these two tests.
@pytest.mark.parametrize("mat", ["elastic", "ladrunoJ2"])
def test_damage_advisory_never_for_non_damage_materials(mat, capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    else:
        K = E / (3 * (1 - 2 * NU))
        G = E / (2 * (1 + NU))
        ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce",
                       250.0, 0.0, 0.0, 2000.0)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    out = capfd.readouterr()
    text = out.out + out.err
    assert _advisory_count(text) == 0, (
        f"advisory must NEVER fire for {mat}; got {text[:400]!r}"
    )


def test_damage_advisory_fires_once_for_asdconcrete(capfd):
    """Process-once semantics: this is the FIRST crack-band material in the
    process, so the advisory fires here — exactly once — and never again."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _asdconcrete(_ASD_TAG)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, _ASD_TAG)   # setDomain fires the probe
    out = capfd.readouterr()
    text = out.out + out.err
    assert _advisory_count(text) == 1, (
        f"damage advisory must fire EXACTLY once at setDomain; got "
        f"{_advisory_count(text)} in {text[:400]!r}"
    )
    assert _ADVISORY_KEY in text
    # advisory only — the element must exist and the run proceeds
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags

    # process-once: a SECOND crack-band element must NOT print again
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _asdconcrete(_ASD_TAG)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, _ASD_TAG)
    out = capfd.readouterr()
    text = out.out + out.err
    assert _advisory_count(text) == 0, (
        f"advisory is process-once; second element must be silent, got "
        f"{text[:400]!r}"
    )


# ==========================================================================
# 9. LadrunoRecorder round-trip — GP rule + GLOBAL_GP_COORDS oracle
# ==========================================================================
# The independent basis / GP tables come from the P0 sympy oracle
# tests/hex20_reference.py (ADR 72 §6 P2 debt c: import, never re-transcribe —
# ordering corrections can't drift). sympy is in the Zone-A env; without it
# only this oracle-backed test skips.
def test_ladruno_recorder_gp_roundtrip(tmp_path):
    h5py = pytest.importorskip("h5py")
    np = pytest.importorskip("numpy")
    ref = pytest.importorskip("hex20_reference")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    out = str(tmp_path / "lb20.ladruno")

    _build_common()
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADED.items():
        ops.load(n, fx, fy, fz)
    ops.recorder("ladruno", out, "-E", "stresses")
    _static_rig()
    assert ops.analyze(1) == 0

    ops.eleResponse(1, "forces")
    live = list(ops.eleResponse(1, "stresses"))
    ops.wipe()                                # flush/close the HDF5

    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no .ladruno output produced"

    with h5py.File(files[0], "r") as h:
        buckets = {}       # RESULTS/ON_ELEMENTS/stresses/<bucket>  (DATA)
        model_groups = {}  # MODEL_STAGE[*]/MODEL/ELEMENTS/<bucket> (GP metadata)

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
        g = buckets[keys[0]]

        # 27 GPs x 6 components, recorded == live
        data = np.asarray(g["DATA"][...], dtype=float)
        assert data.shape[-1] == 6 * 27, data.shape
        rec = data[-1].reshape(-1)
        assert rec == pytest.approx(live, rel=1e-12, abs=1e-14)

        # model-stage element group carries the quadrature metadata
        mkeys = [k for k in model_groups if "LadrunoBrick20" in k]
        assert mkeys, f"no LadrunoBrick20 model group in {sorted(model_groups)}"
        mg = model_groups[mkeys[0]]
        assert int(np.atleast_1d(mg.attrs["NUM_GP"])[0]) == 27

        # GP_PARAM == the brcshl 27-pt table (rule reuse: Hexahedron_GaussLegendre_3)
        gp = np.asarray(mg["QUADRATURE"]["GP_PARAM"][...], dtype=float).reshape(27, 3)
        assert np.allclose(gp, ref.GP27_F[:, :3], atol=1e-14), "GP_PARAM != brcshl table"

        # GLOBAL_GP_COORDS == sum N_i(xi_L) X_i (oracle sympy-derived basis)
        ggp = np.asarray(mg["GLOBAL_GP_COORDS"][...], dtype=float).reshape(27, 3)
        X = np.asarray([_NODES[t] for t in _CONN], dtype=float)
        for L in range(27):
            N, _ = ref.shape_numeric(tuple(ref.GP27_F[L, :3]))
            xg = N @ X
            assert np.allclose(ggp[L], xg, atol=1e-12), (
                f"GLOBAL_GP_COORDS[{L}] {ggp[L]} != {xg}"
            )


# ==========================================================================
# 10. R3 hardening — detJ gate, factory 3D validation, -damp, response guard
# ==========================================================================
def test_nonpositive_jacobian_fails_loudly(capfd):
    """A mid-edge node pushed past the quarter point makes detJ < 0 at the
    near-corner GPs. The setDomain-time gate (F1) must error LOUDLY (element
    tag + GP + detJ value) and the analysis must FAIL (update returns -1) —
    never a silent run over a sign-flipped volume."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    cube = _hex20_nodes([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    cube[9] = (0.05, 0.0, 0.0)   # mid-edge of edge 1-2, well past the 1/4 point
    for t, c in cube.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, 1)   # setDomain runs the detJ gate
    out = capfd.readouterr()
    text = out.out + out.err
    assert "Non-positive Jacobian" in text, (
        f"setDomain must report the non-positive detJ; got {text[:400]!r}"
    )

    # the run fails loudly, not silently
    for n in _BOTTOM:
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    _static_rig()
    assert ops.analyze(1) != 0, (
        "analysis over an inverted-geometry element must FAIL, not run silently"
    )


def test_valid_curved_quarter_point_is_not_flagged(capfd):
    """Control: a mildly curved (but valid) midside placement keeps detJ > 0
    and must NOT trip the gate."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    cube = _hex20_nodes([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    cube[9] = (0.45, 0.02, 0.0)     # mild curvature, well inside the 1/4 limit
    for t, c in cube.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    capfd.readouterr()
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    out = capfd.readouterr()
    text = out.out + out.err
    assert "Non-positive Jacobian" not in text
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags


def test_factory_rejects_non_3d_material(capfd):
    """F2: a material with no 'ThreeDimensional' copy (plane-stress-only
    PlaneStressUserMaterial) is refused at the FACTORY — no element created."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    # nstatevs=40 nprops=7: fc ft fcu epsc0 epscu epstu stc (plane-stress only)
    ops.nDMaterial("PlaneStressUserMaterial", 1, 40, 7,
                   20.7, 2.07, 4.14, -0.002, -0.008, 0.001, 0.08)
    capfd.readouterr()
    try:
        ops.element("LadrunoBrick20", 1, *_CONN, 1)
    except Exception:
        pass  # raising is an acceptable refusal
    out = capfd.readouterr()
    text = out.out + out.err
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags, "non-3D material must be refused at the factory"
    assert "ThreeDimensional" in text, (
        f"factory must explain the 3D-capability refusal; got {text[:400]!r}"
    )


def test_trailing_damp_rejected():
    """F6: a trailing '-damp' with no damping tag is an input error."""
    assert _refused("-damp"), "trailing -damp with no tag must be refused"


def test_material_response_without_gp_number_is_guarded():
    """F5: eleResponse(..., 'material') with no GP number must not crash
    (setResponse argc guard) — it returns an empty response."""
    _build_common()
    ops.element("LadrunoBrick20", 1, *_CONN, 1)
    try:
        res = ops.eleResponse(1, "material")
    except Exception:
        res = None  # a clean interpreter-level refusal is fine; a crash is not
    assert not res
