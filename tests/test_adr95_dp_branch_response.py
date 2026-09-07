"""ADR-95 P0 — the `ladrunoBranch` diagnostic response on UW `DruckerPrager`.

WHAT THIS PINS
--------------
ADR-95 asks which of three events kills the quadratic elements on the
Prandtl-Reissner deck (note `95_prandtl_reissner_quadratic_root_cause_plan.md`).
Hypothesis H1 names the UW Drucker-Prager **two-surface** corner: the model
carries a cone `f1` AND a tension cutoff `f2` at `I1 = T = sqrt(2/3)*sigma_y/rho`,
and at the campaign's apex regulariser `SY = 0.2 kPa` that cutoff sits at
`I1 ~ 1 kPa` — i.e. essentially at zero mean stress, so any Gauss point that
heaves into tension changes CONSTITUTIVE OPERATOR, not just state.  H2 names
loss of ellipticity (Rudnicki-Rice) instead.  Neither can be tested without a
way to read, per Gauss point, *which branch the return map actually took* and
*whether the acoustic tensor is still positive*.

`ladrunoBranch` is that readout.  It is a pure observer: `plastic_integrator`
only WRITES the bookkeeping, never reads it, and the standing collapse gate
`tests/test_r3_prandtl_collapse_gate.py` must stay bit-identical.

THE PAYLOAD (8 slots, and the whole point of asserting the width)
-----------------------------------------------------------------
    [0] branch    0 elastic, 1 f1 only (cone), 2 f2 only (cutoff), 3 corner
                  — the FINAL Jact after the while(!okay) loop
    [1] gamma0    plastic multiplier of f1
    [2] gamma1    plastic multiplier of f2
    [3] f1_trial  f1 at the trial state
    [4] f2_trial  f2 at the trial state
    [5] forced    1 if the `count > 3` forced-accept bailout fired
    [6] I1        first invariant of the RETURNED stress
    [7] detAmin   min over ~200 deterministic directions of det(n.D_ep.n)/(2G)^3

A campaign that read a 7- or 9-wide vector would silently mis-index every
column of its CSV, so the width is asserted before anything else.

HOW THE BRANCHES ARE REACHED
----------------------------
Every leg imposes a KNOWN HOMOGENEOUS STRAIN on a single element by prescribing
all 24 DOF with `sp()` under the **Penalty** handler.  Penalty, not
Transformation, on purpose: a single hex has every node on three faces, so
prescribing a uniform strain field constrains all 24 DOF and the Transformation
handler would hand the solver a zero-size system.  With the penalty ~10 orders
above the material stiffness the achieved strain is the target to ~1e-10, which
is far tighter than any branch predicate here.

The three targets are computed by hand from the model's own algebra, so a leg
that lands on the wrong branch is a real defect and not a tuning accident:

  (a) elastic     mild compression + small shear, f1 and f2 both < 0
  (b) cone        compression + large shear: f1 > 0, f2 << 0  -> branch 1
  (c) cutoff      triaxial TENSION with I1 >> T: f2 > 0, and because the cone
                  apex sits exactly at I1 = T, f1 > 0 there too — so the honest
                  expectation is branch 2 OR 3, and on this deck it is 3.
"""
import math
import os
import sys

import pytest

# ADR-95 rule in force: every probe runs the binary built in THIS worktree.
# `_testbed` does a bare `import opensees`, which an installed Ladruno .pth can
# answer; putting dist/bin first makes the local build win unaided (see the
# memory entry "Installed Ladruno .pth hijacks import opensees").
#
# LADRUNO_DIST_BIN overrides the directory.  It exists because a long campaign
# run holds `dist/bin/opensees.pyd` open for hours, and on Windows that BLOCKS
# the next build's copy step silently (build.bat prints "The process cannot
# access the file" and carries on), leaving a stale .pyd behind a fresh
# OpenSees.exe.  Pointing this at a staged copy lets a rebuild be verified
# without killing another agent's run.
_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# --- deck (ADR-95 §2: phi_txc = 20 deg, rho_bar = 0, SY = 0.2, no hardening) --
K_EL = 1.0e4
G_EL = 4.0e3
SY = 0.2
PHI_TXC = 20.0


def _rho_from_phi_txc(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


RHO = _rho_from_phi_txc(PHI_TXC)
RHO_BAR = 0.0
T_CUT = math.sqrt(2.0 / 3.0) * SY / RHO       # tension cutoff, I1 units

_HEX = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
_HEX20_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7)]
# BezierTet10 mid-edge (vertexA, vertexB) order for control points 5..10
_TET_EDGE_V = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

PENALTY = 1.0e14


def _dp_material(tag=1):
    # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta  -> perfectly plastic
    ops.nDMaterial("DruckerPrager", tag, K_EL, G_EL, SY, RHO, RHO_BAR,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _branch(tag, gp=1):
    r = ops.eleResponse(tag, "material", gp, "ladrunoBranch")
    r = list(r) if r else []
    assert len(r) == 8, (
        f"ladrunoBranch returned width {len(r)} on element {tag} GP {gp}; the "
        "ADR-95 sampler indexes 8 fixed columns and would mis-label every one"
    )
    return r


def _drive_uniform_strain(exx, eyy, ezz, gxy=0.0, gyz=0.0, gzx=0.0):
    """One LadrunoBrick with the linear field u_i = eps_ij x_j on ALL 24 DOF."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.node(i, x, y, z)
    _dp_material(1)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1,
                "-geom", "linear", "-formulation", "bbar")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.sp(i, 1, exx * x + 0.5 * gxy * y + 0.5 * gzx * z)
        ops.sp(i, 2, 0.5 * gxy * x + eyy * y + 0.5 * gyz * z)
        ops.sp(i, 3, 0.5 * gzx * x + 0.5 * gyz * y + ezz * z)

    ops.constraints("Penalty", PENALTY, PENALTY)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0, "uniform-strain driver failed to converge"
    return 1


# ---------------------------------------------------------------------------
def test_build_stamp_is_printed():
    """No ADR-95 number is attributable without the engine hash next to it."""
    assert hasattr(ops, "ladrunoBuild"), "ladrunoBuild() missing — wrong build"
    build = ops.ladrunoBuild()
    print(f"ADR-95 P0 engine build (ladrunoBuild): {build}")
    assert isinstance(build, str) and len(build) >= 7


def test_elastic_leg_is_branch_zero_and_elliptic():
    """(a) mild compression + small shear: f1, f2 < 0 by hand."""
    tag = _drive_uniform_strain(-2.0e-5, -2.0e-5, -2.0e-5, gxy=1.0e-5)
    for gp in range(1, 9):
        b = _branch(tag, gp)
        assert b[0] == 0.0, f"GP {gp}: expected elastic branch 0, got {b[0]}"
        assert b[1] == 0.0 and b[2] == 0.0, f"GP {gp}: elastic gammas non-zero: {b[1:3]}"
        assert b[3] < 0.0 and b[4] < 0.0, f"GP {gp}: trial f1/f2 not both < 0: {b[3:5]}"
        assert b[5] == 0.0, f"GP {gp}: forced-accept fired on an elastic step"
        assert b[7] > 0.0, (
            f"GP {gp}: detAmin = {b[7]} <= 0 with the ELASTIC tangent — the "
            "acoustic tensor of isotropic elasticity is positive definite, so "
            "this is a Voigt/4th-order mapping bug, not a material event"
        )
        # isotropic elasticity: min_n det(n.Ce.n)/(2G)^3 = (K + 4G/3)/(8G)
        expect = (K_EL + 4.0 * G_EL / 3.0) / (8.0 * G_EL)
        assert b[7] == pytest.approx(expect, rel=1e-9), (
            f"GP {gp}: detAmin {b[7]} != closed form {expect}"
        )


def test_cone_leg_is_branch_one_with_positive_gamma0():
    """(b) compression + large shear: f1 > 0, f2 far below the cutoff."""
    tag = _drive_uniform_strain(-1.0e-4, -1.0e-4, -1.0e-4, gxy=8.0e-4)
    for gp in range(1, 9):
        b = _branch(tag, gp)
        assert b[0] == 1.0, f"GP {gp}: expected cone branch 1, got {b[0]}"
        assert b[1] > 0.0, f"GP {gp}: gamma0 = {b[1]} not positive on the cone"
        assert b[3] > 0.0, f"GP {gp}: trial f1 = {b[3]} not positive"
        assert b[4] < 0.0, f"GP {gp}: trial f2 = {b[4]} — cutoff should be inactive"
        assert b[5] == 0.0, f"GP {gp}: forced-accept fired ({b[5]})"
        assert b[6] < 0.0, f"GP {gp}: I1 = {b[6]} should be compressive"


def test_tension_leg_takes_the_cutoff_or_corner_branch():
    """(c) triaxial tension with I1 >> T: the f2 / corner operator."""
    tag = _drive_uniform_strain(3.0e-4, 1.0e-4, 1.0e-4)
    for gp in range(1, 9):
        b = _branch(tag, gp)
        assert b[0] in (2.0, 3.0), (
            f"GP {gp}: expected cutoff (2) or corner (3), got branch {b[0]} — "
            "the tension cutoff is what ADR-95 H1 is about"
        )
        assert b[4] > 0.0, f"GP {gp}: trial f2 = {b[4]} not positive in tension"
        assert b[6] >= T_CUT, (
            f"GP {gp}: I1 = {b[6]} below the cutoff T = {T_CUT} yet the branch "
            f"is {b[0]}"
        )


def test_tension_cutoff_multiplier_is_structurally_zero_in_vanilla():
    """SENTINEL for the vanilla defect P0 found while instrumenting the branch.

    `plastic_integrator`'s residual assembly reads

        for (i = 0; i < 2; i++) {
            if      (Jact(i) == 1) { R(0) = ...; g(0,0) = ...; }
            else if (Jact(i) == 2) { R(1) = ...; g(1,1) = ...; }
        }

    but `Jact` only ever holds 0 or 1 — the `== 2` arm is DEAD.  So `R(1)` is
    never assembled, `g` stays lower-triangular with `g(1,1) = 1`, and
    `gamma(1)` (the tension-cutoff multiplier) comes out EXACTLY zero on every
    branch, corner included.  With `rho_bar = 0` the stress update
    `I1 -= 9K*rho_bar*gamma0 + 9K*gamma1` is then identically zero as well: the
    cutoff SELECTS a different consistent tangent but applies NO stress return,
    so `I1` stays wherever the trial state put it, arbitrarily far above `T`.

    That is exactly the asymmetry ADR-95 H1 is about — the corner is an
    OPERATOR switch, not a stress-return event — so it is pinned here.  If this
    test ever fails, someone has repaired the cutoff and H1's mechanism must be
    re-derived from scratch before any ADR-95 conclusion is carried forward.
    """
    tag = _drive_uniform_strain(3.0e-4, 1.0e-4, 1.0e-4)
    b = _branch(tag, 1)
    assert b[0] == 3.0, f"expected the corner branch here, got {b[0]}"
    assert b[2] == 0.0, (
        f"gamma1 = {b[2]} != 0 — the DruckerPrager tension-cutoff multiplier is "
        "no longer structurally zero; ADR-95 H1 must be re-derived"
    )
    i1_trial = 3.0 * K_EL * (3.0e-4 + 1.0e-4 + 1.0e-4)
    assert b[6] == pytest.approx(i1_trial, rel=1e-8), (
        f"I1 = {b[6]} != trial {i1_trial} — the cutoff now returns stress; "
        "ADR-95 H1 must be re-derived"
    )


def test_token_is_absent_on_a_material_that_does_not_define_it():
    """Guard against the token being answered by some generic fallback."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.node(i, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1)
    r = ops.eleResponse(1, "material", 1, "ladrunoBranch")
    assert not r, f"ElasticIsotropic answered ladrunoBranch with {list(r)}"


# ---------------------------------------------------------------------------
# element-side delegation: the campaign queries the SAME token on the quadratic
# elements, so each must forward an unknown material token to the NDMaterial.
# ---------------------------------------------------------------------------
def _elastic_tip_load(nfix, nload, ndof_nodes):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in nload:
        ops.load(n, 0.0, 0.0, -1.0e-3)
    for n in nfix:
        ops.fix(n, 1, 1, 1)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0


def test_brick20_forwards_the_token():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = {i + 1: _HEX[i] for i in range(8)}
    for k, (a, b) in enumerate(_HEX20_EDGES):
        ca, cb = _HEX[a], _HEX[b]
        coords[9 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
    _dp_material(1)
    ops.element("LadrunoBrick20", 1, *range(1, 21), 1)
    _elastic_tip_load(nfix=(1, 2, 3, 4, 9, 10, 11, 12),
                      nload=(5, 6, 7, 8), ndof_nodes=20)
    b = _branch(1, 1)
    assert b[0] == 0.0, f"LadrunoBrick20 GP1 not elastic under a tiny load: {b[0]}"
    assert b[7] > 0.0, f"LadrunoBrick20 detAmin = {b[7]} <= 0 in the elastic range"


def test_beziertet10_forwards_the_token():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    verts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    coords = {i + 1: verts[i] for i in range(4)}
    for k, (a, b) in enumerate(_TET_EDGE_V):
        ca, cb = verts[a], verts[b]
        coords[5 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
    _dp_material(1)
    ops.element("BezierTet10", 1, *range(1, 11), 1)
    # fix the z = 0 face (vertices 1,2,3 and their three mid-edges 5,6,7)
    _elastic_tip_load(nfix=(1, 2, 3, 5, 6, 7), nload=(4,), ndof_nodes=10)
    b = _branch(1, 1)
    assert b[0] == 0.0, f"BezierTet10 GP1 not elastic under a tiny load: {b[0]}"
    assert b[7] > 0.0, f"BezierTet10 detAmin = {b[7]} <= 0 in the elastic range"


def test_tennodetetrahedron_forwards_the_token():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    verts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    # TenNodeTetrahedron mid-side order: 12, 23, 31, 14, 24, 34 (0-based below)
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    coords = {i + 1: verts[i] for i in range(4)}
    for k, (a, b) in enumerate(edges):
        ca, cb = verts[a], verts[b]
        coords[5 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
    _dp_material(1)
    ops.element("TenNodeTetrahedron", 1, *range(1, 11), 1)
    _elastic_tip_load(nfix=(1, 2, 3, 5, 6, 7), nload=(4,), ndof_nodes=10)
    b = _branch(1, 1)
    assert b[0] == 0.0, f"TenNodeTetrahedron GP1 not elastic: {b[0]}"
    assert b[7] > 0.0, f"TenNodeTetrahedron detAmin = {b[7]} <= 0 in the elastic range"
