"""LadrunoUP ``-geom hypo`` (ADR 79 P2) — the u-p rate-form large-strain battery.

hypo on the Biot u-p element = the ADR-79 P1 brick lane extended to the
coupled blocks:

  * u-rows: spatial total-stress internal force int B'(sigma'_pushed - alpha*p*m) dv
            on the CURRENT configuration; K = int B'cB dv + total-stress
            geometric term (symmetric).
  * p-row coupling: the per-GP alpha*tr(deps)*dv/dt increment — the delta-lnJ
            continuity (tr(deps) sums to ln J in the midpoint limit and is
            EXACTLY zero on a rigid increment: the ADR-78 chord defect cannot
            arise by construction). S and alpha*Htilde go incremental with it
            (the ADR-78 GN22/GN11 pairing).
  * H/S:    rebuilt per evaluation on the current configuration; the storage
            coefficient evolves with the KINEMATIC porosity n(J) = 1-(1-n0)/J
            (always on under hypo); permeability evolution (Kozeny-Carman
            scaling) is a constitutive CHOICE and opt-in via -kozenyCarman.
  * seepage drive: spatial gradients are already global-frame — no pull-back.

Gates:
  1. parser guards (2D / bbar / TH / -kozenyCarman without hypo all fatal).
  2. STATIC rigid rotation of a saturated hex: sigma' = 0, f = 0, p = 0.
  3. TRANSIENT UNDRAINED rigid rotation (per-DOF Path series): excess p ~ 0 —
     the tr(deps) coupling is rigid-exact per step.
  4. consolidation under hypo at small strain reproduces linear
     (refinement-aware — the incremental p-row is a different discretization,
     the ADR-78 gate-5 pattern).
  5. drained-equivalence: saturated fully-drained gravity box + footing push
     under -geom hypo reproduces the dry LadrunoBrick -geom hypo twin with
     buoyant weight (roller sides; summed-reaction weight check catches the
     -b convention trap: UP -b is an ACCELERATION x rho_sat, brick -b a
     FORCE/VOLUME).
  6. porosity/J identity: drained uniaxial-strain compression -> per-GP
     n == 1-(1-n0)/J exactly (hypoState response), J tracks the imposed
     stretch; with -kozenyCarman the reported kcScale equals the KC formula
     exactly and consolidation gets SLOWER (direction gate).
  7. undrained large-strain storage sanity: prescribed uniaxial-strain
     compression at k~0 builds p ~ Qbar*alpha*|ln(lambda)| (the finite-strain
     volume bookkeeping).

Run:  py -3.12 -m pytest tests/test_ladruno_up_element_hypo.py -x -q
"""
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (the corot-battery idiom)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST)

pytestmark = [pytest.mark.zone_b]

_E = 200.0e3
_NU = 0.3
_RHO_SAT = 2.0
_RHO_W = 1.0
_G = 9.81
_KF = 2.2e6
_PORO = 0.4

_HEX_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.05),
    3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_HEX_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

_CUBE_NODES = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}


def _rot(thz, thy):
    cz, sz = math.cos(thz), math.sin(thz)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(thy), math.sin(thy)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return Ry @ Rz


def _up_hex(geom, nodes=_HEX_NODES, k=1.0e-3, body=(0.0, 0.0, 0.0),
            extra=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in nodes.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    args = ["-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
            "-perm", k, k, k, "-stab", "off", "-geom", geom, *extra]
    if body != (0.0, 0.0, 0.0):
        args += ["-body", *body]
    ops.element("LadrunoUP", 1, *_HEX_CONN, 1, *args)


def _static(tol=1.0e-10, iters=50):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


# ========================================================================== #
#  gate 1 — parser guards                                                    #
# ========================================================================== #
def test_parser_rejects_hypo_in_2d():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for tag, (x, y) in {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}.items():
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    with pytest.raises(Exception):
        ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, "-geom", "hypo")


def test_parser_rejects_hypo_with_bbar():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in _HEX_NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    with pytest.raises(Exception):
        ops.element("LadrunoUP", 1, *_HEX_CONN, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, 1e-4,
                    "-formulation", "bbar", "-geom", "hypo")


def test_parser_rejects_kozeny_carman_without_hypo():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in _HEX_NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    with pytest.raises(Exception):
        ops.element("LadrunoUP", 1, *_HEX_CONN, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, 1e-4, "-kozenyCarman")


# ========================================================================== #
#  gate 2 — STATIC rigid rotation: sigma' = 0, f = 0, p = 0                  #
# ========================================================================== #
@pytest.mark.parametrize("thz,thy", [(0.35, 0.25), (1.20, -0.80),
                                     (math.radians(179.0), 0.0)])
def test_static_rigid_rotation_stress_free_and_dry(thz, thy):
    _up_hex("hypo")
    R = _rot(thz, thy)
    ops.fix(1, 0, 0, 0, 1)              # pin one p (unique p solution = 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _HEX_NODES.items():
        d = R @ np.array([x, y, z]) - np.array([x, y, z])
        for i in range(3):
            ops.sp(tag, i + 1, float(d[i]))
    _static()
    assert ops.analyze(1) == 0, "rigid-rotation solve failed"

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    assert s.size == 48
    assert np.abs(s).max() <= 1.0e-7 * _E, (
        f"rigid rotation induced effective stress {np.abs(s).max():.3e}")
    p = np.array([ops.nodeDisp(n, 4) for n in range(1, 9)])
    assert np.abs(p).max() <= 1.0e-6, (
        f"rigid rotation induced static pore pressure {np.abs(p).max():.3e}")
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, (
        f"rigid rotation induced force {np.abs(f).max():.3e}")
    # hypoState sanity: rigid => J == 1 per GP
    hs = np.array(ops.eleResponse(1, "hypoState"), dtype=float).reshape(-1, 3)
    assert np.abs(hs[:, 0] - 1.0).max() <= 1.0e-10, "rigid rotation changed J"


# ========================================================================== #
#  gate 3 — TRANSIENT UNDRAINED rigid rotation: excess p ~ 0                 #
# ========================================================================== #
def test_transient_undrained_rigid_rotation_zero_excess_p():
    nsteps, dt = 30, 0.01
    times = [dt * (kk + 1) for kk in range(nsteps)]
    _up_hex("hypo", k=1.0e-8)           # undrained; no p BC

    ts = 1
    pat = 1
    for tag, (x, y, z) in _HEX_NODES.items():
        X = np.array([x, y, z])
        for dof in range(3):
            vals = []
            for kk in range(nsteps):
                th = (kk + 1) * (0.5 * math.pi) / nsteps
                vals.append(float((_rot(th, 0.0) @ X - X)[dof]))
            ops.timeSeries("Path", ts, "-time", *times, "-values", *vals,
                           "-useLast")
            ops.pattern("Plain", pat, ts)
            ops.sp(tag, dof + 1, 1.0)
            ts += 1
            pat += 1

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    pmax = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "undrained rigid-rotation step failed"
        p = np.array([ops.nodeDisp(n, 4) for n in range(1, 9)])
        pmax = max(pmax, np.abs(p).max())

    # tr(deps) == 0 exactly on a rigid increment (the HW midpoint property) =>
    # zero coupling drive; same discrimination logic as the corot gate.
    assert pmax <= 1.0e-4, f"undrained rigid rotation built excess p = {pmax:.3e}"

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    assert np.abs(s).max() <= 1.0e-6 * _E, (
        f"undrained rigid rotation built effective stress {np.abs(s).max():.3e}")


# ========================================================================== #
#  gate 4 — consolidation column: hypo -> linear at small strain             #
# ========================================================================== #
def _column_mesh(nz=8, lz=4.0):
    tags = {}
    t = 1
    for k in range(nz + 1):
        z = lz * k / nz
        for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1)):
            ops.node(t, float(x), float(y), z)
            tags[(k, (x, y))] = t
            t += 1
    conns = []
    for k in range(nz):
        conns.append([tags[(k, (0, 0))], tags[(k, (1, 0))], tags[(k, (1, 1))],
                      tags[(k, (0, 1))], tags[(k + 1, (0, 0))],
                      tags[(k + 1, (1, 0))], tags[(k + 1, (1, 1))],
                      tags[(k + 1, (0, 1))]])
    return tags, conns


def _consolidation_run(geom, nz=8, nsteps=20, dt=0.08, q=1.0, extra=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _column_mesh(nz=nz)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-5, 1e-5, 1e-5, "-geom", geom, *extra)
    top = nz
    for (k, xy), t in tags.items():
        fx, fy = 1, 1
        fz = 1 if k == 0 else 0
        fp = 1 if k == top else 0
        ops.fix(t, fx, fy, fz, fp)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for xy in ((0, 0), (1, 0), (1, 1), (0, 1)):
        ops.load(tags[(top, xy)], 0.0, 0.0, -q / 4.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        traj.append(ops.nodeDisp(tags[(top, (0, 0))], 3))
        traj.append(ops.nodeDisp(tags[(0, (0, 0))], 4))
    return np.array(traj)


def test_consolidation_hypo_reduces_to_linear():
    lin1 = _consolidation_run("linear", nsteps=20, dt=0.08)
    hyp1 = _consolidation_run("hypo", nsteps=20, dt=0.08)
    lin2 = _consolidation_run("linear", nsteps=40, dt=0.04)
    hyp2 = _consolidation_run("hypo", nsteps=40, dt=0.04)
    scale = np.abs(lin1).max()
    err1 = np.abs(hyp1 - lin1).max()
    err2 = np.abs(hyp2 - lin2).max()
    assert err1 <= 2.0e-2 * scale, (
        f"hypo-vs-linear gap too large at dt=0.08: {err1:.3e} vs {scale:.3e}")
    if err1 > 1.0e-8 * scale:
        assert err2 <= 0.80 * err1, (
            "hypo-vs-linear gap did not shrink under dt-halving "
            f"(err(dt)={err1:.3e}, err(dt/2)={err2:.3e})")


def test_consolidation_kozeny_carman_is_slower():
    # compaction reduces n => KC reduces k => the base excess p must decay
    # SLOWER with -kozenyCarman than without. Direction gate on the same mesh
    # (large q so the strain — and hence the KC effect — is non-trivial).
    plain = _consolidation_run("hypo", nsteps=30, dt=0.08, q=0.02 * _E)
    kc = _consolidation_run("hypo", nsteps=30, dt=0.08, q=0.02 * _E,
                            extra=("-kozenyCarman",))
    p_plain = plain[1::2]
    p_kc = kc[1::2]
    assert p_plain[0] > 0.0 and p_kc[0] > 0.0, "no excess pressure built; vacuous"
    # compare late-time residual pressure: KC must retain more
    assert p_kc[-1] > p_plain[-1] * 1.001, (
        f"KC did not slow consolidation (late p: kc={p_kc[-1]:.4e} vs "
        f"plain={p_plain[-1]:.4e})")


# ========================================================================== #
#  gate 5 — drained-equivalence vs a dry LadrunoBrick (both hypo)            #
# ========================================================================== #
_EQ_N = 4
_EQ_L = 4.0
_EQ_F = 2000.0


def _box_mesh(n, l):
    tags = {}
    t = 1
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(t, l * i / n, l * j / n, l * k / n)
                tags[(i, j, k)] = t
                t += 1
    conns = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                conns.append([tags[(i, j, k)], tags[(i + 1, j, k)],
                              tags[(i + 1, j + 1, k)], tags[(i, j + 1, k)],
                              tags[(i, j, k + 1)], tags[(i + 1, j, k + 1)],
                              tags[(i + 1, j + 1, k + 1)], tags[(i, j + 1, k + 1)]])
    return tags, conns


def _eq_fix_solids(tags, n, has_p):
    for (i, j, k), t in tags.items():
        fx = 1 if (i == 0 or i == n) else 0
        fy = 1 if (j == 0 or j == n) else 0
        fz = 1 if k == 0 else 0
        if k == 0:
            fx = fy = 1
        if has_p:
            fp = 1 if k == n else 0
            ops.fix(t, fx, fy, fz, fp)
        elif fx or fy or fz:
            ops.fix(t, fx, fy, fz)


def _eq_footing_loads(tags, n, ndf):
    c = n // 2
    w = _EQ_F / 16.0
    zeros = [0.0] * (ndf - 3)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            wt = w * (2.0 if di == 0 else 1.0) * (2.0 if dj == 0 else 1.0)
            ops.load(tags[(c + di, c + dj, n)], 0.0, 0.0, -wt, *zeros)


def _eq_solve_static():
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.2)
    ops.analysis("Static")


def test_drained_equivalence_up_hypo_vs_dry_brick_hypo():
    # UP leg: saturated, fully drained, gravity as ACCELERATION
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _box_mesh(_EQ_N, _EQ_L)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1.0, 1.0, 1.0,
                    "-body", 0.0, 0.0, -_G, "-geom", "hypo")
    _eq_fix_solids(tags, _EQ_N, has_p=True)
    _eq_solve_static()
    assert ops.analyze(5) == 0, "UP gravity stage failed"

    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in tags.values())
    W = _RHO_SAT * _G * _EQ_L**3
    assert abs(rz - W) <= 1.0e-6 * W, (
        f"UP total vertical reaction {rz:.6f} != saturated weight {W:.6f}")
    p_base = ops.nodeDisp(tags[(0, 0, 0)], 4)
    assert abs(p_base - _RHO_W * _G * _EQ_L) <= 1.0e-3 * (_RHO_W * _G * _EQ_L)

    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    _eq_footing_loads(tags, _EQ_N, ndf=4)
    c = _EQ_N // 2
    up_backbone = []
    for _ in range(5):
        assert ops.analyze(1) == 0, "UP footing stage failed"
        up_backbone.append(ops.nodeDisp(tags[(c, c, _EQ_N)], 3))

    # dry twin: buoyant weight as FORCE PER VOLUME
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    tags, conns = _box_mesh(_EQ_N, _EQ_L)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    gamma_b = (_RHO_SAT - _RHO_W) * _G
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoBrick", e, *conn, 1, "-formulation", "std",
                    "-geom", "hypo", "-b", 0.0, 0.0, -gamma_b)
    _eq_fix_solids(tags, _EQ_N, has_p=False)
    _eq_solve_static()
    assert ops.analyze(5) == 0, "brick gravity stage failed"

    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in tags.values())
    Wb = gamma_b * _EQ_L**3
    assert abs(rz - Wb) <= 1.0e-6 * Wb

    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    _eq_footing_loads(tags, _EQ_N, ndf=3)
    brick_backbone = []
    for _ in range(5):
        assert ops.analyze(1) == 0, "brick footing stage failed"
        brick_backbone.append(ops.nodeDisp(tags[(c, c, _EQ_N)], 3))

    up_b = np.array(up_backbone)
    br_b = np.array(brick_backbone)
    rel = np.abs(up_b - br_b).max() / np.abs(br_b).max()
    print(f"hypo drained-equivalence backbone rel gap = {rel:.3e}")
    assert rel <= 1.0e-2, (
        f"UP-hypo vs dry-brick-hypo drained backbones differ by {rel:.3e}")


# ========================================================================== #
#  gate 6 — porosity/J identity + Kozeny-Carman formula (hypoState)          #
# ========================================================================== #
def _uniaxial_compress(lam_final, nsteps, k=1.0, extra=()):
    """Drained uniaxial-strain compression of a unit cube: prescribed z-crush,
    roller sides, drained everywhere (all p fixed to 0)."""
    _up_hex("hypo", nodes=_CUBE_NODES, k=k, extra=extra)
    for tag, (x, y, z) in _CUBE_NODES.items():
        ops.fix(tag, 1, 1, 1 if z == 0.0 else 0, 1)   # fully drained: p == 0
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.timeSeries("Constant", 1)
    top = [t for t, (x, y, z) in _CUBE_NODES.items() if z == 1.0]
    for k_step in range(1, nsteps + 1):
        lam = 1.0 + (lam_final - 1.0) * k_step / nsteps
        pat = 100 + k_step
        ops.pattern("Plain", pat, 1)
        for t in top:
            ops.sp(t, 3, lam - 1.0)
        assert ops.analyze(1) == 0, f"compression step {k_step} failed"
        ops.remove("loadPattern", pat)


def test_porosity_tracks_J_and_kozeny_carman_formula():
    lam = 0.80
    _uniaxial_compress(lam, nsteps=25, extra=("-kozenyCarman",))
    hs = np.array(ops.eleResponse(1, "hypoState"), dtype=float).reshape(-1, 3)
    J, n, kc = hs[:, 0], hs[:, 1], hs[:, 2]
    # uniform uniaxial strain: J == lambda at every GP (midpoint-rule exact
    # kinematics: J is det F, computed exactly from the trial config)
    assert np.abs(J - lam).max() <= 1.0e-8, f"J != lambda: {J}"
    n_want = 1.0 - (1.0 - _PORO) / J
    assert np.abs(n - n_want).max() <= 1.0e-10, "porosity != 1-(1-n0)/J"
    kc_want = (n**3 / (1.0 - n)**2) * ((1.0 - _PORO)**2 / _PORO**3)
    assert np.abs(kc - kc_want).max() <= 1.0e-10 * np.abs(kc_want).max(), (
        "kcScale != Kozeny-Carman formula")
    assert kc.max() < 1.0, "compression must REDUCE the KC permeability scale"


# ========================================================================== #
#  gate 7 — undrained large-strain storage: p ~ Qbar*alpha*|ln lambda|       #
# ========================================================================== #
def test_undrained_compression_builds_lnJ_pressure():
    lam, nsteps, dt = 0.90, 40, 0.01
    _up_hex("hypo", nodes=_CUBE_NODES, k=1.0e-12)    # undrained, no p BC
    for tag, (x, y, z) in _CUBE_NODES.items():
        ops.fix(tag, 1, 1, 1 if z == 0.0 else 0, 0)  # p free everywhere
    top = [t for t, (x, y, z) in _CUBE_NODES.items() if z == 1.0]
    times = [dt * (kk + 1) for kk in range(nsteps)]
    ts = 1
    for t in top:
        vals = [(1.0 + (lam - 1.0) * (kk + 1) / nsteps) - 1.0
                for kk in range(nsteps)]
        ops.timeSeries("Path", ts, "-time", *times, "-values", *vals, "-useLast")
        ops.pattern("Plain", 100 + ts, ts)
        ops.sp(t, 3, 1.0)
        ts += 1
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    # NormUnbalance, not NormDispIncr: the undrained p scale (~Qbar*|ln lam| ~
    # 6e5) amplifies displacement-norm roundoff past any tight disp tolerance
    # while the residual sits at ~1e-17.
    ops.test("NormUnbalance", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "undrained compression step failed"

    p = np.mean([ops.nodeDisp(n, 4) for n in range(1, 9)])
    # storage bookkeeping: at n ~ n(J), 1/Qbar = n/Kf (infinite grains) and the
    # accumulated coupling is alpha*ln(J) => p ~ -Qbar*ln(lambda) > 0 in
    # compression. n evolves 0.4 -> ~0.33 over the path; use the mid-path n for
    # the analytic target and gate at 5%.
    n_mid = 1.0 - (1.0 - _PORO) / math.sqrt(lam)     # n at J = sqrt(lam)-ish
    p_want = -(_KF / n_mid) * math.log(lam)
    assert p > 0.0, f"undrained compression built negative p = {p:.4e}"
    assert abs(p - p_want) <= 0.05 * p_want, (
        f"undrained large-strain storage: p = {p:.5e}, expected "
        f"~Qbar*|ln(lambda)| = {p_want:.5e} (5% gate)")
