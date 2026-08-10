"""LadrunoUP ``-geom corot`` (ADR 78) — the u-p corotational battery.

Corot on the Biot u-p element = the SolidTransformation seam (ADR 09/10)
extended to the coupled blocks (ADR 78 §3):

  * K  : core solid K globalized through seam 3 with the TOTAL core force
         (∫B'sigma' dV − Q·p), so K_geo carries the total-stress prestress.
  * Q  : the −Q·p leg rides INSIDE that same core force (frame consistency is
         structural); the tangent coupling is −R̄Q and the damp p-row block is
         (R̄Q)^T — its exact transpose.
  * H,S: reference-configuration (exact under pure rotation for skeleton-
         attached permeability; O(strain) otherwise — the corot premise).
  * seepage drive pulled back: drive = R^T(b_f − u_ddot).

Gates (ADR 78 §6), in priority order:

  1. no-regression: '-geom linear' explicit == flag omitted (same code path).
     Parser guards: corot in 2D fatal; -geom finite fatal (material blocker).
  2. STATIC rigid rotation of a saturated hex: sigma' = 0, force = 0, p = 0 —
     the u-row frame gate (an algebraic identity, held to polar precision).
  3. TRANSIENT UNDRAINED rigid rotation: excess p ~ zero. The killer u-p
     gate: the storage-coupling residual is the INCREMENTAL form
     Q^T(u_d − u_d,committed)/dt (ADR 78 §3.3), which is exactly zero under
     rigid motion. A velocity-form coupling would see the chord's apparent
     compression 2(1−cos dth)/dt per step: cumulative eps_v ~ th_tot^2/N
     = 0.082 here, times Qbar ~ 5.5e6 => p ~ 4e5. Newton tolerance is also
     Qbar-amplified (p_noise ~ Qbar*tol_disp ~ 5e-6 at 1e-12), so the assert
     sits at 1e-4 — ~9 orders below the failure mode, ~1.5 above noise.
  4. seepage-drive frame: a block rotated 90 deg under gravity drive reaches
     the hydrostatic p field measured along the CURRENT (rotated) depth —
     catches a missing R^T pull-back (invisible to gates 2-3, which run
     zero-gravity).
  5. Terzaghi/consolidation under corot at small strain reproduces linear
     (corot -> linear as rotations -> 0).
  6. ** drained-equivalence: with k large (fully drained) and hydrostatic p, a
     LadrunoUP column under gravity must reproduce a DRY LadrunoBrick column
     carrying buoyant weight (grad·sigma + gamma_sat*g = 0  <=>
     grad·sigma' + gamma'*g = 0), same mesh, both -geom corot — measured on
     the TIMs footing model at −0.33 %..−0.53 %, so sub-percent is the gate.
     Includes the summed-reaction weight check that catches the -b convention
     trap (LadrunoUP -b is an ACCELERATION *rho_sat; LadrunoBrick -b is a
     FORCE PER VOLUME) instantly.

Run:  py -3.12 -m pytest tests/test_ladruno_up_element_corot.py -x -q
"""
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

pytestmark = [pytest.mark.zone_b]

# --------------------------------------------------------------------------
# shared parameters
# --------------------------------------------------------------------------
_E = 200.0e3          # kPa-ish scale
_NU = 0.3
_RHO_SAT = 2.0        # t/m^3 saturated mixture (material getRho)
_RHO_W = 1.0
_G = 9.81
_KF = 2.2e6           # fluid bulk
_PORO = 0.4

# distorted positive-Jacobian hex (the LadrunoBrick corot battery geometry) —
# R = polar(H) and the spin map are genuinely exercised, not a trivial cube.
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


def _rot(thz, thy):
    cz, sz = math.cos(thz), math.sin(thz)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(thy), math.sin(thy)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return Ry @ Rz


def _up_hex(geom, k=1.0e-3, body=(0.0, 0.0, 0.0), fluid_body=None, stab="off"):
    """One distorted saturated H8 in a fresh 3D/4 model."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in _HEX_NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    args = ["-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
            "-perm", k, k, k, "-stab", stab, "-geom", geom]
    if body != (0.0, 0.0, 0.0):
        args += ["-body", *body]
    if fluid_body is not None:
        args += ["-fluidBody", *fluid_body]
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
#  gate 1 — no-regression + parser guards                                    #
# ========================================================================== #
def _mini_consolidation(geom_args):
    """3-step consolidation on one hex; returns (u, p) trajectory arrays."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in _HEX_NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    ops.element("LadrunoUP", 1, *_HEX_CONN, 1,
                "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                "-perm", 1e-4, 1e-4, 1e-4, *geom_args)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 0)
    for n in (5, 6, 7, 8):
        ops.fix(n, 1, 1, 0, 1)          # drained top, lateral roller
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    out = []
    for _ in range(3):
        assert ops.analyze(1, 0.05) == 0
        for n in range(1, 9):
            out.extend(ops.nodeDisp(n))
    return np.array(out)


def test_geom_linear_flag_is_the_default_path():
    a = _mini_consolidation([])
    b = _mini_consolidation(["-geom", "linear"])
    assert np.array_equal(a, b), "-geom linear must be BIT-identical to the default"


def test_parser_rejects_corot_in_2d():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for tag, (x, y) in {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}.items():
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    with pytest.raises(Exception):
        ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, "-geom", "corot")


def test_parser_rejects_geom_finite():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for tag, (x, y, z) in _HEX_NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    with pytest.raises(Exception):
        ops.element("LadrunoUP", 1, *_HEX_CONN, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, 1e-4, "-geom", "finite")


# ========================================================================== #
#  gate 2 — STATIC rigid rotation: sigma' = 0, f = 0, p = 0                  #
# ========================================================================== #
@pytest.mark.parametrize("thz,thy", [(0.35, 0.25), (1.20, -0.80),
                                     (math.radians(179.0), 0.0)])
def test_static_rigid_rotation_stress_free_and_dry(thz, thy):
    _up_hex("corot")
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


# ========================================================================== #
#  gate 3 — TRANSIENT UNDRAINED rigid rotation: excess p exactly zero        #
# ========================================================================== #
def test_transient_undrained_rigid_rotation_zero_excess_p():
    # 90 deg about z in 30 steps; every CONVERGED state is a rigid rotation
    # (per-DOF Path series carry the exact cos/sin path — a chord-ramped sp
    # would pass through genuinely strained intermediate states).
    nsteps, dt = 30, 0.01
    times = [dt * (k + 1) for k in range(nsteps)]
    _up_hex("corot", k=1.0e-8)          # undrained: negligible drainage, no p BC

    ts = 1
    pat = 1
    for tag, (x, y, z) in _HEX_NODES.items():
        X = np.array([x, y, z])
        for dof in range(3):
            vals = []
            for k in range(nsteps):
                th = (k + 1) * (0.5 * math.pi) / nsteps
                vals.append(float((_rot(th, 0.0) @ X - X)[dof]))
            # -useLast: at the final step the accumulated domain time can
            # overshoot the last Path point by 1 ulp, and PathSeries then
            # returns factor 0 — snapping the block back to u=0 (review F3).
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

    # exact zero in exact arithmetic (the incremental coupling, ADR 78 §3.3);
    # Newton-tolerance noise is Qbar-amplified to ~5e-6; the velocity-form
    # failure mode is ~4e5. 1e-4 discriminates by ~9 orders.
    assert pmax <= 1.0e-4, f"undrained rigid rotation built excess p = {pmax:.3e}"

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    assert np.abs(s).max() <= 1.0e-6 * _E, (
        f"undrained rigid rotation built effective stress {np.abs(s).max():.3e}")


# ========================================================================== #
#  gate 4 — seepage-drive frame: hydrostatic along the ROTATED depth         #
# ========================================================================== #
def test_seepage_drive_follows_the_rotated_frame():
    # rotate the block 90 deg about y (z-axis -> global x), fluid gravity on:
    # the steady p field must be hydrostatic in the CURRENT configuration,
    # i.e. p = p_ref - rho_w*g*(z_cur - z_ref). A missing R^T pull-back
    # leaves it hydrostatic along the REFERENCE z — off by 90 degrees.
    _up_hex("corot", k=1.0e-3,
            body=(0.0, 0.0, 0.0),               # no solid weight: isolate seepage
            fluid_body=(0.0, 0.0, -_G))
    R = _rot(0.0, 0.5 * math.pi)
    ops.fix(1, 0, 0, 0, 1)                       # p reference at node 1: p=0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _HEX_NODES.items():
        d = R @ np.array([x, y, z]) - np.array([x, y, z])
        for i in range(3):
            ops.sp(tag, i + 1, float(d[i]))
    _static(tol=1.0e-10)
    assert ops.analyze(1) == 0

    zc = {t: float((R @ np.array(X))[2]) for t, X in _HEX_NODES.items()}
    p1 = ops.nodeDisp(1, 4)
    for t in range(2, 9):
        p_expect = p1 - _RHO_W * _G * (zc[t] - zc[1])
        p_got = ops.nodeDisp(t, 4)
        assert abs(p_got - p_expect) <= 1.0e-6 * max(1.0, abs(p_expect)), (
            f"node {t}: p = {p_got:.6e}, expected current-config hydrostatic "
            f"{p_expect:.6e} (drive pull-back missing?)")


# ========================================================================== #
#  gate 5 — consolidation column: corot -> linear at small strain            #
# ========================================================================== #
def _column_mesh(nz=8, lz=4.0):
    """1x1xLz column of H8, nodes[(k)] = layer list bottom->top."""
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
        n = [tags[(k, (0, 0))], tags[(k, (1, 0))], tags[(k, (1, 1))], tags[(k, (0, 1))],
             tags[(k + 1, (0, 0))], tags[(k + 1, (1, 0))], tags[(k + 1, (1, 1))],
             tags[(k + 1, (0, 1))]]
        conns.append(n)
    return tags, conns


def _consolidation_run(geom, nz=8, nsteps=25, dt=0.05, q=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _column_mesh(nz=nz)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-5, 1e-5, 1e-5, "-geom", geom)
    top = nz
    for (k, xy), t in tags.items():
        fx = 1
        fy = 1
        fz = 1 if k == 0 else 0
        fp = 1 if k == top else 0        # drained top only
        ops.fix(t, fx, fy, fz, fp)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for xy in ((0, 0), (1, 0), (1, 1), (0, 1)):
        ops.load(tags[(top, xy)], 0.0, 0.0, -q / 4.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-8, 50, 0)   # 1e-10 sits at the roundoff floor here
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        traj.append(ops.nodeDisp(tags[(top, (0, 0))], 3))     # settlement
        traj.append(ops.nodeDisp(tags[(0, (0, 0))], 4))       # base p
    return np.array(traj)


def test_consolidation_corot_reduces_to_linear():
    # ADR 78 §6 gate 4 (refinement-aware): the corot storage coupling is the
    # INCREMENTAL discretization Q^T*du_d/dt (§3.3), a deliberately different
    # scheme from linear's damp·v — so corot vs linear differ at O(dt) even
    # with zero rotation. Gate: the gap shrinks under dt-halving AND sits
    # under a small absolute cap at the run dt.
    lin1 = _consolidation_run("linear", nsteps=20, dt=0.08)
    cor1 = _consolidation_run("corot", nsteps=20, dt=0.08)
    lin2 = _consolidation_run("linear", nsteps=40, dt=0.04)
    cor2 = _consolidation_run("corot", nsteps=40, dt=0.04)
    scale = np.abs(lin1).max()
    # compare at matching physical times: lin2/cor2 sample every 2nd step
    err1 = np.abs(cor1 - lin1).max()
    err2 = np.abs(cor2 - lin2).max()
    assert err1 <= 2.0e-2 * scale, (
        f"corot-vs-linear gap too large at dt=0.08: {err1:.3e} vs scale {scale:.3e}")
    if err1 > 1.0e-8 * scale:      # ratio only meaningful above noise
        assert err2 <= 0.75 * err1, (
            "corot-vs-linear gap did not shrink under dt-halving "
            f"(err(dt)={err1:.3e}, err(dt/2)={err2:.3e}) — the difference is "
            "not the documented O(dt) coupling-scheme term")


# ========================================================================== #
#  gate 6 — ** drained-equivalence vs a dry LadrunoBrick (both corot)        #
# ========================================================================== #
# Footing-on-a-box, the TIMs configuration in miniature. SIDE FACES ARE
# ROLLERS — that is load-bearing for the equivalence: with free sides the UP
# model's total-traction-free surface implies sigma'·n = p != 0 while the dry
# buoyant twin has sigma'·n = 0 (different natural BCs, measured ~3% here).
# With displacement BCs everywhere the discrete identity
# int(grad(Na)·p) = int(Na·rho_w·g) is Gauss-exact on rectangular H8 with the
# (exactly hydrostatic, nodally linear) p field, so the LINEAR leg matches to
# solver tolerance and the corot leg gates sub-percent.
_EQ_N = 4          # 4x4x4 box of H8
_EQ_L = 4.0        # box edge (m); footing = central 2x2-element patch
_EQ_F = 2000.0     # total footing load


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
            fp = 1 if k == n else 0          # drained top (water table at surface)
            ops.fix(t, fx, fy, fz, fp)
        elif fx or fy or fz:
            ops.fix(t, fx, fy, fz)


def _eq_footing_loads(tags, n, ndf):
    # central 2x2-element patch => 3x3 loaded top nodes, bilinear-consistent
    # weights 1/2/4 (corners/edges/center), total 16w = F
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


def _equiv_up(geom):
    """Saturated fully-drained box: gravity (accel conv) + hydrostatic p."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _box_mesh(_EQ_N, _EQ_L)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    for e, conn in enumerate(conns, start=1):
        # LadrunoUP -body is an ACCELERATION (x rho_sat internally) — ADR 78 §5
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1.0, 1.0, 1.0,          # k large: fully drained
                    "-body", 0.0, 0.0, -_G, "-geom", geom)
    _eq_fix_solids(tags, _EQ_N, has_p=True)
    _eq_solve_static()
    assert ops.analyze(5) == 0, "UP gravity stage failed"

    # weight check (the -b convention catcher): total vertical reaction must
    # equal the TOTAL saturated weight rho_sat*g*V (u-rows carry the mixture).
    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in tags.values())
    W = _RHO_SAT * _G * _EQ_L**3
    assert abs(rz - W) <= 1.0e-6 * W, (
        f"UP total vertical reaction {rz:.6f} != saturated weight {W:.6f} "
        "(-b convention broken?)")

    # hydrostatic-p check: base p = rho_w*g*L
    p_base = ops.nodeDisp(tags[(0, 0, 0)], 4)
    assert abs(p_base - _RHO_W * _G * _EQ_L) <= 1.0e-4 * (_RHO_W * _G * _EQ_L)

    # footing push, 5 steps; record the center-settlement backbone
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    _eq_footing_loads(tags, _EQ_N, ndf=4)
    c = _EQ_N // 2
    backbone = []
    for _ in range(5):
        assert ops.analyze(1) == 0, "UP footing stage failed"
        backbone.append(ops.nodeDisp(tags[(c, c, _EQ_N)], 3))
    return np.array(backbone)


def _equiv_brick(geom):
    """The dry twin: buoyant unit weight as FORCE PER VOLUME (-b convention)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    tags, conns = _box_mesh(_EQ_N, _EQ_L)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
    gamma_b = (_RHO_SAT - _RHO_W) * _G
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoBrick", e, *conn, 1, "-formulation", "std",
                    "-geom", geom, "-b", 0.0, 0.0, -gamma_b)
    _eq_fix_solids(tags, _EQ_N, has_p=False)
    _eq_solve_static()
    assert ops.analyze(5) == 0, "brick gravity stage failed"

    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in tags.values())
    Wb = gamma_b * _EQ_L**3
    assert abs(rz - Wb) <= 1.0e-6 * Wb, (
        f"brick total vertical reaction {rz:.6f} != buoyant weight {Wb:.6f}")

    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    _eq_footing_loads(tags, _EQ_N, ndf=3)
    c = _EQ_N // 2
    backbone = []
    for _ in range(5):
        assert ops.analyze(1) == 0, "brick footing stage failed"
        backbone.append(ops.nodeDisp(tags[(c, c, _EQ_N)], 3))
    return np.array(backbone)


# ========================================================================== #
#  gate 7 — DB roundtrip under corot (setDomain rebuilds the committed        #
#  corot state u_d,n from committed nodal displacements — ADR 78 §3.3)        #
# ========================================================================== #
def test_corot_database_roundtrip(tmp_path):
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 4)
        for tag, (x, y, z) in _HEX_NODES.items():
            ops.node(tag, x, y, z)
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO_SAT)
        ops.element("LadrunoUP", 1, *_HEX_CONN, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, 1e-4, "-geom", "corot")
        for n in (1, 2, 3, 4):
            ops.fix(n, 1, 1, 1, 0)
        for n in (5, 6, 7, 8):
            ops.fix(n, 0, 0, 0, 1)          # drained top
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for n in (5, 6, 7, 8):
            ops.load(n, 0.3, 0.0, -1.0, 0.0)   # oblique: nonzero u_d + rotation
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormUnbalance", 1.0e-8, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")

    db = str(tmp_path / "corot_up.db")

    # uninterrupted reference: 5 steps
    build()
    for _ in range(5):
        assert ops.analyze(1, 0.05) == 0
    ref_disp = np.array([ops.nodeDisp(n) for n in range(1, 9)], dtype=object)

    # interrupted twin: 3 steps -> save -> wipe -> restore -> exactness -> +2
    build()
    for _ in range(3):
        assert ops.analyze(1, 0.05) == 0
    mid_p = [ops.nodeDisp(n, 4) for n in range(1, 9)]
    mid_s = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    ops.database("File", db)
    ops.save(1)
    ops.wipe()
    ops.database("File", db)
    ops.restore(1)
    assert max(abs(ops.nodeDisp(n, 4) - mid_p[n - 1]) for n in range(1, 9)) == 0.0, (
        "nodal pressure drift after corot restore")
    s_after = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    assert np.abs(s_after - mid_s).max() < 1.0e-10, "GP stress drift after restore"

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(2):
        assert ops.analyze(1, 0.05) == 0, "post-restore corot continuation failed"

    # continuation matches the uninterrupted run (uDn_ rebuilt by setDomain
    # from committed disps => the incremental coupling differences correctly)
    for n in range(1, 9):
        got = np.array(ops.nodeDisp(n), dtype=float)
        ref = np.array(ref_disp[n - 1], dtype=float)
        assert np.abs(got - ref).max() <= 1.0e-8 * max(1.0, np.abs(ref).max()), (
            f"node {n}: post-restore trajectory diverged from the "
            f"uninterrupted run ({got} vs {ref})")


@pytest.mark.parametrize("geom,tol", [("linear", 1.0e-6), ("corot", 0.7e-2)])
def test_drained_equivalence_vs_dry_brick(geom, tol):
    bb_up = _equiv_up(geom)
    bb_br = _equiv_brick(geom)
    # linear: discrete-exact (Gauss-exact divergence identity) => solver tol.
    # corot: sub-percent (TIMs measured the pair at -0.33%..-0.53% on the real
    # footing model; the small p-field/kinematics cross terms are O(strain)).
    rel = np.abs(bb_up - bb_br) / np.abs(bb_br)
    assert rel.max() <= tol, (
        f"[{geom}] drained-equivalence backbone broke: max rel diff "
        f"{rel.max():.3e} (UP {bb_up[-1]:.6e} vs brick {bb_br[-1]:.6e})")


def test_drained_equivalence_geometric_shift_matches():
    # Review F1: at this load level the whole geometric effect (~6e-4 rel) is
    # far inside the 0.7% corot-leg tolerance, so a corot flag that silently
    # ran the linear path — or a deleted K_geo — would pass the backbone gate.
    # Gate the corot ACTIVATION itself: the corot-minus-linear shift of the UP
    # twin must match the dry brick's (measured agreement 1.4%; assert 25% of
    # the shift, which a no-op corot (shift -> 0) or a broken K_geo misses).
    up_shift = _equiv_up("corot") - _equiv_up("linear")
    br_shift = _equiv_brick("corot") - _equiv_brick("linear")
    scale = np.abs(br_shift).max()
    assert scale > 1.0e-7, (
        f"brick corot-vs-linear shift ~0 ({scale:.3e}) — the probe load no "
        "longer engages geometric nonlinearity; retune _EQ_F")
    assert np.abs(up_shift - br_shift).max() <= 0.25 * scale, (
        "UP corot geometric shift does not track the dry brick twin: "
        f"max diff {np.abs(up_shift - br_shift).max():.3e} vs shift scale "
        f"{scale:.3e} (corot silently inactive on the u-p element?)")
