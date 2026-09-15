"""ADR-92 F10 — why does an IMPL-EX SANISAND push wall on a SELF-WEIGHT deck?

THE REPORTED FINDING (TIMs act, plane-strain strip, same material class and the
same engine as the ADR-95 slab campaign that reached s/B = 0.15):

    B = 1.5 m strip, 15B x 12B box, gamma' = 9.81 kN/m3 as a body force,
    7.65 kPa surcharge outside the footprint, 18.4 kN/m of footing dead load,
    K0 = 0.455 (Jaky) reached through a temporary nu* in a pressure-independent
    elastic stage, then `updateMaterialStage 1`, then a displacement push.
    LadrunoSANISAND, `-implex -implexControl`, `-maxSubsteps 1000`,
    `-Pmin 0.0101`.  The control's refusal counter climbs from the FIRST push
    step (204 120 at step 3 on 9 720 Gauss points), every refusal cuts the step,
    and the harness step falls below its floor at s/B = 0.0125.  Minimum p' in
    the block is 3.8 kPa, so this is NOT the ADR-93 apex/zero-confinement wall.

This driver reproduces that on the fork's own deck and READS the refusals.

WHAT IS REUSED AND WHAT IS NEW
------------------------------
The deck architecture -- plane-strain one-element-thick slab of
`LadrunoBrick -formulation bbar`, `eleLoad -type -selfWeight` against the
element's own `-b`, `confine -> flip -> push`, `LoadControl(-ds)` on an `sp`
pattern under a `Linear` series with the `Transformation` handler, the adaptive
halve/double controller with a pinned subdivision budget, the resultant identity
and the 1-D geostatic patch control -- is `sanisand_tau0_band.py`'s (ADR 90
WP-A2), which is itself R3's (`tests/test_r3_prandtl_collapse_gate.py`).  The
graded-mesh helpers `_graded` / `_n_graded` are R3's, transcribed here rather
than imported because R3 imports `_testbed`, which binds a different engine.

New here:
  * the TIMs GEOMETRY (B = 1.5, 15B x 12B) and LOADS (gamma' = 9.81 buoyant,
    7.65 kPa surcharge outside the footprint, 18.4 kN/m of footing dead load);
  * `nu` as a LEG KNOB -- the "nu* device".  SANISAND's elastic K0 is
    `nu/(1-nu)` and is depth-independent because `nu` is constant, so setting
    `nu* = K0/(1+K0)` at construction is exactly the device TIMs describe.
    ONE DECLARED DIFFERENCE: the fork's `LadrunoSANISAND` takes `nu` as a
    positional constant and exposes no way to put it back after the K0 stage, so
    `nu*` is held for the whole leg.  Legs B and C therefore bracket the device:
    at K0 = 0.455 `nu* = 0.31271` is the material's OWN calibrated `nu = 0.3129`
    to three decimals (so leg B is simultaneously the "K0 reached natively"
    control), while at K0 = 0.818 `nu* = 0.45` is far from it and near-
    incompressible.
  * the PER-GAUSS-POINT REFUSAL CENSUS (`census`): a few push steps taken with
    the control tolerance set so high that nothing can refuse, then every Gauss
    point's `implexDetail` read back beside its p', eta, psi, M^d, depth and
    distance from the footing edge.  That is the measurement the leg table
    cannot make: a leg that refuses tells you THAT it refused, this tells you
    WHICH points would refuse at which tolerance and what state they are in.

CALIBRATION
-----------
Gorini's calibrated `_PARAMS` (ADR-86 sec.5, `tests/test_ladruno_sanisand.py`,
e_init = 0.6944), i.e. the ADR-92 CP1 set -- the fork's own closest set to the
TIMs calibration.  DECLARED DIFFERENCE: TIMs quote e0 = 0.6944 and psi ~ -0.13;
this deck's psi is MEASURED and reported per leg rather than assumed, and the
elastic constants (G0, nu) are the fork's, not TIMs'.

UNITS: kPa, m, Mg.  Stress from `eleResponse ... stress` is TENSION-positive;
p, q, eta below are in SANISAND's own COMPRESSION-positive convention.

USAGE
-----
    python3.12 -u f10_selfweight_wall.py census --leg B --out <dir>
    python3.12 -u f10_selfweight_wall.py leg    --leg B --out <dir> --wall 420
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

# --------------------------------------------------------------------------
# Engine binding.  F10 runs against the PINNED release build of ladruno tip
# 9c2f964; a refusal census with no engine hash attached cannot be
# re-attributed later.  Override with LADRUNO_DIST_BIN / LADRUNO_F10_EXPECT_BUILD
# (the latter accepts `any`, and then says so loudly).
# --------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_BIN = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "..", "release-build-9c2f964", "dist", "bin"))
_BIN = os.environ.get("LADRUNO_DIST_BIN", _DEFAULT_BIN)
EXPECTED_BUILD = os.environ.get("LADRUNO_F10_EXPECT_BUILD", "9c2f964")

if not os.path.isdir(_BIN):
    raise SystemExit(
        f"\nADR-92 F10 driver: no OpenSees build directory at\n    {_BIN}\n"
        f"Set LADRUNO_DIST_BIN to the directory holding the opensees.pyd you "
        f"mean to test.\n")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import numpy as np                                               # noqa: E402
import opensees as ops                                           # noqa: E402

# --------------------------------------------------------------------------
# Material: the ADR-92 CP1 / ADR-86 sec.5 Gorini set, transcribed from
# tests/test_ladruno_sanisand.py::_PARAMS (kPa).  `nu` (index 1) is the ONE
# constant a leg is allowed to move, and only through the nu* device.
# --------------------------------------------------------------------------
PARAMS = [264.32,   # G0
          0.3129,   # nu          <- the nu* knob
          0.6944,   # e_init
          1.33090,  # Mc
          0.71,     # c
          0.027,    # lambda_c
          0.83,     # e0
          0.45,     # ksi
          101.0,    # P_atm
          0.005,    # m
          1.3,      # h0
          0.968,    # ch
          3.5,      # nb
          0.05,     # A0
          5.75,     # nd
          12.5,     # z_max
          1100.0,   # cz
          2.0]      # Rho
M_C, N_B, N_D = PARAMS[3], PARAMS[12], PARAMS[14]
P_ATM = PARAMS[8]

INT_SCHEME, TAN_TYPE, JACO_TYPE = 1, 2, 1     # ModifiedEuler / consistent / 1
TOL_F = TOL_R = 1.0e-7
OPT_PRESIDUAL = 0.0
OPT_PMIN = 1.0e-4 * P_ATM                     # 0.0101 kPa -- the TIMs value

# --------------------------------------------------------------------------
# Geometry (TIMs): B = 1.5 m strip, 15B wide x 12B deep, one element thick.
# --------------------------------------------------------------------------
B_FOOT = 1.5
XLIM = 0.5 * 15.0 * B_FOOT        # 11.25 m  (15 B total width)
ZBOT = -12.0 * B_FOOT             # -18.0 m  (12 B depth)
THICK = 0.5
R_GRADE = 1.35
FINE_DEPTH = 3.0                  # 2 B of h0-sized elements under the footing

GAMMA_EFF = 9.81                  # kN/m3, buoyant (TIMs); dry rho*g would be 19.62
Q_SURCH = 7.65                    # kPa, OUTSIDE the footprint (TIMs)
Q_FOOT_TOT = 18.4                 # kN/m of footing dead load (TIMs)

# --------------------------------------------------------------------------
# The push controller.  DS_* / GROW_AFTER / SUBDIV_BUDGET are R3's pinned
# constants (see that module's CONSTRAINT 2); DS_MAX is scaled to this deck.
# --------------------------------------------------------------------------
DS_BASE, DS_MIN, DS_MAX = 2.0e-5, 2.0e-7, 1.0e-3
GROW_AFTER = 6
SUBDIV_BUDGET = 80
SFRAC = 0.05                      # push target s/B (the reported wall is ~0.012)
N_GRAV = 10                       # stage-0 ramp steps
GRAV_TOL, GRAV_ITER = 1.0e-9, 40
PUSH_TOL_REL = 1.0e-5             # x the total applied vertical load


# --------------------------------------------------------------------------
# R3's graded-mesh helpers, transcribed verbatim (see the module docstring for
# why they are not imported).  Source: tests/test_r3_prandtl_collapse_gate.py.
# --------------------------------------------------------------------------
def _graded(lo, hi, n, h0, fine_at_lo):
    if abs(h0 * n - (hi - lo)) < 1e-12:
        return np.linspace(lo, hi, n + 1)
    a, b = 1.0, 4.0
    for _ in range(200):
        r = 0.5 * (a + b)
        tot = h0 * n if abs(r - 1) < 1e-12 else h0 * (r ** n - 1) / (r - 1)
        a, b = (r, b) if tot < hi - lo else (a, r)
    r = 0.5 * (a + b)
    e = np.concatenate([[0.0], np.cumsum(h0 * r ** np.arange(n))])
    e *= (hi - lo) / e[-1]
    return lo + e if fine_at_lo else hi - e[::-1]


def _n_graded(length, h0, r=R_GRADE):
    return max(1, int(np.ceil(np.log(1.0 + length * (r - 1.0) / h0) / np.log(r))))


def strip_mesh(h0):
    """Graded plane-strain strip mesh, TIMs geometry.  Same construction as
    R3's `_strip_mesh` with B and the box taken from this module."""
    half = 0.5 * B_FOOT
    nf = max(1, int(round(B_FOOT / h0)))
    no = _n_graded(XLIM - half, h0)
    nzo = _n_graded(-ZBOT - FINE_DEPTH, h0)
    x = np.unique(np.concatenate([_graded(-XLIM, -half, no, h0, False),
                                  np.linspace(-half, half, nf + 1),
                                  _graded(half, XLIM, no, h0, True)]))
    y = np.array([0.0, THICK])
    z = np.unique(np.concatenate([
        _graded(ZBOT, -FINE_DEPTH, nzo, h0, False),
        np.linspace(-FINE_DEPTH, 0.0, int(round(FINE_DEPTH / h0)) + 1)]))
    nx, ny, nz = len(x), len(y), len(z)
    idx = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    nodes = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1).reshape(-1, 3)
    hexes = np.array([[idx[i, j, k], idx[i+1, j, k], idx[i+1, j+1, k], idx[i, j+1, k],
                       idx[i, j, k+1], idx[i+1, j, k+1], idx[i+1, j+1, k+1],
                       idx[i, j+1, k+1]]
                      for i in range(nx - 1) for j in range(ny - 1)
                      for k in range(nz - 1)], dtype=np.int32)
    p = nodes[hexes]
    vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                    np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (vol > 0).all(), "inverted hexes in the strip mesh"
    cen = p.mean(axis=1)

    tol = 1e-9
    top = np.abs(nodes[:, 2]) < tol
    sets = dict(top=np.where(top)[0],
                bottom=np.where(np.abs(nodes[:, 2] - ZBOT) < tol)[0],
                xface=np.where(np.abs(np.abs(nodes[:, 0]) - XLIM) < tol)[0],
                footing=np.where(top & (np.abs(nodes[:, 0]) <= half + tol))[0])
    assert len(sets["footing"]) == 2 * (nf + 1), len(sets["footing"])

    trib = np.zeros(len(nodes))
    for h in hexes:
        face = h[[4, 5, 6, 7]]
        q = nodes[face]
        if np.any(np.abs(q[:, 2]) > tol):
            continue
        d1, d2 = q[2, :2] - q[0, :2], q[3, :2] - q[1, :2]
        trib[face] += 0.5 * abs(d1[0] * d2[1] - d1[1] * d2[0]) / 4.0
    assert abs(trib.sum() - 2 * XLIM * THICK) < 1e-6, trib.sum()
    return nodes, hexes, trib, sets, vol, cen


# --------------------------------------------------------------------------
# state helpers
# --------------------------------------------------------------------------
def dev_and_p(st):
    """(p, q, eta) COMPRESSION-positive from a TENSION-positive stress vector."""
    p = -(st[0] + st[1] + st[2]) / 3.0
    s = (st[0] + p, st[1] + p, st[2] + p, st[3], st[4], st[5])
    j2 = 0.5 * (s[0] ** 2 + s[1] ** 2 + s[2] ** 2) + s[3] ** 2 + s[4] ** 2 + s[5] ** 2
    q = math.sqrt(3.0 * j2)
    return p, q, (q / p if p > 1.0e-9 else float("nan"))


def m_bd(psi):
    """SANISAND's bounding and dilatancy stress ratios in triaxial compression."""
    return M_C * math.exp(-N_B * psi), M_C * math.exp(N_D * psi)


def resp(e, gp, name, n):
    try:
        v = ops.eleResponse(e, "material", gp, name)
    except Exception:
        return None
    return v if (v is not None and len(v) >= n) else None


def refusals():
    v = resp(1, 1, "implexRefusals", 4)
    if v is None:
        return dict(total=0, sign=0, control=0, companion=0)
    return dict(total=int(round(v[0])), sign=int(round(v[1])),
                control=int(round(v[2])), companion=int(round(v[3])))


def guards():
    v = resp(1, 1, "implexGuards", 7)
    if v is None:
        return [0] * 7
    return [int(round(x)) for x in v[:7]]


def k0_to_nu(k0):
    return k0 / (1.0 + k0)


# --------------------------------------------------------------------------
# the legs
# --------------------------------------------------------------------------
LEGS = {
    "A":  ("weightless + uniform surcharge (the ADR-95 campaign condition), "
           "implex + control", dict(gamma=0.0, q_uniform=10.0, q_surch=0.0,
                                    q_foot=0.0, k0=0.455)),
    "B":  ("self-weight, K0 = 0.455 via nu*, implex + control 0.05 (the "
           "reported configuration)", dict(k0=0.455)),
    "C":  ("self-weight, K0 = 0.818 via nu* = 0.45 (the ADR-95 slab's own "
           "Poisson ratio), implex + control 0.05", dict(k0=0.818)),
    "D":  ("B with -implexFactor controlIter", dict(k0=0.455, factor="controlIter")),
    "E":  ("B with a hold at the flip + a tiny first push step, geometric growth",
           dict(k0=0.455, hold=True, ds0=2.0e-6, grow=1.5)),
    "F1": ("B at the C++ default tolerance 0.1", dict(k0=0.455, tol=0.1)),
    "F2": ("B at tol 0.05 with the reduction floor RAISED to 0.5 (the floor "
           "binds early, so refusals become floor fallbacks)",
           dict(k0=0.455, redlim=0.5)),
    "F3": ("B at tol 0.5 -- how far the tolerance alone carries a self-weight "
           "deck", dict(k0=0.455, tol=0.5)),
    "G":  ("B implicit (no -implex) -- the reference reach",
           dict(k0=0.455, implex=False)),
    "H":  ("B with a heavy 100 kPa uniform surcharge: min p' raised two decades, "
           "self-weight and the dilatant K0 state kept",
           dict(k0=0.455, q_uniform=100.0)),
    "I":  ("B with -implexGuard off -- the DECISIVE leg for the P2-2/control "
           "interaction: stop forcing f = 0 on a reversing/softening predecessor "
           "and see whether the control still refuses",
           dict(k0=0.455, guard="off")),
    "J":  ("B with -implexGuard off AND -implexTrialGuard off -- neither f = 0 "
           "path armed", dict(k0=0.455, guard="off", trialguard="off")),
    "K":  ("B with the controller's growth factor pinned at 1.0 (ds never "
           "doubles, so the clock ratio f can never exceed 1)",
           dict(k0=0.455, grow=1.0)),
    "L":  ("B at a CONSTANT ds = 8e-5 m -- the largest step the s/B = 0.0085 "
           "refinement probe measured at 4.4x under tol (max err 0.011)",
           dict(k0=0.455, grow=1.0, ds0=8.0e-5)),
    "M":  ("B with a GENTLE growth factor 1.25 instead of 2, so the clock ratio "
           "f on a growth step is 1.25 and not 2",
           dict(k0=0.455, grow=1.25)),
    "N":  ("B with -implexControl REMOVED -- bare `-implex`, the way the ADR-95 "
           "reference campaign ran it; growth factor still 2.0, everything else "
           "identical. THE DECISIVE ARM.", dict(k0=0.455, control=False)),
    "N1": ("N with the growth factor also pinned at 1.0 -- separates 'no control' "
           "from 'no growth'", dict(k0=0.455, control=False, grow=1.0)),
}


def build(h0, k0, gamma, q_uniform, q_surch, q_foot, implex, tol, redlim,
          factor, maxsubsteps, guard="on", trialguard="on", control=True,
          verbose=True):
    """Build + confine + flip.  Returns the deck dict."""
    nodes, hexes, trib, sets, vol, cen = strip_mesh(h0)
    n_hex, n_nodes = len(hexes), len(nodes)
    nu = k0_to_nu(k0)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))

    par = list(PARAMS)
    par[1] = nu                                   # <- the nu* device
    ops.nDMaterial(
        "LadrunoSANISAND", 1, *par,
        INT_SCHEME, TAN_TYPE, JACO_TYPE, TOL_F, TOL_R,
        "-Presidual", OPT_PRESIDUAL, "-Pmin", OPT_PMIN, "-honorTolR", 0,
        "-maxSubsteps", int(maxsubsteps),
        *(("-implex",) if implex else ()),
        # ADR-92 F10 review round 1, BLOCKER 1: `-implexControl` used to be
        # hard-wired onto `-implex` here, so no leg could run IMPL-EX the way
        # the ADR-95 reference campaign actually ran it (`sanisand_path_diag.py`
        # passes `-implex` alone).  It is now a leg knob, and leg N is that arm.
        *(("-implexControl", float(tol), float(redlim))
          if (implex and control) else ()),
        *(("-implexFactor", factor)
          if (implex and control and factor != "fixed") else ()),
        *(("-implexGuard", guard) if (implex and guard != "on") else ()),
        *(("-implexTrialGuard", trialguard)
          if (implex and trialguard != "on") else ()))

    for e, conn in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                    "-geom", "linear", "-b", 0.0, 0.0, -gamma,
                    "-formulation", "bbar")

    # ONE fix() per node.  u_y = 0 everywhere == plane strain.  ROUGH footing
    # (u_x = 0 on its nodes) -- the ADR-92 CP1 convention, and the rigid-footing
    # kinematics TIMs get from LadrunoKinematicCoupling.
    bt, xf = set(sets["bottom"].tolist()), set(sets["xface"].tolist())
    ft = set(sets["footing"].tolist())
    for n in range(n_nodes):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if (n in xf or n in ft) else 0, 1, 0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso", "-matrixType", 0)   # unsymmetric consistent tangent
        solver = "Pardiso"
    except Exception:
        ops.system("UmfPack")
        solver = "UmfPack"
    ops.test("NormDispIncr", GRAV_TOL, GRAV_ITER, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 0)

    # ---- stage 0a: self weight ------------------------------------------
    patch_err = float("nan")
    want = 0.0
    if gamma > 0.0:
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", *range(1, n_hex + 1), "-type", "-selfWeight",
                    0.0, 0.0, 1.0)
        ops.integrator("LoadControl", 1.0 / N_GRAV)
        for k in range(N_GRAV):
            assert ops.analyze(1) == 0, f"gravity step {k + 1} failed"
        ops.reactions()
        rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
        want = gamma * float(vol.sum())
        assert abs(rz / want - 1.0) < 1e-6, (
            f"gravity resultant identity violated: {rz} vs {want} kN")
        # The 1-D geostatic patch -- exact because nu is constant, so K0 is
        # depth-independent even though SANISAND's moduli go as sqrt(p).
        # b-bar makes the stress piecewise constant at the element CENTROID.
        ezz = exx = 0.0
        for e in range(1, n_hex + 1):
            szz_ex = gamma * float(cen[e - 1, 2])
            sxx_ex = k0 * szz_ex
            st = resp(e, 1, "stress", 6)
            if st is None:
                continue
            ezz = max(ezz, abs(st[2] - szz_ex) / abs(szz_ex))
            exx = max(exx, abs(st[0] - sxx_ex) / abs(sxx_ex))
        patch_err = max(ezz, exx)
        assert patch_err < 1e-6, (
            f"1-D geostatic patch failed at {patch_err:.3e} -- the body-force "
            "or K0 convention is wrong and nothing downstream is valid")
        ops.loadConst("-time", 0.0)

    # ---- stage 0b: surface loads ----------------------------------------
    # q_uniform on the WHOLE top face, q_surch OUTSIDE the footprint only, and
    # the footing's own dead load spread over the footprint.  Applied after the
    # geostatic controls (which assume sigma_zz = gamma z) and before the flip.
    half = 0.5 * B_FOOT
    q_foot_kpa = (q_foot / B_FOOT) if q_foot > 0.0 else 0.0
    surf = 0.0
    if q_uniform > 0.0 or q_surch > 0.0 or q_foot_kpa > 0.0:
        ops.timeSeries("Linear", 3)
        ops.pattern("Plain", 3, 3)
        for n in sets["top"]:
            if trib[n] <= 0:
                continue
            inside = abs(nodes[n, 0]) <= half + 1e-9
            q = q_uniform + (q_foot_kpa if inside else q_surch)
            if q == 0.0:
                continue
            ops.load(int(n) + 1, 0.0, 0.0, -q * float(trib[n]))
            surf += q * float(trib[n])
        ops.integrator("LoadControl", 1.0 / N_GRAV)
        for k in range(N_GRAV):
            assert ops.analyze(1) == 0, f"surface-load step {k + 1} failed"
        ops.reactions()
        rz2 = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
        assert abs(rz2 / (want + surf) - 1.0) < 1e-6, (
            f"surface-load resultant identity violated: {rz2} vs {want + surf}")

    # ---- the at-rest census ---------------------------------------------
    origin = gp_census(n_hex, cen, half, at_origin=True)
    eta_max = max((r["eta"] for r in origin if r["eta"] == r["eta"]), default=0.0)
    assert eta_max < M_C, (
        f"max eta = {eta_max:.5f} at the stage flip is at or above M_c = {M_C} "
        "-- Elastic2Plastic will inflate the calibrated friction constant")

    ops.updateMaterialStage("-material", 1, "-stage", 1)

    if verbose:
        print(f"    mesh {n_hex} hexes / {n_nodes} nodes / {3*n_nodes} DOF, "
              f"solver {solver}, nu* = {nu:.6f} (K0 = {k0:.4f})", flush=True)
        print(f"    geostatic patch max rel err = {patch_err:.3e}, "
              f"eta_max/M_c at the flip = {eta_max / M_C:.4f}", flush=True)
    return dict(nodes=nodes, hexes=hexes, trib=trib, sets=sets, cen=cen,
                n_hex=n_hex, half=half, solver=solver, nu=nu,
                patch_err=patch_err, eta_max=eta_max, origin=origin,
                applied=want + surf)


def gp_census(n_hex, cen, half, at_origin=False):
    """Per-Gauss-point state + (after a push step) the IMPL-EX error.

    b-bar makes p element-constant, so the element CENTROID is the exact
    locator for p'; eta and the error still vary Gauss point to Gauss point.
    """
    out = []
    for e in range(1, n_hex + 1):
        xc, zc = float(cen[e - 1, 0]), float(cen[e - 1, 2])
        dx = max(0.0, abs(xc) - half)      # horizontal distance from the edge
        for gp in range(1, 9):
            st = resp(e, gp, "stress", 6)
            if st is None:
                continue
            p, q, eta = dev_and_p(st)
            ps = resp(e, gp, "psi", 1)
            psi = float(ps[0]) if ps else float("nan")
            mb, md = m_bd(psi) if psi == psi else (float("nan"), float("nan"))
            row = dict(ele=e, gp=gp, x=xc, z=zc, dx=dx, p=p, q=q, eta=eta,
                       psi=psi, Mb=mb, Md=md)
            if not at_origin:
                d = resp(e, gp, "implexDetail", 6)
                if d is not None:
                    row.update(err=float(d[0]), err_dev=float(d[1]),
                               err_vol=float(d[2]), clamp=float(d[3]),
                               f=float(d[5]))
            out.append(row)
    return out


# --------------------------------------------------------------------------
# the push
# --------------------------------------------------------------------------
def push_setup(deck):
    foot = [int(n) + 1 for n in deck["sets"]["footing"]]
    ops.reactions()
    r0 = sum(ops.nodeReaction(t, 3) for t in foot)
    uz0 = ops.nodeDisp(foot[0], 3)
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso", "-matrixType", 0)
    except Exception:
        ops.system("UmfPack")
    ops.analysis("Static")
    return foot, uz0, r0


def q_now(foot, r0):
    ops.reactions()
    return (-(sum(ops.nodeReaction(t, 3) for t in foot) - r0)
            / (B_FOOT * THICK))


def run_leg(args, name, desc, kw):
    tol = kw.get("tol", args.tol)
    redlim = kw.get("redlim", args.redlim)
    implex = kw.get("implex", True)
    factor = kw.get("factor", "fixed")
    gamma = kw.get("gamma", GAMMA_EFF)
    q_uniform = kw.get("q_uniform", 0.0)
    q_surch = kw.get("q_surch", Q_SURCH)
    q_foot = kw.get("q_foot", Q_FOOT_TOT)
    k0 = kw.get("k0", 0.455)
    hold = kw.get("hold", False)
    ds0 = kw.get("ds0", DS_BASE)
    grow = kw.get("grow", 2.0)
    control = kw.get("control", True)

    print(f"=== leg {name}: {desc}", flush=True)
    t0 = time.time()
    deck = build(args.h0, k0, gamma, q_uniform, q_surch, q_foot, implex,
                 tol, redlim, factor, args.maxsubsteps,
                 kw.get("guard", "on"), kw.get("trialguard", "on"), control)
    foot, uz0, r0 = push_setup(deck)
    ptol = PUSH_TOL_REL * max(deck["applied"], 1.0)
    ladder = [("Newton", ptol, 25, 0),
              ("NewtonLineSearch", ptol, 40, 0),
              ("KrylovNewton", 10.0 * ptol, 60, 1)]

    if hold:      # a zero-increment re-equilibration at the flip (P2-3 rule)
        ops.integrator("LoadControl", 0.0)
        ops.test("NormUnbalance", ptol, 25, 0)
        ops.algorithm("Newton")
        rc = ops.analyze(1)
        print(f"    hold at the flip: analyze(1) rc = {rc}, guards = {guards()}",
              flush=True)

    out = os.path.join(args.out, f"f10_{name}.csv")
    # ADR-92 F10 review round 1, nit 13.  `hypo_bearing/README.md` records this
    # exact failure: two processes on the same leg interleave rows into one CSV
    # and leave a torn line, and the summary then reduces a file neither process
    # wrote.  It cost legs L and M a re-run here.  Same guard as the ADR-79
    # runner: refuse a CSV another process touched in the last 180 s
    # (F10_FORCE=1 overrides).
    if (os.path.exists(out) and not os.environ.get("F10_FORCE")
            and time.time() - os.path.getmtime(out) < 180.0):
        raise SystemExit(
            f"\n{out} was written {time.time()-os.path.getmtime(out):.0f} s "
            f"ago -- another process is probably running leg {name}. Two "
            f"writers interleave rows and tear the file. Set F10_FORCE=1 to "
            f"override.\n")
    fh = open(out, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["step", "s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed",
                 "ref_total", "ref_control", "ref_companion", "ref_sign",
                 "g_floor", "g_guard", "g_hold", "g_rev", "g_trial", "g_flip",
                 "g_backoff", "wall_s"])

    smax = SFRAC * B_FOOT
    ds, good, nsub, nfail, nstep = ds0, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    ref_prev = refusals()
    rows = []
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > args.wall:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now/B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(algo)
            if ops.analyze(1) == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good, nsub = 0, nsub + 1
            ds *= 0.5
            if nsub > SUBDIV_BUDGET:
                mode = "BUDGET"
                verdict = f"subdivision budget spent at s/B = {s_now/B_FOOT:.5f}"
                break
            if ds < DS_MIN:
                mode = "FLOOR"
                verdict = (f"step collapsed to the {DS_MIN*1e3:.4g} mm floor at "
                           f"s/B = {s_now/B_FOOT:.5f}")
                break
            continue
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:
            ds, good = min(grow * ds, DS_MAX), 0
        nstep += 1
        s = uz0 - ops.getTime()
        q = q_now(foot, r0)
        r = refusals()
        g = guards()
        rows.append((s / B_FOOT, q))
        wr.writerow([nstep, f"{s:.9g}", f"{s/B_FOOT:.9g}", f"{q:.6g}",
                     f"{ds*1e3:.6g}", relaxed, r["total"], r["control"],
                     r["companion"], r["sign"], *g, f"{time.time()-t0:.1f}"])
        fh.flush()
        if nstep <= 5 or nstep % 25 == 0:
            print(f"    step {nstep:4d}  s/B {s/B_FOOT:.6f}  q {q:8.3f}  "
                  f"ds {ds*1e3:.4g} mm  refusals total {r['total']} "
                  f"(ctl {r['control']}, comp {r['companion']}, sign {r['sign']}; "
                  f"+{r['total']-ref_prev['total']} this step)  "
                  f"guards floor/guard/trial = {g[0]}/{g[1]}/{g[4]}", flush=True)
        ref_prev = r
    fh.close()
    wall = time.time() - t0
    r = refusals()
    g = guards()
    s_end = rows[-1][0] if rows else 0.0
    q_end = rows[-1][1] if rows else 0.0
    q_max = max((x[1] for x in rows), default=0.0)
    print(f"    MODE = {mode}  [{verdict}]  wall = {wall:.1f}s  steps = {nstep} "
          f"nsub = {nsub}  nfail = {nfail}", flush=True)
    print(f"    s/B reached = {s_end:.6f}   q_end = {q_end:.3f}  "
          f"q_max = {q_max:.3f} kPa", flush=True)
    print(f"    refusals: total {r['total']}  control {r['control']}  "
          f"companion {r['companion']}  signChange {r['sign']}", flush=True)
    print(f"    guards: floor {g[0]}  f=0 {g[1]}  hold {g[2]}  rev {g[3]}  "
          f"trialF0 {g[4]}  holdSkip {g[5]}  backoff {g[6]}", flush=True)
    meta = dict(leg=name, desc=desc, build=ops.ladrunoBuild(), h0=args.h0,
                k0=k0, nu=deck["nu"], gamma=gamma, q_uniform=q_uniform,
                q_surch=q_surch, q_foot=q_foot, implex=implex, tol=tol,
                redlim=redlim, factor=factor, maxsubsteps=args.maxsubsteps,
                hold=hold, ds0=ds0, grow=grow, control=control,
                guard=kw.get("guard", "on"),
                trialguard=kw.get("trialguard", "on"),
                mode=mode, verdict=verdict,
                wall=wall, steps=nstep, nsub=nsub, nfail=nfail,
                s_over_B=s_end, q_end=q_end, q_max=q_max, refusals=r,
                guards=g, patch_err=deck["patch_err"],
                eta_max_over_Mc=deck["eta_max"] / M_C,
                n_hex=deck["n_hex"], n_gp=8 * deck["n_hex"])
    with open(os.path.join(args.out, f"f10_{name}.json"), "w") as f:
        json.dump(meta, f, indent=1, default=float)
    return meta


# --------------------------------------------------------------------------
# the census
# --------------------------------------------------------------------------
def run_census(args, name, desc, kw):
    """Take the first few push steps with the control tolerance set so high
    that NOTHING can refuse, then read every Gauss point's extrapolation error
    beside its state.  This is the measurement a walled leg cannot make."""
    gamma = kw.get("gamma", GAMMA_EFF)
    q_uniform = kw.get("q_uniform", 0.0)
    q_surch = kw.get("q_surch", Q_SURCH)
    q_foot = kw.get("q_foot", Q_FOOT_TOT)
    k0 = kw.get("k0", 0.455)

    print(f"=== census {name}: {desc}", flush=True)
    t0 = time.time()
    deck = build(args.h0, k0, gamma, q_uniform, q_surch, q_foot, True,
                 args.census_tol, 1.0e-6, "fixed", args.maxsubsteps,
                 kw.get("guard", "on"), kw.get("trialguard", "on"))

    org = deck["origin"]
    dil = [r for r in org if r["eta"] == r["eta"] and r["eta"] > r["Md"]]
    pmin = min(r["p"] for r in org)
    print(f"    ORIGIN: {len(org)} GPs, p' in [{pmin:.4g}, "
          f"{max(r['p'] for r in org):.4g}] kPa; psi in "
          f"[{min(r['psi'] for r in org):.4f}, {max(r['psi'] for r in org):.4f}]",
          flush=True)
    print(f"    ORIGIN: eta/M^d > 1 (DILATANT at rest) at {len(dil)}/{len(org)} "
          f"= {100.0*len(dil)/len(org):.2f} % of Gauss points", flush=True)
    print(f"    ORIGIN: median eta = "
          f"{np.median([r['eta'] for r in org]):.4f}, median M^d = "
          f"{np.median([r['Md'] for r in org]):.4f}, median M^b = "
          f"{np.median([r['Mb'] for r in org]):.4f}", flush=True)

    foot, uz0, r0 = push_setup(deck)
    ptol = PUSH_TOL_REL * max(deck["applied"], 1.0)
    sfx = f"{name}_ds{args.census_ds:g}"
    rows_out, have = [], []
    for step in range(1, args.census_steps + 1):
        ops.integrator("LoadControl", -args.census_ds)
        ops.test("NormUnbalance", ptol, 60, 0)
        ops.algorithm("Newton")
        rc = ops.analyze(1)
        if rc != 0:
            ops.test("NormUnbalance", 10 * ptol, 80, 0)
            ops.algorithm("KrylovNewton")
            rc = ops.analyze(1)
        s = uz0 - ops.getTime()
        if rc != 0:
            # A failed step reverts the domain, so every `implexDetail` reads the
            # PREVIOUS commit (or zero).  Censusing it would report the guard,
            # not the material: stop instead, and say why.
            print(f"  -- step {step}: rc={rc} -- the analysis did NOT converge at "
                  f"ds = {args.census_ds:g} m; the census stops here rather than "
                  f"report a reverted state.", flush=True)
            break
        c = gp_census(deck["n_hex"], deck["cen"], deck["half"])
        have = [r for r in c if "err" in r]
        errs = np.array([r["err"] for r in have]) if have else np.zeros(0)
        r_led, g_led = refusals(), guards()
        print(f"  -- step {step}: rc={rc}  s/B = {s/B_FOOT:.6f}  "
              f"q = {q_now(foot, r0):.3f} kPa  GPs with an error = {len(have)}",
              flush=True)
        if len(errs):
            for t in (0.05, 0.1, 0.2, 0.5):
                n = int((errs > t).sum())
                print(f"       error > {t:<5}: {n:6d} GPs "
                      f"({100.0*n/len(errs):6.2f} %)", flush=True)
            print(f"       error: median {np.median(errs):.4g}  p90 "
                  f"{np.percentile(errs,90):.4g}  max {errs.max():.4g}", flush=True)
            print(f"       ledger: refusals {r_led}  guards {g_led}", flush=True)
            bad = sorted(have, key=lambda r: -r["err"])[:12]
            print("       worst 12 GPs:  ele/gp     z      dx       p'    eta  "
                  "eta/M^d      psi      err      dev      vol       f", flush=True)
            for r in bad:
                print(f"         {r['ele']:5d}/{r['gp']} {r['z']:7.3f} "
                      f"{r['dx']:7.3f} {r['p']:8.3f} {r['eta']:6.3f} "
                      f"{r['eta']/r['Md']:7.3f} {r['psi']:8.4f} "
                      f"{r['err']:8.4f} {r['err_dev']:8.4f} {r['err_vol']:8.4f} "
                      f"{r['f']:7.4f}", flush=True)
            over = [r for r in have if r["err"] > 0.05]
            if over:
                zz = np.array([r["z"] for r in over])
                dd = np.array([r["dx"] for r in over])
                pp = np.array([r["p"] for r in over])
                ee = np.array([r["eta"] / r["Md"] for r in over])
                print(f"       >0.05 population: depth |z| "
                      f"[{-zz.max():.2f}, {-zz.min():.2f}] m (median "
                      f"{-np.median(zz):.2f}); dx from the footing edge "
                      f"[{dd.min():.2f}, {dd.max():.2f}] m (median "
                      f"{np.median(dd):.2f}); p' [{pp.min():.3f}, {pp.max():.3f}]"
                      f" kPa (median {np.median(pp):.3f}); eta/M^d "
                      f"[{ee.min():.3f}, {ee.max():.3f}] (median "
                      f"{np.median(ee):.3f})", flush=True)
                print(f"       >0.05 population: "
                      f"{100*float((dd > 2*B_FOOT).mean()):.1f} % sit more than "
                      f"2B from the footing edge, "
                      f"{100*float((-zz > 2*B_FOOT).mean()):.1f} % deeper than 2B",
                      flush=True)
        rows_out.append(dict(step=step, rc=rc, s_over_B=s / B_FOOT, n=len(errs),
                             n_over_005=int((errs > 0.05).sum()) if len(errs) else 0,
                             n_over_01=int((errs > 0.1).sum()) if len(errs) else 0,
                             med=float(np.median(errs)) if len(errs) else 0.0,
                             mx=float(errs.max()) if len(errs) else 0.0))
    if have:
        path = os.path.join(args.out, f"f10_census_{sfx}.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ele", "gp", "x", "z", "dx", "p", "q", "eta", "psi",
                        "Mb", "Md", "eta_over_Md", "err", "err_dev",
                        "err_vol", "f"])
            for r in have:
                w.writerow([r["ele"], r["gp"], f"{r['x']:.6g}", f"{r['z']:.6g}",
                            f"{r['dx']:.6g}", f"{r['p']:.6g}", f"{r['q']:.6g}",
                            f"{r['eta']:.6g}", f"{r['psi']:.6g}",
                            f"{r['Mb']:.6g}", f"{r['Md']:.6g}",
                            f"{r['eta']/r['Md']:.6g}", f"{r['err']:.6g}",
                            f"{r['err_dev']:.6g}", f"{r['err_vol']:.6g}",
                            f"{r['f']:.6g}"])
        print(f"    wrote {path}", flush=True)
    meta = dict(leg=name, kind="census", build=ops.ladrunoBuild(), h0=args.h0,
                k0=k0, nu=deck["nu"], gamma=gamma, q_uniform=q_uniform,
                census_tol=args.census_tol, census_ds=args.census_ds,
                n_gp=len(org), p_min=pmin,
                frac_dilatant_at_rest=len(dil) / len(org),
                eta_max_over_Mc=deck["eta_max"] / M_C,
                patch_err=deck["patch_err"], steps=rows_out,
                wall=time.time() - t0)
    with open(os.path.join(args.out, f"f10_census_{sfx}.json"), "w") as f:
        json.dump(meta, f, indent=1, default=float)
    return meta



# --------------------------------------------------------------------------
# the step-size refinement probe
# --------------------------------------------------------------------------
def run_probe(args, name, desc, kw):
    """Walk to a FIXED settlement on a refusal-free path (constant ds, control
    tolerance set so nothing can refuse), then take ONE step of size
    `--census-ds` and census the error.

    THE QUESTION: at the state where leg B seizes, does the IMPL-EX error at a
    Gauss point SCALE with the step?  If it does, the wall is the controller
    walking into an absolute tolerance and a smaller step is always a way out.
    If it does not, the companion is making a jump independent of the increment
    -- the un-primed-step signature (`LadrunoSANISAND.cpp:2895-2900`) reached on
    a PRIMED step, where no exemption applies and no subdivision can help.

    Every ds is a SEPARATE process run to the same settlement on the same
    deterministic path, so the probe steps all start from the same committed
    state.
    """
    gamma = kw.get("gamma", GAMMA_EFF)
    q_uniform = kw.get("q_uniform", 0.0)
    q_surch = kw.get("q_surch", Q_SURCH)
    q_foot = kw.get("q_foot", Q_FOOT_TOT)
    k0 = kw.get("k0", 0.455)

    print(f"=== probe {name} to s/B = {args.probe_s}, then ONE step at "
          f"ds = {args.census_ds:g} m", flush=True)
    t0 = time.time()
    deck = build(args.h0, k0, gamma, q_uniform, q_surch, q_foot, True,
                 args.census_tol, 1.0e-6, "fixed", args.maxsubsteps,
                 kw.get("guard", "on"), kw.get("trialguard", "on"))
    foot, uz0, r0 = push_setup(deck)
    ptol = PUSH_TOL_REL * max(deck["applied"], 1.0)
    target = args.probe_s * B_FOOT
    n = 0
    while (uz0 - ops.getTime()) < target - 1e-12:
        if time.time() - t0 > args.wall:
            print(f"    WALL before the probe settlement "
                  f"(s/B = {(uz0-ops.getTime())/B_FOOT:.6f})", flush=True)
            return None
        ops.integrator("LoadControl", -DS_BASE)
        ops.test("NormUnbalance", ptol, 25, 0)
        ops.algorithm("Newton")
        if ops.analyze(1) != 0:
            ops.test("NormUnbalance", 10 * ptol, 60, 0)
            ops.algorithm("KrylovNewton")
            if ops.analyze(1) != 0:
                print(f"    walk-up FAILED at s/B = "
                      f"{(uz0-ops.getTime())/B_FOOT:.6f}", flush=True)
                return None
        n += 1
    s_walk = uz0 - ops.getTime()
    print(f"    walked {n} steps of {DS_BASE:g} m to s/B = {s_walk/B_FOOT:.6f} "
          f"(q = {q_now(foot, r0):.3f} kPa) in {time.time()-t0:.0f} s", flush=True)

    ops.integrator("LoadControl", -args.census_ds)
    ops.test("NormUnbalance", ptol, 60, 0)
    ops.algorithm("Newton")
    rc = ops.analyze(1)
    if rc != 0:
        ops.test("NormUnbalance", 10 * ptol, 80, 0)
        ops.algorithm("KrylovNewton")
        rc = ops.analyze(1)
    c = gp_census(deck["n_hex"], deck["cen"], deck["half"])
    have = [r for r in c if "err" in r]
    errs = np.array([r["err"] for r in have]) if have else np.zeros(0)
    print(f"    PROBE rc={rc}  ds={args.census_ds:g}  n={len(errs)}  "
          f"median {np.median(errs):.5g}  p99 {np.percentile(errs,99):.5g}  "
          f"max {errs.max():.5g}  n>0.05 {int((errs>0.05).sum())}  "
          f"n>0.1 {int((errs>0.1).sum())}", flush=True)
    bad = sorted(have, key=lambda r: -r["err"])[:6]
    for r in bad:
        print(f"      worst {r['ele']:5d}/{r['gp']} z {r['z']:7.3f} dx "
              f"{r['dx']:7.3f} p' {r['p']:8.3f} eta {r['eta']:6.3f} "
              f"eta/M^b {r['eta']/r['Mb']:6.3f} eta/M^d {r['eta']/r['Md']:6.3f} "
              f"psi {r['psi']:8.4f} err {r['err']:8.5f} dev {r['err_dev']:8.5f} "
              f"vol {r['err_vol']:8.5f} f {r['f']:7.4f}", flush=True)
    meta = dict(leg=name, kind="probe", build=ops.ladrunoBuild(),
                probe_s=args.probe_s, s_walked=s_walk / B_FOOT, rc=rc,
                ds=args.census_ds, n=len(errs),
                median=float(np.median(errs)), max=float(errs.max()),
                n_over_005=int((errs > 0.05).sum()),
                n_over_01=int((errs > 0.1).sum()),
                worst=[{k: float(v) if k not in ("ele", "gp") else int(v)
                        for k, v in r.items()} for r in bad],
                wall=time.time() - t0)
    with open(os.path.join(
            args.out,
            f"f10_probe_{name}_s{args.probe_s:g}_ds{args.census_ds:g}.json"),
            "w") as f:
        json.dump(meta, f, indent=1, default=float)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["leg", "census", "probe"])
    ap.add_argument("--leg", required=True)
    ap.add_argument("--out", default=_HERE)
    ap.add_argument("--h0", type=float, default=0.5)
    ap.add_argument("--wall", type=float, default=420.0)
    ap.add_argument("--tol", type=float, default=0.05)
    ap.add_argument("--redlim", type=float, default=0.01)
    ap.add_argument("--maxsubsteps", type=int, default=1000)
    ap.add_argument("--census-tol", type=float, default=1.0e9)
    ap.add_argument("--census-ds", type=float, default=DS_BASE)
    ap.add_argument("--census-steps", type=int, default=3)
    ap.add_argument("--probe-s", type=float, default=0.0085)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    build_hash = ops.ladrunoBuild()
    print(f"[engine] {ops.__file__}", flush=True)
    print(f"[engine] ladrunoBuild() = {build_hash}", flush=True)
    if EXPECTED_BUILD.lower() != "any":
        assert build_hash.startswith(EXPECTED_BUILD), (
            f"ladrunoBuild() = {build_hash}, expected a build starting "
            f"{EXPECTED_BUILD}; set LADRUNO_F10_EXPECT_BUILD to re-attribute.")

    desc, kw = LEGS[args.leg]
    if args.cmd == "leg":
        run_leg(args, args.leg, desc, kw)
    elif args.cmd == "probe":
        run_probe(args, args.leg, desc, kw)
    else:
        run_census(args, args.leg, desc, kw)


if __name__ == "__main__":
    main()
