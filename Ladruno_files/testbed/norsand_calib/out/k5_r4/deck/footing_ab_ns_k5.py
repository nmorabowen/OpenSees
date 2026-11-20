"""WP-144 P3 NORSAND VARIANT (esmeralda ~/wp144/deck, 2026-10-02).  COPY of ref_campaign/footing_ab.py (md5 98399e0a,
~/ladruno_wp152/deck_campaign_toyoura, 2026-09-30); the originals under ~/ladruno_wp152 are untouched.  What this file adds
(the SANISAND machinery is otherwise unchanged; ref_ladder/ keeps the sibling deck with the seating/RO options):
  * --mat norsand --ns-params <json>: nDMaterial LadrunoNorSand, the material's flags taken from the parameter file
    (energy-agnostic: the energy's own flags are a pass-through list).  ONE MATERIAL PER ELEMENT ROW, because the
    material carries its own initial state: sigma0 = gamma' z_row [-K0, -1, -K0, 0, 0, 0] (11 22 33 12 23 13), v0 = 1 + e,
    pi_i0 from the file's rule.  NorSand has no nu* device and no stages: the K0 state is PRESCRIBED, gravity is one
    full step on an (already) equilibrated field, and updateMaterialStage / the stage-flip poissonRatio are not called.
  * --system Pardiso|UmfPack (Pardiso falls back to UmfPack with a log line if the engine has no Pardiso).
  * deck-side census for NorSand from the element response stressesPlaneStrain: min p', max eta, number of Gauss
    points below the reported p'-floor; psi / state at checkpoints; steps.csv gains three trailing columns
    (n_p_lt_F, min_psi, n_latched); summary.json gains "norsand" and "plateau" (plateau.py, the TIMs WP1 rule).
  * --mesh-perturb is refused for norsand (the per-row initial stress is exact only on a rectangular grid).
  ROUND 3b (2026-10-03, sheet 144a at 428aa1489; see README section 14):
  * pi_i0 by the UNIFIED rule (S.53) at the K0 state: the yield surface through (p_init, max(eta_init, c2 M)); the closed
    form here, or the kernel's own -pi0_auto (pi0 mode 'auto'; also used for rows below the p' floor); the committed
    pi_i is cross-checked at the start.  The 2026-10-02 'on_surface' bisection stays as a mode.
  * K0 = the campaign's nu*/(1-nu*) = 0.4554 (ref_campaign/footing_ab.py:111-112) in the campaign geometry, 0.5 in the
    Kimura mode, unless --k0 or the file says otherwise.
  * --energy BA06|HAR and --pmin <kPa> pass through to the material (-energy, -pmin); the file's "elastic_map" maps the
    DM04 targets onto the chosen energy at the deck's e_init (HAR: -k -g -n; BA06: -p0 -kappa_hat -mu0), p_a = 101 kPa.
  * --kimura V85.6|V75.1|V64.6 (with --sand toyoura): D_r, gamma and B_prototype from kimura1985_tests_meta.csv.
  * --init-only: build the materials and elements, read and check the initial field, stop (no analysis).

WP-138 -- strip-footing A/B: SANISAND ModifiedEuler (IntScheme 1) vs SAS-ME
(IntScheme 129, WP-129), on the fork's OWN copy of the TIMs 2D footing deck.

This file is the DECK and the in-process RUNNER. It is meant to be executed in
a subprocess by `launch.py` (wall-clock timeout + plain log), never imported by
a long-lived interpreter: the engine binary is chosen by FOOTING_BIN.

THE DECK (spec: Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md
section 1, on origin/wp/127-sanisand-replay-counters; every gap filled here is
listed in 138_footing_sas_me_ab.md section 2)

  * plane strain, B = 1.5 m strip, domain 15B x 12B = 22.5 m x 18 m, FULL width
    (no symmetry: the act saw asymmetric bands);
  * LadrunoQuad -formulation bbar, 90 x 27 = 2 430 elements, 9 720 Gauss points,
    4 860 free DOF. The mesh is RECONSTRUCTED from the act's own element tags and
    centroids in ring_points_b8.csv: 64 columns of B/8 over |x| <= 4B (= 6 m),
    13 geometrically graded columns on each side; 12 rows of B/8 in the top
    1.5B, 15 graded rows below; lower block numbered first (tags 1..1350,
    column-major), top block after it (1351..2430). Element 1950 is then the act's
    1950 (x = 0.8438, y = -0.0937), and so on for every tag in the CSV.
  * buoyant self-weight gamma' = 9.81 kN/m3 (body force, ramped by selfWeight);
    7.65 kPa surcharge on the top surface OUTSIDE the footprint; 18.4 kN/m
    footing dead load on the footing's reference node;
  * K0 = nu/(1-nu) = 0.4553 (Jaky, phi = 33 deg) reached through nu* = 0.312885
    in the elastic stage (held for the whole run: the material takes nu as a
    constant);
  * rigid ROUGH footing: LadrunoKinematicCoupling from a 3-DOF reference node at
    (0, 0) to the 9 footprint nodes (-dof 1 2, default penalty 1e12); the
    reference node's u_x and rotation are fixed (a guided, centrally pushed
    footing), its u_y carries the dead load and then the displacement push;
  * base fixed, sides on rollers (u_x = 0);
  * system Pardiso, NormUnbalance at 1e-5 x the total applied vertical load,
    ladder Newton -> NewtonLineSearch -> KrylovNewton (tol x 10);
  * LadrunoSANISAND, the campaign set, IntScheme 1 (or 129), TanType 0,
    JacoType 1, TolF = TolR = 1e-7, -flipAlphaIn init, -Pmin 0.0101,
    -maxSubsteps 2000, -Presidual 0, -honorTolR 0;
  * staging: stage 0 (elastic) self-weight, then stage 0 surface loads, then
    updateMaterialStage 1, then the push (sp on the reference node under
    LoadControl, pseudo-time = footing displacement).

UNITS kPa, m, kN. OpenSees sign (tension positive) for stress and strain in
every file unless a column says otherwise; alpha, alpha_in, z are the
material's own (compression-positive-stress) ratios, exactly as `state` returns
them and as ladrunoSANISANDReplay expects them.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import math
import os
import sys
import time

# ---------------------------------------------------------------- engine ----
_BIN = os.environ.get("FOOTING_BIN")
if not _BIN or not os.path.isdir(_BIN):
    raise SystemExit("FOOTING_BIN must name the dist/bin directory to test")
# ESMERALDA COPY (2026-09-28): Linux-portable engine load. On Windows the
# original -S runner needs the DLL dir + the 3.12 site-packages; on Linux the
# module is <FOOTING_BIN>/opensees.so and numpy comes from the venv's site.
if sys.platform == "win32":
    os.add_dll_directory(_BIN)
    sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
sys.path.insert(0, _BIN)

import numpy as np            # noqa: E402
import opensees as ops        # noqa: E402

if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
        os.path.normcase(os.path.abspath(_BIN)):
    raise SystemExit(f"loaded {ops.__file__}, wanted a pyd in {_BIN}")

# --------------------------------------------------------------- geometry ---
B_FOOT = 1.5
H0 = B_FOOT / 8.0                 # 0.1875 m
XLIM = 7.5 * B_FOOT               # 11.25 m
YBOT = -12.0 * B_FOOT             # -18 m
X_FINE = 4.0 * B_FOOT             # 6 m: B/8 columns over |x| <= 4B
Y_FINE = 1.5 * B_FOOT             # 2.25 m: B/8 rows in the top 1.5B
N_XG = 13                         # graded columns per side
N_YG = 15                         # graded rows below the fine band
N_YF = 12
MESH = "b8"
MESH_PERTURB = ""               # WP-150 R2: "" | jitter:A[:SEED] | shear:DEG (mesh_perturb.py)


def set_mesh(name):
    """ESMERALDA COPY: --mesh b8 (default, the act's 90 x 27 mesh) or b16 (every
    cell halved: B/16 fine band, 2 x the graded cell counts, which gives the
    graded ratio ~ sqrt of b8's, i.e. each b8 graded cell ~ two b16 cells).
    Same domain, same fine-band extents, same loads/BCs. b16 = 180 x 54 =
    9 720 elements, 38 880 Gauss points, 17 footprint nodes; element tags keep
    the same numbering RULE (lower block first, column-major), so they do NOT
    match the act's CSV tags (that mapping only exists for b8)."""
    global H0, N_XG, N_YG, N_YF, MESH
    r = {"b4": 0.5, "b8": 1, "b16": 2}[name]      # WP-150 R2: b4 added
    H0 = B_FOOT / (8.0 * r)
    # b4 (WP-150 R2): B/4 fine band, graded counts ~ half of b8's (7 x-graded, 8 y-graded,
    # 6 fine rows): 46 x 14 = 644 elements, 5 footprint nodes.
    N_XG, N_YG, N_YF = {"b4": (7, 8, 6), "b8": (13, 15, 12), "b16": (26, 30, 24)}[name]
    MESH = name
    # --fine-depth (orchestrator 19099d, 2026-09-30): fine rows follow Y_FINE; the default 1.5B gives
    # 12 / 24 / 6 rows for b8 / b16 / b4, exactly the previous hard-coded counts.
    N_YF = int(round(Y_FINE / H0))

GAMMA = 9.81                      # kN/m3 buoyant
Q_SURCH = 7.65                    # kPa outside the footprint
W_FOOT = 18.4                     # kN/m footing dead load
NU_STAR = 0.312885                # K0 = 0.4553 (Jaky, 33 deg)
K0 = NU_STAR / (1.0 - NU_STAR)

# --------------------------------------------------------------- material ---
# campaign set (intake section 1); nu is the nu* device
SAN = [264.32, NU_STAR, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0,
       0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
G0, E_INIT, MC, CC, LAMC, E0, XI, PATM, MM = (SAN[0], SAN[2], SAN[3], SAN[4],
                                              SAN[5], SAN[6], SAN[7], SAN[8],
                                              SAN[9])
NB, ND = SAN[12], SAN[14]
MAT = 1

# ------------------------------------------------------------- controller ---
DS_BASE, DS_MIN, DS_MAX = 2.0e-5, 2.0e-7, 1.0e-3   # m (F10's transcription)
GROW_AFTER, GROW = 6, 2.0
N_GRAV = 10
PUSH_TOL_REL = 1.0e-5
LADDER = (("Newton", 1.0, 25), ("NewtonLineSearch", 1.0, 40),
          ("KrylovNewton", 10.0, 60))

T0 = time.time()
LOG = None


def log(msg):
    line = time.strftime("[%Y-%m-%d %H:%M:%S]") + f" [{time.time()-T0:9.1f}s] {msg}"
    print(line, flush=True)


# ------------------------------------------------------------------ mesh ----
def _graded(lo, hi, n, h0, fine_at_lo):
    """n cells over [lo, hi], first cell h0 at the fine end, geometric growth."""
    a, b = 1.0, 4.0
    for _ in range(200):
        r = 0.5 * (a + b)
        tot = h0 * (r ** n - 1) / (r - 1)
        a, b = (r, b) if tot < hi - lo else (a, r)
    r = 0.5 * (a + b)
    e = np.concatenate([[0.0], np.cumsum(h0 * r ** np.arange(n))])
    e *= (hi - lo) / e[-1]
    return (lo + e) if fine_at_lo else (hi - e[::-1])


def mesh_coords():
    xf = np.linspace(-X_FINE, X_FINE, int(round(2 * X_FINE / H0)) + 1)
    x = np.unique(np.round(np.concatenate([
        _graded(-XLIM, -X_FINE, N_XG, H0, False), xf,
        _graded(X_FINE, XLIM, N_XG, H0, True)]), 12))
    y = np.unique(np.round(np.concatenate([
        _graded(YBOT, -Y_FINE, N_YG, H0, False),
        np.linspace(-Y_FINE, 0.0, N_YF + 1)]), 12))
    # WP-150 R2: the general count (= 90 rr + 1 and 27 rr + 1 for b8/b16, as before)
    assert len(x) == 2 * N_XG + int(round(2 * X_FINE / H0)) + 1 and len(y) == N_YG + N_YF + 1, (len(x), len(y))
    return x, y


def element_table(x, y):
    """(tag, i, j) with the act's numbering: lower block first, then top."""
    ncol = len(x) - 1
    out = []
    for i in range(ncol):
        for j in range(N_YG):
            out.append((i * N_YG + j + 1, i, j))
    base = ncol * N_YG
    for i in range(ncol):
        for j in range(N_YF):
            out.append((base + i * N_YF + j + 1, i, N_YG + j))
    out.sort()
    return out


GP_XI = ((-1, -1), (1, -1), (1, 1), (-1, 1))


def build(mat, scheme, extra, maxsub, dp_phi, dp_g, deterministic, tolr=1.0e-7,
          tantype=0, presidual=0.0, einit=None, over=None, system="Pardiso", init_only=False):
    x, y = mesh_coords()
    nx, ny = len(x), len(y)
    # WP-150 R2: interior fine-zone nodes may be moved (mesh-orientation variant);
    # surface row, fine-zone outline and boundaries are untouched.
    from mesh_perturb import node_coords, gauss_points, min_jacobian_ratio
    XN, YN = node_coords(x, y, X_FINE, Y_FINE, H0, MESH_PERTURB)
    if MESH_PERTURB:
        log(f"mesh perturbation {MESH_PERTURB!r}: worst min/max detJ over elements "
            f"{min_jacobian_ratio(XN, YN):.3f}")

    # 2026-11-13 (O15 element control): ELEM selects the element technology (dp material only for non-default):
    #   qbbar (default, LadrunoQuad bbar) | qstd | qssp | qeas | cst (LadrunoCST, 2 per cell)
    #   lst (LadrunoLST std, 2 per cell) | bez | bezb (BezierTri6, std / -bbar, 2 per cell)
    ELEM = os.environ.get("ELEM", "lst" if os.environ.get("ELEM_LST") == "1" else "qbbar")
    if ELEM not in ("qbbar", "qstd", "qssp", "qeas", "cst", "lst", "bez", "bezb"):
        raise SystemExit(f"unknown ELEM {ELEM}")
    global NGP
    NGP = {"qssp": 1, "cst": 1, "lst": 3, "bez": 3, "bezb": 3}.get(ELEM, 4)   # Gauss points per element
    LST = ELEM in ("lst", "bez", "bezb")          # 6-node triangles on the refined node grid
    TRI = ELEM in ("cst", "lst", "bez", "bezb")   # triangles: no quad loop
    if ELEM != "qbbar":
        if MESH_PERTURB or mat not in ("dp", "norsand"):
            raise SystemExit("ELEM: dp or norsand material and no --mesh-perturb only")
        log(f"ELEM: {ELEM}")
    if LST:
        xc_, yc_ = np.asarray(x, float), np.asarray(y, float)
        x = np.empty(2 * len(xc_) - 1); x[0::2] = xc_; x[1::2] = 0.5 * (xc_[:-1] + xc_[1:])
        y = np.empty(2 * len(yc_) - 1); y[0::2] = yc_; y[1::2] = 0.5 * (yc_[:-1] + yc_[1:])
        nx, ny = len(x), len(y)
        XN, YN = np.meshgrid(x, y, indexing="ij")
        log(f"ELEM_LST: node grid {nx} x {ny} (corner grid {len(xc_)} x {len(yc_)}), 2 LadrunoLST per cell")

    def nt(i, j):
        return i * ny + j + 1

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i in range(nx):
        for j in range(ny):
            ops.node(nt(i, j), float(XN[i, j]), float(YN[i, j]))
    REF = 900001
    ops.node(REF, 0.0, 0.0, "-ndf", 3)

    if mat == "sanisand":
        san = list(SAN)
        if einit is not None:
            san[2] = float(einit)
        # ablation overrides: SAN index -> value (h0 10, nb 12, A0 13, nd 14, zmax 15)
        for i, v in (over or {}).items():
            san[i] = float(v)
        args = ["LadrunoSANISAND", MAT, *san, int(scheme), int(tantype), 1, 1.0e-7, float(tolr),
                "-flipAlphaIn", "init", "-Pmin", 0.0101,
                "-maxSubsteps", int(maxsub), "-Presidual", float(presidual),
                "-honorTolR", 0, *extra]
        ops.nDMaterial(*args)
        matdesc = " ".join(str(a) for a in args)
    elif mat == "norsand":
        matdesc = "<per-row LadrunoNorSand, see below>"      # created once the element rows are known
    elif mat == "dp":
        # UW DruckerPrager, psi = 0 (rhoBar = 0), plane-strain matched cone at
        # psi = 0: sqrt(J2) = sin(phi) p  ->  rho = sqrt(2) sin(phi) / 3.
        phi = math.radians(dp_phi)
        rho = math.sqrt(2.0) * math.sin(phi) / 3.0
        rho_bar = 0.0
        if os.environ.get("DP_ASSOC") == "1":
            # 2026-11-13 (O15 control): ASSOCIATED, psi = phi; plane-strain matched cone for associated flow:
            # sqrt(J2) = 3 tan(phi) / sqrt(9 + 12 tan^2 phi) p  ->  rho = rho_bar = sqrt(2) x that / 3
            t = math.tan(phi)
            rho = rho_bar = math.sqrt(2.0) * (3.0 * t / math.sqrt(9.0 + 12.0 * t * t)) / 3.0
        kk = 2.0 * dp_g * (1.0 + NU_STAR) / (3.0 * (1.0 - 2.0 * NU_STAR))
        args = ["DruckerPrager", MAT, kk, dp_g, 0.0, rho, rho_bar,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ops.nDMaterial(*args)
        matdesc = " ".join(str(a) for a in args)
    else:
        raise SystemExit(f"unknown mat {mat}")

    info = {}
    elems = element_table(x, y) if not TRI else []
    TRI_PENDING = []
    if ELEM == "cst":
        tg = 0
        for ic in range(nx - 1):
            for jc in range(ny - 1):
                N = lambda di, dj: nt(ic + di, jc + dj)
                if (ic + jc) % 2 == 0:
                    tris = [((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))]
                else:
                    tris = [((0, 0), (1, 0), (0, 1)), ((1, 0), (1, 1), (0, 1))]
                for tri in tris:
                    tg += 1
                    TRI_PENDING.append((tg, ("LadrunoCST", [N(*v) for v in tri], None)))
                    cx = [x[ic + v[0]] for v in tri]; cy = [y[jc + v[1]] for v in tri]
                    xc, yc = float(np.mean(cx)), float(np.mean(cy))
                    area = 0.5 * abs((cx[1] - cx[0]) * (cy[2] - cy[0]) - (cx[2] - cx[0]) * (cy[1] - cy[0]))
                    info[tg] = dict(xc=xc, yc=yc, gps=[(xc, yc)] * 4, area=area)
                    elems.append((tg, ic, jc))
    if LST:
        tg = 0
        for ic in range((nx - 1) // 2):
            for jc in range((ny - 1) // 2):
                i0, j0 = 2 * ic, 2 * jc
                N = lambda di, dj: nt(i0 + di, j0 + dj)
                # union-jack: alternate the diagonal so the mesh has no global bias
                if (ic + jc) % 2 == 0:   # diagonal (0,0)-(2,2)
                    tris = [((0, 0), (2, 0), (2, 2), (1, 0), (2, 1), (1, 1)),
                            ((0, 0), (2, 2), (0, 2), (1, 1), (1, 2), (0, 1))]
                else:                    # diagonal (2,0)-(0,2)
                    tris = [((0, 0), (2, 0), (0, 2), (1, 0), (1, 1), (0, 1)),
                            ((2, 0), (2, 2), (0, 2), (2, 1), (1, 2), (1, 1))]
                for tri in tris:
                    tg += 1
                    TRI_PENDING.append((tg, ("LadrunoLST" if ELEM == "lst" else "BezierTri6", [N(*v) for v in tri], None)))
                    cx = [x[i0 + v[0]] for v in tri[:3]]; cy = [y[j0 + v[1]] for v in tri[:3]]
                    xc, yc = float(np.mean(cx)), float(np.mean(cy))
                    area = 0.5 * abs((cx[1] - cx[0]) * (cy[2] - cy[0]) - (cx[2] - cx[0]) * (cy[1] - cy[0]))
                    info[tg] = dict(xc=xc, yc=yc, gps=[(xc, yc)] * 4, area=area)
                    elems.append((tg, ic, jc))
    ns_tag, ns_rows = {}, None
    if mat == "norsand":
        if MESH_PERTURB and not MESH_PERTURB.startswith("shear"):
            raise SystemExit("norsand deck: only --mesh-perturb shear:DEG is supported (rows keep their depth)")
        ns_tag, ns_rows, matdesc = ns_make_materials(elems, (y[0::2] if LST else y), einit)   # rows = coarse cells
    for tg, (etype, nodes, _) in TRI_PENDING:
        m = ns_tag.get(tg, MAT)
        if etype == "BezierTri6":   # corners then mid12, mid23, mid31 (control points at edge midpoints)
            ops.element("BezierTri6", tg, *nodes, 1.0, "PlaneStrain", m,
                        *(["-bbar"] if ELEM == "bezb" else []), "-bodyForce", 0.0, -GAMMA)
        else:
            ops.element(etype, tg, *nodes, m, "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
    for tag, i, j in (elems if not TRI else []):
        ops.element("LadrunoQuad", tag, nt(i, j), nt(i + 1, j), nt(i + 1, j + 1),
                    nt(i, j + 1), ns_tag.get(tag, MAT), "-formulation", ELEM[1:],
                    "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
        if not MESH_PERTURB:
            xc = 0.25 * (x[i] + x[i + 1]) * 2
            yc = 0.25 * (y[j] + y[j + 1]) * 2
            gps = []
            for xi, et in GP_XI:
                g = 1.0 / math.sqrt(3.0)
                gps.append((0.5 * (x[i] + x[i + 1]) + 0.5 * xi * g * (x[i + 1] - x[i]),
                            0.5 * (y[j] + y[j + 1]) + 0.5 * et * g * (y[j + 1] - y[j])))
            area = (x[i + 1] - x[i]) * (y[j + 1] - y[j])
        else:   # WP-150 R2: the true (distorted) geometry
            cx = [XN[i, j], XN[i + 1, j], XN[i + 1, j + 1], XN[i, j + 1]]
            cy = [YN[i, j], YN[i + 1, j], YN[i + 1, j + 1], YN[i, j + 1]]
            gps, area = gauss_points(cx, cy)
            xc = float(np.mean([g[0] for g in gps]))
            yc = float(np.mean([g[1] for g in gps]))
        info[tag] = dict(xc=float(xc), yc=float(yc), gps=gps, area=float(area))

    for i in range(nx):
        for j in range(ny):
            if j == 0:
                ops.fix(nt(i, j), 1, 1)
            elif i == 0 or i == nx - 1:
                ops.fix(nt(i, j), 1, 0)
    ops.fix(REF, 1, 0, 1)

    half = 0.5 * B_FOOT
    foot = [nt(i, ny - 1) for i in range(nx) if abs(x[i]) <= half + 1e-9]
    assert len(foot) == (2 if LST else 1) * int(round(B_FOOT / H0)) + 1, len(foot)      # 5 / 9 / 17 (b4 / b8 / b16)
    if os.environ.get("FOOT_SMOOTH") == "1":     # 2026-11-13 (O15 control): SMOOTH footing, vertical DOF only coupled
        ops.element("LadrunoKinematicCoupling", 900001, REF, len(foot), *foot, "-dof", 2)
        log("footing: SMOOTH (LadrunoKinematicCoupling -dof 2 only; footprint nodes free in x)")
    else:
        ops.element("LadrunoKinematicCoupling", 900001, REF, len(foot), *foot,
                    "-dof", 1, 2)
    if mat == "norsand":
        # the initial state, read BEFORE any analysis (stage 0b and the push change pi_i and the stresses): the committed pi_i
        # against what the deck asked for, and the kernel's initial-floor census
        _d0 = dict(tags=[t for t, _, _ in elems], ns_rows=ns_rows)
        NS["_check"] = ns_check_pi0(_d0, read_field_ns(_d0, full=True))
        NS["_floor_init"] = ns_floor_census(_d0)
        log(f"p' floor at the initial state: p_min_eff {NS['pmin_eff']:.4f} kPa, rows below it "
            f"{sum(1 for r in ns_rows if r['floored_init'])} of {len(ns_rows)}; kernel census {NS['_floor_init']}")
    if init_only:       # --init-only: materials and elements exist (initial states committed); no analysis objects, no stages
        return dict(x=x, y=y, nx=nx, ny=ny, tags=[t for t, _, _ in elems], info=info, foot=foot, REF=REF,
                    applied=float("nan"), solver="none (--init-only)", matdesc=matdesc, grav_resultant_err=float("nan"),
                    k0_patch_err=float("nan"), mat=mat, ns_rows=ns_rows)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    if deterministic:
        try:
            ops.system("Pardiso", "-deterministic")
            solver = "Pardiso -deterministic"
        except Exception as exc:   # binary without WP-132
            log(f"NOTE Pardiso -deterministic refused ({exc}); plain Pardiso")
            ops.system("Pardiso")
            solver = "Pardiso"
    else:
        solver = system
        try:
            ops.system(system)
        except Exception as exc:
            log(f"NOTE system {system} refused ({exc}); falling back to UmfPack")
            ops.system("UmfPack")
            solver = "UmfPack (fallback)"
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    if mat != "norsand":
        ops.updateMaterialStage("-material", MAT, "-stage", 0)

    # ---- stage 0a: self weight -------------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tags = [t for t, _, _ in elems]
    ops.eleLoad("-ele", *tags, "-type", "-selfWeight", 0.0, 1.0)
    n_grav = 1 if mat == "norsand" else N_GRAV     # norsand: the field is prescribed in equilibrium, so ONE full step
    ops.integrator("LoadControl", 1.0 / n_grav)
    for k in range(n_grav):
        rc = ops.analyze(1)
        if rc != 0:
            raise SystemExit(f"gravity step {k+1} failed rc={rc}")
    if mat == "norsand":
        umax = max(max(abs(v) for v in ops.nodeDisp(n)[:2]) for n in ops.getNodeTags() if n != REF)
        log(f"norsand gravity: max |u| after the one-step self weight = {umax:.3e} m "
            f"(0 = the prescribed K0 field is in equilibrium with gamma' = {GAMMA})")
    ops.reactions()
    base = [nt(i, 0) for i in range(nx)]
    ry = sum(ops.nodeReaction(n, 2) for n in base)
    want = GAMMA * (2 * XLIM) * (-YBOT)
    res_err = abs(ry / want - 1.0)
    # 1-D geostatic patch at the element average of the 4 Gauss points
    perr = 0.0
    for tag in tags:
        st = ops.eleResponse(tag, "stresses")
        syy = float(np.mean(st[1::3]))
        sxx = float(np.mean(st[0::3]))
        ex = GAMMA * info[tag]["yc"]
        perr = max(perr, abs(syy - ex) / abs(ex), abs(sxx - K0 * ex) / abs(K0 * ex))
    log(f"gravity: base reaction {ry:.4f} kN/m vs gamma*A {want:.4f} "
        f"(rel {res_err:.2e}); 1-D K0 patch max rel err {perr:.2e}")
    ops.loadConst("-time", 0.0)

    # ---- stage 0b: surcharge outside + footing dead load ---------------------
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    trib = np.zeros(nx)
    for i in range(nx - 1):
        if 0.5 * (x[i] + x[i + 1]) > half or 0.5 * (x[i] + x[i + 1]) < -half:
            L = x[i + 1] - x[i]
            trib[i] += 0.5 * L
            trib[i + 1] += 0.5 * L
    surf = 0.0
    for i in range(nx):
        if trib[i] > 0:
            ops.load(nt(i, ny - 1), 0.0, -Q_SURCH * trib[i])
            surf += Q_SURCH * trib[i]
    ops.load(REF, 0.0, -W_FOOT, 0.0)
    ops.integrator("LoadControl", 1.0 / N_GRAV)
    for k in range(N_GRAV):
        _t0 = time.time()
        rc = ops.analyze(1)
        if rc != 0:
            raise SystemExit(f"surface-load step {k+1} failed rc={rc}")
        if mat == "norsand":
            log(f"  stage 0b step {k+1}/{N_GRAV}: Newton iterations {ops.testIter()}, wall {time.time()-_t0:.2f} s")
    ops.reactions()
    ry2 = sum(ops.nodeReaction(n, 2) for n in base)
    applied = want + surf + W_FOOT
    log(f"surface loads: surcharge {surf:.4f} kN/m (+{W_FOOT} footing); base "
        f"reaction {ry2:.4f} vs {applied:.4f} (rel {abs(ry2/applied-1):.2e})")
    ops.loadConst("-time", 0.0)
    deck = dict(x=x, y=y, nx=nx, ny=ny, tags=tags, info=info, foot=foot,
                REF=REF, applied=applied, solver=solver, matdesc=matdesc,
                grav_resultant_err=res_err, k0_patch_err=perr, mat=mat, ns_rows=ns_rows)
    return deck


# ------------------------------------------------------------ field reads ---
STAT_NAMES = ["updates", "meCalls", "substeps", "accepted", "rejectedErr",
              "forcedAtDTmin", "forcedClampMc", "rejectedLowP", "abandonedLowP",
              "capHits", "entryPminClamps", "pnResets", "maxSubstepsOneUpdate",
              "lastSubsteps", "lastForcedAtDTmin", "lastAbandonedLowP",
              "lastCapHit"]


SAS_NAMES = ["updates", "elastic", "substeps", "accepted", "rejectedErr",
             "rejectedLowP", "rejectedNonPosH", "rejectedDrift", "rejectedAlpha",
             "elasticStages", "driftCorrections", "alphaInReseats", "hBrackets",
             "alphaProjected", "intersectFail", "refusals", "refStartF",
             "refStartAlpha", "refStartOther", "refDTmin", "refNonPosH",
             "refLowP", "refDrift", "refAlpha", "refCap",
             "maxSubstepsOneUpdate", "lastSubsteps", "lastRefuseCode",
             "maxAlphaRatio", "lastAlphaRatio", "lastF",
             "entryOverKappa", "rejectedReversal",
             "sas_hFloored", "sas_hSoftCapped", "sas_reseatHeld",
             "sepEntriesTension", "sepEntriesLowP", "sepExits", "sepActive",
             "sepLastCode", "sepMaxP0", "sepHeldHighP", "sepHeldCompressing"]   # WP-152 review appended 4   # WP-152 appended 4   # WP-151 (bd93c558d) appended 3   # beb6d8333 appended the last two
# which census the run reads, and where the columns the step line reports sit.
# ME (substepStats): substeps, capHits, rejectedErr, forcedAtDTmin, rejectedLowP,
# forcedClampMc, lastSubsteps, lastCapHit. SAS (sasStats): substeps, refusals,
# rejectedErr, refDTmin, rejectedLowP, refAlpha, lastSubsteps, lastRefuseCode.
CENSUS = dict(resp="substepStats", names=STAT_NAMES, sub=2, cap=9, rej=4, dtmin=5,
              lowp=7, clamp=6, last_sub=13, last_cap=16)
CENSUS_SAS = dict(resp="sasStats", names=SAS_NAMES, sub=2, cap=15, rej=4, dtmin=19,
                  lowp=5, clamp=23, last_sub=26, last_cap=27)
SAS_REF_CODES = {16: "startOutsideYield", 17: "startAlphaOutsideBounding",
                 18: "startInadmissible", 19: "errorAtDTmin", 20: "loadingNonPosH",
                 21: "tensionAtDTmin(lowP)", 22: "driftFailed",
                 23: "alphaOutsideAtDTmin", 24: "maxSubsteps"}


DEADEND = {}     # SAS-ME: Gauss point index -> step first flagged


def read_field(deck, want_stats=True, want_f=False):
    """Every Gauss point's committed state. Arrays of shape (nGP, ...)."""
    if deck["mat"] == "norsand":
        return read_field_ns(deck, full=want_f)
    tags = deck["tags"]
    n = 4 * len(tags)
    sig = np.zeros((n, 3)); eps = np.zeros((n, 3)); st = np.zeros((n, 26))
    psi = np.full(n, np.nan); fy = np.full(n, np.nan)
    stats = np.zeros((n, len(CENSUS["names"]))) if want_stats else None
    k = 0
    for tag in tags:
        for gp in (1, 2, 3, 4):
            _r = ops.eleResponse(tag, "material", gp, "stress")
            if len(_r) == 0:            # ELEM_LST: 3 Gauss points, the 4th slot stays zero
                k += 1
                continue
            sig[k] = _r
            eps[k] = ops.eleResponse(tag, "material", gp, "strain")
            if deck["mat"] == "sanisand":
                st[k] = ops.eleResponse(tag, "material", gp, "state")
                psi[k] = ops.eleResponse(tag, "material", gp, "psi")[0]
                if want_f:
                    fy[k] = ops.eleResponse(tag, "material", gp, "yieldDistance")[0]
                if want_stats:
                    v = ops.eleResponse(tag, "material", gp, CENSUS["resp"])
                    if v:
                        stats[k, :len(v)] = v
            k += 1
    return dict(sig=sig, eps=eps, st=st, psi=psi, f=fy, stats=stats)


def derived(F):
    """p, q, eta, full 6-stress, rho_alpha, e from a SANISAND field.

    sigma_zz is not exposed by the plane-strain wrapper, so it is RECOVERED
    from the committed psi and e (p_r = 0): e_c = e - psi = e0 - lc (p/Pa)^xi.
    `yieldDistance` checks it at every checkpoint (f recomputed vs f read)."""
    sig, st, psi = F["sig"], F["st"], F["psi"]
    e = st[:, 24]
    arg = (E0 - (e - psi)) / LAMC
    p = PATM * np.power(np.clip(arg, 0.0, None), 1.0 / XI)
    sc = np.zeros((len(p), 6))              # compression positive
    sc[:, 0] = -sig[:, 0]; sc[:, 1] = -sig[:, 1]; sc[:, 3] = -sig[:, 2]
    sc[:, 2] = 3.0 * p - sc[:, 0] - sc[:, 1]
    s = sc.copy(); s[:, :3] -= p[:, None]
    nrm = lambda v: np.sqrt(v[:, 0]**2 + v[:, 1]**2 + v[:, 2]**2
                            + 2 * (v[:, 3]**2 + v[:, 4]**2 + v[:, 5]**2))
    qv = math.sqrt(1.5) * nrm(s)
    eta = np.where(p > 1e-12, qv / np.maximum(p, 1e-300), np.nan)
    al = st[:, 6:12]
    an = nrm(al)
    # alpha's own Lode angle: cos3th = sqrt(6) tr(n^3), n = alpha/|alpha|
    cos3 = np.zeros(len(p))
    for k in range(len(p)):
        if an[k] > 0:
            a = al[k] / an[k]
            Mx = np.array([[a[0], a[3], a[5]], [a[3], a[1], a[4]], [a[5], a[4], a[2]]])
            cos3[k] = math.sqrt(6.0) * np.trace(Mx @ Mx @ Mx)
    cos3 = np.clip(cos3, -1, 1)
    g = 2 * CC / ((1 + CC) - (1 - CC) * cos3)
    ab = math.sqrt(2.0 / 3.0) * (g * MC * np.exp(-NB * psi) - MM)
    rho = an / ab
    # f recomputed from the recovered stress, to check against yieldDistance
    sd = s - p[:, None] * al
    frec = nrm(sd) - math.sqrt(2.0 / 3.0) * MM * p
    mbc = MC * np.exp(-NB * psi)
    return dict(p=p, q=qv, eta=eta, s6=-sc, rho=rho, e=e, frec=frec,
                eta_mb=eta / mbc)


def gp_xy(deck):
    out = []
    for tag in deck["tags"]:
        for gx, gy in deck["info"][tag]["gps"]:
            out.append((tag, gx, gy))
    return out


# ------------------------------------------------------------ replay CSVs ---
# SIGN CONVENTION of every replay CSV: COMPRESSION POSITIVE for sigma, dStrain
# and sigma_next (the model's internal convention, and what the TIMs ring CSVs
# actually carry -- WP-127 finding A), so a row feeds
# `ladrunoSANISANDReplay -convention compressionPositive` and WP-134's
# sanisand_reference unchanged. Shear strain is ENGINEERING (gamma_xy) in slot
# 3. alpha, alpha_in, z are the raw internal ratios. `dt_next` is the pseudo-
# time increment of the next step (= ops_Dt, negative: the push runs time
# down); `prevIncrNorm` is GetNorm_Cov of the increment that produced the
# committed state (the -reversalRel reference the material had).
REPLAY_HEAD = (["element", "gp", "x_m", "y_m", "p_kPa", "eta",
                "eta_over_Mb_compression", "e", "psi"]
               + [f"sigma_{i}" for i in range(6)] + [f"alpha_{i}" for i in range(6)]
               + [f"alpha_in_{i}" for i in range(6)] + [f"z_{i}" for i in range(6)]
               + [f"dStrain_{i}" for i in range(6)]
               + ["step", "s_over_B", "gp_x_m", "gp_y_m", "rho_alpha", "f_read",
                  "f_recomputed", "substeps_next", "capHit_next", "dt_next",
                  "prevIncrNorm"]
               + [f"sigma_next_{i}" for i in range(6)]
               + [f"alpha_next_{i}" for i in range(6)]
               + ["e_next", "select"])


def select_worst(D, F, stats_next=None, n_rho=20, n_p=15, n_sub=15):
    sel = {}
    rho = np.nan_to_num(D["rho"], nan=-1)
    for k in np.argsort(-rho)[:n_rho]:
        sel.setdefault(int(k), []).append("rho")
    for k in np.argsort(D["p"])[:n_p]:
        sel.setdefault(int(k), []).append("lowp")
    if stats_next is not None:
        for k in np.argsort(-stats_next)[:n_sub]:
            if stats_next[k] > 0:
                sel.setdefault(int(k), []).append("sub")
    return sel


def write_replay(path, deck, F, D, Fn, Dn, eps_prev, sel, step, sB, sub_next,
                 cap_next, dt_next):
    """Committed state (F, D) at `step` + the increment to (Fn, Dn).

    Fn may be a committed next step or a committed probe iterate; Dn may be
    None (next-state columns then NaN)."""
    xy = gp_xy(deck)
    nan6 = [float("nan")] * 6
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(REPLAY_HEAD)
        for k in sorted(sel, key=lambda k: -np.nan_to_num(D["rho"][k])):
            tag, gx, gy = xy[k]
            inf = deck["info"][tag]
            de = -(Fn["eps"][k] - F["eps"][k])            # compression positive
            de6 = [de[0], de[1], 0.0, de[2], 0.0, 0.0]
            if eps_prev is not None:
                dp = F["eps"][k] - eps_prev[k]
                pin = math.sqrt(dp[0]**2 + dp[1]**2 + 0.5 * dp[2]**2)
            else:
                pin = 0.0
            st = F["st"][k]
            sn = (-Dn["s6"][k]).tolist() if Dn is not None else nan6
            an = Fn["st"][k][6:12].tolist() if Dn is not None else nan6
            en = float(Fn["st"][k][24]) if Dn is not None else float("nan")
            w.writerow([tag, k % 4 + 1, f"{inf['xc']:.4f}", f"{inf['yc']:.4f}",
                        repr(float(D["p"][k])), repr(float(D["eta"][k])),
                        repr(float(D["eta_mb"][k])), repr(float(D["e"][k])),
                        repr(float(F["psi"][k]))]
                       + [repr(float(v)) for v in -D["s6"][k]]
                       + [repr(float(v)) for v in st[6:12]]
                       + [repr(float(v)) for v in st[18:24]]
                       + [repr(float(v)) for v in st[12:18]]
                       + [repr(float(v)) for v in de6]
                       + [step, f"{sB:.8f}", f"{gx:.4f}", f"{gy:.4f}",
                          repr(float(D["rho"][k])), repr(float(F["f"][k])),
                          repr(float(D["frec"][k])),
                          int(sub_next[k]) if sub_next is not None else -1,
                          int(cap_next[k]) if cap_next is not None else -1,
                          repr(float(dt_next)), repr(float(pin))]
                       + [repr(float(v)) for v in sn]
                       + [repr(float(v)) for v in an]
                       + [repr(en), "+".join(sel[k])])


def save_field(path, deck, F, D, step, sB):
    xy = np.array([(t, gx, gy) for t, gx, gy in gp_xy(deck)])
    nd_tag = np.array(ops.getNodeTags(), dtype=int)
    nd_xy = np.array([ops.nodeCoord(int(n))[:2] for n in nd_tag])
    nd_u = np.array([ops.nodeDisp(int(n))[:2] for n in nd_tag])
    np.savez_compressed(path, step=step, s_over_B=sB, tag=xy[:, 0], gx=xy[:, 1],
                        gy=xy[:, 2], nd_tag=nd_tag, nd_xy=nd_xy, nd_u=nd_u, sig=F["sig"], eps=F["eps"], st=F["st"],
                        psi=F["psi"], f=F["f"],
                        stats=F["stats"] if F["stats"] is not None else 0,
                        **({k: v for k, v in D.items()} if D else {}))


def ring_summary(deck, D, F, top=10):
    xy = gp_xy(deck)
    lines = []
    for k in np.argsort(-np.nan_to_num(D["rho"], nan=-1))[:top]:
        t, gx, gy = xy[k]
        lines.append(f"    rho_a {D['rho'][k]:7.4f} eta {D['eta'][k]:7.3f} p' "
                     f"{D['p'][k]:8.3f} ele {t} gp {k%4+1} x {gx:+.3f} y {gy:+.3f}")
    lines.append("    lowest p':")
    for k in np.argsort(D["p"])[:5]:
        t, gx, gy = xy[k]
        lines.append(f"    p' {D['p'][k]:8.4f} eta {D['eta'][k]:7.3f} rho_a "
                     f"{D['rho'][k]:7.4f} ele {t} gp {k%4+1} x {gx:+.3f} y {gy:+.3f}")
    return "\n".join(lines)


# ------------------------------------------------------------ toyoura -------
# GATE 1 (deck_toyoura COPY, 2026-09-28): `--sand toyoura` re-targets this deck
# at the Kimura, Kusakabe & Saitoh (1985, Geotechnique 35(1):33-45) centrifuge
# tests of ROUGH strip footings on DENSE dry Toyoura sand, modelled at
# PROTOTYPE scale at 1 g (B_proto = B_model * N; stresses identical).
#   * B = 1.2 m (B_model 30 mm at 40 g: Figs 3, 5, 6);
#   * domain 16.67B x 10B = the Mark I specimen container (Table 1, p.34:
#     length 0.5 m, depth 0.3 m) at B_model = 30 mm; the sand-layer depth is
#     NOT stated in the paper, so the full container depth is assumed;
#   * NO surcharge, NO footing dead load (neither is reported);
#   * dry unit weight 15.9 kN/m3 (p.36, the vibro-compacted series of Fig. 3);
#   * e_init from D_r (default 0.88: Kimura's dense tests at D_r 86-90 %,
#     Figs 9/17 and the slope series; Fig. 4's dense band is 89-97 %) with
#     e_max 0.977 / e_min 0.597 (DM04 p.631's own D_r labels, "18.5 to 63.7 %
#     (e 0.907 to 0.735)", = Verdugo & Ishihara); Kimura et al. (1985) do NOT
#     report their batch's e_max/e_min, e or dry density (only gamma 15.9);
#   * K0 = 0.5 (Loukidis & Salgado 2011's main series; Jaky with the DM04
#     critical-state angle, 1 - sin 31.2 deg = 0.48, agrees) via nu* = K0/(1+K0)
#     in the elastic stage, then nu -> DM04's 0.05 at the stage flip
#     (setParameter poissonRatio), so the push runs the DM04 set VERBATIM;
#   * DM04 Table 1 Toyoura constants, P_atm = 100 kPa (not tabulated in DM04;
#     Gate 0 / WP-150 memo section 12 choice);
#   * step controller constants (m) scaled by B/1.5 so the ladder is the same
#     in s/B as the campaign deck's.
SAND = "campaign"
NU_PUSH = None
TYR_EMAX, TYR_EMIN = 0.977, 0.597
TYR = dict(G0=125.0, nu=0.05, Mc=1.25, c=0.712, lamc=0.019, e0=0.934, xi=0.7,
           patm=100.0, m=0.01, h0=7.05, ch=0.968, nb=1.1, A0=0.704, nd=3.5,
           zmax=4.0, cz=600.0)
TYR_INFO = {}


def sand_parser():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--sand", choices=("campaign", "toyoura", "toyoura_campaign"), default="campaign",
                   help="GATE 1: toyoura = Kimura et al. (1985) prototype deck, DM04 Toyoura set; "
                        "toyoura_campaign = the campaign's geometry/loads/staging with the DM04 Toyoura set")
    p.add_argument("--fine-depth", type=float, default=None,
                   help="toyoura modes: depth of the B/8 fine band in B, default 1.5 (orchestrator 19099d)")
    p.add_argument("--bproto", type=float, default=None, help="toyoura: prototype B (m), default 1.2")
    p.add_argument("--dr", type=float, default=None, help="toyoura: relative density (0-1), default 0.88")
    p.add_argument("--k0", type=float, default=None, help="toyoura: K0 of the elastic stage, default 0.5")
    p.add_argument("--gamma", type=float, default=None, help="toyoura: dry unit weight kN/m3, default 15.9")
    p.add_argument("--dom-w", type=float, default=None, help="toyoura: domain width in B, default 16.667")
    p.add_argument("--dom-d", type=float, default=None, help="toyoura: domain depth in B, default 10")
    p.add_argument("--nu-push", type=float, default=None,
                   help="toyoura: nu set at the stage flip (default DM04 0.05; < 0 keeps nu*)")
    p.add_argument("--kimura", default=None,
                   help="toyoura (Kimura mode) only: a test_id of kimura1985_tests_meta.csv (V85.6 | V75.1 | V64.6 | H...): "
                        "D_r, gamma (dry unit weight) and B_prototype from the CSV unless --dr/--gamma/--bproto are given")
    p.add_argument("--kimura-meta", default=None,
                   help="path of tests_meta.csv (default: kimura1985_tests_meta.csv next to this file, an unedited copy of "
                        "norsand_calib/data/kimura1985/tests_meta.csv)")
    return p


KIMURA_INFO = {}


def resolve_kimura(pa):
    """--kimura <test_id>: fill pa.dr, pa.gamma, pa.bproto from the data pack's tests_meta.csv (Kimura, Kusakabe & Saitoh
    1985 Fig. 9: D_r from the curve labels; e_assumed = 0.977 - D_r (0.977 - 0.597) and gamma_derived = Gs 9.81/(1 + e),
    Gs = 2.64, are DERIVED there, the paper gives neither: see the CSV header).  An explicit --dr/--gamma/--bproto wins."""
    path = pa.kimura_meta or os.path.join(os.path.dirname(os.path.abspath(__file__)), "kimura1985_tests_meta.csv")
    raw = open(path, "rb").read()
    lines = [ln for ln in raw.decode("utf-8").splitlines() if ln.strip() and not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    row = next((r for r in rows if r["test_id"] == pa.kimura), None)
    if row is None:
        raise SystemExit(f"--kimura {pa.kimura!r}: not in {path} (test ids: {[r['test_id'] for r in rows]})")
    src = {}
    for key, val in (("dr", float(row["D_r_pct"]) / 100.0), ("gamma", float(row["gamma_derived_kN_m3"])),
                     ("bproto", float(row["B_prototype_m"]))):
        if getattr(pa, key) is None:
            setattr(pa, key, val)
            src[key] = f"CSV {val}"
        else:
            src[key] = f"CLI override {getattr(pa, key)} (CSV {val})"
    e_csv = float(row["e_assumed"])
    e_deck = TYR_EMAX - pa.dr * (TYR_EMAX - TYR_EMIN)
    KIMURA_INFO.update(test_id=pa.kimura, csv=os.path.abspath(path), csv_md5=hashlib.md5(raw).hexdigest(), source=src,
                       e_assumed_csv=e_csv, e_deck=e_deck, q_peak_kPa=float(row["q_peak_kPa"]),
                       s_over_B_at_peak=float(row["s_over_B_at_peak"]), case=row["case"], D_r_pct=float(row["D_r_pct"]))


def set_toyoura(pa):
    global SAND, NU_PUSH, B_FOOT, XLIM, YBOT, X_FINE, Y_FINE, GAMMA, Q_SURCH, W_FOOT
    global NU_STAR, K0, SAN, G0, E_INIT, MC, CC, LAMC, E0, XI, PATM, MM, NB, ND
    global DS_BASE, DS_MIN, DS_MAX
    SAND = "toyoura"
    B_FOOT = 1.2 if pa.bproto is None else float(pa.bproto)
    dom_w = 50.0 / 3.0 if pa.dom_w is None else float(pa.dom_w)
    dom_d = 10.0 if pa.dom_d is None else float(pa.dom_d)
    XLIM = 0.5 * dom_w * B_FOOT
    YBOT = -dom_d * B_FOOT
    X_FINE = 4.0 * B_FOOT
    FD = 1.5 if pa.fine_depth is None else float(pa.fine_depth)
    Y_FINE = FD * B_FOOT
    if not (XLIM > X_FINE + 1.0e-9 and -YBOT > Y_FINE + 1.0e-9):
        raise SystemExit(f"--dom-w/--dom-d must exceed the fine band (8B wide, 1.5B deep)")
    GAMMA = 15.9 if pa.gamma is None else float(pa.gamma)
    Q_SURCH = 0.0
    W_FOOT = 0.0
    K0 = 0.5 if pa.k0 is None else float(pa.k0)
    NU_STAR = K0 / (1.0 + K0)
    dr = 0.88 if pa.dr is None else float(pa.dr)
    e_init = TYR_EMAX - dr * (TYR_EMAX - TYR_EMIN)
    NU_PUSH = TYR["nu"] if pa.nu_push is None else (None if pa.nu_push < 0 else float(pa.nu_push))
    t = TYR
    SAN = [t["G0"], NU_STAR, e_init, t["Mc"], t["c"], t["lamc"], t["e0"], t["xi"],
           t["patm"], t["m"], t["h0"], t["ch"], t["nb"], t["A0"], t["nd"], t["zmax"],
           t["cz"], GAMMA / 9.81]
    G0, E_INIT, MC, CC, LAMC, E0, XI, PATM, MM = (SAN[0], SAN[2], SAN[3], SAN[4],
                                                  SAN[5], SAN[6], SAN[7], SAN[8], SAN[9])
    NB, ND = SAN[12], SAN[14]
    sc = B_FOOT / 1.5
    DS_BASE, DS_MIN, DS_MAX = DS_BASE * sc, DS_MIN * sc, DS_MAX * sc
    TYR_INFO.update(B=B_FOOT, dom_w_B=dom_w, dom_d_B=dom_d, XLIM=XLIM, YBOT=YBOT,
                    gamma=GAMMA, K0=K0, nu_star=NU_STAR, nu_push=NU_PUSH, Dr=dr,
                    e_init=e_init, emax=TYR_EMAX, emin=TYR_EMIN, ds_scale=sc,
                    ds=(DS_BASE, DS_MIN, DS_MAX), SAN=list(SAN))
    TYR_INFO.update(fine_depth_B=FD)


# TOYOURA_CAMPAIGN (deck_campaign_toyoura, orchestrator 19099d, 2026-09-30):
# the campaign's geometry, loads and staging UNCHANGED (B 1.5 m, box 15B x 12B,
# B/8 fine block |x| <= 4B over the top 1.5B, 18.4 kN/m footing weight, K0 =
# 0.4553 via nu* = 0.312885 in the elastic stage, DS ladder unscaled), with the
# DM04 Table 1 Toyoura set VERBATIM (P_atm 100) and nu -> DM04's 0.05 at the
# stage flip, exactly as --sand toyoura does. e_init from --dr (e_max 0.977 /
# e_min 0.597). --gamma, when given, REPLACES the body-force unit weight (the
# campaign's buoyant gamma' = 9.81) AND the gamma_d of the embedment surcharge:
# q0 = gamma_d * d_Ef with the campaign's d_Ef = 0.5 m (15.3 x 0.5 = 7.65 kPa);
# without --gamma both stay the campaign's (9.81 body, q0 7.65).
CAMP_D_EF = 0.5
CAMP_GAMMA_D = 15.3


def set_toyoura_campaign(pa):
    global SAND, NU_PUSH, Y_FINE, GAMMA, Q_SURCH, NU_STAR, K0
    global SAN, G0, E_INIT, MC, CC, LAMC, E0, XI, PATM, MM, NB, ND
    SAND = "toyoura_campaign"
    for k in ("bproto", "dom_w", "dom_d"):
        if getattr(pa, k) is not None:
            raise SystemExit("--bproto/--dom-w/--dom-d are fixed by the campaign geometry in toyoura_campaign")
    FD = 1.5 if pa.fine_depth is None else float(pa.fine_depth)
    Y_FINE = FD * B_FOOT
    if not (-YBOT > Y_FINE + 1.0e-9):
        raise SystemExit("--fine-depth must stay above the box base")
    gamma_d = CAMP_GAMMA_D if pa.gamma is None else float(pa.gamma)
    GAMMA = GAMMA if pa.gamma is None else float(pa.gamma)
    Q_SURCH = gamma_d * CAMP_D_EF
    if pa.k0 is not None:
        K0 = float(pa.k0)
        NU_STAR = K0 / (1.0 + K0)
    dr = 0.88 if pa.dr is None else float(pa.dr)
    e_init = TYR_EMAX - dr * (TYR_EMAX - TYR_EMIN)
    NU_PUSH = TYR["nu"] if pa.nu_push is None else (None if pa.nu_push < 0 else float(pa.nu_push))
    t = TYR
    SAN = [t["G0"], NU_STAR, e_init, t["Mc"], t["c"], t["lamc"], t["e0"], t["xi"],
           t["patm"], t["m"], t["h0"], t["ch"], t["nb"], t["A0"], t["nd"], t["zmax"],
           t["cz"], gamma_d / 9.81]
    G0, E_INIT, MC, CC, LAMC, E0, XI, PATM, MM = (SAN[0], SAN[2], SAN[3], SAN[4],
                                                  SAN[5], SAN[6], SAN[7], SAN[8], SAN[9])
    NB, ND = SAN[12], SAN[14]
    TYR_INFO.update(mode="toyoura_campaign", B=B_FOOT, XLIM=XLIM, YBOT=YBOT,
                    dom_w_B=2 * XLIM / B_FOOT, dom_d_B=-YBOT / B_FOOT,
                    X_FINE=X_FINE, Y_FINE=Y_FINE, fine_depth_B=FD,
                    gamma_body=GAMMA, gamma_d_q0=gamma_d, d_Ef=CAMP_D_EF, q0=Q_SURCH,
                    W_FOOT=W_FOOT, K0=K0, nu_star=NU_STAR, nu_push=NU_PUSH, Dr=dr,
                    e_init=e_init, emax=TYR_EMAX, emin=TYR_EMIN,
                    ds=(DS_BASE, DS_MIN, DS_MAX), SAN=list(SAN))



# ------------------------------------------------------------ NORSAND (WP-144 P3) -
# Parameter file (JSON, see ns_params_toyoura.json).  Two formats:
#   * "elastic_map" (round 3b, default): the DM04 elastic targets (G0, nu, n_exp) mapped onto the energy chosen by
#     "energy" / --energy at the deck's own void ratio (the BVP policy 'global' of harness/energy.py: ONE material per
#     row family, e = the deck's e_init, so --dr / --kimura move the stiffness with e).  HAR: k = K(p_a, e)/p_a,
#     g = G(p_a, e)/p_a, n = n_exp (exact on the axis at every p; sheet 2.3); BA06: p0 = -p_rep, kappa_hat = p_rep/K(p_rep, e),
#     mu0 = G(p_rep, e) (sheet 15).  "material_flags" then carries NO energy flag.
#   * legacy (the 2026-10-02 format): "material_flags" is an ordered pass-through list that includes the energy's own flags.
# p_a is ONE number (the -p_a flag): the HAR reference pressure, the fork CSL's and the DM04 target's (TIMs: 101 kPa,
# sheet 2.4 / round 3b A4).  The deck owns what is a property of the BOUNDARY-VALUE PROBLEM: the initial state of each
# row (sigma0, v0, pi_i0).  The one place the deck is energy-aware is the p' floor bookkeeping (p_ref = |p0| under BA06,
# p_a under HAR, sheet 2.4): it decides which rows the floor touches at the initial state.
NS = {}
NS_RESERVED = {"-sigma0", "-pi0", "-pi0_auto", "-v0"}
NS_ENERGY_FLAGS = {"BA06": ("-p0", "-kappa_hat", "-eps_v0", "-mu0", "-alpha0"),
                   "HAR": ("-k", "-g", "-n", "-G0", "-nu", "-e_ref")}
P_MIN_FRAC = 5.0e-3                  # kernel default p_min = 5e-3 p_ref (sheet 9.7)
_NS_PROBE = [100000]
NS_TAG = {}                          # element tag -> row material tag (set by ns_make_materials)
NS_K0_SOURCE = [""]


@contextlib.contextmanager
def silenced():
    """Silence the C++ streams (fd 1 and 2) around constructor probes; Python's own prints are unaffected."""
    sys.stdout.flush()
    sys.stderr.flush()
    fd1, fd2 = os.dup(1), os.dup(2)
    dn = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(dn, 1)
        os.dup2(dn, 2)
        yield
    finally:
        os.dup2(fd1, 1)
        os.dup2(fd2, 2)
        for f in (dn, fd1, fd2):
            os.close(f)


def ns_fv(flags, name, default=None):
    """First value of flag `name` in a [[flag, value, ...], ...] list (default if absent)."""
    for fl in flags:
        if fl[0] == name:
            return fl[1] if len(fl) > 1 else True
    return default


def load_ns_params(path, k0_cli=None, energy_cli=None, pmin_cli=None, prep_cli=None, pi0_cli=None):
    """Read the parameter file, resolve the energy, the p' floor and K0.  Sets NS and (if needed) K0.
    K0 order: --k0 > the file's K0 (when not null) > the deck's mode default (campaign / toyoura_campaign: the
    campaign's nu* = 0.312885, K0 = nu*/(1 - nu*) = 0.4554, ref_campaign/footing_ab.py:111-112; Kimura mode: 0.5)."""
    global NS, K0
    raw = open(path, "rb").read()
    NS = json.loads(raw.decode("utf-8"))
    NS["_path"] = os.path.abspath(path)
    NS["_md5"] = hashlib.md5(raw).hexdigest()
    base = NS["material_flags"]
    for fl in base:
        if not isinstance(fl, list) or not fl or not str(fl[0]).startswith("-"):
            raise SystemExit(f"ns-params: bad material_flags entry {fl!r}")
        if fl[0] in NS_RESERVED:
            raise SystemExit(f"ns-params: {fl[0]} is set per row by the deck, not by the file")
    names = [fl[0] for fl in base]
    if len(set(names)) != len(names):
        raise SystemExit(f"ns-params: duplicate flags in material_flags: {sorted(n for n in set(names) if names.count(n) > 1)}")
    em = NS.get("elastic_map")
    flags_pmin = ns_fv(base, "-pmin")
    if em is not None:
        clash = [n for n in names if n == "-energy" or n in NS_ENERGY_FLAGS["BA06"] or n in NS_ENERGY_FLAGS["HAR"]]
        if clash:
            raise SystemExit(f"ns-params: elastic_map generates the energy flags; remove {clash} from material_flags")
        if ns_fv(base, "-p_a") is None:
            raise SystemExit("ns-params: elastic_map needs -p_a in material_flags (the ONE p_a: HAR energy, CSL, DM04 target)")
        energy = energy_cli or NS.get("energy")
        if energy not in ("BA06", "HAR"):
            raise SystemExit(f"ns-params: energy must be BA06 or HAR (got {energy!r}; --energy / the file's 'energy')")
        NS["_p_rep"] = None
        if energy == "BA06":
            pr = prep_cli if prep_cli is not None else em.get("p_rep")
            if pr is None or not float(pr) > 0.0:
                raise SystemExit("ns-params: BA06 from elastic_map needs p_rep > 0 (file elastic_map.p_rep or --p-rep)")
            NS["_p_rep"] = float(pr)
        NS["_energy_flags"] = None
    else:
        if prep_cli is not None:
            raise SystemExit("--p-rep needs a parameter file with an elastic_map")
        has_har = ("-energy" in names and str(ns_fv(base, "-energy")) == "HAR") or any(n in names for n in ("-k", "-g", "-G0"))
        implied = "HAR" if has_har else "BA06"
        if energy_cli is not None and energy_cli != implied:
            raise SystemExit(f"--energy {energy_cli} contradicts the energy flags of the file's material_flags ({implied}); "
                             "the kernel REFUSES the other energy's constants (sheet 2.4), so the deck stops here")
        energy = implied
        energy_names = set(NS_ENERGY_FLAGS["BA06"]) | set(NS_ENERGY_FLAGS["HAR"]) | {"-energy", "-pmin"}
        NS["_energy_flags"] = [list(fl) for fl in base if fl[0] in energy_names and fl[0] != "-pmin"]
        base = [fl for fl in base if fl[0] not in energy_names]
    if pi0_cli is not None:
        NS["pi0"] = dict(NS.get("pi0") or {}, mode=pi0_cli)
    NS["_base_flags"] = [list(fl) for fl in base if fl[0] != "-pmin"]
    NS["energy_eff"] = energy
    NS["_energy_explicit"] = energy_cli is not None
    pmin = pmin_cli if pmin_cli is not None else (NS.get("pmin") if NS.get("pmin") is not None else flags_pmin)
    if pmin is not None:
        pmin = float(pmin)
        if pmin < 0.0:
            raise SystemExit("p_min < 0 is refused by the kernel (sheet 9.7); 0 switches the floor off")
    NS["_pmin_flag"] = pmin
    if k0_cli is not None:
        K0 = float(k0_cli)
        NS_K0_SOURCE[0] = "--k0"
    elif NS.get("K0") is not None:
        K0 = float(NS["K0"])
        NS_K0_SOURCE[0] = "the parameter file's K0"
    else:
        NS_K0_SOURCE[0] = ("the campaign's nu* = %.6f, K0 = nu*/(1-nu*) (ref_campaign/footing_ab.py:111-112)" % NU_STAR
                           if SAND in ("campaign", "toyoura_campaign") else
                           "the Kimura-mode default 0.5 (Loukidis & Salgado 2011 series, via nu* = K0/(1+K0))")
    return NS


def _dm04_G(em, p_a, p, e):
    return em["G0"] * p_a * (2.97 - e) ** 2 / (1.0 + e) * (p / p_a) ** float(em.get("n_exp", 0.5))


def _dm04_K(em, p_a, p, e):
    nu = em["nu"]
    return 2.0 * (1.0 + nu) * _dm04_G(em, p_a, p, e) / (3.0 * (1.0 - 2.0 * nu))


def ns_final_flags(e_init):
    """The material's flag list for the deck's void ratio: material_flags + the energy flags + -pmin.
    -energy BA06 is emitted only when asked for explicitly (--energy BA06): BA06 is the material's default, so the list
    stays valid on a build that predates -energy."""
    flags = [list(fl) for fl in NS["_base_flags"]]
    en, em = NS["energy_eff"], NS.get("elastic_map")
    if em is not None:
        p_a = float(ns_fv(flags, "-p_a"))
        e_ref = float(em["e_rep"]) if em.get("e_rep") is not None else float(e_init)
        if en == "HAR":
            if not (0.0 <= float(em.get("n_exp", 0.5)) < 1.0):
                raise SystemExit("elastic_map.n_exp must satisfy 0 <= n < 1 for HAR")
            flags += [["-energy", "HAR"], ["-k", _dm04_K(em, p_a, p_a, e_ref) / p_a],
                      ["-g", _dm04_G(em, p_a, p_a, e_ref) / p_a], ["-n", float(em.get("n_exp", 0.5))]]
        else:
            pr = NS["_p_rep"]
            if NS["_energy_explicit"]:
                flags.append(["-energy", "BA06"])
            flags += [["-p0", -pr], ["-kappa_hat", pr / _dm04_K(em, p_a, pr, e_ref)], ["-mu0", _dm04_G(em, p_a, pr, e_ref)]]
    else:
        flags += [list(fl) for fl in NS["_energy_flags"]]
        if NS["_energy_explicit"] and ns_fv(flags, "-energy") is None:
            flags.append(["-energy", en])
    if NS["_pmin_flag"] is not None:
        flags.append(["-pmin", NS["_pmin_flag"]])
    return flags


def ns_pmin_eff(flags):
    """p_min the kernel will use: the -pmin flag, else 5e-3 p_ref (p_ref = |p0| under BA06, p_a under HAR; sheet 2.4, 9.7)."""
    pm = ns_fv(flags, "-pmin")
    if pm is not None:
        return float(pm)
    if NS["energy_eff"] == "HAR":
        pref = ns_fv(flags, "-p_a")
        if pref is None:
            raise SystemExit("HAR needs -p_a")
    else:
        pref = abs(float(ns_fv(flags, "-p0")))
    return P_MIN_FRAC * float(pref)


def _ns_args(tag, sig0, pi0, v0):
    """The nDMaterial command; pi0 None = the kernel's own unified rule (-pi0_auto, sheet 5.4 / S.53)."""
    a = ["LadrunoNorSand", int(tag)]
    for fl in NS["flags"]:
        a += [(float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v) for v in fl]
    a += ["-sigma0", *[float(s) for s in sig0]]
    a += ["-pi0_auto"] if pi0 is None else ["-pi0", float(pi0)]
    return a + ["-v0", float(v0)]


def _ns_try(tag, sig0, pi0, v0):
    try:
        with silenced():
            ops.nDMaterial(*_ns_args(tag, sig0, pi0, v0))
        return True
    except Exception:
        return False


def _ns_create(tag, sig0, pi0, v0):
    """Create the row material; on a refusal re-run it with the kernel's messages visible and stop with the command."""
    args = _ns_args(tag, sig0, pi0, v0)
    try:
        with silenced():
            ops.nDMaterial(*args)
    except Exception:
        log(f"LadrunoNorSand REFUSED the row material: nDMaterial {' '.join(str(a) for a in args)}")
        log("  (a build that predates round 3b has no -energy / -pmin / -pi0_auto: rows below the p' floor need -pi0_auto, "
            "--energy HAR needs -energy -k -g -n; use --pmin 0 on such a build only if the file has no HAR flags)")
        try:
            ops.nDMaterial(*args)           # again, unsilenced: the kernel's own message goes to stderr
        except Exception as exc:            # noqa: BLE001
            raise SystemExit(f"row material refused: {exc}")
        raise SystemExit("row material refused on the first try but accepted on the second: stopping")


def ns_unified_pi0(sig0, flags):
    """pi_i0 (negative) by the unified rule (S.53), sheet 5.4, for a K0 row sigma0 = [-K0 sv, -sv, -K0 sv, 0, 0, 0]:
    the yield surface through (p_init, eta*), eta* = max(eta_init, c2 M).  The row is a compression-corner state
    (theta = pi/3, zeta = 1), so eta_init = q/p = 3(1-K0)/(1+2K0) exactly; c2 = the smooth cap's upper bound
    (planar: c2, else c1; no cap: 0).  pi_i0 = p exp(eta*/M - 1) (N = 0) or p [(1-N)/(1 - eta* N/M)]^((1-N)/N).
    Returns (pi_i0, eta_init, eta_star).  Refuses what the kernel refuses (eta* >= M/N)."""
    sv, sh, sz = -sig0[1], -sig0[0], -sig0[2]
    if not (sv > 0.0 and abs(sh - sz) <= 1e-12 * sv and sh <= sv * (1.0 + 1e-12)):
        raise SystemExit(f"ns unified pi_i0: closed form is for compression-corner K0 rows (K0 <= 1), got sigma0 {sig0}")
    p, q = (sv + sh + sz) / 3.0, sv - sh
    M, N = float(ns_fv(flags, "-M")), float(ns_fv(flags, "-N"))
    cap = str(ns_fv(flags, "-cap", "none"))
    c1, c2 = ns_fv(flags, "-c1"), ns_fv(flags, "-c2")
    c2e = 0.0 if cap == "none" else (float(c2) if c2 is not None else float(c1))
    eta_init = q / p
    eta_star = max(eta_init, c2e * M)
    if N > 0.0 and eta_star * N / M >= 1.0:
        raise SystemExit(f"ns unified pi_i0: eta* = {eta_star:.4f} >= M/N = {M / N:.4f}: no surface passes through the state")
    ratio = math.exp(eta_star / M - 1.0) if N == 0.0 else ((1.0 - N) / (1.0 - eta_star * N / M)) ** ((1.0 - N) / N)
    return -p * ratio, eta_init, eta_star


def _ns_pi0(sig0, v0):
    """pi_i0 (negative) for a row, by the file's rule.  Returns (pi0 or None for -pi0_auto, info).
      unified   : (S.53) closed form; a row below the p' floor goes to the kernel's -pi0_auto (the rule needs the FLOORED state)
      auto      : the kernel's -pi0_auto for every row (the closed form is still computed as a shadow, for the cross-check)
      on_surface: the smallest |pi_i0| the constructor accepts x (1 + margin) (bisection on the refusal; the 2026-10-02 rule)
      ratio     : pi_i0 = -ratio p"""
    flags = NS["flags"]
    p = -(sig0[0] + sig0[1] + sig0[2]) / 3.0
    rule = NS["pi0"]
    mode = rule["mode"]
    pmin = NS["pmin_eff"]
    floored = pmin > 0.0 and p < pmin * (1.0 + 1e-9)
    info = dict(floored=bool(floored), closed=None, eta_init=None, eta_star=None, source="")
    if mode in ("unified", "auto") and not floored:
        info["closed"], info["eta_init"], info["eta_star"] = ns_unified_pi0(sig0, flags)
    if mode == "unified":
        if floored:
            info["source"] = "kernel -pi0_auto (p' below the floor: S.53 needs the floored state)"
            return None, info
        info["source"] = "S.53 closed form"
        return info["closed"], info
    if mode == "auto":
        info["source"] = "kernel -pi0_auto"
        return None, info
    if mode == "ratio":
        info["source"] = "ratio" + (" (of the floored p)" if floored else "")
        return -float(rule["ratio"]) * (max(p, pmin) if floored else p), info   # 2026-11-12: a row below the floor takes the floored p
    if mode != "on_surface":
        raise SystemExit(f"ns-params: unknown pi0 mode {mode!r} (unified | auto | on_surface | ratio)")

    def ok(r):
        _NS_PROBE[0] += 1
        return _ns_try(_NS_PROBE[0], sig0, -r * p, v0)
    lo, hi = 0.02, 50.0
    if not ok(hi):
        raise SystemExit(f"ns-params: even |pi_i0| = {hi} p is refused for sigma0 {sig0}, v0 {v0}")
    if ok(lo):
        r_s = lo
    else:
        for _ in range(45):
            mid = math.sqrt(lo * hi)
            if ok(mid):
                hi = mid
            else:
                lo = mid
        r_s = hi
    info["source"] = "on_surface bisection"
    info["ratio_surface"] = r_s
    return -r_s * (1.0 + float(rule.get("margin", 0.02))) * p, info


NGP = 4


def ns_make_materials(elems, y, einit):
    """One LadrunoNorSand per element ROW (tags 1..nrows), initial state sigma0 = gamma' z_c [-K0 -1 -K0 0 0 0]."""
    global NS_TAG
    e = NS.get("e_init")
    e = float(einit) if e is None else float(e)
    v0 = 1.0 + e
    NS["flags"] = ns_final_flags(e)
    NS["pmin_eff"] = ns_pmin_eff(NS["flags"])
    rows, tag_of_row = [], {}
    for tag, i, j in elems:
        if j not in tag_of_row:
            tag_of_row[j] = len(tag_of_row) + 1
    ns_tag = {tag: tag_of_row[j] for tag, i, j in elems}
    NS_TAG = ns_tag
    for j in sorted(tag_of_row):
        z = -0.5 * (y[j] + y[j + 1])                       # depth of the row centroid (m)
        sv = GAMMA * z                                      # vertical effective stress (kPa, compression positive)
        sig0 = [-K0 * sv, -sv, -K0 * sv, 0.0, 0.0, 0.0]
        pi0, info = _ns_pi0(sig0, v0)
        mt = tag_of_row[j]
        _ns_create(mt, sig0, pi0, v0)
        p = sv * (1.0 + 2.0 * K0) / 3.0
        rows.append(dict(row=j, mat=mt, z=z, sigma_v=sv, p=p, pi_i0=pi0, ratio_pi_p=(-pi0 / p if pi0 is not None else None),
                         ratio_surface=info.get("ratio_surface"), pi_i0_source=info["source"], floored_init=info["floored"],
                         pi_i0_closed=info["closed"], eta_init=info["eta_init"], eta_star=info["eta_star"]))
    nfl = sum(1 for r in rows if r["floored_init"])
    rr = [r["ratio_pi_p"] for r in rows if r["ratio_pi_p"] is not None]
    mode = NS["pi0"]["mode"]
    desc = (f"LadrunoNorSand per row: {len(rows)} rows, e_init {e:.5f} (v0 {v0:.5f}), K0 {K0:.5f} ({NS_K0_SOURCE[0]}), "
            f"energy {NS['energy_eff']}, p_min_eff {NS['pmin_eff']:.4f} kPa ({nfl} rows below it at the initial state), "
            f"flags {NS['flags']}, pi0 rule {mode}"
            + (f" margin {NS['pi0'].get('margin')}" if mode == 'on_surface' else "")
            + f"; row p' {rows[0]['p']:.4f} .. {rows[-1]['p']:.3f} kPa, "
            + (f"pi_i0/p {min(rr):.4f} .. {max(rr):.4f}" if rr else "pi_i0 from the kernel")
            + (f" ({len(rows) - len(rr)} rows from the kernel)" if rr and len(rr) < len(rows) else ""))
    return ns_tag, rows, desc


def ns_check_pi0(deck, Fz):
    """Cross-check of the initial pi_i: the committed state's pi_i against what the deck asked for.
      unified rows : the state's pi_i = the S.53 closed form passed with -pi0 (constructor round trip)
      auto rows    : the kernel's own -pi0_auto against the deck's shadow closed form (the kernel-vs-deck check of S.53)
      floored rows : kernel only (reported as pi_i0/p)
    Returns a dict and logs one line."""
    rows = {r["mat"]: r for r in deck["ns_rows"]}
    st = Fz["st"]
    dev_pass, dev_kernel, fl_ratio = [], [], []
    for k, tag in enumerate(deck["tags"]):
        r = rows[NS_TAG[tag]]
        for gp in range(NGP):
            pi = st[4 * k + gp, 0]
            if r["pi_i0"] is not None:
                dev_pass.append(abs(pi / r["pi_i0"] - 1.0))
            elif r["pi_i0_closed"] is not None:
                dev_kernel.append(abs(pi / r["pi_i0_closed"] - 1.0))
            else:
                fl_ratio.append(-pi / r["p"])
    out = dict(n_rows_passed=sum(1 for r in rows.values() if r["pi_i0"] is not None),
               max_dev_passed=(float(max(dev_pass)) if dev_pass else None),
               n_rows_kernel_auto=sum(1 for r in rows.values() if r["pi_i0"] is None and r["pi_i0_closed"] is not None),
               max_dev_kernel_vs_closed=(float(max(dev_kernel)) if dev_kernel else None),
               n_rows_floored_kernel=sum(1 for r in rows.values() if r["pi_i0"] is None and r["pi_i0_closed"] is None),
               floored_ratio_range=([float(min(fl_ratio)), float(max(fl_ratio))] if fl_ratio else None))
    log(f"pi_i0 check: passed-in rows {out['n_rows_passed']} (committed pi_i vs the value passed: max rel dev "
        f"{out['max_dev_passed']}); kernel -pi0_auto rows {out['n_rows_kernel_auto']} (vs the deck's S.53 closed form: max rel dev "
        f"{out['max_dev_kernel_vs_closed']}); floored kernel-only rows {out['n_rows_floored_kernel']} (pi_i0/p range "
        f"{out['floored_ratio_range']})")
    return out


_FLOOR_OK = [None]


def ns_floor_census(deck):
    """Floor counters summed over the Gauss points: responses 'floor' [at_floor, n_f_tr, n_f_post, eps_f_v, W_f] and
    'floorInit' [n_f_init].  None on a build whose material has no floor responses (pre-round-3)."""
    if _FLOOR_OK[0] is False:
        return None
    tot = np.zeros(5)
    n_init = 0
    try:
        for tag in deck["tags"]:
            for gp in range(1, NGP + 1):
                v = ops.eleResponse(tag, "material", gp, "floor")
                if v is None or len(v) < 5:
                    raise RuntimeError("no floor response")
                tot += np.asarray(v[:5], float)
                n_init += int(ops.eleResponse(tag, "material", gp, "floorInit")[0])
    except Exception:
        _FLOOR_OK[0] = False
        return None
    _FLOOR_OK[0] = True
    return dict(n_at_floor=int(tot[0]), n_f_tr=int(tot[1]), n_f_post=int(tot[2]), eps_f_v=float(tot[3]), W_f=float(tot[4]),
                n_f_init=n_init)


def read_field_ns(deck, full=False):
    """NorSand committed field.  Stress (+ sigma_33) from the element's stressesPlaneStrain every call; strain,
    state [pi_i psi_i v v0 eps_p_v eps_p_s] and psi per Gauss point only when full (checkpoints)."""
    tags = deck["tags"]
    n = 4 * len(tags)
    sig = np.full((n, 3), np.nan)
    s33 = np.full(n, np.nan)
    eps = np.zeros((n, 3))
    st = np.full((n, 6), np.nan)
    psi = np.full(n, np.nan)
    k = 0
    for tag in tags:
        sp = ops.eleResponse(tag, "stressesPlaneStrain")        # [s11 s22 s12 s33] per Gauss point
        if len(sp) != 4 * NGP:
            raise SystemExit(f"stressesPlaneStrain: {len(sp)} values on element {tag}, want {4 * NGP}")
        for gp in range(4):
            if gp >= NGP:          # ELEM with fewer Gauss points: the slot stays NaN
                k += 1
                continue
            sig[k] = (sp[4 * gp], sp[4 * gp + 1], sp[4 * gp + 2])
            s33[k] = sp[4 * gp + 3]
            if full:
                eps[k] = ops.eleResponse(tag, "material", gp + 1, "strain")
                st[k] = ops.eleResponse(tag, "material", gp + 1, "state")
                psi[k] = ops.eleResponse(tag, "material", gp + 1, "psi")[0]
            k += 1
    return dict(sig=sig, s33=s33, eps=eps, st=st, psi=psi, f=np.full(n, np.nan), stats=None)


def derived_ns(F):
    s11, s22, s12, s33 = F["sig"][:, 0], F["sig"][:, 1], F["sig"][:, 2], F["s33"]
    pm = (s11 + s22 + s33) / 3.0
    p = -pm
    q = np.sqrt(1.5 * ((s11 - pm) ** 2 + (s22 - pm) ** 2 + (s33 - pm) ** 2 + 2.0 * s12 ** 2))
    eta = np.where(p > 1e-12, q / np.maximum(p, 1e-300), np.nan)
    return dict(p=p, q=q, eta=eta)


def ns_ring_lines(deck, D, F, top=5):
    xy = gp_xy(deck)
    out = ["    lowest p':"]
    for k in np.argsort(D["p"])[:top]:
        t, gx, gy = xy[k]
        out.append(f"    p' {D['p'][k]:9.4f} eta {D['eta'][k]:7.3f} psi {F['psi'][k]:+.4f} ele {t} gp {k%4+1} "
                   f"x {gx:+.3f} y {gy:+.3f}")
    out.append("    highest eta:")
    for k in np.argsort(-np.nan_to_num(D["eta"], nan=-1))[:top]:
        t, gx, gy = xy[k]
        out.append(f"    eta {D['eta'][k]:7.3f} p' {D['p'][k]:9.4f} psi {F['psi'][k]:+.4f} ele {t} gp {k%4+1} "
                   f"x {gx:+.3f} y {gy:+.3f}")
    return "\n".join(out)


def ns_latched(deck):
    """(Gauss points whose material has latched a refusal, total refused trials), response 'refusal'."""
    nl, nref = 0, 0
    for tag in deck["tags"]:
        for gp in range(1, NGP + 1):
            r = ops.eleResponse(tag, "material", gp, "refusal")
            nl += int(r[2] > 0)
            nref += int(r[1])
    return nl, nref


# ------------------------------------------------------------------- run ----
def main(argv=None):
    _pre = sand_parser()
    _pa, _ = _pre.parse_known_args(argv)
    if _pa.kimura is not None:
        if _pa.sand != "toyoura":
            raise SystemExit("--kimura needs --sand toyoura (the Kimura mode: B 0.9 m, no surcharge, no footing weight)")
        resolve_kimura(_pa)
    if _pa.sand == "toyoura":
        set_toyoura(_pa)
    elif _pa.sand == "toyoura_campaign":
        set_toyoura_campaign(_pa)
    elif any(v is not None for k, v in vars(_pa).items() if k != "sand"):
        raise SystemExit("--bproto/--dr/--k0/--gamma/--dom-w/--dom-d/--nu-push/--fine-depth/--kimura need --sand toyoura")
    ap = argparse.ArgumentParser(parents=[_pre])
    ap.add_argument("--mat", choices=("sanisand", "dp", "norsand"), default="sanisand")
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--extra", default="", help="extra material tokens, space separated")
    ap.add_argument("--maxsub", type=int, default=2000)
    ap.add_argument("--tantype", type=int, default=0,
                    help="positional TanType (0 = elastic Ce, the campaign's; 1 = continuum)")
    ap.add_argument("--tolr", type=float, default=1.0e-7,
                    help="positional TolR (IS the substep tolerance under IntScheme 129)")
    ap.add_argument("--dp-phi", type=float, default=38.0)
    ap.add_argument("--dp-g", type=float, default=30000.0)
    ap.add_argument("--target", type=float, default=0.15, help="s/B target")
    ap.add_argument("--wall", type=float, default=36000.0)
    ap.add_argument("--ckpt-every", type=int, default=25)
    ap.add_argument("--field-every-step", type=int, default=1,
                    help="read the full field every step (1) or only at checkpoints (0)")
    ap.add_argument("--deterministic", type=int, default=0)
    ap.add_argument("--ds-min", type=float, default=DS_MIN)
    ap.add_argument("--probe-test-step", type=int, default=0,
                    help="TEST ONLY: stop after N steps and run the wall post-mortem")
    ap.add_argument("--krylov-tol-mult", type=float, default=10.0,
                    help="tolerance multiplier of the KrylovNewton rung (10 = the "
                         "original ladder; 1 = the same 1e-5 x applied as Newton)")
    ap.add_argument("--mesh", choices=("b4", "b8", "b16"), default="b8",
                    help="element size in the fine band: B/8 (default, the act's mesh), B/16, or B/4 (WP-150 R2)")
    ap.add_argument("--mesh-perturb", default="",
                    help="WP-150 R2 mesh-orientation variant: jitter:A[:SEED] or shear:DEG (interior fine-zone nodes)")
    ap.add_argument("--presidual", type=float, default=0.0,
                    help="LadrunoSANISAND -Presidual (kPa); 0 = the campaign's (cohesionless)")
    ap.add_argument("--einit", type=float, default=E_INIT,
                    help="initial void ratio e_init (positional 3rd SANISAND param); "
                         "default = the campaign's 0.6944")
    ap.add_argument("--zmax", type=float, default=SAN[15], help="fabric z_max (campaign 12.5)")
    ap.add_argument("--nb", type=float, default=SAN[12], help="n^b (campaign 3.5)")
    ap.add_argument("--nd", type=float, default=SAN[14], help="n^d (campaign 5.75)")
    ap.add_argument("--A0", type=float, default=SAN[13], help="A0 dilatancy (campaign 0.05)")
    ap.add_argument("--h0", type=float, default=SAN[10], help="h0 hardening (campaign 1.3)")
    ap.add_argument("--c", type=float, default=SAN[4], help="Lode ratio c (campaign 0.71; 7/9 convexity)")
    ap.add_argument("--ns-params", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                        "ns_params_toyoura.json"),
                    help="norsand: JSON parameter file (material flags pass-through + initial-state rules)")
    ap.add_argument("--energy", choices=("BA06", "HAR"), default=None,
                    help="norsand: the elastic energy, passed to the material as -energy (pass-through). With the file's "
                         "elastic_map the DM04 targets are mapped onto it at the deck's e_init; with literal energy flags it "
                         "must agree with them (the kernel refuses the other energy's constants, sheet 2.4)")
    ap.add_argument("--pmin", type=float, default=None,
                    help="norsand: the p' floor, passed to the material as -pmin (kPa; 0 = off; default: the kernel's "
                         "5e-3 p_ref, p_ref = |p0| for BA06, p_a for HAR; sheet 9.7)")
    ap.add_argument("--p-rep", type=float, default=None,
                    help="norsand BA06 from elastic_map: the representative pressure (kPa) at which the constant mu0, "
                         "kappa_hat match the DM04 targets (default: the file's elastic_map.p_rep)")
    ap.add_argument("--pi0", choices=("unified", "auto", "on_surface", "ratio"), default=None,
                    help="norsand: override the file's pi_i0 rule (unified = S.53 closed form in the deck, auto = the kernel's "
                         "-pi0_auto for every row with the closed form as the cross-check, on_surface = the 2026-10-02 bisection)")
    ap.add_argument("--init-only", action="store_true",
                    help="norsand: build the materials and elements, read and check the initial field, then stop")
    ap.add_argument("--system", choices=("Pardiso", "UmfPack", "SparseGeneral"), default="Pardiso",
                    help="linear system. Pardiso needs the pdl_seq build (falls back to UmfPack if refused). UmfPack was measured to fail the first push step on the plain build (symbolic analysis -8, 2026-10-02), SparseGeneral is the plain-build alternative")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    set_mesh(args.mesh)
    ns = args.mat == "norsand"
    if ns:
        load_ns_params(args.ns_params, k0_cli=_pa.k0, energy_cli=args.energy, pmin_cli=args.pmin, prep_cli=args.p_rep,
                       pi0_cli=args.pi0)
    elif (args.energy is not None or args.pmin is not None or args.p_rep is not None or args.init_only
          or args.pi0 is not None):
        raise SystemExit("--energy/--pmin/--p-rep/--pi0/--init-only need --mat norsand")
    global MESH_PERTURB
    MESH_PERTURB = args.mesh_perturb
    global LADDER
    LADDER = tuple((a, (args.krylov_tol_mult if a == "KrylovNewton" else tf), it)
                   for a, tf, it in LADDER)

    out = args.out
    for d in ("ckpt", "replay"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    log(f"engine {ops.__file__}")
    log(f"build {ops.ladrunoBuild().strip().splitlines()[0]}")
    log(f"env MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS')} "
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')} "
        f"MKL_CBWR={os.environ.get('MKL_CBWR')}")
    log(f"args {vars(args)}")
    if SAND in ("toyoura", "toyoura_campaign"):
        log(f"GATE 1 sand {SAND}: {TYR_INFO}")
    if KIMURA_INFO:
        log(f"kimura test {KIMURA_INFO['test_id']}: {KIMURA_INFO}")
    extra = []
    for tok in args.extra.split():
        try:
            extra.append(float(tok) if any(c in tok for c in ".e") and not tok.startswith("-") else int(tok))
        except ValueError:
            extra.append(tok)
    san = args.mat == "sanisand"
    global CENSUS
    if san and args.scheme == 129:
        CENSUS = CENSUS_SAS
    deck = build(args.mat, args.scheme, extra, args.maxsub, args.dp_phi, args.dp_g,
                 bool(args.deterministic), args.tolr, args.tantype,
                 presidual=args.presidual, einit=args.einit,
                 over={15: args.zmax, 12: args.nb, 14: args.nd, 13: args.A0, 10: args.h0, 4: args.c},
                 system=args.system, init_only=args.init_only)
    global NB, ND
    NB, ND = args.nb, args.nd          # derived() diagnostics follow the ablation
    log(f"ablation knobs: zmax = {args.zmax}, nb = {args.nb}, nd = {args.nd}, "
        f"A0 = {args.A0}, h0 = {args.h0} (campaign 12.5 / 3.5 / 5.75 / 0.05 / 1.3)")
    log(f"lode c = {args.c} (campaign 0.71)")
    try:
        _so = ops.eleResponse(ops.getEleTags()[0], "material", 1, "sasOptions")
        log(f"sasOptions (ele {ops.getEleTags()[0]} gp 1) = {list(_so)}")
    except Exception as _exc:
        log(f"sasOptions unavailable: {_exc}")
    log(f"ladder knobs: Presidual = {args.presidual} kPa, e_init = {args.einit} "
        f"(campaign 0.6944). Initial stress is e_init-INDEPENDENT: gamma' = {GAMMA} "
        f"kN/m3 fixed, K0 = {K0:.4f} (" + (NS_K0_SOURCE[0] if ns else "nu*/(1-nu*), stage-0 elastic") + "); "
        f"e_init changes only the elastic moduli magnitude G(e,p) and psi")
    log(f"material: {deck['matdesc']}")
    log(f"mesh {MESH} (fine h = B/{B_FOOT/H0:.0f} = {H0:.5f} m): {len(deck['tags'])} elements, {4*len(deck['tags'])} Gauss points, "
        f"solver {deck['solver']}, applied vertical load {deck['applied']:.3f} kN/m")

    F = read_field(deck, want_stats=san, want_f=(san or ns))
    D = derived(F) if san else (derived_ns(F) if ns else None)
    if ns:
        save_field(os.path.join(out, "ckpt", "field_flip.npz"), deck, F, D, 0, 0.0)
        log(f"norsand initial field: min p' {D['p'].min():.4f} kPa, max eta {np.nanmax(D['eta']):.4f}, "
            f"min psi {np.nanmin(F['psi']):+.4f}, max psi {np.nanmax(F['psi']):+.4f}; "
            f"Gauss points with p' < {NS['p_floor_report_kPa'][0]} kPa: {int(np.sum(D['p'] < NS['p_floor_report_kPa'][0]))}")
        for r in (deck["ns_rows"][0], deck["ns_rows"][len(deck["ns_rows"]) // 2], deck["ns_rows"][-1]):
            pr_txt = f"pi_i0 {r['pi_i0']:10.4f} (pi_i0/p {r['ratio_pi_p']:.4f})" if r["pi_i0"] is not None else "pi_i0 by the kernel"
            log(f"  row {r['row']:3d} (mat {r['mat']:3d}) z {r['z']:8.4f} m sigma_v {r['sigma_v']:9.4f} p' {r['p']:9.4f} "
                f"{pr_txt} [{r['pi_i0_source']}; eta_init {r['eta_init']}, eta* {r['eta_star']}; on-surface ratio {r['ratio_surface']}]")
        if args.init_only:
            log("DONE (--init-only)")
            return 0
    if san:
        log(f"at the flip (stage 0 end): max eta/Mc = {np.nanmax(D['eta'])/MC:.4f}, "
            f"min p' = {D['p'].min():.4f} kPa, max |f_read - f_recomputed| = "
            f"{np.nanmax(np.abs(F['f'] - D['frec'])):.3e}")
        save_field(os.path.join(out, "ckpt", "field_flip.npz"), deck, F, D, 0, 0.0)

        if SAND in ("toyoura", "toyoura_campaign"):
            # GATE 1: the initial state down the centreline and in the far field
            xyp = gp_xy(deck)
            gx = np.array([g[1] for g in xyp]); gy = np.array([g[2] for g in xyp])
            for xlab, x0 in (("centre", 0.0), ("far", 0.75 * XLIM)):
                for zB in (0.05, 0.25, 0.5, 1.0, 2.0, 4.0, 0.9 * (-YBOT) / B_FOOT):
                    k = int(np.argmin(np.abs(gx - x0) + np.abs(gy + zB * B_FOOT)))
                    sv, sh = -F["sig"][k][1], -F["sig"][k][0]
                    log(f"  initial {xlab:6s} z {-gy[k]:7.3f} m (z/B {-gy[k]/B_FOOT:6.3f}) "
                        f"x {gx[k]:+7.3f}: sig_v' {sv:9.3f} sig_h' {sh:9.3f} K {sh/sv:6.4f} "
                        f"(gamma z {GAMMA*(-gy[k]):9.3f}) p' {D['p'][k]:9.3f} e {D['e'][k]:.5f} "
                        f"psi {F['psi'][k]:+.4f}")
    if not ns:
        ops.updateMaterialStage("-material", MAT, "-stage", 1)
    if SAND in ("toyoura", "toyoura_campaign") and san and NU_PUSH is not None:
        # GATE 1: nu* only served K0 in the elastic stage; the push runs DM04's nu
        ops.setParameter("-val", float(NU_PUSH), "-ele", *deck["tags"], "poissonRatio", str(MAT))
        log(f"GATE 1: poissonRatio {NU_STAR:.6f} -> {NU_PUSH} on {len(deck['tags'])} elements "
            f"after the stage flip (checked from the elastic tangent after push step 1)")

    # ---- push setup --------------------------------------------------------
    REF = deck["REF"]
    uy0 = ops.nodeDisp(REF, 2)
    ops.loadConst("-time", uy0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.sp(REF, 2, 1.0)
    ptol = PUSH_TOL_REL * deck["applied"]
    log(f"push: u_y0(ref) = {uy0:.6e} m, NormUnbalance tol {ptol:.4e} kN "
        f"(1e-5 x applied), ladder {LADDER}")

    def q_now():
        ops.reactions()
        return (W_FOOT - ops.nodeReaction(REF, 2)) / B_FOOT

    fh = open(os.path.join(out, "steps.csv"), "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["step", "s_m", "s_over_B", "q_kPa", "ds_m", "rung", "iters",
                 "fails_before", "wall_step_s", "wall_total_s",
                 "sub_step_total", "sub_step_maxpt", "cap_step", "rejErr_step",
                 "dtmin_step", "lowp_rej_step", "clampMc_step",
                 "min_p", "max_eta", "max_rho_alpha", "n_rho_gt_1",
                 "n_p_lt_F", "min_psi", "n_latched"])      # last three: norsand only (blank otherwise)
    smax = args.target * B_FOOT
    ds, good, nstep, nfail_tot = DS_BASE, 0, 0, 0
    mode = "TARGET"
    prev_stats = F["stats"].copy() if san else None
    Fprev, Dprev = F, D                 # committed state of the previous step
    Fpp = Dpp = None                    # ... and of the one before it
    sub_prev = cap_prev = eps_ppp = None
    sB_prev = sB_pp = 0.0
    dt_prev = 0.0
    pending = []                        # checkpoint (F, D, step, sB) awaiting dS
    t_push = time.time()
    rows = []
    first_try_ds = ds
    q0 = q_now()
    log(f"q at push start = {q0:.4f} kPa (dead load / B = {W_FOOT/B_FOOT:.4f})")
    # WP-152 review: an opt-in EnergyBalance recorder (FOOTING_ENERGY=1) on the push
    if os.environ.get("FOOTING_ENERGY"):
        ops.recorder("EnergyBalance", "-file", os.path.join(out, "energy_balance.txt"), "-time")
        log("EnergyBalance recorder -> energy_balance.txt (KE IE DW ULW RES ERR)")
    while True:
        s_now = uy0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        # a running leg's wall cap can be RAISED (never lowered below elapsed)
        # by writing seconds into <out>/WALL_OVERRIDE
        try:
            args.wall = max(args.wall, float(open(os.path.join(out, "WALL_OVERRIDE")).read()))
        except (OSError, ValueError):
            pass
        if time.time() - T0 > args.wall:
            mode = "WALL"
            break
        if args.probe_test_step and nstep >= args.probe_test_step:
            mode = "PROBETEST"
            break
        ds = min(ds, smax - s_now)
        fails_before = 0
        first_try_ds = ds
        t_step = time.time()
        while True:
            ops.integrator("LoadControl", -ds)
            ok, rung, iters = False, -1, 0
            for r, (algo, tf, it) in enumerate(LADDER):
                ops.test("NormUnbalance", tf * ptol, it, 0)
                ops.algorithm(algo)
                if ops.analyze(1) == 0:
                    ok, rung, iters = True, r, ops.testIter()
                    break
            if ok:
                break
            fails_before += 1
            nfail_tot += 1
            good = 0
            ds *= 0.5
            log(f"  step {nstep+1}: ladder failed, ds -> {ds:.4e} m "
                f"(s/B {s_now/B_FOOT:.6f})")
            if ds < args.ds_min or time.time() - T0 > args.wall:
                break
        if not ok:
            mode = "FLOOR" if ds < args.ds_min else "WALL"
            break
        nstep += 1
        good += 1
        wall_step = time.time() - t_step
        s = uy0 - ops.getTime()
        sB = s / B_FOOT
        q = q_now()
        # ---- field + census ----------------------------------------------
        ck = (nstep % args.ckpt_every == 0)
        cen = [0, 0, 0, 0, 0, 0, 0]
        minp = maxeta = maxrho = float("nan"); nrho = 0
        if san and (args.field_every_step or ck):
            Fn = read_field(deck, want_stats=True, want_f=True)
            Dn = derived(Fn)
            dst = Fn["stats"] - prev_stats
            prev_stats = Fn["stats"].copy()
            C = CENSUS
            sub_pt = dst[:, C["sub"]]
            cen = [int(sub_pt.sum()), int(sub_pt.max()), int(dst[:, C["cap"]].sum()),
                   int(dst[:, C["rej"]].sum()), int(dst[:, C["dtmin"]].sum()),
                   int(dst[:, C["lowp"]].sum()), int(dst[:, C["clamp"]].sum())]
            if C is CENSUS_SAS and cen[2] > 0:
                # refusals this step (all attempts, failed ones included: the
                # census survives revertToLastCommit), by code and by location
                xy = gp_xy(deck)
                bycode = {SAS_REF_CODES[i]: int(dst[:, i].sum()) for i in SAS_REF_CODES
                          if dst[:, i].sum() > 0}
                ks = np.argsort(-dst[:, C["cap"]])[:6]
                where = "; ".join(f"ele {xy[k][0]} gp {k%4+1} ({xy[k][1]:+.2f},"
                                  f"{xy[k][2]:+.2f}) x{int(dst[k, C['cap']])}"
                                  for k in ks if dst[k, C["cap"]] > 0)
                log(f"  refusals step {nstep}: {cen[2]} by code {bycode}; "
                    f"n points {int(np.sum(dst[:, C['cap']] > 0))}; top: {where}")
            if C is CENSUS_SAS:
                # DEAD-END census (WP-129 review finding): a point whose committed
                # rho_alpha > 1 + kappa refuses EVERY later increment with
                # startAlphaOutsideBounding, and no step cut can help. Flag points
                # whose refusals this step are ALL refStartAlpha, and every point
                # whose COMMITTED rho_alpha already exceeds 1.1.
                de = np.where((dst[:, 15] > 0) & (dst[:, 17] == dst[:, 15]))[0]
                over = np.where(Dn["rho"] > 1.1)[0]
                xy = gp_xy(deck)
                for k in list(de) + list(over):
                    if k not in DEADEND:
                        DEADEND[k] = nstep
                        log(f"  DEAD-END candidate first seen step {nstep} (s/B {sB:.6f}): "
                            f"ele {xy[k][0]} gp {k%4+1} ({xy[k][1]:+.3f},{xy[k][2]:+.3f}) "
                            f"p' {Dn['p'][k]:.3f} rho_a {Dn['rho'][k]:.4f} "
                            f"refStartAlpha this step {int(dst[k, 17])} of "
                            f"{int(dst[k, 15])} refusals")
                if len(de) or len(over):
                    log(f"  dead-end census step {nstep}: all-startAlpha refusers {len(de)}, "
                        f"committed rho_a>1.1 {len(over)}, cumulative flagged {len(DEADEND)}")
            minp, maxeta = float(Dn["p"].min()), float(np.nanmax(Dn["eta"]))
            maxrho = float(np.nanmax(Dn["rho"]))
            nrho = int(np.sum(Dn["rho"] > 1.0))
            # replay rows: committed state of the PREVIOUS step + this step's dS
            for (Fc, Dc, stc, sbc, epsc) in pending:
                sel = select_worst(Dc, Fc, sub_pt)
                write_replay(os.path.join(out, "replay", f"replay_step{stc:05d}.csv"),
                             deck, Fc, Dc, Fn, Dn, epsc, sel, stc, sbc, sub_pt,
                             dst[:, CENSUS["cap"]], -ds)
            pending = []
            if ck:
                save_field(os.path.join(out, "ckpt", f"field_step{nstep:05d}.npz"),
                           deck, Fn, Dn, nstep, sB)
                pending.append((Fn, Dn, nstep, sB, Fprev["eps"]))
                log(f"  checkpoint step {nstep} s/B {sB:.6f}: max|f_read - f_rec| = "
                    f"{np.nanmax(np.abs(Fn['f'] - Dn['frec'])):.3e}\n"
                    + ring_summary(deck, Dn, Fn))
            eps_ppp = Fpp["eps"] if Fpp is not None else None
            Fpp, Dpp, sB_pp = Fprev, Dprev, sB_prev
            Fprev, Dprev, sub_prev, cap_prev, dt_prev = (Fn, Dn, sub_pt,
                                                         dst[:, CENSUS["cap"]], -ds)
            sB_prev = sB
        ns_cols = ["", "", ""]
        if ns and (args.field_every_step or ck):
            Fn = read_field_ns(deck, full=ck)
            Dn = derived_ns(Fn)
            minp, maxeta = float(Dn["p"].min()), float(np.nanmax(Dn["eta"]))
            nfl = int(np.sum(Dn["p"] < NS["p_floor_report_kPa"][0]))
            ns_cols = [nfl, "", ""]
            if ck:
                nlat, nref = ns_latched(deck)
                ns_cols = [nfl, f"{np.nanmin(Fn['psi']):.5g}", nlat]
                save_field(os.path.join(out, "ckpt", f"field_step{nstep:05d}.npz"), deck, Fn, Dn, nstep, sB)
                log(f"  checkpoint step {nstep} s/B {sB:.6f}: latched points {nlat}, refused trials {nref}, "
                    f"floor census {ns_floor_census(deck)}\n"
                    + ns_ring_lines(deck, Dn, Fn))
        if san and args.field_every_step and CENSUS is CENSUS_SAS:
            # WP-152 review: the separation census of this COMMITTED step (per GP,
            # cumulative columns), with where the separated points sit
            _nm = CENSUS["names"]
            _S = Fn["stats"]
            _act = _S[:, _nm.index("sepActive")] > 0.5
            _xy = np.array([(gx, gy) for _t, gx, gy in gp_xy(deck)])
            _ytop = _xy[:, 1].max()
            _yrow = np.unique(np.round(_xy[:, 1], 6))
            _ytop2 = _yrow[-2] if len(_yrow) > 1 else _ytop       # the second GP row from the top
            _foot = np.abs(_xy[:, 0]) <= 0.5 * B_FOOT
            _nsep = int(_act.sum())
            _row = [nstep, f"{sB:.9e}", f"{q:.6f}", _nsep, int((_act & _foot).sum()),
                    int((_act & (_xy[:, 1] < _ytop2 - 1e-9)).sum()),
                    f"{(_xy[_act, 1].min() if _nsep else float('nan')):.4f}",
                    f"{(np.abs(_xy[_act, 0]).min() if _nsep else float('nan')):.4f}",
                    f"{(np.abs(_xy[_act, 0]).max() if _nsep else float('nan')):.4f}"]
            _row += [int(_S[:, _nm.index(c)].sum()) for c in
                     ("sepEntriesTension", "sepEntriesLowP", "sepExits",
                      "sepHeldHighP", "sepHeldCompressing")]
            _row += [f"{_S[:, _nm.index('sepMaxP0')].max():.6f}",
                     "+".join(f"{c}:{int(n)}" for c, n in zip(*np.unique(
                         _S[_S[:, _nm.index('sepLastCode')] > 0, _nm.index('sepLastCode')],
                         return_counts=True)))]
            _new = not os.path.exists(os.path.join(out, "sep_census.csv"))
            with open(os.path.join(out, "sep_census.csv"), "a") as _fh:
                if _new:
                    _fh.write("step,s_over_B,q_kPa,n_sep,n_sep_footprint,n_sep_below_2nd_row,"
                              "y_min_sep,absx_min_sep,absx_max_sep,entriesTension,entriesLowP,"
                              "exits,heldHighP,heldCompressing,maxP0_at_entry,lastCodes\n")
                _fh.write(",".join(str(v) for v in _row) + "\n")
            if _nsep:
                log(f"  separation step {nstep}: {_nsep} GP separated ({_row[4]} in the footprint, "
                    f"{_row[5]} below the top two GP rows), y_min {_row[6]}, |x| {_row[7]}..{_row[8]}")
        wr.writerow([nstep, f"{s:.9e}", f"{sB:.9e}", f"{q:.6f}", f"{ds:.4e}", rung,
                     iters, fails_before, f"{wall_step:.3f}",
                     f"{time.time()-t_push:.1f}", *cen,
                     f"{minp:.5g}", f"{maxeta:.5g}", f"{maxrho:.5g}", nrho, *ns_cols])
        fh.flush()
        rows.append((sB, q))
        if SAND in ("toyoura", "toyoura_campaign") and san and nstep == 1:
            try:
                for _tg in (deck["tags"][0], deck["tags"][-1]):
                    _C = np.asarray(ops.eleResponse(_tg, "material", 1, "tangent"), float)
                    _r = _C[1] / _C[0]
                    log(f"GATE 1: ele {_tg} gp 1 tangent C12/C11 = {_r:.6f} -> nu = "
                        f"{_r/(1+_r):.6f} (TanType 0 = elastic Ce; want {NU_PUSH if NU_PUSH is not None else NU_STAR})")
            except Exception as _exc:
                log(f"GATE 1: tangent read failed: {_exc}")
        log(f"step {nstep:5d} s/B {sB:.6f} q {q:9.3f} kPa ds {ds:.3e} rung "
            f"{'NLK'[rung]} it {iters:3d} fails {fails_before} wall {wall_step:8.2f}s"
            + (f" | sub {cen[0]} maxpt {cen[1]} cap {cen[2]} rejErr {cen[3]} "
               f"dtmin {cen[4]} | p'min {minp:.3f} eta_max {maxeta:.3f} "
               f"rho_max {maxrho:.3f} n>1 {nrho}" if san else "")
            + (f" | p'min {minp:.4f} eta_max {maxeta:.3f} n(p'<F) {ns_cols[0]}" if ns and ns_cols[0] != "" else ""))
        if good >= GROW_AFTER and ds < DS_MAX:
            ds, good = min(GROW * ds, DS_MAX), 0
    fh.close()

    s_end = rows[-1][0] if rows else 0.0
    q_end = rows[-1][1] if rows else q0
    qmax = max([r[1] for r in rows], default=q0)
    log(f"MODE = {mode}  s/B reached {s_end:.6f}  q_end {q_end:.3f} kPa  q_max "
        f"{qmax:.3f} kPa  steps {nstep}  failed attempts {nfail_tot}  push wall "
        f"{time.time()-t_push:.1f}s")

    summary = dict(mode=mode, s_over_B=s_end, q_end=q_end, q_max=qmax,
                   steps=nstep, failed_attempts=nfail_tot,
                   push_wall_s=time.time() - t_push, total_wall_s=time.time() - T0,
                   build=ops.ladrunoBuild().strip(), engine=ops.__file__,
                   args=vars(args), matdesc=deck["matdesc"], solver=deck["solver"],
                   k0_patch_err=deck["k0_patch_err"],
                   grav_resultant_err=deck["grav_resultant_err"],
                   applied=deck["applied"], ptol=ptol,
                   dead_end_first_step={int(k): int(v) for k, v in DEADEND.items()})
    if SAND in ("toyoura", "toyoura_campaign"):
        summary["toyoura"] = dict(TYR_INFO)
    summary["seating"] = dict(q0_kPa=q0, W_FOOT_kN_m=W_FOOT, note="no seating stage in this deck: q0 = (stage-0b footing weight - REF reaction)/B at the push start")
    if ns:
        Fz = read_field_ns(deck, full=True)
        Dz = derived_ns(Fz)
        nlat, nref = ns_latched(deck)
        summary["norsand"] = dict(
            params=NS["_path"], params_md5=NS["_md5"], params_name=NS.get("name"), status=NS.get("status"),
            K0=K0, K0_source=NS_K0_SOURCE[0], energy=NS["energy_eff"], material_flags_final=NS["flags"],
            p_min_eff=NS["pmin_eff"], pi0_mode=NS["pi0"]["mode"], pi0_check=NS.get("_check"),
            floor_init=NS.get("_floor_init"), floor_end=ns_floor_census(deck), kimura=dict(KIMURA_INFO),
            rows=deck["ns_rows"], latched_points=nlat, refused_trials=nref,
            p_floor_counts={str(F_): int(np.sum(Dz["p"] < F_)) for F_ in NS["p_floor_report_kPa"]},
            min_p_end=float(Dz["p"].min()), max_eta_end=float(np.nanmax(Dz["eta"])),
            min_psi_end=float(np.nanmin(Fz["psi"])), max_psi_end=float(np.nanmax(Fz["psi"])))
        save_field(os.path.join(out, "ckpt", "field_end.npz"), deck, Fz, Dz, nstep, s_end)
        log(f"norsand end: latched {nlat}, refused trials {nref}, p' floor counts "
            f"{summary['norsand']['p_floor_counts']}, min p' {summary['norsand']['min_p_end']:.4f} kPa, "
            f"kernel floor census {summary['norsand']['floor_end']}")
    if san and CENSUS is CENSUS_SAS:
        try:
            Fz = read_field(deck, want_stats=True, want_f=False)
            tot = Fz["stats"].sum(axis=0)
            summary["sas_totals"] = {CENSUS["names"][i]: float(tot[i])
                                     for i in range(min(len(tot), len(CENSUS["names"])))}
            summary["refusals_by_code"] = {SAS_REF_CODES[i]: int(tot[i])
                                           for i in SAS_REF_CODES if i < len(tot)}
            log(f"cumulative sasStats refusals by code: {summary['refusals_by_code']}")
        except Exception as exc:
            log(f"NOTE sas totals failed: {exc}")

    # ---- post-mortem at the wall (SANISAND): replay rows + one probe --------
    if san and mode in ("FLOOR", "WALL", "PROBETEST") and nstep > 0:
        Fn = read_field(deck, want_stats=True, want_f=True)
        Dn = derived(Fn)
        save_field(os.path.join(out, "ckpt", "field_last_converged.npz"), deck, Fn,
                   Dn, nstep, s_end)
        log("state at the last converged step:\n" + ring_summary(deck, Dn, Fn))
        # dump the census of the last converged step as a CSV (per point)
        xyc = np.array([(t, gx, gy) for t, gx, gy in gp_xy(deck)])
        np.savetxt(os.path.join(out, "census_last_converged.csv"),
                   np.column_stack([np.arange(len(Dn["p"])), xyc, Dn["p"], Dn["rho"],
                                    Fn["stats"]]),
                   delimiter=",", header="k,element,gp_x,gp_y,p,rho_alpha,"
                   + ",".join(CENSUS["names"]), comments="",
                   fmt="%.10g")
        # the last converged PAIR: committed state n-1 + the increment of step n
        if args.field_every_step and Fpp is not None:
            sel = select_worst(Dpp, Fpp, sub_prev)
            write_replay(os.path.join(out, "replay", "replay_wall_last_pair.csv"),
                         deck, Fpp, Dpp, Fprev, Dprev, eps_ppp, sel, nstep - 1,
                         sB_pp, sub_prev, cap_prev, dt_prev)
            log("wrote replay_wall_last_pair.csv (state at step n-1 + dStrain of "
                "the last converged step n)")
        if mode in ("FLOOR", "PROBETEST"):
            # probe: from the last converged state, ONE Newton iteration of the
            # first increment the wall step tried (FixedNumIter 1 commits it; the
            # run is over, so committing a non-equilibrium iterate is harmless).
            dsp = first_try_ds
            got = False
            for _ in range(12):
                ops.integrator("LoadControl", -dsp)
                ops.test("FixedNumIter", 1, 0)
                ops.algorithm("Newton")
                if ops.analyze(1) == 0:
                    got = True
                    break
                log(f"wall probe: iterate 1 at ds = {dsp:.4e} refused; halving")
                dsp *= 0.5
            if got:
                Fp = read_field(deck, want_stats=True, want_f=True)
                Dp = derived(Fp)
                subp = Fp["stats"][:, CENSUS["last_sub"]]
                capp = Fp["stats"][:, CENSUS["last_cap"]]
                sel = select_worst(Dn, Fn, subp)
                write_replay(os.path.join(out, "replay", "replay_wall_probe_iter1.csv"),
                             deck, Fn, Dn, Fp, Dp, Fpp["eps"] if Fpp is not None else None,
                             sel, nstep, s_end, subp, capp, -dsp)
                log(f"wall probe: 1 Newton iterate at ds = {dsp:.4e} m committed "
                    f"(first try was {first_try_ds:.4e}); iterate-1 substeps total "
                    f"{int(subp.sum())}, max/pt {int(subp.max())}, capHits {int(capp.sum())}")
                summary["probe_ds"] = dsp
            else:
                log("wall probe: every probe increment refused at iterate 1")
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import plateau as _plat
        summary["plateau"] = _plat.plateau_from_rows([r[0] for r in rows], [r[1] for r in rows],
                                                     q0, mode)
        log("plateau rule (TIMs WP1): " + json.dumps({k: summary["plateau"].get(k) for k in
                                                       ("verdict", "k_norm", "k_tail_kPa_per_sB", "k0_kPa_per_sB", "trend")}))
    except Exception as exc:
        log(f"NOTE plateau verdict failed: {exc}")
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
