"""WP-138 -- strip-footing A/B: SANISAND ModifiedEuler (IntScheme 1) vs SAS-ME
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
import csv
import json
import math
import os
import sys
import time

# ---------------------------------------------------------------- engine ----
_BIN = os.environ.get("FOOTING_BIN")
if not _BIN or not os.path.isdir(_BIN):
    raise SystemExit("FOOTING_BIN must name the dist/bin directory to test")
os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")

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
    assert len(x) == 91 and len(y) == 28, (len(x), len(y))
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


def build(mat, scheme, extra, maxsub, dp_phi, dp_g, deterministic):
    x, y = mesh_coords()
    nx, ny = len(x), len(y)

    def nt(i, j):
        return i * ny + j + 1

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i in range(nx):
        for j in range(ny):
            ops.node(nt(i, j), float(x[i]), float(y[j]))
    REF = 900001
    ops.node(REF, 0.0, 0.0, "-ndf", 3)

    if mat == "sanisand":
        args = ["LadrunoSANISAND", MAT, *SAN, int(scheme), 0, 1, 1.0e-7, 1.0e-7,
                "-flipAlphaIn", "init", "-Pmin", 0.0101,
                "-maxSubsteps", int(maxsub), "-Presidual", 0.0,
                "-honorTolR", 0, *extra]
        ops.nDMaterial(*args)
        matdesc = " ".join(str(a) for a in args)
    elif mat == "dp":
        # UW DruckerPrager, psi = 0 (rhoBar = 0), plane-strain matched cone at
        # psi = 0: sqrt(J2) = sin(phi) p  ->  rho = sqrt(2) sin(phi) / 3.
        phi = math.radians(dp_phi)
        rho = math.sqrt(2.0) * math.sin(phi) / 3.0
        kk = 2.0 * dp_g * (1.0 + NU_STAR) / (3.0 * (1.0 - 2.0 * NU_STAR))
        args = ["DruckerPrager", MAT, kk, dp_g, 0.0, rho, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ops.nDMaterial(*args)
        matdesc = " ".join(str(a) for a in args)
    else:
        raise SystemExit(f"unknown mat {mat}")

    elems = element_table(x, y)
    info = {}
    for tag, i, j in elems:
        ops.element("LadrunoQuad", tag, nt(i, j), nt(i + 1, j), nt(i + 1, j + 1),
                    nt(i, j + 1), MAT, "-formulation", "bbar",
                    "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
        xc = 0.25 * (x[i] + x[i + 1]) * 2
        yc = 0.25 * (y[j] + y[j + 1]) * 2
        gps = []
        for xi, et in GP_XI:
            g = 1.0 / math.sqrt(3.0)
            gps.append((0.5 * (x[i] + x[i + 1]) + 0.5 * xi * g * (x[i + 1] - x[i]),
                        0.5 * (y[j] + y[j + 1]) + 0.5 * et * g * (y[j + 1] - y[j])))
        info[tag] = dict(xc=float(xc), yc=float(yc), gps=gps,
                         area=float((x[i + 1] - x[i]) * (y[j + 1] - y[j])))

    for i in range(nx):
        for j in range(ny):
            if j == 0:
                ops.fix(nt(i, j), 1, 1)
            elif i == 0 or i == nx - 1:
                ops.fix(nt(i, j), 1, 0)
    ops.fix(REF, 1, 0, 1)

    half = 0.5 * B_FOOT
    foot = [nt(i, ny - 1) for i in range(nx) if abs(x[i]) <= half + 1e-9]
    assert len(foot) == 9, len(foot)
    ops.element("LadrunoKinematicCoupling", 900001, REF, len(foot), *foot,
                "-dof", 1, 2)

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
        ops.system("Pardiso")
        solver = "Pardiso"
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.updateMaterialStage("-material", MAT, "-stage", 0)

    # ---- stage 0a: self weight -------------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tags = [t for t, _, _ in elems]
    ops.eleLoad("-ele", *tags, "-type", "-selfWeight", 0.0, 1.0)
    ops.integrator("LoadControl", 1.0 / N_GRAV)
    for k in range(N_GRAV):
        rc = ops.analyze(1)
        if rc != 0:
            raise SystemExit(f"gravity step {k+1} failed rc={rc}")
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
        rc = ops.analyze(1)
        if rc != 0:
            raise SystemExit(f"surface-load step {k+1} failed rc={rc}")
    ops.reactions()
    ry2 = sum(ops.nodeReaction(n, 2) for n in base)
    applied = want + surf + W_FOOT
    log(f"surface loads: surcharge {surf:.4f} kN/m (+{W_FOOT} footing); base "
        f"reaction {ry2:.4f} vs {applied:.4f} (rel {abs(ry2/applied-1):.2e})")
    ops.loadConst("-time", 0.0)
    deck = dict(x=x, y=y, nx=nx, ny=ny, tags=tags, info=info, foot=foot,
                REF=REF, applied=applied, solver=solver, matdesc=matdesc,
                grav_resultant_err=res_err, k0_patch_err=perr, mat=mat)
    return deck


# ------------------------------------------------------------ field reads ---
STAT_NAMES = ["updates", "meCalls", "substeps", "accepted", "rejectedErr",
              "forcedAtDTmin", "forcedClampMc", "rejectedLowP", "abandonedLowP",
              "capHits", "entryPminClamps", "pnResets", "maxSubstepsOneUpdate",
              "lastSubsteps", "lastForcedAtDTmin", "lastAbandonedLowP",
              "lastCapHit"]


def read_field(deck, want_stats=True, want_f=False):
    """Every Gauss point's committed state. Arrays of shape (nGP, ...)."""
    tags = deck["tags"]
    n = 4 * len(tags)
    sig = np.zeros((n, 3)); eps = np.zeros((n, 3)); st = np.zeros((n, 26))
    psi = np.full(n, np.nan); fy = np.full(n, np.nan)
    stats = np.zeros((n, len(STAT_NAMES))) if want_stats else None
    k = 0
    for tag in tags:
        for gp in (1, 2, 3, 4):
            sig[k] = ops.eleResponse(tag, "material", gp, "stress")
            eps[k] = ops.eleResponse(tag, "material", gp, "strain")
            if deck["mat"] == "sanisand":
                st[k] = ops.eleResponse(tag, "material", gp, "state")
                psi[k] = ops.eleResponse(tag, "material", gp, "psi")[0]
                if want_f:
                    fy[k] = ops.eleResponse(tag, "material", gp, "yieldDistance")[0]
                if want_stats:
                    v = ops.eleResponse(tag, "material", gp, "substepStats")
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
    np.savez_compressed(path, step=step, s_over_B=sB, tag=xy[:, 0], gx=xy[:, 1],
                        gy=xy[:, 2], sig=F["sig"], eps=F["eps"], st=F["st"],
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


# ------------------------------------------------------------------- run ----
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", choices=("sanisand", "dp"), default="sanisand")
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--extra", default="", help="extra material tokens, space separated")
    ap.add_argument("--maxsub", type=int, default=2000)
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
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    out = args.out
    for d in ("ckpt", "replay"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    log(f"engine {ops.__file__}")
    log(f"build {ops.ladrunoBuild().strip().splitlines()[0]}")
    log(f"env MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS')} "
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')} "
        f"MKL_CBWR={os.environ.get('MKL_CBWR')}")
    log(f"args {vars(args)}")
    extra = []
    for tok in args.extra.split():
        try:
            extra.append(float(tok) if any(c in tok for c in ".e") and not tok.startswith("-") else int(tok))
        except ValueError:
            extra.append(tok)
    san = args.mat == "sanisand"
    deck = build(args.mat, args.scheme, extra, args.maxsub, args.dp_phi, args.dp_g,
                 bool(args.deterministic))
    log(f"material: {deck['matdesc']}")
    log(f"mesh: {len(deck['tags'])} elements, {4*len(deck['tags'])} Gauss points, "
        f"solver {deck['solver']}, applied vertical load {deck['applied']:.3f} kN/m")

    F = read_field(deck, want_stats=san, want_f=san)
    D = derived(F) if san else None
    if san:
        log(f"at the flip (stage 0 end): max eta/Mc = {np.nanmax(D['eta'])/MC:.4f}, "
            f"min p' = {D['p'].min():.4f} kPa, max |f_read - f_recomputed| = "
            f"{np.nanmax(np.abs(F['f'] - D['frec'])):.3e}")
        save_field(os.path.join(out, "ckpt", "field_flip.npz"), deck, F, D, 0, 0.0)
    ops.updateMaterialStage("-material", MAT, "-stage", 1)

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
                 "min_p", "max_eta", "max_rho_alpha", "n_rho_gt_1"])
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
    while True:
        s_now = uy0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
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
            sub_pt = dst[:, 2]
            cen = [int(sub_pt.sum()), int(sub_pt.max()), int(dst[:, 9].sum()),
                   int(dst[:, 4].sum()), int(dst[:, 5].sum()), int(dst[:, 7].sum()),
                   int(dst[:, 6].sum())]
            minp, maxeta = float(Dn["p"].min()), float(np.nanmax(Dn["eta"]))
            maxrho = float(np.nanmax(Dn["rho"]))
            nrho = int(np.sum(Dn["rho"] > 1.0))
            # replay rows: committed state of the PREVIOUS step + this step's dS
            for (Fc, Dc, stc, sbc, epsc) in pending:
                sel = select_worst(Dc, Fc, sub_pt)
                write_replay(os.path.join(out, "replay", f"replay_step{stc:05d}.csv"),
                             deck, Fc, Dc, Fn, Dn, epsc, sel, stc, sbc, sub_pt,
                             dst[:, 9], -ds)
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
            Fprev, Dprev, sub_prev, cap_prev, dt_prev = Fn, Dn, sub_pt, dst[:, 9], -ds
            sB_prev = sB
        wr.writerow([nstep, f"{s:.9e}", f"{sB:.9e}", f"{q:.6f}", f"{ds:.4e}", rung,
                     iters, fails_before, f"{wall_step:.3f}",
                     f"{time.time()-t_push:.1f}", *cen,
                     f"{minp:.5g}", f"{maxeta:.5g}", f"{maxrho:.5g}", nrho])
        fh.flush()
        rows.append((sB, q))
        log(f"step {nstep:5d} s/B {sB:.6f} q {q:9.3f} kPa ds {ds:.3e} rung "
            f"{'NLK'[rung]} it {iters:3d} fails {fails_before} wall {wall_step:8.2f}s"
            + (f" | sub {cen[0]} maxpt {cen[1]} cap {cen[2]} rejErr {cen[3]} "
               f"dtmin {cen[4]} | p'min {minp:.3f} eta_max {maxeta:.3f} "
               f"rho_max {maxrho:.3f} n>1 {nrho}" if san else ""))
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
                   applied=deck["applied"], ptol=ptol)

    # ---- post-mortem at the wall (SANISAND): replay rows + one probe --------
    if san and mode in ("FLOOR", "WALL", "PROBETEST") and nstep > 0:
        Fn = read_field(deck, want_stats=True, want_f=True)
        Dn = derived(Fn)
        save_field(os.path.join(out, "ckpt", "field_last_converged.npz"), deck, Fn,
                   Dn, nstep, s_end)
        log("state at the last converged step:\n" + ring_summary(deck, Dn, Fn))
        # dump the census of the last converged step as a CSV (per point)
        np.savetxt(os.path.join(out, "census_last_converged.csv"),
                   np.column_stack([np.arange(len(Dn["p"])), Fn["stats"]]),
                   delimiter=",", header="k," + ",".join(STAT_NAMES), comments="",
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
                subp = Fp["stats"][:, 13]
                capp = Fp["stats"][:, 16]
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
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
