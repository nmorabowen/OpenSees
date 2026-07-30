"""ADR-79 follow-on — the NUMERICAL COLLAPSE LOAD of the benchmark footing on
the cone PDMY actually has.

THE QUESTION. The bearing campaign found the PDMY footing hardening past
5.3x Vesic with no limit point, and `cone_probe.py` explained it: PDMY's
failure surface is a Lode-independent Drucker-Prager cone (alpha = 0.2436)
calibrated in TRIAXIAL COMPRESSION, while Vesic's N_gamma is a Mohr-Coulomb
formula. Matched in PLANE STRAIN (Chen & Han) that cone is phi = 53.7 deg,
whose N_gamma term is 10.8 MPa against the 207 kPa anchor the project quotes —
a factor of 52, which would put the benchmark at ~31 % of its own capacity and
dissolve the "9.3x over-strength" anomaly entirely.

That factor is a HYPOTHESIS, and it rests on three things this script removes:

  1. Chen & Han matching assumes ASSOCIATED flow. This PDMY set is
     non-dilatant (dil1 = dil2 = 0), and non-associated flow lowers a collapse
     load — sometimes a lot.
  2. It treats a SQUARE footing as plane strain and then repairs it with an
     empirical shape factor s_gamma = 0.6.
  3. It is a bearing-capacity FORMULA, evaluated far outside the friction-angle
     range any of its coefficients were ever calibrated in (N_q = 673).

So: run the actual 3D footing, on the actual mesh, with the actual cone, and
MEASURE the collapse load.

METHOD. PDMY cannot deliver a collapse load — 20 nested surfaces approach the
envelope asymptotically, and 11 h of pushing got to 31 % of it. So the material
is replaced by an ELASTIC-PERFECTLY-PLASTIC Drucker-Prager sharing the same
failure surface (`dp_cone_check.py` verifies the surface agrees with PDMY's
measured cone to 0.18 %, and that the probe instrument itself is exact to
0.00 %). By the limit theorems a collapse load is a property of the failure
surface and the flow rule, not of the hardening path that reaches them, so the
substitution is legitimate — and it is the ONLY thing changed:

  * same graded 2816-hex mesh (`bearing_mesh.npz`), same 2 x 2 m smooth footing
    driven on 25 nodes, same 10 kPa surcharge with the same tributary areas and
    the same closed-form reaction correction;
  * DRAINED, so the u-p element is dropped for a plain 3-DOF LadrunoBrick with
    the BUOYANT unit weight gamma' = 9.81 kN/m3 as a body force. The campaign
    validated that the u-p model carries buoyancy in the pore-pressure field
    and its effective stress is gamma'*z + q0, so this reproduces the same
    effective-stress state at a third of the DOF and none of the p-row
    conditioning;
  * `-formulation bbar` is the PRIMARY. Non-associated flow with rho_bar = 0
    is exactly ISOCHORIC plastic flow, which is the locking-prone case, and the
    campaign measured locking at 34.3 % of q on this very mesh. `std` is run as
    the paired locking probe, not as the answer.

LEGS. `nonassoc` (rho_bar = 0) is the physical match to PDMY's non-dilatant
set and is the headline number; `assoc` (rho_bar = rho) is the associated
upper bound Chen & Han's matching assumes. The gap between them is the size of
assumption (1). `*_std` are the locking pair, `*_sy2` the numerical-cohesion
sensitivity, `verify_phi20` the machinery check (below).

VERIFICATION LEG. `verify_phi20` reruns everything at a cone matched to
phi_txc = 20 deg, where N_q = 12 and N_gamma = 6 sit inside the range Vesic's
coefficients were fitted in and the plane-strain correction is mild. If the
numerical collapse load agrees with Vesic-at-the-plane-strain-equivalent there
and disagrees at the real cone, the disagreement is a statement about the
FORMULA at phi = 54 deg, not about this solver.

Run:  py -3.12 dp_collapse.py [leg|all]
"""
import os
import sys
import csv
import math
import time

import numpy as np

_WT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
_BIN = os.path.join(_WT, "dist", "bin")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)

# ADR79_NO_ENGINE lets `dp_analyze.py` import this module for its constants and
# its Vesic/cone helpers WITHOUT loading the engine. That matters because the
# only interpreter here with matplotlib is the apeGmsh venv, and that venv
# carries the installed Ladruno's site .pth — which would hijack `opensees` and
# trip the guard below for a run that never touches the solver.
if os.environ.get("ADR79_NO_ENGINE") == "1":
    ops = None
else:
    import opensees as ops  # noqa: E402

    # see bearing_backbone.py: an installed Ladruno hijacks `import opensees`
    # through a site .pth, and that build refuses the options this needs.
    if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
            os.path.normcase(_BIN):
        raise SystemExit(
            f"WRONG ENGINE: imported {ops.__file__}, expected {_BIN}")

HERE = os.path.dirname(os.path.abspath(__file__))
MESH = os.path.join(HERE, "bearing_mesh.npz")

# --- the model, identical to the campaign's except for the material ---------
G9 = 9.81
GAMMA_B = (2.0 - 1.0) * G9            # buoyant unit weight, kN/m3
K_EL, G_EL = 1.5e5, 5.5e4             # PDMY's reference moduli, kPa
B_FOOT, A_FOOT = 2.0, 4.0
Q_SURCH = 10.0                        # kPa
SFRAC = 0.15                          # push to s/B (collapse arrives earlier)

ALPHA_PDMY = 0.24362                  # cone_probe.py, txc + txe
SY_BASE = 0.2                         # kPa; apex regularizer, NOT a soil property


def alpha_from_phi_txc(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


def mc_from_cone(alpha):
    """Mohr-Coulomb angles matching a cone of slope alpha (Chen & Han)."""
    out = {}
    s = 3.0 * math.sqrt(3.0) * alpha / (2.0 + math.sqrt(3.0) * alpha)
    out["txc"] = math.degrees(math.asin(s)) if abs(s) <= 1 else float("nan")
    s = 3.0 * math.sqrt(3.0) * alpha / (2.0 - math.sqrt(3.0) * alpha)
    out["txe"] = math.degrees(math.asin(s)) if abs(s) <= 1 else float("nan")
    d = 1.0 - 12.0 * alpha ** 2
    out["ps"] = math.degrees(math.atan(math.sqrt(9.0 * alpha ** 2 / d))) \
        if d > 0 else float("nan")
    return out


def vesic(phi_deg, q0=Q_SURCH, b=B_FOOT, gamma=GAMMA_B):
    """Vesic q_ult for a SQUARE footing: q0*Nq*sq + 0.5*gamma*B*Ngamma*sgamma."""
    ph = math.radians(phi_deg)
    nq = math.exp(math.pi * math.tan(ph)) * math.tan(math.pi / 4 + ph / 2) ** 2
    ng = 2.0 * (nq + 1.0) * math.tan(ph)
    sq, sg = 1.0 + math.tan(ph), 0.6
    return dict(phi=phi_deg, Nq=nq, Ngamma=ng, q_term=q0 * nq * sq,
                g_term=0.5 * gamma * b * ng * sg,
                q_u=q0 * nq * sq + 0.5 * gamma * b * ng * sg)


# --- legs -------------------------------------------------------------------
def _leg(alpha, assoc, form="bbar", sy=SY_BASE, geom="linear",
         mesh="bearing_mesh.npz"):
    return dict(alpha=alpha, assoc=assoc, form=form, sy=sy, geom=geom,
                mesh=mesh)


LEGS = {
    # headline pair on the measured PDMY cone
    "nonassoc":   _leg(ALPHA_PDMY, False),
    "assoc":      _leg(ALPHA_PDMY, True),
    # locking pair (same physics, std formulation)
    "nonassoc_std": _leg(ALPHA_PDMY, False, form="std"),
    "assoc_std":    _leg(ALPHA_PDMY, True, form="std"),
    # numerical-cohesion sensitivity: sigma_y is a regularizer, so the collapse
    # load must not care about it. 2.0 kPa is 10x the base.
    "nonassoc_sy2": _leg(ALPHA_PDMY, False, sy=2.0),
    # machinery verification at a friction angle Vesic was actually fitted near
    "verify_phi20":     _leg(alpha_from_phi_txc(20.0), True),
    "verify_phi20_na":  _leg(alpha_from_phi_txc(20.0), False),
    # geometry cross-check: a collapse load is a small-strain statement, but the
    # campaign's backbones were run with richer kinematics, so measure the rung.
    "nonassoc_corot": _leg(ALPHA_PDMY, False, geom="corot"),
    # DOMAIN-SIZE control (build_mesh_big.py): 14.5 B of clearance and 10 B of
    # depth against the campaign mesh's 4.5 B / 4 B, identical near-field. A
    # phi_ps = 54 deg general-shear mechanism is far bigger than the campaign
    # box, so without this pair a measured collapse load could be the boundary's.
    "nonassoc_big": _leg(ALPHA_PDMY, False, mesh="bearing_mesh_big.npz"),
    "assoc_big":    _leg(ALPHA_PDMY, True, mesh="bearing_mesh_big.npz"),
}
GROUPS = {"all": ["nonassoc", "assoc"],
          "lock": ["nonassoc_std", "assoc_std"],
          "verify": ["verify_phi20", "verify_phi20_na"],
          "big": ["nonassoc_big", "assoc_big"],
          "extra": ["nonassoc_sy2", "nonassoc_corot"]}

# --- stepping ---------------------------------------------------------------
DS_BASE, DS_MIN, DS_MAX = 2.0e-4, 5.0e-6, 5.0e-4     # m of settlement per step
GROW_AFTER = 6
# KrylovNewton on this problem either converges in a handful of iterations or
# diverges; a generous iteration cap only makes every RETRY expensive (measured:
# a 60-iteration failed attempt cost 25 s against 1.5 s for a good step).
MAXITER = 25
# Convergence is judged on FORCE, not displacement: the campaign measured
# NormDispIncr reporting success at a state with +0.26 kPa of tension under a
# 100 kPa compressive load once the tangent degenerated (LEDGER_quirks). The
# tolerance is referred to the model's total dead load so it means the same
# thing at every stage of the push.
R_REF = GAMMA_B * (20.0 * 20.0 * 8.0) + Q_SURCH * (20.0 * 20.0)   # kN
TOL = 1.0e-5 * R_REF
MAXSTEP = int(os.environ.get("ADR79_MAXSTEP", 100000))


# An algorithm ladder per increment, before the increment is shrunk. A
# perfectly plastic collapse is where Newton is worst — the tangent goes
# singular in the mechanism mode exactly when the answer arrives — and step
# halving alone is not enough: the `nonassoc_std` leg walled out at s/B =
# 0.0022, 46 % of its initial tangent, with nothing measured. Steps that needed
# the RELAXED tolerance are flagged in the CSV so they cannot pass unnoticed.
# The first four legs (nonassoc, assoc, nonassoc_sy2, nonassoc_big) ran before
# this existed; the ladder changes only how FAR a leg gets, not the path, since
# every accepted step meets the same tolerance unless flagged.
def attempts(tol):
    return [("KrylovNewton", tol, MAXITER, 0),
            ("NewtonLineSearch", tol, 2 * MAXITER, 0),
            ("KrylovNewton", 10.0 * tol, 3 * MAXITER, 1)]


def set_attempt(a):
    algo, tol, it, _ = a
    ops.test("NormUnbalance", tol, it, 0)
    if algo == "NewtonLineSearch":
        ops.algorithm("NewtonLineSearch")
    else:
        ops.algorithm("KrylovNewton")


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def load_mesh(fname):
    p = os.path.join(HERE, fname)
    if not os.path.exists(p):
        raise SystemExit(f"{p} missing — run build_mesh.py / build_mesh_big.py")
    z = np.load(p)
    nodes = z["nodes"]
    # domain extents come from the mesh, not from constants, so the big-domain
    # legs share every check below without a second set of magic numbers.
    ext = dict(x=float(np.abs(nodes[:, 0]).max()),
               y=float(np.abs(nodes[:, 1]).max()),
               zbot=float(nodes[:, 2].min()))
    ext["volume"] = 2 * ext["x"] * 2 * ext["y"] * (-ext["zbot"])
    return (nodes, z["hexes"], z["tributary"],
            {k[4:]: z[k] for k in z.files if k.startswith("set_")}, ext)


def pick_system():
    try:
        ops.system("Pardiso")          # tangent is UNSYMMETRIC when rho_bar != rho
        return "Pardiso"
    except Exception:
        ops.system("UmfPack")
        return "UmfPack"


def build(cfg, nodes, hexes, trib, sets):
    rho_dp = math.sqrt(2.0) * cfg["alpha"]
    rho_bar = rho_dp if cfg["assoc"] else 0.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, zz) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(zz))
    # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta -> perfectly plastic
    ops.nDMaterial("DruckerPrager", 1, K_EL, G_EL, cfg["sy"], rho_dp, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    extra = [] if cfg["form"] == "std" else ["-formulation", cfg["form"]]
    for e, conn in enumerate(hexes, start=1):
        # -b is a body force per unit VOLUME (LadrunoBrick.cpp: bodyForce -=
        # dvol*b*N, no density factor), so gamma' goes in directly. It is also
        # applied unconditionally rather than through a load pattern, which is
        # why gravity below is a single full-amplitude static solve.
        ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                    "-geom", cfg["geom"], "-b", 0.0, 0.0, -GAMMA_B, *extra)

    fix = {}
    for n in sets["xface"]:
        fix.setdefault(int(n), [0, 0, 0])[0] = 1
    for n in sets["yface"]:
        fix.setdefault(int(n), [0, 0, 0])[1] = 1
    for n in sets["bottom"]:
        fix[int(n)] = [1, 1, 1]
    for n, f in fix.items():
        ops.fix(n + 1, *f)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        a = float(trib[n])
        if a > 0.0:
            ops.load(int(n) + 1, 0.0, 0.0, -Q_SURCH * a)
    return rho_dp, rho_bar


def gravity(tag, nodes, sets, trib, ext):
    """One static solve. The K0 state is elastic by a wide margin — the elastic
    K0 = nu/(1-nu) = 0.507 mobilises sin(phi) = 0.327 against the cone's txc
    equivalent 0.523 — so this is a linear solve and it is checked as one."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    sysname = pick_system()
    ops.test("NormUnbalance", TOL, MAXITER, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"{tag} gravity failed"
    # global vertical equilibrium — the campaign's first validation check,
    # repeated here because the body-force convention is the one thing this
    # model does differently from the u-p legs.
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want = GAMMA_B * ext["volume"] + Q_SURCH * float(trib.sum())
    log(tag, f"gravity equilibrium: sum Rz = {rz:.2f} kN vs "
        f"gamma'*V + q0*A = {want:.2f} kN "
        f"({100 * (rz / want - 1):+.4f} %)")
    if abs(rz / want - 1) > 1e-6:
        raise SystemExit(f"{tag}: gravity does not balance — body force or "
                         f"tributary convention is wrong")
    return sysname


def run_leg(name):
    cfg = LEGS[name]
    nodes, hexes, trib, sets, ext = load_mesh(cfg["mesh"])
    t_start = time.time()
    rho_dp, rho_bar = build(cfg, nodes, hexes, trib, sets)
    log(name, f"mesh {cfg['mesh']}: {len(nodes)} nodes, {len(hexes)} hex, "
        f"clearance {(ext['x'] - B_FOOT / 2) / B_FOOT:.1f} B, depth "
        f"{-ext['zbot'] / B_FOOT:.1f} B")
    mc = mc_from_cone(cfg["alpha"])
    log(name, f"alpha={cfg['alpha']:.5f} -> rho={rho_dp:.5f}, "
        f"rho_bar={rho_bar:.5f} ({'associated' if cfg['assoc'] else 'NON-associated'}), "
        f"sigma_y={cfg['sy']} kPa, -formulation {cfg['form']}, "
        f"-geom {cfg['geom']}")
    log(name, f"Chen-Han equivalents: phi_txc={mc['txc']:.2f} deg, "
        f"phi_ps={mc['ps']:.2f} deg; Vesic(phi_ps) q_u="
        f"{vesic(mc['ps'])['q_u']:.0f} kPa, Vesic(phi_txc) q_u="
        f"{vesic(mc['txc'])['q_u']:.0f} kPa")
    gravity(name, nodes, sets, trib, ext)
    log(name, f"gravity done [{time.time() - t_start:.0f}s]")

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = Q_SURCH * float(trib[sets["footing"]].sum())
    uz0 = ops.nodeDisp(foot[0], 3)
    ops.loadConst("-time", 0.0)

    smax = SFRAC * B_FOOT
    # load factor lambda in [0,1] IS s/smax: a 2-point Path on the sp value.
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0,
                   "-values", uz0, uz0 - smax, "-useLast")
    ops.pattern("Plain", 2, 2)
    for tt in foot:
        ops.sp(tt, 3, 1.0)

    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    pick_system()
    ops.test("NormUnbalance", TOL, MAXITER, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", DS_BASE / smax)
    ops.analysis("Static")

    path = os.path.join(HERE, f"dpcollapse_{name}.csv")
    if os.path.exists(path) and time.time() - os.path.getmtime(path) < 180 \
            and os.environ.get("ADR79_FORCE") != "1":
        raise SystemExit(f"REFUSING to write {os.path.basename(path)}: another "
                         f"process is probably still appending to it.")
    f = open(path, "w", newline="")
    w = csv.writer(f)
    w.writerow(["s_m", "s_meas_m", "R_kN", "q_kPa", "ds_mm", "relaxed",
                "wall_s"])
    rows, ds, good, nfail, nrelax = [], DS_BASE, 0, 0, 0
    lad = attempts(TOL)
    t0 = time.time()
    nstep = 0
    truncated_at = None
    while nstep < MAXSTEP:
        s_now = ops.getTime() * smax
        if s_now >= smax - 1e-12:
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", ds / smax)
        ok, relaxed = False, 0
        for a in lad:
            set_attempt(a)
            if ops.analyze(1) == 0:
                ok, relaxed = True, a[3]
                break
            nfail += 1
        if not ok:
            good = 0
            ds *= 0.5
            if ds < DS_MIN:
                truncated_at = s_now
                log(name, f"NO CONVERGENCE at s/B={s_now / B_FOOT:.4f} even at "
                    f"ds={ds * 1000:.4f} mm — stopping")
                break
            continue
        nstep += 1
        nrelax += relaxed
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:
            ds, good = min(2.0 * ds, DS_MAX), 0
        ops.reactions()
        R = -sum(ops.nodeReaction(tt, 3) for tt in foot) + r_corr
        s = ops.getTime() * smax
        row = (s, uz0 - ops.nodeDisp(foot[0], 3), R, R / A_FOOT, ds * 1000,
               relaxed, time.time() - t0)
        rows.append(row)
        w.writerow([f"{v:.9g}" for v in row])
        f.flush()
        if nstep % 50 == 0:
            log(name, f"s/B={s / B_FOOT:.4f} q={row[3]:.1f} kPa "
                f"ds={ds * 1000:.3f} mm [{nstep} steps, {nfail} failed "
                f"attempts, {nrelax} relaxed, {row[6]:.0f}s]")
    f.close()
    dump_field(name, cfg, nodes, hexes, len(rows), ext)
    verdict(name, rows, truncated_at, cfg)
    log(name, f"leg done: {nstep} steps, {nfail} failed attempts, "
        f"{nrelax} relaxed, "
        f"{time.time() - t_start:.0f}s -> {os.path.basename(path)}")
    return rows


def dump_field(name, cfg, nodes, hexes, nrow, ext):
    """Where the soil is actually AT its strength when the run ends.

    This is not decoration. Vesic's N_gamma describes a general-shear mechanism
    that at phi_ps = 54 deg is enormous — the log-spiral radius grows by
    exp(theta*tan phi), a factor of 8.5 over a quarter turn — while this domain
    offers 4.5 B of clearance and 4 B of depth. If the fully-mobilised zone
    reaches the roller boundary then the measured collapse load is the box's,
    not the soil's, and that has to be visible rather than argued about.

    Mobilisation is read straight off the stress: m = sqrt(J2) / (3*alpha*p +
    sigma_y/sqrt(3)), which is 1.0 exactly on the yield surface.
    """
    if nrow == 0:
        return
    alpha, sy = cfg["alpha"], cfg["sy"]
    cen = nodes[hexes].mean(axis=1)
    mob = np.zeros(len(hexes))
    pmean = np.zeros(len(hexes))
    for e in range(len(hexes)):
        s = np.array(ops.eleResponse(e + 1, "stress")).reshape(-1, 6).mean(axis=0)
        xx, yy, zz, xy, yz, zx = s
        I1 = xx + yy + zz
        p = -I1 / 3.0
        dev = np.array([[xx - I1 / 3, xy, zx], [xy, yy - I1 / 3, yz],
                        [zx, yz, zz - I1 / 3]])
        sqJ2 = math.sqrt(max(0.5 * np.tensordot(dev, dev), 0.0))
        cap = 3.0 * alpha * p + sy / math.sqrt(3.0)
        mob[e] = sqJ2 / cap if cap > 1e-9 else np.nan
        pmean[e] = p
    np.savez_compressed(os.path.join(HERE, f"dpcollapse_{name}_field.npz"),
                        centroid=cen, mob=mob, p=pmean)
    yielded = mob > 0.99
    log(name, f"mobilisation: {yielded.sum()} / {len(mob)} elements at "
        f"m > 0.99 ({100 * yielded.mean():.1f} %)")
    if yielded.any():
        c = cen[yielded]
        rx, ry, rz = (np.abs(c[:, 0]).max(), np.abs(c[:, 1]).max(),
                      c[:, 2].min())
        # "Reaches the boundary" must be tested by ELEMENT COLUMN membership,
        # not by a fraction of the domain extent. This mesh is graded: the
        # outermost hex is 3.1 m wide, so its centroid sits at |x| = 8.45 in a
        # 10 m half-domain — a "within 90 % of the extent" test calls that
        # contained when the element is in fact touching the roller. (Measured:
        # nonassoc_sy2 reported "contained" while 16 elements of the
        # boundary-adjacent column were at m > 0.99.)
        cols = np.unique(np.round(np.abs(cen[:, 0]), 4))
        rows = np.unique(np.round(cen[:, 2], 4))
        touch_x = np.isclose(np.round(rx, 4), cols.max())
        touch_z = np.isclose(np.round(rz, 4), rows.min())
        log(name, f"  yielded zone reaches |x|<={rx:.2f} m, |y|<={ry:.2f} m, "
            f"z>={rz:.2f} m (domain {ext['x']:.0f}, {ext['y']:.0f}, "
            f"{ext['zbot']:.0f}; outermost element centroid {cols.max():.2f}, "
            f"deepest {rows.min():.2f}) — "
            + ("TOUCHES THE SIDE BOUNDARY" if touch_x else "clear of the sides")
            + ("; TOUCHES THE BASE" if touch_z else "; clear of the base"))


def verdict(name, rows, truncated_at, cfg):
    """A collapse load must be ASSERTED from the tangent, not from where the
    run stopped — the campaign's whole point. Report q_max, the tangent over
    the last decade of settlement, and whether q is still climbing."""
    if len(rows) < 20:
        log(name, "too few steps for a verdict")
        return
    a = np.array(rows)
    s, q = a[:, 0], a[:, 3]
    qmax, imax = q.max(), int(q.argmax())
    tail = q[int(0.9 * len(q)):]
    stail = s[int(0.9 * len(q)):]
    dqds = np.polyfit(stail, tail, 1)[0]
    dqds0 = np.polyfit(s[:len(s) // 20 + 2], q[:len(s) // 20 + 2], 1)[0]
    mc = mc_from_cone(cfg["alpha"])
    log(name, f"q_max = {qmax:.1f} kPa at s/B = {s[imax] / B_FOOT:.4f}; "
        f"q_end = {q[-1]:.1f} kPa at s/B = {s[-1] / B_FOOT:.4f}")
    log(name, f"tail tangent dq/ds = {dqds:.1f} kPa/m "
        f"({100 * dqds / dqds0:.2f} % of initial {dqds0:.0f}) => "
        f"{'PLATEAU (collapse reached)' if abs(dqds) < 0.02 * dqds0 else 'still hardening'}"
        + ("" if truncated_at is None else "  [run truncated on convergence]"))
    log(name, f"vs Vesic(phi_ps={mc['ps']:.2f}) = {vesic(mc['ps'])['q_u']:.0f} kPa"
        f"  -> numeric/Vesic_ps = {qmax / vesic(mc['ps'])['q_u']:.3f}")
    log(name, f"vs Vesic(phi_txc={mc['txc']:.2f}) = {vesic(mc['txc'])['q_u']:.0f} kPa"
        f"  -> numeric/Vesic_txc = {qmax / vesic(mc['txc'])['q_u']:.3f}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    names = GROUPS.get(which, [which])
    for n in names:
        if n not in LEGS:
            raise SystemExit(f"unknown leg {n!r}; pick from {list(LEGS)} "
                             f"or a group {list(GROUPS)}")
    v33 = vesic(33.0)
    log(f"anchors: Vesic(phi=33 nominal) q_u={v33['q_u']:.1f} kPa "
        f"(q-term {v33['q_term']:.1f} + gamma-term {v33['g_term']:.1f}); "
        f"PDMY backbone at s/B=0.15 was 3384 kPa")
    for n in names:
        log("=== leg:", n)
        run_leg(n)


if __name__ == "__main__":
    main()
