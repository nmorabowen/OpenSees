"""WP-238 demo: a plane-strain biaxial compression test on LadrunoSANISAND with a
correlated void-ratio field (Ladruno_scripts/sanisand_void_ratio_field.py).

    set BIAX_BIN=<worktree>\\dist\\bin
    py -3.12 -S -u biaxial.py --out out/seed1 --seed 1
    py -3.12 -S -u biaxial.py --out out/homog            (no field)

THE TEST
  * specimen W x H = 0.10 m x 0.20 m, nx x ny LadrunoQuad -formulation bbar, plane strain;
  * smooth platens: bottom u_y = 0 with u_x fixed at the middle node only; the top nodes
    share u_y (equalDOF to the top-left node) and are free in x;
  * stage 0 (elastic): sigma_3 on the two sides and on the top, as dead nodal loads;
  * stage 1: the top is pushed down (sp on the master, LoadControl, pseudo-time = u_y)
    with the side pressure held;
  * material: the TIMs campaign set, SAS-ME (IntScheme 129) with the R1 options.

THE FIELD (only with --seed): vf.field at the element centroids, applied with vf.apply
before the first step. --std and --theta are ILLUSTRATIVE inputs, not calibrated numbers.

UNITS kPa, m, kN. Output: <out>/run.npz and <out>/log.txt.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

_BIN = os.environ.get("BIAX_BIN")
if not _BIN or not os.path.isdir(_BIN):
    raise SystemExit("BIAX_BIN must name the dist/bin directory to test")
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
if os.environ.get("BIAX_SITE"):            # `python -S` hides numpy: name its site-packages
    sys.path.append(os.environ["BIAX_SITE"])
elif os.name == "nt":
    sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("BIAX_SCRIPTS") or
                os.path.normpath(os.path.join(HERE, "..", "..", "..", "Ladruno_scripts")))

import numpy as np                              # noqa: E402
import opensees as ops                          # noqa: E402
import sanisand_void_ratio_field as vf          # noqa: E402

SAN = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0,
       0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
OPTS = ["SCH", "TAN", 1, 1.0e-7, "TOLR", "-flipAlphaIn", "init", "-Pmin", 0.0101,
        "-maxSubsteps", "MAXSUB", "-Presidual", 0.0, "-honorTolR", 0,
        "-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5]
MAT = 1
# SAS-ME hands back a tangent on which plain Newton converges only linearly; small steps under
# KrylovNewton are several times cheaper per unit strain than large ones (measured, README).
LADDER = (("KrylovNewton", 1.0, 60), ("NewtonLineSearch", 1.0, 60), ("Newton", 1.0, 60))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=None, help="omit for the homogeneous run")
    ap.add_argument("--nx", type=int, default=20)
    ap.add_argument("--ny", type=int, default=40)
    ap.add_argument("--W", type=float, default=0.10)
    ap.add_argument("--H", type=float, default=0.20)
    ap.add_argument("--sig3", type=float, default=100.0)
    ap.add_argument("--std", type=float, default=0.02)
    ap.add_argument("--theta", type=float, default=0.02)
    ap.add_argument("--emin", type=float, default=0.60)
    ap.add_argument("--emax", type=float, default=0.80)
    ap.add_argument("--white", action="store_true",
                    help="negative control: one independent normal per element instead of the field")
    ap.add_argument("--eps-end", type=float, default=0.14)
    ap.add_argument("--deps", type=float, default=2.0e-5)
    ap.add_argument("--snaps", type=int, default=14)
    ap.add_argument("--scheme", type=int, default=129)
    ap.add_argument("--tolr", type=float, default=1.0e-4)
    ap.add_argument("--maxsub", type=int, default=2000)
    ap.add_argument("--wall", type=float, default=1.0e9, help="stop and save after this many seconds")
    ap.add_argument("--log-every", type=int, default=200)
    ap.add_argument("--tantype", type=int, default=1, help="0 elastic, 1 continuum")
    ap.add_argument("--solver", default="static", choices=["static", "dr"],
                    help="dr: LadrunoDynamicRelaxation, continuous platen ramp (no tangent)")
    ap.add_argument("--rate", type=float, default=2.0e-8, help="dr: platen travel per DR step (m)")
    ap.add_argument("--mass-safety", type=float, default=0.5, help="dr: omega_max*dt <= 2f")
    ap.add_argument("--chunk", type=int, default=500, help="dr: steps per curve sample")
    ap.add_argument("--hold-every", type=float, default=0.01,
                    help="dr: axial-strain interval of the hold checks (0 = only the final one)")
    ap.add_argument("--hold-max", type=int, default=40000, help="dr: step cap of one hold")
    ap.add_argument("--hold-tol", type=float, default=1.0e-3,
                    help="dr: hold gate, residual inf-norm / (sig3 * hy)")
    ap.add_argument("--dr-damping", default="kinetic", choices=["kinetic", "viscous"])
    ap.add_argument("--zeta", type=float, default=1.0, help="dr: viscous damping ratio")
    ap.add_argument("--reversal-tol", type=float, default=None,
                    help="material -reversalTol: no alpha_in reset below this strain increment")
    ap.add_argument("--reversal-rel", type=float, default=None, help="material -reversalRel")
    ap.add_argument("--tension-cutoff", type=float, nargs=2, default=None,
                    metavar=("P_SEP", "P_CONTACT"),
                    help="material -sasTensionCutoff (WP-152): a point that refuses at low p separates")
    ap.add_argument("--tol-relax", type=float, default=1.0,
                    help="static: once the step is cut 16x (64x) accept a step at 3x (this "
                         "factor) the unbalance tolerance; 1 = never. Counted and reported.")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "w", encoding="utf-8")
    t0 = time.time()

    def log(m):
        line = f"[{time.time()-t0:8.1f}s] {m}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    nx, ny, W, H = a.nx, a.ny, a.W, a.H
    hx, hy = W / nx, H / ny

    def nt(i, j):
        return j * (nx + 1) + i + 1

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    xy = np.zeros(((nx + 1) * (ny + 1), 2))
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(nt(i, j), i * hx, j * hy)
            xy[nt(i, j) - 1] = (i * hx, j * hy)
    for i in range(nx + 1):
        ops.fix(nt(i, 0), 1 if i == nx // 2 else 0, 1)
    MASTER = nt(0, ny)
    for i in range(1, nx + 1):
        ops.equalDOF(MASTER, nt(i, ny), 2)

    if os.environ.get("BIAX_NOR1"):
        del OPTS[OPTS.index("-sasHFloor"):]
    if a.tension_cutoff is not None:
        OPTS.extend(["-sasTensionCutoff", *a.tension_cutoff])
    if a.reversal_tol is not None:
        OPTS.extend(["-reversalTol", a.reversal_tol])
    if a.reversal_rel is not None:
        OPTS.extend(["-reversalRel", a.reversal_rel])
    sub = {"TAN": a.tantype, "SCH": a.scheme, "TOLR": a.tolr, "MAXSUB": a.maxsub}
    ops.nDMaterial("LadrunoSANISAND", MAT, *SAN, *[sub.get(o, o) if isinstance(o, str) else o for o in OPTS])
    conn, cen, tags = [], [], []
    for j in range(ny):
        for i in range(nx):
            tag = j * nx + i + 1
            c = (nt(i, j), nt(i + 1, j), nt(i + 1, j + 1), nt(i, j + 1))
            ops.element("LadrunoQuad", tag, *c, MAT, "-formulation", "bbar",
                        "-type", "PlaneStrain", "-thick", 1.0)
            conn.append(c)
            cen.append(((i + 0.5) * hx, (j + 0.5) * hy))
            tags.append(tag)
    conn, cen = np.array(conn), np.array(cen)

    # ---- the field ---------------------------------------------------------
    e = np.full(len(tags), SAN[2])
    if a.seed is not None:
        psi0 = ops.eleResponse(tags[0], "material", 1, "psi")
        if a.white:
            z = np.random.default_rng(a.seed).standard_normal(len(tags))
            e = np.clip(SAN[2] + a.std * z, a.emin, a.emax)
            log(f"WHITE NOISE control: seed {a.seed}, std {a.std}")
        else:
            f = vf.field(cen, mean=SAN[2], std=a.std, theta=a.theta, seed=a.seed,
                         e_min=a.emin, e_max=a.emax)
            e = np.asarray(f.values, float)
            log(f"field: seed {a.seed}, std {a.std}, theta {a.theta} m "
                f"({a.theta/hx:.1f} elements), clipped {f.clipped_low}+{f.clipped_high}")
        vf.apply(ops, tags, e)
        log(f"field applied: e in [{e.min():.4f}, {e.max():.4f}], mean {e.mean():.4f}, "
            f"std {e.std():.4f}")
        # the element must forward the parameter: psi = e - e_c moves by exactly the change of e
        psi1 = ops.eleResponse(tags[0], "material", 1, "psi")
        if psi0 and psi1:
            d = (psi1[0] - psi0[0]) - (e[0] - SAN[2])
            log(f"forwarding check, element 1: d(psi) - d(e) = {d:.3e}")
            if abs(d) > 1.0e-12:
                raise SystemExit("LadrunoQuad did not forward voidRatio to its Gauss points")
        else:
            log("forwarding check: psi not readable before the first step")

    # ---- stage 0: confinement ----------------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        w = a.sig3 * hy * (0.5 if j in (0, ny) else 1.0)
        ops.load(nt(0, j), w, 0.0)
        ops.load(nt(nx, j), -w, 0.0)
    ops.load(MASTER, 0.0, -a.sig3 * W)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.updateMaterialStage("-material", MAT, "-stage", 0)
    ops.integrator("LoadControl", 0.1)
    for k in range(10):
        if ops.analyze(1) != 0:
            raise SystemExit(f"confinement step {k} failed")
    ops.loadConst("-time", 0.0)
    s = ops.eleResponse(tags[len(tags) // 2], "material", 1, "stress")
    log(f"after confinement: mid-element stress {['%.3f' % v for v in s]}")
    ops.updateMaterialStage("-material", MAT, "-stage", 1)

    # ---- stage 1: push -----------------------------------------------------
    uy0 = ops.nodeDisp(MASTER, 2)
    ops.loadConst("-time", uy0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.sp(MASTER, 2, 1.0)
    ptol = float(os.environ.get("BIAX_TOL", "1.0e-3")) * a.sig3 * W
    ladder = LADDER if not os.environ.get("BIAX_ALGO") else ((os.environ["BIAX_ALGO"], 1.0, 60),)
    ds_base = a.deps * H
    ds, ds_min, ds_max = ds_base, ds_base / 1000.0, ds_base
    bottom = [nt(i, 0) for i in range(nx + 1)]

    def sig_a():
        ops.reactions()
        return sum(ops.nodeReaction(n, 2) for n in bottom) / W

    def disp():
        return np.array([ops.nodeDisp(n) for n in range(1, len(xy) + 1)])

    if a.solver == "dr":
        return push_dr(a, log, logf, t0, uy0, sig_a, disp, xy, conn, cen, e, tags,
                       [nt(i, ny) for i in range(nx + 1)], hy)

    curve = [(0.0, sig_a())]
    snap_eps = np.linspace(a.eps_end / a.snaps, a.eps_end, a.snaps)
    snaps, snap_at = [], []
    good, nfail, status = 0, 0, "completed"
    relaxed = {}
    u_end = uy0 - a.eps_end * H
    while ops.getTime() > u_end + 1.0e-12:
        eps_now = (uy0 - ops.getTime()) / H
        if time.time() - t0 > a.wall:
            status = f"stopped: wall budget {a.wall:.0f} s at eps_a = {eps_now:.5f}"
            log(status)
            break
        k = len(snaps)
        target = uy0 - snap_eps[k] * H if k < len(snap_eps) else u_end
        step = min(ds, ops.getTime() - target)
        ok = False
        # the unbalance of a stalled step sits at 1 to 3 times the tolerance whatever the
        # step size (measured), so a deeply cut step may be accepted at a declared looser one
        tfs = [1.0]
        if a.tol_relax > 1.0 and step <= ds_base / 16.0:
            tfs.append(min(3.0, a.tol_relax))
        if a.tol_relax > 3.0 and step <= ds_base / 64.0:
            tfs.append(a.tol_relax)
        for tfx in tfs:
            for algo, tf, it in ladder:
                ops.integrator("LoadControl", -step)
                ops.test("NormUnbalance", tfx * tf * ptol, it, 0)
                ops.algorithm(algo)
                if ops.analyze(1) == 0:
                    ok = True
                    break
            if ok:
                if tfx > 1.0:
                    relaxed[tfx] = relaxed.get(tfx, 0) + 1
                break
        if not ok:
            nfail += 1
            good = 0
            ds = step / 4.0
            if ds < ds_min:
                status = f"stopped: step below {ds_min:.2e} m at eps_a = {eps_now:.5f}"
                log(status)
                break
            continue
        good += 1
        if good >= 6 and ds < ds_max:
            ds, good = min(2.0 * ds, ds_max), 0
        eps_now = (uy0 - ops.getTime()) / H
        curve.append((eps_now, sig_a()))
        if len(curve) % a.log_every == 0:
            log(f"  step {len(curve)-1}: eps_a = {eps_now:.5f}, sigma_a = {curve[-1][1]:.3f}, "
                f"ds = {ds:.2e}, cutbacks {nfail}, algo {algo}, iters {ops.testIter()}")
        if k < len(snap_eps) and eps_now >= snap_eps[k] - 1.0e-10:
            snaps.append(disp())
            snap_at.append(eps_now)
            log(f"eps_a = {eps_now:.4f}  sigma_a = {curve[-1][1]:9.3f} kPa  "
                f"steps {len(curve)-1}  cutbacks {nfail}")
    if len(snap_at) == 0 or abs(snap_at[-1] - curve[-1][0]) > 1.0e-9:
        snaps.append(disp())
        snap_at.append(curve[-1][0])

    # stall diagnostic: which Gauss points refuse the update at the minimum step
    sas0 = np.array([[ops.eleResponse(t, "material", g, "sasStats") for g in (1, 2, 3, 4)]
                     for t in tags])
    stall = np.zeros((0, 9))
    if status.startswith("stopped: step below"):
        ops.integrator("LoadControl", -ds_min)
        ops.test("NormUnbalance", ptol, 60, 0)
        ops.algorithm("Newton")
        ops.analyze(1)
        sas1 = np.array([[ops.eleResponse(t, "material", g, "sasStats") for g in (1, 2, 3, 4)]
                         for t in tags])
        dref = sas1[:, :, 15] - sas0[:, :, 15]
        rows = []
        for i, j in zip(*np.nonzero(dref > 0)):
            sg = ops.eleResponse(tags[i], "material", int(j) + 1, "stress")
            ps = ops.eleResponse(tags[i], "material", int(j) + 1, "psi")[0]
            rows.append((tags[i], j + 1, dref[i, j], sas1[i, j, 27], sas1[i, j, 26],
                         -(sg[0] + sg[1]) / 2.0, ps, cen[i, 0], cen[i, 1]))
        stall = np.array(rows).reshape(-1, 9)
        log(f"stall diagnostic, one more attempt at {ds_min:.1e} m: {len(rows)} Gauss points "
            f"refused (of {4*len(tags)}); columns ele gp refusals code lastSubsteps "
            f"-(sxx+syy)/2 psi x y")
        for r in rows[:12]:
            log("   " + " ".join(f"{v:.4g}" for v in r))
            for g in (1, 2, 3, 4):
                sg = ops.eleResponse(int(r[0]), "material", g, "stress")
                log(f"      ele {int(r[0])} gp {g} stress " + " ".join(f"{v:.4f}" for v in sg))
    tot = sas0.sum(axis=(0, 1))
    log("sasStats totals: updates %.0f, refusals %.0f (startF %.0f startAlpha %.0f other %.0f "
        "dTmin %.0f nonPosH %.0f lowP %.0f drift %.0f alpha %.0f cap %.0f), max substeps in "
        "one update %.0f" % (tot[0], tot[15], tot[16], tot[17], tot[18], tot[19], tot[20],
                             tot[21], tot[22], tot[23], tot[24], sas0[:, :, 25].max()))
    log("separation census: entries tension %.0f, low p %.0f, exits %.0f, points separated "
        "now %.0f" % (tot[36], tot[37], tot[38], tot[39]))
    log("steps accepted at a relaxed unbalance tolerance: " +
        (", ".join(f"{v} at {k:g}x" for k, v in sorted(relaxed.items())) or "none") +
        f" (1x = {ptol:.3e} kN)")
    c = census(tags)
    log(f"alpha_in census: {c[0]:.0f} commits, {c[1]:.0f} with alpha_in moved "
        f"({100*c[1]/max(c[0],1):.3f} %), sum |d alpha_in| = {c[2]:.4f}")
    # per-element state at the end
    e_now = np.array([np.mean([ops.eleResponse(t, "material", g, "psi")[0] for g in (1, 2, 3, 4)])
                      for t in tags])
    np.savez(os.path.join(a.out, "run.npz"), xy=xy, conn=conn, cen=cen, e=e,
             curve=np.array(curve), snaps=np.array(snaps), snap_at=np.array(snap_at),
             psi_end=e_now, W=W, H=H, nx=nx, ny=ny, sig3=a.sig3,
             seed=-1 if a.seed is None else a.seed, white=a.white, std=a.std, theta=a.theta,
             status=status, relaxed=np.array(sorted(relaxed.items()), float).reshape(-1, 2),
             ptol=ptol, stall=stall, sas=sas0.sum(axis=1))
    c = np.array(curve)
    ip = int(np.argmax(c[:, 1]))
    log(f"{status}: {len(curve)-1} steps, {nfail} cutbacks; peak sigma_a = {c[ip,1]:.3f} kPa "
        f"at eps_a = {c[ip,0]:.5f}; end sigma_a = {c[-1,1]:.3f} kPa at eps_a = {c[-1,0]:.5f}")
    logf.close()


def census(tags):
    """Sum of LadrunoSANISAND `commitStats` over the Gauss points: [commits, commits whose
    alpha_in moved (a loading reversal on the committed path), sum |d alpha_in|]."""
    c = np.zeros(3)
    for t in tags:
        for g in (1, 2, 3, 4):
            v = ops.eleResponse(t, "material", g, "commitStats")
            if v:
                c += np.array(v[:3], float)
    return c


def push_dr(a, log, logf, t0, uy0, sig_a, disp, xy, conn, cen, e, tags, top, hy):
    """Stage 1 by dynamic relaxation (ADR 21; recipe of note 83 sec. 4.3): the platen moves by
    `rate` per DR step, continuously, through a Linear series on EVERY top node (no equalDOF:
    the fictitious mass is sized on the unconstrained stiffness). A ramped point carries an
    inertial lag, so at intervals the platen is HELD and the state relaxed to rest; the held
    points are the static ones."""
    W, H = a.W, a.H
    ops.remove("loadPattern", 2)                  # the static push on the master node
    ops.remove("timeSeries", 2)
    for n in top[1:]:
        ops.remove("mp", n)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation", "-mass", "gershgorin", "-dt", a.rate,
                   "-damping", *(["kinetic"] if a.dr_damping == "kinetic" else ["viscous", a.zeta]),
                   "-massSafety", a.mass_safety)
    ops.analysis("Transient")
    res_tol = a.hold_tol * a.sig3 * hy

    def dr(what):
        return float(ops.ladrunoDR(what))

    def arm(u):                                   # u_y = -t on the top nodes, continuous
        ops.setTime(-u)
        ops.timeSeries("Linear", 4)
        ops.pattern("Plain", 2, 4)
        for n in top:
            ops.sp(n, 2, -1.0)

    def hold(u):
        """Freeze the platen at u and relax. Returns (sigma_a, steps, settled, res, ke)."""
        ops.remove("loadPattern", 2)
        ops.remove("timeSeries", 4)
        ops.timeSeries("Constant", 3)
        ops.pattern("Plain", 2, 3)
        for n in top:
            ops.sp(n, 2, float(u))
        n_h, quiet, res, ke, s_prev = 0, 0, float("nan"), float("nan"), None
        while n_h < a.hold_max:
            if ops.analyze(a.chunk, a.rate) != 0:
                break
            n_h += a.chunk
            res, ke = dr("residualNorm"), dr("kineticEnergy")
            s_now = sig_a()
            still = s_prev is not None and abs(s_now - s_prev) < 1.0e-4 * a.sig3
            quiet = quiet + 1 if (res < res_tol and still) else 0
            s_prev = s_now
            if quiet >= 3:
                break
        out = (sig_a(), n_h, quiet >= 3, res, ke)
        ops.remove("loadPattern", 2)
        ops.remove("timeSeries", 3)
        arm(u)
        return out

    arm(uy0)
    s0, n0, ok0, r0, k0 = hold(uy0)               # control: a zero push must not move the state
    log(f"DR control, zero push: sigma_a = {s0:.4f} kPa after {n0} steps, res {r0:.2e} kN, "
        f"KE {k0:.2e}, {'settled' if ok0 else 'NOT SETTLED'}; "
        f"stability margin {dr('stabilityMargin'):.3f}")
    log(f"DR push: {a.rate:.2e} m per step ({a.rate/H:.2e} strain), massSafety {a.mass_safety}, "
        f"hold gate res < {res_tol:.2e} kN")
    curve, holds = [(0.0, s0)], [(0.0, s0, s0, n0, int(ok0))]
    snap_eps = np.linspace(a.eps_end / a.snaps, a.eps_end, a.snaps)
    snaps, snap_at, status = [], [], "completed"
    nsteps = int(round(a.eps_end * H / a.rate))
    next_hold = a.hold_every if a.hold_every > 0 else 2.0
    n = 0
    while n < nsteps:
        if time.time() - t0 > a.wall:
            status = f"stopped: wall budget {a.wall:.0f} s at eps_a = {curve[-1][0]:.5f}"
            break
        m = min(a.chunk, nsteps - n)
        if ops.analyze(m, a.rate) != 0:
            status = f"stopped: DR analyze failed at eps_a = {curve[-1][0]:.5f}"
            break
        n += m
        u = ops.nodeDisp(top[0], 2)
        eps, sg = (uy0 - u) / H, sig_a()
        if not np.isfinite(sg):
            status = f"stopped: non-finite state at eps_a = {eps:.5f}"
            break
        curve.append((eps, sg))
        k = len(snaps)
        if k < len(snap_eps) and eps >= snap_eps[k] - 1.0e-10:
            snaps.append(disp())
            snap_at.append(eps)
        if eps >= next_hold - 1.0e-10:
            sh, nh, okh, rh, kh = hold(u)
            holds.append((eps, sg, sh, nh, int(okh)))
            next_hold += a.hold_every
            log(f"eps_a = {eps:.4f}  ramped {sg:9.3f} -> held {sh:9.3f} kPa "
                f"({100*(sg/sh-1):+.2f} % lag), {nh} steps, "
                f"{'settled' if okh else 'NOT SETTLED'}, res {rh:.2e}, "
                f"margin {dr('stabilityMargin'):.3f}")
        elif (len(curve) - 1) % a.log_every == 0:
            log(f"  step {n}/{nsteps}: eps_a = {eps:.5f}, sigma_a = {sg:.3f}, "
                f"res {dr('residualNorm'):.2e}, KE {dr('kineticEnergy'):.2e}")
    log(status)
    c = census(tags)
    log(f"alpha_in census: {c[0]:.0f} commits, {c[1]:.0f} with alpha_in moved "
        f"({100*c[1]/max(c[0],1):.3f} %), sum |d alpha_in| = {c[2]:.4f}")
    if np.isfinite(curve[-1][1]):
        u = ops.nodeDisp(top[0], 2)
        sh, nh, okh, rh, kh = hold(u)             # the last point is a held one
        holds.append((curve[-1][0], curve[-1][1], sh, nh, int(okh)))
        log(f"final hold at eps_a = {curve[-1][0]:.4f}: ramped {curve[-1][1]:.3f} -> held "
            f"{sh:.3f} kPa, {nh} steps, {'settled' if okh else 'NOT SETTLED'}")
        snaps.append(disp())
        snap_at.append(curve[-1][0])
    psi = np.array([np.mean([ops.eleResponse(t, "material", g, "psi")[0] for g in (1, 2, 3, 4)])
                    for t in tags])
    np.savez(os.path.join(a.out, "run.npz"), xy=xy, conn=conn, cen=cen, e=e,
             curve=np.array(curve), snaps=np.array(snaps), snap_at=np.array(snap_at),
             holds=np.array(holds), psi_end=psi, W=W, H=H, nx=a.nx, ny=a.ny, sig3=a.sig3,
             seed=-1 if a.seed is None else a.seed, white=a.white, std=a.std, theta=a.theta,
             status=status, solver="dr", rate=a.rate)
    h = np.array(holds)
    ip = int(np.argmax(h[:, 2]))
    log(f"{status}: {n} DR steps; held peak sigma_a = {h[ip,2]:.3f} kPa at eps_a = {h[ip,0]:.5f}; "
        f"end held {h[-1,2]:.3f} kPa at eps_a = {h[-1,0]:.5f}")
    logf.close()


if __name__ == "__main__":
    main()
