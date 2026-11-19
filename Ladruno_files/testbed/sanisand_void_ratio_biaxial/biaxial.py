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

    curve = [(0.0, sig_a())]
    snap_eps = np.linspace(a.eps_end / a.snaps, a.eps_end, a.snaps)
    snaps, snap_at = [], []
    good, nfail, status = 0, 0, "completed"
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
        for algo, tf, it in ladder:
            ops.integrator("LoadControl", -step)
            ops.test("NormUnbalance", tf * ptol, it, 0)
            ops.algorithm(algo)
            if ops.analyze(1) == 0:
                ok = True
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

    # per-element state at the end
    e_now = np.array([np.mean([ops.eleResponse(t, "material", g, "psi")[0] for g in (1, 2, 3, 4)])
                      for t in tags])
    np.savez(os.path.join(a.out, "run.npz"), xy=xy, conn=conn, cen=cen, e=e,
             curve=np.array(curve), snaps=np.array(snaps), snap_at=np.array(snap_at),
             psi_end=e_now, W=W, H=H, nx=nx, ny=ny, sig3=a.sig3,
             seed=-1 if a.seed is None else a.seed, white=a.white, std=a.std, theta=a.theta,
             status=status)
    c = np.array(curve)
    ip = int(np.argmax(c[:, 1]))
    log(f"{status}: {len(curve)-1} steps, {nfail} cutbacks; peak sigma_a = {c[ip,1]:.3f} kPa "
        f"at eps_a = {c[ip,0]:.5f}; end sigma_a = {c[-1,1]:.3f} kPa at eps_a = {c[-1,0]:.5f}")
    logf.close()


if __name__ == "__main__":
    main()
