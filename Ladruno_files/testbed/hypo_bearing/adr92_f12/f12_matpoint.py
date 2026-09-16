"""WP-105 / F12 phase (a): IntScheme 2 (BackwardEuler_CPPM) vs IntScheme 1 on
LadrunoSANISAND, at the ADR-92 P0/G0 material-point certificate's own instrument.

Adapted from Ladruno_implementation/adr92_p0_oracle/probe_binary_triaxial.py --
same cube, same D-L constants, same -Pmin/-Presidual/-honorTolR, same
consolidate-then-push. Additions: wall clock per step, the `substeps` response
(a CPPM->explicit_integrator fallback detector that needs no source edit), the
`psi` response (for M^b), and a prescribed volumetric-extension path that drives
p onto the p_min floor (the ADR-93 ring regime).

NO CODE CHANGES ANYWHERE. Read-only use of the 634824e1f binary.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

# WP-105 / F12: was a hardcoded absolute scratchpad path. Default now
# resolves to <fork root>/dist/bin (the build.bat REQUIRED layout, see
# CLAUDE.md), computed from this script's own location four levels up
# (Ladruno_files/testbed/hypo_bearing/adr92_f12/); override with
# F12_DIST_BIN to point at a specific worktree's dist/bin.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DIST = os.environ.get("F12_DIST_BIN", os.path.join(REPO_ROOT, "dist", "bin"))
sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

EXPECT_BUILD = "634824e1f"

CONSTS = dict(
    G0=264.32, nu=0.3129, e_init=0.6944, Mc=1.3309, c=0.71, lambda_c=0.027,
    e0=0.83, ksi=0.45, P_atm=101.0, m=0.005, h0=1.3, ch=0.968, nb=3.5,
    A0=0.05, nd=5.75, z_max=12.5, cz=1100.0, rho=2.0,
)
ORDER = ["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi", "P_atm", "m",
         "h0", "ch", "nb", "A0", "nd", "z_max", "cz", "rho"]
PMIN = 0.0101
PRESIDUAL = 0.0
HONOR_TOLR = 0

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def check_build():
    h = ops.ladrunoBuild()
    if not h.startswith(EXPECT_BUILD):
        raise SystemExit("F12 FAIL-LOUD: wrong build " + repr(h))
    return h


MAX_SUBSTEPS = 0   # set from the CLI; 0 = uncapped = vanilla
PUSH_TOL = 1.0e-9  # the ADR-92 certificate probe's own NormDispIncr tolerance


def build_model(p0, e_init, scheme, tan_type, tol):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in enumerate(NODES, start=1):
        ops.node(tag, float(x), float(y), float(z))
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 0, 0)
    for n in (1, 2, 5, 6):
        ops.fix(n, 0, 1, 0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 0, 1)
    vals = dict(CONSTS)
    vals["e_init"] = e_init
    args = [vals[k] for k in ORDER]
    flags = ["-Presidual", PRESIDUAL, "-Pmin", PMIN, "-honorTolR", HONOR_TOLR]
    if MAX_SUBSTEPS:
        flags += ["-maxSubsteps", int(MAX_SUBSTEPS)]
    ops.nDMaterial("LadrunoSANISAND", 1, *args, scheme, tan_type, 1, tol, tol, *flags)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", "bbar")
    q4 = -p0 / 4.0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, q4, 0.0, 0.0)
    for n in (3, 4, 7, 8):
        ops.load(n, 0.0, q4, 0.0)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, q4)


def consolidate():
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        raise SystemExit("consolidation FAILED")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.integrator("LoadControl", 0.0)
    if ops.analyze(5) != 0:
        raise SystemExit("stage-1 re-equilibration FAILED")
    ops.loadConst("-time", 0.0)


def mat_probe():
    sig = list(ops.eleResponse(1, "material", 1, "stress"))
    eps = list(ops.eleResponse(1, "material", 1, "strain"))
    try:
        sub = list(ops.eleResponse(1, "material", 1, "substeps"))
    except Exception:
        sub = [float("nan"), float("nan")]
    if len(sub) < 2:
        sub = [float("nan"), float("nan")]
    try:
        psi = list(ops.eleResponse(1, "material", 1, "psi"))
    except Exception:
        psi = [float("nan")]
    if not psi:
        psi = [float("nan")]
    return sig, eps, sub, psi


def invariants(sig):
    """OpenSees sign convention is tension-positive; return compression-positive p."""
    sm = (sig[0] + sig[1] + sig[2]) / 3.0
    p = -sm
    dev = [sig[0] - sm, sig[1] - sm, sig[2] - sm, sig[3], sig[4], sig[5]]
    nrm = math.sqrt(dev[0] ** 2 + dev[1] ** 2 + dev[2] ** 2
                    + 2.0 * (dev[3] ** 2 + dev[4] ** 2 + dev[5] ** 2))
    q = math.sqrt(1.5) * nrm
    eta = q / p if p > 1.0e-12 else float("nan")
    return p, q, eta


def Mb_from(p, e):
    """Bounding stress ratio in triaxial compression: M_b = Mc*exp(-nb*psi)."""
    ec = CONSTS["e0"] - CONSTS["lambda_c"] * (max(p, PMIN) / CONSTS["P_atm"]) ** CONSTS["ksi"]
    psi = e - ec
    return CONSTS["Mc"] * math.exp(-CONSTS["nb"] * psi), psi


def _row(i, e_init, dt_, it):
    sig, eps, sub, psi = mat_probe()
    p, q, eta = invariants(sig)
    e_cur = e_init - (1.0 + e_init) * (eps[0] + eps[1] + eps[2])
    Mb, psi_calc = Mb_from(p, e_cur)
    r = dict(step=i, t=ops.getTime(), wall=dt_, iters=it, p=p, q=q, eta=eta,
             Mb=Mb, psi_calc=psi_calc, psi_resp=psi[0], e=e_cur,
             subME=sub[0], capHit=sub[1])
    for k in range(6):
        r["sig%d" % k] = sig[k]
    for k in range(6):
        r["eps%d" % k] = eps[k]
    return r


def _dump(rows, out):
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _loop(nstep, dstep, e_init, out):
    rows = []
    done = 0
    stalled = None
    t_push = 0.0
    for i in range(nstep + 1):
        it = 0
        dt_ = 0.0
        if i > 0:
            t0 = time.perf_counter()
            rc = ops.analyze(1)
            dt_ = time.perf_counter() - t0
            t_push += dt_
            try:
                it = ops.testIter()
            except Exception:
                it = -1
            if rc != 0:
                stalled = i
                break
            done = i
        rows.append(_row(i, e_init, dt_, it))
    _dump(rows, out)
    return rows, done, stalled, t_push


def run_triaxial(p0, e_init, scheme, tan_type, tol, nstep, ez_max, out):
    build_model(p0, e_init, scheme, tan_type, tol)
    consolidate()
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for n in (5, 6, 7, 8):
        ops.sp(n, 3, -1.0)
    du = ez_max / nstep
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", PUSH_TOL, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", du)
    ops.analysis("Static")
    return _loop(nstep, du, e_init, out)


def run_floor(p0, e_init, scheme, tan_type, tol, nstep, out,
              evol_max=2.2e-4, gam_max=6.0e-3):
    """ADR-93 ring regime: prescribed volumetric EXTENSION ramped to evol_max by
    t = 0.5 and then held, with deviatoric shear gam_max*t on throughout, from a
    p0 committed state -- p collapses onto the p_min floor while the point keeps
    flowing. Every face normal is prescribed (zero free DOF)."""
    build_model(p0, e_init, scheme, tan_type, tol)
    consolidate()
    ts = [i / float(nstep) for i in range(nstep + 1)]
    ramp = [min(2.0 * t, 1.0) for t in ts]
    ex = [evol_max / 3.0 * r + gam_max / 3.0 * t for r, t in zip(ramp, ts)]
    ez = [evol_max / 3.0 * r - 2.0 * gam_max / 3.0 * t for r, t in zip(ramp, ts)]
    # A Path series returns 0 past its last time point, and LoadControl's
    # accumulated pseudo-time overshoots 1.0 by an ulp on the last step -- which
    # would silently unload the whole path in one increment. Extend flat.
    ts = ts + [2.0]
    ex = ex + [ex[-1]]
    ez = ez + [ez[-1]]
    ops.timeSeries("Path", 3, "-time", *ts, "-values", *ex)
    ops.timeSeries("Path", 4, "-time", *ts, "-values", *ez)
    ops.pattern("Plain", 3, 3)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, 1.0)
    for n in (3, 4, 7, 8):
        ops.sp(n, 2, 1.0)
    ops.pattern("Plain", 4, 4)
    for n in (5, 6, 7, 8):
        ops.sp(n, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", PUSH_TOL, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.analysis("Static")
    return _loop(nstep, 1.0 / nstep, e_init, out)


def run_replay(p0, e_init, scheme, tan_type, tol, nstep, out, path_csv):
    """Strain-path REPLAY -- the ADR-92 G0 discipline (the oracle replays the
    binary's own recorded strain sequence), here with scheme 2 in the oracle's
    seat. Every face normal is prescribed from the recorded path, so there is NO
    global Newton and NO stall: the only thing that differs between the two arms
    is the constitutive integrator. Wall time is then pure constitutive cost."""
    with open(path_csv, newline="", encoding="utf-8") as fh:
        rec = list(csv.DictReader(fh))
    ex_r = [float(r["eps0"]) for r in rec]
    ey_r = [float(r["eps1"]) for r in rec]
    ez_r = [float(r["eps2"]) for r in rec]
    n_r = len(ex_r) - 1

    def interp(arr, u):           # u in [0,1] along the recorded step index
        x = u * n_r
        i = min(int(x), n_r - 1)
        f = x - i
        return arr[i] * (1.0 - f) + arr[i + 1] * f

    build_model(p0, e_init, scheme, tan_type, tol)
    consolidate()
    ts = [i / float(nstep) for i in range(nstep + 1)]
    ex = [interp(ex_r, u) for u in ts]
    ey = [interp(ey_r, u) for u in ts]
    ez = [interp(ez_r, u) for u in ts]
    ts = ts + [2.0]
    ex, ey, ez = ex + [ex[-1]], ey + [ey[-1]], ez + [ez[-1]]
    ops.timeSeries("Path", 3, "-time", *ts, "-values", *ex)
    ops.timeSeries("Path", 5, "-time", *ts, "-values", *ey)
    ops.timeSeries("Path", 4, "-time", *ts, "-values", *ez)
    ops.pattern("Plain", 3, 3)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, 1.0)
    ops.pattern("Plain", 5, 5)
    for n in (3, 4, 7, 8):
        ops.sp(n, 2, 1.0)
    ops.pattern("Plain", 4, 4)
    for n in (5, 6, 7, 8):
        ops.sp(n, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", PUSH_TOL, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.analysis("Static")
    return _loop(nstep, 1.0 / nstep, e_init, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["tx", "floor", "replay"], default="tx")
    ap.add_argument("--path-csv", default=None)
    ap.add_argument("--p0", type=float, default=100.0)
    ap.add_argument("--e-init", type=float, default=CONSTS["e_init"])
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--tan-type", type=int, default=2)
    ap.add_argument("--tol", type=float, default=1.0e-10)
    ap.add_argument("--nstep", type=int, default=200)
    ap.add_argument("--ez-max", type=float, default=0.02)
    ap.add_argument("--tag", default="")
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data"))
    ap.add_argument("--max-substeps", type=int, default=0)
    ap.add_argument("--push-tol", type=float, default=1.0e-9)
    a = ap.parse_args()
    global MAX_SUBSTEPS, PUSH_TOL
    MAX_SUBSTEPS = a.max_substeps
    PUSH_TOL = a.push_tol
    h = check_build()
    os.makedirs(a.outdir, exist_ok=True)
    name = a.tag or ("%s_p%g_e%g_s%d_T%d_n%d"
                     % (a.kind, a.p0, a.e_init, a.scheme, a.tan_type, a.nstep))
    out = os.path.join(a.outdir, name + ".csv")
    t0 = time.perf_counter()
    if a.kind == "tx":
        rows, done, stalled, t_push = run_triaxial(
            a.p0, a.e_init, a.scheme, a.tan_type, a.tol, a.nstep, a.ez_max, out)
    elif a.kind == "replay":
        rows, done, stalled, t_push = run_replay(
            a.p0, a.e_init, a.scheme, a.tan_type, a.tol, a.nstep, out, a.path_csv)
    else:
        rows, done, stalled, t_push = run_floor(
            a.p0, a.e_init, a.scheme, a.tan_type, a.tol, a.nstep, out)
    tot = time.perf_counter() - t0
    etas = [r["eta"] for r in rows if r["eta"] == r["eta"]]
    its = [r["iters"] for r in rows if r["iters"] > 0]
    summ = dict(build=h, kind=a.kind, p0=a.p0, e_init=a.e_init, scheme=a.scheme,
                tan_type=a.tan_type, tol=a.tol, nstep=a.nstep, ez_max=a.ez_max,
                steps_done=done, stalled_at=stalled, wall_push_s=t_push,
                wall_total_s=tot, ms_per_step=1000.0 * t_push / max(done, 1),
                path_csv=a.path_csv, push_tol=a.push_tol,
                max_substeps=a.max_substeps,
                ms_per_step_ok=1000.0 * sum(r["wall"] for r in rows) / max(done, 1),
                wall_failed_step=(t_push - sum(r["wall"] for r in rows)),
                eta_peak=max(etas) if etas else None,
                eta_end=rows[-1]["eta"], p_end=rows[-1]["p"], q_end=rows[-1]["q"],
                Mb_end=rows[-1]["Mb"], e_end=rows[-1]["e"],
                p_min_along=min(r["p"] for r in rows),
                iters_total=sum(its), iters_mean=(float(sum(its)) / max(len(its), 1)),
                iters_max=(max(its) if its else 0),
                subME_max=max(r["subME"] for r in rows),
                subME_nonzero_steps=sum(1 for r in rows if r["subME"] > 0),
                capHit_steps=sum(1 for r in rows if r["capHit"] > 0),
                out=out)
    with open(out.replace(".csv", ".json"), "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1)
    print("@@SUMMARY " + json.dumps(summ))


if __name__ == "__main__":
    main()
