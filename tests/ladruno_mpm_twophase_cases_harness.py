"""WP-232l route-E cases (T2-T5, the A4 patch checks), from the reference blocks of
tests/mpm_twophase_cases.txt and 232l-V's records under tests/data/mpm/twophase/.

Never reads tests/mpm_twophase_reference.py. The closed forms used here are the ones the blocks
name: Terzaghi's series (the formulas of tests/test_ladruno_up_element_analytic.py:89-102), the
hydrostatic column (block t4_ref_hydro), and for Mandel the record's own p_exact at the particles.

Standalone:
    python tests/ladruno_mpm_twophase_cases_harness.py --bin <dist\\bin> <leg> [<leg> ...]
legs: patch, patch_fix2, patch_fix2_controls, t2, t2_reported, t3_gimp, t3_bspline2, t4,
      t4_reported, t5, diagnostic
Every result is appended to tests/data/mpm/twophase/route_e_results.json with ladrunoBuild.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "mpm" / "twophase"
OUT = DATA / "route_e_results.json"


# ---------------------------------------------------------------------------
#  Closed forms (block t2_ref_terzaghi: test_ladruno_up_element_analytic.py:89-102)
# ---------------------------------------------------------------------------
def terzaghi_series(Z, Tv, nterms=300):
    s = 0.0
    for m in range(nterms):
        M = (2 * m + 1) * math.pi / 2.0
        s += (2.0 / M) * math.sin(M * Z) * math.exp(-M * M * Tv)
    return s


def terzaghi_U(Tv, nterms=300):
    s = 0.0
    for m in range(nterms):
        M = (2 * m + 1) * math.pi / 2.0
        s += (2.0 / M ** 2) * math.exp(-M * M * Tv)
    return 1.0 - s


def as_list(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


# ---------------------------------------------------------------------------
#  One deck
# ---------------------------------------------------------------------------
def build(ops, s):
    """s: ndm 2 deck spec (keys below). Mode A, Newmark, UmfPack, Transformation."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    E, nu, rho = s["mat"]
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", s["h"], "-cells", s["cells"][0], s["cells"][1],
            "-basis", s["basis"])
    ops.mpm("mode", s.get("mode", "implicit"))
    if s.get("kin", "hypo") != "hypo":
        ops.mpm("kinematics", s["kin"])
    lo, hi = s["box"]
    ops.mpm("body", 1, "box", lo[0], lo[1], hi[0], hi[1], "-ppc", s.get("ppc", 2), "-mat", 1, "-rho", rho)
    if s.get("fluid") is not None:
        Kf, n, rf, k, stab = s["fluid"]
        sa = ["-stab", "off"] if stab == "off" else ["-stab", "auto", float(stab[1])]
        fb = ["-fluidBody", *s["fluidBody"]] if s.get("fluidBody") is not None else []
        ops.mpm("twoPhase", 1, "-Kf", Kf, "-poro", n, "-rhoF", rf, "-perm", k, k, *sa, *fb)
    if s.get("gravity"):
        ops.mpm("gravity", 0.0, -s["gravity"])
    for face, comps in s.get("fix", []):
        ops.mpm("fix", face, comps)
    for face, p0 in s.get("pfix", []):
        ops.mpm("pfix", face, p0)
    if s.get("plate"):
        ops.mpm("plate", s["plate"][0], s["plate"][1])
    for name, face, tlo, thi, t in s.get("traction", []):
        ops.mpm("traction", name, face, *tlo, *thi, *t)
    if s.get("path"):
        tt, vv = s["path"]
        ops.timeSeries("Path", 77, "-time", *tt, "-values", *vv)
        ops.mpm("timeSeries", 77)
    if s.get("rayleigh"):
        ops.mpm("rayleigh", *s["rayleigh"])
    if s.get("safety"):
        ops.mpm("safety", s["safety"])
    if s.get("mode") == "explicit":           # route S: the engine's own loop, no FE analysis
        return
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    if s.get("mode", "implicit") == "quasiStatic":
        ops.test("NormDispIncr", 1e-12, 25)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", s["dt"])
        ops.analysis("Static")
    else:
        if s.get("linear"):
            ops.algorithm("Linear")
        else:
            ops.test("NormDispIncr", 1e-10, 25)
            ops.algorithm("Newton")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")


def pdata(ops, f, w=1):
    v = as_list(ops.mpm("particleData", f))
    return [v[i:i + w] for i in range(0, len(v), w)] if w > 1 else v


def plate_tag(ops):
    grid = set(int(t) for t in as_list(ops.mpm("info", "gridNodeTags")))
    rest = [t for t in ops.getNodeTags() if t not in grid]
    return rest[0] if len(rest) == 1 else None


# ---------------------------------------------------------------------------
#  A4 patch checks (P9's p-image fold, P10's plate fold)
# ---------------------------------------------------------------------------
def patch_check(ops, basis, tmp):
    """(a) the plate fold carries an exact uniform state: a single-phase (drained) elastic block,
    quasi-static, small strain, rollers on xmin and ymin, xmax traction-free, plate on ymax: one
    load step must give sigma_yy = -F/a and sigma_xx = 0 at every particle, and every particle's
    u_y is linear in its Y, at rel 1e-10. (b) the folds reproduce the fields the images claim:
    at every particle of a two-phase block with `pfix xmax 0` (a pfix-only face) and the plate,
    sum S^p (X - a) = x - a and sum dS^p (X - a) = 1 (the odd p image), and the plate's u_y
    weights reproduce U + e (y - b) (the affine odd image), from `mpm debug particleBlocks`."""
    a, b, h, F = 4.0, 2.0, 0.5, 10.0
    E, nu = 1.0e4, 0.1
    spec = dict(mat=(E, nu, 2.0), h=h, cells=(8, 4), basis=basis, box=((0.0, 0.0), (a, b)), mode="quasiStatic",
                kin="small", fix=[("xmin", "x"), ("ymin", "y")], plate=("ymax", F), dt=1.0)
    build(ops, spec)
    # the folds' algebra, isolated from renormalize: no support node is dropped (a renormalised
    # support does not reproduce linear fields; that is P4's choice, not a fold property)
    ops.mpm("massCutoff", 1e-14)
    X0 = pdata(ops, "x", 2)
    ops.mpm("analyze", 1, 1.0)
    sig = pdata(ops, "sig", 3)
    x1 = pdata(ops, "x", 2)
    syy = -F / a
    ea = max(abs(sg[1] - syy) for sg in sig) / abs(syy)
    ex = max(abs(sg[0]) for sg in sig) / abs(syy)
    # plane strain, sigma_xx = 0: eps_yy = syy (1 - nu^2)/E, u_y = eps_yy Y (roller at y = 0)
    eyy = syy * (1 - nu * nu) / E
    uy = max(abs((x1[k][1] - X0[k][1]) - eyy * X0[k][1]) for k in range(len(X0))) / abs(eyy * b)
    res = {"basis": basis, "syy_rel": ea, "sxx_rel": ex, "uy_rel": uy}
    # (b) fold reproduction from the element dumps
    spec2 = dict(spec, mode="implicit", fluid=(4444.0, 0.4, 1.0, 1.0e-4, "off"), pfix=[("xmax", 0.0)],
                 dt=0.01)
    build(ops, spec2)
    ops.mpm("massCutoff", 1e-14)
    n = len(X0)
    worst_p = worst_g = 0.0
    for p in range(n):
        # the odd fold at xmax only: xmin is impervious (an even image, which by design does not
        # carry a field with dp/dx != 0 there), so particles whose support reaches xmin are skipped
        if X0[p][0] < 2.0 * h:
            continue
        path = str(Path(tmp) / ("patch_%s_%d.json" % (basis, p)))
        ops.mpm("debug", "particleBlocks", p, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        d = json.loads(Path(path).read_text())
        xp = X0[p][0]
        sp = sum(d["Sp"][k] * (d["coords"][k][0] - a) for k in range(d["nen"]))
        gp = sum(d["dSp"][k][0] * (d["coords"][k][0] - a) for k in range(d["nen"]))
        worst_p = max(worst_p, abs(sp - (xp - a)) / a)
        worst_g = max(worst_g, abs(gp - 1.0))
    res.update({"p_odd_value_rel": worst_p, "p_odd_grad_rel": worst_g})
    res["pass"] = max(ea, ex, uy, worst_p, worst_g) <= 1e-10
    return res


# ---------------------------------------------------------------------------
#  T2 Terzaghi (block t2_ref_terzaghi)
# ---------------------------------------------------------------------------
T2_STEPS = {0.05: 450, 0.1: 900, 0.2: 1800, 0.5: 4500}


def run_t2(ops, basis, stab, kin="hypo", upto=None):
    w = 2.0 if basis.startswith("bspline") else 1.0
    L, q, dt = 10.0, 10.0, 0.02
    E, nu, Kf, n, kbar = 1.0e4, 0.2, 4444.0, 0.4, 1.0e-4
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    Qbar = Kf / n
    cv = kbar / (1.0 / Eoed + n / Kf)
    p0 = q * Qbar / (Qbar + Eoed)
    s0 = q * L / (Eoed + Qbar)
    sinf = q * L / Eoed
    spec = dict(mat=(E, nu, 2.0), h=1.0, cells=(int(w), 10), basis=basis, box=((0.0, 0.0), (w, L)),
                kin=kin, fluid=(Kf, n, 1.0, kbar, stab),
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
                plate=("ymax", q * w))
    t0 = time.time()
    build(ops, spec)
    Y0 = [x[1] for x in pdata(ops, "x", 2)]
    pt = None
    U = 0.0
    rows = {}
    last = max(T2_STEPS.values()) if upto is None else T2_STEPS[upto]
    for step in range(1, last + 1):
        ops.mpm("analyze", 1, dt)
        if pt is None:
            pt = plate_tag(ops)
        U += -ops.nodeDisp(pt, 2)
        for Tv, st in T2_STEPS.items():
            if st == step:
                t = step * dt
                Tva = cv * t / (L * L)
                pp = pdata(ops, "p")
                perr = max(abs(pp[k] - p0 * terzaghi_series((L - Y0[k]) / L, Tva)) for k in range(len(pp))) / p0
                sx = s0 + (sinf - s0) * terzaghi_U(Tva)
                rows[str(Tv)] = {"Tv_actual": Tva, "p_err_max": perr, "U_err": abs(U - sx) / (sinf - s0), "U": U}
    gated = {"0.1", "0.2", "0.5"} & set(rows)
    ok = all(rows[k]["p_err_max"] <= 0.02 and rows[k]["U_err"] <= 0.02 for k in gated)
    return {"basis": basis, "stab": stab, "kin": kin, "rows": rows, "pass_bands_0.02": ok,
            "wall_s": time.time() - t0}


# ---------------------------------------------------------------------------
#  T7 route S time step (U6, A6): T6's Terzaghi column, the drained pencil as the yardstick
# ---------------------------------------------------------------------------
def t6_terzaghi_spec(basis, stab, safety=None):
    w = 2.0 if basis.startswith("bspline") else 1.0
    return dict(mat=(1.0e4, 0.2, 2.0), h=1.0, cells=(int(w), 10), basis=basis, box=((0.0, 0.0), (w, 10.0)),
                mode="explicit", fluid=(4444.0, 0.4, 1.0, 1.0e-4, stab),
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
                plate=("ymax", 10.0 * w), path=((0.0, 5.0, 1.0e6), (0.0, 1.0, 1.0)), safety=safety)


def run_t7(ops, basis, stab=("auto", 0.25), nsteps=2000):
    """Stable at 0.9 x the drained pencil for 2000 steps, unstable or refused at 1.5 x (A6). The
    route-S run-time limit and its cap are lifted (`mpm debug routeSLimit off`, debug only), so
    the legs measure the scheme itself (review of record)."""
    build(ops, t6_terzaghi_spec(basis, stab))
    dt_d = as_list(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
    dt_u = as_list(ops.mpm("criticalTimeStep", "undrainedPencil"))[0]
    dt_est = as_list(ops.mpm("criticalTimeStep"))[0]

    def leg(fac):
        build(ops, t6_terzaghi_spec(basis, stab))
        ops.mpm("debug", "routeSLimit", "off")      # review of record: measure above the cap
        dt = fac * dt_d
        vmax = pmax = 0.0
        for s in range(nsteps):
            try:
                ops.mpm("analyze", 1, dt)
            except Exception:
                return {"factor": fac, "dt": dt, "steps": s, "refused": True, "blew_up": False,
                        "vmax": vmax, "pmax": pmax}
            vm = max(abs(x) for x in pdata(ops, "v"))
            pm = max(abs(x) for x in pdata(ops, "p"))
            if not (vm == vm and pm == pm) or vm > 1.0e6 or pm > 1.0e8:
                return {"factor": fac, "dt": dt, "steps": s + 1, "refused": False, "blew_up": True,
                        "vmax": vm, "pmax": pm}
            vmax, pmax = max(vmax, vm), max(pmax, pm)
        return {"factor": fac, "dt": dt, "steps": nsteps, "refused": False, "blew_up": False,
                "vmax": vmax, "pmax": pmax}

    t0 = time.time()
    a, b = leg(0.9), leg(1.5)
    # stable: 2000 steps, bounded (|v| below 1 m/s, |p| below 10 x the load q = 10 kPa)
    stable = (a["steps"] == nsteps and not a["blew_up"] and not a["refused"] and a["vmax"] < 1.0
              and a["pmax"] < 100.0)
    unstable = b["refused"] or b["blew_up"] or b["vmax"] > 100.0 * max(a["vmax"], 1e-30)
    r = {"basis": basis, "stab": stab, "dt_drained_pencil": dt_d, "dt_undrained_pencil": dt_u,
         "dt_estimator_h_over_c": dt_est, "undrained_over_drained": dt_u / dt_d, "at_0.9": a,
         "at_1.5": b, "pass": bool(stable and unstable)}
    if a["refused"] and a["steps"] == 0 and dt_est < 0.9 * dt_d:
        # the estimator (s = 1) refuses 0.9 x the pencil: diagnostic only, at the estimator's own limit
        r["diagnostic_at_estimator"] = leg(0.999 * dt_est / dt_d)
    r["wall_s"] = time.time() - t0
    return r


# ---------------------------------------------------------------------------
#  T6 route S (block t6_ref_ramp): T2's Terzaghi and T3's Mandel under the 5 s ramp, mode F
# ---------------------------------------------------------------------------
def terzaghi_ramp(Z, t, cv, L, tr, nterms=300):
    """p / (B q) under q(t) = q min(t / tr, 1) (Duhamel on the step series, written here from
    dp/dt = cv p'' + B dq/dt with p = 0 at the drained top, Z = depth / L, impervious base)."""
    te = min(t, tr)
    s = 0.0
    for m in range(nterms):
        M = (2 * m + 1) * math.pi / 2.0
        lam = M * M * cv / (L * L)
        s += (2.0 / M) * math.sin(M * Z) * (math.exp(-lam * (t - te)) - math.exp(-lam * t)) / (lam * tr)
    return s


def terzaghi_settle_ramp(t, cv, L, tr, q, B, Eoed, nterms=300):
    """settlement = (q(t) L - int_0^L p dz) / Eoed under the ramp (compression positive)."""
    te = min(t, tr)
    ip = 0.0
    for m in range(nterms):
        M = (2 * m + 1) * math.pi / 2.0
        lam = M * M * cv / (L * L)
        ip += (2.0 * L / (M * M)) * (math.exp(-lam * (t - te)) - math.exp(-lam * t)) / (lam * tr)
    return (q * min(t / tr, 1.0) * L - B * q * ip) / Eoed


def run_t6_terzaghi(ops, basis, stab, fe, dt=0.005, upto=None):
    """T6 Terzaghi, route S: p at the particles' initial depths against the ramp closed form, and the
    settlement read as the particle quadrature sum_p (F_yy - 1) V0_p / width (mode F has no plate
    node to read U from). Bands: the FE record's band_p, band_U per Tv (t6_fe_dt<dt>.json)."""
    L, q, tr = 10.0, 10.0, 5.0
    E, nu, Kf, n, kbar = 1.0e4, 0.2, 4444.0, 0.4, 1.0e-4
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    Qbar = Kf / n
    cv = kbar / (1.0 / Eoed + n / Kf)
    B = Qbar / (Qbar + Eoed)
    p0 = q * B
    s0 = q * L / (Eoed + Qbar)
    sinf = q * L / Eoed
    w = 2.0 if basis.startswith("bspline") else 1.0
    t0 = time.time()
    build(ops, t6_terzaghi_spec(basis, stab))
    X0 = pdata(ops, "x", 2)
    Y0 = [x[1] for x in X0]
    V0 = pdata(ops, "V")
    steps = {str(r["Tv"]): (int(round(r["t"] / dt)), i) for i, r in enumerate(fe["rows"])}
    rows = {}
    done = 0
    if upto is not None:
        steps = {k: v for k, v in steps.items() if float(k) <= upto}
    for Tv, (st, i) in sorted(steps.items(), key=lambda kv: kv[1][0]):
        ops.mpm("analyze", st - done, dt)
        done = st
        t = st * dt
        pp = pdata(ops, "p")
        F = pdata(ops, "F", 4)
        U = -sum((F[k][3] - 1.0) * V0[k] for k in range(len(pp))) / w
        perr = max(abs(pp[k] - p0 * terzaghi_ramp((L - Y0[k]) / L, t, cv, L, tr)) for k in range(len(pp))) / p0
        sx = terzaghi_settle_ramp(t, cv, L, tr, q, B, Eoed)
        rows[Tv] = {"t": t, "Tv_actual": cv * t / (L * L), "p_err_max": perr, "U": U, "U_exact": sx,
                    "U_err": abs(U - sx) / (sinf - s0), "band_p": fe["band_p"][i], "band_U": fe["band_U"][i],
                    "fe_p_err": fe["rows"][i]["p_err_gimp_points"], "fe_U_err": fe["rows"][i]["U_err"]}
    gated = [k for k in rows if k != "0.05"]
    ok = all(rows[k]["p_err_max"] <= rows[k]["band_p"] and rows[k]["U_err"] <= rows[k]["band_U"] for k in gated)
    return {"case": "terzaghi", "basis": basis, "stab": stab, "dt": dt, "rows": rows, "pass": ok,
            "gated": basis.startswith("bspline"), "wall_s": time.time() - t0}


def run_t6_mandel(ops, basis, stab, fe, peak_gate, dt=0.005):
    a, b, h, sig0 = 10.0, 4.0, 1.0, 10.0
    p0 = fe["p0"]
    spec = dict(mat=(1.0e4, 0.1, 2.0), h=h, cells=(10, 4), basis=basis, box=((0.0, 0.0), (a, b)),
                mode="explicit", fluid=(4444.0, 0.4, 1.0, 1.0e-4, stab), fix=[("xmin", "x"), ("ymin", "y")],
                pfix=[("xmax", 0.0)], plate=("ymax", a * sig0), path=((0.0, 5.0, 1.0e6), (0.0, 1.0, 1.0)))
    t0 = time.time()
    build(ops, spec)
    X0 = pdata(ops, "x", 2)
    col = [k for k, x in enumerate(X0) if x[0] < h / 2.0]
    want = {int(round(r["t"] / dt)): i for i, r in enumerate(fe["rows"])}
    rows = {}
    peak = -1e30
    for step in range(1, max(want) + 1):
        ops.mpm("analyze", 1, dt)
        pp = pdata(ops, "p")
        pc = sum(pp[k] for k in col) / len(col) / p0
        peak = max(peak, pc)
        if step in want:
            r = fe["rows"][want[step]]
            rows[str(r["tv"])] = {"p_center": pc, "p_exact_quarter": r["p_exact_quarter_over_p0"],
                                  "err": abs(pc - r["p_exact_quarter_over_p0"]), "band": fe["band"][want[step]],
                                  "fe_err": r["err_quarter_over_p0"]}
    ok = all(v["err"] <= v["band"] for v in rows.values()) and peak >= peak_gate
    return {"case": "mandel", "basis": basis, "stab": stab, "dt": dt, "rows": rows, "peak": peak,
            "peak_gate": peak_gate, "pass": ok, "gated": basis.startswith("bspline"),
            "wall_s": time.time() - t0}


# ---------------------------------------------------------------------------
#  T13 translation patch (pre-registered in the review of record; block t13_translation)
# ---------------------------------------------------------------------------
def run_t13(ops, mode, basis, plate, steps=200, dt=2.0e-3, vy=-1.0):
    """A free two-phase block in rigid translation (v = (0, vy), no gravity, sealed, p = 0: a
    uniform p would not be in equilibrium with its traction-free faces) keeps v, sigma' and p:
    rel |v - v0| / |v0|, |sigma'| / E and |p| / E after every step.
    Rollers on xmin / xmax; without a plate the block floats inside the grid, with one it is
    flush with ymax under `plate ymax 0` (the plate moves with it)."""
    E = 1.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 8, "-basis", basis)
    ops.mpm("mode", mode)
    top = 8.0 if plate else 6.0
    ops.mpm("body", 1, "box", 0.0, top - 3.0, 2.0, top, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1.0e-4, 1.0e-4)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    if plate:
        ops.mpm("plate", "ymax", 0.0)
    n = len(as_list(ops.mpm("particleData", "m")))
    ops.mpm("particleVelocity", *([0.0, vy] * n))
    if mode == "implicit":
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-12, 25)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
    ev = es = ep = 0.0
    for s in range(steps):
        ops.mpm("analyze", 1, dt)
        v = as_list(ops.mpm("particleData", "v"))
        sg = as_list(ops.mpm("particleData", "sig"))
        pp = as_list(ops.mpm("particleData", "p"))
        ev = max(ev, max(abs(v[2 * k]) + abs(v[2 * k + 1] - vy) for k in range(n)) / abs(vy))
        es = max(es, max(abs(x) for x in sg) / E)
        ep = max(ep, max(abs(x) for x in pp) / E)
    return {"mode": mode, "basis": basis, "plate": plate, "steps": steps, "dt": dt, "v_rel": ev,
            "sig_rel": es, "p_rel": ep, "pass": max(ev, es, ep) <= 1e-10}


def case_bands(name):
    """The `band` lines of case block <name> in tests/mpm_twophase_cases.txt: {key: text}."""
    out, on = {}, False
    for ln in (HERE / "mpm_twophase_cases.txt").read_text(encoding="utf-8").splitlines():
        if ln.startswith("case "):
            on = ln.split()[1] == name
            continue
        if on and ln.strip() == "end":
            break
        if on and ln.startswith("band "):
            parts = ln.split(None, 2)
            out[parts[1]] = parts[2] if len(parts) > 2 else ""
    return out


def case_band_line(name, contains):
    """The text after `band ` of the first band line of case <name> that contains <contains>."""
    on = False
    for ln in (HERE / "mpm_twophase_cases.txt").read_text(encoding="utf-8").splitlines():
        if ln.startswith("case "):
            on = ln.split()[1] == name
            continue
        if on and ln.strip() == "end":
            break
        if on and ln.startswith("band ") and contains in ln:
            return ln[5:]
    raise KeyError((name, contains))


def band_value(text):
    """The first number in a band's text (`0.02`, `<= 0.3`, `>= 1.039...`, `0.01 (of ...)`)."""
    import re as _re
    m = _re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return float(m.group(0))


# ---------------------------------------------------------------------------
#  T3 Mandel (block t3_ref_mandel), 5 s ramp
# ---------------------------------------------------------------------------
def run_t3(ops, basis, stab, rec, upto=None):
    a, b, h, dt, sig0 = 10.0, 4.0, 1.0, 0.02, 10.0
    g = [x for x in rec["gated"] if (x["stab"] == "off") == (stab == "off")][0]
    p0 = g["p0"]
    steps = {str(r["tv"]): int(round(r["t"] / dt)) for r in g["rows"]}
    spec = dict(mat=(1.0e4, 0.1, 2.0), h=h, cells=(10, 4), basis=basis, box=((0.0, 0.0), (a, b)),
                fluid=(4444.0, 0.4, 1.0, 1.0e-4, stab), fix=[("xmin", "x"), ("ymin", "y")],
                pfix=[("xmax", 0.0)], plate=("ymax", a * sig0), path=((0.0, 5.0, 1.0e6), (0.0, 1.0, 1.0)))
    t0 = time.time()
    build(ops, spec)
    X0 = pdata(ops, "x", 2)
    col = [k for k, x in enumerate(X0) if x[0] < h / 2.0]
    rows = {}
    peak = -1e30
    if upto is not None:
        steps = {k: v for k, v in steps.items() if float(k) <= upto}
    last = max(steps.values())
    for step in range(1, last + 1):
        ops.mpm("analyze", 1, dt)
        pp = pdata(ops, "p")
        pc = sum(pp[k] for k in col) / len(col) / p0
        peak = max(peak, pc)
        for r, bnd in zip(g["rows"], g["band"]):
            if steps[str(r["tv"])] == step:
                rows[str(r["tv"])] = {"p_center": pc, "p_exact_quarter": r["p_exact_quarter_over_p0"],
                                      "err": abs(pc - r["p_exact_quarter_over_p0"]), "band": bnd,
                                      "fe_err": r["err_quarter_over_p0"]}
    ok = all(v["err"] <= v["band"] for v in rows.values()) and peak >= rec["peak_gate"]
    return {"basis": basis, "stab": stab, "rows": rows, "peak": peak, "peak_gate": rec["peak_gate"],
            "pass": ok, "wall_s": time.time() - t0}


# ---------------------------------------------------------------------------
#  T4 hydrostatics (block t4_ref_hydro), dt 1 s (the FE's), to t = 900 s (Tv 5)
# ---------------------------------------------------------------------------
def run_t4(ops, basis, stab, h=0.5, t_end=900.0, dt=1.0):
    H, g, rs, rf = 10.0, 9.81, 2.0, 1.0
    w = 2 * h if basis.startswith("bspline") else h
    spec = dict(mat=(1.0e4, 0.2, rs), h=h, cells=(int(round(w / h)), int(round(H / h))), basis=basis,
                box=((0.0, 0.0), (w, H)), fluid=(4444.0, 0.4, rf, 1.0e-4, stab), gravity=g,
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)])
    t0 = time.time()
    build(ops, spec)
    Y0 = [x[1] for x in pdata(ops, "x", 2)]
    n = int(round(t_end / dt))
    ops.mpm("analyze", n, dt)
    y = [x[1] for x in pdata(ops, "x", 2)]
    pp = pdata(ops, "p")
    sg = pdata(ops, "sig", 3)
    perr = max(abs(pp[k] - rf * g * (H - y[k])) for k in range(len(pp))) / (rf * g * H)
    serr = max(abs(sg[k][1] - (-rs * g * (H - Y0[k]) + rf * g * (H - y[k]))) for k in range(len(pp))) / \
        ((rs - rf) * g * H)
    worst = max(range(len(pp)), key=lambda k: abs(pp[k] - rf * g * (H - y[k])))
    return {"basis": basis, "stab": stab, "h": h, "t_end": t_end, "p_err": perr, "syy_err": serr,
            "worst_p_Y": Y0[worst], "pass": perr <= 0.01 and serr <= 0.03, "wall_s": time.time() - t0}


# ---------------------------------------------------------------------------
#  T5 B4 footing (block t5_ref_b4), route E
# ---------------------------------------------------------------------------
def run_t5(ops, stab):
    N, h, q, strip = 10, 3.0, 0.1, 3.0
    spec = dict(mat=(25000.0, 0.3, 2.67), h=h, cells=(N, N), basis="gimp", box=((0.0, 0.0), (30.0, 30.0)),
                fluid=(2.2e12, 0.4, 1.0, 1.0e-7, stab), fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")],
                pfix=[("ymax", 0.0)], traction=[("strip", "ymax", (0.0, 30.0), (strip, 30.0), (0.0, -q))],
                path=((0.0, 0.1, 1.0e6), (0.0, 1.0, 1.0)), rayleigh=(0.05, 0.02))
    t0 = time.time()
    build(ops, spec)
    ops.mpm("analyze", 100, 0.01)
    P = {}
    for t in as_list(ops.mpm("info", "gridNodeTags")):
        c = ops.nodeCoord(int(t))
        i, j = int(round(c[0] / h)), int(round(c[1] / h))
        P[(i, j)] = ops.nodeDisp(int(t), 3)
    inter = [(i, j) for i in range(1, N) for j in range(1, N)]
    p = [P[k] for k in inter]
    chi = [(-1.0) ** (i + j) for (i, j) in inter]
    pn = math.sqrt(sum(x * x for x in p))
    lap = [P[(i, j)] - 0.25 * (P[(i + 1, j)] + P[(i - 1, j)] + P[(i, j + 1)] + P[(i, j - 1)]) for (i, j) in inter]
    cb_lap = math.sqrt(sum(x * x for x in lap)) / max(pn, 1e-30)
    cb_lat = abs(sum(a * b for a, b in zip(p, chi))) / max(pn * math.sqrt(len(chi)), 1e-30)
    gate = (cb_lap <= 0.3) if stab != "off" else (cb_lap >= 0.5)
    return {"stab": stab, "CB_lap": cb_lap, "CB_lat": cb_lat, "pmax": max(abs(x) for x in p), "pass": gate,
            "wall_s": time.time() - t0}


def run_t5_routeS(ops, stab):
    """T5 route S (mode F): B4 undamped (A5: mode F has no Rayleigh); the grid p of the last step from
    `mpm info gridPressure` (ids are extended-lattice ids, ghost layer 1). Gate CB_lap <= 0.3 with the
    default stab; `-stab off` reported."""
    N, h, q, strip = 10, 3.0, 0.1, 3.0
    spec = dict(mat=(25000.0, 0.3, 2.67), h=h, cells=(N, N), basis="gimp", box=((0.0, 0.0), (30.0, 30.0)),
                mode="explicit", fluid=(2.2e12, 0.4, 1.0, 1.0e-7, stab),
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
                traction=[("strip", "ymax", (0.0, 30.0), (strip, 30.0), (0.0, -q))],
                path=((0.0, 0.1, 1.0e6), (0.0, 1.0, 1.0)))
    t0 = time.time()
    build(ops, spec)
    ops.mpm("analyze", 100, 0.01)
    flat = as_list(ops.mpm("info", "gridPressure"))
    ext = N + 3
    P = {}
    for k in range(0, len(flat), 3):
        e = int(flat[k])
        P[(e % ext - 1, e // ext - 1)] = flat[k + 1]
    inter = [(i, j) for i in range(1, N) for j in range(1, N)]
    p = [P[k] for k in inter]
    chi = [(-1.0) ** (i + j) for (i, j) in inter]
    pn = math.sqrt(sum(x * x for x in p))
    lap = [P[(i, j)] - 0.25 * (P[(i + 1, j)] + P[(i - 1, j)] + P[(i, j + 1)] + P[(i, j - 1)]) for (i, j) in inter]
    cb_lap = math.sqrt(sum(x * x for x in lap)) / max(pn, 1e-30)
    cb_lat = abs(sum(a * b for a, b in zip(p, chi))) / max(pn * math.sqrt(len(chi)), 1e-30)
    return {"route": "S", "stab": stab, "CB_lap": cb_lap, "CB_lat": cb_lat, "pmax": max(abs(x) for x in p),
            "gated": stab != "off", "pass": (cb_lap <= 0.3) if stab != "off" else None,
            "wall_s": time.time() - t0}


def _dump_all(ops, tmp, n):
    """Every particle's block dump at the committed state: (pi, {node tag: S^p}) per particle.
    Each dump runs the step preparation, which re-seeds the grid nodes with the projection."""
    out = []
    for p in range(n):
        path = str(Path(tmp) / ("sc_%d.json" % p))
        ops.mpm("debug", "particleBlocks", p, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        d = json.loads(Path(path).read_text())
        out.append((d.get("pi", 0.0), {t: s for t, s in zip(d["nodes"], d["Sp"]) if s != 0.0}, d["V"]))
    return out


def _projection(d, q):
    """c^ = M^-1 P^T V q over the nodes the supports touch (P5a's consistent projection, the
    unconditioned pRef), M = sum_p V S^p S^p^T: dense Gaussian elimination, numpy-free.
    (`nodeDisp` after a dump returns the committed grid p, not the re-seeded projection.)"""
    tags = sorted({t for _, S, _ in d for t in S})
    ix = {t: i for i, t in enumerate(tags)}
    m = len(tags)
    A = [[0.0] * (m + 1) for _ in range(m)]
    for (_, S, V), qq in zip(d, q):
        for ta, sa in S.items():
            ra = A[ix[ta]]
            ra[m] += V * sa * qq
            for tb, sb in S.items():
                ra[ix[tb]] += V * sa * sb
    for c in range(m):
        piv = max(range(c, m), key=lambda r: abs(A[r][c]))
        A[c], A[piv] = A[piv], A[c]
        for r in range(c + 1, m):
            f = A[r][c] / A[c][c]
            if f != 0.0:
                Ar, Ac = A[r], A[c]
                for k in range(c, m + 1):
                    Ar[k] -= f * Ac[k]
    x = [0.0] * m
    for r in range(m - 1, -1, -1):
        x[r] = (A[r][m] - sum(A[r][k] * x[k] for k in range(r + 1, m))) / A[r][r]
    return {t: x[ix[t]] for t in tags}


def patch_fix2(ops, tmp, E=1.0e4, drive=True, selfcheck=True, cap=2000, stop=1e-13):
    """Fix 2's patch check (interface 7b, P9 amended; stopping rule and self-check ruled at
    0fc5a3f13): the bspline2 hydrostatic column, kbar = 1, dt 1 s, stepped until the step-to-step
    change of max|p_p| is <= 1e-13 gamma_w H (cap 2000 steps: reaching it is a miss); the base
    particle's p equals gamma_w (H - y) to 1e-10 gamma_w H. drive=False: the pi = 0 twin
    (-fluidBody 0 0, steady p = 0).
    Self-check, after every step (numpy-free, black box), in gamma_w H:
      bookkeeping  |(p_n - pi_n) - [(p_{n-1} - pi_{n-1}) + sum S^p_{n-1} (c_n - c^_{n-1})]|, gate
                   1e-12: q = p - pi is FLIP-carried (alpha_p = 1); S^p_{n-1}, c^_{n-1} are the step's
                   supports and projection (dumps at n-1), c_n the grid p after the step;
      pi_c = pi    |p after the dumps - p before|: the preparation adds pi(x) - pi_c, so 0 iff
                   pi_c = pi(x_p) at the committed positions (gate 1e-12);
      motion       |(p_n - pi_n) - sum S^p_n c^_n|: q's part the moved supports cannot represent
                   (reported, never gated)."""
    H, g, rs, rf, h = 10.0, 9.81, 2.0, 1.0, 0.5
    s = rf * g * H
    spec = dict(mat=(E, 0.2, rs), h=h, cells=(2, 20), basis="bspline2", box=((0.0, 0.0), (1.0, H)),
                fluid=(4444.0, 0.4, rf, 1.0, "off"), gravity=g,
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
                fluidBody=None if drive else (0.0, 0.0))
    build(ops, spec)
    n = len(pdata(ops, "p"))
    book = picc = motion = 0.0
    hist = {}
    prev = None
    if selfcheck:
        d = _dump_all(ops, tmp, n)
        pp = pdata(ops, "p")
        # q_0 = p_init - (P M^-1 P^T V pi(0))_p (A', the invariant's Init): A' is committed with the
        # first accepted step (Fable audit finding 2), so it is built here from the dumps
        cpi = _projection(d, [d[k][0] for k in range(n)])
        q0 = [pp[k] - sum(sp * cpi.get(t, 0.0) for t, sp in d[k][1].items()) for k in range(n)]
        prev = (q0, d, _projection(d, q0))
    pmax_old, steps, change = None, 0, None
    while steps < cap:
        ops.mpm("analyze", 1, 1.0)
        steps += 1
        pp = pdata(ops, "p")
        pmax = max(abs(v) for v in pp)
        change = abs(pmax - pmax_old) / s if pmax_old is not None else None
        pmax_old = pmax
        if selfcheck:
            c = {t: v[1] for t, v in nodal_pressure(ops).items()}
            q_prev, d_prev, chat_prev = prev
            d = _dump_all(ops, tmp, n)
            pp2 = pdata(ops, "p")
            q = [pp[k] - d[k][0] for k in range(n)]
            chat = _projection(d, q)
            b = max(abs(q[k] - q_prev[k] - sum(sp * (c.get(t, 0.0) - chat_prev.get(t, 0.0))
                                               for t, sp in d_prev[k][1].items())) for k in range(n)) / s
            pc = max(abs(pp2[k] - pp[k]) for k in range(n)) / s
            m = max(abs(q[k] - sum(sp * chat.get(t, 0.0) for t, sp in d[k][1].items())) for k in range(n)) / s
            book, picc, motion = max(book, b), max(picc, pc), max(motion, m)
            if steps in (1, 2, 5, 20, 100, 200, 500, 1000) or steps % 1000 == 0:
                hist[str(steps)] = {"bookkeeping": b, "pi_c_eq_pi": pc, "motion": m}
            prev = (q, d, chat)
        if change is not None and change <= stop:
            break
    if selfcheck:
        hist[str(steps)] = {"bookkeeping": b, "pi_c_eq_pi": pc, "motion": m}
    y = [x[1] for x in pdata(ops, "x", 2)]
    pp = pdata(ops, "p")
    ex = [rf * g * (H - y[k]) if drive else 0.0 for k in range(n)]
    err = [abs(pp[k] - ex[k]) / s for k in range(n)]
    base = min(range(n), key=lambda k: y[k])
    steady = change is not None and change <= stop
    # the owner's decision (d064f5b7e): the gates are the grid coefficients and the field at the
    # particles (pi included), 1e-10; the particle-carried p is recorded (P5a's motion residual,
    # gated by T4p at 0.01)
    nodes = nodal_pressure(ops)
    exc = (lambda Y: rf * g * (H - Y)) if drive else (lambda Y: 0.0)
    coef_err = max(abs(v[1] - exc(v[0][1])) for v in nodes.values() if -1e-9 <= v[0][1] <= H + 1e-9) / s
    fp = field_at_particles(ops, tmp, {t_: v[1] for t_, v in nodes.items()})
    field_err = [abs(fp[k]["field"] - ex[k]) / s for k in range(n)]
    r = {"E": E, "drive": drive, "unit": "rhoF g H", "steps": steps, "steady": steady,
         "last_change": change, "grid_coef_err": coef_err, "field_at_particles_err": max(field_err),
         "field_base_err": field_err[base], "missing_nodes": sum(q["missing"] for q in fp),
         "particle_base_err": err[base], "particle_max_err": max(err),
         "pass": steady and coef_err <= 1e-10 and max(field_err) <= 1e-10,
         "rule": "d064f5b7e: grid coefficients + field at the particles <= 1e-10; particle p recorded"}
    if selfcheck:
        r.update({"selfcheck_bookkeeping_max": book, "selfcheck_pi_c_max": picc, "motion_max": motion,
                  "selfcheck_pass": book <= 1e-12 and picc <= 1e-12, "history": hist})
    return r


#  Diagnostic (orchestrator's ruling B on 25cc9b69e; reported, never gated)
# ---------------------------------------------------------------------------
def field_at_particles(ops, tmp, nodal_p):
    """For every particle: the grid field sum_a S^p_a(x_p) p_a with the supports the next step
    uses (`mpm debug particleBlocks` at the committed state), and its S^p, V for the projection.
    nodal_p: tag -> p read BEFORE this call (the call re-seeds the nodes)."""
    n = len(as_list(ops.mpm("particleData", "m")))
    out = []
    for p in range(n):
        path = str(Path(tmp) / ("diag_%d.json" % p))
        ops.mpm("debug", "particleBlocks", p, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        d = json.loads(Path(path).read_text())
        miss = [t for k, t in enumerate(d["nodes"]) if t not in nodal_p and d["Sp"][k] != 0.0]
        # the field includes fix 2's known affine term pi (the dump's "pi"; 0 off B-spline images)
        f = sum(d["Sp"][k] * nodal_p.get(d["nodes"][k], 0.0) for k in range(d["nen"])) + d.get("pi", 0.0)
        out.append({"field": f, "nodes": d["nodes"], "Sp": d["Sp"], "V": d["V"], "missing": len(miss)})
    return out


def residual_projection(fp, r, skip=()):
    """max_a |sum_p S^p_pa V_p r_p / sum_p S^p_pa V_p| over the grid nodes (the residual as the
    grid sees it); `skip`: node tags left out (the drained face, whose p is pinned)."""
    num, den = {}, {}
    for q, rr in zip(fp, r):
        for t, s in zip(q["nodes"], q["Sp"]):
            num[t] = num.get(t, 0.0) + s * q["V"] * rr
            den[t] = den.get(t, 0.0) + s * q["V"]
    return max((abs(num[t] / den[t]) for t in num if den[t] != 0.0 and t not in skip), default=0.0)


def nodal_pressure(ops):
    out = {}
    for t in as_list(ops.mpm("info", "gridNodeTags")):
        out[int(t)] = (ops.nodeCoord(int(t)), ops.nodeDisp(int(t), 3))
    return out


def diag_t2(ops, basis, stab, tmp):
    w = 2.0 if basis.startswith("bspline") else 1.0
    L, q, dt = 10.0, 10.0, 0.02
    E, nu, Kf, n, kbar = 1.0e4, 0.2, 4444.0, 0.4, 1.0e-4
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    cv = kbar / (1.0 / Eoed + n / Kf)
    p0 = q * (Kf / n) / (Kf / n + Eoed)
    spec = dict(mat=(E, nu, 2.0), h=1.0, cells=(int(w), 10), basis=basis, box=((0.0, 0.0), (w, L)),
                fluid=(Kf, n, 1.0, kbar, stab), fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")],
                pfix=[("ymax", 0.0)], plate=("ymax", q * w))
    build(ops, spec)
    Y0 = [x[1] for x in pdata(ops, "x", 2)]
    rows = {}
    done = 0
    for Tv, st in sorted(T2_STEPS.items(), key=lambda kv: kv[1]):
        ops.mpm("analyze", st - done, dt)
        done = st
        Tva = cv * st * dt / (L * L)
        pp = pdata(ops, "p")
        nodes = nodal_pressure(ops)
        fp = field_at_particles(ops, tmp, {t: v[1] for t, v in nodes.items()})
        ex = [p0 * terzaghi_series((L - Y0[k]) / L, Tva) for k in range(len(pp))]
        r = [pp[k] - fp[k]["field"] for k in range(len(pp))]
        node_err = max(abs(v[1] - p0 * terzaghi_series((L - v[0][1]) / L, Tva)) for v in nodes.values()
                       if v[0][1] <= L + 1e-9) / p0
        rows[str(Tv)] = {
            "Tv_actual": Tva,
            "i_particle_p_err": max(abs(pp[k] - ex[k]) for k in range(len(pp))) / p0,
            "ii_field_at_particles_err": max(abs(fp[k]["field"] - ex[k]) for k in range(len(pp))) / p0,
            "iii_grid_node_err": node_err,
            "residual_amp": max(abs(x) for x in r) / p0,
            "residual_projected": residual_projection(fp, r) / p0,
            "residual_projected_free": residual_projection(
                fp, r, skip={t for t, v in nodes.items() if v[0][1] >= L - 1e-9}) / p0,
            "missing_nodes": sum(q["missing"] for q in fp)}
    return {"basis": basis, "stab": stab, "unit": "p0", "rows": rows,
            "note": "B-spline node values are coefficients, not point values: (iii) is indicative there"}


def diag_t4(ops, basis, stab, tmp, h=0.5):
    H, g, rs, rf = 10.0, 9.81, 2.0, 1.0
    w = 2 * h if basis.startswith("bspline") else h
    spec = dict(mat=(1.0e4, 0.2, rs), h=h, cells=(int(round(w / h)), int(round(H / h))), basis=basis,
                box=((0.0, 0.0), (w, H)), fluid=(4444.0, 0.4, rf, 1.0e-4, stab), gravity=g,
                fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)])
    build(ops, spec)
    ops.mpm("analyze", 900, 1.0)
    y = [x[1] for x in pdata(ops, "x", 2)]
    pp = pdata(ops, "p")
    nodes = nodal_pressure(ops)
    fp = field_at_particles(ops, tmp, {t: v[1] for t, v in nodes.items()})
    ex = [rf * g * (H - y[k]) for k in range(len(pp))]
    r = [pp[k] - fp[k]["field"] for k in range(len(pp))]
    s = rf * g * H
    return {"basis": basis, "stab": stab, "h": h, "unit": "rhoF g H",
            "i_particle_p_err": max(abs(pp[k] - ex[k]) for k in range(len(pp))) / s,
            "ii_field_at_particles_err": max(abs(fp[k]["field"] - ex[k]) for k in range(len(pp))) / s,
            "iii_grid_node_err": max(abs(v[1] - rf * g * (H - v[0][1])) for v in nodes.values()
                                     if -1e-9 <= v[0][1] <= H + 1e-9) / s,
            "residual_amp": max(abs(x) for x in r) / s,
            "residual_projected": residual_projection(fp, r) / s,
            "residual_projected_free": residual_projection(
                fp, r, skip={t for t, v in nodes.items() if v[0][1] >= H - 1e-9}) / s,
            "missing_nodes": sum(q["missing"] for q in fp),
            "note": "B-spline node values are coefficients, not point values: (iii) is indicative there"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(a.bin)
    sys.path.insert(0, a.bin)
    import opensees as ops  # noqa: E402
    import tempfile
    build_id = ops.ladrunoBuild()
    out = Path(a.out)
    res = json.loads(out.read_text()) if out.exists() else {}
    AUTO = ["auto", 0.25]
    for leg in a.legs:
        t0 = time.time()
        if leg == "patch":
            with tempfile.TemporaryDirectory() as tmp:
                r = [patch_check(ops, b, tmp) for b in ("bspline2", "gimp")]
        elif leg == "t2":
            r = [run_t2(ops, b, s) for b in ("gimp", "bspline2") for s in ("off", AUTO)]
        elif leg == "t2_reported":
            r = [run_t2(ops, "linear", s) for s in ("off", AUTO)] + [run_t2(ops, "gimp", AUTO, kin="small")]
        elif leg in ("t3_gimp", "t3_bspline2"):
            rec = json.loads((DATA / "t3_mandel_fe.json").read_text())
            r = [run_t3(ops, leg.split("_")[1], s, rec) for s in ("off", AUTO)]
        elif leg == "t4":
            r = [run_t4(ops, b, s) for b in ("gimp", "bspline2") for s in ("off", AUTO)]
        elif leg == "t4_reported":
            r = [run_t4(ops, b, s, h=1.0) for b in ("gimp", "bspline2") for s in ("off", AUTO)]
        elif leg in ("t6_bspline2", "t6_gimp"):
            fe = json.loads((DATA / "t6_fe_dt0.005.json").read_text())
            basis = leg.split("_")[1]
            stabs = {"off": "off", "auto": AUTO}
            def fst(x):
                return "off" if x == "off" else "auto"
            r = [run_t6_terzaghi(ops, basis, stabs[fst(f["stab"])], f) for f in fe["terzaghi"]]
            r += [run_t6_mandel(ops, basis, stabs[fst(f["stab"])], f, 1.0390679304096313) for f in fe["mandel"]]
        elif leg == "t5_routeS":
            r = [run_t5_routeS(ops, AUTO), run_t5_routeS(ops, "off")]
        elif leg == "t13":
            r = [run_t13(ops, m, b, pl) for m in ("implicit", "explicit") for b in ("gimp", "bspline2")
                 for pl in (False, True)]
        elif leg == "t7":
            r = [run_t7(ops, "bspline2"), run_t7(ops, "gimp")]
        elif leg == "t5":
            r = [run_t5(ops, s) for s in (AUTO, "off")]
        elif leg == "patch_fix2":                   # the gate: E 1e4 with the per-step self-check
            with tempfile.TemporaryDirectory() as tmp:
                r = [patch_fix2(ops, tmp)]
        elif leg == "patch_fix2_controls":          # reported beside it (ruling 0fc5a3f13)
            with tempfile.TemporaryDirectory() as tmp:
                r = [patch_fix2(ops, tmp, 1.0e4, drive=False, selfcheck=False),
                     patch_fix2(ops, tmp, 1.0e6, selfcheck=False),
                     patch_fix2(ops, tmp, 1.0e9)]           # rigid: self-check + motion per step
        elif leg == "diagnostic":
            with tempfile.TemporaryDirectory() as tmp:
                r = {"t2": [diag_t2(ops, b, s, tmp) for b in ("gimp", "bspline2") for s in ("off", AUTO)],
                     "t4": [diag_t4(ops, b, s, tmp) for b in ("gimp", "bspline2") for s in ("off", AUTO)],
                     "reported_never_gated": True}
        else:
            raise SystemExit("unknown leg " + leg)
        res[leg] = {"build": build_id, "date": time.strftime("%Y-%m-%d %H:%M"), "wall_s": time.time() - t0,
                    "results": r}
        out.write_text(json.dumps(res, indent=1))
        print(json.dumps({leg: r}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
