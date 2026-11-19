"""WP-232l route-E cases (T2-T5, the A4 patch checks), from the reference blocks of
tests/mpm_twophase_cases.txt and 232l-V's records under tests/data/mpm/twophase/.

Never reads tests/mpm_twophase_reference.py. The closed forms used here are the ones the blocks
name: Terzaghi's series (the formulas of tests/test_ladruno_up_element_analytic.py:89-102), the
hydrostatic column (block t4_ref_hydro), and for Mandel the record's own p_exact at the particles.

Standalone:
    python tests/ladruno_mpm_twophase_cases_harness.py --bin <dist\\bin> <leg> [<leg> ...]
legs: patch, t2, t2_reported, t3_gimp, t3_bspline2, t4, t4_reported, t5
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
        ops.mpm("twoPhase", 1, "-Kf", Kf, "-poro", n, "-rhoF", rf, "-perm", k, k, *sa)
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


def run_t2(ops, basis, stab, kin="hypo"):
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
    last = max(T2_STEPS.values())
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
    gated = {"0.1", "0.2", "0.5"}
    ok = all(rows[k]["p_err_max"] <= 0.02 and rows[k]["U_err"] <= 0.02 for k in gated)
    return {"basis": basis, "stab": stab, "kin": kin, "rows": rows, "pass_bands_0.02": ok,
            "wall_s": time.time() - t0}


# ---------------------------------------------------------------------------
#  T3 Mandel (block t3_ref_mandel), 5 s ramp
# ---------------------------------------------------------------------------
def run_t3(ops, basis, stab, rec):
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("legs", nargs="+")
    a = ap.parse_args(argv)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(a.bin)
    sys.path.insert(0, a.bin)
    import opensees as ops  # noqa: E402
    import tempfile
    build_id = ops.ladrunoBuild()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
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
        elif leg == "t5":
            r = [run_t5(ops, s) for s in (AUTO, "off")]
        else:
            raise SystemExit("unknown leg " + leg)
        res[leg] = {"build": build_id, "date": time.strftime("%Y-%m-%d %H:%M"), "wall_s": time.time() - t0,
                    "results": r}
        OUT.write_text(json.dumps(res, indent=1))
        print(json.dumps({leg: r}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
