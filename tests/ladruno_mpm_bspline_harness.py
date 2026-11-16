"""WP-232i harness: the B-spline rows of tests/mpm_bspline_cases.txt through the shipped binaries.

Written from tests/mpm_bspline_interface.md, tests/mpm_interface.md, tests/mpm_implicit_interface.md
and the cases files; it never reads tests/mpm_reference.py or tests/mpm_implicit_reference.py.

* `mms(ops, row, ncell, dt)` runs one rung of a manufactured-solution row (tests/mpm_mms.py) in
  mode F (explicit), mode A quasi-static or mode A transient, setting the per-particle body force
  every step (`mpm particleBodyForce`) and the initial velocity (`mpm particleVelocity`).
* `ladder(ops, row)` runs a row's rungs at one dt (the finest rung's rule), applies the brief's
  dt validity check at the finest rung (halving dt until it holds) and returns e per rung and p.
* `g8(ops, V, steps)` is the images-against-full-domain equality (row G8).
* `shape_points(...)`, `column(...)`: G1 / G1-K4 through `mpm debug shape`, and R1's column.
"""
from __future__ import annotations

import copy
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import mpm_mms as M  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "tests" / "mpm_bspline_cases.txt"
HOST = 232
E_MMS, NU_MMS, RHO_MMS, C_TIME = 1.0, 0.3, 1.0, 2.0       # cases-file header
OMEGA = C_TIME * math.pi
PERIOD = 2.0 / C_TIME
FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")


# --------------------------------------------------------------------------- #
#  The cases file                                                              #
# --------------------------------------------------------------------------- #
def parse_rows(path=CASES):
    rows, cur = {}, None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        t = raw.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "row":
            cur = {"name": t[1]}
            continue
        if cur is None:
            continue
        if t[0] == "end":
            rows[cur["name"]] = cur
            cur = None
        elif t[0] in ("gate", "loading"):
            cur[t[0]] = raw.split(None, 1)[1]
        else:
            vals = t[1:]
            conv = []
            for v in vals:
                try:
                    conv.append(int(v))
                except ValueError:
                    try:
                        conv.append(float(v))
                    except ValueError:
                        conv.append(v)
            cur[t[0]] = conv[0] if len(conv) == 1 else conv
    return rows


def as_list(v):
    return v if isinstance(v, list) else [v]


# --------------------------------------------------------------------------- #
#  MMS runs                                                                    #
# --------------------------------------------------------------------------- #
def _pl(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def critical_dt(row, ncell):
    """h / (c_p + max |v|) at the MMS's peak velocity (interface §6.1 with the elastic tangent)."""
    lam, mu = M.lame(E_MMS, NU_MMS)
    cp = math.sqrt((lam + 2 * mu) / RHO_MMS)
    A = M.amplitude(float(row["fmax"]))
    return (1.0 / ncell) / (cp + OMEGA * A * math.sqrt(int(row["ndm"])))


def build_mms(ops, row, ncell):
    nd = int(row["ndm"])
    ops.wipe()
    ops.model("basic", "-ndm", nd, "-ndf", nd)
    ops.nDMaterial("ElasticIsotropic", 1, E_MMS, NU_MMS)
    ops.mpm("grid", HOST, "-origin", *([0.0] * nd), "-h", 1.0 / ncell, "-cells", *([ncell] * nd),
            "-basis", row["basis"])
    mode = row.get("mode", "explicit")
    if mode != "explicit":
        ops.mpm("mode", mode)
    ops.mpm("body", 1, "box", *([0.0] * nd), *([1.0] * nd), "-ppc", int(row["ppc"]), "-mat", 1,
            "-rho", RHO_MMS)
    for d in range(nd):                                    # rollers on every face
        ops.mpm("fix", FACES[2 * d], "xyz"[d])
        ops.mpm("fix", FACES[2 * d + 1], "xyz"[d])
    if mode != "quasiStatic":
        ops.mpm("flip", float(row.get("flip", 1.0)))
    X = np.array(_pl(ops.mpm("particleData", "x"))).reshape(-1, nd)
    return X


def _analysis(ops, mode, dlam):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-13, 30)
    ops.algorithm("Newton")
    if mode == "quasiStatic":
        ops.integrator("LoadControl", dlam)
        ops.analysis("Static")
    else:
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")


def mms(ops, row, ncell, dt=None, record=True):
    """One rung. Returns {e, e_quarter, e_max, crossings, N, dt, steps, wall}."""
    nd = int(row["ndm"])
    mode = row.get("mode", "explicit")
    A = M.amplitude(float(row["fmax"]))
    lam, mu = M.lame(E_MMS, NU_MMS)
    t0 = time.perf_counter()
    X = build_mms(ops, row, ncell)
    N = X.shape[0]
    h = 1.0 / ncell
    cell = np.floor(X / h)
    crossings = 0
    hist = []
    if mode == "quasiStatic":
        nsteps = int(row.get("steps", 4))
        _analysis(ops, mode, 1.0 / nsteps)
        for n in range(nsteps):
            lam1 = (n + 1) / nsteps
            b = M.body_force_static(X, lam1, A, RHO_MMS, lam, mu)
            ops.mpm("particleBodyForce", *b.ravel().tolist())
            ops.mpm("analyze", 1, 1.0 / nsteps)
            x = np.array(_pl(ops.mpm("particleData", "x"))).reshape(-1, nd)
            c1 = np.floor(x / h)
            crossings += int(np.any(c1 != cell, axis=1).sum())
            cell = c1
            hist.append((lam1, M.error(x, X, lam1, A)))
        e = hist[-1][1]
        return {"e": e, "e_quarter": None, "e_max": max(v for _, v in hist), "crossings": crossings,
                "N": N, "dt": 1.0 / nsteps, "steps": nsteps, "wall": time.perf_counter() - t0}
    ops.mpm("particleVelocity", *M.velocity0(X, A, OMEGA).ravel().tolist())
    if mode == "implicit":
        _analysis(ops, mode, None)
    nsteps = int(round(PERIOD / dt))
    t = 0.0
    q_step = int(round(0.25 * PERIOD / dt))
    e_q, e_max = None, 0.0
    for n in range(nsteps):
        tb = t if mode == "explicit" else t + dt           # mode F: t_n; mode A: t_{n+1}
        b = M.body_force(X, tb, A, OMEGA, RHO_MMS, lam, mu)
        ops.mpm("particleBodyForce", *b.ravel().tolist())
        ops.mpm("analyze", 1, dt)
        t += dt
        if record or n == nsteps - 1:
            x = np.array(_pl(ops.mpm("particleData", "x"))).reshape(-1, nd)
            c1 = np.floor(x / h)
            crossings += int(np.any(c1 != cell, axis=1).sum())
            cell = c1
            en = M.error(x, X, math.sin(OMEGA * t), A)
            e_max = max(e_max, en)
            if n + 1 == q_step:
                e_q = en
    return {"e": en, "e_quarter": e_q, "e_max": e_max, "crossings": crossings, "N": N, "dt": dt,
            "steps": nsteps, "wall": time.perf_counter() - t0}


def ladder(ops, row, rungs=None, dt=None, validity=True, log=print):
    """Run a row's ladder at one dt; returns {rungs, h, e, p, nonincr, dt, validity, runs}."""
    rungs = rungs or as_list(row["ladder"])
    mode = row.get("mode", "explicit")
    if dt is None and mode == "explicit":
        dt = float(row.get("dt_rule", 0.1)) * critical_dt(row, rungs[-1])
    elif dt is None and mode == "implicit":
        dt = float(row.get("dt_fixed", 0.005))
    while True:
        runs = []
        for nc in rungs:
            r = mms(ops, row, nc, dt)
            log(f"  {row['name']} 1/h={nc}: e={r['e']:.4e} crossings={r['crossings']} "
                f"wall={r['wall']:.1f}s")
            runs.append(r)
        val = None
        if validity and mode == "explicit" and row.get("validity") is not None:
            rh = mms(ops, row, rungs[-1], dt / 2.0, record=False)
            val = abs(rh["e"] - runs[-1]["e"]) / runs[-1]["e"]
            log(f"  {row['name']} validity at 1/h={rungs[-1]}: e(dt/2)={rh['e']:.4e} change {val:.3%}")
            if val >= float(row["validity"]):
                dt /= 2.0
                log(f"  {row['name']}: validity not met; dt halved to {dt:.4e} and the ladder re-run")
                continue
        break
    hs = [1.0 / n for n in rungs]
    es = [r["e"] for r in runs]
    p = M.order(hs, es) if len(rungs) >= 3 else None
    nonincr = all(es[i + 1] <= es[i] for i in range(len(es) - 1))
    return {"row": row["name"], "rungs": rungs, "h": hs, "e": es, "p": p, "nonincr": nonincr,
            "dt": dt, "validity_change": val,
            "e_quarter": [r["e_quarter"] for r in runs], "e_max": [r["e_max"] for r in runs],
            "crossings": [r["crossings"] for r in runs], "N": [r["N"] for r in runs],
            "wall": [r["wall"] for r in runs]}


def verdict(row, res):
    """(ok, text) of a gated ladder row against its parsed gate."""
    ok = True
    parts = []
    if row.get("p_ge") is not None:
        ok = ok and res["p"] is not None and res["p"] >= float(row["p_ge"])
        parts.append(f"p={res['p']:.3f} (>= {row['p_ge']})")
    if int(row.get("nonincr", 0)):
        ok = ok and res["nonincr"]
        parts.append(f"non-increasing={res['nonincr']}")
    return ok, ", ".join(parts)


# --------------------------------------------------------------------------- #
#  G8 — images against the full domain                                        #
# --------------------------------------------------------------------------- #
def g8(ops, V, steps, h=0.125, ppc=2):
    """Returns (rel dx, rel dsigma, max |F_ii - 1|). Plane x = 0 (cases-file row G8)."""
    E, nu, rho = 1.0e3, 0.3, 1.0
    lam, mu = M.lame(E, nu)
    dt = 0.1 * h / math.sqrt((lam + 2 * mu) / rho)
    a, y0, y1 = 0.5, 0.25, 0.75
    hs = h / ppc
    pts = []
    for iy in range(int(round((y1 - y0) / hs))):
        for ix in range(int(round(a / hs))):
            pts.append((((ix + 0.5) * h) / ppc, y0 + ((iy + 0.5) * h) / ppc))

    def run(half):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("ElasticIsotropic", 1, E, nu)
        n = int(round((1.0 if half else 2.0) / h))
        ops.mpm("grid", HOST, "-origin", 0.0 if half else -1.0, 0.0, "-h", h, "-cells", n,
                int(round(1.0 / h)), "-basis", "bspline2")
        if half:
            ops.mpm("fix", "xmin", "x")
        pid = 0
        for xr, y in pts:
            for side in ((1,) if half else (0, 1)):
                x = xr if side else -xr
                v = (V * math.sin(math.pi * x / a) * (1.0 + y),
                     0.5 * V * math.cos(math.pi * x / (2 * a)) * math.sin(math.pi * y))
                ops.mpm("particle", pid, 1, 1, rho, hs * hs, "-x", x, y, "-lp", hs / 2, hs / 2,
                        "-v", *v)
                pid += 1
        for _ in range(steps):
            ops.mpm("analyze", 1, dt)
        return (np.array(_pl(ops.mpm("particleData", "x"))).reshape(-1, 2),
                np.array(_pl(ops.mpm("particleData", "sig"))).reshape(-1, 3),
                np.array(_pl(ops.mpm("particleData", "F"))).reshape(-1, 4))

    xf, sf, Ff = run(False)
    xh, sh, _ = run(True)
    q = np.arange(len(pts)) * 2 + 1               # the x > 0 particle of each pair in the full run
    rx = np.max(np.abs(xh - xf[q])) / np.max(np.abs(xf[q]))
    rs = np.max(np.abs(sh - sf[q])) / np.max(np.abs(sf[q]))
    strain = float(np.max(np.abs(Ff[:, [0, 3]] - 1.0)))
    return float(rx), float(rs), strain


# --------------------------------------------------------------------------- #
#  G1 / G1-K4 through `mpm debug shape`                                       #
# --------------------------------------------------------------------------- #
def node_coord(nid, x0, h, n, nd):
    ex = n[0] + 3
    i = nid % ex - 1
    rest = nid // ex
    if nd == 2:
        return np.array([x0[0] + i * h, x0[1] + (rest - 1) * h])
    ey = n[1] + 3
    j = rest % ey - 1
    k = rest // ey - 1
    return np.array([x0[0] + i * h, x0[1] + j * h, x0[2] + k * h])


def random_positions(rng, x0, h, n, nd, n_uniform, n_band):
    lo = np.array(x0[:nd])
    hi = lo + np.array(n[:nd]) * h
    out = [lo + rng.uniform(0, 1, nd) * (hi - lo) for _ in range(n_uniform)]
    for f in range(2 * nd):
        d = f // 2
        for k in range(n_band):
            x = lo + rng.uniform(0, 1, nd) * (hi - lo)
            s = 0.0 if k < 8 else rng.uniform(0, 2 * h)
            x[d] = lo[d] + s if f % 2 == 0 else hi[d] - s
            out.append(x)
    return out


def shape(ops, x, nd, folded=False):
    args = ["debug", "shape", *[float(v) for v in x], "-lp", *([0.1] * nd)]
    if folded:
        args.append("-folded")
    return _pl(ops.mpm(*args))


# --------------------------------------------------------------------------- #
#  R1: the slice-1 self-weight column per basis and resolution                #
# --------------------------------------------------------------------------- #
def column(ops, case_name, basis, refine):
    """The slice-1 column (tests/mpm_cases.txt) with `basis`, h / refine, width 1.0 (two cells at
    the base resolution, so a B-spline grid has >= 2 cells across), dt / refine, steps x refine.
    Returns {selfweight_error, cell_crossing_meter, wall}."""
    import ladruno_mpm_harness as H1
    c = copy.deepcopy(H1.case_by_name(case_name))
    c["name"] = f"{case_name}_{basis}_r{refine}"
    c["basis"] = basis
    c["h"] = c["h"] / refine
    width = 1.0
    c["cells"] = [int(round(width / c["h"])), c["cells"][1] * refine]
    for b in c["bodies"]:
        b[5] = repr(width)                      # body <id> box x0 y0 x1 y1 ...
    c["dt"] = c["dt"] / refine
    c["steps"] = c["steps"] * refine
    sw = next(e for e in c["expects"] if e["obs"] == "selfweight_error")
    if not any(e["obs"] == "cell_crossing_meter" for e in c["expects"]):
        e2 = dict(sw)
        e2["obs"] = "cell_crossing_meter"
        c["expects"].append(e2)
    t0 = time.perf_counter()
    r = H1.run_python(c, ops)
    obs = H1.observables(c, r)
    return {"wall": time.perf_counter() - t0, "accepted": r.get("stopped") is None,
            "selfweight_error": obs[("selfweight_error", *sw["args"])][0],
            "cell_crossing_meter": obs[("cell_crossing_meter", *sw["args"])][0]}
