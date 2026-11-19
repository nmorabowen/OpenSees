"""Route S Mandel: the linearised amplification by central differences of TRAJECTORIES.

The setters (`particleVelocity`, `particlePressure`) act only before the first step, so the one-step
Jacobian on the full state (x, v, F/sigma', p) cannot be built column by column mid-run. Instead two
runs start from the deck's initial state +/- delta (a sin pattern on v and p, |delta| 1e-7 m/s,
1e-5 kPa); D(t) = (run+ - run-) / 2 is the linearised response, carried through the FULL state by
route S itself. Its per-step growth over the late window is the spectral radius of the step map
(the dominant mode wins over the neutral ones), and D at the end is that mode: localisation by
particle, energy split (kinetic / strain / fluid storage), grid p, and the plate (top-row v_y).

argv: basis dt t_end [twin] [masscut]   twin: 1 = single phase (no twoPhase), masscut: mode-F u cutoff
Writes amp_<tag>.json in the current directory. Diagnostic only (WP-232l T6 Mandel, ruling on 4b74bd335).
"""
import json
import math
import os
import sys
B = r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\compassionate-mccarthy-18cf89\dist\bin"
os.add_dll_directory(B)
sys.path.insert(0, B)
sys.path.insert(0, r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\agent-adf46e5c8275479c0\tests")
import opensees as ops  # noqa: E402
import ladruno_mpm_twophase_cases_harness as H  # noqa: E402

basis, dt, tend = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
twin = len(sys.argv) > 4 and sys.argv[4] == "1"
mcut = float(sys.argv[5]) if len(sys.argv) > 5 else None
a, b, h, sig0 = 10.0, 4.0, 1.0, 10.0
E, nu, rho, Kf, n = 1.0e4, 0.1, 2.0, 4444.0, 0.4
M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))


def deck():
    spec = dict(mat=(E, nu, rho), h=h, cells=(10, 4), basis=basis, box=((0.0, 0.0), (a, b)),
                mode="explicit", fix=[("xmin", "x"), ("ymin", "y")], pfix=[("xmax", 0.0)],
                plate=("ymax", a * sig0), path=((0.0, 5.0, 1.0e6), (0.0, 1.0, 1.0)))
    if not twin:
        spec["fluid"] = (Kf, n, 1.0, 1.0e-4, "off")
    H.build(ops, spec)
    if mcut is not None:
        ops.mpm("massCutoff", mcut)
        ops.mpm("pressureCutoff", max(mcut, 5e-2))


nsteps = int(round(tend / dt))
every = max(1, nsteps // 30)


def run(sign):
    deck()
    N = len(H.pdata(ops, "m"))
    dv = [sign * 1e-7 * math.sin(1.7 * k + 0.3) for k in range(2 * N)]
    ops.mpm("particleVelocity", *dv)
    if not twin:
        ops.mpm("particlePressure", *[sign * 1e-5 * math.sin(2.3 * k + 0.7) for k in range(N)])
    snaps = {}
    for s in range(1, nsteps + 1):
        try:
            ops.mpm("analyze", 1, dt)
        except Exception:
            snaps["failed_at"] = s
            break
        if s % every == 0 or s == nsteps:
            g = H.as_list(ops.mpm("info", "gridPressure")) if not twin else []
            snaps[s] = {"x": H.pdata(ops, "x"), "v": H.pdata(ops, "v"), "sig": H.pdata(ops, "sig"),
                        "p": H.pdata(ops, "p"), "g": g}
    return snaps, N


sp, N = run(+1.0)
sm, _ = run(-1.0)
m = H.pdata(ops, "m")
V = H.pdata(ops, "V")
X0 = H.pdata(ops, "x", 2)
rows = []
last = None
for s in sorted(k for k in sp if isinstance(k, int) and k in sm):
    A, Bm = sp[s], sm[s]
    D = {f: [(x - y) / 2.0 for x, y in zip(A[f], Bm[f])] for f in ("x", "v", "sig", "p")}
    ek = sum(0.5 * m[k] * (D["v"][2 * k] ** 2 + D["v"][2 * k + 1] ** 2) for k in range(N))
    es = sum(0.5 * V[k] * (D["sig"][3 * k] ** 2 + D["sig"][3 * k + 1] ** 2 + 2 * D["sig"][3 * k + 2] ** 2) / M
             for k in range(N))
    ef = sum(0.5 * V[k] * D["p"][k] ** 2 * n / Kf for k in range(N)) if not twin else 0.0
    rows.append({"step": s, "t": s * dt, "energy": ek + es + ef, "kin": ek, "strain": es, "fluid": ef})
    last = (s, D, A, Bm)
# growth per step over the last third (energy ~ amplitude^2)
tail = rows[len(rows) * 2 // 3:]
rate = (math.log(tail[-1]["energy"] / tail[0]["energy"]) / (2.0 * (tail[-1]["step"] - tail[0]["step"]))
        if len(tail) > 1 and tail[0]["energy"] > 0 else float("nan"))
s, D, A, Bm = last
amp = [math.hypot(D["v"][2 * k], D["v"][2 * k + 1]) for k in range(N)]
order = sorted(range(N), key=lambda k: -amp[k])
tot = sum(x * x for x in amp)
top = [{"x0": X0[k], "dv": [D["v"][2 * k], D["v"][2 * k + 1]], "dp": D["p"][k] if not twin else 0.0,
        "share_v2": amp[k] ** 2 / tot} for k in order[:6]]
toprow = [k for k in range(N) if X0[k][1] > b - h / 2]
plate_vy = sum(D["v"][2 * k + 1] for k in toprow) / len(toprow)
gdiff = []
if not twin and A["g"]:
    ga = {int(A["g"][i]): A["g"][i + 1] for i in range(0, len(A["g"]), 3)}
    gb = {int(Bm["g"][i]): Bm["g"][i + 1] for i in range(0, len(Bm["g"]), 3)}
    ext = 10 + 3
    gd = sorted(((abs(ga[k] - gb.get(k, 0.0)) / 2.0, k) for k in ga), reverse=True)[:5]
    gdiff = [{"node_ij": [k % ext - 1, k // ext - 1], "dp": v} for v, k in gd]
out = {"basis": basis, "dt": dt, "t_end": tend, "twin_single_phase": twin, "masscut": mcut,
       "failed_at": sp.get("failed_at") or sm.get("failed_at"),
       "growth_per_step": math.exp(rate) if rate == rate else None,
       "doubling_time_s": (math.log(2.0) / rate * dt) if rate and rate > 0 else None,
       "energy_split_end": rows[-1], "energy_history": [(r["t"], r["energy"]) for r in rows],
       "top_particles_by_dv": top, "plate_dvy_toprow_mean": plate_vy, "grid_dp_top": gdiff}
print(json.dumps({k: out[k] for k in out if k != "energy_history"}, indent=1))
tag = "%s_dt%s_%s%s" % (basis, sys.argv[2], "twin" if twin else "up", "_mc%s" % sys.argv[5] if mcut is not None else "")
json.dump(out, open(os.path.join(os.getcwd(), "amp_%s.json" % tag), "w"), indent=1)
