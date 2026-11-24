"""WP-253: the MPM benchmark matrix (pre-registration `Ladruno_implementation/253_mpm_performance.md`
section 1-4). Builds on the M0b harness (`mpm_m0b/common.py`): every run is one fresh child process
(`bench.py --child cfg.json`), thread environment pinned before import, affinity mask, profiler HDF5
read by the parent.

    python3.12 Ladruno_scripts/mpm_bench/bench.py --list                         # the whole matrix
    python3.12 Ladruno_scripts/mpm_bench/bench.py --dry-run --modes F --N 10000  # a subset, no run
    python3.12 Ladruno_scripts/mpm_bench/bench.py --smoke                        # N = 10^3, 1 repeat
    python3.12 Ladruno_scripts/mpm_bench/bench.py --modes F --dims 3 --bases gimp --mats elastic \
        --N 100000 --threads 1,8 --out results/my_run

Per row (mode, dim, basis, material, N, threads): `repeats` coarse runs (the throughput T and the
phase split) and, in mode F, ONE deep run (the material share, M2's method); the two are separate
processes. Mode F: 5 warm-up + 20 timed steps in one `mpm analyze` under one profiler run; mode A:
2 warm-up + 6 timed steps, one profiler run per step (M1).
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import traceback
from array import array
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchlib as B       # noqa: E402
import models as Mo        # noqa: E402

C = B.C
NS = (10_000, 100_000)
NS_ESMERALDA = (10_000, 100_000, 1_000_000)
THREADS_NMB = (1, 4, 8)
THREADS_ESMERALDA = (1, 4, 8, 24)
THREADS_A = (1, 8)
F_WARMUP, F_TIMED = 5, 20
A_WARMUP, A_TIMED = 2, 6
SMOKE_F_WARMUP, SMOKE_F_TIMED = 2, 5
SMOKE_A_WARMUP, SMOKE_A_TIMED = 1, 2
SPREAD_RERUN = 0.10
RERUN_REPEATS = 5

# the named children of the step (pre-registration section 2); absent in a build without the
# instrumentation commit, in which case only `mpm.advance` is reported
F_SCOPES = ("mpm.f.dtcrit", "mpm.f.admissible", "mpm.f.support", "mpm.f.p2g", "mpm.f.grid",
            "mpm.f.g2p", "mpm.f.project2", "mpm.f.update", "mpm.f.commit")
DIGEST_FIELDS = ("x", "v", "V", "F", "R", "eps", "epsHat", "sigHat", "sig")


# --------------------------------------------------------------------------- #
#  the matrix (pre-registration section 1)                                    #
# --------------------------------------------------------------------------- #
def matrix(site="nmb"):
    """Every row of section 1 for `site`: nmb = N 10^4, 10^5 at threads 1, 4, 8 plus mode A;
    esmeralda = mode F at N 10^4, 10^5, 10^6 and threads 1, 4, 8, 24 (SANISAND not at 10^6)."""
    out = []
    ns = NS_ESMERALDA if site == "esmeralda" else NS
    ths = THREADS_ESMERALDA if site == "esmeralda" else THREADS_NMB
    for N in ns:
        for dim in (3, 2):
            for basis in Mo.BASES:
                for mat in Mo.MATS:
                    if mat == "sanisand" and (basis == "gimp" or N >= 1_000_000):
                        continue            # 232k K0 (bspline2 only); 10^6 one step = minutes
                    for th in ths:
                        out.append(dict(mode="F", dim=dim, basis=basis, mat=mat, N=N, threads=th))
    if site == "nmb":
        for dim in (3, 2):
            for basis in Mo.BASES:
                for mat in ("elastic", "dp"):
                    for th in THREADS_A:
                        out.append(dict(mode="A", dim=dim, basis=basis, mat=mat, N=10_000, threads=th))
        for mat in ("elastic", "dp"):       # 10^5 only for 3D gimp elastic and DP
            for th in THREADS_A:
                out.append(dict(mode="A", dim=3, basis="gimp", mat=mat, N=100_000, threads=th))
        out.append(dict(mode="A", dim=3, basis="gimp", mat="sanisand", N=10_000, threads=1))
    return out


def smoke_matrix(Ns=(1_000,)):
    """A covering handful at small N: F 3D / 2D, gimp / bspline2, elastic / DP / SANISAND
    (bspline2), threads 1 and 2 on one cell; A 3D and 2D elastic, 3D DP."""
    out = []
    for N in Ns:
        for dim in (3, 2):
            for basis in Mo.BASES:
                for mat in ("elastic", "dp"):
                    out.append(dict(mode="F", dim=dim, basis=basis, mat=mat, N=N, threads=1))
            out.append(dict(mode="F", dim=dim, basis="bspline2", mat="sanisand", N=N, threads=1))
        out.append(dict(mode="F", dim=3, basis="gimp", mat="elastic", N=N, threads=2))
        out.append(dict(mode="A", dim=3, basis="gimp", mat="elastic", N=N, threads=1))
        out.append(dict(mode="A", dim=2, basis="bspline2", mat="elastic", N=N, threads=1))
        out.append(dict(mode="A", dim=3, basis="gimp", mat="dp", N=N, threads=1))
    return out


def add_profiles(rows, deep=True):
    """Coarse run configs plus, in mode F, one deep config per row."""
    out = []
    for r in rows:
        out.append(dict(r, profile="coarse"))
        if deep and r["mode"] == "F":
            out.append(dict(r, profile="deep"))
    return out


def steps_for(c, smoke, args):
    if c["mode"] == "F":
        w, t = (SMOKE_F_WARMUP, SMOKE_F_TIMED) if smoke else (F_WARMUP, F_TIMED)
    else:
        w, t = (SMOKE_A_WARMUP, SMOKE_A_TIMED) if smoke else (A_WARMUP, A_TIMED)
    return (args.warmup if args.warmup is not None else w,
            args.timed if args.timed is not None else t)


def tag_of(c, rep):
    return (f"b_{c['mode']}{c['dim']}d_{c['basis']}_{c['mat']}_N{c['N']}_t{c['threads']}"
            f"_{c['profile']}_r{rep}")


def select(cfgs, args):
    out = []
    for c in C.select(cfgs, args):
        if args.dims and c["dim"] not in args.dims:
            continue
        if args.bases and c["basis"] not in args.bases:
            continue
        if args.threads and c["threads"] not in args.threads:
            continue
        out.append(c)
    return out


# --------------------------------------------------------------------------- #
#  child                                                                      #
# --------------------------------------------------------------------------- #
def _digest(ops):
    """sha256 of every particle array the engine exposes, per field (exact float64 bytes)."""
    out = {}
    for f in DIGEST_FIELDS:
        try:
            v = ops.mpm("particleData", f)
        except Exception as e:                       # noqa: BLE001 -- a field the build lacks
            out[f] = f"unavailable: {type(e).__name__}"
            continue
        vals = [float(x) for x in (v if isinstance(v, (list, tuple)) else [v])]
        out[f] = hashlib.sha256(array("d", vals).tobytes()).hexdigest()[:16]
    return out


def child(cfg):
    res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")}}
    try:
        ops = C.import_ops(cfg["pyd"], cfg.get("_affinity"))
        ops.ladrunoThreads(int(cfg["threads"]))             # the WP-107 element-loop knob
        res["conditions"] = C.child_conditions(ops)          # ladrunoBuild, env pins, readback
        if cfg["mode"] == "probe":
            res["ok"] = True
        elif cfg["mode"] == "F":
            child_F(ops, cfg, res)
        else:
            child_A(ops, cfg, res)
    except BaseException as e:                               # noqa: BLE001 -- data, not a crash
        res["error"] = repr(e)
        res["traceback"] = traceback.format_exc()[-3000:]
    try:
        res["mkl"] = B.mkl_info(Path(cfg["pyd"]).parent)
    except Exception as e:                                   # noqa: BLE001
        res["mkl"] = {"error": repr(e)}
    C.write_child_result(cfg, res)
    return 0


def child_F(ops, cfg, res):
    cmds, params = Mo.column_commands(cfg["dim"], cfg["N"], cfg["mat"], cfg["basis"])
    res["model"] = params
    res["commands"] = [[c, [str(a) for a in args]] for c, args in cmds]
    mem = {"start": C.process_memory()}
    for name, args in cmds:
        getattr(ops, name)(*args)
    mem["built"] = C.process_memory()
    res["numParticles"] = C.info(ops, "numParticles")
    dtc = B.first(ops.mpm("criticalTimeStep"))
    dt = Mo.DT_FACTOR * dtc
    res["dt"], res["dt_crit"] = dt, dtc
    t0 = time.perf_counter()
    rc = ops.mpm("analyze", cfg["warmup"], dt)
    res["warmup_py_wall_s"] = time.perf_counter() - t0
    res["warmup_rc"] = rc
    if cfg["mat"] == "sanisand":            # stage 0 (elastic) for the warm-up, then the plastic stage
        res["stage_flip_particles"] = ops.mpm("stage", Mo.MAT_TAG, 1)
        dtc2 = B.first(ops.mpm("criticalTimeStep"))
        res["dt_crit_after_flip"] = dtc2
        if Mo.DT_FACTOR * dtc2 < 0.999 * dt:        # never step above the build's own 0.5 dt_crit limit
            dt = Mo.DT_FACTOR * dtc2
            res["dt_note"] = "dt lowered after the stage flip"
        res["dt"] = dt
    flags = ["-perStep"] + (["-deep"] if cfg["profile"] == "deep" else [])
    prof = C.StepProfiler(ops, cfg["_h5"], flags)
    rc, wall, run = prof(lambda: ops.mpm("analyze", cfg["timed"], dt))
    res["timed"] = {"run": run, "py_wall_s": wall, "steps": cfg["timed"], "rc": rc}
    res["memory"] = {**mem, "end": C.process_memory()}
    if cfg.get("digest"):
        res["digest"] = _digest(ops)
    res["ok"] = (res["warmup_rc"] in (0, None)) and (rc in (0, None))
    if not res["ok"]:
        res["reason"] = f"analyze returned warmup={res['warmup_rc']} timed={rc}"


def child_A(ops, cfg, res):
    cmds, params = Mo.block_commands(cfg["dim"], cfg["N"], cfg["mat"], cfg["basis"])
    res["model"] = params
    res["commands"] = [[c, [str(a) for a in args]] for c, args in cmds]
    mem = {"start": C.process_memory()}
    for name, args in cmds:
        getattr(ops, name)(*args)
    mem["built"] = C.process_memory()
    res["numParticles"] = C.info(ops, "numParticles")
    n_steps = cfg["warmup"] + cfg["timed"]
    dt = params["dt"]
    prof = C.StepProfiler(ops, cfg["_h5"])
    steps = []
    for k in range(n_steps):
        if cfg["mat"] == "sanisand" and k == cfg["warmup"]:
            res["stage_flip_particles"] = ops.mpm("stage", Mo.MAT_TAG, 1)
        before = {w: C.info(ops, w) for w in ("connectivityChanges", "retries", "acceptedSteps")}
        err, rc = None, None
        try:
            rc, wall, run = prof(lambda: ops.mpm("analyze", 1, dt))
            if rc not in (0, None):
                err = f"analyze returned {rc}"
        except Exception as e:                               # noqa: BLE001 -- a divergence is data
            err, wall, run = repr(e), None, f"s{prof.k - 1:04d}"
        after = {w: C.info(ops, w) for w in ("connectivityChanges", "retries", "acceptedSteps")}
        steps.append({"k": k, "warmup": k < cfg["warmup"], "run": run, "py_wall_s": wall,
                      "testIter": int(ops.testIter()) if err is None else None,
                      "connectivity_changed": after["connectivityChanges"] - before["connectivityChanges"],
                      "retries": after["retries"] - before["retries"],
                      "accepted": after["acceptedSteps"] - before["acceptedSteps"],
                      "error": err})
        if err is not None:
            break
    res["steps"] = steps
    res["memory"] = {**mem, "end": C.process_memory()}
    if cfg.get("digest"):
        res["digest"] = _digest(ops)
    res["ok"] = all(s["error"] is None for s in steps) and len(steps) == n_steps
    if not res["ok"]:
        res["reason"] = "; ".join(s["error"] for s in steps if s["error"]) or "short"


# --------------------------------------------------------------------------- #
#  parent: one run -> one row                                                 #
# --------------------------------------------------------------------------- #
def _find(tree, name, out):
    for k, c in (tree or {}).get("children", {}).items():
        if k == name:
            out.append(c)
        _find(c, name, out)
    return out


def phase_split(step):
    """Shares of the `step` wall per named scope; `other` is the rest. Scopes of the
    instrumentation commit may be absent: then only mpm.advance and commit are reported."""
    W = step["wall_ns"]
    ns = {s: sum(n["wall_ns"] for n in _find(step, s, [])) for s in F_SCOPES}
    present = any(v > 0 for v in ns.values())
    adv = C.wall(C.sub(step, "mpm.advance"))
    com = C.wall(C.sub(step, "commit"))
    out = {"mpm.advance": adv, "commit": com}
    if present:
        out.update({s: ns[s] for s in F_SCOPES if s != "mpm.f.commit" or ns[s] > 0})
        named = sum(ns.values()) + com
    else:
        named = adv + com
    out["other"] = W - named
    return {"scopes_present": present, "ns": out, "share": {k: v / W for k, v in out.items()} if W else {}}


def material_wall(step):
    """Σ elem_by_type wall of every `mpm.f.material` scope in the step (deep run), and its calls."""
    nodes = _find(step, "mpm.f.material", [])
    w = sum(r["wall_ns"] for n in nodes for r in (n.get("elem") or []))
    c = sum(r["count"] for n in nodes for r in (n.get("elem") or []))
    return (w, c) if nodes else (None, None)


def make_row(c, r, prof, rep, builds, machine, steps_cfg):
    """Flatten one run (child result + profile) into a results row."""
    cond = r.get("conditions") or {}
    mkl = r.get("mkl") or {}
    row = {"mode": c["mode"], "dim": c["dim"], "basis": c["basis"], "mat": c["mat"],
           "N_nominal": c["N"], "N": r.get("numParticles"), "threads": c["threads"],
           "profile": c["profile"], "repeat": rep, "tag": r.get("tag"),
           "ok": bool(r.get("ok")), "reason": r.get("reason") or r.get("error"),
           "warmup": steps_cfg[0], "timed": steps_cfg[1],
           "ladrunoBuild_exe": builds["exe"], "ladrunoBuild_pyd": cond.get("ladrunoBuild"),
           "ladrunoThreads": cond.get("ladrunoThreads"),
           "OMP_NUM_THREADS": cond.get("OMP_NUM_THREADS"), "MKL_NUM_THREADS": cond.get("MKL_NUM_THREADS"),
           "mkl": B.mkl_summary(mkl), "mkl_detail": mkl, "machine": machine,
           "child_wall_s": r.get("child_wall_s"), "child_rc": r.get("child_rc"),
           "memory_end": (r.get("memory") or {}).get("end"), "digest": r.get("digest"),
           "affinity": (cond.get("affinity") or {}).get("process_mask")}
    if cond.get("ladrunoBuild") != builds["exe"]:
        row["ok"] = False
        row["reason"] = f"build mismatch: exe {builds['exe']} pyd {cond.get('ladrunoBuild')}"
        return row
    if not r.get("ok") or not prof:
        return row
    N = row["N"]
    if c["mode"] == "F":
        p = prof.get((r.get("timed") or {}).get("run", ""))
        step = C.sub(p["tree"], "step") if p and p["tree"] else None
        if step is None or not step["calls"]:
            row["ok"], row["reason"] = False, "no `step` scope in the profile"
            return row
        n, W = step["calls"], step["wall_ns"]
        row.update({"n_timed": n, "W_ns": W, "T": N * n / (W * 1e-9),
                    "ns_per_particle_step": W / (N * n), "dt": r.get("dt"), "dt_crit": r.get("dt_crit")})
        ser = ((p.get("series") or {}).get("wall_ms")) or []
        if ser:
            ser_s = sorted(ser)
            med = ser_s[len(ser_s) // 2] if len(ser_s) % 2 else 0.5 * (ser_s[len(ser_s) // 2 - 1] + ser_s[len(ser_s) // 2])
            row["T_median_step"] = N / (med * 1e-3)
            row["step_max_over_median"] = max(ser) / med if med else None
        row["phases"] = phase_split(step)
        mw, mc = material_wall(step)
        if mw is not None:
            row["material_share"] = mw / W
            row["material_calls_per_particle_step"] = mc / (N * n)
        elif c["profile"] == "deep":
            row["material_share"] = None
    else:
        steps = [s for s in r.get("steps", []) if not s["warmup"]]
        W, ks, dec, nsteps = 0, 0, [], 0
        for s in steps:
            p = prof.get(s["run"])
            d = C.m1_step(p["tree"]) if p and p["tree"] else {"ok": False}
            if not d.get("ok"):
                continue
            d["connectivity_changed"] = bool(s.get("connectivity_changed"))
            dec.append(d)
            W += d["step_wall_ns"]
            ks += s["testIter"] or 0
            nsteps += 1
        if not nsteps:
            row["ok"], row["reason"] = False, "no mpm.a.step scope in the profile"
            return row
        row.update({"n_timed": nsteps, "W_ns": W, "T": N * nsteps / (W * 1e-9),
                    "ns_per_particle_step": W / (N * nsteps), "newton_k_total": ks,
                    "newton_k_mean": ks / nsteps, "T_per_newton_iter": N * ks / (W * 1e-9) if ks else None,
                    "retries": sum(s["retries"] for s in steps)})
        summ = C.m1_run_summary(dec)
        if summ:
            row["a_decomposition"] = {k: summ[k] for k in ("t_rb_ns", "t_it_ns", "R", "f", "k",
                                                            "step_wall_ns", "largest_rebuild_phase",
                                                            "formTangent_share_of_t_it") if k in summ}
            row["a_decomposition"]["rebuild_phases_ns"] = summ.get("rebuild_phases_ns")
    return row


# --------------------------------------------------------------------------- #
#  parent: the run loop (also used by regress.py)                             #
# --------------------------------------------------------------------------- #
def preflight(args, od):
    """ladrunoBuild of OpenSees.exe and of the pyd (one probe child); refuse on a mismatch."""
    exe, err = B.exe_build(args.exe)
    if exe is None:
        raise SystemExit(f"cannot read ladrunoBuild from {args.exe}: {err}")
    r, _ = C.run_child(__file__, dict(mode="probe", pyd=args.pyd, threads=1), od, "probe", 300,
                       args.affinity_mask)
    pyd = ((r.get("conditions") or {}).get("ladrunoBuild"))
    if not r.get("ok") or pyd is None:
        raise SystemExit(f"cannot import {args.pyd}: {r.get('error') or r.get('reason')}")
    if pyd != exe:
        raise SystemExit(f"REFUSED: ladrunoBuild differs, OpenSees.exe {exe} vs opensees.pyd {pyd} "
                         f"(a partial rebuild? see AGENTS.md rule 2)")
    return {"exe": exe, "pyd": pyd, "exe_path": str(Path(args.exe).resolve()),
            "exe_mtime": Path(args.exe).stat().st_mtime, "pyd_mtime": Path(args.pyd).stat().st_mtime,
            "mkl_probe": B.mkl_summary(r.get("mkl"))}


def run_configs(cfgs, args, od, default_reps, doc, label="bench", smoke=False):
    """Run `cfgs` (coarse and deep configs, in order). Appends rows to doc['rows'] and rewrites
    the aggregate after every configuration. Returns the new rows."""
    builds = doc["builds"]
    machine = doc["machine"]
    new_rows = []
    # one deep run per row key, not per repeat
    for c in cfgs:
        reps = 1 if c["profile"] == "deep" else C.repeats_for(c, args, default_reps)
        warm, timed = steps_for(c, smoke, args)
        for rep in range(reps):
            cfg = dict(c, pyd=args.pyd, warmup=warm, timed=timed, digest=bool(args.digest))
            tag = tag_of(c, rep)
            print(f"[{label}] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout, args.affinity_mask)
            r["tag"] = tag
            prof, perr = C.safe_read_profile(str(od / f"{tag}.h5"))
            if perr:
                r["profile_error"] = perr
            row = make_row(c, r, prof, rep, builds, machine, (warm, timed))
            if perr and row["ok"]:
                row["ok"], row["reason"] = False, f"profile unreadable: {perr}"
            new_rows.append(row)
            doc["rows"].append(row)
            print(f"[{label}]   ok={row['ok']} N={row.get('N')} T={row.get('T', 0):.4g} p-steps/s "
                  f"{row.get('reason') or ''} ({row.get('child_wall_s') or 0:.1f} s)", flush=True)
        C.write_json(od / f"{label}_results.json", doc)
    return new_rows


def new_doc(args, site, od):
    builds = preflight(args, od)
    return {"script": "bench", "site": site, "conditions": C.parent_conditions(args),
            "machine": B.machine_info(), "builds": builds,
            "protocol": {"repeats": args.repeats or 3, "rerun_spread": SPREAD_RERUN,
                         "rerun_repeats": RERUN_REPEATS,
                         "F": f"{F_WARMUP} warm-up + {F_TIMED} timed steps, one mpm analyze, one profiler run",
                         "A": f"{A_WARMUP} warm-up + {A_TIMED} timed steps, one profiler run per step",
                         "T": "N * n_timed / sum of the profiler's step wall (mpm.a.step in mode A)"},
            "rows": []}


def build_parser():
    ap = C.base_parser(__doc__)
    ap.set_defaults(pyd=str(B.DEFAULT_PYD))
    ap.add_argument("--exe", default=str(B.DEFAULT_EXE), help="the OpenSees.exe whose ladrunoBuild "
                    "must equal the pyd's (a mismatch refuses the whole invocation)")
    ap.add_argument("--site", choices=["nmb", "esmeralda"], default="nmb",
                    help="which matrix: nmb (N 10^4, 10^5; threads 1, 4, 8; mode A) or esmeralda "
                         "(mode F, N up to 10^6, threads 1, 4, 8, 24)")
    ap.add_argument("--bases", type=lambda s: [x.strip() for x in s.split(",") if x.strip()], default=None)
    ap.add_argument("--dims", type=lambda s: [int(x) for x in s.split(",") if x.strip()], default=None)
    ap.add_argument("--threads", type=lambda s: [int(x) for x in s.split(",") if x.strip()], default=None)
    ap.add_argument("--no-deep", action="store_true", help="skip the mode-F deep (material share) run")
    ap.add_argument("--digest", action="store_true", help="record a sha256 of every particle array "
                    "after the timed steps (the determinism check across thread counts)")
    ap.add_argument("--dry-run", action="store_true", help="print the selected rows with their "
                    "model sizes and exit")
    ap.add_argument("--rerun-from", default=None, help="a results JSON: re-run only its rows whose "
                    "spread exceeds 10 %% with 5 repeats (pre-registration section 4)")
    ap.add_argument("--warmup", type=int, default=None)
    ap.add_argument("--timed", type=int, default=None)
    return ap


def main():
    args = C.parse_args(build_parser())
    if args.child:
        return child(C.load_cfg(args.child))
    smoke = args.smoke
    default_reps = 1 if smoke else 3
    if args.rerun_from:
        import summarize as S
        base = json.loads(Path(args.rerun_from).read_text(encoding="utf-8"))
        keys = S.flagged_keys(S.aggregate([base["rows"]]))
        rows = [dict(zip(("mode", "dim", "basis", "mat", "N", "threads"), k)) for k in keys]
        args.repeats = RERUN_REPEATS
        cfgs = [dict(r, profile="coarse") for r in rows]
        for c in cfgs:
            c["N"] = int(c["N"])
    else:
        rows = smoke_matrix(tuple(args.N) if args.N else (1_000,)) if smoke else matrix(args.site)
        cfgs = add_profiles(select(rows, args), deep=not args.no_deep)
    if args.dry_run:
        for c in cfgs:
            if c["profile"] == "coarse":
                n = Mo.nominal_particles(c["mode"], c["dim"], c["N"])
                w, t = steps_for(c, smoke, args)
                print(f"{c['mode']} {c['dim']}D {c['basis']:8s} {c['mat']:9s} N={c['N']:>8} "
                      f"(model {n:>8}, {100 * (n - c['N']) / c['N']:+.1f} %) threads={c['threads']:>2} "
                      f"{w}+{t} steps x{C.repeats_for(c, args, default_reps)}"
                      f"{' +deep' if any(d['profile'] == 'deep' and B.row_key(d) == B.row_key(c) for d in cfgs) else ''}")
        print(f"{sum(c['profile'] == 'coarse' for c in cfgs)} rows")
        return 0
    if args.list:
        return C.list_runs(cfgs, args, "bench", tag_of, default_reps)
    od = C.out_dir(args, "bench" + ("_smoke" if smoke else f"_{args.site}"))
    doc = new_doc(args, args.site, od)
    print(f"[bench] exe == pyd == {doc['builds']['exe']}; machine {doc['machine']['cpu_model']}; "
          f"mkl {doc['builds']['mkl_probe']}", flush=True)
    t0 = time.time()
    run_configs(cfgs, args, od, default_reps, doc, "bench", smoke)
    doc["wall_s"] = time.time() - t0
    C.write_json(od / "bench_results.json", doc)
    print(f"[bench] wrote {od / 'bench_results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
