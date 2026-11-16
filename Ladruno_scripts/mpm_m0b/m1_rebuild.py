"""WP-232g M1 (and M1-t): the mode-A rebuild against one Newton iteration.

For every configuration (N, mode qs|transient, material elastic|dp, solver Pardiso|UmfPack,
empty-node treatment renormalize|exclude|ghost) and repeat, one fresh process builds the block
(models.py), takes `warmup` steps unprofiled-for-the-table (still profiled and kept, flagged
warm-up), then `timed` steps, each `mpm analyze 1 dt` under its own coarse profiler run. Per step:
the profiler rollup (raw), the Python wall of the call, Newton iterations (`testIter`), the
`mpm info` counters (connectivityChanges, domainChanges, retries, gridNodes), the private bytes,
the top-layer settlement. The parent decomposes each step (common.m1_step) and aggregates.

    python3.12 Ladruno_scripts/mpm_m0b/m1_rebuild.py [--pyd PATH] [--smoke] [--only N=1000,...]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C       # noqa: E402
import models as Mo      # noqa: E402

NS = (1000, 10000, 100000)
MODES = ("qs", "transient")
MATS = ("elastic", "dp")
SOLVERS = ("Pardiso", "UmfPack")
TREATMENTS = ("renormalize", "exclude", "ghost")
THREADS = 8                       # implicit runs (232g pinned conditions)
WARMUP, TIMED = 2, 6
SMOKE_WARMUP, SMOKE_TIMED = 2, 6      # the full run's steps: N = 10^3 is cheap


def configs(smoke):
    out = []
    if smoke:
        # a covering design at N = 10^3: every solver x mode x material on the default treatment,
        # and the two other treatments on (qs, elastic, Pardiso) -- every code path, 10 children.
        for mode in MODES:
            for mat in MATS:
                for solver in SOLVERS:
                    out.append(dict(N=1000, mode=mode, mat=mat, solver=solver, treatment="renormalize"))
        for tr in ("exclude", "ghost"):
            out.append(dict(N=1000, mode="qs", mat="elastic", solver="Pardiso", treatment=tr))
        return out
    for N in NS:
        for mode in MODES:
            for mat in MATS:
                for solver in SOLVERS:
                    for tr in TREATMENTS:
                        out.append(dict(N=N, mode=mode, mat=mat, solver=solver, treatment=tr))
    return out


# --------------------------------------------------------------------------- #
#  child                                                                      #
# --------------------------------------------------------------------------- #
def child(cfg):
    res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")}}
    try:
        ops = C.import_ops(cfg["pyd"], cfg.get("_affinity"))
        res["conditions"] = C.child_conditions(ops)
        n_steps = cfg["warmup"] + cfg["timed"]
        cmds, params = Mo.block_commands(cfg["N"], cfg["mat"], cfg["mode"], cfg["solver"],
                                         cfg["treatment"])
        res["model"] = params
        res["commands"] = [[c, [str(a) for a in args]] for c, args in cmds]
        mem = {"start": C.process_memory()}
        for name, args in cmds:
            getattr(ops, name)(*args)
        mem["built"] = C.process_memory()
        res["numParticles"] = C.info(ops, "numParticles")
        X = C.particle_field(ops, "x")
        top = Mo.top_layer(X, params["h"])
        prof = C.StepProfiler(ops, cfg["_h5"])
        steps = []
        dt = params["dt"]
        for k in range(n_steps):
            before = {w: C.info(ops, w) for w in ("connectivityChanges", "domainChanges", "retries",
                                                  "acceptedSteps")}
            err = None

            def one():
                ops.mpm("analyze", 1, dt)
            try:
                _, wall, run = prof(one)
            except Exception as e:                       # a divergence is data
                err = repr(e)
                wall, run = None, f"s{prof.k - 1:04d}"
            after = {w: C.info(ops, w) for w in ("connectivityChanges", "domainChanges", "retries",
                                                 "acceptedSteps", "gridNodes")}
            x = C.particle_field(ops, "x")
            settle = sum(x[3 * p + 2] - X[3 * p + 2] for p in top) / len(top) if top else None
            row = {"k": k, "warmup": k < cfg["warmup"], "run": run, "py_wall_s": wall,
                   "testIter": int(ops.testIter()) if err is None else None,
                   "connectivity_changed": after["connectivityChanges"] - before["connectivityChanges"],
                   "domainChanges": after["domainChanges"] - before["domainChanges"],
                   "retries": after["retries"] - before["retries"],
                   "accepted": after["acceptedSteps"] - before["acceptedSteps"],
                   "gridNodes": after["gridNodes"], "settlement_top": settle,
                   "memory": C.process_memory(), "error": err}
            try:
                row["nEqn"] = int(ops.systemSize())
            except Exception:
                row["nEqn"] = None
            steps.append(row)
            if err is not None:
                break
        res["steps"] = steps
        res["memory"] = mem
        if cfg["mat"] == "dp":
            sig = C.particle_field(ops, "sig")
            res["yield_fraction_end"] = Mo.yield_fraction(sig, "dp", {"rho": Mo.DP_RHO,
                                                                       "sigmaY": Mo.DP_SIGY})
        res["ok"] = all(s["error"] is None for s in steps) and len(steps) == n_steps
    except BaseException as e:                            # the parent reads whatever got here
        res["error"] = repr(e)
    C.write_child_result(cfg, res)
    return 0


# --------------------------------------------------------------------------- #
#  parent                                                                     #
# --------------------------------------------------------------------------- #
def tag_of(c, rep):
    return f"m1_N{c['N']}_{c['mode']}_{c['mat']}_{c['solver']}_{c['treatment']}_r{rep}"


MIN_TIMED_AT_1E5 = 4      # the orchestrator's round-2 floor at N = 10^5


def main():
    ap = C.base_parser(__doc__)
    ap.add_argument("--warmup", type=int, default=None)
    ap.add_argument("--timed", type=int, default=None)
    args = C.parse_args(ap)
    if args.child:
        return child(C.load_cfg(args.child))
    default_reps = 1 if args.smoke else 3
    warm = args.warmup if args.warmup is not None else (SMOKE_WARMUP if args.smoke else WARMUP)
    timed = args.timed if args.timed is not None else (SMOKE_TIMED if args.smoke else TIMED)
    cfgs = C.select(configs(args.smoke), args)
    if args.list:
        return C.list_runs(cfgs, args, "m1", tag_of, default_reps, {"warmup": warm, "timed": timed})
    od = C.out_dir(args, "m1" + ("_smoke" if args.smoke else ""))
    doc = {"script": "m1_rebuild", "conditions": C.parent_conditions(args),
           "pinned": {"threads": THREADS, "warmup_steps": warm, "timed_steps": timed,
                      "timed_steps_floor_at_1e5": MIN_TIMED_AT_1E5,
                      "repeats": args.repeats or default_reps, "repeats_by_N": args.repeats_by_N,
                      "statistic": "median of the repeats, spread = (max-min)/median"},
           "definitions": "232g brief: t_rb = prepare + domainChanged + every soe.symbolic; "
                          "t_it = (formTangent + linearSolve - symbolic + formUnbalance + update + "
                          "convTest) / k; R = t_rb / t_it; f = t_rb / step wall (mpm.a.step)",
           "runs": [], "summary": []}
    t0 = time.time()
    for c in cfgs:
        per_rep = []
        tc = max(timed, MIN_TIMED_AT_1E5) if c["N"] >= 100000 else timed
        for rep in range(C.repeats_for(c, args, default_reps)):
            cfg = dict(c, pyd=args.pyd, threads=THREADS, warmup=warm, timed=tc)
            tag = tag_of(c, rep)
            print(f"[m1] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout, args.affinity_mask)
            prof, perr = C.safe_read_profile(str(od / f"{tag}.h5"))
            if perr:
                r["profile_error"] = perr
            for s in r.get("steps", []):
                p = prof.get(s["run"])
                s["profile"] = p["tree"] if p else None
                s["derived"] = C.m1_step(p["tree"]) if p and p["tree"] else {"ok": False}
                s["derived"]["connectivity_changed"] = bool(s.get("connectivity_changed"))
            timed_steps = [s["derived"] for s in r.get("steps", []) if not s["warmup"]
                           and s.get("error") is None]
            r["run_summary"] = C.m1_run_summary(timed_steps)
            r["connectivity_changed_steps"] = [s["k"] for s in r.get("steps", [])
                                               if s.get("connectivity_changed")]
            r["tag"] = tag
            doc["runs"].append(r)
            per_rep.append(r)
            rs = r["run_summary"] or {}
            print(f"[m1]   ok={r.get('ok')} {r.get('reason', '')} R={rs.get('R')} f={rs.get('f')} "
                  f"k={rs.get('k')} ({r.get('child_wall_s', 0):.1f} s)", flush=True)
        doc["summary"].append(summarize(c, per_rep))
        C.write_json(od / "m1_results.json", doc)
    doc["fits"] = fits(doc["summary"])
    doc["wall_s"] = time.time() - t0
    C.write_json(od / "m1_results.json", doc)
    print(f"[m1] wrote {od / 'm1_results.json'}")
    return 0


def summarize(c, reps):
    rs = [r.get("run_summary") or {} for r in reps]
    out = {"config": c, "ok_runs": sum(1 for r in reps if r.get("ok")), "repeats": len(reps)}
    out["failures"] = [{"rep": i, "reason": r.get("reason"), "error": r.get("error"),
                        "child_wall_s": r.get("child_wall_s")} for i, r in enumerate(reps) if not r.get("ok")]
    out["child_wall_s"] = [r.get("child_wall_s") for r in reps]
    for q in ("t_rb_ns", "t_it_ns", "R", "k", "f", "step_wall_ns", "symbolic_share_of_t_rb",
              "formTangent_share_of_t_it", "t_rb_changed_ns", "t_rb_unchanged_ns",
              "symbolic_changed_ns", "symbolic_unchanged_ns", "n_changed", "n_unchanged"):
        out[q] = C.median_spread([x.get(q) for x in rs])
    phases = {}
    for x in rs:
        for p, v in (x.get("rebuild_phases_ns") or {}).items():
            phases.setdefault(p, []).append(v)
    out["rebuild_phases_ns"] = {p: C.median_spread(v)["median"] for p, v in phases.items()}
    if out["rebuild_phases_ns"]:
        out["largest_rebuild_phase"] = max(out["rebuild_phases_ns"],
                                           key=lambda p: out["rebuild_phases_ns"][p] or 0)
    ends = [r["steps"][-1].get("settlement_top") for r in reps if r.get("steps")]
    out["settlement_top_end"] = C.median_spread(ends)
    out["connectivity_changed_timed_steps"] = [
        sum(1 for s in r.get("steps", []) if not s["warmup"] and s.get("connectivity_changed"))
        for r in reps]
    out["retries"] = [sum(s.get("retries") or 0 for s in r.get("steps", [])) for r in reps]
    return out


def fits(summary):
    """log-log slope of R (and t_rb, t_it) against N per (mode, mat, solver, treatment)."""
    groups = {}
    for s in summary:
        c = s["config"]
        key = (c["mode"], c["mat"], c["solver"], c["treatment"])
        groups.setdefault(key, []).append(s)
    out = []
    for key, rows in groups.items():
        rows = sorted(rows, key=lambda s: s["config"]["N"])
        Ns = [s["config"]["N"] for s in rows]
        row = {"mode": key[0], "mat": key[1], "solver": key[2], "treatment": key[3], "N": Ns}
        for q in ("R", "t_rb_ns", "t_it_ns", "f"):
            row["slope_" + q] = C.loglog_slope(Ns, [s[q]["median"] for s in rows])
        out.append(row)
    return out


if __name__ == "__main__":
    sys.exit(main())
