"""WP-232g M4: the 232f data owed on the quasi-static cantilever.

The cases `cantilever_qs_r1..r3` of tests/mpm_implicit_cases.txt, built by the 232f harness's own
translator (tests/ladruno_mpm_implicit_harness.py: model_commands / analysis_commands -- its
Plain / RCM / UmfPack / NormDispIncr 1e-10 30 / Newton / LoadControl analysis, unchanged), under
`mpm emptyNodes renormalize` with `mpm massCutoff` 5e-2 (the default), 1e-6 and 1e-10. Each step
is `mpm analyze 1 dt` under its own coarse profiler run. Recorded per step: Newton iterations
(`testIter` and the profiler's linearSolve calls), retries (`mpm cut` re-tries are data: a
divergence and its step), the tip deflection (the cases file's tip_disp_y: mass-weighted mean
of y - Y over the particles with X >= x_from), and formTangent's share of t_it.

    python3.12 Ladruno_scripts/mpm_m0b/m4_cantilever.py [--pyd PATH] [--smoke] [--only ...]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C       # noqa: E402

CASES = ("cantilever_qs_r1", "cantilever_qs_r2", "cantilever_qs_r3")
CUTOFFS = (5.0e-2, 1.0e-6, 1.0e-10)
THREADS = 8
GATE_TOL = 0.05          # the 232f gate tolerance (relative)


def configs(smoke):
    return [dict(case=case, masscut=mc) for case in (CASES[:1] if smoke else CASES) for mc in CUTOFFS]


def child(cfg):
    res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")}}
    try:
        ops = C.import_ops(cfg["pyd"])
        res["conditions"] = C.child_conditions(ops)
        sys.path.insert(0, str(C.REPO / "tests"))
        import ladruno_mpm_implicit_harness as H
        c = H.case_by_name(cfg["case"])
        xf = H.tip_args(c)[0]
        cmds = list(H.model_commands(c))
        cmds.append(("mpm", ["massCutoff", cfg["masscut"]]))
        cmds.append(("mpm", ["emptyNodes", "renormalize"]))
        cmds += list(H.analysis_commands(c))
        res["commands"] = [[n, [str(a) for a in H._py_args(a)]] for n, a in cmds]
        for name, args in cmds:
            getattr(ops, name)(*H._py_args(args))
        nd = c["ndm"]
        X = C.particle_field(ops, "x")
        m = C.particle_field(ops, "m")
        res["numParticles"] = len(m)
        res["x_from"] = xf
        dt = c["dt"]
        prof = C.StepProfiler(ops, cfg["_h5"])
        steps = []
        for k in range(c["steps"]):
            r0 = C.info(ops, "retries")
            err = None
            try:
                _, wall, run = prof(lambda: ops.mpm("analyze", 1, dt))
            except Exception as e:
                err, wall, run = repr(e), None, f"s{prof.k - 1:04d}"
            x = C.particle_field(ops, "x")
            steps.append({"k": k, "run": run, "py_wall_s": wall, "error": err,
                          "testIter": int(ops.testIter()),
                          "retries": C.info(ops, "retries") - r0,
                          "tip": H._tip(nd, X, x, m, xf),
                          # the accepted sub-steps of this call (a cut step shows several)
                          "dtHistory": H._as_list(ops.mpm("info", "dtHistory"))})
            if err is not None:
                break
        res["steps"] = steps
        res["completed_steps"] = sum(1 for s in steps if s["error"] is None)
        res["converged_all"] = res["completed_steps"] == c["steps"]
        res["first_retry_step"] = next((s["k"] for s in steps if s["retries"]), None)
        res["tip_end"] = steps[-1]["tip"] if steps and steps[-1]["error"] is None else None
        res["ok"] = True
    except BaseException as e:
        res["error"] = repr(e)
    C.write_child_result(cfg, res)
    return 0


def tag_of(c, rep):
    return f"m4_{c['case']}_mc{c['masscut']:.0e}_r{rep}"


def main():
    ap = C.base_parser(__doc__)
    args = ap.parse_args()
    if args.child:
        return child(C.load_cfg(args.child))
    reps = args.repeats or (1 if args.smoke else 3)
    od = C.out_dir(args, "m4" + ("_smoke" if args.smoke else ""))
    doc = {"script": "m4_cantilever", "conditions": C.parent_conditions(args),
           "pinned": {"threads": THREADS, "repeats": reps, "warmup_steps": 0,
                      "note": "every step is data (the ramp's 20 steps); timing medians over repeats"},
           "runs": [], "summary": []}
    t0 = time.time()
    for c in configs(args.smoke):
        if not C.keep(c, args.only):
            continue
        per = []
        for rep in range(reps):
            cfg = dict(c, pyd=args.pyd, threads=THREADS)
            tag = tag_of(c, rep)
            print(f"[m4] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout)
            prof = C.read_profile(str(od / f"{tag}.h5"))
            for s in r.get("steps", []):
                p = prof.get(s["run"])
                s["profile"] = p["tree"] if p else None
                s["derived"] = C.m1_step(p["tree"]) if p and p["tree"] else {"ok": False}
            ok_steps = [s["derived"] for s in r.get("steps", []) if s.get("error") is None]
            r["run_summary"] = C.m1_run_summary(ok_steps)
            r["newton_iterations_per_step"] = [s["derived"].get("k") for s in r.get("steps", [])]
            r["tag"] = tag
            doc["runs"].append(r)
            per.append(r)
            print(f"[m4]   ok={r.get('ok')} steps={r.get('completed_steps')} tip={r.get('tip_end')} "
                  f"first_retry={r.get('first_retry_step')} ({r.get('child_wall_s', 0):.1f} s)", flush=True)
        rs = [r.get("run_summary") or {} for r in per]
        doc["summary"].append({
            "config": c, "repeats": len(per),
            "completed_steps": [r.get("completed_steps") for r in per],
            "first_retry_step": [r.get("first_retry_step") for r in per],
            "tip_end": C.median_spread([r.get("tip_end") for r in per]),
            "newton_iterations_per_step": per[0].get("newton_iterations_per_step") if per else None,
            "formTangent_share_of_t_it": C.median_spread([x.get("formTangent_share_of_t_it") for x in rs]),
            "t_it_ns": C.median_spread([x.get("t_it_ns") for x in rs])})
        C.write_json(od / "m4_results.json", doc)
    doc["tip_comparison"] = compare(doc["summary"])
    doc["wall_s"] = time.time() - t0
    C.write_json(od / "m4_results.json", doc)
    print(f"[m4] wrote {od / 'm4_results.json'}")
    return 0


def compare(summary):
    """Per case: the tip at 1e-6 and 1e-10 against 5e-2's, relative (the 232f gate tolerance)."""
    out = []
    for case in CASES:
        rows = {s["config"]["masscut"]: s for s in summary if s["config"]["case"] == case}
        base = (rows.get(5.0e-2) or {}).get("tip_end", {}).get("median")
        row = {"case": case, "tip_5e-2": base}
        for mc in (1.0e-6, 1.0e-10):
            s = rows.get(mc)
            tip = (s or {}).get("tip_end", {}).get("median")
            done = (s or {}).get("completed_steps") or []
            row[f"{mc:.0e}"] = {"tip": tip, "completed_steps": done,
                                "rel_diff": abs(tip - base) / abs(base) if (tip is not None and base) else None,
                                "exceeds_gate_tol": (abs(tip - base) / abs(base) > GATE_TOL)
                                if (tip is not None and base) else None}
        out.append(row)
    return out


if __name__ == "__main__":
    sys.exit(main())
