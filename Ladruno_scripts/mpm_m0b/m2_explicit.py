"""WP-232g M2: explicit (mode F) throughput, the material update's share, and the SANISAND
projection.

column runs   the 3D quarter-column collapse (models.column_commands) at N = 10^3, 10^4, 10^5,
              ElasticIsotropic / DruckerPrager / LadrunoJ2, threads 1 (and the default for
              elastic). Each configuration runs TWICE, in separate processes: `coarse` (the
              throughput: no per-call clock reads) and `deep` (the material share: the host's
              deep scope "mpm.f.material" with one ElemScope per material call). Every step is
              `mpm analyze 1 dt` under its own profiler run.
driver runs   the single-point driver (models.point_build): one zero-free-DOF stdBrick on the
              fork's confine-first path, deviatoric leg profiled -deep; the per-call cost of
              ElasticIsotropic, DruckerPrager, LadrunoJ2 and LadrunoSANISAND from
              `elem.update` (per class tag) / 8 points and `commit.elements` / 8.

The SANISAND row is a projection (SANISAND is not admitted at particles in v1); it is composed
and validated on DruckerPrager in the parent (`projection`), two ways, both reported.

    python3.12 Ladruno_scripts/mpm_m0b/m2_explicit.py [--pyd PATH] [--smoke] [--only ...]
"""
from __future__ import annotations

import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C       # noqa: E402
import models as Mo      # noqa: E402

NS = (1000, 10000, 100000)
MATS = ("elastic", "dp", "j2")
DRIVER_MATS = ("elastic", "dp", "j2", "sanisand")
WARMUP, TIMED = 5, 20
SMOKE_WARMUP, SMOKE_TIMED = 2, 5
STD_BRICK_GP = 8


DRIVER_SCALES = (0.5, 1.0, 2.0)     # the plastic per-call sensitivity band (round 2, item 6)


def configs(smoke):
    out = []
    for N in ((1000,) if smoke else NS):
        for mat in MATS:
            for prof in ("coarse", "deep"):
                out.append(dict(kind="column", N=N, mat=mat, threads=1, profile=prof))
        for prof in ("coarse", "deep"):
            out.append(dict(kind="column", N=N, mat="elastic", threads=None, profile=prof))
    for mat in DRIVER_MATS:                       # after every column: they set the increment
        for sc in ((1.0,) if mat == "elastic" else DRIVER_SCALES):
            out.append(dict(kind="driver", N=0, mat=mat, scale=sc, threads=1, profile="deep"))
    return out


# --------------------------------------------------------------------------- #
#  child                                                                      #
# --------------------------------------------------------------------------- #
def child(cfg):
    res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")}}
    try:
        ops = C.import_ops(cfg["pyd"], cfg.get("_affinity"))
        res["conditions"] = C.child_conditions(ops)
        if cfg["kind"] == "driver":
            child_driver(ops, cfg, res)
        else:
            child_column(ops, cfg, res)
    except BaseException as e:
        res["error"] = repr(e)
    C.write_child_result(cfg, res)
    return 0


def _yield(ops, mat, params):
    if mat not in ("dp", "j2"):
        return None
    return Mo.yield_fraction(C.particle_field(ops, "sig"), mat, params[mat])


def _median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return None
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def child_column(ops, cfg, res):
    """Warm-up steps in one unprofiled call; the timed steps in ONE `mpm analyze` call under one
    profiler run with -perStep (the per-step wall series), so no report or reset interrupts the
    stepping (a report between steps evicts the caches the next step then pays for). One more,
    unprofiled step afterwards measures the per-step particle strain increment (epsHat)."""
    cmds, params = Mo.column_commands(cfg["N"], cfg["mat"])
    res["model"] = params
    res["commands"] = [[c, [str(a) for a in args]] for c, args in cmds]
    for name, args in cmds:
        getattr(ops, name)(*args)
    res["numParticles"] = C.info(ops, "numParticles")
    dtc = float(ops.mpm("criticalTimeStep"))
    dt = Mo.DT_FACTOR * dtc
    res["dt"], res["dt_crit"] = dt, dtc
    t0 = time.perf_counter()
    ops.mpm("analyze", cfg["warmup"], dt)
    res["warmup_py_wall_s"] = time.perf_counter() - t0
    yf0 = _yield(ops, cfg["mat"], params)
    flags = ["-perStep"] + (["-deep"] if cfg["profile"] == "deep" else [])
    prof = C.StepProfiler(ops, cfg["_h5"], flags)
    _, wall, run = prof(lambda: ops.mpm("analyze", cfg["timed"], dt))
    yf1 = _yield(ops, cfg["mat"], params)
    res["timed"] = {"run": run, "py_wall_s": wall, "steps": cfg["timed"]}
    res["memory_end"] = C.process_memory()
    # the plastic share phi that weights the driver's per-call costs in the projection: the
    # fraction of particles on their yield surface, at the start and the end of the timed window
    res["yield_fraction_timed_start"] = yf0
    res["yield_fraction_end"] = yf1
    res["yield_fraction_timed_mean"] = (0.5 * (yf0 + yf1)) if (yf0 is not None and yf1 is not None) else None
    # the per-step strain increment the driver matches (round 2, item 6): one more step
    e0 = C.particle_field(ops, "epsHat")
    ops.mpm("analyze", 1, dt)
    e1 = C.particle_field(ops, "epsHat")
    nv = 6
    norms = [Mo.strain_norm([e1[nv * p + k] - e0[nv * p + k] for k in range(nv)])
             for p in range(len(e0) // nv)]
    res["strain_increment_norm"] = {"median": _median(norms), "mean": sum(norms) / len(norms) if norms else None,
                                    "max": max(norms) if norms else None,
                                    "definition": "per particle, ||epsHat(n+1) - epsHat(n)|| with "
                                                  "engineering shear, one step after the timed window"}
    if cfg["mat"] in ("dp", "j2"):
        sig = C.particle_field(ops, "sig")
        on = [p for p in range(len(norms))
              if Mo.yield_fraction(sig[6 * p:6 * p + 6], cfg["mat"], params[cfg["mat"]]) == 1.0]
        res["strain_increment_norm"]["median_yielded"] = _median([norms[p] for p in on])
        res["strain_increment_norm"]["n_yielded"] = len(on)
    res["ok"] = True


def child_driver(ops, cfg, res):
    """Two passes over the same deterministic path. Pass 1 (unprofiled) classifies every step
    elastic or plastic; pass 2 rebuilds the point and profiles each contiguous run of
    same-class steps of one leg as one profiler run (no report between its steps)."""
    mat = cfg["mat"]
    sanisand = mat == "sanisand"
    inc = cfg["increment"] * cfg["scale"]
    legs, _ = Mo.point_legs(mat, inc)

    def stage(s):
        if sanisand:      # explicit stage switches, as the fork's deck does (static mElastFlag)
            ops.updateMaterialStage("-material", 1, "-stage", s)

    res["path"] = Mo.point_build(ops, mat, inc)
    res["material_command"] = [str(a) for a in Mo.point_material(mat)]
    res["increment"] = {"matched": cfg["increment"], "scale": cfg["scale"], "used": inc,
                        "source": cfg.get("increment_source")}
    flags = []
    for name, n, st, _ in legs:
        stage(st)
        for i in range(n):
            if ops.analyze(1) != 0:
                raise RuntimeError(f"{name} step {i + 1} failed (pass 1)")
            stress = [float(v) for v in ops.eleResponse(1, "material", 1, "stress")]
            flags.append((st == 1) if sanisand else Mo.point_on_yield(mat, stress))
    res["stress_end"] = [float(v) for v in ops.eleResponse(1, "material", 1, "stress")]
    # pass 2
    Mo.point_build(ops, mat, inc)
    prof = C.StepProfiler(ops, cfg["_h5"], ["-deep"])
    segs, k = [], 0
    for name, n, st, leg_inc in legs:
        stage(st)
        i = 0
        while i < n:
            plastic = flags[k + i]
            j = i
            while j < n and flags[k + j] == plastic:
                j += 1
            m = j - i

            def block(m=m, name=name):
                for _ in range(m):
                    if ops.analyze(1) != 0:
                        raise RuntimeError(f"{name} step failed (pass 2)")
            _, wall, run = prof(block)
            segs.append({"leg": name, "first": i, "steps": m, "plastic": plastic, "run": run,
                         "increment_norm": leg_inc, "py_wall_s": wall})
            i = j
        k += n
    res["steps"] = segs           # one row per profiled segment
    res["ok"] = True


# --------------------------------------------------------------------------- #
#  parent                                                                     #
# --------------------------------------------------------------------------- #
def column_window(tree, N):
    """The timed window (one profiler run over `timed` steps): totals, with the step count."""
    st = C.sub(tree, "step")
    if st is None:
        return {"ok": False}
    adv = C.sub(st, "mpm.advance")
    mat = C.sub(adv, "mpm.f.material")
    elem = (mat or {}).get("elem") or []
    m_wall = sum(r["wall_ns"] for r in elem)
    m_calls = sum(r["count"] for r in elem)
    return {"ok": True, "steps": st["calls"], "step_ns": st["wall_ns"], "advance_ns": C.wall(adv),
            "commit_ns": C.wall(C.sub(st, "commit")), "material_ns": m_wall if mat else None,
            "material_calls": m_calls if mat else None,
            "material_calls_per_particle_step": (m_calls / (N * st["calls"])) if (mat and N and st["calls"]) else None,
            "elem_rows": elem,
            "foreign_top_level": sorted(set(C.top_level_names(tree)) - C.EXPECTED_TOP)}


def driver_percall(tree):
    """One driver step's per-call cost [ns]: trial = elem.update wall / (element updates x 8
    points), commit = commit.elements wall / (commits x 8). The element update includes the
    stdBrick's own strain computation per point (its overhead; see `projection`)."""
    rows = []

    def collect(n):
        for name, c in n["children"].items():
            if name == "elem.update":
                rows.extend(c.get("elem") or [])
            collect(c)
    collect(tree)
    upd_wall = sum(r["wall_ns"] for r in rows)
    upd_cnt = sum(r["count"] for r in rows)
    com_wall = C.total(tree, "commit.elements")
    com_cnt = C.total_calls(tree, "commit.elements")
    return {"trial_wall_ns": upd_wall, "trial_calls": upd_cnt * STD_BRICK_GP,
            "commit_wall_ns": com_wall, "commit_calls": com_cnt * STD_BRICK_GP}


def driver_summary(r):
    """Per-call costs (trial + commit, ns per point-step): elastic over every elastic step,
    plastic over the measure leg's plastic steps (the matched increment) only."""
    acc = {}
    for s in r.get("steps", []):
        d = s.get("derived")
        if not d:
            continue
        if s["plastic"] and s["leg"] != "measure":
            key = "plastic_pre"
        else:
            key = "plastic" if s["plastic"] else "elastic"
        a = acc.setdefault(key, {"tw": 0, "tc": 0, "cw": 0, "cc": 0, "n": 0})
        a["tw"] += d["trial_wall_ns"]
        a["tc"] += d["trial_calls"]
        a["cw"] += d["commit_wall_ns"]
        a["cc"] += d["commit_calls"]
        a["n"] += s["steps"]
    out = {}
    for key, a in acc.items():
        if a["tc"] and a["cc"]:
            out[f"{key}_trial_ns"] = a["tw"] / a["tc"]
            out[f"{key}_commit_ns"] = a["cw"] / a["cc"]
            out[f"{key}_per_point_step_ns"] = a["tw"] / a["tc"] + a["cw"] / a["cc"]
            out[f"{key}_steps"] = a["n"]
    inc = (r.get("increment") or {}).get("used")
    if inc is not None:
        out["measure_increment_norm"] = inc
    return out


def tag_of(c, rep):
    th = "def" if c["threads"] is None else c["threads"]
    sc = f"_x{c['scale']:g}" if c["kind"] == "driver" else ""
    return f"m2_{c['kind']}_N{c['N']}_{c['mat']}{sc}_t{th}_{c['profile']}_r{rep}"


def matched_increment(summary):
    """The DP column's median per-step particle strain-increment norm, at the largest N that
    ran (coarse run, median over repeats). None when no DP column ran."""
    rows = [s for s in summary if s["config"]["kind"] == "column" and s["config"]["mat"] == "dp"
            and s["config"]["profile"] == "coarse" and s.get("strain_increment_norm_median")]
    if not rows:
        return None, None
    s = max(rows, key=lambda s: s["config"]["N"])
    return s["strain_increment_norm_median"]["median"], f"DP column, N = {s['config']['N']}, coarse"


def main():
    ap = C.base_parser(__doc__)
    ap.add_argument("--warmup", type=int, default=None)
    ap.add_argument("--timed", type=int, default=None)
    ap.add_argument("--increment", type=float, default=None,
                    help="the driver's matched strain-increment norm (default: measured in this run's "
                         "DP column)")
    args = C.parse_args(ap)
    if args.child:
        return child(C.load_cfg(args.child))
    default_reps = 1 if args.smoke else 3
    warm = args.warmup if args.warmup is not None else (SMOKE_WARMUP if args.smoke else WARMUP)
    timed = args.timed if args.timed is not None else (SMOKE_TIMED if args.smoke else TIMED)
    cfgs = C.select(configs(args.smoke), args)
    if args.list:
        return C.list_runs(cfgs, args, "m2", tag_of, default_reps, {"warmup": warm, "timed": timed})
    od = C.out_dir(args, "m2" + ("_smoke" if args.smoke else ""))
    doc = {"script": "m2_explicit", "conditions": C.parent_conditions(args),
           "pinned": {"threads": "1, plus the default (unset) row for elastic", "warmup_steps": warm,
                      "timed_steps": timed, "repeats": args.repeats or default_reps,
                      "throughput": "N * timed / sum of the profiler's `step` wall (coarse run)",
                      "material share": "sum of mpm.f.material elem_by_type wall / sum of step wall"
                                        " (deep run)"},
           "runs": [], "summary": []}
    t0 = time.time()
    for c in cfgs:                         # columns first: the driver needs their increment
        per_rep = []
        extra = {}
        if c["kind"] == "driver":
            inc, src = (args.increment, "--increment") if args.increment else matched_increment(doc["summary"])
            if inc is None:
                inc, src = Mo.POINT_DEFAULT_INCREMENT, "default (no DP column in this run)"
            extra = {"increment": inc, "increment_source": src}
            doc["driver_increment"] = {"norm": inc, "source": src}
        for rep in range(C.repeats_for(c, args, default_reps)):
            cfg = dict(c, pyd=args.pyd, warmup=warm, timed=timed, **extra)
            tag = tag_of(c, rep)
            print(f"[m2] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout, args.affinity_mask)
            prof, perr = C.safe_read_profile(str(od / f"{tag}.h5"))
            if perr:
                r["profile_error"] = perr
            N = r.get("numParticles") or c["N"]
            if c["kind"] == "column":
                p = prof.get((r.get("timed") or {}).get("run", ""))
                r["profile"] = p["tree"] if p else None
                r["series"] = p["series"] if p else None          # the per-step wall rows
                r["window"] = column_window(p["tree"], N) if p and p["tree"] else {"ok": False}
                r["run_summary"] = column_summary(r, N)
            else:
                for s in r.get("steps", []):
                    p = prof.get(s["run"])
                    s["profile"] = p["tree"] if p else None
                    s["derived"] = driver_percall(p["tree"]) if p and p["tree"] else None
                r["run_summary"] = driver_summary(r)
            r["tag"] = tag
            doc["runs"].append(r)
            per_rep.append(r)
            print(f"[m2]   ok={r.get('ok')} {r.get('reason', '')} {short(r.get('run_summary'))} "
                  f"({r.get('child_wall_s', 0):.1f} s)", flush=True)
        doc["summary"].append(summarize(c, per_rep))
        C.write_json(od / "m2_results.json", doc)
    doc["projection"] = projection(doc["summary"])
    doc["wall_s"] = time.time() - t0
    C.write_json(od / "m2_results.json", doc)
    print(f"[m2] wrote {od / 'm2_results.json'}")
    return 0


def short(s):
    if not s:
        return "-"
    keys = ("particle_steps_per_s", "material_share", "elastic_per_point_step_ns",
            "plastic_per_point_step_ns")
    return " ".join(f"{k}={s.get(k):.4g}" for k in keys if isinstance(s.get(k), (int, float)))


def column_summary(r, N):
    w = r.get("window") or {}
    if not w.get("ok") or not w.get("steps"):
        return None
    n, step = w["steps"], w["step_ns"]
    py = (r.get("timed") or {}).get("py_wall_s")
    out = {"n_timed_steps": n, "step_ns_mean": step / n,
           "particle_steps_per_s": N * n / (step * 1e-9) if step else None,
           "particle_steps_per_s_python_wall": N * n / py if py else None,
           "ns_per_particle_step": step / (N * n)}
    # robust to a single descheduled step: the median of the per-step series (-perStep wall_ms)
    ser = (r.get("series") or {}).get("wall_ms") or []
    if ser:
        med_ns = statistics.median(ser) * 1.0e6
        out["ns_per_particle_step_median"] = med_ns / N
        out["particle_steps_per_s_median"] = N / (med_ns * 1e-9) if med_ns else None
        out["step_ns_max_over_median"] = max(ser) * 1.0e6 / med_ns if med_ns else None
    if w["material_ns"] is not None:
        mat = w["material_ns"]
        out["material_ns_per_particle_step"] = mat / (N * n)
        out["nonmaterial_ns_per_particle_step"] = (step - mat) / (N * n)
        out["material_share"] = mat / step if step else None
        out["nonmaterial_share"] = 1.0 - out["material_share"] if step else None
        out["material_calls_per_particle_step"] = w["material_calls_per_particle_step"]
        if ser:
            out["nonmaterial_ns_per_particle_step_median"] = out["ns_per_particle_step_median"] *                 out["nonmaterial_share"]
    if r.get("yield_fraction_timed_mean") is not None:
        out["yield_fraction_timed_mean"] = r["yield_fraction_timed_mean"]
    for k, v in (r.get("strain_increment_norm") or {}).items():
        if isinstance(v, (int, float)):
            out["strain_increment_norm_" + k] = v
    return out


def summarize(c, reps):
    rs = [r.get("run_summary") or {} for r in reps]
    keys = set()
    for x in rs:
        keys |= {k for k, v in x.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
    out = {"config": c, "ok_runs": sum(1 for r in reps if r.get("ok")), "repeats": len(reps),
           "failures": [{"rep": i, "reason": r.get("reason"), "child_wall_s": r.get("child_wall_s")}
                        for i, r in enumerate(reps) if not r.get("ok")]}
    for k in sorted(keys):
        out[k] = C.median_spread([x.get(k) for x in rs])
    yf = [r.get("yield_fraction_end") for r in reps if r.get("yield_fraction_end") is not None]
    if yf:
        out["yield_fraction_end"] = C.median_spread(yf)
    return out


def projection(summary):
    """SANISAND, projected; the composition validated on DruckerPrager (the 232g M2 rule: a DP
    projection that misses the DP measurement by > 25 % makes the SANISAND row unreliable).

    The driver gives each material an elastic and a plastic per-call cost c_el(X), c_pl(X)
    (trial + commit per point-step). A particle-step of material X with plastic share phi costs
        c(X, phi) = (1 - phi) c_el(X) + phi c_pl(X)
    phi(DP) is the DP column's measured mean yield fraction over the timed steps (same N);
    SANISAND is plastic whenever it loads (a bounding-surface model), so phi = 1 is the row, and
    phi = phi(DP) is reported beside it as the lower end.

    literal      t(X) = nonmaterial(elastic column, deep run) + c(X, phi)
    difference   t(X) = t(elastic column, coarse run) + c(X, phi) - c_el(elastic)
                 (the stdBrick's own per-point strain work, inside every driver figure, cancels)
    Both forms are reported; the orchestrator picks."""
    def get(kind, N, mat, prof, threads=1, scale=1.0):
        for s in summary:
            c = s["config"]
            if c["kind"] == kind and c["N"] == N and c["mat"] == mat and c["profile"] == prof \
                    and c["threads"] == threads and (kind != "driver" or c.get("scale") == scale):
                return s
        return None

    def med(s, k):
        return ((s or {}).get(k) or {}).get("median") if s else None

    drv = {}
    for m in DRIVER_MATS:
        s = get("driver", 0, m, "deep")
        drv[m] = {"elastic": med(s, "elastic_per_point_step_ns"),
                  "plastic": med(s, "plastic_per_point_step_ns")}

    def cost(m, phi):
        el, pl = drv[m]["elastic"], drv[m]["plastic"]
        if phi is None:
            return None
        if pl is None:
            return el if phi == 0 else None
        if el is None:
            return pl if phi == 1 else None
        return (1.0 - phi) * el + phi * pl

    band = {m: {f"x{sc:g}": med(get("driver", 0, m, "deep", scale=sc), "plastic_per_point_step_ns")
                for sc in DRIVER_SCALES} for m in DRIVER_MATS if m != "elastic"}
    out = {"driver_per_call_ns": drv,
           "driver_plastic_per_call_band_ns": band,
           "sanisand_over_dp_ratio_band": {k: (band["sanisand"][k] / band["dp"][k])
                                           if band["sanisand"][k] and band["dp"][k] else None
                                           for k in band["dp"]},
           "sanisand_over_dp_cost_ratio_plastic":
               (drv["sanisand"]["plastic"] / drv["dp"]["plastic"])
               if drv["sanisand"]["plastic"] and drv["dp"]["plastic"] else None,
           "rows": []}
    c_el_elastic = drv["elastic"]["elastic"]
    for N in sorted({s["config"]["N"] for s in summary if s["config"]["kind"] == "column"}):
        el_deep, el_coarse = get("column", N, "elastic", "deep"), get("column", N, "elastic", "coarse")
        dp_coarse = get("column", N, "dp", "coarse")
        # per-step medians (a descheduled step must not decide the rule); the means are in the summary
        nonmat = med(el_deep, "nonmaterial_ns_per_particle_step_median")
        t_el = med(el_coarse, "ns_per_particle_step_median")
        t_dp = med(dp_coarse, "ns_per_particle_step_median")
        phi_dp = med(dp_coarse, "yield_fraction_timed_mean")
        row = {"N": N, "nonmaterial_elastic_ns": nonmat, "t_elastic_ns": t_el,
               "t_dp_measured_ns": t_dp, "phi_dp": phi_dp,
               "elastic_nonmaterial_share_t1": med(el_deep, "nonmaterial_share")}
        for form in ("literal", "difference"):
            def compose(cx):
                if cx is None:
                    return None
                if form == "literal":
                    return nonmat + cx if nonmat is not None else None
                return t_el + cx - c_el_elastic if (t_el is not None and c_el_elastic is not None) else None
            p_dp = compose(cost("dp", phi_dp))
            miss = abs(p_dp - t_dp) / t_dp if (p_dp is not None and t_dp) else None
            p_s1 = compose(cost("sanisand", 1.0))
            p_sl = compose(cost("sanisand", phi_dp))
            row[form] = {"t_dp_projected_ns": p_dp, "dp_miss": miss,
                         "sanisand_reliable": (miss <= 0.25) if miss is not None else None,
                         "t_sanisand_projected_ns_phi1": p_s1,
                         "t_sanisand_projected_ns_phi_dp": p_sl,
                         "sanisand_particle_steps_per_s_phi1": (1e9 / p_s1) if p_s1 else None}
        out["rows"].append(row)
    return out


if __name__ == "__main__":
    sys.exit(main())
