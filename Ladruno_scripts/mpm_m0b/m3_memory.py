"""WP-232g M3: memory per particle-element.

Mode A (the M1 block, quasi-static, Pardiso, renormalize) and mode F (the M2 column), elastic
and DruckerPrager, at N = 10^4 and 10^5 (smoke: 10^3 and 2 x 10^3, for a slope). One fresh
process per run. Private bytes (GetProcessMemoryInfo PrivateUsage) are read at four stages:

  s0  after `wipe` / `model`              (the interpreter, the binary, MKL)
  s1  after every model command           (particles, their material copies, the analysis objects)
  s2  after the first step                (mode A: elements, grid nodes, FE_Elements, DOF_Groups,
                                           the SOE and the factor; mode F: the dense grid arrays)
  s3  after the second step

and three profiler `-memory` windows (the model commands, step 1, step 2) give what Matrix /
Vector / ID storage each window still holds and the live-object census per class tag. The parent
reports the slope of private bytes against N after the first step (the 232g definition), and
the slope of every stage and every counted item, so the per-particle bytes can be attributed.

    python3.12 Ladruno_scripts/mpm_m0b/m3_memory.py [--pyd PATH] [--smoke] [--only ...]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C       # noqa: E402
import models as Mo      # noqa: E402

NS = (10000, 100000)
SMOKE_NS = (1000, 2000)


def configs(smoke):
    out = []
    for N in (SMOKE_NS if smoke else NS):
        for mat in ("elastic", "dp"):
            out.append(dict(mode="A", N=N, mat=mat, threads=8))
            out.append(dict(mode="F", N=N, mat=mat, threads=1))
    return out


def child(cfg):
    res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")}}
    try:
        ops = C.import_ops(cfg["pyd"])
        res["conditions"] = C.child_conditions(ops)
        if cfg["mode"] == "A":
            cmds, params = Mo.block_commands(cfg["N"], cfg["mat"], "qs", "Pardiso", "renormalize")
        else:
            cmds, params = Mo.column_commands(cfg["N"], cfg["mat"])
        res["model"] = params
        stages = {}
        getattr(ops, cmds[0][0])(*cmds[0][1])           # wipe
        getattr(ops, cmds[1][0])(*cmds[1][1])           # model
        stages["s0_model"] = C.process_memory()
        prof = C.StepProfiler(ops, cfg["_h5"], ["-memory"])

        def build():
            for name, args in cmds[2:]:
                getattr(ops, name)(*args)
        _, _, run_build = prof(build)
        stages["s1_built"] = C.process_memory()
        res["numParticles"] = C.info(ops, "numParticles")
        dt = params["dt"] if cfg["mode"] == "A" else Mo.DT_FACTOR * float(ops.mpm("criticalTimeStep"))
        _, _, run1 = prof(lambda: ops.mpm("analyze", 1, dt))
        stages["s2_step1"] = C.process_memory()
        _, _, run2 = prof(lambda: ops.mpm("analyze", 1, dt))
        stages["s3_step2"] = C.process_memory()
        res["stages"] = stages
        res["runs"] = {"build": run_build, "step1": run1, "step2": run2}
        if cfg["mode"] == "A":
            res["gridNodes"] = C.info(ops, "gridNodes")
            try:
                res["nEqn"] = int(ops.systemSize())
            except Exception:
                res["nEqn"] = None
        res["ok"] = True
    except BaseException as e:
        res["error"] = repr(e)
    C.write_child_result(cfg, res)
    return 0


def tag_of(c, rep):
    return f"m3_{c['mode']}_N{c['N']}_{c['mat']}_r{rep}"


def main():
    ap = C.base_parser(__doc__)
    args = ap.parse_args()
    if args.child:
        return child(C.load_cfg(args.child))
    reps = args.repeats or 1           # memory is deterministic enough; --repeats 3 to check
    od = C.out_dir(args, "m3" + ("_smoke" if args.smoke else ""))
    doc = {"script": "m3_memory", "conditions": C.parent_conditions(args),
           "pinned": {"threads": "mode A 8, mode F 1", "repeats": reps,
                      "definition": "slope of private bytes against N after the first step (s2)"},
           "runs": [], "slopes": []}
    t0 = time.time()
    rows = []
    for c in configs(args.smoke):
        if not C.keep(c, args.only):
            continue
        for rep in range(reps):
            cfg = dict(c, pyd=args.pyd)
            tag = tag_of(c, rep)
            print(f"[m3] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout)
            prof = C.read_profile(str(od / f"{tag}.h5"))
            r["profiler_memory"] = {k: (prof.get(v) or {}).get("memory")
                                    for k, v in (r.get("runs") or {}).items()}
            r["tag"] = tag
            doc["runs"].append(r)
            rows.append(r)
            s2 = ((r.get("stages") or {}).get("s2_step1") or {}).get("private")
            print(f"[m3]   ok={r.get('ok')} private after step 1 = {s2} "
                  f"({r.get('child_wall_s', 0):.1f} s)", flush=True)
        C.write_json(od / "m3_results.json", doc)
    doc["slopes"] = slopes(rows)
    doc["wall_s"] = time.time() - t0
    C.write_json(od / "m3_results.json", doc)
    print(f"[m3] wrote {od / 'm3_results.json'}")
    return 0


def _items(r):
    """Every quantity of one run that may scale with N, as a flat dict."""
    st = r.get("stages") or {}
    out = {}
    for s, m in st.items():
        if m and "private" in m:
            out["private." + s] = m["private"]
    if all(("private." + s) in out for s in ("s0_model", "s1_built", "s2_step1")):
        out["delta.build"] = out["private.s1_built"] - out["private.s0_model"]
        out["delta.step1"] = out["private.s2_step1"] - out["private.s1_built"]
    for win, mem in (r.get("profiler_memory") or {}).items():
        if not mem:
            continue
        for k in ("matrix_live", "vector_live", "id_live", "live_bytes", "peak_bytes"):
            if k in mem:
                out[f"profiler.{win}.{k}"] = mem[k]
        for tag, cnt in (mem.get("components_live") or {}).items():
            out[f"census.{win}.classTag{tag}"] = cnt
    if r.get("gridNodes") is not None:
        out["gridNodes"] = r["gridNodes"]
    if r.get("nEqn") is not None:
        out["nEqn"] = r["nEqn"]
    return out


def slopes(rows):
    """Per (mode, mat): d(quantity)/dN between the two largest N (bytes or counts per
    particle), plus the plain per-particle value at each N."""
    groups = {}
    for r in rows:
        if not r.get("ok"):
            continue
        c = r["config"]
        groups.setdefault((c["mode"], c["mat"]), {}).setdefault(r.get("numParticles") or c["N"], []).append(r)
    out = []
    for (mode, mat), byN in groups.items():
        Ns = sorted(byN)
        if len(Ns) < 2:
            out.append({"mode": mode, "mat": mat, "N": Ns, "note": "one N only: no slope"})
            continue
        n1, n2 = Ns[-2], Ns[-1]
        a = _avg([_items(r) for r in byN[n1]])
        b = _avg([_items(r) for r in byN[n2]])
        sl = {k: (b[k] - a[k]) / (n2 - n1) for k in b if k in a}
        out.append({"mode": mode, "mat": mat, "N": [n1, n2],
                    "bytes_per_particle_after_step1": sl.get("private.s2_step1"),
                    "slopes_per_particle": sl,
                    "per_particle_at_N": {n: {k: v / n for k, v in _avg([_items(r) for r in byN[n]]).items()}
                                          for n in Ns}})
    return out


def _avg(dicts):
    keys = set().union(*dicts) if dicts else set()
    return {k: sum(d[k] for d in dicts if k in d) / sum(1 for d in dicts if k in d) for k in keys}


if __name__ == "__main__":
    sys.exit(main())
