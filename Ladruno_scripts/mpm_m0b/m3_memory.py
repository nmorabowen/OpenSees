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

import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C       # noqa: E402
import models as Mo      # noqa: E402

NS = (1000, 10000, 100000)
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
        ops = C.import_ops(cfg["pyd"], cfg.get("_affinity"))
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
            # the support-size histogram (nodes per particle element): numDOF = 3 nen decides
            # whether an FE_Element owns its tangent (numDOF > 64) or shares the class pool
            hist = {}
            for t in ops.getEleTags():
                n = len(ops.eleNodes(int(t)))
                hist[n] = hist.get(n, 0) + 1
            res["support_histogram"] = {str(k): v for k, v in sorted(hist.items())}
        else:
            # mode F has no element to ask: count each particle's uGIMP support from its position
            # (node i is in the support iff |x - x_i| < h + lp, per axis; origin 0)
            x, lp, h = C.particle_field(ops, "x"), C.particle_field(ops, "lp"), params["h"]
            hist = {}
            for p in range(len(x) // 3):
                n = 1
                for a in range(3):
                    lo = math.floor((x[3 * p + a] - h - lp[3 * p + a]) / h) + 1
                    hi = math.ceil((x[3 * p + a] + h + lp[3 * p + a]) / h) - 1
                    n *= max(hi - lo + 1, 0)
                hist[n] = hist.get(n, 0) + 1
            res["support_histogram"] = {str(k): v for k, v in sorted(hist.items())}
        res["ok"] = True
    except BaseException as e:
        res["error"] = repr(e)
    C.write_child_result(cfg, res)
    return 0


# --------------------------------------------------------------------------- #
#  the analytical per-particle byte count (round 2, item 7)                   #
# --------------------------------------------------------------------------- #
# Estimates for MSVC x64, counted from the headers (sizeof is not measurable from Python):
# a heap block costs its payload + HEAP; a Vector object 24 B (sz, theData, fromFree) plus its
# data block; a Matrix object 32 B plus its data block; an ID 32 B plus 4 B per entry; a
# std::vector 24 B plus its block; a polymorphic base (TaggedObject + MovableObject) 48 B.
HEAP = 16
VEC_OBJ, MAT_OBJ, ID_OBJ, STDV, BASE = 24, 32, 32, 24, 48
ANALYTIC_NOTE = {
    "heap_block_overhead_B": HEAP, "Vector_object_B": VEC_OBJ, "Matrix_object_B": MAT_OBJ,
    "ID_object_B": ID_OBJ, "std_vector_B": STDV, "polymorphic_base_B": BASE,
    "sources": "LadrunoMPMPoint.h members; FE_Element.cpp:188-219 (own tangent and residual "
               "only when numDOF > MAX_NUM_DOF = 64, else the class pool); "
               "ElasticIsotropicThreeDimensional.h (2 Vector(6)); DruckerPrager.h (8 Vector(6), "
               "4 Matrix(6,6), 26 scalars); Node.h (Crd, disp/vel/accel arrays and Vector views, "
               "unbalLoad); LadrunoMPMCore.h Particles (57 doubles + 20 B of ids per particle, held "
               "twice: the state and the `mpm reset` snapshot)",
    "not counted": "the SOE's matrix and PARDISO's factor (MKL-internal); the allocator's page "
                   "slack; Python",
}


def _vec(n):
    return VEC_OBJ + (8 * n + HEAP if n else 0)


def _mat(r, c):
    return MAT_OBJ + 8 * r * c + HEAP


def _idv(n):
    return ID_OBJ + 4 * n + HEAP


def _stdv(nbytes):
    return STDV + (nbytes + HEAP if nbytes else 0)


def material_bytes(mat):
    if mat == "elastic":        # ElasticIsotropicThreeDimensional: E, v, rho, parameterID + 2 Vector(6)
        return BASE + 32 + 2 * (_vec(6) - VEC_OBJ) + 2 * VEC_OBJ + HEAP
    # DruckerPrager3D: 26 scalars, 8 Vector(6) (mState taken as 6), 4 Matrix(6,6)
    return BASE + 26 * 8 + 8 * _vec(6) + 4 * _mat(6, 6) + HEAP


def element_bytes(nen, ndm=3):
    """LadrunoMPMPoint with nen support nodes (its own storage; the tangent is a class pool)."""
    fixed = BASE + 64 + 16 + 8 * (9 * 3 + 1 + 18 + 3 + 3) + 16      # ints, link, State, caches
    per = (_idv(nen) + _stdv(8 * nen)                                # connect, nodes
           + 5 * 0 + _stdv(8 * nen) + _stdv(8 * ndm * nen) * 3 + _stdv(8 * nen)   # S dS vHat aHat ghostK
           + _stdv(8 * ndm * nen) * 3                                # g, fInt, uT
           + 2 * (_vec(ndm * nen) + HEAP))                           # P, Pinc (new Vector)
    return fixed + per + HEAP


def fe_element_bytes(nen, ndf=3):
    n = ndf * nen
    own = (_mat(n, n) + HEAP + _vec(n) + HEAP) if n > 64 else 0
    return 120 + _idv(nen) + _idv(n) + own + HEAP, own


NODE_BYTES = (320 + _vec(3) + HEAP + (8 * 12 + HEAP) + 4 * (VEC_OBJ + HEAP)
              + (8 * 6 + HEAP) + 2 * (VEC_OBJ + HEAP) + (8 * 6 + HEAP) + 2 * (VEC_OBJ + HEAP)
              + _vec(3) + HEAP)
DOF_GROUP_BYTES = 120 + _idv(3) + HEAP
GRID_NODE_HOST_BYTES = 64 + 48 + 48     # liveTag map entry, liveTags set entry, Domain's tag map
PARTICLE_STATE_BYTES = 2 * (57 * 8 + 20) + (8 + 1 + 4 + 24 + 8 + 1 + 40)


def analytic(r):
    """Per-particle bytes by item, from the run's own counts (support histogram, grid nodes)."""
    c, m = r["config"], r.get("model") or {}
    N = r.get("numParticles") or c["N"]
    out = {"material_copy": material_bytes(c["mat"]), "particle_state_and_snapshot": PARTICLE_STATE_BYTES}
    if c["mode"] == "A":
        hist = {int(k): v for k, v in (r.get("support_histogram") or {}).items()}
        tot = sum(hist.values()) or 1
        el = sum(v * element_bytes(k) for k, v in hist.items()) / tot
        fe = sum(v * fe_element_bytes(k)[0] for k, v in hist.items()) / tot
        fe_own = sum(v * fe_element_bytes(k)[1] for k, v in hist.items()) / tot
        g = (r.get("gridNodes") or 0) / N
        out.update({"element": el, "fe_element": fe, "of_which_fe_own_tangent": fe_own,
                    "grid_nodes": g * (NODE_BYTES + DOF_GROUP_BYTES + GRID_NODE_HOST_BYTES),
                    "grid_nodes_per_particle": g,
                    "share_numDOF_gt_64": sum(v for k, v in hist.items() if 3 * k > 64) / tot,
                    "mean_support_nodes": sum(k * v for k, v in hist.items()) / tot})
    else:
        hist = {int(k): v for k, v in (r.get("support_histogram") or {}).items()}
        tot = sum(hist.values()) or 1
        nen = sum(k * v for k, v in hist.items()) / tot if hist else 27   # counted after step 1
        gc = m.get("grid_cells") or [0, 0, 0]
        gnodes = 1
        for n in gc:
            gnodes *= (n + 3)        # real nodes plus the ghost layer on both faces
        out.update({"support_arrays": 8 + nen * (8 + 8 + 24),
                    "trial_arrays": 43 * 8,
                    "dense_grid": gnodes * (19 * 8 + 2) / N, "dense_grid_nodes_per_particle": gnodes / N,
                    "mean_support_nodes": nen})
    out["total"] = sum(v for k, v in out.items()
                       if k not in ("of_which_fe_own_tangent", "grid_nodes_per_particle",
                                    "share_numDOF_gt_64", "mean_support_nodes",
                                    "dense_grid_nodes_per_particle"))
    return out


def tag_of(c, rep):
    return f"m3_{c['mode']}_N{c['N']}_{c['mat']}_r{rep}"


def main():
    ap = C.base_parser(__doc__)
    args = C.parse_args(ap)
    if args.child:
        return child(C.load_cfg(args.child))
    default_reps = 1                   # the plan passes 2
    cfgs = C.select(configs(args.smoke), args)
    if args.list:
        return C.list_runs(cfgs, args, "m3", tag_of, default_reps)
    od = C.out_dir(args, "m3" + ("_smoke" if args.smoke else ""))
    doc = {"script": "m3_memory", "conditions": C.parent_conditions(args),
           "pinned": {"threads": "mode A 8, mode F 1", "repeats": args.repeats or default_reps,
                      "definition": "slope of private bytes against N after the first step (s2), "
                                    "fitted on the two largest N (10^4 -> 10^5)"},
           "analytic_constants": ANALYTIC_NOTE,
           "runs": [], "slopes": []}
    t0 = time.time()
    rows = []
    for c in cfgs:
        for rep in range(C.repeats_for(c, args, default_reps)):
            cfg = dict(c, pyd=args.pyd)
            tag = tag_of(c, rep)
            print(f"[m3] {tag} ...", flush=True)
            r, _ = C.run_child(__file__, cfg, od, tag, args.timeout, args.affinity_mask)
            prof, perr = C.safe_read_profile(str(od / f"{tag}.h5"))
            if perr:
                r["profile_error"] = perr
            if r.get("ok"):
                r["analytic"] = analytic(r)
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
    doc["failures"] = [{"tag": r.get("tag"), "reason": r.get("reason"), "error": r.get("error"),
                        "child_wall_s": r.get("child_wall_s")} for r in rows if not r.get("ok")]
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
        an = byN[n2][0].get("analytic") or {}
        meas = sl.get("private.s2_step1")
        out.append({"mode": mode, "mat": mat, "N": [n1, n2],
                    "bytes_per_particle_after_step1": meas,
                    "analytic_per_particle": an,
                    "measured_minus_analytic": (meas - an["total"]) if (meas is not None and an) else None,
                    "note": "the difference is what the analytic count leaves out: mode A the SOE "
                              "matrix and the factor (MKL), mode F allocator slack; a negative one "
                              "means an over-count (e.g. the 27-node support assumed in mode F)",
                    "slopes_per_particle": sl,
                    "per_particle_at_N": {n: {k: v / n for k, v in _avg([_items(r) for r in byN[n]]).items()}
                                          for n in Ns}})
    return out


def _avg(dicts):
    keys = set().union(*dicts) if dicts else set()
    return {k: sum(d[k] for d in dicts if k in d) / sum(1 for d in dicts if k in d) for k in keys}


if __name__ == "__main__":
    sys.exit(main())
