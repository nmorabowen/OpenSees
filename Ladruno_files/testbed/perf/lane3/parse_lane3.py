"""ADR-75b L3-0 — turn a profiler h5 into the Lane-3 gate table.

Usage:  python parse_lane3.py <file.h5> [--json out.json]

Why this exists (and why parse_prof.py is not enough): ADR-75 §3 gates Lane 3 on a
">40% element fraction", but the shipped profiler reports per-PHASE totals, and a phase
is not a threadable loop. `formTangent` contains BOTH the element kernel (threadable, if
the element is re-entrant) AND the addA scatter (the race §4 has to solve). Threading
buys you the first and not the second, so a phase percentage systematically OVERSTATES
the Amdahl headroom. This parser separates them.

The three threadable loops (ADR-75b §2), and how each is recovered here:

  loop A  Domain::update()                    scope "elem.update"    -- NO FP reduction
  loop B  IncrementalIntegrator::formTangent  scope "elem.tangent"   -- addA reduction
  loop C  ...::formElementResidual            scope "elem.residual"  -- addB reduction

For each, the NAMED deep scope wraps the WHOLE loop (kernel + scatter + iteration), while
the per-classTag `elem_by_type` rows inside it account ONLY for the kernel call
(getTangent / getResidual / update). So:

    kernel  = sum(elem_by_type wall under the loop scope)
    scatter = loop scope wall - kernel          <-- the non-threadable remainder

A loop scope can appear at SEVERAL places in the tree (e.g. `elem.update` under both
`update` and `newStep`; `soe.factor` under both `linearSolve` and `update`). Every
occurrence is summed, and the per-site breakdown is kept, because ADR-40b's lane-E
finding -- a hidden second factorization booked outside `linearSolve` -- was only
visible per-site.
"""
import os, sys, json, collections

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py

LOOPS = {
    "A_update":   "elem.update",
    "B_tangent":  "elem.tangent",
    "C_residual": "elem.residual",
}


def collect(path):
    """Walk the rollup; return (meta, scopes, elem_rows).

    scopes:    full/path -> {"wall_ms", "calls"}
    elem_rows: full/path of the OWNING scope -> [ {classTag, count, wall_ms, us_per_ele} ]
    """
    scopes = {}
    elem_rows = collections.defaultdict(list)

    with h5py.File(path, "r") as f:
        runs = f["runs"]
        run = runs[list(runs)[0]]
        meta = {}
        for k in run.attrs:
            v = run.attrs[k]
            meta[k] = v.item() if hasattr(v, "item") else v

        def walk(g, p=""):
            for k in g:
                it = g[k]
                if isinstance(it, h5py.Group):
                    pp = f"{p}/{k}" if p else k
                    a = dict(it.attrs)
                    if "wall_ns" in a:
                        scopes[pp] = {
                            "wall_ms": int(a["wall_ns"]) / 1e6,
                            "calls": int(a.get("calls", 0)),
                        }
                    walk(it, pp)
                elif k == "elem_by_type":
                    for row in it[()]:
                        elem_rows[p].append({
                            "classTag": int(row["classTag"]),
                            "count": int(row["count"]),
                            "wall_ms": int(row["wall_ns"]) / 1e6,
                            "us_per_ele": float(row["wall_ns_per_elem"]) / 1e3,
                        })

        if "rollup" in run and "root" in run["rollup"]:
            walk(run["rollup"]["root"], "root")

    return meta, scopes, elem_rows


def analyse(meta, scopes, elem_rows):
    # The profiled total. "root/step" is the driver-level scope every lane emits.
    step_ms = scopes.get("root/step", {}).get("wall_ms") or meta.get("wall_ms_total", 0.0)

    out = {
        "meta": {k: meta.get(k) for k in
                 ("solver", "algorithm", "integrator", "nDOF", "nSteps",
                  "threads", "wall_ms_total", "timestamp")},
        "step_ms": step_ms,
        "loops": {},
        "phases": {},
        "classtags": {},
    }

    # ---- coarse phases, for context and for the attribution gap ----
    for name in ("linearSolve", "formTangent", "formUnbalance", "update", "newStep",
                 "commit", "convTest", "soe.symbolic", "soe.factor", "soe.trisolve"):
        tot = sum(v["wall_ms"] for k, v in scopes.items() if k.rsplit("/", 1)[-1] == name)
        if tot:
            out["phases"][name] = {"ms": tot, "pct_step": 100.0 * tot / step_ms}

    # ---- the three threadable loops: kernel vs scatter ----
    for key, scope_name in LOOPS.items():
        sites = {k: v for k, v in scopes.items() if k.rsplit("/", 1)[-1] == scope_name}
        if not sites:
            continue
        loop_ms = sum(v["wall_ms"] for v in sites.values())
        kernel_ms = 0.0
        per_tag = collections.defaultdict(lambda: {"count": 0, "wall_ms": 0.0})
        for site in sites:
            for r in elem_rows.get(site, []):
                kernel_ms += r["wall_ms"]
                per_tag[r["classTag"]]["count"] += r["count"]
                per_tag[r["classTag"]]["wall_ms"] += r["wall_ms"]
        scatter_ms = loop_ms - kernel_ms
        out["loops"][key] = {
            "scope": scope_name,
            "sites": {k: v["wall_ms"] for k, v in sorted(sites.items())},
            "loop_ms": loop_ms,
            "kernel_ms": kernel_ms,
            "scatter_ms": scatter_ms,
            "pct_step_loop": 100.0 * loop_ms / step_ms,
            "pct_step_kernel": 100.0 * kernel_ms / step_ms,
            # The share of the loop that threading does NOT buy you.
            "pct_loop_scatter": (100.0 * scatter_ms / loop_ms) if loop_ms else 0.0,
            "per_classtag": {
                str(t): {
                    "count": d["count"],
                    "wall_ms": d["wall_ms"],
                    "us_per_ele": (1000.0 * d["wall_ms"] / d["count"]) if d["count"] else 0.0,
                } for t, d in sorted(per_tag.items())
            },
        }

    # ---- totals + Amdahl ----
    kernel_total = sum(l["kernel_ms"] for l in out["loops"].values())
    loop_total = sum(l["loop_ms"] for l in out["loops"].values())
    out["totals"] = {
        "kernel_ms": kernel_total,
        "loop_ms": loop_total,
        "pct_step_kernel": 100.0 * kernel_total / step_ms,
        "pct_step_loop": 100.0 * loop_total / step_ms,
    }
    # Amdahl on the KERNEL fraction only: the scatter stays serial in ordered mode and
    # is contended in fast mode, so the kernel share is the honest parallel fraction.
    p = kernel_total / step_ms
    out["amdahl_kernel_only"] = {str(n): 1.0 / ((1.0 - p) + p / n) for n in (2, 4, 8, 16)}
    return out


def fmt(a):
    L = []
    m = a["meta"]
    L.append(f"solver={m.get('solver')}  nDOF={m.get('nDOF')}  steps={m.get('nSteps')}  "
             f"threads={m.get('threads')}  step={a['step_ms']:.1f} ms")
    L.append("")
    L.append("PHASES (% of step)")
    for k, v in sorted(a["phases"].items(), key=lambda kv: -kv[1]["ms"]):
        L.append(f"  {k:<16} {v['ms']:>10.1f} ms  {v['pct_step']:>6.2f}%")
    L.append("")
    L.append("THREADABLE LOOPS -- kernel vs scatter")
    L.append(f"  {'loop':<12} {'loop ms':>10} {'kernel ms':>10} {'scatter ms':>11} "
             f"{'%step loop':>11} {'%step kern':>11} {'%loop scat':>11}")
    for k in ("A_update", "B_tangent", "C_residual"):
        if k not in a["loops"]:
            continue
        d = a["loops"][k]
        L.append(f"  {k:<12} {d['loop_ms']:>10.1f} {d['kernel_ms']:>10.1f} "
                 f"{d['scatter_ms']:>11.1f} {d['pct_step_loop']:>10.2f}% "
                 f"{d['pct_step_kernel']:>10.2f}% {d['pct_loop_scatter']:>10.2f}%")
    t = a["totals"]
    L.append(f"  {'TOTAL':<12} {t['loop_ms']:>10.1f} {t['kernel_ms']:>10.1f} "
             f"{t['loop_ms']-t['kernel_ms']:>11.1f} {t['pct_step_loop']:>10.2f}% "
             f"{t['pct_step_kernel']:>10.2f}%")
    L.append("")
    L.append("PER-CLASSTAG us/ele (threading regression risk: compare to barrier cost)")
    for k in ("A_update", "B_tangent", "C_residual"):
        if k not in a["loops"]:
            continue
        for tag, d in a["loops"][k]["per_classtag"].items():
            L.append(f"  {k:<12} classTag={tag:<7} n={d['count']:<9} "
                     f"{d['wall_ms']:>9.1f} ms  {d['us_per_ele']:>8.2f} us/ele")
    L.append("")
    L.append("AMDAHL ceiling on the KERNEL fraction only "
             f"(p={t['pct_step_kernel']:.1f}%)")
    L.append("  " + "  ".join(f"{n}T={v:.2f}x" for n, v in a["amdahl_kernel_only"].items()))
    L.append("")
    L.append("LOOP SITES (a loop scope can appear in several phases)")
    for k in ("A_update", "B_tangent", "C_residual"):
        if k in a["loops"]:
            for site, ms in a["loops"][k]["sites"].items():
                L.append(f"  {k:<12} {ms:>9.1f} ms  {site}")
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    h5path = sys.argv[1]
    a = analyse(*collect(h5path))
    print(f"=== {os.path.basename(h5path)} ===")
    print(fmt(a))
    if "--json" in sys.argv:
        j = sys.argv[sys.argv.index("--json") + 1]
        with open(j, "w") as f:
            json.dump(a, f, indent=2)
        print(f"\n[json] {j}")
