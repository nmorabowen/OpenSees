"""ADR-75 P1h -- turn a profiler h5 into the PARDISO solver-phase split.

Usage:  python p1h_parse_split.py <file.h5> [--json out.json]
        (needs h5py: C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe)

Why this exists (and why parse_lane3.py is not enough): parse_lane3 knows the three
UmfPack-era scope names (soe.symbolic / soe.factor / soe.trisolve) because those are all
that existed when it was written. PR #667 added four more brackets that are PARDISO-only
or new, and two of them change how the split must be READ:

  soe.cgs      phase 23, the `-krylov` CGS path. NOT a subset of soe.factor: when CGS
               gives up, PARDISO refactorizes INSIDE the phase-23 call (Intel's automatic
               fallback), so that factorization bills here. A -krylov run's soe.factor
               therefore UNDERSTATES total factorization work, and the two must be
               reported as separate columns, never summed into "factor".
  dc.s.fill    setSize CSR pattern build (sorted insertion) -- once per sparsity pattern,
  dc.s.verify  + the P1d symmetry-contract sweep. Not per-solve; do not divide by nSolve.
  soe.addA     DEEP-gated per-element scatter attribution. Absent from a coarse run BY
               DESIGN -- absence is the gating working, not a missing cost.

The multi-site rule from ADR-40b/75b applies to every one of them: a scope appears at
SEVERAL places in the tree (soe.factor under `linearSolve` AND under `update`/`newStep`
when the integrator solves too). Every occurrence is summed and the per-site breakdown is
kept -- that gap is the whole finding of the LEDGER_quirks row this parser services.
"""
import os, sys, json, collections

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py

# The solver-internal brackets. Order matters: it is the report order.
SOLVER_SCOPES = ("soe.symbolic", "soe.factor", "soe.trisolve", "soe.cgs")
SETSIZE_SCOPES = ("dc.s.fill", "dc.s.verify")
DEEP_SCOPES = ("soe.addA",)
# Coarse driver phases, for context and for the attribution gap.
PHASES = ("linearSolve", "formTangent", "formUnbalance", "update", "newStep",
          "commit", "convTest")


def collect(path):
    """Walk the rollup; return (meta, scopes) with scopes: full/path -> {wall_ms, calls}."""
    scopes = {}
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
                        scopes[pp] = {"wall_ms": int(a["wall_ns"]) / 1e6,
                                      "calls": int(a.get("calls", 0))}
                    walk(it, pp)

        if "rollup" in run and "root" in run["rollup"]:
            walk(run["rollup"]["root"], "root")
    return meta, scopes


def _sum_sites(scopes, name):
    """Sum a leaf scope name over EVERY site it appears at; keep the per-site rows."""
    sites = {k: v for k, v in scopes.items() if k.rsplit("/", 1)[-1] == name}
    return (sum(v["wall_ms"] for v in sites.values()),
            sum(v["calls"] for v in sites.values()),
            {k: round(v["wall_ms"], 3) for k, v in sorted(sites.items())})


def analyse(meta, scopes):
    step_ms = scopes.get("root/step", {}).get("wall_ms") or meta.get("wall_ms_total", 0.0)
    out = {"meta": {k: meta.get(k) for k in
                    ("solver", "algorithm", "integrator", "nDOF", "nSteps",
                     "threads", "wall_ms_total", "timestamp")},
           "step_ms": step_ms, "phases": {}, "solver": {}, "setSize": {}, "deep": {}}

    for name in PHASES:
        ms, calls, sites = _sum_sites(scopes, name)
        if ms:
            out["phases"][name] = {"ms": round(ms, 3), "calls": calls,
                                   "pct_step": round(100.0 * ms / step_ms, 3)}

    for group, names in (("solver", SOLVER_SCOPES), ("setSize", SETSIZE_SCOPES),
                         ("deep", DEEP_SCOPES)):
        for name in names:
            ms, calls, sites = _sum_sites(scopes, name)
            if not sites:
                continue
            out[group][name] = {"ms": round(ms, 3), "calls": calls,
                                "pct_step": round(100.0 * ms / step_ms, 3),
                                "sites": sites}

    # ---- the headline: what fraction of solver time is factorization? ----
    fac = out["solver"].get("soe.factor", {}).get("ms", 0.0)
    tri = out["solver"].get("soe.trisolve", {}).get("ms", 0.0)
    sym = out["solver"].get("soe.symbolic", {}).get("ms", 0.0)
    cgs = out["solver"].get("soe.cgs", {}).get("ms", 0.0)
    tot = fac + tri + sym + cgs
    ls = out["phases"].get("linearSolve", {}).get("ms", 0.0)
    out["split"] = {
        "solver_total_ms": round(tot, 3),
        "pct_step": round(100.0 * tot / step_ms, 3) if step_ms else 0.0,
        # Shares OF THE SOLVER, which is the number that bounds every reuse lever.
        "pct_solver_factor":   round(100.0 * fac / tot, 2) if tot else 0.0,
        "pct_solver_trisolve": round(100.0 * tri / tot, 2) if tot else 0.0,
        "pct_solver_symbolic": round(100.0 * sym / tot, 2) if tot else 0.0,
        "pct_solver_cgs":      round(100.0 * cgs / tot, 2) if tot else 0.0,
        # ADR-40b's trap: solver work booked OUTSIDE the linearSolve phase.
        "linearSolve_ms": round(ls, 3),
        "attribution_gap_ms": round(tot - ls, 3),
        "gap_note": ("solver total EXCEEDS linearSolve => the integrator solves too "
                     "(DisplacementControl/ArcLength shape)" if tot > ls else
                     "solver total within linearSolve (expected for LoadControl)"),
        "cgs_caveat": ("soe.cgs > 0: soe.factor UNDERSTATES factorization -- a CGS give-up "
                       "refactorizes inside phase 23 and bills to soe.cgs"
                       if cgs > 0 else None),
    }
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    path = sys.argv[1]
    meta, scopes = collect(path)
    out = analyse(meta, scopes)

    m, s = out["meta"], out["split"]
    print(f"\n=== {os.path.basename(path)} ===")
    print(f"solver={m.get('solver')} nDOF={m.get('nDOF')} steps={m.get('nSteps')} "
          f"step={out['step_ms']:.1f} ms")
    print(f"  (meta 'threads'={m.get('threads')} -- UNRELIABLE, see ADR-75 sec.9 open item; "
          f"record MKL_NUM_THREADS from the environment)")
    print(f"\n  {'scope':<16} {'ms':>10} {'calls':>7} {'%step':>7} {'%solver':>8}")
    for name in SOLVER_SCOPES:
        d = out["solver"].get(name)
        if not d:
            print(f"  {name:<16} {'--':>10} {'--':>7} {'--':>7} {'--':>8}   (absent)")
            continue
        pct_solver = 100.0 * d["ms"] / s["solver_total_ms"] if s["solver_total_ms"] else 0
        print(f"  {name:<16} {d['ms']:>10.2f} {d['calls']:>7} {d['pct_step']:>7.2f} "
              f"{pct_solver:>8.2f}")
    print(f"  {'SOLVER TOTAL':<16} {s['solver_total_ms']:>10.2f} {'':>7} "
          f"{s['pct_step']:>7.2f}")
    print(f"  {'linearSolve':<16} {s['linearSolve_ms']:>10.2f}")
    print(f"  attribution gap  {s['attribution_gap_ms']:>10.2f}  ({s['gap_note']})")
    if s["cgs_caveat"]:
        print(f"  !! {s['cgs_caveat']}")
    for group, label in (("setSize", "setSize"), ("deep", "deep")):
        for name, d in out[group].items():
            print(f"  {name:<16} {d['ms']:>10.2f} {d['calls']:>7} {d['pct_step']:>7.2f}"
                  f"   ({label})")

    if "--json" in sys.argv:
        dst = sys.argv[sys.argv.index("--json") + 1]
        with open(dst, "w") as f:
            json.dump(out, f, indent=2)
        print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
