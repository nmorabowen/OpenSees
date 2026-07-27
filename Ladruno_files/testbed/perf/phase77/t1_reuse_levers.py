"""ADR-77 T1 — the reuse levers as SHIPPED, measured on the transient step.

T0b established the frame this stage is answering inside of: at the 4-thread
operating point, assembly is ~65% of the step and `linearSolve` only ~26%. The
textbook ranking of the Newton family assumes the opposite (that factorization
is what you are trying to skip), so it may invert here. L3-0b warned exactly
this: measure, do not assume.

Arms (all at identical convergence tolerance -- that is what makes the
comparison fair; only maxIter differs from T0, raised so the slow-converging
initial-stiffness arms can actually finish):

  Newton                              re-assemble + re-factor EVERY iteration
  Newton -initial                     cached Ki assembly, re-factor every iter
  ModifiedNewton                      1 assembly + 1 factor per STEP
  ModifiedNewton -initial             cached Ki, 1 assembly + 1 factor per step
  ModifiedNewton -initial -factoronce 1 factorization per ANALYSIS
  KrylovNewton                        accelerated modified Newton

WHY `-initial` IS A LEGITIMATE ARM ON THIS DECK (checked in source 2026-07-26,
and it is the ADR-76 §3 trap NOT firing for once):

  LadrunoJ2::getInitialTangent()   -> buildElasticTangent(), state-INDEPENDENT
  LadrunoBrick::getInitialStiff()  -> computes once, CACHES in Ki, returns *Ki

ADR-76 §3's warning is that `NDMaterial::getInitialTangent()` DEFAULTS to
`getTangent()`, so on most upstream solid models `-initial` is silently full
Newton at full cost. The FORK's own material and element both override it
properly, so here `-initial` is a real initial-stiffness iteration AND its
assembly is nearly free (a cached matrix copy plus the M/C adds). That is a
result worth recording on its own -- it says the ADR-76 trap is an upstream-
material property, not a fork-wide one.

AND IT MAKES `-factoronce` MATHEMATICALLY LEGITIMATE HERE, which is rare:
K_eff = c1*Ki + c2*C + c3*M with Ki constant, Rayleigh coefficients constant,
and dt constant is a genuinely constant matrix across the whole analysis. The
ADR-76 staleness trap needs the tangent to actually change; on this arm it
cannot. The answer check below is what proves the claim rather than assuming
it -- if `-factoronce` is silently wrong, uz will say so.

Run (operating point = 4 threads per T0b):
  OPS_PYD=<worktree>\\dist\\bin MKL_NUM_THREADS=4 python3.12 -S t1_reuse_levers.py

LOAD SCALE IS A FIRST-CLASS KNOB, and finding out why is the reason this stage
has one. MEASURED 2026-07-26 at the T0 deck's default load: peak von Mises is
252.5 MPa against s0 = 250 -- the deck is barely 1% past yield, a hairline
plastic zone. On such a deck K_t ~= K_i, so `ModifiedNewton` and
`ModifiedNewton -initial` return IDENTICAL iteration counts (61 at n=6) and the
initial-stiffness arms look far better than they would on a genuinely nonlinear
model. T0/T0b were pure timing anatomy and were unaffected by this; T1 is about
ITERATION COUNTS, so it is not. Load x4 reaches 611 MPa (2.4x yield) and x10
reaches 1773 MPa, so the sweep spans hairline -> solidly -> heavily plastic and
the deliverable becomes "which lever at what nonlinearity" rather than one
number that silently only holds near the elastic limit.

Env knobs: OPS_T1_N, OPS_T1_STEPS, OPS_T1_REPEATS, OPS_T1_SYSTEM, OPS_T1_TAG,
           OPS_T1_ARMS (comma list of arm keys), OPS_T1_MAXITER,
           OPS_T1_LOAD (load multiplier, default 1.0)
"""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", os.environ["MKL_NUM_THREADS"])
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees"
    r"\.claude\worktrees\pardisio-profiling-0a03b1\dist\bin",
)
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

HERE = Path(__file__).resolve().parent
N = int(os.environ.get("OPS_T1_N", "15"))
NSTEPS = int(os.environ.get("OPS_T1_STEPS", "6"))
REPEATS = int(os.environ.get("OPS_T1_REPEATS", "2"))
SYSTEM = os.environ.get("OPS_T1_SYSTEM", "Pardiso")
MAXITER = int(os.environ.get("OPS_T1_MAXITER", "300"))
THREADS = os.environ["MKL_NUM_THREADS"]
TAG = os.environ.get("OPS_T1_TAG", "")

# identical deck to T0 -- same constants, so T0/T0b numbers stay comparable
E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
RHO = 7.85e-9
L = 100.0
DT = 2.0e-4
TOL = 1.0e-7
LOAD = float(os.environ.get("OPS_T1_LOAD", "1.0"))

ARMS = {
    "newton":            ("Newton",         []),
    "newton_init":       ("Newton",         ["-initial"]),
    "modnewton":         ("ModifiedNewton", []),
    "modnewton_init":    ("ModifiedNewton", ["-initial"]),
    "modnewton_init_fo": ("ModifiedNewton", ["-initial", "-factoronce"]),
    "krylov":            ("KrylovNewton",   []),
}


def build(n, algo_name, algo_opts):
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso,
                   "-rho", RHO)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", e,
                    nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                    nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                    1)
                e += 1
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), 4.0e5 * LOAD, 0.0, -4.0e5 * LOAD)

    ops.constraints("Plain")
    ops.numberer("RCM")          # serial: the Ladruno numberer is parallel-only (T0)
    ops.system(SYSTEM)
    ops.test("NormDispIncr", TOL, MAXITER)
    ops.algorithm(algo_name, *algo_opts)
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ndof = (nnode - (nx + 1) * (ny + 1)) * 3
    return nid(nx // 2, ny // 2, nz), e - 1, ndof


def run_once(n, algo_name, algo_opts, profile_tag=None):
    tip, nele, ndof = build(n, algo_name, algo_opts)
    if profile_tag:
        ops.profiler("start", "-deep", "-perStep")
    iters = []
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1, DT) != 0:
            raise RuntimeError(f"analyze FAILED at step {s+1}")
        iters.append(int(ops.numIter()))
    wall = time.perf_counter() - t0
    if profile_tag:
        ops.profiler("stop")
        ops.profiler("report", str(HERE / f"t1_{profile_tag}_n{n}{TAG}.h5"),
                     "-run", profile_tag)
        ops.profiler("reset")
    return wall, ops.nodeDisp(tip, 3), iters, nele, ndof


def main():
    keys = [k.strip() for k in os.environ.get("OPS_T1_ARMS", ",".join(ARMS)).split(",")
            if k.strip()]
    print(f"ADR-77 T1 reuse levers | n={N} steps={NSTEPS} repeats={REPEATS} "
          f"system={SYSTEM} tol={TOL:g} maxIter={MAXITER} threads={THREADS} "
          f"loadx{LOAD:g}")
    out = {"n": N, "steps": NSTEPS, "system": SYSTEM, "threads": THREADS,
           "tol": TOL, "maxIter": MAXITER, "load_scale": LOAD, "arms": {}}
    ref = None
    for key in keys:
        name, opts = ARMS[key]
        label = " ".join([name] + opts)
        try:
            samples, uz, iters = [], None, None
            for _ in range(REPEATS):
                w, uz, iters, nele, ndof = run_once(N, name, opts)
                samples.append(w)
            med = statistics.median(samples)
            run_once(N, name, opts, profile_tag=key)
        except RuntimeError as exc:
            print(f"  {label:<38} DID NOT CONVERGE: {exc}")
            out["arms"][key] = {"label": label, "failed": str(exc)}
            continue
        if ref is None:
            ref = uz
        rel = abs(uz - ref) / abs(ref) if ref else 0.0
        tot = sum(iters)
        # The answer check is the whole point for -factoronce: a stale-tangent
        # bug shows up HERE, not in the wall clock.
        ok = "OK" if rel <= 1e-6 else "** ANSWER DIFFERS **"
        out["arms"][key] = {"label": label, "wall_s": med, "samples": sorted(samples),
                            "uz": uz, "rel_err_vs_newton": rel, "iters": iters,
                            "total_iters": tot, "iters_per_step": tot / len(iters),
                            "ms_per_iter": med * 1e3 / tot, "answer": ok}
        print(f"  {label:<38} wall={med:7.3f}s  iters={tot:>4} "
              f"({tot/len(iters):5.2f}/step)  ms/iter={med*1e3/tot:7.1f}  "
              f"relerr={rel:.1e} {ok}")
    if "newton" in out["arms"] and out["arms"]["newton"].get("wall_s"):
        base = out["arms"]["newton"]["wall_s"]
        print("\n  speedup vs Newton:")
        for key, a in out["arms"].items():
            if a.get("wall_s"):
                print(f"    {a['label']:<38} {base / a['wall_s']:5.2f}x")
    (HERE / f"t1_reuse_levers{TAG}.json").write_text(json.dumps(out, indent=2))
    print(f"wrote t1_reuse_levers{TAG}.json + per-arm .h5 profiles")


if __name__ == "__main__":
    main()
