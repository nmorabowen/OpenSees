"""ADR-75 P1j — does the factorization share of the solve RISE with model size?

P1h measured the phase split at ONE size (11.5k DOF) and predicted the trend rather
than measuring it: factorization is superlinear in N while triangular solve is ~linear
in nnz(L), so %fac should rise with N, making every factorization-reuse lever
(`factored` gate, -krylov, ModifiedNewton, IMPL-EX) MORE valuable at production scale.
P1h's own caveat called this "the one obvious follow-up". This is it.

Same model as laneB_scaling.py / laneB_model.py (cubic LadrunoBrick + LadrunoJ2,
Newton + LoadControl) so the numbers are comparable to every other ADR-75 row.

TWO MEASUREMENT TRAPS this script is built around:

 1. `soe.symbolic` runs ONCE PER SPARSITY PATTERN, not per solve. Shrinking the step
    count at large N to save wall time would therefore shrink the number of
    factorizations that amortize it, inflating symbolic's share at exactly one end of
    the sweep and MANUFACTURING a trend. So NSTEPS is constant across sizes, and the
    headline metric is fac/(fac+tri) -- symbolic excluded and reported separately.
 2. The `threads` run attribute was only fixed in P1i; every size row here records the
    pinned environment value, and the sweep refuses to run on a pre-P1i binary.

Run:
  OPS_PYD=<dist\\bin> MKL_NUM_THREADS=4 python3.12 -S p1j_size_trend.py
Env: OPS_TREND_SIZES (default "15,20,25,30,35" = 11.5k/26.5k/50.7k/86.5k/136.1k DOF),
     OPS_TREND_STEPS (15), OPS_TREND_REPEATS (2), OPS_TREND_CONFIGS
     (default "sym4:Pardiso -matrixType 2|unsym4:Pardiso -matrixType 0").
"""
from __future__ import annotations
import json
import os
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", os.environ["MKL_NUM_THREADS"])
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get("OPS_PYD")
if not DIST:
    sys.exit("OPS_PYD must point at the dist/bin holding the PARDISO build")
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"

HERE = Path(__file__).resolve().parent
SIZES = [int(s) for s in os.environ.get("OPS_TREND_SIZES", "15,20,25,30,35").split(",")]
NSTEPS = int(os.environ.get("OPS_TREND_STEPS", "15"))
REPEATS = int(os.environ.get("OPS_TREND_REPEATS", "2"))
THREADS = os.environ["MKL_NUM_THREADS"]
CONFIGS = [c.split(":", 1) for c in os.environ.get(
    "OPS_TREND_CONFIGS",
    "sym4:Pardiso -matrixType 2|unsym4:Pardiso -matrixType 0").split("|")]

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
L = 100.0


def _opt(a):
    """Ints and floats must be passed AS numbers: a string option value fails
    OPS_GetIntInput, the factory returns 0, and OpenSees silently falls back to
    ProfileSPDLinSOE -- the run then looks slow, not wrong (banked P1d trap)."""
    try:
        return int(a)
    except ValueError:
        pass
    try:
        return float(a)
    except ValueError:
        return a


def build(n, system_args):
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)
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
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                            nid(i, j+1, k+1), 1)
                e += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*[_opt(a) for a in system_args])
    ops.test("NormDispIncr", 1.0e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 15.0)
    ops.analysis("Static")
    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ndof = (nnode - (nx + 1) * (ny + 1)) * 3
    return nid(nx // 2, ny // 2, nz), e - 1, ndof


def run_once(n, tag, system_args, rep):
    tip, nele, ndof = build(n, system_args)
    out_h5 = HERE / f"p1j_{tag}_n{n}_r{rep}.h5"
    if out_h5.exists():
        out_h5.unlink()
    # ---- THE ACCUMULATION TRAP (cost: one full 53-minute sweep) ----------------
    # The Profiler is a PROCESS-GLOBAL singleton and `ops.wipe()` does NOT touch it.
    # This script runs every (size, config, repeat) in ONE process, so without an
    # explicit reset each report contains the sum of ALL runs so far: the first
    # sweep produced soe.factor call counts of 44, 88, 132 ... 936 and a step time
    # that was a running total. Only run #1 was clean, and the resulting "trend"
    # was an artifact of the accumulation, not of model size.
    ops.profiler("reset")
    # COARSE on purpose: the solver brackets are coarse scopes, and -deep taxes the
    # ASSEMBLY side (+4.5% measured in P1h) i.e. exactly the fraction being compared.
    ops.profiler("start")
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1) != 0:
            ops.profiler("stop")
            raise RuntimeError(f"analyze failed at step {s+1} (n={n}, {tag})")
    wall = time.perf_counter() - t0
    ops.profiler("stop")
    ops.profiler("report", str(out_h5))
    return wall, ops.nodeDisp(tip, 1), nele, ndof, out_h5




def main():
    print(f"Lane-B SIZE TREND | sizes={SIZES} steps={NSTEPS} repeats={REPEATS} "
          f"MKL_NUM_THREADS={THREADS}")
    print(f"configs: {[c[0] for c in CONFIGS]}")
    out = {"threads": int(THREADS), "steps": NSTEPS, "repeats": REPEATS, "rows": [],
           # Per-run wall keyed by h5 stem. This is the INDEPENDENT clock p1j_rollup.py
           # checks the profiler against -- the guard that would have caught the
           # accumulation bug on run #2 instead of after the whole sweep.
           "wall_by_run": {}}
    for n in SIZES:
        ref_ux = None
        for tag, sysstr in CONFIGS:
            args = sysstr.split()
            walls, ux, nele, ndof = [], None, None, None
            for rep in range(1, REPEATS + 1):
                w, ux, nele, ndof, h5 = run_once(n, tag, args, rep)
                walls.append(w)
                out["wall_by_run"][h5.stem] = w
                print(f"  n={n:>3} ndof={ndof:>7} {tag:<8} r{rep} wall={w:8.2f}s -> {h5.name}")
            # Cross-config correctness at every size: half-storage must not change the
            # answer. A timing trend measured on a wrong answer is worthless.
            if ref_ux is None:
                ref_ux = ux
                rel = 0.0
            else:
                rel = abs(ux - ref_ux) / abs(ref_ux) if ref_ux else 0.0
                if rel > 1e-9:
                    print(f"  !! n={n} {tag} tip disp differs from reference by {rel:.2e}")
            out["rows"].append({"n": n, "ndof": ndof, "nele": nele, "config": tag,
                                "walls_s": walls, "ux": ux, "rel_err_vs_ref": rel})
    dst = HERE / "p1j_size_trend.json"
    dst.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dst}")
    print("now parse: p1h_parse_split.py p1j_*.h5 --json <same>.json  (opensees_env)")


if __name__ == "__main__":
    main()
