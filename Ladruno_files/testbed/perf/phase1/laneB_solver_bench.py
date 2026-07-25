"""ADR-75 P1 gate — Lane-B desktop sparse-direct solver bench.

Establishes the **UmfPack baseline** (and serial sparse-direct comparators) that
a wired `system Pardiso` must beat, on the SAME Lane-B model the ADR-40b
dominance report used (3D LadrunoBrick + LadrunoJ2, ~11.5k DOF, implicit static,
15 steps). PARDISO/MUMPS are NOT in the serial build today (both unwired — the
ADR-75 finding), so their rows record `unavailable` until P1/P2 land; the harness
is ready to include them the moment they register.

Metric: wall seconds for the full 15-step analyze() (build excluded), MEDIAN of
REPEATS interleaved rounds (A/B/A/B... to blunt background contamination, per
PHASE0_RECIPE.md). Threads pinned to 1 for the baseline — the PARDISO thread
sweep (1/2/4/8) is a separate driver added with P1.

Correctness cross-check: every solver must reproduce the same tip displacement
to 1e-9 relative, else its timing is meaningless.

Run (this machine):
  C:\\Users\\nmb\\AppData\\Local\\Python\\bin\\python3.12 -S laneB_solver_bench.py
Optional: set OPS_PYD to override the pyd dir; MKL_NUM_THREADS/OMP_NUM_THREADS.
"""
from __future__ import annotations
import json
import os
import statistics
import sys
import time
from pathlib import Path

# stream progress even when piped (Python block-buffers a non-TTY stdout)
sys.stdout.reconfigure(line_buffering=True)

# --- threads pinned BEFORE importing the solver (baseline policy) ----------
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("PMI_RANK", "0")

# --- pyd wiring (this machine; override with OPS_PYD) -----------------------
DIST = os.environ.get("OPS_PYD", r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__)

HERE = Path(__file__).resolve().parent
REPEATS = 7
NSTEPS = 15

# Serial sparse-direct comparators available today + the P1/P2 targets (expected
# unavailable in the serial build until wired). "args" are extra system() tokens.
SOLVERS = [
    ("UmfPack",  ["UmfPack"]),                 # baseline
    ("SuperLU",  ["SuperLU"]),
    ("SparseSYM", ["SparseSYM"]),
    ("Mumps",    ["Mumps", "-ICNTL14", "200"]),  # P2 — unwired in serial today
    ("Pardiso",  ["Pardiso"]),                   # P1 — unwired today
]

# --- Lane-B model (verbatim geometry from phase0/laneB_model.py) ------------
E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
nx = ny = nz = 15
L = 100.0


def _nid(i, j, k):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def build(system_args):
    """Build Lane-B with the given `system` tokens. Returns tip node tag."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(_nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    eTag = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", eTag,
                    _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1), _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                    1)
                eTag += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(_nid(i, j, nz), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    return _nid(nx // 2, ny // 2, nz)


def run_once(system_args):
    """One timed 15-step solve. Returns (wall_s, tip_ux) or raises."""
    tip = build(system_args)
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"analyze failed at step {s + 1}")
    wall = time.perf_counter() - t0
    ux = ops.nodeDisp(tip, 1)
    if not (ux == ux) or abs(ux) == float("inf"):  # NaN/inf guard
        raise RuntimeError("non-finite tip displacement")
    return wall, ux


def probe(system_args):
    """Is this solver usable in this build? One trial run."""
    try:
        return True, run_once(system_args)
    except Exception as e:  # noqa: BLE001
        return False, str(e)


def main():
    print(f"Lane-B solver bench | REPEATS={REPEATS} NSTEPS={NSTEPS} "
          f"MKL_NUM_THREADS={os.environ['MKL_NUM_THREADS']}")
    avail = {}
    for name, args in SOLVERS:
        ok, info = probe(args)
        avail[name] = (args, ok, info)
        print(f"  probe {name:9s}: {'available' if ok else 'UNAVAILABLE — ' + str(info)[:70]}")

    live = [(n, a) for n, (a, ok, _) in avail.items() if ok]
    samples = {n: [] for n, _ in live}
    tips = {}
    for r in range(REPEATS):                       # interleaved rounds
        for name, args in live:
            wall, ux = run_once(args)
            samples[name].append(wall)
            tips[name] = ux

    base = tips.get("UmfPack")
    results = {}
    print(f"\n{'solver':10s} {'median_s':>10s} {'vs UmfPack':>11s}  tip_ux (rel err)")
    for name, _ in live:
        med = statistics.median(samples[name])
        rel = (med / statistics.median(samples["UmfPack"])) if "UmfPack" in samples else float("nan")
        err = abs(tips[name] - base) / abs(base) if base else float("nan")
        results[name] = {"median_s": med, "ratio_vs_umfpack": rel, "tip_ux": tips[name],
                         "tip_rel_err": err, "samples": sorted(samples[name])}
        print(f"{name:10s} {med:10.4f} {rel:10.2f}x  {tips[name]:.9e} ({err:.1e})")

    for name, (args, ok, info) in avail.items():
        if not ok:
            results[name] = {"unavailable": True, "reason": str(info)[:200], "args": args}

    out = HERE / "laneB_solver_baseline.json"
    out.write_text(json.dumps({
        "model": "laneB 15^3 LadrunoBrick+LadrunoJ2 ~11.5k DOF",
        "repeats": REPEATS, "nsteps": NSTEPS,
        "mkl_threads": os.environ["MKL_NUM_THREADS"],
        "pyd": ops.__file__, "results": results,
    }, indent=2))
    print(f"\nwrote {out}")
    print("NOTE: Pardiso/Mumps unavailable = expected (unwired in serial build; "
          "ADR-75). Re-run after P1/P2 to fill those rows + add the thread sweep.")


if __name__ == "__main__":
    main()
