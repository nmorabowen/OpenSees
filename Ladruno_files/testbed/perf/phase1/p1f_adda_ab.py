"""ADR-75 P1f — interleaved A/B for the `PARDISOGenLinSOE::addA` scatter fix.

WHY A DRIVER AND NOT JUST TWO BENCH RUNS. The first attempt measured a baseline,
rebuilt, and compared. It was worthless: the baseline came out 30.87 s where the
SAME configuration had measured 24.46-24.61 s earlier in the session, because CI
polling was running against the machine during the timing window. A ~25% drift
between sessions dwarfs the effect being measured.

Two binaries cannot be imported into one process (`import opensees` binds one
.pyd), so interleaving has to happen at the PROCESS level: this driver alternates
OLD/NEW/OLD/NEW subprocesses inside a single wall-clock window, and takes the
median of each. That is the same discipline laneB_solver_bench.py applies across
solvers, lifted one level up.

  OPS_PYD_OLD  directory holding the pre-fix opensees.pyd (+ its MKL DLLs)
  OPS_PYD_NEW  directory holding the post-fix build (defaults to dist/bin)

Correctness: the fix must be EXACT — same entries, same order, same summation —
so the tip displacements must agree bit-for-bit, not approximately.
"""
from __future__ import annotations
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
HERE = Path(__file__).resolve().parent
PY = sys.executable

ROOT = (r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
        r"\opensees-pr-strategy-9c17b8")
OLD = os.environ.get("OPS_PYD_OLD")
NEW = os.environ.get("OPS_PYD_NEW", os.path.join(ROOT, "dist", "bin"))
REPEATS = int(os.environ.get("OPS_BENCH_REPEATS", "5"))
THREADS = os.environ.get("MKL_NUM_THREADS", "4")
SIZES = [int(s) for s in os.environ.get("OPS_BENCH_SIZES", "20,25").split(",")]
CONFIGS = [s.strip() for s in
           os.environ.get("OPS_BENCH_CONFIGS", "direct,spd").split(",")]

SYS_ARGS = {
    "direct": ["Pardiso"],                       # -matrixType 0, the default path
    "spd":    ["Pardiso", "-matrixType", 1],     # half-storage uses the same scan
    "krylov": ["Pardiso", "-krylov", 6],
}

DRIVER = r'''
import os, sys, time
os.environ["MKL_NUM_THREADS"] = "{threads}"
os.environ["OMP_NUM_THREADS"] = "{threads}"
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
DIST = r"{dist}"
sys.path.insert(0, DIST); os.add_dll_directory(DIST)
import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), "WRONG PYD " + ops.__file__
E, nu = 200000.0, 0.3
Kb = E/(3.0*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu)); LEN=100.0
n = {n}; NS = 15
def nid(i,j,k): return 1 + i + (n+1)*(j + (n+1)*k)
ops.wipe(); ops.model("basic","-ndm",3,"-ndf",3)
ops.nDMaterial("LadrunoJ2",1,Kb,G,"-iso","voce",250.0,0.0,0.0,2000.0)
for k in range(n+1):
    for j in range(n+1):
        for i in range(n+1): ops.node(nid(i,j,k), i*LEN, j*LEN, k*LEN)
for j in range(n+1):
    for i in range(n+1): ops.fix(nid(i,j,0),1,1,1)
e=1
for k in range(n):
    for j in range(n):
        for i in range(n):
            ops.element("LadrunoBrick",e, nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),
                        nid(i,j+1,k), nid(i,j,k+1),nid(i+1,j,k+1),
                        nid(i+1,j+1,k+1),nid(i,j+1,k+1),1); e+=1
ops.timeSeries("Linear",1); ops.pattern("Plain",1,1)
for j in range(n+1):
    for i in range(n+1): ops.load(nid(i,j,n), 4.0e5, 0.0, -4.0e5)
ops.constraints("Plain"); ops.numberer("RCM")
ops.system(*{sysargs})
ops.test("NormDispIncr",1e-7,25); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0/NS); ops.analysis("Static")
t0=time.perf_counter()
for s in range(NS):
    if ops.analyze(1) != 0:
        print("ANALYZE-FAILED", s); raise SystemExit(1)
wall=time.perf_counter()-t0
print("RESULT wall=%.6f ux=%.17e" % (wall, ops.nodeDisp(nid(n//2,n//2,n),1)))
'''


def one(dist, n, cfg):
    src = DRIVER.format(dist=dist, n=n, sysargs=repr(SYS_ARGS[cfg]),
                        threads=THREADS)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
        f.write(src)
        p = f.name
    out = subprocess.run([PY, "-S", p], capture_output=True, text=True)
    os.unlink(p)
    m = re.search(r"RESULT wall=([\d.]+) ux=(\S+)", out.stdout + out.stderr)
    if not m:
        raise RuntimeError((out.stdout + out.stderr)[-500:])
    return float(m.group(1)), m.group(2)


def main():
    if not OLD or not os.path.isdir(OLD):
        sys.exit("set OPS_PYD_OLD to the directory holding the PRE-fix "
                 "opensees.pyd (copy the whole dist/bin, it needs the MKL DLLs)")
    print(f"P1f addA A/B | threads={THREADS} repeats={REPEATS} sizes={SIZES} "
          f"configs={CONFIGS}\n  OLD={OLD}\n  NEW={NEW}\n")
    out = {"threads": THREADS, "repeats": REPEATS, "old": OLD, "new": NEW,
           "sizes": {}}
    for n in SIZES:
        ndof = 3 * (n + 1) ** 2 * n
        print(f"=== grid {n}^3 (~{ndof/1000:.1f}k DOF) ===")
        for cfg in CONFIGS:
            so, sn, uo, un = [], [], None, None
            for _ in range(REPEATS):        # interleave OLD/NEW every round
                w, u = one(OLD, n, cfg); so.append(w); uo = u
                w, u = one(NEW, n, cfg); sn.append(w); un = u
            mo, mn = statistics.median(so), statistics.median(sn)
            exact = (uo == un)
            print(f"  {cfg:7s} old {mo:8.3f}s  new {mn:8.3f}s  "
                  f"speedup {mo/mn:5.3f}x  ux {'BIT-IDENTICAL' if exact else '*** DIFFERS ***'}")
            if not exact:
                print(f"          old ux={uo}\n          new ux={un}")
            out["sizes"].setdefault(str(n), {})[cfg] = {
                "old_median_s": mo, "new_median_s": mn, "speedup": mo / mn,
                "ux_old": uo, "ux_new": un, "bit_identical": exact,
                "old_samples": sorted(so), "new_samples": sorted(sn)}
    dest = HERE / f"p1f_adda_ab_t{THREADS}.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
