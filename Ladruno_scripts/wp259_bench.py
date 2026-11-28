"""WP-259 bench: a stdBrick concrete cantilever at 1 / 4 / 8 threads (reported, no target).

The deck is the gate's `3d` cantilever (tests/test_wp259_concrete_threads.py) refined to
NX x NY x NZ = 40 x 8 x 12 (3840 stdBrick, 30 720 material points), pushed past cracking with
DisplacementControl + FixedNumIter. Two algorithms, so the element share is visible:

  * `mnewton`: ModifiedNewton -initial with UmfPack. The factorization happens once; every
    iteration is the threaded update + residual loops and one triangular solve, so this row
    is close to the element-loop scaling itself.
  * `newton`:  full Newton with `system Pardiso` (unsymmetric; MKL threads as configured).
    The threaded tangent loop is in it, and so is a factorization per iteration.

Each (material, algorithm, threads) runs in its own child process, R repeats, median wall of
the push only (model build excluded). The 1-thread and n-thread nodal displacements are also
compared (digest) so a fast-but-wrong row cannot pass unnoticed.

    set PYTHONPATH=<worktree>\\dist\\bin
    python3.12 Ladruno_scripts/wp259_bench.py --repeats 3 [--out bench.json]
"""
import argparse
import hashlib
import json
import os
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
sys.path.insert(0, str(TESTS))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

from _testbed.subprocess_run import run_python_script   # noqa: E402
import test_wp259_concrete_threads as G                   # noqa: E402

_CHILD = r'''
import json, sys, time
sys.path.insert(0, sys.argv[1])
from _testbed import ops
cfg = json.loads(sys.argv[2])
mat, algo, nthreads = cfg["mat"], cfg["algo"], cfg["threads"]
L, Bw, Hd, NX, NY, NZ = 400.0, 100.0, 100.0, 40, 8, 12
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
dx, dy, dz = L / NX, Bw / NY, Hd / NZ
nid = lambda i, j, k: 1 + i + (NX + 1) * (j + (NY + 1) * k)
for k in range(NZ + 1):
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, j, k), i * dx, j * dy, k * dz)
ops.nDMaterial(mat[0], 1, *mat[1:])
e = 0
for k in range(NZ):
    for j in range(NY):
        for i in range(NX):
            e += 1
            ops.element("stdBrick", e, nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k),
                        nid(i, j + 1, k), nid(i, j, k + 1), nid(i + 1, j, k + 1),
                        nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1), 1)
for k in range(NZ + 1):
    for j in range(NY + 1):
        ops.fix(nid(0, j, k), 1, 1, 1)
master = nid(NX, 0, NZ)
for k in range(NZ + 1):
    for j in range(NY + 1):
        if nid(NX, j, k) != master:
            ops.equalDOF(master, nid(NX, j, k), 3)
ops.ladrunoThreads(nthreads)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(master, 0.0, 0.0, -1.0)
ops.constraints("Transformation")
ops.numberer("RCM")
if algo == "mnewton":
    ops.system("UmfPack")
    ops.algorithm("ModifiedNewton", "-initial")
    iters = 20
else:
    ops.system("Pardiso")
    ops.algorithm("Newton")
    iters = 6
ops.test("FixedNumIter", iters, 0)
ops.integrator("DisplacementControl", master, 3, -0.03)
ops.analysis("Static")
t0 = time.perf_counter()
rcs = [ops.analyze(1) for _ in range(cfg["steps"])]
wall = time.perf_counter() - t0
disp = [x.hex() for n in range(1, nid(NX, NY, NZ) + 1) for x in ops.nodeDisp(n)]
ops.ladrunoThreads(1)
print("RESULT " + json.dumps({"build": ops.ladrunoBuild(), "rcs": rcs, "wall": wall,
                              "disp": disp if algo == "mnewton" else None}))
'''


def run(mat, algo, n, steps):
    cfg = {"mat": G.MATERIALS[mat][0], "algo": algo, "threads": n, "steps": steps}
    code, out = run_python_script(_CHILD, argv=[str(TESTS), json.dumps(cfg)], timeout=3600)
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    if code != 0 or not res:
        raise RuntimeError(out[-3000:])
    return json.loads(res[-1][len("RESULT "):])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--steps", type=int, default=6)
    ap.add_argument("--threads", default="1,4,8")
    ap.add_argument("--mats", default="c3d,rc")
    ap.add_argument("--algos", default="mnewton,newton")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    threads = [int(t) for t in a.threads.split(",")]
    rows = []
    for mat in a.mats.split(","):
        for algo in a.algos.split(","):
            ref_digest, t1 = None, None
            for n in threads:
                walls, digest, build, rcs = [], None, None, None
                for _ in range(a.repeats):
                    r = run(mat, algo, n, a.steps)
                    walls.append(r["wall"])
                    build, rcs = r["build"], r["rcs"]
                    if r["disp"] is not None:
                        digest = hashlib.sha256("".join(r["disp"]).encode()).hexdigest()[:16]
                med = statistics.median(walls)
                if n == threads[0]:
                    t1, ref_digest = med, digest
                row = {"mat": mat, "algo": algo, "threads": n, "median_s": round(med, 3),
                       "spread": round((max(walls) - min(walls)) / med, 3),
                       "speedup": round(t1 / med, 2), "rcs": rcs, "build": build,
                       "digest_equal_1t": None if digest is None else digest == ref_digest}
                rows.append(row)
                print(json.dumps(row), flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
