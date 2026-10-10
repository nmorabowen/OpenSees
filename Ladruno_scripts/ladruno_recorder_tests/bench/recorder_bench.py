"""Ladruno recorder write-path benchmark (WP-164 step 0 and regression gate).

Three cases, each run in a fresh interpreter with the recorder ON and OFF, so the
recorder's own cost is (wall_on - wall_off) for the same analysis:

  small_explicit  2D quad strip, explicit central difference, many steps; the
                  recorder holds only SMALL channels (a 10-node region + -G energy).
                  Exposes per-step HDF5 overhead independent of model size (P2:
                  reopen + flush every step, recompress of a partial chunk).
  envelope_medium 3D brick block, static, `-envelope` on nodal + element stress.
                  Exposes the per-step delete/recreate of every envelope group (P1).
  large_slab      3D brick block, static, streaming element stress at the GPs.
                  Exposes per-step deflate of whole-slab chunks (P3/P4) and, on the
                  read side, the cost of ONE element's time history (P3).

Usage (the build python, matching the opensees.pyd ABI):
    py -3.12 recorder_bench.py <dist\\bin> <out_dir> [case ...] [--repeat N]
Prints one BENCH_ROW json line per case and a summary table; writes
<out_dir>/bench_results.json. Wall times on a shared, loaded box are noisy:
compare medians over --repeat runs, and compare two builds on the SAME box.
"""
from __future__ import annotations

import json
import os
import statistics
import subprocess
import sys
import time

CASES = {
    "small_explicit": dict(nx=40, steps=20000),
    "envelope_medium": dict(nx=24, ny=24, nz=8, steps=40),
    "large_slab": dict(nx=40, ny=40, nz=12, steps=20),
}

# extra recorder tokens for every case, e.g. --extra "-compress 1"
EXTRA = []

RECORDER_ARGS = {
    "small_explicit": ["-R", 1, "-N", "displacement", "-G", "energy"],
    "envelope_medium": ["-N", "displacement", "reactionForce", "-E", "stresses", "-envelope"],
    "large_slab": ["-N", "displacement", "-E", "stresses"],
}

CHILD = r'''
import os, sys, time, json
DIST, OUT, CASE, MODE = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
P = json.loads(sys.argv[5]); RARGS = json.loads(sys.argv[6])
os.add_dll_directory(DIST); sys.path.insert(0, DIST)
import opensees as ops

fname = os.path.join(OUT, f"bench_{CASE}.ladruno")
if MODE == "on" and os.path.exists(fname):
    os.remove(fname)

def brick_block(nx, ny, nz):
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e7, 0.2)
    def nid(i, j, k): return 1 + i + (nx + 1) * (j + (ny + 1) * k)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), float(i), float(j), float(k))
                if k == 0:
                    ops.fix(nid(i, j, k), 1, 1, 1)
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                ops.element("stdBrick", e, nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k),
                            nid(i, j+1, k), nid(i, j, k+1), nid(i+1, j, k+1),
                            nid(i+1, j+1, k+1), nid(i, j+1, k+1), 1)
    return [nid(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)]

ops.wipe()
if CASE == "small_explicit":
    nx = P["nx"]
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e6, 0.25, 1.0)
    for j in range(2):
        for i in range(nx + 1):
            ops.node(1 + i + j * (nx + 1), float(i), float(j))
            ops.mass(1 + i + j * (nx + 1), 1.0, 1.0)
    ops.fix(1, 1, 1); ops.fix(nx + 2, 1, 1)
    for i in range(1, nx + 1):
        ops.element("quad", i, i, i + 1, nx + 2 + i, nx + 1 + i, 1.0, "PlaneStrain", 1)
    ops.region(1, "-node", *[nx + 1 - k for k in range(10)])
    if MODE == "on":
        ops.recorder("ladruno", fname, *RARGS)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1); ops.load(nx + 1, 1.0, 0.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.algorithm("Linear")
    ops.integrator("CentralDifference"); ops.analysis("Transient")
    t0 = time.perf_counter(); c0 = time.process_time()
    rc = ops.analyze(P["steps"], 1.0e-3)
else:
    top = brick_block(P["nx"], P["ny"], P["nz"])
    if MODE == "on":
        ops.recorder("ladruno", fname, *RARGS)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, 1.0, 0.0, -1.0)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 5); ops.algorithm("Linear", "-factorOnce")
    ops.integrator("LoadControl", 1.0 / P["steps"]); ops.analysis("Static")
    t0 = time.perf_counter(); c0 = time.process_time()
    rc = ops.analyze(P["steps"])
if rc != 0:
    raise SystemExit(f"analyze failed rc={rc}")
wall = time.perf_counter() - t0
ops.wipe()   # recorder destructor: final flush/close
wall_total = time.perf_counter() - t0
cpu_total = time.process_time() - c0
size = os.path.getsize(fname) if (MODE == "on" and os.path.exists(fname)) else 0
print("BENCH_JSON " + json.dumps(dict(case=CASE, mode=MODE, wall=wall,
                                      wall_total=wall_total, cpu_total=cpu_total,
                                      bytes=size)))
'''


def run_child(dist, out, case, mode):
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1", HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([sys.executable, "-S", "-c", CHILD, dist, out, case, mode,
                        json.dumps(CASES[case]), json.dumps(RECORDER_ARGS[case] + EXTRA)],
                       env=env, capture_output=True, text=True, timeout=7200,
                       stdin=subprocess.DEVNULL)
    for line in r.stdout.splitlines():
        if line.startswith("BENCH_JSON "):
            return json.loads(line[len("BENCH_JSON "):])
    raise RuntimeError(f"{case}/{mode} failed rc={r.returncode}\n{r.stdout[-2000:]}\n"
                       f"{r.stderr[-2000:]}")


def read_one_history(path, case):
    """Read side (P3): wall time to read ONE entity's full time history."""
    import h5py
    with h5py.File(path, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
        if case == "large_slab":
            grp = f[f"{stage}/RESULTS/ON_ELEMENTS/stresses"]
            data = next(iter(grp.values()))["DATA"]
        else:
            data = f[f"{stage}/RESULTS/ON_NODES/DISPLACEMENT/DATA"]
        k = data.shape[1] // 2
        t0 = time.perf_counter()
        _ = data[:, k, :]
        return time.perf_counter() - t0, list(data.shape), list(data.chunks or [])


def main():
    argv = sys.argv[1:]
    repeat = 1
    if "--repeat" in argv:
        i = argv.index("--repeat")
        repeat = int(argv[i + 1])
        del argv[i:i + 2]
    if "--extra" in argv:
        i = argv.index("--extra")
        EXTRA.extend(int(t) if t.lstrip("-").isdigit() and not t.startswith("-") else t
                     for t in argv[i + 1].split())
        del argv[i:i + 2]
    dist, out = argv[0], argv[1]
    cases = argv[2:] or list(CASES)
    os.makedirs(out, exist_ok=True)
    results = []
    for case in cases:
        on, off = [], []
        for _ in range(repeat):
            off.append(run_child(dist, out, case, "off"))
            on.append(run_child(dist, out, case, "on"))
        row = dict(case=case, params=CASES[case],
                   wall_off=statistics.median(r["wall_total"] for r in off),
                   wall_on=statistics.median(r["wall_total"] for r in on),
                   cpu_off=statistics.median(r["cpu_total"] for r in off),
                   cpu_on=statistics.median(r["cpu_total"] for r in on),
                   extra=" ".join(str(t) for t in EXTRA),
                   bytes=on[-1]["bytes"])
        row["recorder_s"] = row["wall_on"] - row["wall_off"]
        row["recorder_cpu_s"] = row["cpu_on"] - row["cpu_off"]
        fpath = os.path.join(out, f"bench_{case}.ladruno")
        if case != "envelope_medium" and os.path.exists(fpath):
            row["read_one_s"], row["data_shape"], row["chunks"] = read_one_history(fpath, case)
        results.append(row)
        print("BENCH_ROW " + json.dumps(row), flush=True)
    with open(os.path.join(out, "bench_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    print("\ncase              wall_off   wall_on  recorder_s  rec_cpu_s      MB  read_one_s")
    for r in results:
        print(f"{r['case']:<16} {r['wall_off']:9.2f} {r['wall_on']:9.2f} {r['recorder_s']:11.2f}"
              f" {r['recorder_cpu_s']:10.2f} {r['bytes']/1e6:7.1f}"
              f"  {r.get('read_one_s', float('nan')):10.4f}")


if __name__ == "__main__":
    main()
