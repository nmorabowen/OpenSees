"""WP-187 timing bench: H5DRM read cost per step, before and after.

A synthetic .h5drm with the Paper 1 Legacy files' layout -- (3 * Nstations) x Nt
float64 datasets, CHUNKED (192, 125), DEFLATE level 4 -- drives a linear stdBrick
box whose every node is a station (outer shell external, the rest internal), so the
DRM layer is the two outer node shells (about 4 900 of 9 000 nodes at the default
30 x 30 x 10 grid; the 4D deck matches 5 091). The solve is factored once
(`algorithm Linear -factorOnce`), so a step is the H5DRM read plus the DRM boundary
forces plus one triangular solve; a run with no pattern gives that floor.

Arms (classic Tcl OpenSees.exe):
  base       : same model, no load pattern (the non-H5DRM floor)
  old        : --exe-old (a pre-WP-187 binary, e.g. the paper1-single 1a4141ee4 pin)
  legacy     : new binary, -prefetch 0 -chunkCacheMB 0 (the old read path in the new code)
  cache      : new binary, -prefetch 0 (sized chunk cache only)
  default    : new binary, no flags (sized cache + chunk-length time block)

Two station orders bracket the access pattern: "aligned" (data rows in node-tag
order, consecutive reads hit the same chunk) and "shuffled" (random rows, every
node in a different chunk row).

Every arm writes all nodal displacements of the last step with %.17g; the script
checks they are byte-identical to the old binary (or to legacy) and reports the
mean step time over the timed window.

Usage:
  python bench_h5drm_reads.py --exe-new <dist/bin/OpenSees.exe> [--exe-old <pin OpenSees.exe>]
         [--nx 30 --ny 30 --nz 10 --nt 1000 --steps 300] [--work DIR] [--json out.json]
"""
import argparse
import json
import re
import os
import subprocess
import time

import numpy as np
import h5py

CHUNK = (192, 125)
DT = 0.01
H = 0.01          # km
CRD = 1000.0      # km -> m


def write_file(path, nx, ny, nz, nt, order, seed=187):
    rng = np.random.default_rng(seed)
    n = nx * ny * nz
    st = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]   # node-tag order
    if order == "aligned":
        slot = np.arange(n)
    else:
        slot = rng.permutation(n)
    xyz = np.array([(i * H, j * H, k * H) for i, j, k in st])
    internal = np.array([0 < i < nx - 1 and 0 < j < ny - 1 and 0 < k < nz - 1 for i, j, k in st])
    t = np.arange(nt) * DT
    with h5py.File(path, "w") as f:
        d = f.create_group("DRM_Data")
        d["xyz"] = xyz
        d["internal"] = internal
        d["data_location"] = (3 * slot).astype(np.int32)
        dsd = d.create_dataset("displacement", shape=(3 * n, nt), dtype="f8", chunks=CHUNK,
                               compression="gzip", compression_opts=4)
        dsa = d.create_dataset("acceleration", shape=(3 * n, nt), dtype="f8", chunks=CHUNK,
                               compression="gzip", compression_opts=4)
        # smooth band-limited motion per station (compresses about 3:1, like the real files)
        nb = 4096
        for r0 in range(0, 3 * n, nb):
            r1 = min(3 * n, r0 + nb)
            m = r1 - r0
            w = 2 * np.pi * (0.5 + 4.0 * rng.random((m, 1)))
            ph = 2 * np.pi * rng.random((m, 1))
            amp = 1e-2 * (0.5 + rng.random((m, 1)))
            env = np.minimum(1.0, t / 1.0)[None, :]
            dsd[r0:r1, :] = amp * env * np.sin(w * t[None, :] + ph)
            dsa[r0:r1, :] = -amp * env * w * w * np.sin(w * t[None, :] + ph)
        mdat = f.create_group("DRM_Metadata")
        mdat["dt"] = DT
        mdat["tstart"] = 0.0
        mdat["tend"] = (nt - 1) * DT
        mdat["drmbox_x0"] = np.zeros(3)
        for kname, v in (("drmbox_xmin", 0.0), ("drmbox_xmax", (nx - 1) * H),
                         ("drmbox_ymin", 0.0), ("drmbox_ymax", (ny - 1) * H),
                         ("drmbox_zmin", 0.0), ("drmbox_zmax", (nz - 1) * H)):
            mdat[kname] = float(v)
        mdat["h"] = np.array([H, H, H])


def deck(path, nx, ny, nz, steps, flags, out, pattern=True, warm=5):
    tag = lambda i, j, k: 1 + (k * ny + j) * nx + i
    L = ["model basic -ndm 3 -ndf 3"]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                L.append(f"node {tag(i, j, k)} {i * H * CRD!r} {j * H * CRD!r} {k * H * CRD!r}")
    L.append("nDMaterial ElasticIsotropic 1 2.0e8 0.3 2.0")
    e = 1
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                ns = [tag(i, j, k), tag(i + 1, j, k), tag(i + 1, j + 1, k), tag(i, j + 1, k),
                      tag(i, j, k + 1), tag(i + 1, j, k + 1), tag(i + 1, j + 1, k + 1), tag(i, j + 1, k + 1)]
                L.append(f"element stdBrick {e} {' '.join(map(str, ns))} 1")
                e += 1
    for i in range(nx):
        for j in range(ny):
            L.append(f"fix {tag(i, j, 0)} 1 1 1")
    if pattern:
        L.append(f"pattern H5DRM 1 {{{path}}} 1.0 {CRD!r} 1.0e-3 1 1 0 0 0 1 0 0 0 1 0 0 0 {' '.join(map(str, flags))}")
    L += ["constraints Transformation", "numberer RCM", "system UmfPack",
          "test NormDispIncr 1e-10 10", "algorithm Linear -factorOnce", "integrator Newmark 0.5 0.25",
          "analysis Transient",
          f"analyze {warm} {DT!r}",
          "set t0 [clock microseconds]",
          f"analyze {steps} {DT!r}",
          "set t1 [clock microseconds]",
          f"set fo [open {{{out}}} w]",
          f"puts $fo [format \"STEP_US %.1f\" [expr {{($t1 - $t0) / double({steps})}}]]",
          "foreach t [lsort -integer [getNodeTags]] { foreach v [nodeDisp $t] { puts -nonewline $fo [format \"%.17g \" $v] } }",
          "puts $fo \"\"",
          "close $fo",
          "wipe"]   # the pattern's destructor prints the -verbose read timer
    return "\n".join(L) + "\n"


def run(exe, deck_path, out, work, tcl_library=None):
    env = dict(os.environ)
    if tcl_library:
        env["TCL_LIBRARY"] = tcl_library      # a pin bin/ without lib/tcl8.6 next to it
    t = time.perf_counter()
    r = subprocess.run([exe, deck_path], cwd=work, capture_output=True, text=True, env=env)
    wall = time.perf_counter() - t
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"{exe} failed:\n{r.stderr[-3000:]}")
    lines = open(out).read().splitlines()
    us = float(lines[0].split()[1])
    log = r.stdout + r.stderr
    m = re.search(r"WP-187 read time \S+ s over \d+ steps \((\S+) ms/step\)", log)
    read_ms = float(m.group(1)) if m else None
    return us, lines[1], wall, len(log), read_ms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe-new", required=True)
    ap.add_argument("--exe-old")
    ap.add_argument("--nx", type=int, default=30)
    ap.add_argument("--ny", type=int, default=30)
    ap.add_argument("--nz", type=int, default=10)
    ap.add_argument("--nt", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--orders", default="aligned,shuffled")
    ap.add_argument("--work", default=os.path.abspath("wp187_bench_work"))
    ap.add_argument("--json")
    ap.add_argument("--tcl-library", default=None,
                    help="TCL_LIBRARY for every run (default: <exe-new>/../lib/tcl8.6 when it exists)")
    a = ap.parse_args()
    if a.tcl_library is None:
        cand = os.path.normpath(os.path.join(os.path.dirname(a.exe_new), "..", "lib", "tcl8.6"))
        a.tcl_library = cand if os.path.isdir(cand) else None
    os.makedirs(a.work, exist_ok=True)
    res = {"grid": [a.nx, a.ny, a.nz], "nt": a.nt, "steps": a.steps, "chunk": CHUNK, "arms": {}}
    for order in a.orders.split(","):
        h5 = os.path.join(a.work, f"bench_{order}.h5drm").replace("\\", "/")
        if not os.path.exists(h5):
            t = time.perf_counter()
            write_file(h5, a.nx, a.ny, a.nz, a.nt, order)
            print(f"[{order}] wrote {h5} ({os.path.getsize(h5) / 2**20:.0f} MiB) in {time.perf_counter() - t:.0f} s", flush=True)
        # -verbose on the new binary: per-step lines (as the old binary always printed) and the read timer
        arms = [("base", a.exe_new, None), ("legacy", a.exe_new, ["-prefetch", 0, "-chunkCacheMB", 0, "-verbose"]),
                ("cache", a.exe_new, ["-prefetch", 0, "-verbose"]), ("default", a.exe_new, ["-verbose"]),
                ("default_quiet", a.exe_new, [])]
        if a.exe_old:
            arms.insert(1, ("old", a.exe_old, []))
        ref = None
        for name, exe, flags in arms:
            out = os.path.join(a.work, f"{order}_{name}.out").replace("\\", "/")
            dk = os.path.join(a.work, f"{order}_{name}.tcl")
            open(dk, "w").write(deck(h5, a.nx, a.ny, a.nz, a.steps, flags or [], out, pattern=(name != "base")))
            us, disp, wall, logb, read_ms = run(exe, dk, out, a.work, a.tcl_library)
            same = None
            if name != "base":
                if ref is None:
                    ref = disp
                same = (disp == ref)
            res["arms"][f"{order}/{name}"] = {"ms_per_step": us / 1000.0, "identical_to_first": same,
                                             "read_ms_per_step": read_ms, "wall_s": wall, "log_bytes": logb}
            print(f"[{order}] {name:13s} {us / 1000.0:9.2f} ms/step  read {read_ms} ms/step  identical={same}  "
                  f"log={logb} B  wall={wall:.0f} s", flush=True)
        base = res["arms"][f"{order}/base"]["ms_per_step"]
        for name in ("old", "legacy", "cache", "default", "default_quiet"):
            k = f"{order}/{name}"
            if k in res["arms"]:
                res["arms"][k]["h5drm_ms_per_step"] = res["arms"][k]["ms_per_step"] - base
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
