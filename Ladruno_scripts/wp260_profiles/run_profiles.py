"""WP-260 Wave 0 (`260-w0`) profiler runner — parent/child HDF5-isolated harness.

Modes per deck x size (`wp-briefs.md` §2):
  - `1t_deep`   : 1 thread, `-deep -perStep` (phase shares; deep gate refuses threading)
  - `1t_coarse` : 1 thread, `-perStep`       (serial wall time; deep overhead = 1t_deep/1t_coarse)
  - `nt_coarse` : N threads (default 8 on P-cores), `-perStep` (production wall & scaling)

Usage:
    py -3.12 Ladruno_scripts/wp260_profiles/run_profiles.py --smoke --verify-d2-linear
    py -3.12 Ladruno_scripts/wp260_profiles/run_profiles.py --full --export-dir results/<sha9>
"""
from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# NOTE: NEVER import h5py or numpy at module top level — child processes run `_run_child`
# and import `opensees.pyd`, which bundles its own HDF5 DLL. Import h5py/numpy only inside
# parent-only functions (`_ensure_h5drm`, `_extract_h5`).

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PYD = ROOT / "dist" / "bin" / "opensees.pyd"
DEFAULT_EXE = ROOT / "dist" / "bin" / "OpenSees.exe"
DEFAULT_OUT = Path(tempfile.gettempdir()) / "wp260w0_prof"

ENV_PREFIXES = ("OMP_", "KMP_", "MKL_", "LADRUNO_")


def _collect_env(env_map: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    src = env_map if env_map is not None else os.environ
    return {k: src[k] for k in sorted(src) if k.startswith(ENV_PREFIXES)}


def _apply_and_read_affinity(mask: Optional[int]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"requested": hex(mask) if mask is not None else None,
                           "set_ok": None, "process_mask": None, "system_mask": None}
    if sys.platform != "win32":
        return out
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    h = k32.GetCurrentProcess()
    if mask is not None:
        k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        ok = bool(k32.SetProcessAffinityMask(h, ctypes.c_size_t(mask)))
        out["set_ok"] = ok
        if not ok:
            raise OSError(f"SetProcessAffinityMask({hex(mask)}) failed: winerr={ctypes.get_last_error()}")
    pm, sm = ctypes.c_size_t(0), ctypes.c_size_t(0)
    k32.GetProcessAffinityMask.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)
    ]
    if k32.GetProcessAffinityMask(h, ctypes.byref(pm), ctypes.byref(sm)):
        out["process_mask"] = hex(pm.value)
        out["system_mask"] = hex(sm.value)
    return out


def _probe_exe_build(exe_path: Path) -> str:
    tcl = "puts [ladrunoBuild]; exit\n"
    cp = subprocess.run([str(exe_path)], input=tcl, text=True, capture_output=True, check=True)
    combined = cp.stdout + "\n" + cp.stderr
    for ln in reversed(combined.splitlines()):
        s = ln.strip()
        if s and len(s) >= 7 and all(c in "0123456789abcdefABCDEF-dirty" for c in s):
            return s
    return ""


def _parse_threading(threads: int, stderr_txt: str) -> Tuple[str, List[str]]:
    lines = [ln.strip() for ln in stderr_txt.splitlines() if "ladrunoThreads" in ln]
    if threads <= 1:
        return "no (1T)", lines
    joined = "\n".join(lines)
    if "THREADED" in joined and "PARTLY" not in joined and "NOT" not in joined and "SERIAL" not in joined:
        return f"yes ({threads}T)", lines
    if "PARTLY" in joined:
        return f"partial ({threads}T)", lines
    if "NOT" in joined or "SERIAL" in joined or "DEEP profiler gate" in joined:
        return "no (refused -> serial)", lines
    return f"MKL-only/unannounced ({threads}T)", lines


# --------------------------------------------------------------------------- #
#  Parent-side H5DRM generator + HDF5 tree/series reader (h5py parent ONLY)   #
# --------------------------------------------------------------------------- #
def _ensure_h5drm(out_dir: Path, size: str) -> Path:
    import h5py
    import numpy as np

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import decks  # noqa: E402

    p = d2_h5drm_path = out_dir / f"d2_wave_{size}.h5drm"
    cfg = decks.d2_h5drm_params(size)
    dx, nx, ny, nz = cfg["dx"], cfg["nx"], cfg["ny"], cfg["nz"]
    n_pml = cfg["n_pml"]
    # Emit stations on the 2-node shell surrounding the interior/PML boundary
    # Inner shell: i in [n_pml, nx - n_pml], j in [n_pml, ny - n_pml], k in [n_pml, nz]
    # Outer shell adjacent: one step outside (excluding top z=0)
    coords = []
    internal = []
    for iz in range(max(0, n_pml - 1), nz + 1):
        for iy in range(max(0, n_pml - 1), min(ny + 1, ny - n_pml + 2)):
            for ix in range(max(0, n_pml - 1), min(nx + 1, nx - n_pml + 2)):
                is_int = (n_pml <= ix <= nx - n_pml and
                          n_pml <= iy <= ny - n_pml and
                          iz >= n_pml)
                in_core = (n_pml + 1 <= ix <= nx - n_pml - 1 and
                           n_pml + 1 <= iy <= ny - n_pml - 1 and
                           iz >= n_pml + 1)
                if in_core:
                    continue  # interior points well inside need zero DRM load
                coords.append((cfg["x0"] + ix * dx, cfg["y0"] + iy * dx, cfg["z0"] + iz * dx))
                internal.append(is_int)

    xyz = np.array(coords, dtype=np.float64)
    internal_arr = np.array(internal, dtype=bool)
    n = len(xyz)
    dt, nt = cfg["dt"], cfg["nt"]
    t = np.arange(nt, dtype=np.float64) * dt

    # Smooth vertically-propagating SH wavefield u_x(z, t)
    cs = math.sqrt(2.5e8 / (2.0 * (1.0 + 0.25) * 2000.0))
    dis = np.zeros((3 * n, nt), dtype=np.float64)
    vel = np.zeros((3 * n, nt), dtype=np.float64)
    acc = np.zeros((3 * n, nt), dtype=np.float64)
    w = 2.0 * math.pi * 60.0
    amp = 1.0e-3
    for idx in range(n):
        tau = t - (xyz[idx, 2] - cfg["zmin"]) / cs
        mask = tau > 0.0
        s = np.zeros(nt)
        c = np.zeros(nt)
        s[mask] = np.sin(w * tau[mask]) * np.exp(-80.0 * tau[mask])
        ds = np.zeros(nt)
        dds = np.zeros(nt)
        ds[mask] = (w * np.cos(w * tau[mask]) - 80.0 * np.sin(w * tau[mask])) * np.exp(-80.0 * tau[mask])
        dds[mask] = ((-w * w * np.sin(w * tau[mask]) - 160.0 * w * np.cos(w * tau[mask])
                      + 6400.0 * np.sin(w * tau[mask])) * np.exp(-80.0 * tau[mask]))
        row = 3 * idx
        dis[row, :] = amp * s
        vel[row, :] = amp * ds
        acc[row, :] = amp * dds

    if p.exists():
        p.unlink()
    with h5py.File(p, "w") as f:
        d = f.create_group("DRM_Data")
        d["xyz"] = xyz
        d["internal"] = internal_arr
        d["data_location"] = (3 * np.arange(n, dtype=np.int32))
        d.create_dataset("displacement", data=dis)
        d.create_dataset("velocity", data=vel)
        d.create_dataset("acceleration", data=acc)
        m = f.create_group("DRM_Metadata")
        m["dt"] = float(dt)
        m["tstart"] = 0.0
        m["tend"] = float((nt - 1) * dt)
        m["drmbox_x0"] = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        for kname in ("drmbox_xmin", "drmbox_xmax", "drmbox_ymin", "drmbox_ymax", "drmbox_zmin", "drmbox_zmax"):
            m[kname] = float(cfg[kname.split("_")[1]])
        m["h"] = np.array([dx, dx, dx], dtype=np.float64)
    return d2_h5drm_path


def _extract_h5(h5_path: Path, run_id: str) -> Dict[str, Any]:
    import h5py

    def read_grp(g) -> Dict[str, Any]:
        node = {
            "wall_ns": int(g.attrs.get("wall_ns", 0)),
            "calls": int(g.attrs.get("calls", 0)),
            "elem_by_type": [],
            "children": {},
        }
        if "elem_by_type" in g:
            for row in g["elem_by_type"][:]:
                node["elem_by_type"].append([int(row["classTag"]), int(row["count"]), int(row["wall_ns"])])
        for k in g.keys():
            if k == "elem_by_type":
                continue
            obj = g[k]
            if isinstance(obj, h5py.Group):
                node["children"][k] = read_grp(obj)
        return node

    with h5py.File(h5_path, "r") as f:
        grp = f["runs"][run_id]
        tree = read_grp(grp["rollup"]["root"])
        series_ms = (
            [float(x) for x in grp["series"]["wall_ms"][:].reshape(-1)]
            if "series" in grp and "wall_ms" in grp["series"]
            else []
        )
    return {"tree": tree, "series_wall_ms": series_ms}


# --------------------------------------------------------------------------- #
#  Child execution (NEVER imports h5py!)                                      #
# --------------------------------------------------------------------------- #
def _run_child(args) -> int:
    aff_mask = None if args.affinity.lower() == "none" else int(args.affinity, 0)
    aff_info = _apply_and_read_affinity(aff_mask)
    pyd_path = Path(args.pyd).resolve()
    bin_dir = pyd_path.parent
    if sys.platform == "win32":
        os.add_dll_directory(str(bin_dir))
    sys.path.insert(0, str(bin_dir))
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    import opensees as ops  # noqa: E402
    import decks            # noqa: E402

    build_id = str(ops.ladrunoBuild()).strip()
    ops.ladrunoThreads(args.threads)

    if args.verify_d2_linear_child:
        # Compare Linear -factorOnce vs Newton nodal displacements on D2 smoke
        res_vecs = {}
        for algo in ("LinearFactorOnce", "Newton"):
            cfg = decks.run_d2_cpml_wave(ops, size="smoke", h5drm_path=args.h5drm, algo=algo)
            for _ in range(cfg["nsteps"]):
                if ops.analyze(1, cfg["dt"]) != 0:
                    raise RuntimeError(f"D2 verify ({algo}) analyze failed")
            tags = sorted(ops.getNodeTags())
            res_vecs[algo] = [v for t in tags for v in ops.nodeDisp(t)]
        a, b = res_vecs["LinearFactorOnce"], res_vecs["Newton"]
        max_diff = max(abs(x - y) for x, y in zip(a, b))
        max_u = max(abs(x) for x in b)
        print(f"D2 Linear -factorOnce vs Newton: max|u|={max_u:.6e}, max|du|={max_diff:.6e}")
        if not (max_u > 1e-9 and max_diff < 1e-10):
            raise RuntimeError(f"D2 Linear -factorOnce mismatch: max|u|={max_u}, max|du|={max_diff}")
        return 0

    cfg = decks.DECKS[args.deck](ops, size=args.size, h5drm_path=args.h5drm)
    nsteps = cfg["nsteps"]

    ops.profiler("reset")
    flags = ["-deep", "-perStep"] if args.deep else ["-perStep"]
    ops.profiler("start", *flags)
    t0 = time.perf_counter()
    try:
        if cfg["kind"] == "MPM_ModeA":
            for _ in range(nsteps):
                rc = ops.mpm("analyze", 1, cfg["dt"])
                if rc != 0:
                    raise RuntimeError(f"{args.deck}: mpm analyze returned {rc}")
        elif cfg["kind"] == "Transient":
            rc = ops.analyze(nsteps, cfg["dt"])
            if rc != 0:
                raise RuntimeError(f"{args.deck}: analyze returned {rc}")
        else:
            rc = ops.analyze(nsteps)
            if rc != 0:
                raise RuntimeError(f"{args.deck}: analyze returned {rc}")
    finally:
        wall_s = time.perf_counter() - t0
        ops.profiler("stop")

    h5_path = Path(args.h5)
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    run_id = f"{args.deck}_{args.size}_{args.mode}_r{args.rep}"
    ops.profiler("report", str(h5_path), "-run", run_id)
    ops.wipe()

    payload = {
        "ok": True,
        "deck": args.deck,
        "size": args.size,
        "family": cfg["family"],
        "mode": args.mode,
        "rep": args.rep,
        "threads": args.threads,
        "deep": args.deep,
        "nsteps": nsteps,
        "nelm": cfg["nelm"],
        "ndof_est": cfg.get("ndof_est", 0),
        "wall_s": wall_s,
        "h5": str(h5_path),
        "run_id": run_id,
        "build_id": build_id,
        "pyd": str(ops.__file__),
        "affinity": aff_info,
        "env": _collect_env(),
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


# Note: Child process above imports `opensees`, never `h5py`. To keep `h5py` strictly out of child,
# guard top-level h5py/numpy usage inside parent helpers or check `--child` before importing h5py.
# Wait — in this file, let's defer `h5py`/`numpy` imports inside parent-only functions!
def _parent_write_h5drm(out_dir: Path, size: str) -> Path:
    return _ensure_h5drm(out_dir, size)


def _spawn_one(pyd: Path, out_dir: Path, deck: str, size: str, mode: str, threads: int,
               deep: bool, rep: int, affinity: str, h5drm_path: Optional[Path]) -> Dict[str, Any]:
    tag = f"{deck}__{size}__{mode}__r{rep}"
    h5_path = out_dir / f"{tag}.h5"
    json_path = out_dir / f"{tag}.json"
    for p in (h5_path, json_path):
        if p.exists():
            p.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    # Never override OMP_WAIT_POLICY or KMP_BLOCKTIME (fork defaults preserved per brief §2)

    cmd = [
        sys.executable, str(Path(__file__).resolve()),
        "--child",
        "--pyd", str(pyd),
        "--deck", deck,
        "--size", size,
        "--mode", mode,
        "--threads", str(threads),
        "--rep", str(rep),
        "--affinity", affinity,
        "--h5", str(h5_path),
        "--json-out", str(json_path),
    ]
    if deep:
        cmd.append("--deep")
    if h5drm_path is not None:
        cmd.extend(["--h5drm", str(h5drm_path)])

    cp = subprocess.run(cmd, env=env, text=True, capture_output=True)
    stderr_tail = (cp.stderr or "")[-4000:]
    if cp.returncode != 0 or not json_path.exists():
        return {
            "ok": False, "deck": deck, "size": size, "mode": mode, "rep": rep,
            "threads": threads, "returncode": cp.returncode,
            "stdout_tail": (cp.stdout or "")[-2000:],
            "stderr_tail": stderr_tail,
        }
    rec = json.loads(json_path.read_text(encoding="utf-8"))
    rec["stderr_tail"] = stderr_tail
    threaded, anns = _parse_threading(threads, stderr_tail)
    rec["threaded"] = threaded
    rec["threads_announcements"] = anns
    parsed = _extract_h5(h5_path, rec["run_id"])
    series_ms = parsed["series_wall_ms"]
    rec["series_wall_ms"] = series_ms
    rec["step1_wall_ms"] = series_ms[0] if series_ms else None
    rec["steady_wall_s"] = (sum(series_ms[1:]) / 1000.0) if len(series_ms) > 1 else rec["wall_s"]
    rec["_tree"] = parsed["tree"]
    return rec


def main() -> int:
    # Child mode exits BEFORE parent h5py/numpy code executes
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyd", type=Path, default=DEFAULT_PYD)
    ap.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--export-dir", type=Path, default=None)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--only", type=str, default="")
    ap.add_argument("--sizes", type=str, default="")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--affinity", type=str, default="none")
    ap.add_argument("--verify-d2-linear", action="store_true")
    # child-only flags
    ap.add_argument("--child", action="store_true")
    ap.add_argument("--verify-d2-linear-child", action="store_true")
    ap.add_argument("--deck", type=str)
    ap.add_argument("--size", type=str, default="S")
    ap.add_argument("--mode", type=str)
    ap.add_argument("--rep", type=int, default=0)
    ap.add_argument("--h5", type=str)
    ap.add_argument("--json-out", type=str)
    ap.add_argument("--h5drm", type=str, default=None)
    ap.add_argument("--deep", action="store_true")
    args = ap.parse_args()

    if args.child or args.verify_d2_linear_child:
        return _run_child(args)

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import decks  # noqa: E402

    pyd = args.pyd.resolve()
    exe = args.exe.resolve()
    out_dir = args.out.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not pyd.exists() or not exe.exists():
        print(f"ERROR: missing binary pyd={pyd.exists()} exe={exe.exists()}", file=sys.stderr)
        return 2

    exe_build = _probe_exe_build(exe)
    print(f"OpenSees.exe ladrunoBuild: {exe_build}")

    if args.verify_d2_linear:
        drm_smoke = _parent_write_h5drm(out_dir, "smoke")
        cp = subprocess.run([
            sys.executable, str(Path(__file__).resolve()),
            "--verify-d2-linear-child", "--pyd", str(pyd), "--h5drm", str(drm_smoke),
            "--affinity", args.affinity, "--threads", "1",
        ], text=True, capture_output=True)
        print((cp.stdout or "").strip())
        if cp.returncode != 0:
            print(cp.stderr, file=sys.stderr)
            return cp.returncode

    only_decks = set(args.only.split(",")) if args.only else set(decks.DECKS)
    filter_sizes = set(args.sizes.split(",")) if args.sizes else None

    targets: List[Tuple[str, str]] = []
    for d in decks.DECKS:
        if d not in only_decks:
            continue
        if args.smoke:
            szs = ["smoke"]
        else:
            szs = [s for s in decks.DECK_SIZES[d] if (not filter_sizes or s in filter_sizes)]
        for s in szs:
            targets.append((d, s))

    h5drm_files: Dict[str, Path] = {}
    for d, s in targets:
        if d == "d2_cpml_wave" and s not in h5drm_files:
            h5drm_files[s] = _parent_write_h5drm(out_dir, s)

    base_reps = 1 if args.smoke or not args.full else 3
    modes = [("1t_deep", 1, True), ("1t_coarse", 1, False), ("nt_coarse", args.threads, False)]

    records: List[Dict[str, Any]] = []
    trees: Dict[str, Any] = {}

    def run_combo(rep_idx: int, deck: str, size: str, mode: str, nt: int, deep: bool) -> Dict[str, Any]:
        print(f"[rep {rep_idx + 1}] {deck:25s} [{size:5s}] {mode:10s} ...", end="", flush=True)
        rec = _spawn_one(pyd, out_dir, deck, size, mode, nt, deep, rep_idx, args.affinity,
                         h5drm_files.get(size) if deck == "d2_cpml_wave" else None)
        tr = rec.pop("_tree", None)
        if tr is not None:
            trees[rec["run_id"]] = {"tree": tr, "series_wall_ms": rec.get("series_wall_ms", [])}
        records.append(rec)
        if rec.get("ok"):
            print(f" steady={rec['steady_wall_s']:7.3f}s (tot={rec['wall_s']:6.3f}s) {rec['threaded']}", flush=True)
        else:
            print(f" FAILED (rc={rec.get('returncode')})", flush=True)
            print(rec.get("stderr_tail", ""), file=sys.stderr)
        return rec

    for rep in range(base_reps):
        for deck, size in targets:
            for mode, nt, deep in modes:
                run_combo(rep, deck, size, mode, nt, deep)

    # Brief §2: raise from 3 to 5 reps if spread > 5% on steady_wall_s
    if args.full and base_reps == 3:
        for deck, size in targets:
            for mode, nt, deep in modes:
                vals = [r["steady_wall_s"] for r in records
                        if r.get("ok") and r["deck"] == deck and r["size"] == size and r["mode"] == mode]
                if len(vals) == 3:
                    med = statistics.median(vals)
                    spread = (max(vals) - min(vals)) / med if med > 0 else 0.0
                    if spread > 0.05:
                        print(f"  -> spread {100 * spread:.1f}% > 5% on {deck}[{size}] {mode}; running reps 4..5")
                        for rep in (3, 4):
                            run_combo(rep, deck, size, mode, nt, deep)

    mpath = out_dir / "manifest.json"
    if args.only and mpath.exists():
        old = json.loads(mpath.read_text(encoding="utf-8"))
        kept = [r for r in old.get("records", []) if (r.get("deck"), r.get("size")) not in set(targets)]
        records = kept + records

    manifest = {
        "exe_build": exe_build,
        "pyd": str(pyd),
        "smoke": args.smoke,
        "production_threads": args.threads,
        "affinity_arg": args.affinity,
        "parent_env": _collect_env(),
        "base_reps": base_reps,
        "records": records,
    }
    mpath.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    tpath = out_dir / "trees.json"
    tpath.write_text(json.dumps(trees, indent=1), encoding="utf-8")
    print(f"Wrote {mpath} and {tpath}")

    if args.export_dir is not None:
        exp = args.export_dir.resolve()
        exp.mkdir(parents=True, exist_ok=True)
        (exp / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        (exp / "trees.json").write_text(json.dumps(trees, indent=1), encoding="utf-8")
        print(f"Exported JSON artifacts to {exp}")

    failed = [r for r in records if not r.get("ok")]
    return 1 if failed else 0


if __name__ == "__main__":
    # Enforce child isolation: only import h5py/numpy in parent invocations
    if "--child" in sys.argv or "--verify-d2-linear-child" in sys.argv:
        pass
    raise SystemExit(main())
