"""WP-260 Wave 0 (`260-w0`) profiler runner — parent/child HDF5-isolated harness.

Usage:
    python Ladruno_scripts/wp260_profiles/run_profiles.py --smoke
    python Ladruno_scripts/wp260_profiles/run_profiles.py --full
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PYD = ROOT / "dist" / "bin" / "opensees.pyd"
DEFAULT_EXE = ROOT / "dist" / "bin" / "OpenSees.exe"
DEFAULT_OUT = Path(tempfile.gettempdir()) / "wp260w0_prof"

# i7-13700KF P-cores (8 physical cores x 2 HT = logical 0..15)
P_CORES_MASK = 0xFFFF


def _pin_affinity(mask: int) -> None:
    if sys.platform != "win32" or not mask:
        return
    try:
        k32 = ctypes.windll.kernel32
        k32.SetProcessAffinityMask(k32.GetCurrentProcess(), ctypes.c_size_t(mask))
    except Exception:
        pass


def _probe_exe_build(exe_path: Path) -> str:
    tcl = "puts [ladrunoBuild]; exit\n"
    cp = subprocess.run([str(exe_path)], input=tcl, text=True, capture_output=True, check=True)
    combined = cp.stdout + "\n" + cp.stderr
    for ln in reversed(combined.splitlines()):
        s = ln.strip()
        if s and len(s) >= 7 and all(c in "0123456789abcdefABCDEF-dirty" for c in s):
            return s
    return ""


def _run_child(args) -> int:
    """Child process ONLY — NEVER import h5py here (HDF5 DLL clash with opensees.pyd)."""
    _pin_affinity(P_CORES_MASK)
    pyd_path = Path(args.pyd).resolve()
    bin_dir = pyd_path.parent
    if sys.platform == "win32":
        os.add_dll_directory(str(bin_dir))
        # Intel oneAPI MKL / MPI runtime if present
        for env_k in ("MKL_ROOT", "ONEAPI_ROOT", "I_MPI_ROOT"):
            v = os.environ.get(env_k)
            if v:
                for sub in ("bin", "redist/intel64", "lib/intel64"):
                    d = Path(v) / sub
                    if d.is_dir():
                        try:
                            os.add_dll_directory(str(d))
                        except OSError:
                            pass
    sys.path.insert(0, str(bin_dir))
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    import opensees as ops   # noqa: E402
    import decks             # noqa: E402

    build_id = str(ops.ladrunoBuild()).strip()
    ops.ladrunoThreads(args.threads)

    cfg = decks.DECKS[args.deck](ops, smoke=args.smoke)
    nsteps = cfg["nsteps"]

    ops.profiler("reset")
    flags = ["-deep", "-perStep"] if args.deep else ["-perStep"]
    ops.profiler("start", *flags)
    t0 = time.perf_counter()
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
    wall_s = time.perf_counter() - t0
    ops.profiler("stop")

    h5_path = Path(args.h5)
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    run_id = f"{args.deck}_{args.mode}_r{args.rep}"
    ops.profiler("report", str(h5_path), "-run", run_id)

    payload = {
        "ok": True,
        "deck": args.deck,
        "family": cfg["family"],
        "mode": args.mode,
        "rep": args.rep,
        "threads": args.threads,
        "deep": args.deep,
        "nsteps": nsteps,
        "nelm": cfg["nelm"],
        "wall_s": wall_s,
        "h5": str(h5_path),
        "run_id": run_id,
        "build_id": build_id,
        "pyd": str(ops.__file__),
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


def _spawn_one(pyd: Path, out_dir: Path, deck: str, mode: str, threads: int,
               deep: bool, rep: int, smoke: bool) -> Dict:
    tag = f"{deck}__{mode}__r{rep}"
    h5_path = out_dir / f"{tag}.h5"
    json_path = out_dir / f"{tag}.json"
    if h5_path.exists():
        h5_path.unlink()
    if json_path.exists():
        json_path.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    env["OMP_WAIT_POLICY"] = "PASSIVE"

    cmd = [
        sys.executable, str(Path(__file__).resolve()),
        "--child",
        "--pyd", str(pyd),
        "--deck", deck,
        "--mode", mode,
        "--threads", str(threads),
        "--rep", str(rep),
        "--h5", str(h5_path),
        "--json-out", str(json_path),
    ]
    if deep:
        cmd.append("--deep")
    if smoke:
        cmd.append("--smoke")

    cp = subprocess.run(cmd, env=env, text=True, capture_output=True)
    if cp.returncode != 0 or not json_path.exists():
        return {
            "ok": False,
            "deck": deck,
            "mode": mode,
            "rep": rep,
            "returncode": cp.returncode,
            "stdout": cp.stdout[-2000:],
            "stderr": cp.stderr[-2000:],
        }
    return json.loads(json_path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyd", type=Path, default=DEFAULT_PYD)
    ap.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--only", type=str, default="")
    # child-only flags
    ap.add_argument("--child", action="store_true")
    ap.add_argument("--deck", type=str)
    ap.add_argument("--mode", type=str)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--rep", type=int, default=0)
    ap.add_argument("--h5", type=str)
    ap.add_argument("--json-out", type=str)
    ap.add_argument("--deep", action="store_true")
    args = ap.parse_args()

    if args.child:
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

    selected: List[str] = [d for d in decks.DECKS if (not args.only or d in args.only.split(","))]
    reps = 1 if args.smoke or not args.full else 3
    modes = [("1t_deep", 1, True), ("8t_coarse", 8, False)]

    records: List[Dict] = []
    # Interleaved repetitions across decks
    for rep in range(reps):
        for deck in selected:
            for mode, nt, deep in modes:
                print(f"[rep {rep + 1}/{reps}] {deck:26s} {mode:10s} ...", end="", flush=True)
                rec = _spawn_one(pyd, out_dir, deck, mode, nt, deep, rep, args.smoke)
                records.append(rec)
                if rec.get("ok"):
                    print(f" {rec['wall_s']:7.3f} s  (build={rec['build_id'][:9]})", flush=True)
                else:
                    print(f" FAILED (rc={rec.get('returncode')})", flush=True)
                    print(rec.get("stderr", ""), file=sys.stderr)

    mpath = out_dir / "manifest.json"
    if args.only and mpath.exists():
        old = json.loads(mpath.read_text(encoding="utf-8"))
        kept = [r for r in old.get("records", []) if r.get("deck") not in set(selected)]
        records = kept + records
    manifest = {
        "exe_build": exe_build,
        "pyd": str(pyd),
        "smoke": args.smoke,
        "reps": reps,
        "records": records,
    }
    mpath.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {mpath}")
    failed = [r for r in records if not r.get("ok")]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
