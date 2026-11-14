"""Kill stress for the MPM particle writer (WP-232c, 232a section 2.5).

What it shows: whether a run of the REAL binary killed mid-``mpm analyze`` always leaves
a ``.ladruno`` particle file that a crash reader can read. Each trial runs the strict kill
test's child deck (``tests/test_ladruno_mpm_contracts_io.py``, ``_KILL_CHILD``) against the
repo's ``dist/bin/opensees.pyd``, kills it after a wait (half the trials at the kill
test's fixed waits, half at random waits), then checks the file:

- it opens by the reader rule (``tests/mpm_contracts.py::open_particle_file``: a plain
  read-only open, else SWMR-read);
- ``committed_outputs()`` (the crash rule) does not raise;
- the committed rows read, and match the strict kill test's assertions: ids exact,
  volumes > 0, data finite, OFFSETS non-decreasing, STEP increasing, at least 1 output.

It also counts how many killed files a plain open accepts (a killed SWMR file needs the
SWMR-read open). Failing trials keep their file under ``<out_dir>/t<nnnn>/``; all results
go to ``<out_dir>/results.json``.

Measured (2026-11-14): build 7 (barrier, no SWMR) 1 unreadable in 137 kills; build 8
(eviction off) 3 in 400; build 9 (SWMR, cbe1b7149) 0 in 400, every killed file refused by
a plain open and read through ``open_particle_file``.

Run from the repo root, after a full ``build.bat`` (the pyd must be current):

    python3.12 Ladruno_scripts/mpm_crash_tests/kill_stress.py 400 4 <out_dir>

Arguments: number of trials, parallel workers, output directory.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import os
import random
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
PYD = ROOT / "dist" / "bin" / "opensees.pyd"
FIXTURE = TESTS / "data" / "mpm" / "particle_set_v1.h5"
sys.path.insert(0, str(TESTS))
sys.path.insert(0, str(ROOT / "Ladruno_scripts" / "ladruno_recorder_tests"))
import mpm_contracts as mc

CHILD = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_MPM_TESTS"])
from _testbed import ops
out, fixture = sys.argv[1], sys.argv[2]
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
ops.nDMaterial("ElasticIsotropic", 2, 2.0e4, 0.3)
ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "gimp")
ops.mpm("particles", "-file", fixture)
ops.mpm("output", out, "-every", 1, "-fields", "POSITION", "VELOCITY", "STRESS", "VOLUME")
print("CHILD START", flush=True)
ops.mpm("analyze", 50000000, 1e-6)
"""
FIXED = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
FIELDS = ("POSITION", "VELOCITY", "STRESS", "VOLUME")


def check(path, set_ids):
    """Returns (status, committed, T, detail) for one killed file."""
    with mc.open_particle_file(path) as f:
        stage = f["MODEL_STAGE[0]"]
        n = mc.committed_outputs(stage)
        op = stage["RESULTS/ON_PARTICLES"]
        steps, offs = op["STEP"][...], op["OFFSETS"][...]
        nid = op["ID"].shape[0]
        T = min(len(steps), len(op["TIME"]), len(offs) - 1)
        while T > 0 and (offs[T] > nid or any(offs[T] > op[k]["DATA"].shape[0] for k in FIELDS)):
            T -= 1
        if T < 1:
            return "NO_OUTPUT", n, T, f"T={T}"
        ids = op["ID"][: offs[T]]
        for k in FIELDS:
            if not np.all(np.isfinite(op[k]["DATA"][: offs[T]])):
                return "NONFINITE", n, T, k
        if not np.all(np.diff(offs[: T + 1]) >= 0):
            return "OFFS_DECREASE", n, T, ""
        vol = op["VOLUME"]["DATA"][: offs[T], 0]
        for k in range(T):
            if not np.array_equal(ids[offs[k]:offs[k + 1]], set_ids):
                return "ID_MISMATCH", n, T, f"output {k}"
            if not np.all(vol[offs[k]:offs[k + 1]] > 0.0):
                return "VOL_NONPOS", n, T, f"output {k}"
        if not np.all(np.diff(steps[:T]) > 0):
            return "STEP_ORDER", n, T, ""
        return "OK", n, T, ""


def trial(i, wait, out_dir, set_ids):
    d = os.path.join(out_dir, f"t{i:04d}")
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "killed.ladruno")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1", LADRUNO_MPM_TESTS=str(TESTS),
               PYTHONPATH=str(PYD.parent), LADRUNO_EXPECT_OPENSEES=str(PYD))
    env.pop("HDF5_USE_FILE_LOCKING", None)
    p = subprocess.Popen([sys.executable, "-c", CHILD, out, str(FIXTURE)], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, stdin=subprocess.DEVNULL, text=True, env=env)
    try:
        line = p.stdout.readline()
        if "CHILD START" not in line:
            p.kill()
            return {"i": i, "wait": wait, "status": "CHILD_FAILED", "detail": line + p.stderr.read()[-400:]}
        time.sleep(wait)
    finally:
        p.kill()
        p.wait()
    try:
        with h5py.File(out, "r"):
            plain = True
    except Exception:  # noqa: BLE001
        plain = False
    try:
        status, n, T, detail = check(out, set_ids)
    except Exception as e:  # noqa: BLE001
        status, n, T, detail = "UNREADABLE", -1, -1, "".join(
            traceback.format_exception_only(type(e), e)).strip()[-500:]
    size = os.path.getsize(out) if os.path.exists(out) else -1
    if status == "OK":
        shutil.rmtree(d, ignore_errors=True)
    return {"i": i, "wait": round(wait, 3), "status": status, "committed": n, "T": T, "size": size,
                "plain_open": plain, "detail": detail}


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    n_trials, workers, out_dir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    if not PYD.exists():
        sys.exit(f"no binary at {PYD}: run Ladruno_scripts\\build.bat first")
    os.makedirs(out_dir, exist_ok=True)
    with h5py.File(FIXTURE, "r") as f:
        set_ids = f["PARTICLES/ID"][...]
    rng = random.Random(232)
    waits = [FIXED[i % len(FIXED)] if i % 2 == 0 else rng.uniform(0.2, 4.0) for i in range(n_trials)]
    results = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(trial, i, w, out_dir, set_ids) for i, w in enumerate(waits)]
        for k, fu in enumerate(cf.as_completed(futs), 1):
            r = fu.result()
            results.append(r)
            if r["status"] != "OK":
                print("FAIL", json.dumps(r), flush=True)
            if k % 50 == 0:
                bad = sum(1 for x in results if x["status"] != "OK")
                print(f"progress {k}/{n_trials} failures={bad} elapsed={time.time() - t0:.0f}s", flush=True)
    results.sort(key=lambda r: r["i"])
    with open(os.path.join(out_dir, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    c = Counter(r["status"] for r in results)
    ok = [r for r in results if r["status"] == "OK"]
    print("plain read-only open worked on", sum(1 for r in results if r.get("plain_open")),
          "of", len(results), "killed files (the rest need the SWMR-read rule)", flush=True)
    print("SUMMARY", dict(c), "trials", len(results),
          "committed range", (min(r["committed"] for r in ok), max(r["committed"] for r in ok)) if ok else None,
          f"elapsed {time.time() - t0:.0f}s", flush=True)
    sys.exit(0 if c.get("OK", 0) == len(results) else 1)


if __name__ == "__main__":
    main()
