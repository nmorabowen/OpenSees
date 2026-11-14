"""Kill stress across ``mpm reset``: the stage-layout window of the MPM particle writer.

What it shows: whether a run killed while the writer lays out a new stage leaves a file
whose earlier stages can still be read. The child deck alternates ``mpm analyze NSTEPS``
and ``mpm reset`` against the repo's ``dist/bin/opensees.pyd``, so the writer closes the
file, reopens it read-write (no SWMR), lays out the next stage and re-arms SWMR many times
a second. Each killed file is classified by the 232a section 2.5 reader rule:

- ``plain`` / ``swmr``: opened by ``tests/mpm_contracts.py::open_particle_file``; every
  stage's committed outputs are checked (ids exact, volumes > 0, data finite);
- ``refused``: both opens refuse ("not already open for SWMR writing", reader rule step 3,
  ``StageLayoutKilledError``): the run was killed between a reset's read-write reopen and
  ``H5Fstart_swmr_write``, which leaves the superblock flags at WRITE_ACCESS (0x01) only.
  The file is then recovered on a COPY with ``recover_particle_file`` and every recovered
  stage is checked the same way; only the last stage may be "torn while laid out".

A trial is bad when a stage that opens holds bad rows, or a refused file does not recover.
Failing trials keep their files under ``<out_dir>/t<nnnn>/``; all results go to
``<out_dir>/results.json``.

Measured (2026-11-14, build 9 cbe1b7149, a reset every 100 steps, 400 kills): 290 opened
with SWMR-read, 89 plainly, 21 (5.25 %) refused by both opens; 21 / 21 recovered every
earlier stage (in 2, the stage being laid out did not open: no committed output).

Run from the repo root, after a full ``build.bat`` (the pyd must be current):

    python3.12 Ladruno_scripts/mpm_crash_tests/reset_stress.py 400 4 <out_dir> [nsteps]

Arguments: number of trials, parallel workers, output directory, steps per analyze
between resets (default 100).
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
out, fixture, nsteps = sys.argv[1], sys.argv[2], int(sys.argv[3])
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
ops.nDMaterial("ElasticIsotropic", 2, 2.0e4, 0.3)
ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "gimp")
ops.mpm("particles", "-file", fixture)
ops.mpm("output", out, "-every", 1, "-fields", "POSITION", "VELOCITY", "STRESS", "VOLUME")
print("CHILD START", flush=True)
while True:
    ops.mpm("analyze", nsteps, 1e-6)
    ops.mpm("reset")
"""
FIELDS = ("POSITION", "VELOCITY", "STRESS", "VOLUME")


def check_stage(st, set_ids):
    """One stage: (status, committed outputs)."""
    n = mc.committed_outputs(st)
    op = st["RESULTS/ON_PARTICLES"]
    offs = op["OFFSETS"][: n + 1]
    end = int(offs[n]) if n else 0
    ids = op["ID"][:end]
    vol = op["VOLUME"]["DATA"][:end, 0]
    if not all(np.all(np.isfinite(op[k]["DATA"][:end])) for k in FIELDS):
        return "NONFINITE", n
    bad = [k for k in range(n)
           if not np.array_equal(ids[offs[k]:offs[k + 1]], set_ids)
           or not np.all(vol[offs[k]:offs[k + 1]] > 0.0)]
    return ("BAD_ROWS" if bad else "OK"), n


def check_file(f, set_ids):
    """Every particle stage of an open file: [(name, status, committed)]."""
    out = []
    for name in sorted((k for k in f if k.startswith("MODEL_STAGE[")), key=lambda s: int(s[12:-1])):
        st = f[name]
        if "RESULTS/ON_PARTICLES/OFFSETS" not in st:
            out.append((name, "NO_RESULTS", 0))
            continue
        try:
            out.append((name, *check_stage(st, set_ids)))
        except Exception as e:  # noqa: BLE001
            out.append((name, "STAGE_FAIL: " + str(e).splitlines()[0][:100], -1))
    return out


def trial(i, wait, out_dir, nsteps, set_ids):
    d = os.path.join(out_dir, f"t{i:04d}")
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "killed.ladruno")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1", LADRUNO_MPM_TESTS=str(TESTS),
               PYTHONPATH=str(PYD.parent), LADRUNO_EXPECT_OPENSEES=str(PYD))
    env.pop("HDF5_USE_FILE_LOCKING", None)
    p = subprocess.Popen([sys.executable, "-c", CHILD, out, str(FIXTURE), str(nsteps)],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         stdin=subprocess.DEVNULL, text=True, env=env)
    try:
        line = p.stdout.readline()
        if "CHILD START" not in line:
            p.kill()
            return {"i": i, "wait": wait, "mode": "CHILD_FAILED", "bad": True,
                    "detail": line + p.stderr.read()[-400:]}
        time.sleep(wait)
    finally:
        p.kill()
        p.wait()
    res = {"i": i, "wait": round(wait, 3)}
    try:
        with mc.open_particle_file(out) as f:
            res["mode"] = "swmr" if f.swmr_mode else "plain"
            res["stages"] = check_file(f, set_ids)
        bad = any(s[1] not in ("OK", "NO_RESULTS") for s in res["stages"])
    except mc.StageLayoutKilledError:
        res["mode"] = "refused"
        try:
            rep = mc.recover_particle_file(out, os.path.join(d, "recovered.ladruno"))
            res["recovered"] = rep
            with h5py.File(os.path.join(d, "recovered.ladruno"), "r") as f:
                rows = {name: status for name, status, _ in check_file(f, set_ids)}
            done = [r["stage"] for r in rep if "committed" in r]
            bad = any(rows.get(s) != "OK" for s in done)
        except Exception as e:  # noqa: BLE001
            res["recovered"] = "RECOVERY_FAILED: " + str(e).splitlines()[0][:200]
            bad = True
    except Exception as e:  # noqa: BLE001
        res["mode"] = "UNREADABLE"
        res["detail"] = str(e).splitlines()[0][:200]
        bad = True
    if not bad:
        shutil.rmtree(d, ignore_errors=True)
    res["bad"] = bad
    return res


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    n_trials, workers, out_dir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    nsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    if not PYD.exists():
        sys.exit(f"no binary at {PYD}: run Ladruno_scripts\\build.bat first")
    os.makedirs(out_dir, exist_ok=True)
    with h5py.File(FIXTURE, "r") as f:
        set_ids = f["PARTICLES/ID"][...]
    rng = random.Random(1232)
    waits = [rng.uniform(0.2, 3.0) for _ in range(n_trials)]
    results = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(trial, i, w, out_dir, nsteps, set_ids) for i, w in enumerate(waits)]
        for k, fu in enumerate(cf.as_completed(futs), 1):
            r = fu.result()
            results.append(r)
            if r.get("bad"):
                print("BAD", json.dumps(r)[:700], flush=True)
            if k % 50 == 0:
                print(f"progress {k}/{n_trials} bad={sum(1 for x in results if x.get('bad'))} "
                      f"elapsed={time.time() - t0:.0f}s", flush=True)
    results.sort(key=lambda r: r["i"])
    with open(os.path.join(out_dir, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    refused = [r for r in results if r["mode"] == "refused"]
    torn = sum(1 for r in refused if isinstance(r.get("recovered"), list)
               and any(x.get("torn_while_laid_out") for x in r["recovered"]))
    print("open modes:", dict(Counter(r["mode"] for r in results)))
    print("stage statuses (opened files):",
          dict(Counter(s[1] for r in results for s in r.get("stages", []))))
    print(f"refused files: {len(refused)}, recovered "
          f"{sum(1 for r in refused if isinstance(r.get('recovered'), list))}, "
          f"with the last stage torn while laid out: {torn}")
    bad = sum(1 for r in results if r.get("bad"))
    print("SUMMARY bad", bad, "of", len(results), f"elapsed {time.time() - t0:.0f}s", flush=True)
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
