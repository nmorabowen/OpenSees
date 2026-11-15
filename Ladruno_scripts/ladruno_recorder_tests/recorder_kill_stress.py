"""Kill stress of the LadrunoRecorder write pattern, on its h5py mirror (no OpenSees binary).

``recorder_crash_replay.py`` checks every byte-level kill point but cannot say how often a real
kill lands in a dangerous one: that depends on how long each phase takes. This script measures
it: it runs the replay model's writer (the same HDF5 call sequence, through the default sec2
driver onto a real file, with the recorder's WALL-CLOCK ``-flush`` cadence) in a child process,
kills it with TerminateProcess (``Popen.kill``) after a random wait, and reads the file with the
replay model's checker. It is the mirror, not the C++ recorder: the Python overhead per step
dilutes the time spent inside ``H5Fflush``, so the C++ recorder's rate with ``-flush 0`` is
expected to be HIGHER (WP-164 measured the per-step flush at ~99 % of a small explicit run).

    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_kill_stress.py --kills 200 --flush 0
    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_kill_stress.py --kills 200 --flush 10
"""
from __future__ import annotations

import argparse
import collections
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recorder_crash_replay as rcr


def child(a):
    results = rcr.parse_results(a.results)
    state = {"last": None}

    def due(k):
        now = time.perf_counter()
        if a.work_ms > 0:   # stand-in for the analysis step between two records
            end = now + a.work_ms / 1000.0
            while time.perf_counter() < end:
                pass
            now = time.perf_counter()
        if state["last"] is None or a.flush <= 0 or now - state["last"] >= a.flush:
            state["last"] = now
            return True
        return False

    a.steps = 10 ** 9
    a.check_from, a.check_to, a.verbose = 1, 0, False
    rcr.write_run(rcr.NoLog(), a, results, real_path=a.child, flush_due=due)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kills", type=int, default=100)
    ap.add_argument("--flush", type=float, default=0.0, help="-flush seconds (0 = every step)")
    ap.add_argument("--min-wait", type=float, default=1.5, help="earliest kill, s after launch")
    ap.add_argument("--max-wait", type=float, default=6.0)
    ap.add_argument("--nodes", type=int, default=200)
    ap.add_argument("--results", default="DISPLACEMENT:3,VELOCITY:3")
    ap.add_argument("--compress", type=int, default=1)
    ap.add_argument("--precision", choices=("f64", "f32"), default="f64")
    ap.add_argument("--chunk-scale", type=float, default=1.0)
    ap.add_argument("--scale-axis", action="store_true")
    ap.add_argument("--work-ms", type=float, default=0.0, help="busy analysis time per step")
    ap.add_argument("--swmr", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--keep", default="", help="directory to copy failing files into")
    ap.add_argument("--child", default="", help=argparse.SUPPRESS)
    ap.add_argument("--flush-every", type=int, default=0, help=argparse.SUPPRESS)
    a = ap.parse_args()
    if a.child:
        child(a)
        return

    results = rcr.parse_results(a.results)
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp(prefix="rec_kill_")
    stats = collections.Counter()
    details = []
    passthru = ["--flush", str(a.flush), "--nodes", str(a.nodes), "--results", a.results,
                "--compress", str(a.compress), "--precision", a.precision,
                "--chunk-scale", str(a.chunk_scale), "--work-ms", str(a.work_ms)]
    if a.scale_axis:
        passthru.append("--scale-axis")
    if a.swmr:
        passthru.append("--swmr")
    t0 = time.time()
    for i in range(a.kills):
        path = os.path.join(tmp, f"k{i}.ladruno")
        p = subprocess.Popen([sys.executable, __file__, "--child", path] + passthru)
        wait = rng.uniform(a.min_wait, a.max_wait)
        time.sleep(wait)
        if p.poll() is not None:
            stats["CHILD_EXITED"] += 1
            continue
        p.kill()                      # TerminateProcess on Windows
        p.wait()
        status, n, detail = rcr.check(path, results, None, a.swmr, a.precision)
        stats[status] += 1
        if status != "OK":
            details.append((i, round(wait, 2), status, n, detail))
            if a.keep:
                os.makedirs(a.keep, exist_ok=True)
                os.replace(path, os.path.join(a.keep, f"k{i}_{status}.ladruno"))
                continue
        try:
            os.remove(path)
        except OSError:
            pass
    print(f"kills={a.kills} flush={a.flush}s nodes={a.nodes} results={a.results} "
          f"compress={a.compress} work_ms={a.work_ms} swmr={a.swmr} ({time.time() - t0:.0f} s)")
    print(dict(stats))
    for d in details[:30]:
        print("  ", d)


if __name__ == "__main__":
    main()
