"""WP-259 race stress: repeat the concrete determinism decks at 8 threads and compare digests.

A data race shows as a digest that differs on SOME runs. One pass of the pytest gate
(tests/test_wp259_concrete_threads.py) runs each thread count once; this script runs selected
configurations R times at the top thread count and checks every run against the 1-thread digest
(nodal displacements and load factor after every step, every step's return code, every Gauss
point's diagnostics, and the number of return-map warning lines). Same method as
Ladruno_scripts/wp258_race_stress.py (no TSan build: GCC's libgomp is not TSan-instrumented).

    set PYTHONPATH=<worktree>\\dist\\bin
    python3.12 Ladruno_scripts/wp259_race_stress.py --repeats 10 --threads 8 [--out stress.json]
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
sys.path.insert(0, str(TESTS))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import test_wp259_concrete_threads as FE                  # noqa: E402

CASES = [("c3d", "3d"), ("c3d", "3dcr"), ("c3d_implex", "shell"), ("c3d", "ps"),
         ("rc", "3d"), ("rc_implex", "shell"), ("rc_interlock", "shell"), ("rc", "ps")]


def _digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def digest(mat, host, n):
    r, text = FE._run_fe(mat, host, n)
    return _digest([r["steps"], r["pts"], text.count(FE._WARN), text.count(FE._SUPPRESSED)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows, bad = [], 0
    for mat, host in CASES:
        t0 = time.time()
        ref = digest(mat, host, 1)
        mism = sum(digest(mat, host, a.threads) != ref for _ in range(a.repeats))
        bad += mism
        rows.append({"case": "%s-%s" % (mat, host), "threads": a.threads, "repeats": a.repeats,
                     "mismatches": mism, "seconds": round(time.time() - t0, 1)})
        print(json.dumps(rows[-1]), flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps({"rows": rows, "mismatches": bad}, indent=1))
    print("TOTAL mismatches:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
