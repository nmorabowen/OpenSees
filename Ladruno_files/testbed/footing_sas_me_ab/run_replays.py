"""WP-138: batch the C++ replay (and optionally the oracle) over replay CSVs.

    python run_replays.py --bin A|B --label A [--scheme 1] [--extra "..."] [--oracle]
        CSV [CSV ...]

Writes <csv-dir>/../replay_out/<csvname>.<label>.json (and .oracle.json).
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY312 = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
PY311 = r"C:\Users\nmora\AppData\Local\Programs\Python\Python311\python.exe"
sys.path.insert(0, HERE)
BINS = {
    "A": r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\tims-implementation-review-3733c6\dist\bin",
    "B": os.path.join(HERE, "binB_5c8dcd0e0"),   # snapshot, see launch.py
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="A")
    ap.add_argument("--label", default=None)
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--extra", default="")
    ap.add_argument("--maxsub", type=int, default=2000)
    ap.add_argument("--oracle", action="store_true")
    ap.add_argument("--tolr", type=float, default=1.0e-7)
    ap.add_argument("--ref-dir", default=os.environ.get("SANISAND_REF_DIR", ""))
    ap.add_argument("csvs", nargs="+")
    a = ap.parse_args()
    label = a.label or a.bin
    env = dict(os.environ)
    env["FOOTING_BIN"] = BINS.get(a.bin, a.bin)
    env["MKL_NUM_THREADS"] = env["OMP_NUM_THREADS"] = "1"
    for c in a.csvs:
        od = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(c))), "replay_out")
        os.makedirs(od, exist_ok=True)
        base = os.path.splitext(os.path.basename(c))[0]
        oj = os.path.join(od, f"{base}.{label}.json")
        r = subprocess.run([PY312, "-S", os.path.join(HERE, "replay_cxx.py"), "--csv", c,
                            "--out", oj, "--scheme", str(a.scheme), "--extra", a.extra,
                            "--maxsub", str(a.maxsub), "--tolr", str(a.tolr)], env=env, capture_output=True,
                           text=True, timeout=7200)
        last = [l for l in r.stdout.splitlines() if "rows" in l]
        print(f"{base} [{label}]: {last[-1] if last else r.stdout[-500:] + r.stderr[-800:]}",
              flush=True)
        if a.oracle:
            oo = os.path.join(od, f"{base}.oracle.json")
            r = subprocess.run([PY311, os.path.join(HERE, "replay_oracle.py"), "--csv", c,
                                "--out", oo, "--ref-dir", a.ref_dir], capture_output=True,
                               text=True, timeout=36000)
            print(f"{base} [oracle]: rc {r.returncode} {r.stderr[-400:]}", flush=True)


if __name__ == "__main__":
    main()
