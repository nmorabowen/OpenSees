"""WP-138 launcher: run one leg of footing_ab.py in a SUBPROCESS with a hard
wall-clock timeout, a pinned thread budget and a plain slurm-like log at
<out>/logs/log.log.

    python launch.py --leg <name> --bin A|B|<dir> --timeout <s> -- <footing_ab args>

Threads: MKL_NUM_THREADS = OMP_NUM_THREADS = 1 by default. The SANISAND update
is serial anyway (WP-107 refuses to thread it) and PARDISO is ~7 % of a step on
this deck (intake section 1.2), so one thread costs little and buys bit-for-bit
reproducible runs (75c trap 7) -- which is what an A/B comparison needs, since
the intake measured a 30 % run-to-run spread of the wall with threads (1.6).
MKL_CBWR=COMPATIBLE pins the same MKL code path for both binaries.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
BINS = {
    "A": r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\tims-implementation-review-3733c6\dist\bin",
    # snapshot of agent-a6b3c48d38c9ca3d3\dist\bin taken 2026-09-27 15:10 (that
    # folder is being rebuilt); the pyd reports ladrunoBuild cdf43685f, which is
    # SRC-identical to WP-129 head 5c8dcd0e0 (the later commits touch docs only).
    # B_sasme_default was started from the original folder before the snapshot.
    "B": os.path.join(HERE, "binB_beb6d8333"),   # WP-129 review-fixed build (verified ladrunoBuild beb6d8333)
    "B_cdf": os.path.join(HERE, "binB_5c8dcd0e0"),   # pre-review build cdf43685f (provisional)
    "B_orig": r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\agent-a6b3c48d38c9ca3d3\dist\bin",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", required=True)
    ap.add_argument("--bin", default="A")
    ap.add_argument("--timeout", type=float, default=40000.0)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("rest", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    rest = a.rest[1:] if a.rest and a.rest[0] == "--" else a.rest
    out = os.path.join(HERE, "runs", a.leg)
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)
    env = dict(os.environ)
    env["FOOTING_BIN"] = BINS.get(a.bin, a.bin)
    env["MKL_NUM_THREADS"] = env["OMP_NUM_THREADS"] = str(a.threads)
    env.setdefault("MKL_CBWR", "COMPATIBLE")
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    cmd = [PY, "-S", "-u", os.path.join(HERE, "footing_ab.py"), "--out", out] + rest
    logp = os.path.join(out, "logs", "log.log")
    with open(logp, "a", buffering=1) as lg:
        lg.write(f"==== {time.strftime('%Y-%m-%d %H:%M:%S')} leg {a.leg} bin "
                 f"{env['FOOTING_BIN']} timeout {a.timeout}s threads {a.threads}\n")
        lg.write("==== cmd " + " ".join(cmd) + "\n")
        t0 = time.time()
        try:
            rc = subprocess.run(cmd, env=env, stdout=lg, stderr=subprocess.STDOUT,
                                timeout=a.timeout).returncode
        except subprocess.TimeoutExpired:
            rc = "TIMEOUT"
        lg.write(f"==== {time.strftime('%Y-%m-%d %H:%M:%S')} exit {rc} after "
                 f"{time.time()-t0:.1f}s\n")
    print(f"leg {a.leg} exit {rc}")


if __name__ == "__main__":
    main()
