"""WP-144 P3 launcher (COPY of ref_ladder/launch.py; edits: engine dirs -> ~/wp144/bin/{plain,pdl_seq}, TMPDIR -> ~/wp144/tmp, runs footing_ab_ns_k5.py).
WP-138 launcher, ESMERALDA COPY (Linux): run one leg of footing_ab.py in a
SUBPROCESS with a pinned thread budget and a plain slurm-like log at
<HERE>/runs/<leg>/logs/log.log. Meant to be the command of an
`srun --exclusive` step, so it runs ON a compute node.

    python launch.py --leg <name> [--bin E|<dir>] [--timeout s] -- <footing_ab args>

Threads: MKL_NUM_THREADS = OMP_NUM_THREADS = 1 (bit-for-bit reproducible runs).
The build links the SEQUENTIAL MKL 3-layer (spack oneMKL 2024.2.2) for
PARDISO + BLAS/LAPACK; MKL_CBWR=COMPATIBLE pins the MKL code path (on this AMD
Zen3 host only AUTO and COMPATIBLE are accepted anyway).
Python: the conan venv's 3.10 (has numpy), NOT -S (the Windows -S trap does
not exist here); footing_ab.py asserts the module it loaded is <bin>/opensees.so.
"""
import argparse
import os
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.expanduser("~/wp144/venv/bin/python")      # 3.10 + numpy (same ABI as the module)
MKL = ("/mnt/nfshare/software/spack/opt/spack/linux-zen3/intel-oneapi-mkl-2024.2.2-"
       "4olf4bnrh4yzaqbhq742doqvnim6yynl/mkl/2024.2/lib")
BINS = {"plain": os.path.expanduser("~/wp144/bin/plain"),          # build/Release (no Pardiso)
        "pdl": os.path.expanduser("~/wp144/bin/pdl_seq"),           # build/pdl_seq (Pardiso, sequential MKL)
        "E": os.path.expanduser("~/wp144/bin/pdl_seq")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", required=True)
    ap.add_argument("--bin", default="pdl")
    ap.add_argument("--timeout", type=float, default=3.0e6)
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
    env["LD_LIBRARY_PATH"] = MKL + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    env["TMPDIR"] = os.path.expanduser("~/wp144/tmp")   # never write into the other WPs' trees
    cmd = [PY, "-u", os.path.join(HERE, "footing_ab_ns_k5.py"), "--out", out] + rest
    logp = os.path.join(out, "logs", "log.log")
    with open(logp, "a", buffering=1) as lg:
        lg.write(f"==== {time.strftime('%Y-%m-%d %H:%M:%S')} leg {a.leg} host "
                 f"{socket.gethostname()} slurm_job {os.environ.get('SLURM_JOB_ID')} bin "
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
