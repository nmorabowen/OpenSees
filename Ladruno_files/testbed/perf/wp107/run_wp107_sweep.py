"""WP-107 sweep driver -- identity gate + speed-up table.

Runs `wp107_strip_bench.py` as a SEPARATE PROCESS per configuration (a fresh
interpreter per run, so nothing carries over) and reports:

  * BIT-IDENTITY: every thread count's load-settlement curve against the
    1-thread baseline, as an exact string compare of the repr()'d doubles plus
    the max absolute difference. ADR-75b section 7's correctness protocol item 4
    requires N >= 10 repeats per thread count and that each configuration
    reproduce ITSELF, not merely the baseline -- `--repeats` does that.
  * WALL: per-step wall and the speed-up against 1 thread.

Usage (from the worktree root):

    python3.12 -u Ladruno_files/testbed/perf/wp107/run_wp107_sweep.py \\
        --outdir Ladruno_files/testbed/perf/wp107/run1 \\
        --h 0.2 --steps 10 --threads 1 4 8 --repeats 3
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "wp107_strip_bench.py")


def run(outdir, tag, threads, extra, env):
    out = os.path.join(outdir, "curve_%s.csv" % tag)
    log = os.path.join(outdir, "log_%s.txt" % tag)
    cmd = [sys.executable, "-u", BENCH, "--threads", str(threads), "--out", out] + extra
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    txt = open(log, encoding="utf-8", errors="replace").read()
    wall = None
    for line in txt.splitlines():
        if line.startswith("WALLLINE"):
            for kv in line.split():
                if kv.startswith("perstep="):
                    wall = float(kv.split("=", 1)[1])
    return rc, out, log, wall, txt


def curve(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def maxabs(a, b):
    """Max |diff| over the two numeric columns of two curve files."""
    ra = [l.split(",") for l in a.strip().splitlines()[1:]]
    rb = [l.split(",") for l in b.strip().splitlines()[1:]]
    m = 0.0
    for x, y in zip(ra, rb):
        for i in (2, 3):
            m = max(m, abs(float(x[i]) - float(y[i])))
    return m, (len(ra), len(rb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--h", type=float, default=0.2)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 8])
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--mat", default="sanisand")
    ap.add_argument("--tangent", type=int, default=0)
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--system", default="BandGeneral")
    ap.add_argument("--implex", action="store_true")
    ap.add_argument("--ds", type=float, default=0.002)
    ap.add_argument("--load", default="oedometer")
    ap.add_argument("--label", default="")
    ap.add_argument("--pyd", default=None, help="dist/bin to put on PYTHONPATH")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    env = dict(os.environ)
    env["MKL_NUM_THREADS"] = "1"
    if args.pyd:
        env["PYTHONPATH"] = args.pyd

    extra = ["--h", str(args.h), "--steps", str(args.steps),
             "--mat", args.mat, "--tangent", str(args.tangent),
             "--scheme", str(args.scheme), "--system", args.system,
             "--ds", str(args.ds), "--load", args.load]
    if args.implex:
        extra.append("--implex")

    results = []
    base_txt = None
    for t in args.threads:
        for r in range(args.repeats):
            tag = "%s%s_t%d_r%d" % (args.label, args.mat, t, r)
            rc, out, log, wall, txt = run(args.outdir, tag, t, extra, env)
            ok = (rc == 0 and os.path.exists(out))
            c = curve(out) if ok else None
            if base_txt is None and ok:
                base_txt = c
                base_tag = tag
            ident = None
            diff = None
            if ok and base_txt is not None:
                ident = (c == base_txt)
                try:
                    diff, _ = maxabs(base_txt, c)
                except Exception:
                    diff = float("nan")
            refused = ("running the element loop SERIAL" in txt
                       or "running SERIAL" in txt)
            results.append((tag, t, r, rc, wall, ident, diff, refused))
            print("%-28s threads=%d rep=%d rc=%d perstep=%s identical=%s maxdiff=%s refused=%s"
                  % (tag, t, r, rc, ("%.5f" % wall) if wall else "n/a",
                     ident, ("%.3e" % diff) if diff is not None else "n/a", refused))

    # --- summary table
    print()
    print("| threads | per-step wall (s, min over repeats) | speed-up | curve bit-identical | max abs diff | refused? |")
    print("|---|---|---|---|---|---|")
    base = None
    for t in args.threads:
        rows = [x for x in results if x[1] == t and x[4] is not None]
        if not rows:
            print("| %d | FAILED | - | - | - | - |" % t)
            continue
        w = min(x[4] for x in rows)
        if base is None:
            base = w
        allid = all(x[5] for x in rows)
        mx = max((x[6] or 0.0) for x in rows)
        ref = any(x[7] for x in rows)
        print("| %d | %.5f | %.2fx | %s | %.3e | %s |"
              % (t, w, base / w, "YES" if allid else "**NO**", mx, "yes" if ref else "no"))


if __name__ == "__main__":
    main()
