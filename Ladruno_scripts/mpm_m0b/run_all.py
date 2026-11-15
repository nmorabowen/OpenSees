"""WP-232g (M0b): run every measurement script against one build, in sequence.

    python3.12 Ladruno_scripts/mpm_m0b/run_all.py [--pyd PATH] [--smoke] [--only m1,m3]
                                                   [--out DIR] [--repeats N]

Each script runs in its own process and starts a fresh child per run; results land in
<out>/<script>/ (default results/run_all_<stamp>/). The orchestrator launches the full runs
(the 232g brief): nothing else heavy may run on the machine meanwhile.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C       # noqa: E402

SCRIPTS = {"m1": "m1_rebuild.py", "m2": "m2_explicit.py", "m3": "m3_memory.py", "m4": "m4_cantilever.py"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pyd", default=str(C.DEFAULT_PYD))
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--only", default=None, help="comma-separated subset of m1,m2,m3,m4")
    ap.add_argument("--out", default=None)
    ap.add_argument("--repeats", type=int, default=None)
    args = ap.parse_args()
    root = Path(args.out) if args.out else \
        C.RESULTS / f"run_all{'_smoke' if args.smoke else ''}_{_dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
    root.mkdir(parents=True, exist_ok=True)
    which = [w.strip() for w in args.only.split(",")] if args.only else list(SCRIPTS)
    index = {"pyd": str(Path(args.pyd).resolve()), "smoke": args.smoke, "scripts": {}}
    rc_all = 0
    for w in which:
        cmd = [sys.executable, str(HERE / SCRIPTS[w]), "--pyd", args.pyd, "--out", str(root / w)]
        if args.smoke:
            cmd.append("--smoke")
        if args.repeats:
            cmd += ["--repeats", str(args.repeats)]
        print(f"[run_all] {' '.join(cmd)}", flush=True)
        t0 = time.time()
        rc = subprocess.call(cmd, stdin=subprocess.DEVNULL)
        index["scripts"][w] = {"rc": rc, "wall_s": time.time() - t0, "dir": str(root / w)}
        rc_all = rc_all or rc
        (root / "index.json").write_text(json.dumps(index, indent=1), encoding="utf-8")
    print(f"[run_all] done: {root}")
    return rc_all


if __name__ == "__main__":
    sys.exit(main())
