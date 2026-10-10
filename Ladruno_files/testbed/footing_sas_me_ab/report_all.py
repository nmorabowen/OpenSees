"""WP-138: run replay_report.py over every replay CSV of the given legs that
has oracle / ME / SAS outputs in replay_out/.  python report_all.py LEG [LEG ...]"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for leg in sys.argv[1:]:
    rd = os.path.join(HERE, "runs", leg)
    for c in sorted(glob.glob(os.path.join(rd, "replay", "*.csv"))):
        base = os.path.splitext(os.path.basename(c))[0]
        od = os.path.join(rd, "replay_out")
        args = [sys.executable, os.path.join(HERE, "replay_report.py"), "--csv", c]
        o = os.path.join(od, f"{base}.oracle.json")
        if os.path.exists(o):
            args += ["--oracle", o]
        for lab in ("ME", "SAS"):
            j = os.path.join(od, f"{base}.{lab}.json")
            if os.path.exists(j):
                args += ["--cxx", j, "--label", lab]
        if "--cxx" not in args:
            continue
        args += ["--out", os.path.join(od, f"{base}.table.csv")]
        r = subprocess.run(args, capture_output=True, text=True)
        tail = r.stdout.strip().splitlines()[-3:]
        print(f"== {leg}/{base}\n" + "\n".join(tail) + (r.stderr[-400:] if r.returncode else ""))
