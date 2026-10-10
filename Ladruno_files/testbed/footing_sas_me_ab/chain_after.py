"""WP-138: start a leg only after another leg's launcher has written its exit
line AND its footing_ab log shows MODE (a real end, not a kill).

    python chain_after.py WAIT_LEG -- <launch.py args>
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
wait_leg = sys.argv[1]
rest = sys.argv[sys.argv.index("--") + 1:]
log = os.path.join(HERE, "runs", wait_leg, "logs", "log.log")
while True:
    try:
        txt = open(log).read()
    except OSError:
        txt = ""
    if "MODE = " in txt:
        break
    if "==== " in txt and " exit " in txt.splitlines()[-1]:
        print(f"{wait_leg} exited WITHOUT a MODE line (killed?); not starting", flush=True)
        sys.exit(1)
    time.sleep(120)
subprocess.run([sys.executable, "-S", os.path.join(HERE, "launch.py")] + rest)
