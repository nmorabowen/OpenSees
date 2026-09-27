#!/usr/bin/env python3
"""WP-124 regression battery: every tests/*.py that builds one of the eight continuum shells,
each file in its OWN process (a crash cannot hide the rest), against a chosen dist/bin.
Writes {nodeid: outcome} as JSON; --compare prints the tests whose outcome changed.

    <py3.12> battery.py DIST_BIN OUT.json          (DIST_BIN: a dist/bin folder)
    <py3.12> battery.py --compare BEFORE.json AFTER.json
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
ELEMS = re.compile(r"['\"](LadrunoQuad|LadrunoCST|LadrunoLST|LadrunoCSTPair|LadrunoBrick|"
                   r"LadrunoBrick20|BezierTri6|BezierTet10)['\"]")
RUNNER = r'''
import os, sys
DIST, ROOT = sys.argv[1], sys.argv[2]
assert sys.flags.no_site
os.add_dll_directory(DIST); sys.path.insert(0, DIST); sys.path.insert(0, os.path.join(ROOT, "tests"))
sys.path.append(os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages"))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees
assert os.path.normcase(opensees.__file__) == os.path.normcase(os.path.join(DIST, "opensees.pyd")), opensees.__file__
os.chdir(ROOT)
import pytest
sys.exit(pytest.main(["-p", "no:cacheprovider", "-q", "-rA", "--tb=no", *sys.argv[3:]]))
'''


def files():
    return sorted(str(p.relative_to(ROOT)).replace("\\", "/") for p in (ROOT / "tests").glob("test_*.py")
                  if ELEMS.search(p.read_text(encoding="utf-8", errors="replace")))


def run(dist, out):
    runner = ROOT / "build" / "wp124_battery_runner.py"
    runner.write_text(RUNNER, encoding="utf-8")
    res = {}
    for f in files():
        r = subprocess.run([PY, "-S", str(runner), dist, str(ROOT), f], cwd=ROOT,
                           capture_output=True, text=True, errors="replace", timeout=3600)
        got = 0
        for ln in r.stdout.splitlines():
            m = re.match(r"^(PASSED|FAILED|ERROR|SKIPPED|XFAIL|XPASS) (\S+)", ln)
            if m:
                res[m.group(2)] = m.group(1)
                got += 1
        if got == 0 or r.returncode not in (0, 1, 5):
            res[f"{f}::<file>"] = f"RC{r.returncode}" + ("" if got else " (no results)")
        print(f"{f:60s} rc={r.returncode} tests={got}", flush=True)
    json.dump(res, open(out, "w", encoding="utf-8"), indent=0, sort_keys=True)
    counts = {}
    for v in res.values():
        counts[v] = counts.get(v, 0) + 1
    print(out, counts)


def compare(a, b):
    A = json.load(open(a, encoding="utf-8"))
    B = json.load(open(b, encoding="utf-8"))
    changed = [(k, A.get(k, "-"), B.get(k, "-")) for k in sorted(set(A) | set(B)) if A.get(k) != B.get(k)]
    print(f"{len(A)} / {len(B)} results; {len(changed)} changed")
    for k, x, y in changed:
        print(f"  {x:8s} -> {y:8s} {k}")


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        run(sys.argv[1], sys.argv[2])
