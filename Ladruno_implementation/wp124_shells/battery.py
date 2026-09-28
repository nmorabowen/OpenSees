#!/usr/bin/env python3
"""WP-124 regression battery: every tests/*.py that builds one of the eight continuum shells,
each file in its OWN process (a crash cannot hide the rest), against a chosen dist/bin.
Writes {nodeid: outcome} as JSON; --compare prints the tests whose outcome changed.

    <py3.12> battery.py DIST_BIN OUT.json [--files test_a.py test_b.py ...]
    <py3.12> battery.py --compare BEFORE.json AFTER.json

The binary under test reaches EVERY interpreter, not just the parent:
  * the parent runs `python -S` with DIST_BIN wired in-process and asserts opensees.__file__;
  * CHILD interpreters a test spawns itself (`sys.executable ...`, no -S -- e.g.
    test_ladruno_sanisand_flip_determinism.py `_run_child`, importing via tests/_testbed) do NOT
    inherit in-process sys.path / add_dll_directory, so the pin travels in the environment:
    PYTHONPATH = DIST_BIN;tests;site-packages and PATH = DIST_BIN first. Without it the children
    found no `opensees` at all and every such test "failed" identically on both binaries.
  * every file runs with cwd = tests/ (tests such as test_adr94_matrix.py open paths relative
    to it);
  * a tracked file a test REWRITES (test_adr94_matrix.py regenerates
    Ladruno_implementation/_adr94_matrix.md) is restored after that file, if it was clean
    before the run, and reported under "<file>::<restored>".
Node ids are relative to tests/ (`test_x.py::name`), so JSONs from before this fix (keys
`tests\\test_x.py::name`) do not compare key-for-key with new ones.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
SITE = str(Path(PY).parent / "Lib" / "site-packages")
ELEMS = re.compile(r"['\"](LadrunoQuad|LadrunoCST|LadrunoLST|LadrunoCSTPair|LadrunoBrick|"
                   r"LadrunoBrick20|BezierTri6|BezierTet10)['\"]")
RUNNER = r'''
import os, sys
DIST, TESTS = sys.argv[1], sys.argv[2]
assert sys.flags.no_site
os.add_dll_directory(DIST); sys.path.insert(0, DIST); sys.path.insert(0, TESTS)
sys.path.append(os.environ["BATTERY_SITE"])
import opensees
assert os.path.normcase(opensees.__file__) == os.path.normcase(os.path.join(DIST, "opensees.pyd")), opensees.__file__
import pytest
sys.exit(pytest.main(["-p", "no:cacheprovider", "-q", "-rA", "--tb=line", *sys.argv[3:]]))
'''


def files():
    return sorted(p.name for p in TESTS.glob("test_*.py")
                  if ELEMS.search(p.read_text(encoding="utf-8", errors="replace")))


def dirty_tracked():
    r = subprocess.run(["git", "diff", "--name-only"], cwd=ROOT, capture_output=True, text=True)
    return set(r.stdout.split())


def child_env(dist):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([dist, str(TESTS), SITE])
    env["PATH"] = dist + os.pathsep + env.get("PATH", "")
    env["BATTERY_SITE"] = SITE
    env.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    return env


def run(dist, out, only=None):
    dist = str(Path(dist).resolve())
    runner = ROOT / "build" / "wp124_battery_runner.py"
    runner.parent.mkdir(exist_ok=True)
    runner.write_text(RUNNER, encoding="utf-8")
    env = child_env(dist)
    before = dirty_tracked()
    res = {}
    for f in (only or files()):
        r = subprocess.run([PY, "-S", str(runner), dist, str(TESTS), f], cwd=TESTS, env=env,
                           capture_output=True, text=True, errors="replace", timeout=3600)
        got = 0
        for ln in r.stdout.splitlines():
            m = re.match(r"^(PASSED|FAILED|ERROR|SKIPPED|XFAIL|XPASS) (\S+)", ln)
            if m:
                res[m.group(2)] = m.group(1)
                got += 1
        if got == 0 or r.returncode not in (0, 1, 5):
            res[f"{f}::<file>"] = f"RC{r.returncode}" + ("" if got else " (no results)")
        rewritten = sorted(dirty_tracked() - before)
        if rewritten:                                   # a test rewrote a tracked file: put it back
            subprocess.run(["git", "checkout", "--", *rewritten], cwd=ROOT, check=True)
            res[f"{f}::<restored>"] = ",".join(rewritten)
        print(f"{f:60s} rc={r.returncode} tests={got}" + (f" restored={rewritten}" if rewritten else ""),
              flush=True)
    json.dump(res, open(out, "w", encoding="utf-8"), indent=0, sort_keys=True)
    counts = {}
    for k, v in res.items():
        if not k.endswith("::<restored>"):
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
        only = sys.argv[sys.argv.index("--files") + 1:] if "--files" in sys.argv else None
        run(sys.argv[1], sys.argv[2], only)
