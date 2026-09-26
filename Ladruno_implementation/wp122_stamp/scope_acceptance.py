#!/usr/bin/env python3
"""WP-122 acceptance: the stamp is what puts a file in the quirk lint's scope.

Plants the same L2 incident (a process-wide singleton `instance()` with no reset
from Domain::clearAll and no waiver) into one newly stamped file in two trees:
  before = `git archive <before-ref>`  (file unstamped)  -> lint MUST stay silent on it
  after  = `git archive <after-ref>`   (file stamped)    -> lint MUST flag it
Also checks the unplanted after-tree is clean. Trees go to a temp dir; the
worktree is never modified.

    python scope_acceptance.py [--before origin/ladruno] [--after HEAD]
"""
from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGET = "SRC/utility/LadrunoThreads.cpp"
PLANT = (
    "\n// WP-122 scope-acceptance plant (never committed)\n"
    "class Wp122PlantRegistry {\n"
    "public:\n"
    "  static Wp122PlantRegistry &instance() { static Wp122PlantRegistry r; return r; }\n"
    "  int n = 0;\n"
    "};\n"
)


def extract(ref, dest):
    data = subprocess.run(["git", "archive", ref, "SRC", "ci", "Ladruno_implementation",
                           ".claude/skills"], cwd=REPO, capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as tf:
        tf.extractall(dest)


def lint(root):
    r = subprocess.run([sys.executable, str(REPO / "ci" / "check_quirk_patterns.py"), "--root", str(root)],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", default="origin/ladruno")
    ap.add_argument("--after", default="HEAD")
    a = ap.parse_args()
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for label, ref, plant, expect in (("after, unplanted", a.after, False, "clean"),
                                          ("before + plant", a.before, True, "silent on target"),
                                          ("after + plant", a.after, True, "flags target")):
            d = Path(tmp) / label.replace(" ", "_").replace(",", "").replace("+", "p")
            d.mkdir()
            extract(ref, d)
            if plant:
                with open(d / TARGET, "a", encoding="utf-8") as fh:
                    fh.write(PLANT)
            rc, out = lint(d)
            hit = "LadrunoThreads" in out
            if expect == "clean":
                good = rc == 0
            elif expect == "silent on target":
                good = not hit
            else:
                good = hit and rc != 0
            ok &= good
            first = out.strip().splitlines()[-1] if out.strip() else ""
            print(f"[{'PASS' if good else 'FAIL'}] {label:18s} ({ref}): expect {expect}; "
                  f"rc={rc}, target flagged={hit} | {first}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
