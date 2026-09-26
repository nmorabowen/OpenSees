#!/usr/bin/env python3
"""WP-122 acceptance for the CI step `stamp_headers.py --check`.

In `git archive` copies of <ref> (never the worktree):
  unchanged                                  -> --check MUST exit 0
  stamp block removed from a GLOBS file      -> --check MUST exit 1 and name the file
  stamp block edited (stale art/credit)      -> --check MUST exit 1 and name the file
  new Ladruno-named source not in GLOBS      -> --check MUST exit 1 and name the file
The script derives ROOT from its own location, so each tree runs its own copy.

History (the incident): the Ladruno-path rule, evaluated on <before> with that
tree's own GLOBS, must flag exactly the Ladruno-named files among WP-120 R1's 31
unstamped ones (14), plus any stamped-but-not-in-GLOBS Ladruno-named files.

    python stamp_gate_acceptance.py [--ref HEAD] [--before origin/ladruno]
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGET = "SRC/element/LadrunoMassCache.h"      # one of the 31 WP-122 stamped
PLANT = "SRC/element/LadrunoWp122Plant.h"      # a new fork file nobody added to GLOBS
STAMP = "LADRUNO-HEADER-START"


def extract(ref, dest):
    data = subprocess.run(["git", "archive", ref, "SRC", "Ladruno_scripts"], cwd=REPO,
                          capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as tf:
        tf.extractall(dest)


def mutate(root, how):
    if how == "plant":
        (root / PLANT).write_text("#pragma once\n// WP-122 acceptance plant\nint wp122Plant();\n",
                                  encoding="utf-8")
        return PLANT
    p = root / TARGET
    raw = p.read_bytes().decode("utf-8")
    if how == "remove":
        new = re.sub(r"// LADRUNO-HEADER-START.*?// LADRUNO-HEADER-END\r?\n", "", raw, count=1, flags=re.S)
    else:
        new = raw.replace("Ladruno — a research fork of OpenSees", "Ladruno - an old credit line", 1)
    assert new != raw, how
    p.write_bytes(new.encode("utf-8"))
    return TARGET


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD")
    ap.add_argument("--before", default="origin/ladruno")
    a = ap.parse_args()
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for label, how, want_rc in (("unchanged", None, 0), ("stamp removed", "remove", 1),
                                    ("stamp stale", "stale", 1), ("unlisted plant", "plant", 1)):
            d = Path(tmp) / label.replace(" ", "_")
            d.mkdir()
            extract(a.ref, d)
            named_file = mutate(d, how) if how else None
            r = subprocess.run([sys.executable, str(d / "Ladruno_scripts" / "stamp_headers.py"), "--check"],
                               capture_output=True, text=True, encoding="utf-8", errors="replace")
            named = bool(named_file) and named_file in r.stdout
            good = r.returncode == want_rc and (named if want_rc else True)
            ok &= good
            last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-120:]
            print(f"[{'PASS' if good else 'FAIL'}] {label:14s}: rc={r.returncode} (want {want_rc}), "
                  f"file named={named} | {last}")

        # history: the new rule on the pre-WP-122 tree, with that tree's own GLOBS
        d = Path(tmp) / "before"
        d.mkdir()
        extract(a.before, d)
        old = load(d / "Ladruno_scripts" / "stamp_headers.py", "sh_before")
        new = load(REPO / "Ladruno_scripts" / "stamp_headers.py", "sh_after")
        flagged = new.ladruno_named_outside_globs(d, old.authored_files())
        rel = [p.relative_to(d).as_posix() for p in flagged]
        unstamped = [x for x in rel if STAMP not in (d / x).read_text(encoding="utf-8", errors="replace")]
        stamped = [x for x in rel if x not in unstamped]
        good = len(unstamped) == 14
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] history ({a.before}): rule flags {len(rel)} files -- "
              f"{len(unstamped)} unstamped (want 14 of WP-120 R1's 31), {len(stamped)} stamped but "
              f"outside GLOBS")
        for x in unstamped:
            print("       unstamped:", x)
        for x in stamped:
            print("       outside GLOBS:", x)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
