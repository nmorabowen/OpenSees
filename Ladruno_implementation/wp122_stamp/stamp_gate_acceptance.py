#!/usr/bin/env python3
"""WP-122 acceptance for the CI step `stamp_headers.py --check`.

In `git archive` copies of <ref> (never the worktree):
  unchanged                               -> --check MUST exit 0
  stamp block removed from a GLOBS file   -> --check MUST exit 1 and name the file
  stamp block edited (stale art/credit)   -> --check MUST exit 1 and name the file
The script derives ROOT from its own location, so each tree runs its own copy.

    python stamp_gate_acceptance.py [--ref HEAD]
"""
from __future__ import annotations

import argparse
import io
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGET = "SRC/element/LadrunoMassCache.h"      # one of the 31 WP-122 stamped


def extract(ref, dest):
    data = subprocess.run(["git", "archive", ref, "SRC", "Ladruno_scripts"], cwd=REPO,
                          capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as tf:
        tf.extractall(dest)


def mutate(root, how):
    p = root / TARGET
    raw = p.read_bytes().decode("utf-8")
    if how == "remove":
        new = re.sub(r"// LADRUNO-HEADER-START.*?// LADRUNO-HEADER-END\r?\n", "", raw, count=1, flags=re.S)
    else:
        new = raw.replace("Ladruno — a research fork of OpenSees", "Ladruno - an old credit line", 1)
    assert new != raw, how
    p.write_bytes(new.encode("utf-8"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD")
    a = ap.parse_args()
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for label, how, want_rc in (("unchanged", None, 0), ("stamp removed", "remove", 1),
                                    ("stamp stale", "stale", 1)):
            d = Path(tmp) / label.replace(" ", "_")
            d.mkdir()
            extract(a.ref, d)
            if how:
                mutate(d, how)
            r = subprocess.run([sys.executable, str(d / "Ladruno_scripts" / "stamp_headers.py"), "--check"],
                               capture_output=True, text=True, encoding="utf-8", errors="replace")
            named = TARGET in r.stdout
            good = r.returncode == want_rc and (named if want_rc else True)
            ok &= good
            last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-120:]
            print(f"[{'PASS' if good else 'FAIL'}] {label:14s}: rc={r.returncode} (want {want_rc}), "
                  f"target named={named} | {last}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
