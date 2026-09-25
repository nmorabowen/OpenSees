#!/usr/bin/env python3
"""WP-120 scope check: which SRC files are fork-authored, and does the stamp say so?

  fork-added  = C/C++ files under SRC on --ref that exist neither on upstream/master
                nor at the merge-base of the two (i.e. the fork created them)
  stamped     = files carrying LADRUNO-HEADER-START (what the quirk lint scans)
  GLOBS       = Ladruno_scripts/stamp_headers.py GLOBS (what stamping maintains)

Prints the drift between the three and writes the unstamped fork-added list to
--out (default: unstamped_fork_files.txt next to this script).

    python inventory.py [--ref origin/ladruno] [--upstream upstream/master]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import ROOT, SUFFIXES, rel, split_fork_vanilla  # noqa: E402

sys.path.insert(0, str(ROOT / "Ladruno_scripts"))
import stamp_headers  # noqa: E402


def tree(ref):
    out = subprocess.run(["git", "ls-tree", "-r", "--name-only", ref, "SRC"], cwd=ROOT,
                         capture_output=True, text=True).stdout.split()
    return {p for p in out if Path(p).suffix in SUFFIXES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="origin/ladruno")
    ap.add_argument("--upstream", default="upstream/master")
    ap.add_argument("--out", default=str(HERE / "unstamped_fork_files.txt"))
    a = ap.parse_args()
    mb = subprocess.run(["git", "merge-base", a.ref, a.upstream], cwd=ROOT,
                        capture_output=True, text=True).stdout.strip()
    added = tree(a.ref) - tree(a.upstream) - tree(mb)
    fork, _ = split_fork_vanilla()
    stamped = {rel(p) for p in fork}
    globs = {rel(p) for p in stamp_headers.authored_files()}
    dead_globs = [g for g in stamp_headers.GLOBS if not list(ROOT.glob(g))]
    unstamped = sorted(added - stamped)
    lines = {p: sum(1 for _ in open(ROOT / p, encoding="utf-8", errors="replace")) for p in unstamped}
    print(f"merge-base {mb[:9]}; fork-added {len(added)}; stamped {len(stamped)}; GLOBS files {len(globs)}")
    print(f"fork-added but UNSTAMPED (invisible to the quirk lint): {len(unstamped)} files, "
          f"{sum(lines.values())} lines")
    for p in unstamped:
        print(f"   {lines[p]:5d}  {p}")
    print(f"stamped but not fork-added (upstream files the fork stamped): {sorted(stamped - added)}")
    print(f"stamped but missing from GLOBS (a re-stamp would not maintain them): {sorted(stamped - globs)}")
    print(f"GLOBS entries matching no file: {dead_globs}")
    Path(a.out).write_text("\n".join(unstamped) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
