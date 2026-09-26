#!/usr/bin/env python3
"""WP-122 acceptance for the CI step `stamp_headers.py --check`.

Each case runs in its own `git archive` copy of <ref> (never the worktree); the
script derives ROOT from its own location, so each tree runs its own copy.

  unchanged                                  -> exit 0
  stamp block removed from a GLOBS file      -> exit 1, file named
  stamp block edited (stale art/credit)      -> exit 1, file named
  new Ladruno-named source, not in GLOBS     -> exit 1, file named
  new NEUTRAL-named source, not in GLOBS     -> exit 1, file named (upstream-manifest rule)
  stale manifest (an upstream file dropped)  -> exit 1, file named + refresh hint
  NOT_STAMPED exemption for the neutral file -> exit 0
  NOT_STAMPED exemption for a GLOBS file     -> exit 1 (stale exemption)
  manifest file deleted                      -> exit 1

History (the incident): both GLOBS-coverage rules, evaluated on <before> with that
tree's own GLOBS and this branch's manifest, must flag all 31 files WP-120 R1 found
unstamped, and nothing else unstamped.

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
TARGET = "SRC/element/LadrunoMassCache.h"                 # one of the 31 WP-122 stamped
PLANT_L = "SRC/element/LadrunoWp122Plant.h"               # names itself
PLANT_N = "SRC/analysis/integrator/Wp122NeutralPlant.cpp"  # does not (the CriticalTimeStep shape)
UPSTREAM_FILE = "SRC/element/fourNodeQuad/FourNodeQuad.cpp"
MANIFEST = "Ladruno_scripts/upstream_src_manifest.txt"
SCRIPT = "Ladruno_scripts/stamp_headers.py"
STAMP = "LADRUNO-HEADER-START"
EXEMPT_DECL = "NOT_STAMPED: dict[str, str] = {}"


def extract(ref, dest):
    data = subprocess.run(["git", "archive", ref, "SRC", "Ladruno_scripts"], cwd=REPO,
                          capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as tf:
        tf.extractall(dest)


def edit(path, old, new):
    raw = path.read_bytes().decode("utf-8")
    assert raw.count(old) == 1, (path, old[:40])
    path.write_bytes(raw.replace(old, new).encode("utf-8"))


def mutate(root, how):
    """Apply one mutation; return the path --check must name (or None)."""
    if how == "remove":
        p = root / TARGET
        raw = p.read_bytes().decode("utf-8")
        new = re.sub(r"// LADRUNO-HEADER-START.*?// LADRUNO-HEADER-END\r?\n", "", raw, count=1, flags=re.S)
        assert new != raw
        p.write_bytes(new.encode("utf-8"))
        return TARGET
    if how == "stale":
        edit(root / TARGET, "Ladruno — a research fork of OpenSees", "Ladruno - an old credit line")
        return TARGET
    if how in ("plant_l", "plant_n", "exempt_ok"):
        rel = PLANT_L if how == "plant_l" else PLANT_N
        (root / rel).write_text("#pragma once\n// WP-122 acceptance plant\nint wp122Plant();\n", encoding="utf-8")
        if how == "exempt_ok":
            edit(root / SCRIPT, EXEMPT_DECL,
                 'NOT_STAMPED: dict[str, str] = {"%s": "acceptance: vendored third-party file"}' % PLANT_N)
            return None
        return rel
    if how == "exempt_stale":
        edit(root / SCRIPT, EXEMPT_DECL,
             'NOT_STAMPED: dict[str, str] = {"%s": "acceptance: exempt but also in GLOBS"}' % TARGET)
        return TARGET
    if how == "manifest_stale":
        p = root / MANIFEST
        raw = p.read_bytes().decode("utf-8")
        new, n = re.subn(r"^" + re.escape(UPSTREAM_FILE) + r"\r?\n", "", raw, flags=re.M)
        assert n == 1
        p.write_bytes(new.encode("utf-8"))
        return UPSTREAM_FILE
    if how == "manifest_gone":
        (root / MANIFEST).unlink()
        return "upstream_src_manifest.txt"
    raise ValueError(how)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CASES = (
    ("unchanged", None, 0),
    ("stamp removed", "remove", 1),
    ("stamp stale", "stale", 1),
    ("Ladruno plant", "plant_l", 1),
    ("neutral plant", "plant_n", 1),
    ("manifest stale", "manifest_stale", 1),
    ("exempt ok", "exempt_ok", 0),
    ("exempt stale", "exempt_stale", 1),
    ("manifest gone", "manifest_gone", 1),
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD")
    ap.add_argument("--before", default="origin/ladruno")
    a = ap.parse_args()
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for label, how, want_rc in CASES:
            d = Path(tmp) / label.replace(" ", "_")
            d.mkdir()
            extract(a.ref, d)
            named_file = mutate(d, how) if how else None
            r = subprocess.run([sys.executable, str(d / SCRIPT), "--check"],
                               capture_output=True, text=True, encoding="utf-8", errors="replace")
            named = named_file is not None and named_file in r.stdout
            good = r.returncode == want_rc and (named if named_file else True)
            if how == "manifest_stale":
                good = good and "--refresh-upstream-manifest" in r.stdout
            ok &= good
            last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-120:]
            print(f"[{'PASS' if good else 'FAIL'}] {label:15s}: rc={r.returncode} (want {want_rc}), "
                  f"file named={named} | {last[:110]}")

        # history: both coverage rules on the pre-WP-122 tree, with that tree's own GLOBS
        d = Path(tmp) / "before"
        d.mkdir()
        extract(a.before, d)
        old = load(d / SCRIPT, "sh_before")
        new = load(REPO / SCRIPT, "sh_after")
        manifest, _ = new.read_upstream_manifest(REPO)
        authored = old.authored_files()
        flagged = {p.resolve() for p in new.ladruno_named_outside_globs(d, authored)}
        flagged |= {p.resolve() for p in new.fork_sources_outside_globs(d, authored, manifest, {})}
        rel = sorted(p.relative_to(d.resolve()).as_posix() for p in flagged)
        unstamped = [x for x in rel if STAMP not in (d / x).read_text(encoding="utf-8", errors="replace")]
        stamped = [x for x in rel if x not in unstamped]
        # every unstamped fork file in the before tree must be flagged
        missed = [p.relative_to(d).as_posix() for p in new.tracked_sources(d)
                  if p.relative_to(d).as_posix() not in manifest
                  and STAMP not in p.read_text(encoding="utf-8", errors="replace")
                  and p.relative_to(d).as_posix() not in unstamped]
        good = len(unstamped) == 31 and not missed
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] history ({a.before}): rules flag {len(rel)} files -- "
              f"{len(unstamped)} unstamped (want all 31 of WP-120 R1), {len(stamped)} stamped but "
              f"outside GLOBS; unstamped fork files missed: {len(missed)}")
        for x in missed:
            print("       MISSED:", x)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
