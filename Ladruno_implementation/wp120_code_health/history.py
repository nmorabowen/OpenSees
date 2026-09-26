#!/usr/bin/env python3
"""WP-120 step 2: cross-reference clone families with bug-fix history.

For each family (a set of files from clones.py --json), count
  * PRs on `ladruno` (first-parent history; a merge commit's diff vs its first
    parent = the PR's whole change) that touch >= 2 member files, and among them
    the ones whose title/branch reads as a FIX (regex below);
  * LEDGER_quirks.md sections (split on markdown headings) that name >= 2 member
    classes.
A family with replicated-fix history is a refactor candidate; one without is not (yet).

    python history.py exact.json [norm.json ...] [--ref origin/ladruno]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT  # noqa: E402

FIX = re.compile(r"\b(fix\w*|bug\w*|wrong|correct\w*|repair|sign|double|twice|clobber\w*|stale|"
                 r"missing|guard|hole|leak\w*|crash\w*|segfault|nan|refresh|snapshot|frozen|"
                 r"regression|restore\w*|broken|never|false)\b", re.I)


def _git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


def pr_units(ref):
    """One unit per first-parent commit of `ref` (a PR merge or a squash)."""
    meta = {}
    for rec in _git("log", "--first-parent", ref, "--format=%H%x1f%s%x1f%b%x1e").split("\x1e"):
        parts = rec.strip("\n").split("\x1f")
        if len(parts) >= 2:
            body = parts[2].strip().splitlines() if len(parts) > 2 else []
            meta[parts[0]] = (parts[1], body[0] if body else "")
    units, cur = [], None
    for ln in _git("log", "--first-parent", ref, "-m", "--name-only", "--format=@@%H").splitlines():
        if ln.startswith("@@"):
            h = ln[2:].strip()
            subj, title = meta.get(h, ("", ""))
            cur = {"sha": h[:9], "subject": subj, "title": title, "files": set()}
            units.append(cur)
        elif ln.strip() and cur is not None:
            cur["files"].add(ln.strip())
    return units


def quirk_sections():
    text = (ROOT / "Ladruno_implementation" / "LEDGER_quirks.md").read_text(encoding="utf-8", errors="replace")
    secs, head, buf = [], "(preamble)", []
    for ln in text.splitlines():
        if re.match(r"#{2,4}\s", ln):
            secs.append((head, "\n".join(buf)))
            head, buf = ln.strip("# ").strip(), []
        else:
            buf.append(ln)
    secs.append((head, "\n".join(buf)))
    return secs


def cls(path):
    stem = Path(path).stem
    return stem[4:] if stem.startswith("OPS_") else stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--ref", default="origin/ladruno")
    ap.add_argument("--out")
    a = ap.parse_args()

    fams, seen = [], set()
    for j in a.jsons:
        for fm in json.loads(Path(j).read_text(encoding="utf-8"))["families"]:
            key = frozenset(fm["files"])
            if key not in seen:
                seen.add(key)
                fams.append(fm)
    units = pr_units(a.ref)
    secs = quirk_sections()
    rows = []
    for fm in fams:
        members = set(fm["files"])
        classes = {cls(f) for f in members}
        both, fixes = [], []
        for u in units:
            hit = members & u["files"]
            if len(hit) >= 2:
                both.append(u)
                label = u["subject"] + " " + u["title"]
                if FIX.search(label):
                    fixes.append({"sha": u["sha"], "label": (u["title"] or u["subject"])[:110],
                                  "files": sorted(Path(x).name for x in hit)})
        qhits = []
        for head, body in secs:
            named = sorted(c for c in classes if re.search(r"\b" + re.escape(c) + r"\b", head + "\n" + body))
            if len(named) >= 2:
                qhits.append({"heading": head[:110], "classes": named})
        rows.append({"files": sorted(members), "copies": len(members),
                     "shared_lines_sum_pairs": fm["shared_lines_sum_pairs"],
                     "prs_touching_2plus": len(both), "fix_prs": fixes, "quirk_sections": qhits})
    for r in rows:
        print(f"== {r['copies']} files ({', '.join(Path(f).name for f in r['files'][:6])}"
              f"{' ...' if r['copies'] > 6 else ''}) shared={r['shared_lines_sum_pairs']}  "
              f"PRs>=2={r['prs_touching_2plus']} fix-PRs={len(r['fix_prs'])} quirk-secs={len(r['quirk_sections'])}")
        for f in r["fix_prs"][:12]:
            print(f"     PR {f['sha']} {f['label']}  [{', '.join(f['files'][:5])}]")
        for q in r["quirk_sections"][:10]:
            print(f"     Q  {q['heading']}  [{', '.join(q['classes'][:6])}]")
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
