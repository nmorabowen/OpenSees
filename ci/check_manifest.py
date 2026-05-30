#!/usr/bin/env python3
"""Manifest ledger gate (G9). Every Ladruno classTag must be accounted for.

For each Ladruno tag in SRC/classTags.h (comment contains "Ladruno"):
  ERROR (exit 1): no manifest row with that tag_symbol.
  ERROR:          row status == active but its pytest file is missing.
  WARN:           row status in {remediation, planned} with a missing test
                  (work pending — recorded, not yet done).

Requires PyYAML (CI installs it: `pip install pyyaml`).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLASSTAGS = ROOT / "SRC" / "classTags.h"
MANIFEST = ROOT / "Ladruno_implementation" / "testbed" / "manifest.yaml"

DEF = re.compile(r"^\s*#define\s+(\w+)\s+(-?\d+)\b(.*)$")


def ladruno_tags():
    tags = {}
    for line in CLASSTAGS.read_text(errors="replace").splitlines():
        m = DEF.match(line)
        if m and "ladruno" in m.group(3).lower():
            tags[m.group(1)] = int(m.group(2))
    return tags


def main(argv):
    try:
        import yaml
    except ModuleNotFoundError:
        print("check_manifest: SKIP (PyYAML not installed)")
        return 0

    data = yaml.safe_load(MANIFEST.read_text()) or {}
    rows = {r.get("tag_symbol"): r for r in data.get("features", []) if r.get("tag_symbol")}

    errors, warns = [], []
    for sym, val in ladruno_tags().items():
        row = rows.get(sym)
        if row is None:
            errors.append(f"classTag {sym} ({val}) has no manifest row")
            continue
        status = row.get("status", "active")
        pytest_path = row.get("pytest")
        tcl = str(row.get("tcl", ""))
        has_test = bool(pytest_path and (ROOT / pytest_path).exists())
        waived = tcl.startswith("WAIVED")
        if not has_test:
            msg = f"{sym}: manifest pytest '{pytest_path}' missing"
            if status in ("remediation", "planned"):
                warns.append(msg + f" (status={status}, pending)")
            else:
                errors.append(msg + f" (status={status})")
        if not (has_test or waived):
            warns.append(f"{sym}: neither a pytest file nor a WAIVED tcl reason")

    for w in warns:
        print("WARN  " + w)
    for e in errors:
        print("ERROR " + e)
    if errors:
        print(f"check_manifest: FAIL ({len(errors)} errors, {len(warns)} warns)")
        return 1
    print(f"check_manifest: OK ({len(warns)} warns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
