#!/usr/bin/env python3
"""classTag collision + drift gate (Axis 1, G2). Dependency-free.

Three checks across SRC/classTags.h and DEVELOPER/core/classTags.h:

  ERROR (exit 1): a value collision WITHIN a *_TAG family where at least one of
                  the colliding symbols is a Ladruno tag (comment contains
                  "Ladruno"). Protects new work without false-failing on the
                  inherited upstream mess.
  WARN:           any other within-family value collision (inherited upstream);
                  a symbol whose value disagrees between the two header copies
                  (the 205-line drift); a Ladruno tag below the private
                  band >= LADRUNO_BAND.

  --strict promotes all WARN to ERROR.

Family = the text before "_TAG_" (ELE_TAG_BezierTri6 -> family "ELE"). The
broker dispatches by classTag value within a category, so collisions only
matter inside a family.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HEADERS = [ROOT / "SRC" / "classTags.h", ROOT / "DEVELOPER" / "core" / "classTags.h"]
LADRUNO_BAND = 33000

DEF = re.compile(r"^\s*#define\s+(\w+)\s+(-?\d+)\b(.*)$")


def parse(path):
    """-> {symbol: (value, family, is_ladruno)}"""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        m = DEF.match(line)
        if not m:
            continue
        sym, val, tail = m.group(1), int(m.group(2)), m.group(3)
        if "_TAG_" not in sym:
            continue
        # family = everything up to the final _<Name> segment, so sibling
        # sub-namespaces stay distinct (DEG_TAG_STIFF vs DEG_TAG_UNLOAD vs
        # DEG_TAG_STRENGTH all share base "DEG_TAG_*" but dispatch separately).
        family = sym.rsplit("_", 1)[0]
        is_lad = "ladruno" in tail.lower()
        out[sym] = (val, family, is_lad)
    return out


def main(argv):
    strict = "--strict" in argv
    verbose = "--verbose" in argv  # also surface inherited-upstream collisions
    errors, warns = [], []

    parsed = {p: parse(p) for p in HEADERS}
    rel = {p: str(p.relative_to(ROOT)) for p in HEADERS}

    # within-file family collisions: ladruno-involved -> ERROR always; pure
    # inherited-upstream -> WARN only under --verbose (unactionable noise else).
    for p, table in parsed.items():
        byfam = {}
        for sym, (val, fam, lad) in table.items():
            byfam.setdefault((fam, val), []).append((sym, lad))
        for (fam, val), syms in byfam.items():
            if len(syms) < 2:
                continue
            names = ", ".join(s for s, _ in syms)
            involves_lad = any(lad for _, lad in syms)
            msg = f"{rel[p]}: family {fam}_* value {val} shared by: {names}"
            if involves_lad:
                errors.append(msg)
            elif verbose:
                warns.append(msg)

    # cross-file symbol disagreement (the 205-line drift)
    a, b = parsed[HEADERS[0]], parsed[HEADERS[1]]
    for sym in set(a) & set(b):
        if a[sym][0] != b[sym][0]:
            warns.append(
                f"drift: {sym} = {a[sym][0]} in {rel[HEADERS[0]]} but "
                f"{b[sym][0]} in {rel[HEADERS[1]]}"
            )

    # ladruno band policy
    for sym, (val, fam, lad) in parsed[HEADERS[0]].items():
        if lad and val < LADRUNO_BAND:
            warns.append(
                f"ladruno tag {sym} = {val} is below the private band "
                f">= {LADRUNO_BAND} (collision-prone against upstream growth)"
            )

    if strict:
        errors += warns
        warns = []

    for w in warns:
        print("WARN  " + w)
    for e in errors:
        print("ERROR " + e)
    if errors:
        print(f"check_classtags: FAIL ({len(errors)} errors, {len(warns)} warns)")
        return 1
    print(f"check_classtags: OK ({len(warns)} warns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
