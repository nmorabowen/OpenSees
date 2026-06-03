#!/usr/bin/env python3
"""Stamp the LADRUNO author/banner header onto every fork-authored source file.

Source of truth for the ASCII art is `Ladruno_scripts/banner_ASCII.txt` (the same
file that drives the runtime splash via patch_banner.py), so the art never drifts
between the console banner and the file headers.

The block is delimited by `// LADRUNO-HEADER-START` / `// LADRUNO-HEADER-END`
markers and is inserted *after* the leading OpenSees/PEER `/* ... */` comment so
upstream attribution stays first. Running again is idempotent: an existing block
is replaced in place (edit the art or credit here, re-run, done). Original line
endings (LF/CRLF) and any BOM are preserved so the stamp produces a clean diff.

    python Ladruno_scripts/stamp_headers.py            # stamp all authored files
    python Ladruno_scripts/stamp_headers.py --check     # report only, exit 1 if any stale
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# --- credit line (the four authors) ----------------------------------------
CREDIT = [
    "Ladruno — a research fork of OpenSees",
    "Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  "
    "José Abell  ·  Guppi",
]

START = "// LADRUNO-HEADER-START"
END = "// LADRUNO-HEADER-END"
RULE = "// " + "=" * 74

# --- the set of files WE authored (new fork code, not vanilla edits) --------
# Whole directories that are entirely ours, plus specific files in shared dirs.
GLOBS = [
    "SRC/element/ladrunoBrick/*.cpp", "SRC/element/ladrunoBrick/*.h",
    "SRC/element/ladrunoIMKBeam/*.cpp", "SRC/element/ladrunoIMKBeam/*.h",
    "SRC/element/ladrunoEmbeddedRebar/*.cpp", "SRC/element/ladrunoEmbeddedRebar/*.h",
    "SRC/element/bezierTriangle/*.cpp", "SRC/element/bezierTriangle/*.h",
    "SRC/element/bezierTetrahedron/*.cpp", "SRC/element/bezierTetrahedron/*.h",
    "SRC/element/solidTransformation/*.cpp", "SRC/element/solidTransformation/*.h",
    "SRC/utility/profiler/*.cpp", "SRC/utility/profiler/*.h",
    "SRC/material/nD/LadrunoJ2.*", "SRC/material/nD/LadrunoJ2Kernel.h",
    "SRC/material/nD/LadrunoJ2Finite.*",
    "SRC/material/nD/LadrunoHardening.h",
    "SRC/material/uniaxial/LadrunoUniaxialJ2.*",
    "SRC/material/uniaxial/LadrunoRebarBuckling.*",
    "SRC/material/uniaxial/LadrunoBondSlip.*",
    "SRC/material/nD/LogStrainNDMaterial.*", "SRC/material/nD/LogStrainKernel.h",
    "SRC/material/nD/FiniteStrainNDMaterial.h",
    "SRC/material/nD/InitDefGradNDMaterial.*",
    "SRC/material/nD/StagedStrainNDMaterial.*",
    "SRC/analysis/integrator/CentralDifferenceLadruno.*",
    "SRC/analysis/integrator/LadrunoArcLength.*",
    "SRC/analysis/integrator/LadrunoDynamicRelaxation.*",
    "SRC/analysis/integrator/LadrunoFictitiousMass.h",
    "SRC/analysis/integrator/ExplicitBathe.*",
    "SRC/analysis/integrator/ExplicitBatheLNVD.*",
    "SRC/recorder/LadrunoRecorder.*",
    "SRC/recorder/Ladruno_*.cpp", "SRC/recorder/Ladruno_*.h",
    "SRC/recorder/LadrunoMonitor*.cpp", "SRC/recorder/LadrunoMonitor*.h",
    "SRC/recorder/EnergyBalanceRecorder.*", "SRC/recorder/EnergyBalanceKernel.h",
]
SUFFIXES = {".cpp", ".h", ".hpp", ".cc", ".cxx"}


def authored_files() -> list[Path]:
    seen: dict[Path, None] = {}
    for g in GLOBS:
        for p in ROOT.glob(g):
            if p.is_file() and p.suffix in SUFFIXES:
                seen[p.resolve()] = None
    return sorted(seen)


def build_block(eol: str) -> str:
    art = (ROOT / "Ladruno_scripts" / "banner_ASCII.txt").read_text(
        encoding="utf-8").rstrip("\n").split("\n")
    lines = [START, RULE, "//"]
    for a in art:
        lines.append(("//  " + a).rstrip())
    lines.append("//")
    for c in CREDIT:
        lines.append("//  " + c)
    lines.append("//")
    lines.append("// Header auto-stamped by Ladruno_scripts/stamp_headers.py "
                 "(art: banner_ASCII.txt).")
    lines.append("// Do not hand-edit between the markers; edit the script/art "
                 "and re-run instead.")
    lines.append(RULE)
    lines.append(END)
    return eol.join(lines) + eol


_BLOCK_RE = re.compile(
    re.escape(START) + r".*?" + re.escape(END) + r"(?:\r?\n)?", re.S)


def restamp(text: str, block: str, eol: str) -> str:
    """Return text with the header block inserted or replaced (idempotent)."""
    if START in text:
        return _BLOCK_RE.sub(lambda _: block, text, count=1)
    # No existing block: insert after the leading /* ... */ comment if present.
    stripped = text.lstrip()
    if stripped.startswith("/*"):
        close = text.find("*/")
        nl = text.find("\n", close)
        if nl == -1:
            return text + eol + eol + block
        head, tail = text[:nl + 1], text[nl + 1:]
        return head + eol + block + tail
    return block + eol + text


def main() -> int:
    check = "--check" in sys.argv[1:]
    files = authored_files()
    changed: list[Path] = []
    for p in files:
        with open(p, "r", encoding="utf-8", newline="") as fh:
            raw = fh.read()
        bom = ""
        if raw.startswith("﻿"):
            bom, raw = "﻿", raw[1:]
        eol = "\r\n" if "\r\n" in raw else "\n"
        new = bom + restamp(raw, build_block(eol), eol)
        if new != bom + raw:
            changed.append(p)
            if not check:
                with open(p, "w", encoding="utf-8", newline="") as fh:
                    fh.write(new)

    rel = lambda p: p.relative_to(ROOT).as_posix()
    if check:
        if changed:
            print("STALE / unstamped ({}):".format(len(changed)))
            for p in changed:
                print("  " + rel(p))
            return 1
        print("All {} authored files carry a current header.".format(len(files)))
        return 0

    print("Stamped {} of {} authored files (already-current files unchanged):"
          .format(len(changed), len(files)))
    for p in changed:
        print("  " + rel(p))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
