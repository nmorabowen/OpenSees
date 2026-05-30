#!/usr/bin/env python3
"""
patch_banner.py — embed the Ladruno splash banner and feature list into the
OpenSees interpreter startup code.

Two pieces of generated content, two source-of-truth text files:

  1. The ASCII-art banner       <- Ladruno_scripts/banner_ASCII.txt
       spliced into SRC/tcl/tclMain.cpp between
           // BANNER-START  ...  // BANNER-END

  2. The active-feature list    <- Ladruno_scripts/banner_features.txt
       spliced into BOTH
           SRC/tcl/tclMain.cpp            (Tcl  / `OpenSees`)
           SRC/interpreter/PythonModule.cpp (openseespy / openseesmp)
       between
           // FEATURES-START  ...  // FEATURES-END

The .txt files are the source of truth; the C string literals between the
markers are mechanically generated copies. Do NOT hand-edit between markers.

Workflow:
    1. Edit Ladruno_scripts/banner_ASCII.txt    (the art)   and/or
       Ladruno_scripts/banner_features.txt      (the feature list)
    2. python Ladruno_scripts/patch_banner.py
    3. Ladruno_scripts/build.bat OpenSees OpenSeesSP OpenSeesMP

banner_features.txt format: one feature per line. Blank lines and lines
starting with '#' are ignored. Each remaining line becomes a bullet under a
"Ladruno fork — active features:" header. Keep lines short (~70 display
cols) so a default console does not wrap them.
"""
from __future__ import annotations

import os
import re
import shutil
import sys
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Consolidated layout (2026-05-30): the repo root IS the parent of
# Ladruno_scripts; SRC sits directly under it. (The old path assumed an inner
# OpenSees/ directory that no longer exists — keep this in sync if the layout
# moves again.)
ROOT       = os.path.dirname(SCRIPT_DIR)
SRC        = os.path.join(ROOT, "SRC")
TCL_MAIN   = os.path.join(SRC, "tcl", "tclMain.cpp")
PY_MODULE  = os.path.join(SRC, "interpreter", "PythonModule.cpp")
BACKUP_DIR = os.path.join(SCRIPT_DIR, "patch_banner_backups")

BANNER_TXT   = os.path.join(SCRIPT_DIR, "banner_ASCII.txt")
FEATURES_TXT = os.path.join(SCRIPT_DIR, "banner_features.txt")

# Banner markers (banner art lives in tclMain.cpp only).
BANNER_START = "// BANNER-START"
BANNER_END   = "// BANNER-END"
# Feature markers (feature list lives in both tclMain.cpp and PythonModule.cpp).
FEAT_START   = "// FEATURES-START"
FEAT_END     = "// FEATURES-END"

FEATURES_HEADER = "      Ladruno fork — active features:"
BULLET_INDENT   = "        • "


def escape_for_c(line: str) -> str:
    """Escape one display line into a C string literal, incl. trailing \\n."""
    line = line.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{line}\\n"'


def find_marker_line(src: str, marker: str, from_idx: int = 0) -> tuple[int, int, str]:
    """Return (line_start, line_end_incl_newline, indent) for the line that
    contains `marker`. Matches by substring so the caller need not know the
    indentation. Raises ValueError if not found."""
    idx = src.find(marker, from_idx)
    if idx < 0:
        raise ValueError(marker)
    line_start = src.rfind("\n", 0, idx) + 1          # char after prev newline
    line_end   = src.find("\n", idx)
    line_end   = len(src) if line_end < 0 else line_end + 1
    indent     = src[line_start:idx]                  # whitespace before marker
    return line_start, line_end, indent


def splice_block(src: str, start_marker: str, end_marker: str,
                 var_name: str, display_lines: list[str]) -> str:
    """Replace the START..END marker block (inclusive) with a freshly built
    C string-literal declaration. The marker lines keep their original indent;
    the literal lines start at column 0, matching the existing banner style."""
    s_start, _, indent = find_marker_line(src, start_marker)
    _, e_end, _        = find_marker_line(src, end_marker, s_start)

    body = "\n".join(escape_for_c(l) for l in display_lines)
    block = (
        f"{indent}{start_marker}\n"
        f"{indent}static const char *{var_name} =\n"
        f"{body};\n"
        f"{indent}{end_marker}\n"
    )
    return src[:s_start] + block + src[e_end:]


def read_lines(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().splitlines()


def read_features(path: str) -> list[str]:
    """Read banner_features.txt -> display lines (header + bullets + blank)."""
    feats = [l.strip() for l in read_lines(path)]
    feats = [l for l in feats if l and not l.startswith("#")]
    return [FEATURES_HEADER] + [BULLET_INDENT + f for f in feats] + [""]


def backup(path: str, stamp: str) -> str:
    os.makedirs(BACKUP_DIR, exist_ok=True)
    dst = os.path.join(BACKUP_DIR, f"{os.path.basename(path)}.{stamp}.bak")
    shutil.copy2(path, dst)
    return dst


def write_src(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def main(argv: list[str]) -> int:
    for p in (BANNER_TXT, FEATURES_TXT, TCL_MAIN, PY_MODULE):
        if not os.path.exists(p):
            print(f"ERROR: required file not found: {p}", file=sys.stderr)
            return 1

    banner_lines   = read_lines(BANNER_TXT)
    feature_lines  = read_features(FEATURES_TXT)
    stamp          = datetime.now().strftime("%Y%m%d_%H%M%S")

    # --- tclMain.cpp: banner art + feature list ---
    tcl_src = open(TCL_MAIN, encoding="utf-8").read()
    try:
        tcl_src = splice_block(tcl_src, BANNER_START, BANNER_END, "kBanner", banner_lines)
        tcl_src = splice_block(tcl_src, FEAT_START, FEAT_END, "kFeatures", feature_lines)
    except ValueError as e:
        print(f"ERROR: marker {e} not found in {TCL_MAIN}", file=sys.stderr)
        return 1
    b = backup(TCL_MAIN, stamp)
    write_src(TCL_MAIN, tcl_src)
    print(f"patched: {TCL_MAIN}\n  backup: {b}")

    # --- PythonModule.cpp: feature list only ---
    py_src = open(PY_MODULE, encoding="utf-8").read()
    try:
        py_src = splice_block(py_src, FEAT_START, FEAT_END, "kFeatures", feature_lines)
    except ValueError as e:
        print(f"ERROR: marker {e} not found in {PY_MODULE}", file=sys.stderr)
        return 1
    b = backup(PY_MODULE, stamp)
    write_src(PY_MODULE, py_src)
    print(f"patched: {PY_MODULE}\n  backup: {b}")

    n_feat = len(feature_lines) - 2  # minus header and trailing blank
    print(f"\nbanner: {len(banner_lines)} lines   features: {n_feat}")
    print("Now rebuild:")
    print("  Ladruno_scripts\\build.bat OpenSees OpenSeesSP OpenSeesMP")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
