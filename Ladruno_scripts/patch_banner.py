#!/usr/bin/env python3
"""
patch_banner.py — embed an ASCII-art banner from a .txt file into the
splash code in OpenSees/SRC/tcl/tclMain.cpp.

Reads `banner_ASCII.txt` (or any path passed as argv[1]), rebuilds the
C string literal block in `tclMain.cpp` between the
    // BANNER-START
    static const char *kBanner =
"...";
    // BANNER-END
markers, and writes back. Re-run after every edit to the .txt file.

Workflow:
    1. Edit `Ladruno_scripts/banner_ASCII.txt`
    2. python Ladruno_scripts/patch_banner.py
    3. Ladruno_scripts/build.bat OpenSees OpenSeesSP OpenSeesMP

The .txt is the source of truth; the C string in tclMain.cpp is a
mechanically generated copy. Do not hand-edit the C string between the
markers.
"""
from __future__ import annotations

import os
import re
import shutil
import sys
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT       = os.path.dirname(SCRIPT_DIR)
TARGET     = os.path.join(ROOT, "OpenSees", "SRC", "tcl", "tclMain.cpp")
BACKUP_DIR = os.path.join(ROOT, "Ladruno_files", "patch_banner_backups")

DEFAULT_BANNER_TXT = os.path.join(SCRIPT_DIR, "banner_ASCII.txt")

# Markers must match the comments in tclMain.cpp exactly (including indent).
START_MARKER = "    // BANNER-START\n"
END_MARKER   = "    // BANNER-END\n"


def escape_for_c(line: str) -> str:
    """Escape a single source line into a C string literal *body*.

    The banner txt is UTF-8; we keep multi-byte characters as raw bytes
    in the source file (which is also UTF-8). We only escape backslashes
    and double quotes — the two characters that would terminate or
    misparse the C literal.
    """
    line = line.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{line}\\n"'


def build_c_block(banner_lines: list[str]) -> str:
    body = "\n".join(escape_for_c(l) for l in banner_lines)
    return (
        f"{START_MARKER}"
        f"    static const char *kBanner =\n"
        f"{body};\n"
        f"{END_MARKER}"
    )


def main(argv: list[str]) -> int:
    banner_path = argv[1] if len(argv) > 1 else DEFAULT_BANNER_TXT

    if not os.path.exists(banner_path):
        print(f"ERROR: banner file not found: {banner_path}", file=sys.stderr)
        return 1

    if not os.path.exists(TARGET):
        print(f"ERROR: target source file not found: {TARGET}", file=sys.stderr)
        return 1

    # Read banner txt (UTF-8). Strip trailing whitespace per line; keep blank
    # lines (they're meaningful in art). Drop the trailing newline of the file.
    with open(banner_path, "r", encoding="utf-8") as f:
        banner_lines = f.read().splitlines()

    # Read existing tclMain.cpp
    with open(TARGET, "r", encoding="utf-8") as f:
        src = f.read()

    # Locate marker block by plain string search (avoids re.sub's
    # backslash-as-group-reference interpretation eating our literal `\n`
    # escapes inside the C string).
    start_idx = src.find(START_MARKER)
    end_idx   = src.find(END_MARKER, start_idx + len(START_MARKER))
    if start_idx < 0 or end_idx < 0:
        print(
            "ERROR: Could not find BANNER-START / BANNER-END markers in\n"
            f"       {TARGET}\n"
            "       Either the file was hand-edited or this is a fresh upstream\n"
            "       clone. Re-apply the banner-marker patch before running this.",
            file=sys.stderr,
        )
        return 1
    end_idx += len(END_MARKER)  # include the trailing END_MARKER line

    # Backup before overwrite
    os.makedirs(BACKUP_DIR, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = os.path.join(BACKUP_DIR, f"tclMain.cpp.{stamp}.bak")
    shutil.copy2(TARGET, backup_path)

    # Build replacement and splice
    new_block = build_c_block(banner_lines)
    new_src = src[:start_idx] + new_block + src[end_idx:]
    with open(TARGET, "w", encoding="utf-8", newline="\n") as f:
        f.write(new_src)

    # Report
    width = max((len(l) for l in banner_lines), default=0)
    print(f"banner txt: {banner_path}")
    print(f"  lines:      {len(banner_lines)}")
    print(f"  max width:  {width} chars (display width may differ for UTF-8)")
    print(f"target:     {TARGET}")
    print(f"backup:     {backup_path}")
    print(f"")
    print(f"Now rebuild:")
    print(f"  Ladruno_scripts\\build.bat OpenSees OpenSeesSP OpenSeesMP")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
