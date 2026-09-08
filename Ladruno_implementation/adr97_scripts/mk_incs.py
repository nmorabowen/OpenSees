"""Extract the -I flags of the ASDP translation unit from build.ninja into a
GCC response file, so `g++ -fsyntax-only` can pre-flight the header edits with
the same include set MSVC uses (and catch the GCC-only two-phase-lookup /
temporary-binding errors before the 20-minute build)."""
import io
import os
import sys

NINJA = sys.argv[1]
OUT = sys.argv[2]
for line in io.open(NINJA, encoding="utf-8", errors="replace"):
    if line.startswith("  INCLUDES = "):
        toks = [t.replace("\\", "/")
                for t in line[len("  INCLUDES = "):].split() if t.startswith("-I")]
        io.open(OUT, "w", newline="\n").write("\n".join(toks) + "\n")
        print(len(toks), "include dirs ->", os.path.basename(OUT))
        break
else:
    raise SystemExit("no INCLUDES line in " + NINJA)
