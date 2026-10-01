#!/usr/bin/env bash
# WP-144 P1b: g++ -fsyntax-only on the NorSand shell TUs. Run from the repo root (Esmeralda ~/wp144/repo).
# Uses the real kernel header if SRC/material/nD/LadrunoNorSandKernel.h exists, else the throwaway stub.
# Prints ONLY diagnostics located in LadrunoNorSand* / LadrunoNorSand*Kernel files (upstream-header -Wextra noise is dropped).
INC=$(find SRC -type d -printf '-I%p ')
STUB="-I Ladruno_files/testbed/norsand_oracle/kernel_parity/stub"
if [ -f SRC/material/nD/LadrunoNorSandKernel.h ]; then echo "(real kernel header)"; STUB=""; else echo "(STUB kernel header)"; fi
rc=0
for f in "$@"; do
  echo "== $f"
  out=$(g++ -std=c++17 -fsyntax-only -Wall -Wextra $INC $STUB "$f" 2>&1); st=$?
  echo "$out" | grep -E -A4 "^[^ ]*LadrunoNorSand[^ ]*:[0-9]+:[0-9]+: (warning|error)" | head -80
  echo "   g++ exit $st, errors: $(echo "$out" | grep -c ' error:')"
  [ $st -ne 0 ] && rc=1
done
exit $rc
