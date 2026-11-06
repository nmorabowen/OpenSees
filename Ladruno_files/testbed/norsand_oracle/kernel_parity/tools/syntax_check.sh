#!/usr/bin/env bash
# WP-144: g++ -fsyntax-only -Wall -Wextra -Wshadow on the NorSand shell TUs, against the REAL kernel header.
# Run from the repo root (Esmeralda ~/wp144/repo):
#   bash Ladruno_files/testbed/norsand_oracle/kernel_parity/tools/syntax_check.sh SRC/material/nD/LadrunoNorSand*.cpp
# Prints ONLY diagnostics located in LadrunoNorSand* files (upstream-header -Wextra noise is dropped).
INC=$(find SRC -type d -printf '-I%p ')
rc=0
for f in "$@"; do
  echo "== $f"
  out=$(g++ -std=c++17 -fsyntax-only -Wall -Wextra -Wshadow $INC "$f" 2>&1); st=$?
  echo "$out" | grep -E -A4 "^[^ ]*LadrunoNorSand[^ ]*:[0-9]+:[0-9]+: (warning|error)" | head -80
  echo "   g++ exit $st, errors: $(echo "$out" | grep -c ' error:')"
  [ $st -ne 0 ] && rc=1
done
exit $rc
