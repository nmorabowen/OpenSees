#!/usr/bin/env bash
# WP-144: syntax-check the registration TUs that gained a NorSand line (errors only; these include hundreds of
# headers). Run from the repo root (Esmeralda ~/wp144/repo):
#   bash Ladruno_files/testbed/norsand_oracle/kernel_parity/tools/syntax_check_reg.sh <TU> ...
INC=$(find SRC -type d -printf '-I%p ')
for f in "$@"; do
  echo "== $f"
  out=$(timeout 280 g++ -std=c++17 -fsyntax-only -w $INC "$f" 2>&1); st=$?
  echo "$out" | grep -E "error|LadrunoNorSand" | head -20
  echo "   g++ exit $st"
done
