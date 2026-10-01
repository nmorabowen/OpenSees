#!/usr/bin/env bash
# Syntax-check the registration TUs that gained a NorSand line (errors only; these include hundreds of headers).
INC=$(find SRC -type d -printf '-I%p ')
STUB="-I Ladruno_files/testbed/norsand_oracle/kernel_parity/stub"
[ -f SRC/material/nD/LadrunoNorSandKernel.h ] && STUB=""
for f in "$@"; do
  echo "== $f"
  out=$(timeout 280 g++ -std=c++17 -fsyntax-only -w $INC $STUB "$f" 2>&1); st=$?
  echo "$out" | grep -E "error|LadrunoNorSand" | head -20
  echo "   g++ exit $st"
done
