#!/usr/bin/env bash
# WP-130 after review round 1: the RECOMMENDED recipe on F12's bearing deck,
# plus an IntScheme-1 control run back to back on the same (loaded) box, so the
# wall-clock comparison is like for like. SEQUENTIAL. Needs LADRUNO_DIST_BIN, PY.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
W=${WALL:-1200}
A="$HERE/run_arm.sh"
echo "=== rec_s2: IntScheme 2 (fixed default) -cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on ==="
bash $A rec_s2 2 $W -- -cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on
echo "=== rec_s1_T2: IntScheme 1, TanType 2 (same-load control) ==="
bash $A rec_s1_T2 1 $W
echo "=== DONE ==="
