#!/usr/bin/env bash
# concrete3d_replay_ci.sh [outdir] -- the Linux (glibc) lane of the LadrunoConcrete3D record replay (apeGmsh #2621 F1).
# Run from anywhere inside the repo; needs g++ and git. Builds against THIS TREE's kernel header, always.
#
# The committed libmmt records (concrete3d_replay_records.txt) were written by the Windows pyd and replay bit-exactly
# only under Intel's libmmt (concrete3d_replay_libmmt.bat). The glibc fixture is concrete3d_replay_records.linux-glibc.txt:
#   - committed: it is replayed bit-exactly on this tree's kernel, every record, every pass, every value. THAT IS THE
#     GATE. A kernel change that moves any value fails it by design; that PR re-records (the uploaded artifact below,
#     or `record` on a glibc box), commits the new fixture, and the reviewers judge its diff;
#   - not committed: the lane is RECORD-ONLY AND NOT A GATE. It records the libmmt fixture's inputs under glibc into
#     <outdir>/concrete3d_replay_records.linux-glibc.txt for upload, and prints so. The glibc fixture must be committed
#     from that artifact before the PR that adds this lane merges.
# Also checked: the libm guard refuses the libmmt records under glibc (a seam check, not a replay result). Reported,
# never gated: the cross-libm diff of the libmmt records under glibc (a libm difference by construction).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
out="${1:-${RUNNER_TEMP:-/tmp}/c3d_replay}"
mkdir -p "$out"
fix="$here/concrete3d_replay_records.txt"
glibc_fix="$here/concrete3d_replay_records.linux-glibc.txt"
blob="$(git -C "$repo" hash-object "$repo/SRC/material/nD/LadrunoConcrete3DKernel.h")"

g++ -std=c++17 -O2 -DC3D_REPLAY_LIBM='"linux-glibc"' -I "$repo" "$here/concrete3d_replay.cpp" -o "$out/concrete3d_replay"
echo "built against this tree's kernel header, blob $blob"

echo "== libm guard: the libmmt records must be refused under glibc"
set +e; "$out/concrete3d_replay" check "$fix" > "$out/guard.txt" 2>&1; rc=$?; set -e
cat "$out/guard.txt"
[ "$rc" -eq 2 ] && grep -q "do not replay bit-exactly under linux-glibc" "$out/guard.txt" \
  || { echo "ERROR: the libm guard did not refuse (exit $rc)"; exit 1; }

echo "== recording the inputs under glibc (artifact)"
"$out/concrete3d_replay" record "$fix" "$out/concrete3d_replay_records.linux-glibc.txt" --kernel-blob "$blob"

echo "== reported, not gated: the libmmt records replayed under glibc (cross-libm diagnostic)"
"$out/concrete3d_replay" check "$fix" --cross || true

if [ -f "$glibc_fix" ]; then
  echo "== GATE: the committed glibc fixture, bit-exact on this tree's kernel"
  "$out/concrete3d_replay" check "$glibc_fix"
else
  echo "::warning::RECORD-ONLY, NOT A GATE: no glibc fixture is committed (tests/_testbed/concrete3d_replay_records.linux-glibc.txt)."
  echo "Commit this run's artifact concrete3d_replay_records.linux-glibc.txt there before the PR that adds this lane merges."
fi
