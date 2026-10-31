#!/usr/bin/env bash
# concrete3d_replay_ci.sh [outdir] -- the Linux (glibc) lane of the LadrunoConcrete3D record replay (apeGmsh #2621 F1).
# Run from anywhere inside the repo; needs g++ and git (the full history only when the kernel header has moved).
#
# The committed records (concrete3d_replay_records.txt) were written by the Windows pyd, whose libm is Intel's libmmt;
# they replay bit-exactly only there (concrete3d_replay_libmmt.bat). No Linux binary has written records, so on glibc
# this is RECORD-AND-COMPARE: the testbed is built against the kernel the records were written by (the fixture's
# kernel-blob, extracted from git when the tree's header differs), re-records the same inputs under glibc, and a fresh
# process must replay that recording bit-exactly. Gated:
#   1. the fixture parses fail-closed and the libm guard refuses the libmmt records under glibc (exit 2),
#   2. record under glibc, then check the recording: every record, every pass, bit-exact.
# Reported, not gated: the cross-libm diff of the committed libmmt records under glibc (a libm difference by
# construction, never a kernel verdict). The glibc recording is left in <outdir> for upload.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
out="${1:-${RUNNER_TEMP:-/tmp}/c3d_replay}"
mkdir -p "$out"
fix="$here/concrete3d_replay_records.txt"
hdr="SRC/material/nD/LadrunoConcrete3DKernel.h"

want="$(grep -m1 '^#@ ' "$fix" | tr ' ' '\n' | sed -n 's/^kernel-blob=//p')"
[ -n "$want" ] || { echo "ERROR: the fixture names no kernel-blob"; exit 2; }
have="$(git -C "$repo" hash-object "$repo/$hdr")"
inc="$repo"
if [ "$have" != "$want" ]; then
  echo "kernel header $have differs from the recording kernel $want: building the replay against the recorded blob"
  mkdir -p "$out/ref/$(dirname "$hdr")"
  if [ "$(git -C "$repo" rev-parse --is-shallow-repository)" = true ]; then git -C "$repo" fetch --quiet --unshallow; fi
  git -C "$repo" cat-file blob "$want" > "$out/ref/$hdr" || { echo "ERROR: blob $want is not in this clone (fetch the full history)"; exit 2; }
  inc="$out/ref"
fi

g++ -std=c++17 -O2 -DC3D_REPLAY_LIBM='"linux-glibc"' -I "$inc" "$here/concrete3d_replay.cpp" -o "$out/concrete3d_replay"

echo "== 1. libm guard: the libmmt records must be refused under glibc"
set +e; "$out/concrete3d_replay" check "$fix" > "$out/guard.txt" 2>&1; rc=$?; set -e
cat "$out/guard.txt"
[ "$rc" -eq 2 ] && grep -q "do not replay bit-exactly under linux-glibc" "$out/guard.txt" \
  || { echo "ERROR: the libm guard did not refuse (exit $rc)"; exit 1; }

echo "== 2. record-and-compare under glibc"
"$out/concrete3d_replay" record "$fix" "$out/concrete3d_replay_records.linux-glibc.txt"
"$out/concrete3d_replay" check "$out/concrete3d_replay_records.linux-glibc.txt"

echo "== reported, not gated: the committed libmmt records replayed under glibc (cross-libm diagnostic)"
"$out/concrete3d_replay" check "$fix" --cross || true
