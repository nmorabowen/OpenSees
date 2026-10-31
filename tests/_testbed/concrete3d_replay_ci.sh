#!/usr/bin/env bash
# concrete3d_replay_ci.sh [outdir] -- the Linux (glibc) lane of the LadrunoConcrete3D record replay (apeGmsh #2621 F1).
# Run from anywhere inside the repo; needs g++, git and ldd. Builds against THIS TREE's kernel header, always.
#
# The committed libmmt records (concrete3d_replay_records.txt) were written by the Windows pyd and replay bit-exactly
# only under Intel's libmmt (concrete3d_replay_libmmt.bat). The glibc fixture is concrete3d_replay_records.linux-glibc.txt,
# and the job runs on a pinned image (ubuntu-24.04 in ladruno.yml). The libm tag is linux-glibc-<version from
# `ldd --version`>: it is compiled in, stamped into every recording's #@ line, and compared by the libm guard, so a
# glibc update on the image reads as a LIBM MISMATCH (exit 2), never as a kernel change.
#   - glibc fixture committed: it is replayed bit-exactly on this tree's kernel, every record, every pass, every value.
#     THAT IS THE GATE. A kernel change that moves any value fails it by design; that PR re-records (this run's
#     artifact, or `record` on the same image), commits the new fixture, and the reviewers judge its diff.
#   - glibc fixture missing: the job FAILS, after recording. The flow: the first CI run of the PR that adds this lane is
#     red by design and uploads <outdir>/concrete3d_replay_records.linux-glibc.txt; that file is committed as
#     tests/_testbed/concrete3d_replay_records.linux-glibc.txt, and the next run replays it as the gate.
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

# e.g. "ldd (Ubuntu GLIBC 2.39-0ubuntu8.4) 2.39" -> 2.39-0ubuntu8.4 (the distro build, not just the upstream version)
ldd1="$(ldd --version | head -1)"
glibc="$(printf '%s\n' "$ldd1" | sed -n 's/.*GLIBC \([^) ]*\)).*/\1/p')"
[ -n "$glibc" ] || glibc="$(printf '%s\n' "$ldd1" | awk '{print $NF}')"
case "$glibc" in ''|*[!A-Za-z0-9.+~_-]*) echo "ERROR: cannot read the glibc version from: $ldd1"; exit 2;; esac
tag="linux-glibc-$glibc"
echo "libm tag $tag (from: $ldd1)"

g++ -std=c++17 -O2 -DC3D_REPLAY_LIBM="\"$tag\"" -I "$repo" "$here/concrete3d_replay.cpp" -o "$out/concrete3d_replay"
echo "built against this tree's kernel header, blob $blob"

echo "== libm guard: the libmmt records must be refused under glibc"
set +e; "$out/concrete3d_replay" check "$fix" > "$out/guard.txt" 2>&1; rc=$?; set -e
cat "$out/guard.txt"
[ "$rc" -eq 2 ] && grep -q "do not replay bit-exactly under $tag" "$out/guard.txt" \
  || { echo "ERROR: the libm guard did not refuse (exit $rc)"; exit 1; }

echo "== recording the inputs under glibc (artifact)"
"$out/concrete3d_replay" record "$fix" "$out/concrete3d_replay_records.linux-glibc.txt" --kernel-blob "$blob"

echo "== reported, not gated: the libmmt records replayed under glibc (cross-libm diagnostic)"
"$out/concrete3d_replay" check "$fix" --cross || true

if [ -f "$glibc_fix" ]; then
  echo "== GATE: the committed glibc fixture, bit-exact on this tree's kernel (refused if its glibc is not $glibc)"
  "$out/concrete3d_replay" check "$glibc_fix"
else
  echo "::error::No glibc fixture is committed (tests/_testbed/concrete3d_replay_records.linux-glibc.txt): FAILING BY DESIGN."
  echo "Commit this run's artifact concrete3d_replay_records.linux-glibc.txt ($tag) at that path; the next run gates on it."
  exit 1
fi
