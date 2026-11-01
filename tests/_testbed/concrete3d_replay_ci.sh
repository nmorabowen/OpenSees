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
ldd1="$(ldd --version | sed -n 1p)"
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

echo "== fixture guard: a corrupted fixture must be refused, never truncated, defaulted or ignored"
# name | awk edit of the libmmt fixture | exit | the message it must produce. Every edit is one an earlier loader
# accepted (REPLAY PASS 28/28): a value that truncates back, trailing junk, an unchecked index, key, line or domain, a
# field that was parsed and never compared. Exit 2 = refused at load; exit 1 = a compared value differs.
while IFS='|' read -r m prog want_rc want; do
  awk "$prog" "$fix" > "$out/$m.txt"
  cmp -s "$fix" "$out/$m.txt" && { echo "ERROR: mutant $m changed nothing"; exit 1; }
  set +e; "$out/concrete3d_replay" check "$out/$m.txt" --cross > "$out/$m.log" 2>&1; rc=$?; set -e
  echo "$m: exit $rc: $(grep -E 'REPLAY (ERROR|FAIL)|MISMATCH' "$out/$m.log" | head -1)"
  [ "$rc" -eq "$want_rc" ] && grep -q "$want" "$out/$m.log" || { echo "ERROR: mutant $m did not give exit $want_rc with '$want' (exit $rc)"; exit 1; }
done <<'MUTANTS'
O_status_frac|!d && /^O /{$2=$2+0.5; d=1} 1|2|O status is not an integer
Z_status_frac|!d && /^Z /{$6=$6+0.5; d=1} 1|2|Z status is not an integer
Z_pass_index|!d && /^Z /{$2=$2+1; d=1} 1|2|is not its position
Z_pass_range|!d && /^Z /{$2="4294967296"; d=1} 1|2|Z pass is out of int range
P_int_frac|!d && /^P /{sub(/ maxSubIncr=[0-9]+/, "&.5"); d=1} 1|2|P maxSubIncr is not an integer
P_num_junk|!d && /^P /{sub(/ E=[^ ]+/, "&zz"); d=1} 1|2|not a number
P_ctTemper_domain|/^P /{sub(/ ctTemper=0/, " ctTemper=3")} 1|2|P ctTemper=3 is not a value the wrapper writes
P_implex_domain|!d && /^P /{sub(/ implex=0/, " implex=2"); d=1} 1|2|P implex=2 is not a value the wrapper writes
P_maxSubIncr_domain|!d && /^P /{sub(/ maxSubIncr=10/, " maxSubIncr=9"); d=1} 1|2|P maxSubIncr=9 is not a value the wrapper writes
K_extra_key|!d && /^K /{$0=$0 " junk=1"; d=1} 1|2|unexpected key junk
K_dup_key|!d && /^K /{$0=$0 " seq=1"; d=1} 1|2|duplicate key seq
W_twice|!d && /^W /{print; d=1} 1|2|two W lines
G_drop|!d && /^G /{d=1; next} 1|2|lacks its G line
G_class|!d && /^G /{sub(/ class=[^ ]+/, " class=other"); d=1} 1|2|G class= is not a record class
C0_not_E2|!d && /^C /{$2="1e-3"; d=1} 1|2|is not E\[2\] bit for bit
meta_twice|!d && /^#@ /{print "#@ "; d=1} 1|2|the #@ line must come once
meta_key|/^#@ /{$0=$0 " junk=1"} 1|2|the #@ line has an unexpected key junk
meta_blob|/^#@ /{sub(/kernel-blob=[0-9a-f]+/, "kernel-blob=x")} 1|2|kernel-blob is not a 40-hex
T2_end_R|!d && /^T /{$4="0.5"; d=1} 1|1|T2end
MUTANTS

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
