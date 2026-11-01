#!/usr/bin/env bash
# concrete3d_replay_ci.sh [outdir]            -- the Linux (glibc) lane of the LadrunoConcrete3D record replay
# concrete3d_replay_ci.sh --contract [outdir] -- its field-classification contract (a workflow step of its own)
# (apeGmsh #2621 F1). Run from anywhere inside the repo; needs g++, git, awk and ldd. Builds against THIS TREE's kernel
# header, always.
#
# The committed libmmt records (concrete3d_replay_records.txt) were written by the Windows pyd and replay bit-exactly
# only under Intel's libmmt (concrete3d_replay_libmmt.bat). The glibc fixture is concrete3d_replay_records.linux-glibc.txt,
# and the job runs on a pinned image (ubuntu-24.04 in ladruno.yml). The libm tag is linux-glibc-<version from
# `ldd --version`>: it is compiled in, stamped into every recording's #@ line, and compared by the libm guard, so a
# glibc update on the image reads as a LIBM MISMATCH (exit 2), never as a kernel change. `check` also compares the
# fixture's #@ kernel-blob with this tree's kernel header blob, so any kernel edit needs a re-record.
#   - glibc fixture committed: it is replayed bit-exactly on this tree's kernel, every record, every pass, every value,
#     the full returned State (N lines) included. THAT IS THE GATE. A kernel change that moves any value fails it by
#     design; that PR re-records (this run's artifact, or `record` on the same image), commits the new fixture, and
#     the reviewers judge its diff. After the gate: the glibc fixture must hold the libmmt fixture's records, order and
#     inputs (`inputs-match`); the libmmt fixture must be recorded on this tree's kernel blob (the blob tie, proven by
#     its own mutant); corrupted outputs must be caught; and two mutated kernels must FAIL the replay on their N lines.
#   - --contract (the next workflow step, after the gate): the testbed's field table, the single source of the FIELD
#     CLASSIFICATION, must match the header block it generates (`fields --check`), and every row `contract` generates
#     from it must give its class's outcome on the glibc fixture: each IN row FAILS (or is refused where its entry says
#     so), each OUT row fails naming its compare, each PIN row is refused, and each UNREAD and INFO row still PASSES
#     (the coverage limit, enforced). Each failure line names the field and its class.
#   - glibc fixture missing: the job FAILS, after recording. The flow: the first CI run of the PR that adds this lane is
#     red by design and uploads <outdir>/concrete3d_replay_records.linux-glibc.txt; that file is committed as
#     tests/_testbed/concrete3d_replay_records.linux-glibc.txt, and the next run replays it as the gate.
# Also checked: the libm guard refuses the libmmt records under glibc (a seam check, not a replay result), and corrupted
# libmmt fixtures are refused at load. Reported, never gated: the cross-libm diff of the libmmt records under glibc.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
contract=0
if [ "${1:-}" = --contract ]; then contract=1; shift; fi
out="${1:-${RUNNER_TEMP:-/tmp}/c3d_replay}"
mkdir -p "$out"
fix="$here/concrete3d_replay_records.txt"
glibc_fix="$here/concrete3d_replay_records.linux-glibc.txt"
kernel="SRC/material/nD/LadrunoConcrete3DKernel.h"
blob="$(git -C "$repo" hash-object "$repo/$kernel")"

# e.g. "ldd (Ubuntu GLIBC 2.39-0ubuntu8.4) 2.39" -> 2.39-0ubuntu8.4 (the distro build, not just the upstream version)
ldd1="$(ldd --version | sed -n 1p)"
glibc="$(printf '%s\n' "$ldd1" | sed -n 's/.*GLIBC \([^) ]*\)).*/\1/p')"
[ -n "$glibc" ] || glibc="$(printf '%s\n' "$ldd1" | awk '{print $NF}')"
case "$glibc" in ''|*[!A-Za-z0-9.+~_-]*) echo "ERROR: cannot read the glibc version from: $ldd1"; exit 2;; esac
tag="linux-glibc-$glibc"
echo "libm tag $tag (from: $ldd1)"

# build <exe> <include root>: the testbed against <root>/SRC/material/nD/LadrunoConcrete3DKernel.h
build() { g++ -std=c++17 -O2 -DC3D_REPLAY_LIBM="\"$tag\"" -I "$2" "$here/concrete3d_replay.cpp" -o "$1"; }
exe="$out/concrete3d_replay"

if [ "$contract" = 1 ]; then
  # the replay step's testbed when it is newer than both of its sources (one build per job), else a fresh build
  if [ -x "$exe" ] && [ "$exe" -nt "$here/concrete3d_replay.cpp" ] && [ "$exe" -nt "$repo/$kernel" ]; then
    echo "reusing $exe (built by the replay step from this tree's kernel header, blob $blob)"
  else
    build "$exe" "$repo"; echo "built against this tree's kernel header, blob $blob"
  fi
  [ -f "$glibc_fix" ] || { echo "::error::No glibc fixture is committed: the contract runs on it (see the replay step)."; exit 1; }
  echo "== the header's FIELD TABLE block must be the compiled field table (its single source)"
  "$exe" fields --check "$here/concrete3d_replay.cpp"
  echo "== contract: every row generated from the field table, on the glibc fixture (this libm's own records)"
  set +e; "$exe" contract "$glibc_fix" --kernel-blob "$blob" > "$out/contract.txt" 2>&1; rc=$?; set -e
  cat "$out/contract.txt"
  grep '^CONTRACT FAIL: ' "$out/contract.txt" | sed 's/^CONTRACT FAIL: /::error::field-classification contract: /' || true
  [ "$rc" -eq 0 ] || { echo "ERROR: the field-classification contract failed (exit $rc): each CONTRACT FAIL line names the field and its class"; exit 1; }
  exit 0
fi

build "$exe" "$repo"
echo "built against this tree's kernel header, blob $blob"

echo "== libm guard: the libmmt records must be refused under glibc"
set +e; "$exe" check "$fix" --kernel-blob "$blob" > "$out/guard.txt" 2>&1; rc=$?; set -e
cat "$out/guard.txt"
[ "$rc" -eq 2 ] && grep -q "do not replay bit-exactly under $tag" "$out/guard.txt" \
  || { echo "ERROR: the libm guard did not refuse (exit $rc)"; exit 1; }

echo "== recording the inputs under glibc (artifact)"
"$exe" record "$fix" "$out/concrete3d_replay_records.linux-glibc.txt" --kernel-blob "$blob"

echo "== reported, not gated: the libmmt records replayed under glibc (cross-libm diagnostic)"
"$exe" check "$fix" --cross --kernel-blob "$blob" || true

# mutants <fixture> <command...> : each row of stdin is name | awk edit (applied to <fixture>) | exit | message (a grep
# pattern the run's output must contain). The command runs with the word MUT replaced by the mutated fixture. A
# mutant that changes nothing is an error. Exit 2 = refused at load; exit 1 = a compared value differs; exit 0 = the
# replay passes (the UNREAD-input rows, which must NOT be read).
mutants() {
  local src="$1" m prog want_rc want rc a; shift
  local -a cmd
  while IFS='|' read -r m prog want_rc want; do
    awk "$prog" "$src" > "$out/$m.txt"
    cmp -s "$src" "$out/$m.txt" && { echo "ERROR: mutant $m changed nothing"; exit 1; }
    cmd=(); for a in "$@"; do if [ "$a" = MUT ]; then cmd+=("$out/$m.txt"); else cmd+=("$a"); fi; done
    set +e; "${cmd[@]}" > "$out/$m.log" 2>&1; rc=$?; set -e
    echo "$m: exit $rc: $(grep -E 'REPLAY (ERROR|FAIL|PASS)|MISMATCH|INPUTS (DIFFER|MATCH)' "$out/$m.log" | head -1)"
    [ "$rc" -eq "$want_rc" ] && grep -q "$want" "$out/$m.log" || { echo "ERROR: mutant $m did not give exit $want_rc with '$want' (exit $rc)"; exit 1; }
  done
}

# Every row is an edit an earlier loader accepted (REPLAY PASS): a value that truncates back, trailing junk, a
# non-canonical spelling, an unchecked index, key, line, parameter domain or dump rule, a field that was parsed and
# never compared, a # line that hides a record.
echo "== load guard: corrupted libmmt fixtures must be refused at load (libm-independent: load runs before any replay)"
mutants "$fix" "$exe" check MUT --cross --kernel-blob "$blob" <<'MUTANTS'
O_status_frac|!d && /^O /{$2=$2+0.5; d=1} 1|2|O\[0\] is not an integer
O_status_negzero|!d && /^O 0 /{$2="-0"; d=1} 1|2|O\[0\] is not a canonical integer
Z_status_frac|!d && /^Z /{$6=$6+0.5; d=1} 1|2|Z\[4\] is not an integer
Z_pass_index|!d && /^Z /{$2=$2+1; d=1} 1|2|is not its position
Z_pass_range|!d && /^Z /{$2="4294967296"; d=1} 1|2|Z\[0\] is out of int range
P_int_frac|!d && /^P /{sub(/ maxSubIncr=[0-9]+/, "&.5"); d=1} 1|2|P maxSubIncr is not an integer
P_int_noncanon|!d && /^P /{sub(/ ctTemper=0 /, " ctTemper=00 "); d=1} 1|2|P ctTemper is not a canonical integer
P_num_junk|!d && /^P /{sub(/ E=[^ ]+/, "&zz"); d=1} 1|2|not a number
P_num_noncanon|!d && /^P /{sub(/ rho=0 /, " rho=0.0 "); d=1} 1|2|not canonical
E_num_noncanon|!d && /^E / && $6=="0"{$6="0.0"; d=1} 1|2|not canonical
P_ctTemper_domain|/^P /{sub(/ ctTemper=0/, " ctTemper=3")} 1|2|P ctTemper=3 is not a value the wrapper writes
P_implex_domain|!d && /^P /{sub(/ implex=0/, " implex=2"); d=1} 1|2|P implex=2 is not a value the wrapper writes
P_maxSubIncr_domain|!d && /^P /{sub(/ maxSubIncr=10/, " maxSubIncr=9"); d=1} 1|2|P maxSubIncr=9 is not a value the wrapper writes
P_implexRmax_rule|!d && /^P /{sub(/ implexRmax=2 /, " implexRmax=2.5 "); d=1} 1|2|is not the dump's constant 2
P_lch_ref_rule|!d && /^P /{sub(/ lch_ref=[^ ]+/, " lch_ref=1"); d=1} 1|2|P lch_ref is not lch bit for bit
K_extra_key|!d && /^K /{$0=$0 " junk=1"; d=1} 1|2|unexpected key junk
K_dup_key|!d && /^K /{$0=$0 " seq=1"; d=1} 1|2|duplicate key seq
K_int_noncanon|!d && /^K /{sub(/ tag=1 /, " tag=+1 "); d=1} 1|2|K tag is not a canonical integer
K_counter_noncanon|!d && /^K /{sub(/ seq=/, " seq=0"); d=1} 1|2|K seq is not a canonical whole number
K_cls_formula|!d && /^K kind=C/ && / subInfo=0 /{sub(/ cls=0 /, " cls=1 "); d=1} 1|2|is not what the dump writes
W_twice|!d && /^W /{print; d=1} 1|2|two W lines
G_drop|!d && /^G /{d=1; next} 1|2|lacks its G line
G_class|!d && /^G /{sub(/ class=[^ ]+/, " class=other"); d=1} 1|2|G class= is not a record class
N_drop|!d && /^N /{d=1; next} 1|2|lacks its N line
C0_not_E2|!d && /^C /{$2="0.001"; d=1} 1|2|is not E\[2\] bit for bit
C1_pass0_pin|!d && /^K /{f=($0 ~ /^K kind=[CF] / && $0 ~ / pass=0 /)} !d && f && /^C /{$3="0.001"; d=1} 1|2|loop start C\[1\] is not its trial T\[2\]
comment_record|/^K /{n++} n==2{$0="# " $0} 1|2|a # line after the first record
meta_twice|!d && /^#@ /{print "#@ "; d=1} 1|2|the #@ line must come once
meta_key|/^#@ /{$0=$0 " junk=1"} 1|2|the #@ line has an unexpected key junk
meta_blob|/^#@ /{sub(/kernel-blob=[0-9a-f]+/, "kernel-blob=x")} 1|2|kernel-blob is not a 40-hex
MUTANTS

if [ ! -f "$glibc_fix" ]; then
  echo "::error::No glibc fixture is committed (tests/_testbed/concrete3d_replay_records.linux-glibc.txt): FAILING BY DESIGN."
  echo "Commit this run's artifact concrete3d_replay_records.linux-glibc.txt ($tag) at that path; the next run gates on it."
  exit 1
fi

echo "== GATE: the committed glibc fixture, bit-exact on this tree's kernel (refused if its glibc is not $glibc)"
"$exe" check "$glibc_fix" --kernel-blob "$blob"

echo "== the glibc fixture holds the libmmt fixture's records, in its order, with its inputs"
"$exe" inputs-match "$glibc_fix" "$fix"
# inputs-match compares inputs only: the libmmt fixture (the pyd's values) must also be re-recorded on THIS kernel.
# Its values are proven only by the self-hosted libmmt replay; this ties its blob to the tree on every PR run.
# blob_tie <fixture>: returns 1, with the ERROR line, unless <fixture> is recorded on this tree's kernel blob
blob_tie() {
  local fb; fb="$(sed -n 's/^#@ .*kernel-blob=\([0-9a-f]*\).*/\1/p' "$1")"
  [ "$fb" = "$blob" ] || { echo "ERROR: the libmmt fixture is recorded on kernel-blob=$fb, this tree's kernel is $blob: re-record it under libmmt (concrete3d_replay_libmmt.bat record) in this PR"; return 1; }
  echo "libmmt fixture kernel-blob matches this tree ($blob)"
}
blob_tie "$fix" || exit 1
# the tie must fire: the same function on a copy of the libmmt fixture recorded on another kernel blob
other="$(printf '%040d' 0)"; [ "$other" != "$blob" ] || other="$(printf '1%039d' 0)"
awk -v b="$other" '/^#@ /{sub(/kernel-blob=[0-9a-f]+/, "kernel-blob=" b)} 1' "$fix" > "$out/blob_tie_mut.txt"
cmp -s "$fix" "$out/blob_tie_mut.txt" && { echo "ERROR: the blob-tie mutant changed nothing"; exit 1; }
set +e; blob_tie "$out/blob_tie_mut.txt" > "$out/blob_tie_mut.log" 2>&1; rc=$?; set -e
echo "blob-tie mutant: exit $rc: $(head -1 "$out/blob_tie_mut.log")"
[ "$rc" -eq 1 ] && grep -q "^ERROR: the libmmt fixture is recorded on kernel-blob=$other, this tree's kernel is $blob" "$out/blob_tie_mut.log" \
  || { echo "ERROR: the blob tie did not fire on a fixture recorded on another kernel blob (exit $rc)"; exit 1; }
mutants "$glibc_fix" "$exe" inputs-match MUT "$fix" <<'MUTANTS'
rec_drop|/^K /{n++} n!=2|1|INPUTS DIFFER: 27 records vs 28
S_edit|!d && /^S /{$2=sprintf("%.17g", $2*1.25); d=1} 1|1|INPUTS DIFFER at record 0 (seq=4523): S line
K_ele_edit|!d && /^K /{sub(/ ele=[0-9]+/, "&0"); d=1} 1|1|INPUTS DIFFER at record 0 (seq=4523): K ele=
T2_R_exempt|!d && /^T /{$4="0.5"; d=1} 1|0|INPUTS MATCH: 28 records
MUTANTS

echo "== compare guard: an output corrupted in the glibc fixture must be caught (on this libm's own records, after the gate)"
mutants "$glibc_fix" "$exe" check MUT --kernel-blob "$blob" <<'MUTANTS'
T2_end_R|!d && /^T /{$4="0.5"; d=1} 1|1|1 values differ: T2end
N_state_value|!d && /^N /{$28=sprintf("%.17g", $28 * 1.25 + 1); d=1} 1|1|1 values differ: N.sigcMax
MUTANTS

# The coverage limit (the UNREAD rows) is enforced by --contract, the next workflow step: its rows are generated from
# the field table in concrete3d_replay.cpp, never written here.

# Mutated kernels: a copy of this tree's header with an exact, pinned edit (each anchor must occur exactly once, or
# the step fails loudly), built as a second testbed and replayed on the glibc fixture. Each moves only returned State
# fields that O, D, W and Z never show, so only the N line catches it (the 968a01f8d testbed, which had no N line,
# replays both 28/28). kmut <name> <field the failure must name first> <anchor> <replacement> [<anchor> <replacement>...]
kmut() {
  local name="$1" field="$2" root="$out/kmut_$1" rc kb; shift 2
  mkdir -p "$root/$(dirname "$kernel")"
  cp "$repo/$kernel" "$root/$kernel"
  while [ "$#" -ge 2 ]; do
    awk -v a="$1" -v b="$2" '{ i = index($0, a); if (i) { $0 = substr($0, 1, i - 1) b substr($0, i + length(a)); n++ } print }
                             END { if (n != 1) { print "anchor found " n + 0 " times, not once: " a > "/dev/stderr"; exit 3 } }' \
      "$root/$kernel" > "$root/k.tmp" || { echo "ERROR: kernel mutant $name: its anchor is not in the header exactly once: $1"; exit 1; }
    mv "$root/k.tmp" "$root/$kernel"; shift 2
  done
  cmp -s "$repo/$kernel" "$root/$kernel" && { echo "ERROR: kernel mutant $name changed nothing"; exit 1; }
  build "$root/concrete3d_replay" "$root"
  kb="$(git hash-object "$root/$kernel")"
  set +e; "$root/concrete3d_replay" check "$glibc_fix" --kernel-blob "$kb" > "$root/check.log" 2>&1; rc=$?; set -e
  echo "kernel mutant $name: exit $rc: $(grep -E '^  MISMATCH' "$root/check.log" | head -1)"
  [ "$rc" -eq 1 ] && grep -q "MISMATCH .*values differ: N\.$field(" "$root/check.log" \
    || { echo "ERROR: kernel mutant $name did not fail the replay on N.$field (exit $rc)"; tail -5 "$root/check.log"; exit 1; }
}
echo "== kernel mutants: a kernel that changes only the returned State must fail the replay, on its N line"
kmut sigcMax_etPrev sigcMax "out.sigcMax = sigcMax;" "out.sigcMax = 2.0 * sigcMax + 1.0;" \
                            "out.etPrev = etc;" "out.etPrev = 1.5 * etc + 1e-3;"
kmut dwt dwt "    out.dwt = wt_impl - in.wt;" "    out.dwt = wt_impl + in.wt;"
echo "== all glibc-lane checks passed"
