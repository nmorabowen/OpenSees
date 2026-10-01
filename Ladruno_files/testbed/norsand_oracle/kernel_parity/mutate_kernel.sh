#!/usr/bin/env bash
# WP-144 P1a: break the kernel-parity gate on purpose (material checklist, "A test can be GREEN because
# of the very bug"). Each mutant is a one-line edit of a COPY of SRC/material/nD/LadrunoNorSandKernel.h,
# built through NS_KERNEL_INCLUDE; every mutant must make kernel_parity FAIL. Linux / Esmeralda:
#   bash Ladruno_files/testbed/norsand_oracle/kernel_parity/mutate_kernel.sh      (from the repo root)
# 2026-10-01: 9 / 9 mutants fail.
set -u
ROOT=$(cd "$(dirname "$0")/../../../.." && pwd)
H=$ROOT/SRC/material/nD/LadrunoNorSandKernel.h
cd "$ROOT/Ladruno_files/testbed/norsand_oracle"
run() {
  local tag="$1" from="$2" to="$3"
  local d=${TMPDIR:-/tmp}/nsmut_$USER/$tag; mkdir -p $d; cp $H $d/
  python - "$d/LadrunoNorSandKernel.h" "$from" "$to" <<'EOF'
import sys
p, a, b = sys.argv[1:]
s = open(p).read()
n = s.count(a)
assert n >= 1, ("pattern not found", a)
s = s.replace(a, b, 1)
open(p, "w").write(s)
EOF
  local out
  out=$(NS_KERNEL_INCLUDE=$d python -m pytest kernel_parity -q -p no:cacheprovider -x -k "not self_check" 2>&1 | tail -4)
  echo "== $tag: $(echo "$out" | grep -E "passed|failed|error" | tail -1)"
  echo "$out" | grep -E "^E |AssertionError" | head -2
}
run spin_half_to_one 'spectral(V, atilde, g, 0.5, C4)' 'spectral(V, atilde, g, 1.0, C4)'
run shear_col_no_sum 'C4[i][j][k][l] + C4[i][j][l][k]' 'C4[i][j][k][l]'
run vfac_v_not_v0 'const double vfac = n.v0;' 'const double vfac = v;'
run cap_drop_weta_qab 'fl.q_ab[a][b] = w * qu_ab[a][b] + w_eta * Y.eta_p / 3.0 * g_a[a];' 'fl.q_ab[a][b] = w * qu_ab[a][b];'
run scan_step_2e-3 'constexpr double PI_SCAN_REL = 1.0e-3;' 'constexpr double PI_SCAN_REL = 2.0e-3;'
run corner_branch_off 'if (std::fabs(s3) < CORNER_SIN3T) {' 'if (std::fabs(s3) < 0.0) {'
run trial_tol_zero 'if (fl0.F <= F_TRIAL_TOL_REL * std::fabs(P.p0)) {' 'if (fl0.F <= 0.0) {'
run pistar_a_drop_POm 'const double pistar_a = (ps / (3.0 * p)) + pe.POm * pe.fl.Om_a[a];' 'const double pistar_a = (ps / (3.0 * p));'
run substep_halvings_7 'constexpr int    MAX_SUBSTEP_HALVINGS = 8;' 'constexpr int    MAX_SUBSTEP_HALVINGS = 7;'
