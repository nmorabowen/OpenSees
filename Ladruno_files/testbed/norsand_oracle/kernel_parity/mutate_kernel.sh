#!/usr/bin/env bash
# WP-144 P1a: break the kernel-parity gate on purpose (material checklist, "A test can be GREEN because
# of the very bug"). Each mutant is a one-line edit of a COPY of SRC/material/nD/LadrunoNorSandKernel.h,
# built through NS_KERNEL_INCLUDE; every mutant must make kernel_parity FAIL. Linux / Esmeralda:
#   bash Ladruno_files/testbed/norsand_oracle/kernel_parity/mutate_kernel.sh      (from the repo root)
# Options (environment): NJOBS=<n> runs n mutants in parallel (default 1; each has its own build and pytest
# basetemp); MUT_ONLY="<tag> <tag> ..." runs only those mutants. The summary line "KILLED k / n" closes the output;
# a mutant that SURVIVES prints "== <tag>: SURVIVED".
# 2026-10-01: 13 / 13 mutants fail: the 9 P1a kernel mutants + the 4 chain mutants of the chained substep
# tangent (sheet 9.6): last_substep_tangent, chain_drop_Spi_carry, chain_drop_v_column, chain_assemble_1e-6_off.
# 2026-10-01 (G2 owner decision, exponential specific-volume update v_{n+1} = v_n exp(tr deps), vfac = v_{n+1}):
# vfac_v_not_v0 (now the correct code) is replaced by its reverse vfac_v0_not_v, and two v-law mutants are
# added: linear_v_update (the superseded v += v0 tr deps) and chain_Sv_v0 (S^v_{k+1} with v0, not v_{k+1}).
# 2026-10-01 after the v-law port: 15 / 15 mutants fail (Esmeralda).
# 2026-10-03 round 3b (energy option, p' floor, unified pi_i0; sheet §2.4, §5.4, §9.7, §10.2): the HAR->BA06
# swaps (killed by the K1.1h / K1.11 closed forms and by HAR parity, NOT by self-FD: A4), the floor mutants
# M-F1 ... M-F9 of sheet §9.7 that are kernel-level (M-F3a in its three forms, M-F3d in two), the (S.56) gate
# (ungated, W_ramp inverted), the (S.53) rule, and E_f. trial_tol_zero follows the p_ref scaling.
set -u
ROOT=$(cd "$(dirname "$0")/../../../.." && pwd)
H=$ROOT/SRC/material/nD/LadrunoNorSandKernel.h
cd "$ROOT/Ladruno_files/testbed/norsand_oracle"
NJOBS=${NJOBS:-1}
WORK=${TMPDIR:-/tmp}/nsmut_$USER
mkdir -p "$WORK"
MUTANTS=()
mut() { MUTANTS+=("$1"$'\x1f'"$2"$'\x1f'"$3"$'\x1f'"${4:-1}"); }

run_one() {
  local tag="$1" from="$2" to="$3" occ="$4"
  local d=$WORK/$tag
  rm -rf "$d"; mkdir -p "$d"; cp "$H" "$d/"
  python - "$d/LadrunoNorSandKernel.h" "$from" "$to" "$occ" <<'EOF'
import sys
p, a, b, occ = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
s = open(p).read()
i = -1
for _ in range(occ):
    i = s.find(a, i + 1)
    assert i >= 0, ("pattern not found (occurrence %d)" % occ, a)
s = s[:i] + b + s[i + len(a):]
open(p, "w").write(s)
EOF
  [ $? -eq 0 ] || { echo "== $tag: PATTERN NOT FOUND"; return; }
  local out
  out=$(NS_KERNEL_INCLUDE=$d python -m pytest kernel_parity -q -p no:cacheprovider -x -k "not self_check" \
        --basetemp="$d/pt" 2>&1 | tail -4)
  if echo "$out" | grep -qE "[0-9]+ failed|error"; then
    echo "== $tag: KILLED $(echo "$out" | grep -E "passed|failed|error" | tail -1)"
  else
    echo "== $tag: SURVIVED $(echo "$out" | grep -E "passed|failed|error" | tail -1)"
  fi
  echo "$out" | grep -E "^E |AssertionError|^FAILED" | head -2
}
export -f run_one
export H WORK

# ---- P1a / chain / v-law mutants ----------------------------------------------------------------------------------
mut spin_half_to_one 'spectral(V, atilde, g, 0.5, C4)' 'spectral(V, atilde, g, 1.0, C4)'
mut shear_col_no_sum 'C4[i][j][k][l] + C4[i][j][l][k]' 'C4[i][j][k][l]'
mut vfac_v0_not_v 'const double vfac = v;' 'const double vfac = n.v0;'
mut cap_drop_weta_qab 'fl.q_ab[a][b] = w * qu_ab[a][b] + w_eta * Y.eta_p / 3.0 * g_a[a];' 'fl.q_ab[a][b] = w * qu_ab[a][b];'
mut scan_step_2e-3 'constexpr double PI_SCAN_REL = 1.0e-3;' 'constexpr double PI_SCAN_REL = 2.0e-3;'
mut corner_branch_off 'if (std::fabs(s3) < CORNER_SIN3T) {' 'if (std::fabs(s3) < 0.0) {'
mut trial_tol_zero 'if (fl0.F <= F_TRIAL_TOL_REL * pRef(P)) {' 'if (fl0.F <= 0.0) {'
mut pistar_a_drop_POm 'const double pistar_a = (ps / (3.0 * p)) + pe.POm * pe.fl.Om_a[a];' 'const double pistar_a = (ps / (3.0 * p));'
mut substep_halvings_7 'constexpr int    MAX_SUBSTEP_HALVINGS = 8;' 'constexpr int    MAX_SUBSTEP_HALVINGS = 7;'
# sheet 9.6 chained substep tangent: the ladder returning the LAST sub-increment's CTO (the pre-2026-10-01
# contract) must be killed, and so must a chain that drops a term of (S.46) or is off by 1e-6 in the (S.47) assembly
mut last_substep_tangent 'run_fractions(P, n, deps, fr, m, m > 1, ro);' 'run_fractions(P, n, deps, fr, m, false, ro);'
mut chain_drop_Spi_carry '+ ((1.0 - ch.kappa) / ch.c) * c.S_pi[J]' '+ 0.0 * c.S_pi[J]'
mut chain_drop_v_column 'const double S_v = v_new * cum * chain_trE(J);' 'const double S_v = 0.0 * v_new * cum * chain_trE(J);'
# G2 exponential v-law: the superseded linear update, and the chain's S^v with v0 instead of v_{k+1}
mut linear_v_update 'const double v = n.v * std::exp(tr);' 'const double v = n.v + n.v0 * tr;'
mut chain_Sv_v0 'chain_propagate(cs, fr[k], ro.cur.st.v,' 'chain_propagate(cs, fr[k], n.v0,'
mut chain_assemble_1e-6_off 'tangent_small4(res.ae, res.sig, res.eps_e, V, Ae);' 'tangent_small4(res.ae, res.sig, res.eps_e, V, Ae); for (int i_ = 0; i_ < 3; ++i_) for (int j_ = 0; j_ < 3; ++j_) for (int k_ = 0; k_ < 3; ++k_) Ae[i_][j_][k_][k_] *= 1.0 + 1e-6;'

# ---- round 3b: the energy option (sheet 2.4; A4: killed by K1.1h/K1.11 and HAR parity, not by self-FD) -----------
# occurrences of 'if (P.energy == 1) {': 1 elastic(), 2 energy_psi(), 3 floor_ev(), 4 invert_elastic()
mut HAR_to_BA06_silently_elastic 'if (P.energy == 1) {' 'if (P.energy == 2) {' 1
mut HAR_to_BA06_energy_psi 'if (P.energy == 1) {' 'if (P.energy == 2) {' 2
mut HAR_to_BA06_invert_elastic 'if (P.energy == 1) {' 'if (P.energy == 2) {' 4
mut HAR_t4_term_as_AB06_eq64 'ratio = 3.0 * g * pa * w;' 'ratio = D22;'
mut MF8_pref_p0_under_HAR 'return P.energy == 1 ? P.p_a : std::fabs(P.p0);' 'return P.energy == 1 ? 1.0 : std::fabs(P.p0);'
# ---- round 3b: the p' floor, sheet 9.7 M-F1 ... M-F9 (kernel-level) ------------------------------------------------
mut MF1_floor_not_applied 'if (P.p_min <= 0.0) return;' 'if (P.p_min <= 0.0 || true) return;'
mut MF2_trial_events_not_counted 'o.st.eps_f_v = n.eps_f_v + o.res.dfv_tr + o.res.dfv_post;' 'o.st.eps_f_v = n.eps_f_v + o.res.dfv_post;'
mut MF2_n_f_tr_not_counted 'o.st.n_f_tr = n.n_f_tr + (o.res.floor_tr ? 1 : 0);' 'o.st.n_f_tr = n.n_f_tr;'
mut MF3a_Phi_tr_dropped 'fl_tr.active ? fl_tr.Phi : nullptr' 'nullptr'
mut MF3a_Phi_post_dropped 'fl_post.active ? fl_post.Phi : nullptr' 'nullptr'
mut MF3a_ae_not_at_floored_state 'fl_post.active ? el_c.ae : nullptr)' 'nullptr)'
mut MF3a_elastic_unprojected 'if (fl_tr.active) {                    // (S.32f) elastic' 'if (false) {                    // (S.32f) elastic'
mut MF3b_epsp_dropped '(1.0 / 3.0) * epsp * SQ23() * nh[b]' '0.0 * epsp * SQ23() * nh[b]'
mut MF3c_vcol_on_floored_trial ') + b[i][2] * Phi_tr[2][j]) - bs[i];' ') + b[i][2] * Phi_tr[2][j]) - bs[i] * ((Phi_tr[0][j] + Phi_tr[1][j]) + Phi_tr[2][j]);'
mut MF3d_chain_Phi_tr_omitted 'apply_op4(Op, T);' '(void)Op;'
mut MF3d_chain_Phi_post_omitted 'apply_op4(Op, c.S_eps);' '(void)Op;'
mut MF4_HAR_floor_with_BA06_inverse 'if (P.energy == 1) {' 'if (P.energy == 2) {' 3
mut MF5_trial_floor_skipped 'floor_project(P, eps_raw, fl_tr);' '{ Params P0_ = P; P0_.p_min = 0.0; floor_project(P0_, eps_raw, fl_tr); }'
mut MF6_projection_moves_deviator 'fr.eps_f[a] = eps_e[a] - dfv / 3.0;' 'fr.eps_f[a] = eps_e[a] - dfv / 3.0 - 1e-3 * (eps_e[a] - ev / 3.0);'
mut MF7_pi_altered_by_projection 'R.pi = pe.pi; R.dlam = pe.dlam;' 'R.pi = fl_post.active ? pe.pi * (1.0 + fl_post.dfv) : pe.pi; R.dlam = pe.dlam;'
mut MF9_floor_inside_local_Newton 'int e = elastic(P, eps_e, pe.el);' 'FloorResult ff_; floor_project(P, eps_e, ff_); int e = elastic(P, ff_.eps_f, pe.el);'
mut floor_E_f_bound_not_Psi 'R.ef_post = fl_post.active ? floor_energy(P, pe.eps_e, fl_post.eps_f, fl_post.dfv) : 0.0;' 'R.ef_post = fl_post.active ? P.p_min * fl_post.dfv : 0.0;'
mut floor_always_active 'if (elastic(P, eps_e, el) == EE_NONE                     // in dom Psi here (always EE_NONE)' 'if (elastic(P, eps_e, el) == EE_NONE && false'
# ---- round 3b: (S.56) gated to cap = smooth (A3), W_ramp's ratio (A2); the unified pi_i0 rule (S.53) ---------------
mut S56_ungated 'if (P.cap == 2 && detail::PI_SCAN_REL > W / 10.0 + 1e-15) {' 'if (detail::PI_SCAN_REL > W / 10.0 + 1e-15) {'
mut S56_W_ramp_inverted '1.0 - std::pow((1.0 - P.c2 * P.N) / (1.0 - P.c1 * P.N), (1.0 - P.N) / P.N);' '1.0 - std::pow((1.0 - P.c1 * P.N) / (1.0 - P.c2 * P.N), (1.0 - P.N) / P.N);'
mut S53_unified_rule_dropped 'eta = std::fmax(eta, c2 * P.M);' '(void)c2;'
mut S53_c2_for_cap_none 'const double c2 = (P.cap == 0) ? 0.0 : P.c2;' 'const double c2 = P.c2;'

if [ -n "${MUT_ONLY:-}" ]; then
  SEL=()
  for m in "${MUTANTS[@]}"; do
    t=${m%%$'\x1f'*}
    for o in $MUT_ONLY; do [ "$t" = "$o" ] && SEL+=("$m"); done
  done
  MUTANTS=("${SEL[@]}")
fi
OUT=$WORK/results.txt
: > "$OUT"
printf '%s\0' "${MUTANTS[@]}" | xargs -0 -P "$NJOBS" -I{} bash -c 'IFS=$'"'"'\x1f'"'"' read -r t f to oc <<< "$1"; run_one "$t" "$f" "$to" "$oc"' _ {} | tee -a "$OUT"
n=${#MUTANTS[@]}
k=$(grep -c ": KILLED" "$OUT")
echo "KILLED $k / $n"
grep -E ": SURVIVED|PATTERN NOT FOUND" "$OUT" || true
