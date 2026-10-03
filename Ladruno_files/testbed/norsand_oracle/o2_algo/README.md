# O2 — algorithmic oracle for LadrunoNORSAND (WP-144)

The backward-Euler return map of sheet `Ladruno_implementation/144a_norsand_equation_sheet.md`
§9 (AB06 Box 2) in Python/numpy: elastic predictor, trial yield check, 4 local unknowns
x = (ε₁ᵉ, ε₂ᵉ, ε₃ᵉ, Δλ) with the 4×4 Jacobian (S.30), nested scalar solve for π_i (S.27/S.37),
spectral decomposition, and the closed-form consistent tangent (S.31)–(S.34). This is the
algorithm the C++ kernel must reproduce to ~1e-10; every branch and tolerance here is the
kernel's contract. Written from the sheet only (independent of O1 and of the Deriver's scripts).

Sheet sections used: §1–§10 (algebra), §13 (K1 closed forms), §14 (K2 protocol, S.43/S.44).
Derivatives are the sheet's closed forms transcribed by hand (no sympy in the code path).

## Files
| file | content |
|---|---|
| `params.py` | `Params` (sheet §1.3 names), `validate()` with the owner-approved refusals |
| `kernel.py` | invariants (S.1, S.6, S.7), vertex rule §3.2, ζ WW/GA in y-form (S.8–S.11), BA06 energy (S.3–S.5), F and Q derivatives (S.12–S.21), cap (S.35–S.37), CSL (S.22), π_i* (S.23–S.24), nested π_i solve, residual/Jacobian (S.29–S.30), ã^ep (S.31–S.32), spectral tangents (S.33, S.34), chained substep tangent §9.6 (`chain_data` (S.45), `chain_propagate` (S.46), `chain_assemble` (S.47)) |
| `api.py` | `State`, `initial_state`, `step` (ladder + chained tangent), `step_fractions` (prescribed α_k, no ladder), `run_path`, `tangent`, `tangent_last_substep`, `tangent_finite`, drivers `triaxial`, `k2_path` |
| `acoustic.py` | `acoustic_min_det` (coarse (θ,φ) sweep + Nelder-Mead), (S.44) transcription for self-check |
| `selfcheck.py` | the self-checks (not the gate suite) |

## Run
From `Ladruno_files/testbed/norsand_oracle/`:
```
python -m o2_algo.selfcheck              # all groups
python -m o2_algo.selfcheck jac cto      # groups: jac | cto | newton | k1 | k2 | cap | chain  (k2 --table for the §14 sensitivity table)
```
Interface (identical in O1): `Params`, `State(sigma, eps_e, pi_i, v, D, eps_p_v, eps_p_s, flags)`,
`initial_state(params, sigma0, v0, pi_i0)`, `run_path(params, state0, deps[n,3,3])`,
`tangent(params, state)`, `acoustic_min_det(params, state, tangent4)`,
`triaxial(params, state0, kind, axial_strain_total, n_incr)`, `k2_path(params, state0, n_max)`.
Extra: `tangent_finite` (S.34) and `State.finite` for the K2 finite-strain protocol (diagonal,
fixed-direction paths only: ε̃ = ε^e_n + ln f; v = v₀J is then the same v-law as small strain, see below).

**Specific volume (sheet §1.2, G2 owner decision 2026-10-01, decision 1 = option b): exponential
update in both modes.** `api._step_once` sets v_{n+1} = v_n exp(tr Δε) (⇔ v = v₀ exp(tr ε)) and calls
`kernel.return_map(…, v, vfac)` with vfac = v = v_{n+1}: the trial-strain derivative in (S.31) is
∂v_{n+1}/∂ε̃_b = v_{n+1}, the converged specific volume of the step (not v_n, not v₀). The chain
(S.46) uses S^v_{k+1} = v_{k+1} (Σ_{j≤k} α_j) tr E_J with v_{k+1} the sub-increment's own converged v
(`chain_propagate(…, v_new=cur.v, …)`). v₀ is carried only as the committed datum (`initial_state`,
`State.v0`) and enters no derivative. This supersedes the G0/G1 linear rule v = v₀(1 + tr ε),
vfac = v₀, S^v = v₀ Σα tr E_J (BA06 Box 2 step 6b / 2.71): the two agree to first order in tr ε, and the
exponential form makes the LogStrain wrapper exact (v = v₀J). The small-strain and finite modes now
differ only in the tangent assembly ((S.33) vs (S.34)).

## Numerical contract (module constants in `kernel.py`)
| constant | value | meaning |
|---|---|---|
| `RES_TOL` | 1e-12 | local Newton converged when ‖(r₁,r₂,r₃, r₄/\|p₀\|)‖₂ ≤ RES_TOL |
| `MAX_LOCAL_ITERS` | 30 | then refusal `local_noconv` |
| `MAX_LINESEARCH` | 10 | halvings of the Newton step (Δλ included) until the scaled residual decreases; a nested failure or an inadmissible iterate at the trial point counts as "does not decrease"; else refusal `local_linesearch[:<nested reason>]` |
| `PI_TOL_REL` | 1e-12 | nested: \|r(π_i)\| ≤ PI_TOL_REL·\|π_{i,n}\| |
| `PI_SCAN_REL`, `PI_SCAN_MAX` | 1e-3, 1000 | nested solve selects the root **continuous with π_{i,n}**: scan from π_{i,n} in fixed steps of PI_SCAN_REL·\|π_{i,n}\| (away from 0 if r(π_{i,n}) > 0, toward 0 if r < 0) and stop at the FIRST sign change; if \|r\| grows between two scan points before a sign change the near root has annihilated in a fold (cap ramp) → `pi_fold`; PI_SCAN_MAX steps (travel \|π_{i,n}\|) without a sign change → `pi_nobracket` |
| `MAX_PI_ITERS` | 50 | safeguarded Newton inside that first bracket (start at the end with the smaller \|r\|; Newton step if strictly inside, else bisection; bracket updated by sign); else `pi_noconv` |
| `MAX_SUBSTEP_HALVINGS` | 8 | `api.step`: a refused increment is retried as 2, 4, …, 2^8 = 256 equal sub-increments (each a full BE step chained on the previous one); `flags['substeps']` = count used (1 = none); still refused at 2^8 → the whole-increment trial-elastic state is returned with reason `<finest reason> (substeps exhausted at 2^8)`; a substepped increment returns the CHAINED tangent of sheet §9.6 (below) |
| `R_TOL_REL` | 1e-8 | vertex rule: R < R_TOL_REL·\|p\| ⇒ n̂ = y_a = n̂_ab = y_ab = 0, Ω = Ω_a = 0, F := pη, f_a = F_p/3 |
| `CORNER_SIN3T` | 1e-8 | \|sin 3θ\| below ⇒ corner branch (S.9), ζ_yy := 0 |
| `F_TRIAL_TOL_REL` | 1e-10 | trial is plastic iff F(σ^tr, π_{i,n}) > F_TRIAL_TOL_REL·\|p₀\| |
| `EPS_S_TOL` | 1e-14 | ε_s below ⇒ n̂ᵉ := 0, q/ε_s := D₂₂ in (S.3) |
| `REPEATED_EIG_TOL` / `REPEATED_STRETCH_TOL` | 1e-10 | repeated-eigenvalue limits of g_ab (S.33) / γ̃_ab (S.34) |
| `PI_MAX_NEG` | 0 | any π_i iterate ≥ 0 is an evaluation failure |

Other contract points: p ≥ 0 or the §7 guard B ≤ 0 at an iterate is an evaluation failure
(`p_or_pi_nonneg`, `B_nonpos`) → line-search backtrack, refusal if at the trial state. A
converged Δλ < 0 is refused (`negative_dlambda`). On refusal `State.flags.refused = True` with
`reason`; the state returned is the trial-elastic one and `run_path` repeats it for the remaining
increments (nothing is raised mid-path). D = Δλ Σ σ_a q_a with converged, capped q_a (S.38).
Planar cap: w = 1 for η ≥ c₁M (cap for η < χ_cap M strictly, BA06 2.76). Params refusals:
N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄) else ValueError; warn on ρ > ρ̄; GA needs ρ, ρ̄ ∈ [7/9, 1]; WW needs
ρ, ρ̄ ∈ (½, 1] — **ρ = ½ exactly is refused** (owner decision 2026-10-01: at ½ the compression
corner of the Willam–Warnke section is a vertex).
The tangent is non-symmetric (major asymmetry 1–5 % measured): never a symmetric solver.

**Tangent on a substepped increment: the chained consistent tangent (sheet §9.6; owner decision
2026-10-01, supersedes "CTO of the last sub-increment").** When `api.step` had to substep
(`flags['substeps'] > 1`), `tangent()` returns the exact derivative of the increment's FINAL
stress with respect to the TOTAL strain increment Δε, chained through every accepted
sub-increment: `kernel.chain_data` takes the (S.45) blocks at each converged plastic sub-step
(b = J⁻¹, u = b t with t = (Δλ q_{a,π}, F_π), w = Π_xᵀ b[:, :3], κ = Π_x·u, c = r'(π_i), Π_v),
`kernel.chain_propagate` runs the recursion (S.46) on the six Δε columns (state per column: the
full 3×3 ∂ε^e/∂Δε_J, ∂π_i/∂Δε_J, and the shared cumulative fraction; ∂v_{k+1}/∂Δε_J = v_{k+1} Σα tr E_J
in closed form (v_{k+1} the sub-increment's converged v; was v₀ Σα tr E_J before G2); an elastic
sub-increment takes the elastic line), and `kernel.chain_assemble` forms
C = a^e(ε^e_m) : S^ε_m (S.47) with a^e in the (S.33) form on the eigen-data of the final ε^e_m.
Columns are the kernel's six tensor slots {00,11,22,01,12,02} with E_J = e_k⊗e_l + e_l⊗e_k on the
shear slots; the 3×3×3×3 returned has C[:,:,k,l] = C[:,:,l,k] = column/2 there, so `C:E = column`
and the parity reduction `C4_ijkl + C4_ijlk` is the kernel's 6×6 column exactly. Each operator
keeps the (S.33) row convention and symmetrises its own input (sheet §9.4/§9.6 contract). On a
NON-substepped increment (`substeps = 1`) the tangent is the closed-form CTO (S.31)–(S.33),
bit-identical to the oracle before §9.6; the chain for m = 1 reduces to it (8.6e-16 at distinct
trial eigenvalues; ≤ ~1e-8 inside the 1e-10 repeated-eigenvalue band, where the two limits
differ by the averaged-Φ-rows effect the sheet documents). A refused ladder level discards its
sensitivities and the finer level restarts from S₀ = 0. `api.step_fractions(P, st, deps,
fractions, chain=True)` takes the increment with prescribed fractions (any α_k > 0, Σ = 1, e.g.
a recursive-halving shape) and no ladder, chaining even m = 1 (that is how the reduction and the
FD checks below are measured); `api.tangent_last_substep` is the old last-sub-increment CTO,
kept only to quantify its error. `flags['pattern']` records the branch pattern ("PPEP"). Not
chained: finite mode (`State.finite`, the K2 diagonal protocol only) keeps the (S.34) tangent of
the last sub-increment — §9.6 is small strain (the LogStrain wrapper wraps the small-strain
kernel), and the finite-mode chain is not derived in the sheet.

**Why the nested solve scans (G1 defect, 2026-10-01).** With the smooth cap Ω = w(η)Ω^u depends
on π_i, and across the ramp entry η = c₁M the nested residual r(π_i) folds: the loop gain
√(2/3) h Δλ |π_i* − π_i| Ω^u w_η |η_π| reaches ≈ 0.96 and r is nearly flat over ~0.01|π_{i,n}|,
with a near root (w ≈ 6e-4, continuous with π_{i,n}) and a far root (w = 1). The earlier
factor-2 bracket from π_{i,n} enclosed the far root and the local Newton then failed its line
search (`test_g1_cap_smooth_runs_to_completion[O2-dev2e-3]`). The scan selects the near root;
when it has annihilated (|r| grows before any sign change) the iterate is rejected (`pi_fold`),
the local step is backtracked, and if the increment is still refused it is substepped.

## Self-check results (Esmeralda, numpy 2.2.6 / scipy 1.15.3; re-run 2026-10-01 after the scan/substep change)
- 4×4 Jacobian vs central FD (paper, fork, fork+smooth cap; 3 plastic states each): ≤ 3.9e-8 at
  h = 1e-8 (h = 1e-6 gives ~5e-7 purely from O(h²) truncation of r₄ in Δλ, J₄₄ ≈ 1e5).
- CTO vs FD of converged stress, non-coaxial shear, off corners: 4.4e-9 (paper WW), 4.2e-9
  (fork WW), 3.3e-9 (GA), 4.5e-8 (cap active), 1.3e-9 (elastic, repeated-eigenvalue start);
  all O(h²). Exact TXC corner, WW: 6.0e-4 / 6.3e-5 / 6.4e-6 at h = 1e-4/1e-5/1e-6 (O(h), §4.3);
  GA at the corner 2.3e-9.
- Local Newton, hard non-coaxial step: 2.9e0, 1.2e0, 3.1e-1, 2.9e-2, 2.3e-4, 1.5e-8, 1e-15
  (Δλ = 1.279e-3, 6 iterations; nested work 4690 scalar evaluations, the price of the fixed
  1e-3 scan on a step whose π_i root is far from π_{i,n}).
- K1: p(−0.01) = −271.828183; closed non-coaxial loop W/Σ|W| = 6e-17, state return 3e-14;
  ζ corners exact; refusals fire (incl. WW ρ = ½ and ρ̄ = ½); image point F = 9e-14; flow-rule
  identity 1.3e-14 at every plastic step; min D = +1.4e-2 (drained) / +1.1e-2 (undrained); H
  sign change: D − χψ_i = 5.8e-4 (fork) / 1.8e-4 (paper) with Δε_a = 5e-4; drained asymptote
  η → 1.408 (M 1.331) at 25 % strain, ψ_i → −2.7e-2; undrained CS p = −398.401 vs closed form
  −398.471 at ε_a = −2.0 (rel 1.8e-4; 1.1e-8 at −4.0), q/|p| = 1.33091 = M; paper set
  p = −780.85 vs −785.77 (rel 6.3e-3, still approaching: ψ_i = −8.3e-5), q/|p| = 1.20015.
- Finite mode: ã^ep vs FD 2.7e-9; (S.44) vs n·a·n contraction 2e-16.
- K2 nominal (π_{i,0} = −60.4, χ = −3.5, v_c0 = 1.81): ρ = 0.7/ρ̄ = 0.8 → n_first = 23
  (interpolated 22.40); ρ = ρ̄ = 1 → n_first = 27 (26.47). Ordering holds, gap 4 (paper 22/26).
  Unchanged by the nested-solve rewrite (no cap ⇒ r monotone ⇒ the same unique root).
- Smooth cap (G1.cap path, c₁ = 0.05, c₂ = 0.15, π_{i,0} = −80, ψ_{i,0} = −0.05), `selfcheck cap`:
  amp 1e-4 and 5e-4 complete at n = 40…640 with no substepping (5e-4: 169–937 nested evaluations
  per run, end w ≈ 5.9e-3); amp 2e-3 completes at every n: n = 40 substeps 30 of 34 plastic
  increments (max 4 sub-increments), n = 80 substeps 58 (max 2), n ≥ 160 none; min D = 0 (every
  D ≥ 0); endpoint π_i = −117.592 (n = 40) → −117.550 (n = 640), η plateau 0.0848, w = 0.063.
  Ramp entry at n = 40: step 9 Δλ = 7.8560e-4, π_i = −80.0131, η = 0.06471, w = 5.713e-4
  (identical to the monolithic 5-unknown solve of the triage); step 10 Δλ = 7.6155e-4,
  π_i = −80.1671; first substepped increment is step 11 (2), steps 12–40 need 4.
  Gate file `tests/test_g1_convergence_tangents.py`: 41/41 pass, including `[O2-dev2e-3]` and
  the n = 320/640 refinement (σ error 2.473e-4 → 1.237e-4, π_i 2.471e-4 → 1.236e-4: order 1.00).
  Planar / no-cap near-isotropic paths: refused (planar `local_linesearch`, no cap `local_noconv` since the G2
  exponential v-update, both `(substeps exhausted at 2^8)`), reported, as the gate requires.
- Chained substep tangent (§9.6), `selfcheck chain` (2026-10-01): max over the six kernel columns
  of ‖C_J − FD_J‖/‖C_J‖, central FD of the whole increment with the fractions held fixed.
  (A) AMP_STOP smooth cap n = 40, all 30 substepped increments of the ladder (step 11 m = 2 "PP",
  steps 12–40 m = 4 "PPPP"): 3.4e-6…1.1e-5 / 3.5e-8…8.6e-8 / 7.6e-10…1.2e-8 at h = 1e-6 / 1e-7 /
  1e-8 (O(h²) to the round-off floor); the last-sub-increment CTO is 0.51 (m = 2) to 0.79–0.90
  (m = 4) off. (B) generic non-coaxial plastic increment (fork WW, three shears, θ = 0.82), forced
  m = 8: 2.7e-8 / 1.3e-10 / 1.4e-9; m = 2: 2.6e-8 / 4.2e-10 / 1.7e-9 (last-sub CTO 2.1e-2 / 1.2e-2
  off). (C) m = 1: chain vs (S.33) 8.6e-16 (plastic), 0 (elastic); the ladder's m = 1 state carries
  no chain and `tangent()` is bit-identical to (S.33). (E) non-uniform α: (½,¼,⅛,⅛) on (B)
  2.7e-8 / 4.1e-10 / 1.0e-9; (¼,¼,¼,⅛,⅛) on AMP step 20 4.8e-6 / 5.2e-8 / 2.7e-9; (½,¼,¼) on AMP
  step 11 8.8e-6 / 6.6e-8 / 8.9e-9. (F) vertex (hydrostatic from the apex, π_i = −46.4758 frozen,
  exactly repeated eigenvalues): chain m = 1 vs (S.33) 4.6e-17; max|C:1|/max|C| = 3.8e-16 (m = 1),
  3.2e-16 (α = ½,¼,¼) with the FD along 1 at 0 and 2.0e-11.
  The full G1 suite is unchanged (see the gate record in `tests/`); the non-substepped tangent
  path is untouched by construction.

## Self-check re-run after the exponential v-update (Esmeralda, 2026-10-01, G2 owner decision; before → after)
All seven groups (jac, cto, newton, k1, k2, cap, chain) re-run with v_{n+1} = v_n exp(tr Δε), vfac = v_{n+1},
S^v = v_{k+1} Σα tr E_J. Every FD agreement stays at its O(h²)/round-off floor, which is the check that the
new v-factors are the consistent derivatives (the sheet's vfac = v₀ form would sit at ~1e-5…1e-4 here).
- Jacobian vs FD, max at h = 1e-8: 3.9e-8 → 4.5e-8 (paper state 4; the other eight states 5e-11…8e-9).
- CTO (S.33) vs FD at h = 1e-6: paper WW 4.44e-9 → 4.44e-9, fork WW 4.19e-9 → 4.19e-9, GA 3.31e-9 → 3.31e-9,
  cap active 4.49e-8 → 4.41e-8, elastic 1.30e-9 → 1.30e-9, TXC corner GA 2.30e-9 → 2.31e-9; WW corner O(h)
  unchanged (6.0e-4 / 6.3e-5 / 6.4e-6).
- Local Newton, hard step: same 7-iterate history to 1.14e-15 → 8.53e-16; Δλ 1.2792e-3 → 1.2794e-3,
  nested 4690 → 4690.
- K1: closed forms and refusals identical; flow-rule identity 1.3e-14 → 1.1e-14; min D unchanged
  (1.376e-2 / 1.589e-2); drained TXC path physics moves at the expected second order in tr ε:
  K1.8 η at 25 % strain 1.40836 → 1.40617 (fork), 1.22732 → 1.22727 (paper); K1.6 D − χψ_i at the H sign
  change 5.84e-4 → 5.14e-4 (fork), 1.84e-4 → 1.50e-4 (paper), same step (104 / 99). Undrained K1.7
  unchanged to all printed digits (isochoric: exp(0) = 1, v const 0.0 exactly).
- Finite mode: ã^ep vs FD 2.70e-9 → 2.70e-9 (unchanged by construction: finite mode already used this law);
  (S.44) 2.06e-16; K2 n_first 23 / 27, n_interp 22.396 / 26.470 unchanged.
- Smooth cap: substep counts, iteration counts and min D identical; endpoint π_i −117.5922 → −117.5881 (n = 40),
  −117.5495 → −117.5455 (n = 640): a 3.5e-5 relative shift, the v₀x²/2 effect at tr ε ≈ 8e-3.
- Chain (A) AMP_STOP 30 substepped increments: 3.4e-6…1.1e-5 / 3.5e-8…8.6e-8 / 2.0e-9…1.5e-8 at h = 1e-6 / 1e-7
  / 1e-8 (was …/ 7.6e-10…1.2e-8); last-sub CTO 0.51–0.90 unchanged. (B) m = 8: 2.74e-8 / 1.57e-10 / 9.8e-10
  (was 1.33e-10 / 1.44e-9); m = 2: 2.56e-8 / 3.5e-10 / 1.3e-9. (C) m = 1 vs (S.33) 1.0e-15 (was 8.6e-16),
  ladder m = 1 bit-identical: True. (E) (½,¼,⅛,⅛) 2.71e-8 / 3.4e-10 / 1.6e-9; AMP step 20 4.76e-6 / 5.1e-8 /
  3.9e-9; AMP step 11 8.77e-6 / 6.5e-8 / 9.7e-9. (F) vertex identical (4.6e-17, 3.8e-16, 3.2e-16, 2.0e-11).
  Logs: session scratchpad `o2_before.log` / `o2_after.log`.
