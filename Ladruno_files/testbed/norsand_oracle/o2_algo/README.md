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
| `kernel.py` | invariants (S.1, S.6, S.7), vertex rule §3.2, ζ WW/GA in y-form (S.8–S.11), BA06 energy (S.3–S.5), F and Q derivatives (S.12–S.21), cap (S.35–S.37), CSL (S.22), π_i* (S.23–S.24), nested π_i solve, residual/Jacobian (S.29–S.30), ã^ep (S.31–S.32), spectral tangents (S.33, S.34) |
| `api.py` | `State`, `initial_state`, `step`, `run_path`, `tangent`, `tangent_finite`, drivers `triaxial`, `k2_path` |
| `acoustic.py` | `acoustic_min_det` (coarse (θ,φ) sweep + Nelder-Mead), (S.44) transcription for self-check |
| `selfcheck.py` | the self-checks (not the gate suite) |

## Run
From `Ladruno_files/testbed/norsand_oracle/`:
```
python -m o2_algo.selfcheck              # all groups
python -m o2_algo.selfcheck jac cto      # groups: jac | cto | newton | k1 | k2 | cap  (k2 --table for the §14 sensitivity table)
```
Interface (identical in O1): `Params`, `State(sigma, eps_e, pi_i, v, D, eps_p_v, eps_p_s, flags)`,
`initial_state(params, sigma0, v0, pi_i0)`, `run_path(params, state0, deps[n,3,3])`,
`tangent(params, state)`, `acoustic_min_det(params, state, tangent4)`,
`triaxial(params, state0, kind, axial_strain_total, n_incr)`, `k2_path(params, state0, n_max)`.
Extra: `tangent_finite` (S.34) and `State.finite` for the K2 finite-strain protocol (diagonal,
fixed-direction paths only: ε̃ = ε^e_n + ln f, v = v₀J, ∂v/∂ε̃ = v).

## Numerical contract (module constants in `kernel.py`)
| constant | value | meaning |
|---|---|---|
| `RES_TOL` | 1e-12 | local Newton converged when ‖(r₁,r₂,r₃, r₄/\|p₀\|)‖₂ ≤ RES_TOL |
| `MAX_LOCAL_ITERS` | 30 | then refusal `local_noconv` |
| `MAX_LINESEARCH` | 10 | halvings of the Newton step (Δλ included) until the scaled residual decreases; a nested failure or an inadmissible iterate at the trial point counts as "does not decrease"; else refusal `local_linesearch[:<nested reason>]` |
| `PI_TOL_REL` | 1e-12 | nested: \|r(π_i)\| ≤ PI_TOL_REL·\|π_{i,n}\| |
| `PI_SCAN_REL`, `PI_SCAN_MAX` | 1e-3, 1000 | nested solve selects the root **continuous with π_{i,n}**: scan from π_{i,n} in fixed steps of PI_SCAN_REL·\|π_{i,n}\| (away from 0 if r(π_{i,n}) > 0, toward 0 if r < 0) and stop at the FIRST sign change; if \|r\| grows between two scan points before a sign change the near root has annihilated in a fold (cap ramp) → `pi_fold`; PI_SCAN_MAX steps (travel \|π_{i,n}\|) without a sign change → `pi_nobracket` |
| `MAX_PI_ITERS` | 50 | safeguarded Newton inside that first bracket (start at the end with the smaller \|r\|; Newton step if strictly inside, else bisection; bracket updated by sign); else `pi_noconv` |
| `MAX_SUBSTEP_HALVINGS` | 8 | `api.step`: a refused increment is retried as 2, 4, …, 2^8 = 256 equal sub-increments (each a full BE step chained on the previous one); `flags['substeps']` = count used (1 = none); still refused at 2^8 → the whole-increment trial-elastic state is returned with reason `<finest reason> (substeps exhausted at 2^8)` |
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

**Tangent on a substepped increment.** The consistent tangent is always the closed-form CTO
(S.31)–(S.34) at the converged root of the last backward-Euler solve. When `api.step` had to
substep, that is the CTO of the LAST sub-increment (its own trial strain, converged stress and
Δλ); it is a consistent linearisation of that sub-step, not of the whole increment (the chain
rule through the earlier sub-increments is not assembled). The C++ kernel does the same.

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
  Planar / no-cap near-isotropic paths: refused `local_linesearch (substeps exhausted at 2^8)`,
  reported, as the gate requires.
