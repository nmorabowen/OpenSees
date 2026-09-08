---
title: ADR 97 verdict — ASDPlasticMaterial3D closest-point return map
project: Ladruno
status: verdict shipped (opt-in) — default flip is a separate PR
owner: nmora
tags: [implementation, material, review, verdict]
---

# ADR 97 verdict — `ASDPlasticMaterial3D` closest-point return map

Build at closeout: `7e93e4381` for every measurement in §1-§4 below (HEAD
`5e89b8a12` — the merge of P4 with P6 — is SRC-identical to it,
`git diff 7e93e4381 HEAD -- SRC` empty); the banner line (§ deliverable 2) and
the final full-battery re-run were verified on the subsequent build
`33c670bf5` (banner-only source change, `ops.ladrunoBuild()` verified on both).
Plan:
[[97_ladruno_asdp_closest_point_adr]]. Stack: #817 (`wp/97a-plan-oracles`, P0) →
#819 (`wp/97b-cp-smooth`, P1) → #824 (`wp/97c-cp-principal`, P2) → #825
(`wp/97d-cp-hoekbrown`, P3) → #827 (`wp/97f-explicit-gate`, P5) → #829
(`wp/97e-numalg-repoint`, P4), plus #826 (`wp/97g-measure`, P6, based on #824).
This document is P7 (`wp/97h-closeout`), based on #829 with #826 merged in.
Reports: [[reviews/adr97_p1_report]]…[[reviews/adr97_p5_report]],
[[reviews/adr97_p6_measurement]]. Mutation records:
[[reviews/adr97_p1_mutation]]…[[reviews/adr97_p4_mutation]]. Oracles:
`Ladruno_implementation/adr97_oracle/` (`README.md` is the P0 findings ledger).

## §0 Verdict

**Warrant-grade for the 23 of 46 specializations it supports, shipped strictly
opt-in.** `Closest_Point` + `Algorithmic` answer ADR-94's M3 finding — no
shipped tangent was the tangent of the committed map, and `Backward_Euler` is
an Ortiz–Simo cutting plane, not a closest-point return — for VonMises,
Drucker–Prager (including the apex), MohrCoulomb, and MohrCoulombTensionCutoff,
each verified against an independent numpy oracle (complex-step Jacobians,
exact to machine precision) *and* against a central difference of the binary's
own assembled residual on a free-DOF finite-element rig, with no shared code
between the two checks. `Backward_Euler` + `Secant` — every result the fork has
published on this material was obtained on that pair — is untouched: 23 decks /
282 committed-stress rows byte-identical in fresh subprocesses at every phase,
with exactly ONE deliberate, root-caused exception (below).

Four mutation gates (drop one term of the consistent tangent, or the
Numerical_Algorithmic repoint itself, on a scratch build never committed) each
turned the load-bearing tests red at the historically-recorded magnitude of the
defect they exist to catch, and left every test that should be untouched green
— the sharpest evidence available that the test suite measures what it claims
to. P6's mesh-scale measurement (§3) found `Closest_Point` converges strictly
more of a bearing-capacity load history than `Backward_Euler` under identical
criteria, at equal-or-fewer total Newton iterations, on two material families
at two mesh scales — and recommends the default flip. That flip is explicitly
**not** part of this warrant: D1 keeps `Backward_Euler` + `Secant` as the
shipped default, and the recommendation is the owner's decision in a separate
PR (§3).

**The one documented exception to byte-identity.** `test_adr97_p4_inertness.py`'s
`cube/vm/BE/Numerical_Algorithmic_FirstOrder/plastic` baseline entry moved by up
to 5.428e-09 absolute in 40/60 committed-stress components when P4 re-pointed
`Numerical_Algorithmic_*` at the actual committed map. Root-caused, not assumed:
that deck is load-controlled with free DOFs, and `Backward_Euler` calls
`ComputeTangentStiffness()` at the tail of *every* `setTrialStrainIncr()` — every
outer-Newton trial, not only the converged one — so a materially different
*returned tangent* (now a real FD of `Backward_Euler` itself, instead of the old
FD of an unrelated map) changes which of several `NormDispIncr`-tolerance-
equivalent states the outer Newton lands on. Two independent signatures confirm
this is a Newton-path effect, not a code-path regression: (1) the same deck's
purely elastic leg is bit-identical (an elastic trial is exactly linear, so old
and new FD reproduce the same analytical `E` regardless of which map is
differentiated); (2) the magnitude (5.428e-09) sits at the same order as the
file's own documented cross-platform compiler-noise floor (5.8e-09, MSVC vs
GCC/libm). All 22 other decks in the same 23-deck baseline — every other
`cube/vm/BE/*` combination, every fully-prescribed tet deck, all four
explicit-integrator decks — regenerated bit-identical. `Backward_Euler`'s own
source has zero calls to `compute_local_stress()` remaining anywhere (grep-
confirmed): D1's byte-identity promise holds for the *map*; this one baseline
entry's *numbers* moved because they exercise the (correctly re-pointed)
`Numerical_Algorithmic_*` tangent option, not `Backward_Euler`'s return.

Scope: **23 of 46 registered specializations** (VonMises, Drucker–Prager,
MohrCoulomb, MohrCoulombTensionCutoff — smooth + principal-space families).
StiffSoil (two-surface) and RoundedMohrCoulomb stay refused at parse time (no
runtime coverage existed before this ADR either, per ADR-94 §5); the 23 mixed
YF/PF pairings that do not match a supported family are refused at parse time
and pinned by `static_assert` in a g++ pre-flight. Every refusal names the
family and cites this ADR.

## §1 Verification manifest

One row per gate per WP. "Test" is the file that carries the assertion;
"oracle" is the independent reference; run ids are the **final** Zone-A
(Ubuntu) conclusion for that branch (an earlier failing run exists for #819,
#825 and #826 — see the gate-7 note below the table).

| WP | Gate | Measured | Test | Oracle | Zone-A run (final) |
|---|---|---|---|---|---|
| P1 (#819) | 1 — correctness | committed stress vs oracle 7.4e-14…9.2e-9 (VM/DP, 6 cases); `\|f\|` ≤ 7.1e-15; local Newton ≤ 3 iters (gate ≤5) | `tests/test_adr97_p1_smooth.py` | `adr97_oracle/cppm_vm.py`, `cppm_dp.py` | [34177905021](https://github.com/nmorabowen/OpenSees/actions/runs/34177905021) |
| P1 | 2 — tangent | `Algorithmic` vs FD of the binary's own residual 1.51e-11 (4 DOF) / 2.15e-10 (12 DOF), vs pinned `Continuum` negative control 0.573447; two-cube `testIter` 16 vs 113 (7.1×) | same file | `adr97_oracle/fd_tangent_driver.py` | same run |
| P1 | 3 — path independence | AF step-refinement error ratio CP:BE = 1:4.07 (predicted 4.3×); linear hardening agrees to 1.45e-16 | same file | `adr97_oracle/path_independence.py` | same run |
| P1 | 4 — BE inertness | 23 decks / 282 rows byte-identical (fresh subprocess); CP≡BE 0.0…3.8e-17 non-rotating, CP≠BE 2.3e-2 on AF | `tests/test_adr97_p4_inertness.py` | pre-change dumps, same binary | same run |
| P1 | 5 — mutation | drop `dl·dm/dσ` term of Ξ → `Algorithmic` becomes exactly `Continuum` (0.573447/0.670133 reproduced to the digit); AF local Newton 3→13 iters; 4 killed / 39 survivors | scratch build, [[reviews/adr97_p1_mutation]] | — | — |
| P1 | 6 — fail-loud | 13/13: starved refusal, `Algorithmic` cross-refusal (D2), unconverted-family refusal, unknown-token refusal, `strict_convergence` inertness, B4 hydrostatic-tension no-NaN | `tests/test_adr97_p6_failloud.py` | — | same run |
| P1 | 7 — portability | first run [34174010743](https://github.com/nmorabowen/OpenSees/actions/runs/34174010743) FAILED on 3 decks (platform-dependent float pins); fixed by `94abcdbb7` (gate-4 byte-identity scoped exact-on-win32 / 1e-6-relative-elsewhere) | — | — | **34177905021 (success)** |
| P2 (#824) | 1a/1b/1c | 6 oracle regions to 1.5e-16…9.1e-16; path `\|f\|` ≤ 2.0e-14; MCTC bit-identical to BE (gap 0.0) on hydrostatic-tension and Rankine-face | `tests/test_adr97_p2_principal.py` | `adr97_oracle/cppm_mc.py`, ADR-84 MCTC decks | [34177906457](https://github.com/nmorabowen/OpenSees/actions/runs/34177906457) |
| P2 | 2 — tangent | `Algorithmic` vs FD 2.88e-11 (degenerate edge), 8.40e-9/1.08e-8 (face, free-node rig); BE fails to converge on the same rig; global Newton 16 vs 62 (MC), 22 vs 74 (MCTC) | same file | `adr97_oracle/fd_tangent_driver.py` (new free-node rig) | same run |
| P2 | 4 — BE inertness + finding | 23/282 byte-identical; CP step-size-independent (N=1..40, 1e-16 floor); CP≡BE at apex 2.1e-16; **finding:** BE's analytic MC gradient (`MC_ds=0`) is 2.9e-1 off the exact return, only its own FD (`MC_ds>0`) agrees — pinned both directions, not fixed (D1) | same file | same | same run |
| P2 | 5 — mutation | eigenprojection rotation term → 1.0 on the non-degenerate branch: sheared-face FD error 1.08e-08→1.29 (129%); axis-aligned face rig stops converging; 4 killed / 34 survivors | scratch build, [[reviews/adr97_p2_mutation]] | — | — |
| P2 | 6 — fail-loud | 12 tests: 6 mixed-pairing refusals, `MC_phi==0` refusal, 3 hydrostatic/near-degenerate no-NaN commits, `Algorithmic` cross-refusal, `strict_convergence` inertness | same file | — | same run |
| P2 | 7 — portability | first run failed with the same gate-4 platform issue as P1 (pre-`94abcdbb7`); final green | — | — | **34177906457 (success)** |
| P3 (#825) | 1a/1b/1c | 7 oracle regions (incl. the header's mis-called "apex") to 7.3e-16…1.85e-13; `\|f\|` ≤ 7.3e-11; exact closed-form tension plateau 244.4854419 kPa (1.16e-16 rel.), sharper than ADR-94's own 5e-2 pin on `T` | `tests/test_adr97_p3_hoekbrown.py` | `adr97_oracle/cppm_hb.py` | [34183565868](https://github.com/nmorabowen/OpenSees/actions/runs/34183565868) |
| P3 | 2 — tangent | `Algorithmic` vs FD 3.34e-9/3.87e-9 (face, curvature+rotation live), 1.49e-10 (degenerate edge); BE fails to converge (both `Continuum` and `Secant`) where CP converges in 16 iters | same file | `adr97_oracle/fd_tangent_driver.py` | same run |
| P3 | 4 — BE inertness + gap pinned | 23/282 byte-identical; potential-gap pinned both directions: CP plastic `eps_vol` +2.208e-05 vs BE +2.09e-13 (5.91e-2 relative stress gap = 26.2% of strength scale) | same file | same | same run |
| P3 | 5 — mutation | drop curvature term of `dm/dy` → face rig stops converging, edge-rig FD error 1.49e-10→8.31e-3 (5.6e7×); local Newton 4/4/3/5→7/8/7/7; 11 killed / 21 survivors | scratch build, [[reviews/adr97_p3_mutation]] | — | — |
| P3 | 6 — fail-loud | all 7 deck-reachable mixed HB pairings refused (both directions), 4 more by `static_assert`; matched pairing accepted; `strict_convergence` inertness; no-NaN at hydrostatic/near-degenerate tension | same file | — | same run |
| P3 | 7 — portability | green | — | — | **34183565868 (success)** |
| P4 (#829) | 2/6 — repoint correctness | `Numerical_Algorithmic_{First,Second}Order` now differentiate the actual committed map: BE vs FD 2.06e-8/3.39e-8; CP `Algorithmic` vs `Numerical_Algorithmic_SecondOrder` 1.21e-8 (all ≤ 1e-6 gate); two-cube `testIter` sum 16 (was 113 pre-fix); H6 4.6%→~2-3e-8 | `tests/test_adr97_p4_numalg.py`, `test_adr94_hlist_numerics.py::test_H6...` | `adr97_oracle/fd_tangent_driver.py` | [34188231854](https://github.com/nmorabowen/OpenSees/actions/runs/34188231854) |
| P4 | 4 — BE inertness | 22/23 decks byte-identical; the one documented exception is analyzed in §0 | `tests/test_adr97_p4_inertness.py` | pre-change dumps | same run |
| P4 | 5 — mutation | revert to `compute_local_stress()` FD → error returns to 4.574e-2 (historical 4.6%); 4 of 6 tests + H6 killed, 2 refusal-propagation tests correctly survive | scratch build, [[reviews/adr97_p4_mutation]] | — | — |
| P4 | 6 — fail-loud (refusal propagation) | `Backward_Euler` P2a exhaustion deck refuses identically under `Secant`/`Numerical_Algorithmic_SecondOrder`; `Closest_Point` P1 starved-VM reproducer still refuses under `Numerical_Algorithmic_SecondOrder` | same file | — | same run |
| P4 | 7 — portability | green (resolved after the in-progress check during drafting; re-confirmed) | — | — | **34188231854 (success)** |
| P5 (#827) | 6 — D5 gate | 4 explicit integrators gated behind `experimental_integrator`, checked both directions + typo rejection + `Algorithmic` cross-refusal on an opted-in explicit method | `tests/test_adr94_hlist_mechanical.py`, `test_adr94a_fail_loud.py`, `test_adr97_p6_failloud.py` (appended) | — | [34184110377](https://github.com/nmorabowen/OpenSees/actions/runs/34184110377) |
| P5 | M9 fix | `StiffSoilShear` `cot(phi)` NaN at `phi==0` fixed algebraically (~1e-14 rel. match for `phi!=0`); PF 0/0 at hydrostatic axis fixed with a zero-guard; both regression-tested (finite at step 1) | `tests/test_adr97_p5_stiffsoil.py` | standalone g++ probe, ADR-94 R3a re-run | same run |
| P5 | 4 — BE inertness (incl. explicit decks) | 23/282 byte-identical, including the 4 explicit-integrator decks now opted in via the gate-4 dumper | `tests/test_adr97_p4_inertness.py` | pre-change dumps | same run |
| P5 | 7 — portability | green | — | — | **34184110377 (success)** |
| P6 (#826) | measurement | `Closest_Point` converges strictly more of the load history than `Backward_Euler` at equal-or-fewer total Newton iterations, on MC and MCTC, at 3993 and 24000 DOF (§3 table) | `tests/test_adr97_p6_measure_smoke.py` (smoke only; full sweep is a driver, not a pytest gate) | `adr97_oracle/measure_p6.py` | [34183306667](https://github.com/nmorabowen/OpenSees/actions/runs/34183306667) |
| P6 | 7 — portability | green (smoke test only; the mesh-scale sweep itself is not a CI gate — minutes-to-hours per configuration) | — | — | **34183306667 (success)** |

**Gate-7 caveat, all branches.** Zone-A (Ubuntu) is the required check; the
self-hosted Zone-B / full-suite jobs on several of these runs stay `queued`
indefinitely on this fork's current runner capacity (documented behavior, not a
failure) — the table above reports the Zone-A (Ubuntu) job's conclusion only,
per the orchestrator's instruction. P4 (#829, run 34188231854) was still
`in_progress` at the time this table was first drafted; re-checked after the
P7 build/battery — **Zone-A (Ubuntu) resolved `success`** — so all seven
branches in the stack are green.

## §2 Mutation-gate record

Four independent mutations, each on a scratch build never committed, each
reverted and re-measured identical before the WP's own PR was pushed.

| WP | Term dropped | Tests turned red | Magnitude | Survivors | Interpretation |
|---|---|---|---|---|---|
| P1 | `dl·dm/dσ` term of the algorithmic elastic modulus Ξ, in `cp_assemble`'s `J_ss` block | `test_gate2_algorithmic_matches_a_finite_difference_of_the_committed_map`, `..._sheared_twelve_dof_rig`, `..._two_cube_newton_cost`, `test_gate1_newton_converges_in_at_most_five_iterations` (4) | FD error 1.51e-11→0.573447, 2.15e-10→0.670133 — **exactly** the pinned `Continuum` values, to the digit; AF Newton 3→13 iters | 39/43 | Proves Ξ, not just "a" term, is what separates `Algorithmic` from `Continuum` — Ξ→E is the whole content of ADR-94 M3. Committed stress untouched (residual not mutated) on all 39. |
| P2 | eigenprojection rotation term `T[3+s,3+s] = (y_i−y_j)/(x_i−x_j)` → 1.0 (non-degenerate branch only) | 4 of 38: 2 face FD tests, the degenerate-edge FD test (unpredicted — see below), the BE-cannot-do-face negative control | Sheared face FD error 1.08e-8→1.29 (129%); axis-aligned face rig stops converging (`analyze→-3`); degenerate-edge 2.88e-11→8.74e-2 | 34/38 | Rotation term is most of the shear block of `C_alg`, not a correction. Degenerate-edge test also dies because the oedometric rig's `s1==s2≠s3` state has only ONE of its three shear slots on the l'Hôpital branch — the other two hit the mutated branch; the pure l'Hôpital path is therefore verified only in combination (recorded, not a gap in coverage that can be closed without a hydrostatic rig that has nothing to measure). |
| P3 | curvature term of `dm/dy` (`A.col(i) += ... contrib`) in `hb_assemble`'s Jacobian | 11 of 32: the 3 face/edge FD tests plus 8 others reading `cp_iterations` or convergence | Face rig stops converging entirely; edge-rig FD error 1.49e-10→8.31e-3 (5.6e7×); local Newton 4/4/3/5→7/8/7/7 | 21/32 | Confirms the P3 tests measure surface *curvature*, not just that a return lands on the surface — committed stress unchanged (1.6e-15…1.9e-13 vs oracle) because the residual is untouched. |
| P4 | reverted `numerical_tangent_of_committed_map()`-based wrappers to pre-P4 `compute_local_stress()` FD | 4 of 6 in `test_adr97_p4_numalg.py` + `test_H6_no_tangent_option_reproduces_the_consistent_tangent` (5 total) | FD error returns to 4.574e-2, matching the historically-recorded pre-fix 4.6% (wp/94c) almost to the digit | 2/8 (the 2 refusal-propagation tests, correctly orthogonal) | Confirms the re-point, not merely "some" change to the wrapper, is what the suite measures — the mutated number reproduces the exact historical defect. |

No mutation gate was run for P5 (D5 is a parser-level gate and M9 is two
narrowly-scoped arithmetic fixes — not new algorithmic machinery of the kind
P1–P4's gates target; recorded as scoped out, not an omission).

## §3 P6 measurement and the D1 default-flip recommendation

Full data: [[reviews/adr97_p6_measurement]]. Summary (MC deck, E=30000 kPa,
φ=32°, ψ=8°, c=15 kPa; MCTC deck is the Cerro-Lindo-like EDZ material from
`tests/test_asdplastic_mctc.py`; both non-associated, `system UmfPack`,
`numberer RCM`, `NormDispIncr 1e-8`):

| Family | Scale (DOF) | `BE_Secant` (shipped default) | `CP_Algorithmic` | `CP_Algorithmic_Krylov` |
|---|---|---|---|---|
| MC | 3993 | 7/10 steps, 60 iters | **10/10**, 39 iters | **10/10**, 53 iters, fastest wall-clock |
| MC | 24000 | 7/12 steps, 56 iters | **12/12**, 48 iters | **12/12**, 66 iters, fastest wall-clock |
| MCTC | 3993 | 6/10 steps, 112 iters | **9/10**, 90 iters | **9/10**, 70 iters, fastest wall-clock |
| MCTC | 24000 | 6/12 steps, 113 iters | **11/12**, 140 iters | **11/12**, 117 iters, fastest wall-clock |

On every configuration and scale measured, `Closest_Point` (with either tangent)
converges strictly more of the load history than `Backward_Euler`/`Secant`
under identical convergence criteria, and the committed-stress gap on steps
both maps complete stays under 1% (0.02%–0.68%, shrinking with mesh refinement
— consistent with the known small divergence when the flow direction rotates
mildly over a step, not a correctness gap). `BE_Continuum` never recovers the
steps `BE_Secant` misses either — the failure is the cutting-plane *map*, not
the tangent quality feeding Newton.

**Recommendation (the report's own §6, restated here for the record): flip the
default to `integration_method Closest_Point` + `tangent_type Algorithmic`,
leaving `algorithm KrylovNewton` as the (already-default) algorithm.**

**This is stated plainly as the owner's decision, not this PR's action, and it
is scoped:**

1. **Scope to the 23/46 supported specializations.** VonMises, Drucker–Prager,
   MohrCoulomb, MohrCoulombTensionCutoff are covered; the flip cannot be
   unconditional while StiffSoil/RoundedMohrCoulomb/the 23 mixed pairings are
   refused at parse time under `Closest_Point` — either those families need
   their own coverage first, or the default flip needs to be conditional on
   family support (a parser change, not just an option-default change).
2. **One geometry, one mesh regularity, one solver.** The measurement used a
   strip-footing bearing-capacity model on a regular hex mesh with `UmfPack`;
   the ADR-75 PARDISO path and MPI decomposition (`sendSelf`/`recvSelf`, ADR-94
   D3) were not exercised for `Closest_Point`'s per-tag option maps.
3. **Not a correctness re-verification of every downstream result.** The flip
   changes the *default* a deck gets with no explicit `integration_method`/
   `tangent_type` tokens; every published fork result (Cerro Lindo, ADR-84,
   the D3 finite-strain oracle) was obtained on explicit `Backward_Euler` +
   whatever tangent those decks name, and stays reproducible either way — but
   any *future* deck that omits the tokens would silently change behavior the
   day the flip ships, which is exactly why D1 gates it on a separate PR.

## §4 Owner items found on the way, not fixed under D1

D1 forbids changing `Backward_Euler`'s answer in this ADR (every published
fork result depends on it staying byte-identical). Each item below was found,
quantified, and pinned as a test that will turn red the day the shipped code
changes — they are handed to the owner as their own separate PRs, each already
scoped by the report/test that found it.

| # | Item | Where it lives | Quantified | Pinned by |
|---|---|---|---|---|
| 1 | `HoekBrown_PF::g` is evaluated in the un-negated (tension-positive) frame while `HoekBrown_YF` negates first; its branch always falls into the compressive `else`, so `g` collapses to a **Tresca** potential — `HB_mb_psi` is inert, flow is exactly non-dilatant, 32.86° off the frame-consistent normal even at `mb_psi==mb`, and no shipped flow direction has a return to the apex past the tensile corner (the mechanism behind the ADR-94 H10a residual) | `PlasticFlowDirections/HoekBrown_PF.h:66` | plastic `eps_vol` +2.208e-05 (frame-consistent) vs BE's +2.09e-13 (zero dilation); 5.91e-2 relative stress gap = 26.2% of strength scale | `tests/test_adr97_p3_hoekbrown.py` (gate-4 potential-gap test, both directions); [[reviews/adr97_p3_report]] §"Also fixed"/finding; `adr97_oracle/README.md` findings 5–7; `LEDGER_quirks` "Hoek-Brown (ADR-97 P3)" |
| 2 | `MohrCoulomb_YF`/`_PF`'s shipped **analytic** Lode-angle gradient (`MC_ds=0`, the branch every deck in the repo actually uses) is wrong: it agrees with the exact closest-point return only through the header's OWN finite difference (`MC_ds>0`) | `YieldFunctions/MohrCoulomb_YF.h`, `PlasticFlowDirections/MohrCoulomb_PF.h` (the `c1/c2/c3` branch) | `MC_ds=0`: 2.9e-1 off, step-size dependent; `MC_ds=1e-4`: 2.1e-14 | `tests/test_adr97_p2_principal.py::test_gate4_backward_euler_agrees_only_through_its_own_finite_difference` (both directions); `LEDGER_quirks` "Mohr-Coulomb principal-space return (ADR-97 P2)" |
| 3 | `DruckerPrager_YF`'s cohesion internal variable is commented out of `f` itself (`DruckerPrager_YF.h:26`) while `yf_hardening` still contributes `df/dk = -1` for it — a hardening term matching no term of the actual yield function | `YieldFunctions/DruckerPrager_YF.h:26` | DP-with-cohesion-hardening is perfectly plastic under `Closest_Point` (uses the true `df/dq=0`) and hardening under `Backward_Euler` (uses the phantom term); IV itself still grows (`k=0.190036` in the P1 probe) and is read by nothing | [[reviews/adr97_p1_report]] §3 "P0 findings honoured, not fixed"; `adr97_oracle/README.md` finding 2 |
| 4 | `ArmstrongFrederick`'s hardening policy adds the ENGINEERING-shear (doubled) flow direction to a STRESS-like back stress, and has no `2/3` factor on `h_a` (commented out) | `AllASDHardeningFunctions.h:157-160` | rotating-normal path: header convention gives `σ12=23.78/σ23=37.77`, tensor-consistent gives `13.49/29.23` — 43%/23% difference on a sheared path | `adr97_oracle/README.md` finding 1; both `Closest_Point` and `Backward_Euler` share the policy, so P1 mirrors the shipped (wrong) convention under D1 |
| 5 | `Backward_Euler`'s DP-apex and MC-apex region tests (`check_apex_region`) are **Euclidean** normal-cone tests where the exact boundary is in the **elastic metric** — the header says so in its own comments | `YieldFunctions/DruckerPrager_YF.h`, `MohrCoulomb_YF.h` `check_apex_region` | DP: exact slope `(K·η̄+h_k/η)/G` vs header `η` — misclassifies in BOTH directions (e.g. inadmissible negative √J2 committed at `(p-p_apex)/q=0.36` under `η̄=0.2`); HB: header always over-claims, 32/400 scanned trials, up to 122.6 kPa lost strength | `adr97_oracle/README.md` DP block + finding 3, HB block + finding 7; `LEDGER_quirks` "Euclidean vs elastic-metric apex classification" (wp/94c, re-quantified here) |
| 6 | `StiffSoilCap`'s uninitialized cap internal variable `pc0==0` gives an anomalous first response (wrong physics, not NaN) | `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/StiffSoilCap_YF.h` (via the pre-investigation's `stiffsoilcap_drive.py` probe) | not quantified beyond "anomalous," out of scope for this ADR's M9 pass | [[reviews/adr97_p5_report]] §4 "Not done / left for the owner" |
| 7 | (carried from ADR-94, unchanged by this ADR) `Runge_Kutta_45_Error_Control_old` still contains a bare `exit(-1)` on a NaN guard — if this integrator is ever reachable again it terminates the whole process rather than returning a material failure code | `ASDPlasticMaterial3D.h:5385,5576` | n/a — the integrator is refused at parse time (ADR-94 precedent, reaffirmed by ADR-97 D5's own text) so the line is currently unreachable, not exercised | [[reviews/adr94_verdict]] §1 M2/row 109; not re-tested here because ADR-97 does not touch `RK45_old`'s refusal |

## §5 Residual risks

- **Support is 23 of 46 (50%).** The 23 unsupported specializations — 11 mixed
  Hoek-Brown pairings, 6 mixed Mohr-Coulomb pairings, StiffSoil's 2, and the
  remaining crosses — are refused at parse time with a message naming the
  family and this ADR; none can be silently mis-integrated. Widening support
  is new algorithm work (a composite active set for StiffSoil's two surfaces;
  RoundedMohrCoulomb is not even registered in `AllYieldFunctions.h`), not a
  closeout task.
- **Apex tangent has no finite-difference gate, by construction.** The DP, MC
  and HB apex tangents are rank 0 (perfectly plastic) or rank 1 (hardening,
  bulk-only) because the vertex does not move under a perturbation small enough
  to stay in the apex region — there is no FD check that can exercise it
  without leaving the region it is meant to test. Verified instead by the
  closed-form algebra matching the Jacobian-derived tangent to the oracle's own
  floor (P1: DP apex exact 0.0/1.32e-16; P2: MC apex exact zero for
  non-hardening; P3: HB apex exact zero) — a different, but not a weaker, form
  of evidence.
- **Wall-clock numbers in §3 were measured on a shared dev machine** (per
  `ladruno-concurrent-worktrees`) with other background activity at various
  points; the iteration counts (`ops.testIter()` sums), which are
  load-independent, are the primary evidence and are robust — the wall-clock
  ordering (`CP_Algorithmic_Krylov` fastest, `Backward_Euler` slowest) held at
  every family/scale measured, but the exact seconds are not laboratory-clean.
- **The refused-pairing set relies on compile-time family markers staying in
  sync with the runtime generator.** A future new YF/PF specialization that
  reuses an existing family marker incorrectly would be silently accepted or
  refused for the wrong reason; the `static_assert` pre-flight (§ "found while
  implementing" entries across P2/P3) is the safety net, not a runtime check.
- **All seven Zone-A (Ubuntu) runs are green as of this writing** (§1). P4's
  run (#829, 34188231854) was in progress when the manifest table was first
  drafted and resolved `success` on re-check; no branch in the stack has an
  outstanding required-check failure at closeout time.

## §6 Merge guide

Order (base chain, oldest first): **#817 → #819 → #824 → #825 → #827 → #829 →
#826 → this PR (P7, `wp/97h-closeout`)**. #826 is based on #824 (`wp/97c`), not
on #827/#829 — merge it after #829 in this order, retargeting its base to
`ladruno` at that point, same as every other PR in the chain.

For each PR, in order:

1. `gh pr edit <n> --base ladruno` once every PR *below* it in the stack has
   already been merged (a PR's base can only safely become `ladruno` after its
   own parent branch is gone).
2. Re-dispatch Zone-A on the retargeted branch:
   `gh workflow run ladruno.yml --ref <branch>` — a base retarget does not
   automatically re-run CI, and `gh pr checks` after a `ready` flip can lie
   (`ladruno-pr-ready-flip-no-zone-a`) — dispatch and watch the *run id*, not
   the PR checks tab.
3. **Before merging into `ladruno` locally**, if GitHub reports the PR as
   "dirty" or shows an unexpected diff, merge `origin/ladruno` into the branch
   locally first and re-push — the union-merge trap (GitHub's own merge
   preview can disagree with a local three-way merge when the branch has
   moved).
4. Flip to ready (`gh pr ready <n>`) only once Zone-A (Ubuntu) is green on the
   retargeted branch's own run.
5. **Merge commits, never squash, never `--auto`.** The owner merges; agents do
   not (ADR-87 D9/D10).
6. Repeat for the next PR in the stack.

After #826 merges, do the same retarget/dispatch/merge sequence for this
closeout PR itself, based at that point directly on `ladruno`.

The default-flip PR (§3) is explicitly **not** in this chain — it is a new,
separate work package the owner opens after reviewing §3, scoped to the 23
supported specializations.

## §7 Orchestration lessons (for `LEDGER_quirks` and future multi-WP ADRs)

1. **A subagent that backgrounds its own build is never re-woken by the build
   finishing.** Three sessions across this ADR parked on a build or a
   measurement sweep running in the background and were only revived by an
   unrelated notification landing or a later finisher session picking up the
   worktree. Watch with a bounded foreground `until`-loop on the artifact's own
   mtime (`dist/bin/opensees.pyd`), never a detached background watcher with no
   one polling it.
2. **`SendMessage` between concurrent agent sessions was unavailable** during
   this ADR's execution window, which is why (1) above could not be worked
   around by a wake-up call — the finisher pattern (a later session resumes a
   stalled worktree by checking its state, not by being messaged) is the
   fallback, and should be assumed as the default, not the exception, when
   planning a multi-WP ADR with background builds.
3. **A duplicate P6 orchestrator ran concurrently** with the one whose report
   is cited here — two sessions independently began the mesh-scale measurement
   lane. Per `ladruno-duplicate-lane-work`, check `git worktree list` and
   `gh pr list` for an existing branch/PR on the same phase before starting a
   new one; a stray uncommitted duplicate should be rescued as its own commit,
   not discarded, in case it measured something the surviving lane did not.
4. **`pytest --ignore=<bare filename>` inside a glob-expanded invocation does
   not reliably exclude that file.** `pytest test_adr94*.py
   --ignore=test_adr94_matrix.py` still collected and executed
   `test_adr94_matrix.py` in this ADR's P5 session, silently overwriting a
   hand-annotated tracked doc (`_adr94_matrix.md`) with the sweep's own
   regenerated version — caught by `git status` before commit, recovered with
   `git checkout --`. Build the file list explicitly instead
   (`ls test_adr94*.py | grep -v matrix`) whenever a glob and an `--ignore`
   need to disagree in the same invocation; run `test_adr94_matrix.py` only via
   `--collect-only` (safe — fixtures are lazy) or standalone.
5. **`HB_sigma_ci` vs `HB_sigci` produced a vacuous refusal test, twice.** Both
   P1's and P2's Hoek-Brown refusal assertions used the non-existent parameter
   name `HB_sigma_ci` (the real one is `HB_sigci`); under the ADR-94 contract a
   missing parameter is itself a construction error, so the refusal "passed"
   for the wrong reason and never exercised the family gate at all. The general
   fix, applied from P3 onward: every refusal test must also be checked in the
   *positive* direction (the identical deck must construct under an integrator
   that does support it) — a test that only ever refuses cannot distinguish
   "refused because unsupported" from "refused because the deck itself is
   broken."
6. **`ops.ladrunoBuild()` cannot see an uncommitted mutation.** During mutation
   testing, the provenance stamp still reports the last *committed* hash even
   when the binary was built from a scratch, never-committed source edit —
   verifying `ladrunoBuild() == HEAD` after a mutation-gate build proves
   nothing about which source the binary contains; the only reliable check is
   the build log's timestamp plus `git status`/`git diff` on the working tree
   at build time, recorded in the mutation report, not re-derived from the
   binary after the fact.

Items 4–6 are also filed in [[LEDGER_quirks]] under their respective P2/P3/P4/P5
sections (cross-referenced, not duplicated in full there).

## §8 Test totals

Full explicit-list battery, run once on `7e93e4381` and re-confirmed identical
on the post-banner build `33c670bf5`: `test_adr84_p2a_strict_convergence.py`,
`test_adr84_p3_confined_corner.py`, `test_adr94_components.py`,
`test_adr94_contract.py`, `test_adr94_hlist_hb.py`,
`test_adr94_hlist_mechanical.py`, `test_adr94_hlist_numerics.py`,
`test_adr94_redblue_blue.py`, `test_adr94_redblue_cpp.py`,
`test_adr94_redblue_numerics.py`, `test_adr94a_fail_loud.py`,
`test_adr94b_statics.py`, `test_adr94c_numerics.py`, `test_asdplastic_mctc.py`,
`test_asdplastic_response_tags.py`, `test_adr97_p1_smooth.py`,
`test_adr97_p2_principal.py`, `test_adr97_p3_hoekbrown.py`,
`test_adr97_p4_inertness.py`, `test_adr97_p4_numalg.py`,
`test_adr97_p5_stiffsoil.py`, `test_adr97_p6_failloud.py`,
`test_adr97_p6_measure_smoke.py` (23 files; `test_adr94_matrix.py` deliberately
excluded from execution — its module-scoped fixture rewrites a hand-annotated
tracked doc, §7 item 4; `--collect-only` confirms 1 test collected, unrun):

**237 passed, 2 skipped, 0 failed** (20.35 s). The 2 skips are the pre-existing
ADR-94 R3a time-boxed skips (unrelated to this ADR). `git status` confirmed
clean before and after the run — `_adr94_matrix.md` untouched.

Progression across the stack, same explicit-list shape, growing as files were
added: P1 134/2 skip (18 files) → P2 172/2 (19) → P3 204/2 (19) → P5 230/2 (22,
D5+M9 tests added) → P4 236/2/0 (22, final on `7e93e4381`) → **P7 237/2/0 (23,
P6's smoke test added)**.

## See also

[[97_ladruno_asdp_closest_point_adr]] · [[reviews/adr94_verdict]] (§1 M2/M3/M9,
§6 D2, §7 — the finding this ADR closes) · [[reviews/adr97_p1_report]] ·
[[reviews/adr97_p2_report]] · [[reviews/adr97_p3_report]] ·
[[reviews/adr97_p4_report]] · [[reviews/adr97_p5_report]] ·
[[reviews/adr97_p6_measurement]] · [[LEDGER_vanilla_files]] ·
[[LEDGER_implementations]] · [[LEDGER_quirks]]
