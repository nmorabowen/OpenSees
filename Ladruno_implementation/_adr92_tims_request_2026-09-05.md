---
title: "ADR 91 — IMPL-EX integration for LadrunoSANISAND: a solver that cannot fail at the corner, inside the coupled element"
project: Ladruno
type: ADR
status: PROPOSED (request from the TIMs response-curve-matrix act, 2026-09-05); P0 not started
priority: high
owner: nmora
requested_by: TIMs Workbench, act `work/ape/response-curve-matrix` (gupi claude-code)
related:
  - "[[86_ladruno_sanisand_adr]]"
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[31_ladruno_concrete3d_adr]]"
  - "[[71_ladruno_up_adr]]"
  - "[[84_ladruno_mc_tension_cutoff_adr]]"
  - "[[LEDGER_implementations]]"
tags: [adr, material, sanisand, implex, integration, footing, tims]
updated: 2026-09-05
---
> [!note] PRESERVED AS RECEIVED — do not edit
> This is the TIMs request of 2026-09-05, copied verbatim from the untracked
> `Ladruno_implementation/91_ladruno_sanisand_implex_adr.md` in the main checkout at
> `3f003d110`. It is kept unedited so the fork's decisions can be read against what was
> actually asked for. **It contains six statements about this code base that the fork-side
> review corrected** — see `_adr92_sanisand_implex_scoping.md` §3. The ADR that supersedes
> it is `92_ladruno_sanisand_implex_adr.md`. Renumbered 91 -> 92 (ADR 91 is the shell
> stiffness modifiers).


# ADR 91 — IMPL-EX for `LadrunoSANISAND`

Numbered 91 because 89 is not on disk and may be reserved; renumber if 89 is free.

## 1. What / Why

The TIMs response-curve-matrix act has spent five days (2026-09-01 to 05) pushing the Gorini footing (1.5 m square, 12 × 6 × 7.5 m half-domain, `LadrunoUP` H8 bbar, coarse mesh 21 058 DOF) to failure on `ManzariDafalias` and on `LadrunoSANISAND -Presidual 0.0 -Pmin 0.0101`, drained and undrained, on Esmeralda. The record is `Tries/response-curve-matrix/model/harness/ESMERALDA.md` §§1–33 and the session notes of 2026-09-04/05 in that act. What it established about the solver, in the order it was learned:

1. **Every drained leg ends in a convergence wall, not a limit point.** Thirty production legs; no leg met the WP1 criterion (tail tangent ≤ 2 % of the first-step tangent over Δ(s/B) ≥ 0.01). Two read as limit points met by walls by a tail-tangent rule; the rest are walls with the load still rising.
2. **Part of the wall was the harness's ladder.** With the third rung (KrylovNewton at ten times the tolerance) removed and a fixed 1e-4 m increment, bare `LadrunoSANISAND` reaches s/B 0.0206 where the full ladder walled at 0.0068 (jobs 146299 vs 146175, ESMERALDA §§30–31). The relaxed rung commits states converged to 1e-4 that nothing afterwards converges from.
3. **The remainder of the wall has an address.** With a clean solver history the bare column walls at s/B 0.0206 with Newton and line search failing at 1e-4 m and at every halving to 5e-5. The plates at the last checkpoint (S8f) put the failing skin at and diagonally beyond the footing's geometric corner on both flanks, with the minimum mean stress pinned at the material's own floor 0.0101 kPa; the fine mesh (175 290 DOF) resolves the same thing as a stress concentration two element rows thick at the corner with p′ < 0, sharpened by refinement, and fails there as a load drop at 54 kPa. A cohesionless material at zero confinement has zero strength and, for SANISAND with G ∝ √p′, zero stiffness: the consistent tangent is singular in that ring and Newton has no descent direction. Every regularisation that gives the ring strength or confinement (10 kPa surcharge, a 2 kPa Drucker–Prager crust in the top row, vanilla's 1.01 kPa residual pressure scoped to the top row) carries the leg past it. The unregularised problem is the one the act must also be able to solve, because the regularisations each add a capacity term.
4. **The matrix-free routes were scoped and one failed for a material reason.** `LadrunoDynamicRelaxation` (Gershgorin mass, kinetic damping) relaxes every 1 mm hold to a residual below 1e-6 kN on the ndf-3 `LadrunoBrick` twin and yet produces a convex curve at twice Newton's load by s/B 0.02 (job 146300, ESMERALDA §33). On a perfectly plastic cone the same scheme agrees with Newton to 3 %. Measured on the desktop patch with SANISAND itself (`Session Notes/2026-09-05 — Relaxation and Memory.md`, 2026-09-05): every relaxation setting gives a different "equilibrium", from 0.10 to 3.08 times Newton's load depending on the damping law, the mass safety and the hold size; at a platen held at 3.000 mm the back-stress norm at the footing-edge Gauss point runs a full cycle (0.646 → 0.536 → 0.631) and the fabric norm grows by 35 % with nothing moving; by s/B 0.01 relaxation has densified the sand two to seven times more than Newton (e = 0.6774–0.6874 against Newton's 0.6920 and the explicit transient's 0.6906). The fork commits every sweep (`DirectIntegrationAnalysis.cpp:259`) and offers no staged commit, so the cure would be C++ in the integrator, not a flag. The explicit central-difference transient with ρ×100 mass scaling lands at 1.17/1.22/0.93 times Newton at s/B 0.005/0.01/0.015 and rings by a quarter. `CentralDifferenceLadruno` avoids that by construction (monotonic path in time) at 90 000 steps to s/B 0.10 on the coarse mesh, single-threaded, and neither scheme can run on `LadrunoUP`, whose pressure rows carry no mass (ADR 71 §6).

IMPL-EX (Oliver, Huespe & Cante 2008) is the one candidate that addresses 3 and 4 at once and keeps the coupled element: the internal variables at step n+1 are extrapolated from the converged step n, the stress handed to the element is computed from that extrapolated state, and the algorithmic tangent is the elastic operator, symmetric and positive definite at every Gauss point regardless of confinement. The global system is linear in the step. There is no fictitious mass, no oscillation, no memory artefact; the load path is the one the footing imposes. The fork already carries an IMPL-EX implementation with error control in `ASDConcrete3DMaterial` (`-implex`, `-implexControl`, `-implexAlpha`, `implex_error_tolerance`, `implex_time_redution_limit`, recorder `implexError`), written for softening concrete for the same reason, and `LadrunoConcrete3D` inherits it. This ADR asks for the same on `LadrunoSANISAND`.

## 2. Where it goes in the material

`LadrunoSANISAND` subclasses `ManzariDafalias` (ADR 86). The base carries four local integrators selected by `mScheme` (`ManzariDafalias.h:156`): 0 forward Euler, 1 backward Euler implicit (`BackwardEuler_CPPM`), 2 backward Euler with stability considerations, 3 constrained-strain-increment explicit; the OpenSees default for the 3D class is scheme 1 with tangent type 0 and the modified-Euler substepping (`ModifiedEuler`) behind scheme 0/3. IMPL-EX is meaningful only on top of an implicit local return (schemes 1 and 2): the committed state at step n comes from the implicit return, and the extrapolation acts on the increment of the plastic multiplier.

The plastic multiplier in ManzariDafalias is the scalar `dGamma` (`mDGamma`, position [25] of the recorded `gp_state`). The internal variables it drives are the back-stress `mAlpha`, the fabric tensor `mFabric`, the void ratio `mVoidRatio` (kinematic, from the strain), `mAlpha_in` (the last reversal). The IMPL-EX step:

- Extrapolate: `dGamma_tilde(n+1) = dGamma(n) · Δt(n+1)/Δt(n)` (with the `implex_alpha` scale of the concrete implementation; under displacement control Δt is the pseudo-time of the step and the ratio is the ratio of push increments, 1 at a fixed increment).
- Stress: elastic predictor from the strain increment, minus the plastic correction evaluated with `dGamma_tilde` and the flow direction and dilatancy of the **committed** state n (`n`, `R`, `D` at σ_n, α_n, z_n, e_n), so that the update is explicit in the state and implicit in nothing.
- Tangent: the elastic operator at the committed state (G(p_n), K(p_n)). Constant within the step, symmetric, positive definite wherever p_n ≥ p_min.
- Commit: at commit, run the implicit return (`BackwardEuler_CPPM`) from the committed state n with the actual strain increment to obtain the true state n+1 and the true `dGamma(n+1)` for the next extrapolation; store `implex_error` as the norm of the difference between the extrapolated stress and the implicit stress, as the concrete class does.
- Control (`-implexControl`): if `implex_error` exceeds `implex_error_tolerance`, request a step reduction through the same mechanism the concrete class uses (the element/analysis reads the error and the integrator's step controller halves), down to `implex_time_reduction_limit`.

Reversals: `mAlpha_in` updates at load reversal on the committed path only; the extrapolated step must not trigger a reversal update. The fabric tensor evolves with the dilatancy sign on the committed path only.

Low confinement: at Gauss points where p_n is at the floor `p_min`, the elastic operator is the floor's operator (small, positive) and the extrapolated plastic correction scales with `dGamma(n)`, which is bounded because the committed return bounded it; the extrapolated stress cannot leave the floor by more than the elastic predictor allows in one increment. This is the property that removes the corner as a solver event. It does not remove it as a mechanics event: the ring still flows, and the load that the flow carries is what the regularisation study prices.

## 3. Decision requested

Add to `LadrunoSANISAND` (3D and plane-strain) the options

```
nDMaterial LadrunoSANISAND tag ... -Presidual pr -Pmin pmin -honorTolR h \
    -implex [-implexControl tol reductionLimit] [-implexAlpha a]
```

with the semantics of `ASDConcrete3DMaterial`: `-implex` on selects the IMPL-EX stress and tangent for `setTrialStrain`/`getStress`/`getTangent`, with the implicit return at `commitState`; `-implexControl` enables the error-driven step request; `-implexAlpha` scales the extrapolation (1.0 = Oliver's). Recorders `implexError` and `avgImplexError` as in the concrete class. `revertToLastCommit` restores the committed state and the stored `dGamma(n)`. `getCopy` carries the flags (the FSPM lesson of ADR 86 / PR #779). Refuse `-implex` on `mScheme` 0 and 3 with a sentence.

Scope: `LadrunoSANISAND` only, not vanilla `ManzariDafalias`; the base class is not touched except where a virtual hook is needed to reach `dGamma` and the committed flow direction, which ADR 86's tripwire seam (`ManzariDafalias.h:161`) shows how to add without altering vanilla behaviour.

## 4. Usage, in the act's harness

`drive.py --cell D-L --integrator newton --ladder strict --ds 1e-4 --ds-max 1e-4 --implex [--implex-control 0.05 0.01]`: the harness's S0 material line gains `-implex`; S5 keeps the strict ladder (Newton at the declared tolerance converges in one iteration on a linear step and the ladder never fires); the checkpoint dump adds `implex_error` per Gauss point; `run.json` records the flags under `declared.integrator`.

## 5. Acceptance / tests (`tests/test_ladruno_sanisand_implex.py`)

1. **Material point, the Level 0 anchors.** The act's drained triaxial instrument (`Tries/response-curve-matrix/model/probe_fspm_sanisand.tcl` and its `LadrunoSANISAND` branch, Level 0 of the act) at the declared strain increment: `-implex` against the implicit return within 1 % in η_end and ε_v at the anchor increment, and the difference halving when the increment halves (first-order convergence, the IMPL-EX signature). At the G1 increment the two must agree to the G1 tolerance.
2. **Reversal.** A cyclic triaxial leg: `mAlpha_in` updates only at committed reversals; the extrapolated stress never triggers one; the implex error spikes at the reversal step and decays.
3. **Corner patch.** `Tries/response-curve-matrix/model/harness/probe_mc_tangent.py --probe patch` rebuilt with `LadrunoSANISAND` (the scoping note `2026-09-05 — Scoping the Explicit Route.md` has the patch and the harness's step controller): with `-implex` the push reaches s/B 0.05 at a fixed 1e-4 m with the implex error below tolerance, where Newton on the implicit material walls at 0.011 under the controller and passes at a fixed 1e-4; the curve within 3 % of Newton's fixed-step curve where both exist.
4. **The coarse bare leg.** `D-L` coarse, strict ladder, 1e-4, `-implex`: passes s/B 0.0206 (job 146299's wall) and reaches s/B 0.10 or a WP1 limit point; the curve to 0.0206 within 3 % of 146299's; the implex error field at the corner reported.
5. **The coupled row.** `U-L` coarse (transient, k = 1e-10 m/s, Newmark 0.6/0.3025, the S7b increment 3.125e-6 m): `-implex` runs inside `LadrunoUP` without change to the element; the pore pressure under the footprint compared with the implicit leg at matched settlement.
6. **`getCopy`, `sendSelf`/`recvSelf`** carry the flags and the stored `dGamma(n)`; the FSPM wrapper over `LadrunoSANISAND -implex` constructs (the ADR 86 acceptance harness).
7. **Vanilla untouched**: `ManzariDafalias` decks bit-identical before and after (the ADR 86 bit-identity gate, `LEDGER_vanilla_files`).

## 6. Phasing

- **P0** — the flag on `LadrunoSANISAND3D`, extrapolated stress and elastic tangent, implicit commit, `implexError` recorder, tests 1, 3, 6, 7. Build by the human; desktop smoke by the act.
- **P1** — `-implexControl` with the step request, test 2, and the plane-strain class (the 2D exploration act, `work/ape/2d-exploration`, T6, needs it on the strip).
- **P2** — tests 4 and 5 on Esmeralda under the act's harness; the reading of the corner's implex error field against the strict leg's plates.
- **Not in scope**: IMPL-EX on vanilla `ManzariDafalias`; any change to `LadrunoUP`; the rate regularisation of ADR 90, which addresses localisation width and not the solver, and which remains the alternative if the flow at the corner turns out to need a length scale rather than a solvable step.

## 7. Risks

- **First-order step error.** The extrapolated path lags the implicit one by O(Δt); at the act's 1e-4 m increment the material-point test (1) says how much. If the error at the corner is dominated by the ring's own flow, `-implexControl` will halve the step there and the cost rises where the wall was; that is the correct behaviour and it is bounded by the reduction limit.
- **Overshoot at the floor.** Where p_n is at `p_min` and `dGamma(n)` is large from the previous step's flow, the extrapolated correction can push the trial stress to the floor in one increment; the elastic tangent at the floor is small and the global step remains solvable, but the implex error there will be large. Test 3 measures it; the fabric and reversal guards of §2 keep it from contaminating the committed state.
- **Two tangents in one mesh.** Under `-implex` the SANISAND elements return an elastic tangent while a crust or a footing returns its own; the global matrix stays positive definite as long as every material's tangent is. The crust's Drucker–Prager consistent tangent is not always so; the crust leg under `-implex` should use the crust's elastic tangent too, or the crust should be run implicit and the sand implex with the ladder kept at rung 0 and 1.
- **Reading.** An IMPL-EX curve is not a converged-equilibrium curve in the implicit sense; every row satisfies equilibrium with the extrapolated stress. The act must report `implexError` beside every WP1 verdict, and a limit point called on an IMPL-EX leg must be confirmed by the implicit material up to the last settlement the implicit solver reaches.

## 8. Ledger / obligations

- `LEDGER_implementations`: one row, `LadrunoSANISAND -implex`, ADR 91, files `LadrunoSANISAND.{h,cpp}`, `LadrunoSANISAND3D.{h,cpp}`, `LadrunoSANISANDPlaneStrain.{h,cpp}`, a virtual seam in `ManzariDafalias.h` if needed (vanilla behaviour unchanged, the bit-identity gate).
- The act records the build hash on first use in `APE/information/environment.md` (the human's, to `main`) and every IMPL-EX leg names `-implex` and its control values in `run.json`.
- ADR 90's D2 (the reopened non-associated softening width question) is not closed by this ADR and is not touched by it.

## 9. Staffing

P0: Opus, high, in the fork, with the concrete class as the template and ADR 86's seam as the pattern; the author does not write the tests. Tests: Sonnet from the acceptance list above, with the act's Level 0 instrument as the oracle. Review: Fable, once, on the reversal and floor guards of §2 and on test 3's reading. Build: the human.

## Log

- 2026-09-05 — Requested by the act after the strict-ladder diagnostic (ESMERALDA §31) located the remaining wall at the footing's corner and the relaxation route was found to load SANISAND's memory (§33). Written on `ladruno` at `3f003d110`, uncommitted.
