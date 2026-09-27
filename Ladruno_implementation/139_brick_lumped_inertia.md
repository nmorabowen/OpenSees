# WP-139 — LadrunoBrick `-lumped`: one mass model for the residual AND the tangent

Revision 0 (scoping). Opened 2026-09-27 by the owner from WP-124 gap **C8**
([[124_continuum_shell_helpers]]). Branch `wp/139-brick-lumped-inertia`, cut from `ladruno` @ `64a0341a6`
(WP-124 merged). Status: **scoping; draft PR.** No production code changed yet.

## Problem

Under `element LadrunoBrick … -lumped` (`massType == 1`) the element uses TWO different mass matrices:

| Consumer | Mass it uses | Code |
|---|---|---|
| inertia RESIDUAL (`getResistingForceIncInertia` → `formInertiaTerms(0)`) | **consistent** `ρ Σ_gp N_j N_k dV` (always, whatever `massType`) | `LadrunoBrick.cpp` `formInertiaTerms`, the `resid += temp * momentum` loop |
| mass TANGENT (`getMass` → `formInertiaTerms(1)`), αM Rayleigh, ground-motion load (`addInertiaLoadToUnbalance`) | **row-sum lumped** `diag(ρ Σ_gp N_j dV)` | same loop, `massType == 1` branch |

So implicit dynamics integrates a hybrid: consistent inertia with a lumped Jacobian, lumped αM damping and a
lumped ground load. Newton's Jacobian is not the derivative of the residual, so it converges only LINEARLY.
Observed in WP-124: an ELASTIC `-lumped` Newmark step does not reach `NormDispIncr 1e-10` in 40 iterations
(the WP-124 fingerprint's `Brick/lumped` case needs `1e-7 / 300`).

**Where it comes from.** Upstream `SRC/element/brick/Brick.cpp:882-890` has the same code (lumped only in the
tangent branch); LadrunoBrick inherited it. LadrunoBrick's "documented as intentional" note (the ADR-77 mass-cache
comment) records the fact, it does not justify the physics. **LadrunoBrick20 already fixed this pattern (F-1):**
under `-lumped` its residual inertia uses the same cached lumped diagonal as the tangent
(`LadrunoBrick20::formInertiaResidual`, `resid(c) += M0(c,c) a(c)`), so the two cannot disagree.

**Where it does NOT bite.** Explicit `CentralDifferenceLadruno` forms the residual at trial acceleration 0 (the
ADR-68 T7 skip returns before the residual inertia loop), so only the lumped M is ever used — the `-lumped`
flag's main purpose (explicit Δt) is consistent today. To be verified per explicit integrator (vanilla
`CentralDifference`, `ExplicitBathe`, SMS) in P0.

## Options

| | What | Consequence |
|---|---|---|
| **A (recommended)** | Brick20's F-1: under `-lumped` the residual inertia is `M_L a` with the SAME lumped diagonal as `getMass` | `-lumped` means one lumped mass model everywhere; Newton quadratic again; implicit `-lumped` results CHANGE (they become proper lumped-mass dynamics); explicit runs expected bit-identical (P0 checks) |
| B | Keep the consistent residual, make the tangent consistent too | `-lumped` would no longer lump the explicit Δt mass → defeats the flag |
| C | Keep the hybrid, document it | Status quo: linearly convergent implicit `-lumped`, a mass model no textbook describes |

## Plan

- **P0 — measure** (no code): the WP-124 fingerprint suite + an explicit lane per integrator on `Brick/lumped`
  (CDL, vanilla CentralDifference, ExplicitBathe, SMS) to record which paths read the residual inertia at a
  nonzero trial acceleration. Newton iteration counts for an elastic `-lumped` Newmark step (baseline).
- **P1 — fix (option A)**: under `massType == 1`, `formInertiaTerms(0)` adds `M_L(c,c) a(c)` from the same
  `getMass()` diagonal (the per-instance `LadrunoMassCache`), in the style of Brick20 F-1. Consistent mass
  (`massType == 0`) untouched: its path must stay bit-identical.
- **Proof**: fingerprint — ONLY `Brick/lumped` implicit series may change; every other variant, and every
  explicit lane, byte-identical. Absolute oracles: (i) elastic `-lumped` Newmark step converges in ≤ 3 Newton
  iterations (quadratic); (ii) free rigid body under a nodal force: `a = F / m` exactly on every node (lumped);
  (iii) lumped SDOF-like patch: period from `M_L`, not `M_c`. Mutation row: revert to the consistent residual
  → (i) fails. Upstream `Brick` left alone (vanilla; out of scope unless the owner asks).
- Ledgers: LEDGER_implementations row, LEDGER_quirks entry (the hybrid and why it hid: Windows/explicit runs
  never see it), the `ladruno-new-element` guide item "one mass model for residual AND tangent".

## Open questions (owner)

- Confirm option A. It changes implicit `-lumped` LadrunoBrick results by design.
- Upstream `Brick` carries the same hybrid: leave it (vanilla), or fix it too as a separate vanilla edit?
