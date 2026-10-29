---
wp: WP-200
title: "LadrunoForceBeamColumn re-solved every update from the committed state: Newton floor and a stiff predictor tangent"
date: 2026-10-29
---
### LadrunoForceBeamColumn Newton stalls near |R| ~ 1e-2, or the next step's predictor tangent is 20-45 % stiffer than forceBeamColumn's

- **Bites:** (1) plain Newton with `NormDispIncr 1e-8` stalls with a constant
  residual (|R| ≈ 0.014 N on a 5 m RC cantilever) and fails steps that vanilla
  `forceBeamColumn` converges in 3 iterations. (2) On frames whose beams carry
  `-beamUniform` loads, each step's first Newton residual was ~1e6 N against
  ~1e3 N for vanilla. On frame-ladruno-validation slice 04a this was 6-9
  iterations per step against 3-5, and a hard failure at 2 % drift where
  vanilla completes all 141 steps under the same ladder.
- **Why:** `evalBulkState` started every call from the COMMITTED state
  (`Secommit`, `kvbulkCommit`, `vscommit`/`fscommit`/`Ssrcommit`, increment
  `v - vCommitBulk`). Vanilla iterates from the last trial state. (1) Each
  update re-solved the whole step and stopped as soon as `|dW| < tol`, so the
  force carried an error of the size of the element tolerance that no global
  iteration removed. (2) `DisplacementControl::domainChanged` (every `analyze`
  after `ops.integrator`) calls `applyLoadDomain` and `formUnbalance`. That
  sets `eleLoadDirty` on span-loaded members, and `getResistingForce` then
  re-solves them at zero increment. The re-solve fed the sub-tolerance section
  residual (~1e-12 strain) into the fibres. Concrete02/Steel02 read a tiny
  reversed increment as unloading and returned their unloading tangent, and
  the predictor was then formed on that tangent.
- **Not a symptom by itself:** "K after a zero-increment `Linear` step equals
  K at the converged state" does NOT hold in general. That step applies the
  residual correction, and a ~1e-15 strain change of the wrong sign flips a
  fibre on vanilla too (13 % on a 3D RC cantilever). Compare the Newton
  history and the tangent across a load re-application, not that step.
- **Fix/status:** fixed 2026-10-29. `evalBulkState` iterates from the trial
  state and `vTrialBulk` (reset to `vCommitBulk` by `revertToLastCommit` and
  `recvSelf`). `update()` skips the re-solve when the basic deformation and the
  element loads (pointer, factor, data) are unchanged since the last solve.
  2D and 3D. Regression: `tests/test_ladruno_force_beam_column_trial_state.py`.
