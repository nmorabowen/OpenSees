---
wp: WP-215
title: "A relative residual ||Kx-b||/||b|| refuses good solves under Penalty; use the backward error"
date: 2026-10-26
pr: "#1027"
---
### Under Penalty a backward-stable K solve shows ‖Kx − b‖ / ‖b‖ ≈ 1e-4

- **Bites:** any acceptance check of a linear solve written as
  ‖Kx − b‖ / ‖b‖ < tol. WP-213's `LadrunoRitzKernel` used tol = 1e-4: on the
  tied-floor Penalty model (α = 1e14, kip-in) `eigen -fullGenLapack 1` then
  `residualModes -dir 1` was refused with "did not solve K x = f (relative
  residual 0.00016)" on UmfPack; with 2-3 modes, or after a re-run of ARPACK,
  it passed (#1025).
- **Why:** a backward-stable solve leaves a residual of order ε‖K‖‖x‖. Under
  Penalty ‖K‖ ≈ α while ‖b‖ ≈ k‖x‖, so the ratio is ≈ ε·α/k ≈ 1e-4: round-off,
  not a failed solve. Which right-hand side the check sees (it depends on the
  mode count) decides the outcome.
- **Fix/status:** measure the normwise backward error
  ‖Kx − b‖ / (‖ |K||x| ‖ + ‖b‖), with |K||x| from the same element sweep as
  Kx (`LadrunoRitzKernel::applyK` `absKX`). It is never above the old ratio,
  sits at ~ε for a good solve and at O(0.1) for `system Diagonal`, which is
  still refused. Same family as the Rayleigh-quotient floor
  (WP-213 penalty-eigenvalue-roundoff). Learned 2026-10-26 (WP-215, #1025).
