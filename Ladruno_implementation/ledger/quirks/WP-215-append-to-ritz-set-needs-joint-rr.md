---
wp: WP-215
title: "Vectors M-orthogonal to a Ritz set are not K-orthogonal to it; appending needs a joint Rayleigh-Ritz"
date: 2026-10-26
pr: "#1027"
---
### Appending residual vectors to Ritz vectors couples them through K

- **Bites:** extending a mode set by deflating new vectors against it in the
  M-inner product and running Rayleigh–Ritz on the new vectors only (what
  `residualModes` does after `eigen`). After `ritzVectors` the set is not an
  eigenbasis: V3, `ritzVectors 3 -dir 1` then `residualModes -dir 1` gave
  ΨᵀKR up to −442 against λ_r = 355, 634, the CQC base shear 1.137 × instead
  of 1.004 × the reference, and the rigid static response 42 % off. Every
  number looked plausible.
- **Why:** M-orthogonality gives K-orthogonality only for eigenvectors
  (φᵀKr = λφᵀMr = 0); for a Ritz vector ψ, Kψ ≠ λMψ, so K⁻¹MΨ leaves span Ψ.
- **Fix/status:** `residualModes` detects a kept RITZ set and re-solves
  span[Ψ, R] by one joint Rayleigh–Ritz step, re-installing all vectors as
  RITZ (test asserts ΨᵀKΨ diagonal to 1e-8). Same rule for any future
  "append to a non-eigen basis". The WP-213 oracle had mirrored the kernel and
  could not see it: oracles now build the subspace independently. Learned
  2026-10-26 (WP-215, review of #1027).
