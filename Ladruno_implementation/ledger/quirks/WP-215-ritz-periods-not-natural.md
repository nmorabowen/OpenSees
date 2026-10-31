---
wp: WP-215
title: "Ritz periods beyond the first few are not natural periods"
date: 2026-10-26
pr: "#1027"
---
### `ritzVectors` periods above the first few are not natural periods (and a Ritz set is not a frequency-domain basis)

- **Bites:** after `ritzVectors $n`, `modalProperties`, `nodeEigenvector`,
  recorders and `eigen`-style period tables all read the Ritz pairs as if they
  were modes. Reporting T_k for k beyond the first few as "the structure's
  periods", checking period separation with them, or feeding the set to
  `frequencyResponse` / `steadyStateDynamics` / `randomResponse` gives wrong
  answers silently in stock OpenSees terms.
- **Why:** a Ritz value is the eigenvalue of K restricted to the Krylov
  subspace of the starting loads. The lowest ones are upper bounds on the
  lowest eigenvalues and converge fast; the highest ones bound nothing useful
  and their shapes mix many modes. The basis is tuned to the loads (it
  reproduces their static response exactly), not to a frequency band.
- **Fix/status:** `modalProperties -print` / `-file` tags every Ritz row
  `# Z` with a NOTE that says so, `ritzVectors` prints the same caution, and
  the three frequency-domain commands warn on a Ritz set. Use `eigen` for
  natural periods and frequency-domain work; use `ritzVectors` for RSA and
  modal time history. Learned 2026-10-26 (WP-215).
