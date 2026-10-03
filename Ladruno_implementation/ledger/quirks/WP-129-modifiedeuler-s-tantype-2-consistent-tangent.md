---
wp: WP-129
title: "ModifiedEuler's TanType 2 \"consistent\" tangent chain accumulates T where the recurrence needs dT (WP-129, code-read)"
legacy_seq: 490
---
### `ModifiedEuler`'s `TanType 2` "consistent" tangent chain accumulates `T` where the recurrence needs `dT` (WP-129, code-read)
- **Bites:** after each accepted substep, `T += dT; aCep_Consistent = ½(aCep1 + aCep2)·(aD·aCep_Consistent + T·I)`. For purely elastic stages (`aCep = Ce`) this gives `Σ T_k·Ce`, i.e. `(N+1)/2 · Ce` for N equal substeps, where the derivative of the update is `Ce` (the recurrence with `dT` gives exactly `Ce`). A drift correction resets it (`Stress_Correction` sets `aCep_Consistent = aCep` when it acts), so the delivered TanType-2 operator is a history-dependent mix. Not measured on a deck here; the elastic arithmetic is exact.
- **Workaround/status:** SAS-ME returns ONE continuum tangent at the end state for TanType 1 and 2 (SAS practice). ModifiedEuler unchanged. Candidate for a separate fix (it moves TanType-2 decks).
