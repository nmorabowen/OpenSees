---
wp: ADR-94
title: "A bit-identity gate cannot certify a fix that removes shared static state (ADR-94 wp/94b)"
legacy_seq: 401
---
### A bit-identity gate cannot certify a fix that removes shared static state (ADR-94 wp/94b)
- **Bites:** wp/94b's gate compared committed-stress histories of 23 single-element decks dumped IN ONE PROCESS on the pre-fix build (`229842f7f`) against the per-instance build (`11e3a1283`): 14/23 decks deviated, mostly 1e-9..1e-16 absolute, and the Hoek-Brown deck by 0.12 kPa at step 1 (2e-5 relative). Neither is a wrong answer: (a) the old build's statics carried state between decks run sequentially in one process, so the BASELINE was the contaminated side; (b) the fix changes the assembled tangent, hence the global Newton path, so converged stresses differ at the global-tolerance level — on rock with E ~ 1e7 kPa a `NormDispIncr 1e-8` tolerance is Δσ ≈ E·1e-8 ≈ 0.1 kPa, which is what was measured.
- **Why:** "byte-identical before/after" presumes the change is inert on the converged path; a tangent fix is not inert on the PATH, only on the RESULT, and a static-state fix is not even inert on the baseline.
- **Rule:** gate a static-state or tangent fix on (1) per-deck runs in FRESH subprocesses on both builds, and (2) converged-result agreement within the global tolerance (the R1 blue two-cube test does this at 1e-6), plus (3) iteration counts (the two-cube model went from +62 % to −12 % versus the separate sum). Keep bit-identity gates for changes that do not touch the tangent. Also: `tests/test_adr94_matrix.py` regenerates the tracked `_adr94_matrix.md` on every run — revert it before committing unless the regeneration is intended.
