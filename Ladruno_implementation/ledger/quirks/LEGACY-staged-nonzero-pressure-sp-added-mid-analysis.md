---
wp: LEGACY
title: "Staged NONZERO pressure sp added mid-analysis under constraints('Transformation') converges cleanly to a WRONG steady state on LadrunoUP (interior p wildly off…"
legacy_seq: 172
---
### Staged NONZERO pressure `sp` added mid-analysis under `constraints('Transformation')` converges cleanly to a WRONG steady state on LadrunoUP (interior p wildly off; Penalty/Lagrange correct)
- **Bites:** the ADR-71 §3.2 initialization recipe — run a stage, then add `sp p=<head>` and continue. Under Transformation the model converges (rc=0) to interior p ≈ −73 for a top head of +1 on a sealed 4-element column (measured); the identical `sp` present from step 1 is handled correctly by all three handlers, and Penalty/Lagrange are correct in both sequences.
- **Why (suspected):** Transformation condenses constrained DOFs at analysis-setup time; a mid-analysis `sp` after `wipeAnalysis` re-setup interacts with the committed-but-unconstrained p state; exact mechanism not chased (P4 revisit alongside the gravity/hydrostatic init recipe).
- **Rule for the guide:** staged prescribed-head sequences use `constraints('Penalty', ...)` (or Lagrange); mirrors the existing fully-prescribed-rig Transformation trap. Repro pinned in `tests/test_ladruno_up_element_analytic.py` (ADR-71 P1, 2026-07-11).
