---
wp: ADR-95
title: "ADR-95: a \"control that must fail\" can be a bug in disguise (2026-09-07)"
date: 2026-09-07
legacy_seq: 396
---
## ADR-95: a "control that must fail" can be a bug in disguise (2026-09-07)

The R3 Prandtl gate's associated-flow control asserted that the ψ = φ leg must NOT produce a
capacity, because it had only ever been observed seizing on the step floor while hardening. That
ending was the vanilla UW `DruckerPrager` dead tension-cutoff/corner branch (ADR-95, PR #803):
dilatant flow reaches I1 = T earlier, so the associated leg hit the defect first. On the repaired
material it plateaus at 1.6026 of the non-associated exact (h0 = 0.5). Rule: a control whose
expected outcome is "the solver fails" must state WHY it fails and be re-checked whenever the
material or solver changes; assert the discrimination you need (here: distinct answers), not the
failure mode you happened to observe.
