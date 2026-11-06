---
wp: ADR-84
title: "The same Backward_Euler calls a step \"elastic\" whenever f merely DECREASED — perpetuating an existing violation (FIXED opt-in, ADR-84 P2a)"
legacy_seq: 305
---
### The same `Backward_Euler` calls a step "elastic" whenever f merely DECREASED — perpetuating an existing violation (FIXED opt-in, ADR-84 P2a)

- **Bites:** the elastic early-exit reads
  `if ((yf_val_start <= 0 && yf_val_end <= 0) || (yf_val_start - yf_val_end > tol_yf)) { Stiffness = Eelastic; return 0; }`.
  The second disjunct accepts the step as elastic on the sole grounds that `f`
  **went down by more than tol** — with no requirement that the end state be
  admissible. From an admissible commit (`f_start <= 0`) it cannot manufacture a
  violation, so it is invisible in any clean-history test; but from an
  already-inadmissible commit — exactly what the exhaustion-accept above
  produces — it will happily carry `f_end > 0` forward step after step as long
  as the value keeps shrinking, never engaging the plastic corrector that would
  actually pull the point back onto the surface. The two defects compound: one
  creates the inadmissible state, the other protects it from correction.
- **Why:** the disjunct exists for elastic UNLOADING from a plastic state, where
  `f_start ≈ 0` and `f_end < 0` — a legitimate case that the first disjunct's
  `yf_val_start <= 0.0` misses on the boundary. But "f decreased" is a much
  weaker test than "the end state is admissible", and the weaker test is what
  got written.
- **Workaround/status (2026-08-13, ADR-84 P2a, PR):** under the same
  `strict_convergence 1` flag the second disjunct additionally requires
  `yf_val_end <= tol_yf`; when it does not hold the step falls through to the
  plastic corrector instead of being accepted. The unloading case is untouched
  (`f_end < 0` satisfies the added condition trivially). Default 0 keeps the
  upstream disjunct exactly as written. Note the ordering: the MCTC
  `special_return` hook (ADR-84 P0) sits BELOW this exit, so a step wrongly
  classified as elastic here never reaches the hook — fixing the classification
  is what lets the hook see the states it was built for.
