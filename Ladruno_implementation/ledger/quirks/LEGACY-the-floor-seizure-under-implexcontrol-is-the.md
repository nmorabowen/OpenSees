---
wp: LEGACY
title: "The FLOOR seizure under -implexControl is the implexPrimed bare > 0.0 sign test: a 1e-12 plastic history forfeits the un-primed exemption and is refused on an…"
legacy_seq: 433
---
### The FLOOR seizure under `-implexControl` is the `implexPrimed` bare `> 0.0` sign test: a 1e-12 plastic history forfeits the un-primed exemption and is refused on an error that does not decay with `dt`
- **Bites:** the run does not spend its subdivision budget — it reaches the step
  FLOOR. Halving stops helping: the same Gauss point refuses at every rung all
  the way down to `DS_MIN`.
- **Why:** `LadrunoSANISAND.cpp:3021` is
  `const bool implexPrimed = (this->GetNorm_Cov(mImplexDEpsP) > 0.0);`. The
  un-primed exemption at `:3023` exists precisely because the companion's
  drift-correction jump does not scale with `d_eps` — the source says so at
  `:3005`-`:3014` ("it ASYMPTOTED at ~0.076 instead of decaying ... the signature
  of a companion jump that is independent of the increment ... a dead analysis").
  But ANY non-zero committed plastic increment, however small, forfeits the
  exemption. A Gauss point that took essentially no plastic strain last step and
  yields in this one is therefore refused on exactly the quantity the exemption
  exists to tolerate.
- **Measured (ADR-92 F10 leg B, `out/refusal_warnings_B.txt`):** refusals at
  `|d_eps_p(n)| =` 9.09e-13, 6.66e-12, 9.58e-21 and 2.40e-21. At 9.09e-13 the
  step halves `4e-5 -> 2e-5` and the error moves **0.2243 -> 0.2143, i.e. 4.5 %**
  — it does not decay. The other 39 of 49 throttled lines
  (`|d_eps_p(n)| >= 5e-8`) DO decay first order (0.299 / 0.149 / 0.084 as `|dt|`
  goes `4e-5 -> 2e-5 -> 1e-5`), so the two families are cleanly separable in the
  warning text itself.
- **Diagnostic that does NOT see it:** a step-size refinement probe of the bulk
  error field at the same settlement reports **zero** Gauss points over tolerance
  at every `ds` from 8e-5 down to 5e-7, first order throughout. The seizure is a
  handful of points reached only along a path that has already been refusing.
  Read the throttled warning lines, not the field.
- **Workaround/status (2026-09-14):** unfixed, and deliberately — ADR-92 F10 is a
  diagnosis. The matching fix is the shape P2-5b already used for the
  loading-reversal reset: a RELATIVE test (`||d_eps_p(n)|| > rel * ||d_eps||`)
  instead of an absolute/sign one. It needs its own gate and mutation score.
  Until then: an `-implexControl` refusal whose warning line reports a
  `|d_eps_p(n)|` many decades below the strain increment is a first-yield event,
  not an extrapolation failure, and no amount of subdivision will clear it —
  turn the control off instead.
