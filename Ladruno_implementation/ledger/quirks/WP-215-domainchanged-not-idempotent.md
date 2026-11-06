---
wp: WP-215
title: "Analysis::domainChanged() is not idempotent mid-analysis; Linear -factorOnce ignored it"
date: 2026-10-26
pr: "#1027"
---
### Calling the analysis' `domainChanged()` "to be safe" cuts a CentralDifference history; `-factorOnce` kept a stale factor

- **Bites:** a command that needs the model numbered (`ritzVectors`, any
  eigen-like helper) and calls `StaticAnalysis/DirectIntegrationAnalysis::
  domainChanged()` unconditionally. Mid-transient, `CentralDifference::
  domainChanged()` sets U(t − Δt) = U(t), so the next steps run from a lost
  velocity. Separately, `Linear -factorOnce` (factorOnce = 2 after the first
  step) skipped `formTangent` forever, even after a renumbering or after
  another command factored K into the analysis' LinearSOE.
- **Why:** `domainChanged()` re-reads the nodes into the integrator (W1 in
  ADR 98: integrators cache U/U̇/Ü and read nodes only there); `Linear` had
  no `domainChanged()` override.
- **Fix/status:** renumber only when `Domain::hasDomainChanged()` differs
  from the analysis' stamp (`ladrunoDomainStamp()`, eigen's rule), and
  `Linear::domainChanged()` re-arms `-factorOnce`. Tested by comparing a
  transient history with and without a `ritzVectors` call in the middle.
  Learned 2026-10-26 (WP-215, review of #1027).
