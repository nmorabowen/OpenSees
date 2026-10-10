---
wp: ADR-69
title: "ADR-69 P2 -- 1 vanilla row(s)"
files: ["`SRC/analysis/integrator/IncrementalIntegrator.{h,cpp}`"]
table: "main"
legacy_seq: [237]
---
| `SRC/analysis/integrator/IncrementalIntegrator.{h,cpp}` | `// Ladruno` ADR-69 P2: modal-damping energy publisher. `addModalDampingForce()` additionally captures the dissipation rate `-(dampingForces . v) = v'C_modal v` at each iterate (OVERWRITE, never accumulate — the last iterate of a step is the converged one); `commit()` trapezoid-integrates over the domain-time delta and publishes to `MODAL_WORK` BEFORE `commitDomain()` (recorders read the registry inside), declaring the channel at the first valid commit — which precedes the recorder's first-record column fixing. First valid commit seeds only (the step's start time predates the first modal `formUnbalance`). `mdRateValid` cleared per commit so a stale rate can never re-publish; failed/reverted steps never publish. 5 in-class-initialized scratch members in the header. COVERAGE: integrators that override `commit()` without chaining to the base (HHT family, `*_TP` explicit) do not publish — their modal dissipation stays in RES (documented; the canonical modal-damping pairing, Newmark, uses the base commit). | ADR-69 P2 |
