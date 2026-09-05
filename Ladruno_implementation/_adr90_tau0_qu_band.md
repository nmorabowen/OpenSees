---
title: "ADR 90 / WP-A2 — the tau = 0 collapse-load band on softening SANISAND"
project: Ladruno
type: measurement report (ADR-90 un-descope gate)
status: "RESULTS 2026-09-05 — verdict in §1"
owner: nmora
adr: ADR 90
engine_build: 548fe911427e90a2edfead05cb3672a738d25b6d
driver: Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py
related:
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[_adr90_regularization_planning_brief]]"
  - "[[_adr90_a0_results]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[86_ladruno_sanisand_apegmsh_emitter_guide]]"
  - "[[testbed/00_canonical_testbed]]"
  - "[[LEDGER_quirks]]"
tags: [adr90, regularization, sanisand, localization, collapse-load, measurement]
updated: 2026-09-05
---

# ADR 90 / WP-A2 — the tau = 0 collapse-load band on softening SANISAND

> [!abstract] The question, and why it is the un-descope gate
> ADR 90's warrant rests on a premise nobody had measured: that the
> **unregularized** collapse load of the campaign's own problem is *not*
> already mesh-converged inside the campaign's own tolerance. The fork has that
> measurement for `DruckerPrager` — `tests/test_r3_prandtl_collapse_gate.py`
> reads 1.0849 / 0.9977 / 0.9514 of exact Prandtl–Reissner at
> h0 = 1.0 / 0.5 / 0.25, every leg a genuine plateau — and did **not** have it
> for `LadrunoSANISAND`, whose non-associated **softening** is the reason
> ADR 90 exists.
>
> **Nothing in `SRC/` was changed and no build was made.** This measures ADR 90
> §3.1 option (d) — *"do nothing to the lane; report the band"* — which that
> section itself calls *"the negative control of the whole ADR, not an option"*.

PLACEHOLDER
