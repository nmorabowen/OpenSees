---
wp: ADR-69
title: "Overlay energy accounting: the +Q.p coupling work rides ADR-69's ULW (external-load-work) channel -- closure holds, attribution is merged"
legacy_seq: 192
---
### Overlay energy accounting: the `+Q.p` coupling work rides ADR-69's ULW (external-load-work) channel -- closure holds, attribution is merged
- **Bites:** reading an ADR-69 `EnergyBalanceRecorder` breakdown on an overlay run, there is no separate "pore-coupling work" channel and the external-load work looks inflated -- you cannot tell coupling work from genuine external load work in the per-channel numbers.
- **Why:** the overlay injects its forces through `Node::addUnbalancedLoad`, and the ADR-69 kernel's ULW = integral of v^T P_ext dt reads `Node::getUnbalancedLoad` -- so the `+Q.p` forces are INSIDE the external-work channel by construction (verified at P3 pin 3.A). The closure residual ERR therefore stays within the ADR-69 bound (measured by battery gate (g)); only the ATTRIBUTION is merged, not the balance.
- **Workaround/status (2026-07-18, ADR-73 P3):** documented, not silently absent -- energy closure on overlay runs is trustworthy; per-channel attribution of coupling work is not separable. P4 may split it into its own channel; no recorder code change at P3.
