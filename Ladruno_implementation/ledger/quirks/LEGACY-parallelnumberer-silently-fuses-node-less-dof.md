---
wp: LEGACY
title: "ParallelNumberer silently FUSES node-less DOF groups (Lagrange multipliers) across ranks — getRef() returns −1 for all of them"
legacy_seq: 204
---
### `ParallelNumberer` silently FUSES node-less DOF groups (Lagrange multipliers) across ranks — `getRef()` returns −1 for all of them
- **Bites:** under MP with the Lagrange handler, every rank's multiplier DOF groups share ref=−1; the stock gather-merge dedups by ref ⇒ all of them collapse into ONE merged vertex ⇒ silently wrong numbering. Never observed in production only because the fork's MP lanes use Transformation/Plain handlers.
- **Workaround/status:** `LadrunoParallelNumberer` hard-errors on ref<0 with a message naming the gap; stock still fuses. Use Transformation/Penalty handlers under MP. *2026-07-22 (found in the ADR-74 N2 review).*
