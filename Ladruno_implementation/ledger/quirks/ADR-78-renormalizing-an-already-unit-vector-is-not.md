---
wp: ADR-78
title: "Renormalizing an already-unit vector is NOT bitwise-idempotent — serialization round trips through a validating constructor drift by 1 ulp (ADR 78 P2 review F2)"
legacy_seq: 300
---
### Renormalizing an already-unit vector is NOT bitwise-idempotent — serialization round trips through a validating constructor drift by 1 ulp (ADR 78 P2 review F2)
- **Bites:** any unpack/restore path that rebuilds state through the same validating entry point that normalized it originally. `sqrt(n·n)` of a stored unit normal is 1±1ulp for generic directions (measured: re-dividing changed bits on 3368 of 10000 random unit vectors), so `save → restore → re-verify` reported "definitions differ" on a model nobody touched — the verify instrument poisoned by its own rebuild path.
- **Why it evades the usual guards:** axis-aligned test vectors ((0,0,1) etc.) have EXACT norms, so every convenient test normal passes bit-exact; only a tilted normal exposes it. And behavioral gates can't see 1 ulp — only a bit-compare can.
- **Rule:** make normalization idempotent at the choke point (`|nrm−1| < 1e-12 ⇒ passthrough`), and gate serialization with (a) a tilted/irrational-normed vector in the test model and (b) a restore→re-save BYTE compare of the packed payload, not just response parity.
