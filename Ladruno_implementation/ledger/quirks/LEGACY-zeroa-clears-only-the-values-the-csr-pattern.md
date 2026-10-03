---
wp: LEGACY
title: "zeroA() clears only the VALUES — the CSR pattern (colA/rowStartA) is built once in setSize() and survives every assembly, so OpenSees already satisfies \"freeze…"
legacy_seq: 217
---
### `zeroA()` clears only the VALUES — the CSR pattern (`colA`/`rowStartA`) is built once in `setSize()` and survives every assembly, so OpenSees already satisfies "freeze the sparsity"
- **Bites (in the good direction):** a plan that budgets real work for building and maintaining a frozen sparsity graph before threaded assembly (the expensive half of the Kratos atomic-scatter pattern). It is already paid. In `PARDISOGenLinSOE`: `colA`/`rowStartA` are filled once from the DOF graph (`:221-297`); `zeroA()` zeros `A[]` and clears `factored` and touches neither; `addA` (`:370-478`) then locates each target by a **read-only** linear search over the frozen row and does exactly one read-modify-write, `A[k] += m(i,j)`.
- **Consequences:** (1) the entire implicit assembly race is *one* `+=` on a shared `double` at an index computed from immutable data — an atomic on one statement per SOE, no coloring, no thread-private matrices; (2) but the inner search makes the scatter **O(idSize² × rowlen)**, a pre-existing *serial* inefficiency worth fixing independently of threading (this is the "check is O(nnz) against an O(nnz × rowlen) loop" shape); (3) anything that can resize the SOE mid-run (`domainChanged` from ADR-51 element removal / ADR-60 contact re-emission) invalidates the freeze and must force a threaded phase **off**, not race with it — same family as the `-factorOnce` staleness caveat.
- **Workaround/status:** recorded as settled evidence in [[75b_ladruno_threaded_assembly_adr]] §2.2/§4.1. *2026-07-25 (ADR-75b L3-0).*
