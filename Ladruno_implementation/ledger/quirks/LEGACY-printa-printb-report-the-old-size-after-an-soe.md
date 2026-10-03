---
wp: LEGACY
title: "printA / printB report the OLD size after an SOE SHRINKS — a silent wrong answer, no crash"
legacy_seq: 162
---
### `printA` / `printB` report the OLD size after an SOE SHRINKS — a silent wrong answer, no crash

- **Bites:** any model whose equation count goes DOWN between analyses on the same
  SOE under `system FullGeneral` — a second `analyze` with more DOFs constrained,
  staged construction, `remove element`, a contact set that retires. `printA`/
  `printB` then hand back the PREVIOUS, larger system: measured on a 12 -> 6
  shrink, `systemSize()` correctly reported 6 while `printA("-ret")` returned 144
  values and `printB("-ret")` 12. 108 stale entries presented as the tangent.
- **Why:** `FullGenLinSOE::setSize()` built `vectX`/`vectB`/`matA` with **`Bsize`**
  (and `Bsize x Bsize`), which is a grow-only high-water CAPACITY, not `size`, the
  live equation count. `Asize`/`Bsize` only ever grow, so after a shrink the
  wrappers describe the old system over the same storage. It never faulted because
  `Bsize**2 == Asize` exactly — the read stays in bounds, it is just wrong.
- **Why it hid for so long:** it is not a crash, and it is invisible unless you
  compare the returned LENGTH against an independent oracle. `len(B)` is NOT that
  oracle — it is the very quantity that goes stale, so a probe that checks
  `len(A) == len(B)**2` passes happily on the bug. Use `systemSize()`.
- **Status (2026-08-18): FIXED** — wrappers now dimensioned by `size`. `size` is
  correct by construction: `addA` indexes `A + col*size`, `addB`/`setB` over
  `[0,size)`. **FullGen was the lone outlier**: BandGen, BandSPD, ProfileSPD,
  SProfileSPD, Diagonal, SparseGenCol and SymSparse already used `size` (surveyed,
  so this quirk is FullGeneral-only). Allocation still tracks capacity — only the
  wrapper extent changed, no reallocation behaviour moved. Gate:
  `tests/test_printa_unsized_soe.py::test_printa_after_shrink_reports_the_live_size`.
- **Sibling of the entry below**, and the reason to read both: that one is `printA`
  KILLING the process on an SOE that was never sized; this one is `printA` LYING
  about an SOE that was sized and then shrank. Same command, same file, opposite
  symptom — one is impossible to miss, the other impossible to notice.
