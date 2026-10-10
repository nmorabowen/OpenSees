---
wp: LEGACY
title: "MKL dfeast_srci SHRINKS its in/out m0 argument in place — don't feed the shrunken value back into your own subspace bookkeeping"
legacy_seq: 165
---
### MKL `dfeast_srci` SHRINKS its in/out `m0` argument in place — don't feed the shrunken value back into your own subspace bookkeeping

- **Bites:** any FEAST RCI driver that passes its subspace-size variable by pointer
  into `dfeast_srci` and afterwards compares the found count `m` against that same
  variable (saturation / auto-enlarge logic). Surfaced building ADR 43 P3c-serial:
  a band holding 5 modes with a seeded `m0 = 15` came back with `m0` REWRITTEN to 5,
  so `m == m0` looked like subspace saturation, the enlargement re-ran from the
  shrunken value, plateaued (`1.5*5+8 → shrunk to 5 → again`), and after
  `maxEnlarge` retries the solver REFUSED a perfectly complete result.
- **Why:** the RCI's `m0` is in/out session state — when FEAST's reduced Rayleigh–Ritz
  detects the band holds fewer modes than the subspace, it adapts the working
  subspace DOWN and records that in `m0`. That is convergence bookkeeping, not a
  saturation verdict. (The packed driver `dfeast_scsrgv` does NOT exhibit this —
  its `m0` comes back unchanged, which is why the P1/P2 battery never saw it.)
- **Fix:** give the RCI session its OWN `m0` copy (`int m0Rci = m0;` …
  `dfeast_srci(..., &m0Rci, ...)`) and keep the caller's requested `m0` for the
  `m >= m0` saturation compare. Post-shrink, the live right-hand sides at
  `ijob=11` are the first `m0Rci` columns of `workc` — solve with `m0Rci` RHS,
  not the original `m0`. See `runFeastRci` in `FeastEigenSolver.cpp`.
