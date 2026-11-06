---
wp: WP-163
title: "Ordered comparisons with NaN are always false: a min/max accumulator silently keeps its pre-NaN extremes (WP-163)"
date: 2026-10-03
---
### Ordered comparisons with NaN are always false: a min/max accumulator silently keeps its pre-NaN extremes (WP-163)
- **Bites:** the `-envelope` sink updated MIN/MAX/ABSMAX with `<`/`>`; a run that went NaN at step k (explicit
  runs commit NaN — CDL does not trap it) kept the finite pre-divergence extremes, so the envelope looked healthy;
  a NaN FIRST sample stuck forever.
- **Workaround/status:** NaN is now deliberately sticky (first NaN poisons MIN/MAX/ABSMAX, ARG_STEP = that step).
  Do not build with `/fp:fast` / `-ffast-math`: the `v != v` test (and `std::isnan`) is not reliable under it.
  The size-mismatch rows of an element bucket are NaN-filled for the same reason (a zero row looked like data).
