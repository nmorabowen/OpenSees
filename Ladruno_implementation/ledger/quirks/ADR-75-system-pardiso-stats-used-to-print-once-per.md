---
wp: ADR-75
title: "system Pardiso -stats used to print once per PATTERN, not once per FACTORIZATION (ADR-75 P1k)"
legacy_seq: 427
---
## `system Pardiso -stats` used to print once per PATTERN, not once per FACTORIZATION (ADR-75 P1k)

The P1d `-stats` implementation gated its print on a `statsDone` flag reset only
when the symbolic phase (11) re-ran — i.e. once per sparsity PATTERN. That reads
as "once per factorization" until you actually refactorize the same pattern:
`ModifiedNewton` holding `A` fixed never re-triggers it (correct — nothing new
to report), but plain `Newton` under `LoadControl`/`StaticIntegrator` DOES
reassemble+refactorize every step against the SAME sparsity pattern (`zeroA()`
sets `factored = false` on every tangent assembly, `setSize()`/phase 11 is not
re-run), and the old `-stats` printed its block exactly ONCE across the whole
analysis — silently dropping every later refactorization's numbers. Not a bug
exactly (the counters ARE mostly pattern-invariant — peak/permanent symbolic
memory really doesn't change step to step), but it means "grep the log for
`PARDISO stats:`" answered "did the model fit at all", not "how many times did
this thing refactorize and what did each cost" — which is what TIMs PM-01 D26
actually wanted. Fixed by gating the print on a per-`solve()`-call local
(`didFactorNow`, true only when THAT call executed phase 22) instead of a
member latch. Lesson for the next "once per X" cache-style gate on a `-stats`
style flag: ask explicitly whether X is "this analysis" or "this event", because
OpenSees's own re-solve/re-factor cadence (one pattern, many refactorizations)
makes those two very different answers.
