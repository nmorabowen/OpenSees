---
wp: LEGACY
title: "A material cannot refuse at commit: Domain::commit() discards commitState()'s return and returns 0 unconditionally"
legacy_seq: 379
---
## A material cannot refuse at commit: `Domain::commit()` discards `commitState()`'s return and returns 0 unconditionally

Found by the ADR-92 red/blue review (B3), diagnosing why the `-implexControl`
commit-time companion refusal in `LadrunoSANISAND::commitState()` never
propagated to the analysis.

`LadrunoSANISAND.cpp:1660` returns `LADRUNO_MATERIAL_REFUSED` from
`commitState()` when the IMPL-EX companion (`-maxSubsteps` scheme 1) exhausts
its substep budget. But `Domain.cpp:2244` calls elements at commit as a bare
`elePtr->commitState();` inside a `while` loop with the return value discarded,
and `Domain::commit()` itself ends `return 0;` **unconditionally** — there is
no code path in it that can propagate a non-zero child return. Whatever
`AnalysisModel::commitDomain()` does with a non-zero `commit()` result, it can
never see this refusal: the element and material returns are thrown away two
frames before that check. The step is silently accepted, `committedTime`
advances, recorders write against the un-refused state, and — because the
early `return` in `commitState()` sits before the history-variable update at
`:1673-1682` — `mEpsilon_n`, `mImplexDtCommit`, `mImplexStepArmed` and
`mImplexDEpsP` are all left stale for every subsequent step.

**The rule that follows:** commit-time is not a refusal site by construction
in this codebase — the object model has no channel for it. A material (or
element) that needs to refuse work discovered only at commit must make that
refusal *observable* some other way: a throttled `opserr`, a running counter
exposed as an `eleResponse`/`setParameter` value a driver or recorder can poll
(the model this fork already uses for `n_material_refused`-shaped counts), or
a material-level flag a recorder reads. Returning a sentinel from
`commitState()` and trusting it to reach the analysis does nothing — the same
`Domain::commit()` discard applies to every element, not just
`LadrunoSANISAND`'s.
