---
wp: LEGACY
title: "Domain::revertToLastCommit() sets dT = 0 and re-applies the load — so a rate-dependent material runs the FIRST evaluation of every retried step UNREGULARIZED"
legacy_seq: 353
---
### `Domain::revertToLastCommit()` sets `dT = 0` and re-applies the load — so a rate-dependent material runs the FIRST evaluation of every retried step UNREGULARIZED
- **Bites:** a run with step cutbacks silently mixes regularized and unregularized steps. With the fork's `dt <= 0 => beta = 1` convention (deliberate: a missing time increment must not turn a material elastic) the first evaluation after every cutback takes the **inviscid** branch, which **dumps the entire accumulated overstress in one committed step** — a finite stress drop with no strain increment. Nothing in the output says it happened.
- **Why:** `Domain::revertToLastCommit()` does `currentTime = committedTime; dT = 0.0; this->applyLoad(currentTime);` (`SRC/domain/domain/Domain.cpp:2334-2339`). The same `dt <= 0` path is reached by `loadConst` (which makes the increment negative). A **held-load** step is the mirror image: there `dt > 0` but the strain rate is zero, so the model relaxes fully toward the inviscid backbone with time constant `tau` — **staged geostatic steps do exactly this**.
- **Workaround:** count the steps that commit with `tau > 0` but `beta == 1`, expose the count in the provenance block, and **fail acceptance on a non-zero count**. Keep `tau = 0` during staging and gate post-gravity byte-identity against a run with no wrapper at all. "Inert without a positive `ops_Dt`" is the wrong mental model — it is not inert, it discharges. Learned 2026-09-05, ADR-90 3-lens review, [[90_ladruno_viscoplastic_regularization_adr]] §4.2.
