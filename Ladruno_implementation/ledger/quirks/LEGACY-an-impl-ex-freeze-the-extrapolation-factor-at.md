---
wp: LEGACY
title: "An IMPL-EX \"freeze the extrapolation factor at the first trial call of the step\" rule is consumed by Domain::revertToLastCommit()'s own trailing update()"
legacy_seq: 378
---
## An IMPL-EX "freeze the extrapolation factor at the first trial call of the step" rule is consumed by `Domain::revertToLastCommit()`'s own trailing `update()`

Found while diagnosing the above; it is a defect in ADR-92 P1's own code and it
hit the **default** `-implexDt pseudo` source.

`Domain::revertToLastCommit()` does not stop at reverting. It sets `dT = 0.0`,
re-applies the committed load, and ends `return this->update();` — one
zero-strain-increment state determination through every material at
`ops_Dt == 0`. A material that freezes a per-step quantity "on the first trial
call after a commit or revert" will therefore freeze it **from that spurious
call**: `dt = 0`, extrapolation factor `f = 0`, arm consumed. The retried step
then runs its entire ladder rung with the plastic extrapolation switched off,
and the reported `implexError` is the error of an operator the deck never asked
for. Nothing crashes and the committed answer stays correct (the companion
return at commit is untouched), so this is invisible without looking for it.

**The rule that works:** arm from the first trial call whose strain increment is
**non-zero**; a zero increment takes `f = 0` for that evaluation only (which is
what a `LoadControl 0.0` hold requires anyway — no strain advanced, no plastic
flow predicted) and leaves the step armed. Under a pseudo-time source this
cannot perturb any step that actually moves, because `ops_Dt` is identical on
every call of a step, so *which* call arms it is unobservable.

This is the same `dT = 0` trap already recorded for rate-dependent materials
("`Domain::revertToLastCommit()` sets `dT = 0` and re-applies the load — so a
rate-dependent material runs the FIRST evaluation of every retried step
UNREGULARIZED"). **Anything that keys off "the first call of a step" must
assume that call is a zero-increment ghost.**
