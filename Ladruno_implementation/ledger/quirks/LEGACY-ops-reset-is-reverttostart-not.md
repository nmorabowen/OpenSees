---
wp: LEGACY
title: "ops.reset() is revertToStart, NOT revertToLastCommit — it ZEROES the committed state, and a probe that uses it to \"return to the base state\" measures a pristin…"
legacy_seq: 377
---
## `ops.reset()` is `revertToStart`, NOT `revertToLastCommit` — it ZEROES the committed state, and a probe that uses it to "return to the base state" measures a pristine material

Cost half a day on ADR-92 P1's tangent-identity gate, which failed on its own
precondition with a "before" of `-100.2 kPa` and an "after" of
`[-0.0101, -0.0101, -0.0101, -0.0, -0.0, -0.0]`.

**The chain, at source.** `ops.reset()` → `OPS_resetModel()`
(`SRC/interpreter/OpenSeesCommands.cpp:2534`) → `Domain::revertToStart()` →
every material's `revertToStart()` → `ManzariDafalias::initialize()`, which
does `mSigma_n.Zero(); mEpsilon_n.Zero(); mAlpha_n.Zero(); mFabric_n.Zero()`.
The entire committed state is gone. `mElastFlag` is a **static** and is NOT
touched, so the material is still on the plastic stage and the next state
determination runs the full plastic path against a zero committed state.

**Then `Domain::revertToStart()` ends `return this->update();`** — like
`revertToLastCommit()`, it pushes one state determination through every element
before returning. So the value a probe reads *after* `ops.reset()` is not the
committed stress; it is whatever the material returns from a zero-strain
increment against a zeroed committed state. Under `LadrunoSANISAND -implex`
that is `sigma~ = 0`, which trips the ADR-92 `p_min` floor clamp and comes back
as `p_min*I1` — `[-0.0101]*3` at the element face, with three IEEE `-0.0`
shears from `getStress()`'s `-1.0 *` flip. **The `-0.0` triple is the
fingerprint**: it can only come from a stress tensor whose deviator is exactly
zero, which no real load path produces.

**What to use instead.** A `StaticAnalysis` step that fails already calls
`the_Domain->revertToLastCommit()` and `theIntegrator->revertToLastStep()`
itself (`StaticAnalysis.cpp:185` and four sibling sites). So a probe built on a
deliberately-failed `analyze(1)` needs **no** reset call at all — the correct
revert has already happened, and adding `ops.reset()` destroys it. There is no
Python verb that calls `Domain::revertToLastCommit()` on its own; a failed step
is the way to reach it.
