---
wp: LEGACY
title: "Newton cannot reach a tight NormDispIncr on SANISAND — the substepped ModifiedEuler return makes the discrete map only piecewise smooth and the residual STALLS…"
legacy_seq: 356
---
## Newton cannot reach a tight `NormDispIncr` on SANISAND — the substepped `ModifiedEuler` return makes the discrete map only piecewise smooth and the residual STALLS around 1e-6 m

**Found 2026-09-05, ADR-90 WP-A2.**

`ManzariDafalias`'s default integration scheme is a **substepped** `ModifiedEuler` return with a
**hardcoded** substep tolerance `TolE = 1e-4` (exactly what the fork's `-honorTolR` flag exists to
expose — `86_ladruno_sanisand_apegmsh_emitter_guide` §1). The stress a Gauss point returns is
therefore a piecewise-smooth function of the strain increment, the assembled residual inherits
that, and **Newton stops converging quadratically and stalls.**

Measured on a strip-footing leg (h0 = 0.5, 390 hexes): over 47 failed convergence attempts the
residual displacement norm stalls at a **median of 6.6e-7 m** (min 3.4e-8, max 1.3e-4). A
`test NormDispIncr 1e-8` target is therefore not merely tight, it is **unreachable** — and the run
that nominally used it was in fact carried by the relaxed third rung of its algorithm ladder on
**18 of its 26 steps**, i.e. 65 of every 125 state-determination passes were spent failing two
rungs that could not succeed. A study in that state is measuring its own convergence test: the
ADR-63 note-71 failure mode.

- **Use `NormUnbalance`** — which is also what `tests/test_r3_prandtl_collapse_gate.py` actually
  uses (`tol = 1.0e-5 * want`), despite "NormDispIncr per the R3 gate" appearing in more than one
  downstream brief.
- **Scale the force tolerance to a deck-intrinsic force** (the model's own weight `gamma*V`, the
  applied load, ...). Measured, same deck and wall budget: `NormDispIncr 1e-8 m` reached
  s/B = 0.00106; `NormUnbalance 1e-6 gamma*V` reached 0.00218; `NormUnbalance 1e-5 gamma*V`
  reached 0.00442 — with the answer moving by a median of 0.3-0.75 % between all three at matched
  settlement.
- **A displacement-norm tolerance is not mesh-neutral, which is disqualifying inside a
  mesh-convergence study.** The norm runs over the free DOFs, and there are 3.6x more of them at
  h0 = 0.25 than at h0 = 1.0, so the same nominal number is a different physical requirement on
  each mesh of the sequence. A force tolerance scaled by `gamma*V` is identical on all three by
  construction.

> **DOCUMENTED, NOT CHANGED, in WP-86b (ADR-86b T3, PR pending).** Written up as a deck rule in
> `86_ladruno_sanisand_apegmsh_emitter_guide.md` §6 with the measured table. **Deliberately NO
> runtime warning:** a material cannot see which convergence test the deck installed, so any check
> would have to live in the analysis layer and would fire on every non-SANISAND deck that ever uses
> a displacement norm. The stall is a property of the substepped return, not of the deck, so the
> guidance is the fix.
