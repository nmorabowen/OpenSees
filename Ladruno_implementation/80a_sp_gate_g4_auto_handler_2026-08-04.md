---
title: ADR-80 gate G4 — AutoConstraintHandler missing updateElement(): reproduced, fixed, gated
project: Ladruno
status: complete
priority: medium
owner: nmora
tags:
  - implementation
  - solver
  - constraints
  - validation-gate
aliases:
  - ADR-80 G4
updated: 2026-08-04
parent: "[[80_ladruno_sp_imposition_strengthening_adr]]"
---

# ADR-80 gate G4 — `AutoConstraintHandler` silent wrong answer

> Gate G4 of [[80_ladruno_sp_imposition_strengthening_adr]] (S2). The other P0
> gates (G1 mesh scaling, G2 step scaling, G3 tangent-vs-residual dominance)
> need the Cerro Lindo fuse decks and are **not** covered here. G4 was
> self-contained and ran first.

## The question

[[80_sp_prescribed_displacement_findings]] §4 candidate B claimed, **from source
reading only**, that `AutoConstraintHandler::applyLoad()` omits the
`updateElement()` loop that `TransformationConstraintHandler::enforceSPs` has at
`:518-521`, and that with a non-homogeneous `sp` this yields *a wrong answer, not
an error*. The ADR required a reproducer before any fix. Refuting the inference
was an equally valid outcome.

**Verdict: CONFIRMED, and the inference understated the blast radius.**

## Method

Two models, both driven to `δ = 3.0` with a free interior DOF whose exact answer
is `δ/2 = 1.5` by symmetry:

- **truss** — two equal `Truss` elements in series; node 1 fixed, node 3 driven,
  node 2 free.
- **brick** — two stacked `stdBrick`s with `ux`,`uy` fixed everywhere (a 1D chain
  of two identical brick springs); bottom face `uz` fixed, top face `uz` driven,
  mid plane free.

Crossed over `constraints` ∈ {`Transformation`, `Auto`, `Plain`} × `test` ∈
{`NormDispIncr`, `EnergyIncr`, `NormUnbalance`}, all pinned to `Newton` +
`LoadControl 1.0`, tol `1e-10`. `Transformation` is the reference; the gate is
"`Auto` must match it".

Deck and artifacts (JSON, per the ADR's saved-artifact rule):

- `Ladruno_implementation/sp_gates/g4_auto_handler.tcl`
- `Ladruno_implementation/sp_gates/g4_auto_handler_PREFIX_2026-08-04.json`
- `Ladruno_implementation/sp_gates/g4_auto_handler_POSTFIX_2026-08-04.json`

## Result

Interior displacement (exact = 1.500000) and Newton iteration count.

| model | test | Transformation | **Auto — PRE-fix** | Auto — POST-fix |
|---|---|---|---|---|
| truss | `NormDispIncr` | 1.5 (2 it) | **0.0 (1 it)** ✗ | 1.5 (2 it) ✓ |
| truss | `EnergyIncr` | 1.5 (1 it) | **0.0 (1 it)** ✗ | 1.5 (1 it) ✓ |
| truss | `NormUnbalance` | 1.5 (1 it) | 1.5 (2 it) | 1.5 (1 it) ✓ |
| brick | `NormDispIncr` | 1.5 (2 it) | **0.0 (1 it)** ✗ | 1.5 (2 it) ✓ |
| brick | `EnergyIncr` | 1.5 (2 it) | **0.0 (1 it)** ✗ | 1.5 (2 it) ✓ |
| brick | `NormUnbalance` | 1.5 (2 it) | 1.5 (3 it) | 1.5 (2 it) ✓ |

**`analyze()` returned 0 in every pre-fix wrong-answer cell.** No warning, no
error, no diagnostic — the interior DOF simply never moved.

Post-fix `Auto` matches `Transformation` **on displacement and on iteration
count, in every cell** — including the two `NormUnbalance` cells where the
pre-fix path had been reaching the right answer by a longer route (2→1 and 3→2
iterations). That exact-march agreement is the expected outcome, since `Auto`
routes SP nodes through the same `TransformationDOF_Group`.

## Three things the gate changed about the ADR

1. **`EnergyIncr` is fooled too.** The findings note named only
   `CTestNormDispIncr`. Energy is `dU·R`, so every convergence measure built on
   `dU` collapses with it. The exposure is "dU-based tests", not one class.
2. **Test-dependence *is* the mechanism, and it is why this hid.**
   `NormUnbalance` sees the large *stale* residual, refuses to converge at
   iteration 1, and the next `LoadControl::update` → `updateDomain()` silently
   repairs the state. **The same deck is right or wrong purely by choice of
   convergence test.** When auditing an old `constraints Auto` result the
   question is not "did it converge" but "what did `test` say".
3. **The one-iteration convergence is the diagnostic tell.** A dU-based test
   converging in a single Newton iteration means the first residual carried no
   information about the prescribed motion. Kept as an explicit assertion in the
   regression test.

**Negative control:** `constraints Plain` also gives 0.0 — but it *warns*
(`non-homogeneos constraint for node N homogeneous constraint assumed`) and its
behaviour is unchanged by the fix. Plain is wrong and loud, by design; `Auto` was
wrong and **silent**. That contrast is the whole defect.

## Fix

`SRC/analysis/handler/AutoConstraintHandler.{h,cpp}` — a `theFEs` member holding
the FE_Elements that touch an SP-constrained node (membership computed in
`handle()` from the existing `sp_map`, O(1) per node), the `updateElement()` loop
in `applyLoad()`, and a `clearAll()` reset. Mirrors the Transformation handler.
`domainChanged()` always calls `clearAll()` immediately before `handle()`
(`StaticAnalysis.cpp:386`), so `theFEs` cannot accumulate duplicates.

Ledgered in [[LEDGER_vanilla_files]]; quirk row in [[LEDGER_quirks]] updated from
"source-level inference, NOT yet reproduced" to CONFIRMED + FIXED.

**Blast radius:** answer-changing *by design* for the broken case. Homogeneous-only
(`fix`-only) models are numerically unaffected — `enforceSPs` rewrites an unchanged
value, so the added update is idempotent; asserted by
`test_homogeneous_only_model_is_unaffected`.

## Gate: `tests/test_auto_handler_sp_update.py`

15 Zone-A tests (upstream `Truss`/`stdBrick` + `Elastic`, no MKL, no fork
elements). Asserts Auto-vs-Transformation agreement across all three test types
on both models, that the interior actually moves, the ≥2-iteration tell, the
Transformation premise, and the homogeneous-only blast-radius guard.

**Verified non-vacuous:** run against the stale pre-fix `opensees.pyd`, **10 of
15 fail** with exactly the documented signature (`u = 0.0`, `assert 1 >= 2`);
all 15 pass on the fixed module.

## What this does and does not settle for ADR-80

Settles **S2 only**. S2 was always independent of the predictor question — it is
a plain missing-update bug, not a conditioning problem.

It says **nothing** about S1 (`LadrunoLoadControl -extrapolate`): the ×28.6
prescribed-displacement penalty is a `Transformation` + `LoadControl`
*conditioning* pathology on a correct answer, whereas G4 was `Auto` returning an
*incorrect* one. G1–G3 remain the gates that decide S1, and G3 still bounds what
a predictor can buy. Do not read this result as evidence for the predictor.
