---
title: ADR-80 P3 — LadrunoLoadControl -tangentPredictor MEETS its acceptance gate exactly
project: Ladruno
status: complete (positive result)
priority: high
owner: nmora
tags:
  - implementation
  - solver
  - integrator
  - constraints
  - validation-gate
aliases:
  - ADR-80 P3 verdict
  - tangentPredictor
updated: 2026-08-04
parent: "[[80_ladruno_sp_imposition_strengthening_adr]]"
---

# ADR-80 P3 verdict — candidate C recovers the whole excess

> The S1 report ([[80c_s1_extrapolate_verdict_2026-08-04]]) closed by saying
> *"open D6/Candidate C as the next step rather than tuning `frac`."* This is
> that step, and it is the first ADR-80 remedy to pass its own gate.

## Result, against the gate S1 failed

Same model, same marches, same binary. Everything comes from
`sp_gates/p3_tangent_predictor_acceptance.tcl` →
`p3_tangent_predictor_acceptance_2026-08-04.json`.

### The gate that matters — adaptive march with cutbacks, `maxIter` = 5

| variant | increments | **cutbacks** | iterations |
|---|---|---|---|
| elastic control (the ceiling) | 6 | 0 | 12 |
| stock `LoadControl` | 43 | **23** | 224 |
| `LadrunoLoadControl` (no flag) | 43 | 23 | 224 |
| `-extrapolate 1.0` | 43 | 23 | 224 |
| **`-tangentPredictor`** | **6** | **0** | **12** |

**Cutbacks 23 → 0. Iterations 224 → 12, ×18.7 → ×1.00.** The tangent route does
not approach the elastic control, it *equals* it — every digit.

### Fixed march, both caller idioms

| idiom | variant | iterations |
|---|---|---|
| — | elastic control | 20 |
| reissue | stock / no flag / `-extrapolate 1.0` | 60 / 60 / **60** |
| reissue | **`-tangentPredictor`** | **20** |
| persist | stock / no flag / `-extrapolate 1.0` | 60 / 60 / 35 |
| persist | **`-tangentPredictor`** | **20** |

The `-extrapolate` row reproduces 80c exactly, including its defect: **35 under
`persist`, 60 (completely inert) under `reissue`**. The tangent route scores 20
under both, because it carries no cross-step state at all.

### Controls that had to not move, and did not

- **Answer:** `u_top` = 0.150000 and — the probe that actually discriminates —
  **`u_mid` = 0.075000 exactly in all 48 rows.**
- **Traction drive:** 20 (fixed) / 12 (adaptive) with and without the flag.
- **Elastic material + `-tangentPredictor`:** 20 / 12, i.e. identical to the
  elastic control. No cost, no harm where there is nothing to fix.
- **`-extrapolate 0` is still bit-identical to stock** in every row.

⚠ **`u_top` is the DRIVEN node and is not a discriminator** — it reports the
prescribed value back to itself and would read 0.150000 even for a step that
moved nothing but the boundary layer. The first version of both the harness and
the pytest probed only `u_top`; that was a vacuous check and is corrected. Every
figure above is confirmed at the midpoint.

## What it does

`-tangentPredictor` is the Kratos `b −= K·Δu_D` route — the D6 candidate,
Kratos's `use_old_stiffness_in_first_iteration`.

```
newStep()      Domain::applyLoad(lambda)   <- sets every sp VALUE
                                              but does NOT reach the handler,
                                              so enforceSPs never runs and the
                                              elements stay at the COMMITTED state
iteration 1    formTangent  -> K at the committed state
               formUnbalance-> P_ext - R_int(committed),  then  b -= K*du_D
               solve        -> dU = the interior's linearized response
update(dU)     incrDisp(dU); NOW applyLoadDomain(lambda) -> enforceSPs writes the
                                              driven face and pre-updates its
                                              layer, whose neighbours have
                                              already moved
iterations 2+  stock path, unchanged
```

The constitutive law is **never evaluated in the lagging state**. That is the
whole difference from `-extrapolate`, which moved the starting guess closer but
left an overstrain that still tripped the spurious yield.

### Why it works where the history predictor could not

Three properties, all of them consequences of being **stateless**:

1. **It fires on the first increment.** `-extrapolate` cannot — there is no
   previous increment to extrapolate.
2. **It fires after a cutback.** `-extrapolate` deliberately skips there
   (the failed attempt's partial ΔU is garbage). Cutbacks are exactly where the
   cost lives, so S1 was structurally blind to its own headline metric.
3. **It survives `integrator LadrunoLoadControl ...` being re-issued every
   step** — the fork's robust-solve driver and g5's ladder. That idiom destroyed
   `-extrapolate`'s stored increment and made it measure *zero* (80c §2, armed =
   0). Measured again here: `reissue`/`ext1` = 60 against `persist`/`ext1` = 35.

## Implementation, and one finding

`LadrunoLoadControl` gains `-tangentPredictor` (no new class, no new class tag).
The route needed one thing that did not exist, and getting it produced a
correction to the ADR:

> ⚠ **ADR-80 D6 named `TransformationFE::getTangForce()` as "the hook" for this
> route. It cannot be.** Its only input is a GLOBAL vector indexed by `myID`,
> and a Transformation-eliminated dof has `myID == -1` — the prescribed
> increment is *not representable in its argument*. The same reason the route is
> needed at all (no equation ⇒ no column) is the reason that signature cannot
> carry it. The stub is left untouched; a new `getSPTangentForce(Integrator*)`
> was added instead, which reads the increment off the `SP_Constraint` through
> the constrained node's `TransformationDOF_Group`.

Four additive vanilla touch-points, all no-ops with no stock caller
(see [[LEDGER_vanilla_files]]): `DOF_Group::getSPDispIncr` (default 0) +
`TransformationDOF_Group`'s override; `FE_Element::getSPTangentForce`
(default 0) + `TransformationFE`'s override, which builds `du` in original
element-dof space, forms `f = K du` at the committed state, and returns `T^t f`
so the caller can post it with the element's ordinary `getID()` — whose `-1`
entries drop the eliminated rows for free.

### The hazard this route creates, and the guard for it

Iteration 1 deliberately runs **before** the sps are enforced. If the `K·Δu_D`
term ever fails to assemble, `dU` comes back ≈ 0 and a dU-based convergence test
declares victory on a step that never moved — **precisely the
`AutoConstraintHandler` silent wrong answer that ADR-80 S2 just fixed**. So
`formUnbalance()` counts the contributing elements and, on zero, enforces the sps
immediately, re-assembles, and **disables the route for the rest of the run**
with a note. Reached on (a) a model with no non-homogeneous `sp` — the common
case, verified in the traction rows — and (b) a constraint handler whose sp'd
elements are not `TransformationFE`s. Gated by
`tests/test_adr80_tangent_predictor.py` (11 cases, Zone-A), including the
interior-did-move assertion stated positively so it cannot pass vacuously.

`-extrapolate` and `-tangentPredictor` **refuse to compose** (they measure their
increments from different origins and would double-count); the flag combination
is a loud warning that disables `-extrapolate`.

## Scope — what this measurement does NOT license

- **The synthetic ceiling is ×18.7 here, not the field's ×28.6.** Cerro Lindo is
  a different model with 52 cutbacks. Field validation stays open and is still
  gated on those decks being available (ADR §"Field validation").
- **The model is built so that nothing yields physically** (converged σ = 300 MPa
  against `f_y` = 379.5). The G3 SCOPE LIMIT applies unchanged: **in a model that
  genuinely yields, part of the tangent collapse is real** and no predictor can
  recover it. Do not quote "recovers the full excess" for a plastically-hinging
  model.
- **Contact and `-geom` are untested here.** ADR §Risks flagged both; nothing in
  P3 addresses them.
- **Cost when it fires:** one extra `applyLoadDomain` per increment, plus one
  element-tangent formation per sp-touching element per increment. No extra
  solve — the predictor IS iteration 1, which the algorithm was going to run
  anyway. On these decks it is a net saving of 212 of 224 iterations.

## Re-verification on landing — 2026-09-05

The 2026-08-04 snapshot sat uncommitted in a dead worktree for a month and was
rescued in the ADR-87 D10 sweep 2 (`branch_graveyard.md`). Landed on
`wp/80-p3-tangent-predictor` (PR #786) as one commit on top of `ladruno`
`f71b9196b`, 229 commits past the snapshot's base, and re-verified there:

- Build `4410bd6a0` (`ladrunoBuild` == HEAD on both `OpenSees.exe` and
  `opensees.pyd`).
- `sp_gates/p3_tangent_predictor_acceptance.tcl` re-run →
  `p3_tangent_predictor_acceptance_2026-09-05.json`, **byte-identical** to the
  2026-08-04 JSON: adaptive `maxIter` 5, J2, stock 43 inc / 23 cutbacks / 224
  iters; `-tangentPredictor` 6 / 0 / 12; `u_mid` = 0.075000 on every row.
- `tests/test_adr80_tangent_predictor.py`: 11 passed.

So the verdict above survives ADR-76 (tangent reuse) and ADR-77 (implicit
transient) having moved the analysis subsystem underneath it in the meantime.
