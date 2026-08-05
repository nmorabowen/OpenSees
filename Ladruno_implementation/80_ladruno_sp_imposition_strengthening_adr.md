---
title: Non-homogeneous SP imposition strengthening — a static predictor, not a new handler
project: Ladruno
status: S1 BUILT but its acceptance gate FAILED — D6 (getTangForce) triggered
priority: medium
owner: nmora
tags:
  - implementation
  - solver
  - constraints
  - integrator
aliases:
  - ADR-80
  - SP imposition strengthening
  - LadrunoLoadControl
updated: 2026-08-04
parent: "[[80_sp_prescribed_displacement_findings]]"
---

# ADR-80 — Non-homogeneous SP imposition strengthening

> Written on top of [[80_sp_prescribed_displacement_findings]] (PR #695), which
> records the mechanism, the ×28.6 measurement, the cross-code survey and the
> candidate ranking. This ADR takes the decisions and sequences the work.
> Every source reference below was re-verified in this tree on 2026-08-04.

## What

Three decisions, phased so that **no C++ is written before the validation
gates have run**:

1. **S1 (primary): `LadrunoLoadControl`** — a fork-owned static integrator,
   strict superset of stock `LoadControl`, adding a **`-extrapolate <frac>`
   displacement predictor** (Abaqus `EXTRAPOLATION=LINEAR` analogue).
   `-extrapolate 0.0` is the default and must be **bit-identical to stock**.
   Integrator classTag **33015** (next free in the ladruno integrator band;
   the ND-registry 33015 = `LadrunoRCConcrete` is a different registry, not a
   collision).
2. **S2 (independent bug): `AutoConstraintHandler::applyLoad()` missing the
   `updateElement()` loop** — reproduce first, then the ~5-line ledgered
   vanilla fix + regression test. **DONE — see [[80a_sp_gate_g4_auto_handler_2026-08-04]].**
   **Fork-only (decided 2026-08-04): NOT offered upstream.** The fix is
   self-contained and would apply cleanly upstream, but upstreaming is its own
   workstream (their review cadence, their regression bar, a behaviour-changing
   fix to a handler other people's decks depend on) and this ADR is not
   authorized to open it. Revisit deliberately, not as a side effect.
3. **S3 (gated diagnostics): iterate-level visibility** — an opt-in
   per-iterate dump (residual norm, ‖dU‖, count of Gauss points plastic in
   the iterate but elastic in committed state). Deferred until S1's gates
   are done; scoped here so it is not forgotten.

**Explicitly NOT in scope** (findings §4, upheld): a new SP constraint
handler; `Penalty`/`Lagrange` as workarounds; `updateDomain()` in
`LoadControl::newStep`; implementing `TransformationFE::getTangForce()`
(Kratos-style `b −= K·Δu_D`) — **demand-driven only, if S1 proves
insufficient**; any arc-length/relaxation device (the elastic control has 0
cutbacks and no limit point — nothing for them to do).

## Why

A non-homogeneous `sp` under `constraints Transformation` +
`integrator LoadControl` forces the **first constitutive evaluation of every
increment inside `newStep`**, with the driven face advanced by the full
increment and the interior frozen at the last commit — an overstrain of order
`L/h` in one element layer
([TransformationConstraintHandler.cpp:518-521](../SRC/analysis/handler/TransformationConstraintHandler.cpp)
pre-updates every constrained element;
[LoadControl.cpp:130](../SRC/analysis/integrator/LoadControl.cpp) has nothing
in front of `applyLoadDomain`). Path-independent laws recover in 2–3
iterations; a path-dependent law spuriously yields, its consistent tangent
collapses (`2G·H/(3G+H)` ≈ 1 % of elastic at `H/E` = 1 %), and the active set
oscillates. Measured on the Cerro Lindo rung-4 fuse, all runs pinned to
`KrylovNewton`: **×28.6 Newton iterations and 52 cutbacks under prescribed
displacement, ×0.95 under traction, identical answer to 0.00 %**.

The SP elimination itself is coherent — the constrained DOF has no equation
(`setID(dof,-1)`), so pre-updating the elements is the only way the
prescribed motion reaches the RHS. **The missing piece is a static
predictor**, and the correct ordering already exists in the tree:
`DisplacementControl::newStep` does `incrDisp → applyLoadDomain →
updateDomain`
([DisplacementControl.cpp:287-289](../SRC/analysis/integrator/DisplacementControl.cpp)).
Abaqus defaults to 100 % linear extrapolation of the previous increment;
Kratos ships `use_old_stiffness_in_first_iteration` and enables it precisely
for its J2 / plastic-damage / fatigue tests. One principle, three codes: *the
first linear solve of an increment distributes the boundary motion; the
constitutive law should not be evaluated before that distribution happens.*

## Decisions

- **D1 — the handler is not the defect.** No new SP handler, no edit to
  `TransformationConstraintHandler`'s algorithm. The fix is a predictor in
  the static integrator.
- **D2 — fork class, not a vanilla edit.** Per the fork rule (never edit an
  upstream integrator's algorithm), S1 is a new `LadrunoLoadControl` class
  mirroring how `LadrunoArcLength` relates to stock `ArcLength`. Stock
  `LoadControl.cpp` stays vanilla.
- **D3 — Abaqus-shaped guards, conservative by default.** Extrapolation is
  OFF: on the first increment of the analysis; when
  `Δλ/Δλ_prev ∉ [0.5, 2.0]`; on the first `newStep` after a failed step
  (cutback); and after `domainChanged` or `recvSelf` (state invalidated).
  The predictor moves only the starting guess — it **cannot change a
  converged answer** (trial states are rebuilt from committed state each
  iteration).
- **D4 — gates before C++.** Phase 0 (below) runs the findings' §5 gates
  first. ~~Gate 3 (`Newton -initial` discrimination) **bounds what a predictor
  can buy** — if the collapsed-tangent half dominates, the expected win
  shrinks~~ and this ADR must say so rather than tune toward the preferred
  conclusion.
  > ⚠ **CORRECTED 2026-08-04 (P0 complete).** The struck clause was wrong and
  > is repeated in the G3 gate description below — see that note and
  > [[80b_sp_gates_g1g2g3_2026-08-04]] §G3 for the argument. In short: the
  > collapsed tangent is *downstream of the same overstrain* the predictor
  > removes, so G3's split (measured: tangent 58 %, residual 42 %) describes
  > how the cost decomposes, **not** a cap on the remedy. The decision itself
  > — gates before C++ — stands and was honoured; only its stated rationale
  > was faulty. **P0 is now complete and S1 is GO.**
- **D5 — S2 is fixed only after reproduction.** ✅ **HONOURED AND CLOSED.** The
  [AutoConstraintHandler.cpp:573-584](../SRC/analysis/handler/AutoConstraintHandler.cpp)
  omission was a source-level inference of a **silent wrong answer** (stale
  first residual + `CTestNormDispIncr` with no minimum-iteration guard ⇒
  converged-at-iteration-1). Reproduced at G4, then fixed, then gated by
  `tests/test_auto_handler_sp_update.py` — see
  [[80a_sp_gate_g4_auto_handler_2026-08-04]]. **The inference understated it:
  `EnergyIncr` is fooled too (anything built on `dU` is), and test-dependence
  turned out to *be* the mechanism — the same deck is right or wrong purely by
  choice of convergence test.**
- **D6 — measure, then decide on C/D.** `TransformationFE::getTangForce()`
  (the stub at
  [TransformationFE.cpp:447-453](../SRC/analysis/fe_ele/transformation/TransformationFE.cpp))
  is the hook for the principled Kratos-style route, but it changes behaviour
  for every model with a non-homogeneous `sp` — high blast radius for the
  same benefit S1 gets cheaply. Deferred unless S1's acceptance gate fails.

## S1 design — `LadrunoLoadControl`

### Command

```tcl
integrator LadrunoLoadControl $dLambda <$numIter $minLambda $maxLambda> <-extrapolate $frac>
```

```python
ops.integrator('LadrunoLoadControl', dLambda, numIter, minLambda, maxLambda, '-extrapolate', frac)
```

`frac = 0.0` (default) ⇒ bit-identical to stock `LoadControl` (same Δλ
adaptation by iteration count, same `applyLoadDomain`-only `newStep`).
`frac = 1.0` ⇒ full Abaqus-style linear extrapolation. Crosswalk to document
in the integrators guide: **`EXTRAPOLATION=LINEAR` → `-extrapolate 1.0`**, in
the same style as `*STATIC, RIKS` → `LadrunoArcLength`.

### Predictor mechanics (`newStep`, before `applyLoadDomain`)

```
if (frac > 0 && havePrev && !firstStep && !justCutBack
    && ratio = Δλ/Δλ_prev ∈ [0.5, 2.0]) {
    theModel->incrDisp(frac * ratio * ΔU_prev);   // tangent-consistent guess
    theModel->applyLoadDomain(λ);                 // enforceSPs now pre-updates a layer
                                                  //   whose neighbours already moved
    theModel->updateDomain();                     // physical-strain first evaluation
} else {
    theModel->applyLoadDomain(λ);                 // stock path, byte-identical
}
```

The ordering is copied from `DisplacementControl::newStep` — already proven
in the tree. Cost when active: one `incrDisp` + one extra `updateDomain` per
increment. Cost when `frac = 0`: zero (branch not taken).

### State machine (the part that must not be improvised)

- `update(deltaU)` accumulates `ΔU_accum += deltaU` (the integrator already
  sees every iterative correction).
- `commit()` (override) promotes `ΔU_prev ← ΔU_accum`, `Δλ_prev ← Δλ`,
  sets `havePrev`, clears `justCutBack`, then calls the base `commit()`.
- `newStep()` zeroes `ΔU_accum` at entry. If entered twice with no
  intervening `commit()` (= the caller cut back after a failed
  `analyze`), set `justCutBack` — the failed attempt's partial `ΔU` is
  garbage and Abaqus skips extrapolation after a sharp Δt change anyway.
- `domainChanged()` resizes/zeroes both vectors and clears `havePrev`
  (equation numbering changed — the stored vector is meaningless).
- `sendSelf/recvSelf` transmit **parameters only** (`frac` + stock params);
  the receiving side starts with `havePrev = false`. Predictor state is
  advisory, so losing it on migration is safe by construction.

### Registration footprint (all per `LadrunoArcLength` precedent)

| touch-point | edit |
|---|---|
| `SRC/analysis/integrator/LadrunoLoadControl.{h,cpp}` | **new fork files** |
| `SRC/classTags.h` | `INTEGRATOR_TAGS_LadrunoLoadControl 33015` + band comment |
| `SRC/interpreter/OpenSeesCommands.{cpp,h}` | `OPS_LadrunoLoadControl()` + dispatch (`// Ladruno` marked) |
| `SRC/interpreter/PythonWrapper.cpp`, `TclWrapper.cpp` | expose the command |
| `SRC/analysis/integrator/CMakeLists.txt` | add source |
| `Ladruno_scripts/banner_features.txt` → `patch_banner.py` | banner line (at ship) |

Ledger plan: one `LEDGER_implementations` row (feature, tag 33015, files,
status); `LEDGER_vanilla_files` rows for `classTags.h` / interpreter files
(already-ledgered files get their rows extended); S2 adds a
`LEDGER_vanilla_files` row for `AutoConstraintHandler.cpp`.

## Validation gates

### Phase 0 gates (findings §5 — run on the existing binary, no C++)

All hex8-class (~22.6 k DOF), algorithm pinned, **JSON artifacts saved**
(the precursor diagnostic wrote to stdout only and its figures are now
uncitable — do not repeat that).

- **G1 — overstrain scales as Δδ/h.** Vary only the driven cover's
  through-thickness element size, covers plastic; penalty should fall with
  coarsening, grow with refinement; elastic-cover control at each size.
  Sharpest prediction of the mechanism — a refutation here stops S1.
- **G2 — penalty scales with step size.** 10/20/40 steps, measured per unit
  of *progress*. A prior contrary observation ("halving the step made it
  worse") was taken unpinned under `Newton -initial`; settle it pinned.
- **G3 — which half dominates.** Re-run covers-plastic-displacement pinned to
  `Newton -initial`: substitutes the elastic tangent for the collapsed one
  while keeping the spurious residual. Recovery ⇒ tangent half dominates;
  no recovery ⇒ residual half dominates (Abaqus doctrine predicts the
  latter). ~~**Bounds the predictor's ceiling.**~~
  > ⚠ **CORRECTED 2026-08-04 after running it** — see
  > [[80b_sp_gates_g1g2g3_2026-08-04]] §G3. **G3 does NOT bound the predictor's
  > ceiling, and the original sentence ("a predictor attacks the residual half
  > directly") is wrong.** The two halves are not independent remedies to be
  > split: the tangent half is *downstream of the same overstrain* — the
  > consistent tangent collapses **because** the layer spuriously yields, and it
  > yields **because** it was overstrained. A predictor that removes the
  > overstrain removes the spurious yielding, hence the collapsed tangent too.
  > G3 measures how the cost **decomposes given the overstrain happens**; the
  > predictor's ceiling is the **full** excess over the elastic control.
  >
  > **SCOPE LIMIT (added 2026-08-04, second pass).** That last sentence holds
  > **only where the yielding is entirely spurious** — which is how the G0/G3
  > model is deliberately built (converged σ = 300 MPa against `f_y` = 379.5,
  > nothing yields physically). **In a model that genuinely yields, part of the
  > tangent collapse is real** — the material *is* plastic at the converged
  > state — and no predictor can recover that part. There the predictor's
  > ceiling is the excess attributable to the *spurious* fraction only, which
  > is not what G3 measures and which this ADR has **not** measured. Do not
  > quote "the predictor recovers the full excess" for a plastically-hinging
  > model. The Cerro Lindo case happens to sit on the safe side of this line
  > (the cover never yields, §1), so S1's field target is unaffected.
  > Measured synthetically: tangent 58 %, residual 42 % — so scoping S1 to the
  > residual half would have set a target wrong by >2×.
- **G4 — reproduce the `AutoConstraintHandler` defect.** One brick,
  prescribed face displacement, `constraints Auto` vs `Transformation`,
  `test NormDispIncr`. Expected signature: iteration-1 "convergence" with
  only the boundary layer moved. Gates S2.

### S1 acceptance gates

**Primary (fast, deterministic, no external decks) — added 2026-08-04 after P0.**
Use the synthetic harness `sp_gates/g123_predictor_gates.tcl` as S1's inner loop:

- `-extrapolate 0` reproduces the G0 disp/J2 row **exactly** (60 iterations at
  `N` = 20, 10 steps) — the bit-identical-to-stock requirement.
- `-extrapolate 1.0` drives that row from **60 toward the elastic 20**
  (×3.00 → ×1.00). Judge against the **full** excess, not the 42 % residual
  share — see the G3 correction above.
- Run the acceptance model **in or above** G1's saturation plateau (`N` ≥ 20),
  or the gate understates what the predictor achieved.
- **Do not quote ×28.6 as the target.** That figure is a different model *with
  52 cutbacks*; the synthetic ceiling is ×3.00.

**Field validation (deferred until the Cerro Lindo decks are available):**

- `-extrapolate 0` reproduces the recorded stock march increment for
  increment (δ = 0.0312 → … → 1.0726 mm ladder; 70 increments, 0 cutbacks).
- `-extrapolate 1.0` collapses the covers-plastic displacement case from
  **324 inc / 52 cutbacks / 4 255 iterations** toward the elastic
  **30 / 0 / 149**. *Anything less than a large improvement means the
  mechanism is wrong or incomplete — report that, do not tune.*
- The converged answer does not move (920.0 kN coarse / 1 219.6 kN rated).
- Zone-A pytest green (C++ no-regression gate).

### S2 acceptance gates

- Reproducer demonstrates the wrong answer (or refutes the inference — also
  a valid outcome; record it and close S2 as no-defect).
- Post-fix: `Auto` matches `Transformation` on the reproducer; regression
  test added. **Stays in the fork — not offered upstream** (decided
  2026-08-04; see S2 above). The adjacent `sp -subtractInit` no-op fix on the
  same code path also lives fork-side.

## Implementation plan

| phase | deliverable | gate to proceed | est. size |
|---|---|---|---|
| ~~**P0**~~ **DONE** | G4 → [[80a_sp_gate_g4_auto_handler_2026-08-04]]; G1/G2/G3 → [[80b_sp_gates_g1g2g3_2026-08-04]] (synthetic re-derivation + JSON artifacts) | G1 confirmed (monotone, **saturating**); G2 settled; G3 run **and its stated inference corrected** | done |
| **P1** | S2: reproducer → fix → regression test (own small PR; **fork-only**) | G4 reproduces | ~5 lines C++ + 1 test |
| ~~**P2**~~ **DONE — GATE FAILED** | S1 `LadrunoLoadControl` (tag 33015) + `ladrunoLoadControl` runtime cmd shipped; `-extrapolate 0` bit-identical to stock | **Acceptance gate FAILED: cutbacks 23 → 23, iterations −10 %. Predictor demonstrably fired (armed=63) — a capability failure, not a wiring one.** See [[80c_s1_extrapolate_verdict_2026-08-04]] | done |
| **P3** **← next** | **C (`getTangForce`) — NOW LIVE, D6 triggered by P2's failure**; S3 diagnostics still deferred | P2's gate failed | the principled route: eliminates the overstrain instead of reducing it |

P0 needs the Cerro Lindo fuse decks
(`C:\nmb\My Libraries\Cerro Lindo\Informe No3\Models\Fuse FEM\04_solid_fuse\`)
or a synthetic re-derivation of the same geometry; G4 is self-contained
(one brick).

## Risks and open questions (carried from findings §6)

- **Caller-driven adaptive marches.** The `[0.5, 2]` ratio guard is copied
  from Abaqus's abandonment criterion, but Abaqus owns its stepper and ours
  does not — an aggressive external ladder will trip the guard often and
  silently get stock behaviour. Acceptable (safe direction), but the P2
  report must count how often the predictor actually fired.
- **Contact and `-geom`.** Unknown whether the predictor helps or hurts under
  `LadrunoContact` active-set changes or corotational kinematics — Abaqus
  warns extrapolation can iterate excessively across abrupt BC changes.
  Out of scope here; flagged for the campaign rungs that hit it.
- **Overclaim guard.** Abaqus rates plastic active-set switching "usually
  less severe" than contact, and for submodeling documents the *opposite*
  drive preference (traction is the fragile option near limit loads). Do not
  generalize S1's result to "Neumann is safer for path-dependent materials".
- **Unreproducible precursor figures.** The 895 239 N vs 3.9e-9 N iteration-1
  residual pair predates the corrected mechanism reading and is uncitable;
  regenerate under P0 with artifacts before either number is used.
