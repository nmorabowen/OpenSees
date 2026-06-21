---
title: "ADR-30 P5 — ExplicitBathe / LNVD adoption of LadrunoProjectionHandler"
project: Ladruno
type: design
status: design (pre-code gate)
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"
  - "[[projection_handler_handoff]]"
  - "[[_adr30_loop_state]]"
tags: [adr, constraints, explicit-dynamics, projection, noh-bathe, sms]
---

# P5 — ExplicitBathe / LNVD adoption of the projection consumer

## Goal
Extend M-orthogonal acceleration projection (the `LadrunoProjectionHandler`, classTag
33001) from the single `CentralDifferenceLadruno` integrator to the **Noh–Bathe** explicit
family, so `constraints LadrunoProjection` works under those integrators too. Backlog item
#2 in the handoff.

## Scope — which classes actually need code

The consumer is a 1-method interface (`LadrunoProjectionConsumer::setConstraintProjector`).
The handler `dynamic_cast`s the integrator to it. So adoption = "implement the consumer +
place the leap-frog hooks." Inheritance map (verified from the headers):

| class | base | overrides newStep/update? | action |
|---|---|---|---|
| `CentralDifferenceLadruno` | TransientIntegrator | — | **already a consumer** (P0–P3) |
| `CentralDifferenceSMS` | CentralDifferenceLadruno | no (only domainChanged) | **inherits projection — already works** |
| `CentralDifferenceSMSConsistent` | CentralDifferenceLadruno | no (domainChanged + refineAccel) | **inherits — already works** |
| `ExplicitBathe` | TransientIntegrator | yes | **ADD consumer + hooks** |
| `ExplicitBatheSMS` | ExplicitBathe | no (domainChanged + send/recv) | inherits |
| `ExplicitBatheSMSConsistent` | ExplicitBathe | no (domainChanged + refineAccel) | inherits |
| `ExplicitBatheLNVD` | TransientIntegrator | yes | **ADD consumer + hooks** |
| `ExplicitBatheLNVDSMS` | ExplicitBatheLNVD | no (domainChanged) | inherits |
| `ExplicitBatheLNVDSMSConsistent` | ExplicitBatheLNVD | no (domainChanged + refineAccel) | inherits |

**Net: only TWO base classes get edited — `ExplicitBathe` and `ExplicitBatheLNVD`.** The
4 SMS/Consistent subclasses inherit projection for free, and all SMS `domainChanged()`
overrides already chain to their base `domainChanged()` (so the IC check propagates). The
entire CentralDifference family already had projection via inheritance.

## The Noh–Bathe step (both bases are structurally identical; LNVD = EB + FLAC damping)

One time step = two sub-steps, with TWO diagonal solves:
- **newStep()**: predict `U_tpdt, V_fake` at `t+p·dt`; setResponse; updateDomain. The
  solution algorithm then forms M, formUnbalance, solves → `A_tpdt`.
- **update(U)** (called once; `U == A_tpdt`):
  1. `A_tpdt = U` (+ optional `refineAccel` for Consistent);
  2. compute `V_tpdt`, predict `U_tdt` at `t+dt`; setResponse; updateDomain;
  3. **second solve done manually inside update()**: `formUnbalance(); solve(); A_tdt = X`
     (+ optional `refineAccel`);
  4. compute `V_tdt`; final setResponse.
- **commit()**: `A_t = A_tdt` (carries to next step's predictor → manifold invariance).

## Hook placement (mirrors the CDL contract — project AFTER refine, before the accel is used)

1. **`setConstraintProjector(p)`** — store `theProjector=p`; `massBuilt=false`.
2. **`domainChanged()`** — after seeding committed state into `U_t/V_t/A_t`: set
   `massBuilt=false`; if projector, `enforceIC(*U_t,*V_t)` (no mass needed; uses L+delta).
   Abort on IC violation (`-projectICs` snaps instead).
3. **`newStep()` first build** — gated on `theProjector && !massBuilt` (NOT `firstStep`, so it
   re-arms after every domainChanged for staged construction):
   `formTangent(CURRENT_TANGENT)` to assemble raw diag(M) into the Diagonal SOE →
   `buildMass(theLinSOE)` (must read M before the algorithm's solve factors it to 1/m) →
   `massBuilt=true` → `project(*A_t)` (the committed a0 onto the manifold). The algorithm
   re-forms M before its own solve, so the extra assembly is harmless (cost: 1 mass assembly
   on the first step of each stage).
4. **`update()` sub-step 1** — after `*A_tpdt=U` (+refine): `project(*A_tpdt)`.
5. **`update()` sub-step 2** — after `*A_tdt=getX()` (+refine): `project(*A_tdt)`.
6. **`commit()`** — scatter `tieForceAtEqn` onto nodes (guarded by `isMassReady()`), so the
   `constraintTieForce` recorder (P4a) works under these integrators too. The cached tie
   force is from the LAST `project()` of the step (= sub-step 2, `A_tdt`).

`A_tpdt`/`A_tdt` are member Vectors → projected in place, no `Aproj` scratch needed (unlike
CDL whose `update(const Vector&)` had to copy U first).

## Why this is correct (reuses the proven P0–P3 algebra)
- The projector is general-C and integrator-agnostic (Gate-A). It only needs (a) diag(M) once,
  (b) `project()` called on every solved accel before that accel advances state, (c) the
  carried accel (`A_t`) to be the projected one. All three are satisfied above.
- SMS lumped injection happens at domainChanged → buildMass (at first newStep) reads the
  SCALED mass → projection is M-orthogonal w.r.t. the actual integration mass. Correct.
- Consistent variants: project with the lumped diag(M) AFTER `refineAccel`, exactly as
  `CentralDifferenceSMSConsistent` already does (shipped). M-orthogonality is w.r.t. the
  lumped mass by convention.
- Prescribed motion (P4b/P4c) is handler-side (`applyLoad`) and integrator-agnostic — it
  already works with any integrator whose `setResponse` skips `eqn<0` DOFs (these do).

## Open questions for the pre-code gate
- **Q1.** Is `formTangent` at the top of the first newStep safe before the predictor mutates
  `U_tpdt`? (Domain is at committed state post-domainChanged; M is state-independent. Expect
  YES, matches CDL's starter forming M at committed config.)
- **Q2.** Tie-force semantics: caching from sub-step 2 only (not a p-weighted blend of both
  sub-steps). Acceptable as "the step's constraint reaction"? CDL caches its single update's
  projection; this is the analogue. Expect YES.
- **Q3.** Does any SMS `domainChanged` override run BEFORE seeding state such that the base
  `enforceIC` sees stale `U_t/V_t`? (They chain `Base::domainChanged()` FIRST, then inject
  mass — so enforceIC runs on freshly-seeded committed state. Expect SAFE.)

## Tests — tests/test_adr30_projection_p5.py
- **TP5-1** ExplicitBathe + equalDOF: momentum conservation + tie exact vs CDL.
- **TP5-2** ExplicitBathe vs Transformation+FullGeneral reference (free-DOF response match).
- **TP5-3** ExplicitBatheLNVD + rigidDiaphragm: slave rz tie exact.
- **TP5-4** ExplicitBatheSMS (inherited projection) + equalDOF: ties hold under mass scaling.
- **TP5-5** ExplicitBatheLNVDSMSConsistent (inherited + refineAccel) + equalDOF: ties hold.
- **TP5-6** no-projection regression: plain ExplicitBathe (no LadrunoProjection) byte-unchanged.
- **TP5-7** constraintTieForce recorder readback under ExplicitBathe.

## Ledger / bookkeeping
- No new classTag (editing existing fork integrators).
- `LEDGER_implementations.md`: note projection-adoption on the EB/LNVD/SMS rows.
- `LEDGER_vanilla_files.md`: none (all six are fork-authored).
- Banner: no new feature line (projection already listed).
