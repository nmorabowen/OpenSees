---
title: "ADR 236 — Material history transfer and a commit primitive: restart on a new mesh (WP-236)"
project: Ladruno
type: ADR
status: "Proposed (2026-11-18) — design only, no C++"
priority: medium
owner: nmora
related:
  - "apeGmsh ADR 0136 (Proposed): adaptive mesh convergence. This ADR is the fork half of its Level B."
  - "[[232_ladruno_mpm_general_driver_adr]] §1 and §9 item 9: the FE → MPM handoff waits on this primitive"
  - "[[90_ladruno_viscoplastic_regularization_adr]]: regularization of softening soil is decided there, not here"
  - "[[60a_pfem_crosspollination_amendment]]: apeGmsh is a seeder, not a runtime remesher"
  - "[[79_ladruno_geom_hypo_adr]] §8–§9: the bearing campaign and its unconverged control sequence"
  - "[[157_mortar_friction_pair_state]], [[39_ladruno_contact_domain_adr]]: contact path state"
  - "WP-99 (commit-time refusal seam), WP-104 (`wipe` zeroes process-wide state), WP-168 (`LadrunoCommandTable.h`), WP-192 and WP-224 (threaded commit)"
tags: [adr, material, state-transfer, adaptivity, remeshing, restart, soil, contact, wp-236]
updated: 2026-11-18
---

# ADR 236 — Material history transfer and a commit primitive

> **Lineage.** This is the state-transfer draft that
> [[232_ladruno_mpm_general_driver_adr]] cites by file name
> (`88_ladruno_state_transfer_adaptivity_adr.md`, 2026-09-02, uncommitted in
> the owner's main checkout). Number 88 was taken by the H5DRM higher-order
> work (#778, #800) while the draft sat. It is filed here as 236, the next
> number above every WP in use on 2026-11-18. Its apeGmsh sibling, which the
> draft called "0103", is apeGmsh ADR 0136.
>
> **Checked against** `origin/ladruno` `cfe478254` on 2026-11-18. Every code
> reference below was re-read on that commit and is cited as `path::symbol`.
> Under ADR-87, a claim below that has not been measured is marked as
> proposed.

## 1. Context

### 1.1 Two consumers need the same primitive

- **Mesh adaptivity (apeGmsh ADR 0136).** The loop solves on a mesh,
  estimates the error, builds a new mesh, and either reruns the whole
  history or continues from the state it had. Continuing needs the nodal
  state and the material history written onto a model that did not compute
  them.
- **The FE → MPM handoff (WP-232, phase 9).** ADR 232 §1 says it "needs a
  material state-transfer primitive (`getHistory` / `setHistory`)" and that
  "Nothing has landed."

In both cases the receiving model is a **new model**. The fork does not
remesh a running domain:

- `PartitionedDomain::repartition` prints one line and returns 0.
- `OPS_removeObject` (`SRC/interpreter/OpenSeesMiscCommands.cpp`) has its
  elemental-load cleanup commented out.
- `FE_Datastore::commitState` and `::restoreState` send and receive the
  whole domain by tag, so a restore cannot land on a different topology.
- ADR-60a already settled that apeGmsh seeds models and does not remesh at
  run time.

### 1.2 What exists

| Need | State at `cfe478254` |
|---|---|
| Write nodal state | `OPS_setNodeDisp/Vel/Accel` take `-commit` but set one DOF per call, each starting from the committed vector. `OPS_LadrunoSetNodeTrial` sets full trial vectors and never commits. |
| Make a transient integrator read it | `Newmark::domainChanged` re-reads committed u, v, a. It runs only if `domainChange()` was raised (quirk `WP-232-transient-integrators-ignore-node-writes-between-steps`). |
| Measure the unbalance of an injected state | `OPS_LadrunoTrialResidualNorm` (`SRC/interpreter/OpenSeesCommands.cpp`). It is read-only, needs an active analysis, and its header calls the residual approximate for path-dependent materials. |
| Read material history | No base-class method. `NDMaterial.h` has no `getHistory`, `setHistory` or `getState`. |
| Write material history | Nothing. `LadrunoJ2::setStateToCommitted` is private. `InitStressNDMaterial` injects a stress only. `StagedStrainNDMaterial` captures a birth strain; it cannot be given one. |
| Reach one integration point | `SRC/element/LadrunoElementShell.h::forwardToMaterialPoint` routes `setParameter("material", pt, …)`. `LadrunoBrick::updateParameter` broadcasts to all eight points, so a per-point value does not pass through the element's update. |
| Register a command once | `SRC/interpreter/LadrunoCommandTable.h`, one `LADRUNO_COMMAND` row for all three interpreters (WP-168). |
| Refuse at commit | `SRC/material/LadrunoMaterialStatus.h::ladrunoNoteCommitRefusal`, read in `Domain::commit` (WP-99). |
| Bulk Gauss-point output | `LadrunoRecorder` writes `GLOBAL_GP_COORDS` beside `QUADRATURE/GP_PARAM` and `GP_WEIGHT` under `MODEL_STAGE[*]/MODEL/ELEMENTS/<classTag>-<Class>[g]/`. |

### 1.3 What "the state" is, per material

| Material | Committed history | Note |
|---|---|---|
| `ElasticIsotropic` | none | the identity case |
| `LadrunoJ2` | `epsP_n[6]`, `ebarP_n`, `alpha_n[8][6]`, `D_n`, `cEps22`, plus IMPL-EX increments `dGamma_n`, `dD_n` | increments are not transferred |
| `DruckerPrager` (post-`31322a47a`) | `mEpsilon_n_p`, `mBeta_n`, `mAlpha1_n`, `mAlpha2_n`, `mSigma`, `mElastFlag` | has a protected `Vector getState()`; `materialState` is its stage flip |
| `ManzariDafalias`, `LadrunoSANISAND`, `SAniSandMS` | strain, stress, elastic strain, α, α_in, fabric, void ratio | `const Vector getState()` exists, non-virtual, and is inconsistent across views (quirk `ADR-86-manzaridafalias-getstate-is-shadowed-by-all`); the stage is class-wide static (quirk `WP-232k-sanisand-class-wide-stage-at-particles`) |
| `PM4Sand` | same family, plus `initialize(Vector initStress)` | |
| `LadrunoConcrete3D` | 40 kernel doubles, `cEps22`, and a `trialRefused` flag | a vector, but not in v1 |
| `ASDConcrete3D` | `CrackPlanes` objects and a 6×6 `PT_commit` | not a vector |
| `PressureDependMultiYield` | an array of nested yield surfaces, a PPZ block, and process-static parameter arrays | not in v1 |

The September draft assumed that stress and equivalent plastic strain were
enough. The table is why they are not.

## 2. Decision

### D1 — `NDMaterial::getHistory` and `setHistory`

Three virtuals on the base class (proposed):

```cpp
virtual int getHistorySize(void) const { return -1; }
virtual int getHistory(Vector &h) const { return -1; }
virtual int setHistory(const Vector &h) { return -1; }
```

- **The default refuses.** A return of −1 anywhere in a model aborts the
  command and names the element and the point. A material that cannot
  serialize its history is not a restart material.
- **The layout belongs to the class.** Callers treat it as opaque. Its first
  entry is a layout version, so a stale vector fails on write.
  `LadrunoSANISAND.h` (`LWIRE_TAG`, WP-151) already versions its wire block
  this way.
- **It holds committed, independent history only.** No trial values, no
  tangents. IMPL-EX increments are written as zero, so extrapolation
  restarts.
- **`setHistory` sets the committed state and makes the trial state equal
  to it.**
- **Stress** is recomputed from the injected displacements where the law
  allows it. Where stress is itself state (`mSigma`, the SANISAND family),
  the layout carries it.

**The names avoid three clashes.** `getState()` already exists, non-virtual,
on the SANISAND family and protected on `DruckerPrager`; a base virtual of
that name would be hidden. The parameter word `"materialState"` is already
the stage flip, and `isMaterialPointToken` excludes it. ADR 232 already
calls this primitive `getHistory` / `setHistory`.

**Version 1 materials:** `ElasticIsotropic`, `LadrunoJ2` and
`DruckerPrager`. These are the three ADR 232 §6 admitted first, so the
handoff and the restart share one list. The SANISAND family and `PM4Sand`
follow (S2). `PressureDependMultiYield`, `LadrunoConcrete3D`,
`ASDConcrete3D`, fiber sections and beams are out.

### D2 — Four commands, one row each in `LadrunoCommandTable.h`

| Command | Does |
|---|---|
| `ladrunoGetHistory eleTag ip` | returns the history vector of one point |
| `ladrunoSetHistory eleTag {h_0} {h_1} …` | writes every point of one element |
| `ladrunoSetNodeTrial … -commit` | the existing command, with a flag that also commits the node |
| `ladrunoCommitState` | makes the injected state the committed state |

- **Reaching a point.** The commands resolve the material pointer through a
  fork-side interface on the Ladruno solids, not through `Parameter`
  (`LadrunoBrick::updateParameter` broadcasts, and `Parameter` carries a
  scalar). Stock elements are not in v1. The exact interface is settled in
  S1.
- **`ladrunoCommitState`** calls `Domain::update`, then `Domain::commit`,
  then raises `domainChange()`.
  - Going through `Domain::commit` inherits the WP-99 refusal seam.
  - Raising `domainChange()` is what makes a transient integrator re-read
    u, v and a.
  - `Domain::commit` also fires recorders and advances committed time. The
    command therefore runs before any recorder is declared, and refuses if
    one exists.
- **`ladrunoSetNodeTrial` without `-commit`** stays byte-identical. The
  ADR-52 half-increment gate uses it
  (`tests/test_adr52_w1i1b_halfres_gate.py`).

### D3 — Re-equilibration uses the stage flip, where there is one

The transferred state is not in equilibrium on the new mesh. A zero-load
Newton step restores it, and in doing so return-maps the injected history
and adds plastic strain.

- For a material with an elastic stage (`DruckerPrager` `materialState`,
  the SANISAND family), the balancing step runs in the elastic stage and
  the flip follows it. No new code.
- For `LadrunoJ2`, which has no stage, S3 measures the drift first. A
  freeze flag is built only if the drift fails the gate. If it is built, it
  is read inside `commitState`, written only between steps, and reset by
  `wipe` (WP-104). That is safe under the threaded commit of WP-192 and
  WP-224, where anything written inside commit must be atomic.

The plastic strain added by the balancing step is recorded per point and
reported. It is a gate (INV-3).

### D4 — Contact state is reset, not transferred

On the new mesh the contact surfaces are declared again and the engine
resolves fresh pairs. `FrictionState`, `MortarNormalState` and
`MortarFrictionState` (ADR 157) start from zero, as after
`LadrunoContactDomain::revertToStart`. The normal law is one-sided
(`LadrunoFrictionKernel.h`), so uplift and re-closure need nothing. The
first step pays some extra augmentations to recover the multiplier. Slip
history is lost; transferring it through the keyed maps is a later slice.

### D5 — No regularization is decided here

The September draft proposed a Duvaut–Lions wrapper so that softening soil
could be remeshed. [[90_ladruno_viscoplastic_regularization_adr]] has since
measured that design. Its status reads "D2 REOPENED; close-out PROPOSED at
CP2", its D2 recommends "the generic wrapper not at all", and its §8.1 is
"PROPOSED, not decided". This ADR adds nothing to that question. For a
softening soil, a restart transfers the history it is given and makes no
claim about the band width.

### D6 — The mesh is built on the reference geometry

Remeshing a deformed boundary is ADR 232 §7.5 and phase 10. This ADR does
not need it and does not own it.

### D7 — Bulk history in `.ladruno` (proposed, later slice)

One command per element is enough for S1. A model-sized transfer wants one
file. The route is a new additive group for material history, following the
`ON_PARTICLES` precedent of WP-232a: a minor-version bump, a section in
`ladruno_schema_v1.md`, and a `validate()` arm in `ladruno_format.py`. The
group name and the version number are fixed when the slice opens.

## 3. Invariants

- **INV-1 — the identity round trip.** Same mesh, same deck: read the
  history, start a new process, write it back, commit. The trial residual
  is zero to 1e-10 and the continued load path matches the uninterrupted
  run to 1e-8. A material is admitted only after it passes.
- **INV-2 — no silent partial transfer.** Any −1 aborts, naming the element
  and the point.
- **INV-3 — the balancing step is measured.** The plastic strain or damage
  it adds is reported as a maximum and an L2 norm.
- **INV-4 — the layout is versioned.** A mismatch fails in `setHistory`.
- **INV-5 — `wipe` leaves nothing behind.** Any process-wide state this ADR
  adds is reset in `OPS_clearAllNDMaterial`.
- **INV-6 — the flagless paths do not move.** `ladrunoSetNodeTrial` without
  `-commit`, and every material that does not override the new virtuals,
  give byte-identical results.

## 4. Slices

| # | Content | Gate | Size |
|---|---|---|---|
| S0 | This ADR | Owner sign-off; ADR-87 D4 verdict if the owner treats it as a family ADR | S |
| S1 | The three virtuals; `ElasticIsotropic`, `LadrunoJ2`, `DruckerPrager`; the four commands; the point interface on the Ladruno solids | INV-1 on an elastic block and on a J2 and a DP notched bar past first yield; INV-2 on a model with one material that does not implement; INV-6; the mutation lane fails when `setHistory` is a no-op | M |
| S2 | SANISAND family and `PM4Sand` | INV-1 on a SANISAND deck, with the class-wide stage handled explicitly | M |
| S3 | The balancing step: stage flip on DP, measured drift on J2 from a deliberately perturbed state | Both numbers in this ADR; the freeze flag built or declined, with the reason | S |
| S4 | Bulk history group in `.ladruno` (D7) | The oracle and the format validator pass; apeGmsh reads it | M |
| S5 | Footing restart with a contact reset: monotonic moment–rotation, restarted mid-path on the same mesh | The path continues within the augmentation tolerance; the overturning moment is unchanged to 1e-6 | M |
| Later | Contact slip history; `LadrunoConcrete3D`; nested-surface and crack-plane state | per feature | L |

S1 is also what ADR 232's phase 9 waits for.

**The footing benchmark is a target, not an oracle.** ADR 79 §9 gives the
mild-cone control as 1.0020, 0.9517 and 0.9260 of the exact
Prandtl–Reissner load at 4, 8 and 16 elements across the footing. That is
a sequence still moving by 2.6 % at its last step, and ADR 232 calls it "a
monotone sequence, not a converged value". The series at the real cone
terminates earlier at every refinement. Closing the control sequence is
work for apeGmsh ADR 0136's Level A. It is not a gate here.

## 5. What this ADR does not do

- It does not remesh a running domain, repartition, or add a background
  grid.
- It does not regularize anything (D5).
- It does not transfer fiber-section or beam state.
- It does not estimate error. That is client-side in apeGmsh, on the
  Gauss-point data the recorder already writes.

## 6. Open questions for the owner

1. **Is S1 worth building for adaptivity alone**, or only when the MPM
   handoff (ADR 232 phase 9) needs it? One implementation serves both.
2. **Should the virtuals go on `NDMaterial`** (an upstream file, a full
   rebuild, one `vanilla/` ledger row) or on a fork-side mix-in that the
   admitted materials inherit? The proposal is the base class: the MPM
   handoff reaches materials through `NDMaterial*`.
3. **Does `ladrunoCommitState` refuse when recorders exist**, or bypass
   them? The proposal is to refuse.
4. **Is this a family ADR under ADR-87 D4**, needing an out-of-family
   verdict before acceptance?

## 7. Review log

- **2026-09-02.** A red and a blue team reviewed the brainstorm across both
  repositories. Falsified: stress and equivalent plastic strain as the
  state; "re-equilibrate at frozen load" as a defined operation; mutating
  the running domain. Confirmed: the commit-boundary cadence of ADR-60a and
  nearest-point copy as the default transfer on the apeGmsh side.
- **2026-11-18, before filing.** Every reference was re-read on
  `cfe478254`. Changes from the September text:
  - the number (88 → 236) and the sibling (0103 → 0136);
  - D5 replaces a proposed Duvaut–Lions wrapper with a pointer to ADR 90;
  - the "1.0020 oracle" is restated as an unconverged sequence;
  - `LadrunoConcrete3D` has 41 committed values, not 38, and its relaxation
    is a linear blend `dt/(eta+dt)`, not an exponential;
  - a third name clash (`DruckerPrager::getState`) and the SANISAND
    class-wide stage are recorded;
  - `MortarFrictionState` joins the contact state that is reset;
  - commands register through `LadrunoCommandTable.h`;
  - the freeze flag is demoted from a decision to a measured fallback (D3).

## 8. Implementation log

Nothing built.
