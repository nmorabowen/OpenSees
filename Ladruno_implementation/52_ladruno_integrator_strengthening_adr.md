---
title: Integrator Strengthening — fork-scoped roadmap (implicit + explicit)
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - integrator
  - explicit-dynamics
  - adaptive-time-stepping
  - bulk-viscosity
  - roadmap
aliases:
  - ADR-52
  - Integrator Strengthening
updated: 2026-06-23
parent: "[[49_ladruno_integrator_study_workflow_adr]]"
---

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

# ADR-52 — Integrator strengthening (fork-scoped)

> **Rev 2 (2026-06-23):** revised after a multi-agent adversarial review
> (27 agents, 21/22 findings confirmed against source). The strategic
> direction was upheld; this revision corrects the class-tag band, the W3-I2
> footprint (it *does* need a ledgered header edit), the registration
> touch-points, the W1-E2 footprint bucket, the bulk-viscosity energy
> bookkeeping, and several provenance/numerics overstatements. See
> §"Adversarial-review corrections" for the audit trail.

## What

A roadmap to strengthen OpenSees time/path integration in **both lanes** —
implicit and explicit — acting **on fork-owned code wherever possible**, and
paying only small, *ledgered* vanilla edits where unavoidable. It turns the gap
verdict of [[49a_integrator_scorecard_2026-06-23]] into a sequenced program.

The guiding principle: **stop adding point-solution integrator classes; add
reusable *capability layers* (adaptivity, stabilization, sensitivity) and
consolidate fork-owned classes.** Spend new effort where the fork actually
lives — robust nonlinear statics, contact, and explicit quasi-statics.

**In scope:** an error-gated adaptive-Δt layer hosted in the (fork-owned,
Python) robust-solve driver; explicit bulk viscosity on the 3 fork continuum
elements; collapsing the fork's own `ExplicitBathe*` class explosion;
SMS reporting + parallel-safety fixes; (deferred) sensitivity-carrying
subclasses of HHT / generalized-α; (gated) a tunable implicit composite built
on TRBDF2.

**Not in scope:** editing the *algorithm* of any upstream integrator;
deleting/merging upstream classes; hourglass control (element technology —
separate ADR); explicit adaptive Δt (out of scope by design — see §"Explicitly
dropped"). The cleanup items become **documentation**, not code.

## Why

The [[49_ladruno_integrator_study_workflow_adr]] study (run
[[49a_integrator_scorecard_2026-06-23]]) found:

- One real capability gap vs Abaqus/LS-DYNA: **error-gated automatic time
  stepping**. (OpenSees *does* ship iteration-count variable stepping —
  `VariableTimeStepDirectIntegrationAnalysis` via `analyze numIncr dt dtMin
  dtMax Jd` — so the gap is the *accuracy/error* gate, not adaptivity per se;
  49a:191 overstated this as "none built-in".)
- A second gap: **explicit bulk viscosity** (b1/b2) — absent everywhere in
  OpenSees incl. the fork.
- A fork-owned maintainability problem: **6 `ExplicitBathe*` classTags** from
  parallel sibling classes, with untested cross-layer combinations.
- A sensitivity blind spot: DDM exists only on `Newmark`/`DisplacementControl`/
  `ArcLength`/`MinUnbalDispNorm` — **not** on the modern dampers `HHT` /
  `GeneralizedAlpha`, nor any explicit integrator.

The fork constraint (CLAUDE.md): **do not touch vanilla OpenSees if avoidable.**
Operationally:

1. **Never edit an upstream integrator's *algorithm*** (`.cpp` math).
2. **Minimize new integrator classes** — prefer extending the fork-owned
   robust-solve driver and refactoring fork-owned classes.

A *new C++ integrator class* costs a fixed set of small, ledgered hooks (see
§Where for the verified list). The robust-solve driver is
**`Ladruno_scripts/robust_drive.py`** (OpenSeesPy, not C++) — so the adaptive-Δt
work lands with **zero C++/vanilla footprint** (the convergence-driven tier
re-hosts an existing built-in; only the error-gated tier may need a small
fork-owned read-only helper).

## Where

**Verified registration touch-points** for a *new* C++ integrator (grep-checked
during the review — corrects Rev 1's list):

- `SRC/tcl/commands.cpp` (Tcl factory strcmp chain, ~:231-239, :5630-5686)
- `SRC/interpreter/OpenSeesCommands.cpp` (+`.h:603-611`) (Python factory chain)
- `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` (broker switch)
- `SRC/runtime/runtime/TclPackageClassBroker.cpp` (second broker)
- `SRC/classTags.h` (≥33000 ladruno band)
- `SRC/analysis/integrator/CMakeLists.txt` + `Makefile`
- *Not* per-integrator hooks: `DirectIntegrationAnalysis.cpp`, `TclWrapper.cpp`,
  `PythonWrapper.cpp`/`PythonModule.cpp` (generic one-time wiring / banner only).

Per wave:

- **Modify (fork-owned, zero vanilla):**
  - `Ladruno_scripts/robust_drive.py` — adaptive-Δt layer (W1-I1).
  - `SRC/analysis/integrator/CentralDifferenceSMS*.{h,cpp}`,
    `ExplicitBatheSMS*.{h,cpp}` — Δt_cr reporting + parallel audit (W1-E3).
- **Modify (fork-owned classes, but *non-zero* ledgered vanilla — registration
  deletions):**
  - `SRC/analysis/integrator/ExplicitBathe*.{h,cpp}` (6 fork classes) — collapse
    to one + flags (W1-E2). Touches the registration files above (deletions) →
    ledger in [[LEDGER_vanilla_files]].
- **New code (fork-owned):**
  - Bulk-viscosity on the 3 fork continuum elements:
    `SRC/element/ladrunoBrick/`, `SRC/element/ladrunoPlane/` (LadrunoBrick,
    LadrunoQuad, LadrunoCST) + a flag on the fork explicit integrators (W2-E1).
  - `LadrunoHHT`, `LadrunoGeneralizedAlpha` (sensitivity subclasses, W3-I2),
    `LadrunoBatheImplicit` *as a subclass of `TRBDF2`* (gated, W3-I3).
- **Vanilla edits (accepted, ledgered):**
  - W1-E2 registration deletions (above).
  - **W3-I2 base-header promotion:** `HHT.h` and `GeneralizedAlpha.h` declare all
    members `private:`; sensitivity subclassing is impossible without promoting
    the needed members to `protected:` (no `.cpp`/algorithm change). Mark
    `// Ladruno: protected for sensitivity subclass`; record in
    [[LEDGER_vanilla_files]]. (Newmark works only because `Newmark.h` is already
    `protected`.)
- **Reference / similar fork impls:** `CentralDifferenceLadruno.{h,cpp}`
  (CFL/projector/first-step pattern), `CriticalTimeStep.{h,cpp}`,
  `LadrunoMassScaling.{h,cpp}`, `LadrunoBrick.cpp:1659-1665` (element c_d/L_e/
  ε̇_vol for viscous hourglass — the pattern W2-E1 reuses),
  `Newmark.{h,cpp}` (the 5-method sensitivity seam), [[ladruno_integrators_guide]],
  [[31_ladruno_robust_solve_driver_adr]].
- **Class-tag reservations:** integrator band **occupied 33000–33012**
  (`classTags.h:1143-1155`). Next free → `LadrunoBatheImplicit` **33013**,
  `LadrunoHHT` **33014**, `LadrunoGeneralizedAlpha` **33015**. Assign from the
  top of the band at build time; **do NOT pre-reserve tags W1-E2 may free**
  (33002/33009/33011) — re-check `classTags.h` immediately before each PR.

## How

Sequenced by vanilla footprint × value. Each item is a self-contained sub-task.

### Wave 1 — lowest footprint (do first)

**W1-I1 — Adaptive implicit Δt (Python; convergence tier re-hosts a built-in,
error tier is net-new).**
Add a step-control layer to `robust_drive.py` (today a *statics* driver:
LoadControl/DisplacementControl + a quasi-static DR rung). Two tiers:
- *W1-I1a — convergence/cost-driven (now, zero C++):* a Python re-host of the
  existing `VariableTimeStepDirectIntegrationAnalysis` (`analyze numIncr dt
  dtMin dtMax Jd`) so the robust driver can adapt Δt across the transient
  integrators. **Caveat:** iteration-count control is a *cost*/convergence
  criterion, **not an accuracy** one — it can grow Δt on a stable-but-inaccurate
  stiff response. This is a new transient lane bolted onto the statics spine;
  it does **not** by itself close the Abaqus/LS-DYNA accuracy gap.
- *W1-I1b — error-gated (net-new):* an Abaqus-style **half-increment residual**
  accuracy gate. Not exposed by today's converged-step getters; if not faithfully
  buildable in Python (`printB`/`setNode*` approximations aside), add a
  *fork-owned read-only* residual helper rather than editing the analysis loop.
  **Abaqus accuracy parity is reached only once W1-I1b lands.**
- *Test:* gate the adaptive trajectory against the `DisplacementControl
  [reference]` path-follower already in `torture_snapthrough.py`/
  `torture_softening.py` within a stated error bound; `torture_*.py` as
  regression + wall-clock baseline.

**W1-E2 — Collapse `ExplicitBathe*` 6 → 1 + flags.**
Fold `ExplicitBathe{,LNVD,SMS,SMSConsistent,LNVDSMS,LNVDSMSConsistent}` into a
single `ExplicitBathe` with orthogonal flags `-lnvd`, `-sms`, `-consistent`.
- *Structure (corrected):* these are **sibling** classes (all `: public
  TransientIntegrator`), not an inheritance chain. `-sms`/`-consistent` select
  the mass-scaling path in `domainChanged` (`ExplicitBatheSMS.cpp:150-154`);
  `-lnvd` routes through a **unified `formUnbalance()` override** on the hot path
  (`ExplicitBatheLNVD.cpp:558-570`) plus `alpha_flac`/`addLocalDamping()` state.
- *Footprint (corrected):* **net-negative line count, but NON-zero vanilla** —
  it *deletes* registration across the standard files (rebase-fragile, ledger
  it). The net-negative win is on **command/wrapper entries**, not the broker.
- *Serialization (corrected):* each class bakes its tag into the base ctor, so
  `getClassTag()`/`sendSelf` use e.g. 33009. **Keep the 6 broker cases** (route
  each to `new ExplicitBathe()` and extend `sendSelf/recvSelf` to round-trip the
  `{lnvd,sms,consistent}` flag set decoded from the incoming tag) so saved-DB /
  parallel `recvSelf` of retired tags still reconstructs (precedent:
  [[LEDGER_vanilla_files]] "no type for class tag 33000/33001").
- *Aliasing (corrected):* there is **no alias facility** — selection is a
  hand-written strcmp chain (`OpenSeesCommands.cpp:1878-1903`, `commands.cpp:
  5630-5686`). During the one-release deprecation window the ~12 dispatch
  branches are **retained** (each forced to `OPS_ExplicitBathe()` with fixed
  flags); the command/broker surface shrinks only after the window closes. The
  class `.cpp/.h` are removable in-PR.
- Then add the cross-layer validation that is missing today (LNVD × consistent
  PCG). Keep classTag 33000.

**W1-E3 — SMS honesty + parallel audit (fork-internal, zero vanilla).**
(a) When SMS is active, report Δt_cr as a *"pre-scaling estimate"* (still
computed at `CentralDifferenceSMS.cpp:166`) instead of silently rejecting
`-cflAbort`/`-recompute` (MF-1). (b) Resolve the consistent-SMS parallel
contradiction (`CentralDifferenceSMSConsistent` "not parallel-safe" vs
`ExplicitBatheSMSConsistent` claims distributed PCG) — audit
`LadrunoConsistentRefine` and document the real story.

### Wave 2 — fork-only, scoped

**W2-E1 — Explicit bulk viscosity (b1 linear / b2 quadratic).**
The standard explicit stabilizer absent everywhere in OpenSees. Element-level
pressure term:
- linear `p₁ = b1·ρ·c_d·L_e·ε̇_vol` (damps highest element dilatational mode),
- quadratic `p₂ = ρ·(b2·L_e)²·|ε̇_vol|·ε̇_vol` (**compression only**; smears shocks),
- added to element internal pressure, **excluded from `getStress()`**.
- *Provenance (corrected):* `c_d=√(D₀₀/ρ)`, `L_e`, `ε̇_vol` are computed **in the
  element** (the pattern `LadrunoBrick.cpp:1659-1665` already uses for viscous
  hourglass). `CriticalTimeStep` returns only λ_max/Δt scalars — **it is not a
  data source** for this.
- *Scope (corrected):* available only on the **3 fork continuum elements**
  (LadrunoBrick, LadrunoQuad, LadrunoCST). **Not** on vanilla stdBrick/SSPbrick/
  quad — document this availability gap prominently.
- *Δt_cr (corrected):* adds only a **bounded** ξ≈b1 fraction, reducing Δt_cr by
  `(√(1+ξ²)−ξ)` — NOT the unbounded `betaK·ω_max/2` collapse. (Any ξ>0 *does*
  reduce the realized damped step; "does not collapse Δt_cr" is true only vs
  unbounded betaK.)
- *Energy (corrected — see MUST-FIX #4):* bulk viscosity is non-conservative
  (Abaqus ALLVD). **Decision:** include the viscous pressure in
  `getResistingForce()` so internal energy IE closes the fork balance
  `RES = ULW−(KE+IE+DW)` (`EnergyBalanceKernel.h:80-138`), while **excluding** it
  from `getStress()`. Ship an **energy-closure / ALLVD gate** as a required test.
- *Defaults:* Abaqus `b1=0.06, b2=1.2`.
- *Relation to contact:* this damps the highest *element dilatational* mode and
  is **complementary to**, not a substitute for, the already-shipped contact
  `-visc` stabilizer (#385/#387).
- *Tests:* numpy FD operator oracle (~1e-10); compression-only branch test;
  energy-closure gate (≤0.01%); Zone-A byte-identity on the b1=b2=0 default
  path; 1-D bar shock — verify the **undamped** 2/ω estimate is unchanged AND
  the damped Δt_cr drops by ≤(1−b1)≈6% at b1=0.06; reported stress excludes the
  viscous pressure.

### Wave 3 — standard new-class cost (deferred)

**W3-I2 — Sensitivity-carrying `LadrunoHHT` / `LadrunoGeneralizedAlpha`.**
Subclass `HHT`/`GeneralizedAlpha` and add DDM. **Footprint (corrected):** this
needs **2 ledgered base-header edits** — promote the required members in
`HHT.h`/`GeneralizedAlpha.h` from `private:` to `protected:` (no algorithm
change) — *plus* the standard registration hooks.
The sensitivity seam is **five virtual overrides, not three** (the Newmark
pattern): `formSensitivityRHS`, `saveSensitivity`, `commitSensitivity`, **and**
`formEleResidual`/`formNodUnbalance` branching on `sensitivityFlag`
(`Newmark.cpp:577-747`) — re-deriving the α-weighted M/C sensitivity terms (the
base `TransientIntegrator` paths have no sensitivity branch). Unblocks
reliability/fragility/FORM on the numerically-damped integrators.
*Tests:* FD gradient-check oracle (~1e-6); Zone-A byte-identity on the
sensitivity-off path.

**W3-I3 — Tunable implicit composite (Bathe-Noh 2012) — GATED → NO-GO (2026-06-24).**
`TRBDF2` **already IS** the Bathe (2007) Trap+BDF2 composite
(`TRBDF2.cpp:29-30`), with constants hard-coded (`:114-144`). The gap would be the
**user-tunable γ/ρ∞ variant (Bathe & Noh 2012)** — a fork subclass of `TRBDF2`
overriding the constant-setup block.

**Gate result — DO NOT BUILD (benchmark `Ladruno_scripts/w3i3_bathe_gate_benchmark.py`).**
Benchmarked the existing integrators as proxies (fixed `TRBDF2` vs tunable monolithic
`HHT(α)` + trapezoidal) on (1) a mass→ENT-stop impact, dt swept well-resolved→2·Tc;
(2) a rigid-limit stiffness stress k=1e6→1e12 at fixed dt; (3) a 1-D bar step-load wave
at Courant 1 and 4. Findings:
- **No robustness gap.** Every scheme converged everywhere, including the rigid-impact
  limit (k=1e12, ~3e-3 steps/contact). The composite's headline robustness edge never
  manifested as an HHT-fails-where-TRBDF2-survives case on these contact problems.
- **Tunability already covered.** Trapezoidal (ρ∞=1) injects spurious energy on under-
  resolved impact (E_retain → 7–9×); `HHT(α)` dials that down monotonically and can be
  tuned to near-ideal energy at each dt — i.e. the fork already spans ρ∞.
- **The composite's real edge is MAX dissipation, which `TRBDF2` already delivers.** On
  the under-resolved wave (Courant 4) `TRBDF2` suppressed wavefront ringing best
  (ripple 0.99 vs HHT's ~1.20 at any α, trapezoidal 1.30). But that is the ρ∞=0 corner —
  a *tunable* composite (ρ∞>0 = less dissipation) would be *worse* there, not better.
- So the "composite + partial dissipation" quadrant a Bathe-Noh variant would occupy is
  **not demonstrated to be needed**: when you want the composite you want its full
  dissipation (= `TRBDF2`, set-and-forget); when you want partial/tunable dissipation
  `HHT`/`GeneralizedAlpha` suffice with no robustness penalty.

**Guidance (→ `ladruno_integrators_guide`), not code:** use `TRBDF2` for set-and-forget
maximal high-frequency dissipation on stiff/contact/wave problems; use `HHT`/
`GeneralizedAlpha` for dialed/partial dissipation. **Revisit trigger:** a concrete stiff
multi-DOF contact problem where `HHT`/`GeneralizedAlpha` fail to converge but a composite
survives *and* full `TRBDF2` dissipation is too much — this gate found no such case.
(Caveat: the benchmark is SDOF impact + a 1-D bar; it did not exercise large 3-D contact
with material softening, the composite's most-favorable regime.)

### Explicitly dropped (documentation, not code)

These would require editing vanilla algorithm files → guidance in
[[ladruno_integrators_guide]] / [[LEDGER_quirks]], **not** edits:

- ~~Merge `Newmark1` → `Newmark`~~ → document "known dup, prefer `Newmark`".
- ~~Merge `ArcLength1` → `ArcLength`~~ → document "use `ArcLength`/`LadrunoArcLength`".
- ~~Fix `ExplicitDifference` / remove `CentralDifference{Alternative,NoDamping}`~~
  → document "avoid; use `CentralDifferenceLadruno` / `ExplicitBathe`".
  `ExplicitDifference` is a **lower-quality duplicate** of
  `CentralDifferenceLadruno` (no CFL abort; crude `Ut₋₁=Ut` startup
  `:265`; permissive `updateCount>2` `:274`; non-standard output reconstruction
  `:303-310`) — same core leap-frog, so the gap is code quality, not capability.
  A genuine fix could be **upstreamed as a real-OpenSees PR**, separate from the fork.
- ~~Validate/clean `EQPath`, `HSConstraint`~~ → mark experimental in docs.
- **Explicit adaptive Δt** — out of scope **by design**: SMS is the competing
  strategy under explicit (fixed-target mass scaling; SMS hard-rejects
  `-cflAbort`/`-recompute` at `CentralDifferenceSMS.cpp:91-96`). A contact-aware
  Δt controller on `CentralDifferenceLadruno` (non-SMS baseline) is a possible
  **future ADR**, not this one.

### Vanilla-footprint summary (corrected)

| Item | Footprint | Lane |
|---|---|---|
| W1-I1a adaptive Δt (convergence tier, robust_drive.py) | **zero** (Python; re-hosts a built-in) | implicit |
| W1-I1b error-gated tier | **conditional** — zero if Python-buildable, else low (fork read-only helper) | implicit |
| W1-E2 ExplicitBathe collapse | **net-negative line count, NON-zero vanilla** (ledgered registration deletions) | explicit |
| W1-E3 SMS honesty / parallel audit | **zero** (fork) | explicit |
| W2-E1 bulk viscosity | **low** (3 fork elements + integrator flag) | explicit |
| W3-I2 sensitivity subclasses | standard hooks **+ 2 base-header private→protected edits** (ledgered) | implicit |
| W3-I3 implicit Bathe (gated) | standard hooks (subclass of TRBDF2) | implicit |

## Risks / open questions

> [!done] RESOLVED (#407)
> W1-I1b: can the half-increment-residual error gate be computed from OpenSeesPy
> alone? **No** — `setNodeDisp` triggers no `Element::update()`, so the assembled
> residual ignores the displacement-dependent internal force at an injected state.
> Shipped two small read-only fork commands (`ladrunoTrialResidualNorm`,
> `ladrunoSetNodeTrial`); registration-only vanilla, no classTag/header promotion.

> [!question]
> W1-E2: keep the deprecated dispatch branches for one release (no alias
> facility exists — they are retained strcmp branches), or hard-cut and break
> existing fork models/banners?

> [!question]
> W3-I3: subclass `TRBDF2` (override constant-setup for tunable γ/ρ∞) vs a
> from-scratch class — confirm the subclass can reach the needed members.

**Risks surfaced by the adversarial review (decided / to-track):**

1. **Energy-recorder closure for non-conservative element dissipation.** Any
   element dissipative term that is neither Rayleigh `C` (DW) nor in
   `getResistingForce()` (IE) silently breaks `RES = ULW−(KE+IE+DW)`
   (`EnergyBalanceKernel.h:80-138`). *Decided* (MUST-FIX #4): viscous pressure
   goes into `getResistingForce()`, excluded from `getStress()`; ship an ALLVD
   closure gate.
2. **Class-tag band contention with in-flight waves.** Band occupied through
   33012; W1-E2 may free 33002/33009/33011 mid-program. Assign new tags only
   from the confirmed top (≥33013), re-checked against `classTags.h` immediately
   before each PR.
3. **Serialization compatibility of the ExplicitBathe collapse.** Must preserve
   `recvSelf`/DB/parallel reconstruction under all retired tags (keep broker
   cases; extend `sendSelf/recvSelf` to carry the flag set). Risk of silently
   breaking saved models / MPI runs if broker cases are deleted.
4. **W3-I2 forces a base-header edit the constraint discourages.** DDM on HHT/
   GeneralizedAlpha cannot be pure subclassing — it needs a `private→protected`
   header promotion (lightest path) or a full copy-fork. Recorded as an explicit
   exception to rule #1, ledgered.
5. **Bulk-viscosity availability gap.** Off-by-default and confined to 3 fork
   continuum elements → silently unavailable on the vanilla meshes most explicit
   users run. Document prominently.

- Numerical: `b2` (quadratic) must act in compression only, else spurious
  tensile damping; verify the bounded Δt_cr reduction empirically.
- Parallel: W1-E3 must settle whether consistent-SMS PCG is rank-safe before any
  new consistent path is added.
- Banner/ledger: W1-E2 changes shipped feature rows — update
  `Ladruno_scripts/banner_features.txt` + run `patch_banner.py`, and reconcile
  [[LEDGER_implementations]].

## Adversarial-review corrections (audit trail)

Rev 2 folded in 4 MUST-FIX + 15 SHOULD-FIX confirmed findings from the
2026-06-23 multi-agent review (workflow `wf_309e2697-eed`, 21/22 findings
confirmed). Headline corrections: class-tag band → 33013+ (was wrongly 33008+);
W3-I2 needs ledgered header promotion (was "no base edit"); W3-I2 = 5 overrides
(was 3); bulk-viscosity energy bookkeeping decided + gated (was open question);
registration touch-point list verified; W1-E2 footprint = non-zero vanilla (was
"net-negative vanilla"); W1-I1 re-hosts an existing built-in + iteration-count ≠
accuracy; bulk-viscosity provenance is the element, not CriticalTimeStep; TRBDF2
already is Bathe-2007 (gap = tunable 2012 variant). Only refuted finding:
the bulk-viscosity *formula itself* is sound and needed no change.

## Implementation log

Suggested order: W1-E3 → W1-E2 → W1-I1a → W2-E1 → W3-I2 → W3-I3-gate. Each wave =
its own PR on `ladruno`.

- **2026-06-24 — W1-E3 shipped (#394).** SMS dt_cr honesty + parallel-safety doc
  fix. `-cflAbort`/`-recompute` under the 4 SMS integrators
  (`CentralDifferenceSMS{,Consistent}`, `ExplicitBatheSMS{,Consistent}`) now
  downgrade to report-only (the pre-scaling dt_cr is surfaced) instead of refusing
  the run; the stale "consistent variant is NOT parallel-safe" comments were
  corrected after auditing `LadrunoConsistentRefine.h` (the consistent PCG *is*
  parallel-safe — global inner products + shared-DOF assembly, ADR-38 V5). Fork-only,
  no serialization/numerical change. Zone-A green.
- **2026-06-24 — W1-I1a shipped (#396).** Transient adaptive-Δt lane
  (`robust_transient()` + `TransientResult`) added to `Ladruno_scripts/robust_drive.py`:
  a pure-Python re-host of the built-in variable-step transient analysis
  (iteration-count dt sizing) + the algorithm ladder + honest `integrated`
  (not-accuracy-certified) verdict. Zero vanilla footprint. Self-tested live against
  the dist build (linear SDOF + Newmark; static self-test unaffected). The
  accuracy-grade follow-up is **W1-I1b** (half-increment-residual error gate).
- **2026-06-24 — W2-E1 shipped for LadrunoBrick (#399).** Explicit bulk viscosity
  (`-bulkViscosity b1 b2`, alias `-bv`): viscous volumetric artificial-pressure
  stress `s=c_bulk·ε̇_vol` into the resisting force only (excluded from reported
  stress; `c_d` from the initial elastic tangent, `L_e=vol^(1/3)`), dissipative by
  construction, off-by-default (bit-identical). Guarded to `-geom linear` + std/bbar
  (warn+zero otherwise). Passed an 18-agent adversarial review (merge-with-fixes:
  applied corot guard B1, formulation guard B2, initial-tangent S1, hoist S2, doc S5).
  Validated at runtime by a new `zone_a` test (off-path bit-identity + dissipation/sign)
  that CI builds-and-runs. Fork-only.
- **2026-06-24 — W2-E1 completed for LadrunoQuad + LadrunoCST (#403).** Same viscous
  volumetric stress ported to the 2D fork continuum elements: normal comps xx,yy;
  `L_e=getCharacteristicLength()` (Quad √area, CST √(2·area)); 2D is `-geom linear`
  only (no geom guard); Quad wired through std/bbar (SSP/EAS warn+zero at parse), CST
  single-GP std. Off-by-default bit-identical; coeffs threaded through ctors (no
  body-force collision), sendSelf/recvSelf, Print. Runtime `zone_a` test
  (`..._2d.py`) covers off-path identity + dissipation/sign for both, under explicit
  central difference (the runtime gate caught — and rejected — an earlier implicit-
  Newmark test that blew up at large b1). **Bulk viscosity now on all 3 fork continuum
  elements.**
- *W2-E1 follow-ups (deferred, non-blocking):* **S3** `bvDissipated`/ALLVD recorder
  channel (energy balance already closes); **S4** one-time warning when material
  `rho==0`; extend bulk viscosity to the uri/ssp/eas single-point Brick/Quad paths.
- **2026-06-24 — W1-I1b shipped (#407).** Half-increment-residual ACCURACY gate —
  upgrades W1-I1a from convergence-driven to accuracy-grade (the one true gap vs
  Abaqus/LS-DYNA). **Open question RESOLVED: not Python-only.** Probed the dist build:
  `setNodeDisp` sets only the node trial vector and triggers no `Element::update()`, so
  `reactions()`/`printB()` (even after `updateElementDomain`) report a residual that
  ignores the displacement-dependent internal force (`|b|`=0 at an injected state). So
  two small fork commands (registration-only vanilla, **no classTag, no header
  promotion** — much lighter than W3-I2's footprint): **`ladrunoTrialResidualNorm
  <loadTime>`** (drives the element `update()` loop with a forced POSITIVE half-step dt,
  then the active integrator's `formUnbalance()` → inf-norm of the free-DOF dynamic
  unbalance; optional `loadTime` re-applies loads at the midpoint; no commit) and
  **`ladrunoSetNodeTrial`** (full-vector trial setter — the per-dof `setNodeDisp` cannot
  build a multi-dof trial state). `robust_transient(error_gate=True, haftol=…)` builds
  the constant-avg-accel midpoint state, reads the residual, and sizes the next Δt from
  `min(iter-count, (haftol/r_half)^(1/order))`. **FEED-FORWARD** (OpenSees commits on
  success → no mid-step rejection); verdict `accuracy_gated` iff every committed step met
  `haftol`, else `integrated` (`n_overtol`/`halfres_max` carry the evidence). **Fidelity:
  EXACT for rate-/path-independent (elastic) materials; APPROXIMATE for inelastic/rate-
  dependent (post-commit reference is t_{n+1}, representative +dt imposed).** Adversarial
  review (23 agents, 6/17 confirmed) caught a MAJOR the elastic-only test masked: the
  midpoint `applyLoad` left the global `ops_Dt` NEGATIVE (`t_mid−t_{n+1}`), corrupting
  rate-dependent materials' relaxation (`exp(-ops_Dt/tR)`→growing) — fixed by the forced
  positive dt + save/restore. Tests (`tests/test_adr52_w1i1b_halfres_gate.py`, zone_a):
  numpy oracle incl. load-at-midpoint; Newmark-consistency check; full-vector-setter vs
  per-dof loss; committed-state-untouched; gate-improves-accuracy vs fine-Δt ref; Maxwell
  rate-dependent regression (guards the `ops_Dt` fix); gate-off == W1-I1a.
- **2026-06-24 — W3-I3 GATED → NO-GO.** Benchmarked the existing integrators as proxies
  (fixed `TRBDF2` vs tunable `HHT(α)`/trapezoidal) on impact + rigid-limit + 1-D wave
  problems (`Ladruno_scripts/w3i3_bathe_gate_benchmark.py`). No robustness gap (all
  converged to the rigid limit); `HHT` already spans ρ∞; the composite's edge is MAX
  dissipation which `TRBDF2` already gives. A *tunable* composite's quadrant (composite +
  partial dissipation) is not demonstrated to be needed → don't build; documented the
  use-the-right-tool guidance + a revisit trigger (stiff multi-DOF contact where HHT
  diverges but a composite survives). See the W3-I3 section above. No code shipped (the
  benchmark script is the evidence artifact).
- **2026-06-24 — W3-I2 PR1 shipped (#413): `LadrunoHHT`.** Sensitivity-carrying (DDM)
  subclass of `HHT`, classTag **33013** (the lowest free integrator tag — W3-I3's reserved
  33013 was NO-GO #410, never built; assigned from the top of the free band per the
  documented rule). Unblocks reliability/fragility/FORM on the numerically-damped HHT (DDM
  shipped in vanilla only on Newmark/DisplacementControl/ArcLength/MinUnbalDispNorm).
  **Vanilla footprint was bigger than the ADR predicted but stayed header-only:** besides
  the `private:`→`protected:` promotion, `HHT`'s ctors hardcode `INTEGRATOR_TAGS_HHT` (no
  classTag param, unlike `Newmark`), so an extra protected **inline** classTag ctor was
  added to `HHT.h` — keeping `HHT.cpp` byte-identical. Derivation: HHT's `U/Udot/Udotdot`
  follow the Newmark recurrence ⇒ `saveSensitivity`/`commitSensitivity`/`formSensitivityRHS`/
  `formIndependentSensitivityRHS`/`computeSensitivities` are copies of Newmark; only
  `formEleResidual`/`formNodUnbalance` differ — α-weighted damping multiplicator, `∂C/∂h` on
  `Ualphadot`, and an EXTRA element-only `−K·(1−α)·dUₙ` term (via `addK_Force`, consistent
  tangent; vanishes at α=1). FD-vs-DDM oracle (undamped + mass-prop damped, param=E) passed
  in CI on the first run; 22-agent adversarial review (15 raw → 3 confirmed, 0 blockers).
  **New CI gate discovered:** `ci/check_manifest.py` (G9) requires every Ladruno classTag to
  have a row in `Ladruno_implementation/testbed/manifest.yaml`. **Lesson:** `Element` base
  `getTangentStiffSensitivity`/`getMassSensitivity` return zero + warn (betaK·K DDM not
  implemented) and Truss doesn't override them ⇒ the `∂C/∂h`-on-`Ualphadot` term can't be
  FD-tested on a Truss; it's derivation-validated + pinned by the α=1→Newmark reduction.
- **2026-06-24 — W3-I2 PR2 (#415): `LadrunoGeneralizedAlpha`.** DDM subclass of
  `GeneralizedAlpha` (Chung-Hulbert), classTag **33014** — strict superset of `LadrunoHHT`
  with TWO spectral params (αF on K/C at `Ualpha`/`Ualphadot`, αM on M at `Ualphadotdot`).
  Same header-only `GeneralizedAlpha.h` edit (promotion + inline classTag ctor). **Crux
  discovered by the FD oracle + adversarial review:** OpenSees `GeneralizedAlpha`'s tangent
  (`αM·c3·M`) is **inconsistent with its own primal residual**, which integrates inertia at
  the full step `Udotdot` (`update()` sets accel=`Udotdot`) ⇒ effective Jacobian M-coef `c3`.
  A first cut built the DDM tangent-consistent (`αM`, inertia at `Ualphadotdot`) and the
  Zone-A FD-vs-DDM oracle FAILED ~2e-3 at αM=0.9. Fix: build the sensitivity residual
  **primal-consistent** (M at `Udotdot`, no αM — like Newmark) AND **re-form the
  sensitivity-solve tangent with `c3·M`** (a `sensTangentFlag` branch in
  `formEleTangent`/`formNodTangent` + `formTangent()` in `computeSensitivities`) rather than
  reuse the inconsistent factored primal tangent. K/C (αF) terms need no fix. Primal path
  untouched ⇒ byte-identical. Reduces to Newmark DDM at αM=αF=1. Base-class quirk logged in
  [[LEDGER_quirks]]. Same test battery as PR1. **W3-I2 complete; ADR-52 remaining: W1-E2 only.**
- **2026-06-24 — W1-E2 shipped (#419): `ExplicitBathe*` 6→1 collapse.** The last wave.
  Folded the six explicit-Bathe class tags into ONE `ExplicitBathe` (33000) selected by
  orthogonal flags `-lnvd <alpha>` / `-sms <dtTarget>` / `-consistent`. **Architecture
  correction:** the six were a **2-base × 3-SMS-mode lattice** (`ExplicitBathe` /
  `ExplicitBatheLNVD` bases, each × {none, sms-lumped, sms-consistent}), NOT 6 flat siblings
  as the handoff said. The two bases' `newStep`/`update`/`commit` stepping was **byte-
  identical**; only LNVD's `formUnbalance()`+`addLocalDamping()`+`alpha_flac` and the SMS
  `domainChanged()` injection differed → both now flag-gated on the one class (LNVD via the
  unified `formUnbalance` override that reduces to the base when off; SMS via the
  `domainChanged` branch + the already-present `refineAccel` hook for `-consistent`).
  **Serialization (the crux):** the flag combo DERIVES the classTag (`tagForFlags`) so a
  serialized object reports its matching legacy tag; both brokers route all six tags →
  `ExplicitBathe::makeForBroker(tag)` (flags decoded from the tag) → `recvSelf` fills a
  fixed-size param superset. **No new tag; no retired tag freed** — 33002/09/10/11/12 stay as
  deprecated-but-recognized aliases. The five retired command names keep working (each OPS_
  parser preserves the exact historical positional grammar, forcing a fixed flag set).
  Footprint (all ledgered): delete 5 `.cpp` + 5 `.h`; rewrite `ExplicitBathe.{h,cpp}`; collapse
  the broker cases + drop the 5 retired includes (both brokers); drop the 5 build rows
  (CMakeLists + Makefile); annotate the 5 tags deprecated in `classTags.h`; consolidate the
  banner. Tests (`tests/test_adr52_w1e2_explicitbathe_collapse.py`, zone_a): each new flag form
  byte-identical to its legacy alias (5 combos); the **LNVD×consistent cross-layer** (the
  previously-untested combo) composes — stable, no nodal-mass mutation, relaxes a loaded bar to
  PL/EA, PCG active; `-lnvd 0.0` / `-consistent`-without-`-sms` reduce to base. The existing
  per-alias batteries still run under the deprecated names (free byte-identity regression).
  **6-lens adversarial Workflow review (9 agents): 2 confirmed MAJORs, 0 blockers — both fixed
  before merge, both defects ONLY on the NEW unified command surface (legacy aliases unaffected):**
  (1) the new `-lnvd <alpha>` optional-value peek used `strtod(OPS_GetString())`, but under
  openseespy `OPS_GetString` returns `"Invalid String Input!"` for a numeric PyFloat arg → the
  alpha was silently dropped to the 0.8 default (the test used 0.8 so it passed by coincidence;
  `-lnvd 0.0` reduce-to-base would have caught it). Fixed by adopting the proven contact `-soft`
  idiom (classify the peek by leading `-`, read the value with `OPS_GetDoubleInput`); the test now
  uses alpha=0.6 to expose a dropped alpha. (2) the unified parser did not downgrade
  `-cflAbort`/`-recompute` to report-only under `-sms` (the SMS aliases do) → `-sms … -cflAbort`
  would hard-abort at the pre-scaling Noh-Bathe limit, the MF-1 hazard `-sms` exists to avoid +
  a docstring contradiction. Fixed with a post-loop downgrade; new test asserts `-sms -cflAbort`
  == plain `-sms`. (Minor, documented intentional: the two LNVD-SMS aliases historically REFUSED
  `-cflAbort`/`-recompute`; routed through the shared impl they now DOWNGRADE like the other SMS
  forms — the consistent W1-E3 #394 behavior.) Serialization-collapse pattern logged in
  [[LEDGER_quirks]]. **ADR-52 COMPLETE** (deferred W2-E1 S3/S4/uri-ssp follow-ups optional,
  non-blocking).
- *Remaining waves:* **NONE — ADR-52 complete.** W1-E2 #419, W3-I2 #413+#415, W3-I3 NO-GO #410,
  W1-E3 #394, W1-I1a #396, W2-E1 #399+#403, W1-I1b #407.
