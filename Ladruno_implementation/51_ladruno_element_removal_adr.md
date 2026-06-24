---
title: Element Removal — runtime topology change for progressive collapse (pillar 1 of the AEM-like collapse program)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - subsystem
  - collapse
  - element-removal
  - contact
  - explicit
  - implicit
---

# Element Removal — runtime topology change for progressive collapse

> **ADR 51.** Pillar **1 of 3** of the AEM-like progressive-collapse program scoped
> in [[50_aem_opensees_scoping_report]] (§3): **runtime element removal** →
> **broad-phase recontact discovery** → **debris collision**. This ADR owns pillar
> 1 only; pillars 2–3 are the runtime-discovery generalization of
> [[39_ladruno_contact_domain_adr]] (see [[48_ladruno_contact_capstone_adr]]).
>
> **Revised 2026-06-23 after a 4-lens adversarial review** (dynamics/mechanics ·
> OpenSees architecture/code-reality · contact seam · scope/phasing/fidelity). The
> gate **confirmed the core mechanism** but **falsified the execution-model premise,
> two correctness claims, the debris-handoff seam, and the "generality/immediately
> useful" framing**. All corrections are folded in below and flagged `[GATE]`.
> Per-lens verdicts at the bottom under *Adversarial review log*.
>
> **Amended 2026-06-23 (project requirement):** **OpenSeesMP support is REQUIRED** —
> the fork must run partitioned collapse models. Q-PARALLEL is reframed from a v1
> *limitation* to a v1 *requirement*; the design is grounded in the fact that
> `PartitionedDomain` already routes removal to the owning `Subdomain` (see **H7**).
>
> Companion skills (synthesis + primary sources, installed for all sessions):
> `applied-element-method` (AEM/ELS, the target behaviour) and `opensees-contact`
> (`element_removal.md` = the full `RemoveRecorder` reading; `broad_phase_collision.md`;
> `separation_robust_integration.md`). PEER 2007/10 (Talaat & Mosalam) + the
> `RemoveRecorder` source are the primary references.

## What

A **runtime element-removal subsystem** that lets a converged nonlinear dynamic
analysis **delete a failed element from the `Domain` mid-run** — with the correct
transient physics and clean topology. This is the FEM-side of "element separation":
the structure is followed through the **continuum stage** (intact FEM) into the
**discrete stage** (members gone) without remeshing.

OpenSees already ships exactly this mechanism — **`RemoveRecorder`** (Talaat &
Mosalam, `SRC/recorder/RemoveRecorder.{cpp,h}`, `recorder Collapse`).

> `[GATE]` **"Modernized port, not green-field" — reframed to match where the risk
> actually is.** Decompose stock into (a) release + topology cleanup + mass
> transfer, (b) the criterion `switch`, (c) collision (which stock *lacks*). v1
> **ports (a)** (low-risk mechanics), **green-fields (b)** as a virtual-dispatch
> criterion interface, and the contact seam (c) **has no stock analog at all**. So
> "port" describes only the low-risk half; the two things v1 "adds over stock" are
> precisely the parts that cannot be ported because stock lacks them. The novelty —
> and the risk — lives in the criterion interface and (deferred, see H6) the contact
> seam, not in the ported mechanics.

**v1 ships:**

> `[GATE]` **Execution model — the original "post-convergence, no re-handle" framing
> was FACTUALLY WRONG and is corrected.** `Domain::removeElement`/`removeNode` set
> the `domainChange` flag (`Domain.cpp:1161/1181`). `record()` runs *inside*
> `Domain::commit()`, so removal lands at end-of-step *N*; the **next** `analyze()`
> sees `hasDomainChanged()` and runs a **full `domainChanged()` →
> `AnalysisModel::clearAll()` + `ConstraintHandler::handle()` + renumber + SOE
> rebuild** before evaluating step *N+1* (`DirectIntegrationAnalysis.cpp:206/423`).
> **Consequence (good):** there is **no stale graph**; the analysis layer re-forms
> itself — the "no analysis-loop surgery" claim survives, but *because of* the forced
> re-handle, not because nothing changes. **Consequence (cost):** every removal
> triggers a re-handle + renumber + symbolic re-factor on the next step — a real
> per-event cost, and the reason removal cadence matters. The entire Q-RECONTACT /
> debris framing is rebuilt on this corrected model (see H6).

- A **removal hook** on the `Recorder::record()` contract — called inside
  `Domain::commit()`, i.e. **only after a converged step** (a non-converged step
  returns before `commit()`, so a criterion cannot trip on an unconverged state —
  verified, `DirectIntegrationAnalysis.cpp:226-264`).
- A **removal-criterion interface** any element/material can expose (generalizing
  the `getRemCriteria` response that **only `ForceBeamColumn3d` (+ its thermal
  variant) implements** today — grep-verified). Ships with the **stock criterion
  set**: strain min/max, damage-index, Elwood shear/axial limit-state, plus a
  **`UserDefined`** callback.
  > `[GATE]` **The "generality" win is an *interface* claim, not a *capability*
  > claim, in v1.** The shipped criteria are exactly stock's, and the Elwood limit
  > states are intrinsically RC-column constructs on the limit-state material. So a
  > brick/shell becomes removable *in principle* but **v1 ships no new
  > physically-meaningful trip criterion for any non-beam element** beyond generic
  > strain min/max + `UserDefined`. Stated plainly so the generality sell isn't
  > vapor.
- **Per-instance / per-Domain state — NOT stock's process-global statics.**
  > `[GATE]` **Blocker fix.** Stock keeps `remEleList`/`remEles`/`numRemEles`/
  > `fileName`/`theFile` as **class statics** never reset between analyses (freed
  > only when the last recorder dies). Under OpenSeesPy multi-run / `wipe`+rebuild /
  > MP this leaks and *silently treats stale tags from a previous model as
  > already-removed*. The port **must** carry this state per-recorder (or on the
  > Domain), with a P1 gate: *two sequential analyses in one process give identical
  > results.*
- **Correct dynamic-equilibrium force release** — remove the element (true rank
  change via `Domain::removeElement`); the *missing resisting force* drives the
  transient in the next residual. Re-apply the removed element's **self-weight** and
  **decrement lumped nodal mass** — see H3 for the corrected (momentum-faithful)
  formulation; stock's `0.5·mₑ` rule is exact only for a 2-node prismatic member.
- **Full topology cleanup** — elemental loads; on node removal **nodal loads, SP-
  *and MP-* constraints** (`[GATE]` stock + the original ADR omitted MP_Constraints —
  diaphragm/`equalDOF` retained nodes would dangle and crash the handler; RC frames
  have diaphragms everywhere); region-membership purge of dead tags.
- **Integrator-agnostic + contact-safe**: validated under implicit (generalized-α /
  HHT) *and* explicit (central difference), and with a `LadrunoContactDomain`
  attached (removal must not silently corrupt committed contact **pair state** — see
  H6, the re-keying hazard).

**NOT in v1 (later phases / sibling ADRs):**
- **Automatic recontact / collision of removed bodies** — pillar 2 (the ADR-39
  generalization). `[GATE]` v1 **does not** ship a "debris record"; see H6.
- **Gradual / cohesive separation** (spring-by-spring, à la AEM). v1 removal is
  **whole-element / binary**.
- **The AEM discretization itself** ([[50_aem_opensees_scoping_report]] §3.1).
- General column **end-moment release** (`D_M=1` strut conversion).
- **Cross-partition floating-island AUTO-detection** under OpenSeesMP — a
  disconnected sub-assembly that *spans ranks* can't be seen by a per-partition
  incidence map; the analyst pre-list still works, auto-detection is post-v1. (NB:
  this is the *only* parallel piece deferred — partition-interior removal and
  boundary-node coordination are **in scope**; MP is a requirement, see **H7**.)

> `[GATE]` **NOT in v1 — *fidelity* (the honesty non-goal, imported from
> [[50_aem_opensees_scoping_report]] §2.4).** Binary FEM-element removal is a
> **lower-fidelity, mesh-path-dependent** approximation of AEM separation: the
> collapse path is **biased to element boundaries** (erosion-class artifact), and
> binary removal **over-stiffens until the trip, then over-releases** (vs AEM's
> gradual spring rupture). Per §2.4 we **do not claim accuracy parity with AEM** —
> the win is *workflow/automation*, a path through separation on the OpenSees base,
> **not** proven collapse-accuracy supremacy. Closing the fidelity gap (cohesive
> separation) is explicitly out of scope.

## Why

**[[50_aem_opensees_scoping_report]] §3.5 named element removal as pillar 1** — the
substrate every later pillar mutates (they all change `Domain` topology), so
de-risking the removal/cleanup/release path first is correct sequencing.

> `[GATE]` **"Immediately useful now" — scoped down to what v1 actually delivers.**
> The standalone near-term capability is **force-based RC-frame, column-removal
> alternate-load-path studies** — essentially the **stock `RemoveRecorder`
> envelope, modernized** (per-instance state, MP cleanup, integrator-correct). It is
> **not** "remove any element," and the named use case (the RC-frame study) lands at
> **P2**, not P1 (P1 only reproduces stock — see Phasing). Anything past
> column-removal-on-a-frame is gated behind **pillar 2** (recontact). The honest
> sell: a clean, modern, integrator-correct element-removal substrate now; general
> collapse later.

**Why not just use stock `RemoveRecorder`?** It works but is narrow and dated:
criteria only for `ForceBeamColumn3d`; **process-global static state** (broken under
re-run/MP); **partial cleanup** (no MP_Constraints, no region purge); parks removed
elements with no path to recontact; predates the fork's generalized-α / explicit /
contact substrate and was never checked against them.

**The dynamic-equilibrium release (the one correctness invariant).** Deleting a
load-bearing element drops its internal resisting force `f_k = aₖᵀrₖ(u⁻)` from the
assembled residual; displacement can't jump (inertia), so the structure is out of
equilibrium by `f_k`, applied as a step → the dynamic overshoot collapse codes
require. The right implementation: remove the element, let the next dynamic residual
carry `−f_k` (only the lost *self-weight* is re-applied, separately, as nodal loads;
the non-gravity part of `f_k` is correctly *not* re-applied — its removal *is* the
forcing).

> `[GATE]` **"Never quasi-static / never ε-soft-kill" — softened from absolute to
> bounded.** Both dismissals were over-stated:
> - **Static pushdown is the standard alternate-path baseline** (Marjanishvili &
>   Agnew 2006; the DAF≈2 is applied *because* you ran it statically) — indeed the
>   P2 gate compares against it. Quasi-static is wrong only when removal creates a
>   **kinematic mechanism** (then `K_T` is singular while the dynamic `K_eff = K_T +
>   (1+α)/(βΔt²)M` stays PD); in a redundant frame a static re-solve is well-posed.
> - **ε-soft-kill is wrong only in its *naïve* form** (no force bookkeeping → `κ ∼
>   k_max/(ε k_min) → ∞`, a floppy full-mass mode, slow wrong-frequency spring-back).
>   Soft-kill *with explicit stored-force release* is a legitimate erosion-avoidance
>   technique. We choose true removal for cleanliness, not because the alternatives
>   are universally invalid.

> `[GATE]` **Hook-timing lag (acknowledged).** Because `record()` is post-commit,
> removal at end-of-step *N* first perturbs the residual at *N+1* — fixed-step
> "through-the-event" integration with an O(Δt) impulse-smearing error. The
> event-aware Δt-shrink in H5 is what bounds it; they are now tied together rather
> than the release being sold as a clean Heaviside.

## Where

- **Port base:** `SRC/recorder/RemoveRecorder.{cpp,h}` (`RECORDER_TAGS_RemoveRecorder
  = 15`). Read for the `record()` hook, `elimElem`/`elimNode`/`elimSecondaries`/
  `updateNodalMasses`, the **static lists** (the B1 liability), and the
  `revertToStart`-not-`delete` park (the H6 ownership hazard).
- **Domain topology API:** `Domain.h:123-124` `removeElement`/`removeNode`. `[GATE]`
  their bodies do only the container removal + `domainChange()`
  (`Domain.cpp:1150/1170`); **all dependent-object cleanup is the caller's
  responsibility** — including `Domain::removeMP_Constraints(nodeTag)`
  (`Domain.cpp:1310`), which the cleanup path must call.
- **Re-form path:** `DirectIntegrationAnalysis.cpp:206` (`hasDomainChanged`) → `:423`
  (`handle`) — the mechanism that keeps the analysis layer consistent.
- **The criteria gap:** `getRemCriteria` only in `ForceBeamColumn3d.cpp` (+ thermal);
  Elwood limit states on `SRC/material/uniaxial/limitState/limitCurve/AxialCurve.cpp`.
- **Contact seam:** `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — `[GATE]` note
  its data model is surfaces + pair definitions + **tag/ordinal-keyed** path state
  (`FrictionState`, `MortarNormalState`); it has **no struct for a free rigid body**
  and its commit/revert hooks iterate only those state maps (the H6 reality check).
- **Integrator context:** generalized-α / HHT; `CentralDifference*`; the fork's
  `LadrunoMassScaling` / `CriticalTimeStep` / `LadrunoDynamicRelaxation`.
- **New code (proposed):** `SRC/domain/collapse/LadrunoElementRemoval.{h,cpp}` +
  `LadrunoRemovalCriterion`. **Class tag TBD — assign at P1 against
  `SRC/classTags.h`, record in [[LEDGER_implementations]]** (+ banner line). Touching
  `Domain` (MP cleanup, optional incidence map) → [[LEDGER_vanilla_files]] rows +
  `// Ladruno` markers.

## How

### H1 — the hook
Keep the post-commit `Recorder` pattern; removal triggers the forced re-handle next
step (corrected model above). Generalize the criterion check to virtual dispatch,
not a hard-coded `switch`.

### H2 — the removal-criterion interface
```
class LadrunoRemovalCriterion { virtual bool tripped(Element*, double t) = 0; };
```
Ship the stock set + `UserDefined` (Tcl/Python callback). `[GATE]` This decouples
removal from the element class at the *dispatch* level; it does **not** by itself
make non-beam elements physically removable (no new trip criterion for them in v1).

### H3 — release + cleanup (corrected for momentum + completeness)
> `[GATE]` **Mass/weight transfer — stock's formula is wrong off the 2-node-prismatic
> path; fix it.** Stock `updateNodalMasses` (a) subtracts the *translational*
> `0.5·mₑ` from **every** mass diagonal **including rotational DOFs**
> (`RemoveRecorder.cpp:1251-1254`) — dimensionally wrong, **does not conserve angular
> momentum**; (b) re-applies weight as a point **force** with **no moment couple**
> (`:1264`) — drops the fixed-end moment of distributed self-weight; (c) does the
> transfer **only if the user passes `-mass`**. v1 must instead **pull the
> element's consistent/lumped mass and body-force from the element** and re-apply the
> full nodal **force+moment** vector. **Until that lands, validated removability is
> restricted to 2-node prismatic members** (resolves the generality↔Q-MASS
> contradiction the gate flagged).
- `removeElement`; delete its `ElementalLoad`s from every `LoadPattern`.
- Node removal: delete `NodalLoad`s and `SP_Constraint`s, **and `MP_Constraint`s on
  BOTH roles** — `[GATE]` `Domain::removeMP_Constraints(int)` keys on the
  **constrained** node only (`Domain.h:130`), so it **misses** a removed node that is
  the **retained** (master) node of a rigid diaphragm / `equalDOF` / rigid link
  (`MP_Constraint` exposes `getNodeRetained()`/`getNodeConstrained()`,
  `MP_Constraint.h:75-76`) → a dangling constraint that crashes the handler on the
  next re-handle. The cleanup must **sweep `getMPs()` and delete any constraint where
  *either* role matches the removed node** — this is the realistic case for
  column-removal-with-a-diaphragm (the master is the retained node). Then purge the
  node/element tags from any `MeshRegion`.
- **Do not delete** the element/node object mid-analysis (`revertToStart` + park) —
  other recorders / contact adapters may hold the pointer. `[GATE]` Ownership of the
  parked object is **undefined** between the removal subsystem and any consumer;
  pin it (see H6) — do not hand a `revertToStart`'d, recorder-owned object to the
  contact engine as a "rigid body."
- **Rebuild the Rayleigh damping matrix** after removal. `[GATE]` `C = a₀M + a₁K`;
  removal changes M and K, so a stale C injects spurious damping at the event — the
  source skill flags this; the original ADR dropped it.

### H4 — topology cleanup automation `[GATE: demoted]`
A Domain-side node→element incidence map *could* auto-detect dangling/floating
sub-assemblies (dropping stock's analyst pre-list). But `[GATE]` an element graph
**already exists** and is **rebuilt from scratch on every `domainChanged`**
(`eleGraphBuiltFlag`); a hand-maintained incremental map is **new upstream `Domain`
state** with `currentGeoTag`/parallel-send implications — **not** a cheap free-rider.
This is a **convenience** phase (stock shipped useful studies for 15 years with
analyst pre-listing); it is the **first thing to cut under schedule pressure** and is
**re-ordered after** the strategically load-bearing work (see Phasing).

### H5 — integrator robustness
- Default the collapse lane to **generalized-α** with a **documented `ρ∞` knob**
  (`[GATE]` ρ∞≈0.7 is a *starting* value, not "the" default: it damps
  high-frequency content that in a fragmenting collapse **is** physical, and the
  dissipation-free explicit lane will **not** agree with it on peak demands — so it's
  a stated modeling decision).
- **Event-aware step control**: proactively shrink Δt for a few steps after a
  removal; `test EnergyIncr`; relax tol at the event; restore after. (Also bounds the
  H1 timing lag.)
- The dynamic `K_eff` mass term keeps a post-removal mechanism solvable — never fall
  back to a quasi-static step *in the mechanism case* (cf. the softened Why).
- **Explicit lane:** recompute `Δt_crit` after removal (the removed element or its
  neighbours may have set `ω_max`).
- **Energy accounting** `[GATE]`: binary removal `revertToStart`s the element, so its
  stored elastic strain energy is **deleted** from the system (not physically
  dissipated). State this and bound it; it matters for collapse energetics and is the
  honest cost of binary (vs cohesive) release.
- Reuse **viscous contact stabilization** (`-visc μc`) as the event regularizer when
  contact is engaged.

### H6 — the contact seam `[GATE: de-scoped to an inert stub or cut]`
> The original "register a debris record (nodes + velocity, rigid flag) with
> `LadrunoContactDomain`" seam **does not survive review** and is removed from v1's
> shipping scope:
> - **No home:** `LadrunoContactDomain` has no struct/API for a free rigid body; its
>   state is pair-definitions + tag-keyed path state.
> - **No owner:** once `removeNode` runs, nothing integrates the body's motion (the
>   contact engine is not in the analysis loop except as a stateless adapter source);
>   Q-RECONTACT's kinematics-ownership question is unresolved.
> - **No consumer:** pillar 2 will define what the record must carry (current-config
>   AABBs, per-body velocity for CCD, friction-history carry). Building the producer
>   against an undecided contract invites a retrofit *of the seam itself*.
>
> **v1 instead** either (a) ships **no debris record** (honest pillar boundary —
> removal stands alone; pillar 2 designs the record with its consumer), or (b) ships
> an **inert, default-empty `DebrisBody` stub asserted byte-identical/no-effect**,
> mirroring the contact lane's "store-the-definition, inert" P1b pattern
> ([[48_ladruno_contact_capstone_adr]] contract — the shared-hook guarantee is a
> *requirement on the PR, not inherited*).
>
> **The real contact interaction v1 MUST get right** `[GATE]` is **committed
> pair-state integrity under removal**. `FrictionState` is keyed by `(contactTag,
> slaveTag, segIndex)` and `MortarNormalState` by `(contactTag, slaveNodeTag)`, where
> `segIndex` is a **global master-segment ordinal**. If removal **renumbers a contact
> surface** (drops a master segment from the middle), surviving segments' ordinals
> shift and committed friction/λ_N **silently alias onto the wrong physical pair** —
> a wrong-answer, not a crash. v1 contract: *committed contact path-state must be
> stable across a removal-induced re-handle* — prove surfaces aren't re-enumerated,
> or re-key by surviving node-tag pairs, or document that removing a contact-surface
> member invalidates that contact's history.

### H7 — parallelism (OpenSeesMP) is a REQUIREMENT, not a deferral `[GATE: amended]`
Element removal **must run under OpenSeesMP** (partitioned collapse models). Source
reality is favourable: the **container plumbing already exists** —
`PartitionedDomain::removeElement` searches the subdomains and routes removal to the
owning `Subdomain` (`PartitionedDomain.cpp:654/678`); `removeNode` sweeps the local
domain **and every subdomain** (`:692/706`); `removeMP_Constraints` accumulates
across them (`:814/828`). So the work is **coordination, not capability**:
- **Per-instance state (B-static) is the hard precondition.** Process-global static
  removed-lists aren't merely leaky — they are **incoherent across MPI ranks**. Each
  rank/`Subdomain` must own its removal state. This is the single gating fix; nothing
  parallel works until it lands.
- **Boundary/shared-node removal coordination — the crux.** A partition-boundary node
  lives in several `Subdomain`s (+ as an external node on the master). A criterion
  that trips on a *remote actor's local element* makes a **local decision** whose
  **action may be global** (the shared node + its mass/load share must be removed
  consistently on every owner). `PartitionedDomain::removeNode` already does the
  multi-subdomain sweep *when driven from the master*; the missing piece is a
  **rank-consistent removal protocol** so a removal decided on one rank propagates.
  This is its own phase.
- **The forced re-handle is a parallel re-partition / re-number** under
  `PartitionedDomain` (the ADR-30 parallel numberer) — expensive, so removal cadence
  matters even more than in serial.
- **Floating-island AUTO-detection across ranks is deferred** (per-partition
  incidence map can't see a span-the-ranks island; analyst pre-list still works).
- **Contact-under-MP is NOT pillar 1's problem.** Pillar 1 (removal) has **no
  contact**, so it reaches MP well before recontact does — ADR-39 explicitly defers
  *MPI contact decomposition* to v2/v3. Clean story: **removal can be MP-capable now;
  recontact-under-MP (pillars 2–3) is a separate, later parallel-contact problem.**

**Honest gradient:** partition-*interior* removal under OpenSeesMP is tractable once
B-static lands (route + per-instance state); **boundary-node coordination** is the
genuine new work (its own phase); cross-rank floating-island auto-detection is
deferred.

### Testing
- `[GATE]` **Oracle, corrected.** The ≤2× DAF "vs analytic" check is well-posed
  **only for a linear, undamped SDOF** step release (the `1−cos` bound). Use it as a
  **linear verification** case only. For the **nonlinear** column-removal lane, the
  gate is an **energy + momentum balance** check, *not* a scalar DAF (yielding +
  damping + MDOF modal superposition break the 2× bound — it can be exceeded or
  undercut).
- Port the **trussed-canopy** benchmark (PEER §6.5) — cascade + secondary/floating.
- **Column-removal** RC frame **with a rigid diaphragm** under generalized-α —
  exercises MP-constraint cleanup (B2) and alternate-load-path redistribution.
- **Static-state gate** `[GATE]`: two sequential `analyze()` runs in one process give
  identical results (catches the static-list liability).
- **Contact pair-state gate** `[GATE]`: remove a non-contacting element adjacent to a
  contacting pair → the surviving pair's committed `FrictionState`/`λ_N` is
  **byte-identical** before/after (catches the re-keying aliasing — the old
  "active set uncorrupted" clause was a tautology the forced re-handle passes
  trivially).
- **No-regression**: with no removal recorder, byte-identical to stock (Zone-A).
- **Integrator matrix**: implicit + explicit; with/without an attached
  `LadrunoContactDomain`.

## Phasing

> `[GATE]` **Re-labeled and re-ordered after review.** P1 is a **regression floor**,
> not a value rung (it reproduces a 2007 capability — necessary discipline, zero new
> user capability). The first user-facing collapse study lands at **P2**. P3
> (incidence-map automation) is **demoted below** P4 and flagged cut-first. The
> debris seam is an **inert stub** (or deferred), not a shipped record.

| Phase | Scope | Gate |
|---|---|---|
| **P1** (regression floor) | Generalized hook + `LadrunoRemovalCriterion` (stock criteria) + release/cleanup **with MP-constraint + region purge** + **per-instance state**; byte-identical no-op when unused | Trussed-canopy reproduces stock; **two sequential analyses identical**; Zone-A green |
| **P2** (value) | Dynamic-release correctness (momentum-faithful mass/weight + Rayleigh rebuild) + event-aware stepping under generalized-α/HHT **and** explicit; **the column-removal RC-frame study** | Linear-SDOF overshoot oracle ≤2%; nonlinear lane passes energy+momentum balance; diaphragm column-removal converges |
| **P3** (convenience, cut-first) | Auto dangling/floating detection via the Domain incidence map | Dangling/floating auto-detected on a meshed-column removal |
| **P4** (seam, inert) | Contact-safe removal: **committed pair-state integrity** under removal; optional inert `DebrisBody` stub | Remove-near-contact → surviving pair friction/λ_N byte-identical; stub (if present) asserted no-effect |
| **P-MP** (required) | OpenSeesMP: per-instance state coherent across ranks (rides on P1); **rank-consistent boundary-node removal protocol**; parallel re-handle tolerates per-step removal | Partition-interior column-removal gives bit-identical results to serial; a boundary-node removal is consistent across ranks; no static-state leak across ranks |
| **→ pillar 2** | Recontact: runtime body-discovery broad phase **(new code, NOT a generalization of the bucket sort — see Risks)** designs + consumes the debris record; *contact-under-MP is its own v2/v3 problem* | (sibling ADR — ADR-39 generalization) |

## Risks / open questions

> [!question] Q-RECONTACT — design the debris record WITH pillar 2, not before
> Who owns a removed body's kinematics, and what does the record carry? `[GATE]`
> Resolved direction: **do not build the producer before pillar 2's broad-phase
> consumer pins the contract** (current-config AABBs, per-body velocity, friction
> carry). v1 ships at most an inert stub.

> [!question] Q-BUCKETSORT-REWRITE — pillar 2 is a rewrite, not a generalization
> `[GATE]` `LadrunoContactBucketSort` is a **uniform grid sized to the median master-
> segment diagonal**, built from **reference coords** over **predefined master/slave
> surfaces**. Every pillar-2 need breaks a structural assumption: current-config
> re-emit (rebuild every step), open-universe body discovery (no predefined
> master/slave), size-dispersed debris (a single median cell size is **fatal** —
> `broad_phase_collision.md §3`), edge-edge + self-contact (a different narrow phase).
> The skill prescribes a **BVH/loose-octree/hierarchical-hash replacement**. Call
> pillar 2 a **head start on the static case + new code**, not a refactor.

> [!question] Q-COHESIVE — binary vs gradual separation
> Binary removal injects the full `f_k` transient that a cohesive `G_f` law would
> *bound*; generalized-α then **dissipates** it (so the "physics" becomes
> ρ∞/Δt-tuning-dependent — it *hides* an artifact). Is a cohesive traction–separation
> law worth it, or is binary + dissipation enough for the target cases? See
> `separation_robust_integration.md §7`.

> [!question] Q-MASS — consistent mass/body-force for non-beam elements
> `[GATE]` Promoted from "deferred" to a **P2 in-scope requirement** for any element
> v1 claims to support beyond 2-node prismatic (the generality↔mass-model
> contradiction). Pull consistent mass + full force+moment from the element.

> [!question] Q-CRITERIA — only *externally-cached* criterion state is at risk
> `[GATE]` Reworded: the post-commit hook placement already guarantees a criterion
> can't trip on an unconverged state, and element-committed state reverts with the
> Domain. The risk is **only** a criterion that caches path-state *outside* the
> element's own commit cycle — that must commit/revert with the Domain.

> [!question] Q-PARALLEL — OpenSeesMP is a REQUIREMENT (see H7)
> `[GATE]` Reframed from limitation to requirement. The container plumbing exists
> (`PartitionedDomain` routes removal to subdomains); the open design is the
> **rank-consistent boundary-node removal protocol** and confirming the parallel
> re-handle (ADR-30 numberer) tolerates per-step topology change at scale. Hard
> precondition: per-instance state (B-static) — static globals are incoherent across
> ranks. Cross-rank floating-island auto-detection stays deferred. Contact-under-MP
> is pillars 2–3, not this ADR.

## Ledger / banner

On implementation: [[LEDGER_implementations]] row (feature, kind, class tag, files,
status, PR) + `Ladruno_scripts/banner_features.txt` line (`patch_banner.py`);
[[LEDGER_vanilla_files]] rows + `// Ladruno` markers for the `Domain` MP-cleanup /
incidence-map touches and any `ForceBeamColumn3d` change; [[LEDGER_quirks]] for the
static-state and pair-state-re-keying gotchas.

## Adversarial review log

4 independent lenses, 2026-06-23, each tasked to **falsify** against the source.
Net: **core mechanism CONFIRMED; execution-model premise + 2 correctness claims +
the debris seam + the generality framing FALSIFIED**; all folded in as `[GATE]` above.

**Lens A — dynamics/mechanics.** REFUTED: the ≤2× DAF as a *nonlinear* acceptance
oracle (linear-SDOF bound only) → B-oracle; "conserves linear *and* angular
momentum" (code subtracts translational mass from rotational diagonals + no weight
moment) → B-mass. WEAKENED: hook-timing exactness (O(Δt) lag); "never quasi-static /
soft-kill" (over-stated); binary "is enough" (dissipation *hides* the artifact).
MISSING: Rayleigh rebuild, energy accounting, explicit Δt_crit recompute. CONFIRMED:
the dynamic-release invariant is the strongest part.

**Lens B — OpenSees architecture.** CONFIRMED: the analysis layer re-forms via
`domainChanged` → no stale graph (the load-bearing correctness claim holds);
criterion runs post-convergence only; `getRemCriteria` is `ForceBeamColumn3d`-only.
REFUTED: "safe to opt in" — **process-global static state** breaks re-run/MP →
B-static. LANDMINES: **MP_Constraint cleanup omitted** → B-MP; region staleness;
"topology removal supported" is true only for the container (cleanup is the caller's
job); the incidence map is new upstream state, not free reuse.

**Lens C — contact seam.** The decisive lens. REFUTED the original "post-commit, no
re-handle" premise (removal *forces* a re-handle that rebuilds the stateless adapters
— so the active set is safe, but for the opposite reason stated) → B-execmodel.
REFUTED the **debris-handoff seam** (no home/owner/consumer in
`LadrunoContactDomain`) → de-scoped to a stub. WEAKENED "generalizes the bucket sort"
→ it's a rewrite. Surfaced the **deepest unexamined bug**: committed
`FrictionState`/`MortarNormalState` **re-keying / aliasing** when removal renumbers a
contact surface (wrong-answer) → B-pairstate.

**Lens D — scope/phasing/fidelity.** "Right bet, right build order, wrong sales
pitch." WEAKENED: "immediately useful now" (the study lands at P2; ≈ stock envelope);
the "generality" win (interface, not capability — no non-beam criterion); "port, not
green-field" (the *new* parts can't be ported). REFUTED-by-omission: the parent §2.4
**fidelity caveat** (no AEM accuracy parity; mesh-path-dependence; erosion artifacts)
→ B-fidelity. WEAKENED phasing: P1 is a regression floor; P3 is gold-plating ahead of
P4; the consumer-less debris record is YAGNI. CONCEDED: pillar order, positioning vs
ADR 39/48/50, and the dynamic-release invariant all survive.

**Blocker ledger (all folded in):** B-execmodel (re-handle premise) · B-static
(per-instance state) · B-MP (constraint cleanup) · B-mass (momentum-faithful
transfer) · B-oracle (linear-only DAF) · B-pairstate (contact re-keying) ·
B-debris (seam → stub) · B-fidelity (import §2.4). Should-fixes: Rayleigh rebuild,
energy accounting, soften the dismissals, ρ∞ knob, region purge, P3↔P4 re-order, P1
relabel, "bucket-sort = rewrite." Survived unchallenged: the rank-change +
dynamic-residual mechanism, the pillar order, the no-duplication positioning.

## See also
- [[50_aem_opensees_scoping_report]] — the parent program (3 pillars; this is #1).
- [[39_ladruno_contact_domain_adr]] · [[48_ladruno_contact_capstone_adr]] — the
  contact engine pillar 2 rewrites for recontact.
- [[41_ladruno_mortar_alm_contact_adr]] — mortar/ALM lane (narrow phase).
- Skills: `applied-element-method` · `opensees-contact`
  (`element_removal.md`, `broad_phase_collision.md`,
  `separation_robust_integration.md`).
- Primary: Talaat & Mosalam, PEER 2007/10 + EQE 38(5):609–634; Elwood, CJCE
  31(5):846–859; the `RemoveRecorder` source.
