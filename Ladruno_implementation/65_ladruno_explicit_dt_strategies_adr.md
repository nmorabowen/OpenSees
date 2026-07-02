---
title: Explicit Δt strategies beyond SMS — adaptive-Δt driver, subcycling, IMEX partitioning
project: Ladruno
status: scoping (decision-capture; Route B accept-on-trigger, Routes A/C deferred with triggers)
priority: medium
owner: nmora
tags:
  - adr
  - integrator
  - explicit
  - mass-scaling
  - subcycling
---

# ADR-65 — Explicit Δt strategies beyond SMS

> ADR-65. Status: **SCOPING — decision-capture, no code.** Next free ADR slot after
> ADR-64 (shell-to-solid tie) → this is **65**.
> Family: ADR-35/36/38 (HRZ lumping + selective/consistent mass scaling, the shipped
> answer this ADR looks past) · ADR-52 (integrator strengthening / ExplicitBathe flag
> class) · ADR-49/49a (integrator study + scorecard) · ADR-30 (explicit constraint
> projection — interface machinery any partition scheme must respect) · ADR-39/41/60
> (contact engine — re-emission cadence interacts with sub-stepping) · ADR-29/58/62
> (rigidification substitutes: RBE2, RigidBody scoping, LadrunoTie).
> **Revision v2** (2026-07-02) — post 3-lens adversarial gate. v1 findings fixed:
> `criticalTimeStep()` returns the PRE-scaling pencil under SMS (v1 claimed post-scaling
> — verified false, twice independently); §Why argued against *lumped* SMS while the
> fork ships *consistent* SMS (ADR-38) which already removes the period-shift failure
> mode; leap-frog stagger must be re-seeded on every Δt change; the P0 oracle is an
> optimistic upper bound and its histogram is not yet dumpable; GC stability and the
> energy-audit epistemics were overstated; triggers made operational; "Route 0 —
> shipped substitutes" row added; Route B downgraded ACCEPT → **accept-on-trigger**;
> Olovsson citation corrected (IJNME 63, +Unosson).
>
> **Decisions (v2, 2026-07-02):**
> - **Route 0 (shipped substitutes): the default.** Consistent SMS (ADR-38),
>   rigidification (RBE2 33012 / LadrunoTie / ADR-58), mesh repair, whole-model
>   implicit. Every other route must first defeat this row.
> - **Route B (adaptive-Δt driver): ACCEPT-ON-TRIGGER.** Small, analysis-layer, the
>   shrink-on-failure hooks exist. Build after a build-free headroom oracle (§Route B
>   P0) shows ≥2× wasted steps on a named campaign/battery case — not opportunistically.
> - **Route A (subcycling / multi-time-step): DEFER with trigger.** The real
>   alternative to SMS in its residual niche, major undertaking. Activate per
>   §Triggers, and run the (corrected) P0 payoff oracle before any C++.
> - **Route C (implicit–explicit domain partitioning): DEFER, long-horizon.**
>   Analysis-architecture project, not an integrator.

---

## What

Strategies for the regime where the explicit critical time step Δt_cr is throttled by
a small subset of the mesh, evaluated against what the fork already ships:

- **Route 0 — shipped substitutes (baseline, not a build).** Lumped SMS (ADR-36),
  consistent SMS (ADR-38), rigidify the tiny stiff cluster (RBE2 / LadrunoTie /
  future ADR-58 RigidBody), mesh repair, or run the whole model implicit.
- **Route A — Subcycling / multi-time-step (MTS), Belytschko–Yen–Mullen nodal
  partitions.** Small elements step at Δt/k instead of getting heavier. No added
  inertia; trades it for interface-transmission artifacts (see §Route A).
- **Route B — Adaptive-Δt driver.** An analysis-layer runner that re-queries
  `criticalTimeStep()` as stiffness evolves and grows/shrinks the analysis Δt with
  hysteresis, sub-dividing on step failure via `revertToLastStep()`.
- **Route C — Implicit–explicit (IMEX) domain partitioning.** Stiff region implicit,
  rest explicit, coupled at the interface (GC-style dual Schur).

This ADR records the decision on each and concrete activation triggers. No code, no
classTag reservation, no ledger rows (doc-only).

## Why — the escalation ladder

The failure mode this family addresses: a few small/stiff elements set Δt_cr for the
whole mesh. The fork's shipped answer is mass scaling, and it is a *ladder*, not one
tool:

1. **Lumped SMS (ADR-36)** adds real translational inertia at the offending nodes.
   Cheap, robust; its artifact is global — ADR-38's gate measured **f₁ shift −53.4%**
   at an aggressive dtTarget. The `-maxAddedMass` diagnostic (default 5% of total
   translational element mass, warn-and-proceed at `CentralDifferenceSMS.cpp:261–265`
   — nothing is actually capped) is the lumped-mode tripwire.
2. **Consistent SMS (ADR-38, shipped, MPI-validated)** builds the Olovsson-type
   artificial mass to annihilate rigid-body-mode inertia: same gate measured
   **f₁ shift −0.17%** at the same dtTarget, at roughly half the added mass. The
   *global period-shift* objection to mass scaling is therefore already solved in
   this fork. What consistent SMS does **not** fix: added inertia still loads the
   local/deformation modes — wave transmission *through* the scaled bin, impact
   response *of* the scaled interface hardware.
3. **Rigidification** — if the small stiff cluster is effectively rigid hardware
   (couplers, anchors, platens), RBE2 / LadrunoTie remove it from the pencil entirely,
   zero inertia artifact, zero new machinery.

**Route A exists for what survives the ladder**: the fast bin is load-bearing,
deformable, and its *local dynamics matter* (wave transmission / impact through the
interface), so any added inertia — lumped or consistent — pollutes exactly the
response being studied, and rigidification is wrong by definition. In that residual
niche, subcycling's work scales with Σ_p N_p/Δt_p instead of N/Δt_min with **no added
inertia** (its own artifact is interface wave reflection from the Δt/dispersion
mismatch between partitions — same mechanism as an abrupt mesh-size change — plus GC
interface dissipation if GC is used; "no inertia pollution" ≠ "exact").

Route B is orthogonal: on softening/contact problems the safe Δt is time-varying and
today the user hand-picks a globally safe value or babysits restarts. The comment at
`ExplicitBathe.h:111–113` already prescribes the driver pattern.

Route C is the "why not both" endpoint, but it is an analysis-architecture project and
its headline use case has Route 0 substitutes today.

## Route A — Subcycling / multi-time-step (DEFERRED)

### Sketch

Nodal-partition explicit–explicit subcycling (Belytschko–Yen–Mullen 1979 — the paper
that introduced explicit–explicit partitions with distinct time steps; Neal &
Belytschko 1989 for non-integer ratios):

1. From the per-element pencil Δt_e, bin nodes into partitions with steps Δt, Δt/k₁ …
   (integer ratios, few bins).
2. Advance the small-Δt partition k times per major step; interface nodes see
   interpolated (constant-velocity) kinematics from the large-Δt side.
3. Commit/collision/recorder semantics stay on the **major** step.

Modern variants replace the interpolated interface with a Lagrange-multiplier /
dual-Schur interface: **GC method** (Gravouil & Combescure 2001) and **PH method**
(Prakash & Hjelmstad 2004 — energy-preserving, removes GC's interface dissipation).
Characterize GC honestly: the coupling is *provably non-destabilizing* — the coupled
scheme is stable provided **each subdomain satisfies its own scheme's stability
limit** (the fast bin keeps its own CFL; no Δt limit disappears) — and it is
dissipative at the interface for step ratio > 1 (conservative at ratio 1).

If we build, the interface treatment is the load-bearing choice. Classic
constant-velocity (Belytschko-type) subcycling carries Daniel's (1998) result:
**narrow bands of unstable Δt-ratio/frequency combinations** ("statistical
instability"). Two consequences for the validation plan:

- The instability is caught **a priori** by a linear two-partition
  amplification-matrix / spectral-radius sweep over step ratio and frequency — that
  sweep is a *mandatory* gate if the interpolated interface is chosen.
- An energy audit alone is **not** a reliable catcher: a battery case at a benign
  ratio passes while production sits in a bad band. `EnergyBalanceRecorder` stays in
  the gate as a runtime tripwire, not as the proof.

### Why it's a major undertaking (the honest cost list)

- **Interface stability** — the core math risk (Daniel bands for interpolation; GC's
  interface dissipation must be measured, not assumed).
- **Framework cadence** — OpenSees has one global time loop; `Domain::commit()` fires
  per analysis step and downstream machinery keys off it (ADR-60's contact re-sort
  trigger is commit-anchored, `Domain.cpp:2201–2208`; recorders record per commit).
  Sub-steps must be invisible to all of it: materials in the small partition commit k×
  per major step, contact re-emission and constraint projection (ADR-30) need a
  defined cadence for pairs/constraints that *span* partitions.
- **MPI lockstep** — OpenSeesMP exchanges interface state per step; the small
  partition multiplies communication by k unless processor boundaries are aligned with
  Δt partitions (a partitioner-aware change).
- **Recorder/output semantics** — sub-step states are internal; energy bookkeeping
  must integrate across sub-steps.
- **SMS composition** — the natural production shape is hybrid (subcycle the worst
  bin, SMS the mild tail), doubling the interaction surface.

Effort: **L** (multi-month; thesis-chapter scale). New integrator class (classTag from
the ≥33000 pool, `LEDGER_implementations.md` row) + partitioner + interface kernel +
MPI work + validation battery incl. the amplification sweep.

### P0 payoff oracle (run before any code)

Build-free in spirit, with one plumbing caveat. Evaluate:

```
speedup_raw ≈ (N_total / Δt_min) / Σ_bins (N_bin / Δt_bin)
```

**This is an optimistic upper bound**, and the corrections bite hardest exactly in the
motivating case (scattered small stiff elements):

- **Interface halo**: every element touching a fast node evaluates at the minor
  cadence — each fast region drags a one-element ring into the fast bin. Bin on
  halo-expanded element sets, not the raw histogram.
- **Integer-ratio rounding**: achievable Δt_bin = Δt_major / integer — up to ~2×
  below the histogram value per bin.
- **Minor-cadence fixed costs**: contact broad-phase, ADR-30 projection, and MPI
  interface exchange scale with sub-step *count*, not fast-bin size.

Decision rule: as a **veto** the raw formula is sound (< 2× raw ⇒ don't build,
guaranteed). As a **go** signal, require ≥2× on the halo-expanded, integer-rounded
estimate (or ≥4× raw as the lazy haircut).

**Plumbing caveat**: the pencil is computed per element inside
`buildMassScaling`/`computeCriticalTimeStep` (`LadrunoMassScaling.h`) but only
aggregates are exposed (`MassScalingReport`, min + governing tag). The oracle needs a
small dump hook (e.g. a `-dumpDt <file>` option) or an offline script re-deriving
Δt_e from element K/M — budget that before calling the oracle "build-free".

## Route B — Adaptive-Δt driver (ACCEPT-ON-TRIGGER)

### P0 headroom oracle + trigger (before writing the driver)

Consistent with the fork's decision culture (ADR-61 shelved a sound design on a P0
oracle; ADR-63 Q-MULTISHELL not built for zero consumers), Route B does not ship on
"opportunistically". Build-free oracle: on a named Zone-B softening/contact battery
case (or the next explicit campaign's model), log `criticalTimeStep()` per step at a
coarse cadence and integrate the headroom ratio `dt_cr(t) / dt_fixed`. **Trigger:
oracle shows a variable-Δt run would save ≥2× the steps (or the fixed-Δt run dies and
restart-babysitting is the current workaround) on at least one named consumer.**

### Sketch

An analysis-layer runner, not a new scheme. Hooks that exist and are verified:

- `criticalTimeStep()` is valid before the first step (`CentralDifferenceLadruno.cpp:335`,
  computed unconditionally in `domainChanged`) and refreshable via `-recompute N`
  (`:111`, `:466`);
- `revertToLastStep()` reseeds the leap-frog starter and preserves the recompute
  cadence (`:420–436`) — and the reseed is valid standalone, not just after failure.

Driver loop (v1, pseudo):

```
dt = safety * criticalTimeStep()
while t < t_end:
    ok = analyze(1, dt)   # KE-breaker trip surfaces as a failed analyze — one path
    if not ok:  revertToLastStep(); dt *= 0.5; retry (bounded); flag the run
    every N steps (or after domainChange): dt_cr = criticalTimeStep()
        shrink immediately if dt > safety*dt_cr
        grow with hysteresis (≤ ~5%/step, only after M clean steps)
    on ANY dt change (grow or shrink): revertToLastStep() after the last
        commit, then newStep(dt_new)      # stagger reseed — see below
```

**Stagger reseed is mandatory, not hygiene.** `newStep()` advances
`Vhalf += Δt_new·a_n` with no previous-Δt memory (`CentralDifferenceLadruno.cpp:574`);
after a committed step at Δt_old the half-step velocity is centered for Δt_old, and
the correct advance is `((Δt_old+Δt_new)/2)·a_n`. A bare `newStep(dt_new)` injects a
`(Δt_new−Δt_old)/2·a` velocity error at **every** ratchet — small per change under 5%
hysteresis but systematic during ramps. `revertToLastStep()` re-arms the starter and
rebuilds `v_{−1/2}` at the new Δt from committed state, which makes the loop exact.
(The shrink-on-failure path was already correct for this reason.)

Deliverable: `Ladruno_scripts/adaptive_explicit_driver.py` (+ thin Tcl proc), tests in
the Zone-B battery on the oracle's named case. While in there, fix the stale pointer
at `ExplicitBathe.h:113` — it cites an "adaptive_substep example" that doesn't exist;
name the script to match or update the comment. Effort: **S** (days).

### Constraints and open edges

- **SMS interop — the hard truth (v1 verified):** under SMS, `criticalTimeStep()`
  returns the **PRE-scaling** element pencil (`CentralDifferenceLadruno.cpp:742–746`;
  same in `ExplicitBathe.cpp:1435–1441`). The post-scaling effective limit
  (`setSMSEffectiveLimit` ← `minDtSelfReport`, P3 #475) is consumed *only* by the
  `newStep` report and **has no getter**. A naive driver under SMS would collapse Δt
  to the un-scaled sliver limit and defeat SMS entirely. Therefore **v1 refuses to
  run under SMS** (assert no SMS integrator). A v2 SMS-aware driver has two explicit
  prerequisites: (a) a `getSMSEffectiveLimit()` accessor + command plumbing; (b) an
  answer to re-scaling: raising Δt above the effective limit would need a fresh
  scaling pass at a new dtTarget, and **no such path exists today** — `dtTarget` is
  fixed at construction; `-recompute` under SMS is parse-and-discard report-only
  (MF-1; the P3 "consume" is argument hygiene, `CentralDifferenceSMS.cpp:104–111`,
  `:141–148`). Suspect the answer is "SMS runs stay fixed-Δt" — an added-mass jump
  mid-run is a momentum discontinuity.
- **Variable-Δt output**: recorders emit non-uniform time series; post-processing
  keys on the stored time (LadrunoRecorder writes a per-step `TIME` attribute) —
  document, don't "fix".
- **Re-query cost**: the tangent-based refresh is O(N_elements) eigen work per call
  (the `ExplicitBathe.h` caveat block says "use a large N" for exactly this reason);
  the driver's N needs cost guidance in the script's defaults, not a magic number.
- **MPI**: Δt decisions and revert/retry must be rank-collective — a rank-local
  failure with local halving forks Δt across ranks or deadlocks. v1: allreduce the
  (ok, dt_cr) pair before deciding.
- **KE-breaker semantics**: the breaker trips pre-commit, so single-step revert is
  state-valid — but divergence detected at step n was typically seeded over several
  already-committed steps at too-large Δt. Retry-halving stops the compounding; it
  does not repair committed pollution. The driver flags the run (results before the
  trip are suspect), never just silently retries.
- **Restart/checkpoint** under a non-uniform Δt history: out of v1 scope; document
  that a restarted run re-derives Δt from `criticalTimeStep()`, not from the history.
- **apeGmsh exposure**: none in v1 (script-level tool); revisit only if the driver
  becomes the default explicit runner.

## Route C — IMEX domain partitioning (DEFERRED, long-horizon)

Element/nodal implicit–explicit partitions (Hughes & Liu 1978, *element* partitions;
Belytschko & Mullen 1978, *nodal* partitions), or GC-style multi-time-step with the
stiff subdomain implicit. In OpenSees terms this is **not an integrator**: one
`Analysis` owns one `Integrator`, one solve path. IMEX needs either a compound
integrator (two child integrators + interface Lagrange multipliers) or a
partitioned-analysis layer — both analysis-architecture projects on the scale of
ADR-30's handler work or larger.

Route 0 substitutes for the headline use case (stiff structure on soft soil, DRM
input): whole-model implicit HHT (fine at seismic frequencies), or whole-model
explicit with consistent SMS confined to the structure.

Build-order note, stated conditionally: **if** Route A is built with a GC dual-Schur
interface, Route C reuses most of that interface machinery and A-before-C is the
cheap order. If A ships the interpolated interface (or never ships), C inherits
nothing and is evaluated standalone on its own trigger.

Effort: **L–XL**. No trigger currently in sight; revisit on a concrete project where
both Route 0 substitutes measurably fail.

## Decision matrix

| Route | Verdict | Effort | Dynamic artifact | Framework risk | Blocked on |
|---|---|:---:|---|:---:|---|
| 0 — shipped substitutes (consistent SMS / rigidify / mesh repair / implicit) | **DEFAULT** | — | consistent SMS: deformation-mode inertia only (f₁ −0.17% measured) | none | nothing — shipped |
| B — adaptive-Δt driver | **ACCEPT-ON-TRIGGER** | S | none (with mandatory stagger reseed) | low (reseed discipline; MPI-collective decisions) | headroom oracle ≥2× on a named consumer |
| A — subcycling / MTS | **DEFER + trigger** | L | interface wave reflection (Δt/dispersion mismatch) + GC dissipation; no added inertia | high (commit cadence, MPI, recorders, contact) | §Triggers + corrected P0 oracle |
| C — IMEX partitioning | **DEFER** | L–XL | interface coupling artifacts | very high (analysis architecture) | a real project; GC machinery from A *if* A goes GC |

## Triggers (operational — what re-opens the deferred routes)

Route A activates when **all** of, evidence attached to the reopened ADR:
1. **≥2 distinct production models within a year** where, *running consistent SMS*
   (ADR-38 — not lumped), the added-mass diagnostic still exceeds the admissible
   fraction for the study, **and** a documented before/after modal or
   wave-transmission diff shows the residual deformation-mode inertia distorts the
   response being studied (>1% on the metric of record). Note the tooling gap: ADR-36
   left the modal check "optional"; a before/after-scaling modal-diff utility is part
   of the evidence cost, and "trips the `-maxAddedMass` warning" alone is **not** an
   event — it's warn-and-proceed.
2. Route 0 exhausted: rigidification inapplicable (the fast bin is load-bearing and
   deformable) and mesh repair unavailable (geometry-driven refinement).
3. The **corrected** P0 payoff oracle (halo-expanded bins, integer-rounded ratios)
   predicts ≥2× over the best admissible consistent-SMS configuration — after the
   dump-hook plumbing is budgeted.

Route B activates on its own §P0 headroom oracle (≥2× wasted steps or
death-and-babysit on a named consumer).

Route C activates only on a named project where both Route 0 substitutes measurably
fail — and is evaluated standalone unless Route A already shipped a GC interface.

## Risks / open questions

> [!question]
> Route A interface: constant-velocity interpolation (cheap; Daniel band-instability
> ⇒ mandatory amplification-matrix sweep) vs GC dual Schur (non-destabilizing,
> dissipative for ratio > 1, heavier)? Leaning GC if built — but the choice gates the
> Route C reuse claim, so decide it *first* if C's trigger looks live.

> [!question]
> Route B v2 under SMS: is "SMS runs stay fixed-Δt" the final answer, or is a
> re-scaling pass at a new dtTarget acceptable despite the added-mass momentum jump?
> Blocked on the `getSMSEffectiveLimit()` accessor either way.

> [!question]
> Route A × contact: pairs spanning partitions — minor or major cadence? Minor is
> correct physics (impact happens in the fast bin), but ADR-60's re-emission trigger
> is commit-anchored. Needs a design note before any build.

## References

- Belytschko, Yen & Mullen (1979), *Mixed methods for time integration*, CMAME 17/18:259–275.
- Neal & Belytschko (1989), *Explicit–explicit subcycling with non-integer time step
  ratios*, Computers & Structures 31.
- Daniel (1998), *A study of the stability of subcycling algorithms in structural
  dynamics*, CMAME 156:1–13 — the band ("statistical") instability result, found via
  amplification-matrix analysis.
- Gravouil & Combescure (2001), *Multi-time-step explicit–implicit method for
  non-linear structural dynamics*, IJNME 50:199–225 — the GC method.
- Prakash & Hjelmstad (2004), *A FETI-based multi-time-step coupling method for
  Newmark schemes in structural dynamics*, IJNME 61(13):2183–2204 — the PH method.
- Hughes & Liu (1978), *Implicit–explicit finite elements in transient analysis*
  (Parts I–II), J. Appl. Mech. 45:371–378 — element partitions; Belytschko & Mullen
  (1978), IJNME 12 — nodal partitions.
- Olovsson, Simonsson & Unosson (2005), *Selective mass scaling for explicit finite
  element analyses*, IJNME 63(10):1436–1445 — SMS theory baseline (ADR-36/38).
- Fork anchors: `SRC/analysis/integrator/CentralDifferenceSMS.cpp` (cap warning
  :261–265 warn-and-proceed, self-report consume :232–234, report-only downgrade
  :141–148), `CentralDifferenceLadruno.cpp` (`-recompute` :111, pre-step dt_cr
  validity :335, `revertToLastStep` :420–436, uniform-Δt `Vhalf` update :574,
  pre-scaling `getCriticalTimeStep` :742–746), `ExplicitBathe.h:102–113` (driver
  pattern pointer + stale example name), `LadrunoMassScaling.h` (`minDtSelfReport`;
  per-element pencil computed but not dumpable).

## Implementation log

*(Route B only, when its oracle fires; Routes A/C stay here until a trigger fires.)*
