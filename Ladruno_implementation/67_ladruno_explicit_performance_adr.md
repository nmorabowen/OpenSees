---
title: Explicit-integrator performance — cheaper steps & cheaper Δt setup
project: Ladruno
status: scoping (design-capture, no code yet; P1 accepted, P2 conditional, P3 deferred)
priority: medium
owner: nmora
tags:
  - adr
  - integrator
  - explicit
  - performance
  - mass-scaling
---

# ADR-67 — Explicit-integrator performance

> ADR-67. Status: **SCOPING — design-capture, no code.** Next free slot after ADR-66
> (solid-shell). Companion axis to **ADR-65** (explicit Δt *strategies* = **fewer steps**);
> this ADR is the orthogonal axis = **cheaper steps + cheaper Δt setup**, i.e. the
> integrator-side work that is *not* the element force loop.
> Family: ADR-05 (CentralDifferenceLadruno) · ADR-35/36/38 (HRZ + selective/consistent
> mass scaling — the consistent PCG this ADR optimizes) · ADR-52 (ExplicitBathe flag
> class) · ADR-40 (general performance/profiler) · ADR-49/49a (integrator study) ·
> ADR-30 (explicit constraint projection, in the per-step path) · ADR-65 (Δt strategies —
> the fewer-steps lever this ADR defers to).
> Provenance: 2026-07-04/05 skill-guided review (the `explicit-dynamics` skill's
> `advanced_methods_landscape.md` "what to pursue" ranking + `stability_and_timestep.md`
> + `mass_scaling.md`), every claim cross-checked against current source.
> Revision 2026-07-05: adversarial re-review (F1–F6) corrected the P2.1 closed-form,
> re-scoped P1.1 to the whole `refineAccel` body, fixed the P3.1 warm-start contract,
> moved P2.3 off the byte-identical gate, and widened the P1.2 cache-invalidation trigger.
> Revision 2026-07-05b: independent (Fable) cross-review **DROPPED the original P1.2** — its
> premise (a periodic initial-stiffness eigensolve) is unreachable, `-recompute` force-upgrades
> to `-tangent` (`CentralDifferenceLadruno.cpp:119`) — and replaced it with **P1.2′** (massless-
> element short-circuit, byte-identical, reachable); corrected P2.1 (the element pencil is
> **fixity-blind**, closed form is `k(1/m₁+1/m₂)`, its massless-spring target belongs to P1.2′);
> added a required mid-run-topology-change test and the `updateParameter`/`updateMaterialStage`
> staleness channel.

## Context — the performance model

Explicit wall-time = **(number of steps) × (cost per step)**. Two independent levers:

- **Fewer steps** = larger stable Δt → **ADR-65** (SMS / subcycle / IMEX / KR-α). *Out of
  scope here* except to point at it.
- **Cheaper step** + **cheaper Δt setup** = *this ADR*.

**The dominant per-step cost is the element internal-force assembly** (`formUnbalance` — the
element loop), NOT the integrator. That is `opensees-performance` territory (SoA batching /
vectorization / GPU, ADR-40). The *integrator* can only (a) cut the step count [ADR-65], or
(b) remove its own redundant work. This ADR is strictly (b).

## What is already optimal (do NOT touch)

- **Lumped step** (`CentralDifferenceLadruno::update`, `CentralDifferenceLadruno.cpp:646-687`):
  O(N) diagonal `M⁻¹` (in the `Linear` algorithm's SOE solve), O(N) velocity/`setResponse`,
  O(N) KE-breaker dot. No hidden O(N²). Leave it.
- **The eigensolve already requests the largest eigenvalue only** — `DSYGVX RANGE='I',
  IL=IU=n` (`CriticalTimeStep.cpp:107-110`), the ~2× cheaper symmetric-definite path with a
  `DGGEV` fallback. The obvious "don't compute all eigenvalues" win is already banked.
- **Element self-report short-circuit** already skips the eigensolve for elements that carry
  their own bound (`CriticalTimeStep.cpp:316-323`).
- Both SMS paths are MPI bit-exact (ADR-38 V5).

## Findings & phases (integrator-side, ranked by value × (1/risk))

### P1 — accepted (byte-identical or strictly-less-work, low risk)

**P1.1 — Pre-allocate the consistent-SMS PCG scratch (whole `refineAccel` body, not just
`consistentPCG`).**
`consistentPCG` allocates **five full `neq`-length `Vector`s every step** (`Ax, z, p, Ap,
res` — `LadrunoMassScaling.h:559-562`), called once per step from `refineAccel`. But the
churn is larger than those five: the caller `applyConsistentRefine` **also** allocates
`Vector r(neq)` every step (`LadrunoConsistentRefine.h:149`), and the **MPI path
additionally rebuilds** `std::vector<double> w` and `Vector ones(neq)` **every step**
(`LadrunoConsistentRefine.h:87-91`). A hoist scoped to `consistentPCG` alone leaves ≥1
(serial) / ≥3 (MPI) `neq`-allocations per step on the table. For a 100k-DOF model the
five-vector slice alone is ~4 MB malloc/free *per step*, thousands of times — pure
allocator churn. **Fix:** hoist the *entire* `refineAccel` scratch set (the 5 PCG vectors
**plus** the caller's `r`, and the MPI `w`/`ones`) to reusable buffers owned by the
consistent integrator (sized in `domainChanged`, reused every step). The MPI `w` rebuild is
the one judgment call — the current code rebuilds it each step for robustness to
re-partition; cache it but invalidate on `domainChanged` (topology/partition change).
**Result: byte-identical.** Highest-ROI, lowest-risk item; the single cleanest.

**P1.2 — ~~Cache per-element λ_max in initial-stiffness mode~~ — DROPPED (premise unreachable).**
The original P1.2 assumed `-cfl` **without** `-tangent` still runs a *periodic* `-recompute N`
eigensolve worth caching away. **It does not.** All three parsers **force
`cflUseTangent = true` the moment `-recompute` is given** (`CentralDifferenceLadruno.cpp:119` —
source comment: *"recomputing the initial stiffness every N steps is pointless"*;
`ExplicitBathe.cpp:210,375`), and the SMS parsers downgrade `-recompute` to report-only. In
pure initial-stiffness `-cfl` the eigensolve already runs **exactly once**, in `domainChanged`
(`CentralDifferenceLadruno.cpp:338`); there is **no periodic initial-stiffness recompute to
skip**, and "-recompute without -tangent" cannot be produced from the script interface (the
flag auto-upgrades). The only residual — reusing `λ_max` across successive `domainChanged`
events — is both tiny and exactly where the stale-cache instability hazard lives, and the
safety-required full invalidation on `domainChanged` makes it a literal no-op. **Superseded by
P1.2′.**

**P1.2′ — Massless-element short-circuit before `DSYGVX`/`DGGEV` (the reachable byte-identical
win).** An element with all-zero lumped mass (`mdiag[i] ≤ 0 ∀i`) — e.g. a zeroLength
interface/contact spring whose inertia is carried by its nodes' *other* elements — fails the
`mPositive` check (`CriticalTimeStep.cpp:92-94`) and falls into the **most expensive** path,
`DGGEV` (`:140-188`), only to have **every** eigenpair rejected (`β=0` ⇒ `betaTol=0` ⇒ all
skipped) and return `-1` ("no bound"). **Fix:** pre-check `all mdiag ≤ 0 → return -1`, skipping
the full general eigensolve. **Byte-identical** (same `-1`, same "contributes no `Δt` bound"
semantics — a massless element has no inertial timescale), and it targets *exactly* the tiny
interface springs P2.1 was chasing — but correctly (they have **no finite** `λ_max`, so P2.1's
"closed form" was never applicable to them; see P2.1). Reachable, safe, and needs no config the
parser forbids.

### P2 — conditional (accuracy-neutral fast-paths; small numerical surface, gate on batteries)

**P2.1 — Closed-form `λ_max` fast-path for simple elements.**
Truss / zeroLength / spring / 1-DOF-pencil elements have `λ_max` analytically — skip
`DSYGVX` and its ~9 per-element `new[]/delete[]` (`elementLambdaMax`, `CriticalTimeStep.cpp:
256-269`; `maxGeneralizedEigenvalue` allocs `CriticalTimeStep.cpp:97-134`).

**Two facts reshape this item (2026-07-05 re-review) — read before implementing.**

*(a) The element pencil is fixity-blind, so the closed form is NOT `λ_max = k/m` and NOT
fixity-dependent.* `elementLambdaMax` builds the pencil from **only** `ele->getInitialStiff()`
and the lumped element mass (`CriticalTimeStep.cpp:252-265`) — boundary conditions never enter
it (grounding is applied at the DOF/assembly layer, not the element). So a grounded spring and a
free–free truss produce the **same** element number. For a 2-node axial pencil
`K=k[[1,-1],[-1,1]]`, `M=diag(m₁,m₂)` the eigenvalues are `{0, k(1/m₁+1/m₂)}`, i.e.
**`λ_max = k(1/m₁+1/m₂)`** (= `4k/m_tot` for equal split; `Δt_cr = 2/ω_max = ℓ/c`, the
`stability_and_timestep.md` element check). Coding the naive `λ_max = k/m_tot` overestimates
the step by 2× → **unconditionally unstable**. The real variation axes are **unequal nodal
lumping** and **zero-mass DOFs**, plus which DOF *directions* the zeroLength material couples
(arbitrary per assigned material) — *not* grounding. The equivalence-gate mesh must vary
**nodal-mass asymmetry and include zero-mass DOFs**, not fixity (a grounded-vs-free-free mesh
tests the same number twice — vacuous).

*(b) P2.1's own headline population — tiny interface/contact springs — is mostly massless at
the element level, so it has NO finite `λ_max`.* Those springs carry no element mass (inertia
lives on their nodes via other elements) ⇒ they correctly return `-1` (no bound), handled far
more safely by **P1.2′** (the massless short-circuit) than by any closed form. So P2.1's
finite-`λ` closed form applies only to the **mass-bearing** simple elements (a truss with its
own `ρAL`), a much smaller population than the ADR first claimed.

Net: **P2.1 is downgraded** — P1.2′ takes the massless springs (byte-identical, safe); the
finite closed form is a narrow, must-gate-hard optimization for mass-bearing 1-DOF-pencil
elements only, worth it only if a profile shows their eigensolve setup actually costs. The
equivalence gate (closed form == eigensolve, to rounding, over a mass-asymmetric + zero-mass
mesh) is load-bearing, not a formality — the wrong formula is a silent stability hazard.

**P2.2 — Reuse eigensolve scratch under `-tangent`.**
`elementLambdaMax` mallocs `K_data`, `M_data` + 7 LAPACK arrays *per element per call*. Under
`-tangent` (re-eigensolved every step for state-dependent tangent) reuse **a single
max-`n`-sized scratch set across the element loop** — the loop is serial (`while` at
`CriticalTimeStep.cpp:307`), so this is one buffer sized to the largest element, not
thread-local storage (there is no threading here to be local to). Byte-identical.

**P2.3 — `consistentMatVec`: precompute lumped diagonal.**
`consistentMatVec` divides `x(i)/diagMinv[i]` (`LadrunoMassScaling.h:528-529`) on *every*
matvec of *every* PCG iteration. Precompute `mlump[i] = 1/diagMinv[i]` once per step and
multiply. **NOT byte-identical** — `x*(1/d)` is two roundings vs one `x/d`, so the PCG
follows a slightly different floating-point trajectory and the last bits of the recorder
output move; this item is gated by the **equivalence/energy-band** gate (with P2.1), **not**
the byte-identical gate. Scope note: the division sits in the O(`neq`) diagonal-scatter
loop, **not** the O(Σ`nl²`) element-block inner loop that dominates a matvec — so the saving
is a per-`neq` division, real but modest; the block matmul is the hot path. Low priority
relative to P1.

**P2.4 — Optional cheap `-cfl fast` sizing (ℓ/c bound).**
Belytschko's element-size/wave-speed bound `Δt_e ≈ ℓ_e/c_e` (`stability_and_timestep.md`)
is a conservative estimate that skips the eigensolve for continuum elements. Offer as an
opt-in `-cfl fast` mode: far cheaper setup, slightly conservative Δt. **Not** byte-identical
(different, conservative Δt) → opt-in flag, clearly documented as an estimate.

### P3 — deferred (touch accuracy; need the ADR-38 validation battery to move)

**P3.1 — Warm-start / adaptive-tolerance the consistent PCG.**
`consistentPCG` re-seeds from the diagonal solve each step and drives to `tol=1e-10`
(`LadrunoMassScaling.h:557`), 3–21 iters (ADR-38). A **warm-start from the previous step's
converged `a`** and/or a looser explicit-appropriate `tol` could cut iterations — but both
move the solution, so gate behind the ADR-38 f₁-preservation + energy-closure batteries and
prove the drift is below the existing gate bands before shipping.

**Warm-start requires an API change, not just a re-seed — the current signature forbids it.**
`consistentPCG` overloads its `a` argument as three things at once
(`LadrunoMassScaling.h:550`): the RHS *source*, the initial *guess*, and the *output*. The
caller recovers the right-hand side from the **entry** value of `a`:
`r(i) = a(i) / Ainv[i]` (= `M_lump · a_entry`, `LadrunoConsistentRefine.h:151` serial,
`:96-98` MPI), relying on `a` entering as the fresh diagonal solve `M_lump⁻¹ r`. If you seed
`a` with the previous step's `a_prev = M̃⁻¹ r_prev`, the caller computes
`r = M_lump · a_prev ≠ r_current` and **corrupts the RHS**, not just the guess. A real
warm-start must decouple the three roles: pass `r` explicitly and add a *separate* guess
vector — a change to the shared serial+MPI `refineAccel` body (`applyConsistentRefine`,
`LadrunoConsistentRefine.h:57`) and its two call sites. This raises the effort/risk above a
tuning knob: **medium value, real risk, and a signature/contract change** — do it only with
the ADR-38 battery in hand and re-verify the np=1 == np=2 bit-exact after touching the shared
body. (Adaptive `tol` alone does *not* need the signature change — it is the cheaper half of
this item and can be tried first.)

## Fewer-steps levers (cross-ref only — belong to ADR-65 / the skill landscape)

Recorded so this ADR is self-contained; **not** owned here:
- **Consistent (Olovsson) SMS is the shipped "bigger Δt without period corruption" lever**
  (−0.17% f₁, ADR-38) — should be the *default* recommendation over lumped SMS.
- **Expose `HHTGeneralizedExplicit`** (vanilla, dissipative explicit, single ρ∞ knob) — the
  skill's #1 "cheapest real win": kills contact/hourglass chatter at ~1× cost, enabling a
  stable run at a larger Δt. Documentation + validation, *not* a build. (Candidate for an
  ADR-52 addendum or ADR-65 Route 0.)
- **Subcycling** = ADR-65 Route A (deferred; AVI / hardened clustered LTS, ADR-scale).
- **KR-α** = vanilla `KRAlphaExplicit`, reachable but unwired; dense O(N²) + softening-only;
  niche (moderate-DOF softening / RTHS), not large SSI/contact. ADR-65 Route (model-based).

## Decision

- **Do P1.1 + P1.2′ first** — byte-identical, self-contained, reachable. **P1.1** targets the
  path that actually costs integrator time (consistent-SMS per-step allocation); **P1.2′** (the
  massless-element short-circuit) skips the wasted `DGGEV` on interface springs. (The original
  **P1.2 λ_max-cache is DROPPED** — its premise, a periodic initial-stiffness eigensolve, does
  not exist on any reachable code path: `-recompute` force-upgrades to `-tangent`,
  `CentralDifferenceLadruno.cpp:119`.) Verify against the existing
  `test_centralDifferenceSMSConsistent_*`, `test_massScaling_validation*`, and CDL batteries
  with a **byte-identical assertion** (the no-op-preserving discipline of ADR-05/36/38: default
  paths must not move) — **plus a new mid-run `remove element` + continue case** so the P1.1
  buffer-lifetime invariant is actually exercised (see Validation).
- **P2 conditional** — **P2.2** (reuse eigensolve scratch, byte-identical) is the cheapest and
  can ride alongside P1. **P2.1 is downgraded** — its massless-spring target is taken (correctly)
  by P1.2′; the residual finite closed form (`k(1/m₁+1/m₂)`, fixity-*independent*) applies only
  to mass-bearing 1-DOF-pencil elements, and the wrong form is a silent stability hazard, so it
  needs the load-bearing equivalence gate — pick it up only if a profile shows those elements'
  eigensolve setup actually costs. **P2.3** is a modest per-`neq` win and is **not**
  byte-identical (equivalence gate). P2.4 only if a user wants the fast estimate.
- **P3 deferred** — only with the ADR-38 accuracy battery in hand.
- **Fewer-steps work stays in ADR-65**; the one cheap cross-over worth pulling forward is
  exposing `HHTGeneralizedExplicit` (separate small effort).

## Validation

- **Byte-identical gate** for P1.1/P1.2′/P2.2: default lumped + consistent runs must be
  bit-for-bit unchanged (compare recorder output pre/post on the SMS + CDL batteries).
  (P2.3 is **excluded** — `x*(1/d)` ≠ `x/d` bit-for-bit; see the equivalence gate.)
- **Lifetime-invariant gate (NEW, required for P1.1)**: a **mid-run `remove element` +
  continue** case on a consistent-SMS run — the existing SMS/CDL batteries are all static-mesh
  (`test_centralDifferenceSMSConsistent_*`, `test_massScaling_*` do **no** mid-run topology
  change), so without this the P1.1 buffer-resize / MPI-`w`-staleness invariant — the item's
  *only* real risk — ships untested. The byte-identical gate on a fixed mesh cannot see it.
- **Equivalence / energy-band gate** for P2.1 **and P2.3**: for P2.1, closed-form `λ_max`
  == eigensolve `λ_max` (to rounding) on a mesh that varies **nodal-mass asymmetry
  (`m₁≠m₂`) and includes zero-mass DOFs** (*not* grounding — the pencil is fixity-blind, so a
  grounded-vs-free-free mesh tests the same number twice), plus a brick; the governing
  element/tag must not change. For P2.3, the run must stay within the ADR-38 f₁-preservation +
  energy-closure bands (not bit-exact).
- **Perf gate**: profile (ADR-40) a large consistent-SMS run pre/post P1.1 (per-step alloc
  bytes → ~0) and a setup-heavy interface-spring model pre/post P1.2′/P2.2 (per-element
  eigensolve time; count of elements taking the massless short-circuit).
- MPI: re-confirm np=1 == np=2 bit-exact after touching `consistentPCG`/`consistentMatVec`.

## Gotchas (fold to LEDGER_quirks on ship)

- `consistentPCG`/`consistentMatVec` scratch is shared with the **distributed** PCG
  (`LadrunoConsistentRefine.h`, `consistentParPCG`) — any buffer hoist must keep both the
  serial and MPI paths correct (the serial path must stay the `w==null` reduction of the
  distributed one).
- Any per-element `λ_max` caching (P2.1, or a future cross-`domainChanged` scheme) must be
  **invalidated on `domainChanged`** (element add/remove, ADR-51) — a stale cache after a
  topology change is a silent **wrong-Δt → unstable-run** hazard, strictly worse than the
  recompute it saves.
- **`λ_max` staleness has a second channel beyond `domainChanged` and Rayleigh:**
  `updateParameter` on `E`/stiffness and `updateMaterialStage` (soil-stage switch — squarely
  this fork's SSI lane) change the **initial** stiffness *without* firing `domainChanged`. Any
  initial-stiffness `λ_max` cache (P2.1) can therefore go silently stale even when
  undamped-only. Either don't cache across these events or hook them explicitly.
- `-recompute` **cannot** be requested without `-tangent` (the parser force-upgrades,
  `CentralDifferenceLadruno.cpp:119`) — so there is no "initial-stiffness `-recompute`" path to
  special-case (this is why the original P1.2 was dropped). Any future skip that *does* touch
  the tangent-recompute path must gate on `!cflUseTangent`.

## Implementation log

- **2026-07-06 — lane-D drill-down closes the 40b Finding-3 items and adds ONE new design item
  ([[40b_phase0_dominance_report]] §lane-D addendum has the numbers).**
  - **P-NEW-1 (recipe, SHIPPED as documentation): `algorithm Linear -factorOnce` on plain
    constant-M fixed-Δt explicit runs.** The per-step tangent reassembly (16.2% of explicit
    wall) is vanilla `Linear` re-forming the constant lumped-mass effective tangent; the vanilla
    `-factorOnce` flag skips it. Measured live-wave A/B: **−17.9% wall, md5-bit-identical**
    (G-BYTE). HAZARD (quirks ledger): no domainChanged reset — never combine with `-reemit`
    contact / element removal / staged activation / variable Δt. A robust integrator-side
    tangent cache with domainChanged+Δt invalidation is the follow-up if those combinations
    ever need the win.
  - **P-NEW-2 (design item, NOT authorized): skip the second constitutive pass in
    `CentralDifferenceLadruno::update()`.** The step runs `updateDomain()` twice —
    `newStep` (trial advance; the real state determination) and `update()`
    (`:654-655`, full-step (u,v,a) snapshot push at **unchanged displacement**). For
    rate-independent materials the second pass recomputes identical constitutive state:
    **~22% of explicit wall**. A skip must gate on rate-dependence / velocity consumers
    (element damping reads v, LNVD publishes from the residual path, recorders read committed
    state) — G-BYTE only under a proven "no consumer of the second pass" condition. Scope as
    its own item with an adversarial gate before any code; combined ceiling with P-NEW-1
    ≈35-40% of the explicit step.
  - **newStep is exonerated:** its 22.5% is the leap-frog trial-state constitutive pass
    (elem.update 5.63 µs/ele), not integrator overhead. With full attribution the explicit lane
    is **~89% element work** — rank-7 threading (ADR-40) owns the big lever, this ADR owns the
    two items above.

- **2026-07-06 — P-NEW-1 PLANNED + ADVERSARIALLY GATED + BUILT: constant-mass tangent cache,
  default ON.** Design: `massTangentValid` latch owned by `CentralDifferenceLadruno`;
  `formTangent(statFlag)` override skips the whole assembly loop when latched (Diagonal SOE
  keeps its factored-in-place state — exactly the measured `-factorOnce` behavior); cleared in
  `domainChanged()` (the choke point for contact re-emission / removal / staged activation /
  SMS injection / every SOE resize), and — per the gate — in `setConstraintProjector()` and
  `revertToLastStep()`, keeping it inseparable from the ADR-30 `massBuilt` latch. Escape hatch
  `-noMassCache`. **Adversarial gate verdict (Opus, 10 attack items, all source-verified):
  default-ON G-BYTE defensible.** Findings folded: (1) latch-coupling guard (projector could
  read factored 1/M as mass under future reordering — the two extra clear sites); (5)
  `statusFlag` hygiene mirrored in the skip path; (6) "pure mass" rests on the inherited
  `getCFactor()==0` zero-scaling the base modal-damping branch — documented in the header, do
  not give CDL a nonzero cFactor without revisiting; (8) `updateParameter` on DENSITY is the
  one invisible event — documented exclusion + one-time runtime note when the deck declares
  parameters (`updateMaterialStage` dropped from the exclusion list: stiffness-only, irrelevant
  to a mass tangent); (9) flag transient/not serialized, all ctors init false; MPI path
  verified SAFER than serial (MPIDiagonalSOE is already factor-once by construction; no
  collective mismatch possible — collectives live in solve(), which every rank calls every
  step). SMS subclasses inherit; every SMS mass mutation verified bracketed by domainChanged.
  Gate also added the missing test: **revert-and-retry** (the only path exercising the
  top-priority guard). Zone-A gate file: `tests/test_cdl_mass_cache.py` (MC-1 quiet-path
  bit-identity, MC-2 remove-element invalidation, MC-3 failed-step revert+retry — all EXACT
  float equality vs `-noMassCache`). With this, the `-factorOnce` recipe advice is obsolete:
  every explicit CDL run gets the ~18% by default, contact/removal models included.

- **2026-07-06 — P-NEW-2 GATED + BUILT: `-commitSolveState` (OPT-IN, default OFF).** Skips the
  second element constitutive pass in `update()` (~22% of explicit wall, 40b lane-D) while
  KEEPING `setResponse(u, v_full, a)` — the **load-bearing invariant**: committed nodal output,
  recorders, KE, contact `-visc`, and the brick viscous-hourglass commit all read node state,
  which still gets v_{n+1}. Element state is committed as the solve actually used it (same
  strain u_{n+1}, lagged v_{n+1/2}). **Adversarial gate verdict (Opus, 10 items + inventory,
  all source-verified): SOUND as opt-in; default-ON NOT defensible** — the G-EQUIV class is
  non-empty in practice. **The exact G-EQUIV trigger set (gate inventory):** elements whose
  `update()` passes a velocity-derived rate to `setTrialStrain` — ZeroLength,
  CoupledZeroLength, ZeroLengthVG_HG, TwoNodeLink, Truss, Truss2, CorotTruss, CorotTruss2,
  MultipleShearSpring, MultipleNormalSpring — feeding rate-consuming materials
  (Viscous/ViscousDamper/ViscoelasticGap-type). Their committed state under the flag reflects
  v_{n+1/2} (arguably MORE consistent — the step's forces were computed there — but a results
  change). Safe-by-construction findings: TwoNodeLinkSection/TrussSection/ZeroLengthSection
  pass deformation only; ALL fork materials rate-independent (LadrunoJ2/Finite, Concrete3D —
  incl. `-eta` Duvaut-Lions: reads the SCALAR ops_Dt, identical in both passes, NOT node
  velocity — LogStrain, Staged*, Lemaitre); LadrunoBrick `update()` is pure-displacement
  (bulk viscosity + viscous hourglass read velocity at FORCE time / commit from nodes);
  contact adapters are FE_Elements outside `Domain::update()` entirely; RemoveRecorder
  criteria read trial displacement only; recorders run post-commit on committed node state;
  no-arg `updateDomain()` applies no loads; orthogonal to the P-NEW-1 mass-cache latch; DDM
  unreachable under explicit. One-time semantic note at first skip (gate P0 — never silent).
  Zone-A gate file: `tests/test_cdl_commit_solve_state.py` (CS-1 rate-independent
  bit-identity; CS-2 ZeroLength+Viscous AND Truss+Viscous characterization — delta REQUIRED
  nonzero and bounded; CS-3 brick `-bulkViscosity` + uri `-hourglass viscous` bit-identity —
  proves the setResponse invariant on the fork's own velocity-reading element; CS-4
  remove-element + revert/retry composition with P-NEW-1). Combined with P-NEW-1 the explicit
  step sheds ~35-40% on rate-independent decks (P-NEW-1 automatic; P-NEW-2 opt-in per deck).
