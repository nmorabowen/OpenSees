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

**P1.1 — Pre-allocate the consistent-SMS PCG scratch.**
`consistentPCG` allocates **five full `neq`-length `Vector`s every step** (`Ax, z, p, Ap,
res` — `LadrunoMassScaling.h:559-562`), called once per step from `refineAccel`. For a
100k-DOF model that is ~4 MB malloc/free *per step*, thousands of times — pure allocator
churn. **Fix:** hoist to reusable buffers owned by the consistent integrator (sized in
`domainChanged`, reused every step). **Result: byte-identical.** Highest-ROI, lowest-risk
item; the single cleanest.

**P1.2 — Cache per-element λ_max in initial-stiffness mode; make `-recompute` a no-op there.**
In `-cfl` **without** `-tangent`, each element's `λ_max` (initial `K`, lumped `M`) is
**constant** — but `-recompute N` re-runs the full `DSYGVX` on *every* element every N steps
(`computeCriticalTimeStep`, `CriticalTimeStep.cpp:288-367`), recomputing an unchanged number
at O(N_elem·n³). **Fix:** in initial-stiffness mode, compute once in `domainChanged`, cache,
and skip the periodic recompute entirely (or at minimum **warn that `-recompute` without
`-tangent` is a no-op**). **Result: byte-identical**, eliminates the entire `-recompute`
cost in the common case. (Under `-tangent` the recompute is genuinely needed — see P2.)

### P2 — conditional (accuracy-neutral fast-paths; small numerical surface, gate on batteries)

**P2.1 — Closed-form `λ_max` fast-path for simple elements.**
Truss / zeroLength / spring / 1-DOF-pencil elements have `λ_max = k/m` analytically — skip
`DSYGVX` and its ~9 per-element `new[]/delete[]` (`elementLambdaMax`, `CriticalTimeStep.cpp:
256-269`; `maxGeneralizedEigenvalue` allocs `CriticalTimeStep.cpp:97-134`). These are
*exactly* the tiny stiff interface/contact springs that dominate the models people mass-scale
and the biggest `-tangent`-mode setup cost. Detect by element class / `n==small` and validate
the closed form equals the eigensolve on a mixed mesh. Low risk (same number, cheaper route).

**P2.2 — Reuse eigensolve scratch under `-tangent`.**
`elementLambdaMax` mallocs `K_data`, `M_data` + 7 LAPACK arrays *per element per call*. Under
`-tangent` (re-eigensolved every step for state-dependent tangent) reuse thread-local scratch.
Byte-identical.

**P2.3 — `consistentMatVec` inner-loop: precompute lumped diagonal.**
`consistentMatVec` divides `x(i)/diagMinv[i]` in the hot loop (`LadrunoMassScaling.h:528-529`)
on *every* matvec of *every* PCG iteration. Precompute `mlump[i] = 1/diagMinv[i]` once per
step and multiply. Byte-identical (to rounding), removes a division from the innermost loop.

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
prove the drift is below the existing gate bands before shipping. Medium value, real risk.

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

- **Do P1.1 + P1.2 first** — byte-identical, self-contained, and they target the two paths
  that actually cost integrator time (consistent-SMS per-step allocation; `-recompute` /
  `-tangent` setup). Verify against the existing `test_centralDifferenceSMSConsistent_*`,
  `test_massScaling_validation*`, and CDL batteries with a **byte-identical assertion** (the
  no-op-preserving discipline of ADR-05/36/38: default paths must not move).
- **P2 conditional** — pick up P2.1/P2.3 alongside P1 if cheap; P2.4 only if a user wants the
  fast estimate.
- **P3 deferred** — only with the ADR-38 accuracy battery in hand.
- **Fewer-steps work stays in ADR-65**; the one cheap cross-over worth pulling forward is
  exposing `HHTGeneralizedExplicit` (separate small effort).

## Validation

- **Byte-identical gate** for P1.1/P1.2/P2.2/P2.3: default lumped + consistent runs must be
  bit-for-bit unchanged (compare recorder output pre/post on the SMS + CDL batteries).
- **Equivalence gate** for P2.1: closed-form `λ_max` == eigensolve `λ_max` (to rounding) on a
  mixed truss+brick mesh; the governing element/tag must not change.
- **Perf gate**: profile (ADR-40) a large consistent-SMS run pre/post P1.1 (per-step alloc
  bytes → ~0) and a `-recompute`/`-tangent` run pre/post P1.2/P2.2 (setup time).
- MPI: re-confirm np=1 == np=2 bit-exact after touching `consistentPCG`/`consistentMatVec`.

## Gotchas (fold to LEDGER_quirks on ship)

- `consistentPCG`/`consistentMatVec` scratch is shared with the **distributed** PCG
  (`LadrunoConsistentRefine.h`, `consistentParPCG`) — any buffer hoist must keep both the
  serial and MPI paths correct (the serial path must stay the `w==null` reduction of the
  distributed one).
- Per-element `λ_max` caching must be **invalidated on `domainChanged`** (element
  add/remove, ADR-51) — a stale cache after topology change is a silent wrong-Δt hazard.
- `-recompute` no-op-in-initial-stiffness must NOT silently change behavior under `-tangent`
  (where recompute is load-bearing) — gate the skip on `!cflUseTangent`.
