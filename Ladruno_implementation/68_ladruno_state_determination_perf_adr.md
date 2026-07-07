---
title: Element / state-determination performance — the kernel work-removal sub-ADR
project: Ladruno
status: scoping (design-capture, no code; Phase-0-gated per ADR-40 P1/P8)
priority: medium
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - adr
  - performance
  - state-determination
  - element
  - material
  - kernel
  - sub-adr
---

# ADR-68 — Element / state-determination performance

> ADR-68. Status: **SCOPING — design-capture, no code.** This is the **OPTIMIZE** sub-ADR that
> [[40_ladruno_performance_adr]] **P8** ("each cleared item spawns its own sub-ADR") calls for on
> the element / state-determination lane — the ranks 4/5/6 kernel work-removal items plus the six
> hidden per-step cost centers [[40a_kratos_crosspollination_amendment]] §3 surfaced. It is
> **strictly downstream of ADR-40 Phase 0**: no OPTIMIZE item here may start before the profiler
> attributes state-determination dominance *on the lane that item serves* (ADR-40 P1).
> Family: ADR-40 (perf program / Phase-0 gate + profiler) · ADR-40a (instrument scopes + KRATOS
> crosswalk + the six cost centers) · ADR-67 (integrator-side perf — the *other* axis; the
> element loop is out of its scope and is the whole of this one) · ADR-09 (LadrunoBrick) ·
> ADR-10/13/16 (J2 / uniaxial-J2 / finite-J2) · ADR-14 (IMK beam) · ADR-19/20 (brick EAS) ·
> ADR-31 (Concrete3D) · ADR-06 (profiler — the instrument).
> Provenance: 2026-07-05 review runway (ADR-67 adversarial review + Fable cross-review) that
> established the gate taxonomy below and confirmed the element loop — not the integrator — is
> where the time is; ADR-40a §3 for the cost-center file:line (38-agent source-checked).

## Context — why a dedicated sub-ADR, and why now

ADR-40 concluded, and ADR-67 conceded, that **the dominant per-step cost is the element
internal-force / state-determination loop**, not the integrator and (for the fork's lane) not the
solver. That loop is where the real speed is — but ADR-40 only carries it as **ranked table rows**
(ranks 4-6) and ADR-40a only **instruments** it (profiler sub-scopes + benches). Neither is a
*design*: neither pins, per kernel, the correctness gate the optimization must pass, the oracle
that enforces it, or the lane whose Phase-0 number authorizes it. This ADR is that design.

The organizing insight from the ADR-67 review: **these kernels are all the same shape** — redundant
per-Gauss-point / per-element iterative work (a tangent filled when unused, a geometry rebuilt every
form pass, an inner Newton run to tolerance every global iteration, an eigendecomposition redone
per call). The lever is uniformly *lazy / cached / warm-started per-GP work*. But they split sharply
on one axis that decides everything about how they ship — **whether the fix changes the numbers**.

## The gate taxonomy (the core contribution — apply per item, do not blanket)

ADR-40's backwards-compat note lumps all kernel refactors as "byte-identical to the 1e-12 oracle."
The ADR-67 review proved that is **wrong for a subset** and the distinction is load-bearing:

| Gate | Definition | Applies to | Enforcement |
|---|---|---|---|
| **G-BYTE** (byte-identical) | recorder output bit-for-bit unchanged | work-*skipping* that computes the same values (lazy tangent, cached geometry, short-circuits) | pre/post bit-compare on the item's lane bench |
| **G-EQUIV** (equivalence-band) | a *different* numerical result within a stated band | algorithm *replacement* (closed-form return vs nested Newton) — two roundings or a different converged point | a **kernel-specific oracle** + band; the 3D oracle cannot cover a plane-stress path |
| **G-OPT** (opt-in, documented) | intentionally different / conservative | approximations offered as a flag (`ℓ/c` sizing) | flag + doc; never a default |

**Mis-filing an item is a defect.** ADR-67 P2.3 was filed G-BYTE but is G-EQUIV (`x*(1/d)≠x/d`);
rank 5 is genuinely G-BYTE but rank 6 is G-EQUIV and needs an oracle the fork does not yet have.
Every item below carries its gate explicitly.

## Optimization targets (ranked within the lane; each Phase-0-gated on *its* lane)

### T1 — J2 stress-core + lazy `consistentTangent` via `tang_flag`  (= ADR-40 rank 5) · **G-BYTE**
`returnMap` fills `Dtan[6][6]` **unconditionally** — there is no `tang_flag` in the signature
(`LadrunoJ2Kernel.h:160-167`), so every call assembles the 6×6 even on the residual-only passes that
never read it; `returnMapDamaged` adds a *second* 6×6 (`:366,:479`). **Fix:** thread a per-call
`tang_flag` (KRATOS pattern `small_displacement.cpp:148-152` / `constitutive_law.h:98`, ADR-40a §1)
so the stress core runs alone when the tangent is unused; drop the redundant `M/Mp/normM` recompute.
**Gate G-BYTE** — when the tangent *is* requested it is bit-identical; skipping an unused output
cannot move results. Lane: **fiber frames + any J2 solid**. Guard: the stress path's use of `Dtan`
as elastic-tangent scratch (`LadrunoJ2Kernel.h:213`) must be preserved. **Ceiling caveat:** rank 5
optimizes only the inner 6×6 — it cannot touch the LogStrain per-GP tax (T5), so judge its payoff
against a cost it cannot remove (ADR-40a §3.5).

### T2 — J2 plane-stress projected return map  (= ADR-40 rank 6) · **G-EQUIV (needs a new oracle)**
The plane-stress / plate-fiber path enforces `σ₂₂=0` by an **outer Newton on ε₂₂**, `maxIt=25`
(`LadrunoJ2.cpp:304`) to rel-tol `1e-10·‖σ‖` (`:311`), each sweep running the **inner** 3D scalar
return `maxIter=50` (`LadrunoJ2Kernel.h:224`) and rebuilding a full 6×6 — worst case ~1250 scalar
steps + condensation per material point per iteration. **Fix:** replace with a plane-stress-projected
closed-form return (Simo & Hughes Box 3.3 + Box 3.4/§3.4.5; de Souza Neto §9.2.3/§9.4 **[acquire]**).
**Gate G-EQUIV** — the closed form converges to a *different point at ~1e-10* than the nested Newton;
it **cannot** be byte-identical, and the existing 3D `LadrunoJ2` 1e-12 numpy oracle
(`tests/ladrunoj2_reference.py`) does **not** cover the condensed path. **Blocker:** this item needs
a **new plane-stress/plate-fiber oracle** derived independently, plus a stated equivalence band,
*before* any code. Lane: **fiber shells / plate-fiber sections**. Highest value-per-cost of the J2
work, highest ceremony.

### T3 — LadrunoBrick shape-gradient cache + de-static  (= ADR-40 rank 4) · **G-BYTE**
`computeBasis()/shp3d()` rebuild the reference-config `Shape[...]` gradients **every** `formTangent`
(`LadrunoBrick.cpp:534` static scratch; recompute per call) instead of caching them at the reference
configuration; the same block is file-scope `static` (`:86`, `:524-536`) — a data race (T-thread
predecessor). **Fix:** cache reference-config gradients once; move the scratch off file-scope.
**Gate G-BYTE** (same gradients, computed once). **Carve-out (ADR-40a C16):** the de-static may
proceed **now on thread-safety/correctness grounds** (it is the hard predecessor of any element-loop
threading); its **performance** payoff stays Phase-0-gated. Lane: **3D solid / moderate-3D concrete**.

### T4 — LadrunoBrick `-formulation eas` inner Newton + `Kaa` reuse  (ADR-40a §3.1) · **G-BYTE / G-EQUIV split**
`formEAStrue` runs a per-element nested Newton (`LadrunoBrick.cpp:2660-2698`) re-evaluating all 8 GP
constitutive models and solving a 9×9 `Kaa.Solve` (`:2695`) every inner iteration, then condenses
`K* = Kdd − Kda·Kaa⁻¹·Kad` (`:2733-2734`) — reached per element **per global iteration** inside
`formTangent` (`:447`, `:803`, `:962`). **Fix options, different gates:** (a) reuse the `Kaa`
factorization across the condensation solves and warm-start the EAS iterate from the last global
iteration — **G-BYTE** if the converged α is unchanged; (b) cap/relax the inner tolerance —
**G-EQUIV**, needs a band. Lane: **EAS solid**. Note: the data-dependent iteration count makes this
the worst case for any static-partition element-loop threading (ADR-40a §3.1) — relevant to ADR-40
rank 7 scoping.

### T5 — LogStrain finite-strain per-GP tax  (ADR-40a §3.5) · **G-BYTE (cache) / needs its own study**
Every `LogStrain` call runs ≥3 cyclic-Jacobi 3×3 eigensolves (`LogStrainKernel.h:96-130`, 50 sweeps),
degeneracy branching (`:169-209`), and a 4th-order push-forward (`:279-305`) wrapping the inner 6×6.
This is **fixed per-GP overhead T1's `tang_flag` cannot remove.** **Fix candidates:** cache the
eigendecomposition between the stress and tangent calls of the same step; cheaper degeneracy path;
push-forward only when the tangent is requested (composes with T1). **Gate G-BYTE** for the caching
variants (same spectral data reused). Lane: **any finite-strain material**. Needs a small dedicated
study — flagged here so T1's measured payoff is not misattributed.

### T6 — LadrunoIMKBeam hinge Newton  (ADR-40a §3.2) · **G-BYTE (warm-start) / measure first**
`ladrunoIMKSolveAxis` (`LadrunoIMKHinge.h:52-148`) runs a 2×2 internal Newton up to **25 iterations**
per bending axis with softening resolution + flexibility condensation, **every form pass** — twice
per 3D element (`LadrunoIMKBeam.cpp:226,233`), once for 2D. **For a bare concentrated-plasticity
frame this — not brick geometry — is the dominant element cost, and ADR-40's Phase-0 bench (a) uses
vanilla `forceBeamColumn`, which never exercises it.** **Prerequisite:** ADR-40 needs an **IMK
moment-frame bench** (ADR-40a §3.2) or this lane is invisible. **Fix:** warm-start the hinge Newton
from the last committed state (**G-BYTE** if converged to the same tol). Lane: **concentrated-plasticity
frames** — the fork's IMK lane, entirely unscoped by ADR-40's benches.

## Phase-0 verdicts (2026-07-06 first measured pass — see [[40b_phase0_dominance_report]])

| Target | Verdict | One-line basis |
|---|---|---|
| T1 | **NOT authorized** | no lane is J2-tangent-bound (solid-lane J2 residual = 3.9% of step) |
| T2 | **DEMOTED** — ceiling ≈2% of wall | elastic-baseline diff on the shell lane; iteration-count inflation, not the material return, owns the cost. The G-EQUIV oracle ceremony is not worth it |
| T3 | **Cache KILLED as designed; de-static proceeds (2026-07-07 micro-drill)** | measured split of the 40.6 µs implicit tangent: geometry 2-3%, **B·D·Bᵀ ≈95%**, material copies ≈2% — cache ceiling ~0.7% of implicit wall vs 2.1 KB/ele G-BYTE. Explicit-only opt-in variant PARKED (~9-11% of explicit element work) behind T7 + rank 7. Solid-lane axis shifts to the B·D·Bᵀ kernel (upper-triangle symmetry ~2× ceiling, blocking/SIMD) + SOE factor reuse (rank 8/10) |
| **T7 (NEW, 2026-07-07)** — CDL inertia no-op skip · **G-BYTE** | **Authorized for design + gate** | under CDL the residual's `formInertiaTerms` runs at trial accel ≡ 0 (`Azero` solve): momentum = Σ N_j·a_j = 0 exactly, yet the block (incl. its own `computeBasis`+8×`shp3d` recompute) costs ~22-30% of the residual call ≈ **~11% of explicit element work — for exact zeros**. Fix: skip when all 24 nodal trial-accel components are exactly 0.0 (≈20 ns check). Gate note: skipped adds are `x += temp·0.0` — exact except the IEEE `-0.0 + (+0.0) = +0.0` sign edge; gate test must use exact-float comparison. Needs its own Opus adversarial pass before code (implicit-integrator passthrough, Rayleigh interaction, `mass.Zero()` side effect when `tangFlag==0`) |
| T4/T5 | Unmeasured | need the EAS + finite-strain benches; Lemaitre bbar-vs-eas dump is prior art |
| T6 | **DEMOTED — refuted by the `elem.update` instrument (2026-07-06 re-measure)** | hinge Newton = 0.62 µs/ele, ~14% of the update phase, total 33004 work ≈12% of step; the 35.4% band was integrator/domain update MACHINERY, not element work. Warm-start ceiling <1% — dropped. Lane E's real profile = solve ~30% + update machinery ~30% + per-step glue ~25% (small-model fixed-overhead regime) |

## T7 — SHIPPED (2026-07-07): CDL inertia no-op skip · G-BYTE

`LadrunoBrick::formInertiaTerms(tangFlag)` now returns immediately when
`tangFlag==0 && inertiaSkip` and **every** nodal trial-accel component is exactly
`0.0` (8 nodes × 3 dof, ~24 exact-float compares). Under CDL the residual is posed
at trial accel ≡ 0, so the skipped block contributes only zeros; the mass matrix
(`tangFlag==1`, from `getMass`/`addInertiaLoadToUnbalance`) is never skipped.
Escape: `element LadrunoBrick … -noInertiaSkip` (a transient, **unserialized**
per-element flag set via `setInertiaSkip()` — both settings are G-BYTE-identical,
so a broker-reconstructed element defaulting to skip-ON is correct).

**Adversarial gate (Opus) → G-BYTE DEFENSIBLE, and STRONGER than ==-identical.**
The reviewer refuted the design's own conceded `-0.0` edge: `resid` is
`Zero()`-primed and its **last** write in every formulation is `resid += bodyForce`
(or the finite inline body subtraction), which renormalizes any `-0.0` to `+0.0`
*before* `getResistingForceIncInertia` reads it — so the skipped `resid += temp·(+0.0)`
is a provable no-op with **no** sign-of-zero divergence. It is bit-identical, not
merely value-equal.
- **LOAD-BEARING INVARIANT (recorded per the gate):** the `resid.Zero()`-prime +
  trailing `resid += bodyForce` normalization is what makes the skip bit-identical.
  Do **not** optimize away the trailing bodyForce add when `bodyForce==0`, nor
  reorder `globalizeForce` after it — either weakens bit-identity to `==`-identity.
- `mass.Zero()` side effect: benign — `mass` is read only by `getMass`/
  `addInertiaLoadToUnbalance`, both of which call `formInertiaTerms(1)` first (which
  re-zeros); `Element::getRayleighDampingForces` uses the virtual `getMass()`, not
  the member. Implicit integrators pass through untouched (nonzero trial accel).

**Validation:** Zone-A `tests/test_brick_inertia_skip.py` — IS-1 CDL-wave
displacement bit-identity, **IS-2 signbit lock** (recorded nodal accel compared
byte-for-byte via `struct.pack`, catching any `-0.0/+0.0` a value-`==` gate would
miss — the gate-recommended assertion), IS-3 Newmark implicit passthrough, IS-4
Rayleigh `-alphaM` under CDL — **4/4 PASS**. Perf: lane-D A/B (`laneD_ab_bench.py`
+ new `ELE_FLAGS`/`NSTEPS` knobs, median-of-5 interleaved, 10³ LadrunoBrick+J2,
1000 CDL steps): **−8.2%** wall (21.93 s skip vs 23.89 s `-noInertiaSkip`), in the
predicted ~5–10% band; **uz bit-identical to 17 digits** across configs. Fork-owned
file (no vanilla footprint, no class tag).

## Phasing (inherits ADR-40 Phase 0 as a hard gate)

- **Phase 0 (ADR-40, prerequisite — not owned here):** publish the attributed per-phase breakdown
  **per lane**. This ADR adds one requirement to it: the bench set must cover **all four lanes** the
  targets serve — fiber frame (T1), fiber-shell/plate-fiber (T2), 3D solid/EAS (T3/T4),
  **concentrated-plasticity frame (T6)** and **finite-strain (T5)** — or those targets have no
  authorizing number. (T6 and T5 benches are the gap ADR-40a §3.2/§3.5 already flagged.)
- **Phase 1 — G-BYTE work, lane-gated (T1, T3-cache, T5-cache, T6-warmstart).** Ship the
  work-skipping / caching items whose lane Phase 0 shows state-determination-bound. Each rides its
  lane bench with a **bit-for-bit** assertion. Lowest risk; no new oracle.
- **Phase 2 — G-EQUIV work (T2, T4-tolerance).** Only after the **new plane-stress oracle** (T2) and
  a stated band exist. These change the numbers; they need the equivalence gate, not the byte gate.
- **T3 de-static** proceeds out-of-band **now** on thread-safety grounds (ADR-40a C16 carve-out);
  its speed credit is booked only under Phase 1.

## Decision

- **This ADR does not authorize any code.** It is the design that lets a spawned item ship correctly
  the moment its lane clears Phase 0. The measure-first discipline (ADR-40 P1) is binding.
- **Gate per item, never blanket** — T1/T3/T5/T6 are G-BYTE; T2 and the T4-tolerance variant are
  G-EQUIV and blocked on a kernel-specific oracle. Do not file a G-EQUIV item under the byte gate
  (the ADR-67 P2.3 defect).
- **Sequence by lane, not by rank number** — the dominant kernel differs per lane (fiber → T1;
  concentrated plasticity → T6; finite strain → T5; 3D solid → T3/T4). Pick the target whose lane
  Phase 0 says is state-determination-bound; drop the rest for that lane.
- **T6 and T5 are the coverage risk** — ADR-40's current benches (vanilla forceBeamColumn, small-strain
  J2) make the two lanes where an IMK/LogStrain kernel dominates *invisible*. Adding those two benches
  is Phase-0 work, prerequisite to this ADR delivering anything on those lanes.

## Validation

- **Per-item gate** as tabled: G-BYTE ⇒ bit-for-bit recorder compare on the item's lane bench;
  G-EQUIV ⇒ kernel oracle + published band; G-OPT ⇒ flag + doc.
- **New oracles required before code:** a plane-stress/plate-fiber J2 oracle (T2) — the 3D
  `tests/ladrunoj2_reference.py` 1e-12 oracle does not cover the condensed path.
- **Perf gate:** each item profiled (ADR-06 / ADR-40 rank 1) pre/post on its lane bench, reporting the
  attributed sub-scope ADR-40a §3 defined (EAS inner Newton, IMK hinge Newton, LogStrain eigensolve
  tax) — the payoff is the movement of *that* sub-scope, not end-to-end wall-clock alone.
- **Regression:** the existing `LadrunoJ2` 1e-12 oracle and `LadrunoConcrete3D` `run_tangent_gate`
  fixture stay byte-identical after every G-BYTE refactor (ADR-40 backwards-compat).

## Gotchas / risks (fold to LEDGER_quirks on ship)

- **The 3D oracle does not cover plane stress** (T2) — a G-EQUIV item validated against the wrong
  oracle ships an unverified numerical change. Build the plane-stress oracle *first*.
- **`tang_flag` must be per-call, never per-instance memoization** (ADR-40 anti-goal) — dirty-flag
  tangent caching fights the commit/setTrial cycle and the condensation path.
- **T1 cannot remove the LogStrain tax (T5)** and **T3 (geometry) is not the IMK/EAS cost** — do not
  let a payoff measured on one lane justify the work on another; each lane needs its own Phase-0 number.
- **T3 de-static is the thread-safety predecessor for ADR-40 rank 7, but not sufficient** — the
  assembly-layer `FE_Element`/`DOF_Group` class-wide static pools (`FE_Element.h:132-133`,
  `DOF_Group.h:150-151`) and the vanilla element library's file-scope statics race independently
  (ADR-40 rank-7 scope correction, 2026-07-05). Threading is a separate, larger audit.
- **Backwards-compat with LogStrain reuse** — T1's `tang_flag` must keep the default `returnMap`
  signature usable by `LogStrainNDMaterial` (finite-strain lift) or add an overload; do not break the
  existing caller.
