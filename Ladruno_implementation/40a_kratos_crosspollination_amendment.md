---
title: ADR-40 amendment — KRATOS cross-pollination + scope-completeness (adversarially reviewed)
project: Ladruno
status: proposed
priority: medium
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - implementation
  - adr
  - performance
  - amendment
  - kratos
  - cross-pollination
  - adversarial-review
---

# ADR-40 amendment — KRATOS cross-pollination + scope-completeness

> **Amendment, not replacement.** This sharpens [[40_ladruno_performance_adr]] with (a) a
> source-verified cross-pollination crosswalk against KRATOS Multiphysics, (b) **six** hidden
> per-step cost centers the original program does not scope, and (c) corrections to several
> decisions that an adversarial review found over-claimed. Every claim below was source-checked
> in both trees during a 38-agent adversarial-review workflow (5 hostile lenses → independent
> verify → synthesize), run to completion across two resumes past a session limit. **One critique
> was refuted (C14 — the amendment already subordinated it to the anti-goal gate); 27 findings
> upheld/partial with corrections.**
> The governing principles of ADR-40 are unchanged and binding here: **P1 measure-before-optimize**,
> **P2 reuse > reinvent**, the anti-goals, and **no gating on cross-tool numbers**.

## Provenance & verification status

- KRATOS repo audited: `C:\Users\nmb\Documents\Github\Kratos` (builder&solver, linear solvers,
  constitutive-law, geometry, explicit path).
- Fork claims re-verified at file:line in `SRC/`.
- **Verification complete.** All five lenses including *decision-soundness* landed (the workflow
  was resumed twice past a session limit). The decision-soundness pass added three corrections
  folded in below: the **C9 over-claim** (rank-11 implicit is *not* "largely done"; `MumpsParallelSOE`
  is inherited stock, and the ICNTL(7) premise is an open gate), the **C5/C6 complacency strip**
  ("parity" is not the bar — reframed as a reuse-prior negative result), and **two further hidden
  cost centers** (§3.5/§3.6). The only refuted critique was against C14 — *because* this amendment
  already handled it correctly.

---

## §1 — Corrected KRATOS reference citations (the crosswalk, now precise)

These are the cross-pollination references for ADR-40 ranks 4/5/8/10. All confirmed; line numbers
corrected where the first pass was wrong.

| Ref | ADR-40 rank | KRATOS reference (verified) | Correction made |
|---|---|---|---|
| Lazy tangent via per-call flag | rank 5 (`tang_flag`) | `small_displacement.cpp:148-152` sets `ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR` true/false from `CalculateStiffnessMatrixFlag`; flag defined `constitutive_law.h:98`. Stateless, per-call — the exact "tang_flag, NOT per-instance memoization" pattern. | **Citation repointed** from `base_solid_element.cpp:147-152` (that's `InitializeNonLinearIteration`; `BaseSolidElement::CalculateAll` is an abstract stub that throws) to the derived `small_displacement.cpp`. |
| Geometry/shape-fn cache | rank 4 (de-static) | Shape-function values & local gradients stored per `IntegrationMethod` in immutable `std::array` containers (`geometry_shape_function_container.h:82,86`), returned by `const&` (`:245-253`); surveyed displacement elements carry no function-local statics. | **Quantifier scoped** to surveyed elements (`small_displacement`, `base_solid_element`) — not asserted over all ~hundreds of KRATOS element classes. |
| Numeric-factorization-reuse gap | rank 8 | `eigen_direct_solver.h` **already splits** `Compute` (factorize, `:99`) from `Solve` (tri-solve, `:116`) and documents reuse intent (`:88`), but `Solve()` calls `Compute()` **unconditionally**, with no factored-flag skip; re-factorization fires **every `Solve()` = every Newton iteration**. `KeepSystemConstantDuringIterations==true` saves *assembly* (`BuildRHSAndSolve`) but **not** factorization. | **Strengthened**: KRATOS is a *tighter* mirror of fork rank-8 than first stated — the solve-side split already exists; only the skip-guard is missing. **No shippable reuse target exists in KRATOS** — rank 8 must be *built* (in-tree `MumpsSolver.cpp` job=2 is the reference), not ported. |
| Factor+solve lumped in profiler | rank 1 | `residualbased_block_builder_and_solver.h:530-534,642,683-687` wraps the whole solve in one opaque `Timer("Solve")` scope; the Eigen wrapper does `Compute`+`Solve` with no internal split. | Confirmed. **Reframed as a reuse-prior *negative result*:** KRATOS offers no finer timer to copy, so rank 1 must build its own scope — this removes an external pressure source, it does **not** retire rank 1, which stays gated on the fork's own Phase-0 data (not on KRATOS parity). |
| GeoMechanics defaults to AMGCL | (informs C7/C10) | `geomechanics_solver.py:115-125` base default = `amgcl`/`amg`/`ilu0`/`gmres`; the **U-Pw subclass** `geomechanics_U_Pw_solver.py:79-89` overrides to the **same** amgcl default. A stage's `ProjectParameters` may override. | Confirmed; the subclass override **strengthens** the claim (it resolves to amgcl, not away from it). |

---

## §2 — Instrument-vs-optimize discipline (keeps the amendment inside P1)

Every amendment item is one of two kinds. **INSTRUMENT** work is Phase-0-safe now; **OPTIMIZE**
work inherits ADR-40's Phase-0 gate and may not start ahead of a measured dominance result.

| INSTRUMENT-now (Phase-0-safe) | OPTIMIZE (Phase-0-gated) |
|---|---|
| Reference citations (§1) | Atomic-scatter for rank-7 implicit half (§4, C2) |
| Profiler **scopes** for the six hidden costs (§3) | Rank-4/7 *performance* promotion (§4, C16) |
| Contact profiler **pre-wiring** (§4, C11) | Regime-split SIMD (§4, C14) |
| KRATOS anchor wiring, one matched model (§4, C7) | Iterative-AMG evaluation (§4, C10) |
| rank-4 de-static **on thread-safety grounds only** (§4) | rank-4 de-static *for its speed payoff* |

---

## §3 — Six hidden per-step cost centers ADR-40 does not scope (NEW)

The scope-completeness lens found six real per-step costs buried inside scopes ADR-40 treats as
cheap. Each needs a **rank-1 profiler sub-scope** and a **rank-3 bench**. These are the most
material additions in this amendment.

1. **LadrunoBrick `-formulation eas` inner Newton + static condensation.**
   `formEAStrue` runs a per-element nested Newton (`while` loop `LadrunoBrick.cpp:2660-2698`,
   re-evaluating all 8 GP constitutive models and solving a 9×9 `Kaa.Solve` `:2695` each
   iteration to a data-dependent count), then static condensation `K*=Kdd−Kda·Kaa⁻¹·Kad`
   (`:2733-2734`). It is reached from both `getTangentStiff` (`:447`) and the residual path
   (`:803,:962`) — i.e. **per element per global iteration, hidden inside `formTangent`**.
   → rank-1 sub-scope inside `formEAStrue` (inner-iteration count + the two `Kaa` solves);
   rank-3 bench `eas` vs `ssp`/`std`. **Also the worst case for §4-C2 atomic-scatter:** the
   data-dependent EAS iteration count makes a static partition straggle.

2. **LadrunoIMKBeam / IMKBeam2d hinge Newton.** `ladrunoIMKSolveAxis`
   (`LadrunoIMKHinge.h:52-148`) runs a 2×2 internal Newton up to **25 iterations**
   per bending axis with softening (`k_t<0`) resolution + flexibility condensation, every form
   pass — twice per 3D element (Mz/My, `LadrunoIMKBeam.cpp:226,233`), once for 2D
   (`LadrunoIMKBeam2d.cpp:225`). **For a bare concentrated-plasticity frame this, not brick
   geometry, is the dominant element cost.** ADR-40 has zero IMK mention, and Phase-0 bench (a)
   is *forceBeamColumn* (vanilla), which never exercises it.
   → rank-1 scope around `ladrunoIMKSolveAxis`; add an **IMK moment-frame bench**.

3. **LadrunoArcLength double solve — and the factor-vs-solve scope must move to the SOE layer.**
   `LadrunoArcLength` does a predictor solve in `newStep` (`:316-322`, **unscoped** — the
   `NewtonRaphson("linearSolve")` scope ADR-40 relies on, `:66`, never wraps it) and a second
   corrector solve in `update` (`:378-380`, **mis-bucketed** under the `"update"` scope). An
   arc-length run (snap-through; the von Mises truss the fork validates against) thus
   **under-reports solve cost ~2×**. → the rank-1 factor-vs-solve scope ADR-40 plans (`:70-72`)
   must live in the **SOE layer** (`SRC/system_of_eqn/…solve()`), not only in the equiSolnAlgo,
   so integrator-driven solves are captured. **rank-8 pays double here** — both solves reuse the
   one factorization from `newStep:316`.

4. **LadrunoDynamicRelaxation fictitious-mass rebuild.** `buildGershgorinDiagonal`
   (`LadrunoFictitiousMass.h:94`) calls `getTangentStiff()` over **all** FE_Elements — a full
   element-tangent assembly — and it is **not one-time**: `-recompute N` (`:325-331`) and,
   on the **default path** (`autoRefresh=true`, `dampMode=0`), an auto-rebuild at **every
   kinetic-energy peak** (`:423-424`). So a default quasi-static DR run carries recurring
   `formTangent`-equivalent passes inside the "explicit" driver. (The in-code "built once"
   comment `:197` and LEDGER row 58 are stale.)
   → rank-1 scope around `buildFictitiousMass`; add a **DR quasi-static bench**.

5. **LogStrain finite-strain kinematic tax.** Every `LogStrain` constitutive call runs ≥3 cyclic-
   Jacobi 3×3 eigensolves (`LogStrainKernel.h:96-130`, 50 sweeps) + degeneracy branching
   (`:169-209`) + a full 4th-order spatial-tangent push-forward (`spatial_tangent_full :279-305`)
   wrapping the inner 6×6 (`LogStrainNDMaterial.cpp:141-208`). **rank-5's `tang_flag` optimizes
   only the inner 6×6 and cannot touch this fixed per-GP overhead** — so rank-5's measured payoff
   must be judged against a cost it can't remove. → rank-1 scope around the eigensolve + push-forward;
   rank-3 finite-vs-small-strain bench on the same model.

6. **LadrunoRecorder per-step HDF5 I/O.** `LadrunoRecorder` grows a chunked dataset one slab
   **per step** through `shuffle`+`deflate(level 4)` (`Ladruno_Hdf5.h:338-340`, `appendSlab3d`),
   with per-step in-place envelope rewrites + H5 flush (`LadrunoRecorder.cpp:401-405`) and
   one-file-per-partition-per-rank fan-out (`:489-503`). At MPI scale this can make the **C7
   cross-tool anchor silently I/O/compression-bound**. → rank-1 scope separating deflate+shuffle
   from the envelope rewrite/flush; rank-3 recorders-on-vs-off (and envelope-mode) bench variant;
   treat the hard-coded deflate level + flush cadence as *measured* levers evaluated before ranks 4-7.

> Together with the SMS-PCG (C12) and constraint-projection (C13) costs already noted, the
> "explicit step" is **not** the textbook element-loop + nodal-divide ADR-40 rank-1 assumes.
> Phase-0's explicit profiler must attribute: element loop, SMS-PCG, projection, **EAS inner
> Newton, IMK hinge Newton, DR mass-rebuild, LogStrain eigensolve tax**, contact, recorder I/O.

---

## §4 — Corrected / narrowed decisions

- **C2 — atomic-scatter for rank-7 implicit half → demote to "feasibility evidence + phase-2.5
  benchmark candidate."** Corrections: (i) the KRATOS assembly loop is a **raw OpenMP
  `#pragma omp parallel firstprivate(...)` + `omp for schedule(guided,512)`**
  (`residualbased_block_builder_and_solver.h:228-263`), **not** `block_for_each` (that's used for
  DOF/graph/sparsity setup only); the per-entry `AtomicAdd`-into-CSR
  (`Assemble`→`AssembleRowContribution` `:1543-1637`, `atomic_utilities.h`) with **no coloring**
  is the real, load-bearing reference. (ii) **Precondition:** atomic-scatter needs a frozen
  symbolic pattern + a directly-addressable value array → **couple to rank-10**; it does **not**
  apply to band/profile/skyline SOEs (already an anti-goal lane). (iii) **It needs no new SOE
  class** — OpenSees' default `SparseGenColLinSOE` already exposes a stable CSC value buffer `A[k]`
  indexed by `colStartA/rowA` (`SparseGenColLinSOE.cpp:308-363`); the adaptation is a thread-safe
  atomic-add path replacing the serial `addA()` linear search. (iv) Respect **P5/P1**: this stays
  a phase-2.5 item, gated on Phase-0 showing assembly >40% and after explicit v1 lands; **atomic
  contention vs colored-scatter throughput is itself a measurement** (EAS/IMK load imbalance is
  the hazard), not a settled win.

- **C10 — iterative AMG "escape hatch" → REWRITE.** The original anti-goal excludes AMG for
  **convergence** ("fails on the fork's ill-conditioned softening tangents", ADR-40 `:190-192`),
  **not memory**. The "memory-bound escape hatch" framing misattributed the rationale and is
  internally incoherent (AMG can't converge on the very softening tangents a fill-wall would fire
  on). Restate: *AMG remains an anti-goal for the softening-tangent lane. IF a future scaling
  bench MEASURES a MUMPS-fill memory wall on a **well-conditioned large-elastic** model (the only
  regime where AMG converges), THEN KRATOS AMGCL is the reference for a **demand-gated sub-ADR**
  (per P8) — not an in-scope item.*

- **C11 — contact → reframe as PRE-WIRING, correct the cost model.** Shipped contact is
  **P1 plumbing + P2a only**: a single rigid analytical **frictionless** penalty plane
  (`gap=n·(X+u−p0)`, `resid=−kₙ⟨−g⟩₊n`; `LadrunoContactFE.cpp:60-131`). The broad-phase
  spatial hash (ADR-39 **P2.5**), deformable narrow-phase closest-point (**P2/P3**), and IMPL-EX
  Coulomb friction (**P3**) are **deferred**; mortar+ALM is the unstarted ADR-41 draft. So C11's
  cost model is the **future** design. Keep the recommendation but reframe: **pre-wire** rank-1
  scopes and scaffold a rank-3 contact bench *ahead of* those phases landing — instrument, then
  rank when ADR-39 lands the data structures.

- **C16 — rank-4/7 "promote to primary" → promote the MEASUREMENT, not the optimization.**
  "CFL-bound" (number of steps) does **not** entail "element-loop-bound" (where time goes within
  a step) — a diagonal-mass step can be dominated by SMS-PCG (C12), projection (C13), contact
  (C11), or the §3 hidden costs. ADR-40 already gates rank 4 (`:121`) and rank 7
  (`:124`, `>40% ele`). Restate: *the explicit-wave bench (item (c), `:158`) must publish the
  per-step split; that split — not the CFL argument — decides whether ranks 4/7 pay.*
  **Carve-out:** rank-4 de-static may proceed **now on thread-safety/correctness grounds** (the
  `LadrunoBrick` data race, `:90`, hard predecessor of rank 7) — its **performance** payoff stays
  gated. Also drop the unqualified "no solve" for the SMS path (C12 establishes a PCG).

- **C7 — KRATOS anchor → matched-physics, supplements (not replaces) the Abaqus arm.** The
  anchor is valid but must be a **dry large dynamic SOLID** model
  (StructuralMechanics / explicit central-difference), **not** GeoMechanics U-Pw (different
  physics per C8, parked per C15). It fills the cross-tool slot the license-blocked Abaqus arm
  can't (`:163-167`), **inheriting `:162` "recorded once, never a gate" and the `:203-205`
  anti-goal by reference.** Effort cap: **one** matched model, recorded once — not a CI number.

- **C14 — regime-split SIMD → hypothesis to TEST, not a standing exception.** ADR-40 v1.2 already
  conditionally opens GPU constitutive-batching gated on Phase-0 dominance + roofline. Keep C14 as
  *"test with the Phase-0 roofline (low-divergence elastic-wave kernel vs branchy softening)"* —
  not a pre-authorization of SIMD work the ADR defers.

- **C9 — rank-11 implicit is NOT "largely done" → downgrade (decision-soundness, HIGH).** The
  distributed direct-solve *primitive* exists, but `MumpsParallelSOE`/`MumpsParallelSolver` are
  **inherited stock OpenSees** (upstream CVS header; plain registration at
  `OpenSeesCommands.cpp:4172-4184`), **not fork-authored**. Remaining implicit work is **four** open
  items, not two: (a) ParMETIS partition quality vs serial-METIS-on-rank-0; (b) the explicit-diagonal
  parallel numberer (ADR-30 v1); (c) **verify ADR-30's ICNTL(7) ordering-irrelevance premise** —
  still open per ADR-40 Q3 (`:229-232`); *gate, do not assume*; (d) **measure** distributed
  correctness / load-balance / scaling of the inherited MUMPS path on the fork's lane (P1). Inherit
  ADR-30's gate by reference; do not assert completion.

- **C5/C6 — strip the complacency coda (decision-soundness, MEDIUM).** "Parity with KRATOS" is not
  ADR-40's success criterion (P1: gated on the fork's own Phase-0 dominance). Frame C5/C6 as a
  **reuse-prior *negative result***: KRATOS also re-factorizes every solve and also lumps factor +
  tri-solve into one timer, so there is **no shippable reference to port** for rank 1 or rank 8 —
  both must be *built* (in-tree `MumpsSolver.cpp` job=2 is the rank-8 reference). Informational,
  **not** a reason to deprioritize either.

---

## §5 — Caveats folded into existing claims

- **C8** — implicit Newmark/HHT integrators are vanilla and un-reimplemented, **but the implicit
  path is exercised**: the coupling elements (`LadrunoEmbeddedRebar/Node 33005/33006`,
  `RBE2/RBE3 33011/33012`) ship `getDamp`/`getRayleighDampingForces` overrides + a Newmark smoke
  test to fix a first-implicit-step OOB crash. Say "vanilla integrators; touched only via
  element-side damping overrides," not "inert."
- **C8/C15** — the fork is dry at the **continuum** level (`LadrunoBrick getNumDOF()==24`), but
  `LadrunoEmbeddedNode` ships an **experimental, flag-gated UP pressure-DOF tie** (`-pressure`,
  `ndf≥ndm+1`; LEDGER row 89 / PR #205). Note it so "u-p out of scope" isn't read as "zero
  pressure-DOF handling." Does not change any rank/anti-goal.

---

## §6 — Net edits to ADR-40

1. **Rank 1 (profiler):** add sub-scopes for EAS inner Newton, IMK hinge Newton, DR mass-rebuild,
   **LogStrain eigensolve tax, recorder HDF5 compression/flush**, SMS-PCG, constraint projection,
   contact (pre-wired); **move the factor-vs-solve scope to the SOE layer** so integrator-driven
   solves (arc-length, DR) are captured. Add the explicit-iterative/projection axis to the Phase-0
   dominance question (`:216-218`).
2. **Rank 3 (benches):** add an **IMK moment-frame** bench, an **arc-length snap-through** bench,
   a **DR quasi-static** bench, an **`eas`-vs-`ssp`** bench, a **finite-vs-small-strain** bench, a
   **recorders-on-vs-off** variant, and a **dry large-dynamic-solid** model matched to KRATOS as
   the cross-tool anchor.
3. **Rank 7 / P5:** record atomic-scatter-into-`SparseGenColLinSOE` as the phase-2.5 implicit
   candidate (couple to rank-10; benchmark vs colored-scatter; EAS/IMK load-imbalance caveat).
4. **Decisions:** rewrite the C10 AMG carve-out (convergence, demand-gated sub-ADR); restate C16
   as measurement-promotion with the rank-4 thread-safety carve-out; bound the C7 anchor;
   **downgrade the C9 rank-11 "largely done" claim to a four-item open list** (inherit ADR-30's
   gate by reference); **strip the C5/C6 parity coda** (reuse-prior negative result).
5. **Anti-goals:** unchanged in force; C10/C14 clarified to sit *inside* them, not beside them.

## Implementation log

- **2026-06-22 — adversarial-review pass (complete).** 38-agent workflow (5 hostile lenses →
  independent verify → synthesize), run to completion across two resumes past a session limit.
  **27 findings upheld/partial, 1 refuted (the C14 critique — already handled).** All KRATOS + fork
  claims re-checked at file:line. The decision-soundness lens added the C9 downgrade, the C5/C6
  reuse-prior reframe, and two further hidden cost centers (LogStrain tax, recorder I/O). Ready for
  review toward `accepted`.
