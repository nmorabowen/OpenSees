---
title: Performance program for the Ladruno fork (measurement-first)
project: Ladruno
status: proposed
priority: medium
owner: nmora
tags:
  - implementation
  - adr
  - performance
  - profiling
  - solver
  - parallel
  - explicit
  - program
---

# Performance program for the Ladruno fork (measurement-first)

> **Program / roadmap ADR.** Not a single feature — a **prioritized, phased performance
> program** for the fork, produced by an 11-agent panel discussion (5 lenses propose →
> cross-critique → synthesis; see Implementation log). The governing principle is
> **measure before you optimize**: there is today **no attributed phase breakdown and no
> fork-vs-Abaqus / fork-vs-stock-OpenSees number anywhere**, so every speed claim below is
> a *hypothesis gated on Phase 0*. The second principle is **reuse > reinvent** (the
> profiler, the benchmark harness, MUMPS/UMFPACK reuse patterns, and OpenMP are all already
> in-tree or one CMake line away). Spawns per-item sub-ADRs as each clears its Phase-0 gate.

## What

A four-phase program of **11 ranked interventions**, each scoped to the fork's actual lane
(earthquake / structural: fiber frames, moderate 3D `LadrunoConcrete3D`, explicit dynamics)
and **away** from Abaqus's lane (large 3D continuum + general contact). Two items
(`profiler-publish`, `UMFPACK STRATEGY_AUTO`) are independently source-verified and need no
measurement; the remaining nine are gated on the Phase-0 profiler/benchmark output.

**In scope:** an attributed intra-step profiler breakdown (reusing the shipped profiler); a
tiered benchmark suite with one cross-tool anchor and a parallel-scaling pair; the
correctness-adjacent UMFPACK strategy fix; in-lane kernel/assembly work-removal (brick
geometry caching, J2 stress-core/lazy-tangent split, plane-stress condensation collapse);
shared-memory OpenMP threading of the element loop (explicit diagonal path first); solver
factorization reuse (numeric, then symbolic); BLAS/MKL link verification; negative-pivot
export for the robust driver. **Out of scope (anti-goals, see Risks):** iterative
Krylov+AMG / out-of-core / GPU for large-3D continuum; hand-rolled solvers/preconditioners;
SIMD-batching the return map ahead of algorithmic work-removal; the ParMETIS/HPC stack
ahead of a measured production-scale setup bottleneck.

## Why

The fork is correctness-first and well-verified, but **un-profiled**. Optimizing now is
flying blind: the `GROUNDING` prior (fiber frames are *not* solver-bound) is plausible but
unmeasured, and the wrong guess wastes the scarcest resource — single-author weeks. The
program therefore front-loads two cheap instruments (a published profiler breakdown + an
extended benchmark suite with a cross-tool anchor) whose output **gates** which of the
downstream items actually pays. It also fixes one near-free, correctness-adjacent solver
default that today silently mis-orders the matrix the flagship concrete material explicitly
requires. Everything else is deferred behind data.

## Where

Grounded touch-points per item (file:line independently verified for ranks 1–2; ranks 3–11
are panel-sourced, to be confirmed as each clears Phase 0):

- **Profiler (rank 1):** `SRC/utility/profiler/` (`PerfClock`, `Profiler`, `ProfilerMacros`,
  `ProfilerHDF5Writer`) — *already ships*; `profiler` command registered in
  `SRC/interpreter/OpenSeesCommands.cpp`; live per-classTag scopes in
  `IncrementalIntegrator::formTangent/formElementResidual`, `TransientIntegrator`. Verify
  `theSOE->solve()` is bracketed distinctly from `formTangent`; confirm the explicit
  `CentralDifferenceLadruno::formUnbalance` path is scoped. **No new profiler.**
- **UMFPACK strategy (rank 2):** `SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.cpp`
  **lines 206 and 234** hardcode `UMFPACK_STRATEGY_SYMMETRIC` (PIVOT_TOLERANCE 1.0); parser
  `OPS_UmfpackGenLinSolver` in `OpenSeesCommands.cpp`.
- **Benchmark suite (rank 3):** `Ladruno_files/testbed/perf/runner.py` (extend, keep policy).
- **Brick geometry (rank 4):** `SRC/element/ladrunoBrick/LadrunoBrick.cpp`
  (`computeBasis`/`shp3d`/`computeB`; function-local `static` scratch → per-element members).
- **J2 kernel (ranks 5–6):** `SRC/material/nD/LadrunoJ2Kernel.h` (`returnMap`/`returnMapDamaged`
  recompute), `LadrunoJ2.cpp` (`setTrialStrain` plane-stress/plate-fiber 25-iter
  condensation, `condenseTangent`).
- **OpenMP (rank 7):** `CMakeLists.txt` (`find_package(OpenMP)`), `SRC/analysis/integrator/`
  element residual loop, `SRC/graph` (coloring, phase 2).
- **Solver reuse (ranks 8, 10):** `UmfpackGenLinSolver` (Numeric member), `MumpsSolver.cpp`
  (job=3/job=2 reference pattern already in-tree), `SuperLU.cpp`, the SOE `factored` flag.
- **Negative pivot (rank 9):** PARDISO `iparm[21/22]`, MUMPS `INFOG(12)` → R-INSTAB classifier
  in [[31_ladruno_robust_solve_driver_adr]].
- **HPC (rank 11):** `DomainPartitioner` (serial METIS on rank 0), [[30_ladruno_parallel_numberer_adr]].

## How

### The ranked program

| # | Item | Effort | Risk | Gated on |
|---|------|--------|------|----------|
| 1 | Publish attributed phase breakdown from the **existing** profiler (fiber-frame + 3D-RC) | low | low | — |
| 2 | Default UMFPACK to `STRATEGY_AUTO` (correctness-adjacent, regression-free) | low | low | — |
| 3 | Extend benchmark suite: tiered + 1 cross-tool anchor + small scaling pair | low | low | — (better after 1) |
| 4 | De-static brick scratch + cache reference-config shape gradients (thread-safety enabler) | low | low | 1 |
| 5 | J2 kernel: stress-core + lazy `consistentTangent`; drop redundant M/Mp/normM recompute | med | low | 1 |
| 6 | Collapse LadrunoJ2 plane-stress/plate-fiber condensation (skip full tangent in inner sweeps) | med | med | 5, 1 |
| 7 | ONE OpenMP element-loop effort — v1 explicit diagonal private-buffer reduction | med | med | 1 (>40% ele), 4 |
| 8 | Verify BLAS link is threaded MKL; port MUMPS factored-flag (solve-only) reuse to UMFPACK/SuperLU | med | med | 1 (solve material), 2 |
| 9 | Export negative-pivot/inertia count (MUMPS INFOG(12) / PARDISO) to the R-INSTAB classifier | low | low | — |
| 10 | Reuse symbolic factorization across `setSize` when sparsity is bit-identical | med | med | 8, 1 |
| 11 | (Conditional, demand-gated) ParMETIS partitioning + LadrunoParallelNumberer | high | med | 3 scaling benches |

### Phasing

- **Phase 0 — MEASURE (ranks 1–3).** Publish the per-phase split (formTangent/assembly,
  state-determination per classTag, geometry, addA/addB scatter, factor-vs-solve, update,
  recorders) on a fiber-frame and a 3D-RC model; land the `STRATEGY_AUTO` fix; extend the
  benchmark suite with the cross-tool anchor and a 2–8-rank scaling pair. **Phase 0 output
  gates everything after it.**
- **Phase 1 — IN-LANE, LOW-RISK, REUSE-HEAVY (ranks 4–6).** Brick de-static + shape-gradient
  cache; consolidated J2 kernel refactor (one shared converged-state POD struct); plane-stress
  condensation collapse. Attacks the fork's *own* self-imposed costs; gated on Phase 0 showing
  state-determination/assembly dominance.
- **Phase 2 — PARALLELISM + SOLVER REUSE (ranks 7–9).** Single OpenMP element-loop effort
  (explicit diagonal v1 only); BLAS-link check + UMFPACK numeric reuse; negative-pivot export.
- **Phase 3 — SECOND-ORDER + CONDITIONAL HPC (ranks 10–11).** Symbolic-reuse-across-`setSize`;
  and *only if* scaling benches prove a production-scale setup Amdahl term, the co-designed
  ParMETIS partitioner + parallel numberer. Implicit-path OpenMP colored-scatter is a deferred
  phase-2.5 follow-on, not v1.

### Benchmark plan (measurement-first)

Two instruments, both reusing what ships:
1. **Intra-step profiler** — in-tree (`SRC/utility/profiler`). Remaining work: verify
   `solve()` is scoped distinctly from `formTangent`, scope the explicit `formUnbalance`
   path, and **publish** a per-phase breakdown.
2. **Benchmark suite** — extend `runner.py` keeping its locked policy (threads pinned to 1,
   warmup discarded, median-of-7, +10% WARN / +25% FAIL, explicit re-baseline; populate the
   empty `baselines/`). Matched problems: (a) fiber-frame forceBeamColumn chain [exists];
   (b) moderate-3D `LadrunoConcrete3D` RC block (unsymmetric UMFPACK path — ranks 2/8 target);
   (c) explicit `CentralDifferenceLadruno` + SMS/HRZ step (ranks 4/7 target); (d) `LadrunoJ2`
   plane-stress/plate-fiber material-point loop (ranks 5/6 target). Metrics: median
   s/iteration end-to-end **plus** the profiler phase split; factor-count and factor-vs-solve
   for solver items; speedup vs core count for threading. Baselines: (i) self-regression per
   bench; (ii) **one cross-tool anchor** — a shared model in stock OpenSees vs fork vs (if
   licensed) Abaqus, recorded once as a reference ratio (not a gate) to kill the "no measured
   comparison" blind spot; (iii) a minimal weak+strong scaling pair (2–8 ranks, one node) to
   give the HPC items a gate. **Correctness gates ride alongside every kernel change:** the
   `LadrunoJ2` 1e-12 numpy oracle (`tests/ladrunoj2_reference.py`) and the `LadrunoConcrete3D`
   `run_tangent_gate` fixture must stay byte-identical after refactors.

## Decisions

| # | Decision | Rationale | Consequence / extension point |
|---|----------|-----------|-------------------------------|
| P1 | **Measurement-first**: ranks 1–3 land before any optimization | No attributed breakdown and no cross-tool number exist; the dominance question (state-determination vs assembly vs solve) decides which items pay | Phase 0 output is the hard gate for ranks 4–11; items can be dropped if the data says so |
| P2 | **Reuse the shipped profiler**, do not build one | `SRC/utility/profiler` already provides `PerfClock`/`Profiler`/`ProfilerMacros` + registered `profiler -deep`; the "build an ADR-06 profiler" prerequisite is stale | Only gap-filling scopes + publishing data |
| P3 | `UMFPACK_STRATEGY_AUTO` as the default (expose `-strategy`/`-pivotTol`) | Verified hardcoded `SYMMETRIC` (lines 206/234) mis-orders the unconditionally non-symmetric `LadrunoConcrete3D` tangent; AUTO reproduces today's path on SPD assemblies (zero regression) | Lands in Phase 0; record delta on the 3D-RC bench |
| P4 | Attack the fork's **own** self-imposed costs before chasing Abaqus | Brick geometry recompute, J2 redundant recompute, 25-iter plane-stress condensation are in-lane and verified; large-3D continuum is Abaqus's lane | Ranks 4–6; gated on Phase 0 dominance |
| P5 | **One** OpenMP element-loop work item, explicit diagonal path first | Three lenses proposed the same work; the explicit diagonal SOE needs no graph coloring (exact per-thread reduction) | Implicit colored-scatter deferred to phase-2.5; auto-disable below an element-count threshold |
| P6 | Numeric-factor reuse drives off the **SOE `factored` flag**, not a new algorithm-level flag | The in-tree MUMPS job=3 pattern is the reference; an invented flag risks solving against a stale LU | Scope to ModifiedNewton/InitialStiffness first; guard that A was untouched |
| P7 | **Defer the ParMETIS/HPC stack** until scaling benches prove a production-scale setup bottleneck | ADR 30's own gate: nothing needs to change for today's 66k–1M-DOF runs; ParMETIS is a heavy cross-platform (Windows) bet that does not touch per-step cost | Rank 11 stays demand-gated; co-design with [[30_ladruno_parallel_numberer_adr]] + apeGmsh per-rank shared-node emission |
| P8 | Each cleared item spawns its **own sub-ADR** with a validation plan | Matches the fork's per-feature ADR + tiered-test discipline | This ADR is the umbrella; sub-ADRs carry the implementation detail |

## Risks / open questions

> [!warning] **Anti-goals — do NOT chase these.**
> - Do **not** rebuild a profiler (it ships) — verify coverage and publish data.
> - Do **not** optimize the large-3D-continuum (+contact) regime (iterative Krylov+AMG/ILU,
>   out-of-core direct, GPU dense-plastic offload). It is Abaqus's lane and fails on the
>   fork's ill-conditioned softening tangents.
> - Do **not** hand-roll a Krylov/preconditioner or leave the unpreconditioned
>   `ConjugateGradientSolver` as a pretend option — if ever justified, wrap Eigen or PETSc KSP/PC.
> - Do **not** SIMD-batch the J2 return map before the algorithmic work-removal (ranks 5–6)
>   and a roofline showing a compute-bound dense-plastic regime. Try `/O2 /fp:fast`
>   autovectorization on the POD `double[6]` loops first.
> - Do **not** enable OpenMP element threading by default or on the implicit colored-scatter
>   path first — gate on a >40% element-fraction profiler result (fiber frames can regress).
> - Do **not** add per-instance dirty-flag tangent memoization to J2 (fragile under the
>   commit/setTrial cycle; fights the condensation path) — use the kernel split + `tang_flag`.
> - Do **not** build ParMETIS before scaling benches prove a production-scale bottleneck.

> [!question] **Phase-0 dominance question (gates the whole program).** On real production
> models, is step time state-determination-bound, assembly/geometry-bound, or solve-bound?
> The profiler output decides which of ranks 4–10 actually pay.

> [!question] **Which algorithms run in production** — full Newton vs modified-Newton /
> Krylov / IMPL-EX / dynamic-relaxation? Ranks 5 (lazy tangent) and 8 (numeric reuse) pay
> ONLY when the tangent is reused across correctors; full-Newton-dominated production
> deprioritizes both together.

> [!question] **Is the default BLAS/LAPACK threaded MKL or reference Netlib?** A
> `dumpbin /dependents` + `dgbsv_` symbol trace answers it; if reference, relinking is a
> near-free win and reorders rank 8.

> [!question] **ADR-30 Q3:** does MUMPS's internal `ICNTL7` ordering make equation-numbering
> order irrelevant on the implicit path? If yes, the parallel numberer benefit is confined to
> the diagonal/explicit path — verify before any implicit-path claim.

> [!question] **OpenMP safety:** are all Ladruno materials' `setTrialStrain`/commit state
> strictly per-Gauss-point (specifically `LadrunoConcrete3D` history variables), or is there
> shared mutable state beyond the identified brick function-local static scratch?

> [!question] **Cross-tool ratio:** what is the actual fork-vs-Abaqus (and fork-vs-stock)
> ratio per model class? Within ~1.5–2x on fiber frames (optimization is polish) or 5–10x
> somewhere it should not compete (effort misdirected)? The cross-tool anchor is the only way
> to know.

- **Topology-change signal (rank 10):** the symbolic-reuse fingerprint must use an exact
  bit-identical `Ap/Ai` compare (not a lossy hash) and a reliable constraint-handler
  topology-changed signal covering contact status, staged activation, element birth/death,
  and `ASDEmbeddedNode`, or it silently reuses a wrong ordering.
- **Backwards compat:** every kernel refactor keeps the default `returnMap` signature
  byte-identical (LogStrain finite-strain reuse + the 1e-12 oracle); P3 reproduces today's
  path on symmetric assemblies.
- **Ledger/banner debt:** each merging sub-ADR carries its `LEDGER_implementations.md` row
  and (if a new classTag) the `classTags.h` + broker + banner work per the CLAUDE.md workflow.

## Implementation log

- **2026-06-21 — origin: 11-agent performance panel (5 lenses propose → cross-critique →
  synthesis judge).** Lenses: linear-solver/sparse-LA, parallelism/HPC, explicit-dynamics
  vectorization, element/material kernel efficiency, scope/dependency/measurement pragmatist.
  The panel was grounded with the post-adversarial-review facts (inherited MUMPS/PETSc/
  OpenSeesSP-MP stack; Ladruno self-costs: non-symmetric concrete tangent, plane-stress
  condensation; existing levers IMPL-EX/SMS/HRZ; and the hard caveat that **no measured
  benchmark exists**). The synthesis corrected one stale panel premise — it proposed "build a
  profiler first," but the profiler already ships — and merged three duplicate OpenMP
  proposals into one gated item.
- **Independently verified before authoring (do not re-litigate):** (1) `SRC/utility/profiler/`
  exists with `PerfClock/Profiler/ProfilerMacros/ProfilerHDF5Writer` and the `profiler` command
  is registered in `OpenSeesCommands.cpp` (rank 1 / P2 stand). (2)
  `UmfpackGenLinSolver.cpp:206` and `:234` hardcode `UMFPACK_STRATEGY_SYMMETRIC` (rank 2 / P3
  stand). All other file:line citations are panel-sourced and to be confirmed as each item
  clears its Phase-0 gate.
- **Status:** proposed. Next action = Phase 0 (publish the profiler breakdown + extend the
  benchmark suite + land the `STRATEGY_AUTO` fix). Each downstream item promotes to its own
  sub-ADR once Phase 0 data justifies it.

*(filled in as items execute; per-item detail moves to its sub-ADR and to
`Ladruno_internal/` on completion.)*
