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
> program** for the fork, produced by an 11-agent panel discussion and refined by a v1.1
> adversarial review (see Implementation log). The governing principle is **measure before you
> optimize**: there is today **no attributed phase breakdown and no fork-vs-Abaqus /
> fork-vs-stock-OpenSees number anywhere**, so every speed claim below is a *hypothesis gated on
> Phase 0*. The second principle is **reuse > reinvent** (the profiler, the benchmark harness,
> MUMPS/UMFPACK reuse patterns, and OpenMP are all already in-tree or one CMake line away).
> Per-rank technical grounding is curated in the `opensees-performance` skill (`reading_map.md`).
> Spawns per-item sub-ADRs as each clears its Phase-0 gate.

## What

A four-phase program of **11 ranked interventions** (+ one bench-hygiene item), each scoped to the
fork's actual lane (earthquake / structural: fiber frames, moderate 3D `LadrunoConcrete3D`, explicit
dynamics) and **away** from Abaqus's lane (large 3D continuum + general contact). Two items
(`profiler-publish`, `UMFPACK STRATEGY_AUTO`) are independently source-verified and need no
measurement; the remaining nine are gated on the Phase-0 profiler/benchmark output.

**In scope:** an attributed intra-step profiler breakdown (reusing the shipped profiler); a
tiered benchmark suite with one cross-tool anchor and a parallel-scaling pair; DOF-numberer
bench hygiene; the correctness-adjacent UMFPACK strategy fix; in-lane kernel/assembly work-removal
(brick geometry caching, J2 stress-core/lazy-tangent split, plane-stress condensation collapse);
shared-memory OpenMP threading of the element loop (explicit diagonal path first); solver
factorization reuse (numeric, then symbolic); BLAS/MKL link verification; negative-pivot export.
**Out of scope (anti-goals, see Risks):** iterative Krylov+AMG / out-of-core / GPU for large-3D
continuum; hand-rolled solvers/preconditioners; SIMD-batching the return map ahead of algorithmic
work-removal; the ParMETIS/HPC stack ahead of a measured production-scale setup bottleneck.

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

Grounded touch-points per item (file:line verified for ranks 1–2 and the v1.1-confirmed items;
remaining panel-sourced lines to be confirmed as each clears Phase 0):

- **Profiler (rank 1):** `SRC/utility/profiler/` (`PerfClock`, `Profiler`, `ProfilerMacros`,
  `ProfilerHDF5Writer`) — *already ships*; `profiler` command registered in `OpenSeesCommands.cpp`.
  **`NewtonRaphson` already scopes formUnbalance / formTangent / linearSolve(`theSOE->solve()`) /
  update distinctly** (`NewtonRaphson.cpp:144/165/192/200`) — that path is DONE. The Phase-0 gaps:
  `ModifiedNewton.cpp` (the factorization-reuse path for ranks 8/10) has **zero** `OPS_PROFILE`
  scopes; the explicit `CentralDifferenceLadruno` path is unscoped; and **no factor-vs-triangular-
  solve scope exists** anywhere under `SRC/system_of_eqn/` (today one coarse `linearSolve` scope
  wraps factor+solve together) — yet that split is the gating metric for ranks 8/10, so Phase 0 must
  build it. **No new profiler.**
- **Serial DOF numbering (item 3b):** `numberer Plain|RCM|AMD` (`SRC/graph/numberer` +
  `SRC/analysis/numberer`; dispatch `OpenSeesCommands.cpp:1559-1572`). OpenSees' **default is already
  RCM** (`OpenSeesCommands.cpp:290-291` + 4 more sites). For `BandGeneral/BandSPD/ProfileSPD`, factor
  + storage cost is derived straight from DOF-index differences
  (`BandGenLinSOE.cpp:146-168`: `numSubD/numSuperD` from `diff=vertexNum-otherNum`,
  `newSize=size*(2*numSubD+numSuperD+1)`), so RCM/AMD cut it; it is **~neutral** for UMFPACK/MUMPS,
  which reorder internally (`UmfpackGenLinSolver.cpp:217/245`; `MumpsSolver.cpp:75` `ICNTL(7)`). The
  fork bench leaves it on the table: `Ladruno_files/testbed/perf/runner.py:69-70` does
  `system("BandGeneral"); numberer("Plain")`. **Action = bench hygiene, not a code change.**
- **UMFPACK strategy (rank 2):** `UmfpackGenLinSolver.cpp` **205/206 (dl)** and **233/234 (di)** set
  `Control[UMFPACK_PIVOT_TOLERANCE]=1.0` then `Control[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_SYMMETRIC`
  (the AUTO fix must touch **both** branches). Parser `OPS_UmfpackGenLinSolver` (`:43-63`) currently
  accepts only `useLongIndices` — extend it for `-strategy` and `-pivotTol`.
- **Benchmark suite (rank 3):** `Ladruno_files/testbed/perf/runner.py` (extend, keep policy).
- **Brick geometry (rank 4):** `LadrunoBrick.cpp:86` file-scope `static Matrix B(6,3)` + function-local
  static scratch (`:526-536`, `:600`, `:710-713`) — a genuine OpenMP **data race** (the rank-7
  enabler) AND a per-call geometry recompute (`Shape[...]` rebuilt every `formTangent` via
  `computeBasis()/shp3d()` instead of cached at the reference config). **Hard predecessor of rank 7.**
- **J2 kernel (ranks 5–6):** `LadrunoJ2Kernel.h` — `returnMap` fills `Dtan[6][6]` **unconditionally**
  (no `tang_flag`), `returnMapDamaged` adds a second 6×6 assembly (rank-5 seam = a `tang_flag`
  parameter, **not** dirty-flag memoization). `LadrunoJ2.cpp` plane-stress/plate-fiber: outer `eps_22`
  Newton `maxIt=25` (`:304`) × inner 3D scalar return `maxIter=50` (`LadrunoJ2Kernel.h:224`) →
  worst-case ~1250 scalar steps + a full 6×6 rebuilt every outer sweep; `condenseTangent`.
- **OpenMP (rank 7):** `CMakeLists.txt` (`find_package(OpenMP)`), the element residual loop, `SRC/graph`
  (coloring, phase 2). Gated on rank 4 (race removed) + the OpenMP-safety open question. **Scope
  correction (2026-07-05):** rank 4 (de-static `LadrunoBrick`) is **necessary but far from
  sufficient**. The *assembly layer itself* shares class-wide static scratch: `FE_Element`
  (`FE_Element.h:132-133`, `static Matrix **theMatrices / static Vector **theVectors`) and
  `DOF_Group` (`DOF_Group.h:150-151`) hold pools keyed by DOF size, **shared across all
  instances of the same size** — so `FE_Element::getResidual/getTangent` races under a threaded
  element loop even on a pure-`LadrunoBrick` mesh with rank 4 done. On top of that, file-scope
  static work matrices are the norm across the vanilla element library (`FourNodeQuad.h:121`
  `static Matrix K`, etc.). A private-buffer reduction fixes only the global `addB` scatter, not
  these per-class internal workspaces. Rank 7 must therefore either (a) restrict to a verified
  race-free element set, or (b) carry an assembly-layer + per-element-type de-static audit —
  that audit, not rank 4, is the true predecessor.
- **Solver reuse (ranks 8, 10):** UMFPACK **already retains** the `Symbolic` handle across solves
  (member, ctor `Symbolic(0)` at `UmfpackGenLinSolver.cpp:69`, freed only in setSize/dtor) — symbolic
  analysis is **not** the per-solve waste. The waste: the `Numeric` handle is allocated AND freed on
  **every** solve (`umfpack_di_numeric` `:160` then `umfpack_di_free_numeric` `:176`; dl path
  `:129/:145`). Rank 8 = retain `Numeric`, skip `umfpack_*_numeric` when the SOE `factored` flag says
  `A` is untouched (a **numeric re-factorization avoided**); the `umfpack_*_refactor` primitive (reuse
  pattern + pivot ordering, new values) is the cheaper middle path. Rank 10 = reuse the `Symbolic`
  handle across `setSize`. `MumpsSolver.cpp` (job=3/job=2) is the in-tree reference; the SOE `factored`
  flag is the single source of truth.
- **Negative pivot (rank 9):** PARDISO `iparm[21/22]`, MUMPS `INFOG(12)` → R-INSTAB classifier
  in [[31_ladruno_robust_solve_driver_adr]].
- **HPC (rank 11):** `DomainPartitioner` (serial METIS on rank 0), [[30_ladruno_parallel_numberer_adr]].

## How

### The ranked program

| # | Item | Effort | Risk | Gated on |
|---|------|--------|------|----------|
| 1 | Publish attributed phase breakdown from the **existing** profiler (fiber-frame + 3D-RC); fill the ModifiedNewton/explicit + factor-vs-solve scope gaps | low | low | — |
| 2 | Default UMFPACK to `STRATEGY_AUTO` + expose `-strategy`/`-pivotTol` (correctness-adjacent, regression-free) | low | low | — |
| 3 | Extend benchmark suite: tiered + 1 cross-tool anchor + small scaling pair | low | low | — (better after 1) |
| 3b | DOF-numberer hygiene: take the perf bench off `numberer Plain`; record Plain vs RCM vs AMD as a free data point on band/profile paths | low | low | 1 |
| 4 | De-static brick scratch + cache reference-config shape gradients (thread-safety enabler; **hard predecessor of 7**) | low | low | 1 |
| 5 | J2 kernel: stress-core + lazy `consistentTangent` via a `tang_flag`; drop redundant M/Mp/normM recompute | med | low | 1 |
| 6 | Replace the nested `eps_22` condensation (25×50 Newton, full 6×6 every sweep) with a plane-stress-**projected** return map (closed-form quartic) | med | med | 5, 1 |
| 7 | ONE OpenMP element-loop effort — v1 explicit diagonal private-buffer reduction | med | med | 1 (>40% ele), 4 |
| 8 | Retain the UMFPACK `Numeric` handle (skip refactor when `factored`); port MUMPS solve-only reuse; verify BLAS link is threaded MKL | med | med | 1 (solve material), 2 |
| 9 | Export negative-pivot/inertia count (MUMPS INFOG(12) / PARDISO) to the R-INSTAB classifier | low | low | — |
| 10 | Reuse the `Symbolic` handle across `setSize` when sparsity is bit-identical | med | med | 8, 1 |
| 11 | (Conditional, demand-gated) ParMETIS partitioning + LadrunoParallelNumberer (distinct triggers; see P7) | high | med | 3 scaling benches |

### Phasing

- **Phase 0 — MEASURE (ranks 1–3 + 3b).** Publish the per-phase split (formTangent/assembly,
  state-determination per classTag, geometry, addA/addB scatter, **factor-vs-solve**, update,
  recorders) on a fiber-frame and a 3D-RC model; land the `STRATEGY_AUTO` fix; extend the
  benchmark suite (cross-tool anchor + 2–8-rank scaling pair); fix the bench numberer downgrade.
  **Phase 0 output gates everything after it.**
- **Phase 1 — IN-LANE, LOW-RISK, REUSE-HEAVY (ranks 4–6).** Brick de-static + shape-gradient
  cache; consolidated J2 kernel refactor (one shared converged-state POD struct + `tang_flag`);
  plane-stress projected return map. Gated on Phase 0 showing state-determination/assembly dominance.
- **Phase 2 — PARALLELISM + SOLVER REUSE (ranks 7–9).** Single OpenMP element-loop effort
  (explicit diagonal v1 only); UMFPACK numeric reuse + BLAS-link check; negative-pivot export.
- **Phase 3 — SECOND-ORDER + CONDITIONAL HPC (ranks 10–11).** Symbolic-reuse-across-`setSize`;
  and *only if* scaling benches prove a production-scale setup Amdahl term, the co-designed
  ParMETIS partitioner + parallel numberer. Implicit-path OpenMP colored-scatter is a deferred
  phase-2.5 follow-on, not v1.

### Benchmark plan (measurement-first)

Two instruments, both reusing what ships:
1. **Intra-step profiler** — in-tree (`SRC/utility/profiler`). Remaining work: bring `ModifiedNewton`
   and the explicit `CentralDifferenceLadruno` path up to `NewtonRaphson`'s scope coverage, **build the
   factor-vs-triangular-solve scope** (it does not exist yet, and it is the rank-8/10 gating metric),
   and **publish** the per-phase breakdown.
2. **Benchmark suite** — extend `runner.py` keeping its locked policy (threads pinned to 1,
   warmup discarded, median-of-7, +10% WARN / +25% FAIL, explicit re-baseline; populate the
   empty `baselines/`). Matched problems: (a) fiber-frame forceBeamColumn chain [exists];
   (b) moderate-3D `LadrunoConcrete3D` RC block (unsymmetric UMFPACK path — ranks 2/8 target);
   (c) explicit `CentralDifferenceLadruno` + SMS/HRZ step (ranks 4/7 target); (d) `LadrunoJ2`
   plane-stress/plate-fiber material-point loop (ranks 5/6 target). Metrics: median
   s/iteration end-to-end **plus** the profiler phase split; factor-count and factor-vs-solve
   for solver items; speedup vs core count for threading. Baselines: (i) self-regression per
   bench; (ii) **one cross-tool anchor** — recorded once as a reference ratio, **never a gate**.
   **Status:** the fork-vs-stock-OpenSees arm runs today (fork binary + openseespy on PATH) and
   answers the dominance and self-regression questions; the **Abaqus arm is license-blocked** —
   Abaqus 2025 is installed on the dev machine but **unlicensed (no DSLS server configured)**, so the
   "is the fork ~1.5× or ~10× off an industrial code" question is **deferred** until a DSLS entitlement
   is configured. (iii) a minimal weak+strong scaling pair (2–8 ranks, one node) to give the HPC items
   a gate. **Correctness gates ride alongside every kernel change:** the `LadrunoJ2` 1e-12 numpy oracle
   (`tests/ladrunoj2_reference.py`) and the `LadrunoConcrete3D` `run_tangent_gate` fixture must stay
   byte-identical after refactors.

## Decisions

| # | Decision | Rationale | Consequence / extension point |
|---|----------|-----------|-------------------------------|
| P1 | **Measurement-first**: ranks 1–3 land before any optimization | No attributed breakdown and no cross-tool number exist; the dominance question (state-determination vs assembly vs solve) decides which items pay | Phase 0 output is the hard gate for ranks 4–11; items can be dropped if the data says so |
| P2 | **Reuse the shipped profiler**, do not build one | `SRC/utility/profiler` already provides `PerfClock`/`Profiler`/`ProfilerMacros` + registered `profiler -deep`; NewtonRaphson is already scoped | Phase 0 = fill the ModifiedNewton/explicit + factor-vs-solve scope gaps and publish |
| P3 | `UMFPACK_STRATEGY_AUTO` default (expose `-strategy`/`-pivotTol`) | Verified hardcoded `SYMMETRIC` **+ `PIVOT_TOLERANCE=1.0`** (lines 205/206 & 233/234) both mis-orders AND over-restricts pivoting on the unconditionally non-symmetric `LadrunoConcrete3D` tangent (stability risk on indefinite blocks); AUTO reproduces today's path on SPD assemblies (zero regression) | Lands in Phase 0; record delta on the 3D-RC bench |
| P4 | Attack the fork's **own** self-imposed costs before chasing Abaqus | Brick geometry recompute, J2 redundant recompute, 25×50 plane-stress condensation are in-lane and verified; large-3D continuum is Abaqus's lane | Ranks 4–6; gated on Phase 0 dominance |
| P5 | **One** OpenMP element-loop work item, explicit diagonal path first | Three lenses proposed the same work; the explicit diagonal SOE needs no graph coloring (exact per-thread reduction) | Implicit colored-scatter deferred to phase-2.5; auto-disable below an element-count threshold |
| P6 | Numeric-factor reuse drives off the **SOE `factored` flag**, not a new algorithm-level flag | The in-tree MUMPS job=3 pattern is the reference; an invented flag risks solving against a stale LU | Scope to ModifiedNewton/InitialStiffness first; guard that A was untouched |
| P7 | **Defer the ParMETIS/HPC stack** until scaling benches prove a production-scale bottleneck | ADR 30's own gate: nothing changes for today's 66k–1M-DOF runs; ParMETIS is a heavy cross-platform bet that does not touch per-step cost | Rank 11's two halves have **distinct triggers** (ParMETIS partitioning gated on `DomainPartitioner` serial-METIS-on-rank-0; the parallel numberer on ADR-30's measured numbering term) and can land independently — inherit the gate **by reference** to [[30_ladruno_parallel_numberer_adr]], do not copy its numbers here |
| P8 | Each cleared item spawns its **own sub-ADR** with a validation plan | Matches the fork's per-feature ADR + tiered-test discipline | This ADR is the umbrella; sub-ADRs carry the implementation detail. **First spawned: [[68_ladruno_state_determination_perf_adr]]** — the element/state-determination OPTIMIZE design (ranks 4/5/6 + the §3-40a cost centers), Phase-0-gated, carrying the per-kernel byte-identical-vs-equivalence gate taxonomy |
| P9 | **RCM/AMD numbering** cuts `BandGeneral/BandSPD/ProfileSPD` bandwidth/profile cost; **~neutral** for UMFPACK/MUMPS (internal reorder) | The numberer choice and the rank-2 UMFPACK STRATEGY fix address **different solver families** and must not be conflated; orthogonal to the *parallel* numberer in [[30_ladruno_parallel_numberer_adr]] (distributed-setup Amdahl term) — do not merge | Item 3b is bench hygiene (default is already RCM; the bench downgrades to Plain), not a code lever; whether AMD beats the RCM default on the fiber-frame bench is an open measurement |

## Risks / open questions

> [!warning] **Anti-goals — do NOT chase these.**
> - Do **not** rebuild a profiler (it ships) — fill the scope gaps and publish data.
> - Do **not** optimize the large-3D-continuum (+contact) regime (iterative Krylov+AMG/ILU,
>   out-of-core direct, **GPU *solver* offload**). It is Abaqus's lane and fails on the
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
> - Do **not** treat the cross-tool anchor as a pass/fail gate, and do **not** let the
>   license-blocked Abaqus arm stall any Phase-1+ item — no optimization item may be gated on a
>   fork-vs-Abaqus number.
> - **GPU — split the flavor (it is not a monolith).** GPU *linear-solver* offload
>   (cuSOLVER/cuSPARSE/AmgX) stays an anti-goal (the large-3D softening-tangent regime above). GPU
>   *constitutive-kernel batching* (the per-Gauss-point return map — the genuinely SIMT-amenable part) is a
>   **conditional far-future extension of the vectorization track (ranks 4–7)**, gated on Phase-0
>   state-determination dominance + a roofline (compute-bound) + ranks 5–6 (branch-heavy
>   condensation/IMPL-EX removed) done first; sparse, uncorrelated seismic yielding causes warp divergence
>   that can erase the gain. GPU *contact* is **out** (the fork has no general contact, and that regime is
>   the anti-goal). Path, if ever: matrix-free batching (libCEED/MFEM/Kokkos/deal.II-MF), **not** a GPU
>   solver.

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
> the diagonal/explicit path — verify before any implicit-path claim. (Same mechanism makes the
> serial numberer, item 3b, ~neutral for the sparse solvers.)

> [!question] **OpenMP safety:** are all Ladruno materials' `setTrialStrain`/commit state
> strictly per-Gauss-point (specifically `LadrunoConcrete3D` history variables), or is there
> shared mutable state beyond the identified brick function-local static scratch (rank 4)?
> **This question is broader than first framed** (2026-07-05): the blocking shared state is not
> only in Ladruno materials but in the **assembly layer** (`FE_Element`/`DOF_Group` class-wide
> static matrix/vector pools, `FE_Element.h:132-133`, `DOF_Group.h:150-151`) and in the
> **vanilla element library**'s file-scope static work matrices — both race under a threaded
> element loop independent of any Ladruno code. Answering "is the loop thread-safe" requires an
> audit of the assembly layer + every element type present, not just the Ladruno materials.

> [!question] **Cross-tool ratio:** what is the actual fork-vs-Abaqus (and fork-vs-stock)
> ratio per model class? The fork-vs-stock arm runs now; the Abaqus arm is license-deferred.

- **Topology-change signal (rank 10):** the symbolic-reuse fingerprint must use an exact
  bit-identical `Ap/Ai` compare (not a lossy hash) and a reliable constraint-handler
  topology-changed signal covering contact status, staged activation, element birth/death,
  and `ASDEmbeddedNode`, or it silently reuses a wrong ordering.
- **Backwards compat:** every kernel refactor keeps the default `returnMap` signature
  byte-identical (LogStrain finite-strain reuse + the 1e-12 oracle); P3 reproduces today's
  path on symmetric assemblies.
- **Ledger/banner debt:** each merging sub-ADR carries its `LEDGER_implementations.md` row
  and (if a new classTag) the `classTags.h` + broker + banner work per the CLAUDE.md workflow.

## References (on-hand library)

Per-rank grounding is curated in the `opensees-performance` skill →
`references/reading_map.md`; PDFs in Seafile `C:/nmb/My Libraries/nmb/skills/opensees-performance-manuals/`.

- **Ranks 2 / 8 / 10:** Davis, *Direct Methods for Sparse Linear Systems* (SIAM 2006)
  [.djvu, cited by chapter — the symbolic/numeric split, elimination trees, ordering] + the
  SuiteSparse/UMFPACK source; the MUMPS solve-only (job=3) reuse analogue is Liu, *The Multifrontal
  Method for Sparse Matrix Solution* (SIAM Review 34(1), 1992).
- **Rank 5 (3D stress-core + consistent tangent):** Simo & Hughes, *Computational Inelasticity*
  (1998), Box 3.1 (consistency condition) + Box 3.2 (radial return) + §3.3.2 (exact linearization).
- **Rank 6 (closed-form plane-stress return):** Simo & Hughes Box 3.3 (plane-stress return mapping) +
  Box 3.4 / §3.4.5 (closed-form quartic — explicitly ~4–6 Newton iterations and vectorizable). Secondary
  (matches the code's own `LadrunoJ2.cpp:303` comment): de Souza Neto §9.2.3 / §9.4 — **[acquire]**, not
  on hand.
- **Ranks 4 / 7 (vectorization / roofline):** Williams, Roofline (Berkeley TR EECS-2008-134); Agner Fog.
- **Rank 11:** METIS/ParMETIS (Karypis) + McKenna-Scott-Fenves 2010 (the OpenSees parallel object model)
  are the **implementation** references; Smith 1996 / Toselli & Widlund 2005 are **DD-solver-theory
  background only** (rank 11 builds no Schwarz/FETI/BDDC preconditioner — that lane is an anti-goal).
- **Iterative anti-goal:** Saad, *Iterative Methods*; Barrett, *Templates*.
- **IMPL-EX:** the on-hand PDF is the modified variant (Prazeres et al. 2016); the canonical Oliver et
  al. 2008 IMPL-EX paper is **[acquire]**.

## Implementation log

- **2026-06-21 — origin: 11-agent performance panel (5 lenses propose → cross-critique →
  synthesis judge).** Grounded with the post-adversarial-review facts (inherited MUMPS/PETSc/
  OpenSeesSP-MP stack; Ladruno self-costs: non-symmetric concrete tangent, plane-stress
  condensation; existing levers IMPL-EX/SMS/HRZ; no measured benchmark). Corrected the stale "build a
  profiler" premise (it ships) and merged three duplicate OpenMP proposals into one gated item.
- **Independently verified before authoring (do not re-litigate):** `SRC/utility/profiler/` exists +
  `profiler` registered; `UmfpackGenLinSolver.cpp:206/234` hardcode `UMFPACK_STRATEGY_SYMMETRIC`.
- **2026-06-21 — v1.1 review pass (11-agent discussion + adversarial review, source-re-verified).**
  Added item **3b** DOF-numberer bench hygiene (default is already RCM; the bench downgrades to
  `numberer Plain`; band/profile only, ~neutral for sparse; orthogonal to [[30_ladruno_parallel_numberer_adr|ADR 30 (parallel numberer)]]). Sharpened **rank 8**
  (UMFPACK `Symbolic` already retained — only `Numeric` is freed per solve, so this avoids a numeric
  *re-factorization*; added `umfpack_*_refactor` + Davis vocabulary). Added **rank 2** `PIVOT_TOLERANCE
  1.0` detail + expose `-strategy`/`-pivotTol`. Re-characterized **rank 6** as a nested 25×50 Newton and
  **CORRECTED its citation**: the plane-stress closed-form is Simo & Hughes **Box 3.3 + Box 3.4/§3.4.5**
  — **Box 3.1/3.2 were mis-applied and are RE-HOMED to rank 5** (Box 3.1 is the general 3D consistency
  condition, not plane stress); de Souza Neto marked **[acquire]**. Corrected the **profiler-coverage**
  task (NewtonRaphson is done; `ModifiedNewton` + the explicit path + a not-yet-existing factor-vs-solve
  scope are the Phase-0 gaps). Stated the **cross-tool Abaqus arm is license-blocked** (installed,
  unlicensed, no DSLS) + added a non-gating anti-goal. Pinned **rank-11** refs to METIS/McKenna
  (Smith/Toselli = DD-theory background only). Added this **References** block + a bidirectional
  `opensees-performance` skill cross-link. Status stays: **proposed**.
- **2026-06-21 — v1.2:** scoped GPU explicitly (it had been conflated with solver-offload). Split into
  GPU *solver* offload (anti-goal, Abaqus's lane) vs GPU *constitutive-kernel batching* (the SIMT-amenable
  per-GP return map) — the latter a conditional far-future vectorization-track extension, gated on Phase-0
  dominance + roofline + ranks 5–6, with GPU contact out (no general contact; anti-goal regime). Raised by
  the maintainer's question about GPU for material iterations / contact.

- **2026-07-06 — Phase-0 FIRST MEASURED PASS → [[40b_phase0_dominance_report]].** Five lanes
  profiled with the shipped profiler on the existing binary (zero rebuild). Dominance is per-lane:
  fiber frame = `update` 64.5% (force-based interior — INVISIBLE to per-classTag buckets);
  3D solid J2 = **UmfPack solve 66.4%** (rank 2 PROMOTED); plate-fiber shell = element 55.9% but
  J2 condensation ceiling ≈2% (ADR-68 T2 DEMOTED); explicit CDL = element 48.9% (**rank-7 gate
  PASSED**; "explicit path unscoped" REFUTED — 99.9% closure; two NEW integrator costs found:
  per-step tangent forms 15.6% + newStep 23.9%); IMK = hinge Newton hides in opaque `update`
  (ADR-40a §3.2 confirmed by measurement). **Revised #1 scope gap: `Element::update()` needs a
  per-classTag `elem.update` bucket** — ahead of the ModifiedNewton + factor-vs-solve gaps listed
  above. Verdict table + next actions in 40b.

- **2026-07-25/26 — [[76_ladruno_tangent_reuse_adr]] (assembly/solve-lane sub-ADR, amends this
  program) OPENED AND CLOSED.** From an external issue report (`Newton -initial` re-assembles +
  re-factorizes every iteration, measured 1.61× `ModifiedNewton -initial` at 39k DOF for an
  identical answer). Shipped: the documentation (R1), the `OPS_ModifiedNewton` multi-option parser
  (R4, 11/11 smoke), and a spun-off LAPACK **singular-matrix-reported-SUCCESS** fix in
  BandGen/FullGen/BandSPD (CI-gated regression deck). The engine-side tangent-version counter (R2)
  was **WITHDRAWN after adversarial review** — the load-bearing finding, worth reading from the
  program level: `NDMaterial::getInitialTangent()` *defaults to* `getTangent()`, so "initial
  stiffness" is state-dependent for most solid element/material pairs and `-initial` on such models
  silently IS full Newton. Any future reuse predicate must default false (its Appendix A.4).

*(filled in as items execute; per-item detail moves to its sub-ADR and to
`Ladruno_internal/` on completion.)*
