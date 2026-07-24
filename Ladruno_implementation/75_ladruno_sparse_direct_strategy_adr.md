---
title: Sparse-direct solver strategy — desktop (PARDISO) + cluster (MUMPS) + threaded assembly
project: Ladruno
status: proposed — scoping (code-verified inventory; measure-gates defined; cross-framework precedent folded)
priority: medium
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - adr
  - performance
  - solver
  - sparse-direct
  - pardiso
  - mumps
  - openmp
  - threads
  - mpi
  - desktop
  - cluster
  - sub-adr
---

# ADR-75 — Sparse-direct solver strategy: PARDISO (desktop) + MUMPS (cluster) + threaded assembly

> ADR-75. The **solve/parallelism-lane** perf sub-ADR that [[40_ladruno_performance_adr]] levers
> #1 (sparse direct solvers) and #3 (parallelism) call for. Scope set by a user goal —
> **"we want both desktop and cluster"** performance — and gated by the measured dominance in
> [[40b_phase0_dominance_report]]: only the **3D-solid lane (Lane B, 66% linearSolve)** is
> solver-bound; the fork's primary **fiber-frame lanes are `update`/element-bound**, where *no*
> linear solver helps. This ADR therefore spans three lanes: two solver lanes (one per hardware
> regime) plus the shared-memory **element-assembly** lane that is the only thing that speeds up
> the frame models. Family: ADR-40 (perf program) · ADR-40b (dominance) · ADR-74 (setup lane).

## 1. Context — two hardware regimes, two axes of parallelism

The user wants strong performance on **both** a single desktop and a multi-node cluster. Those
map onto the two orthogonal axes of parallelism, which OpenSees exposes very differently:

- **Threads (shared memory)** — one process, many threads sharing one address space; scales to the
  cores of **one node**; no message passing. Launched as a normal `import openseespy` run.
- **MPI (distributed memory)** — many processes ("ranks"), private memory, explicit messages;
  scales across **many nodes**; launched with `mpiexec`. This is what OpenSees calls **SP**
  (parallel solver, one domain) and **MP** (partitioned domain) — *both are MPI, neither is threads.*

Today OpenSees has **no shared-memory build axis at all**: `find_package(OpenMP)` is absent, the
only `#pragma omp` in `SRC/` are 7 dead lines in PFEM code, and the `omp_set_num_threads` command
is compiled out. The *only* multicore parallelism available in a serial build is whatever **MKL's
threaded BLAS** does *inside* a solver — invisible to OpenSees.

### Dominance recap (why the three lanes, [[40b_phase0_dominance_report]])

| Lane | Dominant cost | Solver (PARDISO/MUMPS) helps? | Threaded assembly helps? |
|---|---|---|---|
| A fiber frame | `update` 64.5% | ❌ (solve ≈3%) | ✅ |
| B 3D solid (LadrunoBrick+J2) | **solve 66%** / element 30% | ✅ | ✅ |
| C plate-fiber shell | element 56% | ❌ | ✅ |
| D explicit CDL | element 49% | ❌ (diagonal, no factor) | ✅ |
| E IMK frame | `update` 35% + solve 31% | partial | ✅ |

A direct solver wins **one** lane. Threaded element assembly is the dominant cost in **four of
five** — and is the *only* lever for the primary frame lanes. This is the measure-first spine of
the whole ADR.

## 2. Code-verified inventory (what exists today)

Registered `system` verbs (from `OpenSeesCommands.cpp`): `UmfPack`, `SuperLU`, `SparseSYM`,
`SparseGEN`, `FullGeneral`, `BandGeneral`, `Diagonal`, `MPIDiagonal`, `Mumps`, `Petsc`, … .

| Solver | Files | State | Threaded? | MPI? |
|---|---|---|---|---|
| **UmfPack** | `linearSOE/umfGEN/` | wired; Lane-B baseline | BLAS only | no |
| **MUMPS** | `linearSOE/mumps/` | wired **MP/SP only** (`_MUMPS` gated by `if(MPI_FOUND)`) | BLAS + `ICNTL16` | **yes** |
| **MKL PARDISO** | `linearSOE/pardiso/PARDISOGenLin{SOE,Solver}` | **in-tree, UNWIRED** (no `system Pardiso`, not in build) | **native OpenMP** | no |
| **PETSc** | `linearSOE/petsc/` | wired but unbuilt on this toolchain | yes | optional |

Key facts established during scoping:
- `MumpsSolver.cpp:50-52` already `#include <libseq\mpi.h>` for a non-MPI build, and
  `OPS_MumpsSolver()` has a serial branch (`:4997`) — a **serial MUMPS path exists in source** but
  is never built (and `libseq`/`libmpiseq` is **not bundled**; the Windows build via
  `scivision/mumps` v5.5.1.5 is always **Intel-MPI, LP64** — nnz capped at 2³¹).
- `PARDISOGenLinSolver.cpp:23` uses `<mkl_pardiso.h>` (`mtype=11`, METIS ordering `iparm[1]=2`) —
  it is **MKL PARDISO**. `FeastEigenSolver.cpp:132` already calls MKL `pardiso()` and factors with
  it, **proving the symbols link and run on the exact oneAPI/Windows toolchain.**
- Factorization reuse is already correct in MUMPS (`job=3` solve-only gated on the SOE `factored`
  flag; `job=5` factor+solve otherwise) — unlike UmfPack, which re-factors every solve.
- `system Mumps` already parses `-matrixType 0|1|2` (unsym / SPD / general symmetric) and
  `-ICNTL7` (ordering) / `-ICNTL14` (workspace) / `-commSplit` (ADR-43 sub-communicators) — but the
  symmetric path is **unexercised** by the perf scripts.

## 3. The three lanes

### Lane 1 — PARDISO owns the desktop (shared-memory direct)
Native-OpenMP MKL PARDISO in the plain serial `opensees.pyd`: all cores on the factor/solve, **no
MPI**, **zero new dependency** (MKL already linked), ~90% of the code already in-tree and proven to
link.
- Finish/verify `PARDISOGenLinSOE`/`Solver`; register `system Pardiso`; lift a `_PARDISO`/link into
  the **serial** `OpenSees`/`OpenSeesPy` targets (mirror the `_MUMPS` wiring, minus the MPI gate).
- Symmetric `mtype` (`2`/`-2`) for symmetric tangents; phase-33 solve-only reuse gated on the SOE
  `factored` flag (match MUMPS); expose thread count (`MKL_NUM_THREADS`/`mkl_set_num_threads`).
- **Only the solve is threaded** — Amdahl: real win on Lane B, ~nil on frame lanes (that's Lane 3).

### Lane 2 — MUMPS owns the cluster (distributed MPI direct)
Already shipped and ADR-74-hardened at 18 M-hex/np240. This lane is **tuning, not building**:
- Exercise `-matrixType 2` (symmetric) in real decks — ~2× time+memory where the tangent is
  symmetric.
- Add a `-BLR`/`ICNTL35` knob (block low-rank compression; MUMPS 5.5.1 supports it) — memory/time
  at large 3D scale, a genuine cluster-only edge PARDISO lacks.
- Tune **hybrid** ranks×threads (`ICNTL16` + threaded MKL): fewer ranks × more threads per node.
- Housekeeping: strip the serial-path `std::cerr` debug chatter (`MumpsSolver.cpp:59,82`).
- **Descoped:** serial/`libseq` MUMPS. PARDISO owns desktop, so the fragile sequential-lib build
  (and the LP64/ILP64 question) is **not** pursued. MUMPS stays cluster-only — its strength.

### Lane 3 — Threaded element assembly (the missing axis; helps BOTH regimes)
The only lever for the frame lanes, and the **threads-per-node half of hybrid MPI+threads** on the
cluster — so it strengthens desktop *and* cluster. Net-new, highest effort/risk.
- **Two hazards:** (a) `static`/shared mutable scratch in element/material kernels → data race
  (the fork already flags "de-static the brick scratch"); (b) the `addA` **scatter collision**
  (elements sharing a DOF write the same entries).
- **Remedies (revised by §4 precedent, in preference order):** (a) **private-buffer reduction**
  (explicit/diagonal SOE — exact, no coloring; *start here*, fewest moving parts); (b) **atomic-
  scatter on a frozen sparsity** (implicit — build the CSR graph once, then a threaded parallel-for
  with per-entry `AtomicAdd`; Kratos's proven answer, simpler than coloring and memory-cheaper than
  thread-private matrices); (c) **element ordering into conflict-free groups** (LS-DYNA — removes the
  write conflict at the source *and* enables SIMD; the stronger but larger option); (d) **graph
  coloring** as a fallback only if atomics contend badly. Note this **supersedes the perf-skill's
  coloring-first sketch** — two production codes went atomic/ordering, not coloring. Hard rule:
  never reallocate during the threaded phase.
- **Gate:** only where the profiler shows **>40% element/`update` fraction**; trivially-cheap
  elements (forceBeamColumn `getTangent` ≈0.1 µs) can *regress* from fork-join overhead.
  "OpenMP-by-default / implicit colored-scatter first" is an **anti-goal** (inherited from ADR-40).
- **First move is instrumentation, not code:** ADR-40b's #1 gap is that the `elem.update` loop
  (force-based / IMK interior iteration — the frame lanes' real cost) is unscoped. Measure it first.

## 4. Cross-framework precedent (Abaqus · Kratos · LS-DYNA)

Consulted per the "reuse mature precedent" discipline (`abaqus-theory` + `kratos` skills; LS-DYNA
from vendor docs). All three independently implement the same desktop-threads / cluster-MPI split we
propose — and two of them (Kratos, LS-DYNA) point away from my initial "graph-coloring" instinct for
the assembly race toward a simpler, better-proven remedy.

| Axis | Abaqus/Standard | Kratos | LS-DYNA | → ADR-75 choice |
|---|---|---|---|---|
| **Desktop (shared-mem)** | `mp_mode=THREADS`, threaded multifrontal | OpenMP `block_for_each` + MKL PARDISO | **SMP** (`ncpu`), threaded multifrontal | **PARDISO** (native OpenMP) |
| **Cluster (distributed)** | `mp_mode=MPI` / **DMP**, distributed direct | Trilinos/Epetra + **Amesos→MUMPS** | **MPP**, distributed multifrontal; **MUMPS** (`LSOLVR=30`) | **MUMPS** (MPI) |
| **Symmetric policy** | auto-detect + symmetric-approx fallback | **explicit** (`pardiso_ldlt/llt` vs `pardiso_lu`) | symmetric default + unsym path | default-symmetric **+ explicit flag** |
| **Assembly-race remedy** | assembly parallel on both axes; solver-owned | **freeze CSR sparsity, then per-entry `AtomicAdd`** — *no coloring, no thread-private matrices* | **element ordering into conflict-free groups** (also enables SIMD) + ordered reduction | explicit private-buffer → **atomic-scatter on frozen sparsity** (coloring only if it contends) |
| **Determinism** | (order-dependent) | atomic order varies | **explicit consistency switch** `ncpu=-N`, ~10–15% cost | **ordered-reduction flag**, documented cost |
| **Solver interface** | one `mp_mode` knob, two backends | registered **factory string**, B&S agnostic | independent `LSOLVR` knob | single `system` verb / factory seam |
| **Big-3D escape hatch** | iterative PCG (narrow: well-conditioned solids only) | AMGCL / ML-AMG (default *distributed*) | **BLR → direct-quality PCG preconditioner**; out-of-core | MUMPS **BLR**; iterative deferred/data-justified |

**Three synthesis lessons this ADR adopts:**

1. **Solver tier and assembly-threading tier are orthogonal, independently-selected layers — all
   three codes keep them separate.** PARDISO/MUMPS own the *solve*; the OpenMP assembly layer helps
   *assembly* regardless of solver. Don't couple them. This is exactly why Lane 3 is its own effort,
   and why PARDISO+MUMPS (which already thread their own factorization) should own the solve while we
   spend engineering on the element/`update` loop.

2. **For the scatter race, prefer freeze-sparsity + atomic-scatter (Kratos) over graph coloring.**
   Kratos's *proven* answer: build the CSR graph once (per-row locks in the cheap symbolic pass),
   then thread the numeric assemble with a plain parallel-for and a per-entry `AtomicAdd` into the
   shared, pre-sized matrix — simpler than coloring and cheaper than thread-private matrices
   (memory ∝ nnz × threads). LS-DYNA's *element-ordering into conflict-free groups* is the stronger
   variant (kills the write conflict at the source **and** enables SIMD). Both beat "color then
   atomic the leftovers." **Hard rule (Kratos): never reallocate/resize during the threaded phase —
   only atomic-increment existing scalars.** Start with the explicit diagonal-SOE private-buffer
   reduction (fewest moving parts, no factorization interaction), design an *ordered* variant from
   day one, and expose **determinism as a flag with a documented cost** (LS-DYNA `ncpu=-N`).

3. **Symmetric-default + explicit override; BLR as the cheap direct→preconditioner escape hatch;
   keep element code parallelism-oblivious behind one SOE/B&S seam.** Default to symmetric
   factorization but *expose* the unsymmetric choice (the fork's contact/non-associated tangents are
   genuinely unsymmetric — don't over-auto-sniff, per Kratos). All three warn that the **distributed
   direct solve is the cluster scaling ceiling** (it does *not* scale like assembly or explicit) —
   so budget for good decomposition/ordering (METIS), reach for **MUMPS BLR** as the first
   memory/scaling relief, and keep an iterative+AMG path on the roadmap as a *data-justified*
   escape hatch, not a default. Portability (Kratos): the element `getTangent`/`getResisting` layer
   stays identical serial vs MPI; only the sparse space + assembler + solver swap.

## 5. Decision

1. **Portfolio, not one solver.** Desktop → **MKL PARDISO** (shared-memory, no MPI). Cluster →
   **MUMPS** (MPI, + BLR + hybrid). Each regime uses its strongest tool; both are ~done.
2. **Single portable verb.** A `system` alias / `-auto` policy resolves to PARDISO in a serial
   build and MUMPS in an MPI build, so the same model script runs on laptop and cluster unchanged.
3. **Threaded assembly is the real prize** but a separate, staged, measurement-gated effort — it is
   the only thing that helps the primary frame lanes and the threads-per-node cluster half.
4. **Do NOT** build serial/`libseq` MUMPS, a GPU solver offload, a hand-rolled Krylov/precond, or
   OpenMP-by-default (ADR-40 anti-goals stand).

### Cross-cutting levers (pay in both regimes)
- **Symmetric factorization, default-on + explicit override** (PARDISO `mtype ±2` / MUMPS
  `-matrixType 2` — ~2× time+memory everywhere the tangent is symmetric). Per §4: default symmetric,
  but *expose* the unsymmetric path rather than over-auto-sniffing — the fork's contact/non-
  associated tangents are genuinely unsymmetric and must be able to select `pardiso_lu`/unsym MUMPS.
- **Factorization reuse** driven off the SOE `factored` flag (present in MUMPS; add PARDISO phase-33)
  — pays under ModifiedNewton/Initial/IMPL-EX, not full Newton.
- **BLR** (MUMPS `ICNTL35`) as the first memory/scaling relief on large 3D, and the cheap
  direct→preconditioner bridge (LS-DYNA precedent) if an iterative path is ever needed.

## 6. Sequencing & gates

- **P1 — PARDISO desktop (small, mostly done).** Wire `system Pardiso` + serial-build link;
  symmetric + reuse. **Gate:** beats UmfPack on Lane B at 4 threads; byte/1e-12 oracle vs UmfPack.
- **P2 — MUMPS cluster tuning.** `-matrixType 2`, `-BLR`, hybrid ranks×threads. **Gate:** the
  ADR-74 rung harness shows a symmetric/BLR win at fixed np; no accuracy regression.
- **P3 — single verb / `-auto`.** Portability polish once P1/P2 land.
- **P4 — threaded assembly (own effort, staged).** (a) scope `elem.update`; (b) de-static kernels;
  (c) explicit private-buffer reduction (ordered variant from day one); (d) **atomic-scatter on the
  frozen `formTangent` sparsity** (Kratos pattern; coloring/element-ordering only if it contends);
  (e) thread the `update` loop. **Gate each stage** on >40% element fraction + oracle correctness.
  Likely its own sub-ADR.

### Bench matrix (the decider)
- **Desktop:** Lane B + one larger 3D solid × {UmfPack baseline, PARDISO @ 1/2/4/8 threads,
  ±symmetric} — median-of-7, threads pinned, via `Ladruno_files/testbed/perf/runner.py`.
- **Cluster:** same models × {MUMPS np-sweep, ±BLR, ranks×threads} — reuse the ADR-74 two-sweep
  method (fixed-np rung + fixed-V np-sweep).

## 7. Correctness & constraints
- **FP determinism:** threaded reduction changes summation order → last-bit drift. Precedent (§4):
  LS-DYNA ships this as an **explicit consistency switch** (`ncpu=-N`, thread-count-independent
  accumulation order, ~10–15% cost) defaulting to the fast/non-deterministic mode. Adopt the same
  shape — an **ordered-reduction flag** with a documented cost — and design the Lane-3 reduction so
  an ordered variant exists from day one. Against the fork's byte-identical/1e-12 oracle discipline,
  the ordered mode is the QA/regression path; the fast mode is opt-in for production speed.
- **LP64 ceiling:** the MUMPS build is LP64 (`-Dintsize64=OFF`), nnz < 2³¹; ILP64 is flagged
  non-trivial in `BUILD_GOTCHAS.md §3` and stays deferred. PARDISO (MKL LP64) shares the ceiling.

## 8. Anti-goals (inherited + lane-specific)
Serial/`libseq` MUMPS · GPU *solver* offload · hand-rolled Krylov/preconditioner · SIMD before
algorithmic work-removal · OpenMP-by-default or implicit colored-scatter before a measured >40%
element fraction · ParMETIS/`cluster_sparse_solver` before a measured deck justifies it.

## 9. Open questions
- `cluster_sparse_solver` (MKL distributed PARDISO) as a *future* unification of both regimes on one
  dependency — attractive, but would replace a proven MUMPS-parallel path; deferred, not chosen.
- Does METIS ordering link into the bundled MUMPS build (fill quality on large 3D)? Confirm.
- Iterative (PETSc AMG-CG) for large *well-conditioned* 3D — data-justified only; softening tangents
  fight preconditioners. Out of scope here.

## 10. Ledger / banner
No source touched yet (scoping ADR). When P1 lands: `LEDGER_implementations.md` row for the PARDISO
SOE/solver + `system Pardiso`; a `banner_features.txt` line; class tag reuse of
`LinSOE_TAGS_PARDISOGenLinSOE 99990` (already in `classTags.h`).
