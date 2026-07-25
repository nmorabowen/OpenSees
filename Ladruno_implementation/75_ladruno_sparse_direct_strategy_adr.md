---
title: Sparse-direct solver strategy — desktop (PARDISO) + cluster (MUMPS) + threaded assembly
project: Ladruno
status: proposed — scoping (code-verified inventory; measure-gates defined; cross-framework precedent folded; revised post-adversarial-review §12)
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

A direct solver wins **one** lane *of this table*. Threaded element assembly is the dominant cost in
four of five. This is the measure-first spine of the whole ADR — **but read the production-regime
correction immediately below before using it to prioritize.**

### ⚠ Production-regime correction (2026-07-24, from the fork owner)

ADR-40b's lane table is a *breadth* survey, and earlier drafts of this ADR wrongly inferred from it
that the fork's **primary** workload is fiber frames. It is not. The production regime is
**huge solid nonlinear models** — i.e. **Lane B is the primary lane**, not a side lane. Consequences,
all of which *raise* the value of the solver lanes:

1. **The solver lanes matter more, not less.** Lane B is 66% `linearSolve` at a mere 11.5k DOF. A 3D
   sparse-direct factorization scales ~O(N^1.5–2) in flops (worse than the ~O(N) assembly), so at
   10⁵–10⁶⁺ DOF the solve fraction **grows** — the 66% is a *floor* for production sizes, and
   PARDISO/MUMPS work compounds with model size.
2. **The measurement basis is under-powered.** Every gate in this ADR was set on an 11.5k-DOF model
   (and ADR-40b spans 0.7–11.5k). That is small by two-plus orders of magnitude versus production.
   **✅ NOW MEASURED — P1c (`phase1/RESULTS_p1c_scaling.md`) confirms the concern was right:** the
   PARDISO win **compounds** with size — **1.61× (11.5k) → 2.15× (26k) → 3.40× (51k) DOF** — because
   UmfPack scales ~O(N²) while PARDISO scales ~O(N^1.45). The P1 headline understated the win for
   this fork's real regime by ~2×. **And a capability wall appeared: UmfPack ran OUT OF MEMORY at
   86,490 DOF while PARDISO solved it in 30.4 s and 136,080 DOF in 68.6 s** — so PARDISO raised the
   largest solvable single-machine model by ≥1.6× in DOF, ceiling untested. That is worth more than
   any speed ratio here: models that previously forced the cluster may now fit on a workstation.
3. **Memory becomes the binding constraint, so BLR is promoted.** At huge 3D scale the direct
   factor's memory — not its time — is what stops a run. MUMPS **BLR** (`ICNTL35`) and out-of-core
   move from "nice tuning" to a primary P2 item, with the accuracy caveat in §5 handled explicitly.
4. **The LP64 ceiling stops being theoretical.** `-Dintsize64=OFF` caps nnz at 2³¹. A large 3D
   factor can genuinely exceed that, so the ILP64 question (§7, currently deferred as "non-trivial")
   needs a *measured* nnz headroom check on a real production deck.
5. **Nonlinear ⇒ factorization reuse pays.** Many Newton iterations per step is exactly the regime
   where the P1 reuse work (and ModifiedNewton/IMPL-EX) converts into wall-clock.
6. **Lane 3 is still the other half.** Even in Lane B, element state determination was 30% — and the
   fork's expensive solid materials (`LadrunoConcrete3D` CDPM2, `LadrunoJ2`) make that fraction
   heavier on real decks than on this elastic-ish benchmark.

## 2. Code-verified inventory (what exists today)

Registered `system` verbs (from `OpenSeesCommands.cpp`): `UmfPack`, `SuperLU`, `SparseSYM`,
`SparseGEN`, `FullGeneral`, `BandGeneral`, `Diagonal`, `MPIDiagonal`, `Mumps`, `Petsc`, … .

| Solver | Files | State | Threaded? | MPI? |
|---|---|---|---|---|
| **UmfPack** | `linearSOE/umfGEN/` | wired; Lane-B baseline | BLAS only | no |
| **MUMPS** | `linearSOE/mumps/` | wired **MP/SP only** (`_MUMPS` gated by `if(MPI_FOUND)`) | BLAS + `ICNTL16` | **yes** |
| **MKL PARDISO** | `linearSOE/pardiso/PARDISOGenLin{SOE,Solver}` | **2019 prototype, UNWIRED & never compiled here**; re-factors every solve | native OpenMP (prototype hints `iparm[2]=1`) | no |
| **PETSc** | `linearSOE/petsc/` | wired but unbuilt on this toolchain | yes | optional |

Key facts established during scoping:
- `MumpsSolver.cpp:50-52` already `#include <libseq\mpi.h>` for a non-MPI build, and
  `OPS_MumpsSolver()` has a serial branch (`:4997`) — a **serial MUMPS path exists in source** but
  is never built (and `libseq`/`libmpiseq` is **not bundled**; the Windows build via
  `scivision/mumps` v5.5.1.5 is always **Intel-MPI, LP64** — nnz capped at 2³¹).
- `PARDISOGenLinSolver.cpp:23` uses `<mkl_pardiso.h>` — it is **MKL PARDISO**. But it is a **2019
  contributed prototype (M. Salehi), not in any build target**, so its compile status against the
  current MKL headers / `LinearSOE` interface is **unverified**. What is actually proven is that
  FEAST's *own* `pardiso()` wrapper links (`FeastEigenSolver.cpp:132`) — **not** this class. The
  prototype (revised by adversarial review, §12): `solve()` runs phase 11→22→**33→−1 every call**
  (`PARDISOGenLinSolver.cpp:156-203`) — it re-does the **METIS symbolic reorder + numeric factor and
  then frees all memory on every solve**; `pt[64]` is a stack local so it *cannot* persist; `mtype`
  is **hardcoded 11 (unsymmetric)** and the SOE stores **full unsymmetric CSR** (no upper-triangle
  half-storage); and it leaks `iparm` each solve (`:208` delete commented out). So it is a *starting
  point* (CSR conversion + phase skeleton), **not** "90% done."
- Factorization reuse is already correct in MUMPS (`job=3` solve-only gated on the SOE `factored`
  flag; `job=5` factor+solve otherwise) — unlike UmfPack, which re-factors every solve.
- `system Mumps` already parses `-matrixType 0|1|2` (unsym / SPD / general symmetric) and
  `-ICNTL7` (ordering) / `-ICNTL14` (workspace) / `-commSplit` (ADR-43 sub-communicators) — but the
  symmetric path is **unexercised** by the perf scripts.

## 3. The three lanes

### Lane 1 — PARDISO owns the desktop (shared-memory direct)
Native-OpenMP MKL PARDISO in the plain serial `opensees.pyd`: all cores on the factor/solve, **no
MPI**, **zero new dependency** (MKL already linked). The in-tree prototype is a head-start, **not**
finished (§2, §12) — realistic P1 work, in order:
- **Prove it compiles** against current MKL headers / `LinearSOE`; register `system Pardiso`; lift a
  `_PARDISO` flag + link into the **serial** `OpenSees`/`OpenSeesPy` targets (mirror `_MUMPS`, minus
  the MPI gate). Fix the per-solve `iparm` leak.
- **Re-architect for factorization reuse** — persist `pt[]` as a member, split symbolic (11) /
  numeric (22) / solve (33), stop the per-solve phase −1 release, and gate solve-only on the SOE
  `factored` flag (match MUMPS). This is a *restructure*, not a flag.
- **Symmetric path is a SOE change, not just `mtype`** — PARDISO `mtype ±2` needs upper-triangle
  input, so `PARDISOGenLinSOE` must gain half-storage (the path `MumpsSOE` already has). **Measure,
  never assume, the symmetric win:** the one symmetric solver we *can* measure (`SparseSYM`) is
  **2.10× SLOWER** than unsymmetric UmfPack on Lane B (`phase1/RESULTS_laneB_baseline.md`) —
  implementation quality dominates the storage-format advantage. Treat "~2× from symmetric" as an
  unproven hypothesis, not a lever (§12).
- Confirm threading actually engages (`MKL_NUM_THREADS`; the prototype's `iparm[2]=1` is a red flag).
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
2. **Explicit solver verbs, not `-auto` magic.** (Revised — §12.) The author writes `system Pardiso`
   or `system Mumps` explicitly; portability comes from *documentation + a thin build-time guard*
   that errors clearly if you ask for a solver this build lacks — **not** from silent build-dependent
   resolution. In a research code where solver choice changes convergence and last-bit results,
   implicit resolution is a footgun (and contradicts the Kratos explicit-factory precedent in §4).
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
- **BLR** (MUMPS `ICNTL35`) as memory/scaling relief on large 3D — **but it is an *approximate*
  factorization** (low-rank truncation with an accuracy tolerance), i.e. a direct→preconditioner
  bridge, **not** an exact drop-in. It is off-limits for the byte-identical/1e-12 oracle paths and
  must be opt-in with a documented accuracy knob (§12).

## 6. Sequencing & gates

- **P0 — portfolio vs unify trade study. ✅ DONE** → portfolio confirmed
  ([[75a_p0_portfolio_vs_unify_trade_study]]); P1 unblocked, scope fixed to serial/shared-memory.
- **P1 — PARDISO desktop. ✅ DONE, GATE PASSED** (`phase1/RESULTS_p1_pardiso.md`). `system Pardiso`
  is built and working in the serial module (477 pardiso symbols, was 0). Lane B, same binary:
  **1.71× faster than UmfPack at 4 threads** (10.396 s vs 17.775 s), 1.76× at 8, and **1.19× even
  single-threaded**. Tip displacement **bit-identical to UmfPack at every thread count (rel err
  0.0)** — so threading introduced no FP drift and the §7 determinism concern is Lane-3-only.
  Scaling flattens past 4 threads (1.50×→1.58×, memory-bandwidth-bound) ⇒ **4 threads is the
  recommended desktop default**. The measured 1.76× sits just under the predicted ~2.2× Amdahl
  ceiling, confirming the residual ~34% is Lane-3 territory. *(Original scoping, for the record:)* Compile-verify the prototype; wire
  `system Pardiso` + serial-build link; **re-architect factorization reuse** (persist `pt`, drop the
  per-solve release); add symmetric SOE half-storage; fix the `iparm` leak.
  **Gate (now a measured number): beat UmfPack's 22.711 s on Lane B** at 4 threads *and* threading
  verified engaged; bit-identical/1e-12 tip displacement vs the locked baseline
  (`phase1/RESULTS_laneB_baseline.md`).
- **P2 — MUMPS cluster tuning.** `-matrixType 2`, `-BLR`, hybrid ranks×threads. **Gate:** the
  ADR-74 rung harness shows a symmetric/BLR win at fixed np; no accuracy regression.
- **P3 — explicit-verb portability polish.** Clear build-time errors + docs (no `-auto`; §12) once
  P1/P2 land. **Preceded by the P0 trade study below** if unify-on-MKL is chosen.
- **P4 — threaded assembly (own effort, staged).** (a) scope `elem.update`; (b) de-static kernels;
  (c) explicit private-buffer reduction (ordered variant from day one); (d) **atomic-scatter on the
  frozen `formTangent` sparsity** (Kratos pattern; coloring/element-ordering only if it contends);
  (e) thread the `update` loop. **Gate each stage** on >40% element fraction + oracle correctness.
  Likely its own sub-ADR.

### Bench matrix (the decider)
- **Desktop:** Lane B + one larger 3D solid × {UmfPack baseline, PARDISO @ 1/2/4/8 threads,
  ±symmetric} — median-of-7, threads pinned. **Harness written and ready:**
  `Ladruno_files/testbed/perf/phase1/laneB_solver_bench.py` (same Lane-B model as ADR-40b;
  interleaved rounds; probes each solver and records `unavailable` for the unwired ones; asserts a
  1e-9 tip-displacement cross-check so a timing is never reported for a wrong answer).
  **✅ EXECUTED 2026-07-24 — baseline locked** (`phase1/RESULTS_laneB_baseline.md`):
  **UmfPack 22.711 s = the P1 gate**; SparseSYM 2.10× slower; SuperLU 3.46× slower; all three
  **bit-identical** tip displacement (rel err 0.0). `Mumps` and `Pardiso` both fail at *runtime*
  with `WARNING unknown system type` — **empirical confirmation** of the §2 static finding that the
  desktop regime has no threaded sparse-direct solver today.
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
element fraction · ParMETIS before a measured deck justifies it · **`cluster_sparse_solver` /
distributed PARDISO** (P0-decided: MUMPS is mandatory for CMS/FEAST/PFEM regardless, so this only
*adds* a third solver family — [[75a_p0_portfolio_vs_unify_trade_study]]).

## 9. Open questions
- ~~**P0 DECISION — portfolio vs. unify-on-MKL.**~~ **CLOSED → portfolio confirmed**
  ([[75a_p0_portfolio_vs_unify_trade_study]]). Decisive finding: **MUMPS cannot be removed** — it is
  load-bearing for **CMS** (`LadrunoCMSMumps`, a `FATAL_ERROR` build gate at `CMakeLists.txt:645`),
  **distributed FEAST/modal** (`LadrunoDistBlockZKernel`), and **PFEM** (3 solvers). So unify-on-MKL
  would *still* link MUMPS and merely **add a third** solver family — its "one dependency, one test
  surface" premise is false. `cluster_sparse_solver` moved to §8 anti-goals; P1 is unblocked with
  scope fixed to serial/shared-memory only.
- Does METIS ordering link into the bundled MUMPS build (fill quality on large 3D)? Confirm.
- **Non-MKL desktop gap:** with serial-MUMPS descoped and PARDISO MKL-only, a non-MKL build
  (Zone-A Ubuntu / OpenBLAS) has **no threaded desktop solver** — UmfPack is the only fallback.
  Accept, or keep a non-MKL threaded option on the table.
- Iterative (PETSc AMG-CG) for large *well-conditioned* 3D — data-justified only; softening tangents
  fight preconditioners. Out of scope here.

## 10. Ledger / banner
No source touched yet (scoping ADR). When P1 lands: `LEDGER_implementations.md` row for the PARDISO
SOE/solver + `system Pardiso`; a `banner_features.txt` line; class tag
`LinSOE_TAGS_PARDISOGenLinSOE 99990` (already in `classTags.h`) — **but 99990 is off the fork's
33xxx Ladruno convention (upstream-prototype value); re-tag into range and LEDGER-check for collision
before shipping.**

## 11. Architectural risks (register)
The risk profile is **bimodal**: Lanes 1–2 are low-risk (self-contained `LinearSOE` back-ends);
essentially all architectural risk sits in Lane 3.
1. **Threaded assembly is a whole-codebase re-entrancy invariant** — every element/material kernel
   (many vanilla-upstream) must lose its `static`/shared scratch; one miss = silent, thread-only,
   nondeterministic wrong answers. *Highest.* Mitigate: explicit-diagonal path first, ThreadSanitizer,
   per-classTag gating.
2. **Threading vs the byte-identical QA discipline** — threaded FP reduction breaks byte-identical by
   construction; needs the ordered-reduction CI policy decided *before* any threaded code lands.
3. **Assembly loop is a central chokepoint** — the FE_Element scatter is shared by static/transient/
   eigen/sensitivity; keep the threaded path behind a default-off flag so serial stays byte-identical.
4. **Nested threading / oversubscription** — MKL solver threads × OpenMP assembly threads × MPI ranks;
   needs one coordinated thread-count policy or benches mislead.
5. **The unify seam isn't clean** — PARDISO CSR vs serial-MUMPS COO vs distributed `a_loc`; no
   templated sparse-space, so numberer/graph coupling leaks. "One verb" is a build-aware SOE factory.
6. **PARDISO ties desktop-perf to MKL** — non-MKL builds get no threaded desktop solver (see §9).
7. **Reward back-loaded onto the riskiest lane** — Lanes 1–2 only help 3D-solid; the frame-lane payoff
   lives entirely in Lane 3 (highest risk). Real failure mode: ship the easy lanes, stall before the
   one the primary models need.
8. **Stale-LU reuse trap (both solvers)** — reuse gated on `factored` assumes `A` untouched; any path
   that refills `A` without clearing `factored` silently solves a stale factor.

## 12. Revision log — adversarial review (post-merge, evidence-backed)
Corrections folded after reading `PARDISOGenLinSolver.cpp` (the review caught claims I'd made from a
header + grep, not the implementation):
- **PARDISO maturity overclaim → corrected.** It is a 2019 prototype never compiled in this build; it
  re-factors + frees memory **every solve** (`:156-203`), hardcodes `mtype=11`, has no symmetric SOE
  storage, and leaks `iparm`. "~90% done/proven to link" was false ("proven" was FEAST's own wrapper,
  not this class). P1 re-scoped small→**medium**; §2/§3-Lane1/§6-P1 rewritten.
- **`-auto` magic → dropped** for explicit `system Pardiso`/`Mumps` (contradicted the Kratos precedent
  it cited; hurts research-code reproducibility). §5.2, §6-P3.
- **Portfolio-vs-unify → elevated** from footnote to a **P0 decision** (the "proven MUMPS" objection
  weakened once PARDISO isn't nearly-free). §9.
- **BLR → caveated** as an *approximate* factorization, off-limits for oracle paths. §5, §9.
- **"~2× symmetric" → hedged** (sparse-symmetric time savings often <2×; measure). §3-Lane1, §5.
- **Added:** the §11 risk register, the non-MKL-desktop gap (§9), and the classTag-convention fix (§10).
- **Strategy direction survives** — desktop-threads / cluster-MPI / assembly-threading is sound and
  cross-framework-validated; the corrections are to effort estimates, one contradiction, and one
  under-argued decision, not to the architecture.
