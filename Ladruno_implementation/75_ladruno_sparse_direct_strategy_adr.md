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

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

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
  input, so `PARDISOGenLinSOE` must gain half-storage (the path `MumpsSOE` already has).
  **DONE in P1d, and the hedge is resolved IN FAVOUR of symmetric:** 1.94-1.96x faster than UmfPack
  (~1.25x vs unsymmetric PARDISO) and **-41.8% peak memory**, bit-identical answers (`phase1/RESULTS_p1d_symmetric.md`). The
  caution was well-placed but pointed at the wrong culprit — `SparseSYM` being **2.10× SLOWER** than
  unsymmetric UmfPack (`phase1/RESULTS_laneB_baseline.md`) was *its* implementation quality, not a
  property of symmetric storage. Kept as the standing lesson: the measurement was cheap and the prior
  was wrong in both directions at once.
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
- **➜ NOW ITS OWN SUB-ADR: [[75b_ladruno_threaded_assembly_adr]]** (opened 2026-07-25). It settles the
  two policy questions that constrain everything downstream — the **ordered-reduction / determinism CI
  policy** (§3 there: per-loop-class; ordered mode is the CI gate and the default when threading is on,
  fast/atomic mode opt-in and forbidden on oracle paths — LS-DYNA's *shape* with the default
  **inverted**, because this fork's QA is exact rather than tolerance-based) and the **scatter remedy**
  (§4 there: freeze-sparsity + atomic scatter as settled, with the new checked fact that OpenSees
  *already* freezes the sparsity — `zeroA` never touches `colA`/`rowStartA` — so the race is one
  `A[k] +=`). Stage L3-0 is measured.
- ~~**First move is instrumentation, not code:** ADR-40b's #1 gap is that the `elem.update` loop
  (force-based / IMK interior iteration — the frame lanes' real cost) is unscoped. Measure it first.~~
  **CORRECTED — this premise was stale.** The `elem.update` instrument was built and shipped
  **2026-07-06** (`Domain.cpp:2394`/`:2403`, macro `ProfilerMacros.h:114`), and ADR-40b's own addenda
  already used it. The real gap was that the >40% gate had never been evaluated **per threadable
  loop** — a profiler *phase* like `formTangent` bundles the threadable kernel with the
  non-threadable `addA` scatter. **Closed by L3-0**
  (`Ladruno_files/testbed/perf/lane3/RESULTS_l3a_update_scope.md`).
- **Measured Lane-3 gate (L3-0):** the >40% element gate **fails on Lane B under UmfPack (35.8%
  kernel) and passes under `Pardiso -matrixType 2` @4T (74.9%)** — solve falls 55.9% → 16.9%, so
  **Lane 1's win is what creates Lane 3's business case**. Lanes A and D are element-bound regardless
  (81.9% / 86.8%). `Domain::update` alone is 80.8% (lane A) / 54.4% (lane D) and has **no FP reduction
  at all**, making it both the biggest frame-lane lever and the only one with zero determinism cost.

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

### 5.2-bis — Symmetric is opt-in, NOT default (revised by P1d measurement)

§5's original wording followed the Kratos precedent to *"default symmetric, expose the override"*.
**P1d reverses the default while confirming the lever.** `-matrixType 2` is now measured at 1.35×
faster and −41.8% peak memory (`RESULTS_p1d_symmetric.md`), yet the default remains `0`:

- half-storage reads only the `col >= row` half of each element matrix — **no averaging, no
  detection** — so on a genuinely unsymmetric tangent it silently solves a *different system*;
- this fork has such tangents in production: `LadrunoContact`, non-associated flow, follower loads,
  `LadrunoUP` (whose ledger row already says "needs a general solver");
- the Kratos precedent assumes an element library that is symmetric by construction. This one isn't.

So: **explicit token, loudly documented, with `-stats` to make the payoff self-evident.** Authors of
symmetric models should pass `-matrixType 2`; prefer `2` over `1` (SPD/Cholesky) because `1` fails
outright on the indefinite tangent any softening or buckling model eventually has — the `-4`
zero-pivot message now says exactly that.

### Cross-cutting levers (pay in both regimes)
- **Symmetric factorization — SHIPPED for PARDISO (P1d), opt-in** (PARDISO `-matrixType 1|2` →
  `mtype ±2` / MUMPS `-matrixType 2`). **No longer a hypothesis: ~1.25x vs unsymmetric PARDISO (1.94-1.96x vs
  UmfPack) and -41.8% peak memory, measured over two sweeps.** The old "~2× time+memory, but `SparseSYM` was 2.10× slower so measure it" hedge is
  resolved — that was `SparseSYM`'s implementation, not symmetric storage. MUMPS `-matrixType 2`
  remains **unexercised** on the cluster and inherits a strong prior from this result.
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
  0.0)**. ~~so threading introduced no FP drift and the §7 determinism concern is Lane-3-only.~~
  **⚠ CORRECTED by P1f:** the *accuracy* claim stands, but "no FP drift / Lane-3-only" does not — it
  rested on one run per thread count. With **10 runs of one binary at a fixed thread count**,
  threaded PARDISO returns **2 distinct results in a 5/5 split** (~1 ULP) at 4 threads while being
  10/10 identical at 1 thread; reproduced on a pre-P1f binary, so it is MKL's, not ours, and it is
  size-dependent (an 8³ model is stable even at 4T, which is why the single-sample 15³ check
  passed). **So §7's determinism concern is NOT Lane-3-only** — a byte-identical gate must pin
  `MKL_NUM_THREADS=1` or use `iparm[33]` CNR (legal here since we set `iparm[1]=2`; not yet exposed).
  Scaling flattens past 4 threads (1.50×→1.58×, memory-bandwidth-bound) ⇒ **4 threads is the
  recommended desktop default**. The measured 1.76× sits just under the predicted ~2.2× Amdahl
  ceiling, confirming the residual ~34% is Lane-3 territory. *(Original scoping, for the record:)* Compile-verify the prototype; wire
  `system Pardiso` + serial-build link; **re-architect factorization reuse** (persist `pt`, drop the
  per-solve release); add symmetric SOE half-storage; fix the `iparm` leak.
  **Gate (now a measured number): beat UmfPack's 22.711 s on Lane B** at 4 threads *and* threading
  verified engaged; bit-identical/1e-12 tip displacement vs the locked baseline
  (`phase1/RESULTS_laneB_baseline.md`).
- **P1d — symmetric PARDISO (`-matrixType`). ✅ DONE, WINS ON BOTH AXES**
  (`phase1/RESULTS_p1d_symmetric.md`). `PARDISOGenLinSOE` gained upper-triangle half-storage and the
  solver now *derives* `mtype` from it (11 / 2 / −2), so storage and factorization mode cannot
  disagree. `system Pardiso -matrixType 0|1|2` (+ `-symmetric`/`-spd`), in **both** interpreters —
  P1b had registered the verb only in `OpenSeesCommands.cpp`, so `OpenSees.exe` never had it; the Tcl
  chain is a separate if-ladder and is wired here.
  **Measured, on Lane B, same binary, interleaved, TWO independent sweeps:** `-matrixType 2` is
  **1.94-1.96x faster than UmfPack** at 4 threads (1.62-1.63x single-threaded), reproducible to
  +-1%; versus *unsymmetric PARDISO* it is **~1.25x** (range 1.24-1.35x, the unsym anchor being the
  noisy term). Tip displacement is **bit-identical** in every configuration of both runs. **The `SparseSYM`-based worry (§3, "symmetric
  ≠ automatically better", 2.10× SLOWER) is refuted for PARDISO** — it was an implementation-quality
  artifact of `SparseSYM`, not a property of symmetric storage.
  **And the memory result is the bigger one - `-stats` (new, `iparm[14]/[15]/[16]`, mirroring the
  MUMPS `-stats`) shows peak memory -41.8%** (105.60 -> 61.48 MB; SPD -46.4%), stored nnz −49.3%
  (exact by construction), factor nnz −47.3%. **Directly contrast P2b: BLR cut the stored factors
  21.8% but peak only 4.6%.** Symmetric shrinks the *fronts themselves*, so the fit/no-fit number
  actually moves — and it is **exact**, so unlike BLR it is legal on byte-identical/oracle paths.
  Against the P1c capability wall this is the largest memory lever ADR-75 has produced.
  **Default stays `0` (unsymmetric)** — see §5.2-bis.
- **P1e — factorization-preconditioned CGS (`-krylov <L>`). ✅ DONE, AND IT COVERS THE CASE P1a CANNOT**
  (`phase1/RESULTS_p1e_krylov.md`). P1a's reuse gate keys on the SOE `factored` flag, which answers
  "is A unchanged" — so it pays under ModifiedNewton/Initial/Krylov/IMPL-EX and **nothing under full
  Newton**. `iparm[3]` is the complementary axis: it reuses the retained L/U as a *preconditioner*
  for a tangent that HAS changed. `system Pardiso -krylov <digits>` in **both** interpreters.
  **Measured on Lane B under full Newton:** **1.51× vs direct PARDISO at 50.7k DOF / 4 threads**
  (1.22× @11.5k, 1.30× @26.5k — *rising with N*), 1.51-1.94× single-threaded, and **1.57× stacked
  with P1d `-matrixType 1`** (vs 1.34× for half-storage alone; the levers compose sublinearly since
  both attack the same factorization cost). Tip displacement bit-identical in every row.
  **⚠ Two size/thread traps recorded:** threading *erodes* the win (factorization parallelizes, the
  triangular solves in each CGS iteration do not) so 11.5k/4T reads a misleading 1.22×, and the
  trend only reasserts itself with N — quoting the small-model number would repeat the ADR-40b
  mistake. **⚠ Scope limit: Intel documents `K=2` for symmetric POSITIVE DEFINITE only, so
  `-matrixType 2` (`mtype -2`) — the right choice for a softening/buckling tangent — cannot use this
  lever at all** (warns, falls back to direct). The nonlinear-softening class that most needs solver
  speed is therefore the least served; the 1.57× composition figure is SPD-only.
  Three implementation subtleties, all banked in `LEDGER_quirks.md`: phase **23** is mandatory (the
  automatic direct fallback is documented for 23 only; under 33 the same failure is just `error=-4`);
  a CGS *win* leaves the stored factors **stale**, so the phase-33 shortcut must be forbidden
  afterwards (new `factorsCurrent`, distinct from `haveFactors`) or a later same-A solve silently
  answers the previous tangent; and `iparm[3]` must stay 0 for a pattern's first numeric pass.
  The stale-factor guard was verified by **deliberately removing it and re-measuring** — it is
  load-bearing, but Lane B's tip displacement (1.2e-11) and iteration count (1.01×, and only 1.03×
  even swept to perfect plasticity) do **not** detect it, because Newton's residual test launders a
  stale-factor solve into a slower quasi-Newton. The check that does discriminate (`-krylov` vs
  *direct* under the same ModifiedNewton: 0.0 vs 1.2e-11) is what `p1e_smoke.py` now asserts.
  **Adversarial review also refuted a second hypothesis of mine — preconditioner aging.** CGS never
  refactors while it wins, so the preconditioner stays the step-1 factorization and the `-stats`
  iteration drift (1→3→4→5) looked like a short-benchmark bias. Measured: 60 steps gives **1.34×**
  vs 15 steps' 1.30×, and over 150 steps / 340 solves the iteration count **plateaus at 4–5**. No
  refresh policy needed. **⚠ Top residual risk: the FALLBACK branch has never executed** — 340
  solves, `L` up to 9, perfect plasticity at double load (tip 55.5, 18 CGS iterations): zero
  fallbacks. Its decode is unverified; mitigation is that the branch's only correctness-relevant
  action (`factorsCurrent = true`) is unconditional and right for both `iparm[19] < 0` and `== 0`,
  so a bad decode corrupts a diagnostic string, not an answer.
  **FOLLOW-UP (softening / limit-point sweep, `p1e_softening_probe.py` + `p1e_prepost_probe.py`):**
  the "cost when CGS fails" risk is **CLOSED — the failure never happens.** On softening tangents
  (`Hiso` −2000…−20000) driven to a limit point under `-matrixType 0`, CGS works ~4× harder (15–19
  iterations vs 4–5) and `-krylov` is **still 1.08–1.16× faster**; zero fallbacks again. So the
  fallback branch is **effectively unreachable**, not merely untested — leave it as documented dead
  code rather than engineer a synthetic trigger. **⚠ What replaces that risk is narrower and real:
  `-krylov` can change the POST-PEAK branch.** Stopping at increasing step counts shows the two are
  **bit-identical through the whole physically meaningful range** (ux 0.394→4.725, steps 5–32) and
  separate only *after* the limit point, where `LoadControl` on a softening structure is ill-posed
  and both answers are non-physical (direct 6435 vs krylov −5502 at step 35). Post-limit-point path
  chaos amplifying a 1-ULP difference, not a defect — but it is a **reproducibility** hazard: keep
  `-krylov` off when the deliverable is a post-peak path (progressive-collapse / AEM lane).
- **P1f — the `addA` scatter quadratic, and the determinism claim it broke. ✅ DONE**
  (`phase1/p1f_adda_ab.py`, `p1f_determinism_probe.py`). `PARDISOGenLinSOE::addA` rescanned the
  whole CSR row for each of `idSize²` element entries — `O(idSize² × rowlen)`, which ADR-75b's L3-0
  profile measured at **1699 ms of an ~11.9 s step and 1.28× slower than UmfPack's** scatter.
  Replaced by a binary search, legal because `setSize` **enforces** ascending CSR rather than
  assuming it. **Measured 1.098× at 26.5k DOF** (`-matrixType 1`: 1.029×, less because half-storage
  already scatters fewer entries into shorter rows). **This was the right next lever precisely
  because P1a/P1d/P1e worked:** every solver win raises the assembly share (L3-0 put the element
  fraction at 35.8% under UmfPack → 74.9% under PARDISO), and threading a quadratic would only buy
  a parallel quadratic — so this lands *before* Lane 3's L3-1.
  **⚠ And it broke a published claim.** The A/B reported "ux DIFFERS" at 4 threads. The tempting
  reading was "the change is not exact"; the correct question was whether each binary reproduces
  *itself*. It does not: **10 runs of ONE binary at 4 threads give 2 distinct results in a 5/5
  split** (~1 ULP), on the *pre-P1f* binary too, while 1 thread is 10/10 identical and old==new.
  So the change is exact, and **threaded MKL PARDISO is not byte-reproducible run-to-run** — which
  refutes P1's "bit-identical at every thread count ⇒ the §7 determinism concern is Lane-3-only".
  That conclusion came from **one run per thread count**, a design that cannot distinguish
  deterministic from lucky. Corrections landed in `RESULTS_p1_pardiso.md` and §P1 above;
  `LEDGER_quirks.md` carries the full entry. Mitigation: pin `MKL_NUM_THREADS=1` for byte-identical
  gates, or expose `iparm[33]` CNR — available to us since we set `iparm[1]=2` (Intel forbids CNR
  only for `iparm[1]=3`), **not currently wired**.
  **Opt-in, off by default** — byte-identical to P1d when absent.
- **P2 — MUMPS cluster tuning. 🟡 `-BLR` SHIPPED, effect NOT yet validated at scale**
  (`phase1/RESULTS_p2_blr.md`). `system Mumps -BLR <eps>` (+ raw `-ICNTL35`/`-CNTL7`) is wired,
  propagates to subordinate ranks via `sendSelf`/`recvSelf` (without which rank 0 would factor BLR
  while the others factored full-rank), and **demonstrably engages** above a size threshold
  (rel diff 1.86e-12 at `eps=1e-4`; a *smaller* model showed BLR silently not engaging at all).
  **⚠ Measured counterintuitive result: at ~30k DOF / np2 BLR is 1.70× SLOWER at `1e-9` and 3.17×
  slower at `1e-4`** — compression overhead exceeds flop savings on small fronts, AND in a
  *nonlinear* loop a looser tolerance returns a less accurate correction so Newton needs more
  iterations (more solves). **So BLR is a MEMORY lever, not a speed lever**, its justification is the
  P1c capability wall, and it stays opt-in/off-by-default. **P2b `-stats` SHIPPED** (`system Mumps ... -stats` dumps
  `INFOG(9)/(21)/(22)`, `RINFOG(3)`, BLR `RINFOG(14)/(15)`), which finally makes compression
  observable — and it produced a **non-obvious measured result: BLR shrinks the FACTORS but barely
  moves PEAK MEMORY.** At `eps=1e-4`: factor entries **−21.8%**, BLR flops **−45%**, but peak
  **MB/proc only −4.6%**; at `eps=1e-9` it is strictly worse (**+8.4% MB/proc**, no flop saving).
  Peak factorization memory is dominated by the **active frontal/working space**, not the stored
  factors — so "BLR saves memory" is true of factor storage and largely false of the allocation that
  actually decides whether a model fits. **At ~32k DOF/np2 BLR is a win on no axis.** This bounds the
  small end only and does *not* refute BLR on production-size fronts. **Still open:** the crossover —
  run a production deck with `-stats`, BLR on/off, comparing `INFOG(21)` and wall time.
  `-matrixType 2` (symmetric) and hybrid ranks×threads remain untouched.
- **P3 — explicit-verb portability polish.** Clear build-time errors + docs (no `-auto`; §12) once
  P1/P2 land. **Preceded by the P0 trade study below** if unify-on-MKL is chosen.
- **P4 — threaded assembly. ➜ SPLIT OUT to [[75b_ladruno_threaded_assembly_adr]]** (2026-07-25); read
  that for the staging (L3-0…L3-5), the determinism policy, and the risk register. **L3-0 (measure) is
  DONE** — `lane3/RESULTS_l3a_update_scope.md`. Status of the original sketch: (a) "scope
  `elem.update`" was **already shipped 2026-07-06**, so L3-0 instead produced the per-loop gate table
  the >40% gate actually needs; (d) the frozen-sparsity precondition turns out to **already hold** in
  OpenSees; (e) threading the `update` loop is **promoted from last to first** — it is the largest
  frame-lane fraction *and* the only loop with no FP reduction, hence bit-identical when threaded.
  Two hazards L3-0's review surfaced that the sketch did not anticipate: `ops_TheActiveElement` is a
  mutable global written inside that loop, and ~13 element files (incl. `LadrunoRigidBody`) write
  **shared node trial state** there — an *ordering* race no reduction policy can fix.
  *(Original sketch, for the record:)* (a) scope `elem.update`; (b) de-static kernels;
  (c) explicit private-buffer reduction (ordered variant from day one); (d) **atomic-scatter on the
  frozen `formTangent` sparsity** (Kratos pattern; coloring/element-ordering only if it contends);
  (e) thread the `update` loop. **Gate each stage** on >40% element fraction + oracle correctness.

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
- ~~**PARDISO instrumentation parity**~~ ✅ **CLOSED 2026-07-27,
  [#667](https://github.com/nmorabowen/OpenSees/pull/667) (`62768d1f1`).** Raised by the
  [[76_ladruno_tangent_reuse_adr]] issue report, tracked here. `PARDISOGenLinSolver` carried no
  `OPS_PROFILE` scopes where `UmfpackGenLinSolver` had `soe.symbolic`/`soe.factor`/`soe.trisolve`,
  so on a PARDISO run `linearSolve` was one opaque block. Now 7 brackets: the three UmfPack names
  (reused verbatim so a cross-solver profile lines the phases up), a PARDISO-only `soe.cgs` for
  phase 23, `dc.s.fill`/`dc.s.verify` on the `setSize` CSR build, and a DEEP-gated `soe.addA`.
  Verified by `tests/test_pardiso_solver.py::test_profiler_brackets_present`; measured split
  published in `Ladruno_files/testbed/perf/phase1/RESULTS_p1h_phase_split.md`, and its size
  trend over 11.5k-136k DOF in `RESULTS_p1j_size_trend.md` (P1j): `fac/(fac+tri)` rises
  **71.4% -> 93.7%**, factorization **15.0% -> 56.2% of step**, driven by a measured exponent
  gap of ~0.64-0.72 between `soe.factor` and `soe.trisolve`. Every factorization-reuse lever
  is therefore worth materially more at production scale than the 11.5k gate suggested.
- ~~**the other half of that item:** the profiler HDF5 run attributes record `threads=1`/`nElem=0`
  regardless of configuration~~ ✅ **CLOSED 2026-07-27 (P1i, this PR).** `Profiler.cpp` set
  `m.threads = threads_.size()` — profiler-*registered* threads, 1 on any single-threaded command
  layer regardless of `MKL_NUM_THREADS`; `nElem`/`nNode` were promised by a comment and filled by
  nobody. Now `resolveRunThreads()` (env, with the registered count as a floor) plus a
  `mkl_get_max_threads()` override under `-D_PARDISO`, and `nElem`/`nNode` filled at all four
  `buildMeta()` call sites — `profiler report` and `checkpoint`, each existing twice across the
  Python and Tcl ladders. **`nnz` stays 0 by decision** (no size-agnostic `LinearSOE` accessor;
  filling it means a virtual on an upstream base class). **`nSteps` was never a bug** — it derives
  from the per-step series, so 0 without `-perStep` is correct.
  ⇒ **ADR-75 §9 has no remaining instrumentation items.**

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
   that refills `A` without clearing `factored` silently solves a stale factor. The assembly-side
   mirror of this axis (skip the re-assembly so `factored` stays set) was examined and **WITHDRAWN**
   in [[76_ladruno_tangent_reuse_adr]] — its §4.1 invalidator inventory is the definitive list of
   paths that change `A` without any observable event, and its §4.3 records that combining `-krylov`
   with any assembly-skip (`ModifiedNewton -factoronce`) is a pessimization (a CGS win leaves
   `factorsCurrent == false`, so every subsequent solve of the unchanged matrix re-enters phase 23).

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
