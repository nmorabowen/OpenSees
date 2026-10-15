---
title: "WP-192 — threaded element loops for the Paper 1 element set (per-type opt-in, MPI ranks)"
project: Ladruno
type: performance work package (element lane)
status: "implemented, draft PR #979 — not merged; default off (1 thread)"
priority: medium
owner: nmora
related:
  - "[[107_ladruno_openmp_element_loop]] (WP-107: loop A, all-or-nothing allowlist, MPI refused)"
  - "[[75b_ladruno_threaded_assembly_adr]] (loops A-E, hazards H1-H7, the ordered-gather design)"
  - "[[188_ladruno_factor_reuse_adr]] (why the element loops become the floor)"
  - "paper repo: review/engstruct-2026/30_campaigns/timehistory/perf/PROFILE_4D.md (§2 phase table)"
tags: [performance, openmp, element, assembly, mpi, hybrid, wp-192]
updated: 2026-10-14
---

# WP-192 — threaded element loops for the Paper 1 element set

> [!summary] The short version
> With METIS, fewer ranks, factor reuse (WP-188) and possibly a faster factorization, the
> MUMPS share of a Paper 1 step shrinks and the **element loops become the floor**. On the
> rev3 1D deck at 2 ranks they are already **~85 % of the step** on the critical rank, and
> **> 99 % of that is `ASDShellQ4`** (layered sections, ASDConcrete3D IMPL-EX). WP-107 could
> not help: one loop, all-or-nothing allowlist with only `LadrunoQuad` on it, MPI refused.
>
> WP-192 threads the **element state determination (`Domain::update`), the integrators'
> tangent and residual loops and the element commit** for the Paper 1 types, with a
> **per-type opt-in** (un-audited elements run serially *in their place*) and an explicit
> **hybrid MPI + threads** acknowledgement:
>
> ```tcl
> ladrunoThreads 4 -hybrid      ;# or LADRUNO_THREADS=4 LADRUNO_THREADS_HYBRID=1
> ```
>
> The assembly is a **chunked compute + ordered serial scatter**: element matrices are formed
> in parallel into per-chunk buffers and scattered by one thread in the serial order, so the
> assembled system — and every result — is **bit-identical to serial at every thread count**
> (measured on the rev3 1D and 3C windows at 2 / 4 / 8 threads: zero difference). On the
> profiled 1D window (2 ranks) the element loops run **1.7× / 2.7× / 3.4× faster at 2 / 4 / 8
> threads per rank** and the step **1.37× / 1.70× / 1.80×**. Default unchanged: 1 thread,
> byte-identical to the base binary, no measurable overhead.

## 1. Why — the measured split (baseline, before any change)

Base binary = `ladruno` 0c4c6a1d8 (+ this ADR stub), Windows/MSVC, local i7-13700KF, P-cores
only (ranks pinned), a shared box (another agent's 4-rank job ran throughout). Decks: the
apeGmsh port's rev3 emissions (`tools.campaign_deck emit 1D --rev rev3 --ranks 2`,
`cases.case_t3 --cell 3C --rev3 --ranks 2`) — the rev3 STKO decks are partitioned for 8
ranks (1D, 3C) and 32 (4D), too many for a local profile.

**1D, 2 ranks** (12 968 ASDShellQ4 + 192 forceBeamColumn; MUMPS, Transformation,
KrylovNewton, Newmark), full gravity, then 200 transient steps (quiet start, 1 iteration a
step), `profiler start -deep -perStep -memory`, ms/step:

| Phase | Rank 0 | Rank 1 (critical) |
|---|---|---|
| newStep (incl. `elem.update`) | 252 (197) | 443 (311) |
| formTangent (`elem.tangent`) | 611 (603) | 915 (903) |
| formUnbalance, 2/step (`elem.residual`) | 659 (656) | 1098 (1088) |
| update (`elem.update`) | 211 (209) | 330 (324) |
| commit (element commit + recorders + EnergyBalance) | 829 | 1189 |
| linearSolve (MUMPS factor + solve, incl. wait) | 1940 | 531 |
| **step** | **4505** | **4508** |

By type (rank 1, all element scopes): **ASDShellQ4 2 523 ms (99.8 %)** — 31–49 µs per update,
82–126 µs per tangent, 50–83 µs per residual call — ForceBeamColumn3d 6 ms (0.2 %).

**3C, 2 ranks** (33 200 stdBrick soil, 4 000 ASDAbsorbingBoundary3D, 6 620
ASDEmbeddedNodeElement, 3 512 ASDShellQ4, 192 forceBeamColumn), shortened gravity, 100
transient steps: step 29.8 s, **linearSolve 27.1–28.3 s (91–95 %)**, element loops ~2.6 s on
rank 1 (8.7 %). By type (rank 1, update + tangent + residual):

| Type | ms/step | share of element time |
|---|---|---|
| stdBrick | 1 340 (tangent 1 053: 63 µs per element per call) | 54 % |
| ASDShellQ4 | 517 | 21 % |
| ASDAbsorbingBoundary3D | 101 | 4 % |
| ForceBeamColumn3d | 31 | 1.2 % |
| ASDEmbeddedNodeElement | 22 | 0.9 % |
| commit (elements + recorders, not attributed) | 159 | 6 % |

So: 1D (and Tier 1–2) is element-bound already; Tier 3–4 are MUMPS-bound today (PROFILE_4D:
97.6 %) and become element-bound once reuse/METIS/fewer ranks shrink the factorization —
which is the case this WP prepares for.

## 2. What ships

| Piece | Where | Notes |
|---|---|---|
| `ladruno_parallelFor`, the thread-region flag `ladruno_tlInRegion`, `ladruno_loopThreads`, the hybrid acknowledgement, `LadrunoOutputGuard`, the `LADRUNO_THREADS_LOOPS` bisection mask | `SRC/utility/LadrunoThreads.{h,cpp}` | OPS_Utilities is compiled with OpenMP; OPS_Analysis is not (PFEM's dormant pragmas), so the integrators call the parallel-for through a function pointer |
| Per-type update loop | `Domain::ladrunoThreadedUpdate` | maximal runs of consecutive audited elements are threaded; every other element runs alone on the calling thread between them, so WP-107's "an un-audited element may write shared node state" stays fenced (it still sees exactly its serial state) |
| Commit loop | `Domain::ladrunoThreadedCommit` | same shape; opt-in `Element::ladrunoThreadSafeCommit` (default = the form answer) |
| Tangent / residual loops | `SRC/analysis/integrator/LadrunoThreadedAssembly.{h,cpp}`, hooked in `IncrementalIntegrator::formTangent`, `::formElementResidual`, `TransientIntegrator::formTangent` | §3 |
| Per-thread twins of the class-wide pools | `FE_Element` (tangent/residual pool → `ladrunoTan/ladrunoRes`), `TransformationFE` (`modMatrices`, `dataBuffer`, `localKbuffer`, `dofData`, `theTransformations`), `Element` (Rayleigh pool), `Matrix` (`matrixWork/intWork`) | ADR-75b 5.4-H1, "the 100 % collision no element allowlist can fix". Used only inside a region: every serial path is the vanilla one |
| Hybrid MPI + threads | `ladrunoThreads <n> -hybrid` (classic Tcl and interpreter), `LADRUNO_THREADS_HYBRID=1` | WP-107's MPI fence stays the default; an inherited `LADRUNO_THREADS` alone still cannot oversubscribe a job |
| Opt-in predicates | `Element::ladrunoThreadSafeForm`, `::ladrunoThreadSafeCommit`, `Material::ladrunoThreadSafeForm` (+ the existing `ladrunoThreadSafeUpdate`) | default false everywhere |

Refusals, in order: OpenMP compiled out → 1 thread → MPI binary without the acknowledgement →
PartitionedDomain/Subdomain → deep profiler armed (H6, unchanged) → fewer than 2 audited
elements. Every outcome is announced once per domain per outcome (WP-107 S3), e.g.
`element update loop PARTLY threaded on 4 threads: 46 of 47 elements threaded, 1 run SERIAL
in their place; the first of those, element 990000 (classTag 19), is not on the thread-safe
allowlist` and `element tangent loop threaded on 4 threads: 13160 of 13160 FE_Elements
computed in parallel, scattered in serial order`.

## 3. The assembly strategy: chunked compute, ordered serial scatter

ADR-75b §4 settled "freeze-sparsity + atomic scatter" for a fast mode and a row-wise
**ordered gather** for the reproducible mode, and priced the gather at Σ idSize² of live
storage (~1.5 GB at 1 M DOF, §4.2 (iii)). WP-192 takes neither — and neither graph colouring
nor per-thread copies of A:

1. The FE_Element list is snapshotted (H3) and walked in serial order.
2. A maximal run of consecutive thread-safe FEs, capped at `max(512, 64·threads)`, is
   computed in parallel: each FE's `getTangent()` / `getResidual()` is formed on a worker
   (into that worker's private pools) and **copied into its own slot** of a master-owned
   buffer.
3. One thread then calls `addA` / `addB` for the chunk **in serial order with the same
   matrices and IDs** the serial loop passes. A non-thread-safe FE (constraint FEs, contact,
   an un-audited element, a TransformationFE on a time-varying MP) is computed and scattered
   alone in its serial place.

Consequences:

* The SOE receives the identical sequence of `(FE index, i, j)` contributions (ADR-75b
  4.2 (i)), so **A and b are bit-identical to serial**, at every thread count, by
  construction. The DOF_Group loops keep their serial position (4.2 (ii)).
* Memory is bounded by the chunk, not the model: 512 × 24² × 8 B ≈ 2.4 MB for shells.
* The scatter stays serial (the Amdahl term visible in §6.3).
* Not colouring: it threads the scatter too, but changes the summation order (not
  bit-identical), needs a graph build per `domainChanged`, and buys only the scatter, which
  is not where the time is. Not per-thread SOE copies: memory ∝ nnz × threads (a MUMPS SOE
  of a 200 k-DOF rank is hundreds of MB) plus an FP reduction.

**Determinism: same thread count → same answer; every thread count → the serial answer, bit
for bit**, on a deterministic solver (serial Pardiso at `MKL_NUM_THREADS=1`, UmfPack,
BandGeneral). The update and commit loops have no reduction (WP-107's argument). Under
OpenSeesMP + MUMPS the *base* binary is itself not run-to-run reproducible (WP-188 quirk,
measured again in §6.1), so there threads add no scatter of their own.

## 4. The thread-safety audit

Rule set = WP-107 `Element.h` items 1–5, over the update, tangent/residual and commit call
graphs. Method: every function- and class-scope `static` scratch on those paths became
`thread_local` (zero bytes per element; identical arithmetic on one thread); every
process-wide accumulator was found and either made per-thread or left on the serial path.

| Type (deck role) | What was shared | Fix | update | form | commit |
|---|---|---|---|---|---|
| ASDShellQ4 (walls, slabs) | `ASDShellQ4Globals::instance()` (all work matrices, LHS/RHS), 20 function statics, both transformations, `ASDEICR` | `thread_local` singleton + statics | ✔ | ✔ | ✔ |
| LayeredShellFiberSection | class `stressResultant/tangent/array`, 5 statics | thread_local; answer = every layer | ✔ | ✔ | ✔ |
| PlateFromPlaneStress, PlaneStressMaterial, PlateFiberMaterial, PlateRebar, ElasticMembranePlateSection | class `stress/tangent(/array)`, the condensation work matrices | thread_local; wrappers delegate the answer | ✔ | ✔ | ✔ |
| ASDConcrete3D / ASDConcrete1D (IMPL-EX) | 30 / 10 statics; the `GlobalParameters` IMPL-EX error max/average accumulated in `commitState` | thread_local; the accumulator becomes per-thread partials folded (under a lock) before any read — the max is order-free, `avgImplexError` (a diagnostic no deck reads) may differ in its last bits | ✔ | ✔ | ✔ |
| Hysteretic, Elastic (uniaxial) | nothing on the state path | — | ✔ | ✔ | ✔ |
| ForceBeamColumn3d + FiberSection3d + LinearCrdTransf3d + Lobatto/Legendre | class `theMatrix/theVector/workArea/*Subdivide`, 56 statics; `Tlg/kg` + 57 transformation statics; the fiber location buffers | thread_local; failure-path `opserr` under `LadrunoOutputGuard`; opt-in only with LinearCrdTransf3d, Lobatto/Legendre, explicit fibers, no Damping | ✔ | ✔ | ✔ |
| stdBrick + ElasticIsotropic 3D | class `stiff/resid/mass/damping/xl`, file-scope `B`, 61 statics, `shp3d`; material `sigma/D` | thread_local; no Damping | ✔ | ✔ | ✔ |
| ASDEmbeddedNodeElement | 85 kernel statics | thread_local | ✔ | ✔ | ✔ |
| ASDAbsorbingBoundary3D | 68 kernel statics; commit books `ABSORB_LEAK` (process-wide FP sum) | thread_local; commit stays serial | ✔ | ✔ | serial |
| zeroLength (2C springs, dashpots), everything else | not audited | serial in place | serial | serial | serial |

Found on the way:

* **`Matrix::Solve(Matrix&, Matrix&)`**: the first cut of the per-thread work area left
  `x = b` inside the serial-only branch, so a threaded ForceBeamColumn3d inverted garbage
  (NaN on the second step). The threaded-vs-serial hash caught it; `LADRUNO_THREADS_LOOPS=1`
  isolated it to the update loop in one run — why the mask ships.
* **ASDAbsorbingBoundary3D latent defect**: one kernel sized its static `BB(6, m_num_dofs)`
  once, by the *first* caller, and never resized it (its three siblings resize per call).
  Fixed with `BB.resize`; serial behaviour changes only where it read past the end.
* **Matrix work area**: the triple products choose their code path by the work-area size and
  `Invert` passes it to LAPACK as `lwork`; a worker answers with the SHARED logical size
  (never grown in a region), so it takes the serial path bit for bit.

## 5. MPI: hybrid ranks × threads

OpenSeesMP gives every rank a plain `Domain`, so WP-107 had to fence the whole binary.
WP-192 keeps the fence as the default and lifts it on an explicit acknowledgement. Sizing
is the user's: ranks × threads ≤ cores, each rank pinned to a domain of cores (Intel MPI
`I_MPI_PIN_DOMAIN=[mask0,mask1]`, as below; Slurm `--cpus-per-task=t`, `OMP_PROC_BIND=close`).
MUMPS stays MPI-only (mkl_sequential, the ADR-75 threading-layer rule).

## 6. Verification

### 6.1 Byte-identity at the default (1 thread) and N threads = serial

| Configuration | Comparison | Result |
|---|---|---|
| 1D rev3 serial emission (`--ranks 1 --system Pardiso`, `MKL_NUM_THREADS=1`), 30 gravity increments + 5 transient steps | base vs WP-192 (`LADRUNO_THREADS` unset) | `disp` (every node, 6 dof, 17 digits, every step), `energy_balance`, `th_steps` byte-identical; the four `recorder ladruno` files (f64) identical in every float dataset |
| same | WP-192 at **2, 4, 8 threads** vs base | **identical** (max difference 0) |
| 3C serial emission (bricks, absorbing skin, embedded ties, shells, beams; 47 524 FEs), 30 gravity increments + 5 steps | base vs WP-192 unset; WP-192 at **4 threads** vs base | **identical** (all nodes, every step) |
| 1D 2 ranks, MUMPS | base vs base (two runs) | **differ already**: 9.9e-5 mm on 14.3 mm after 5 increments |
| same, 50 transient steps | base vs WP-192 at 1 thread (serial path) / at 2, 4, 8 | max |Δu| 1.6e-2 mm on 10.7 mm (1.5e-3) for the *unthreaded* run, same order for the threaded ones; WP-192 1 vs 4 threads happened to be byte-identical — the scatter is MUMPS's |

pytest `tests/test_wp192_threaded_element_loops.py` (zone_a): RC wall of layered ASDShellQ4
(ASDConcrete3D with and without IMPL-EX, PlateRebar/Hysteretic, a PlateFiber layer,
`-drillingNL`) + fiber forceBeamColumns (ASDConcrete1D) + rigid diaphragm (TransformationFE),
static push + Newmark/KrylovNewton shake; stdBrick block + ASDEmbeddedNode ties; a mixed
domain with an un-audited zeroLength; full-field md5 per step and the return codes equal at
1/2/4 threads; a repeated-run race probe (5 runs at the most threads = 1 thread).
`test_wp107_threaded_update.py` 18/18 unchanged. **Mutation gate**: (M1)
`ASDShellQ4Globals` back to `static` → 3 tests red; (M2) the FE_Element tangent pool shared
→ 4 red; (M3) `Brick::stiff/resid` shared → the brick test red.

### 6.2 Race testing

ThreadSanitizer is not available for the MSVC build and is not wired into Linux CI. The
race evidence is: per-step full-field hashes at 1/2/4 (pytest) and 2/4/8 (deck windows)
threads equal to serial, the repeated-run probe, the three mutations, and Zone-A (gcc,
Ubuntu) running the same tests green. A TSan job over these tests is a follow-up.

### 6.3 Speedup on the profiled window

1D, 2 ranks, `LADRUNO_THREADS=t LADRUNO_THREADS_HYBRID=1`, ranks pinned
(`I_MPI_PIN_DOMAIN`: t = 1 → cpu 0 | 8; 2 → {0,2} | {8,10}; 4 → {0,2,4,6} | {8,…,14}
(one thread per physical core); 8 → cores 0–3 | 4–7 with Hyper-Threading), 30 gravity
increments + 50 transient steps, coarse profiler (`-perStep`), rank 1 (critical), ms/step:

| threads/rank | step | ×  | element loops¹ | × | tangent | residual (2/step) | update | commit² | window wall |
|---|---|---|---|---|---|---|---|---|---|
| base binary | 3 768 | — | 2 296 | — | 771 | 913 | 268 | 1 030 | 348 s |
| 1 | 3 753 | 1.00 | 2 278 | 1.00 | 768 | 910 | 263 | 1 031 | 339 s |
| 2 | 2 736 | 1.37 | 1 335 | 1.71 | 441 | 498 | 159 | 961 | 250 s |
| 4 | 2 214 | 1.70 | 834 | 2.73 | 270 | 279 | 101 | 937 | 201 s |
| 8 (HT) | 2 080 | 1.80 | 673 | 3.38 | 212 | 206 | 86 | 951 | 186 s |

¹ newStep + formTangent + formUnbalance + update. ² Commit = element commit + node commit +
recorders + EnergyBalance; threading its element part saves ~140 ms at 4 threads
(`LADRUNO_THREADS_LOOPS=7`: 1 073 ms), so **~850 ms of it is not element work** — the four
`recorder ladruno` (f64 in these windows) and `EnergyBalance -v2` with 10 regions. At
8 threads that serial commit is 46 % of the step.

The 3C windows were not timed threaded: at 2 ranks 3C is 91–95 % MUMPS, so t threads can
only remove part of its 8.7 % element share today; the type shares above say what it buys
once the factorization shrinks.

## 7. Not done here (follow-ups, with their measured shares)

| Item | Share where measured | Note |
|---|---|---|
| Commit's non-element part (ladruno recorders, EnergyBalance -v2 regions) | 1D: ~850 ms of 3 753 (23 %) at 1 thread, 46 % of the step at 8 threads | the next floor; per-region energy sums are reductions (ordered-gather candidates) |
| stdBrick tangent recomputed every iteration for a linear soil | 3C: 54 % of element time (63 µs/element/tangent) | a kernel item, not threading: cache the elastic stiffness |
| zeroLength springs / dashpots (2C, 2 487) | not profiled here (2C) | audit ZeroLength + Elastic/Viscous uniaxials; WP-184 dashpot energy |
| Serial scatter (`addA`/`addB`) | part of the residual tangent/residual time at 8 threads (≈ 200 ms each) | coloured or row-owned scatter would break bit-identity; measure first |
| ForceBeamColumn3d with PDelta/Corotational, HingeRadau, a SectionIntegration | 0.2–1.2 % | opt-in currently refuses these configurations |
| ASDShellQ4 / stdBrick with a `Damping` object | — | opt-in refuses; audit `Damping` |
| TSan job (Linux) over the WP-107/192 tests | — | §6.2 |
| Cluster measurement (16×1 vs 16×2 vs 8×4, 2×8 layouts) | — | after the v4 pin; no esmeralda runs in this WP |

## 8. Caveats

* `opserr` from inside a threaded loop is serialised only where WP-192 put a guard
  (ForceBeamColumn3d failure messages); other failure-path prints of audited types
  (PlateFiberMaterial's non-convergence warning) may interleave. Failure paths only.
* The deep profiler still refuses threading (H6): per-type numbers come from serial runs;
  use `profiler start -perStep` to time a threaded run.
* Hyper-Threading: the 8-thread row shares 4 physical cores per rank; most of its gain over
  4 threads is in the tangent/residual loops.
