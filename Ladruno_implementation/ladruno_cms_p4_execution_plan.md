# LadrunoCMS P4 — performance regime + memory scalability: plan & interim verdict

**Status:** P4 NOT met. ADR §10 P4 ([[1000_ladruno_cms_adr]]) requires a
*documented* size/parallelism regime where CMS reduces **time OR memory** vs the
standard OpenSees solver, without violating P2, before `remediation → shipped`.
This document (a) gives the honest interim verdict from the existing one-run
baseline, (b) specifies the campaign to establish (or refute) a winning regime,
and (c) designs the memory-scalability work that is the actual gate. Runtime
items need a working `LADRUNO_CMS=ON` build + the Building-1A harness (see
[[ladruno_cms_p3_execution_plan]] prerequisites).

## Interim verdict (from `scalability/RESULTS.md`, ONE run — not yet statistical)

| Method | Ranks | Solve [s] | Peak RSS | vs CMS |
|---|---:|---:|---|---|
| ARPACK (standard) | 1 | 30.1 | 2.0 GiB | baseline |
| CMS physical | 4 | 311.6 | 10.3 GiB (rank 4.0 GiB) | 7.0x slower, 5x more mem |
| CMS physical | 6 | 208.7 | 6.7 GiB | 6.9x slower |
| FEAST | 4 | 392.0 | 4.5 GiB | CMS 1.26x faster |
| FEAST | 6 | 448.8 | 6.6 GiB | CMS 2.15x faster |

**Honest reading:** at Building-1A scale CMS does **NOT** beat the standard
sequential solver on either axis — it is ~7x slower AND uses more memory. Its
only win is against FEAST (another *parallel* eigensolver), which is not the
ADR's bar ("the standard OpenSees solver"). **So the winning regime the ADR
requires for `shipped` does not currently exist.** Do not claim CMS is faster
than ARPACK on Building 1A (RESULTS.md decision 2).

**Where a real win is plausible — and it is MEMORY, not speed.** CMS cannot
out-run a direct sparse ARPACK at these sizes. Its structural advantage is that
no rank holds the whole problem, so at large `n` (a model that does NOT fit in
one node's RAM) CMS could be the only solver that runs at all — a **memory-
capacity** win, which §10 explicitly allows ("reduces ... memory"). BUT today
memory does NOT scale: the local Craig-Bampton construction materializes DENSE
rank-local workspaces (largest Building-1A rank 4.0 GiB). Until that is fixed the
memory-capacity regime cannot even be demonstrated. **The memory refactor (§4
below) is THE gate to any P4 win.**

## 1. Instrumentation — mostly present; close these gaps

> **DONE 2026-07-26 (ADR §27).** Finer hierarchy phases, per-rank peak RSS and
> the n/r2, r2/r* ratios are wired and reported under `-verbose`; the first
> profile is recorded in §27. Comm volume is still **not** instrumented. The
> profile's headline: T2 is 97% of the hierarchy, refinement is 58% of total,
> and r2/r* = 1.06 (level 1 buys almost nothing on a 1-D chain). It also
> surfaced a hard-coded `maximumRestarts = 20` with no user control, now exposed
> as `-maxRestarts` (§28).

Already reported (`LadrunoCMSEigenSolver.cpp`, `-verbose`): assembly / hierarchy
/ refinement seconds; refinement sub-phases (factorization, solve,
inverseRefinement, massAction, RayleighRitz); dimensions rRaw, rD (=r*),
originalZ. Gaps vs the §10 P4 list:

- **Finer hierarchy phases** — split `hierarchySeconds` into partition, fine
  modes (T2), compatibility (S2), level-1 (T1), global solve (S1),
  back-substitution, publication (the ADR names these individually).
- **Peak RSS per rank** — add a platform query (getrusage `ru_maxrss` on Linux;
  the CI lane is Linux via mpiexec) reported per rank at the end.
- **Communication volume** — count bytes moved by the Allreduce/Bcast/Allgatherv
  calls (a thread-through counter in the refiner + hierarchy merges).
- **Dimension ratios** — the solver already has n (globalDimension), r2 (dim
  after S2) and r* (rD after S1); just emit `n/r2` and `r2/r*` (r2/r* measures
  the marginal contribution of level 1). Trivial arithmetic on existing values.

These are additive reporting only — no algorithm change — and belong in the same
`-verbose` block. Compile-check with the nightly lane.

## 2. The ablation comparison (identical hardware + model)

> **RE-MEASURED 2026-07-26 (ADR §30), and the news is bad.** Not on Building 1A
> — that deck is still missing — but on a like-for-like sheet model where both
> paths solve the identical mesh (`cms_compare_arpack.tcl`; the six eigenvalues
> agree to ~1e-12). Standard = `eigen` with Arpack+UmfPack on 1 rank; CMS = 4
> ranks.
>
> | mesh | n | standard | CMS | ratio |
> |---|---:|---:|---:|---:|
> | 20x20 | 3 200 | 0.0294 s | 0.2335 s | **7.9x** |
> | 40x40 | 12 800 | 0.2061 s | **107.90 s** | **523x** |
> | 60x60 | 28 800 | 0.6480 s | did not finish | — |
>
> The 7.9x at n=3200 reproduces the Building-1A 7.0x almost exactly (§29's work
> moved it down from ~11x). **The second row is the problem:** 4x the DOFs costs
> the standard solver 7x and CMS **462x**. CMS scales far WORSE, which undercuts
> the memory-capacity argument — a capacity win is worthless if the answer never
> arrives. Caveat: CMS parameters were held fixed (k2=12, k1=24) at every size,
> so a campaign that scales k2 with n would do better; 523x is not a tuning gap
> though. Prime suspect: the fixed-interface Lanczos cost growing
> super-linearly in m at constant k2 (consistent with the restart exhaustion in
> §27.4).

Run all four, same Building-1A deck, same box:

- (a) standard solver, no CMS (ARPACK) — the bar;
- (b) CMS **level-2 only** ablation (omit T1) — benchmark-only, NEVER an accepted
  config or fallback (guard it so it cannot be selected as a solver);
- (c) full CMS L2+L1;
- (d) [after P3] CMS with the distributed global solve.

Report per-phase time, peak RSS/rank, comm volume, and n, r2, r*, n/r2, r2/r*
for each. This isolates the reduction each level contributes from the cost of the
final backend.

> **VERDICT DEFERRED TO AN ESMERALDA CAMPAIGN (ADR §31).** All the negative
> evidence was taken at 4 ranks on a desktop — precisely where the source paper's
> own data predicts CMS cannot win (Yu et al. measure 6.2M–13.2M DOF on
> 4 096–65 536 cores, against another *parallel* CMS, not against a sequential
> direct solver; at their smallest configuration the advantage is 1.11x). Closing
> the lane without testing it in the regime it was designed for would mirror the
> original mistake of building it without checking that regime was reachable.
> §31 defines the campaign and its exit criterion, including the competitor that
> was in the tree all along: `eigen` over `MumpsParallelSOE`.
>
> **NEXT SESSION: [[_adr1000_cms_p4_handoff]]** — the two experiments that
> decide whether the scaling collapse is a `k2` tuning artefact or an
> algorithmic property, with acceptance criteria for each outcome (including
> closing P4 honestly if it is algorithmic).

## 3. Statistical discipline

Repeat each accepted point **≥3×** before any performance claim (RESULTS.md
decision 5). The current baseline is explicitly one run; treat its numbers as
indicative, not published.

## 4. Memory scalability — the gate (design, not yet implemented)

The debt: `reduceCraigBampton` and the S2/S1 merges build DENSE rank-local
matrices (unique×unique group pencils, dense-as-CSR fed to MUMPS, the local CB
`phi`/`psi` blocks). For a partition with `m` interior + `b` interface DOFs this
is O((m+b)²) dense per rank — the multi-GB workspaces.

**Design direction (verify each step preserves P2c to 1e-8):**
- Keep the local interior pencils SPARSE end-to-end; drive the fixed-interface
  modes through the existing MUMPS inverse-action + block Lanczos (already
  sparse) without ever forming a dense (m×m) matrix. **NOT DONE** —
  `reduceCraigBampton` still builds dense `phi`/`psi`.
- Assemble the reduced group/leader pencils as SPARSE (they are small but
  currently dense-materialized); feed MUMPS the sparse form. **HALF DONE
  2026-07-26** — `assembleCompatible` is sparse now; `gatherCompatiblePencil`
  still materialises a dense block per participant when unpacking the MPI upper
  triangle.
- Stream the compatibility merges (`assembleCompatible`, `gatherCompatiblePencil`)
  instead of dense unique×unique buffers; replace the O(u²) linear-search key
  matching with a hash/sorted join. **DONE 2026-07-26** for `assembleCompatible`
  and the hash join.

> **Status 2026-07-26 (ADR §26).** `assembleCompatible` no longer materialises a
> `dimension × dimension` dense buffer, and no longer emits a structurally dense
> CSR — the old path kept every upper entry including zeros, which is also what
> handed MUMPS the dense pattern behind the Part-0 ordering failure. Cost goes
> from O(dim²) to O(nnz): measured on an F-subdomain chain the saving grows
> 3.4x → 43x as F goes 4 → 64, an asymptotic improvement rather than a constant.
> Equivalence is exact, not merely in-tolerance — the physical smoke returns
> `maxResidual = 5.80175e-09` with identical dimensions, the same value as before
> the refactor. **Not measured on Building 1A** (the deck is not in the repo).

**The irreducible floor (ADR §13) — state it honestly, do not fight it here.**
Under the current `getEigenvector`/`AnalysisModel` contract each rank stores at
least `8·n·numModes` bytes of *replicated full* eigenvectors (n=500k, 50 modes ≈
190 MiB/rank), plus Domain structures. v1 uses per-mode `MPI_Allgatherv` to bound
the transient buffer but cannot remove the persistent replica. Publishing only
owned DOFs (true distributed output) requires changing `AnalysisModel` and every
eigenvector consumer — **a SEPARATE ADR; do not silently change AnalysisModel.**
So the achievable P4 memory win is in the *workspace* term (the dense CB/merge
temporaries), not the eigenvector-replica floor.

## P4 exit gate

`shipped` requires ONE documented (n, ranks) regime where full CMS reduces time
OR peak-RSS-per-rank vs the standard solver, ≥3 repeats, P2c intact. Given the
speed data, the realistic target is a **memory-capacity** regime at large n after
the §4 workspace refactor: a model that exceeds single-node RAM which CMS solves
and ARPACK cannot. If, after the refactor, no such regime exists at reachable
scales, the honest outcome is to keep CMS at `remediation` as a verified-but-not-
promoted research solver — the Shenwei speedups are not a local acceptance
criterion (§10).
