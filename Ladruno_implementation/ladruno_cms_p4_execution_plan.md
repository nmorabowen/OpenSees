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

Run all four, same Building-1A deck, same box:

- (a) standard solver, no CMS (ARPACK) — the bar;
- (b) CMS **level-2 only** ablation (omit T1) — benchmark-only, NEVER an accepted
  config or fallback (guard it so it cannot be selected as a solver);
- (c) full CMS L2+L1;
- (d) [after P3] CMS with the distributed global solve.

Report per-phase time, peak RSS/rank, comm volume, and n, r2, r*, n/r2, r2/r*
for each. This isolates the reduction each level contributes from the cost of the
final backend.

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
  sparse) without ever forming a dense (m×m) matrix.
- Assemble the reduced group/leader pencils as SPARSE (they are small but
  currently dense-materialized); feed MUMPS the sparse form.
- Stream the compatibility merges (`assembleCompatible`, `gatherCompatiblePencil`)
  instead of dense unique×unique buffers; replace the O(u²) linear-search key
  matching with a hash/sorted join.

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
