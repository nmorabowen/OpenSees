# LadrunoCMS P4 — performance regime + memory scalability: plan & interim verdict

**Status:** P4 NOT met. ADR §10 P4 ([[1000_ladruno_cms_adr]]) requires a
*documented* size/parallelism regime where CMS reduces **time OR memory** vs the
standard OpenSees solver, without violating P2, before `remediation → shipped`.
This document (a) gives the honest interim verdict from the existing one-run
baseline, (b) specifies the campaign to establish (or refute) a winning regime,
and (c) designs the memory-scalability work that is the actual gate. Runtime
items need a working `LADRUNO_CMS=ON` build + the Building-1A harness (see
[[ladruno_cms_p3_execution_plan]] prerequisites).

> **2026-07-26 — §2 is SETTLED and §30's scaling alarm is WITHDRAWN (ADR §32).**
> The 523x of §30 was a `k2` tuning artefact, not a property of the formulation.
> With `k2` chosen per size the CMS/standard ratio **flattens at ~4x** over a 16x
> range in `n`, and CMS fits `n^1.12` against the standard solver's `n^1.34` — it
> scales slightly *better*, not worse. The `O(n^4)` reading is retracted.
>
> **The gate is unchanged and still unmet on the time axis** (CMS never wins on
> speed; ~4x slower is the floor measured). What changes is that the memory-
> capacity rationale — which §30.3 had called into question — **is restored**: a
> constant ~4x is compatible with "the answer arrives", `O(n^4)` was not. So §4
> below (the workspace memory refactor) is again THE gate, exactly as originally
> planned. See ADR §32 for the ladder, the `k2` cliff and the scope limits.

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
>
> **SUPERSEDED 2026-07-26 (ADR §32) — the suspect was right and the caveat was
> understated: 523x *was* a tuning gap.** At n=12800 with k2=12 the local Lanczos
> thrashes (`restarts` 1 → 239, `opApplications` 120 → 22 032 going from m=760 to
> m=3120); raising k2 to 16 on the same binary takes that point from 105.6 s to
> 0.96 s. Final ladder, 5 repeats per arm, standard and CMS back-to-back, minimum
> quoted (the box is shared — see ADR §32.4 on why minimum):
>
> | mesh | n | standard | CMS k2=12 | CMS best k2 | ratio (best) |
> |---|---:|---:|---:|---:|---:|
> | 20x20 | 3 200 | 0.0288 s | 0.2328 s | 0.2328 s (k2=12) | 8.1x |
> | 40x40 | 12 800 | 0.2021 s | 105.63 s | 0.9592 s (k2=16) | **4.7x** |
> | 60x60 | 28 800 | 0.5914 s | >200 s | 2.3163 s (k2=16) | **3.9x** |
> | 80x80 | 51 200 | 1.1902 s | >300 s | 5.1306 s (k2=16) | **4.3x** |
>
> Fits: standard `n^1.34`, CMS-at-best-k2 `n^1.12`. **No automatic k2 heuristic
> was shipped** — the obvious `sqrt(m)` rule is measurably 1.4–8.7x off the
> optimum and its own ratio degrades (11.3x → 37.5x, `n^1.78`), and the optimum
> sits next to a two-order-of-magnitude cliff that is *non-monotone from below*
> (at n=3200: k2=8 fast, 9/10/11 slow, 12 fast again). The evidence-backed next
> step is a **diagnostic**, not a rule — see the new item 2b.

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
> **[[_adr1000_cms_p4_handoff]] is CONSUMED 2026-07-26 (ADR §32).** Both
> experiments ran; verdict **`tuning`, not algorithmic** — the 523x was `k2`, the
> `O(n^4)` reading is retracted, and CMS fits `n^1.12` against the standard
> solver's `n^1.34`. That does **not** reopen or close P4: the §31 deferral
> stands, for the reason §31 gives. What it removes is the argument that would
> have closed P4 without a campaign. **ADR §32.9 lists the three things this
> changes in the campaign design** — notably that campaign requirement 4 ("scale
> `k2` with the subdomain") is not sufficient as written, because the natural rule
> for satisfying it is the wrong one.

## 2b. `k2` convergence diagnostic (the actionable follow-up from ADR §32)

**Problem it solves.** A `k2` below the convergence threshold does not fail — it
*grinds*, silently, for minutes to forever, and returns the right answer if you
wait. The user has no way to tell "this model is big" from "this run is
thrashing". The two look identical from outside.

**What to do — a warning, not a heuristic.** The counters already exist
(`lanczosRestarts`, `lanczosOperatorApplications`, ADR §27) but are only printed
under `-verbose`. Emit a warning from the T2 path when `restarts` crosses a
threshold, naming the lever: *"LadrunoCMS: the fixed-interface basis is not
converging (restarts=N); raise -modesL2."* Cheap, no algorithm change, and it
turns a 300 s silent stall into an actionable message.

**Do this BEFORE the Esmeralda campaign, not after** (ADR §32.9 item 2). In a
cluster queue a `k2`-thrashing job looks exactly like "the model is too big": it
burns its walltime and gets killed. Without the warning the campaign can collect
false negatives and never know.

**What NOT to do yet.** Do not ship an automatic `k2` rule. ADR §32.7 gives the
three measured reasons: `sqrt(m)` is 1.4–8.7x off the optimum, the cliff is
non-monotone from below so "smallest converging `k2`" is unsafe, and four sizes
on one mesh topology cannot calibrate a parameter that decides between 1 s and
>300 s.

**The campaign that would justify a rule** (not started): sweep `k2` across ≥3
interface topologies and ≥2 requested-mode counts, so the dependence on `m` can
be separated from the dependence on `b` and on the number of modes requested.

## 3. Statistical discipline

Repeat each accepted point **≥3×** before any performance claim (RESULTS.md
decision 5). The current baseline is explicitly one run; treat its numbers as
indicative, not published.

> **Reinforced 2026-07-26 (ADR §32.4).** This box is shared and was measured at
> 76% CPU from unrelated processes. A 3-repeat 80x80 point came out at 8.80 s
> with ±11% spread; the identical configuration re-run on a quiet box gave
> 5.25 s at ±1%. Two consequences, now house rules for this lane:
> - **Report the minimum, not the mean.** Contention noise is one-sided, so the
>   minimum is the honest estimator of uncontended cost. Keep the mean and the
>   spread in the log.
> - **The spread is the alarm.** A tight ±1% band is evidence the run was clean;
>   a wide one means re-measure before believing the point. Run the two arms
>   back-to-back in one pass so both see the same machine state.

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

> **Status of the gate, 2026-07-26 (ADR §32).** Still **NOT met**, and the target
> is unchanged: **memory, not speed.** On the time axis CMS is ~4x slower than
> the standard solver at every measured size and there is no crossover in sight
> — that is now a *measured floor* rather than a suspicion, and it is not going
> to be closed by tuning.
>
> The change from §30 is that the gate is **reachable again**. §30 read CMS as
> `O(n^4)`, which would have killed the memory-capacity regime outright (a model
> you can fit but cannot finish is not a win). §32 retracts that: cost is a
> constant ~4x with `n^1.12` scaling, so a capacity win would be a real win.
>
> Order of work from here: (1) the §2b `k2` diagnostic — cheap, and it stops the
> next person losing a session to a silent stall; (2) §4, the workspace memory
> refactor, which is the gate itself; (3) Building 1A, still blocked on the deck
> being absent from the repo, which remains the only measurement that can close
> P4 for real. Do **not** spend further sessions optimising CMS wall clock: §26,
> §29 and §32 together say the speed axis is not where this is decided.
