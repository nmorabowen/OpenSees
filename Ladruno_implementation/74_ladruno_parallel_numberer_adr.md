---
title: ParallelNumberer O(V²) merge — the 16.35 h first-step stall (root cause + fix ladder)
project: Ladruno
status: proposed — root cause CONFIRMED (G0 measured dc.numberDOF ~ N^1.98 + closing arithmetic, two independent derivations); 3-lens adversarial review folded; T0/T1 ready to implement
priority: high
owner: nmora
amends: 40_ladruno_performance_adr
supersedes: 30_ladruno_parallel_numberer_adr
tags:
  - adr
  - performance
  - parallel
  - numberer
  - mpi
  - setup
  - explicit
  - sub-adr
---

# ADR-74 — ParallelNumberer O(V²) merge: the 16.35 h first-step stall

> ADR-74. The **setup-lane** perf sub-ADR that [[40_ladruno_performance_adr]] P8 calls for,
> spawned not by a Phase-0 profiler number but by a **production incident**: slurm 144090
> (18.56 M hex / 19.18 M node / 57.5 M DOF, np240 explicit) spent **58,876 s (16.35 h) on a
> single analysis step** while the other 999 steps averaged 1.42 s. ADR-40 listed "the
> ParMETIS/HPC stack ahead of a *measured production-scale setup bottleneck*" as an
> anti-goal — **that gate is now cleared by measurement** for the numbering path (and only
> that path). **Supersedes [[30_ladruno_parallel_numberer_adr]]** (see §Supersedes).
> Family: ADR-40 (perf program) · ADR-40a (instrument scopes) · ADR-67 (integrator lane) ·
> ADR-68 (element/state-determination lane — the *other* per-step axis; this ADR is the
> one-time-setup axis) · ADR-06 (profiler — the instrument that both missed and then caught
> this, see §Instrumentation) · ADR-51/55/73 (domain-change multipliers, see §Sequencing).
> *Family-link note:* ADR-60–73 landed on `ladruno` with the post-#588 sync (this ADR was
> numbered against that tree, so 74 = next free); all family links now resolve locally.
> External record: Obsidian `Cluster Scalability/Max Model — 18.6M P6 Cluster-Emit
> (20x20km).md` §8 (the incident), `Max Model — 14M-Hex Hero Run.md` §7 (retro-explained
> below). Adversarial reviews: Fable cross-review 2026-07-21 (pre-ADR) + 3-lens gate
> 2026-07-21 (§Adversarial review log — one claim refuted and corrected, two weakened and
> reframed, root cause **stands**).

## What

`ParallelNumberer::numberDOF()` — the engine under `numberer ParallelRCM` /
`ParallelPlain` in OpenSeesMP — gathers every rank's DOF-group graph onto rank 0 and
merges them with per-vertex **linear scans** (`ID::getLocation`), making first-step setup
**O(V²_global) serial on one core** while every other rank busy-polls in `recvID`. At
19.18 M nodes that is ~2×10¹⁴ integer comparisons ≈ 16-19 h — the observed stall, to
order of magnitude, with the same comparison model reproducing the np8 G0 rungs to ≤3%
(§Gates). **G0 measured the isolated `dc.numberDOF` bracket at N^1.98.**

**Delivery vehicle (owner decision 2026-07-21): a fork-native class, not an in-place
edit.** `LadrunoParallelNumberer : public ParallelNumberer` (new files, ladruno numberer
classTag in the ≥33000 band), following the fork's additive discipline: upstream
`ParallelNumberer` stays behavior-untouched; the only upstream-file diffs are (a) a
**ledgered private→protected promotion** in `ParallelNumberer.h` (`processID` /
`numChannels` / `theChannels` / `theNumberer` — same pattern as the HHT.h and
GeneralizedAlpha.h promotions, classTags.h 33013/33014), and (b) additive registration.
`mergeSubGraph` is already `protected`. The tiers:

- **T0 — kill the quadratic (the stall fix):** override `numberDOF(int)` with a dense
  direct-address lookup (node tags are 1..N from gmsh/apeGmsh) + hash-map fallback,
  applied to **both** merge passes *and* the plain-ordering branch. O(V²)→O(V).
  Ordering-stable by construction ⇒ **bit-identical equation numbering vs stock
  `ParallelRCM`** (gate G1). The plain branch is a deliberate exception: fixing its
  latent tag-0 quirk changes the (still-valid) permutation, so the plain verb is gated by
  bijection + same-physics (G1b), not bit-identity.
- **T1 — kill the pointer-chasing residual:** after T0 the wall is `std::map` finds —
  2 per directed adjacency entry in `Graph::addEdge` **and ~1 per entry in RCM's BFS
  reads** (~1.6×10⁹ finds ≈ 25-45 min at 19 M). **Redesign per adversarial review** (the
  originally-specified `startAddEdge`/`addEdgeFast` borrow fails twice: its dense
  snapshot goes stale every subgraph because the merge interleaves vertex-adds and
  edge-adds, costing an O(V·P) rebuild; and it never touches RCM's *reads*): the subclass
  **owns a dense `tag → Vertex*` vector end-to-end** through the merge and calls
  `Vertex::addEdge` directly (`ID::insert` is sorted ⇒ adjacency order canonical ⇒
  bit-identity safe), and builds the merged graph on **dense array-backed storage**
  (merged tags are `getFreeTag`-contiguous, so array storage is exact) so RCM's BFS reads
  are O(1) too. Numbering drops to low minutes at 19 M.
- **T2 — stop centralizing (the architecture):** owner-computes distributed numbering
  (owner rule on shared nodes, per-rank owned-DOF counts, `MPI_Exscan` start offsets, one
  halo exchange) — imported from superseded ADR-30 (§Supersedes) with its open questions.
  O(V/P) time and memory, no rank-0 merged graph (~5 GB at this scale), milliseconds at
  any np. Verified compatible with `MPIDiagonalSOE::setSize`, which discovers neighbors
  by intersecting *global equation numbers* (`MPIDiagonalSOE.cpp:271`) — any consistent
  global bijection works. T2 changes the numbering (it abandons global RCM), so it is
  **opt-in per analysis**, not a silent default — see §Fence.

## Why — the incident, and the arithmetic that closes

**Observed** (slurm 144090, np240 / 15 nodes, `-deep -perStep -memory` profiler): one step
at 58,876 s wall, synchronized to the second across all 240 ranks; CPU/wall = 0.999
(MPI busy-poll); the timed sub-phases (`newStep`+`solveCurrentStep`+`commit`) covered only
2.6% of it — the cost sat in the *untimed* region of the step loop, which is exactly where
`domainChanged()` → `numberDOF()` runs on step 1.

**Predicted from source** (owner derivation + refutation-lens re-derivation, agreeing).
Pass 1: `mergeSubGraph` dedups each incoming vertex against the global ref list with
`vertexRefs.getLocation(ref)` (`ParallelNumberer.cpp:381`); `ID::getLocation` is a linear
scan (`ID.cpp:219-229`). ~90% of a 3D partition's vertices are interior — *not yet
present* — and a linear scan proves absence only by walking the **full** current array.
With 239 remote subgraphs of ~80 k vertices merged sequentially:
Σⱼ j·V²_sub ≈ **1.83×10¹⁴ comparisons**. Pass 2 (adjacency, `:406`/`:413`) adds
~2×10¹³ against the L2-resident `theSubdomainMap`, plus ~25-45 min of `std::map` graph
operations. An early-exit compare loop does not auto-vectorize and retires ~1
element/cycle ≈ 3×10⁹ cmp/s (compute-limited; DRAM bandwidth is near but not binding for
a sequential 4 B-stride scan). Honest total: **~2.05×10¹⁴ cmp ⇒ ~68-72 k s predicted vs
58,876 s observed** — an order-of-magnitude closure over-predicting by 15-25% (boundary
duplication and turbo effects absorb part of it), with **no competing mechanism within
10×**: the fabric gather is bulk-bounded to ~10-30 s (§Where m-wire), map operations cap
at ~45 min, `MPIDiagonalSOE::setSize`'s per-pair intersection is a linear two-pointer
merge. *(An earlier draft claimed "±5%"; the adversarial review showed that figure
dropped the pass-2 term — the corrected framing above is weaker in precision and
stronger in honesty; the decisive evidence is G0.)*

**The strongest evidence — the model reproduces the measured rungs.** The same
comparison model, evaluated at np8, predicts 47.6 s and 711 s for the G0 rungs; measured
first-step walls are 48.95 s and 693.71 s (≤3%), in a regime where **pass 2 dominates**
(pass 2 scales ∝ V²/P — large at np8, negligible at np240) — confirming both quadratic
passes at once. And the isolated `dc.numberDOF` bracket grew ×14.5 for ×3.87 nodes:
**N^1.98** (§Gates G0).

**Why it hid for ~20 years** (the code is stock 2007 upstream): O(V²) with a tiny
constant — ~3 s at 100 k nodes, ~5 min at 1 M — absorbed into "setup" at every scale
anyone ran before the apeGmsh cluster-emit pipeline (ADR 0065) made ≥10 M-node MP decks
routine. And the profiler could not see it: the cost lands *between* the timed brackets,
and the `-perStep` wall column wrote zeros (both fixed, §Instrumentation).

**The 14 M run retro-explained.** The hero run (slurm 144087, 11.09 M nodes, np240,
`ANALYZE_MS` = 20,180 s) reported a puzzling 449 µs/hex·step ≈ 25× the light-rank floor.
If its steady state actually sat at the ~18 µs/hex·step floor (0.82 s/step × 999 ≈ 822 s),
the implied hidden first-step cost is 19,358 s; the quadratic scaled from the 18.6 M stall
predicts 58,876 × (11.09/19.18)² ≈ 19,684 s — **1.7% apart**. The "449 µs" was an average
inflated by the same invisible stall. (Caveat: element changed between runs — stdBrick vs
`LadrunoBrick -lumped` — so the floor assumption is soft; even at 2× the floor the
retro-fit stays consistent. Corroboration, not proof; G0/G3 are the proof.)

**The irony that sharpens the fix:** this model runs `system MPIDiagonal` under explicit
central difference. The LHS is diagonal — no factorization, no fill — so the **RCM
ordering the entire gather-merge exists to feed buys nothing** on this lane. For the
explicit lane the Ladruno **plain** verb (merge + subdomain-order scatter, no RCM pass at
all) is the fastest pre-T2 option; the T2 endgame is to not centralize at all.

## Supersedes — ADR-30, and what is imported from it

[[30_ladruno_parallel_numberer_adr]] ("LadrunoParallelNumberer — distributed DOF
numbering", draft, 2026-06-12) already identified the P0 gather/merge serial term and
sketched what this ADR calls T2 — **under the same class name and file path**. It is
superseded, not merely amended, because its central quantitative premise is now measured
wrong: it modeled setup as "roughly linear in total DOF" and projected **~3 min at 18 M
DOF**; the measured reality is **quadratic and 16.35 h** — a ~300× under-projection,
because the O(V²) merge scan was not visible to a linear model. Concretely:

- **The class name and path belong to this ADR** (T0/T1 subclass of `ParallelNumberer`).
  ADR-30's distributed design becomes **T2** of this ladder; its open questions Q1-Q4
  (owner rule, halo protocol, constraint interaction, staged renumber) carry over as T2's
  design backlog, unresolved.
- **Imported lessons:** ADR-30's staged gate ("renumbers correctly every stage") becomes
  part of G4; its fork-PR #233 registration lesson (classic-Tcl must extend the flat
  factory-table chain, NOT add a nested else-if — MSVC C1061) governs the registration
  work here.
- **Correction propagated:** ADR-30's "roughly linear" claim and 3-min projection are
  struck; ADR-40's references to ADR-30 (P7/P9) should be read as pointing here.
- *Housekeeping:* two `30_*` files exist (`30_ladruno_explicit_constraint_projection_adr`
  owns the `tests/test_adr30_*` battery) — a pre-existing numbering collision this ADR
  sidesteps by taking 74.

## Where

| what | file:line |
|---|---|
| the linear scan | `SRC/matrix/ID.cpp:219-229` (`getLocation`); binary-search sibling `getLocationOrdered` at `:232-253` (unusable here — see Decision log) |
| vertex-dedup pass (pass 1) | `SRC/analysis/numberer/ParallelNumberer.cpp:381` (`vertexRefs.getLocation`) |
| adjacency pass (pass 2) | `ParallelNumberer.cpp:406` (1/vertex) + `:413` (1/directed adjacency entry) vs `theSubdomainMap` — ~2×10¹³ comparisons at np240, **dominant at np8** (∝ V²/P) |
| `ParallelPlain` extra quadratic | `ParallelNumberer.cpp:222` (`ID theOrderedRefs(V)` ⇒ sz=V from call one) + `:241`,`:249` ⇒ V² ≈ 3.7×10¹⁴ **additional**; tag-0 quirk is *accidentally benign* (one valid-permutation shift — vertex 0 numbered last), not corruption |
| the **second** quadratic (survives T0 unless indexed) | `ParallelNumberer.cpp:318-350` — per DOF group holding a `-4`, a full scan of **all** MP_Constraints: O(#constrained-groups × #MPs). Invisible in the MP-free plane-wave gates; bites equalDOF/tie-heavy decks (interop). T0 folds in a one-pass `constrainedNode → MPs` index (cheap, O(#MPs)) |
| broken-upstream variant (scoped OUT) | `ParallelNumberer.cpp:507-599` — `numberDOF(ID&)`: numbering call commented out (`:582`), recv loop strides `i += 2` against a differently-packed buffer (`:532`), start-DOFs never filled. Reachable only from the SP/DomainDecomposition lane (`DomainDecompositionAnalysis.cpp:229/247`). The Ladruno override **hard-errors** with a clear message; no equivalence obligation to code that never worked |
| ID append growth (verified benign) | `ID.cpp:335` (doubling ⇒ amortized O(1); the appends are NOT the problem) |
| post-T0 residual | `Graph.cpp:189-190` (`addEdge`: 2 map finds per directed entry) + RCM BFS `getVertexPtr` reads (`RCM.cpp:145,155,221,231` → `Graph.cpp:279-285`) ⇒ ~1.6×10⁹ finds ≈ 25-45 min — the T1 target, both halves |
| T1 building blocks | `Vertex::addEdge` (sorted `ID::insert` ⇒ canonical adjacency order); `getFreeTag` contiguity (`Graph.cpp:310`); dense array-backed graph storage for O(1) BFS reads. (`Graph::startAddEdge`/`addEdgeFast` exist at `Graph.cpp:223-276` but are NOT used — snapshot-staleness, see Decision log) |
| wire volume (verified-negative for the fabric hypothesis) | `Graph::sendSelf` ships (5V+2E) ints + V doubles ⇒ ~2.6-2.8 GB gathered at 19.18 M ≈ **10-30 s** into P0 over 2.5 GbE; ~6 GB / <1 min at 43 M. Bulk-bounded — hours are impossible here |
| T2 compatibility | `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp:271` (neighbor discovery by global-eq-number intersection — numbering-scheme-agnostic) |
| the setSize quadratic (RESOLVED, §Resolved) | `MPIDiagonalSOE.cpp` `quickSort`/`q_sort` (leftmost-pivot) fed ascending input by `AnalysisModel::getDOFGraph`'s `MapOfTaggedObjects` iteration ⇒ O((N/P)²/2) per rank; fixed by `std::sort` on the raw array |
| registration — the verbs (MP) | `SRC/tcl/commands.cpp:3919-3929` (`_PARALLEL_INTERPRETERS`: construct + `setProcessID(OPS_rank)` + `setChannels` at parse — a new verb slots in cleanly); `SRC/interpreter/OpenSeesCommands.cpp:1575-1581`, `:4146-4168` (needs a new `OPS_LadrunoParallelRCM` factory + `OpenSeesCommands.h` declaration) |
| registration — full surface (per the HHT ledger precedent, LEDGER_vanilla_files.md) | classTags.h · commands.cpp · OpenSeesCommands.{h,cpp} · `FEM_ObjectBrokerAllClasses.cpp` (the real numberer factory — `FEM_ObjectBroker::getNewNumberer` is a stub returning 0) · TclPackageClassBroker.cpp · CMakeLists + Makefile · banner regen (`patch_banner.py`) |

## How — the class, and the T0/T1 design

**`LadrunoParallelNumberer`** (`SRC/analysis/numberer/LadrunoParallelNumberer.{h,cpp}`):
subclass of `ParallelNumberer`; overrides `numberDOF(int)`; overrides `numberDOF(ID&)`
to **error** (broken upstream, §Where); inherits the channel/processID plumbing via the
ledgered promotion (which must include **`numChannels`** — the rank-0 loop needs it);
keeps RCM (or none) as the graph numberer on top — same algorithm, sane data structures.
Two verbs: `numberer LadrunoParallelRCM` (bit-identical to stock `ParallelRCM`, G1) and
`numberer LadrunoParallelPlain` (bijection + same-physics, G1b; help text advertises the
tag-0 fix — for the explicit/diagonal lane this is the recommended verb: no RCM pass).
apeGmsh side: a frozen dataclass in `opensees/analysis/numberer.py` + `_NumbererNS`
method + export (three small edits mirroring `ParallelRCM`); stock primitive untouched;
fork-gated with the emit-time "non-fork MP build" warning pattern.

**T0 — the lookup (hoisted to `numberDOF`, threaded through all merge calls, NOT rebuilt
per subgraph):**

- **Primary:** `std::vector<int32_t> refToLoc` sized maxTag+1, sentinel −1 — O(1),
  76 MB at 19.18 M, perfectly prefetchable. **Size math in 64-bit** (`maxTag` near 2³¹
  must not overflow). Guard: if maxTag > κ·numVertex (κ≈4; sparse/hand-authored tag
  spaces), fall back to `std::unordered_map<int,int>` with `reserve(V)`.
- **`ref < 0` hard-guard** (adversarial find): `Vertex::getRef()` is
  `DOF_Group::getNodeTag()`, which returns **−1 for node-less DOF groups** (Lagrange
  handler multipliers). Stock code *silently fuses* all such vertices across ranks into
  one merged vertex — a pre-existing wrong-numbering landmine (→ LEDGER_quirks); a dense
  table would turn it into a negative-index UB. The Ladruno class **errors loudly** on
  `ref < 0` with a message naming the Lagrange-handler-under-MP gap.
- **Pass 1** (`:381`): table probe; on miss, append exactly as today (same `getFreeTag()`
  sequence, same append order) and record `refToLoc[ref] = loc`.
- **Pass 2** (`:406`,`:413`): per-subgraph `subTagToMerged` map filled during pass 1;
  pass 2 becomes O(E_sub).
- **Plain branch** (`:241`,`:249`): seen-bitmap keyed by merged tag; the tag-0 quirk is
  fixed (behavior change, scoped to the plain verb, G1b-gated).
- **`-4` fixup** (`:318-350`): one pass over MPs building `constrainedNode → MPs` before
  the loop — kills the second quadratic for tie-heavy decks.
- **Do not mirror the 3-arg ctor** (adversarial find): it leaves `theNumberer`
  uninitialized while the dtor `delete`s it — UB grenade with no current callers. The
  subclass exposes only the `GraphNumberer&` and default ctors, and never deletes
  `theNumberer` / never touches `theChannels` (base-owned, two wiring regimes).

**T1 — owned dense structures end-to-end:** the subclass maintains its own
`tag → Vertex*` vector across the whole merge (updated as vertices are created — no
snapshot to go stale), calls `Vertex::addEdge` on both endpoints directly (sorted-set
insert ⇒ adjacency order canonical regardless of insertion order ⇒ RCM input identical),
and hosts the merged graph on dense array-backed storage so RCM's BFS `getVertexPtr`
reads are O(1). This replaces the earlier `startAddEdge`/`addEdgeFast` plan, which the
adversarial review showed cannot deliver: the dense snapshot staleness costs O(V·P)
rebuilds (invisible at np8, minutes at np240 — a G2-passes/G3-stalls trap), and edge-add
acceleration never reaches RCM's *reads*.

**Ordering stability is the load-bearing correctness claim:** the structures above are
pure lookup accelerators — iteration order over subgraph vertices, `getFreeTag` tag
assignment, and append order are untouched ⇒ identical merged graph ⇒ identical RCM input
⇒ **bit-identical equation numbering** (RCM verb). G1 enforces this, not the argument.

## Performance ceiling & sequencing — why T0+T1 first, T2 deferred

T0+T1 is the right *first* fix, **not** the performance optimum. The optimum for the lane
that hit the wall is T2: under `MPIDiagonal` + explicit there is no factorization, so the
ideal numberer does no gather, no merge, no ordering at all. Ranked, with projections:

| design | 19.2 M nodes | 43 M (the P=8 retry) | numbering | correctness risk |
|---|---|---|---|---|
| today | 16.35 h | ~3.4 **days** (quadratic) | RCM | — |
| T0 only | ~25-45 min (map residual) | ~1-2 h | bit-identical | ~none (G1 diff) |
| **T0+T1** | **~1-3 min** | **~5-10 min** | bit-identical | ~none (G1 diff) |
| T0+T1, plain verb | **~1 min** | **~2-5 min** | valid, no RCM | G1b (bijection + physics) |
| T2 | ~ms-s | ~ms-s | different (valid) | own correctness campaign |

*Measured update (N2): T0 alone already lands ~16 s at 2.1 M np8 (N^1.10) — the
"T0 only" row above was conservative; at 19.18 M the numberer term extrapolates to
~2-4 min pre-T1. The table's remaining binder is `dc.setSize` (§Risks).*

*Measured update (N3): T1 lands `dc.numberDOF` 0.61 / 2.52 / 6.36 s on the rungs
(vs T0's 1.61 / 7.24 / 15.48); P0 compute at 2.0 M is merge 2.15 + order 0.45 —
the bracket is now **communication-dominated** (gather wire + non-root scatter
wait). merge+order alone scale ~N^1.14 ⇒ ~33 s at 19.18 M + the ~10-30 s np240
gather ≈ **~1 min total** — beating the "~1-3 min" projection above. With the
setSize quicksort dead too (§Resolved), the whole measured first-step setup at
2.0 M np8 is 9.2 s vs ~2,978 s pre-ADR (~320×).*

T0-before-T2 is deliberate, on three grounds:
1. **Risk asymmetry.** T0+T1's oracle is a bit-diff (G1). T2's oracle is "same physics to
   solver tolerance" — weaker, laborious, on a path (parallel numbering consistency)
   where a subtle bug is *silently wrong physics*, the worst failure class in this code.
2. **Coverage.** T0+T1 fixes the merge for every configuration that reaches it —
   including implicit MP and future distributed factorizing solvers. T2 serves only
   diagonal/explicit.
3. **Sufficiency at K=1.** Post-T0+T1 the 43 M numbering costs ~5-10 min against a
   multi-hour solve, and the rank-0 merged graph (~11-12 GB at 43 M) fits the 60 GB nodes.

**The K× multiplier (adversarial find, changes the T2 trigger):** `domainChanged()`
re-runs on **every** domain change, and the fork's production patterns multiply K:
apeGmsh emits an unconditional per-stage `domainChange` (recorder MODEL_STAGE splitting,
apeGmsh PR #633) ⇒ staged decks pay numbering **per stage**; ADR-51 element removal bumps
the domain stamp per removal event; ADR-55 contact re-discovery likewise. A 20-event
removal history at 19 M re-crosses the hour line even post-T1. The G3 deck (single-stage
plane wave) is the **K=1 best case** and is flagged as such.

**Decision rule (ADR-40's measurement-first principle applied to our own fix):** T0+T1
now; **T2 is scheduled when a measured model makes K × the T0+T1 residual
non-negligible** — single-shot beyond ~100 M nodes, **or K ≫ 1 production decks
(staged / element-removal / contact-rediscovery) at ≥10 M**, or rank-0 memory binding.
**Flip conditions, recorded:** were the explicit lane the only lane and >43 M imminent,
skip T0 and build T2 directly; were upstream acceptance the near-term goal, write the
in-place edit first.

## No *correct* deck-level escape (all routes verified in source)

*(Corrected by the adversarial review: an earlier revision claimed "every numberer verb
in MP crosses the merge" — **false**. The `commands.cpp:3894-3898` Plain/RCM →
`ParallelNumberer` mapping is `_PARALLEL_PROCESSING` = the **SP** build
(CMakeLists.txt:884); OpenSeesMP builds with `_PARALLEL_INTERPRETERS` and maps `Plain` →
per-rank `PlainNumberer`, `RCM` → per-rank serial `DOF_Numberer(RCM)` (`:3907-3911`),
which never touch the merge. The escape verbs exist — they are just **wrong**:)*

- **Per-rank `Plain`/`RCM` in MP is a wrong-answer generator** with any distributed SOE:
  `MPIDiagonalSOE::setSize` identifies shared DOFs by *equality of global equation
  numbers across ranks* (`:216-230`, `:262-276`) and then **sums** shared entries
  (`:150-155`); its only error paths are size-sanity exits — nothing validates that the
  overlap is physically meaningful. Independent per-rank 0..n numbering ⇒ near-total
  spurious overlap ⇒ silently summed unrelated DOFs.
- **`ParallelPlain`** runs the *same* merge, then adds the `:222`/`:241` full-width
  scans: ≈ 3.7×10¹⁴ extra comparisons ⇒ ~2-3 **days** at 19 M. Strictly worse.

The code fix is mandatory. It is also stock-upstream code ⇒ **upstreamable** (T0/T1 are
behavior-identical on the RCM verb — the cleanest possible upstream PR shape).

## Gates

| Gate | Delivers | Pass criterion | Status |
|---|---|---|---|
| **G0** — scaling sweep (pre-fix) | measured `dc.numberDOF` exponent vs N at fixed np | super-linear growth across 0.25/1.0/2.0 M-hex rungs, np8 | **PASSED — 3 points: `dc.numberDOF` 39.90 / 578.87 / 2477.50 s ⇒ N^2.01**; comparison model reproduces rungs 1-2 ≤3%; step-1 isolation via the per-step wall column (48.95 s vs 0.36 s median); the merge additionally verified *live* (ADR74DBG trace: dedup + scatter correct) |
| **G1** — numbering byte-identity (RCM verb) | the new class changes no numbers | two-deck oracle (§Discovered en route): Deck A (`system Mumps`, global ids) + Deck B (`system MPIDiagonal`, end-state), np2+np8 | **PASSED at N2**: chain decks bit-identical with the O(V) engine; **production-scale artifact**: the real 0.25 M-hex 3D rung (269,314 nodes, np8, MPIDiagonal end-state) — stock `numbering_new_*` vs T0 `numbering_t0_*` in `~/ladruno_nsweep/rung-0.25M`, **all 8 ranks bit-identical**; plus the MP/`-4` fixup gate (Plain handler + double-equalDOF-on-one-node, np2+np8, bit-identical) |
| **G1b** — plain verb | valid numbering, tag-0 fix documented | bijection + same-physics + the tag-0 micro-case | **PASSED at N2**: bijection + same-physics (rel 1e-10; ordering-noise floor ~1e-15, so 5 orders of margin) on the chain; the tag-0 divergence test asserts stock-plain and Ladruno-plain dumps DIFFER while both remain valid bijections (the fix is real and visible); *scope note:* rung-scale same-physics vs stock-plain is infeasible — stock `ParallelPlain` at 0.25 M is the ~3.7×10¹⁴-comparison path (days) |
| **G1c** — fallback path | the κ-guard hash path is gated, not just written | tag-offset clone crossing κ, np2 **and** np8, with an **engagement observable** | **PASSED at N2**: strided deck (maxTag ~10⁷) bit-identical to stock at np2+np8, and the harness asserts the `RefIndex -> hash fallback` observable actually fired (an oracle that cannot fail is not an oracle) |
| **G2** — the fix, measured | T0(+T1) timing on the G0 rungs | `dc.n.merge` exponent ≈ 1; per-bracket attribution | **PASSED — three points**: `dc.numberDOF` 39.90→**1.73 s** (23×), 578.87→**7.60 s** (76×), 2477.50→**16.48 s** (**150×**); exponent **N^1.10 over 3 points** (1→2: 1.09, 2→3: 1.10); residual attributed by sub-scopes to the T1-scoped map ops (2.0 M: merge 7.91 s, order 5.18 s, gather 1.01 s, scatter 0.19 s). `dc.setSize` 3-point confirm: 10.1/101.4/501.1 s ≈ N^1.9 — the pre-G3 blocker, now with its own clean curve |
| **G3** — cluster confirmation | the production-scale kill | 18.6 M deck re-run on Esmeralda with the patched binary, **per-bracket budgets**: `dc.numberDOF` ≤ ~3 min, `dc.setSize` reported separately with its `dc.s.*` sub-brackets (quadratic sort fixed 2026-07-22, §Resolved — expect seconds: graph build + np-sequential bcast rounds), first step ≤ minutes, total ≈ 25 min; per-step wall column confirms step 1 carried setup. Flagged K=1 best case | pending cluster rebuild (G1/G2 + the setSize fix all green) |
| **G4** (T2 only) | distributed numbering correctness | same-physics vs T0 numbering (solver-tolerance, NOT bit-identical); np-sweep shows O(V/P); **ADR-30's staged gate imported**: renumbers correctly every stage | deferred until T2 is scheduled |

**Instrumentation for G2/G3:** the single `dc.numberDOF` bracket cannot prove *which*
term died — the Ladruno class adds P0-side sub-scopes **`dc.n.gather` / `dc.n.merge` /
`dc.n.order` / `dc.n.scatter`**, so G2/G3 attribute the kill and quantitatively retire
the fabric hypothesis rather than retiring it by elimination.

**SP scope: OUT for v1.** The `_PARALLEL_PROCESSING` (OpenSeesSP) lane and the
`FEM_ObjectBrokerAllClasses` case are not exercised by any gate here; the fork's
production lane is OpenSeesMP. The broker case ships only if/when an SP smoke test joins
the battery — until then the ADR does not claim SP support.

## Implementation plan — phases N0-N5 (oracle-first, the fork discipline)

Each phase is one PR-sized unit with its gate named up front; no phase starts before the
previous phase's gate is green. The harness comes **before** the class (N0) so every
later phase regresses against an instrument already proven able to fail.

| Phase | Delivers | Gate |
|---|---|---|
| **N0 — the oracle** | The G1 instrument, *before any numberer code*: a fork-only debug verb `ladrunoNumbering <file>` (additive registration) dumping per-rank sorted `(nodeTag, ndf, eqn IDs)`; `tests/test_adr74_numberer_1.py` with (a) the **bijection assertion** (every eqn 0..neqn−1 exactly once across ranks), (b) a stock-vs-stock identity run (`ParallelRCM` vs itself, np2+np8, 0.25 M rung) proving the harness reads true, and (c) a **mutation test** — perturb one equation ID in a dump copy and assert the harness FAILS (an oracle that cannot fail is not an oracle) | harness detects identity AND detects the planted mutation |
| **N1 — the class, delegate mode** | `LadrunoParallelNumberer.{h,cpp}` skeleton: ledgered promotion in `ParallelNumberer.h` (incl. `numChannels`); classTag (≥33000 numberer band); registration in both interpreters per the PR #233 flat-factory lesson; `numberDOF(ID&)` → hard error; `numberDOF(int)` **delegates to the base** (zero new behavior); apeGmsh primitive `ops.numberer.LadrunoParallelRCM()` (+ plain) with the emit-time fork gate | G1 harness: delegate vs stock **bit-identical** np2+np8; full MP battery + apeGmsh plane-wave gates green; build byte-identical for decks not using the new verbs |
| **N2 — T0 (kill the quadratics)** | Own-body `numberDOF(int)`: dense `refToLoc` (64-bit size math, κ-guard hash fallback, `ref < 0` hard error) in pass 1; per-subgraph `subTagToMerged` in pass 2; seen-bitmap plain branch (tag-0 fixed); `constrainedNode → MPs` index for the `-4` fixup; `dc.n.gather/merge/order/scatter` sub-scopes | **G1** (RCM verb bit-identical, np2+np8) + **G1b** (plain: bijection + same-physics + the tag-0 micro-case) + **G1c** (tag-offset clone crosses κ, fallback engages, still bit-identical) + G0 rungs re-run: `dc.n.merge` exponent ≈ 1 |
| **N3 — T1 (kill the map residual)** | Owned `tag → Vertex*` vector maintained through the merge; `Vertex::addEdge` called directly on both endpoints; merged graph on dense array-backed storage; RCM runs against it | **G1 again** (the critical regression — bit-identity must survive T1) + **G2 full**: 3-rung sweep on the fixed binary, per-bracket table, projections vs the §Sequencing table |
| **N4 — G3, the production kill** | Cherry-pick/rebase instrumentation + N1-N3 onto the cluster tree (`ladruno-p5-build` — the three touched files are hash-verified identical); rebuild Esmeralda binary; re-run the 18.6 M deck with `numberer LadrunoParallelRCM` (and once with `LadrunoParallelPlain`) under per-bracket budgets | **G3**: `dc.numberDOF` ≤ ~3 min; `dc.setSize` reported separately (its ~N² anomaly investigated **before** this run, see Risks); first step ≤ minutes; total ≈ 25 min; per-step column shows step 1 |
| **N5 — ship hygiene** | LEDGER_vanilla_files row (the promotion); LEDGER_quirks entries (Lagrange `ref=−1` fusion, tag-0, broken `numberDOF(ID&)`); banner regen; CHANGELOG; the §Risks external-record corrections applied to the two Obsidian notes + ADR-30 header; apeGmsh plane-wave family decks switched to the Ladruno verb; upstream-PR extraction decision recorded **+ the apeGmsh typed primitives `ops.numberer.LadrunoParallelRCM/Plain` (slipped from N1 — tracked here so it cannot be lost)** | docs land in the same PR as the switch (the build-control rule); vault corrections verifiable |

**Validation principles pinned:** "same-physics" (G1b/G4) means the apeGmsh plane-wave
gate metric (response match ≤1.5e-15 on the 0.25 M rung); every gate that can be run at
np2 runs at np2 *and* np8 (two channel topologies); the `dc.setSize` anomaly is a
**blocking pre-G3 task** so G3's budget table attributes honestly; and no phase may
widen its scope — a discovered defect goes to the Risks list, not into the running PR
(this ADR was under-scoped twice during review; the antidote is scope discipline).

## Discovered en route (N0, 2026-07-22): MPIDiagonalSOE localizes the numbering

The N0 oracle's first real run caught a semantic this ADR (and both prior reviews) did
not know: **under `system MPIDiagonal`, the globally-consistent numbering the
ParallelNumberer computes exists only transiently.** `MPIDiagonalSOE::setSize` first
builds its shared-DOF exchange map by intersecting the *global* equation numbers
(`:216-276`), then deliberately **rewrites every DOF group's ids to rank-local
0..n_local−1** (`MPIDiagonalSOE.cpp:408-429`, "renumber DOFs 0 through size",
`newDOF = myDOFs.getLocationOrdered(dof); dofPtr->setID(i, newDOF)`) and has the FE
elements re-cache. Established by runtime trace (temporary `ADR74DBG` prints):
`numberDOF` demonstrably assigns correct, consistent globals (merge + dedup + scatter
all verified live at np2 — `mergedVertex=11` from 6+6, scatter pairs correct), and the
post-analysis ids are those globals compacted per rank (exact −offset match). MUMPS
SOEs do **not** renumber (no `setID` in `Mumps*SOE`). Consequences:

- **No physics defect** — shared-DOF exchange is built from consistent globals before
  localization; past MP runs (cluster and local) are untouched.
- **G1's oracle is two-deck** (implemented in `tests/test_adr74_numberer_1.py`):
  Deck A (`system Mumps`, global ids survive) carries the *strict* numberer-identity
  oracle — shared-node cross-rank equality + global bijection + byte-identity; Deck B
  (`system MPIDiagonal`, the production stack) carries the *end-state* oracle —
  per-rank dense local bijection + byte-identity. Post-localization identity is
  necessary but not sufficient for numberer identity (order-preserving compaction can
  alias distinct global numberings), hence Deck A is the load-bearing G1 gate.
- **The irony sharpens again**: on the explicit lane the merge's RCM output is not only
  never factorized — it is *discarded within the same domainChanged call* by the
  localization. T2 is further strengthened: any consistent global bijection suffices,
  because MPIDiagonal re-localizes regardless.
- The stall analysis is untouched: the merge provably runs (live trace + measured
  N^1.98-2.01 + model reproduction ≤3%) before the localization discards its ordering.

## Resolved (2026-07-22): the `dc.setSize` ~N^1.9 anomaly — quicksort on sorted input

The pre-G3 blocker (§Risks) is root-caused, fixed, and measured. **None of the ADR's
original suspects were it** — the np-round bcast+intersections loop is O(P·N/P) with a
linear two-pointer merge, `ID::insert` is not on the path, and the localization pass is
O(N log N) binary searches. The quadratic is one line earlier:
`MPIDiagonalSOE::setSize` copies its local DOF list out of the DOF graph and sorts it
with a **hand-rolled leftmost-pivot quicksort** (`quickSort`/`q_sort`,
`MPIDiagonalSOE.cpp`). The DOF graph is `MapOfTaggedObjects`-backed (`std::map`,
`AnalysisModel::getDOFGraph`), so its vertex iterator hands the equation numbers over in
**ascending order — the array arrives already sorted**, and leftmost-pivot quicksort on
sorted input is the textbook deterministic O(n²/2)-comparison worst case, per rank,
every domain change. (Bonus defect: O(n) recursion depth on the same input.)

**Why it eluded the np240 incident analysis:** the cost is O((N/P)²) — *per-rank*
quadratic. At np240/19.18 M it is ~30-60 s, invisible inside the 16.35 h numberer stall;
at np8 rungs it is the dominant post-T0 term. Same reason it never hurt small-np decks
for 20 years. Like the numberer term it multiplies by K (per domainChanged).

**Attribution (dc.s.\* sub-brackets, np8, max over ranks; pre-fix = instrumented stock
sort, this session — the G2-era 10.1/101.4/501.1 s bracket totals reproduced as
6.54/93.54/(not rerun) under warmer machine state):**

| bracket | 0.25 M pre → post | 1.0 M pre → post | 2.0 M pre* → post |
|---|---|---|---|
| `dc.setSize` | 6.54 → **0.37 s** (17.7×) | 93.54 → **1.61 s** (58×) | 501.1* → **3.11 s** (**161×**) |
| `dc.s.sort` (the kill) | 6.23 → 0.00 | **91.93 → 0.00** | — → 0.01 |
| `dc.s.graph` (residual #1) | 0.31 | 1.39 | 2.50 |
| `dc.s.bcast` | 0.37 → 0.07 | 1.86 → 0.35 | — → 0.58 |
| everything else | ≤0.01 | ≤0.03 | ≤0.08 |

\*2.0 M pre-fix from the N2 T0 profiles (bracket total; predates the sub-brackets).
Pre-fix `dc.s.sort` measured exponent **N^1.99** (6.23→91.93 for ×3.87 nodes) with
95-98% of the bracket; the n²/2-comparison model closes all three rungs (5.7×10⁹ /
8.6×10¹⁰ / 3.4×10¹¹ compares ≈ 6-10 / 86-143 / 341-568 s predicted). Post-fix the
bracket is **linear** (N^1.08 / N^0.97) and dominated by `dc.s.graph` — the
`getDOFGraph` build, ~1.2 s/M-node/rank-share, fine until measured otherwise.

**The fix** (in-place ledgered vanilla edit, NOT a subclass — one call site, output
provably identical): `std::sort(myDOFsArray, myDOFsArray+size)` replaces
`quickSort(myDOFs, size)`; `q_sort` is kept for provenance. Tags are unique, so any
correct sort yields the identical ascending array ⇒ identical `myDOFs`, identical
shared-DOF intersections, identical localization. **Gates:** (a) end-state numbering
dumps (`ladrunoNumbering`) byte-identical across all 8 ranks on the real 0.25 M rung,
fixed vs stock-sort binary — and the stock-sort instrumented dumps byte-match the N2
`numbering_t0_*` artifacts (oracle-chain null check); (b) `tests/test_adr74_numberer_1.py`
full suite green on the fixed binary; (c) the 3-rung post-fix curve above.

**G3 implication:** at 19.18 M/np240 the setSize budget drops from ~30-60 s of sort to
seconds (graph build + bcast rounds); with T0 the whole first-step setup at cluster
scale is now projected ≤ ~5 min pre-T1. The `dc.s.*` brackets ship, so G3 attributes
setSize honestly per its gate row. `DistributedDiagonalSOE` (the SP sibling) does not
sort — grep confirms `MPIDiagonalSOE` held the only copy of `q_sort` in the tree.

## Resolved (2026-07-22): `dc.handle` ~N^2 — the getLocation-in-loop family, twice over

The third setup quadratic, flagged by the owner from the collector table
(0.26 / 3.53 / 14.40 s over the rungs ⇒ N^1.93-2.00) and confirmed per-rank by a
fixed-V np-sweep (falls ~1/P²: 3.25/0.85/0.29/0.14 s at np2/4/8/16 —
`~/ladruno_nsweep/psweep-np*`). Two mechanisms, attributed by new `dc.h.*`
sub-brackets (`splists/nodes/eleclass/fecreate`), fixed in two vanilla files:

1. **`TransformationConstraintHandler::handle()`** — `ID::getLocation` linear
   scans inside the node and element loops: the element-classification loop
   (O(#ele × nodes/ele × #SP) — measured dominant, `dc.h.eleclass` 2.79 s of
   3.62 at 1.0 M, N^1.94), the per-node `constrainedNodesSP/MP` scans + O(#SP)
   forward re-scans, the SP/MP list-build membership scans, and the FE-creation
   `transformedEle.getLocation` (O(#ele × #transformed)). On slab decks
   (fixed-depth: the production 20×20 km meshes and every rung) the constrained
   base area ∝ N, so #SP ∝ N and all of them go quadratic. Fix: one
   `unordered_set` of constrained nodes + first-index maps + per-node SP index
   lists (ascending = the stock loc+rescan visit order) + a `transformedEle`
   set twin — **searches only**; every mutation, creation order, and branch
   quirk (incl. the stock `numSPConstraints != 0` MP-path oddity) preserved.
2. **`TransformationDOF_Group` SP-only ctor** — exposed only after fix 1: the
   ctor swept **every SP in the domain** per constrained node
   (`dc.h.nodes` residual 1.27 s at 2.0 M, ~N^2.3). The sweep is provably
   redundant: both in-tree callers (TransformationConstraintHandler and
   AutoConstraintHandler) follow the ctor with `addSP_Constraint` for every SP
   of the node from `getDomainAndLoadPatternSPs` ⊇ `getSPs`, and
   `addSP_Constraint` sets the same `theSPs[dof]` pointer — identical end
   state. Upstream itself had already commented out the equivalent sweep in
   the MP ctor. Removed (original kept in a comment).

**Measured (np8, max over ranks):** `dc.handle` 0.28 / 3.62 / 14.40 →
**0.04 / 0.17 / 0.34 s** (7× / 21× / **42×**), post-fix exponent N^1.0-1.07;
all sub-brackets ≤ 0.22 s at 2.0 M. **Gates:** end-state numbering dumps
byte-identical to the N2 `numbering_t0` artifacts on all 8 ranks (the fix is
upstream of numbering — DOF_Group creation order feeds it — so this gate is
load-bearing, and it was run twice: after fix 1 and after fix 1+2);
`test_adr74_numberer_1.py` 18/18 (twice, same points). With all three
quadratics dead, the measured first-step `domainChanged` at 2.0 M np8 totals
**9.6 s vs ~2,993 s pre-ADR (~313×)**, every remaining term ~linear.
Cluster note: per-rank quadratic ⇒ the np240 impact was small (~1-2 s inside
the incident); the scenarios this fix rescues are the 43 M P=8 retry (~1.6 h
avoided) and every np8-class local sweep. The K× domainChanged multiplier
(§Sequencing) applied to this term too and is now retired with it.

## Instrumentation (shipped with this investigation, branch `perf/step-stall-instrumentation`)

Commit `e0113e53a` — two profiler gaps closed, prerequisites for G0-G3:
1. the `-perStep` wall column now records real per-step wall (was: `0.0f` pushed,
   `phases` never seeded ⇒ zero-width column ⇒ the old profiles physically cannot say
   *which* step stalled). **First use isolated the stall to step 1 by measurement**
   (48.95 s vs 0.36 s median on the 0.25 M rung).
2. `domainChanged` is bracketed (`dc.handle` / `dc.numberDOF` / `dc.setSize` at
   `DirectIntegrationAnalysis.cpp:427/440/451`) inside the `step` scope.

*Harness note (lesson recorded):* the sweep collector initially validated against a
synthetic HDF5 it generated itself and missed the writer's real `rollup/root/...` level —
a self-consistent-fake-test trap. Fixed by verifying against a real `profile_0.h5`;
readers of profiler output must anchor paths at `runs/<id>/rollup/root/step/...`.

## Fence — what this ADR is and is NOT

- **Not the element-kernel lane.** Steady state was already at the ~18 µs/hex·step floor;
  per-step optimization is ADR-68 (cluster tree). This ADR is the **per-domain-change**
  axis — priced at K× per §Sequencing, not once.
- **Not a SIMD project.** Vectorizing the scan (AVX2 compare+movemask) buys a realistic
  2-3× against ~10⁵× from T0. Rejected — see Decision log.
- **Not a general HPC stack adoption.** ADR-40's anti-goal stands *except* where this
  incident measured the bottleneck. ParMETIS/PT-Scotch ordering enters only if a
  distributed *fill-reducing* ordering is ever needed — and MUMPS/SuperLU_DIST prefer
  their own internal orderings, so the likely T2 shape never needs it.
- **RCM survives.** T0/T1 keep global RCM byte-identical for every lane that factorizes.
  T2's skip-the-ordering mode is scoped to diagonal/explicit configurations and opt-in.
  *Staged-mix note:* numberer choice is per-analysis; a gravity stage on a factorizing
  solver may want `LadrunoParallelRCM` while the dynamic stage opts into T2 — supported
  by the machinery (each stage renumbers anyway), stated here so nobody treats it as a
  conflict.
- **The 2.5 GbE fabric hypothesis is now arithmetically cornered, retired at G3.** The
  gather is bulk-bounded to ~10-30 s (§Where), G0 measured the quadratic by name, and the
  comparison model reproduces the measured rungs — but only the G3 cluster run
  extinguishes it formally.

## Decision log

- **Direct-address vector over `unordered_map` as primary** — tags are dense 1..N in the
  apeGmsh pipeline; 76 MB flat array beats ~1-1.5 GB of hash nodes and malloc churn.
  Fallback keeps correctness for sparse tag spaces (G1c-gated). (Fable review, folded.)
- **`getLocationOrdered` rejected** — `vertexRefs` is append-ordered and its order is
  load-bearing (positional correspondence with `vertexTags`); sorting breaks the mapping.
- **SIMD rejected as the stall fix** — compute-limited scan, 2-3× available vs ~10⁵× from
  T0; SIMD effort belongs on the per-step kernels (ADR-68 lane).
- **Fix covers pass 2, the plain branch, and the `-4` fixup — not just pass 1** — pass 2
  dominates at low np (G0 measured it), the plain branch is ~2× worse than the merge,
  and the `-4` loop is a second quadratic on tie-heavy decks. (Fable + 3-lens reviews;
  the original proposal under-scoped twice.)
- **T1 redesigned to owned dense structures** (3-lens review, two lenses converging) —
  the `startAddEdge`/`addEdgeFast` borrow fails on snapshot staleness (O(V·P) rebuilds,
  invisible at np8, a G2-passes/G3-stalls trap) and never accelerates RCM's reads. The
  owned `tag → Vertex*` + `Vertex::addEdge`-direct + dense-storage-graph design is less
  code, keeps bit-identity (sorted adjacency is insertion-order-canonical), and removes
  the only load-bearing dependency on Graph internals.
- **`numberDOF(ID&)` scoped out** — broken upstream (numbering call commented out,
  mismatched recv layout); the override hard-errors; the defect is ledgered. No
  bit-identity obligation to code that never worked.
- **Fork-native class over in-place edit** (owner, 2026-07-21) — honors the additive
  contract; keeps stock `ParallelRCM` alive as the in-binary G1 reference; the upstream
  PR is later extracted as the in-place form. **Alternatives weighed, and the named
  smell:** the subclass is *inheritance for plumbing reuse, not substitution*. Accepted
  eyes-open: (a) in-place is upstream-shaped but rejected on governance; (b) a standalone
  `DOF_Numberer` sibling avoids the promotion but duplicates the channel wiring — though
  the review noted the MP-only duplication is just ~35 trivial lines, so this trade flips
  if SP stays permanently out of scope; (c) a merge-strategy object is over-engineering.
  Precedent: LadrunoHHT / LadrunoGeneralizedAlpha (classTags.h 33013/33014).
- **T2 is opt-in, not a default swap** — it abandons global RCM; diagonal/explicit lanes
  lose nothing; every other lane keeps today's path until measured.
- **Plain-verb non-identity is a feature, documented** — the tag-0 quirk fix and the
  no-RCM fast path make `LadrunoParallelPlain` the recommended explicit-lane verb, gated
  by bijection + same-physics rather than bit-identity.

## Risks / open questions

- **`dc.setSize` ~N^1.9 — RESOLVED 2026-07-22** (see §Resolved above): the mechanism
  was none of the suspects listed here (intersections / `ID::insert` / localization)
  but a leftmost-pivot quicksort fed already-sorted input by the map-backed DOF graph
  — O((N/P)²/2) per rank. Fixed with `std::sort` (output-identical, ledgered);
  501.1 → 3.11 s at the 2.0 M rung, post-fix curve linear, byte-identity + 18/18
  suite green. The "~9 h at 19.18 M" naive extrapolation was wrong in the reassuring
  direction: the term is per-rank quadratic, so np240 pays only ~30-60 s — but it
  would have poisoned every np8-class local sweep and the G3 attribution.
- **maxTag density assumption** (T0 primary path): κ-guard + hash fallback, G1c-gated;
  negative-ref micro-case included in the same test file.
- **Rank-0 memory after T0**: the merged graph (~5 GB at 19 M vertices) lives on rank 0
  until T2. Fine at 60 GB nodes; a printed diagnostic (vertices + est. bytes) fires when
  V > 5 M.
- **G1 dump mechanics**: per-rank sorted `(nodeTag, eqn IDs)` file — decide debug verb vs
  test-only fixture at implementation.
- **Upstream PR**: T0/T1 (RCM verb) are byte-identical and upstream-shaped; follow the
  fork's upstream-divergence ledger process; coordinate before opening.
- **External-record corrections owed** (so the vault stays truthful): the 18.6 M note §8
  "first-collective over the fabric" hypothesis (superseded by this mechanism, pending
  G3) and §10 to-do list; the 14 M note §7 449 µs interpretation; the Jun-11→Jul-18
  binary date in both notes' §6; ADR-30's "roughly linear / ~3 min" projection (struck
  here).
- **Edge cases verified, recorded as non-issues:** empty partition dies earlier
  (`AnalysisModel::getDOFGroupGraph` exits on 0 vertices); the direct-address table is
  P0-only by construction (only the `processID == 0` branch merges).

## Adversarial review log

**2026-07-21 — Fable cross-review (pre-ADR):** mechanism + magnitude confirmed
independently; four findings folded (compute-limited not bandwidth-limited; post-T0 map
residual ⇒ T1 added; pass 2 + plain branch added to T0's scope; per-rank Plain confirmed
*incorrect*, not merely slow).

**2026-07-21 — 3-lens gate (refutation / design / completeness), all findings folded:**
- **REFUTED (1):** "every numberer verb in MP crosses the merge" — wrong build macro
  (`_PARALLEL_PROCESSING` is SP; MP's Plain/RCM are per-rank). Corrected to "no *correct*
  escape"; conclusion survives via the verified silent-wrong-answer mechanism.
- **WEAKENED (2):** the "±5% closure" (achieved only by dropping the ADR's own pass-2
  term; reframed as order-of-magnitude + ≤3% np8-rung reproduction — stronger evidence,
  honest precision); T1's "low minutes" via `startAddEdge`/`addEdgeFast` (snapshot
  staleness + untouched RCM reads; redesigned to owned dense structures).
- **SURVIVED (verified at cited lines):** the O(V²) mechanism and count; ID append
  amortization; `ParallelPlain` magnitudes; the MPIDiagonalSOE silent-summation
  mechanism; the post-T0 residual; `getFreeTag` contiguity; the 14 M retro-fit; subclass
  feasibility (promotion genuinely required; both `numberDOF` virtuals dispatch).
- **NEW findings folded:** ADR-30 supersession (name/path collision + invalidated linear
  projection); the K× domainChanged multiplier (new T2 trigger); the `ref = −1`
  Lagrange UB hazard; the `-4`×MPs second quadratic; `numberDOF(ID&)` broken upstream;
  broker factory correction (`FEM_ObjectBrokerAllClasses`, not the stub); G1
  dump/bijection definition + G1b/G1c; `dc.n.*` sub-scopes + per-bracket G3 budgets;
  registration surface per the HHT ledger row + PR #233 flat-factory lesson; wire-volume
  verified-negative; tag-0 quirk accidentally benign; ctor-grenade + ownership notes;
  `numChannels` added to the promotion list; prose fixes (239 subgraphs; exponent-
  dilution direction).

## Implementation log

- **2026-07-21** — Root cause identified from source (owner) after the 144090 profiler
  showed 58,876 s in the untimed step region.
- **2026-07-21** — Profiler gaps closed: commit `e0113e53a` (per-step wall column +
  `dc.*` attribution). Verified: 10/10 standalone harness; scope strings present in the
  rebuilt Windows binary.
- **2026-07-21** — Fable adversarial cross-review; findings folded (see log above).
- **2026-07-21** — G0 sweep, np8 local (0.25/1.0/2.0 M-hex rungs, harness
  `Dropbox/UANDES EC/Cluster Max Models/numberer_sweep/`). **Rungs 1-2 measured:
  `dc.numberDOF` 39.90 s → 578.87 s (×14.5) for ×3.87 nodes ⇒ N^1.98**; step-1 wall
  isolated by the new per-step column (48.95 s vs 0.36 s median; 693.71 s vs 1.53 s);
  comparison model reproduces both rungs ≤3%. `dc.setSize` anomaly logged (Risks).
  Rung 3 (2.0 M) in flight; final table on completion. *(Raw-wall exponent reads ~1.9 —
  slightly below the bracket's 1.98 because the additive linear solve term dilutes the
  raw ratio downward.)*
- **2026-07-21** — 3-lens adversarial gate run; one claim refuted, two weakened, root
  cause stands; this revision folds all findings (supersedes-ADR-30 section, T1
  redesign, corrected escape claim, honest closure framing, expanded gates G1/G1b/G1c,
  `dc.n.*` sub-scopes, K× sequencing, ref<0 / ctor / -4 / ID-variant guards).
- **2026-07-22** — **N0 shipped**: `ladrunoNumbering` verb (classic Tcl, mirrors the
  `profiler` verb pattern) + `tests/test_adr74_numberer_1.py` (two-deck oracle +
  mutation tests). G0 completed: third rung landed — `dc.numberDOF` 39.90 / 578.87 /
  2477.50 s over 0.25/1.0/2.0 M hex ⇒ **N^2.01 over three points**. N0's first real
  run caught the MPIDiagonal id-localization semantic (§Discovered en route),
  established by live `ADR74DBG` trace of `numberDOF` (merge/dedup/scatter verified
  correct at np2; prints reverted after use). Harness lessons banked: mpiexec exits 0
  on Tcl-init failure and on script error — dump-count/content assertions are the real
  failure signal, rc is not; build-tree exes need TCL_LIBRARY (the openseesmp.sh gap,
  now handled in-harness).
- **2026-07-22** — **N1 shipped**: `LadrunoParallelNumberer` (classTag 33000) in
  delegate mode; ledgered promotion in `ParallelNumberer.{h,cpp}` (private→protected +
  two tag-forwarding ctors, incl. the no-numberer form the plain verb needs); verbs
  `LadrunoParallelRCM`/`LadrunoParallelPlain` registered in classic Tcl + the
  interpreter; ledger rows appended. **N1 gate PASSED 13/13**: delegate-vs-stock
  bit-identical at np2+np8 on Mumps-RCM (strict), MPIDiag-RCM (production end-state),
  and MPIDiag-Plain — simultaneously the null-test of the N0 oracle. Next: N2 (T0).
- **2026-07-22** — **N2/T0 shipped**: own-body `numberDOF(int)` — RefIndex dense/hash
  dedup (64-bit sizing, κ-guard), per-subgraph pass-2 map, seen-flag plain branch
  (tag-0 fixed), `ref<0` Lagrange hard-guard, `constrainedNode→MPs` index for the -4
  fixup, `dc.n.gather/merge/order/scatter` sub-scopes. **Gates 14/14**: G1 RCM verb
  bit-identical to stock at np2+np8 on both decks *with the O(V) engine underneath*;
  G1b plain bijection + same-physics; G1c strided-tag hash fallback bit-identical.
  **G2 measured**: 23× / 76× on the rungs, exponent N^2.01 → **N^1.09**; residual =
  the T1-scoped map ops, attributed by the new sub-scopes. `dc.setSize` (~N^2.06) is
  now the dominant first-step term — the pre-G3 blocker (resolved below, 2026-07-22).
- **2026-07-22** — **setSize quadratic killed** (§Resolved; the pre-G3 blocker):
  read-the-source found the leftmost-pivot `q_sort` on the already-sorted DOF list
  (map-iteration order) — O((N/P)²/2)/rank, closing all three rung measurements to
  ~20% via the n²/2-compare model; `dc.s.*` sub-brackets
  (`graph/fill/sort/bcast/shared/localize/eleid`) added inside `dc.setSize` attributed
  95-98% to `dc.s.sort` at measured **N^1.99** pre-fix; `std::sort` fix ⇒
  `dc.setSize` 6.54/93.54/501.1 → **0.37/1.61/3.11 s** (17.7×/58×/**161×**),
  post-fix linear, residual = the `getDOFGraph` build. Gates: 8-rank byte-identity
  (fixed vs stock vs the N2 `numbering_t0` artifacts) + `test_adr74_numberer_1.py`
  **18/18** on the fixed binary. In-place ledgered vanilla edit (one call site);
  `q_sort` retained; sole copy in tree (SP's `DistributedDiagonalSOE` does not sort).
  Shipped as PR #593 (with the missing 33000 manifest row, a pre-existing N1 gap the
  G9 gate surfaced).
- **2026-07-22** — **N3/T1 shipped**: the merged graph now lives on an OWNED
  dense `ArrayOfTaggedObjects`-backed Graph with a `tag → Vertex*` vector
  maintained through the merge — P0's vertices are copied in (tag/ref/weight/
  color + adjacency; contiguity live-checked), `getFreeTag` sequence preserved
  (nextFreeTag == numVertexP0 after the copy), pass-2 edges go in via direct
  `Vertex::addEdge` on both endpoints (sorted `ID::insert` ⇒ order-canonical),
  and RCM runs on the dense graph so every BFS `getVertexPtr` (including the
  GPS trial passes — a full BFS per candidate start vertex) is an array index.
  The model's map-backed group graph is left unmerged. **Gates: G1 SURVIVES T1**
  — 18/18 np2+np8 both decks + the 0.25 M rung end-state dumps byte-identical
  to the N2 `numbering_t0` artifacts on all 8 ranks. **G2 full**: `dc.numberDOF`
  1.61/7.24/15.48 → **0.61/2.52/6.36 s**; P0 merge 7.91→2.15, order 5.18→0.45
  at 2.0 M — the bracket is now gather/scatter (communication) dominated;
  merge+order ~N^1.14 ⇒ ~33 s at 19.18 M ⇒ numbering ≈ ~1 min at cluster scale.
  First-step setup at 2.0 M np8: ~2,978 s pre-ADR → **9.2 s** (~320×).
  Fork-file-only change (`LadrunoParallelNumberer.cpp`) — no ledger row needed.
  Next: N4/G3, the Esmeralda kill-run.
- **2026-07-22** — **`dc.handle` quadratics killed** (§Resolved above; owner-flagged
  from the collector table + fixed-V np-sweep): `dc.h.*` sub-brackets attributed
  the N^1.94 element-classification scans and, once those died, the hidden
  `TransformationDOF_Group` ctor domain-SP sweep (N^2.3 residual). Both fixed
  (searches-only handler rewrite + provably-redundant sweep removal);
  0.28/3.62/14.40 → **0.04/0.17/0.34 s** (42× at 2.0 M), byte-identity ×2 +
  18/18 ×2. First-step `domainChanged` at 2.0 M np8: ~2,993 → **9.6 s (~313×)**,
  all terms now ~linear. Also this session: the fixed-V np-sweep
  (`psweep-np{2,4,16}`) established `dc.setSize` post-fix FALLS ~1/P (bcast term
  0.04→0.07 s np2→16 — no structural np-wall) and the interconnect-study
  caution (the pre-T0 merge is np-invariant global-V² ⇒ mimics an Amdahl
  fraction in strong scaling; the np128 verdict needs the dc.*-bracketed re-run,
  folded into G3).
- **2026-07-22** — **N2 adversarial gate (3 lenses) folded**: transcription-fidelity
  verdict "bit-identity HOLDS for every in-tree path" (mutation-for-mutation verified;
  κ paths answer-identical; -4 order preserved; two stock defects incidentally fixed —
  a leaked V-sized ID per domainChanged on the RCM path and a null-deref on the
  non-root error path). Converged MAJORs applied: fail-STOP (exit −1) at both
  mid-collective error sites (return = distributed deadlock + ignored result on the
  production caller), pass-2 `find()` guard (operator[] would silently wire spurious
  edges to vertex 0), `nOrdered` completeness assertion (the tag-contiguity theorem is
  now live-checked), int64 sizing-hint math. Gates added per review: the MP/`-4` fixup
  identity gate (previously ZERO coverage of the rewritten second quadratic), the
  tag-0 divergence micro-case, the Lagrange `ref<0` loud-failure test, G1c np8 + the
  hash-engagement observable. *Provenance*: all rung artifacts stamp the N1-merge
  banner hash (`d04077d3c`) + the T0 working tree = branch `feat/adr74-n2-t0`; the
  review-fix binary re-ran the full suite green.
- **Cluster note**: the Esmeralda tree (`ladruno-p5-build` @ #580) carries the same
  unfixed code — all three instrumentation-touched files byte-identical to the patch
  base — so the instrumentation + T0/T1 rebase cleanly there. The live cluster binary is
  Jul-18 (not the Jun-11 build the Obsidian notes record).
