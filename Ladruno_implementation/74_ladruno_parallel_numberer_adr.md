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
> *Family-link note:* ADR-60–73 live on the cluster tree (`ladruno-p5-build`), not on this
> branch — the 59→74 numbering jump is deliberate (74 = next free on the newest tree).
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
| **G0** — scaling sweep (pre-fix) | measured `dc.numberDOF` exponent vs N at fixed np | super-linear growth across 0.25/1.0/2.0 M-hex rungs, np8 | **2 of 3 rungs in: `dc.numberDOF` 39.90 s → 578.87 s for ×3.87 nodes ⇒ N^1.98**; comparison model reproduces both rungs ≤3% (47.6/711 s predicted vs 48.95/693.71 s step-1 walls); rung 3 in flight |
| **G1** — numbering byte-identity (RCM verb) | the new class changes no numbers | per-rank sorted `(nodeTag → eqn IDs)` dump, `LadrunoParallelRCM` vs stock `ParallelRCM`, same partition, 0.25 M rung, np2 + np8: **bit-identical**; plus a **global bijection assertion** (every eqn 0..neqn−1 exactly once); vehicle `tests/test_adr74_numberer_*.py` + the apeGmsh plane-wave MP gates (≤1.5e-15) named as the battery | pending T0 |
| **G1b** — plain verb | valid numbering, tag-0 fix documented | bijection assertion (stock plain *fails or double-orders* the constructed tag-0 micro-case — documented); same-physics vs stock on the 0.25 M rung | pending T0 |
| **G1c** — fallback path | the κ-guard hash path is gated, not just written | tag-offset clone of the 0.25 M deck (node tags × large stride) provably crossing κ: fallback engages **and** numbering still bit-identical to stock | pending T0 |
| **G2** — the fix, measured | T0(+T1) timing on the G0 rungs | `dc.n.merge` exponent ≈ 1; per-bracket attribution via the sub-scopes below. *Validity note:* pass 1 is P-independent (Σ ≈ V²/2) and pass 2 ∝ V²/P, so np8 exponents transfer to np240 conservatively; the genuinely P-dependent hazards (snapshot rebuilds, rank-0 memory) are addressed by the owned-structure T1 design | pending T0 |
| **G3** — cluster confirmation | the production-scale kill | 18.6 M deck re-run on Esmeralda with the patched binary, **per-bracket budgets**: `dc.numberDOF` ≤ ~3 min, `dc.setSize` reported separately (np-sequential bcast rounds at np240 — and see the open question below), first step ≤ minutes, total ≈ 25 min; per-step wall column confirms step 1 carried setup. Flagged K=1 best case | pending G1/G2 + cluster rebuild |
| **G4** (T2 only) | distributed numbering correctness | same-physics vs T0 numbering (solver-tolerance, NOT bit-identical); np-sweep shows O(V/P); **ADR-30's staged gate imported**: renumbers correctly every stage | deferred until T2 is scheduled |

**Instrumentation for G2/G3:** the single `dc.numberDOF` bracket cannot prove *which*
term died — the Ladruno class adds P0-side sub-scopes **`dc.n.gather` / `dc.n.merge` /
`dc.n.order` / `dc.n.scatter`**, so G2/G3 attribute the kill and quantitatively retire
the fabric hypothesis rather than retiring it by elimination.

**SP scope: OUT for v1.** The `_PARALLEL_PROCESSING` (OpenSeesSP) lane and the
`FEM_ObjectBrokerAllClasses` case are not exercised by any gate here; the fork's
production lane is OpenSeesMP. The broker case ships only if/when an SP smoke test joins
the battery — until then the ADR does not claim SP support.

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

- **`dc.setSize` shows its own ~N² at np8** (G0 data: 6.23 s → 100.23 s, ×16.1 for ×3.87
  nodes ⇒ ~N^2.06) — mechanism unidentified (`MPIDiagonalSOE::setSize`'s per-pair
  intersection is linear; suspicion falls on the shared-eq list construction or
  `ID::insert` behavior at scale). **Separate defect, needs its own look before G3** —
  budgeted separately in G3 precisely so it cannot masquerade as numberer residue.
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
- **Cluster note**: the Esmeralda tree (`ladruno-p5-build` @ #580) carries the same
  unfixed code — all three instrumentation-touched files byte-identical to the patch
  base — so the instrumentation + T0/T1 rebase cleanly there. The live cluster binary is
  Jul-18 (not the Jun-11 build the Obsidian notes record).
