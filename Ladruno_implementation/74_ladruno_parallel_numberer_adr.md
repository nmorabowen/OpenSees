---
title: ParallelNumberer O(V²) merge — the 16.35 h first-step stall (root cause + fix ladder)
project: Ladruno
status: proposed — root cause CONFIRMED by closing arithmetic (two independent derivations, ±5%); T0 ready to implement; G0 scaling sweep in flight
priority: high
owner: nmora
amends: 40_ladruno_performance_adr
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
> that path). Family: ADR-40 (perf program) · ADR-40a (instrument scopes) · ADR-67
> (integrator lane) · ADR-68 (element/state-determination lane — the *other* per-step axis;
> this ADR is the one-time-setup axis) · ADR-06 (profiler — the instrument that both missed
> and then caught this, see §Instrumentation).
> External record: Obsidian `Cluster Scalability/Max Model — 18.6M P6 Cluster-Emit
> (20x20km).md` §8 (the incident), `Max Model — 14M-Hex Hero Run.md` §7 (retro-explained
> below). Adversarial cross-review: Fable agent, 2026-07-21 (findings folded, §Decision log).

## What

`ParallelNumberer::numberDOF()` — the engine under **both** `numberer ParallelRCM` and
`numberer ParallelPlain` in OpenSeesMP — gathers every rank's DOF-group graph onto rank 0
and merges them with per-vertex **linear scans** (`ID::getLocation`), making first-step
setup **O(V²_global) serial on one core** while every other rank busy-polls in `recvID`.
At 19.18 M nodes that is ~1.8×10¹⁴ integer comparisons ≈ 16-17 h — the observed stall, to
within 5%. This ADR ships a three-tier fix, each independently landable and gated:

**Delivery vehicle (owner decision 2026-07-21): a fork-native class, not an in-place
edit.** `LadrunoParallelNumberer : public ParallelNumberer` (new files, ladruno numberer
classTag in the ≥33000 band), following the fork's additive discipline: upstream
`ParallelNumberer` stays behavior-untouched; the only upstream-file diffs are (a) a
**ledgered private→protected promotion** in `ParallelNumberer.h` (`processID` /
`theChannels` / `theNumberer` — same pattern as the HHT.h and GeneralizedAlpha.h
promotions for LadrunoHHT/LadrunoGeneralizedAlpha, classTags.h 33013/33014), and (b)
additive registration verbs. `mergeSubGraph` is already `protected`. The tiers:

- **T0 — kill the quadratic (the stall fix):** override `numberDOF` with a dense
  direct-address lookup (node tags are 1..N from gmsh/apeGmsh) + hash-map fallback,
  applied to **both** merge passes *and* the plain-ordering branch. O(V²)→O(V).
  Ordering-stable by construction ⇒ **bit-identical equation numbering vs stock
  `ParallelRCM`**, enforced by gate G1.
- **T1 — kill the tree lookups (the residual):** after T0 the wall is 20-45 min of
  `std::map` finds (2 per directed adjacency entry in `Graph::addEdge`, 1 per entry in
  RCM's BFS). Build the merged graph through the in-tree dense path
  `Graph::startAddEdge`/`addEdgeFast` (`Graph.cpp:223-276`, virtual on Graph — no Graph
  edits needed) — merged tags are contiguous from `getFreeTag`. Numbering drops to low
  minutes.
- **T2 — stop centralizing (the architecture):** owner-computes distributed numbering
  (owner rule on shared nodes, per-rank owned-DOF counts, `MPI_Exscan` start offsets, one
  halo exchange). O(V/P) time and memory, no rank-0 merged graph (~5 GB at this scale),
  milliseconds at any np. Verified compatible with `MPIDiagonalSOE::setSize`, which
  discovers neighbors by intersecting *global equation numbers* (`MPIDiagonalSOE.cpp:271`)
  — any consistent global bijection works. T2 changes the numbering (it abandons global
  RCM), so it is **opt-in per analysis**, not a silent default — see §Fence.

## Why — the incident, and the arithmetic that closes

**Observed** (slurm 144090, np240 / 15 nodes, `-deep -perStep -memory` profiler): one step
at 58,876 s wall, synchronized to the second across all 240 ranks; CPU/wall = 0.999
(MPI busy-poll); the timed sub-phases (`newStep`+`solveCurrentStep`+`commit`) covered only
2.6% of it — the cost sat in the *untimed* region of the step loop, which is exactly where
`domainChanged()` → `numberDOF()` runs on step 1.

**Predicted from source.** `mergeSubGraph` dedups each incoming vertex against the global
ref list with `vertexRefs.getLocation(ref)` (`ParallelNumberer.cpp:381`); `ID::getLocation`
is a linear scan (`ID.cpp:219-229`). ~90% of a 3D partition's vertices are interior — *not
yet present* — and a linear scan proves absence only by walking the **full** current array.
With 240 subgraphs of ~80 k vertices merged sequentially: Σ ≈ V_sub² × (239·240/2) ≈
**1.83×10¹⁴ comparisons**. An early-exit compare loop does not auto-vectorize and retires
~1 element/cycle ≈ 3×10⁹ cmp/s (compute-limited: the hardware prefetcher keeps a
sequential 4 B-stride scan fed — DRAM bandwidth is *not* the binding constraint).
**1.83×10¹⁴ / 3×10⁹ ≈ 61,000 s predicted vs 58,876 s observed.** Two independent
derivations (owner + Fable cross-review) land within 5% of the measurement. No second
mechanism is required.

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
retro-fit stays consistent. This is corroboration, not proof; G0/G3 are the proof.)

**The irony that sharpens the fix:** this model runs `system MPIDiagonal` under explicit
central difference. The LHS is diagonal — no factorization, no fill — so the **RCM
ordering the entire gather-merge exists to feed buys nothing** on this lane. The T2
endgame is to not compute it at all where it is provably unused.

## Where

| what | file:line |
|---|---|
| the linear scan | `SRC/matrix/ID.cpp:219-229` (`getLocation`); binary-search sibling `getLocationOrdered` at `:232-253` (unusable here — see Decision log) |
| vertex-dedup pass (pass 1) | `SRC/analysis/numberer/ParallelNumberer.cpp:381` (`vertexRefs.getLocation`) |
| adjacency pass (pass 2) | `ParallelNumberer.cpp:406` (1/vertex) + `:413` (1/directed adjacency entry) vs `theSubdomainMap` — **2-4×10¹³ comparisons on its own** (hours); any fix that only touches pass 1 leaves this |
| `ParallelPlain` extra quadratic | `ParallelNumberer.cpp:222` (`ID theOrderedRefs(V)` ⇒ sz=V from call one) + `:241`,`:249` (`getLocation` per vertex ⇒ V² ≈ 3.7×10¹⁴ **additional**) |
| ID append growth (verified benign) | `ID.cpp:335` (doubling ⇒ amortized O(1); the appends are NOT the problem) |
| post-T0 residual | `Graph.cpp` `addEdge` (2 `std::map` finds per directed entry, ~5.2×10⁸ at this scale) + `RCM.cpp` BFS `getVertexPtr` per adjacency entry ⇒ 20-45 min |
| the T1 cure, already in-tree | `Graph.cpp:223-276` (`startAddEdge`/`addEdgeFast` dense tag-indexed vector); merged tags contiguous via `getFreeTag` (`Graph.cpp:310`) |
| T2 compatibility | `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp:271` (neighbor discovery by global-eq-number intersection — numbering-scheme-agnostic) |
| registration (both verbs → same engine) | `SRC/tcl/commands.cpp:3919-3929`, `SRC/interpreter/OpenSeesCommands.cpp:1575-1581` |

## How — the class, and the T0 design (the stall fix)

**`LadrunoParallelNumberer`** (`SRC/analysis/numberer/LadrunoParallelNumberer.{h,cpp}`):
subclass of `ParallelNumberer`; overrides both `numberDOF` variants; inherits the
channel/processID plumbing via the ledgered promotion; keeps RCM (or none) as the
graph numberer on top — same algorithm, sane data structures. Registration: additive
`numberer LadrunoParallelRCM` (+ plain variant) branches in `commands.cpp` /
`OpenSeesCommands.cpp`; additive FEM_ObjectBroker case for the new classTag. apeGmsh
side: one typed primitive (`ops.numberer.LadrunoParallelRCM()`), fork-gated at run like
every Ladruno surface; the plane-wave family decks switch to it once G1 passes.

One lookup structure, hoisted to `numberDOF` and threaded through all merge calls (NOT
rebuilt per subgraph):

- **Primary:** `std::vector<int32_t> refToLoc` sized maxTag+1, sentinel −1 — node tags from
  gmsh/apeGmsh are dense 1..N, so this is O(1), 76 MB at 19.18 M, perfectly prefetchable.
  Guard: if maxTag > κ·numVertex (sparse/hand-authored tag spaces, κ≈4), fall back to
  `std::unordered_map<int,int>` with `reserve(V)` (~1-1.5 GB, still O(V) total).
- **Pass 1** (`:381`): replace `vertexRefs.getLocation(ref)` with the table probe; on
  miss, append exactly as today (same `getFreeTag()` sequence, same append order) and
  record `refToLoc[ref] = loc`.
- **Pass 2** (`:406`,`:413`): during pass 1 also fill a per-subgraph
  `subTagToMerged` map (V_sub entries, cleared per subgraph); pass 2 becomes O(E_sub).
- **`ParallelPlain` branch** (`:241`,`:249`): replace `theOrderedRefs->getLocation` with a
  seen-bitmap keyed by merged tag. Also fixes the latent tag-0 bug: the ID is zero-filled
  and `START_VERTEX_NUM = 0`, so vertex tag 0 currently reads "already ordered" against
  untouched zeros.

**Ordering stability is the load-bearing correctness claim:** the structures above are
pure lookup accelerators — iteration order over subgraph vertices, `getFreeTag` tag
assignment, and append order are untouched ⇒ identical merged graph ⇒ identical RCM input
⇒ **bit-identical equation numbering**. G1 enforces this, not the argument.

## No deck-level escape (all routes verified in source)

- **Every numberer verb hits the merge in MP builds.** Under `_PARALLEL_PROCESSING`,
  even `numberer Plain` and `numberer RCM` construct `ParallelNumberer`
  (`commands.cpp:3894-3898`). There is no verb in OpenSeesMP that avoids
  `mergeSubGraph`.
- **`ParallelPlain`** runs the *same* merge, then adds the `:222`/`:241` full-width scans:
  ≈ 3.7×10¹⁴ extra comparisons ⇒ ~2-3 **days** at this scale. Strictly worse.
- **Per-rank `Plain`/`RCM`** is a **wrong-answer generator** with any distributed SOE:
  `MPIDiagonalSOE::setSize` identifies shared DOFs by *equality of global equation
  numbers across ranks* — independent per-rank 0..n numbering makes every rank spuriously
  "share" its low-numbered equations and silently sum physically unrelated DOFs.

The code fix is mandatory. It is also stock-upstream code ⇒ **upstreamable** (T0/T1 are
behavior-identical, the cleanest possible upstream PR shape).

## Gates

| Gate | Delivers | Pass criterion | Status |
|---|---|---|---|
| **G0** — scaling sweep (pre-fix) | measured `dc.numberDOF` exponent vs N at fixed np | super-linear (~N²) growth across 0.25/1.0/2.0 M-hex rungs, np8, `numberer_sweep/` harness | **in flight** (local, instrumented binary e0113e53a) |
| **G1** — numbering byte-identity | the new class changes no numbers | DOF-ID dump diff, `LadrunoParallelRCM` vs stock `ParallelRCM`, same partition, 0.25 M rung, np2 + np8: **bit-identical**; full existing MP battery green | pending T0 |
| **G2** — the fix, measured | T0(+T1) timing on the same rungs | `dc.numberDOF` drops from the G0 curve to O(V) (T0) / low-minutes-at-19M-extrapolation (T1); exponent ≈ 1 | pending T0 |
| **G3** — cluster confirmation | the production-scale kill | 18.6 M deck re-run on Esmeralda with the patched binary: first step ≤ minutes, total ≈ 25 min; per-step wall column confirms *which* step carried setup | pending G1/G2 + cluster rebuild |
| **G4** (T2 only) | distributed numbering correctness | same-physics gate vs T0 numbering (response identical to solver tolerance, NOT bit-identical — ordering differs); np-sweep shows O(V/P) | deferred until T2 is scheduled |

## Instrumentation (shipped with this investigation, branch `perf/step-stall-instrumentation`)

Commit `e0113e53a` — two profiler gaps closed, prerequisites for G0-G3:
1. the `-perStep` wall column now records real per-step wall (was: `0.0f` pushed,
   `phases` never seeded ⇒ zero-width column ⇒ the old profiles physically cannot say
   *which* step stalled);
2. `domainChanged` is bracketed (`dc.handle` / `dc.numberDOF` / `dc.setSize`) inside the
   `step` scope, so setup cost is attributed instead of landing in the untimed remainder.

## Fence — what this ADR is and is NOT

- **Not the element-kernel lane.** Steady state was already at the ~18 µs/hex·step floor;
  per-step optimization is [[68_ladruno_state_determination_perf_adr]]. This ADR is the
  **once-per-domain-change** axis only.
- **Not a SIMD project.** Vectorizing the scan (AVX2 compare+movemask) buys a realistic
  2-3× (16 h → 6-8 h) against ~10⁵× from T0. Rejected — see Decision log.
- **Not a general HPC stack adoption.** ADR-40's anti-goal stands *except* where this
  incident measured the bottleneck: the numbering path. ParMETIS/PT-Scotch ordering is
  pulled in by T2 **only if** a distributed *fill-reducing* ordering is ever actually
  needed (implicit + distributed direct solve lane) — and even there, MUMPS/SuperLU_DIST
  prefer their own internal orderings, so the likely T2 shape never needs it.
- **RCM survives.** T0/T1 keep global RCM byte-identical for every lane that factorizes
  (banded/profile serial solvers under OpenSeesSP-style use). T2's skip-the-ordering mode
  is scoped to diagonal/explicit configurations and is opt-in.
- **The 2.5 GbE fabric hypothesis is retired only if G3 passes.** Until the cluster
  confirmation run shows the stall gone, the interconnect explanation (Obsidian
  Interconnect-vs-RAM note) stays formally alive. G0 exponent + the closing arithmetic
  make it unlikely; G3 makes it dead.

## Decision log

- **Direct-address vector over `unordered_map` as primary** — tags are dense 1..N in the
  apeGmsh pipeline; 76 MB flat array beats ~1-1.5 GB of hash nodes and malloc churn.
  Fallback keeps correctness for sparse tag spaces. (Fable review, folded.)
- **`getLocationOrdered` rejected** — `vertexRefs` is append-ordered and its order is
  load-bearing (positional correspondence with `vertexTags`); sorting it breaks the
  mapping, and a separate sorted index is strictly more code than the table for no gain.
- **SIMD rejected as the stall fix** — the scan is compute-limited at ~1 elem/cycle;
  vectorization lifts it 2-3× where T0 lifts it ~10⁵×. SIMD effort belongs on the
  per-step kernels (ADR-68 lane), which run 1000× per analysis; the numberer runs once.
- **Fix must cover pass 2 and the `ParallelPlain` branch, not just pass 1** — pass 2 is
  2-4×10¹³ comparisons (hours) on its own; `ParallelPlain` is ~2× worse than the merge
  itself. Under-scoping T0 to the single `:381` line would leave a still-unusable path.
  (Fable review — the original owner proposal under-scoped exactly this way.)
- **T1 uses the existing `startAddEdge`/`addEdgeFast`** rather than a new structure — the
  dense path already in `Graph.cpp` fits because merged tags are `getFreeTag`-contiguous;
  reuse > reinvent (ADR-40 principle).
- **T2 is opt-in, not a default swap** — it abandons global RCM, so numbering (and any
  ordering-sensitive downstream) changes. Diagonal/explicit lanes lose nothing; every
  other lane keeps today's path until measured.
- **Fork-native class over in-place edit** (owner, 2026-07-21) — `LadrunoParallelNumberer`
  in the ≥33000 band rather than editing `ParallelNumberer::numberDOF` in place. Honors
  the fork's additive contract (upstream classes behavior-untouched; only a ledgered
  visibility promotion + additive registration), gives the deck an explicit opt-in verb,
  and keeps stock `ParallelRCM` alive as the G1 byte-identity reference *inside the same
  binary*. The upstream PR is later extracted as the in-place form of T0/T1.

## Risks / open questions

- **maxTag density assumption** (T0 primary path): compose/renumber flows or hand decks
  could present sparse tags; the κ-guard + hash fallback covers it, but the guard needs a
  unit test with a deliberately sparse tag space.
- **Rank-0 memory after T0**: the merged graph (~5 GB at 19 M vertices: Vertex objects +
  map nodes + adjacency IDs) still lives on rank 0 until T2. Fine at 60 GB nodes;
  worth a printed diagnostic (vertices + est. bytes) when V > 5 M.
- **Wire cost of the gather** stays O(E_global) data volume even after T0 — acceptable
  (the 18.6 M graph is ~GB-scale over the fabric, minutes not hours), eliminated by T2.
- **G1 dump mechanics**: needs a small DOF-ID dump hook (rank, nodeTag, eqn numbers) —
  decide whether that ships as a debug verb or a test-only fixture.
- **Upstream PR**: T0/T1 are byte-identical and upstream-shaped; coordinate with the
  fork's upstream-divergence policy before opening.

## Implementation log

- **2026-07-21** — Root cause identified from source (owner) after the 144090 profiler
  showed 58,876 s in the untimed step region. Closing arithmetic: 1.83×10¹⁴ comparisons
  @ ~3×10⁹/s ≈ 61,000 s vs 58,876 s observed.
- **2026-07-21** — Profiler gaps closed: commit `e0113e53a`
  (`perf/step-stall-instrumentation`) — per-step wall column + `dc.*` attribution.
  Verified: 10/10 standalone harness on `recordStep` semantics; both TUs syntax-clean;
  scope strings present in the rebuilt Windows binary.
- **2026-07-21** — Adversarial cross-review (Fable agent): mechanism + magnitude
  **confirmed independently (±5%)**; four findings folded — (1) compute-limited not
  bandwidth-limited ⇒ SIMD worth 2-3× not <2×; (2) post-T0 residual is 20-45 min of
  `std::map` finds ⇒ T1 added; (3) original fix under-scoped (pass 2 + `ParallelPlain`
  branch added to T0); (4) per-rank `Plain` confirmed *incorrect* (not merely slow) with
  distributed SOEs ⇒ "no deck-level escape" is a hard claim, verified at
  `MPIDiagonalSOE.cpp:271`.
- **2026-07-21** — G0 sweep launched locally (np8; 0.25/1.0/2.0 M-hex rungs; instrumented
  binary; harness `Dropbox/UANDES EC/Cluster Max Models/numberer_sweep/`). Result: *fill
  in on completion.*
- **Cluster note**: the Esmeralda tree (`ladruno-p5-build` @ #580) carries the same
  unfixed code — all three touched files byte-identical to the patch base — so the
  instrumentation + T0 rebase cleanly there. The live cluster binary is Jul-18 (not the
  Jun-11 build the Obsidian note records).
