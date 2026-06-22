# ADR-39 P2.5 — bucket-sort broad phase (drop-in for the brute-force pairing)

> Replaces the handler's O(nSlave·nSeg) brute-force (slave, master-segment) pairing
> with the §26.11 bucket sort, validated oracle-first in
> `contact_prototypes/proto_bucket_sort.py` (P0a, 4/4: superset contract, jittered
> mesh, runaway guard, complexity). Gate = **verify == brute force** (identical
> contact result, fewer adapters). Parent ADR `39_..._adr.md`; driver `_adr39_loop_state.md`.

## What it changes (and what it does NOT)

TODAY (P2b-1/2b), `LadrunoContactHandler::handle()` builds ONE `LadrunoContactFE`
adapter per (slave node × master segment) pair — every pair, brute force. Each
adapter self-deactivates when its slave is not penetrating its segment. Correct but
O(nSlave·nSeg) adapters (most permanently inactive).

P2.5 inserts a **broad phase**: build a spatial bucket grid over the master segments,
then for each slave node emit adapters ONLY for the segments in its 27-neighbour
bucket cell. The narrow phase (projection/gap/penalty in the FE) is UNCHANGED. The
contact RESULT is identical for any pair the broad phase keeps; the gate proves the
kept set is a **superset of every geometrically-near pair**, so nothing active is lost.

NON-goal (scoping, documented): the grid is built at `handle()` time (a domainChanged
event), NOT every step. For STATIC and SMALL-MOTION dynamic contact (the validated
regime + the gate meshes) the handle()-time candidate set covers every pair that
activates → bitwise-identical to brute force. **Large sliding** (a slave migrating
more than the search band from its handle()-time cell) needs the broad phase
RE-EMITTED — that is the epoch-re-emit follow-on (ADR §active-set churn; LS-DYNA
re-sorts every 10–15 cycles). P2.5 ships the sort; epoch re-emit is the next rung
(or rides P3). The handler already rebuilds all adapters on every domainChanged, so
a manual re-`analyze` after remeshing/large motion re-sorts for free.

## Algorithm (from the P0a oracle — transcribe verbatim)

`proto_bucket_sort.py::BucketSort`:
1. **cell size** `LMAX = lmaxFrac · median(segment longest-diagonal)` (§26.11.1).
   lmaxFrac default 1.0 (1.5 for distorted meshes — expose as the cell knob).
2. **bbox** from segment centroids; optional **runaway-node percentile clip**
   (e.g. 1–99%) so one fly-off node can't inflate the bbox and collapse all real
   segments into one cell (P0a Test 3 proved this — un-clipped → real surface
   collapses to ~1 bucket; clipped → recovers).
3. `NX=(xmax−xmin)/LMAX+1`, etc.; **cap** `NX·NY·NZ ≤ min(nSeg, 5000)` (26.37),
   shrinking the largest dim until under cap.
4. **segment-SPAN registration** (LS-DYNA type-13, eq 26.44): register each segment
   in EVERY bucket its CORNER pointers span, not just its centroid bucket. Centroid-
   only (type-4) is deficient (§26.8.1) — a large/poorly-shaped segment's nearest
   centroid ≠ nearest segment → true-near pairs missed (P0a caught 30% misses).
5. **candidates(slavePt)** = union of segments registered in the 27 neighbour buckets
   `[PX−1..PX+1]×…` (26.42), de-duplicated.

**Superset contract (the gate):** for search radius `r ≤ cell` with ±1-neighbour
search, candidates ⊇ {segments within r of the slave}. The realistic contact band
(gap tol + per-epoch slip) is a fraction of a cell, so the in-guarantee contract is
`r ≤ 0.5·cell` (P0a uses r=0.5·cell; r==cell is the boundary edge, not guaranteed).

**Why the C++ cell = `lmaxFrac·median(seg diag)` (no explicit `r_search` term) is
sufficient** (P2.5 code-gate MINOR-1, reconciled): for PENETRATING contact the slave
lies on the master surface (gap ≈ 0), and any point of a bilinear-quad / tri segment
is a convex combination of its corner nodes ⇒ the slave is inside the segment's
corner-AABB ⇒ it is in a cell where SPAN-registration already placed the segment (the
±1 neighbour is pure extra margin). So `r_search ≈ 0 ≪ cell` automatically and the
superset holds without a search-band term. The only thing this does NOT cover is a
slave APPROACHING from more than ~1 cell away pre-penetration between two handle()s —
that is exactly the large-sliding / epoch-re-emit scope limit (below), not a static-
contact gap. (Gate verified 0 misses over >250k probes incl. skewed/ramped/mixed
meshes.) If pre-penetration approach-band tracking is ever needed, re-introduce
`cell = max(median diag, r_search)·margin` then.

## C++ — NEW header-only `SRC/domain/contact/LadrunoContactBucketSort.h`

Mirror the `LadrunoContactKernel.h` discipline: self-contained, OpenSees-free
(raw `double[]` + STL only), so the handler (OPS_Analysis) includes it with no link
inversion. A small class (grid state is non-trivial) rather than pure fns:

```
namespace LadrunoContactBucketSort {
  class Grid {
   public:
    // segCoords: flat nSeg*nps*3 (reference coords); lmaxFrac default 1.0;
    // clipPct: 0 disables the runaway guard, else e.g. 1.0 (clip to [1%,99%]).
    Grid(int nSeg, int nps, const double *segCoords,
         double lmaxFrac = 1.0, double clipPct = 1.0);
    // fill out[] (caller-sized >= nSeg) with the de-duplicated candidate segment
    // indices for slavePt; returns the count. De-dup via an internal stamp array.
    int candidates(const double slavePt[3], int *out) const;
    int numSegments(void) const { return nSeg_; }
    bool runawayGuardFired(void) const { return guardFired_; }
   private:
    int nSeg_, nps_, NX_, NY_, NZ_;
    double xmin_,xmax_,ymin_,ymax_,zmin_,zmax_,lmax_;
    bool guardFired_;
    std::vector<std::vector<int>> grid_;          // size NX*NY*NZ, linear bucket id
    mutable std::vector<int> stamp_; mutable int stampTick_;  // O(1) candidate de-dup
    int px(double x,double lo,double hi,int N) const;
    int linId(int px,int py,int pz) const { return (px-1)+(py-1)*NX_+(pz-1)*NX_*NY_; }
  };
}
```

Notes:
- Linear bucket id (capped ≤ min(nSeg,5000)) → `std::vector<std::vector<int>>` not a
  hash map (simpler, cache-friendly, the cap bounds memory).
- De-dup with a `stamp_` array + monotonically-increasing `stampTick_` per
  candidates() call (no per-call alloc/clear): a seg id is added to `out` only if
  `stamp_[id] != tick`, then stamped. Standard sparse-set de-dup.
- Reference coords (undeformed) are sufficient for the broad phase under the
  small-motion contract; the search band absorbs the displacement. (A future
  epoch re-emit would rebuild from current coords.)

## Handler integration — `LadrunoContactHandler::handle()` segment loop

Per contact definition (a MASTER_SEGMENTS surface vs a SLAVE_NODES surface):
1. Gather the master segment reference coords into a flat `nSeg*nps*3` buffer.
2. Build a `Grid` over them.
3. For each slave node: `grid.candidates(slaveRefCoords, cand)`; loop ONLY over the
   `cand[0..nCand)` segments instead of all `nSeg`. Inside, the existing per-pair
   logic is UNCHANGED (segNodes gather, ndf/degenerate guards, orientDir, P2b-2b
   `-kn auto` resolve, build adapter).

The rigid-plane (P2a) path is untouched (no segments → no broad phase).

## Cell/search knob (parser, optional)

Default lmaxFrac=1.0 + clip 1% works for the gate meshes. Optionally expose
`contact … -cell <lmaxFrac>` later; v1 P2.5 can hard-default (no parser change →
smaller surface, faster gate). DECIDE while coding: keep the default-only path
for P2.5 (zero parser/vanilla churn) unless a test needs the knob.

## Tests — `tests/test_adr39_contact_p25_bucketsort.py`

1. **verify == brute force (THE gate):** a multi-segment master (e.g. a 2×2 or 3×3
   grid of quad faces = 4/9 segments, from a meshed LadrunoBrick block top or a raw
   node grid) + several slaves, some penetrating different segments. Run the SAME
   model twice — once on the bucket-sort handler (default), once forcing brute force
   — assert the per-node contact forces / reactions are bitwise (rel-1e-12) identical.
   (Mechanism for "force brute force": a hidden `-cell` huge value → 1 bucket → all
   segments candidate = brute force; OR an env/static toggle. DECIDE while coding —
   simplest is a cell so large the whole surface is one bucket.)
2. **superset / no missed contact:** a slave directly over an INTERIOR segment of a
   3×3 master grid penetrates → contact force is correct (the right segment was a
   candidate), matching the analytic P/kn.
3. **complexity:** on a larger grid (e.g. 6×6=36 segments) assert the emitted adapter
   count << nSlave·nSeg (query a count via `ladrunoContactInfo` or a new probe; or
   assert the analysis still converges + correct while the candidate set is small).
4. **regression:** the existing single-segment P2b tests still pass (1 segment → grid
   trivially returns it).

## Files (P2.5)
- NEW `SRC/domain/contact/LadrunoContactBucketSort.h` (header-only; stamp_headers + ledger).
- `SRC/analysis/handler/LadrunoContactHandler.cpp` — broad-phase in the segment loop.
- (maybe) parser `-cell` knob — only if a test needs it (prefer none for P2.5).
- `tests/test_adr39_contact_p25_bucketsort.py`.

## Open items to resolve while coding
- Q1: how to "force brute force" in the equivalence test — a huge cell (1 bucket) is
  the cleanest (no toggle plumbing); confirm 1-bucket == all-segments-candidate.
- Q2: dedup stamp array lifetime (mutable members on a const method — fine, it's a
  cache; or make candidates() non-const). Keep it const + mutable.
- Q3: degenerate/collinear master segment diagonal → median diag still > 0 if any
  segment is non-degenerate; guard lmax_ ≥ tiny floor (mirror the P0a 1e-12 floor).
- Q4: a slave exactly on a bucket boundary — px() clamps to [1,N]; ±1 neighbour
  covers it. Confirmed by the P0a 2000-probe superset test (0 misses).
