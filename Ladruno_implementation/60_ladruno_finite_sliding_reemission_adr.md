# ADR-60 — Finite-Sliding NTS Contact Re-Emission (deformed-configuration broad phase)

- **Status:** DESIGN — pre-code adversarial gate COMPLETE 2026-06-29 (5-lens; dispositions folded below). Ready for P0.
- **Owner:** Mora Bowen · Palacios · Abell · Guppi
- **Priority:** high (silent-correctness bug — NTS contact force vanishes with no warning under large sliding)
- **Parents / prior art:** [[39_ladruno_contact_domain_adr]] (NTS broad+narrow phase, P2.5 bucket sort),
  [[41_ladruno_mortar_alm_contact_adr]] (mortar = brute force + node-keyed λ_N), [[57_ladruno_edge_edge_contact_adr]],
  [[48_ladruno_contact_capstone_adr]] (status-of-record).
- **Successor / fence:** [[55_ladruno_contact_runtime_discovery_adr]] (open-universe discovery — this ADR is its
  CLOSED-universe NTS precursor), [[51_ladruno_element_removal_adr]] (newly-exposed faces; the segIndex/removal seam).
- **Class tags:** none new. Behavior of `LadrunoContactHandler` (`HANDLER_TAG_LadrunoContactHandler = 33002`) +
  `LadrunoContactDomain`.

---

## What

Make **node-to-segment (NTS)** contact survive **finite sliding** within declared surfaces: when a slave's
**deformed** position migrates off the candidate segments the NTS bucket sort picked at `handle()` time, **re-emit**
the candidate set from the current (committed) configuration — automatically, at a committed config, with friction
re-engaged cleanly at the crossing.

**Scope narrowed by the gate (2026-06-29):** the silent pass-through is **NTS-only**. The mortar lane is **brute
force** — every (slave facet, master facet) pair is enumerated and the deformed-config clip activates overlaps
(`LadrunoContactHandler.cpp:637-644`), so mortar already handles finite sliding (at O(nSf·nMf) cost; a slave-aware
"mortar P2.5" broad phase is a separate *cost* optimization, not a correctness gap). Edge-edge rides the mortar
brute-force enumeration (`:755`). **So ADR-60 touches the NTS lane only.** One pre-existing mortar friction issue
(stale `gT0` engagement origin when a slave's dominant master facet changes) is noted but is NOT introduced by
re-emission and is out of scope (see Q-MORTAR-GT0).

OUT: runtime self-contact / free-body discovery / runtime surface registration / size-dispersity BVH (→ ADR-55);
newly-exposed faces after element removal (→ ADR-51); parallel/MP re-emit (refused, see D7); mortar re-emission
(unneeded).

---

## Why

The NTS broad phase (`LadrunoContactBucketSort`, ADR-39 P2.5) builds its 27-bucket grid **once per `handle()`** from
**reference** coords — segment nodes via `mn->getCrds()` (`LadrunoContactHandler.cpp:510`), slave query via
`sn->getCrds()` (`:533`). `Node::getCrds()` never changes during analysis, so the candidate `(slave, segment)` set is
**frozen undeformed**. A slave that slides more than the search band (`~cellFrac × median seg diagonal`) off its
reference cell lands on a **non-candidate** segment → no adapter is built → the deformed-aware narrow phase
(`segmentActive`) never sees it → **silent pass-through**: contact force drops to zero, no warning, no crash.

The shipped comment "a manual re-analyze re-sorts for free" (`LadrunoContactBucketSort.h:40`) is **wrong**:
re-`handle()` re-reads `getCrds()` → identical candidate set. This is the rung ADR-39 P2.5 deferred (`:36-40`).
→ `LEDGER_quirks.md`.

---

## Where

- **New** `SRC/domain/contact/LadrunoContactReemit.h` — header-only, OpenSees-free (mirrors `LadrunoContactKernel.h`):
  the migration metric only (`max |x_committed − anchor|` vs `f·band`, with the reference-median band, D6/Finding-6).
  No slip-transfer math (dropped — see D4). Oracle-testable in isolation.
- **Modified** `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — per-contact re-emit anchor store + a **persistent
  per-(contactTag,slaveTag) orientDir** captured at first handle (D3/Finding-4); `needsResort(Domain*)` with the
  augmentation + parallel + empty-anchor guards; a **surface-membership epoch** per contact (H6); config fields with
  **default member initializers** (Finding-4-byte). NTS friction uses plain fresh-slot re-engagement — no transfer code.
- **Modified** `SRC/analysis/handler/LadrunoContactHandler.cpp` — NTS loop only: when re-emit ON, feed the grid
  `getCrds()+getDisp()` (committed); register anchors + the held orientDir; disable the runaway clip on the deformed
  feed and surface `runawayGuardFired()` (Finding-3); exclude self-contact `(slave ∈ segment.nodes)` candidates
  (Finding-4b); detect slide-off-surface (Finding-4a).
- **Modified (1 line, guarded)** `SRC/domain/domain/Domain.cpp` — extend the existing `// Ladruno (ADR-39)` `commit()`
  hook (`~:2199`): `if (theContactDomain != 0 && !contactAugmenting && theContactDomain->needsResort(this)) this->domainChange();`.
  Placed **after** the node-commit loop (`:2185-2187`, ordering locked, Finding-2) and gated on `!contactAugmenting`
  (Finding-GA1). `domainChange()` is the public flag-raiser (`Domain.h:224`). **No new vanilla file.**
- **Reference sources:** LS-DYNA Theory §26.11 (pp.545–550, deformed bbox, 5–15-cycle re-sort + 2–3 stored segments +
  surrounding-segment local tracking); Abaqus TG §5.1.1 (small-sliding = our shipped mode) vs §5.1.2 (finite-sliding =
  re-projection + slide-line smoothing).

---

## How

### Mechanism (cache-free, rides a tested rail)

1. **Anchor + orientDir registration (`handle()`, ON only).** After the NTS grid is built, the handler hands the engine,
   per contact: `{slaveTag → anchor[3]}` (the coords it fed the grid) and the search `band` (reference-median based,
   D6). On the **first** handle the engine also captures, per (contactTag, slaveTag), the **orientDir** from the
   reference config (slave on the correct side) and holds it verbatim thereafter (Finding-4 / BLOCKER-6).
2. **Migration watch (`commit()`).** `Domain::commit()` already calls `theContactDomain->commit()` after committing all
   node states. We add the one guarded line above. `needsResort(Domain*)` first-returns `false` if anchors are empty,
   if `contactAugmenting`, or under a partitioned host; else loops anchors, looks up each slave's CURRENT committed
   position **by tag through the passed `Domain*`** (no cached pointer), and returns `dmax > f·band` OR a cycle-counter
   trip — with **hysteresis** (re-arm only after `dmax < 0.5 f·band`) and a **min-steps floor (default 10, ≈ LS-DYNA
   cadence)** to defeat threshold chatter (Finding-5).
3. **Re-handle (next substep).** `domainChange()` raises `hasDomainChangedFlag`; the driver re-handles at the top of
   the next step (`DirectIntegrationAnalysis::analyzeStep:206`, `StaticAnalysis::analyze:150` inside the loop;
   `VariableTimeStep`/`PFEM`/DomainDecomposition all poll — Finding-3 swept). This `commit()`→`domainChange()`→re-handle
   path is the **already-exercised load-balancing rail** (`StaticAnalysis.cpp:147-149`).
4. **Re-emit (`handle()`).** NTS grid rebuilt from `getCrds()+getDisp()`; new pairs reflect the deformed config;
   anchors reset; friction re-engages on fresh slots (D4).

### Decisions

- **D1 — Trigger/cadence:** migration metric (default `f = 0.5 × band`) + hysteresis + min-steps floor (default 10),
  watched in `commit()`. Hybrid `-resortEvery N` cycle override (= LS-DYNA BSORT; the natural unit under explicit).
- **D2 — Plumbing:** one guarded line in the existing ADR-39 `Domain::commit()` hook; `needsResort(Domain*)`.
- **D3 — Coordinate source:** deformed `getCrds()+getDisp()` when ON, reference `getCrds()` when OFF (the OFF branch is
  the **byte-verbatim** current code, not a `+0.0` ternary — Finding-3). `orientDir` **persisted** per (contactTag,
  slaveTag), captured at first reference handle, held across re-emits (Finding-4). auto-`kn` stays reference
  (config-independent, verified Finding-5).
- **D4 — History at a crossing: RE-ENGAGE, do NOT transfer.** (Revised by the gate — Reviewer-3 HIGH.) The friction
  return map consumes a **displacement-anchored** slip `gTeff = gTvec − gT0` (`LadrunoContactFE.cpp:227-234,655`;
  `LadrunoFrictionKernel.h:102`), NOT a position-space in-plane vector. When a slave re-pairs onto a new `segIndex`,
  `getOrCreateFrictionState` returns a **fresh zeroed slot** (`engaged=false`) which self-captures `gT0 = current
  gTvec`, `gpT = 0` ⇒ `gTeff − gpT = 0` ⇒ **zero stick force at the crossing instant ⇒ already traction-continuous.**
  Rotating `gpT` (the dropped D4-v1) would have injected an `O(traction)` spurious force and is *unsupported by the
  math*; the old slot is GC'd. **Sustained sticking-drag continuity across a crossing is the job of P4 (patch
  adapter)**, exactly as LS-DYNA (stored-neighbor segments) and Abaqus (slide line) do it. Mortar **λ_N** (scalar,
  node-keyed, frame-free) survives finite sliding for free; mortar **friction** `gT0` staleness is a pre-existing,
  out-of-scope issue (Q-MORTAR-GT0).
- **D5 — Scope:** NTS only (see What).
- **D6 — Cost & band:** re-emit rides a full `domainChanged` (renumber+graph+SOE), frequency `≈ 2D/h`. The search
  `band` used by the trigger is keyed off the **reference** median diagonal (deformation-invariant, cheap O(nSeg) from
  `getCrds()`), so it cannot drift away from the grid the next re-emit builds (Finding-6). Default OFF ⇒ zero churn.
- **D7 — Parallel:** re-emit is **serial-only**. `needsResort()` returns `false` + skips anchor registration under a
  `Subdomain`/partitioned host (one-time warning); a slave whose deformed position would re-pair across a partition
  boundary inherits ADR-39's partition-crossing refusal (Finding-MP). MP re-emit → ADR-55 / ADR-51-H7.

### The patch / local-tracking alternative (P4)

Re-emission is unavoidable for unbounded sliding (the `FE_Element`-only-in-`handle()` wall). LS-DYNA and Abaqus both
pair global re-search with **local tracking**. P4 builds each NTS slave adapter connected to its neighborhood **patch**
of segments (+margin); the narrow phase picks the active segment within the patch each step ⇒ within-patch sliding
needs no `domainChanged` **and is friction-continuous by construction** (one adapter owns the patch — the real fix for
sustained-drag continuity that D4 deliberately does NOT attempt). Re-emit then fires only on patch-exit. Staged second
so the connectivity refactor is decoupled from the core fix.

---

## Design-gate BLOCKERs (status after the gate)

1. **BLOCKER-REENTRY** — ✅ RESOLVED. `commit()`→`domainChange()`→re-handle is consumed at each substep and is the
   load-balancing rail. Residual (held-load Uzawa) handled by GA-1.
2. **BLOCKER-IDENTITY** — strengthened (Finding-2/4/5): OFF must be bitwise identical AND assert (a) zero extra
   `domainChange`/`currentGeoTag` advance, (b) byte-identical filled `segCoords`/`slavePt`, (c) empty anchor map, (d)
   `needsResort()` O(1)-false. Struct config fields get **default member initializers** so stack garbage cannot flip
   `enableReemit` true (Finding-4-byte). P0 gate.
3. **BLOCKER-GA1 (augmentation)** — ✅ promoted to P0. The line is gated `!contactAugmenting` AND the `-resortEvery`
   counter must NOT advance during augmentation (`Domain::isContactAugmenting()`, `Domain.h:251`).
4. **BLOCKER-ORIENTDIR (was -6)** — real. Persist per-(contactTag,slaveTag) orientDir from the first reference handle;
   never re-derive from a penetrating current config. P1 gate = convex-ridge crossing keeps a repulsive, non-flipped
   normal.
5. **BLOCKER-CLIP** — disable the percentile runaway clip on the deformed feed (else it collapses the grid and drops
   real pairs during instability); surface `runawayGuardFired()` as a one-time warning. P1.
6. **BLOCKER-SELF-CONTACT** — exclude `(slave ∈ segment.node-tags)` candidates at re-emit (a folded declared surface
   can pair a node with a segment containing it). 3-line filter. P1.
7. **BLOCKER-SLIDE-OFF** — distinguish "migrated to a new in-bounds segment" from "departed the surface"; on departure
   drop the friction slot (clean with re-engagement) and force → 0. P1 oracle includes a slide-off-the-edge case.
8. **BLOCKER-MEMBERSHIP (H6 seam)** — re-emit is **refused** on a contact whose surface membership changed since the
   last handle (a per-contact membership epoch/fingerprint; named error). This is the honest H6 disposition — replaces
   the struck "partial mitigation" claim. With re-engagement (no transfer) the friction-aliasing risk is already much
   reduced, but the refusal is the guarantee.

*Removed by the D4 re-engagement decision:* BLOCKER-TRANSFER-FRAME, the 180°/90° rotation degeneracy, and
BLOCKER-TRANSFER-AMBIGUITY (multi-slot) — there is no slip rotation to get wrong.

---

## Phased rollout (oracle-first)

| Phase | Deliverable | Validation oracle |
|---|---|---|
| **P0** | Plumbing, ZERO behavior change: `LadrunoContactReemit.h` (metric), `needsResort(Domain*)` (empty/augment/parallel guards, hysteresis, floor), anchor + persistent-orientDir store, membership epoch, struct default initializers, the 1-line guarded `Domain::commit()` extension, `-reemit off`/`-resortFrac`/`-resortEvery` parser. | numpy: migration-metric + hysteresis unit tests. C++ build-free: `needsResort` fires iff `dmax>f·band` and is inert when augmenting/parallel/empty. OpenSees: full contact battery byte-identical + currentGeoTag-sequence + segCoords-bytes identical (BLOCKER-IDENTITY). |
| **P1** | Deformed coord source (D3, `getDisp()`) + auto re-emit (D2) wired; **frictionless NTS**; clip-off + runaway warn; self-contact filter; slide-off detection; persistent orientDir. | numpy: slave sliding along piecewise-flat master — active segment ALWAYS in the re-emitted candidate set (**no-pass-through**, tight). OpenSees: sliding-block — (a) no-pass-through across ≥3 crossings; (b) **resultant** slave force continuous to a LOOSE tol that accounts for the documented faceted-normal jump (cite ADR-41 Q-NORMAL — NOT bit-continuity vs a smooth reference); (c) BLOCKER-ORIENTDIR convex-ridge; (d) BLOCKER-CLIP slave-stays-in-contact-while-sibling-diverges; (e) slide-off → force 0 cleanly. |
| **P2** | NTS **friction** under re-emit via fresh-slot re-engagement (no transfer code); GA-1 augmentation guard verified live. | numpy: traction at a crossing = 0 stick-force instant then re-accumulates (re-engagement continuity). OpenSees: dragged block w/ friction across crossings — traction has no spike at re-engagement; energy balance closed; held-load augmentation unaffected by the trigger. |
| **P3** | Metric tuning, `-resortEvery` hybrid, recorder-continuity doc, quirks + ledgers + banner, ADR-55/48 re-fence note. | OpenSees: `f`-sweep cost/accuracy; per-pair recorder discontinuity disclosed (summed `getNtsForce` continuous). |
| **P4** *(patch, optional)* | Patch / local-tracking adapter; within-patch crossing = narrow-phase only, friction-continuous. | numpy: patch-exit detection. OpenSees: same sliding-block, `domainChanged` count P4 ≪ P1; sustained-drag friction continuous (the continuity D4 does not attempt). |

---

## Fence (scope vs the collapse program)

- **ADR-55 (runtime discovery):** owns the OPEN universe (BVH for size dispersity, runtime registration, self-contact,
  free bodies). ADR-60 is the closed-universe NTS precursor; ADR-55 P1/P2 reuse its trigger verbatim and swap the grid
  for a BVH. ADR-55 P1's "friction history survives re-pairing" gate is the SAME deliverable as ADR-60 P2 — a one-line
  forward-fence is added to ADR-55 P1 so one becomes a no-op (P3 doc deliverable).
- **ADR-51 (element removal):** newly-exposed faces OUT. The H6 segIndex-aliasing seam is handled by BLOCKER-MEMBERSHIP
  (refuse re-emit on membership change), not by the (struck) "node-tag anchor partially mitigates" claim.
- **ADR-57 (edge-edge):** rides the mortar brute-force lane; no re-emission work.

---

## Ledger / classTag / banner bookkeeping

- **classTags:** none new (handler stays 33002).
- **`LEDGER_implementations.md`:** row — *Finite-sliding NTS contact re-emission*, files = `LadrunoContactReemit.h`
  (new) + `LadrunoContactDomain.{h,cpp}` + `LadrunoContactHandler.cpp`, status per phase, PR.
- **`LEDGER_vanilla_files.md`:** `Domain.cpp` — extend the existing ADR-39 `commit()` row with the ADR-60 guarded
  `needsResort()→domainChange()` line (one `// Ladruno (ADR-60)` comment). No other vanilla touch.
- **`LEDGER_quirks.md`:** (a) "re-handle re-sorts for free" myth; (b) mortar is brute-force so finite-sliding-correct
  already (NTS is the culled lane); (c) re-emit rides a full `domainChanged` (cost); (d) friction re-engages (not
  transfers) at a crossing — sustained-drag continuity needs P4; (e) `needsResort()` skips held-load augmentation
  commits (GA-1); (f) re-emit is serial-only; (g) the runaway clip is off on the deformed feed.
- **Banner:** `shipped`-gated line in `Ladruno_scripts/banner_features.txt` when P1 lands, then `patch_banner.py`.
- Fix the stale `LadrunoContactBucketSort.h:36-40` SCOPING comment once shipped.

---

## Risks / open questions

- **Q-MORTAR-GT0:** pre-existing — a mortar slave's `gT0` engagement origin is stale if its dominant master facet
  changes under sliding (frame-dependent, captured once). Not introduced by ADR-60; candidate follow-up (re-anchor
  `gT0` on detected dominant-facet change). Out of ADR-60 scope; logged.
- **Q-IMPLICIT-NEWTON:** does implicit Newton converge at a re-emit step (active-set jump + rebuilt SOE)? P2/P3 gate;
  if not, re-emit is declared explicit-first (matches ADR-39 posture / ADR-55 Q-EXPLICIT-CHURN).
- **Q-RESTART:** the anchor + orientDir store is **reconstructed at the next handle()** (re-derivable), so a
  database/sendSelf round-trip that drops it costs at most one conservative re-emit. Documented; no serialization added.
- **Q-PATCH-SCOPE:** confirm P4 belongs in ADR-60 vs its own ADR (connectivity refactor + own gate).

---

## Decision log

- 2026-06-29 — Planning sign-off: NEW ADR-60 (precursor to ADR-55); baseline re-emit P0–P3 + patch P4; `f=0.5` default
  + `-resortEvery N`.
- 2026-06-29 — Grounding: LS-DYNA §26.11 (deformed bbox, 5–15-cycle re-sort + local tracking) + Abaqus §5.1.1/§5.1.2 —
  patch is the standard finite-sliding algorithm, baseline the cheaper subset.
- 2026-06-29 — BLOCKER-REENTRY resolved (load-balancing rail).
- 2026-06-29 — **5-lens adversarial gate (below). Scope narrowed to NTS (mortar is brute-force); D4 changed from
  slip-rotation to fresh-slot re-engagement; orientDir persisted; getDisp() not getCommittedDisp(); GA1/MP/struct-defaults
  promoted to P0; H6 "partial mitigation" struck for an explicit membership-change refusal.**

## Adversarial review log (5-lens, 2026-06-29 — findings + dispositions)

**Lens A — trigger/control-flow.** (A1, major) GA-1 guard must wrap BOTH the metric and the `-resortEvery` counter and
sit below the `contactAugmenting` guard → **adopted, promoted to P0** (D2/BLOCKER-GA1). (A2, major) a rejected step can
leave a re-emit armed off a thrown-away config → mitigated: `commit()` runs only on accepted steps, and with D4
re-engagement the re-handle mutates no committed state → benign; documented. (A3) every shipped driver polls
`hasDomainChanged()` (incl. VariableTimeStep/PFEM/DomainDecomposition) → no missing-driver hole. (A4, major) trigger
fires uncoordinated on every MPI rank → **adopted: serial-only refusal (D7).** (A5, major) floor-default-1 + no
hysteresis is an explicit landmine → **adopted: floor default 10 + hysteresis (D1).**

**Lens B — byte-identity OFF.** (B1) reuse the `theContactDomain!=0` guard; gate the call so contact-present/reemit-OFF
pays nothing → adopted. (B2) `needsResort` early-returns on empty anchors; anchor store/GC strictly under
`if(enableReemit)` → adopted. (B3, med) keep the coord-source OFF branch byte-verbatim; gate must byte-diff
`segCoords` → adopted into BLOCKER-IDENTITY. (B4, HIGH) new struct fields are uninitialized POD → garbage could flip
`enableReemit` true → **adopted: default member initializers** (P0). (B5) one spurious re-emit is observable in
`currentGeoTag` → gate asserts the geo-tag sequence.

**Lens C — history transfer (the decisive lens).** (C1, HIGH) rotating `gpT` is inconsistent with the
displacement-anchored, `gT0`-framed return map and injects `O(traction)` spurious force; **fresh-slot re-engagement is
already traction-continuous** → **adopted: D4 rewritten to re-engagement; slip-rotation dropped.** (C2, med) the
minimal rotation is non-unique at 180° and the energy bound is force-scale not slip-scale → **moot (rotation removed).**
(C3, HIGH) **"mortar friction survives for free" is FALSE** — only scalar λ_N is frame-free; mortar `gT0`/`gpT`/`λ_T`
are frame-dependent + stale on facet change → **adopted: claim corrected; λ_N keep, mortar-friction `gT0` staleness is
pre-existing/out-of-scope (Q-MORTAR-GT0).** (C4) segIndex stable across re-emit (true) but not across removal → see H6;
recommend a runtime `mTags` fingerprint → **adopted as BLOCKER-MEMBERSHIP.** (C5) multi-slot ambiguity → **moot
(re-engagement; old slots GC'd).**

**Lens D — coords & geometry.** (D1, HIGH) `Node::getCommittedDisp()` does not exist → **adopted: use `getDisp()`**
(returns the committed buffer post-commit). (D2, med) committed disp valid at both read sites only because
`commitState()` precedes `theContactDomain->commit()` → **adopted: ordering locked + asserted.** (D3, HIGH) percentile
runaway clip collapses the deformed grid and silently drops pairs; `runawayGuardFired()` never surfaced → **adopted:
BLOCKER-CLIP.** (D4, HIGH) orientDir is stack-local, not persisted, and re-deriving from a penetrating current config
inverts the normal → **adopted: BLOCKER-ORIENTDIR (persist from first reference handle).** (D5) auto-kn genuinely
config-independent → verified, no change. (D6, med) band from deformed `medianSegDiag` drifts vs the next grid → **adopted:
band keyed off the reference median (D6).**

**Lens E — scope & completeness.** (E1, major) "rebuild the bucket grid" is wrong for mortar/edge-edge (brute force,
not the grid) → **adopted: scope narrowed to NTS; mortar already finite-sliding-correct.** (E2, major) GA-1 is a
correctness guard on the added line → promoted to P0. (E3, BLOCKER) H6 "partial mitigation" is false (anchor=node-keyed,
FrictionState=segIndex-keyed) → **adopted: struck; replaced by BLOCKER-MEMBERSHIP refusal.** (E4a, BLOCKER) slide-off
the surface entirely → **adopted: BLOCKER-SLIDE-OFF.** (E4b, major) self-contact via shared nodes after folding →
**adopted: BLOCKER-SELF-CONTACT filter.** (E4c, major) implicit Newton convergence at a re-emit step → Q-IMPLICIT-NEWTON
gate. (E4e, major) restart/sendSelf of the anchor store → **adopted: reconstructed at next handle (Q-RESTART).** (E4f,
major) MP re-emit → **adopted: D7 serial-only refusal.** (E5, major) the P1 continuity gate over-promised vs
faceted-normal chatter → **adopted: split into no-pass-through (tight) + resultant (loose, cite Q-NORMAL).**

## Post-code review + remediation backlog (R1–R8, 2026-06-30)

After P0+P1+P2 shipped (#443/#444/#446/#447) a post-code 5-lens review (recorded in the session memory)
verified the **core is sound vs source** — OFF byte-identical; D4 fresh-slot re-engagement traction-continuous
(`gT0` self-capture ⇒ `gTeff−gpT=0`); `segIndex` a stable global `mTags` ordinal (no aliasing under pure re-emit);
GA-1 holds; `getDisp()` ordering correct; rejected steps don't arm — but found several ADR-promised P0/P1
dispositions that **did not actually land**. Remediation backlog, severity-ranked (recommended order R2→R1+R6→R3→R5→R4/R7/R8):

| ID | Sev | Gap | Status |
|---|---|---|---|
| **R2** | HIGH | Search band was fed the **deformed** `segCoords` (handler), not reference ⇒ not deformation-invariant (contradicts D6 + the `referenceBand()` contract). | **FIXED (this PR):** fill `segCoords` reference-first → band from reference → shift to deformed in place for the grid. |
| R1 | HIGH | BLOCKER-MEMBERSHIP absent — no `mTags` fingerprint, so a friction key aliases if `mTags` reorders via element removal / re-mesh. | open |
| R3 | HIGH (curved) | `orientDir` not persisted; re-derived from reference each handle. Planar-safe (+`-outward` safe), but a sharp convex ridge can flip the normal → silent pass-through. The convex-ridge gate test is absent. **Full fix needs consistent-winding normals or nodal smoothing (ADR-47), not just persistence.** | open — **use `-outward` on curved/non-planar masters** (documented limitation) |
| R5 | MED | BLOCKER-SLIDE-OFF absent (relies on GC + projection; a stale adapter can persist ≤ floor; an edge clamp can hold spurious traction). | open |
| R6 | MED | D7 serial-only refusal absent in `needsResort` (the trigger runs on every rank under MP). | open |
| R4 | MED | `Trigger` reborn every handle (`clearReemit` drops it) ⇒ hysteresis/floor vestigial; the anchor-rebuild is the de-facto rate-limiter (so practical impact is low). | open |
| R7 | MED | `-resortEvery` shipped as a min-floor, not the D1 LS-DYNA-BSORT forced cycle-cadence. | open |
| R8 | LOW | `runawayGuardFired()` never surfaced (dead warning). | open |

**Exposed combo today:** `-reemit` + `-mu` on a curved / re-meshable master (no guard yet gates friction-under-re-emit
for those). Flat / `-outward` masters — the shipped + tested path — are unaffected.
