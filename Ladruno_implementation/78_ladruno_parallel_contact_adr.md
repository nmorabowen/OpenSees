---
title: Parallel (MPI) contact — owner-rank interactions on replicated-boundary-node decks
project: Ladruno
status: proposed — design of record; 2026-08-11 adversarial gate folded; **P0 RUN AND PASSED same day (implicit + explicit, 3-mutation battery)**, surfacing one new blocker (P0.a — the Tcl MPI interpreter has no contact verbs) and one new invariant (INV-6, massless ghosts)
priority: high
owner: nmora
amends: 39_ladruno_contact_domain_adr, 41_ladruno_mortar_alm_contact_adr
tags:
  - adr
  - contact
  - parallel
  - mpi
  - explicit
  - implicit
---

# ADR-78 — Parallel (MPI) contact

> ADR-78. Closes the **"Parallel (MPI): OUT of v1 — refuse partition-crossing pairs with a
> named error"** line in [[39_ladruno_contact_domain_adr]] §Risks, and the matching
> `[SERIAL-ONLY]` entries in the [[48_ladruno_contact_capstone_adr]] status table. Detailed
> design of record for the **parallel lane**; the capstone stays authoritative for global
> contact status and sequencing.
>
> Family: ADR-39 (NTS penalty) · ADR-41 (mortar/ALM) · ADR-47 (deferral ledger) ·
> ADR-48 (capstone) · ADR-57 (edge-edge) · ADR-60 (finite-sliding re-emit) ·
> ADR-74 (ParallelNumberer — the MP setup lane this rides on).
>
> **Cross-library.** The engine half lives here; the deck-authoring half is
> **apeGmsh ADR 0092** (`src/apeGmsh/opensees/architecture/decisions/0092-partitioned-contact-emit.md`),
> which owns owner-rank selection, ghost declaration, and the removal of today's blanket
> `BridgeError`. Neither half ships alone.

## What

Make `contact` / `contactPlane` work under **OpenSeesMP partitioned decks**, where today
every contact verb is refused at emit time and would silently produce a *partial* interface
if it were not.

The scope is deliberately narrow: **one owner rank per contact interaction**. The rank that
owns an interaction runs the entire broad phase, narrow phase, and path state for it, and
reaches the other ranks' nodes through **ghost node declarations** — the same replication
mechanism apeGmsh already uses for cross-partition MP constraints. Contact contributions
land in the distributed SOE and are summed with the native rank's internal forces by the
solver.

**Not in scope** (fenced, see §Non-goals): OpenSeesSP, the `Subdomain`/actor architecture, a
distributed broad phase, dynamic contact repartitioning, cross-rank self-contact, and
contact-aware load balancing (P5 sketch only).

## Why this is far cheaper than it looks

The 2026-08-11 scoping found that the premise behind ADR-39's deferral — that parallel
contact needs OpenSees' `PartitionedDomain` / `Subdomain` actor machinery, dynamic ghost
nodes, and an in-Newton allreduce — **does not apply to how the fork is actually driven in
parallel.**

apeGmsh partitioned decks are **manual domain decomposition**: one plain `Domain` per rank
inside `if {[getPID] == K}` guards, with boundary nodes **replicated by tag** across every
rank that needs them, and global assembly performed by the distributed SOE under
`numberer ParallelPlain` / `ParallelRCM` (apeGmsh ADR 0027 / ADR 0061; the auto-emitted
`system Mumps` chain at `apesees.py:5595`, and — this is the one to confirm in P0 — the
`DistributedDiagonalSOE` / `MPIDiagonalSOE` path for the explicit lane).

Three consequences, and they are the whole design:

1. **Ghost nodes already exist as a shipped, exercised mechanism.** ADR 0027 declares foreign
   nodes on the non-owning rank via `node(tag, *xyz, ndf=…)` and replays the owner's ordered
   `fix` / `remove sp` stream (INV-2, amended twice for exactly the singular-matrix failure
   naive replication causes). Contact needs the same thing, applied to the interface node set.
2. **Contributions to a shared DOF are summed by the SOE.** An adapter assembling into a
   ghost DOF on rank K reaches the physical DOF on its native rank through the solver, not
   through new MPI code. The contact engine stays MPI-free.
3. **Keeping an interaction whole on one rank makes mortar's global reduction rank-local by
   construction.** The mortar pressure `p_I = min(0, λ_I + εN·ḡ_I)` needs the *global*
   weighted gap over every facet touching slave node `I` (`accumulateMortarGap`,
   `LadrunoContactDomain.h:325`, order-independence contract at `:295`). Split the interface
   across ranks and that becomes an allreduce **inside every residual evaluation**. Don't
   split it, and there is nothing to reduce. This single fact is why the MVP is
   owner-rank-per-interaction and not slave-node-partitioned contact.

**Why now.** The MP setup lane is healthy (ADR-74 root-caused and fixed the `ParallelNumberer`
O(V²) stall; production explicit runs at 18.56 M hex / np240 are on record). Contact is the
last flagship capability that forces a user back to single-process — which caps contact models
at whatever one node's RAM holds, exactly the wrong constraint for pounding and
progressive-collapse work.

## Where

**Fork — modify:**
- `SRC/analysis/handler/LadrunoContactHandler.cpp` — the missing-node paths (`:625–648`,
  `:732`, `:972`, `:1157`) must become a **hard, named error** instead of the current
  "malformed model" coordinate backfill; the `-kn auto` and SOFT `m_eff` sizing paths need a
  deck-supplied route (§How D4).
- `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — `sendSelf`/`recvSelf` (absent today);
  a `partitionOwner` field per interaction if P2 needs a run-time assertion.
- `SRC/analysis/handler/LadrunoContactHandler.cpp:1306` — the `sendSelf`/`recvSelf` stubs whose
  own comment defers to `LadrunoContactDomain::sendSelf`.

**Fork — reference (patterns to copy):**
- `SRC/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSOE.{h,cpp}` — the explicit-lane
  shared-DOF summation contract P0 must verify.
- `SRC/domain/constraints/` MP replication behaviour under `ParallelPlain`.

**apeGmsh — modify** (owned by its ADR 0092): `opensees/apesees.py:2309` (the blanket
`BridgeError`), the partition-graph weighting in `mesh/partitioning`, and
`opensees/_internal/build.py:emit_contacts` / `emit_contact_planes` (rank scoping).

**Build:** no new target, no new dependency. MPI is already linked in the MP build.

## How

### D1 — one owner rank per contact interaction

Every `contact` / `contactPlane` verb is emitted inside **exactly one** rank's `if {[getPID]
== K}` block. That rank builds all adapters for the interaction, runs its broad phase, and
holds all path state (`FrictionState`, `MortarNormalState`, `EdgeEdgeState`). Emitting the
same verb on two ranks would double-count the contact force — forbidden by construction
because the emitter emits once.

**Owner = the rank owning the MASTER surface's backing solid elements** `[REVIEW: corrected]`,
ties broken by lowest rank index for determinism. The first draft said "majority of slave
nodes"; the review found that `-kn auto` resolves the owning solid of the **master segment**
(D4), so a master-side owner keeps auto-sizing native and reduces the ghost set to
geometry-only slave nodes. Nothing else in the narrow phase reads elements: mortar integrates
the slave facet from coordinates, and all path state is engine-side.

### D2 — ghost the whole interface on the owner rank

Every node of both surfaces the owner does not natively own is declared in the owner's block
with `node(tag, *xyz, ndf=…)`, plus the owner's ordered SP stream (ADR 0027 INV-2, in both
directions under staging).

**Whole surface, not a geometric halo.** A halo sized from the initial configuration is unsafe
under ADR-60 finite-sliding re-emit and under any large-sliding run: the candidate set churns
every step and the engine has no way to acquire a node mid-analysis. The whole-interface set
is bounded by *interface* size, not model size, and is immune to active-set churn.

### D3 — no MPI in the contact engine

The adapter assembles `F_c` / `K_c` into ghost DOFs; the distributed SOE sums them into the
native rank's equation. `LadrunoContactDomain` gains serialization (so the definitions survive
whatever the MP path does with analysis objects) but **no communication**. This is the single
most important scope boundary in this ADR: the moment contact needs its own MPI calls, the
design has failed and should be re-gated.

### D4 — sizing that reads *assembled* state must move into the deck `[GATE]`

Two shipped features derive their value at `handle()` time from state a ghost node does not
have on the owner rank:

- **`-kn auto`** sizes the penalty from the owning solid element's stiffness. `[REVIEW]`
  `ladrunoResolveAutoKn` (`LadrunoContactHandler.cpp:123–172`) scans the Domain for *the first
  element containing **all** of the **master segment's** nodes* and reduces
  `kn = f_si·mean(nᵀK_block n)`. So it needs the **master** side's backing solid on the rank —
  which the D-REVIEW-1 owner rule now guarantees natively. **And its failure is silent:** no
  owning element ⇒ `return -1.0` ⇒ *the caller skips the pair with a warning*. Under
  partitioning that is a **silently smaller interface** — the same failure class as D5, on a
  different path.
- **`-soft <SOFSCL>`** replaces `kn` with `SOFSCL·4·m_eff/dt²`, where `m_eff` comes from the
  handler's **assembled** nodal-mass cache — rebuilt in *one pass over the rank's elements*
  (`LadrunoContactHandler.cpp:174–181`), because `Node::getMass()` holds only the nodal piece.
  `[REVIEW]` `m_eff = 1/(invM_slave + Σ Nᵢ²·invM_segᵢ)` needs **both** sides' assembled mass,
  so no owner rule can make SOFT native: whichever surface is ghosted contributes zero mass,
  `m_eff` comes out too small, `k_soft` too soft — **silently**. SOFT is therefore refused
  under partitioning until the fork accepts a deck-supplied effective mass.

**The two split apart under the corrected owner rule** `[REVIEW]`:

- **`-kn auto` survives.** Master-side ownership (D1) puts the segment's backing solids on the
  owner rank natively, so `ladrunoResolveAutoKn` resolves exactly as it does in serial —
  *provided the master surface is uncut*, which apeGmsh ADR 0092 INV-4 now enforces by
  weighting the master side as uncuttable. No deck-supplied `kn` needed. The first draft's
  "compute it in apeGmsh and emit a literal" is **withdrawn**; a literal would have silently
  diverged from the serial value.
- **`-soft` does not.** It needs both sides' assembled mass, so it is refused under
  partitioning with a named error until the fork accepts a deck-supplied `m_eff`.

What *must* still change is the **failure behaviour**: auto-kn resolving to `-1.0` must abort
instead of skipping the pair (D5.2), because under partitioning that is the difference between
"your mesh is malformed" and "rank 7 quietly dropped a third of your interface".

### D5 — every silent degradation becomes loud `[REVIEW: widened]`

Three paths currently degrade an interface without failing. All three become hard errors
naming the contact tag, the node/segment, and the rank:

1. **Missing node** — `LadrunoContactHandler.cpp:625–648` marks it `HUGE_VAL` and backfills
   with the segment's first valid coordinate ("malformed model"); the narrow loop skips it.
2. **Auto-kn with no owning solid** — `ladrunoResolveAutoKn` returns `-1.0` and the caller
   *skips the pair with a warning* (`:123–172`). Found by the adversarial review; it was not
   in the first draft of this ADR.
3. **SOFT with an incomplete mass cache** — silently under-stiff (D4).

Under partitioning each of these is *a silently partial interface*, the worst failure mode
available to this design, and a warning on rank 7 of 240 is not a defence.

### D6 — MP-constraint composition is a *conditional* dependency `[REVIEW: corrected]`

The first draft claimed partitioned decks carry cross-rank MP constraints "by construction".
**That is wrong.** apeGmsh replicates a partition-boundary node as a plain `node` line plus a
replayed `fix` / `remove sp` stream (`emit_ghost_sp_ops`, `ghost_tags_by_rank`,
`_bucket_fix_targets_by_rank` — `apesees.py:2501–2870`); shared DOFs are reconciled by the
numberer and the SOE, **not** by MP constraints. A partitioned deck carries MP constraints
only when the *model* declares them.

`LadrunoContactHandler` still replicates `PlainHandler` and does not enforce them
(`LadrunoContactHandler.h:26–34`), so contact + MP stays refused (`apesees.py:5257`). But
ADR-39's **Q-CONSTR is not a blocker for this lane** — it bounds which models the first ship
covers. Two-body pounding, foundation rocking against a rigid plane, and mortar mesh-ties
across ranks all declare no MP constraints.

### Invariants

- **INV-1 (single assembly).** Each contact interaction is assembled by exactly one rank.
- **INV-2 (complete interface).** The owner rank can resolve **every** node of both surfaces;
  a missing node aborts with a named error, never a backfill.
- **INV-3 (rank-local state).** All path state for an interaction lives on one engine;
  `commit()` / `revertToLastCommit()` need no communication.
- **INV-4 (no engine MPI).** `SRC/domain/contact/` contains no MPI call.
- **INV-5 (serial byte-identity).** A non-partitioned deck is byte-identical to pre-ADR-78.
- **INV-6 (massless ghosts).** A ghost node declares geometry + SP state only — never mass,
  never elements, never loads. The distributed solver sums `A` as well as `B`, so a ghost that
  carries mass double-counts it on every shared DOF and both ranks converge on the *same wrong
  answer*. Added after P0 measured it (see §Implementation log, P0.e).

## Phases

| P | Scope | Gate |
|---|---|---|
| **P0** ✅ | **No SRC change.** 2-rank deck: contact wholly on the master-owning rank, slave interface ghosted, under `system Mumps` (implicit) and `MPIDiagonal` (explicit), vs the identical serial model. | **PASSED 2026-08-11** — implicit within 1.6e−14, explicit bit-identical, 3 mutations each fail as designed. `-kn auto` resolves natively on the master-owning rank, confirming D1. See §Implementation log. |
| **P0.5** ✅ | Register the contact family in `SRC/tcl/commands.cpp` (+ declarations in `commands.h`) so the classic engine — and therefore `OpenSeesMP` — can see it. Registration only; no parser duplicated. | **DONE 2026-08-11.** A 2-rank contact deck **in Tcl** matches its serial twin (1.4e−14) and both match the Python-lane reference to every digit. See §Implementation log P0.5. |
| **P1** | D5 fail-loud (missing node → named error) + a `contact`-verb-on-wrong-rank abort. | Serial battery byte-identical; a deliberately incomplete interface aborts instead of running. |
| **P2** | D4 as corrected: `-soft` refused under partitioning with a named error (`-kn auto` needs no change). `LadrunoContactDomain::sendSelf`/`recvSelf`. | The refusal fires; a partitioned `-kn auto` run resolves the **same** penalty as its serial twin. |
| **P3** | apeGmsh ADR 0092: partition-graph weighting, owner selection, ghost + SP replay, `BridgeError` relaxed to a locality check. | A `g.constraints.contact` model emits a runnable 2-rank deck; partition-straddling cases still refuse with the *specific* reason. |
| **P4** | Validation battery: 2-body pounding on 2 ranks (explicit, NTS + `-soft` if P2 shipped it), mortar mesh-tie across ranks (implicit + ALM), 4-rank strong-scaling sanity. | Each matches its serial reference; contact energy balance closes. |
| **P5** | *Deferred sketch only.* Slave-node-partitioned contact + facet ghosting for load balance — requires the in-Newton allreduce for `ḡ_I` that D1 exists to avoid. | — |

## Risks / open questions — **all five settled 2026-08-11, see §Adversarial review log**

> [!success] Q-EXPLICIT — **RESOLVED, mechanism confirmed in source**
> `DistributedDiagonalSolver::solve()` allocates `dataShared[2*numShared]` — *"2 times for A &
> B"* (`DistributedDiagonalSOE.cpp:261`) — packs both the diagonal **mass** `A[loc]` and the
> **residual** `B[loc]` for every shared DOF, gathers to P0, **sums** (`*vectShared +=
> otherShared`, `DistributedDiagonalSolver.cpp:131`), broadcasts back, and writes the summed
> values into each sharing rank's local `A` and `B`. A ghost DOF on the owner rank therefore
> contributes `A=0, B=F_contact`; the native rank contributes `A=m, B=F_int`; after the
> exchange both hold `A=m, B=F_int+F_contact` and `X=B/A` is correct on both. Shared-DOF
> detection is bilateral through the same P0 gather (`:146`, `:213`, `:231`). **The "no engine
> MPI" invariant holds for the explicit lane.** P0 now *confirms* a known mechanism rather
> than discovering one.

> [!question] Q-P0GATHER — **new, raised by the review**
> The exchange above is a **gather-to-rank-0 every solve**, O(numShared) through one process.
> Ghosting a whole interface adds that interface's DOFs to the shared set on the owner rank.
> Negligible for one interface; worth measuring on the many-interaction case at P4. This
> replaces Q-GHOSTCOST as the real cost of ghosting.

> [!success] Q-GHOSTCOST — **RESOLVED, not a constraint**
> Deck cost is ~50 B per ghost `node` line. A 100³-hex model has a ~10⁴-node interface ⇒
> ~0.5 MB against ADR 0061's measured ~250 B/hex ⇒ ~250 MB deck: **~0.2%.** No ghost budget,
> no refusal threshold, no halo fallback. ADR 0092's sign-off question 2 is withdrawn.

> [!success] Q-RECORDER — **RESOLVED, already shipped**
> apeGmsh `f5358318` (2026-08-10) settled exactly this contract: across partitions a read is
> loud only when **every** rank is missing the group — *"a rank owning none of the recorded
> elements legitimately writes none, so `_per_partition` drops it from the stitch and
> re-raises only when nothing survived."* Contact-on-one-rank is the same shape. No new work;
> ADR 0092 INV-6 downgrades to a regression test.

> [!success] Q-HANDLER — **RESOLVED, not a blocker** (see D6)
> The "partitioned decks carry MP constraints by construction" premise was false. Q-CONSTR
> bounds the first ship's model coverage; it does not gate the lane.

> [!question] Q-BALANCE — **refined, still open**
> The review sharpened the failure mode: many interactions **balance naturally** (each gets
> its own owner rank), so debris-field / AEM contact (ADR-51/54) is *fine*. The pathological
> case is a **single dominant interface** — one soil-structure interface on a np=240 run puts
> all of its contact work on one rank. Measured at P4; the fix is ADR-78 P5 and costs the
> in-Newton allreduce D1 exists to avoid.

## Non-goals

- **OpenSeesSP.** Explicitly out (decided 2026-08-11). SP replicates the whole domain per
  process, so contact is already correct there and buys nothing this ADR is about.
- **`Subdomain` / actor-based partitioning.** The fork's parallel decks do not use it; contact
  will not be the reason to start.
- **Distributed broad phase**, dynamic contact repartitioning, cross-rank self-contact.
- **Contact-aware load balancing** (P5 sketch only).

## Adversarial review log (2026-08-11, source-grounded, pre-P0)

**Verdict: SALVAGEABLE-WITH-CHANGES → now SOUND.** The architecture survives; **two of its own
claims were falsified** and one risk was upgraded from "accepted" to "must be caught". All
corrections are folded in above and flagged `[REVIEW]`.

**Confirmed sound:**
1. **The enabling premise holds.** Replicated boundary nodes + a distributed SOE that sums
   shared-DOF contributions is real, for the explicit lane as well as the implicit one —
   `DistributedDiagonalSolver.cpp:100–151`. The engine needs no MPI.
2. **Owner-rank-per-interaction keeps mortar's `ḡ_I` reduction rank-local.** The alternative
   costs an allreduce per residual evaluation. The trade is correctly identified.
3. **Whole-interface ghosting over a halo.** Justified by ADR-60 re-emit *and* cheap
   (Q-GHOSTCOST: ~0.2% of deck).

**Falsified (BLOCKER / MAJOR):**
- **B1 [D1, owner rule] — "majority of slave nodes" is the wrong side.** `ladrunoResolveAutoKn`
  (`LadrunoContactHandler.cpp:123–172`) resolves the owning solid of the **master segment**.
  Master-side ownership keeps `-kn auto` native and shrinks the ghost set to geometry-only
  slave nodes. → D1 corrected; apeGmsh ADR 0092 INV-1 corrected; INV-4 sharpened to *cut on
  the slave side, never the master*.
- **B2 [D4/D5] — auto-kn failure is silent, and was missed.** No owning solid ⇒ `return -1.0`
  ⇒ *the caller skips the pair with a warning*. A partitioned run would lose contact pairs
  without failing. → D5 widened from one path to three.
- **B3 [D6] — "partitioned decks carry MP constraints by construction" is false.** apeGmsh
  ghosts are `node` + replayed `fix` (`emit_ghost_sp_ops`), not MP constraints. → Q-CONSTR
  demoted from hard dependency to a coverage bound; the first ship is genuinely useful
  without it.

**Reframed:**
- **Q-BALANCE** — many interactions balance naturally (one owner each); the pathological case
  is a *single dominant* interface. The original framing had this backwards.
- **Q-GHOSTCOST → Q-P0GATHER** — the cost of ghosting is not deck bytes, it is the
  gather-to-rank-0 shared exchange volume in `DistributedDiagonalSolver::solve()`.

**Unchanged by review:** INV-1..INV-5, the phase ladder, and the non-goals.

## Implementation log

### P0 — RUN 2026-08-11. **PASS on both lanes, with one new blocker.**

Model: two stacked unit `stdBrick`s (E=2e7, ν=0), NTS penalty contact, `-kn auto`,
`-outward 0 0 1`, lateral DOFs fixed ⇒ a 1-D column. Bottom block = master (rank 0),
top block = slave (rank 1), slave interface nodes 11–14 ghosted onto rank 0 with their
`fix` replayed. Serial reference = the same model in one domain. Build pinned to the
installed Ladruno tree (Aug-10, `4a6aeec5`) via `LADRUNO_OPENSEES_BIN` /
`LADRUNO_OPENSEESMP_BIN` — see the P0.b finding below.

**Implicit lane** (`numberer ParallelPlain` + `system Mumps` vs serial `RCM`/`UmfPack`):

| quantity | serial | 2-rank | rel. Δ |
|---|---|---|---|
| `w15` (top, rank 1) | −5.624999999999997e−3 | −5.624999999999919e−3 | 1.4e−14 |
| `w11` (interface) | −5.124999999999998e−3 | −5.124999999999914e−3 | 1.6e−14 |
| `w5` (master face, rank 0) | −5.000000000000000e−4 | −5.000000000000034e−4 | 6.8e−15 |
| `ΣR_base` | 10000.0 | 10000.0 | 0 |

Serial is analytically exact (each block compresses `PL/EA` = 5e−4; `ΣR` = the 10 kN
applied). `w11` on rank 1 and the **ghost** `w11` on rank 0 agree bit-for-bit.

**Explicit lane** (`CentralDifferenceLadruno`, 200 × 5e−5 s, `MPIDiagonal` vs serial
`Diagonal`) — **bit-identical to every printed digit**:
`w15` 0.0037458245314017234 · `w11` 0.003692142684617696 · `w5` −4.5278662184168605e−05,
serial and 2-rank alike, ghost included. **Q-EXPLICIT is now measured, not just read.**

**Mutation battery** (a green smoke proves nothing until a broken one fails):

| mutation | result |
|---|---|
| ghosts not declared on rank 0 | `WARNING … contact 1: node 11 of contactSurface 2 does not exist; contact SKIPPED`, then `analyze -3` |
| contact verb emitted twice | converges to a **plausible wrong answer**: `w15` = −2.8125e−3, exactly half the penetration (double stiffness). `ΣR` still 10000 — equilibrium holds, compliance is wrong |
| ghost nodes given mass | converges, both ranks **agree**, answer wrong by ~20% (`w11` 3.69e−3 → 2.96e−3) |

### P0 findings that change the ADR

- **P0.a — NEW BLOCKER: the Tcl MPI interpreter cannot see the contact verbs.**
  `OpenSeesMP` links `${OPS_SRC_DIR}/tcl/commands.cpp` (CMakeLists §"OpenSeesMP Tcl MPI
  Interpreter"), but `contactSurface`/`contact`/`contactPlane` are registered **only** in
  `TclWrapper.cpp:1955` and `PythonWrapper.cpp:3334`. Measured: `OpenSeesMP.exe` (and
  `OpenSees.exe`, and the installed Aug-10 build) answer `invalid command name
  "contactSurface"`. **apeGmsh's partitioned emit target is Tcl** (ADR 0061 per-rank
  `.tcl`), so a correctly partitioned contact deck could not run today for a reason that
  has nothing to do with partitioning. P0 ran on `openseesmp` (the Python MPI module,
  which links `PythonWrapper.cpp`). ⇒ **New phase P0.5, ahead of everything else:** either
  register the contact verbs in `commands.cpp`, or make partitioned contact emit a
  per-rank **Python** deck. This is a decision for apeGmsh ADR 0092 to fence.
- **P0.b — build resolution is a live hazard.** A bare `import opensees` / `import
  openseesmp` in the dev venv resolved to a copy left in *another session's scratchpad*
  (`ladruno_opensees.pth` → `_ladruno_opensees_boot`). Serial and MP must be pinned to the
  same build or the comparison is meaningless. The `LADRUNO_OPENSEES_BIN` /
  `LADRUNO_OPENSEESMP_BIN` escape hatch exists for exactly this; P0 uses it.
- **P0.c — D5 is partly overstated, and the real gap is narrower.** The missing-**slave**-node
  path is *not* silent: the handler names the contact, the surface, and the node, and skips
  the contact. It is still a **warning, not an abort**, so a model whose bodies are otherwise
  restrained converges with the contact quietly absent — the substance of D5 stands, but its
  first bullet should say "warns and skips the whole contact", not "backfills a coordinate".
  (The coordinate backfill at `:625–648` is the **master**-segment path and was not exercised.)
- **P0.d — INV-1 cannot be enforced by the engine.** The duplicate-contact mutation *was*
  caught by a handler warning — but that detector is **rank-local** (each rank's handler sees
  only its own contacts), so a contact verb emitted on two ranks would produce the same
  plausible-wrong answer with **no** warning at all. Single assembly must be guaranteed at
  emit time. This upgrades apeGmsh ADR 0092 INV-1 from "structural convenience" to "the only
  available guard".
- **P0.e — NEW INV-6: ghosts must carry no mass.** Measured above. Because the distributed
  solver sums `A` as well as `B`, a ghost that declares mass double-counts it on every shared
  DOF — and both ranks then agree on the *same wrong answer*. apeGmsh's ADR 0027 ghosts are
  `node` + replayed `fix` only, so this is already true; it now has to stay true deliberately.

**Not covered by P0:** mortar/ALM across ranks, friction, `-soft` (refused by design),
`contactPlane`, >2 ranks, and a cut master surface. Those are P4.

### P0.5 — DONE 2026-08-11. Contact family registered in the classic Tcl engine.

`SRC/tcl/commands.cpp` (+147) and `SRC/tcl/commands.h` (+40), both purely additive:
ten `extern` declarations, ten two-line bridges
(`OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, &theDomain)` → the shared
`OPS_Ladruno*()`), ten `Tcl_CreateCommand` registrations, and ten forward declarations in
the header. **No parser is duplicated** — the `OPS_*` entry points live in
`OpenSeesOutputCommands.cpp`, which is in the `OPS_INTERPRETER` object library folded into
`OpenSeesLIB`, so they were already *linked* into every target and only the registration was
missing. Mirrors the ADR-44 modal-family bridges verbatim, including *why* the header
declarations are needed: the registrations sit ~5600 lines before the definitions.

The whole family is registered, not just the three definition verbs — the ALM workflow
(`ladrunoBeginAugment` / `ladrunoMortarPenetration` / `ladrunoEndAugment`) is how a mortar
contact is driven, so a deck that could DEFINE a contact but not augment it would be a trap.
Verified that the query commands survive the bridge: `OPS_SetIntOutput` /
`OPS_SetDoubleOutput` in `elementAPI_TCL.cpp` append to the interpreter result.

**All ten pairings and all ten runtime paths were checked, not just the two the P0
deck exercises.** The name→bridge→`OPS_` chain was composed for each command and
diffed against the `TclWrapper.cpp` reference mapping — a swapped pair (e.g.
`ladrunoMortarPenetration` wired to `OPS_LadrunoEdgePenetration`) compiles cleanly
and is silently wrong, so it cannot be left to the compiler. All ten match. The
eight commands the P0 deck never calls were then smoke-run against a domain with a
live contact: `ladrunoContactInfo` → `1 0 0 0`, the three penetration queries →
`0.00000000000000000000` (the `%35.20f` of `OPS_SetDoubleOutput`), the augment pair
→ empty, `ladrunoContactForce` → its own arity error, `contactPlane` → its usage
path. None answered `invalid command name`.

**Measured on the freshly built engine (`Ladruno_scripts\build.bat OpenSees OpenSeesMP`):**

| | serial Tcl | 2-rank Tcl | Python reference (installed build) |
|---|---|---|---|
| `w15` | −5.624999999999997e−3 | −5.624999999999919e−3 | −5.624999999999997e−3 / −…919e−3 |
| `w11` | −5.124999999999998e−3 | −5.124999999999914e−3 | −5.124999999999998e−3 / −…914e−3 |
| `w5` | −5.000000000000000e−4 | −5.000000000000034e−4 | idem |
| `ΣR` | 10000.0 | 10000.0 | 10000.0 |

`contactSurface` now answers on `OpenSees.exe` **and** `OpenSeesMP.exe` (was
`invalid command name`). Two results worth separating:

1. **P0.5's own gate** — a partitioned contact deck in **Tcl**, the language apeGmsh's
   per-rank emit actually speaks (ADR 0061), runs and matches its serial twin to 1.4e−14,
   with the rank-0 **ghost** `w11` bit-identical to rank 1's native value.
2. **A free build-provenance check** — the serial Tcl numbers reproduce the Python-lane
   reference *to every digit*, and the 2-rank Tcl numbers reproduce the Python 2-rank
   numbers to every digit. Two interpreters, two builds, one answer. That retires the open
   worry about whether the rebuilt engine was configured equivalently (lp64 / sequential /
   static MKL, MUMPS from `mumps-install`).

### P1 — DONE 2026-08-11. Every silent degradation aborts, and it tears the job down.

**Fifteen paths converted** from warn-and-continue to abort: the three
`ladrunoSurfaceNodesOk` call sites (NTS / mortar / edge-edge), `-kn auto` failing to size a
penalty, SOFT with no assembled mass at an interface node (D4), and eleven configuration
skips (surface undefined, wrong surface kinds, unsupported `nodesPerSeg`, `kn <= 0`, slave
`ndf != 3`, and the `contactPlane` equivalents). `grep "; skipped"` in the handler now returns
**zero**. Two were worse than the rest: `contactPlane` with an undefined slave surface was a
bare `continue` with *no message at all*, and the rigid-plane `ndf < ndm` path dropped
individual slaves, leaving a plane that restrained only some of them.

**The abort had to leave the handler's translation unit.** `LadrunoContactHandler.cpp` is
compiled into the `OPS_Analysis` **OBJECT library** — built once, no parallel define, folded
into every target — so the first version of this fix, an `#ifdef _PARALLEL_*` block written
inside the handler, **compiled to nothing in every build, OpenSeesMP included**. It compiled
clean and left both happy paths bit-identical; only the mutation deck caught it, by hanging
exactly as it had before. Same trap ADR-02 calls "the core blocker" for `OpenSeesCommands.cpp`.

Fixed the way Patch 8 / Patch 9 fixed theirs: the one parallel-sensitive function lives in
`SRC/analysis/handler/LadrunoContactAbort.cpp`, listed in `OPS_MPI_PER_TARGET_SOURCES` and
compiled **per target** with that target's defines. Verified in the generated `build.ninja` —
`-D_PARALLEL_INTERPRETERS` present on `OpenSeesMP`'s copy, absent on the serial one. Isolating
one function rather than moving the 1400-line handler keeps the duplicated compilation
trivial; the CMake edit surface is identical either way.

**Why teardown and not just `return -1`:** measured twice. A failing rank returns from
`handle()` while its peers walk into the numberer's gather and block on a rank that never
arrives. The job hung until an external timeout, twice leaving orphaned `OpenSeesMP.exe`
processes alive for hours — which on the second occasion held the binary open and failed the
next link. On a 240-rank run that burns the allocation and buries the FATAL line in one log.

| | before P1 | after P1 |
|---|---|---|
| serial happy path | reference values | **unchanged** |
| 2-rank happy path | reference values | **unchanged** |
| `noghost` mutation | warned, limped to `analyze=-3`, **hung 7 min**, orphans left | FATAL + teardown, **terminated in 1.1 s**, zero orphans |

**Found in passing, not fixed:** `AutoConstraintHandler.cpp` is in the same object library and
is *not* listed per-target, so its `MPI_Allreduce` — the one computing the global penalty for
`constraints Auto` — is compiled out of every target. Under MPI each rank would size its
penalty from rank-local stiffness. Verified from its `DEFINES` line; pre-existing and unrelated
to contact, so it deserves its own change rather than riding along here.

### P1 amendments — adversarial review, same day. Two of P1's own claims were wrong.

**A1 — the SOFT guard rejected correct models, on a false premise.** It scanned *both* surfaces
and justified itself with "k_soft collapses to 0 and the contact silently does nothing". Both
halves are refuted by the adapter it claims to protect:

* `gapModeInvMass` (`LadrunoContactFE.cpp:448`) treats a DOF with `m ≤ 0` as **infinite** mass —
  *"the correct gap-mode treatment of a FIXED master node … a fixed/rigid master contributes 0,
  leaving m_eff = m_slave"*. A rigid indenter or drum meshed as fixed, element-free master nodes
  is a standard pattern, and P1 aborted on it.
* `softKn` (`LadrunoContactFE.cpp:516`) already warns once and **falls back to the configured
  kn**. The pre-P1 behaviour was loud-warn-plus-safe-fallback, never silence.

Narrowed to **slave-surface nodes only** (a massless slave really does defeat soft-sizing, and
is the ghosted-node case ADR 0092 INV-3 refuses), and the message corrected.

**A2 — the conversion was incomplete, and P1's own verification hid it.** P1 reported
`grep "; skipped"` returning zero. Literally true, substantively misleading: the pattern matched
what had been *changed* rather than what was *claimed*. Four paths survived —
`"; facet skipped"` (mortar slave ndf), `"; pair skipped"` (mortar auto-penalty, edge-edge
auto-penalty), and one with **no message at all** (`if (!ok) continue;` on a master segment node
with `ndf != 3`, which `ladrunoSurfaceNodesOk` cannot catch because it checks existence and 3-D
coordinates but not ndf). All four now abort. The mortar auto-penalty one mattered most: it is
the same resolver and same failure as the NTS `-kn auto` abort that D5.2 named fatal, and
apeGmsh ADR 0092 INV-1 leans on it — the node-majority owner proxy is only safe because a
mis-pick fails loudly. That backstop did not exist for mortar mesh-ties, which P4 plans to run
partitioned.

### P1 KNOWN LIMITATION — incompatible with runtime element removal

`handle()` re-runs on **every** `domainChanged()`, not once, and `LadrunoContactDomain` has **no
API to retire a contact or prune a surface**. `RemoveRecorder` (the Elwood / Talaat-Mosalam
collapse lane) removes elements *and dangling nodes* at runtime. So a collapse run with a
declared contact, where a removed node belonged to a contact surface, now **aborts mid-run**
where it previously limped on with a warning — and under MPI the teardown skips `H5Fclose`, so
earlier stages' MPCO output can be lost.

The handler cannot presently distinguish "your deck has a typo" from "the analysis legitimately
removed this node", and the abort is right for the first and wrong for the second. Closing this
needs surface-membership pruning on `removeNode`/`removeElement`, or an explicit
contact-removal command. **Until then, do not combine `recorder Collapse` with a declared
contact.** Recorded here rather than fixed because it is a design gap, not a defect in P1.

### FOLLOW-UP — `AutoConstraintHandler`'s allreduce had never been compiled

Open item 1 from the P1 handoff, and the third instance of the same trap in this ADR.

`AutoConstraintHandler.cpp` sizes the `constraints Auto` penalty from the order of magnitude of
the **mean element-diagonal stiffness**, `PVAL = 10^(round(log10(KAVG)) + oom)`. Upstream wrote
four `MPI_Allreduce` calls to make `KAVG` a *global* quantity. That file compiles into the
`OPS_Analysis` OBJECT library — built **once, with no parallel define**, and folded into every
target — so its `#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)` block
compiled to nothing in **every** binary, `OpenSeesMP` included. It has been dead since the
handler was contributed.

Under MPI each rank therefore averaged only its **own** elements. `PVAL` is the value
`getPenaltyValue()` falls back to for any constrained node with no locally attached elements —
which is precisely a partition-straddling interface node — so the ranks sharing a constraint
size it from disjoint stiffness samples, with nothing printed.

Fixed the way P1 fixed its `MPI_Abort`: the reduction moved to
`SRC/analysis/handler/LadrunoAutoPenaltyReduce.cpp`, listed in `OPS_MPI_PER_TARGET_SOURCES`
(the P1 list, renamed from `OPS_CONTACT_PER_TARGET_SOURCES` now that it holds two unrelated
files) and compiled per target with that target's defines. The behaviour is a verbatim move —
same four collectives in the same order, same `OPS_PARTITIONED` quick-return, same warning
strings; only the compilation unit changed, which was the entire bug.

`AutoConstraintHandler.cpp` is deliberately left with **no `#ifdef _PARALLEL_*` at all** and an
unconditional call. That is the part that matters going forward: the trap is not that someone
wrote the wrong guard, it is that a correct-looking guard in an OPS_Analysis TU is invisible.
Removing the guard entirely is what stops it reopening silently.

**Scope, stated honestly.** This restores the *global* statistic (`KMIN`/`KMAX`/`KAVG` and the
`m_global_penalty` derived from it). It does **not** make the *per-node* penalties global —
those accumulate from elements attached to the node on that rank, so a replicated interface node
still gets a partial sum on each rank. That is inherent to an element-loop accumulation and would
need element ghosting rather than a reduction; it is also far less severe, since both ranks'
contributions are summed by the distributed SOE into a single, merely larger, effective penalty.
The fallback path fixed here is the one that can be orders of magnitude wrong.

**Gate:** `Ladruno_files/testbed/auto_penalty_mpi/`. The element populations are chosen so the
three possible answers are four orders of magnitude apart — 100 trusses at `k=1e2` give
`PVAL=1e5`, one truss at `k=1e6` gives `1e9`, and together they give exactly `1e7`. A 2-rank run
must print `1e7` on **both** ranks; the pre-fix signature is `1e5` on one and `1e9` on the other.
The two serial single-partition runs are what make that a fix rather than a coincidence: they
establish that `1e5` and `1e9` are what the two partitions produce alone.
