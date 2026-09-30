---
title: "WP-157 — Mortar friction path state per (slave node, facet pair)"
project: Ladruno
type: ADR (amends the ADR-41 C3 mortar friction lane; status row in the ADR-48 capstone)
status: "PR open (not merged) — resolves LEDGER_quirks MAJOR-1 (the C3.1 gate, #377) for friction"
owner: nmora
related:
  - "[[41_ladruno_mortar_alm_contact_adr]] (C3.1–C3.3: the lane this amends)"
  - "[[48_ladruno_contact_capstone_adr]] (status-of-record row; contract #1: path state on the Domain)"
  - "[[155_pile_contact_r05]] (§5–§6: the pile cohesion failure and the recommendation for this slice)"
  - "[[LEDGER_quirks]] (Mortar friction committed slip is last-writer-wins …)"
  - "apeGmsh piles-validation ladder/R1_plumbing/VERDICT.md, ladder/R0_5_fork_slice/REVIEW_fable.md finding 2"
tags: [adr, contact, mortar, friction, shared-node, pile, wp-157]
updated: 2026-09-30
---

# WP-157 — Mortar friction path state per (slave node, facet pair)

> [!summary] The short version
> The mortar friction state (committed slip `gpT`, engagement origin `gT0`/`engaged`, tangential
> multiplier `λ_T`, and the committed double-buffer) now lives in a `MortarFrictionState` keyed
> **(contactTag, slave node, slave-facet ordinal, master-facet ordinal)**, not in the per-node
> `MortarNormalState`. Every facet pair runs its own return map on its own local slip in its own
> tangent plane, so it must own the state it reads back. `λ_N` stays per node. The kernel math is
> unchanged. With one facet pair per slave node the arithmetic is the same, so the result is
> bit-identical; see §4 for the battery and the byte-identity dump.

## 1. The defect

ADR-41 C3.1 (#377) folded the friction fields into `MortarNormalState`, keyed (contactTag, slave
node). A slave node is integrated by every (slave facet, master facet) adapter that touches it. Each
adapter computes its own local weighted slip

    gbarT_I^pair = P_n(pair) [ Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K ] / a_I^pair

runs `frictionReturnMap`, and wrote `gpTtrial`/`lambdaTtrial` into the node's slot; the first pair
to evaluate captured `gT0`. The last writer won. The C3.1 gate fenced this as MAJOR-1 ("matched
meshes only"); C3.3 and WP-155 re-confirmed it. It has two effects:

1. **Order dependence where the pairs disagree.** With a non-uniform slip field the pairs sharing a
   node have different local slips, so the committed slip is that of whichever FE was evaluated last.
   Oracle T2: a linear slip field on a flat non-matching interface, forward against reverse sweep,
   differs by 27 % of the cap at step 2.
2. **A slip or traction read in the wrong plane.** On a curved or creased interface the pairs
   sharing a node have different normals. `gpT` and `λ_T` are 3-vectors in the writer's tangent
   plane, and the reader never projects them. From step 2 on, a pair read a neighbour's `λ_T`/`gpT`
   with a component along its own normal. Oracle T3 (a 24-facet polygon, the slave rotated half a
   facet, as on the R1 pile): step 1 is exact; at step 2 the normal leak `|t·n|/c` is 1.0 and the
   lateral friction force is 57 % low.

A uniform field (the matched battery, the flat patch under uniform shear) makes all pairs agree,
so the defect never showed there.

## 2. Decisions

### D1 — key the friction state per (slave node, facet pair)

`LadrunoContactDomain::MortarFrictionState` holds `gpT, gpTtrial, gT0, engaged, lambdaT,
lambdaTtrial, gT0committed, engagedCommitted` (the fields removed from `MortarNormalState`). The key
is `MortarPairKey {c, n, sf, mf}`: contact tag, slave node tag, and the GLOBAL slave- and
master-facet ordinals the handler loops over. Both ordinals are rebuild-stable, the same argument
as the NTS `FrictionState` key (contactTag, slaveTag, segIndex). The FE learns `mf` through
`LadrunoContactFE::setMortarMasterFacet`, called by the handler right after construction (3D and
2D loops). A 4-int composite key, never a hash (the EdgeKey rule).

This is the design the kernel already implies. The force of a pair is `D^pair t^pair`, computed from
the pair's local slip and local pressure, the C2.2 rule that keeps `R(u)` deterministic per facet.
The state is an output of that same local return map. Giving the pair its own slot makes it a
self-consistent integration cell: it re-reads exactly what it wrote. That gives commit invariance:
in full slip, re-evaluating at the converged `u` after `commit()` returns the same traction.

**Alternatives rejected.**

- *Key by `feTag`.* FE tags are reassigned at every `handle()` (the reason `theMortarFacetContribs`
  is cleared there).
- *One nodal state, blended by area over the pairs.* The slip is a return-map output in each pair's
  own plane. A blend needs a projection rule. It also hands each pair back a state it did not
  produce, so commit invariance fails: the blend moves at every commit. The ledger listed it as an
  option for the flat case only.
- *A nodal return map on the global weighted slip* (Σ over pairs, like `λ_N`'s accumulator). This is
  the consistent nodal-multiplier mortar. The residual sweep would need the global slip before
  every facet is evaluated, which the C2.2 rule forbids; it is the dual/nodal-LM redesign that
  ADR-47 defers.

### D2 — `λ_N` stays per node

The normal multiplier is Uzawa'd from the LINEAR global gap accumulator (`gtGlobal/aGlobal`), which
is order independent. Only the friction fields move. The tangents read `λ_N` from the node slot
(`nst`) and the friction state from the pair slot.

### D3 — lifecycle (capstone contract #1)

`commit()` promotes `gpT`, `λ_T` (gated by the WP-155 `-augment` mode, per contact) and the
engagement double-buffer. `revertToLastCommit()` drops the trials and restores `gT0`/`engaged`.
`revertToStart()` clears the store. GC: the handler marks every live frictional (node, pair) between
`mortarNormalGCBegin()` and `mortarNormalGCEnd()`, and End erases the unmarked slots. Only
frictional contacts are marked, so frictionless mortar and ties allocate nothing.

### D4 — a new pair engages fresh

A pair that first overlaps later (finite sliding onto a new master facet, or a re-pair) captures
`gT0` at its first in-contact evaluation, so it starts with zero stick traction. Its weight is its
overlap area, which grows from zero, and its traction reaches the cap after a slip of `cap/ε_T`.
The transient is of that order. The NTS per-segment key behaves the same way (ADR-60 D4). The
per-node layout inherited the old pair's state instead, which is right only on a flat interface.

### D5 — byte identity

When each slave node sees exactly one facet pair, the pair slot receives exactly the reads and
writes the node slot did, so the result is bit-identical. At nodes with several pairs the result
changes by design.

## 3. Oracle

`contact_prototypes/proto_adr157_mortar_pair_friction.py` (numpy only) mirrors
`addMortarFriction` line by line on a 1-D mortar interface embedded in 3-D and runs both layouts.
The kernel math does not change, so this oracle pins the storage invariants.

| Gate | Result |
|---|---|
| T1 flat non-matching, uniform slip, 3 alignments × 2 sweep orders × 3 cones: stick `−k_t s L`, slip `−cap L` | both layouts exact (≤ 1.6e-16) |
| T2 flat, linear slip, step 2: order independence and `f_K = −cap a_K` | pair: 3.6e-17 / 3.6e-17; node: order diff 0.27 cap, error 0.27 cap |
| T3 24-gon, lateral slip, cohesion: `|t·n| = 0` and `F = Σ_f c L_f dir_f` | pair: exact at both steps; node: step 2 leak 1.0, Fx −43.43 vs −100.86 |
| T4 one pair per node | node == pair bit for bit |

## 4. Verification

**New tests** (`tests/test_adr157_mortar_pair_friction.py`), 12 cases:

- (a) A flat non-matching patch (3×3 on 2×2, 2×2 on 3×3, and a graded mesh) under consistent
  shear and normal loads, stick then slip, for cohesion, Coulomb, and Coulomb + c capped by τmax.
  The slide equals `(Q − cap)/k` to 1e-6. It passes on the base too, because the pairs agree under a
  uniform field; it pins the analytic patch and alignment independence.
- (b) A creased roof, cohesion only, pressed down in three displacement-driven steps. The ridge
  nodes are shared across the two flanks. The contact force equals the analytic value
  `Fz = 2A(p cos α + c sin α)` with `Fx = 0`: the new binary is within 1e-7, while d63f49750 gives
  `Fx = 0.43` against a cohesion force of 7.3, with Fz off by 3e-4 from step 2 on. A split-ridge
  model agrees with the shared one.
- (c) The roof under force control, cohesion past the cap, in 4 steps. The new binary converges in
  every step and the ridge stays on the symmetry plane to 1e-11. d63f49750 drifts off the plane at
  step 2 and fails to converge at step 3.

Result: 12/12 on the new binary; 10/12 on d63f49750 (the two tests above fail).

**Contact battery** (`tests/*{contact,mortar,adr39,adr41,adr57,adr85,adr96,ladrunoTie}*`, 60
files): see the PR body for the counts on d63f49750 and on this head.

**Pile reproducers (R1):** see the PR body and the piles-validation `ladder/R0_6_fork_slice/BUILDER.md`.

## 5. Descoped

- The geometric tangent terms (∂D/∂u, ∂M/∂u, ∂n/∂u) stay deferred, as in C2/C3.
- No transfer of path state between pairs in finite sliding (D4).
- The 2D mortar lane uses the same key; no new 2D test was added (its battery still passes).
