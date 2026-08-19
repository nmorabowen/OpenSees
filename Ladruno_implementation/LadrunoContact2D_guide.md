---
title: LadrunoContact — 2D plane-model contact user guide
project: Ladruno
tags:
  - user-guide
  - contact
  - 2d
---

# LadrunoContact — 2D plane-model contact user guide

Companion to [[85_ladruno_contact_2d_adr]] (the design record) and the
existing 3D `LadrunoContact` documentation. This page covers only what is
specific to `-ndm 2` decks. All lanes shipped across ADR-85 T0–T4: rigid
plane, NTS penalty (with friction), mortar/ALM (with friction and tie), and
the D4 radial end-cap. The 2D contact lane is **complete**. ADR-85 **F1**
added one follow-on orientation mode on the NTS lane, `-outward winding`
(see "Orientation" below).

## Command surface (2D forms)

```tcl
contactSurface <masterTag> -master 2 <n0> <n1> <n1> <n2> ...
contactSurface <slaveTag>  -slave  <n0> <n1> ...
contactSurface <slaveTag>  -slave-segments 2 <n0> <n1> <n1> <n2> ...

contact <tag> <masterTag> <slaveTag> <kn>|auto <kt> <mu> \
    [-outward <ox> <oy> | -outward winding] [-consistanttan] [-geomtan] [-visc <muc>]

contact <tag> <masterTag> <slaveTag> "-mortar" -epsN <val>|auto \
    [-thickness <h>] [-tie] [-outward <ox> <oy>] ...

contactPlane <tag> <slaveSurfTag> <nx> <ny> <px> <py> <kn> [opts]   # 7-arg 2D form
contactPlane <tag> <slaveSurfTag> <nx> <ny> 0 <px> <py> 0 <kn>      # legacy 9-arg form, still valid
```

`-master`/`-slave-segments` take `nodesPerSeg = 2` (a 2-node line segment,
the 2D analogue of the 3D 3/4-node facet). `-slave` (a plain node set) is
also accepted for the NTS lane, exactly as in 3D. `ndf == ndm == 2` exactly
is required on every NTS/mortar/tie surface (the rigid-plane lane alone
keeps `ndf >= ndm`, for 3D-shell-on-2D-plane style decks); a mismatch is a
named, loud abort — never a silent DOF-shift.

## **THE CHAINED STRIDE-2 PAIR-LIST CONVENTION — READ THIS BEFORE DECLARING A
## MULTI-SEGMENT MASTER OR SLAVE-SEGMENTS SURFACE**

**`-master 2` and `-slave-segments 2` take a FLAT STRIDE-2 PAIR LIST, chained
head-to-tail — NOT a node chain.** Three segments need **six** tags, listed
as consecutive pairs:

```tcl
contactSurface 10 -master 2  101 102  102 103  103 104     ;# CORRECT: 3 segments
```

The natural-looking shorthand `-master 2 101 102 103 104` (four tags) is
**silently legal** and declares **two disjoint segments** `(101,102)` and
`(103,104)` with a **hole** where `(102,103)` should be. The parser cannot
object — an even tag count with no repeated node is indistinguishable from
a genuinely disjoint declaration (which is also legitimate). The result:
the deck **converges, balances its reactions, and transmits the load
through the wrong distribution** — the ADR-78 "converged, balanced, wrong"
failure shape, not a crash. This was hit during T3 gate authoring (measured:
master row forces `[5/18, 8/18, 5/18]*P` instead of the intended
`[1/4, 1/2, 1/4]*P` — the *exact* discrete solution of the *holed* surface,
not a numerical error).

Build multi-segment lists programmatically, e.g.:

```python
def seg_pairs(chain):        # chain = [n0, n1, n2, ..., nk] in surface order
    tags = []
    for a, b in zip(chain[:-1], chain[1:]):
        tags += [a, b]
    return tags

ops.contactSurface(10, "-master", 2, *seg_pairs([101, 102, 103, 104]))
```

Surface order is also the adjacency contract the concave-vertex primitive
and the D4 end-cap key off: `handle()` runs an O(nSeg) chain-integrity scan
that FATALs by name on a permuted or reversed listing (re-order/re-wind
guidance included in the message) — you cannot silently reach the holed
shape above through THAT specific mistake; the hole above is legal precisely
because it declares two genuinely separate segments, which the scan cannot
(and must not) refuse.

**One exception, and it is the ONLY place the holed shape is refused:**
under `-outward winding` (ADR-85 F1, below) the master must form one
connected chain, so the four-tag declaration above draws a named FATAL. That
is not the scan getting smarter — it is a separate requirement that exists
because winding bypasses the orientation vote, which is what otherwise
catches a second run wound the wrong way.

## Vertex / corner policy

At an **interior** shared vertex between two consecutive (chained) segments,
the 2D lane classifies the corner from the reference-configuration
orientation:

- **Convex** corners (material fans outward) are covered twice by their two
  adjacent segments and **never carry force** — no vertex pair is
  instantiated.
- **Concave** (re-entrant) corners are a real coverage hole — both adjacent
  segments refuse a penetrating slave there — and DO get a vertex pair
  (radial normal from the corner point, committed side sign at first
  capture, never re-derived while interpenetrating).

At an **open terminal** (a chain end with no neighbouring segment on that
side), T4 instantiates a **radial end-cap** vertex pair — see "The D4
end-cap" below.

## The D4 end-cap (ADR-85 T4)

Every open terminal of a declared 2D master chain (both ends of a
stand-alone segment; each unmatched end of a longer chain) is automatically
capped by a radial vertex pair, replacing the retired `NTS2D_END_SLACK`
parametric window. No user action is required to enable it. Practical
notes:

- The cap activates **only within a bounded tangential reach** of the
  terminal (a small, fixed fraction of the adjacent segment's own
  parametric length) — it does **not** reach out and claim slaves that are
  genuinely elsewhere on the model (e.g. a slave that used to be covered by
  a segment that was later removed via `remove element`/`remove node`).
- The cap is **C0** with the adjoining segment at the terminal itself: no
  window, no force-discontinuity cliff for a boundary slave that drifts
  slightly past the nominal end under load (the scenario this replaces
  `NTS2D_END_SLACK` for).
- Friction on an end-cap pair, if `-mu > 0` is declared on the contact, uses
  the same scalar return map as every other 2D NTS pair.

## Orientation: the centroid vote, and `-outward winding`

The 2D NTS and mortar lanes derive their orientation from an
**interface-level reference vote** computed once at `handle()` from the
reference-configuration surface centroids (**not** the shipped 3D default's
per-pair, ridge-flip-prone datum). If the master and slave surfaces are
**flush** (coincident or nearly coincident centroids — e.g. a zero-gap
masonry joint or a footing seated directly on soil with no initial gap),
the centroid vote is genuinely ambiguous and the deck draws a **named
abort** rather than guessing. A vote that resolves for *some* master
segments but not others (a strongly curved or inconsistently wound master —
a closed loop is the extreme case) is also a named refusal; the engine will
not guess a per-segment sign.

There are therefore **two** ways to orient a 2D interface. Pick one — giving
both is a named refusal.

### 1. `-outward ox oy` — declare a direction

A unit-ish 2-component direction toward the slave's allowed half-space. The
ONLY vector form accepted in 2D; the legacy 3-component 3D form is rejected
on a 2D surface. Works on **both** the NTS and the mortar lanes.

```tcl
contact 1 10 20 $kn 0.0 0.0 -outward 0.0 1.0
```

### 2. `-outward winding` — let the chain declare the side (**NTS only**, ADR-85 F1)

```tcl
contactSurface 10 -master 2  101 102  102 103  103 104
contact 1 10 20 $kn 0.0 0.0 -outward winding
```

The centroid vote is **bypassed**: `sigma` is fixed at `+1`, so every master
segment's normal is `perp(t) = (-t_y, t_x)` of **its own travel direction**.
Read it as: **the slave lies to the LEFT of the chain's traversal.**

Because the sign is per-segment-tangent rather than one global direction,
this is what unlocks the two decks the vote refuses:

- **flush interfaces** (masonry joint, footing on soil, any zero-gap deck) —
  no centroid separation is needed, because none is consulted;
- **closed-loop and strongly curved masters** (a ring, a full indenter
  profile) — every segment gets its own perpendicular, so there is nothing
  to split. List a closed loop **clockwise** for outward-facing normals,
  **counter-clockwise** for inward-facing ones.

**The cost, stated plainly: winding turns a wrong master into a converged
wrong answer.** The centroid vote is currently what catches a master
boundary that faces *away* from its slave; bypassing it, that deck runs and
transmits nothing (or transmits on the wrong side) instead of aborting. The
engine keeps the one check it can still make — see the connectivity
requirement below — but "is this the face I meant?" is now the caller's.

**Winding requires ONE CONNECTED CHAIN.** A **disjoint** master (the
four-tag `-master 2 101 102 103 104` shape from the pair-list section above:
two separate runs, legal, and invisible to the chain-integrity scan because
every node is used exactly once) is a **named FATAL** under `-outward
winding`. It has to be: with the vote bypassed, a second run wound the other
way would be a silent wrong-side normal rather than the split-vote refusal
that catches it today. Declare one `contact` per connected run, or re-wind
the runs into a single head-to-tail chain. Open chains, single segments and
closed loops all satisfy it.

**Winding is refused on the `-mortar` lane, deliberately.** The
chain-integrity scan is an NTS-lane construct; the 2D mortar lane runs no
such scan, so the connectivity invariant winding leans on is not established
there, and hoisting the scan would newly FATAL permuted or reversed mortar
masters that ship accepted today. The consequence is an asymmetry worth
stating twice: **flush 2D MORTAR interfaces still abort, and mortar users
still pass `-outward ox oy`.**

### Closed loops: give the ring ENOUGH SEGMENTS, or it transmits nothing

Winding makes a closed-loop master declarable for the first time, which also
makes a pre-existing property of the NTS ownership rule reachable for the
first time. **A ring whose facets are coarse relative to its own diameter
converges, balances its reactions, and transmits EXACTLY ZERO.** Measured for
a regular n-gon of circumradius 1 with a slave seeded just inside one facet,
at the default `-cell 1.0`: **n ≤ 8 inert, n ≥ 9 exact.** Forcing a single
broad-phase bucket (a very large `-cell`) reproduces the inert result at any
n, which is the tell.

Why, in one sentence each: the broad-phase grid is capped at `nSeg` cells
with a ±1-cell search, so a coarse ring cannot be separated from its own far
side; the far-side facet then arms with a body-diameter "penetration" and
kicks the slave toward the ring's interior; and from the interior every facet
projects the slave in bounds, so the ordered-ownership rule ("stand down if
your predecessor is in bounds") has no starting point on a closed chain and
the deferrals cycle all the way around. An open chain cannot do this — its
first segment has no predecessor and always arms.

Practical rule: **size a closed loop's facets from the loop's own extent,
not from the surrounding elastic mesh** — the same rule as the curved-master
section below, now with a second reason. If a ring transmits nothing, refine
it (or check `printA`: a spring-only tangent with no `kn` in it is the
signature). Full write-up in `LEDGER_quirks.md`; the limitation is pinned by
`tests/test_adr85_contact2d_t5_winding.py::test_coarse_closed_loop_transmits_nothing`.

## 2D units / thickness table

The 2D lane has **two independent thickness conventions** — do not conflate
them:

| Quantity | Where it lives | Thickness handling |
|---|---|---|
| Plane element thickness (`"quad"`/`"tri"` etc., the `THICK` argument) | The **element** declaration | Baked into the element's own stiffness; contact never re-reads or re-derives it |
| NTS `kn`/`kt` (penalty) | `contact ... <kn> <kt> <mu>` | Force **per unit length of gap** — as in 3D, `kn` is *never* a pressure. `-kn auto` resolves from the owning element's `getInitialStiff()`, which already folds the element's own thickness — **no `-thickness` parameter exists on the NTS lane at all.** |
| Mortar `-thickness h` (default **1.0**) | `contact ... "-mortar" ... -thickness <h>` | Interval integrals on the mortar lane naturally produce force **per unit thickness**. `h` is applied **once**, at the 2D injection site, to: `epsN`, `epsT`, the viscous `-visc` coefficient, the friction clamps (`-tauMax`/`-cohesion`), and the tie stiffness. |
| `-epsN auto` under mortar `-thickness h != 1` | mortar, auto-resolved | The **auto-resolved value is NOT h-scaled** (it already absorbs the owning element's thickness via `getInitialStiff()`) — re-scaling it by `h` would be an h² error. Regression-gated: `tests/test_adr85_contact2d_t3_thickness.py::test_2d_epsN_auto_thickness_no_double_count`. |
| `lambda_N` (mortar ALM multiplier) | internal, queried via `ladrunoMortarPenetration` etc. | **Never separately scaled** — it inherits `h` through `epsN*gbar` in the Uzawa update. Scaling it again is the same h² error class. |
| `ladrunoContactForce(node)` | query | **NTS lane only**, fed exclusively from the SEGMENT/end-cap branch. Reports the SCALAR normal traction magnitude `tn`, not a force *vector* — if the active pair's normal is not purely axis-aligned (as it generally is not near a corner or the D4 end-cap), do not assume this equals any single global force component. |
| `ladrunoMortarPenetration()` | query | Dimension-blind (a length); unaffected by `-thickness`. |

**A unit-thickness convention for benchmarking:** when comparing a 2D
plane-strain deck against a force-per-unit-thickness analytic reference
(e.g. the T4 Hertz benchmark), declare elements with thickness `1.0` — the
summed FE contact force then equals the reference force-per-unit-thickness
numerically, with no conversion factor anywhere in the comparison.

## Curved masters: facet length must track the expected penetration, not the mesh

On a **curved** 2D master (a faceted arc, e.g. a rigid indenter profile),
the NTS narrow phase stops arming a pair once its penetration exceeds
roughly **twice the local master facet length** (an interior-seam
consequence of the ordered-ownership rule, not a bug). A deck that seeds a
slave with a large initial penetration relative to a fine facet mesh will
therefore silently **disarm the interior of the contact patch**, collapsing
the transmitted load onto a rim of a few surviving nodes — this looks like
ordinary discretization error and is not. If your master arc's facet length
is driven by the SAME mesh parameter as the surrounding elastic mesh (the
natural-looking choice), refining that mesh can make this *worse*, not
better. Size the facet length from the **expected penetration depth**
instead, independent of the elastic mesh's own resolution — see
`tests/test_adr85_contact2d_t4_hertz.py`'s `hm` sizing for a worked example.

## Known coverage gaps (demand-driven, not blocking)

- **`WARN_VTX2D_SIDE_FLIP`** — REPRODUCED for the first time in T4 (four
  prior attempts across T1b-T3 on the CONCAVE-vertex path failed; the
  latch is dimension-generic, keyed by vertex node tag, and fires
  identically for the D4 end-cap's single-tangent commit). Rig: an
  open-terminal end-cap slave offset both tangentially and radially from
  the terminal (so the radial normal starts with a non-degenerate y
  component), pushed to commit a side sign, then load-reversed hard enough
  to cross the adjoining segment's extended line while staying inside the
  reach bound -- `committed -1, live 1` fires reliably on the very first
  reversal iterate. The rig as constructed does not yet reach a converged
  equilibrium past the flip (Newton stalls after the warning fires, which
  is itself informative -- a committed-vs-live disagreement is exactly the
  ADR-57 "attracts when it should repel" hazard the latch exists to catch
  loudly, not a rig bug to engineer around blindly). **Left as a T4+
  follow-up:** tune the rig (softer approach, smaller reversal increments,
  or a stiffer stabilizing spring) to a converged, permanent regression
  test; the reproduction itself is the debt this item existed to close.
- **2D axisymmetric contact** and **parallel/DDM 2D contact** are explicitly
  out of scope (see the ADR's "Not in scope" section) — both are
  formulation changes, not parameterizations of the shipped 2D lane.
- **Coarse CLOSED LOOPS transmit nothing** (ADR-85 F1) — see "Closed loops:
  give the ring ENOUGH SEGMENTS" under Orientation above. Disclosed and
  measured, not fixed: the fix touches the shared 2D ordered-ownership
  precedence rule, which would move every existing 2D deck's force
  distribution and needs its own gate.
- **`-outward winding` does not check that the master FACES its slave**
  (ADR-85 F1) — the centroid vote was silently doing that job as well, and
  winding bypasses the vote. A backwards-wound master runs and gives a
  converged wrong answer. Whatever generates the deck should own that check.
