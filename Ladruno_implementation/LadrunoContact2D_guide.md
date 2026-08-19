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
the D4 radial end-cap. The 2D contact lane is **complete**.

## Command surface (2D forms)

```tcl
contactSurface <masterTag> -master 2 <n0> <n1> <n1> <n2> ...
contactSurface <slaveTag>  -slave  <n0> <n1> ...
contactSurface <slaveTag>  -slave-segments 2 <n0> <n1> <n1> <n2> ...

contact <tag> <masterTag> <slaveTag> <kn>|auto <kt> <mu> \
    [-outward <ox> <oy>] [-consistanttan] [-geomtan] [-visc <muc>]

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

## Flush interfaces require `-outward`

The 2D NTS and mortar lanes derive their orientation from an
**interface-level reference vote** computed once at `handle()` from the
reference-configuration surface centroids (**not** the shipped 3D default's
per-pair, ridge-flip-prone datum). If the master and slave surfaces are
**flush** (coincident or nearly coincident centroids — e.g. a zero-gap
masonry joint or a footing seated directly on soil with no initial gap),
the centroid vote is genuinely ambiguous and the deck draws a **named
abort** rather than guessing. **Flush 2D interfaces must supply `-outward
ox oy`** (a unit-ish 2-component direction, the ONLY form accepted in 2D —
the legacy 3-component 3D form is rejected on a 2D surface):

```tcl
contact 1 10 20 $kn 0.0 0.0 -outward 0.0 1.0
```

A vote that resolves for *some* master segments but not others (a strongly
curved or inconsistently wound master) is also a named refusal — split the
surface into separate contacts or re-wind it; the engine will not guess a
per-segment sign.

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
