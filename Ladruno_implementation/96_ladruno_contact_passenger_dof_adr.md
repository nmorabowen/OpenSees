---
title: ADR 96 — Contact and ZeroLength on ndf ≥ 3 nodes (the pressure DOF as a passenger)
project: Ladruno
status: accepted
priority: high
owner: nmora
tags:
  - implementation
  - contact
  - element
  - up
related:
  - "[[78_ladruno_up_corot_adr]]"
  - "[[47_ladruno_contact_deferrals_adr]]"
  - "[[71_ladruno_up_family_adr]]"
  - "[[39_ladruno_contact_domain_adr]]"
  - "[[_tims_proposed_model_requests_2026-09-07]]"
---

# ADR 96 — Contact and ZeroLength on ndf ≥ 3 nodes

## What

The 3-D contact lanes (rigid plane, NTS, mortar/ALM, edge-edge) and the upstream
`ZeroLength` element accept nodes with **more than three DOFs**. Each node's
**first three DOFs are the translations** (the fork's ndf convention, ADR 71) and are
the only ones the contact adapter or the spring reads, writes or couples; the
fourth DOF of a `LadrunoUP` node (pore pressure) and the rotations of an ndf-6
beam or shell node **ride as passengers** — untouched. **This is NOT ADR 47
deferral 9.** There is no gap flow, no pressure penetration, no thermal or
hydraulic coupling across the interface and no coupled contact physics of any
kind; deferral 9 stays deferred and must not be read as having covered, or as
being covered by, this plumbing item. It is exactly the sibling item ADR 78
§4 recorded (`78_ladruno_up_corot_adr.md:277-290`), sequenced now because the
TIMs proposed model needs an openable interface between an ndf-3 footing skin
and ndf-4 saturated sand (request `_tims_proposed_model_requests_2026-09-07.md`,
F1). The 2-D contact lane is **not touched**: its `ndf == ndm == 2` equality gate
(ADR 85, `LadrunoContactHandler.cpp:153-181`) stays as it is; a 2-D u-p node is
a separate decision with its own adoption record in apeGmsh.

## Why

Two independent guards refused every mixed-ndf interface:

1. **Contact.** `FE_Element::setID()` (`SRC/analysis/fe_ele/FE_Element.cpp`) copies
   *every* equation of every connected `DOF_Group` into `myID` and returns `-3`
   once it runs past `numDOF`. On an ndf-4 node the fork's 3-per-node adapter
   layout (`LadrunoContactFE`, `ndof = 3·(1+n_ps)`) therefore received the pore
   pressure's equation in a translation slot and aborted. The handler pre-empted
   that with six `getNumberDOF() != 3` FATALs (ADR 78 P1) so the failure was at
   least loud — but it made every u-p contact model impossible.
2. **ZeroLength.** `ZeroLength::setDomain()` refuses nodes of differing ndf
   ("differing dof at ends", warn and return: the element is silently absent),
   and dispatches its element size on `(dimension, ndf)` pairs, so a (3,4) pair
   had no formulation at all.

The gate on TIMs' side (PM-01): an openable base under pure vertical compression
must be inert — the summed normal traction equals the applied load, no negative
contact pressure, and the pore-pressure field identical to the same column bonded
by `equalDOF` on DOFs 1–3.

## Where

- Modify (fork): `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — `setID()` override.
- Modify (fork): `SRC/analysis/handler/LadrunoContactHandler.cpp` — the six 3-D
  guards (`:1825, 1857, 2335, 2358, 2563, 2575` at `bc63a388e`).
- Modify (**upstream**, `// Ladruno (ADR-96)` marked, vanilla ledger rows):
  `SRC/element/zeroLength/ZeroLength.{h,cpp}` — passenger mode.
- Not touched: `LadrunoContactProjection.h` (kernel-local `3·(1+n_ps)` scratch
  sizes remain correct: the adapter still has exactly three slots per node),
  `LadrunoContactHandler.cpp:153-181` (2-D gate), `ZeroLength.cpp` parser.
- Tests: `tests/test_adr96_passenger_dof.py`; the byte-identity harness
  `tests/_testbed/contact_dump.py`; the existing ADR-39/41/57/85 battery.
- Docs: this ADR; `ndf_and_mixed_models_guide.md`;
  `zerolength_and_link_springs_guide.md`; `LadrunoContact2D_guide.md` (one note).

## How — decisions

**D1. Translations first, the rest passengers.** On every node a contact
adapter or a passenger-mode ZeroLength touches, DOFs 1–3 are the translations.
That is the fork's ndf convention (ADR 71: `LadrunoUP` is `[u_x u_y u_z p]`,
beams/shells `[u | θ]`); nothing here re-derives it from the element type. A
node with ndf < 3 in a 3-D lane is still refused.

**D2. The contact fix is a `setID()` override, not a handler-side ID map.**
`LadrunoContactFE::setID()` walks the adapter's own ordered node list (the
layout each constructor declares: `[slave | seg_1..n]`, `[slave facets | master
facets]`, `[sa sb ma mb]`, `[slave]`) and copies the first `ndm` equations of each
node's `DOF_Group`. `FE_Element::numDOF` and `theModel` are private, so the
override sizes itself from `myID.Size()` and reaches the groups through
`Node::getDOF_GroupPtr()`. Layout mismatches and groups with fewer than `ndm`
equations are FATAL with a message; `EMPTY` mode returns 0 (no connectivity);
any other mode falls through to the base. The residual/tangent/B-operator code
in the adapter is untouched — it never assumed anything about the node beyond
"its first three DOFs are the translations".

**D3. Handler guards relax from `!= 3` to `< 3`.** The two FATAL sites (NTS slave,
mortar slave) keep their messages, reworded to "< 3"; the four silent
`ok = false; continue` master-side skips (NTS master, mortar master, edge-edge
slave and master) become `< 3` as well, so an ndf-4 master facet is no longer
dropped without a word. No new message where there was none — that is an ADR 78
P1 item, not this one.

**D4. ZeroLength keeps its 6-slot core and scatters.** In 3-D, when both nodes
have ndf ≥ 3 and the pair is not a vanilla `(3,3)` or `(6,6)` pair, `setDomain`
enters passenger mode: `numDOF = 6`, `elemType = D3N6`, the vanilla `t1d`
transformation and material loop are untouched, and every public accessor
(`getTangentStiff`, `getInitialStiff`, `getDamp`, `getMass`, `getResistingForce`,
`getDampingForce`, `getResistingForceIncInertia`, the three sensitivity
accessors) scatters the 6×6 / 6-vector core into an element of
`dofNd1 + dofNd2` slots with node 2's translations at offset `dofNd1`. Node
displacement and velocity differences are taken over the first three entries
only. `getNumDOF()` reports the element size so `FE_Element` sizes its ID map
from the nodes' own groups — the layout `[node1 all | node2 all]` is exactly what
`Element::getRayleighDampingForces()` assembles, so `-doRayleigh` stays
consistent (the Rayleigh path returns the element-sized base matrix and adds the
element-sized damping forces to the scattered vector). The vanilla path is
byte-identical: passenger mode is a boolean decided in `setDomain`, every other
branch is `numDOFPassenger == 0 ⇒ the old expression`.

**D5. Rotational springs are refused in passenger mode.** A `-dir` of 4..6 has no
slot in the translational core; `setDomain` prints a message naming the element,
the two ndf values and the offending direction, and disables the element
(vanilla's own behaviour for an unsupported pair). A `(6,6)` pair is still the
vanilla `D3N12` element with rotations.

**D6. Wire format unchanged.** `numDOF` is already sent; `numDOFPassenger` and
the scatter targets are re-derived in `setDomain` on the receiving side.

**D7. Element responses stay core-sized.** `force`, `deformation`,
`dampingForces`, `material` responses on a passenger-mode ZeroLength are the
6-slot core (three per node), as on a `(3,3)` pair; a recorder sees the same
columns whether the soil node is ndf-3 or ndf-4.

**D8. What is not claimed.** No coupling, no mass, no pressure term; no 2-D lane;
no `ZeroLengthSection` / `ZeroLengthND` / `TwoNodeLink` (they keep their equal-ends
rule, see the guide); no parallel-contact re-verification beyond the serial
battery (the adapter's `setID` runs on every rank identically).

## Validation gates

- **G1 — bit-identical on ndf-3 models.** `tests/_testbed/contact_dump.py`
  artifact before (`bc63a388e`+F4 build) and after: identical as bytes. The
  ADR-39/41/57/85 pytest battery unchanged.
- **G2 — the TIMs gate.** A `LadrunoUP` column (H8, `-formulation bbar -pOrder
  equal`, drained static path) under an ndf-3 rigid platen in NTS contact on its
  top face, pushed statically: Σ contact normal traction = applied load, no
  negative contact pressure, and the pore-pressure field identical to the same
  column bonded by `equalDOF` on DOFs 1–3.
- **G3 — the spring.** A `zeroLength` with `ENTMaterial` between an ndf-3 and an
  ndf-4 node opens under tension (zero force) and carries compression, with the
  ndf-4 node's fourth DOF unchanged; the same deck on a (4,4) pair and a (3,3)
  pair agrees.

## Risks / open questions

> [!question]
> Should the mortar and edge-edge master-side `< 3` skips become FATAL like the
> slave side? ADR 78 P1 left them silent; this ADR only relaxes the comparison.
> Deferred to the next contact review.

> [!question]
> The passenger-mode ZeroLength returns `Element::getDamp()` directly under
> `-doRayleigh 1`; the vanilla path copies it into the shared static `ZeroLengthM6`.
> Both are correct; only the vanilla copy would have resized the shared static on a
> size mismatch, which is why the passenger path does not copy.

- A `(6,4)` pair (beam on saturated soil) enters passenger mode with a
  translational core; there is no rotational spring on it by D5.
- Parallel contact (ADR 78) builds adapters per rank through the same
  constructors; `setID` is rank-local. Not re-run here (serial battery only).

## Implementation log

- 2026-09-07 — branch `wp/96-contact-passenger-dof` cut from `bc63a388e`.
  Request note citations corrected in the PR (not in the note): the six handler
  guards are at `:1825, 1857, 2335, 2358, 2563, 2575`, not the ADR-78 line numbers
  the note repeats; the `3·(1+n_ps)` sizing lives in `LadrunoContactFE.cpp`, not
  the handler; `ZeroLength.cpp:76, 153-181, 217-260` are parser code — the gate is
  `setDomain()` `:610-672`. The request's "ADR 94" number was taken (ASDPlastic
  review) and 95 is reserved (Prandtl); this is 96.
