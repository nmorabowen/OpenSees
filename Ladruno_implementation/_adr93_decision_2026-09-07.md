---
title: "ADR 93 — decision memo (TIMs proposed-model request F2, 2026-09-07)"
project: Ladruno
type: decision-memo
status: PROPOSED — for the owner to confirm and fold into 93_ladruno_sanisand_zero_confinement_adr.md (on wp/93-sanisand-zero-confinement, #801)
owner: nmora
related:
  - "[[93_ladruno_sanisand_zero_confinement_adr]]"
  - "[[_tims_proposed_model_requests_2026-09-07]]"
  - "[[86_ladruno_sanisand_adr]]"
tags: [sanisand, zero-confinement, decision, tims]
---

# ADR 93 — decision memo

> [!note] Why a memo and not an edit of the ADR. ADR 93 exists only on `wp/93-sanisand-zero-confinement`
> (draft #801), a live work package with commits four hours old; writing into another lane's branch is the
> duplicate-lane trap (`LEDGER_quirks`). This memo is the decision text, ready to paste as the ADR's §6 and
> to flip its `status:` from BRAINSTORM, once the owner says so. Nothing here is implemented.

## 6. Decision (2026-09-07, TIMs request F2 — proposed, awaiting the owner's word)

**Context that closed the question.** The TIMs proposed model (PM-01, `work/ape/proposed-model`)
cures the zero-confinement ring **by physics on its side**: a minimum embedment
`d_Ef ≥ 0.4 m` as a dry surcharge outside the footprint (7.7 kPa at the 0.5 m default, against
the 5–6 kPa the response-curve act measured the wall to clear at) plus the footing's own weight
(12.3 kPa) — PM-01 D20 / D8. Their order if that is not enough: raise `d_Ef` first, then this
ADR's II.1. A strength floor, a residual pressure (II.4) and a material crust (III.1) are
**excluded on their side with measured reasons** (PM-01 §17.5). Regularisation is closed
(`90_…_tims_report.md` §6, accepted 2026-09-07). That answers §5's first question: the campaign's
question is `N_γ` on an idealised half-space; the ring is *paid for*, not *stiffened*, and the
BVP-side candidates (III) are off the table for the campaign.

**D1 — Candidate II.1 (elastic-only floor `p_r,e`) is the fork's fallback of record.** It is the
one candidate that puts a floor under the stiffness without adding strength or cohesion, it is
calibratable separately from everything the strength calibration used, and it can be gated for
capacity-neutrality. Every other material-side candidate is closed for the campaign: II.4 stays a
control arm only; II.3 waits on I.1's clamp census; II.2 (D5a) is folded into the P0 whose owner is
the material's author.

**D2 — Not implemented until the owner says go.** TIMs do not need it today; building it now
would put an un-measured parameter on the wire of a calibrated material (constraint 1). When the
owner says go, the change is exactly the one §II.1 already specifies:

- the three `GetElasticModuli` overloads (`ManzariDafalias.cpp:4834-4835, 4876-4877, 4896-4897`
  as of `bc63a388e`): `sqrt(pn / P_atm)` → `sqrt((pn + p_r,e) / P_atm)` in `G` (and through it `K`),
  **nowhere else** — the plastic-side `m_Presidual` untouched, `p_min` untouched;
- one parser flag on `LadrunoSANISAND` (`-Pelastic p_r,e`, default 0 ⇒ byte-identical), echoed at
  construction like `-Presidual`, carried on the fork wire (`Vector(4)` → `Vector(5)`), copied by
  `getCopy`, survives `revertToStart` (the three fingerprints `test_ladruno_sanisand.py` already
  has for `p_r`);
- the capacity-neutrality gate: on a one-element drained triaxial at `p_0 = 0.01 P_atm` and on the
  ADR's strip control (`sanisand_tau0_band.py`), the peak resultant moves by **less than 1 %**
  (stated tolerance; the owner may tighten it) when `p_r,e` is halved from its chosen value —
  the floor is allowed to change the *stiffness* of the ring and forbidden to change its
  *strength*;
- `Elastic2Plastic` (`:5113-5136`) re-checked with the floor present: the sub-`p_min` reset and
  the `M_c` raise both test the **stress**, which the floor does not touch, so the stage switch
  behaves as today — to be confirmed by the existing zero-count gate on the confined deck, not
  assumed.

**D3 — Vanilla footprint.** The three sites are in `ManzariDafalias.cpp`; either a
`// Ladruno` marked edit with a vanilla ledger row (the ADR-86 PR-2 pattern already on those
lines) or a virtual override in the subclass — the former, because the three overloads are not
virtual today and making them so is itself a vanilla edit.

**Open for the owner (unchanged from §5):** whether `p_r,e` is acceptable to the material's author
as a fork parameter (default 0, echoed) or must stay a deck-level flag; who owns D5a.

**Status after this section: DECIDED (II.1 as fallback, not built) — pending the owner's
confirmation to flip `status:` above.**
