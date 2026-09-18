# WP-114 — BezierTri6 `-bbar` under non-associated flow (TIMs F17)

Status: **scoping** (2026-09-18). Branch `wp/114-beziertri6-bbar-nonassoc`, cut from
`ladruno` @ `48c0e99bc`. Source request: TIMs Workbench `2d-model` act,
`fork_request_bezier_2026-09-18.md`. We reproduce on our own test beds and never run
the Workbench models.

## Claims to test (from the act, not yet verified here)

- Strip footing, band-refined tri6 field, `-bbar`. UW DruckerPrager φ=33°: associated flow
  reaches a clean limit of 465.9 kPa. At ψ=0 it walls at s/B=0.0022, with 123 GPs at the
  cutoff or apex. On the same field, `LadrunoQuad -bbar` reaches a peak and a plateau with 0
  cutoff/apex points.
- LadrunoSANISAND walls at s/B ≤ 0.0006 on the Bézier mesh, but not on the quad mesh.
- **Every Pardiso factorisation perturbs 1 816 pivots on the Bézier field**, from the
  elastic K0 stage onward. The quad field shows none.

## Asks → phases

| # | Ask | Phase | Depends on |
|---|-----|-------|-----------|
| 1 | Identify the DOFs behind the 1 816 perturbed pivots (bisect: tie, surface mid-edge CPs, corner elements, bbar off) | A | — |
| 2 | ψ=0 vs ψ=φ punch patch; numerical vs analytic element tangent at a plastic, non-associated state | A (tangent probe) → B (punch) | — |
| 3 | Why cutoff/apex points appear only on the Bézier elements | B | 1, 2 |
| 4 | Verdict: a fix behind `-bbar`, or a documented limitation in the element guide | C | 1–3 |

## Leads found during scoping (hypotheses)

- `BezierTri6::computeBBarMatrix` uses the **3D** B-bar split (÷3 terms) on a 3-row
  plane-strain B. The ε_zz row, (B̄−B)/3, is dropped. So at a GP, ε_xx+ε_yy =
  div u + (2/3)(avg − local): the dilatation is only partly projected, and the ε_zz
  the material sees is 0, not the projected value. Compare this with LadrunoQuad's
  plane-strain B-bar.

## Deliverables

- A reproducible pivot-bisection script.
- A numerical-vs-analytic tangent test at a plastic, non-associated state.
- A non-associated punch test on BezierTri6 `-bbar`. It either reaches the quad plateau,
  or the limitation is documented.
- Ledgers: a quirks row, an implementations row if we fix something, and the guide
  (`04_bezier_elements.md`).
