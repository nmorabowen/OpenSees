---
wp: LEGACY
title: "ops.printA('-ret') returns the raw OpenSees Matrix buffer, which is COLUMN-major — read row-major, an unsymmetric tangent looks transposed and FD-vs-tangent te…"
legacy_seq: 171
---
### `ops.printA('-ret')` returns the raw OpenSees Matrix buffer, which is COLUMN-major — read row-major, an unsymmetric tangent looks transposed and FD-vs-tangent tests false-fail at exactly the asymmetry magnitude
- **Bites:** any test/tool that reshapes `printA('-ret')` into `(neq, neq)` C-order and compares against an oracle or FD residual. On symmetric tangents the bug is invisible; on LadrunoUP's unsymmetric [K,−Q;0,H] the −Q block appears in the p-row/u-col slot and the check fails at |Q|/|K| (measured 2.6e-5 — small enough to chase as a "tolerance problem" for hours).
- **Why:** `OpenSeesCommands.cpp:2590` hands back `&A(0,0)` flat; OpenSees `Matrix` storage is column-major (Fortran order).
- **Workaround:** reshape Fortran-order or transpose after reshape (`np.array(ret).reshape(neq, neq, order='F')`). Pinned in `tests/test_ladruno_up_element_equiv.py` helper (ADR-71 P1, 2026-07-11).
