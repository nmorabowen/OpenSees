---
wp: LEGACY
title: "analyze() returns rc=0 on a NaN-poisoned system — \"rc==0\" does NOT mean \"numbers\""
legacy_seq: 200
---
### `analyze()` returns rc=0 on a NaN-poisoned system — "rc==0" does NOT mean "numbers"
- **Bites:** trusting the analyze return code as a health signal. `FullGeneral` + `algorithm Linear` solved a stiffness matrix full of NaN (degenerate-eas 1/det ≈ 1e17 blowup) and reported SUCCESS; `nodeDisp` was NaN with no error printed at any level.
- **Why:** LAPACK `dgesv` propagates NaN without setting its info flag, and no layer above it (SOE / algorithm / analysis) checks the solution for finiteness.
- **Consequence:** element-level degeneracy/finiteness guards are the ONLY defense against silent-NaN results in linear analyses. Tests asserting "the analysis must fail" must either construct a genuinely SINGULAR system (zeroed row → dgesv info>0 → rc!=0) or assert output finiteness explicitly — never rely on NaN tripping the solver.
- **Status (2026-07-20):** found while probing the eas degeneracy-guard axis-collapse hole; the guard fix makes eas elements refuse loudly (their zeroed block → singular SOE → honest rc!=0), but the generic solver blind spot remains (upstream-class, unfixed).
