---
wp: LEGACY
title: "Truss (and most elements) default to a LUMPED mass — -lump diagonal ≡ -lump rowsum on them"
legacy_seq: 20
---
### `Truss` (and most elements) default to a LUMPED mass — `-lump diagonal` ≡ `-lump rowsum` on them
- **Bites:** the ExplicitBathe / CriticalTimeStep `-lump <rowsum|diagonal>` option
  only changes `dt_cr` when the element's `getMass()` returns a **consistent**
  (non-diagonal) matrix. The `Truss` element defaults to `cMass=0` (lumped:
  `0.5*rho*L` on the diagonal, zero off-diagonals). For a diagonal matrix,
  "diagonal-of-consistent" (take `M_ii`) and "row-sum" (sum each row onto the
  diagonal) are identical, so both lumpings return the same `dt_cr` (= `le/c` for
  a 2-node bar). A test that asserts `dt_diag < dt_rowsum` on a default lumped bar
  is **vacuous and fails** (the two tie).
- **Why:** `Truss::getMass()` returns lumped mass unless built with `-cMass 1`;
  the consistent branch is `rho*A*Le/6 * [[2,1],[1,2]]`. Only then do the lumpings
  diverge: diagonal-of-consistent = `rho*A*Le/3`, row-sum = `rho*A*Le/2`, giving
  `dt_diag/dt_rowsum = sqrt((Le/3)/(Le/2)) = sqrt(2/3) ≈ 0.816` (per-element pencil).
- **How it surfaced:** `tests/test_explicitBathe_integrator.py::test_eb_lump_diagonal`
  failed on Zone-A (Ubuntu) CI with `diagonal=0.02000 vs rowsum=0.02000` (a tie at
  `le/c`). The Linux result is **correct**; the source `_verify_explicit.py` test 11
  had the same faulty premise (built a lumped bar yet asserted the inequality).
- **Workaround/status:** build the bar with `-cMass 1` (consistent mass) for any
  test that means to exercise the diagonal-vs-rowsum *difference*. Fixed in the
  Zone-A port (PR #40) and mirrored back into `Ladruno_scripts/_verify_explicit.py`.
  *(Learned 2026-05-31, Zone-A ExplicitBathe battery.)*
