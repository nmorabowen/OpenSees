---
wp: WP-207
pr: "#1015"
title: "responseSpectrumAnalysis -mode k wipes the gravity state (replace + commit)"
date: 2026-10-22
---
### `responseSpectrumAnalysis -mode k` wipes the gravity state (axial forces go to 0, later `eigen` loses P-Delta)

- **Bites:** gravity, `loadConst`, `eigen`, `modalProperties`, then per-mode
  `responseSpectrumAnalysis ... -mode k`. After the first mode the committed
  displacement is exactly the modal one, the column axial force is 0, a P-Delta
  element's mode force lacks its `N_g·u_k` term, and a second `eigen` returns the
  unsoftened eigenvalues (3-storey PDelta column: λ 32.28 → 39.60).
- **Why:** stock `solveMode()` calls `node->setTrialDisp(u_modal, i)` with the
  *absolute* modal displacement and `endMode()` commits it, so there is nothing
  for `revertToLastCommit` to return to. Per-mode displacements (from the
  eigenvectors) are still right; per-mode forces of geometric-stiffness elements
  are first-order only.
- **Fix/status:** use `-increment` (WP-207): trial = committed + sign·u_k, update,
  record, no commit; `responseSpectrumAnalysis -revert` to return to gravity.
  Without it, re-run gravity before any later `eigen`. Learned 2026-10-22.
