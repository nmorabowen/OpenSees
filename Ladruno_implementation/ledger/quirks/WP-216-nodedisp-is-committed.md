---
wp: WP-216
title: "nodeDisp returns the committed displacement, so -increment trial states read zero"
date: 2026-10-27
pr: "#TBD"
---
### `nodeDisp` is the committed displacement: after `responseSpectrumAnalysis -increment` it does not show the mode

- **Bites:** a per-mode Python loop that forms modal inertia (λ M u) or drifts
  from `nodeDisp` after `responseSpectrumAnalysis ... -mode k -increment`
  reads the committed state (zero for a plain RSA), with no error.
  `eleForce` and `reactions` do see the mode, so the numbers disagree silently.
- **Why:** `Node::getResponse(Disp)` returns `getDisp()` (committed);
  `-increment` (WP-207) sets a trial and never commits it.
- **Fix/status:** read trial minus committed as the `IncrDisp` node response
  (`nodeResponse $node $dof 4`), or use `responseSpectrumPerMode -inertia`
  (WP-216), which forms M φ s Γ Sa itself. Learned 2026-10-27 (WP-216).
