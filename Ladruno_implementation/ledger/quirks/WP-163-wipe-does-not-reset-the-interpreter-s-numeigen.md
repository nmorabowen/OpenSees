---
wp: WP-163
title: "wipe does not reset the interpreter's numEigen — gate any modal path on Domain::getNumEigenvalues(), never on OPS_GetNumEigen() (WP-163)"
date: 2026-10-03
---
### `wipe` does not reset the interpreter's `numEigen` — gate any modal path on `Domain::getNumEigenvalues()`, never on `OPS_GetNumEigen()` (WP-163)
- **Bites:** `eigen 3; wipe; <new model>; recorder ladruno f -N modesOfVibration; analyze 1` killed the process:
  the recorder gated its modal write on `*OPS_GetNumEigen()` (still 3 — the Tcl global `numEigen` and
  `OpenSeesCommands::numEigen` survive `wipe`/`wipeAnalysis`), then called `Domain::getEigenvalues()` on a domain
  whose spectrum `clearAll()` had deleted → `exit(-1)` ("Eigenvalues were never set"). `Node::getEigenvectors()`
  on a node with no eigenvectors is the same trap (`exit(0)` — exit code 0, so a harness sees "success").
- **Why:** the interpreter count and the domain spectrum are separate state with separate lifetimes; only the
  domain's is cleared by `wipe`.
- **Workaround/status:** fixed in the Ladruno recorder (WP-163 R1): gate on `Domain::getNumEigenvalues()` and
  probe `Node::getNumEigenvectors()` first (both ADR46, non-exiting). Any other fork code that reads
  `OPS_GetNumEigen()` to decide whether to touch eigen data has the same bug. Resetting `numEigen` in `wipe`
  would need a vanilla edit in both interpreters (not done).
