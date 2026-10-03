---
wp: LEGACY
title: "Domain::commit() discards element commit returns — a material cannot refuse at commit; LadrunoQuad propagates any nonzero update code, LadrunoBrick only the se…"
legacy_seq: 446
---
### `Domain::commit()` discards element commit returns — a material cannot refuse at commit; `LadrunoQuad` propagates any nonzero update code, `LadrunoBrick` only the sentinel
- **Bites:** you write a material that detects, at `commitState()`, that the step
  it is being asked to commit is not integrable — and you return a failure code.
  Nothing happens. The step is committed, the analysis reports it converged, and
  the run walks on. Measured instance: `LadrunoSANISAND` under `-implex` without
  `-implexControl` on the TIMs plane-strain strip — **25.9 million** commit-time
  companion cap hits, a straight-line load–settlement curve to 2 674 kPa, every
  step "converged".
- **Why:** `Domain::commit()` (`SRC/domain/domain/Domain.cpp`, the element loop)
  is a bare `elePtr->commitState();`. The return value is not captured, not
  summed, not tested. Nothing downstream of it exists to propagate: by the time
  `commitState()` runs the algorithm has already declared convergence and
  `StaticAnalysis::analyze()` is past its failure branch
  (`StaticAnalysis.cpp`, the `revertToLastCommit` + `return -3` path belongs to a
  failed `solveCurrentStep`, not to a failed commit). **A refusal is only
  actionable at the TRIAL** (`setTrialStrain`).
- **And even at the trial the elements disagree** — 26 forward, 1 sentinel-only,
  25 discard, out of 52 NDMaterial hosts. The full audited table is the entry
  **"Element refusal roster"** below; it is the only authoritative copy. The
  short of it: the same refusal cuts the step on a `LadrunoQuad` mesh and is
  invisible on an `SSPquad` one, and a material that returns "some nonzero
  value" rather than the sentinel is silently swallowed by `LadrunoBrick`
  specifically. Four shipped `opserr` texts stated this wrongly ("today
  LadrunoBrick", and `QuadUP` listed as a discarder when it is in fact a
  propagator); WP-99 corrected them and then, after review round 1 found the
  REPLACEMENT list was still a wrong closed list, made them non-exhaustive and
  pointed them here.
- **Workaround/status (2026-09-14, revised after review round 1):** WP-99 makes
  the refusal leave the element path entirely. A material calls
  `ladrunoNoteCommitRefusal()` (`SRC/material/LadrunoMaterialStatus.h`) from its
  `commitState()`; `Domain::commit()` checks that counter after its element loop
  and **returns a failure**, which `AnalysisModel::commitDomain()` turns into -2
  and every analysis class turns into `-4`. That is element-independent, so a
  DISCARD element cannot swallow it. `LadrunoSANISAND` additionally keeps a
  sticky per-instance latch as a second line of defence, for a driver that
  ignores the analysis return code. **Why the latch alone was not enough, and
  this is the measurement that decided the design:** two stacked `stdBrick` with
  the lower element starved (`-maxSubsteps 2`) under `algorithm Linear` ran **20
  further accepted steps**, `analyze() == 0` throughout, with the refusing
  element frozen as a rigid inclusion — and *more quietly than before*, because a
  latched `commitState()` returns early and so the old 10-per-process cap
  warnings stopped firing too.
  **Raw element commit codes are still NOT propagated** and must not be: ADR-33/34
  requires that a negative "best-state" code (ASDConcrete3D and friends) not fail
  a step, and at commit there is no sentinel-filtering element in the path to
  tell a declaration from a diagnostic — which is exactly why the declaration
  goes out of band instead. If you want a *recoverable* refusal, refuse at the
  trial (`-implexControl` is the SANISAND example): a commit-time refusal is
  fatal by construction, because the nodes and the sibling integration points
  have already committed by the time it happens.
  Cross-links: [[LEDGER_implementations]] "IMPL-EX commit-time refusal latch",
  [[LadrunoSANISAND_implex_guide]] §9.
