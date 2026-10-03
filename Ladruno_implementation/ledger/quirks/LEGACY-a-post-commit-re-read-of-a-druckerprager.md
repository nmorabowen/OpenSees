---
wp: LEGACY
title: "A post-commit re-read of a DruckerPrager tangent (printA, eleResponse 'stiffness', a recorder) is the ELASTIC branch roughly half the time — it is not the stat…"
legacy_seq: 429
---
### A post-commit re-read of a DruckerPrager tangent (`printA`, `eleResponse 'stiffness'`, a recorder) is the ELASTIC branch roughly half the time — it is not the state the step actually converged with

- **Bites:** any diagnostic or test that reads a DruckerPrager-family
  material's tangent AFTER a step has committed — `ops.printA('-ret')`,
  `ops.eleResponse(tag, 'stiffness')`, a stiffness recorder sampled between
  steps — gets the ELASTIC tangent (`mCep == mCe`, exactly symmetric even
  under non-associated flow) on a roughly-coin-flip fraction of plastic
  steps, even though the step just converged plastically. Reading only the
  LAST step of a ramp for an "is this run plastic" fingerprint (asymmetry,
  branch id, whatever) is a knife edge: a *relative* `1e-12` nudge to a
  material constant (e.g. `sig_y`) is enough to flip the read.
  `tests/test_upstream_symmetrize_fixes.py::test_zerolengthnd_unsym_tangent`
  hit exactly this — its final-step-only asymmetry read measured `0.0` on
  `ladruno` head (a false regression signal; the mirror-assembly bug it
  guards was NOT back), fixed by sampling the max over every converged step
  instead of the last one.
- **Why:** the element's `getTangentStiff()` (e.g.
  `SRC/element/zeroLength/ZeroLengthND.cpp:377`) calls
  `setTrialStrain(committed strain)` to re-form the tangent for a caller
  like `printA`, which re-runs `DruckerPrager::plastic_integrator()`
  (`SRC/material/nD/UWmaterials/DruckerPrager.cpp:446`) from a stress state
  sitting EXACTLY on the yield surface (it's the state the step just
  converged to and committed). That integrator tests `f1 <= fTOL` with
  upstream's `fTOL = 0.0`; on this re-formation `f1_trial` is pure roundoff
  (`+-1e-15`-ish), so which branch it takes — elastic (`mCep = mCe`,
  symmetric) or plastic (`mCep` possibly unsymmetric under non-associated
  flow) — is decided by floating-point noise, not by physics.
- **Not a real problem:** analysis RESULTS are unaffected. The next step's
  Newton iterations re-form the tangent at a genuine trial state away from
  the surface, so the roundoff-branch tangent is only ever used for one
  diagnostic read, never to compute an increment that matters. Do NOT "fix"
  this by changing `fTOL` off `0.0` — it's vanilla upstream, has a wide
  blast radius, and isn't the defect; the fix belongs in the diagnostic/test
  code that reads the tangent.
- **Workaround/status (2026-09-08):** any test asserting on a post-commit
  DruckerPrager tangent must sample over multiple steps (max, or count of
  nonzero reads) rather than trusting the final one. Multi-Gauss-point
  volume/area elements (BezierTet10, BezierTri6) mask this in practice —
  it's unlikely EVERY Gauss point lands on the roundoff boundary
  simultaneously — but that's luck, not a guarantee; `test_beziertet10_unsym_tangent.py`
  and `test_beziertri6_unsym_tangent.py` were hardened the same way as a
  precaution. `test_pardiso_asym_rearm.py` was left alone: its guard samples
  the SOE's assembly-time (mid-Newton-iteration) tangent across ~135
  assemblies via a console-log latch, not a single post-commit read, so it
  doesn't share this mechanism.
