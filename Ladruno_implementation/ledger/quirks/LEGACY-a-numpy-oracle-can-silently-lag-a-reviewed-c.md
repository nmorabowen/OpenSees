---
wp: LEGACY
title: "A numpy ORACLE can silently lag a reviewed C++ fix — a red zone_a gate may indict the reference, not the shipped code"
legacy_seq: 261
---
### A numpy ORACLE can silently lag a reviewed C++ fix — a red zone_a gate may indict the reference, not the shipped code
- **Bites:** any subsystem gated by a hand-written numpy oracle that mirrors a
  C++ kernel (`concrete3d_ref.py`↔`LadrunoConcrete3DKernel.h`, and the same
  pattern in the J2/logstrain/up/hypo `*_reference.py` pairs). A fix applied to
  the kernel during PR review does NOT propagate to the oracle, and nothing in
  CI notices: the kernel-vs-oracle test compares them on a COMMITTED FIXTURE of
  paths that may never exercise the diverging branch.
- **The case:** `test_p2i_multiaxial_apportioning_gate` asserts
  `I3_pure_compression_wt < 1e-9` (no spurious compression→tension damage) and
  measured **0.997**. It looks like a formulation bug in the tensile damage
  gate. It is not — the gate `sig_t_drive = E*et if max(w) > 1e-6*ft else 0`
  is CORRECT and opens legitimately, because the effective stress handed to it
  genuinely IS tensile: under uniaxial-STRAIN compression the hardening Newton
  overshoots to `rho<0`, the apex branch teleports to the hydrostatic-TENSION
  vertex, and a trial with max principal **−23.76** returns **[+2.94,+2.94,+2.94]**
  with `conv=True`. `f==0` holds at the apex BY CONSTRUCTION, so the oracle's
  `(converged or apex) and |f_indep|<tol` reports success for a sign-flipped,
  inadmissible state. `kp` also jumps 0.034→0.182 in one step.
- **The kernel was already right.** PR #249's adversarial review added to
  `returnMapHardening` an admissibility test (`dlam>=0 && kp>=kp_n`), refusal to
  report converged, and a fallback to the ELASTIC PREDICTOR so the caller cuts
  the step — with the comment "this deliberately diverges from the numpy
  oracle's (equally-arbitrary) apex teleport — the kernel is the safe reference
  here." The oracle received only the HONEST-f-recompute half of that fix and
  kept the literal pre-fix expression the C++ comment calls out as lying. So
  **the shipped material never had the bug**; only the reference did.
- **It shipped red and stayed red.** The gate produces the byte-identical
  0.9971183764898133 at `c349e8763` (#336), the very commit that introduced it —
  #336 landed AFTER #249, so the gate was authored against an oracle that
  already lagged the kernel. It has never passed.
- **Rule:** when a zone_a gate goes red on a subsystem that has BOTH a numpy
  oracle and a C++ kernel, diff the two implementations of the disputed branch
  BEFORE touching the assertion — the oracle is as likely to be stale as the
  code. And when a PR-review fix lands in a kernel, port it to the oracle in the
  SAME PR, because the fixture-based agreement test will not catch the drift.
  Do not "fix" a gate by weakening its assertion until you have established
  which side is wrong. Cross-ref [[ladruno-adr79-geom-hypo]].
  *2026-08-04 (P2i apex-teleport hunt).*
