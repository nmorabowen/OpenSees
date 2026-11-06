---
wp: LEGACY
title: "A refuse/halve/regrow path does not just cost reach — it STIFFENS the curve, by ~19 % at the settlement where it dies"
legacy_seq: 435
---
### A refuse/halve/regrow path does not just cost reach — it STIFFENS the curve, by ~19 % at the settlement where it dies
- **Bites:** two `-implex` legs of the same deck are compared at matched
  settlement and the refusing one reads higher. It is easy to mistake for a real
  difference between the flags being compared.
- **Measured (ADR-92 F10):** leg B (`tol 0.05`, growth x2, 724 refusals) against
  leg N (control off, 0 refusals), at matched `s/B`: **+1.72 % at 0.0030,
  +2.56 % at 0.0050, +5.99 % at 0.0070, +19.40 % at 0.0085** — monotone. At leg
  B's terminal settlement the whole leg set splits by whether it refused
  (B +19.4 %, L +20.2 %, I +18.5 %, F2 +12.1 %; F1 +2.4 %, M +1.5 %, J +0.7 %,
  K +0.3 %, F3 -0.7 %). A refusal-free constant-`ds` walk agrees with N to 0.4 %.
- **Why:** equilibrium is found on the EXTRAPOLATED stress, so a leg that
  refused, halved and re-grew arrives at the same settlement along a different
  strain path. The material state is the implicit companion's at every commit;
  the divergence lives in the structure's kinematics.
- **Workaround/status (2026-09-14):** never quote a `q` from a leg that refused
  without saying so, and never compare two `-implex` legs on load unless both
  ran refusal-free. Compare reach and refusal counts; compare load only against
  a refusal-free arm.
