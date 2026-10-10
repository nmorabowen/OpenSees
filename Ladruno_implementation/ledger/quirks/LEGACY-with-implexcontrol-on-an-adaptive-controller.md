---
wp: LEGACY
title: "With -implexControl ON, an adaptive controller that DOUBLES the step is a refusal generator — but the growth rule is innocent with the control off"
legacy_seq: 434
---
### With `-implexControl` ON, an adaptive controller that DOUBLES the step is a refusal generator — but the growth rule is innocent with the control off
- **Bites:** with the control on, the refusal count tracks the controller's
  growth factor and nothing else: x2 -> 724 refusals and FLOOR at `s/B = 0.0085`;
  x1.25 -> 253 and 0.0113; x1.0 -> **6** and 0.0265, zero subdivisions in 1 984
  converged steps. It reads as a material effect and is not one.
- **Why:** `implexError = ||sigma~ - sigma_impl|| / (||sigma_impl|| + P_atm*||eps||)`
  (`LadrunoSANISAND.cpp:2629`-`:2646`) is first order in the step — measured on a
  controlled refinement at a fixed committed state (max error 1.126e-2 /
  5.778e-3 / 3.039e-3 / 1.672e-3 for `ds` = 8e-5 / 4e-5 / 2e-5 / 1e-5, ratios
  1.95 / 1.90 / 1.82, onto a `dt`-independent floor of ~4e-4). `-implexControl`
  bounds it ABSOLUTELY, so a halve-on-failure / double-after-N controller can
  only find the bound by crossing it — with the clock ratio `f = dt_{n+1}/dt_n`
  sitting at 2 on exactly the crossing step, so the extrapolated plastic
  increment is doubled on top of the doubled strain increment.
- **But it is a CO-FACTOR, not a cause:** with the control OFF, growth x2 and
  growth x1.0 both reach the target and agree to **0.383 % mean / 0.779 % max**,
  and x2 is **6x faster** (58 s vs 364 s). Do not "fix" a growth rule that is
  only pathological under a flag you may not need.
- **Workaround/status (2026-09-14):** if you keep `-implexControl`, pin the
  growth factor at 1.0 (x1.25 is not enough) or drive the controller off
  `avgImplexError` so the bound is approached from below. Also read the
  `q`-column warning below.
