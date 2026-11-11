---
wp: WP-230
title: "\"Clamp p at a floor and keep the deviator\" is a property of the discrete update, not of the model: WHICH deviator is kept decides whether the state is admissible"
pr: "#1072"
date: 2026-11-11
---
### "Clamp p at a floor and keep the deviator" is a property of the discrete update, not of the model: WHICH deviator is kept decides whether the state is admissible (WP-230)
- **Bites:** porting a low-pressure rule of the form "if p < p_floor, set p = p_floor and keep the shear" into an integrator that REFUSES instead of returning a trial stress, or comparing two codes that both "clamp".
  - In the exact model the trajectory reaches p → 0 on the yield surface, with a deviator that vanishes with p (q ≤ M·p). There is no deviator to keep at the floor.
  - A rule that keeps the deviator of an elastic-predictor TRIAL stress (computed from a higher pressure) and resets p lands outside the bounding surface by an amount set by the step size.
  - SAS-ME has no trial stress at a refusal. Keeping the last COMMITTED deviator gives: after a tension refusal (E1) the committed p is already below the floor, clamping RAISES p, and the clamped stress is INSIDE the bounding surface (oracle, five paths: stress ratio / bounding ratio ≤ 0.34). Only after an accuracy or cost refusal at p0 > p_floor (E2 with p_sep > p_floor) does the clamp lower p and push the stress ratio out, by p0/p_floor.
  - The yield cone is thin (m = 0.005 in the campaign set), so f/(√(2/3)·m·p) of a clamped state is 100–200 cone radii whenever α is merely stale. That number alone does not say the STRESS is inadmissible.
- **Rule:** when a clamp is asked for, state which deviator is kept and measure two things separately: f in cone radii (α stale or not) and the stress ratio over the bounding ratio (`sas_clampMaxBound`; above 1 no α makes the state admissible). Never add the stress change to α: it is a stress ratio.
