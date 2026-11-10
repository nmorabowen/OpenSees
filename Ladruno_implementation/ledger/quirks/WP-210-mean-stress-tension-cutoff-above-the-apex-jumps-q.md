---
wp: WP-210
title: "A mean-stress tension cutoff ABOVE the cone apex that collapses the stress to p_cut*I makes sigma(eps) jump at p_tr = p_cut"
date: 2026-11-09
---
### Footing BVPs stall with "material refused" on a mesh that converged when coarser

- **Bites:** a (p, q) cone model with cohesion has its apex at p = -c cot(phi). A tension
  cutoff p >= p_cut = -sigma_t with the default sigma_t = 0 sits ABOVE that apex. LadrunoHS
  sent every trial with p_tr <= p_cut to sigma = p_cut*I (q = 0), which is only the exact
  return when p_cut IS the apex. Just above p_tr = p_cut the ordinary shear return keeps
  q ~ M_c a_apex (about 20 kPa for c = 10, phi = 25); just below, q = 0. Surface Gauss points
  beside a footing edge sit near p = 0, so Newton chatters across the jump. Finer meshes put
  points closer to the edge singularity and stall earlier.
- **Misread as:** "return-map refusal in the shear+cap corner". The refusals were the
  |deps| > 5 runaway-increment guard tripping on Newton iterates that had ALREADY diverged;
  the oracle solves every refused state. Read a refusal dump (st_n, eps_trial) before
  blaming the solver it names.
- **Why the gates missed it:** the tension FD row tested a trial deep in tension, where the
  collapsed map is constant, and the round-3 "apex exactness" assert pinned the wrong
  behaviour. Nothing scanned a strain path ACROSS p_tr = p_cut.
- **Fix (WP-210 round 4):** `solve_tension_corner` holds p = p_cut and keeps q continuous
  (q_tr if admissible, else the shear root at p_cut, the failure line or the non-hardening
  cap), with exact sensitivities; it steps aside when dilatancy alone lifts p above p_cut.
  The apex collapse stays only for sigma_t >= c cot(phi). Gates: kernel check cases
  `tension_corner` (FD tangent, rc 18) and `tension_continuity` (no jump along a path
  through p_cut, rc 19); oracle T2-C row and continuity scan. Mutation (old collapse) fails
  rc 18 with FD error 157.
- **Generalises to:** any multi-surface return whose surfaces intersect away from the
  model's natural corner. Test continuity of sigma(eps) along a strain path that crosses
  every surface boundary, not only FD tangents inside each regime.
