---
wp: WP-199
title: "A lumped PML must lump its inertial memory with the mass rule: HRZ mass + consistent C-PML damping grows without bound"
date: 2026-11-16
---
### A lumped PML must lump its inertial memory with the mass rule: HRZ mass + consistent C-PML damping grows without bound

- **Bites:** `LadrunoCPML ... -lumped` (the default for explicit and for the
  apeGmsh emitter) under a long run. Absorption is excellent for ~1 s, then the
  energy grows exponentially (x10 per ~0.05-0.1 s) in a stretched layer, as a
  mode polarised along the layer (u_x in a base layer) that alternates sign node
  to node. The growth rate does not depend on dt (semi-discrete, not Newmark),
  and kappaMax, alphaMax, R0, a structured layer mesh or a finer h do not cure
  it. Found by apeGmsh_Examples example 04 (Lamb + tunnel, 2026-11-16).
- **Why:** the stretched inertial operator `J_s(w)(iw)^2 u` splits into an
  acceleration piece (`kappa0 kappa1 kappa2 * a`), a velocity memory piece
  (`g_m[v]`) and corner displacement-inertial pieces. With `-lumped` the first
  used the HRZ diagonal while the memory pieces kept the consistent
  `rho N_a N_b` coupling. The layer was then not a complex-stretched copy of
  any one lattice. The Q4/T3 lumped lattice has grid modes that are backward
  in the layer normal (the P sheet near the tangential Nyquist: at kx = pi/h,
  omega falls 462 -> 308 rad/s as ky grows, for the ex04 soil), and a PML
  amplifies backward modes instead of damping them. Consistent mass (one rule
  for every piece) was stable all along.
- **Fix (`LadrunoCPMLKernel.h`):** under `-lumped` every node carries its own
  copy of each GP's inertial cascades, fed with that node's own v / u and
  weighted by the per-GP HRZ factor `N_a^2 / sum_b N_b^2`. For every supported
  quadrature (symmetric GP sets) `sum_b N_b^2` is the same at all GPs, so the
  weights reproduce the element-level HRZ mass exactly: the whole inertial
  operator is now lumped one way. The state size grows by
  nGP x nNodes x (inertial cascades); pack/unpack carries it.
- **Cost:** short-run reflection in the 2-D lumped anchors rises from the
  (unstable) mixed scheme's values to the consistent-mass level: Q4 0.9-2.1 %,
  Lamb 1.07 % (still 36x below `ASDAbsorbingBoundary2D`), unstructured T3 at the
  45-deg corner 3.6 % (gate 3 % -> 4 %; consistent layer mass gives 1.2 %).
- **Guard:** `tests/test_ladruno_cpml_late_time_stability.py`, a 2000-step
  free-surface half-space over a base + side layer (lumped, lumped kappaMax=1,
  consistent). Before: interior KE back to the source peak by 1.1 s, u_max
  6e13 m at 3 s. After: KE ~1e-10 of the peak at 3 s. The old checks missed it
  because the Lamb/Green anchors run <= 190 steps and the 2500-step Lyapunov
  case is a 2x2 toy without a free surface. Rule: any new PML/absorber option
  needs a >= 2000-step half-space run, not only a short reflection anchor.
