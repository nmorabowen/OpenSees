---
wp: WP-232l
title: "A lumped particle-to-grid projection makes FLIP particles carry a residual that undoes the lumping (two-phase MPM pore pressure)"
date: 2026-11-18
---
### FLIP particle pressure drifts from an accurate grid field next to a drained face

- **Bites:** any MPM scalar field carried on particles with FLIP increments and projected to the
  grid with the row-sum-lumped matrix, p^ = M_L^-1 P^T V p. Projection then interpolation is not
  the identity, so FLIP's fixed point (p^ = the grid field g, a pinned node included) is
  p_part = P g + r with r != 0. The particles carry a residual that compensates the lumping. PIC
  is worse: the lumped smoothing is applied every step.
- **Seen (WP-232l, build `fd81cbdd3`):** Terzaghi on gimp, Tv 0.1. The grid field at the particle
  points is within 0.0016 p0 and the nodes within 0.0015, but the particle p misses by 0.089 p0;
  the top particle is at -0.044 p0 where the field is +0.045. The hydrostatic column has an exact
  field (1e-4) and a particle error of 0.029 gamma_w H. PIC gives -0.83 p0. Records:
  `tests/data/mpm/twophase/route_e_results.json` (blocks t2, t4, diagnostic).
- **Fix (P5a, interface 7b):** project with the consistent matrix M = P^T V P (a least-squares
  fit), so M^-1 P^T V P = I. Guard the solve against singular M; one particle per cell on a
  B-spline grid is refused.

#### An even B-spline image is not zero flux under a normal fluid drive

- **Bites:** an impervious B-spline face with gravity (or `-fluidBody`) along its normal. An even
  image imposes dp/dn = 0, while zero Darcy flux needs dp/dn = rhoF b.n. Seen: the hydrostatic
  column's base particle on bspline2, 0.055 gamma_w H.
- **Fix (P9 amended, fix 2):** an affine even image, c_g = c_m + rhoF lambda b_d (X_g - X_m).

#### A GIMP tail of ~1e-8 makes a near-singular pore-pressure row

- **Bites:** a two-phase particle that moves a little toward a node only single-phase particles
  hold. With p-activity "any S^p > 0", the node's p row carries a weight of ~1e-8; seen at
  p = -6.7e3 after two steps.
- **Fix (P4 amended, fix 3):** p-active iff the two-phase nodal mass > `pressureCutoff` x m_max,
  and S^p is renormalised over the p-active nodes (route S's rule, `3d3699a50`).
