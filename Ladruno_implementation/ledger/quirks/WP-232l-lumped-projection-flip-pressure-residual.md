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
- **Second bite:** the known term pi = rhoF lambda b . (affine image part) is outside the range of
  P, and it changes as the particles move. Seen on `d06695367` (no compensation, base 1.0e-3) and
  `89718fdad` (A' once at the start, but G2P did not carry pi's change with motion: base 1.2e-3;
  a rigid skeleton gave 1.9e-11).
- **Fix (interface 7b, ruling on `89718fdad`):** the particle carries q = p - pi_c; the
  projection fits q, G2P is FLIP on q, and after every accepted step the host sets
  p <- p + pi(x_new, drive) - pi_c, pi_c <- pi. A' is the invariant's Init.
- **Third bite:** the u side's mass-cutoff renormalisation rescaled every component of the kept
  entries, the odd-folded p weights and pi's coefficients included. When the soil settles away
  from a drained B-spline face, the row beyond it falls below the u cutoff and is dropped, and the
  extra division breaks linear reproduction there. The steady field was offset by a uniform
  1.24e-3 gamma_w H (build `0fc5a3f13`, E 1e4 column).
- **Fix (ruling `bb1421559`):** the u pass rescales the u components only; the p weights and the
  affine coefficients follow fix 3 on the raw weights, a folded-away entry retained even when
  u-dropped. Unit test: `test_drained_face_p_weights_reproduce_linear_after_u_drop`.

#### A GIMP tail of ~1e-8 makes a near-singular pore-pressure row

- **Bites:** a two-phase particle that moves a little toward a node only single-phase particles
  hold. With p-activity "any S^p > 0", the node's p row carries a weight of ~1e-8; seen at
  p = -6.7e3 after two steps.
- **Fix (P4 amended, fix 3):** p-active iff the two-phase nodal mass > `pressureCutoff` x m_max,
  and S^p is renormalised over the p-active nodes (route S's rule, `3d3699a50`).

#### MUSL's re-projected velocity breaks the IMEX u-p elimination (route S)

- **Bites:** a semi-implicit u-p MPM step (explicit skeleton, implicit p, v1 eliminated) whose
  skeleton strain comes from MUSL's second projection v** while the fluid balance uses the grid
  v1. The difference is O(dt) FLIP null-space content, so part of the undrained coupling acts
  explicitly. Seen on Mandel (bspline2, build `4b74bd335`): rho = 1.00093 per step at 0.30 x the
  drained pencil, in the particle column beside the drained free face; the single-phase twin is
  stable.
- **Fix (interface 7c, owner's decision):** two-phase particles take L from v1 (USL for the
  skeleton). Light u nodes (m_a <= max(massCutoff, pressureCutoff) m_max) are dropped whole in a
  two-phase step; their correction dt M^-1 Q p grows like 1/delta on a GIMP tail (gimp Mandel).
