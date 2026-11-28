---
wp: WP-263
title: "A Dual pushed through the RC eigen formulas loses the principal-axis rotation exactly at gxy = 0 and at repeated principals"
date: 2026-11-28
status: "open"
---
### Dual-number RC tangent misses the shear column at gxy = 0, or is NaN / basis-dependent at equal principals

- **Bites:** templating `LadrunoRCKernel.h` (or any spectral split: Concrete3D, ASDConcrete3D
  clones) by pushing the Dual through the existing eigen code.
  - **`eps1FromMembrane`:** p1 = (gxy/2, lam1 - exx)/nrm. At gxy = 0 (uniaxial tension along x,
    the commonest RC state) the vector cancels to zero, the scale-relative fallback sets
    p1 = (1, 0) as a CONSTANT, and dp1/dgamma, the tension-stiffening and crack-capture
    rotation, is lost. Near gxy = 0 a nonzero `nrm` of round-off size gives a huge, wrong
    derivative instead.
  - **`eig3x3sym` (cyclic Jacobi):** the iteration stops on the VALUE of the off-diagonals, so
    the eigenvector derivatives are not converged. At repeated eigenvalues they are undefined.
    The eigenvalue derivatives `v^T dA v` are right, but inside a cluster they depend on the
    arbitrary basis Jacobi returned.
  - **`lublinerCriterion`:** `sqrt(3 J2)` at J2 = 0 (hydrostatic stress, including the zero state
    through the compressive measure) gives 0.5/0 * 0 = NaN.
- **Fix (WP-263):** none of these sees a Dual. Each takes its values from the double function,
  so they are bit-identical, and its derivatives from a closed form:
  - `decompose`: the isotropic-function derivative of ST in the eigenbasis,
    `G_ij = (<li>+ - <lj>+)/(li - lj)`, `(H(li) + H(lj))/2` when the two are equal, `H(0) = 1/2`
    (no 1/(li - lj) singularity), with dSC = dS - dST and principal-stress derivatives averaged
    over a cluster (|gap| < 1e-9 of the stress scale);
  - `eps1FromMembrane`: dp1 = dtheta (-p1y, p1x), dtheta = (D dg - g dD)/(2(D^2 + g^2));
  - `eps1Spatial`: the derivative of the top-cluster eigenprojector;
  - `sqrtKink`: zero derivative at a non-positive argument.
  Oracle: `tests/test_ladrunoRCConcrete_dual_cpp.py`. On an equal TENSILE cluster the map has a
  genuine kink (max eigenvalue), and the central FD converges to the dual at first order only.
