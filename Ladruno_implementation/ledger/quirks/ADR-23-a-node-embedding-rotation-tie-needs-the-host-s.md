---
wp: ADR-23
title: "A node-embedding ROTATION tie needs the host's ∂N/∂x, not its weights N_i (ADR 23 Phase 2 / UR)"
legacy_seq: 67
---
### A node-embedding ROTATION tie needs the host's `∂N/∂x`, not its weights `N_i` (ADR 23 Phase 2 / UR)
- **Why it bites:** the translational (U) and pressure (UP) ties only need the host
  shape-function WEIGHTS `N_i(ξ)` (`getInterpolationWeights`, ADR 20). The rotation
  (UR) tie ties the constrained node's rotations to the host CONTINUUM rotation
  `θ = ½ curl(u) = skew(∇u)`, which is built from the host displacement GRADIENT — so
  it needs `∂N_i/∂x` (cartesian shape derivatives), a DIFFERENT host query that
  weights cannot supply. Hence the new vanilla `Element::getInterpolationGradients(ξ,dNdx)`
  (default −1; overridden on `LadrunoBrick` via `shp3d`, `BezierTet10` via
  `computeJacobian`). The translational rows of the UR `B`-matrix still use `N_i`; only
  the rotation rows use `∂N/∂x`.
- **Volume host ⇒ PURE `skew(∇u)`, not ASD's mixed convention.** `ASDEmbeddedNodeElement`
  embeds into a planar tri/tet *surface*, so it builds a 2D local frame and uses the
  surface SLOPE (factor 1) for the two bending rotations + `½ curl` (factor ½) only for
  the drilling — a mixed convention forced by the missing out-of-plane derivative.
  `LadrunoEmbeddedNode` embeds into a 3D VOLUME host (hex/tet) where all 9 gradient
  components are available, so it uses the dimensionally-clean **pure continuum rotation**
  `θ = ½ Σ_i (∇N_i × u_i)` (½ on all three, no local frame, frame-objective) — de Souza
  Neto §3. The host operator block is `½·skew(∇N_i)`; the gradient virtual returns global
  cartesian `∂N/∂x` directly, so NO `R`-rotation of the block is needed.
- **UR is mesh-limited (UR-4):** on a CST (3-node tri) / TET4 (4-node tet) host `∂N/∂x`
  is element-CONSTANT ⇒ the UR constraint collapses to a single element-wide RIGID-SPIN
  tie (no intra-element rotational gradient). Moment-critical embeds (anchors, headed
  studs) need a higher-order host (`BezierTet10`) where `∂N/∂x` varies with ξ. Document,
  don't silently sell as exact. 2026-06-04.
