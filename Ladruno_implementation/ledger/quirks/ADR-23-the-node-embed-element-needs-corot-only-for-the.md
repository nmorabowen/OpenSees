---
wp: ADR-23
title: "The node-embed element needs -corot ONLY for the D9 MATERIAL frame, never for the penalty tie (ADR 23 Phase 2b v2)"
legacy_seq: 68
---
### The node-embed element needs `-corot` ONLY for the D9 MATERIAL frame, never for the penalty tie (ADR 23 Phase 2b v2)
- **The split:** the isotropic/penalty U/UP/UR tie is already frame-objective (`D=K·I` has no
  preferred axis; gap + weights both transform with the host), so it needs **no** co-rotation —
  unlike the anisotropic `LadrunoEmbeddedRebar`, whose frozen bar `dir` goes stale. But the D9
  **material interface** reintroduces a preferred axis (the `-normal`/tangent frame carrying
  per-direction uniaxials), so a directional contact normal DOES go stale under host rotation.
  `-corot` co-rotates that frame with the host CONTINUUM rotation `θ_host = skew(∇u)|_ξ`,
  `frameCur = R(θ_host)·frame` — **reusing the UR `∂N/∂x`/`rotOper` machinery verbatim** (not the
  rebar's secant-to-point-B trick; the node element has a normal+tangents frame, not a bar axis,
  so the natural rotation source is the host continuum spin, factored into `hostContinuumRotation`).
  3D = Rodrigues exp-map of the axial vector `θ_host`; 2D = the single drilling planar rotation.
- **Mechanically distinct from UR:** UR ties the cNode's rotation DOFs to `θ_host` (a constraint
  on a DOF); `-corot` rotates the *material frame* used by the translational interface (no rotation
  DOF needed — material mode runs at ndf=ndm). They share only `θ_host`/`gradN`; `-corot` is
  material-mode-only (parse-time reject otherwise) and independent of `-rot`.
- **The dropped `∂e_d/∂u` tangent term is EXACT, not approximate, when the host DOFs are
  prescribed.** `-corot` inherits the rebar's dropped consistent-tangent term (R7/D9-5: residual
  exact, tangent inexact ⇒ frame-objective for explicit, step-halving for implicit, may slow
  Newton on stiff-normal large-per-step-rotation contact). But note `frameCur` depends on
  `θ_host`, which is a function of the HOST translations **only** (`∂frameCur/∂u_cNode = 0`). So in
  a Zone-A mechanics test where the host is prescribed (`sp`/`fix`), the cNode tangent is *exact*
  and Newton converges quadratically — the inexactness only bites when the host continuum is free
  and spinning per-step. 2026-06-07.
