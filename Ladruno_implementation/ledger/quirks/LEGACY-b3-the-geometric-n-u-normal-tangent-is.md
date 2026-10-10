---
wp: LEGACY
title: "B3: the geometric ∂n/∂u normal tangent is SYMMETRIC but INDEFINITE — a convergence-basin tradeoff"
legacy_seq: 132
---
### B3: the geometric ∂n/∂u normal tangent is SYMMETRIC but INDEFINITE — a convergence-basin tradeoff
- **Bites:** `contact … -geomtan` on large, soft-master contact patches. `K_geom = kn·gN·H` with the gap
  `gN < 0` in contact, so the geometric block SUBTRACTS from the main `kn·BᵀB` (PSD) ⇒ the contact tangent
  can be indefinite (still symmetric — it's the Hessian of the scalar gap — so ProfileSPD factors it, but
  Newton can leave the convergence basin far from the solution). **Observed:** a 749-slave-node fixed-sphere
  Hertz patch on a soft deformable master DIVERGES with `-geomtan` but CONVERGES without it; a single
  warped-quad slide (few DOFs) converges FASTER with it (quadratic, 4 vs 7 iters). So the geometric tangent
  improves LOCAL (near-solution) convergence but can shrink the GLOBAL basin on big soft patches.
- **Why it's GATED off-default:** exactly this tradeoff. `-geomtan` is an opt-in refinement (like
  `-consistanttan`), default OFF ⇒ the robust `kn·BᵀB` + byte-identity. Turn it on for curved /
  large-sliding interfaces where quadratic Newton is wanted and the patch is well-conditioned; pair with a
  line search / load stepping for robustness on large soft patches. Found by the B3 Hertz study.
