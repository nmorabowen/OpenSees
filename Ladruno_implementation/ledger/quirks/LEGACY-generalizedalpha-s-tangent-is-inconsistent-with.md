---
wp: LEGACY
title: "GeneralizedAlpha's tangent is INCONSISTENT with its own residual for αM≠1 (matters for DDM/Newton)"
legacy_seq: 134
---
### GeneralizedAlpha's tangent is INCONSISTENT with its own residual for αM≠1 (matters for DDM/Newton)
- **Bites:** anything that relies on `GeneralizedAlpha::formEleTangent` being the true ∂R/∂U — DDM
  sensitivity (ADR-52 W3-I2 `LadrunoGeneralizedAlpha`) most sharply, and Newton convergence rate generally.
- **The quirk:** `GeneralizedAlpha::update()` calls `setResponse(*Ualpha, *Ualphadot, *Udotdot)` — it sets
  the model acceleration to the **full-step `Udotdot`**, not the αM-intermediate `Ualphadotdot`. So the
  PRIMAL dynamic residual the elements assemble (`getResistingForceIncInertia`) is
  `R = F(t+αF·dt) − P(Ualpha) − C·Ualphadot − M·Udotdot`, whose consistent Jacobian has **M-coef `c3`**.
  But `formEleTangent` emits `αF·K + αF·c2·C + **αM·c3·M**` — M-coef `αM·c3`. For αM≠1 the tangent ≠
  ∂R/∂U. (Strictly this is also a non-textbook generalized-α: inertia "should" act at `Ualphadotdot`. HHT
  is FINE — it has no αM, M acts at the full step, tangent `c3·M` is consistent.)
- **Why it usually goes unnoticed:** for a LINEAR system Newton still converges to the residual-defined
  fixed point regardless of the (inexact) tangent — just not quadratically. Only when you DIFFERENTIATE the
  converged solution (DDM) does the tangent inconsistency bite: a tangent-reuse DDM gives a biased gradient
  (the W3-I2 PR2 FD-vs-DDM oracle failed ~2e-3 at αM=0.9 before the fix).
- **Fix used by `LadrunoGeneralizedAlpha` (#415):** build the sensitivity RESIDUAL primal-consistent (M at
  `Udotdot`, no αM — like Newmark) AND **re-form the sensitivity-solve tangent with `c3·M`** (a
  `sensTangentFlag` branch in `formEleTangent`/`formNodTangent` + a `formTangent(CURRENT_TANGENT)` in
  `computeSensitivities`) instead of reusing the inconsistent factored primal tangent. Primal path
  untouched ⇒ byte-identical. The K/C terms (αF) need no such fix — K acts at `Ualpha`, C at `Ualphadot`,
  consistent with the tangent. Found by the W3-I2 PR2 adversarial review + the Zone-A FD oracle.
