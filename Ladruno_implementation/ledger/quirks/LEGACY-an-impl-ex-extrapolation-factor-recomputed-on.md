---
wp: LEGACY
title: "An IMPL-EX extrapolation factor recomputed on every setTrialStrain destroys the tangent identity — and the p_min clamp destroys it again, at the floor"
legacy_seq: 376
---
## An IMPL-EX extrapolation factor recomputed on every `setTrialStrain` destroys the tangent identity — and the `p_min` clamp destroys it again, at the floor

Two traps found while writing ADR-92 P1's `-implex` on `LadrunoSANISAND`. Both
are general to IMPL-EX in this fork, not to SANISAND.

**1. `f` must be frozen for the whole step.** The extrapolated response is

```
sigma~ = sigma_n + Ce(p_n) : ( d_eps - f * d_eps_p(n) ),   f = (dt_{n+1}/dt_n)*alpha
```

and the operator handed to the element is `Ce(p_n)`, on the claim that it *is*
`d(sigma~)/d(eps)`. That claim holds only while `f` does not depend on the trial
strain. With a pseudo-time source (`ops_Dt`, which is what `ASDConcrete3DMaterial`
samples) `f` is constant within a step anyway and the trap never fires — which is
exactly why it is easy to write a strain-based or increment-based time source
later and not notice. Under `-implexDt strain`, `dt_{n+1} = ||d_eps||`, so an `f`
recomputed on every pass carries a `d f/d eps` term and the delivered tangent is
wrong by it — silently, on a run whose global iteration count *looks* like the
one-iteration IMPL-EX advertises. `LadrunoSANISAND` therefore computes `f` ONCE,
at the first trial call after a commit or a `revertToLastCommit`, and freezes it
for the rest of the step. Anything else that grows an IMPL-EX path must do the
same, and its tangent-identity test must be run on a NON-pseudo time source or it
proves nothing.

**2. A floor clamp on the extrapolated stress is a projection, so it breaks the
identity where it fires.** P0 measured `sigma~` crossing `p_min` into tension
(`-1.37 kPa` on a committed `+0.0101 kPa` state), and the repair is the code's own
device: `sigma~ = dev(sigma~) + p_min*I1`. But that is a projection, so wherever
it acts `d(sigma~)/d(eps)` is the deviatoric projection of `Ce`, not `Ce`. The
material keeps returning `Ce` (symmetric, positive definite — the reason IMPL-EX
was asked for), which makes the tangent knowingly inexact at exactly the Gauss
points sitting on the floor, i.e. the free-surface ring. **A tangent-identity gate
must therefore be run on a path where the clamp is idle, and a test that reports
"identity to machine precision" on a clamped path is measuring nothing.** The
deeper fix — bounding `f` so the whole `sigma~` stays admissible without a
projection — keeps the identity exact everywhere and is ADR-92 P1 section 4's
open alternative.

**3. A P0 numpy oracle written before the clamp existed is not a parity reference
on a clamped path.** `adr92_p0_oracle/sanisand_implex_oracle.py`'s
`Implex.extrapolate` has no clamp (the clamp is a P1 item the oracle's own G3 gate
*discovered*). So the 1e-8 parity gate is meaningful only where the C++ clamp
stays idle — on the G3 corner path the two MUST differ, and a parity run there
that passes is evidence the C++ clamp is not firing, not evidence of parity.
