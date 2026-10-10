---
wp: LEGACY
title: "LadrunoSANISAND at the fork-default -Presidual 0 clamps a free-surface Gauss point on a DILATING deck, and says so -- but only in opserr"
legacy_seq: 361
---
## `LadrunoSANISAND` at the fork-default `-Presidual 0` clamps a free-surface Gauss point on a DILATING deck, and says so -- but only in `opserr`

**Found 2026-09-05, ADR-90 WP-A2.**

With `-Presidual 0.0` (the fork default -- a cohesionless sand has no cohesion) the low-p floor is
`-Pmin` alone, `1e-3 * P_atm = 0.101` kPa. On a strip footing on a DENSE (`e_init = 0.60`,
`psi = -0.22`) sand, the dilating soil under the footing edge drives a shallow Gauss point onto
that floor and the material logs

    WARNING ManzariDafalias::ModifiedEuler() - material tag 1: mean stress p = 0.100973 is below
    the floor m_Pmin + m_Presidual = 0.101; CLAMPING the stress to p = 0.101 (deviator
    preserved). The result at this integration point is set by the clamp, not by the model.

- **Measured onset:** `s/B ~ 0.0153` on the COARSEST mesh (h0 = 1.0) and the DENSE density only.
  The loose (`e_init = 0.6944`) legs and both finer meshes fired **zero** clamps over the same and
  deeper settlements -- so it is a dense-dilatant / coarse-element / large-settlement effect, not a
  property of the deck.
- **It is silent in Python.** The message goes to `opserr`, so a driver that does not capture it
  (see the `ops.logFile` quirk) will keep integrating a model whose answer at that point is the
  clamp's. Count it, do not hope to read it.
- **Any deep push on dense sand owes a decision here** -- a declared non-zero `-Presidual`, a small
  surcharge to keep the free surface confined, or an explicitly accepted and disclosed clamp.
  The gravity state is no warning at all: the shallowest Gauss point sat at 1.56-6.25 kPa, 15-60x
  the floor, on every mesh, before the push ever started.

> **DOCUMENTED, NOT CHANGED, in WP-86b (ADR-86b T4, PR pending). The default STAYS `-Presidual 0`.**
> The three options — a declared non-zero `-Presidual`, a small surcharge, an explicitly accepted
> clamp — with what each buys and costs, and a fork-side (non-binding) recommendation, are written
> up in `86_ladruno_sanisand_handoff.md` §5b and, in consumer language, in
> `86_ladruno_sanisand_apegmsh_emitter_guide.md` §7.
> **It is a modelling decision on a CALIBRATED soil, not a numerical tidy-up, and it is not one
> parameter:** per the ADR-86 PR-3 tripwire memo, `p_residual` ALSO bounds the `D_factor` dilatancy
> sigmoid from below — vanilla's 1.01 kPa held `D_factor >= 0.4278` while `p_r = 0` drops that floor
> to 4.83e-4, a factor of **886** — so restoring a non-zero `p_r` re-engages ADR-86 **D5a**, which
> is still open. Whichever option is taken, COUNT the clamp events and report the number.
> **It gets closer, not further, after WP-86b:** the substep cap and the `TanType` default exist to
> let a leg reach deeper settlements, and this clamp is the next thing waiting there.
