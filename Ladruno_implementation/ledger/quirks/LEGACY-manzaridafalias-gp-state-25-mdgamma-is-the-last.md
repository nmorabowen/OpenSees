---
wp: LEGACY
title: "ManzariDafalias gp_state[25] (mDGamma) is the LAST SUBSTEP's plastic multiplier, not the step total, under every substepped scheme"
legacy_seq: 364
---
## `ManzariDafalias` `gp_state[25]` (`mDGamma`) is the LAST SUBSTEP's plastic multiplier, not the step total, under every substepped scheme

**Found 2026-09-05, ADR-92 scoping; measured by the P0 oracle (`_adr92_p0_oracle_results` §6).**

`BackwardEuler_CPPM` zeroes `NextDGamma` and solves for it (`ManzariDafalias.cpp:2220`, `:2274`), so under
`IntScheme 2` it is the increment's plastic multiplier. Under `ModifiedEuler` (scheme 1, the
deck default), `RungeKutta4/45` and the `MaxStrainInc`/`MaxEnergyInc` family, every substep
**overwrites** it (`ModifiedEuler` `:1498`, `:1570`; `RungeKutta4` `:1744-1807`; `RungeKutta45`
`:1973-2016`; `ForwardEuler` `:1342`, which the `MaxStrainInc`/`MaxEnergyInc` FE variants call)
and nothing sums it — the recorder sees whatever
the last substep computed, which depends on how many substeps the controller chose.

- **Bites:** anything that reads `state[25]` as "how much plastic flow this step" — a recorder,
  a dilatancy post-processor, and any IMPL-EX that extrapolates `dGamma`. The TIMs request for
  ADR-92 proposed exactly that; the oracle measured the two extrapolation forms **identical to
  2e-10 kPa under scheme 2 and 37 kPa apart on a 277 kPa stress under scheme 1**.
- **Workaround:** integrate the plastic strain yourself from the committed `eps - eps_e`
  (`getState` entries 0-5 are `eps_e`; `strain` is `eps`) — that IS a step total. ADR-92 D1.
