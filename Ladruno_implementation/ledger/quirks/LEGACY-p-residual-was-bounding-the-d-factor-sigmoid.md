---
wp: LEGACY
title: "p_residual was BOUNDING the D_factor sigmoid — the two low-p devices are coupled"
legacy_seq: 342
---
## `p_residual` was BOUNDING the `D_factor` sigmoid — the two low-`p` devices are coupled

**Measured 2026-08-27, ADR-86 PR-3.** `GetStateDependent` computes one `p` that
**includes `m_Presidual`** (`:4914`) and hands it to both `GetPSI` and the `D_factor` sigmoid.
The sigmoid's half-suppression point is `7.6349/7.2713 = 1.0500 kPa`; vanilla's `m_Presidual` is
**1.01 kPa**. So in vanilla the sigmoid can never suppress dilatancy below `D_factor = 0.4278`,
however low the true confinement. With `p_residual = 0` — which is `LadrunoSANISAND`'s default —
the floor drops to `4.830e-4`, a factor of **886**.

Measured on the confine-first deck, same prescribed strain path both legs: min `D_factor` along
the path is **0.7227** at `p_r = 1.01` and **0.0016821** at `p_r = 0`.

**Consequence for anyone running `p_r = 0`:** you have not only removed an apparent cohesion,
you have un-masked a dilatancy suppressor that vanilla kept bounded. See
[[86_ladruno_sanisand_pr3_tripwire_memo]] §1; the modelling call (ADR 86 D5a) is open and, per
D8, is Prof. Gorini's.
