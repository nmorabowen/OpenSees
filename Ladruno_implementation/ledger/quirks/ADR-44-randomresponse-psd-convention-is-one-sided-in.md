---
wp: ADR-44
title: "randomResponse PSD convention is ONE-SIDED in Hz — the factor bugs have unmistakable signatures (ADR 44 P3)"
legacy_seq: 175
---
## `randomResponse` PSD convention is ONE-SIDED in Hz — the factor bugs have unmistakable signatures (ADR 44 P3)

The `-inputPSD` time series is sampled at **f in Hz** and read as the **one-sided**
PSD `G(f)` of the base acceleration (`σ_üg² = ∫₀^∞ G df` — the wind / floor-vibration
/ equipment-spec convention). Against the random-vibration-textbook **two-sided
rad/s** PSD `S(Ω)`: `G(f) = 4π·S(Ω=2πf)`, and the white-noise SDOF anchor becomes
`σ_x² = G0/(8ξω³)` (NOT the textbook `πS0/(2ξω³)`). If a future edit scrambles the
convention, the Monte-Carlo gate reads it immediately: a one-sided/two-sided mixup
shows as a ~41 % (√2) RMS error, an Hz/rad mixup as ~150 % (√2π) — both pinned in
`modal_response_p3_spike/psd_rms_oracle.py` (0.6 % agreement when correct). Related
trap in the same spike: a synthetic realization `Σ√(2G·df)·cos(2πf_k t+φ_k)` has
EXACT variance only over a full period `T = 1/df` — validate over whole periods or
the input-variance check itself wobbles.
