---
wp: ADR-44
title: "ADR-44 LadrunoModalResponse — modal-transient (P1a) + FRF/SSD (P2) gotchas"
legacy_seq: 168
---
### ADR-44 LadrunoModalResponse — modal-transient (P1a) + FRF/SSD (P2) gotchas

- **Element-level stiffness-proportional (`betaK`) Rayleigh damping ≠ assembled
  `a1·K` (Truss).** A direct `Newmark` run with `rayleigh 0 0 0 a1` (or `betaKcurr`
  / `betaKinit` — all three identical) on a linear Truss chain differs from the
  EXACT classical-modal solution of `M ü + a1 K u̇ + K u = −MR ü_g` by **several
  percent, dt-INVARIANT** (measured 4.4% on a 2-DOF chain; does NOT shrink as
  dt→0, so it is NOT Newmark truncation). The exact modal answer was pinned three
  ways (numpy modal-superposition 1e-18, numpy full-matrix Newmark of (M,a1K,K)
  ~9e-5, and `modalResponseHistory` itself). Mass-proportional (`alphaM`) Rayleigh
  DOES match OpenSees-Newmark to truncation (~1e-4). **Consequence for validation:**
  never oracle an exact modal-damping feature against an OpenSees direct-`Newmark`
  run that uses `betaK` — use `alphaM`-only, `modalDamping`, or an explicit
  full-matrix Newmark reference. (Root cause: element `getRayleighDampingForces`
  builds the damping force from a per-element stiffness that is not the same as the
  globally-assembled `K` for the truss basic-system mapping — not chased further.)

- **`Path` timeSeries `getFactor(t)` returns 0 at exactly the record end.** For a
  transient whose last station time `t = nsteps·dt` coincides with the final
  abscissa, floating-point round-off of `nsteps·dt` lands just past the last
  sample ⇒ `getFactor` returns 0, not the last value. Symptom: the FINAL committed
  station is slightly off (all earlier stations exact). Fix in models/tests: pad
  the record with ≥1 trailing sample so the analysis window is strictly interior.
  (This is inherent to any getFactor-sampling integrator, incl. UniformExcitation.)

- **`eigen` on a tiny model needs `-fullGenLapack`.** The default ARPACK path
  requires `NEV < N` (NCV must be > NEV and ≤ N); a 1- or 2-DOF model with
  `eigen 1`/`eigen 2` fails `_saupd info = -3`. Use `eigen -fullGenLapack N`
  (dense) for small verification models.

- **Recorder text precision.** Node recorders default to ~6 significant figures;
  a modal-transient vs analytic assert at 1e-8 is dominated by that rounding.
  Add `-precision 15`/`16` to the recorder (already noted for ADR-46; re-bitten here).

- **`responseSpectrumAnalysis -scale` is a dead knob (stock).** The Petracca
  `ResponseSpectrumAnalysis::solveMode` computes `u = V*Vscale*MPF*Sa/λ` and never
  multiplies by `m_scale`, so the `-scale` factor the command parses is silently
  ignored. Not our bug (upstream); the ADR-44 P1b `-combine` path matches this
  behavior (no scale) for consistency rather than silently diverging. Flagged
  during the P1b wiring.

- **FRF sign convention (P2) is `e^{+iΩt}` — pin it or the phase silently flips.**
  `frequencyResponse`/`steadyStateDynamics` use `H_a(Ω)=1/(ω_a²−Ω²+iΩd_a)` with a
  `+iΩd_a` imaginary part, i.e. the `e^{+iΩt}` time convention → the response LAGS
  90° at resonance (`u=+i/(ωd)` for the mass-normalized SDOF, `angle=+90°` because
  the base-accel modal load carries the extra `−Γ`). The opposite (`e^{−iΩt}`,
  `−iΩd_a`) magnitude is IDENTICAL and only the phase sign differs, so a flipped
  convention passes every |·|/RMS test and is invisible until someone reads phase —
  the classic frequency-domain bug. Pinned two ways: the resonance-phase assert in
  `test_sdof_frf_closed_form`, and the end-to-end match to the direct complex solve
  `(K−Ω²M+iΩC)⁻¹(−MR)` (same convention on both sides) in `test_mdof_frf_vs_direct`
  and `modal_response_p2_spike/frf_oracle.py`. Magnitude gates alone would NOT catch
  a sign flip.

- **P2 frequencies are Hz, not rad/s.** `-freq fmin fmax nf` bounds and the output
  `f` column are Hz; internally `Ω=2πf`. So the FRF magnitude peaks at `f_a=ω_a/2π`,
  not at `ω_a`. (A `-radps` option was deliberately NOT added in v1 to keep the knob
  count down; documented in the guide.) Reusing the P1a `LadrunoModalDamping.coeff`
  gives `d_a=2ξ_aω_a` with `ω_a` in rad/s, consistent with the internal `Ω`.

- **P2 reuses P1a's normalization by construction — do NOT re-derive it.** The FRF
  recovery weight is `ψ_a(node,dof)·(−Γ_a)` with `ψ_a=eigenvector·Vscale` and `Γ_a`
  the `modalProperties` participation factor — byte-for-byte the P1a per-mode
  recovery ingredients, so the frequency-domain result is the analytic steady state
  of the exact same modal ODE P1a integrates in time. This is WHY the numpy oracle
  can use clean mass-normalized modes (`m_a=1`, `ψ=φ`, `Γ=φᵀMR`) while the C++ uses
  OpenSees' Vscale/Γ: P1a already proved the two normalizations agree physically
  (it matched direct Newmark), so P2 inherits that instead of re-validating it.

- **An UNDAMPED mode sampled exactly at `Ω=ω_a` gives an infinite FRF.** With
  `d_a=0` the denominator `ω_a²−Ω²+iΩd_a` is real and hits 0 at resonance → `inf`.
  Harmless with any real damping (imag part ≠ 0). **CORRECTION (ADR-44 P2 review):**
  the `-biased` grid clusters points in a ±5% window around each in-band `ω_a/2π`
  AND its `k=0` sample lands **EXACTLY ON** `ω_a/2π` (`f = fa + halfw·(k/NCLUST)`,
  k=0 ⇒ f=fa — `LadrunoModalResponse.cpp:632-634`). So a `-biased` sweep with
  **zero** damping DOES manufacture this singularity (inf/NaN row at every in-band
  mode) — the parser now emits a one-line WARNING when `-biased` is combined with an
  exactly-zero `d_a`. A `-lin` grid landing exactly on an undamped resonance would
  do the same. The same singularity appears
  at `Ω=0` (the `f=0` sample of a `fmin=0` sweep) if a RIGID mode is retained
  (`ω_a=0` → denom `0+0i`) — an unrestrained base-excited structure genuinely has an
  undefined DC steady displacement. Both are honest physical singularities (NaN/inf
  row), documented not guarded; fully-supported base-excitation models never carry a
  rigid mode.
