# norsand_calib/harness: the WP-144 P3 calibration harness (Zone B)

Fits LadrunoNORSAND's (χ, h, N, N̄, ρ, ρ̄) to element-test curves through **O2**, the algorithmic oracle the C++
kernel matches to ~1e-10 ([`../../norsand_oracle/o2_algo`](../../norsand_oracle/o2_algo/README.md)). Plan:
[`144_ladruno_norsand_plan.md`](../../../../Ladruno_implementation/144_ladruno_norsand_plan.md) §3 P3 and §5.2 K3
(re-aimed 2026-10-02: Toyoura, Tatsuoka et al. 1986 drained plane strain, DM04 Table 1 CSL). Equations:
[`144a_norsand_equation_sheet.md`](../../../../Ladruno_implementation/144a_norsand_equation_sheet.md).
numpy/scipy, so it runs on Esmeralda (`~/wp144/venv`), from `norsand_calib/`.

| File | What it does |
|---|---|
| `sand.py` | What the fit does not touch: the fork CSL (S.22), M, p_a, the DM04 elastic targets G(p, e), ν. `toyoura_dm04()` reads the data pack's DM04 Table 1 transcription (`../data/dm04/toyoura_table1.csv`, p_a 100 kPa assumed). `TOYOURA_PLACEHOLDER` (as recalled, p_a 101) is the smoke's set only |
| `energy.py` | The energy plug. `get("BA06")` maps the targets to (p₀, κ̂, μ₀, α₀ = 0) at a representative p (sheet §15; `ElasticPolicy` 'per_test' or 'global'). `get("HAR")` is the §2.3 slot: `available()` is False until an oracle has an `energy` field; binding it is `HAR.O2_FIELDS` (and `HAR.constants` if the sheet defines HAR's constants differently). No other file names an energy |
| `model.py` | Fitted vector + `Setup` (sand, energy, policy, WW, smooth cap c₁ 0.05 / c₂ 0.15, π_i0 rule) → oracle `Params` and the isotropic initial state. π_i0 rules: 'ramp_end' (default: the surface through (p_init, c₂M), see Findings), 'on_surface' (the apex through the isotropic stress, π_i0 = p (1−N)^((1−N)/N)), 'ratio' |
| `drivers.py` | O2 element drivers PS, TC, TE (drained, mixed control on O2's step: Newton with its consistent tangent, then bracket + Brent), PSU, TCU (constant volume, TIMs rung C6); `run_o1` is the O1 counterpart |
| `data.py` | Curve CSV schema `eps_a_pct, sr, eps_v_pct` + metadata (docstring), point-test CSV for φ_peak(e), ε_peak(e) |
| `objective.py` | Residuals in sr and ε_v up to the data peak + 2 % (post-peak weight 0.5), plus peak φ and ε_peak; σ's stated in its docstring (`Weights`) |
| `fit.py` | The refusal rules as a box reparametrisation; `multistart` (least_squares TRF, scrambled Sobol starts, processes); `identifiability` (Jacobian in ln θ: singular values, condition, correlation); `profile` |
| `validate_drivers.py` | Drivers vs O1, first order (G1 criteria); K1.7 undrained endpoint; the start census (apex / ramp_end) |
| `smoke_recovery.py` | The harness gate: recover known parameters from an O2-generated synthetic set |
| `fit_tatsuoka.py` | The K3 fit entry point on the data pack |
| `run_smoke.sbatch` | The smoke on a compute node |

Outputs land in `../out/` (`validate_drivers.json`, `smoke_recovery.json`, `synthetic/*.csv`).

## Run

```bash
python -m harness.validate_drivers            # ~1 min
python -m harness.smoke_recovery --starts 8 --workers 8      # or: sbatch harness/run_smoke.sbatch
```

The K3 fit is `python -m harness.fit_tatsuoka` (data pack `../data/`, DM04 sand via `sand.toyoura_dm04()`, ρ pinned at
c = 0.712, free χ, h, N, N̄, ρ̄; `--eval name=value,…` evaluates once). One evaluation of the full Tatsuoka objective
(4 curves + 22 point tests) took 85 s with the apex start (Esmeralda login node, 2026-10-02). A 5-parameter start
runs about 6 evaluations per TRF iteration, so the fit belongs on a compute node with one process per start.

## Results (Esmeralda login node, numpy 2.2.6 / scipy 1.15.3, 2026-10-02; logs and JSON in `../out/`)

**Driver validation** (`validate_drivers.log/.json`). Every driver passes. It converges to O1 at first order,
n = 25…400, at σ3′ 49 kPa, e 0.716, π_i0 = 0.8p. Orders are the least-squares log-log slopes:

| Driver | e_σ | e_π | e_ε |
|---|---|---|---|
| PS | 0.972 | 0.971 | 0.955 |
| TC | 0.988 | 0.987 | 0.940 |
| TE (to 1.5 %) | 0.996 | 0.994 | 0.980 |
| PSU | 1.007 | 1.006 | – |
| TCU | 0.999 | 0.999 | – |

K1.7 undrained endpoint at |ε_a| = 200 %:
- PSU: p within 6.8e-6 of p_cs, ζq/|p| within 8.7e-8 of M.
- TCU: 9.5e-5 and 1.0e-6.

**Recovery smoke** (`smoke_recovery.log/.json`). **PASS.** Synthetic set: two O2 PS curves, 188 residuals, ramp_end
start, BA06 'per_test', placeholder Toyoura. Truth: χ −3, h 150, N 0.3, N̄ 0.2, ρ 0.712, ρ̄ 0.75. Runs used 4 Sobol
starts and `max_nfev` 60.

| Config | Best cost | Max rel. error | Starts that recovered | Jacobian (d r/d ln θ) condition |
|---|---|---|---|---|
| `six` | 1.8e-22 | 2.6e-13 | 2 of 4 | 2.0e2 |
| `rho_pinned` | 6.4e-22 | 1.2e-12 | 1 of 4 | 1.4e2 |

- The other starts stopped on xtol in local minima: costs 7.6 and 39.8 (`six`); 10.7, 1.1e4 and 1.6e4 (`rho_pinned`).
  Use many starts.
- Weakest direction: mostly N (0.93), with h and N̄.
- Unit-residual relative SDs (`six`): N 0.21, N̄ 0.16, h 0.060, ρ̄ 0.046, ρ 0.029, χ 0.017.
- Profiles at ×0.95 / ×1.05 are upper bounds after 15 TRF iterations:
  - `rho_pinned`: χ 7.3 / 7.9, h 0.37 / 0.64.
  - `six`: χ 7.1 / 5.6, h 0.36 / 0.33.
  - So χ is sharply defined and h much less so; the other parameters compensate for it.

**Cost.** One smoke residual evaluation (2 curves) takes 2.4–2.7 s. A start takes 44–294 evaluations, 9–37 min. The
smoke took 5339 s wall with 4 niced processes. One evaluation of the full Tatsuoka objective (4 curves + 22 point
tests, ramp_end start) takes 63 s. A 16-start K3 fit is about 3–8 h on 16 cores.

## Findings that bind the fit (Esmeralda, 2026-10-02)

- **The default start is 'ramp_end', not the apex.** From the apex (`pi0_rule='on_surface'`), O2 substeps every
  increment that crosses the smooth-cap ramp (η < c₂M; nested π_i fold). The substep level of an increment changes
  with θ and with the trial lateral strain, and O2's response jumps where it does. The mixed solve saw a jump of
  1.3e-3 p_init at σ3′ 4.9 kPa. The first smoke (apex start, plain Newton lateral solve, killed after its first
  configuration) stopped all 4 starts on xtol at costs 6.0 / 79 / 207 / 3738, recovering none. Its Jacobian had
  σ_max 4.2e6 and condition 1.9e6, and the line scans showed a jump of norm 379 in the residual under a 1e-6 change
  of N or N̄: the perturbed curve's lateral Newton failed at increment 4.
  'ramp_end' puts π_i0 on the surface through (p_init, c₂M), which has no free parameter. On it the residual is
  smooth: |Δr|/δ is constant over δ = 1e-6…1e-3 in every parameter, at the truth and at the stalled point. The
  apex rule is kept as an option.
- **The lateral solve is safeguarded** (`drivers._solve_lateral`): Newton, then bracket + Brent. O2 refuses some
  trial increments (an over-large lateral strain). The earlier "TC from the apex refuses at increment 2"
  (`local_linesearch:pi_fold`, substeps exhausted at 2^8) was such a trial, not the path: with the safeguarded
  solve, TC completes from both starts (census in `validate_drivers.json`). A C++ global Newton trial can hit the
  same refusal.
- **BA06 cannot carry one μ₀ across 4.9 and 49 kPa.** The 'per_test' policy gives every test the stiffness of its
  own confinement. A BVP needs 'global'; that choice belongs to the energy decision (§2.5), not to the harness.
- **e is e_0.05** for every Tatsuoka test (`../data/README.md` §3.1). The 49 kPa and ≥ 98 kPa states are slightly
  denser than modelled. This is stated, not corrected.
