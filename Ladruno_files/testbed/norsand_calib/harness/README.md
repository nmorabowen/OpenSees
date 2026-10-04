# norsand_calib/harness: the WP-144 P3 calibration harness (Zone B)

Fits LadrunoNORSAND's (χ, h, N, N̄, ρ, ρ̄) to element-test curves through **O2**, the algorithmic oracle the C++
kernel matches to ~1e-10 ([`../../norsand_oracle/o2_algo`](../../norsand_oracle/o2_algo/README.md)). Plan:
[`144_ladruno_norsand_plan.md`](../../../../Ladruno_implementation/144_ladruno_norsand_plan.md) §3 P3 and §5.2 K3
(re-aimed 2026-10-02: Toyoura, Tatsuoka et al. 1986 drained plane strain, DM04 Table 1 CSL). Equations:
[`144a_norsand_equation_sheet.md`](../../../../Ladruno_implementation/144a_norsand_equation_sheet.md).
numpy/scipy, so it runs on Esmeralda (`~/wp144/venv`), from `norsand_calib/`.

| File | What it does |
|---|---|
| `sand.py` | What the fit does not touch: the fork CSL (S.22), M, p_a, the DM04 elastic targets G(p, e), ν. `toyoura_dm04()` reads the data pack's DM04 Table 1 transcription (`../data/dm04/toyoura_table1.csv`) at **p_a = `TIMS_P_A` = 101 kPa** by default (round 3b A4: the campaign's `Patm 101`, ONE flag for the HAR energy and the fork CSL; was the pack's 100 kPa assumption). `TOYOURA_PLACEHOLDER` (as recalled, p_a 101) is the smoke's set only |
| `energy.py` | The energy plug. `get("BA06")` maps the targets to (p₀, κ̂, μ₀, α₀ = 0) at a representative p (sheet §15). `get("HAR")` (round 3b, bound to O2's and O1's `energy='HAR'`, `k, g, n_e, p_a`): n = the targets' exponent (0.5), **k = K(p_a, e)/p_a, g = G(p_a, e)/p_a**, so the oracle's own G(p), K(p) equal the DM04 targets at every p (checked to 1e-15) with no representative pressure and the axis ν unchanged; none of the five BA06 values is passed (both oracles refuse them under HAR, sheet §2.4). `ElasticPolicy`: 'per_test', or **'global' for a BVP** (`ElasticPolicy.bvp(p_rep, e_rep)`: one set of parameters for every test of the body; HAR needs only `e_rep`); `tcl_flags(sand, policy, e)` writes the same energy as `nDMaterial LadrunoNorSand` flags (`-energy HAR -k -g -n`, or `-p0 -kappa_hat -mu0`) for a deck. No other file names an energy |
| `model.py` | Fitted vector + `Setup` (sand, energy, policy, WW, smooth cap c₁ 0.05 / c₂ 0.15, π_i0 rule, `p_min`) → oracle `Params` and the isotropic initial state. π_i0 rules: **'unified' (default, round 3b: sheet §5.4 (S.53), computed by the ORACLE's own helper `initial_state(pi_i0=None)`; the surface through (p_init, max(η_init, c₂M)), which for these isotropic starts is the old 'ramp_end' surface to 2e-16)**, 'ramp_end' (the harness's own closed form of the same surface), 'on_surface' (the apex through the isotropic stress, π_i0 = p (1−N)^((1−N)/N)), 'ratio'. An initial state the oracle refuses (§7 guard B ≤ 0, η* ≥ M/N) becomes an incomplete curve with the refusal penalty, not an exception. `Setup.p_min` = None leaves the oracle's floor (O2 5e-3 p_ref, O1 off) |
| `drivers.py` | O2 element drivers PS, TC, TE (drained, mixed control on O2's step: Newton with its consistent tangent, then bracket + Brent), PSU, TCU (constant volume, TIMs rung C6); `run_o1` is the O1 counterpart |
| `data.py` | Curve CSV schema `eps_a_pct, sr, eps_v_pct` + metadata (docstring), point-test CSV for φ_peak(e), ε_peak(e) |
| `objective.py` | Residuals in sr and ε_v up to the data peak + 2 % (post-peak weight 0.5), plus peak φ and ε_peak; σ's stated in its docstring (`Weights`) |
| `fit.py` | The refusal rules as a box reparametrisation; `multistart` (least_squares TRF, scrambled Sobol starts, processes); `identifiability` (Jacobian in ln θ: singular values, condition, correlation); `profile` |
| `validate_drivers.py` | Drivers vs O1, first order (G1 criteria); K1.7 undrained endpoint; the start census (apex / ramp_end / unified); `--energy HAR` runs the same gate with the HAR energy on both oracles |
| `check_plugs.py` | Seconds-long wiring checks (round 3b): HAR plug vs the DM04 targets through the oracle's own elastic law, BA06 plug numbers, the unified π_i0 helper vs the closed form on both oracles (72 cases), the 'global' BVP policy, a short PS run under both energies vs O1 |
| `smoke_recovery.py` | The harness gate: recover known parameters from an O2-generated synthetic set. `--variant clean` (the gate), `noisy` (**round 3b: the inverse-crime caveat removed**: generated at a different step and with Gaussian noise of the objective's σ), `both` |
| `fit_tatsuoka.py` | The K3 fit entry point on the data pack |
| `run_smoke.sbatch` | The smoke on a compute node (`--exclude=node5`; `SMOKE_ARGS` picks the variant) |

Outputs land in `../out/` (`validate_drivers.json`, `smoke_recovery.json`, `synthetic/*.csv`).

## Run

```bash
python -m harness.check_plugs                 # seconds
python -m harness.validate_drivers [--energy HAR]            # ~1 min
python -m harness.smoke_recovery --starts 8 --workers 8      # clean gate; or: sbatch harness/run_smoke.sbatch
python -m harness.smoke_recovery --variant noisy --no-profiles --starts 8 --workers 8   # noisy + different step
# a Windows working tree has CRLF line endings: strip them from run_smoke.sbatch on the cluster (sed -i 's/$//') before sbatch
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

## Round 3b (Esmeralda login node, nice 19, 2026-10-03; logs and JSON in `../out/*r3b*`)

Changes: π_i0 default = the unified rule (S.53) through the oracle's own helper; HAR plug bound to O2/O1 (`energy='HAR'`, k, g, n_e, p_a = 101 kPa); `ElasticPolicy.bvp()`;
`toyoura_dm04()` at p_a 101; `p_min` pass-through; a noisy recovery variant. The compute-node submission was not available to this session, so the smokes ran niced on the login
node (4 workers), as on 2026-10-02.

- **`check_plugs`** PASS (seconds): the oracle's own G(p), K(p) equal the DM04 targets at p = 1…300 kPa to 1.1e-15 under HAR for both policies and two void ratios, axis ν = 0.05 to 3e-17;
  the BA06 numbers are unchanged; the unified π_i0 of both oracles equals the closed form to 2.2e-16 (72 cases: O1/O2, BA06/HAR, three caps, N = 0 and 0.3); the 'global' policy gives one
  parameter set to every test; a 3 % PS run is first order against O1 under both energies (σ error halves from n = 60 to 120: 4.0e-3 → 2.1e-3 BA06, 4.5e-3 → 2.3e-3 HAR).
- **`validate_drivers`**: BA06 reproduces the 2026-10-02 table digit for digit (PS e_σ 9.799e-4 at n = 400, order 0.972, …); **HAR passes the same gate** (PS 0.972 / 0.972 / 0.953, TC 0.986 / 0.986 / 0.937,
  TE 0.994 / 0.992 / 0.976, PSU 1.004 / 1.003, TCU 0.997 / 0.997; K1.7 endpoint p within 8e-5 (PSU) and 5.5e-4 (TCU) of p_cs, ζq/|p| within 1.7e-6 and 9.7e-6 of M). Start census: the
  'unified' start substeps exactly where 'ramp_end' did (PS and TC at increments 0 and 1).
- **Clean recovery smoke PASS, unchanged by the round**: `six` best cost 1.8e-22, max rel error 2.6e-13, 2 of 4 starts recovered (the other two stopped at costs 7.6 and 39.8); `rho_pinned`
  6.4e-22, 1.2e-12, 1 of 4 (others 10.7, 1.1e4, 1.6e4); Jacobian conditions 2.0e2 / 1.4e2; the synthetic data rows are byte-identical to 2026-10-02 (only the `# source` header differs); 5237 s.
  This is the unified π_i0 start: it is the same surface as 'ramp_end', so nothing moved.
- **Noisy recovery smoke PASS** (`smoke_recovery_noisy_r3b.json`, 3335 s; **the inverse-crime caveat is removed**). The set is generated at 0.025 % axial steps while the model fits at 0.05 %, and carries Gaussian
  noise of the objective's own σ (σ_sr 0.10, σ_ev 0.10 %, σ_φ 0.5°, σ_ε,peak 0.25 %; seed 1; the x = 0 row exact). The step difference alone puts |r(truth)| = 2.1 (rms 0.155σ over 188 residuals) under a
  noise-free set, small against the noise. Gate: Mahalanobis² of ln(θ_fit/θ_true) under the Gauss–Newton sandwich covariance at the truth ≤ χ²(0.999, k), and cost ≤ cost at the truth.

  | Config | Mahalanobis² (limit) | χ² fit / at truth (expected Σvar) | max rel error | worst z | Starts consistent | Jacobian condition |
  |---|---|---|---|---|---|---|
  | `six` | 16.3 (22.5) | 117.4 / 120.6 (141.8) | 10.6 % (N) | −0.88 (χ) | 2 of 4 | 1.7e2 |
  | `rho_pinned` | 6.2 (20.5) | 117.8 / 120.6 (141.8) | 10.0 % (N) | +0.49 (N) | 1 of 4 | 1.1e2 |

  Relative SDs of ln θ at this noise level (`six`): χ 1.2 %, h 5.7 %, ρ 2.6 %, ρ̄ 4.0 %, N̄ 14 %, N 20 %. The noisy fit is therefore good to a few % in χ, h, ρ, ρ̄ and not better than the noise in N (10 % off, z 0.5).
  The other starts again stopped in local minima (costs 60–61 and 1.5e4–2.0e4). χ² at the truth sits below its expectation (120.6 vs 141.8) because the realisation is mild; the fit is below the truth, as it must be.
- **One full-data evaluation under each energy** (`fit_eval_r3b_BA06/HAR.log`, χ = −3, h 150, N 0.3, N̄ 0.2, ρ̄ 0.75, ρ 0.712; 4 curves + 22 point tests, p_a 101): BA06 cost 8830, HAR 8684; every curve and point test completes under both. These are not fits.

## Findings that bind the fit (Esmeralda, 2026-10-02)

- **The default start is the unified rule (S.53) = 'ramp_end' for isotropic starts, not the apex.** From the apex (`pi0_rule='on_surface'`), O2 substeps every
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
