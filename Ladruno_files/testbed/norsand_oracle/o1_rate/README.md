# O1: continuum rate oracle for LadrunoNORSAND (WP-144)

This is the truth reference (plan §5.1). It takes the rate equations of equation sheet 144a §12
and integrates them with `scipy.integrate.solve_ivp(method="Radau", rtol=1e-10)`, with atol
scaled per component. It has no return map and no projection:

- The closed-form multiplier (S.41) keeps F = 0 to the ODE tolerance. `flags['max_F_rel']`
  reports the drift, which is about 1e-14.
- The elastic-to-plastic switch is a terminal event at F = 0.
- Unloading, a zero denominator (a stop, never mapped to elastic), the B > 0 guard, cap
  switches and the vertex are events too.

It was written from the sheet alone, independently of O2.

## Files

| File | What it holds |
|---|---|
| `params.py` | `Params` (sheet §1.3) and `validate()`, which applies the refusals of §4 and §11 (condition A). Pass `check=False` to bypass them. |
| `model.py` | The energy (S.2)–(S.5), in coordinate-free tensor form, and ζ (S.8)–(S.11). Also F, Q_A, the cap (S.35)–(S.36), the vertex rule (§3.2), the CSL, π_i* and H (S.22)–(S.25), and a^ep (S.42). |
| `integrator.py` | `State` and `integrate_increment`: one increment, strain control or mixed control. Specific volume is algebraic, v = v₀ exp(tr ε) (G2 owner decision 2026-10-01, sheet §1.2/(S.40)); under `kin="log"` (log strain) that is v = v₀J. |
| `api.py` | `initial_state`, `run_path`, `tangent` (the continuum tangent: loading branch if the last segment was plastic), and `triaxial` (axial = x, signed strain). |
| `localization.py` | `acoustic_min_det`, built from (S.34) + (S.44) with `finite=True`, or from the raw tangent with `finite=False`. Also `k2_path` (S.43). |
| `selfcheck.py` | The self-checks (not the gate suite). |

## Running the self-checks

From `Ladruno_files/testbed/norsand_oracle`:

    python -u -m o1_rate.selfcheck [quick] [basic tangents dissipation rtol caps har floor pi0 undrained_cs drained_cs k2 k2_sensitivity]

With no group named, every group runs. The whole set takes under 1 minute on 24 cores.

## Choices and findings against the sheet

1. **§12, a^e.** a^e is the tensor derivative of σ = p1 + (2/3)(q/ε_s)e. It matches (S.3)+(S.33) to 1e-15.
2. **§3.2, the vertex.** With no cap, near-isotropic compression drives R → 0 in finite time, after which the vertex rule has no consistent rate solution. O1 stops with `vertex_reached` / `vertex_nonisotropic`. A pure hydrostatic path is fine.
3. **§10.1, the planar cap.** The state can be attracted to the corner at η = χ_cap·M (Filippov sliding), where there is no classical rate solution. O1 stops with `cap_sliding`. The smooth cap runs through.
4. **§4.2, WW at ρ = ½ exactly.** ζ = 2cosθ, so ζ'(π/3) = −√3 ≠ 0: the section has a vertex at the compression corner. Owner decision 2026-10-01: O1 **refuses** it (`ValueError`). The admissible WW range is (½, 1] for both ρ and ρ̄; GA stays [7/9, 1]. Before this decision O1 only warned. The ρ = ½ branch of `model.zeta_fun` is kept for direct shape evaluation (the K1.3 corner print), but `Params` can no longer reach it.
5. **§14, the interpolated criterion.** The step index is round(n*); `n_interp` is also returned.
6. **§11.** The ρ ranges are applied to both ρ and ρ̄.
7. **§1.2 / (S.40), G2 owner decision 2026-10-01.** The specific volume is the exponential update v = v₀ exp(tr ε), kept algebraic from the integrated total strain (anchored at the increment start, v = v_n exp(tr ε − tr ε_n)); the ODE copy v̇ = v tr ε̇ is carried in y[7] as a diagnostic (`flags['v_ode']`). The pre-G2 small-strain rule v̇ = v₀ tr ε̇ is gone; `kin` no longer changes v. Self-check K1.10 (in the `drained_cs` group) prints the identity.
8. **§2.3–§2.4, `energy='BA06'|'HAR'` (round 3b).** BA06 stays the default. HAR (`k`, `g`, `n_e`, `p_a`) follows (S.4h)–(S.5h''), with (S.3) used in full. Under HAR, any BA06 value given (`p0 kappa_hat eps_v0 mu0 alpha0`) is refused, and so are `k g n_e` under BA06. `p_a` is the one parameter shared with the fork CSL. The TIMs value is 101 kPa. The dataclass default stays 101.325 only so that no existing fork-CSL result moves, so a HAR run must pass `p_a` explicitly. `p_ref` is |p₀| (BA06) or p_a (HAR). A HAR state leaving dom Ψ is a floor event when the floor is on, and the `p_to_zero` stop when it is off.
9. **§9.7, the floor `p_min` (rate form (S.55)).** `p_min=0` (the default) means off. `"default"` means 5·10⁻³ p_ref. A float is used as given; a negative value is refused. On the floor, O1 integrates the second mechanism ε̇^f = λ̇_f·1/3 with ṗ = 0. The modes are `floor` and `plastic_floor` (Koiter, 2×2). The modes are entered at the `floor_hit` event and left at `floor_release` (λ̇_f → 0⁻) or `unload`. The tie-break of §12 is extended to the second multiplier. Counters are `flags['eps_f_v']`, `W_f = p_min·ε^f_v`, `at_floor`, `floor_segments` and `n_floor_increments`, plus `n_f_init` from `initial_state`. In the rate form the floor's stored-energy gain is exactly W_f, because |p| = p_min on the floor. `model.floor_project`/`floor_target`/`floor_phi`/`floor_energy` are the split operator Π_f (S.48)–(S.51a), with (S.52) E_f available on demand. They are used only by `initial_state` (the init floor) and by the self-checks. On the integration path, O1 never projects. `tangent()` adds the floor mechanism when the last segment ended on the floor. Owner decision (c) holds: the tangent is exact and 1:C = 0, with no regularisation.
10. **§5.4, (S.53).** `initial_state(..., pi_rule="S53")` builds the surface through (p_init, max(η_init, c₂M)). c₂ is the smooth cap's c₂, c₁ for the planar cap (χ_cap), and 0 for no cap. The state is refused if η* ≥ M/N, and the B > 0 guard is checked. The default stays `pi_rule="surface"` (the pre-round-3 O1 rule), so that existing results do not move. With cap = none, the two rules coincide off the axis. An explicit `pi_i0` overrides both.
11. **Regression.** With the defaults (BA06, `p_min=0`, `pi_rule="surface"`), O1 is bit-identical to the pre-round-3 O1 (428aa1489). This was checked on 17 paths: K1.1, K1.2, eight triaxial cases (paper/fork, α₀ = 5, GA), the three cap modes, the surface-rule start and the K2 path, plus a^ep and a^e.
