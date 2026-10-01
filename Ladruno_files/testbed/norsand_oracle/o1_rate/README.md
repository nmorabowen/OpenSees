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
| `integrator.py` | `State` and `integrate_increment`: one increment, strain control or mixed control. `kin="log"` switches to log strain with v = v₀J. |
| `api.py` | `initial_state`, `run_path`, `tangent` (the continuum tangent: loading branch if the last segment was plastic), and `triaxial` (axial = x, signed strain). |
| `localization.py` | `acoustic_min_det`, built from (S.34) + (S.44) with `finite=True`, or from the raw tangent with `finite=False`. Also `k2_path` (S.43). |
| `selfcheck.py` | The self-checks (not the gate suite). |

## Running the self-checks

From `Ladruno_files/testbed/norsand_oracle`:

    python -u -m o1_rate.selfcheck [quick] [basic tangents dissipation rtol caps undrained_cs drained_cs k2 k2_sensitivity]

With no group named, every group runs. The whole set takes under 1 minute on 24 cores.

## Choices and findings against the sheet

1. **§12, a^e.** a^e is the tensor derivative of σ = p1 + (2/3)(q/ε_s)e. It matches (S.3)+(S.33) to 1e-15.
2. **§3.2, the vertex.** With no cap, near-isotropic compression drives R → 0 in finite time, after which the vertex rule has no consistent rate solution. O1 stops with `vertex_reached` / `vertex_nonisotropic`. A pure hydrostatic path is fine.
3. **§10.1, the planar cap.** The state can be attracted to the corner at η = χ_cap·M (Filippov sliding), where there is no classical rate solution. O1 stops with `cap_sliding`. The smooth cap runs through.
4. **§4.2, WW at ρ = ½ exactly.** ζ = 2cosθ, so ζ'(π/3) = −√3 ≠ 0: the section has a vertex at the compression corner. Owner decision 2026-10-01: O1 **refuses** it (`ValueError`). The admissible WW range is (½, 1] for both ρ and ρ̄; GA stays [7/9, 1]. Before this decision O1 only warned. The ρ = ½ branch of `model.zeta_fun` is kept for direct shape evaluation (the K1.3 corner print), but `Params` can no longer reach it.
5. **§14, the interpolated criterion.** The step index is round(n*); `n_interp` is also returned.
6. **§11.** The ρ ranges are applied to both ρ and ρ̄.
