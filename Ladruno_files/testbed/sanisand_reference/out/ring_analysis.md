- **paper**: elastic samples with f > 1e-6 kPa: 0; max |f| at exit over ok runs: 1.7e-07 kPa.
- **paper**: 480 runs; statuses {'ok': 476, 'solver_failed': 3, 'p_floor': 1}; runs with a plastic sample at (α−α_in):n < −1e-12 (h < 0): **0**; admissible starts 468, of which max ρ_α > 1: 0 (α itself moved in 0; the rest are elastic paths on which α^b(ψ) contracted past a fixed α), max ρ_α 0.983.
- **uw_model**: elastic samples with f > 1e-6 kPa: 0; max |f| at exit over ok runs: 1.7e-07 kPa.
- **uw_model**: 480 runs; statuses {'ok': 477, 'solver_failed': 3}; runs with a plastic sample at (α−α_in):n < −1e-12 (h < 0): **0**; admissible starts 468, of which max ρ_α > 1: 0 (α itself moved in 0; the rest are elastic paths on which α^b(ψ) contracted past a fixed α), max ρ_α 0.983.
- **uw_rule**: elastic samples with f > 1e-6 kPa: 0; max |f| at exit over ok runs: 2.9e+02 kPa.
- **uw_rule**: 480 runs; statuses {'ok': 426, 'solver_failed': 27, 'H_nonpositive': 27}; runs with a plastic sample at (α−α_in):n < −1e-12 (h < 0): **146**; admissible starts 468, of which max ρ_α > 1: 25 (α itself moved in 25; the rest are elastic paths on which α^b(ψ) contracted past a fixed α), max ρ_α 1.965.

C++ vs `uw_model` on admissible starts with a reference status ok (‖Δσ_C++ − Δσ_ref‖/‖Δσ_ref‖):

| C++ | substeps ≤ 7 (err = 0 path): n / median / p95 | substeps > 7: n / median / p95 | escapes ρ_α > 1 (≤ 7 / > 7) |
|---|---|---|---|
| ME | 248 / 4.3e-01 / 4.1e+00 | 220 / 6.9e-03 / 8.3e-01 | 25 / 0 |
| ME8 | 192 / 3.9e-01 / 1.3e+00 | 276 / 1.1e-02 / 8.3e-01 | 7 / 1 |

Inadmissible start rows (ρ_α > 1):

| row | ρ_α0 (ρ_b0) | f0 | probe | δ | ref uw_model: status, max ρ_α, end ρ_α, η_end, p_end, f_end | C++ ME: rc, substeps, ρ_α end, η_end, f_after |
|---|---|---|---|---|---|---|
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | isoComp | 1e-05 | ok, 7.28, 6.17, 10.94, 0.537, 1.2e-09 | 0, 1, 4.40, 7.86, 1.2e-09 |
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | shear | 1e-05 | solver_failed, 7.28, 7.28, 12.86, 0.352, 7.6e-10 — solver failed at t=0.0544221, p=0.3519 kPa, mode plastic, (alpha-alpha_in):n=1.386e-10, b:n=3.728e-07, Hs=1.666e-04, rho_b=6.321 | 0, 1, 7.14, 12.66, 1.4e-02 |
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | isoComp | 1e-04 | ok, 7.28, 5.40, 9.59, 1.08, 2.3e-09 | 0, 1, 0.92, 1.79, 5.5e-09 |
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | shear | 1e-04 | solver_failed, 7.28, 7.28, 12.86, 0.352, 7.7e-10 — solver failed at t=0.00544221, p=0.3519 kPa, mode plastic, (alpha-alpha_in):n=1.387e-10, b:n=3.895e-07, Hs=1.740e-04, rho_b=6.321 | 0, 1, 5.80, 10.33, 4.4e-10 |
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | isoComp | 1e-03 | ok, 7.28, 4.35, 7.84, 3.71, 8.1e-09 | 0, 1, 0.29, 0.47, 4.9e-08 |
| b8 1950/3 | 7.28 (6.26) | 7.5e-10 | shear | 1e-03 | solver_failed, 7.28, 7.28, 12.86, 0.352, 7.5e-10 — solver failed at t=0.000544221, p=0.3519 kPa, mode plastic, (alpha-alpha_in):n=1.387e-10, b:n=3.919e-07, Hs=1.751e-04, rho_b=6.321 | 0, 1, 9.99, 17.70, -1.4e-09 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | isoComp | 1e-05 | ok, 6.81, 5.76, 10.22, 0.545, 1.4e-10 | 0, 1, 4.12, 7.36, 1.3e-10 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | shear | 1e-05 | ok, 6.81, 6.70, 11.86, 0.353, 9.2e-11 | 0, 14, 6.70, 11.85, 3.1e-09 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | isoComp | 1e-04 | ok, 6.81, 5.02, 8.94, 1.11, 2.9e-10 | 0, 1, 0.86, 1.69, 5.8e-10 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | shear | 1e-04 | ok, 6.81, 6.16, 10.93, 0.346, 9.0e-11 | 0, 1, 5.45, 9.72, 3.0e-02 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | isoComp | 1e-03 | ok, 6.81, 4.01, 7.25, 4.28, 1.1e-09 | 0, 1, 0.29, 0.47, 5.1e-09 |
| b8 1950/2 | 6.81 (5.86) | 7.8e-11 | shear | 1e-03 | ok, 6.81, 4.29, 7.71, 0.276, 7.2e-11 | 0, 1, 6.79, 12.18, 1.4e-10 |
