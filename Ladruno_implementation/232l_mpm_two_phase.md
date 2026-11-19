# WP-232l — MPM two-phase particles: results

ADR 232 phase 6, PR #1116, branch `wp/232l-mpm-two-phase`. Interface `tests/mpm_twophase_interface.md`
(P1–P13, §7a, §7b, §7c and its amendment). User guide: [[232_ladruno_mpm_guide]], "Two-phase particles".
Records: `tests/data/mpm/twophase/` (`route_e_p5a_results.json`, `route_s_results.json`,
`t0_t1_mpm*.json`, `t0_iii_mpm.json`, `t10_single_phase_bitwise.json`, `r_rows_route_e.json`). Every
record carries `ladrunoBuild`.

## What was built

- **Route E (mode A):** a monolithic u-p step. Each two-phase `LadrunoMPMPoint` assembles
  LadrunoUP's blocks (Q, H, S, H̃, f_seep) at its own weights under Newmark.
- **Route S (mode F):** the explicit skeleton is USL for two-phase particles, and the pressure is
  implicit: [S* + Δt H + Δt² QᵀM⁻¹Q] p by Jacobi-PCG, then the correction. The run-time limit is
  s·Δt_d with s ≤ 0.5 (≤ 0.7 warned).
- **Shared:**
  - The particle pressure carrier P5a: FLIP on a consistent L2 projection, probe-guarded.
  - p-activity on the two-phase mass against `pressureCutoff`.
  - B-spline images per component (P9). Fix 2, the affine p image on impervious faces, is route E
    only (A7).
  - The rigid plate (P10), pfix (P8), traction (P10b), and the q + π invariant.

## Gate rows (final numbers)

| Row | Gate | Final | Build |
|---|---|---|---|
| T0 (i) particle blocks | rel 1e-9 | 4.6e-14 | 844d84212 |
| T0 (ii) three route-E steps | rel 1e-9 | 7.4e-13 | 844d84212 |
| T0 (iii) ten route-S steps (oracle re-derived for §7c, e43e0f58a) | rel 1e-9 | 8.1e-12 | 69ffff490 |
| T1 FE identity (Q4 1 cell, 2×2, H8) | A 1e-12, b 1e-12 | A ≤ 5e-16, b ≤ 1.1e-15 | 844d84212 |
| T2p Terzaghi, route E (gimp, bspline2 × off, auto) | 0.02 | worst p 0.0065, U 0.0049 (Tv ≥ 0.1) | 844d84212 |
| T3 Mandel, route E (4 legs) | bands per tv; peak ≥ 1.039 | pass; peaks 1.080–1.086 | 844d84212 |
| T4p hydrostatic (h = H/20) | 0.01 (p) | gimp 1.0e-4, bspline2 1.3e-5 | 844d84212 |
| Fix-2 patch check (d064f5b7e rule: grid coefficients + field at the particles) | 1e-10 | 5.3e-14 / 4.2e-14 (particle p 6.0e-6, recorded) | 844d84212 |
| T5 checkerboard, route E | auto ≤ 0.3; off ≥ 0.5 | 0.134 / 0.621 | 844d84212 |
| T5 checkerboard, route S (A5: undamped) | ≤ 0.3 | 0.137 (off reported 0.614) | 69ffff490 |
| T6 Terzaghi, route S, Δt 0.005 (bspline2 gated, gimp reported) | 0.02 | bspline2 p ≤ 0.0037, U ≤ 0.0006; gimp p ≤ 0.0029, U ≤ 0.0016 | 69ffff490 |
| T6 Mandel, route S, Δt 0.005 | off 0.0214/0.0263/0.02…, auto 0.0227/0.02…; peak ≥ 1.039 | bspline2 worst 0.0162 / 0.0166, peaks 1.089 / 1.087; gimp 0.0135 / 0.0141 | 69ffff490 |
| T7 route-S Δt limit | stable at 0.9 × drained pencil, unstable or refused at 1.5 × | **recorded miss** (strict xfail); see below | 5d64922a4, 69ffff490 |
| T8 refusals (incl. rank-deficient, stray p, massCutoff > pressureCutoff, route-S cap) | pass | pass | 69ffff490 |
| T9 restart bit for bit | pass | pass | 69ffff490 |
| T10 single-phase bit for bit against cfe478254 | 24 cases | 24/24 | 69ffff490 |
| T11 mutation (MPM ZERO ≥ 0.500, SCALE ≥ 0.398; IDENT turns T0 (i) / T1 red) | — | pending: CI lane on the final head | — |
| T12 engines agree (Tcl, Python) | pass | pass | 69ffff490 |
| Two-phase test files | — | 69 passed, 2 xfailed (T7) | 69ffff490 |

Mandel route S, amplification of the step map at Δt 0.005 (bspline2): ρ = 0.999999998 after §7c,
against 1.00093 before (4b74bd335).

## Recorded misses and owner decisions

- **2026-11-18, route E on fd81cbdd3:** the row-sum-lumped projection biased the particle p (quirk).
  The owner fixed the carrier in this WP (P5a, consistent projection); fixes 2 and 3 were ruled
  as defects.
- **2026-11-18, fix-2 patch check:** three rounds.
  - π's invisible part: A′ and the q + π invariant.
  - u-side renormalisation rescaled the p weights (bb1421559).
  - P5a's motion residual (6.0e-6). **Owner (d064f5b7e):** the check gates the grid coefficients
    and the field at the particles; the particle-carried p is recorded and gated by T4p at 0.01.
- **2026-11-18, Fable carrier audit:** four findings, fixed with tests.
  - The dangling `else` (T2p bspline2 was void before the fix).
  - A′ and π_c commit only on an accepted step.
  - π is evaluated at the step's drive, with A′ at λ(t₀).
  - Mode A refuses massCutoff > pressureCutoff.
- **2026-11-18, A7 (accepted, b5234a259):** route S refuses fix 2's term (an impervious B-spline face
  under a normal fluid drive). A follow-up issue opens at landing.
- **2026-11-19, T6 Mandel (route S):**
  - The diagnosis: on bspline2, a linear instability (ρ 1.00093 at 0.30 × the drained pencil),
    because MUSL's v** fed the skeleton while the fluid balance saw v₁. On gimp, the light-node
    correction.
  - **Owner:** the consistency fix (USL skeleton strain from v₁ for two-phase particles).
  - The orchestrator withdrew its 3d3699a50 ruling: light u nodes are dropped whole in two-phase
    steps.
- **2026-11-19, T7 — recorded miss (owner: "record miss, cap at 0.5"):**
  - On T6's column with the USL skeleton, route S is stable to 0.7 × the drained pencil on gimp and
    0.8 × on bspline2. It blows up at 0.8 × / 0.9 ×.
  - K_f = 1 gives the same thresholds, so the u-p coupling is not the cause.
  - Route S's limit is s·Δt_d, with s ≤ 0.5 by default, up to 0.7 warned, and refused above.
  - Δt_d is computed per epoch and on a cell change. The Mandel step dropped from 12.3 ms to 2.95 ms.

## R rows — route E (reported, never gated; `r_rows_route_e.json`, build 5d64922a4)

- **R1, first settlement** u(0⁺)/u_undrained, MPM against LadrunoUP at the same Δt.
  - T2's soil (K_f 4444), gimp and bspline2: 1.008–1.046 across α₀ off / 0.1 / 0.25 and
    Δt 1e-3 / 0.02 / 0.2. LadrunoUP gives 1.011–1.047, and every pair is within 0.004.
  - Stiff fluid (K_f 2.2e6): gimp tracks LadrunoUP to 0.05 (2.01 against 2.06 at off, Δt 1e-3;
    25.76 against 25.76 at auto 0.25). bspline2 sits 1–4 % lower (1.57 at off, Δt 1e-3).
- **R2, projection diffusion** (P5a, α_p 1, no motion): particle against exact ≤ 1e-3 of the
  amplitude over 2 s. Particle against field: 4.9e-4.
- **R3, cell-crossing meter** (T4's column translated 4 cells on 232k's carrier block): p error up to
  0.018 γ_wH during the translation.
- **R4, cost per step:**

  | Configuration | ms per step |
  |---|---|
  | Route E gimp (40 particles), Linear | 0.60 |
  | Route E gimp, Newton | 1.24 |
  | Route E bspline2 (80 particles), Linear | 2.10 |
  | Route E bspline2, Newton | 4.39 |
  | LadrunoUP Q4 1 × 10, Linear | 0.063 |

- **R5, undrained simple shear:** not run. Particles have no drainage switch (no way to release a
  `pfix` mid-run), and per the protocol none was invented.

## R rows — route S

*(Placeholder: 232l-V's route-S reported rows, run on a dist copy, go here.)*

## Known limits

- **A7:** route S has no fix-2 term.
- **T7:** the cap at s 0.5, up to 0.7 warned.
- There is no drainage switch, and porosity is constant.
- Serial: no new threaded code.
- The particle-pressure motion residual of P5a is recorded, not damped.
