# WP-232l — MPM two-phase particles: results

ADR 232 phase 6, PR #1116, branch `wp/232l-mpm-two-phase`. Interface `tests/mpm_twophase_interface.md`
(P1–P13, §7a, §7b, §7c and its amendment). User guide: [[232_ladruno_mpm_guide]], "Two-phase particles".
Warrant build: **d77caf8fdfbf37a244ab7e2a58e0ccef926b3d3b** (release candidate, 2026-11-19). Records: `tests/data/mpm/twophase/` (`route_e_p5a_results.json`, `route_s_results.json`,
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
| T0 (i) particle blocks | rel 1e-9 | 4.6e-14 | d77caf8fd |
| T0 (ii) three route-E steps | rel 1e-9 | 7.4e-13 | d77caf8fd |
| T0 (iii) ten route-S steps (oracle re-derived for §7c, e43e0f58a) | rel 1e-9 | 8.1e-12 | d77caf8fd |
| T1 FE identity (Q4 1 cell, 2×2, H8) | A 1e-12, b 1e-12 | A ≤ 5e-16, b ≤ 1.1e-15 | d77caf8fd |
| T2p Terzaghi, route E (gimp, bspline2 × off, auto) | 0.02 | worst p 0.0065, U 0.0049 (Tv ≥ 0.1) | d77caf8fd |
| T3 Mandel, route E (4 legs) | bands per tv; peak ≥ 1.039 | pass; peaks 1.080–1.086 | d77caf8fd |
| T4p hydrostatic (h = H/20) | 0.01 (p) | gimp 1.0e-4, bspline2 1.3e-5 | d77caf8fd |
| Fix-2 patch check (d064f5b7e rule: grid coefficients + field at the particles) | 1e-10 | 5.3e-14 / 4.2e-14 (particle p 6.0e-6, recorded) | d77caf8fd |
| T5 checkerboard, route E | auto ≤ 0.3; off ≥ 0.5 | 0.134 / 0.621 | d77caf8fd |
| T5 checkerboard, route S (A5: undamped) | ≤ 0.3 | 0.137 (off reported 0.614) | d77caf8fd |
| T6 Terzaghi, route S, Δt 0.005 (bspline2 gated, gimp reported) | 0.02 | bspline2 p ≤ 0.0037, U ≤ 0.0006; gimp p ≤ 0.0029, U ≤ 0.0016 | d77caf8fd |
| T6 Mandel, route S, Δt 0.005 | off 0.0214/0.0263/0.02…, auto 0.0227/0.02…; peak ≥ 1.039 | bspline2 worst 0.0162 / 0.0166, peaks 1.089 / 1.087; gimp 0.0135 / 0.0141 | d77caf8fd |
| T7 route-S Δt limit | stable at 0.9 × drained pencil, unstable or refused at 1.5 × | **recorded miss** (strict xfail); see below | d77caf8fd |
| T8 refusals (incl. rank-deficient, stray p, massCutoff > pressureCutoff, route-S cap) | pass | pass | d77caf8fd |
| T9 restart bit for bit | pass | pass | d77caf8fd |
| T10 single-phase bit for bit against cfe478254 | 24 cases | 24/24 | d77caf8fd |
| T11 mutation: MPM IDENT (joint Q + projection cut) | the two-phase oracle and T1 go red | 14/14 red (`_twophase_oracle` + `_fe_identity` on `dist/mut_MPM_IDENT`; green on the normal build) | d77caf8fd |
| T11 mutation: MPM ZERO ≥ 0.500, SCALE ≥ 0.398 | CI lane | pending: run 44532817830 | d77caf8fd |
| T13 translation patch (free block, both routes, gimp / bspline2, with / without a plate) | v, σ′, p rel 1e-10 | pending: the review-of-record build | — |
| T12 engines agree (Tcl, Python) | pass | pass | d77caf8fd |
| Two-phase test files | — | 69 passed, 2 xfailed (T7) | d77caf8fd |
| The whole MPM suite (`tests/test_ladruno_mpm_*.py` + the oracle's self-test) | — | 509 passed, 40 skipped (opt-in slow tiers, MUMPS matrix types absent), 14 xfailed (12 pre-registered from 232f–232k + T7 × 2), 0 failed | d77caf8fd |

Mandel route S, amplification of the step map at Δt 0.005 (bspline2): ρ = 0.999999998 after §7c,
against 1.00093 before (4b74bd335); measured on 69ffff490, whose route S equals d77caf8fd's.

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

- **2026-11-19, ready-flip review of record on d77caf8fd (Fable BLOCK + /code-review high):**
  - **Blocker:** route S could not translate a two-phase body. §7c (1) dropped light u nodes
    without renormalising. Corrected: mode A's `renormalize` is applied over each particle's
    u-active nodes in every route-S step.
  - **Plate finding (confirmed with probe_plate.py):** with `plate ymax 0` on bspline2 and uniform
    v_y = −1, both modes were off by 7 % in v and 14 kPa in σ′. The lumped mass of the plate's
    normal component used the unsigned W while the momentum used the affine image's signed weights.
    Fixed: that component's mass is Σ S′_c m_p on the mirror rows and the plate.
  - New gate T13. The gates are now asserted by pytest (`tests/test_ladruno_mpm_twophase_gates.py`).

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

`r_rows_route_s.json` (232l-V, 2d4a95129, build 69ffff490; route S is unchanged since then).

- **R1, first settlement** u(T_r)/s0 at the end of a 20-transit ramp. Route S runs at
  0.49 × Δt_d and reads the settlement as the particle quadrature.
  - T2's soil (closed form under the ramp: 1.077):

    | Leg | LadrunoUP | Route E gimp / bspline2 | Route S gimp / bspline2 |
    |---|---|---|---|
    | off | 1.088 | 1.088 / 1.085 | 1.083 / 1.081 |
    | auto 0.1 | 1.091 | 1.091 / 1.088 | 1.085 / 1.084 |
    | auto 0.25 | 1.095 | 1.095 / 1.092 | 1.090 / 1.088 |

  - Stiff fluid, K_f 2.2e6 (closed form 14.61):

    | Leg | LadrunoUP | Route E gimp / bspline2 | Route S gimp / bspline2 |
    |---|---|---|---|
    | off | 14.60 | 14.42 / 13.92 | 15.37 / 14.63 |
    | auto 0.1 | 21.37 | 21.20 / 20.59 | 22.18 / 21.43 |
    | auto 0.25 | 28.89 | 28.83 / 28.29 | 29.82 / 29.11 |

  - Stabilization inflates the first settlement on the stiff fluid in every code, as TIMs reported.
    The routes track LadrunoUP.
- **R2, cost and Δt** on T2's column:

  | Code | Δt run | ms per step | ms per simulated s |
  |---|---|---|---|
  | Route S gimp (s 0.5 admits 0.00673; Δt_d 0.01346, undrained 0.00952) | 0.0066 | 0.19 | 29.1 |
  | Route E gimp | 0.02 | 0.59 | 29.7 |
  | Route S bspline2 (Δt_d 0.01775) | 0.0087 | 0.72 | 83.2 |
  | Route E bspline2 | 0.02 | 2.06 | 103.1 |
  | LadrunoUP Q4 1 × 10 | 0.02 | 0.064 | 3.2 |

- **R3, route comparison:**
  - First settlement on T2's soil (auto 0.25, ramp): closed form 1.077, LadrunoUP 1.095,
    route E 1.095 / 1.092, route S 1.090 / 1.088.
  - B4 CB_lap (auto / off): route E 0.134 / 0.621, route S 0.137 / 0.614, LadrunoUP 0.126 / 0.603.

## Known limits

- **A7:** route S has no fix-2 term.
- **T7:** the cap at s 0.5, up to 0.7 warned.
- There is no drainage switch, and porosity is constant.
- Serial: no new threaded code.
- The particle-pressure motion residual of P5a is recorded, not damped.
