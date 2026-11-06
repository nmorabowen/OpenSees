---
wp: LEGACY
title: "LadrunoArcLength -stabilize (33004): what viscous regularization can and cannot pass"
legacy_seq: 75
---
### `LadrunoArcLength -stabilize` (33004): what viscous regularization can and cannot pass
Measured on the live build (2026-06-16) while building the ADR-31 rung-4 seam. Four
non-obvious behaviours, all relevant to anyone wiring `-stabilize` into a driver:
- **Pure monotone softening is NOT passable by `-stabilize`.** In stabilize mode the
  integrator IS load control (mutually exclusive with the arc-length quadratic), and a
  softening branch has no equilibrium above the peak. Stabilized load control stalls at
  the strength peak *exactly* like plain load control (softening Concrete02 truss: both
  stop at `λ=30.0, ε=−0.002`). Softening is rung-3's job (switch to displacement control),
  NOT rung-4's. `-stabilize` only helps **snap-through-style limit points** where a
  continuing branch exists.
- **`-adaptStab` PREVENTS crossing a hard limit point.** It rescales `cVisc` each commit
  to hold `dissipVisc/Estrain0 ≈ fTarget`, which keeps the viscous force too weak to push
  through. On the von Mises snap-through truss, `-stabilize -adaptStab` (any `f`) stalls at
  the limit; only `-stabilize` *without* adaptStab and at an **elevated** `f` (≈1e-3…1e-2)
  crosses — and then by a slow **diffusive crawl** (ran to the 2000-step budget, cumulative
  ratio ~3e3), i.e. the "dynamic jump across the unstable branch" ADR-31 flags as R-LOG-MASK,
  not a traced path. Crossing is also **non-monotonic in `f`**: `f`=1e-3,1e-2 cross but
  `f`=5e-2,1e-1,5e-1 stall again (too much damping re-freezes the step).
- **`-cVisc N` is silently overwritten by the first-commit calibration.** `commit()` runs
  the `fTarget` calibration under `if (!cCalibrated)` with NO guard for a user-supplied
  `cVisc`, so `-cVisc 1`, `10`, `100` all produce the identical run (ratio 0.4495). The code
  comment "one-shot calibration … (skipped if -cVisc given)" is aspirational — the skip is
  not implemented. If you need an explicit coefficient today, the `-cVisc` path does not
  deliver it. (Out of scope to fix under ADR-31; flagged for a future ADR-20 follow-up.)
- **openseespy turns a command's stderr WARNING into a raised `OpenSeesError`.** The
  `scaleCVisc(factor)` guard `factor>0` writes a warning to `opserr` and returns −1, but in
  the Python module any `opserr` output during a command raises `opensees.OpenSeesError`
  rather than returning the −1. A driver must pre-validate `factor>0` itself (never rely on
  catching the −1) or wrap the call in try/except. (Tcl sees the −1 return; Python does not.)
  Verified by `torture_stabilize.py` + `test_robust_battery.py::test_stabilize_*` (4 cases).
- **RC-shell Phase 2a interlock — the MCFT `v_ci,max` formula is UNIT-DEPENDENT (SI: N, mm).**
  `LadrunoRCConcrete -interlock` bounds the crack-plane shear at the Vecchio–Collins limit
  `v_ci,max = 0.18·√fc' / (0.31 + 24·w/(a_g+16))` with crack width `w = eps_n·s_theta`. The
  numeric constants (0.18, 0.31, 24, 16) are empirical in **MPa and mm** — `√fc'` is `√(MPa)`,
  `w` and `a_g` are in **mm**. The rest of the kernel is unit-agnostic, but this one law is not:
  use it on an N–mm–tonne–s model (fc' in MPa = N/mm², lengths in mm) or rescale the constants.
  `s_theta` defaults to `-crackSpacing`, else `lch`, else 1.0; `a_g` (`-agg`) defaults to 16 mm.
  NB Phase 2a **CLIPS the smeared (damage-reduced) crack-plane shear `τ_sm=m_σ·sig_ip` to ±v_ci,max**
  (a bound), it does NOT substitute bare-elastic `G·γ`; below the cap the stress is unchanged.
- **Fixed-crack interlock only engages under NON-PROPORTIONAL loading.** The crack normal freezes to
  the principal-tensile direction at cracking; under a *proportional* path stress/strain stay coaxial,
  so the crack-plane shear is ~0 and `v_ci,max` never binds — interlock looks inert. It engages only
  once the principal direction ROTATES off the frozen normal (tension-then-shear, any non-radial path).
  TEST consequence: an off-axis interlock test MUST be two-stage (freeze oblique, then rotate shear onto
  it); a single proportional ramp gives `τ_nt≈0` and tests nothing.
- **rc_shell_ref.py oracle uses RAW `(x,y,q)` backbones; the C++ kernel ADJUSTS them
  (`buildBackbone` E-consistency).** So oracle-vs-C++ **absolute** stress only matches for
  quantities INDEPENDENT of the backbone q-adjustment. The Phase-1 β gate dodges this by being
  a RATIO; the Phase-2a interlock cap test dodges it by asserting only the **CAPPED** crack-plane shear
  `|τ_nt| == v_ci,max` (backbone-free); sub-cap `σ_xy` is the damaged smeared shear and differs
  C++-vs-oracle. To compare absolute normal stresses oracle-vs-C++ you must first port
  `buildBackbone`'s adjust() into the numpy Backbone. (Two-stage driver `_path_tension_then_shear`:
  tension on dof-1∝X, shear on dof-2∝X — disjoint DOFs let one SP per DOF realize the path under Penalty;
  pass `gamma0>0` to freeze an OBLIQUE crack for the off-axis rotation test.)
