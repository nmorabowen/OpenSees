# WP-134 — an independent SANISAND reference integrator (oracle for TIMs F18, WP-128/129)

Pure Python, no C++ change, no build. The code is in
`Ladruno_scripts/sanisand_reference/` (module, CLI, README). The validation runners and
every number below are in `Ladruno_files/testbed/sanisand_reference/` (scripts, with
outputs under `out/`). The tests are `tests/test_sanisand_reference.py` (pure) and
`tests/test_sanisand_reference_crosscheck.py` (skips without the WP-127 pyd).

## 0. Answers in one page

@@SUMMARY@@

## 1. Independence, and what was read

- **Written from the paper:** Dafalias & Manzari (2004), *Simple plasticity sand model
  accounting for fabric change effects*, J. Eng. Mech. 130(6):622–634.
- **Not read:** WP-128's `md_port.py`.
- **`ManzariDafalias.cpp`, read for two purposes only.** Each is recorded in §3:
  1. to list the U. Washington additions to the paper (`GetStateDependent`,
     `GetElasticModuli`, `GetF`, `GetLodeAngle`, `integrate()`'s α_in test,
     `explicit_integrator`'s branch logic);
  2. to find out why ModifiedEuler and RungeKutta45 disagree. Only the header of
     `ModifiedEuler()` was read for this, which is where U9 was found.
- **Not transcribed:** no line of the C++ integration schemes.

**Caveat on equation numbers.** This session could not open the paper's full text, so
the paper's own equation numbers are **not verified here**.
- §2 labels the equations E1–E14 and cites each one by the quantity the paper's
  multiaxial summary defines. These are the standard, widely reproduced DM04 forms.
- Mapping E1–E14 to the paper's numbering needs a check against the PDF (open item O1).

## 2. The equations (compression positive; tensors deviatoric where stated)

Notation:
- p = tr σ/3 and s = σ − pI.
- r = s/p.
- α is the back-stress ratio (deviatoric), z the fabric tensor (deviatoric).
- e is the void ratio.
- A:B = A_ij B_ij and ‖A‖ = √(A:A).

| # | quantity | form used |
|---|---|---|
| E1 | hypo-elasticity | dσ = 2G de + K dε_v I, with de = dev dε and dε_v = tr dε |
| E2 | shear modulus | G = G0 p_at (2.97 − e)²/(1 + e) (p/p_at)^½ |
| E3 | bulk modulus | K = 2(1+ν)/(3(1−2ν)) G |
| E4 | critical state line | e_c = e_c0 − λ_c (p/p_at)^ξ |
| E5 | state parameter | ψ = e − e_c |
| E6 | yield surface | f = ‖s − pα‖ − √(2/3) m p = 0 (m constant) |
| E7 | loading direction | n = (s − pα)/‖s − pα‖ = (r − α)/‖r − α‖ on f = 0 |
| E8 | Lode angle, interpolation | cos3θ = √6 tr n³; g(θ,c) = 2c/((1+c) − (1−c) cos3θ), c = M_e/M_c |
| E9 | bounding and dilatancy images | α^b_θ = √(2/3)[g M_c e^(−n_b ψ) − m] n; α^d_θ = √(2/3)[g M_c e^(n_d ψ) − m] n |
| E10 | hardening | dα = ⟨L⟩ (2/3) h (α^b_θ − α); h = b0/((α − α_in):n); b0 = G0 h0 (1 − c_h e)(p/p_at)^(−½) |
| E11 | plastic modulus | K_p = (2/3) p h (α^b_θ − α):n |
| E12 | flow | dε^p = ⟨L⟩ R, where R = B n − C(n² − I/3) + (D/3) I, B = 1 + (3/2)((1−c)/c) g cos3θ, C = 3√(3/2)((1−c)/c) g |
| E13 | dilatancy, fabric | D = A_d (α^d_θ − α):n, A_d = A0(1 + ⟨z:n⟩); dz = −c_z ⟨−dε^p_v⟩ (z_max n + z) |
| E14 | loading index, strain-driven | L = N/H, N = 2G n:de − K (n:r) dε_v, H = K_p + 2G(B − C tr n³) − K D (n:r) |

Derived facts that the code relies on:
- (a) ∂f/∂σ = Q = n − (1/3)(n:r) I. So consistency, df = Q:dσ − p n:dα = 0, gives
  E14 exactly.
- (b) **B − C tr n³ ≡ 1** for every unit deviatoric n. The two terms cancel through
  cos3θ = √6 tr n³ (pinned by a test), so H = K_p + 2G − K D (n:r).
- (c) With a = (α − α_in):n and H_s = a·H = (2/3) p b0 (b:n) + a (2G − K D n:r), two
  expressions stay finite at a = 0, where h = ∞:
  - L = a N/H_s
  - h L = b0 N/H_s

  This is the exact continuous extension. The reference uses it and needs no cap.
- Void ratio: de = −(1 + e) dε_v (small-strain kinematics).

### Ambiguities in the paper and how they were resolved

| # | ambiguity | resolution here |
|---|---|---|
| A1 | Sign convention, effective stress | compression positive, effective stress; tensors as the fork stores them (§4) |
| A2 | Whether m evolves | constant, as in DM04's final model |
| A3 | Hypo- vs hyper-elasticity | hypo-elastic (E1), moduli at the current state (U8/U9 are the UW variants) |
| A4 | cos3θ from n or r; round-off beyond ±1 | from n (the flow and bounding images live on n); clamped to [−1, 1] |
| A5 | When α_in is updated. The paper resets it at the "initiation of a new loading process", identified by (α − α_in):n < 0 | at every **plastic onset** (the on-surface mode decision) with a < 0, α_in := α; and inside a plastic segment an **event** on a → 0⁻ re-seats it. Result: h ≥ 0 everywhere on the paper branch (§6.4: 0 negative-h samples in 960 ring runs) |
| A6 | Loading/unloading in strain-driven form; H ≤ 0 | plastic iff N > 0 and H > 0; elastic iff N ≤ 0. For **N > 0 with H ≤ 0 the continuous problem has no admissible rate solution**: the elastic branch has df = N > 0 and the plastic branch has L < 0. The integrator **stops** with `H_nonpositive`. It never treats this case as elastic |
| A7 | h at a = 0 | the exact regularised form of fact (c); no cap (UW's 1e10 cap is U7) |
| A8 | Void-ratio evolution | de = −(1 + e) dε_v (UW: −(1 + e_init), U5) |
| A9 | "Inside the bounding surface" for a multiaxial α | the bounding surface in α-space is {√(2/3) α^b(θ(u)) u : u unit deviatoric}, so α is inside iff **ρ_α** = ‖α‖/(√(2/3) α^b(θ_α)) < 1, with θ_α the Lode angle of α itself. WP-128's ρ_b uses θ of n; both are reported. ρ_b can exceed 1 while α is inside, when n and α point differently (§6.4) |
| A10 | e in G | current e (UW: e_init, U4) |
| A11 | p_at | 100 kPa for the Toyoura set; 101 kPa for the campaign set (as run) |
| A12 | ⟨−dε^p_v⟩ in E13 | = ⟨−L D⟩ with L > 0 on the plastic branch |
| A13 | p → 0 | the model is singular: G → 0 and b0 → ∞. Integration stops at `p_floor` (1e-6 kPa default) with status `p_floor` |
| A14 | A start state with f > 0 | undefined in the paper. Refused with `start_outside_yield`; `start_outside="plastic"` is an explicit opt-in |
| A15 | A start state with α outside the bounding surface (ρ_α > 1) | defined by the equations (K_p < 0 pulls α back). Integrated and reported, not refused (§6.5) |

## 3. UW additions (each one a separate switch, paper by default)

| id | UW addition / convention | where (read-only) | switch |
|---|---|---|---|
| U1 | low-p dilatancy sigmoid, D ×= 1/(1 + exp(7.6349 − 7.2713 p)) for p < 0.05 P_atm | `GetStateDependent` | `d_factor` |
| U2 | p → p + p_r in f, n, r, ψ, b0, D, K_p (not in G) | `GetF`, `GetStateDependent`, ... | `p_residual` |
| U3 | G, K use max(p, p_min) (fork `-Pmin`) | `GetElasticModuli` | `p_min` |
| U4 | G with **e_init**, not the current e (the fork's `mUseCurrentVoidRatioInG` seam, false by default) | `GetElasticModuli` | `g_void_ratio="initial"` |
| U5 | e = e_init − (1 + e_init) tr ε, i.e. de = −(1 + e_init) dε_v | `explicit_integrator` | `void_ratio_law="initial"` |
| U6 | α_in := α_n **once per increment**, when (α_n − α_in,n):(Ce:Δε) < 0. There is no reseat inside the increment, so a plastic onset or a zero crossing with (α − α_in):n < 0 runs with **h < 0** (WP-128's **mechanism G**) | `integrate()` | `alpha_in_rule="uw"` |
| U7 | h = 1e10 when \|(α − α_in):n\| < 1e-10, else b0/a (negative allowed) | `GetStateDependent` | `h_cap=1e10` |
| U8 | the elastic predictor and intersection use Ce of the committed state | `explicit_integrator` | `elastic_moduli="frozen"` |
| U9 | **ModifiedEuler never recomputes K, G inside the increment**: all substeps use the committed moduli (`ModifiedEuler()` header: `aC = GetStiffness(K, G)`, no `GetElasticModuli` call). RungeKutta45 does recompute them | `ModifiedEuler` | `elastic_moduli="frozen_increment"` |
| U10 | the loading test at a start on the surface uses n:Δσ_trial/‖Δσ_trial‖ > −√TolF, i.e. n, **not** Q = ∂f/∂σ. It ignores −(n:r)dp, so a path that unloads f (an isotropic compression that lowers η) is taken as "plastic" | `explicit_integrator` | not a model option. Discrete; attributed, not reproduced |
| — | discrete ModifiedEuler behaviour: stress-only error (E); Λ < 0 → elastic with α re-derived from r (F); no cap on step growth; forced accept at dT_min with an M_c clamp | `ModifiedEuler` | not reproduced. The reference has no discretisation; these are what it measures |

Presets:
- `Options()`: the paper.
- `Options.uw()`: U1–U8. The comparator for the C++ RungeKutta45.
- `Options.uw_me()`: U1–U9. The comparator for the C++ ModifiedEuler.
- `ring_variants()["uw_model"]`: the UW **constitutive** additions U1–U5 with the
  paper's α_in rule and continuous moduli. **This is the oracle a corrected C++
  integrator should reproduce**, since WP-129 fixes the integration, not the model.
- `uw_rule`: `uw_model` + U6 + U7.

## 4. Conventions (to exchange states with the C++ unchanged)

- Compression positive: the internal `mSigma`, and the ring CSVs despite their README
  (WP-127 finding A).
- σ, α, α_in and z are 6-vectors of **tensor** components in the order xx yy zz xy yz zx.
- Strain is Voigt with **engineering** shear, as in `ladrunoSANISANDReplay -dStrain`.
- Units are kPa.
- The C++ side runs through `ladrunoSANISANDReplay -convention compressionPositive
  -type 3D`, one call per increment, `-primed 1 -prevIncrNorm 0`. This happens in a
  CPython 3.12 subprocess (`_cxx_runner.py`, `-S`, manual paths, asserting
  `opensees.__file__`) on the WP-127 build in
  `.claude/worktrees/tims-implementation-review-3733c6/dist/bin`, used read-only.
- The replay reports the *committed* e, so the C++ end e is taken as
  e − (1 + e_init) Δε_v.

## 5. Integrator

- The state is (σ, α, z, e, Σ⟨L⟩, ε), 26 unknowns, integrated over pseudo-time
  t ∈ [0, 1] of one prescribed increment.
- Mixed control: each Voigt component is strain-controlled or stress-rate-controlled. The
  unknown strain rates are solved from the current (elastic or elastoplastic) Mandel
  tangent at every right-hand-side call.
- Solver: SciPy `solve_ivp`, Radau IIA (BDF selectable), rtol 1e-10 by default, atol
  scaled per block (σ by p0, α, z by z_max, e). The finite-difference Jacobian uses a
  fixed relative step. SciPy's own `num_jac` overflows on the columns the rates do not
  depend on.

Events end a smooth segment. The mode is decided again at each one:

| mode | event | trigger |
|---|---|---|
| elastic | yield | f − ftol ↑ (ftol = 1e-8 × the cone radius) |
| plastic | unload | N ↓ |
| plastic | H sign | H_s through 0 |
| plastic | α_in reversal | a ↓ 0 (paper rule: reseat) |
| plastic | kinks | ⟨z:n⟩ and ⟨−LD⟩ |
| any | p_floor | p reaches p_floor |

Back-to-back zero-length kink segments drop the kink restarts for one segment, to stop
chatter.

Reported per increment:
- status and t_end;
- the end state;
- f at exit, and max |f| over the plastic segments (the consistency drift);
- max ρ_b and max ρ_α along the path, sampled at every accepted solver step;
- the minimum (α − α_in):n over the plastic samples (the h < 0 detector);
- the segments with their events, and the reseats.

@@RESULTS@@
