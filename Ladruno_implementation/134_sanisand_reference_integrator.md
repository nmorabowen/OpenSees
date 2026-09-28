# WP-134 — an independent SANISAND reference integrator (oracle for TIMs F18, WP-128/129)

Pure Python, no C++ change, no build. The code is in
`Ladruno_scripts/sanisand_reference/` (module, CLI, README). The validation runners and
every number below are in `Ladruno_files/testbed/sanisand_reference/` (scripts, with
outputs under `out/`). The tests are `tests/test_sanisand_reference.py` (pure) and
`tests/test_sanisand_reference_crosscheck.py` (skips without the WP-127 pyd).

## 0. Answers in one page

1. **The reference reproduces DM04 (§6.1).**
   - Elastic isotropic compression matches the closed form to 1.7e-15.
   - On the Toyoura set, the undrained family at e = 0.833 behaves as DM04 shows:
     - dilative at p0 = 100 kPa;
     - phase transformation at 1000 and 2000 kPa (p dips to 659 and 937 kPa, then
       recovers);
     - contractive at 3000 kPa, with a q peak of 1409 and quasi-steady softening.

     All four p0 end on the **same** critical state, p 1087 kPa and q 1359 kPa, with
     η/M − 1 ≤ 3e-4 and |ψ| ≤ 1e-4.
   - Drained, dense sand dilates (ε_v −8.4 %) with a peak η above M; loose sand
     contracts to the CSL (|ψ| 2e-3).
   - Consistency drift is ≤ 3e-7 kPa.
2. **Against the fork's C++ (§6.2), on benign states.** These are on-surface, α
   admissible, at p0 = 20/50/100 kPa, with 4 probe directions × δ 1e-5/1e-4/1e-3.
   - On monotonic plastic increments (class A), the ModifiedEuler-faithful continuum
     `Options.uw_me()` equals the C++ ModifiedEuler at TolR 1e-8 (`-honorTolR 1`): the
     stress increment to a **median 3e-8 relative (max 1.2e-6)**, α to 3e-7.
   - On increments that first cross the m-cone (class B), the agreement is the same at
     δ = 1e-5 (max 7e-7).
   - Every difference between the **paper** and the C++ is attributed:
     - **U4**: G with e_init, 2–4 %.
     - **U9**: ModifiedEuler never updates K, G inside the increment. This costs 0.6 %
       / 6 % / 24 % at δ = 1e-5 / 1e-4 / 1e-3, and **does not shrink with TolR**.
       This is a new finding.
   - **RK45 is not a usable C++ reference**, as WP-128 said: its dT_min is hard-coded
     at 1e-3. Its median disagreement is 1–45 %.
3. **The remaining C++ differences are one discrete mechanism, not tolerance (§6.3).**
   - Every larger disagreement comes from an increment the C++ ModifiedEuler took in
     **≤ 7 substeps with err = 0 exactly**, even at TolE 1e-8. The C++ trace shows
     accepted substeps with `err = 0.0`, and the next one swallowing the rest of the
     increment. These are classes C and D, and the class-B outliers at δ ≥ 1e-4.
   - The mechanism: the Λ < 0 branch is taken as elastic, so both Heun stages agree,
     the error is 0 and the step factor is unbounded. This is **mechanism F**, entered
     through U10: the start-on-surface loading test uses n:Δσ, not ∂f/∂σ:Δσ.
   - Size of the damage: 20–65 % errors in the stress increment on **admissible,
     benign states at 20–100 kPa**.
4. **The WP-128 smallest reproducer (§6.4).** Start at σ = 0.0101 I, α = α_in = z = 0,
   and apply one dε_yy = 1e-4.
   - The continuum gives **η = 0.534, ρ_b = ρ_α = 0.251**. The paper model, the UW
     constitution and the UW α_in rule all agree, and ρ_α ≤ 0.29 on the whole sweep.
   - C++ ME (campaign) gives η 10.83 and ρ 5.10 in one substep. ME at TolE 1e-8 gives
     0.27 in 2394 substeps, off by U9 alone.
   - The escape (ρ_α > 1) sets in between Δp/p = 4.8 (0.81) and 7.3 (1.05) at
     p_s = 1 kPa, and between 4.1 (0.75) and 8.7 (1.19) at p_s = 5 kPa. That is, an
     increment that multiplies p by about 6–8.
5. **Ring replay (§6.5), 80 rows × {isoComp, shear} × δ {1e-5, 1e-4, 1e-3}.**
   - For the oracle (`uw_model`: UW constitution, paper integration), on the 78
     admissible rows:
     - **0 escapes**: max ρ_α 0.983, which is the largest start;
     - **0 samples with h < 0**;
     - |f| at exit ≤ 1.7e-7 kPa.
   - The C++, from the same admissible starts:
     - ModifiedEuler (campaign) escapes **25 times, up to ρ_α 18.5** (shear, 1e-3);
     - ModifiedEuler at TolE 1e-8 escapes **8 times, up to 4.7**.
     - Every ME escape, and 7 of 8 ME8 escapes, sits in the ≤ 7-substep err = 0
       bucket.
6. **The two inadmissible rows (b8 1950 gp 2/3, ρ_α 6.8/7.3, f ≈ 0) (§6.6).** No
   integrator of the model can "take" them in the sense of returning an admissible
   state:
   - Under isoComp, and under shear for gp 2, the reference integrates them honestly.
     α is pulled back by K_p < 0 but stays far outside: ρ_α 4.0–6.7 at the end of the
     increment, f ≈ 0.
   - Under shear, **gp 3's rate equations are singular**. After the in-plastic α_in
     reseat, (α−α_in):n → 0 and b:n → 0 together, so the loading index is 0/0
     (H_s ≈ 2e-4). The reference **stops** at t = 0.054 and says so.
   - The C++ returns `rc 0` for all of them in **one substep**:
     - under isoComp at 1e-3 it teleports to ρ_α 0.29 (the M_c clamp);
     - under shear 1e-5 it commits f = +1.4e-2 kPa.
7. **Mechanisms (§7).**
   - **G is confirmed as a trigger at the model level.** The UW once-per-increment α_in
     rule, integrated exactly (`uw_rule`), puts h < 0 in 146 of 480 ring runs. It
     makes the rate problem inadmissible (`H_nonpositive`) in 27 runs and drives α
     out of the bounding surface from 25 admissible starts (up to ρ_α 1.97). With
     h < 0 the α equation is a *repelling* relaxation: all 43 "ok" `uw_rule` runs
     with |f| > 1e-5 at exit (up to 287 kPa) have h < 0.
   - The paper's event-driven reseat makes G impossible by construction: 0 of 960
     runs.
   - **E is confirmed as the enabler.** The escape needs a stress-only error.
   - **F is more than a compounder.** It produces err = 0 exactly, so no tolerance
     sees it. It alone accounts for the 20–65 % benign-state errors at TolE 1e-8, and
     for every ring escape of the campaign ME.
   - In the smallest reproducer G does not act in the continuum (the UW rule and the
     paper rule give identical answers). There, the trigger is the *discrete*
     first-stage overshoot at h ≈ 1e10, which manufactures the (α−α_in):n < 0 of
     stage 2. That is G created by the discretisation, then passed by E and F.
   - **WP-129 must fix E, F, U9 and the α_in rule together.**

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
| A5 | When α_in is updated. The paper resets it at the "initiation of a new loading process", identified by (α − α_in):n < 0 | at every **plastic onset** (the on-surface mode decision) with a < 0, α_in := α; and inside a plastic segment an **event** on a → 0⁻ re-seats it. Result: h ≥ 0 everywhere on the paper branch (§6.5: 0 negative-h samples in 960 ring runs) |
| A6 | Loading/unloading in strain-driven form; H ≤ 0 | plastic iff N > 0 and H > 0; elastic iff N ≤ 0. For **N > 0 with H ≤ 0 the continuous problem has no admissible rate solution**: the elastic branch has df = N > 0 and the plastic branch has L < 0. The integrator **stops** with `H_nonpositive`. It never treats this case as elastic |
| A7 | h at a = 0 | the exact regularised form of fact (c); no cap (UW's 1e10 cap is U7) |
| A8 | Void-ratio evolution | de = −(1 + e) dε_v (UW: −(1 + e_init), U5) |
| A9 | "Inside the bounding surface" for a multiaxial α | the bounding surface in α-space is {√(2/3) α^b(θ(u)) u : u unit deviatoric}, so α is inside iff **ρ_α** = ‖α‖/(√(2/3) α^b(θ_α)) < 1, with θ_α the Lode angle of α itself. WP-128's ρ_b uses θ of n; both are reported. ρ_b can exceed 1 while α is inside, when n and α point differently |
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

## 6. Validation

### 6.1 Elastic sanity, CSL limit, DM04 behaviour (`run_paper.py` → `out/paper.md`)

Elastic isotropic compression, campaign set: p0 = 20 kPa, ε_v = 3e-3, with G at e_init
so that K = k√p.
- Reference: p = 1085.906892 kPa.
- Closed form (√p0 + k ε_v/2)²: 1085.906892 kPa.
- Relative error 1.7e-15. One elastic segment; α untouched.

Toyoura (DM04 set: G0 125, ν 0.05, M 1.25, c 0.712, λ_c 0.019, e_c0 0.934, ξ 0.7,
m 0.01, h0 7.05, c_h 0.968, n_b 1.1, A0 0.704, n_d 3.5, z_max 4, c_z 600, p_at 100).
Triaxial compression to 30 % axial strain, rtol 1e-8:

| test | ψ0 | q_peak | p_min | p_end | q_end | η_end/M | ψ_end | ε_v end | max ρ_α | max \|f\| plastic |
|---|---|---|---|---|---|---|---|---|---|---|
| dense drained (e 0.735, 100 kPa) | −0.180 | 289.8 | 100.0 | 174.3 | 223.0 | 1.0233 | −0.0188 | −0.0841 | 1.008 | 6.0e-08 |
| loose drained (e 0.95, 100 kPa) | +0.035 | 213.3 | 100.0 | 171.1 | 213.3 | 0.9973 | +0.0021 | +0.0215 | 1.000 | 4.9e-08 |
| e 0.833 undrained, p0 100 | −0.082 | 1358.5 | 85.7 | 1086.5 | 1358.5 | 1.0003 | −0.0001 | 0 | 1.013 | 5.1e-08 |
| e 0.833 undrained, p0 1000 | −0.006 | 1359.0 | 659.4 | 1087.0 | 1359.0 | 1.0002 | −0.0000 | 0 | 1.004 | 1.4e-07 |
| e 0.833 undrained, p0 2000 | +0.054 | 1359.1 | 937.4 | 1087.1 | 1359.1 | 1.0001 | −0.0000 | 0 | 1.001 | 2.0e-07 |
| e 0.833 undrained, p0 3000 | +0.104 | 1408.6 | 1026.1 | 1087.2 | 1359.1 | 1.0001 | −0.0000 | 0 | 1.000 | 2.6e-07 |
| e 0.907 undrained, p0 100 | −0.008 | 201.3 | 48.1 | 161.0 | 201.3 | 1.0007 | −0.0005 | 0 | 1.001 | 3.0e-08 |
| e 0.735 undrained, p0 100 | −0.180 | 3608.1 | 96.9 | 2867.3 | 3583.7 | 0.9999 | +0.0001 | 0 | 1.040 | 4.7e-08 |

Reading the table:
- Every undrained test reaches the CSL.
- The dense drained test is still approaching it at 30 % (ψ −0.019, η/M 1.023),
  post-peak.
- These checks are qualitative. No digitised DM04 curves were available, so the paper's
  figures are matched in kind, not point by point (open item O2).
- **ρ_α slightly above 1 (up to 1.04) is the model, not an error.** In DM04, softening
  happens *because* α^b(ψ) contracts inside α as ψ rises. Then K_p < 0 and α is pulled
  back, so ρ_α ≈ 1 + small on softening branches is admissible behaviour. The "6×" of
  the ring is not.

### 6.2 Cross-check against the C++ on benign states (`run_crosscheck.py`, `analyse_crosscheck.py` → `out/crosscheck*.md`)

States:
- campaign set;
- on the yield surface, e 0.72, α_in = z = 0;
- p0 20 / 50 / 100 kPa;
- directions: TC at η = 0.5 M_c, TC + xy shear at 0.8 M_c, TE at 0.4 c M_c.

Probes: txLoad (+δ xx), txUnload (−δ xx), shear (γ_xy = δ), isoComp (δ, δ, 0), at
δ = 1e-5 / 1e-4 / 1e-3.

C++ runs:
- ME8: ModifiedEuler at TolR 1e-8, `-honorTolR 1`;
- ME: the campaign setting;
- RK45: RungeKutta45 at 1e-10.

Metric: ‖Δσ_ref − Δσ_C++‖/‖Δσ_C++‖, on the increment.

Classes, from the reference `uw_me` path:
- A: monotonic plastic.
- B: elastic transit through the m-cone, then plastic, with h ≥ 0.
- C: B, but the UW α_in rule leaves h < 0 at the onset (mechanism G).
- D: the reference stops. Either `H_nonpositive` under the UW rule, or `p_floor`: the
  increment empties p; the C++ ME8 then hits its 200 000-substep cap.

| δ | class | n | uw_me vs ME8 Δσ (median / max) | Δα (median / max) | uw vs RK45 Δσ (median / max) | paper vs ME8 Δσ (median) | ME8 substeps (median) |
|---|---|---|---|---|---|---|---|
| 1e-5 | A | 12 | 3.6e-8 / 1.2e-6 | 7.2e-7 / 4.5e-6 | 8.5e-3 / 0.80 | 3.7e-2 | 204 |
| 1e-5 | B | 15 | 6.3e-9 / 6.7e-7 | 2.6e-7 / 1.4e-6 | 6.9e-3 / 0.65 | 2.9e-2 | 612 |
| 1e-5 | C | 7 | 0.29 / 0.44 | 1.1 / 1.3 | 0.29 / 0.44 | 4.3e-2 | **1** |
| 1e-5 | D | 2 | 0.52 / 0.63 | 0.92 / 0.94 | 0.22 / 0.41 | 4.2e-2 | **1** |
| 1e-4 | A | 12 | 2.3e-8 / 6.1e-8 | 3.1e-7 / 4.9e-6 | 4.0e-2 / 0.90 | 7.2e-2 | 798 |
| 1e-4 | B | 15 | 4.0e-8 / 0.23 | 1.2e-7 / 0.47 | 3.5e-2 / 0.73 | 8.1e-2 | 1632 |
| 1e-4 | C | 9 | 0.36 / 0.52 | 1.2 / 1.5 | 0.48 / 0.97 | 0.13 | **1** |
| 1e-3 | A | 6 | 4.1e-8 / 9.9e-8 | 2.1e-6 / 7.5e-6 | 0.45 / 0.91 | 0.20 | 1419 |
| 1e-3 | B | 12 | 1.0e-7 / 0.29 | 4.7e-7 / 0.54 | 0.36 / 0.83 | 0.22 | 3974 |
| 1e-3 | C | 4 | 0.53 / 0.65 | 1.8 / 2.2 | 0.70 / 0.94 | 1.4 | **1** |
| 1e-3 | D | 14 | 3.2e-4 / 1.0 | 1.1 / 1.7 | 2.3e-4 / 1.0 | 6.6e-3 | 7706 |

Effect of each UW addition alone on class A: the median relative change of Δσ when
that one addition is switched back to the paper inside `uw_me`.

| addition | δ 1e-5 | δ 1e-4 | δ 1e-3 |
|---|---|---|---|
| U1 D_factor (inactive at p ≥ 20) | 0 | 0 | 0 |
| U3 p_min (inactive) | 0 | 0 | 0 |
| **U4 G(e_init)** | 3.7e-2 | 3.7e-2 | 2.1e-2 |
| U5 e-law(e_init) | 2e-8 | 3e-7 | 6e-7 |
| U6 α_in rule (no reversal in class A) | 0 | 0 | 0 |
| U7 h cap | 2e-14 | 8e-13 | 1e-13 |
| **U9 frozen K, G over the increment** | 5.8e-3 | 5.8e-2 | 0.24 |

**Agreement statement.** On monotonic increments and on δ = 1e-5 cone transits, the
reference with U1–U9 reproduces the C++ ModifiedEuler at TolR 1e-8 to ≤ 1.2e-6 in σ
and ≤ 7.5e-6 in α. Every paper-vs-C++ difference in those classes is U4 + U9.

The class-B outliers at δ ≥ 1e-4 (p20/50/100 TE txLoad, and TE isoComp at 1e-3) and
all of classes C/D are the discrete path of §6.3.

**U9 is a finding for WP-129.** ModifiedEuler converges, as TolR → 0, to a *different*
ODE: hypo-elasticity with the committed moduli held for the whole increment. Its error
is first order in the increment. It is 24 % of the stress increment at δ = 1e-3 and
p = 20–100 kPa, and it is invisible to the error test.

### 6.3 The discrete path behind the remaining differences (`probe_cpp_trace.py` → `out/cpp_trace_probe.json`)

The C++ ModifiedEuler per-substep trace at TolE 1e-8 (`-honorTolR 1`):

| case | path | substeps | trace (T, dT, err, code) |
|---|---|---|---|
| p20 TE txLoad 1e-4 | 4 (unload-then-plastic) | 5 | (0, 1, 3.5e-8, rej), (0, 0.43, 1.5e-8, rej), (0, 0.28, 9.8e-9, acc), (0.28, 0.23, **0.0**, acc), (0.51, 0.49, **0.0**, acc) |
| p20 TC isoComp 1e-5 | 3 ("plastic" per U10) | 1 | (0, 1, **0.0**, acc) |
| p100 TC+shear txLoad 1e-5 | 3 | 1 | (0, 1, **0.0**, acc) |

`err = 0.0` exactly means both Heun stages returned the same stress increment. That is
the Λ < 0 → "elastic" branch in both stages (mechanism F). With the step factor
uncapped, the next substep takes the whole remainder, so no tolerance can see it.
These increments differ from the continuum by 20–65 % in Δσ and 100–200 % in Δα
(classes C/D above).

On the isoComp probes, U10 makes it worse. The start-on-surface test calls the
increment "plastic" (n:Δσ > 0), while ∂f/∂σ:Δσ < 0: isotropic compression lowers η, so
the path is genuinely elastic at first.

### 6.4 WP-128's smallest reproducer (`run_reproducer.py` → `out/reproducer.md`)

Start σ = p_s I, α = α_in = z = 0, e = e_init. Apply one plane-strain dε_yy = +δ.

| p_s | δ | Δp/p (ref) | reference `uw_model`: η / ρ_b / ρ_α | paper η / ρ_b | `uw_rule` η / ρ_b | C++ ME (campaign): substeps, η / ρ_b / ρ_α / f | C++ ME8: substeps, η / ρ_b |
|---|---|---|---|---|---|---|---|
| 0.0101 | 1e-5 | 2.78 | 0.317 / 0.147 / 0.147 | 0.318 / 0.147 | 0.317 / 0.147 | 1, 1.19 / 0.79 / 0.56 / 1.6e-12 | 463, 0.32 / 0.15 |
| 0.0101 | 3e-5 | 13.7 | 0.433 / 0.202 / 0.202 | 0.434 / 0.202 | 0.433 / 0.202 | 1, 3.21 / 2.13 / 1.51 / 1.7e-9 | 1329, 0.44 / 0.20 |
| **0.0101** | **1e-4** | 108 | **0.534 / 0.251 / 0.251** | 0.535 / 0.251 | 0.534 / 0.251 | **1, 10.83 / 5.10 / 5.10** / 2.1e-8 | 2394, 0.57 / 0.27 |
| 0.0101 | 3e-4 | 859 | 0.601 / 0.287 / 0.287 | 0.601 / 0.288 | 0.601 / 0.287 | 1, 0.84 / 16.25 / 11.52 / **11** | 4277, 0.67 / 0.31 |
| 0.1 | 1e-4 | 15.0 | 0.438 / 0.206 / 0.206 | 0.439 / 0.206 | 0.438 / 0.206 | 1, 3.39 / 1.60 / 1.60 / 8e-8 | 3840, 0.44 / 0.21 |
| 1 | 1e-4 | 2.78 | 0.317 / 0.149 / 0.149 | 0.317 / 0.149 | 0.317 / 0.149 | 1, 1.20 / 0.57 / 0.57 | 3022, 0.31 / 0.15 |
| 1 | 1.5e-4 | 4.84 | 0.361 / 0.171 / 0.171 | 0.361 / 0.171 | 0.361 / 0.171 | 1, 1.70 / 0.81 / 0.81 | 3569, 0.36 / 0.17 |
| 1 | 2e-4 | 7.34 | 0.392 / 0.187 / 0.187 | 0.392 / 0.187 | 0.392 / 0.187 | 1, 2.20 / **1.05 / 1.05** | 3997, 0.39 / 0.18 |
| 1 | 2.5e-4 | 10.3 | 0.414 / 0.198 / 0.198 | 0.414 / 0.198 | 0.414 / 0.198 | 1, 2.67 / 1.81 / 1.28 | 4346, 0.42 / 0.20 |
| 1 | 3e-4 | 13.7 | 0.432 / 0.208 / 0.208 | 0.432 / 0.208 | 0.432 / 0.208 | 1, 3.15 / 2.14 / 1.51 | 4639, 0.44 / 0.21 |
| 5 | 3e-4 | 4.13 | 0.347 / 0.168 / 0.168 | 0.347 / 0.168 | 0.347 / 0.168 | 1, 1.54 / 0.75 / 0.75 | 3411, 0.34 / 0.16 |
| 5 | 5e-4 | 8.68 | 0.400 / 0.197 / 0.197 | 0.400 / 0.197 | 0.400 / 0.197 | 1, 2.42 / **1.19 / 1.19** | 4164, 0.40 / 0.19 |

**The "true" answer** (dε_yy = 1e-4) is **η = 0.534, ρ = 0.251**, and ρ stays below 0.29
on every case. C++ ME at TolE 1e-8 lands within 0.02 of it, which is U9.

**Escape onset** (the campaign ME's ρ_α crossing 1): between Δp/p = 4.8 and 7.3 at
p_s = 1 kPa, and between 4.1 and 8.7 at p_s = 5 kPa. That is one increment that
multiplies p by about 6–8, not about 4. WP-128's 3.3–5.2 used a different Δp/p (from
the C++ end state, not the continuum).

The UW α_in rule changes nothing here: α_in = α = 0 and (α−α_in):n ≥ 0 throughout the
continuous path. **In the reproducer, G appears only in the discrete scheme.** The
first Heun stage at h = 1e10 overshoots α, and the second stage then sees
(α−α_in):n < 0.

### 6.5 Ring replay (`run_ring.py`, `analyse_ring.py` → `out/ring_rows.md` (per row), `out/ring_summary.md`, `out/ring_analysis.md`)

Setup:
- 80 rows (b8 40, b16 40), read compression-positive from the WP-127 attachments.
- Campaign set, ν 0.312885, p_min 0.0101, p_r 0.
- The WP-127 probes (isoComp (δ, δ, 0), shear γ_xy = δ) at δ = 1e-5 / 1e-4 / 1e-3.
- Admissible = ρ_α ≤ 1 at the start: 78 rows. b8 1950/2 and 1950/3 are not.
- `ring_rows.md` carries per row: status, f at exit, max ρ_α, p_end, η_end, and the C++
  ME and ME8 results with ‖Δσ_C++ − Δσ_ref‖.

| δ | probe | oracle `uw_model` status | oracle max ρ_α (adm) / escapes | `uw_rule` status | `uw_rule` escapes | C++ ME escapes / max ρ_α | ME rel. diff median / p95 | C++ ME8 escapes / max ρ_α | ME8 rel. diff median / p95 |
|---|---|---|---|---|---|---|---|---|---|
| 1e-5 | isoComp | 80 ok | 0.983 / 0 | 73 ok, 7 H≤0 | 2 | 0 / 0.83 | 0.15 / 0.34 | 0 / 0.83 | 0.15 / 0.34 |
| 1e-5 | shear | 79 ok, 1 stop (1950/3) | 0.983 / 0 | 77 ok, 2 H≤0, 1 stop | 0 | 0 / 0.93 | 4.9e-3 / 1.1 | 0 / 0.94 | 1.3e-4 / 0.66 |
| 1e-4 | isoComp | 80 ok | 0.983 / 0 | 62 ok, 7 H≤0, 11 stop | 9 | 0 / 0.82 | 0.38 / 0.67 | 0 / 0.73 | 0.37 / 0.53 |
| 1e-4 | shear | 79 ok, 1 stop | 0.983 / 0 | 77 ok, 2 H≤0, 1 stop | 0 | 1 / 1.80 | 4.4e-3 / 1.5 | 0 / 0.89 | 2.3e-3 / 1.1 |
| 1e-3 | isoComp | 80 ok | 0.983 / 0 | 60 ok, 7 H≤0, 13 stop | 14 | 6 / 7.35 | 0.82 / 0.96 | 0 / 0.71 | 0.82 / 0.92 |
| 1e-3 | shear | 79 ok, 1 stop | 0.983 / 0 | 77 ok, 2 H≤0, 1 stop | 0 | **18 / 18.5** | 2.3e-2 / 5.8 | **8 / 4.7** | 1.3e-2 / 4.6 |

Over all 480 runs per variant:

- **Paper model and `uw_model` (the oracle):**
  - 0 runs with h < 0 anywhere;
  - 0 escapes from admissible starts;
  - |f| at exit ≤ 1.7e-7 kPa;
  - no elastic sample above the surface.

  The only non-ok statuses are 1950/3 under shear (§6.6). The paper model has one extra
  `p_floor`.
- **`uw_rule`** (the UW α_in rule, integrated exactly):
  - **146 of 480 runs have h < 0**;
  - 27 runs are `H_nonpositive`: there is no admissible rate solution;
  - 27 runs stop as `solver_failed`: stiff and unstable;
  - **25 admissible starts escape** (α itself moves out, up to ρ_α 1.97);
  - 43 "ok" runs exit with |f| > 1e-5 kPa (up to 287 kPa), **all of them with h < 0**.
    With h < 0, dα ∝ h(α^b − α) *repels* α from the bounding surface. The ODE is
    unstable and consistency drifts even at rtol 1e-10.
  - The 333 ok runs with h ≥ 0 keep |f| ≤ 1.6e-7.
- **The C++ versus the oracle, split by the C++ substep count:**

| C++ | ≤ 7 substeps (err = 0 path): n / median / p95 | > 7 substeps: n / median / p95 | escapes (≤ 7 / > 7) |
|---|---|---|---|
| ME (campaign) | 248 / 0.43 / 4.1 | 220 / 6.9e-3 / 0.83 | **25 / 0** |
| ME8 (TolE 1e-8) | 192 / 0.39 / 1.3 | 276 / 1.1e-2 / 0.83 | **7 / 1** |

Where the C++ actually substeps, it sits within about 1 % of the oracle; U9 explains
that level (§6.2). The large differences and every ME escape come from the err = 0
path.

### 6.6 The two inadmissible rows, said plainly

b8 element 1950, Gauss points 2 and 3:
- p' 0.35 kPa, η 12.1 / 12.9;
- ρ_α = 6.81 / 7.28 (ρ_b 5.86 / 6.26);
- f ≈ 0 (7.8e-11 / 7.5e-10): the stress rides the m-cone around a bad α.

| row | probe | δ | reference (`uw_model`) | C++ ME (campaign) |
|---|---|---|---|---|
| 1950/3 | isoComp | 1e-5 / 1e-4 / 1e-3 | ok. α pulled back by K_p < 0 but still outside: ρ_α end 6.17 / 5.40 / 4.35, η 10.9 / 9.6 / 7.8, f ≤ 8e-9 | rc 0, 1 substep: ρ_α 4.40 / **0.92** / **0.29** (η 7.86 / 1.79 / 0.47). The Mc-clamp "teleport" |
| 1950/3 | shear | all | **stops** at t = 0.054 / 0.0054 / 0.00054 (the same strain, 5.4e-7). After an in-plastic α_in reseat, (α−α_in):n → 1e-10 and b:n → 4e-7 together; the loading index is 0/0 (H_s ≈ 2e-4); rate equations singular | rc 0, 1 substep: ρ_α 7.14 / 5.80 / 9.99, and at 1e-5 **f = +1.4e-2 kPa** committed |
| 1950/2 | isoComp | 1e-5 / 1e-4 / 1e-3 | ok: ρ_α end 5.76 / 5.02 / 4.01 | rc 0, 1 substep: ρ_α 4.12 / 0.86 / 0.29 |
| 1950/2 | shear | 1e-5 / 1e-4 / 1e-3 | ok: ρ_α end 6.70 / 6.16 / 4.29 | 14 / 1 / 1 substeps: 6.70 / 5.45 (f +3e-2) / 6.79 |

**What the reference does with them.**
- It does not refuse them. They are on the yield surface, and the equations are defined
  with α outside α^b (A15).
- It integrates what the model says: K_p < 0 pulls α back towards the bounding
  surface. **No increment up to 1e-3 brings it back inside.** Recovery would take tens
  of per-mille of strain at this confinement.
- For 1950/3 under shear the equations themselves are singular, and the reference
  stops.

**Consequences.**
- Neither answer ("still 4–7× outside", or "singular") is a state an integrator should
  hand to a global Newton.
- The C++ answers are worse: a one-substep teleport to ρ 0.29, or an f > 0 commit.
  Neither is an integration of anything.
- As WP-128 said, these states should be **refused with a named code**, and they must
  never be produced. §6.5 shows the continuum never produces them from an admissible
  start.

## 7. Verdict on the leading hypothesis (WP-128 / WP-129)

| mechanism | statement | the reference says |
|---|---|---|
| **G** (trigger) | α_in re-seated once per increment, so (α−α_in):n → 0 then < 0 inside the increment, and h → 1e10 → negative | **Supported, at two levels.** (i) Model level: the UW rule alone, integrated exactly, gives h < 0 in 146/480 ring runs, 27 inadmissible rate problems, 25 escapes (≤ 1.97) and unstable, drifting consistency. The paper's rule, event-driven, gives 0/960. (ii) Discrete level: in the smallest reproducer the continuum never has a < 0, but ME's first stage at h = 1e10 overshoots and *creates* a < 0 in stage 2 |
| **E** (enabler) | stress-only substep error | **Supported.** Every escape passes the stress-only test. At TolE 1e-8, ME converges to the oracle (up to U9) wherever it actually substeps |
| **F** | Λ < 0 (negative denominator or h < 0) treated as elastic, with the step factor uncapped | **Supported, and more than a compounder.** It makes err **exactly 0** (C++ trace, §6.3), so no tolerance catches it. It accounts for all 25 campaign-ME ring escapes and 7 of 8 at TolE 1e-8, and for 20–65 % stress-increment errors on *benign admissible* states at 20–100 kPa (§6.2 class C) |
| **U9** (new) | K, G frozen over the ModifiedEuler increment | 0.6 / 6 / 24 % of Δσ at δ 1e-5 / 1e-4 / 1e-3. It does not shrink with TolR |

**What WP-129 must do, as the oracle sees it.**
- Paper-style α_in reseating: at every plastic onset, and at (α−α_in):n → 0 inside the
  increment.
- α and z in the error norm.
- Classify each stage by the loading numerator N with ∂f/∂σ, and treat N > 0 with H ≤ 0
  as a refusal, not "elastic".
- Cap the step growth.
- Recompute K, G per substep, which removes U9.
- Refuse inadmissible states.

The oracle for all of this is `ring_variants()["uw_model"]`.

## 8. Open items

- **O1** Map E1–E14 to DM04's own equation numbers against the PDF. The full text was
  not available in this session.
- **O2** Match the Toyoura figures of DM04 point by point. The data is not digitised,
  so §6.1 is qualitative.
- **O3** `solver_failed` on 1950/3 under shear is the singular point itself. A
  dedicated `singular_loading_index` status would read better than the SciPy message.
- **O4** Chains of increments (drive the ring states along a path), and a
  global-Newton-like iterate replay. Only single increments were compared with the C++.
- **O5** The cross-check test gates class A only. Classes C/D are documented, not gated,
  because they measure a C++ defect that WP-129 will change.
- **O6** Speed: 0.2–5 s per increment at rtol 1e-10 in pure Python. That is fine for
  hundreds of points, not for a mesh.

## 9. How to reproduce

From the worktree root, with any CPython that has numpy + scipy:

```
python Ladruno_files/testbed/sanisand_reference/run_paper.py
python Ladruno_files/testbed/sanisand_reference/run_crosscheck.py && python Ladruno_files/testbed/sanisand_reference/analyse_crosscheck.py
python Ladruno_files/testbed/sanisand_reference/run_reproducer.py
python Ladruno_files/testbed/sanisand_reference/run_ring.py && python Ladruno_files/testbed/sanisand_reference/analyse_ring.py
<py3.12> -S Ladruno_files/testbed/sanisand_reference/probe_cpp_trace.py
python -m pytest tests/test_sanisand_reference.py tests/test_sanisand_reference_crosscheck.py
```

The C++ side needs the WP-127 build (`ladrunoSANISANDReplay`). Its location is set in
`Ladruno_scripts/sanisand_reference/cxx.py` and can be overridden by environment
variable.
