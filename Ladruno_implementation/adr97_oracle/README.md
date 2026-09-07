# ADR-97 P0 — closest-point return map oracles

Reference implementations, in numpy, of the **closest-point (CPPM)** return map
and its **consistent (algorithmic) tangent** for the three `ASDPlasticMaterial3D`
families ADR-97 targets, plus the element-level finite-difference driver that
measures the shipped binary against its own residual.

Everything here mirrors the **shipped headers**, not a textbook: tension-positive
stress (`meanStress() == trace/3`), Voigt order `[11 22 33 12 23 13]` with
engineering shear, VOIGT (doubled-shear) YF/PF derivatives (ADR-94 wp/94c B5),
`LinearIsotropic3D_EL` elasticity, and the truncated
`SQRT_2_over_3 = 0.816496580928` literal.  Each file's docstring names the header
lines it reproduces and records every convention wart found there.

Regenerate the transcript:

```
LADRUNO_OPENSEES_QUIET=1 \
PYTHONPATH=<worktree>/dist/bin \
python3.12 Ladruno_implementation/adr97_oracle/run_all.py    # -> reference_output.txt
```

Every gate is an `assert` inside the script, so a clean exit **is** the pass.
The current transcript is `reference_output.txt` (all five exit 0), taken against
build `3622d6214ef4cdeb8cf65a102ee35f6cd9973337`.

---

## What each file pins

| file | what it pins | printed reference numbers | C++ test that will consume them (P1–P3) |
|---|---|---|---|
| `asd_common.py` | shared conventions + the Newton / consistent-tangent machinery. Jacobians come from a **complex-step** derivative (exact to machine precision) and are cross-validated once per family against a central difference. | complex-step vs FD Jacobian: `1.391e-13` (VM+AF), `7.843e-14` (DP cone) | — (support module) |
| `cppm_vm.py` | von Mises closest point: (a) perfect, (b) linear isotropic `H=7000`, (c) **Armstrong–Frederick solved implicitly at n+1**; closed-form tangent for (a)/(b), Jacobian-derived tangent for (c); three paths (triaxial, simple shear, rotating normal). | see **VM block** below | `tests/test_adr97_cppm_vm.py` |
| `cppm_dp.py` | Drucker–Prager cone **and** apex closest point, associated and non-associated (`etabar != eta` → unsymmetric tangent), perfect and linear hardening; the apex tangent is rank 0 (perfect) / rank-1 volumetric (hardening); and the **exact vs header apex-region test**. | see **DP block** below | `tests/test_adr97_cppm_dp.py` |
| `cppm_mc.py` | Mohr–Coulomb in principal-stress space (Clausen 2006/2007): return to face / edge σ1=σ2 / edge σ2=σ3 / apex, selected by the boundary-plane test on the elastic predictor, Koiter tangent per region transformed to 6D Voigt with the principal-rotation ("T") term and its degenerate limit. ψ=10° < φ=30°. | see **MC block** below | `tests/test_adr97_cppm_mc.py` |
| `path_independence.py` | the gate that separates CPPM from the shipped cutting plane: step-halving of a 3-leg path (both O(h), same limit) **and** the per-step uniqueness of CPPM vs the iterate-path dependence of `Backward_Euler` for AF. | see **path block** below | `tests/test_adr97_path_independence.py` |
| `fd_tangent_driver.py` | element-level free-DOF FD tangent of the **binary's own** assembled residual vs `printA('-sparse','-ret')`. Importable as `fd_check(mat_fn, ...) -> dict(rel_err, rel_fro, K_fd, K_asm, ...)`. | see **element block** below | `tests/test_adr97_tangent_element.py` |
| `run_all.py` | regenerates `reference_output.txt` (stderr discarded — the C++ side prints a page per material construction on `opserr`). | — | — |

---

### VM block — `cppm_vm.py`

E = 70000, ν = 0.3, k₀ = 30, 10 steps per leg.  Final state after each path:

| case | path | σ | α | k | dΛ (last step) | FD rel err |
|---|---|---|---|---|---|---|
| (a) perfect | triaxial | `[60, 60, 90, 0,0,0]` | 0 | 30 | `3.1843366656e-04` | `1.531e-11` |
| (a) perfect | simple shear | `[0,0,0, 17.3205080757, 0,0]` | 0 | 30 | `2.8284271247e-04` | `1.684e-11` |
| (a) perfect | rotating normal | `[0,0,0, 0.6597367945, 17.3079388537, 0]` | 0 | 30 | `2.1190632791e-04` | `2.273e-11` |
| (b) H=7000 | triaxial | `[55.2147239264, 55.2147239264, 99.5705521472, 0,0,0]` | 0 | `44.3558282209` | `2.9303711647e-04` | `9.622e-12` |
| (b) H=7000 | simple shear | `[0,0,0, 24.5280749163, 0,0]` | 0 | `42.4838719668` | `2.6028470473e-04` | `1.158e-11` |
| (b) H=7000 | rotating normal | `[0,0,0, 2.7663534077, 27.3507956089, 0]` | 0 | `47.6146636537` | `1.9392521616e-04` | `1.248e-11` |
| (c) AF ha=15000 cr=300 | triaxial | `[49.6232363748, 49.6232363748, 110.7535272504, 0,0,0]` | `[-10.3767636252, -10.3767636252, 20.7535272504, 0,0,0]` | 30 | `2.7381670398e-04` | `1.190e-10` |
| (c) AF | simple shear | `[0,0,0, 45.3906185722, 0,0]` | `[0,0,0, 28.0701104965, 0,0]` | 30 | `2.0547089463e-04` | `7.050e-11` |
| (c) AF | rotating normal | `[0,0,0, 23.7842127217, 37.7654987099, 0]` | `[0,0,0, 21.6973739054, 20.5711652319, 0]` | 30 | `1.4942296298e-04` | `1.841e-11` |

Closed form vs Newton-Jacobian tangent for (a)/(b): `max|Δσ| <= 3.4e-13`,
`rel_fro(C) <= 7.9e-15`.  Newtons: 1 iteration (a)/(b), 3 iterations (c);
`|f| <= 7.1e-15` on every committed state.

**Convention contrast** (rotating-normal path, AF): the header's `ha * mdev`
(engineering/doubled shear) gives `σ12 = 23.7842127217`, `σ23 = 37.7654987099`;
the tensor-consistent variant gives `13.4865193`, `29.2286435` — a 43 % / 23 %
difference on a sheared path.

### DP block — `cppm_dp.py`

E = 30000, ν = 0.25, ξ_c = 20, η = 0.4 ⇒ apex at p = **+50** (tension positive);
K = 20000, G = 12000.

| case | region | σ | k | dΛ (last) | ‖C−Cᵀ‖/‖C‖ | rank C | FD rel err |
|---|---|---|---|---|---|---|---|
| associated, perfect | cone | `[-26.1416983071, -26.1416983071, -94.7352064897, 0,0,0]` | 0 | `1.9415646030e-04` | `2.13e-17` | 5 | `1.381e-11` |
| non-assoc `etabar=0.2`, perfect | cone | `[-26.6815746982, -26.6815746982, -89.9603702958, 13.9585578524, 0,0]` | 0 | `2.4142148763e-04` | `1.854e-01` | 5 | `1.893e-11` |
| non-assoc, H = 500 | cone | `[-26.4735236551, -26.4735236551, -90.1380619914, 14.0436481624, 0,0]` | `0.2701984426` | `2.3633849288e-04` | `1.814e-01` | 6 | `2.226e-11` |
| associated, perfect | apex | `[50,50,50,0,0,0]` | 0 | `9.0000000000e-04` | 0 | **0** | `0.000e+00` |
| non-assoc, H = 500 | apex | `[50.6296305482]*3` | `0.2518522193` | `1.7484847733e-03` | `1.32e-16` | **1** | `1.701e-10` |

Closed form vs Newton-Jacobian: `max|Δσ| <= 7.1e-15`, `rel_fro(C) <= 2.3e-16`
in all five.

**Region-test finding.**  The exact (elastic-metric) apex boundary slope is
`(K*etabar + h_k/eta)/G`; the header uses the Euclidean `eta`.  They disagree in
BOTH directions and the disagreement is not benign:

* `etabar = 0.2` (exact slope 0.3333 < header 0.4): a trial at
  `(p−p_apex)/q = 0.36` is classified CONE by the header, and the cone return it
  then runs gives `sqrt(J2)_{n+1} = −0.470588` — an **inadmissible negative
  sqrt(J2)**.  At 0.37: `−0.647059`.
* `etabar = eta = 0.4` (exact slope 0.6667 > header 0.4): trials at
  `(p−p_apex)/q = 0.45` and `0.60` are classified APEX by the header although the
  correct return is to the cone (`sqrt(J2)_{n+1} = 3.42` and `1.05`).

### MC block — `cppm_mc.py`

E = 30000, ν = 0.25, φ = 30°, ψ = 10° (non-associated), c = 10; apex at
`17.3205080757`.  Warrant: the header's invariant `f` equals the principal-stress
`f = 0.5[(s1−s3)+(s1+s3) sinφ] − c cosφ` to `2.132e-14` over 2000 random states.

| trial σ_tr | region | principal return y | dΛ | rank C | FD rel err |
|---|---|---|---|---|---|
| `[-10,-40,-100,0,0,0]` | face | `[-20.4812370658, -39.6838547781, -96.0847273489]` | `[5.998546e-04]` | 5 | `2.022e-11` |
| `[-10,-40,-100,15,-8,5]` | face | `[-20.2611835701, -44.1530665087, -95.4245668617]` | `[9.443248e-04]` | 5 | `1.414e-10` |
| `[-25,-25,-140,0,0,0]` | **edge σ1=σ2** | `[-33.0520613099, -33.0520613099, -133.7972000812]` | `[4.75162e-04, 4.75162e-04]` | 3 | `2.869e-11` |
| `[10,-95,-95,0,0,0]` | **edge σ2=σ3** | `[-18.2208073609, -89.303438234, -89.303438234]` | `[8.075564e-04, 8.075564e-04]` | 3 | `2.499e-11` |
| `[45,45,45,0,0,0]` | apex | `[17.3205080757]*3` | least-squares, sign-free | **0** | 0 (abs) |
| `[44,40,36,0,0,0]` | apex | `[17.3205080757]*3` | least-squares, sign-free | **0** | 0 (abs) |

FD stencil h = 1e-8 relative on the total strain; the printed boundary-plane
margins show every trial sits well inside its region.  Associated (ψ=φ) face
tangent asymmetry `4.420e-17`; non-associated `2.972e-01`.

Two facts the C++ must respect: at the apex under NON-associated flow the three
Koiter multipliers are **not** all positive (the cone of return directions no
longer contains the hydrostatic direction), so the apex region must be selected
by the boundary planes and never by an active-set search on `dΛ >= 0`; and the
apex tangent is the exact **zero** matrix for a non-hardening MC.

### Path block — `path_independence.py`

Step-halving of a 3-leg path, `|σ(N) − σ(N=640)|_inf`:

| | N=5 | N=10 | N=20 | N=40 | limits agree to |
|---|---|---|---|---|---|
| VM+AF, CPPM | `4.5217e-01` | `2.2993e-01` | `1.1415e-01` | `5.5453e-02` | `1.5050e-02` |
| VM+AF, cutting plane | `1.9550e+00` | `9.8499e-01` | `4.8964e-01` | `2.3846e-01` | ” |
| DP+linear, CPPM | `2.6247e-01` | `1.3683e-01` | `6.9039e-02` | `3.3829e-02` | `9.9476e-14` |
| DP+linear, cutting plane | identical to CPPM to 2.8e-14 | | | | ” |

Both maps are O(h) (ratios 1.92–2.06) and converge to the same limit; at finite
step the cutting plane is **4.3× less accurate** on VM+AF.

Same step, same start, different iterate path:

* **CPPM invariance** over three Newton start guesses: `1.137e-12` (gate ≤ 1e-9).
* **Cutting-plane spread** over four equally valid iterate paths (first corrector
  × 1.0, 0.5, 0.25, 1.5), all committing `|Φ| < 1e-10`:
  **σ `9.147082e+00`, α `5.683938e+00`** — on a stress of ~46.
* VM + **linear** hardening: cutting-plane spread `7.105e-15` (path INdependent).
  DP + linear: `7.105e-15`.  **This is why ADR-94 H6 could not see the defect**:
  the shipped cutting plane is only path dependent when `h` is not constant along
  the iterates, i.e. for Armstrong–Frederick (22 of 46 specializations carry AF).

### Element block — `fd_tangent_driver.py`

`stdBrick`, one step, `dEps_zz = -3.6210526316e-03`,
`σ = [-196.842105, -196.842105, -240, 0,0,0]`, 4 free DOFs.

| material | rel_err (max) | rel_fro |
|---|---|---|
| `ElasticIsotropic` (self-test) | `4.852e-12` | — |
| `Backward_Euler` / `Continuum` | **`0.573447`** | `0.589373` |
| `Backward_Euler` / `Secant` (the default) | `0.799355` | `0.730083` |
| `Backward_Euler` / `Elastic` | `1.025263` | `0.917170` |
| `Backward_Euler` / `Numerical_Algorithmic_FirstOrder` | `0.045735` | `0.047005` |
| `Backward_Euler` / `Numerical_Algorithmic_SecondOrder` | `0.045735` | `0.047005` |
| 12-DOF sheared rig / `Continuum` | `0.670133` | `0.472726` |
| 12-DOF sheared rig / `Numerical_Algorithmic_FirstOrder` | `0.053446` | `0.037702` |
| `Closest_Point` / `Algorithmic` | *not accepted by this build — P1 adds it; gate is ≤ 1e-6* | |

These reproduce ADR-94 M3's 57.3 / 79.9 / 102.5 / 4.6 / 4.6 % **without any
numpy reference**, from the binary's own assembled residual — which independently
confirms that ADR-94's reference tangent was the consistent tangent.
`0.573447` is the **negative control** the ADR-97 tests contrast against.

---

## Header findings recorded by these oracles

1. **AF hardening mixes conventions** (`AllASDHardeningFunctions.h:157`):
   `derivative = ha*mdev - cr*mdev_eq*alpha_dev` adds the ENGINEERING-shear
   (doubled) `mdev` to a STRESS-like back stress that is then consumed by
   `dev(sigma) - alpha` and `tensor_dot_stress_like`.  On sheared paths α's shear
   slots grow 2× relative to its normal slots.  There is also no `2/3` on `ha`
   (the `(2./3.)*ha` variants are commented out at lines 158–160).  The same
   policy still contains `cout` debug spam (lines 142–146, 152, 163) on every
   hardening evaluation.
2. **Drucker–Prager's cohesion IV is not in its own `f`**
   (`DruckerPrager_YF.h:26` is commented out) while `yf_hardening` contributes
   `df/dk = -1` times that IV's rate — a hardening term matching no term of `f`.
   The oracle pins the self-consistent `f = sqrt(J2) + eta*p - (xi_c + k)`.
3. **The DP / MC apex region tests are Euclidean, not elastic-metric** — the
   headers say so in their own comments; the numbers above quantify the cost.
4. **`MohrCoulomb_PF.h:83` converts the cohesion to radians**
   (`GET_PARAMETER_VALUE(MC_c)*M_PI/180`).  Harmless today only because `c`
   enters `g` as an additive constant that differentiation kills.
5. **Mohr–Coulomb has no edge/apex algebra at all**: for `|theta| >= 29°` the
   header silently substitutes a Drucker–Prager gradient, and the default path is
   a numerical central difference of `f`.  `cppm_mc.py` is the reference for what
   `Closest_Point` must do instead.
