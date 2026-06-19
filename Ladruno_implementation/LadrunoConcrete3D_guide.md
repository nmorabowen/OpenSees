---
title: "LadrunoConcrete3D — CDPM2 solid-concrete implementation guide"
project: Ladruno
type: implementation guide
status: oracle COMPLETE through P2 damage (P0 surface · P1 return-map/hardening/tangent · P2a-e dual damage + analytic tangent) + adversarial-review-hardened; C++ KERNEL has the return map, the analytic effective tangent, AND the damage stress update (all g++ byte-verified); STILL OWED = the analytic damaged tangent in the kernel (P3b) + the nDMaterial wrapper (classTag 33017). The model is NOT yet callable from Tcl/Python.
related:
  - "[[31_ladruno_concrete3d_adr]]"          # the decision record
  - "[[LadrunoConcrete3D_dev_handoff]]"      # the implementer's brief (build-PR checklist)
  - "[[LadrunoConcrete3D_user_guide]]"       # the short API spec
  - "[[19_ladruno_rc_shell_adr|LadrunoRCConcrete]]"  # the shell/MCFT sibling (33015)
  - "[[project_ladruno_concrete3d]]"          # agent-memory pointer
updated: 2026-06-19
---

# LadrunoConcrete3D — CDPM2 solid-concrete implementation guide

> [!warning] Build status (2026-06-19)
> The **constitutive law is fully verified** by a numpy oracle through P2 dual damage + the analytic
> damaged tangent, hardened by two adversarial-review rounds. The **C++ kernel**
> (`LadrunoConcrete3DKernel.h`) implements the effective-stress return map, the analytic effective
> tangent, **and the dual-damage nominal-stress update** — all g++ byte-verified against the oracle.
> **Not yet shipped:** the analytic *damaged* tangent in the kernel (P3b) and the
> `nDMaterial LadrunoConcrete3D` wrapper (classTag **33017**). So the model is **not callable from
> Tcl/Python yet** — Part III is the *intended* interface plus how to exercise the verified core today.

`LadrunoConcrete3D` is a **CDPM2-grade** (Grassl, Xenos, Nyström, Rempling, Gylltoft 2013;
arXiv:1307.6998) 3D **solid** concrete `nDMaterial`: effective-stress **Menétrey–Willam plasticity**
(a smooth three-invariant surface with Lode dependence and an ellipsoidal compression cap) +
**non-associated flow** + a **confinement-dependent ductility measure** + **dual scalar damage**
(tension `ω_t` / compression `ω_c`) on the spectral tension/compression split, **crack-band
regularized**. It is the solid/triaxial sibling of [[19_ladruno_rc_shell_adr|`LadrunoRCConcrete`]]
(33015, the shell/MCFT-membrane axis): same ASDConcrete3D damage lineage, different physics axis.

## Contents

- **Part I — Theory** — §1 the gap · §2 the MW surface · §3 non-associated flow + semi-implicit return ·
  §4 hardening + ductility · §5 the non-symmetric tangent · §6 dual damage + spectral split ·
  §7 crack-band `G_f`/`G_c` · §8 unilateral crack-closure · §9 the compression→tension coupling
- **Part II — Architecture** — §10 oracle-first · §11 one kernel · §12 the damaged update · §13 the
  analytic damaged tangent · §14 state + the effective/nominal split · §15 the g++ byte-check ·
  §16 the gate ladder
- **Part III — Usage** — §17 command grammar (intended) · §18 parameters + calibration · §19 the
  unsymmetric-solver requirement · §20 exercising the core today · §21 phases + recorders · §22 units
- **Part IV — Quirks** — §23 consolidated quirks
- **Appendices** — A phase/PR history + V&V · B references

---

# Part I — Theory

## 1. The constitutive gap (vs ASDConcrete3D)

`ASDConcrete3D` (the upstream model `LadrunoConcrete3D` is built to surpass) is **damage-only**: its
`compute()` is a C0-predictor + a hardening/softening damage curve, with a **Lubliner** envelope that
has **no Lode dependence and no cap**. Confinement is *emergent* — a calibrated band, ≈Mander ±5% only
out to `p/fc ≈ 0.2` — not constitutive. There is one `lch`, no dilation knob, no unilateral recovery.

CDPM2 closes the gap with a genuine **plasticity + damage** split:

- a **three-invariant Menétrey–Willam surface** — Lode angle `θ` *and* an ellipsoidal compression cap;
- **non-associated flow** (a dilatancy knob `D_f`) — concrete dilates less than associated flow predicts;
- a **ductility measure** `x(σ̄_V)` that makes peak strain grow with confinement *by mechanism* — the
  generalization of Mander's pressure→ductility law, so it extrapolates across triaxiality;
- **dual scalar damage** `ω_t`/`ω_c` on the spectral split → the peak + softening + unilateral recovery.

The **headline decision** (ADR 31): build a **material**, header-only and OpenSees-free
(`LadrunoConcrete3DKernel.h`, the LadrunoJ2Kernel "one core, many views" doctrine). The IMPL-EX /
serialization / crack-band *scalar machinery* is reused from ASDConcrete3D; the **return map donor is
LadrunoJ2Kernel, NOT ASD** — ASD has no return map to port. (A 15-agent ADR-scoping pass caught the
"reuse ASD's return map" premise as false.)

## 2. The Menétrey–Willam effective-stress surface

All plasticity is in **effective** (undamaged) stress `σ̄`. With `ξ = I₁/√3`, `ρ = √(2J₂)`, Lode angle
`θ ∈ [0, π/3]`, and **compression-negative** sign convention (`fc`, `ft` entered as positive
magnitudes), the yield function is Grassl 2013 **Eq.18**:

```
f(ξ,ρ,θ;κ_p) = { [1−q_h1(κ_p)]·A_V² + (√1.5 ρ/fc) }²  +  m0·q_h1²·q_h2·R  −  q_h1²·q_h2²
  A_V = ρ/(√6 fc) + σ̄_V/fc                     (the ellipsoidal-cap base; σ̄_V = I₁/3)
  R   = ρ·r(θ,e)/(√6 fc) + σ̄_V/fc              (the m0-friction bracket; carries the Lode r)
  m0  = 3(fc²−ft²)/(fc·ft) · e/(e+1)            (Eq.20 — friction, fixed by the tensile meridian)
```

At `q_h1 = q_h2 = 1` this reduces to the **failure surface Eq.21** = Menétrey–Willam. The eccentricity
`e ∈ (0.5, 1]` (0.5 is the convexity limit) is the **deviatoric out-of-roundness**; the Willam–Warnke
Lode function `r(θ,e)` (Eq.19) satisfies `r(0)=1/e` (tensile meridian), `r(π/3)=1` (compressive),
so the deviatoric-section ratio is exactly `e`. `e` is fixed from the **Kupfer** equibiaxial strength
`fcc/fc` (default 1.16 → the canonical `e ≈ 0.52`), NOT a free fit knob.

**Machine-exact normalization (P0 gate):** uniaxial compression `(−fc,0,0)` and uniaxial tension
`(+ft,0,0)` land on `f=0` to ~1e-16, the tension condition *fixing* `m0`. `e→1` degenerates the
deviatoric section to a circle (Drucker-Prager-like meridian). **Key identity:** `ξ/(√3 fc) = σ̄_V/fc`
exactly, so the code uses `ξ`/`ρ` throughout.

## 3. Non-associated flow + the semi-implicit return map

The return map works in **invariant `(ξ,ρ)` space with the Lode angle `θ` FROZEN at the trial value**
— the deviatoric return is then radial, which **sidesteps the `∂θ/∂σ` Lode-corner singularity**. This
"semi-implicit" freeze is **exact for a single radial return**; the `O(Δθ)` tax is multi-step only.

- **Perfect plastic (P1):** a 3-unknown Newton on `(ξ, ρ, Δλ)` with an **analytic 3×3 Jacobian**:
  ```
  ξ  = ξ_tr  − 3K·Δλ·m_v       m_v = D_f·m0/(√3 fc)         (volumetric flow, dilatancy-scaled)
  ρ  = ρ_tr  − 2G·Δλ·m_s       m_s = 3ρ/fc² + m0/(√6 fc)    (deviatoric flow — Lode-INDEPENDENT)
  f(ξ,ρ,θ_tr) = 0
  ```
  The **flow potential is Lode-independent** (Eq.22-25 — the yield function has `r`, the flow does
  not), and `D_f ∈ (0,1]` scales only the *volumetric* (dilatant) flow (`D_f=1` ⇒ associated). If the
  deviatoric return drives `ρ<0`, project to the **hydrostatic-tension apex** `ξ = √3 fc/m0`, `ρ=0`.
- **Effective-stress plasticity is MONOTONIC** — it hardens up to the failure surface (reached exactly
  at `κ_p = 1`, where uniaxial compression sits at `σ̄₁₁ = −fc`) and **never peaks or softens**. The
  **peak and softening are entirely the DAMAGE part** (§6). Do not gate "peak == fc"; gate "`σ̄` at
  `κ_p = 1`".
- **Full tensor via the spectral split:** eigendecompose the trial stress → return the principal
  stresses (radial return preserves eigenvectors) → recompose `σ̄ = V·diag(σ̄ₚ)·Vᵀ`.

## 4. Hardening + the confinement-ductility measure

Hardening is the 4-unknown extension `(ξ, ρ, Δλ, κ_p)`, with `f` evaluated at Eq.18 and:

```
q_h1(κ_p)  Eq.30   pre-peak size: q_h0 → 1 over κ_p ∈ [0,1], slope H_p at κ_p=1⁻
q_h2(κ_p)  Eq.31   post-peak: 1 (κ_p<1), then 1 + H_p(κ_p−1)
κ_p evol   Eq.32   κ̇_p = (Δλ·‖m‖/x_h(σ̄_V))·(2 cos θ)²
x_h(σ̄_V)  Eq.33-36  R_h = −σ̄_V/fc − 1/3 (Eq.34) ⇒ more ductile in compression (σ̄_V<0 ⇒ R_h>0 ⇒ larger x_h)
```

`x_h` is the **confinement-ductility measure** — the strain at `κ_p=1` grows ~0.0012 (unconfined) →
~0.0032 (`p/fc=0.1`) purely through `x_h`, **"Mander by mechanism"** (not a backbone). The `Ah/Bh/Ch/Dh`
(Eq.33) set the absolute peak strain and its confinement growth — **recalibrate to your concrete's
measured peak strains**; keep `H_p` small (large `H_p` makes hardening non-monotone). The hardening
return shares the `Voce+linear` `σ_y` contract with `LadrunoJ2` via `LadrunoHardening.h`.

## 5. The consistent tangent — non-symmetric, always

The algorithmic tangent `dσ̄/dε` is the **spectral lift** (de Souza Neto Box A.6) of the principal
Jacobian (implicit-function theorem on the inner Newton residual; the Lode directional gradient by an
isolated scalar micro-FD, the rest closed-form). It is **NON-symmetric** for two independent reasons:

1. **non-associated flow** — `‖C−Cᵀ‖/‖C‖ ≈ 0.46` at `D_f=0.3`;
2. **the frozen-eigenvector spectral recompose** — even in the *associated* limit (`e=1, D_f=1`) the
   asymmetry is `≈0.024` (~20× smaller but nonzero), because `V` held at the trial drops the
   eigenprojection-spin `dV/dε`. This was *falsified* not to be the Lode θ-freeze: a principal-space
   associated off-meridian state is machine-symmetric (~2e-10), and the full-tensor asymmetry scales
   **linearly with shear** (→0 as shear→0).

⇒ **Tier-1 requires an unsymmetric solver UNCONDITIONALLY** (`UmfPack`/`FullGeneral`), not only for
`D_f<1`. There is no runtime guard — see §19.

## 6. Dual scalar damage + the spectral T/C split (where the peak comes from)

Damage knocks the monotonic effective stress down to the **nominal** stress (Grassl 2013 §2.3,
**Eq.1**):

```
σ = (1−ω_t)·⟨σ̄⟩₊ + (1−ω_c)·⟨σ̄⟩₋        (Macaulay split of the EFFECTIVE stress eigenvalues)
```

Both damage variables are driven by **one equivalent strain `ε̃` (Eq.37)** and the
**compressive weight `α_c` (Eq.46)** (0 in tension, 1 in compression):

| | tension | compression |
|---|---|---|
| history rate | `κ̇_dt = ε̃̇` (FULL, Eq.43) | `κ̇_dc = α_c·ε̃̇` (Eq.47) |
| plastic part | `κ̇_dt1 = (1/x_s)‖ε̇_p‖` (Eq.44, **no `α_c`**) | `κ̇_dc1 = (α_c·β_c/x_s)‖ε̇_p‖` (Eq.48) |
| damage-scaled | `κ̇_dt2 = κ̇_dt/x_s` (Eq.45) | `κ̇_dc2 = κ̇_dc/x_s` (Eq.49) |

with the inelastic strain `ε_i = κ_d1 + ω·κ_d2` (Eq.52) and softening ductility
`x_s = 1 + (A_s−1)·R_s`, `R_s = −√6 σ̄_V/ρ̄` for `σ̄_V≤0` else 0 (Eq.56-57; `x_s = A_s` in uniaxial
compression = the confinement-ductility hook). **`β_c = 1` in this monotonic slice** (Eq.50 is the
cyclic damage↔plasticity transition — deferred, §9). The equations were **re-pinned from the arXiv
full text** before coding (the recurring "fetch the actual eqs, don't guess" lesson — the asymmetry
`κ̇_dt = ε̃̇` vs `κ̇_dc = α_c·ε̃̇` is literal).

**Onset is at the failure surface:** `ε̃ = ε0 = ft/E` ⟺ `q_h2 = 1` ⟺ **`κ_p = 1`** — so pre-peak is
pure plasticity and post-peak is damage, no double-counting. Each `ω` is solved implicitly
(`(1−ω)·D = f·exp(−ε_i/ε_f)`) by a **bracketed (bisection-safeguarded) root** — a raw clamped Newton
spuriously *healed* the crack on steep softening (an adversarial-review CRITICAL).

**The per-channel softening drive `D`** is the **extreme effective principal of that sign**
(`max⟨σ̄ᵢ⟩₊` for tension, `max⟨−σ̄ᵢ⟩₊` for compression). This reduces the unified tensor update to the
validated 1D drivers **exactly** in tension and to the onset-crossing step in compression. (Do NOT
feed it `E·κ_d`: right for tension, gives the wrong peak `ft≠fc` for compression.) The genuinely
multiaxial apportioning — extreme-principal vs `‖σ̄_t‖` norm — is a documented refinement.

## 7. Crack-band regularization (`G_f`/`G_c`)

Exponential softening `σ_nom = f·exp(−ε_i/ε_f)` with the **crack-band length** `ε_f = G_f/(ft·lch)`
(tension), `ε_fc = G_c/(fc·lch)` (compression). By construction `∫σ_nom dε_i = f·ε_f = G/lch`, so the
**damage dissipation per crack area is `lch`-objective** — the softening law's `ε_f` wiring is gated at
3 `lch` values (relative err ~1e-4). **Honest caveat (D3/C3, reported not gated):** CDPM2 regularizes
the *damage* softening only; the effective-**plasticity** dissipation is per-volume and *un*-regularized,
so the FE-visible *total* fracture energy is ~30% `lch`-dependent. That is a CDPM2 damage-only-
regularization characteristic, small in realistic regimes; full plastic-dissipation regularization is a
documented follow-on.

> **Gate-honesty note (review):** integrating `σ_nom` over `ε_i` (the softening law's *own* abscissa)
> returns `G/lch` *tautologically* — it catches an `ε_f`-wiring error but is NOT an independent
> objectivity proof. The honest FE-objectivity number is the D3/C3 total-energy spread, reported.

## 8. Unilateral crack-closure — automatic, zero extra state

Because the spectral split is **recomputed from the converged effective stress every step**, a
principal that flips negative (a crack closing) is routed into the `(1−ω_c)` channel and is **no longer
multiplied by `(1−ω_t)`** — full stiffness recovers automatically, tier-independently, with **no stored
projector and no `s_rec` knob** (ADR §4.3 BLOCKING). The trap to avoid (ASDConcrete3D's IMPL-EX
frozen-`PT`): split once and reuse a *frozen* projector — then a cracked-then-closed material carries
almost no compressive stress, and the cyclic value proposition dies exactly where it's sold. Freeze
*only* the damage-driver extrapolation in IMPL-EX, **never the recomposition split**.

## 9. The compression→tension damage coupling (known limitation, P2f)

Per literal CDPM2 **Eq.43** (`κ̇_dt = ε̃̇`, the *full* equivalent strain, **no `(1−α_c)` factor** —
re-confirmed from the arXiv source), a compression excursion accumulates `κ_dt` and so **pre-damages a
subsequent tension reload**. The loss is **graded by the compression damage**: a specimen compressed
only into the hardening range (`κ_p < 1`, no compression damage) retains *full* tensile strength
(peak ≈ `ft`); past compression-damage onset it degrades progressively (`κ_p=1.5` → ~0.6`ft`;
fully crushed → ~0). The **monotonic** tension and compression responses are correct; this **cyclic
T/C coupling** (the dropped `β_c` Eq.50 + the open `α_t`-weighting question: literal-CDPM2 full-`ε̃`
vs a tensile-plastic-strain projection) is the **P2f** increment. It is tracked by the **DT5 diagnostic**
(reported, not gated).

---

# Part II — Architecture

## 10. Oracle-first methodology — the numpy spec is the source of truth

Every increment is pinned in a **numpy oracle** (`tests/_testbed/concrete3d_ref.py`) *before* any C++
is written; the C++ is then byte-verified against the oracle's dumped numbers. This is the
LadrunoJ2/Lemaitre pattern. The oracle is the **verified specification** — the gates (`run_p0_gate` …
`run_p2e_gate`, Part-II §16) live there, the pytest file `tests/test_ladrunoConcrete3D_material.py`
calls them (with the module-level `pytest.mark.zone_a` so CI actually runs them — a missing marker once
silently deselected the whole file, the #249 "CI false-green" lesson). Each increment is its **own PR
off `ladruno`** (the fork auto-merges fast, so old branches strand — always fresh-branch).

## 11. One kernel, the C++ port — `LadrunoConcrete3DKernel.h`

Header-only, OpenSees-free, **all-inline** (no `.cpp`), `{00,11,22,01,12,02}` Voigt with **true-tensor**
off-diagonals. `Params` holds `E,ν,fc,ft,e,m0,Gf,Gc,Df,As,qh0,Hp,Ah..Dh,η,lch`. What is implemented:

| Layer | Functions | Verified |
|---|---|---|
| surface | `invariants`, `lodeR`, `m0Of`, `yieldF`/`yfInv`/`yfInvHard`, `eccentricityFromKupfer` | algebraic identities |
| effective return | `returnMapPrincipal` (3×3 Jac), `returnMapHardening` (4×4 Jac), `eig3sym`, `returnMapTensor` | g++ vs oracle ~1e-13 (pp) / ~1e-8 (hard) |
| effective tangent | `principalJacobian`, `consistentTangent` | g++ ~2e-10 (pp) / ~7e-7 (hard) |
| **damage kinematics** | `equivStrainGeneral` (Eq.37), `alphaCompression` (Eq.46), `damageDrivers` (`ε̃,α_c,x_s`), `solveOmegaBracketed` | (see §12) |
| **damage update** | `plasticStrain6`, `damagedUpdate` (Eq.1 nominal), public `returnMap` | **g++ vs oracle ~1e-14** |

**Still owed:** the **analytic damaged tangent** in the kernel (§13, P3b) and the `nDMaterial` wrapper.

## 12. The damaged constitutive update — effective → nominal

`returnMap(strain, in_State) → out_State` does two stages, mirroring the oracle `damaged_step_tensor`:

1. **Effective-stress return** from the committed *effective* state `in.sigEff` (NOT the nominal
   `in.sig`): `returnMapTensor` → `σ̄`, `κ_p`, and the P1 effective tangent. `sigEffImplicit = σ̄`
   (kept separate so a future IMPL-EX / finite-strain choice never corrupts the recovery — the
   LogStrain `bᵉ` fix, ADR R3).
2. **`damagedUpdate`** → the nominal stress: re-eig `σ̄`, compute `ε̃`/`α_c`/`x_s`, accumulate the six
   damage histories (clamped at the onset `ε0`; the plastic-strain increment is the **tensor Frobenius
   norm** of `ε_p = ε − C⁻¹:σ̄`, via `plasticStrain6` = the isotropic closed form `ε_p,ij = ε_ij −
   [σ̄_ij(1+ν) − ν·tr(σ̄)δ_ij]/E`, no linear solve), solve `ω_t`/`ω_c` (bracketed, **physical-floored**
   §23), and recompose `σ = V·diag[(1−ω_t)⟨σ̄ₚ⟩₊ + (1−ω_c)⟨σ̄ₚ⟩₋]·Vᵀ`.

The damage drivers are **frame-invariant**, so eigenvalue order/sign do not matter; the recompose is
byte-stable across eig conventions.

## 13. The analytic damaged tangent — `C = D_dam : C_eff − σ̄_t⊗∂ω_t/∂ε − σ̄_c⊗∂ω_c/∂ε`

Derived + FD-verified in the **oracle** (`damaged_tangent_analytic`, P2e) to rel ~1e-10 against the P2d
numerical reference; the **kernel port is P3b** (the kernel currently returns the *effective* tangent —
over-stiff on softening). The three pieces (the C++ must replicate exactly):

- **`C_eff`** — the P1 effective consistent tangent (numerical in the oracle; *analytic* in the kernel,
  already shipped). The C++ assembles its own `C_eff` and adds the damage linearization.
- **`D_dam`** — `isotropic_tangent`: the spectral derivative of the per-principal damaged stress with
  `ω` **frozen** (de Souza Neto Box A.6 operational form `dY:S = Σ_a y'_a(E_a:S)E_a + Σ_{a≠b}G_ab E_a S E_b`,
  `G_ab = (y_a−y_b)/(λ_a−λ_b)` → `y'_a` by l'Hôpital at eigenvalue coalescence — the P1 spectral-tangent
  machinery).
- **`∂ω/∂ε`** — the implicit-function theorem on the bracketed `ω`-solve (`H = ∂F/∂ω = D[(1−ω)κ_d2/ε_f − 1]`),
  chained through the histories. The scalar sub-gradients `∂ε̃/∂σ̄`, `∂x_s/∂σ̄`, `∂α_c/∂σ̄` are isolated
  micro-FDs; `∂λ_extreme/∂σ̄` is the analytic eigenprojection (with the **Voigt `[1,1,1,2,2,2]`
  double-contraction weight** on the shear off-diagonals — §23); `∂‖Δε_p‖/∂ε` is closed form.

The damaged tangent is **degraded + INDEFINITE on the softening branch** (`C[0,0]<0`, `λ_min(symC)<0`)
— the concrete **Tier-2 IMPL-EX motivation** — and stays finite across a load reversal and (as a valid
subgradient) at the `σ̄_lat=0` Macaulay kink.

## 14. State + the effective/nominal split

`State` carries: `eps[6]`, `sig[6]` (committed **nominal**), `sigEff[6]` (committed **effective** —
drives the next return *and* the damage plastic strain), `kp`, and the **six damage histories**
`et_max, kdt1, kdt2, kdc, kdc1, kdc2`. The caller commits by copying `out → in`. The wrapper (pending)
serializes `Params`' fixed block + this `State` + the IMPL-EX scalars.

## 15. The g++ byte-verification harness

`gen_concrete3d_fixture.py` dumps the oracle's reference numbers to a committed text fixture
(`concrete3d_oracle_fixture.txt`); the standalone `concrete3d_kernel_check.cpp` (no OpenSees build)
compiles the header, reproduces every case, and diffs at cross-platform floors. Four blocks:
**PATH** (return-map stress, ~1e-13 pp / ~1e-8 hard), **TAN** (analytic effective tangent, per-entry +
asymmetry-norm, ~1e-6), **DMG** (the dual-damage nominal stress — tension/reversal/confined-compression/
shear — to **~1e-14**), and a **robustness fuzzer** (0 inadmissible converged returns over 180k trials;
the NaN-on-degenerate-params net). The pytest regenerates to a *tmp* path and asserts the committed
fixture is current (never self-referential). **Confined-compression damage states are built
STRESS-controlled** (a small confining `σ3` in `drive_damaged_unified`) — confined-*strain* paths sit
in the apex-fragile regime and diverge across platforms (a CI-red incident, §23).

## 16. The gate ladder (what each `run_*_gate` proves)

| Gate | Proves |
|---|---|
| `run_p0_gate` | MW surface fc-normalization: uniaxial C/T on `f=0` (~1e-16), `e` Lode identity, Kupfer `fcc/fc` |
| `run_p1_gate` | return-to-surface (never outside), uniaxial `fc`/`ft`, confined `fcc(σ3)`, dilatancy ordering, apex |
| `run_hardening_gate` | `q_h1/q_h2/x_h` identities, reduce-to-P1, `σ̄=fc` at `κ_p=1`, confinement ↑ ductility |
| `run_tangent_gate` | reduce-to-elastic, non-symmetric (non-assoc), the recompose-spin falsification, Taylor-quadratic, objectivity |
| `run_p2_gate` (P2a) | tensile `ω_t`: onset at `κ_p=1`, peak `ft` while effective monotonic, softens, `ε_f` wiring |
| `run_p2b_gate` | compressive `ω_c` + `α_c` split + `ε̃==ε0` on the surface + `G_c` wiring |
| `run_p2c_gate` | the unified tensor update: split partition, **direct unilateral routing** (compressive entries `ω_t`-invariant), the `ω`-floor (pure-compression `ω_t≈0`), reduce-to-P2a/P2b, end-to-end reversal recovery, frame objectivity; **DT5** diagnostic |
| `run_p2d_gate` | the single-step tensor update == path driver; the numerical damaged tangent: degraded+indefinite, non-symmetric, finite across reversal/eigenvalue-crossing |
| `run_p2e_gate` | the **analytic** damaged tangent == numerical (~1e-10) in tension/confined-compression/shear/reversal; reduces to elastic/effective pre-onset; the Macaulay-kink subgradient |

---

# Part III — Usage

> Until the `nDMaterial` wrapper ships, the commands below are the **intended** interface (so models can
> be planned and reviewed against it); §20 is how to exercise the verified core today.

## 17. Command grammar (intended)

```tcl
nDMaterial LadrunoConcrete3D $tag $E $nu $fc $ft  \
    <-rho $rho>                                   \
    <-e $e | -fcc $fccRatio>                      \
    <-Df $Df>                                     \
    <-hardening $qh0 $Hp>                         \
    <-ductility $Ah $Bh $Ch $Dh>                  \
    <-Gf $Gf -Gc $Gc>  <-As $As>                  \
    <-eta $eta> <-implex>
```
`fc`, `ft` are **positive magnitudes**; the model is compression-negative internally. Consumed by any
3D solid element (`LadrunoBrick`, `stdBrick`, `SSPbrick`, …); the planned `-geom finite` lift rides the
existing `LogStrain`/finite path.

## 18. Parameters, defaults & calibration

| Param | Meaning | Default |
|---|---|---|
| `E, nu` | Young's modulus, Poisson ratio | required |
| `fc, ft` | uniaxial compressive / tensile strength (positive) | required |
| `-rho` | mass density | 0 |
| `-e` | deviatoric eccentricity, **∈ (0.5, 1] (hard convexity bound)** | derived from `-fcc` |
| `-fcc` | equibiaxial/uniaxial ratio `fcc/fc` (sets `e`) | 1.16 (Kupfer → `e≈0.52`) |
| `-Df` | dilatancy (non-associated flow); `<1` realistic, `=1` associated | set `<1` (e.g. 0.7) |
| `-hardening qh0 Hp` | initial yield fraction; hardening modulus | 0.3, 0.5 |
| `-ductility Ah Bh Ch Dh` | confinement-ductility (Eq.33) | CDPM2 literature defaults |
| `-Gf -Gc` | tensile / compressive fracture energy (crack-band) | from `fc` |
| `-As` | softening ductility (`x_s` in compression = the confinement-ductility hook) | 2.0 |
| `-eta` | Duvaut–Lions viscosity (0 = inviscid) | 0 |
| `-implex` | engage Tier-2 IMPL-EX | off |

- **`e` is a validation target, not a fit knob** — leave it derived from `-fcc` unless you have biaxial
  data; it lands at the canonical `e ≈ 0.52`.
- **`-ductility` sets the absolute peak strain + its confinement growth** — the defaults give the right
  *trend*; recalibrate `Ah/Bh/Ch/Dh` to measured peak strains. Keep `Hp` small.
- **`-Gf`/`-Gc`** are crack-band-regularized ⇒ damage dissipation is mesh-objective (the FE-visible
  *total* carries the §7 plastic-dissipation caveat).

## 19. Solver requirement — READ THIS

The consistent tangent is **non-symmetric** (non-associated flow **and** the spectral recompose, §5),
unconditionally — even at `Df=1`. You **must** use an unsymmetric solver:

```tcl
system UmfPack          ;# or FullGeneral
```

A symmetric solver (`ProfileSPD`, `BandSPD`, symmetric `Mumps`) **silently mis-solves** — there is no
runtime guard (the wrapper will warn). The IMPL-EX tier (`-implex`) is SPD ⇒ symmetric is fine there.

## 20. Exercising the verified core today (no wrapper)

```bash
# the full constitutive law, all gates (numpy only):
python tests/_testbed/concrete3d_ref.py
# the pytest gate battery + the g++ kernel byte-check:
python -m pytest tests/test_ladrunoConcrete3D_material.py -q
```
```python
# single-point response from the oracle (the exact constitutive update the wrapper will call):
import sys; sys.path.insert(0, "tests/_testbed"); import numpy as np, concrete3d_ref as R
mp = R.make_material(E=30000., nu=0.2, fc=30., ft=3., Df=0.7)
st = R.make_damage_state(mp)
sig, st, info = R.damaged_step_tensor(st, deps6, mp, Gf=0.1, Gc=5.0, lch=50., As=2.0)
C   = R.damaged_tangent_analytic(st, deps6, mp, 0.1, 5.0, 50., 2.0)   # the analytic damaged tangent
```

## 21. What you get at each phase + recorders

- **Now (oracle + kernel through the damage stress update):** the full **monotonic** stress response —
  pre-peak plasticity, the triaxial strength envelope (Kupfer biaxial, confined `fcc(σ3)`),
  confinement-dependent ductility, **the peak + tension/compression softening + automatic unilateral
  recovery**, crack-band-regularized.
- **P3b:** the analytic damaged tangent in the kernel (Tier-1 implicit accuracy on softening).
- **Wrapper:** callable `nDMaterial`, classTag 33017, recorders — stress/strain/tangent + `damage`
  (tension/compression), `kappa_p`, `equivalentPlasticStrain`, `implexError`.
- **P2f / later:** cyclic (`β_c` + the compression→tension temper), multiaxial-damage apportioning,
  robustness tiers (`-implex`, explicit), finite strain (`-geom finite`), the confined-fiber 1D view.

## 22. Units

Strength/modulus in stress units (e.g. MPa), `Gf`/`Gc` as energy/area (e.g. N/mm), `lch` as length
(mm) consistent with the element — `ε_f = Gf/(ft·lch)` must be a small strain, so keep `lch` at
mm-scale. Compression-negative internally; enter `fc`, `ft` positive.

---

# Part IV — Quirks

## 23. Consolidated quirks

1. **The `ω`-solve needs a physical stress FLOOR.** Solve `ω_t` only when `sig_t_drive > 1e-6·ft`
   (mirror `ω_c`/`fc`). In a uniaxial-*stress* state the lateral effective stress is only driven to ~0,
   leaving a ~1e-10 MPa residual whose **sign is numerical noise**; gating on `> 0` lets it (plus the
   tensile history CDPM2 accumulates even in compression) drive `ω_t → 1` in pure compression. It is
   mask-hidden (the compressive axial principal routes through the `(1−ω_c)` channel ⇒ `σ11` unaffected)
   but it poisons the state and is sign-of-noise fragile. **Apply the floor at every site that
   recomputes `ω`** — the stress update *and* the analytic-tangent recompute, or the tangent decouples
   from the update. (Adversarial-review fix.)

2. **Unilateral recovery is AUTOMATIC if you re-split the converged effective stress every step.** Do
   NOT freeze the tension/compression projectors (the ASD IMPL-EX trap) — a closing crack would keep
   the `(1−ω_t)≈0` factor and carry no compression. Re-split is cheap (one symmetric 3×3 eig you
   already do) and tier-independent; freeze only the damage-driver extrapolation.

3. **Drive each damage channel by the EXTREME effective principal of its sign**, not `E·κ_d`. `E·κ_dt`
   is correct for tension (and reduces to the 1D driver) but gives the wrong peak (`ft≠fc`) for
   compression.

4. **Voigt double-contraction weight `[1,1,1,2,2,2]` on analytic tensor gradients, NOT on micro-FD
   gradients.** An eigenprojection `∂λ/∂σ̄ = E_a` contracted with the stress increment is a tensor
   double contraction (off-diagonals ×2); a per-Voigt-component micro-FD is already the per-component
   partial (no weight). Omitting the ×2 is invisible on coaxial paths and gives a 1e-2 tangent error
   the moment shear appears.

5. **The `σ̄_lat=0` Macaulay kink is a valid-subgradient point, not a bug.** In uniaxial-stress
   compression the lateral eigenvalues sit on the tension/compression boundary; the analytic tangent
   agrees with the central difference on the loaded axial component and differs only on the ~zero-stress
   lateral + coupled-shear directions. Gate the analytic tangent at smooth states (eigenvalues bounded
   from 0); build confined-compression test states **stress-controlled** (small `σ3`), never
   confined-*strain* (apex-fragile, platform-divergent — a CI-red incident).

6. **Effective-stress plasticity is MONOTONIC.** The failure surface is reached at `κ_p=1`; the peak +
   softening are the DAMAGE part. Gate "`σ̄` at `κ_p=1`", never "peak == fc".

7. **Crack-band gates that integrate over `ε_i` are tautological** (return `G/lch` by construction) —
   they verify the `ε_f` wiring, not FE objectivity. The honest objectivity number is the D3/C3
   total-energy spread (reported).

8. **Compression pre-damages subsequent tension** (graded by compression damage) — literal CDPM2 Eq.43,
   not a bug; the cyclic temper is P2f. Tracked by the DT5 diagnostic.

9. **Hardening apex honesty (kernel).** A converged plastic hardening return is gated on admissibility
   (`κ_p≥κ_p,n`, `Δλ≥0`, finite, on-surface); deep-compression near-apex trials honestly report
   `converged=false` + elastic fallback rather than teleporting to the tension vertex. The kernel is
   the *safe* reference and deliberately diverges from the oracle's arbitrary apex teleport there.

10. **Every fork pytest file needs module-level `pytestmark = [pytest.mark.zone_a]`** or CI's
    `pytest -m zone_a` silently deselects it (the #249 CI-false-green incident — the whole file ran
    *nothing* while reporting SUCCESS).

11. **`m0` is derived, never user-set.** `yieldF` trusts `mp.m0`; the wrapper must enforce
    `m0 == m0Of(fc,ft,e)` in one place and reject `e ≤ 0.5` (convexity).

---

# Appendices

## A. Phase + PR history + V&V

| Phase | What | PR |
|---|---|---|
| P0 | MW surface fc-normalization (oracle) | [#240](https://github.com/nmorabowen/OpenSees/pull/240) |
| P1 | perfect-plastic return map + hardening + consistent tangent (oracle) | #240, [#244](https://github.com/nmorabowen/OpenSees/pull/244), [#247](https://github.com/nmorabowen/OpenSees/pull/247) |
| P1 review | review fixes + dev handoff | [#248](https://github.com/nmorabowen/OpenSees/pull/248) |
| P1 C++ | kernel return map + analytic effective tangent + g++ dump gate (5-lens review) | [#249](https://github.com/nmorabowen/OpenSees/pull/249) |
| P2a | tensile `ω_t` + crack-band `G_f` (oracle) | [#259](https://github.com/nmorabowen/OpenSees/pull/259) |
| P2b | compressive `ω_c` + `α_c` split + `G_c` (oracle) | [#261](https://github.com/nmorabowen/OpenSees/pull/261) |
| P2c | unified tensor split + automatic unilateral (oracle) | [#284](https://github.com/nmorabowen/OpenSees/pull/284) |
| P2d | single-step update + numerical damaged tangent (oracle) | [#285](https://github.com/nmorabowen/OpenSees/pull/285) |
| P2e | analytic dual-projector damaged tangent (oracle) | [#286](https://github.com/nmorabowen/OpenSees/pull/286) |
| fix | PE2 cross-platform robustness | [#287](https://github.com/nmorabowen/OpenSees/pull/287) |
| review | 4-dim adversarial review fixes (`ω` floor + de-tautologize DT3 + DT5 diagnostic) | [#288](https://github.com/nmorabowen/OpenSees/pull/288) |
| **P3a** | **C++ kernel damage stress update (g++ byte-verified ~1e-14)** | [#289](https://github.com/nmorabowen/OpenSees/pull/289) |
| P3b | analytic damaged tangent in the kernel | *next* |
| wrapper | `nDMaterial LadrunoConcrete3D` + classTag 33017 | *pending* |
| P2f | cyclic `β_c` + compression→tension temper + multiaxial apportioning | *deferred* |

**V&V status:** oracle gates **18/18** (Zone-A); the g++ kernel byte-check (PATH/TAN/DMG + fuzzer)
green on Linux CI. classTag 33017 **RESERVED** (lands with the wrapper; LogStrain2D-33016 convention).
Two adversarial-review rounds (a 5-lens panel on #249, a 4-dimension workflow on P2c-P2e) — math core
held both times; findings folded in (CI false-green, apex honesty, the `ω` floor, the tautological DT3).

## B. References

- Grassl, Xenos, Nyström, Rempling, Gylltoft (2013), *CDPM2: A damage-plasticity approach to modelling
  the failure of concrete*, Int. J. Solids Struct. 50, 3805-3816 (arXiv:1307.6998) — pinned by Eq. №.
- Menétrey & Willam (1995), *Triaxial failure criterion for concrete and its generalization*, ACI SJ.
- de Souza Neto, Perić, Owen (2008), *Computational Methods for Plasticity* — Box A.6 (isotropic
  tensor-function derivative), the spectral consistent tangent.
- Decision record: [[31_ladruno_concrete3d_adr]]. Implementer brief: [[LadrunoConcrete3D_dev_handoff]].
  The shell/MCFT sibling: [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] (33015).
