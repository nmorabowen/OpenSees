---
title: ADR — Lemaitre coupled ductile damage on the J2–Chaboche family
project: Ladruno
status: accepted (decisions locked 2026-06-02)
priority: high
owner: nmora
tags:
  - implementation
  - material
  - plasticity
  - damage
  - fracture
  - adr
---

# ADR — Lemaitre coupled ductile damage (`LadrunoJ2` + `LadrunoUniaxialJ2`)

**Status:** accepted — decisions D1–D4 locked 2026-06-02 · **Parent / spine:** the combined-hardening
J2 family — 3D [[10_ladruno_j2_plasticity]] (`nDMaterial LadrunoJ2`, classTag
33011) and its uniaxial twin [[13_ladruno_uniaxial_j2_adr]]
(`uniaxialMaterial LadrunoUniaxialJ2`, classTag 33000) · **Book frame:** dSNPO
(de Souza Neto–Perić–Owen 2008) **Ch. 12 §12.3** (Lemaitre's ductile damage
coupled to AF-kinematic + isotropic hardening) and **§12.4** (the integration
algorithm) · **Shared code:** `LadrunoHardening.h`.

> One-line verdict: this is the **natural, book-grounded continuation** of exactly
> the AF + isotropic family we already shipped. dSNPO §12.3 *is* "AF-kinematic +
> isotropic hardening + Lemaitre damage D in one return map." It is a **coupled
> core extension, not a wrapper**: damage degrades stiffness *inside* the flow
> (effective stress `σ̃ = σ/(1−D)`), and the damage rate is **triaxiality-driven**
> in 3D — neither of which a strain-watching `Fatigue`/`MinMax` wrapper can
> express. It lights up ductile fracture on **both** the continuum element
> (`LadrunoBrick` + `LadrunoJ2`) and the fiber (`LadrunoUniaxialJ2`) because the
> isotropic backbone (and now the damage law) is shared code.

---

## 1. Context

### 1.1 Where this sits

The fork now has a clean, **undamaged** combined-hardening J2 spine:

- `LadrunoJ2` (3D, classTag 33011): Voce+linear isotropic + Chaboche AF kinematic
  (`α = Σₖ αₖ`), scalar-Newton return map (`n = M(Δγ)/‖M(Δγ)‖`, Kobayashi–Ohno),
  analytic consistent tangent (reduces term-by-term to `J2Plasticity` at N=0), all
  5 dimensional views, full state recording. IMPL-EX is a **structure-only hook**:
  `dGamma_n` is already committed and serialized.
- `LadrunoUniaxialJ2` (1D, classTag 33000): the uniaxial reduction, sharing the
  isotropic law verbatim via `LadrunoHardening.h` (the V7 1e-12 oracle contract).
  Already carries a committed **plastic-work accumulator** `Cwp` and the same
  `dGamma`/IMPL-EX hook.

The J2 ADR explicitly foresaw this extension. Its IMPL-EX note reads: *"It earns
its keep later for explicit dynamics and any future softening/damage extension
(where implicit Newton stalls on indefinite tangents)."* The uniaxial ADR's §6.1
parked a `LadrunoFatigue` **wrapper** for low-cycle fatigue — Lemaitre is the
**mechanistic, coupled** sibling of that idea, and the two are complementary (see
§1.3).

### 1.2 Why Lemaitre is the strong pick (vs the parked wrapper)

| | `LadrunoFatigue` wrapper (parked, [[13_ladruno_uniaxial_j2_adr]] §6.1) | **Lemaitre coupled damage (this ADR)** |
|---|---|---|
| coupling | post-hoc: watches `Δεᵖ`, scales/zeros stress at `D=1` | **in the return map**: `σ̃ = σ/(1−D)` degrades stiffness *and* accelerates plastic flow |
| driver | plastic-strain range / dissipated energy (rainflow) | **damage energy release rate `Y`** — explicitly **triaxiality-dependent** in 3D |
| 3D fracture | no — a wrapper cannot see `σ_H/σ_eq` | **yes** — notch/triaxiality-driven ductile fracture, the headline |
| applies to | uniaxial fibers only | **both** `LadrunoJ2` (3D continuum) and `LadrunoUniaxialJ2` (fiber) — shared law |
| book home | phenomenological (Coffin–Manson / Miner) | **dSNPO §12.3** — sitting right there, AF + iso + Lemaitre |

The continuum element (`LadrunoBrick`) gets ductile fracture **for free** because
the damaged material plugs into the exact same `setTrialStrain`/`getTangent`
contract — no element change.

### 1.3 The two damage families are complementary, not competing

- **Lemaitre (this ADR):** continuum damage mechanics, triaxiality-driven,
  monotonic-ductile-tearing oriented. The *physics of necking / void growth /
  ductile rupture under multiaxial stress*.
- **`LadrunoFatigue` (still parked):** rainflow low-cycle-fatigue *index*, cheap,
  symmetric-cyclic-seismic oriented, wraps the undamaged core.

Both can coexist: Lemaitre is a **core mode** (`-damage lemaitre …`); the fatigue
index stays a **wrapper**. v1 of this ADR ships Lemaitre only.

---

## 2. Decision (locked 2026-06-02)

1. **[D1 — LOCKED] Damage as a *mode* on the existing classes**, not new sibling
   classes. Add `-damage lemaitre $r $s $pD $Dc` to **both** `LadrunoJ2` and
   `LadrunoUniaxialJ2`. **Default `-damage none` ⇒ `D ≡ 0`**, taking the *exact
   current code path* with zero overhead so **every existing test stays
   bit-identical**. *No new classTag* (damage is a parameter, not a registry class).
2. **Effective-stress coupling, strain-equivalence.** `σ = (1−D) Dᵉ : εᵉ`; yield
   evaluated on the effective relative stress `ξ̃ = dev(σ̃) − α`. The AF return-map
   *direction machinery is unchanged* — damage scales the deviator magnitude, not
   its direction, so the `M(Δγ)` radial structure survives (§3).
3. **[D2 — LOCKED] Numerical scheme = fully-coupled backward-Euler.** Solve the
   `(Δγ, D)` system *simultaneously* — `D` is live inside every Newton iterate, not
   frozen. Because damage preserves the radial direction (item 2), this is **the
   same scalar-Newton-on-Δγ skeleton the spine already ships**, augmented to a
   2-unknown `(Δγ, D)` solve (dSNPO §12.4); it is *not* a from-scratch tensor solve
   (§3.4). The operator-split (frozen-`D`) scheme is kept only as a **robust
   fallback / seed**, not the default.
4. **[D3 — LOCKED] 3D + uniaxial in ONE PR.** Both `LadrunoJ2` and
   `LadrunoUniaxialJ2` get the mode together via the shared `LadrunoDamage.h`, so
   the **3D↔1D 1e-12 damage-reduction oracle (V4 test) lands complete from day
   one**. Note the triaxiality caveat: a `uniaxialMaterial` has no triaxiality
   (`σ_H/σ_eq ≡ 1/3`) ⇒ `R_v = const`, so the uniaxial Lemaitre is the valid
   **degenerate/energy-based** point and oracle; the multiaxial *fracture* story
   still lives in 3D — the single PR just makes both consistent at once.
5. **[D4 — LOCKED] Softening regularization = characteristic-length + IMPL-EX.**
   Local Lemaitre with **characteristic-length regularization** of the damage
   evolution, reusing the `getCharacteristicLength()` + `LadrunoBrick` `lch`
   handshake **already built for ASDConcrete3D**
   ([[11_brick_asdconcrete_integration.md]]) so dissipated fracture energy is
   mesh-objective; **plus a real IMPL-EX code path** (turn the existing
   `dGamma_n` hook on — both classes already commit it) as the robustness route for
   the indefinite softening tangent. Honest limit: `lch` makes *dissipated energy*
   objective, **not** the localization width — nonlocal/gradient damage is a future
   ADR.
6. **Failure handling = cap + flag, no erosion.** Clamp `D` at `D_c`, degrade
   stress, expose a `"damage"`/`"failed"` response. Element erosion/deletion is an
   element concern, explicitly out of scope.
7. **No `exit()`, real `revert*`/`setParameter`/`Print`/`sendSelf`** — Ladruno
   house hygiene, same as the spine. State grows by **one committed scalar** `D_n`
   (3D) / `D_n` (1D) plus the four damage params.

---

## 3. Constitutive model — Lemaitre coupled to AF + isotropic (FEM-expert)

Reference: dSNPO **§12.3** (the coupled model), **§12.4** (integration / Box
12.x), Lemaitre & Chaboche (1990), Lemaitre (1985). The point of this section is
that the model is the **existing `LadrunoJ2` with one scalar field `D` threaded
through the effective stress** — the return-map skeleton is reused, not rebuilt.

### 3.1 State, free energy, effective stress

Internal variables: elastic strain `εᵉ`, accumulated plastic strain `p` (≡ `ε̄ᵖ`),
backstresses `αₖ`, **damage `D ∈ [0, D_c]`** (`D=0` virgin, `D_c` ≈ 0.2–0.5 at
rupture). Helmholtz free energy splits elastic (damaged) + plastic (hardening):

```
ρψ = (1−D)·½ εᵉ : Dᵉ : εᵉ  +  ρψₚ(p, αₖ)
```

**Elastic law** (strain-equivalence ⇒ the `(1−D)` multiplies the elastic stress):

- Direct: `σ = (1−D) Dᵉ : εᵉ`,  **effective** `σ̃ = σ/(1−D) = Dᵉ : εᵉ`
- Voigt: `{σ} = (1−D)·[ K 1⊗1 + 2G I_dev ]{εᵉ}`

### 3.2 Yield surface in *effective* stress space

`ξ̃ = dev(σ̃) − α = s/(1−D) − α` (relative effective deviatoric stress):

- Direct: `Φ = ‖ξ̃‖ − √(2/3)·σ_y(p)`,  `σ_y(p) = σ₀ + Q∞(1−e^{−bp}) + H_iso p`
- Index: `Φ = √(ξ̃_{ij} ξ̃_{ij}) − √(2/3)·σ_y`,  `ξ̃_{ij} = s_{ij}/(1−D) − α_{ij}`
- 1D reduction: `Φ = |σ/(1−D) − X| − σ_y(p)`  (the `√(2/3)` divides out, exactly as
  in [[13_ladruno_uniaxial_j2_adr]] §3.1)

### 3.3 Flow, hardening, damage evolution

| | form |
|---|---|
| plastic flow | `ε̇ᵖ = γ̇ ∂Φ/∂σ = (γ̇/(1−D))·N`,  `N = ξ̃/‖ξ̃‖` — the `1/(1−D)` **accelerates** plastic flow as damage grows |
| accumulated | `ṗ = √(2/3)‖ε̇ᵖ‖ = γ̇/(1−D)` |
| isotropic | `σ_y(p)` (Voce+linear, shared `LadrunoHardening.h`) |
| kinematic (AF) | `α̇ₖ = (2/3)Cₖ ε̇ᵖ_dir − γₖ αₖ ṗ` (unchanged form; written in effective space) |
| **damage** | `Ḋ = (Y/r)^s · ṗ`  for `p ≥ p_D`, else `Ḋ = 0` (damage threshold) |

**Damage energy release rate `Y`** (the triaxiality carrier — dSNPO eq 12.x):

- Direct: `Y = ½ εᵉ : Dᵉ : εᵉ = σ̃_eq²/(2E) · R_v`
- **Triaxiality function** `R_v = (2/3)(1+ν) + 3(1−2ν)·(σ_H/σ_eq)²`,
  `σ_H = ⅓ tr σ`, `σ_eq = √(3/2)‖s‖` (von Mises)
- **1D:** `σ_H/σ_eq = 1/3` always ⇒ `R_v = (2/3)(1+ν) + ⅓(1−2ν) = const`, so
  `Y = σ̃²/(2E)` up to that constant — the uniaxial degenerate point (§2.4).

Four damage constants: `r` (energy denominator, sometimes `S`), `s` (exponent),
`p_D` (plastic-strain threshold below which no damage), `D_c` (critical/rupture).

### 3.4 Return map — fully-coupled `(Δγ, D)` backward-Euler (v1, **D2 locked**)

Backward-Euler, elastic predictor / plastic corrector, with `D` **live** in the
solve. The key structural fact: the **effective trial stress is
damage-independent** (`σ̃ᵗʳ = Dᵉ:(ε_{n+1} − εᵖₙ)` — no `D`), and damage scales the
deviator *magnitude* not its *direction*, so the shipped `M(Δγ)` direction machinery
carries over **unchanged**. Damage adds exactly **one extra unknown** `D` and
**one extra equation** (the integrated evolution law).

1. **Trial:** `σ̃ᵗʳ = Dᵉ:(ε_{n+1} − εᵖₙ)`, `ξ̃ᵗʳ = dev(σ̃ᵗʳ) − αₙ`,
   `Φᵗʳ = ‖ξ̃ᵗʳ‖ − √(2/3)σ_y(pₙ)`.
2. If `Φᵗʳ ≤ 0`: elastic, `σ = (1−Dₙ)σ̃ᵗʳ`, `D^alg = (1−Dₙ)Dᵉ`. Done.
3. Else **plastic corrector — 2-unknown Newton in `(Δγ, D)`:**
   - **Direction (unchanged from the spine):** `n = M(Δγ)/‖M(Δγ)‖`,
     `M(Δγ) = s̃ᵗʳ − Σₖ αₖⁿ/(1+√⅔γₖΔγ)` (Kobayashi–Ohno), `s̃ᵗʳ = dev(σ̃ᵗʳ)`.
   - **Effective deviator / backstress** as functions of `Δγ` (exact spine forms);
     plastic-strain increment carries the damage acceleration `Δp = √(2/3)Δγ/(1−D)`.
   - **R₁ — consistency (effective space):** `‖ξ̃(Δγ)‖ − √(2/3)·σ_y(pₙ + Δp) = 0`.
   - **R₂ — integrated damage:** `D − Dₙ − (Y(Δγ)/r)^s · Δp = 0` for `pₙ+Δp ≥ p_D`,
     else `D − Dₙ = 0`; with `Y(Δγ) = σ̃_eq(Δγ)²/(2E)·R_v(σ_H/σ_eq)`.
   - Solve `[R₁,R₂]ᵀ = 0` for `(Δγ, D)` by Newton (analytic 2×2 Jacobian; all
     `∂/∂Δγ` terms are the spine's `dM/dΔγ` extended by the `1/(1−D)` factors,
     `∂/∂D` couples through `Δp` and the degradation). Seed `Δγ` from the radial /
     frozen-`D` (operator-split) estimate, `D` from `Dₙ`; bracket `D ≤ D_c`.
4. **Degrade:** `σ_{n+1} = (1−D_{n+1}) σ̃_{n+1}`, clamp `D_{n+1}` at `D_c`.

> **Why this is not much more than the spine.** It is the *same scalar-Newton-on-Δγ
> skeleton*, promoted to a 2×2 by adding `D` and the damage residual `R₂` (dSNPO
> §12.4). The radial reduction the spine already exploits is exactly what keeps the
> coupled system small — there is **no return to a full tensor solve**. The
> **operator-split** (freeze `D=Dₙ` over the return, then update `D` explicitly) is
> retained only as the **Newton seed and a robustness fallback** when the coupled
> 2×2 fails to converge near rupture, not as the default path.

### 3.5 Consistent tangent

Differentiate the **converged coupled `(Δγ, D)` state** w.r.t. `ε_{n+1}`. The
algorithmic tangent is the degraded elastoplastic operator plus the
damage-sensitivity rank-one term:

```
D^alg = (1−D_{n+1})·D^alg_plastic  −  σ̃_{n+1} ⊗ ∂D_{n+1}/∂ε_{n+1}
```

where `∂D_{n+1}/∂ε` and `∂Δγ/∂ε` come from inverting the 2×2 Newton Jacobian of
§3.4 (the implicit-function theorem on `[R₁,R₂]=0`) — so the tangent is **consistent
with the fully-coupled solve**, not a frozen-`D` approximation. **Oracle checks:**
with damage off (`r→∞` or `p_D→∞` ⇒ `R₂` collapses to `D=Dₙ=0`) this reduces
**bit-identically** to the shipped `LadrunoJ2` tangent (the cheapest regression,
V0); FD-check `‖D^alg − D^FD‖` at several damaged states (V5).

---

## 4. Architecture — reuse the spine, add one scalar field

### 4.1 Shared damage law header (mirror the `LadrunoHardening.h` pattern)

Factor the damage kinematics into a tiny header-only **`LadrunoDamage.h`** next to
`LadrunoHardening.h`:

- `damageReleaseRate(sigEff_eq, sigH, nu, E) -> Y` (with the `R_v` triaxiality
  function — 3D passes the real ratio, 1D passes the frozen `1/3`).
- `damageIncrement(Y, dp, r, s, pD, p) -> ΔD` and its `∂(ΔD)/∂(·)` derivatives.

Both `LadrunoJ2.cpp` and `LadrunoUniaxialJ2.cpp` `#include` it, so the 3D↔1D
**damage** reduction is byte-identical at triaxiality 1/3 — extending the existing
oracle contract from "isotropic backbone" to "isotropic backbone + damage."

### 4.2 Touch list (all additive; mode is OFF by default)

- **`SRC/material/nD/LadrunoJ2.{h,cpp}`** — add `damageModel` enum (default `NONE`),
  `r,s,pD,Dc` params, committed `D_n` + trial `D`, the §3.4 branch in `integrate()`
  (guarded: `NONE` ⇒ current path untouched), `R_v` from the (already-computed)
  stress, `D` into `sendSelf`/`recvSelf`/`revert*`, `setResponse("damage")`.
- **`SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}`** — same, with `R_v=const`;
  reuse the existing `Cwp` plastic-work accumulator for the energy-based diagnostic.
- **`SRC/material/nD/LadrunoDamage.h`** — new shared header (above).
- **`OPS_LadrunoJ2.cpp` / `OPS_LadrunoUniaxialJ2.cpp`** — parse `-damage lemaitre
  $r $s $pD $Dc` (+ optional `-Dc`, `-implex`); default `none`.
- **Vanilla (→ [[LEDGER_vanilla_files]]):** **no `classTags.h` change** — D1 locked
  to mode-not-class, so damage is a parameter on the existing registered classes; the
  registry/broker entries are already in place. (The rejected sibling-class route
  would have needed `ND_TAG_LadrunoJ2Damage 33012` + `MAT_TAG_…Damage 33001` + broker
  + registry — that churn is exactly what the mode decision avoids.)
- **Build-control:** `LEDGER_implementations.md` (annotate the two existing
  `LadrunoJ2`/`LadrunoUniaxialJ2` rows with the damage mode), `LEDGER_quirks.md`
  (softening/`lch`/IMPL-EX gotcha + the localization-width disclaimer), banner line
  via `Ladruno_scripts/banner_features.txt` → `patch_banner.py`.

### 4.3 IMPL-EX (turn the existing hook on)

Both classes already commit `dGamma_n`. The IMPL-EX path: extrapolate
`Δγ̃ = Δγₙ·(Δtₙ₊₁/Δtₙ)`, freeze the flow direction and `D` from the extrapolation,
and form a **constant SPD tangent** — exactly the robustness escape for the
indefinite softening tangent. Ship it as `-implex` here (its first real consumer),
not just a structural hook. Pairs with explicit dynamics (`CentralDifferenceLadruno`
/ `ExplicitBathe`) where softening + implicit Newton would otherwise stall.

---

## 5. Public API (signature draft)

```tcl
# 3D: combined hardening + Lemaitre ductile damage
nDMaterial LadrunoJ2 $tag  $K $G  -iso voce $sig0 $Qinf $b $Hiso \
                                  -kin $N $C1 $g1 $C2 $g2 $C3 $g3 \
                                  -damage lemaitre $r $s $pD $Dc <-implex>

# uniaxial fiber: same, triaxiality fixed at 1/3
uniaxialMaterial LadrunoUniaxialJ2 $tag  $E  -iso voce $sig0 $Qinf $b $Hiso \
                                        -kin 1 $C1 $g1  -damage lemaitre $r $s $pD $Dc
```

- **`-damage none` (default, or simply omitted) ⇒ identical to today.**
- Responses add `"damage"` (scalar `D`), `"failed"` (bool `D≥D_c`), and the energy
  release rate `"Y"`; the spine's stress/strain/tangent/backStress/plasticStrain/
  equivalentPlasticStrain are unchanged. The uniaxial `"plasticWork"` (`Cwp`) stays.

---

## 6. Testing / oracle matrix (Zone-A pytest)

(Test IDs are `V*` to avoid colliding with the decision labels `D1`–`D4` in §2/§7.)

| # | Test | Reference / oracle |
|---|---|---|
| V0 | **damage-off regression**: `-damage none` (or omitted) | **bit-identical** to shipped `LadrunoJ2`/`LadrunoUniaxialJ2` batteries (9/9 each) — the cheap tripwire |
| V1 | uniaxial monotone tension, damage on | semi-analytic Lemaitre coupon `σ(ε)` with softening onset at `p_D` |
| V2 | damage-off via `p_D = ∞` mid-run | proves the threshold gate + that pre-threshold response == undamaged |
| V3 | **triaxiality sweep (3D)**: uniaxial vs notched/confined at fixed `p` | damage-rate ratio matches `R_v(σ_H/σ_eq)` closed form — THE differentiator |
| V4 | **3D↔1D damage reduction** at triaxiality 1/3, identical `(r,s,pD,Dc)` | `LadrunoJ2`-condensed-to-uniaxial vs `LadrunoUniaxialJ2` ~1e-12 (the shared `LadrunoDamage.h` contract; extends the spine's V7 oracle) |
| V5 | consistent-tangent FD at several damaged states | `‖D^alg − D^FD‖ < tol` |
| V6 | **mesh-objectivity**: notched bar, 2 mesh sizes, `LadrunoBrick`+`LadrunoJ2` | same global force–displacement with `lch` regularization (mirrors the ASDConcrete size-objectivity test, [[11_brick_asdconcrete_integration.md]]) |
| V7 | dissipated energy at rupture (`D→D_c`) | == prescribed fracture energy `G_f` (mesh-objective via `lch`) |
| V8 | IMPL-EX vs implicit on a softening step | bounded lag; SPD tangent; no Newton stall |
| V9 | `sendSelf`/`recvSelf` round-trip incl. `D` | parallel/serialization hygiene |

Smoke (L0): single `LadrunoBrick`+`LadrunoJ2` and single truss/fiber
`LadrunoUniaxialJ2` push-to-rupture, Zone-A CPython 3.12 env
([[../.claude/.../memory/project_opensees_test_env]]).

---

## 7. Risks / open questions

> [!note] **D1 — RESOLVED (mode, not class).** `-damage lemaitre …` on the existing
> classes, default off. No classTag churn, no duplication, current tests
> bit-identical (V0). Residual risk: `integrate()` branches and softening code share
> a TU with the pristine spine — contained by the V0 bit-identical gate.

> [!note] **D2 — RESOLVED (fully-coupled now).** `(Δγ, D)` solved simultaneously
> (dSNPO §12.4), `D` live in each iterate. Feasible without a tensor solve because
> damage preserves the radial `M(Δγ)` direction (§3.4) — it is the spine's scalar
> Newton promoted to a 2×2. **Residual risk:** the coupled 2×2 can stiffen near
> rupture (`D→D_c`, `Y` large); mitigations = operator-split seed, bracket `D≤D_c`,
> and the IMPL-EX escape (D4). The frozen-`D` operator-split stays in the code as
> the fallback, so we keep both robustness paths.

> [!note] **D3 — RESOLVED (3D + uniaxial in one PR).** Both ship together on the
> shared `LadrunoDamage.h`; the V4 1e-12 reduction oracle is complete from day one.
> **Caveat to keep visible:** the uniaxial twin is the *degenerate* (`R_v=const`,
> triaxiality 1/3) point — do not market it as multiaxial fracture; that lives in
> 3D. Larger single PR is the accepted cost.

> [!note] **D4 — RESOLVED (char-length + IMPL-EX).** Local Lemaitre +
> `lch`-regularized damage (mesh-objective `G_f`, reusing the ASDConcrete handshake)
> + a real IMPL-EX path (first consumer of the committed `dGamma_n` hook).
> **Honest disclaimer (must go in the banner/docs, not be silently dropped):** `lch`
> makes the **dissipated energy** mesh-objective but **not** the localization width;
> true width-objectivity needs nonlocal/gradient damage (future ADR).

> [!note] **Triaxiality needs the stress state, not richer strain.** `R_v` is
> computed from `σ_H/σ_eq` of the *converged* stress — already available in both
> classes' return maps. No interface change (the `setTrialStrain` contract stays
> scalar/tensor as-is), mirroring the buckling-length lesson in
> [[13_ladruno_uniaxial_j2_adr]] §6.2.

> [!note] **Backstress-in-effective-space variants.** Lemaitre–Chaboche has a few
> formulations of how `αₖ` couples to `D` (effective vs nominal). v1 follows dSNPO
> §12.3 (yield on `ξ̃ = s/(1−D) − α`, `α` in effective space). Pin the choice in
> code comments; it is the difference that the D4 1e-12 reduction will catch if
> 3D and 1D ever diverge.

- **No `exit()`** anywhere; real `getInitialTangent` (`=(1−Dₙ)Dᵉ` after damage, but
  the *initial* tangent for modal/eigen is the **virgin** `Dᵉ` — `D=0` at start);
  real `revert*`/`setParameter`/`Print`. Same hygiene as the spine.
- **Backwards compat:** additive, default-off; no existing model affected; V0 is
  the proof.
- **Oracle drift tripwire:** the D4 1e-12 gate fails the instant a 3D damage edit
  bypasses the shared `LadrunoDamage.h`.

---

## 8. Implementation plan (step order — one PR, 3D + uniaxial together)

1. **`LadrunoDamage.h`** (shared, header-only, no OpenSees deps — mirror
   `LadrunoHardening.h`): `Y(σ̃_eq, σ_H/σ_eq, ν, E)` with the `R_v` triaxiality
   function, `ΔD(Y, Δp, r, s, p_D, p)` + the `∂/∂Δγ`, `∂/∂D` derivatives the §3.4
   2×2 Jacobian needs.
2. **3D `LadrunoJ2`**: `damageModel` enum (default `NONE`), `r,s,p_D,D_c`, committed
   `D_n`, the **§3.4 fully-coupled `(Δγ,D)` Newton** in `integrate()` (guarded so
   `NONE` is the *exact* current scalar-Newton path — operator-split kept as the
   seed/fallback), §3.5 consistent tangent, `sendSelf`/`recvSelf`/`revert*`,
   `setResponse("damage"|"failed"|"Y")`. **Run the existing 9/9 battery — must be
   byte/behaviour-unchanged (V0)** before adding damage behaviour.
3. **Uniaxial `LadrunoUniaxialJ2`** (same PR): same mode via `LadrunoDamage.h` with
   `R_v=const` (triaxiality 1/3); reuse the existing `Cwp` plastic-work accumulator
   as the energy diagnostic. **Run the existing 9/9 — unchanged (V0).**
4. **Parsers** `OPS_LadrunoJ2.cpp` / `OPS_LadrunoUniaxialJ2.cpp`:
   `-damage lemaitre r s pD Dc <-implex>`, default `none`.
5. **Regularization:** wire `getCharacteristicLength()` into the damage evolution
   (scale the post-threshold softening so `G_f` is mesh-objective); reuse the
   `LadrunoBrick` `lch` handshake.
6. **IMPL-EX code path** (`-implex`): extrapolate `Δγ`, freeze `D`/direction,
   constant SPD tangent — first real consumer of the committed `dGamma_n` hook; the
   softening-robustness escape.
7. **Zone-A battery** §6 (V0–V9), with the **V3 triaxiality sweep**, **V4 3D↔1D
   1e-12 reduction**, and **V6 mesh-objectivity** as the load-bearing tests (they
   prove this is real triaxiality-driven fracture, not a tuning knob).
8. **Build-control:** `LEDGER_implementations.md` (note the mode on the two existing
   rows), `LEDGER_quirks.md` (softening/`lch`/IMPL-EX gotcha + the localization-width
   disclaimer), banner line via `banner_features.txt` → `patch_banner.py`,
   `// Ladruno` comments on any vanilla touch. **PR off `ladruno`** (`--base
   ladruno`, one logical PR, verify OPEN before any follow-up push — the
   stranded-commit trap, [[../.claude/.../memory/feedback_stranded_commits_after_automerge]]).

**Deferred to follow-up ADRs:** nonlocal / gradient-enhanced damage (rigorous
*width*-objective regularization); Gurson void-growth (dSNPO §12.5 — a *different*,
pressure-dependent surface, not this family); anisotropic / tensorial damage; the
`LadrunoFatigue` rainflow wrapper (still parked, complementary to this).
