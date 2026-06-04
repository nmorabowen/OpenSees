---
title: ADR — Uniaxial combined-hardening J2 steel material (LadrunoUniaxialJ2 / "LadrunoJ2 1D")
project: Ladruno
status: proposed
priority: medium
owner: nmora
tags:
  - implementation
  - material
  - uniaxial
  - plasticity
  - steel
  - adr
---

# ADR — Uniaxial combined-hardening J2 (`LadrunoUniaxialJ2`)

**Status:** proposed (2026-06-02) · **Parent / oracle:** the 3D
[[10_ladruno_j2_plasticity]] (`nDMaterial LadrunoJ2`, classTag 33011) ·
**Registry:** `UniaxialMaterial` (NOT `nDMaterial`) · **Consumes:** fiber
sections, truss/zeroLength, any `uniaxialMaterial` slot.

> One-line verdict from the scoping discussion: **do not ship this as a "better
> Steel02".** That pitch loses to two mature incumbents (`Steel02`,
> `Steel4`). Ship it as (1) the **verification oracle** that pins the 3D
> LadrunoJ2 Chaboche calibration under uniaxial stress, (2) a fiber material
> with **true multi-backstress Chaboche ratcheting** — the one thing
> Menegotto–Pinto cannot do — and (3) the clean **core** onto which the genuinely
> novel steel features (fracture / low-cycle-fatigue index, local-buckling-aware
> degradation, weld/HAZ zones) bolt later as *separable layers*.

---

## 1. Context

### 1.1 The incumbents this must respect, not fight

| material | model | cyclic mechanism | ratcheting | cost | role |
|---|---|---|---|---|---|
| `Steel02` (tag 26) | Giuffré–Menegotto–Pinto | smooth `R(R0,cR1,cR2)` asymptote transition → **Bauschinger by geometry** | **no** | explicit, free | the seismic-steel default |
| `Steel4` (tag 87) | M–P + separated iso/kin, asym, ult, fracture | iso + nonlinear-kin (M–P-style) | weak/approx | explicit-ish | the "advanced" uniaxial slot |
| **`LadrunoUniaxialJ2`** (this) | uniaxial von Mises + Voce/tabulated iso + **Chaboche AF (N backstresses)** | mechanistic backstress evolution | **yes (multi-AF)** | closed-form / scalar-Newton return map | oracle + ratcheting + steel-feature substrate |

Earthquake demand on structural steel is broadly **symmetric low-cycle**:
Bauschinger + low-cycle fatigue dominate, and `Steel02`/`Steel4` already cover
that, calibrated to death. The features that survive a uniaxial reduction of J2
— **ratcheting, mean-stress relaxation** — are *not* what governs seismic steel
fibers; they govern pressure vessels, piping, rails, rolling contact, asymmetric
mechanical demand. So a standalone "uniaxial J2 beats Steel02" claim does not
hold and is **explicitly not** the goal here.

### 1.2 The two things that *do* justify building it

1. **Verification oracle for the 3D model.** Under a uniaxial stress state the 3D
   `LadrunoJ2` (its `PlaneStress`/condensed path driven to a 1-D stress) must
   reproduce the closed-form uniaxial law **to ~1e-12**, *using the identical
   `(C_k, γ_k)` and Voce/tabulated parameters*. That cross-check pins the 3D
   Chaboche calibration and the `(2/3)` bookkeeping far more sharply than any 3D
   self-test can — a uniaxial coupon has an analytic answer; a 3D return map does
   not. This is the cheapest high-value regression we can own.

2. **True Chaboche ratcheting in a fiber.** Multi-backstress AF with distinct
   recovery rates `γ_k` is *the* classical ratcheting model. No OpenSees uniaxial
   material reproduces progressive strain accumulation under asymmetric cycling
   mechanistically. Niche, but real (LCF studies, ratcheting-driven demand).

### 1.3 Future steel features — kept in mind, designed-for, NOT in v1

The exciting roadmap (parked here so the v1 architecture doesn't foreclose it):

- **Fracture / low-cycle-fatigue index.** A cumulative damage `D ∈ [0,1]` driven
  by plastic-strain range (Coffin–Manson / Miner) or dissipated energy, that
  degrades stress and/or triggers fiber rupture. Precedent: OpenSees `Fatigue`
  and `MinMax` are *wrapper* `UniaxialMaterial`s — strong argument to implement
  this as **composition over the J2 core**, keeping the core a pure oracle.
- **Local-buckling-aware fibers.** Asymmetric compression-side strength/stiffness
  degradation (à la Dhakal–Maekawa rebar buckling, or `Steel4`'s fracture/ult
  envelope). A compression-degradation layer — flag or wrapper, not core math.
- **Weld / HAZ behaviour.** Largely a **parameter-set** story (reduced `σ_y0`,
  reduced toughness/ductility per fiber) plus a toughness-driven fracture trigger
  — reuses the fracture layer, no new constitutive core.

**Design consequence:** the v1 core stays a clean, undamaged, energy-consistent
combined-hardening law. Damage/buckling/weld enter as **separable layers**
(wrapper or explicit flag), so the oracle property (§1.2.1) is never polluted by
a degradation term.

---

## 2. Decision

1. **Implement as a `UniaxialMaterial`** named **`LadrunoUniaxialJ2`**
   (`uniaxialMaterial LadrunoUniaxialJ2 ...`). Fiber sections, trusses, and
   zeroLength consume `UniaxialMaterial` — an `nDMaterial` reduction cannot feed
   a fiber. (Naming: see §2a.)
2. **Constitutive content = the exact uniaxial reduction of the 3D `LadrunoJ2`**:
   `σ = E(ε−εᵖ)`, yield `|σ−X| − σ_y(ε̄ᵖ)`, associative flow, Voce+linear *or*
   tabulated isotropic, Chaboche superposed AF backstress `X = Σ_{k=1}^N X_k`
   (ship `N=3`, arch for arbitrary `N`). Same hardening laws, same `(C_k, γ_k)`,
   same `σ_y(·)` as the 3D class — **shared code, §4**.
3. **Return map = closed-form direction + scalar consistency solve.** In 1D the
   flow direction `s = sign(σ_tr − X_n)` is *fixed by the trial state* (relative
   stress cannot rotate on a line), so there is **no `M(Δγ)` direction iteration**
   like the 3D kernel. Consistency reduces to **one scalar equation in `Δγ`**:
   closed-form for linear iso, a short scalar Newton for Voce/tabulated. Cheaper
   and more robust than both the 3D kernel and (arguably) `Steel02`'s implicit
   `R`-curvature solve.
4. **Analytic consistent tangent** `E_alg = E·h/(E+h)`, `h = H'_iso + Σ_k H^{kin}_k`
   (§3.4). `getInitialTangent() = E` (elastic — correct for modal/eigen, the bug
   `SimplifiedJ2` has). No `exit()`, real `revert*`, real `setParameter`/`Print`
   — Ladruno house hygiene.
5. **classTag `MAT_TAG_LadrunoUniaxialJ2 = 33000`** (first Ladruno *uniaxial* material;
   the 33000-band is per-registry, so this does not collide with ELE/ND/INTEGRATOR
   33000). `// Ladruno` comment in `classTags.h`, broker + registry entries.
6. **v1 ships the pure core only.** No damage, no buckling, no weld layer — those
   are follow-up ADRs (§6) built as separable layers.

### 2a. Sub-decision — the name (LOCKED 2026-06-02)

**Chosen: `LadrunoUniaxialJ2`** (user call, 2026-06-02). Lineage-explicit: the
J2/oracle parentage is visible in the command string itself, which suits a
material whose *raison d'être* is being the uniaxial twin of `nDMaterial
LadrunoJ2`. It avoids an exact collision (the `uniaxialMaterial` /`nDMaterial`
prefixes + the `Uniaxial` qualifier disambiguate), at the cost of a longer name.
The steel-feature roadmap (§6) layers on top of this core regardless of its name.
*(Rejected `LadrunoUniaxialJ2` — cleaner product name but hides the J2 lineage that is
the whole point here.)*

---

## 3. Constitutive model — the uniaxial reduction (FEM-expert section)

Reference frame: dSNPO (de Souza Neto–Perić–Owen) **Ch. 3 §3.5 / Box 3.1** (1-D
elastoplasticity), **§6.6.4** (Armstrong–Frederick), Simo–Hughes **Box 1.5**,
Ibrahimbegovic **Ch. 3**. The point of this section is to show that the 3-D von
Mises model of [[10_ladruno_j2_plasticity]] *collapses exactly* to the classical
1-D combined-hardening model — which is **why** the oracle test can hit 1e-12 and
**why** the constants `(C_k, γ_k, σ_y)` carry over unchanged.

### 3.1 From 3D von Mises to 1D — the reduction

Take a uniaxial stress state `σ = σ·(e₁⊗e₁)` (only `σ₁₁ ≠ 0`). The deviator is

- **Direct:** `s = dev(σ) = σ·(e₁⊗e₁ − ⅓𝟏)`
- **Index:**  `s_{ij} = σ(δ_{i1}δ_{j1} − ⅓δ_{ij})`
- **Voigt {11,22,33,12,23,13}:** `{s} = σ·{⅔, −⅓, −⅓, 0,0,0}ᵀ`

so `‖s‖ = √(s:s) = |σ|·√(⅔)`. In a uniaxial *process* the backstress inherits the
same structure, `α = X·(e₁⊗e₁ − ⅓𝟏)`, with `X` the **axial backstress** (units of
stress). The relative (shifted) stress is then collinear:

```
ξ = s − α = (σ − X)·(e₁⊗e₁ − ⅓𝟏)        ⇒    ‖ξ‖ = |σ − X|·√(⅔)
```

Substitute into the 3-D yield `f = ‖ξ‖ − √(⅔)·σ_y(ε̄ᵖ)`:

```
f = √(⅔)·( |σ − X| − σ_y(ε̄ᵖ) )   =  0   ⟺   |σ − X| = σ_y(ε̄ᵖ)
```

The `√(⅔)` is common and divides out. **The 3-D model becomes the textbook 1-D
combined-hardening model**, with `E = 9KG/(3K+G)` the uniaxial Young's modulus
recovered from the bulk/shear pair the 3-D class stores.

### 3.2 Governing equations (1D)

| | Direct / scalar form |
|---|---|
| elastic law | `σ = E·(ε − εᵖ)` |
| yield | `f(σ,X,ε̄ᵖ) = |σ − X| − σ_y(ε̄ᵖ)`,  `σ_y = σ₀ + Q∞(1−e^{−bε̄ᵖ}) + H·ε̄ᵖ`  (or tabulated) |
| flow | `ε̇ᵖ = γ̇·s`,  `s = sign(σ − X)` |
| accumulated | `ε̄̇ᵖ = |ε̇ᵖ| = γ̇` |
| backstress (Chaboche AF) | `Ẋ = Σ_k Ẋ_k`,  `Ẋ_k = C_k·ε̇ᵖ − γ_k·X_k·ε̄̇ᵖ = (C_k·s − γ_k·X_k)·γ̇` |
| KKT | `γ̇ ≥ 0`, `f ≤ 0`, `γ̇·f = 0` |

Each AF term **self-saturates** to `X_{k,sat} = C_k/γ_k` (monotonic) and reverses
on unloading (the `−γ_k X_k` *recall* term) — this is the dynamic recovery that
yields Bauschinger, stable hysteresis, and ratcheting. `γ_k = 0` ⇒ linear Prager
term (`= SimplifiedJ2`); `N=1` ⇒ single AF; `N≥2` ⇒ full Chaboche. The `(C_k,γ_k)`
are **identical** to the 3-D class — that is the shared-hardening contract (§4).

### 3.3 Return map — closed-form direction, scalar consistency

Backward-Euler, trial-elastic / plastic-corrector. Trial freezes plastic state:
`σ_tr = E(ε_{n+1} − εᵖ_n)`, `X_{k} = X_{k,n}`, `η_tr = σ_tr − X_n`,
`f_tr = |η_tr| − σ_y(ε̄ᵖ_n)`.

- If `f_tr ≤ 0`: **elastic**, `σ_{n+1}=σ_tr`, `E_alg=E`. Done.
- Else **plastic**. The corrector direction is **fixed**: `s = sign(η_tr)` (the
  relative stress slides toward `X` along the axis — it cannot change sign). With
  `Δγ ≥ 0`:

```
σ_{n+1}   = σ_tr − E·Δγ·s
X_{k,n+1} = ( X_{k,n} + C_k·Δγ·s ) / ( 1 + γ_k·Δγ )          (per-term closed form)
```

Consistency `s·(σ_{n+1} − X_{n+1}) = σ_y(ε̄ᵖ_n + Δγ)` is **one scalar equation**:

```
R(Δγ) = |η_tr| − E·Δγ − Σ_k [ (s·X_{k,n} + C_k·Δγ)/(1+γ_k·Δγ) − s·X_{k,n} ]
                       − σ_y(ε̄ᵖ_n + Δγ)  = 0
```

`dR/dΔγ` is analytic (`d/dΔγ` of each rational term + `−H'_iso − E`). Solve by
Newton seeded at the perfectly-plastic estimate `Δγ₀ = f_tr/(E + h_0)`. For
**linear** iso + linear kin, `R` is rational and a couple of Newton steps are
exact; for **Voce/tabulated**, the scalar Newton converges in 2–4 steps. Crucially
there is **no tensor `M(Δγ)` direction solve** (the 3-D kernel's complication) —
`s` is constant — so this is strictly simpler and more robust than the 3-D path,
and than `Steel02`'s implicit `R`-curvature evaluation.

### 3.4 Consistent (algorithmic) tangent

Differentiate the converged consistency condition. The result is the classical
elastoplastic tangent

```
E_alg = E·h / (E + h),     h = H'_iso(ε̄ᵖ_{n+1}) + Σ_k H^{kin}_k
H'_iso = b·Q∞·e^{−b ε̄ᵖ} + H            (or the tabulated-spline slope)
H^{kin}_k = ( C_k − γ_k·s·X_{k,n+1} ) / ( 1 + γ_k·Δγ )      (AF effective modulus)
```

`h` is the total plastic hardening modulus (iso + all backstresses). At `Δγ→0`
this is the continuum elastoplastic tangent; for `H^{kin}_k→C_k`, `γ_k→0` it
reduces to linear hardening `E_alg = E·h/(E+h)`. **Oracle check:** `N=0` (pure
Voce iso) must reproduce the uniaxial `J2Plasticity` tangent, and the full model
must match the 3-D `LadrunoJ2` condensed-to-uniaxial tangent — both to ~1e-12 on
the *closed-form* law (no `γ*` fudge here, unlike the 3-D internal `(1−1e-8)`).

### 3.5 Implementation sketch

```cpp
// committed state per material point
struct LadrunoUniaxialJ2State {
  double epsP;        // plastic strain
  double ebarP;       // accumulated plastic strain
  double X[NMAX];     // N backstresses (ship N=3)
  double dGamma;      // last Δγ — IMPL-EX hook only, see 10_ladruno_j2 §IMPL-EX
};

// setTrialStrain(eps):
sigTr = E*(eps - st.epsP);
eta   = sigTr - sum(st.X);
ftr   = fabs(eta) - sigmaY(st.ebarP);
if (ftr <= 0) { sig = sigTr; Etan = E; return; }   // elastic
s = (eta >= 0) ? 1.0 : -1.0;
dG = ftr / (E + h0);                                // seed
for (it=0; it<maxit; ++it) {                        // SCALAR Newton on Δγ
   // assemble R(dG), dRdG from the rational backstress terms + sigmaY(ebarP+dG)
   dG -= R/dRdG;  if (fabs(R) < tol*sigma_scale) break;
}
sig    = sigTr - E*dG*s;
ebarP += dG;  epsP += dG*s;
for k: X[k] = (X[k] + C[k]*dG*s)/(1+g[k]*dG);
Etan   = E*h/(E+h);                                 // h from converged state
```

---

## 4. Architecture — code sharing with the 3D class (the oracle contract)

The oracle property is only real if the **hardening physics is literally the same
code**. Decision:

- Factor the hardening laws out of the 3-D `LadrunoJ2Kernel.h` into a tiny
  **header-only `LadrunoHardening.h`**:
  - `sigmaY(ebarP, isoParams)` + `dSigmaY/dEbarP` — Voce+linear and tabulated
    (C¹ monotone spline) in one place.
  - per-term AF backstress update `(X_n + C·Δγ·dir)/(1+γ·Δγ)` and its `H^{kin}`
    derivative, written generically (scalar `dir`/`X` for 1-D; the 3-D kernel
    feeds the same routine its component values).
- The **return-map structure differs and stays separate**: 3-D uses the
  `M(Δγ)`-direction scalar Newton (Kobayashi–Ohno); 1-D uses the fixed-`s` scalar
  Newton above. Do **not** force the 1-D law through the 6-component tensor
  machinery — that would discard the closed-form-direction win and *create* drift
  risk instead of removing it.
- `LadrunoUniaxialJ2.{h,cpp}` (UniaxialMaterial) `#include`s `LadrunoHardening.h`. The
  3-D class is refactored to consume the same header in the *same* PR (or a
  precursor) so the two cannot diverge.

Files:
- **New:** `SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}`,
  `SRC/material/nD/ladrunoJ2/LadrunoHardening.h` (shared),
  `OPS_LadrunoUniaxialJ2.cpp` parser, `CMakeLists.txt` entry.
- **Modify (vanilla → [[LEDGER_vanilla_files]]):** `SRC/classTags.h`
  (`MAT_TAG_LadrunoUniaxialJ2 33000`), `FEM_ObjectBrokerAllClasses` (broker entry),
  the uniaxial `OPS_*` registry, `SRC/material/uniaxial/CMakeLists.txt`.
- **Reference:** `Steel02.cpp`/`Steel4.cpp` (UniaxialMaterial idioms,
  `sendSelf`/`recvSelf`, fiber consumption), `J2Plasticity.cpp` (the `N=0` oracle
  target), `Fatigue.cpp`/`MinMaxMaterial.cpp` (the wrapper pattern for §6).

---

## 5. Public API

```tcl
# analytic Voce+linear iso, N Chaboche backstresses
uniaxialMaterial LadrunoUniaxialJ2 $tag  $E  -iso voce $sig0 $Qinf $b $Hiso \
                                        -kin $N  $C1 $g1  $C2 $g2  $C3 $g3
# tabulated isotropic curve, single AF
uniaxialMaterial LadrunoUniaxialJ2 $tag  $E  -iso table {e0 s0  e1 s1 ...} \
                                        -kin 1  $C1 $g1
# pure isotropic (no Bauschinger) ⇒ -kin 0 ;  linear Prager ⇒ g_k = 0
```

- `E` directly (uniaxial), *not* `K,G` — a fiber material is parameterized in
  Young's modulus; the oracle mapping `E=9KG/(3K+G)` is documented for the
  cross-check, not exposed.
- `getInitialTangent()=E`; `getTangent()=E_alg`; `getStress`/`getStrain` standard.

---

## 6. Future-feature roadmap (separate ADRs — design hooks only in v1)

Kept explicitly so v1 doesn't foreclose them. **None are v1.**

1. **`LadrunoFatigue` / fracture-LCF index.** Prefer a **wrapper**
   `UniaxialMaterial` (compose over `LadrunoUniaxialJ2`, mirror OpenSees `Fatigue`):
   accumulate `D` from plastic-strain half-cycles (rainflow / Coffin–Manson
   `D += Δε_p^{...}/ε_f`) or dissipated energy; at `D=1` zero the stress
   (rupture). Keeps the J2 core a pure oracle. *Hook needed in v1:* expose
   `ε̄ᵖ`, per-step `Δεᵖ`, and dissipated plastic work via `getResponse`
   ("plasticStrain", "plasticWork") so a wrapper can read them.
2. **Local-buckling-aware degradation.** Compression-side strength/stiffness
   knockdown (Dhakal–Maekawa / Gomes–Appleton style — the OpenSees `ReinforcingSteel`
   `-DMBuckling`/`-GABuckling` precedent). **Core flag, NOT a wrapper** — buckling
   *rewrites the compression-side return* asymmetrically (a buckled bar is not
   symmetric on the next tension excursion), so a pure strain-watching wrapper
   (the §6.1 fatigue pattern) cannot express it; it must see the *signed* stress
   state, which the J2 core carries.
   - **Buckling needs geometry, not richer strain.** Strain is dimensionless;
     buckling scales with slenderness `λ = L/D` (Euler `∝ 1/λ²`). The `setTrialStrain`
     interface stays **scalar** — the buckling length enters as a **parameter**:
     user-given `λ` at construction, *or* pushed via `setParameter` from the
     section/element (which knows the physical tie spacing). This is exactly the
     `ReinforcingSteel` `lsr = unsupported length / bar diameter` channel. **Do
     not widen the strain interface.**
   - **v1 forward-compat hooks (cheap, do now):** (a) design the `setParameter`
     dispatch so a buckling-length/slenderness key can be added without reshaping
     it (the core already ships real `setParameter` for `E,σ₀,Q∞,b,Cₖ,γₖ`); (b)
     keep signed stress + plastic-strain history reachable (already true). No
     buckling code in v1 — just don't foreclose it.
   - **Scope boundary:** this is the *phenomenological* route (a uniaxial σ–ε
     envelope that knows its `λ`). The *geometric* route — discretizing the bar as
     a laterally-deflecting beam-column sub-element — is mechanistic but is **not a
     uniaxial material**; explicitly out of this family.
3. **Weld/HAZ.** Mostly a parameter-set + the fracture layer (reduced toughness).
   No new core. Possibly a convenience constructor / material-zone tagging at the
   section level.

---

## 7. Testing / oracle matrix (Zone-A pytest)

| # | Test | Reference / oracle |
|---|---|---|
| V0 | elastic round-trip (`f<0`) | `σ=Eε`, `E_alg=E` |
| V1 | **reduce-to-`J2Plasticity`**, `N=0`, Voce monotone | analytic uniaxial backbone, ~1e-12 |
| V2 | monotone, linear-kin ≡ iso (`N=1`,`γ=0`) | pins `C_k` slope scaling vs analytic |
| V3 | reversed cycle, single AF | Bauschinger offset, saturation span `2·X_sat` |
| V4 | strain-controlled cyclic, `N=3` | stable hysteresis loop (Chaboche calib) |
| V5 | **stress-controlled ratcheting** | progressive `εᵖ` accumulation — the headline differentiator vs `Steel02` |
| V6 | consistent-tangent FD | `|E_alg − dσ/dε_FD| < tol` at several plastic states |
| **V7** | **3-D `LadrunoJ2`-condensed-to-uniaxial vs `LadrunoUniaxialJ2`** | **~1e-12** with identical `(C_k,γ_k,iso)` — THE reason this material exists |
| V8 | `getInitialTangent`==`E` after deep plastic loading | proves the `SimplifiedJ2` modal bug is absent |
| V9 | `sendSelf`/`recvSelf` round-trip of full state (`N=3`) | parallel/serialization hygiene |

Smoke (L0): single truss / 1-fiber section push-pull, `uniaxialMaterial
LadrunoUniaxialJ2`, Zone-A pytest (CPython 3.12 env per
[[../.claude/.../memory/project_opensees_test_env]]).

---

## 8. Risks / open questions

> [!note]
> **Naming (§2a) — RESOLVED:** `LadrunoUniaxialJ2` (user, 2026-06-02).

> [!question]
> **Tabulated-iso C¹ spline reuse.** Share the *exact* spline the 3-D class uses
> (so V7 holds) — confirm the 3-D tabulated path is actually factored into
> `LadrunoHardening.h`, not still inline in the kernel. If the 3-D tabulated mode
> is still deferred (it is, per [[10_ladruno_j2_plasticity]] "still deferred"),
> v1 ships **Voce+linear only** and tabulated lands when the 3-D one does — keeps
> the shared-header contract honest.

> [!question]
> **Ratcheting fidelity disclaimer.** Single-AF over-predicts ratcheting (the V5
> "signature"); multi-AF with split `γ_k` improves it but Ohno–Wang is the real
> fix. Document the limit; do **not** claim quantitative ratcheting accuracy in
> the banner.

- **No `exit()`**, real `getInitialTangent`, real `revert*`/`setParameter`/`Print`
  — the `SimplifiedJ2` anti-patterns, corrected (same rule as the 3-D class).
- **Backwards compat:** new material, additive; no existing model affected.
- **Oracle drift:** the V7 1e-12 gate is the tripwire — if a 3-D hardening edit
  ever breaks it, the shared header has been bypassed.

---

## 9. Implementation plan (step order)

1. **Factor `LadrunoHardening.h`** out of `LadrunoJ2Kernel.h` (Voce+linear `σ_y`
   + `dσ_y`; per-term AF update + `H^{kin}`). Refactor 3-D `LadrunoJ2` to consume
   it; **re-run the 3-D battery (8/8) — byte/behaviour unchanged** before adding
   any 1-D code. *(This is the precondition for the oracle property.)*
2. **`LadrunoUniaxialJ2.{h,cpp}`** (UniaxialMaterial): state `{εᵖ, ε̄ᵖ, X[N], Δγ}`,
   `setTrialStrain` per §3.3/§3.5, `E_alg` per §3.4, `getInitialTangent=E`,
   clean `commit/revert/sendSelf/recvSelf`, `getResponse("plasticStrain"|
   "plasticWork"|"backStress")` (the §6.1 fatigue hooks), real `Print` +
   **extensible `setParameter` dispatch** (the §6.2 buckling-length channel must
   be addable later without reshaping it). `setTrialStrain` stays scalar.
3. **`OPS_LadrunoUniaxialJ2.cpp`** parser (`-iso voce|table`, `-kin N ...`); register
   in the uniaxial `OPS_*` table + Python/Tcl.
4. **classTag + broker + build** — `MAT_TAG_LadrunoUniaxialJ2 33000` (`// Ladruno`),
   `FEM_ObjectBrokerAllClasses`, `material/uniaxial/CMakeLists.txt`.
5. **Zone-A battery** §7 (V0–V9), incl. the **V7 3-D cross-check** (drive 3-D
   `LadrunoJ2` `PlaneStress`→uniaxial and compare).
6. **Build-control:** `LEDGER_implementations.md` row (classTag 33000, uniaxial,
   shipped); `LEDGER_vanilla_files.md` rows for every `classTags.h`/broker/registry
   edit; `LEDGER_quirks.md` note (uniaxial-band 33000 reuse across registries);
   banner line in `Ladruno_scripts/banner_features.txt` → `patch_banner.py`.
7. **PR off `ladruno`** (`--base ladruno`, per CLAUDE.md). One logical PR; verify
   OPEN before any follow-up push (the stranded-commit trap).

**Deferred to follow-up ADRs:** tabulated iso (gated on the 3-D one, §8);
`LadrunoFatigue` fracture/LCF wrapper (§6.1); local-buckling degradation (§6.2);
weld/HAZ zones (§6.3); IMPL-EX code path (structure-only hook stored in `Δγ`).
