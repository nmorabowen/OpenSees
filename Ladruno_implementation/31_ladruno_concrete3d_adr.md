---
title: "ADR 31 — LadrunoConcrete3D (CDPM2-grade solid concrete): implementation spec"
project: Ladruno
type: ADR / implementation spec
status: PROPOSED — formulation scoped; 15-agent adversarial sweep (5 design dims × red-team) folded in; no code
related:
  - "[[19_ladruno_rc_shell_adr]]"        # SIBLING: LadrunoRCConcrete (33015) — shell/membrane/MCFT axis
  - "[[10_ladruno_j2_plasticity]]"       # the return-map IMPL-EX donor (LadrunoJ2Kernel returnMapDamaged + dScaleOverride)
  - "[[09_finite_strain_material_wrapper]]" # LogStrainNDMaterial (33010) finite-strain lift
  - "[[09_ladruno_brick]]"               # LadrunoBrick -geom finite host + gmsh SENB harness
  - "[[15_lemaitre_ductile_damage_adr]]" # shared LadrunoDamage.h; IMPL-EX-freezes-damage precedent
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
tags: [adr, material, concrete, plastic-damage, cdpm2, menetrey-willam, triaxial, confinement, implex, finite-strain]
updated: 2026-06-16
---

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

# ADR 31 — `LadrunoConcrete3D` (CDPM2-grade solid concrete)

**Status:** PROPOSED. Formulation scoped and **adversarially reviewed** (15-agent workflow:
4 ground-truth recon readers + 5 design panels + 5 red-team critics + synthesis; 15 blocking
+ 22 major findings triaged and folded in below). This is the implementation spec; **no code
has landed.** classTag **33017** (ND band; 33016 is reserved for `LogStrain2D` per
[[25_ladruno_plane_elements_adr|ADR 25]]).

> [!info] What LadrunoConcrete3D is, in one line
> A **3D solid** concrete constitutive model: effective-stress **Menétrey–Willam plasticity**
> (smooth 3-invariant surface, non-associated flow, confinement-aware ductility) + **dual
> scalar damage** (tension `ω_t` / compression `ω_c`), integrated with a real semi-implicit
> return map and an IMPL-EX robustness tier. ≈ **CDPM2** (Grassl et al. 2013) with the fork's
> robustness stack. It is the **triaxial/confinement** sibling of [[19_ladruno_rc_shell_adr|`LadrunoRCConcrete`]]
> (the shell/membrane/MCFT model), **not** a replacement for it.

---

## 1. Driver & goal — why ASDConcrete3D is insufficient

ASDConcrete3D (Petracca) is the fork's robust general-purpose concrete model and the donor
for much of this design. But the recon pass **verified in code** two structural limits that
no amount of input tuning fixes:

1. **No real plasticity / no triaxial surface.** ASDConcrete3D is **damage-only**: `compute()`
   (`ASDConcrete3DMaterial.cpp:2322–2481`) takes the elastic predictor `σ̃ = C₀:(εₙ−εₙ₋₁)`
   (cpp:2344–2348), spectrally splits it, and fabricates "plastic damage" by back-solving
   `d_plastic = 1 − q/σ̃_eq` off a 1-D hardening backbone (cpp:2429–2433). **There is no yield
   function, no flow rule, no plastic multiplier, no return map.** The `lublinerCriterion`
   (cpp:2483–2497) is only a *scalar equivalent-strain measure* — a Drucker-Prager-like
   envelope with **no Lode-angle (θ) dependence and no compression cap**. Confinement is
   therefore **emergent, not constitutive**: the multiaxial state shifts position on a *fixed*
   compression backbone via the envelope's hydrostatic (`α·I₁`) term and the `Kc` meridian knob.
   This works in the **calibrated band** — [[21_rc3d_validation_gates|Gate 2]] verified
   `fcc/Mander` within **5%** for `p/fc∈[0,0.20]`, including the ε_peak ductility growth — but
   it is a **single-meridian, no-cap, fixed-ductility, no-dilation-knob** approximation: it does
   not generalize across triaxiality/loading paths (no Lode), bounds nothing at high hydrostatic
   compression (no cap), cannot *control* confinement-ductility (a byproduct of `Kc`+backbone,
   not a designed function of `σ`), and is a robustness sink under high `p` (Gate 2 needed
   `KrylovNewton`+step-halving; plain Newton diverged). §4.1–§4.2's MW cap, Lode `r(θ,e)`,
   ductility measure `x(σ)`, and dilatancy knob `Df` are the constitutive versions of all four.

2. **Compression under-regularized; no unilateral effect.** One scalar `lch` feeds both
   backbones, the compressive measure uses only negative principals, and damage is monotonic
   with no crack-closure stiffness recovery.

The state-of-the-art single-surface model that fixes both while staying robust+fast is
**CDPM2** (effective-stress plasticity + two-scalar damage). Its anatomy maps almost 1:1 onto
ASDConcrete3D's machinery — so this is **re-aiming a proven architecture, not a rewrite.**

---

## 2. Decision summary

Build a **header-only, OpenSees-free `LadrunoConcrete3DKernel.h`** (the
[[10_ladruno_j2_plasticity|`LadrunoJ2Kernel`]] "one core, many views" doctrine) implementing:

- **Plasticity:** effective-stress **Menétrey–Willam** 3-invariant yield surface (friction
  `m0`, eccentricity `e`, Willam–Warnke Lode `r(θ,e)`), **non-associated** flow with explicit
  dilatancy, hardening `qh1`(pre-peak)/`qh2`(post-peak) driven by a **ductility measure
  `x(σ)`** so post-peak compressive ductility grows with confinement. A **real semi-implicit
  return map** with a dedicated **apex/vertex sub-algorithm**.
- **Damage:** dual scalar `ω_t`/`ω_c` on the spectral tension/compression split (ASDConcrete3D
  machinery, reused), strain-driven, crack-band-regularized in **both** `Gf` and `Gc`, with a
  **unilateral** crack-closure recovery.
- **Robustness (three tiers, one kernel):** Tier-1 implicit-accurate return map (default);
  Tier-2 error-controlled **IMPL-EX** (`-implex`, user-selected); Tier-3 **explicit-dynamic**
  (no tangent needed). Duvaut–Lions viscous (`-eta`) available under any tier.
- **Finite strain:** free via `nDMaterial LogStrain $t $c3d` + `LadrunoBrick -geom finite`,
  with a kernel out-contract fix so IMPL-EX does not corrupt the LogStrain `bᵉ` recovery (§4.4).

**Reuse provenance — corrected after red-team (this was the #1 false premise in the draft):**

| Reuse from | What | NOT from |
|---|---|---|
| **ASDConcrete3D** | `StressDecomposition` spectral projectors `PT`/`PC`; `HardeningLaw`+`regularize()` crack-band; `CrackPlanes`; **per-sub-object serialization**; IMPL-EX *scalar bookkeeping* (`time_factor`, `dtime_0`, EC error monitor at cpp:1749–1762, two-step `svt_commit_old` rotation) | — |
| **LadrunoJ2Kernel** | the **return-map IMPL-EX seam** (`returnMapDamaged` + `dScaleOverride`); analytic consistent tangent scaffold | — |
| **NEW code** | the MW surface, non-associated flow, ductility `x(σ)`, the semi-implicit return map + apex handler, dual-projector tangent | ❌ there is **no plastic update to "port wholesale" from ASDConcrete3D** — it has none |

---

## 3. Scope-fence & classTag

- **`LadrunoConcrete3D` (33017)** owns the **solid / triaxial / confinement** axis: 3D continuum
  (`LadrunoBrick`, `stdBrick`, tets), confined columns, joints, deep beams, bearing.
- **`LadrunoRCConcrete` (33015)** owns the **shell / cracked-membrane** axis: in-plane RC walls
  via MCFT compression softening on `ASDShellQ4`+`LayeredShellFiberSection`.
- **Siblings, not competitors.** Different classTags, different element hosts, documented
  relationship. They share the ASDConcrete3D spectral/hardening lineage but diverge on physics
  (return-map plasticity here; β-softening membrane there). No tag or scope collision.

`#define ND_TAG_LadrunoConcrete3D 33017` (ND band; after `LogStrain2D` 33016-RESERVED).

---

## 4. Formulation (each dimension folds its red-team fixes)

### 4.1 Yield surface — Menétrey–Willam, normalized

**PINNED against Grassl et al. 2013 (IJSS 50:3805, arXiv:1307.6998) — equation numbers cited.**
fc-normalized three-invariant form. With `σ̄_V = I₁/3` (Eq.12; in code stored as `ξ = I₁/√3`, so
`ξ/(√3·fc) = σ̄_V/fc`), `ρ̄ = √(2J₂)` (Eq.13), Lode angle `θ̄` (Eq.14):

```
                ⎧                                      ⎫²
f_p(σ̄_V,ρ̄,θ̄;κp)=⎨ [1−qh1(κp)]·(ρ̄/(√6 fc) + σ̄_V/fc)² + √(3/2)·ρ̄/fc ⎬          ← Eq.(18)
                ⎩                                      ⎭
              + m0·qh1²(κp)·qh2(κp)·[ ρ̄/(√6 fc)·r(θ̄,e) + σ̄_V/fc ] − qh1²(κp)·qh2²(κp)

qh1=qh2=1  ⇒  (3/2)ρ̄²/fc² + m0[ρ̄·r/(√6 fc) + σ̄_V/fc] − 1 = 0   (failure surface Eq.21 = Menétrey-Willam)
m0 = 3(fc²−ft²)/(fc·ft) · e/(e+1)                               ← Eq.(20)
r(θ̄,e) = Willam-Warnke elliptic Lode fn (Eq.19), e ∈ (0.5,1]   (0.5 = convexity limit, EXCLUSIVE)
```

The `[1−qh1]` term is an **ellipsoidal hardening cap** that closes the surface during hardening
(`qh1<1`) and vanishes at peak (`qh1=1` → the open MW cone). My earlier draft's `qh1·qh2`/`qh1²qh2²`
guess was **wrong** for `qh≠1`; both reduce to Eq.21 at `qh=1` (why P0/P1 passed), but only Eq.(18)
is correct under hardening — now implemented and gated (HA reduce-to-P1 = 2.6e-13).

**Fixes folded in:**
- **[DONE] Pinned against the literature.** `f`, `m0`, `r(θ,e)` carry Eq. numbers (18/19/20/21);
  P0 oracle asserts the un-hardened surface reproduces the meridian/eccentricity identity **and**
  Kupfer `fcc/fc=1.16` (recovering `e≈0.52`) before any return-map code — both PASS.
- **[MAJOR] `m0` and `e` are coupled — do not freeze both independently.** Choose explicitly:
  **(a)** freeze `e≈0.525` and treat equibiaxial strength as a *derived output* (set the Kupfer
  tolerance to accommodate), **or (b)** derive `e` from an equibiaxial-ratio input
  (`fcc/fc`, default 1.16) so the deviatoric trace passes through the user's strengths by
  construction. **Recommend (b)** for "usable out of the box." Either way, `e` is a *validation
  target, never fit to Kupfer*.

### 4.2 Return map — hardening, ductility, non-associated flow, apex (PINNED to Grassl 2013)

**Hardening laws (Eqs. 30–31), `κp`-driven:**
```
qh1(κp) = qh0 + (1−qh0)(κp³−3κp²+3κp) − Hp(κp³−3κp²+2κp)   if κp<1, else 1      ← Eq.(30)
qh2(κp) = 1 if κp<1, else 1 + Hp(κp−1)                                            ← Eq.(31)
```
`qh1` ramps the surface from `qh0` to the failure cone over `κp∈[0,1]`; `qh2` adds post-peak
plastic hardening for `κp>1`. **Effective-stress plasticity is MONOTONIC (no peak)** — the
failure surface is reached exactly at `κp=1` (gated: `σ11(κp=1)=fc` to 7e-4); the *peak/softening
is the DAMAGE part (P2)*, not plasticity.

**Hardening-variable evolution + ductility measure (Eqs. 32–36):**
```
κ̇p = (λ̇‖m‖ / xh(σ̄_V)) · (2cosθ̄)²                                              ← Eq.(32)
xh(σ̄_V) = Ah − (Ah−Bh)exp(−Rh/Ch)        if Rh≥0   ;   Eh exp(Rh/Fh)+Dh  if Rh<0  ← Eq.(33)
Rh = −σ̄_V/fc − 1/3   (Eq.34) ;  Eh=Bh−Dh (Eq.35) ;  Fh=(Bh−Dh)Ch/(Ah−Bh) (Eq.36)
```
`Rh` is **the confinement term**: compression (`σ̄_V<0`) → `Rh>0` → larger `xh` → slower `κp` →
more plastic strain to reach the failure surface ⇒ **confinement-dependent ductility, by
mechanism** (gated HD: strain at `κp=1` grows 0.0012→0.0032 from unconfined to `p/fc=0.1`). Ship
CDPM2's published `Aₕ/Bₕ/Cₕ/Dₕ` (calibrated from peak strains — recalibrate for fork data, §6).

- **Non-associated potential** `g_p` (Eqs. 22–29) — **Lode-INDEPENDENT** (no `r(θ̄)`; the yield
  surface has `r`, the potential does not). The deviatoric flow `∂g/∂ρ̄ = 3ρ̄/fc² + m0/(√6 fc)`
  carries **no `r`** (corrected from the P0 draft). Dilatancy via `Df = −m₂/m₁` (Eq.27) sets the
  volumetric flow through `mg(σ̄_V)` (Eq.23) which decreases with confinement. v1 oracle uses a
  simplified constant-`Df` volumetric flow; the full `mg` potential is a follow-on (does not
  change peak strength — flow-independent). Non-associated ⇒ **non-symmetric** tangent (§4.4/§4.5).
- **Semi-implicit return — TWO maps (load-bearing for the C++ implementer):** the **perfect-plastic
  limit is a 3-unknown `(ξ,ρ,Δλ)` Newton with a closed-form analytic Jacobian**; the **hardening map
  adds `κp` as a 4th unknown with a numerical Jacobian** (so its stress is only ~1e-8 — the C++ build
  PR owes the **analytic 4×4 Jacobian**, incl. piecewise `∂xh/∂σ̄_V` across `Rh=0` and cubic
  `∂qh1/∂κp`, FD-checked against the oracle). `θ̄` frozen → radial deviatoric return. The **consistent
  tangent is taken on the perfect-plastic map** (analytic inner Jacobian → clean FD). **Reduces to the
  verified perfect-plastic map at qh0=1,Hp=0 (HA = 2.6e-13).**
- **[KNOWN GAPS — build-PR scope, must be in the handoff]** (a) the **hardening apex is unhandled**
  beyond a hydrostatic-vertex projection; near the tensile apex the frozen-θ residual can read ~0
  while the point is off-surface — the **honest independent-θ `f` flag (added) reports
  `converged=False`** there, but the dedicated apex/Koiter sub-algorithm + drift correction are owed.
  (b) **Off-meridian first-yield drift** (κp~0, small surface + steep `qh1` ramp): the semi-implicit
  return leaves an O(0.1) off-surface drift that does NOT vanish under refinement — needs sub-stepping
  near first yield. Verified only on the compressive meridian (the axisymmetric driver can't make
  off-meridian states); diagnostic `HE` records it.
- **[MAJOR] Order of operations for `x(σ)` vs `Gc`.** lch-regularize the *intrinsic* softening
  law **first** (so `Gc` is size-objective at reference confinement), **then** apply `x(σ)` to
  the plastic increment in a way that **preserves total compressive fracture energy per unit
  area at each confinement**. Gate this with a `Gc` mesh-objectivity test **at multiple
  confinements** (not just unconfined).
- **[BLOCKING/MAJOR] Apex & Lode-corner handling.** The MW surface closes to a hydrostatic-
  tension **apex** (singular flow direction) and has **Lode corners at θ=0/60°**. The Gf SENB
  tension battery drives stress straight to the apex, so a naïve single-surface return diverges
  or returns to the wrong point. Ship a **dedicated apex/vertex return** (hydrostatic-axis
  return when the trial projects beyond the cone vertex; Koiter multi-surface or C¹-blended
  flow at corners). Tier-1 "unconditional material stability" holds **only with** this handler.
- **[MAJOR] Semi-implicit freeze set is load-bearing.** Pin which quantities (`θ`, `x(σ)`) are
  held at trial vs updated, against Grassl 2013, and add an oracle **drift check**: after the
  return, evaluate `f` at the converged stress with the *updated* `θ` and assert `|f|/fc < tol`
  on triaxial paths; document the O(Δθ) semi-implicit tax as **distinct from** the IMPL-EX
  O(Δt) tax and prove it vanishes under step refinement.

### 4.3 Damage & regularization — dual scalar, with the energy contract pinned

- Reuse ASDConcrete3D's spectral split (`PT`/`PC`, `computePjj`) to apportion effective stress
  into tension/compression; `ω_t` acts on `ST`, `ω_c` on `SC`; nominal
  `σ = (1−ω̄_t)ST + (1−ω̄_c)SC`.
- **[BLOCKING] Backbone abscissa contract — keep it STRAIN-like.** `regularize()` and
  `computeFractureEnergy()` (`ASDConcrete3DMaterial.cpp:947–1001, 1183–1255`) are **hard-wired
  to a strain abscissa** (the energy integral `∫y dx`, the `E=y/x` secant, the `gmin` triangle).
  Do **not** reinterpret the abscissa as a dimensionless `κ_p`. Keep the damage abscissa an
  **equivalent inelastic strain**; apply `x(σ)` as a *rate* multiplier inside the kernel, never
  as the backbone abscissa. (Resolves the draft's internal strain-vs-`κ_p` contradiction.)
- **[BLOCKING] Prove Gf/Gc, don't assume "same container ⇒ free objectivity."** CDPM2 carries a
  *real* plastic strain `κ_p` separate from the damage driver `κ_d`, so ASDConcrete3D's
  `q`-based plastic bookkeeping in `computeFractureEnergy` no longer holds. The **only honest
  gate** is a numpy-oracle test that integrates the **actual dissipated energy (plastic +
  damage)** over a fully-softened uniaxial path and asserts it equals `Gf` (tension) / `Gc`
  (compression) **after `regularize()` at 2–3 `lch` values**.
- **[MAJOR] `gmin` snap-back floor breaks exact energy on coarse meshes.** `gmin = 1.01·½·peak.y·peak.x`
  (cpp:961–965) clamps dissipation up when `lch` is large, so "area = Gf·lch_ref/lch" and the
  floor cannot both hold. Assert the objectivity tests stay in the regime where the floor does
  **not** bind, and **emit a one-time warning** when it does (energy is being clamped).
- **[BLOCKING] Unilateral recovery must be tier-independent.** Re-splitting the *converged
  effective stress* each step is what follows load reversal — but in IMPL-EX the projector is
  **frozen** (`PT=PT_commit`, cpp:2352–2363), so a crack closing during a Tier-2 step keeps
  multiplying by `(1−ω_t)` and the cyclic value proposition silently dies exactly where it's
  sold (near peak). **Fix:** recompute the spectral split of the converged effective stress for
  the nominal recomposition in **both** tiers (cheap, no tangent); freeze **only** the
  damage-driver extrapolation. Put the recovery knob on the correct channel —
  `ω_t,eff = ω_t·(1 − s_rec·g_close(σ))` with `g_close` a principal-sign closure indicator —
  **not** as a `PC`-component multiplier. Gate with a tension→compression reload stiffness-
  recovery test.
- **[MAJOR] Cross-damage coupling default.** `cdf` defaults to 0 (coupling OFF) and `alpha` is
  hard-coded 1.0 in `mix_dam` (cpp:2436–2448). Decide CDPM2's actual `ω_c←ω_t` coupling, map to
  `(cdf, alpha)` with **non-trivial fc/ft-derived defaults**, expose `alpha`, and test the
  default-parameter combined-T/C descending branch — or rewrite the rationale if `cdf=0` ships.
- **[MAJOR] Full dual-projector tangent.** The damaged consistent tangent is **not** J2's
  single-scalar cross-term: it needs `−ST⊗∂ω_t/∂ε − SC⊗∂ω_c/∂ε`, the coupling derivative if on,
  and a documented decision on the `∂PT/∂σ` term (include it, or justify neglect and accept
  more Tier-2 drops near eigenvalue coalescence). FD-verify across a reversal **and** an
  eigenvalue-crossing path, not just monotonic loading.

### 4.4 Robustness — three tiers, one kernel

- **Tier-1 (default): implicit-accurate** semi-implicit return map. The consistent tangent is
  **non-symmetric** ⇒ **requires an unsymmetric solver** (`UmfPack`/`FullGeneral`). **Verified
  in the oracle (P1-tangent gate):** strongly non-symmetric for non-associated flow
  (`‖C−Cᵀ‖/‖C‖ ≈ 0.46` at `Df=0.3`), and **non-symmetric even in the associated limit** (`e=1,
  Df=1` → `≈0.024`, ~20× smaller but nonzero). **The associated-limit residual is the frozen-
  eigenvector spectral recompose** (`σ=V·diag(σₚ)·Vᵀ` with `V` held at the trial drops the
  eigenprojection/spin terms `dV/dε`) — falsified that it is the Lode θ-freeze: a principal-space
  off-meridian associated state is machine-symmetric (~2e-10), and the full-tensor asymmetry scales
  **linearly with shear** (→0 as shear→0), FD-step-independent (gate T3c). So the unsymmetric solver
  is a **hard requirement UNCONDITIONALLY**, not only for `Df<1` (OpenSees has no runtime guard —
  document in the user guide + banner note; warn on a symmetric solver).
- **[BLOCKING] Tier-2 (IMPL-EX) must freeze the PLASTIC state too, not just damage.** Freezing
  only `ω` and still solving the non-associated softening return implicitly leaves the tangent
  non-symmetric and **indefinite on the softening branch** — the SPD promise is false. Following
  Oliver/Huespe, **extrapolate the plastic-strain increment + the dual damage** so the explicit
  effective stress is linear in `Δε` and the Tier-2 secant is `D_dam(ω̃):C₀`. **[CORRECTED — oracle
  P3 IMPL-EX adversarial review, #301 review/NUM-1]:** this secant is symmetric-part SPD **only in
  SINGLE-SIGN principal regimes** (all-tensile / all-compressive `σ̄_x`, where `D_dam` is one positive
  scaling — the original `(1−ω̃)·C₀` claim). On a **MIXED-SIGN, high-`ω` direction-contrast** state (a
  tensile-damaged principal beside an undamaged compressive one, `ω_t > ~0.97`) the two branch slopes
  `(1−ω_t) ≠ (1−ω_c)` make `D_dam:C₀` non-commuting and its symmetric part **indefinite** (`λ_min ≈
  −5e2` for `σ̄=[1,−2,−2]`, `ω_t=0.99`) — the intrinsic dual-damage IMPL-EX limitation, not a fixable
  bug (genuine dual-damage *consistency* and unconditional SPD are mutually exclusive). The secant
  nonetheless remains the **exact consistent tangent of the reported explicit stress** (FD-verified)
  and is far better-conditioned than Tier-1, and the **committed physics is exact**. Falsification
  gates (oracle PI1/PI5/PI6): single-sign snap-back SPD; `D_dam(ω̃):C₀ == d(σ_rep)/d(Δε)`; the
  mixed-sign indefiniteness **pinned** (not claimed away); and the extrapolation **time-ratio CLAMPED
  to `[0,R_max]`** (a step-cutting solver's `dt`-growth otherwise over-extrapolates `ω̃` past `[0,1)`
  ⇒ stress collapse + injects unbounded plastic strain via `Δε̃ᵖ`; negative `dt` ⇒ `r=0`).
- **Tier-3 (explicit-dynamic):** same kernel, `do_tangent=false`; no global tangent ⇒ softening
  is a non-issue (`CentralDifferenceLadruno`/`ExplicitBathe`; LS-DYNA CSCM/KCC philosophy).
- **[BLOCKING] Duvaut–Lions at the PLASTIC level.** Relax the plastic multiplier toward the
  inviscid return with factor `Δt/(η+Δt)` so `η→0` recovers the inviscid return **exactly**
  (and the *same* instructions, for the byte gate). Do **not** claim to inherit ASDConcrete3D's
  `rate_coeff` — that blends a *damage driver*, not a plastic stress. Validate with the
  closed-form 1-D overstress oracle.
- **[BLOCKING] Finite-strain IMPL-EX must not corrupt LogStrain.** `LogStrainNDMaterial.cpp:190–204`
  recovers `εᴱ = Cᵉ:τ` from the **returned** stress assuming a linear-elastic inner law. Under
  IMPL-EX the returned stress is frozen/extrapolated, so the committed `bᵉ = exp[2εᴱ]` is
  polluted and compounds step-over-step. **Fix the kernel out-contract:** `returnMap` must
  expose the **implicit effective stress** separately, and LogStrain recovers `εᴱ` from
  `σ_eff,implicit` regardless of which tier degraded the reported nominal stress. Regression:
  rigid-rotation + softening under `-geom finite` with forced demotes, assert objectivity AND
  committed `bᵉ` matches the all-Tier-1 trajectory to O(Δt).
- **[BLOCKING/MAJOR] The automatic implicit↔IMPL-EX hybrid is DESCOPED from v1.** It is the
  single biggest robustness risk (chattering, energy jumps at the switch, per-GP tier
  divergence) and has no gate. **v1 ships the three tiers as user-selected modes** — Tier-1
  default, `-implex` for Tier-2, explicit solver for Tier-3 — exactly ASDConcrete3D's proven UX,
  with **no automatic switching** and therefore none of the switch bug surface. The hybrid is a
  **Phase-5 research item** with its own co-designed gate: if ever built, the switch fires
  **only at a step boundary** (full rollback + re-integrate, never mid-Newton), latches ≥2 steps
  (anti-chatter), **re-seeds** IMPL-EX history (`svt_commit_old`, `dtime_0`) on every transition,
  and commits the **implicit** stress/damage while using the frozen-SPD tangent **only as the
  iteration matrix** ("IMPL-EX for the Jacobian, implicit for the residual"), gated by
  `EnergyBalanceRecorder`. The demote trigger must **not** be "symmetric-part loses positive-
  definiteness" — that is a category error for non-associated flow; use local return-map residual
  stagnation / local Newton solvability / global `ConvergenceTest` stall.

### 4.5 Architecture

- **Kernel:** header-only OpenSees-free `LadrunoConcrete3DKernel.h`, de-static'd spectral
  helpers. **Tensor convention (LadrunoJ2 lineage):** the kernel stores symmetric tensors in
  `{00,11,22,01,12,02}` with **true tensor** off-diagonals (J2 squares them directly); the
  engineering↔tensor conversion is the wrapper's responsibility at the OpenSees boundary, exactly
  as in `LadrunoJ2.cpp`. **NB for P2:** ASDConcrete3D's spectral split helpers use the engineering
  shear convention — reconcile (convert at the borrow seam) before reusing them for the dual-damage
  projectors, or the off-diagonal damage apportioning is wrong by a factor. The `nDMaterial
  LadrunoConcrete3D` wrapper (33017) includes the kernel.
- **[MAJOR] v1 views = 3D + finite only.** The whole validation heart is the 3D triaxial
  battery; `LogStrain`-finite is free. **Defer PlaneStrain/AxiSymmetric to Phase 2** with their
  own reduce-to-3D verification — and note AxiSymmetric is **not** a pure `vmap` of a 3D run
  (hoop stress couples through the radius); if kept it must be verified against a thick-cylinder
  closed form, not a 3D probe.
- **[BLOCKING] Serialization = hybrid.** Use ASDConcrete3D's per-sub-object
  `serializationDataSize()/serialize()/deserialize()` for the **variable-length** borrowed
  objects (HardeningLaw backbones, CrackPlanes) — they are variable-length and cannot be
  flattened into the J2/RC fixed-index `Vector`. Use the flat-index `Vector` idiom **only** for
  the fixed scalar block (`K,G,fc,ft,Gf,Gc,e,Df,η`, flags, the 6-tensor committed history, and
  the IMPL-EX `svt_commit/svt_commit_old/dtime_n_commit`).
- **[BLOCKING] IMPL-EX is ported from ASDConcrete3D, NOT inherited from the RC sibling** —
  `LadrunoRCConcrete`/`LadrunoRCKernel` carry **zero** IMPL-EX (recon: grep count 0). Budget it
  as a fresh port (no RC precedent for its serialization/restart correctness).

### 4.6 Confined-fiber view — hoop-spring condensation (P5 view, "Mander by mechanism")

A distinguished reduced view that brings the triaxial model into a **1-D fiber** beam-column
without a pre-baked Mander backbone: the *same* kernel, condensed against the lateral confinement
supplied by the transverse steel. The whole confinement-aware stack (3D solid → 2D → confined
fiber) then runs **one constitutive law** with consistent physics at every scale.

**Mechanism = static condensation with a non-zero lateral residual.** Every reduced view drives
the axial strain `ε11` (from the fiber section) and Newton-solves the lateral strains so a chosen
lateral residual vanishes. The dimension is set by the residual:

```
free surface  (plain 2D/1D reduction):  R = σ_lat(ε_lat)                 = 0
ACTIVE confinement (prescribed):         R = σ_lat(ε_lat) + p_hoop         = 0   (constant p input)
PASSIVE confinement (real hoops):        R = σ_lat(ε_lat) + σ_hoop(ε_lat)  = 0   (hoop spring/law)
```

This is the **same nested-Newton condensation** the plate-fiber view already uses for `σ33=0`
([[19_ladruno_rc_shell_adr|`LadrunoRCConcrete`]] / `LadrunoJ2` PlaneStress); only the residual's
target moves from `0` to the hoop term. Output: `σ11` + the **condensed 1-D tangent**
`dσ11/dε11` (static condensation of the 6×6 with the hoop spring `∂σ_hoop/∂ε_lat` added to the
lateral diagonal).

**Why passive (spring), not just active (prescribed `p`):** hoop confinement is *generated* by
the concrete dilating against the hoops, not imposed. The **non-associated dilatancy `Df`** sets
the lateral expansion → mobilizes hoop tension → supplies the confining pressure, self-consistently;
the strength **and** ductility gains then emerge from the MW cap + Lode + ductility measure `x(σ)`.
That is the constitutive generalization of what Mander encodes empirically (pressure→ductility),
and it rides the compressive meridian where CDPM2 is correct and ASD's single-meridian envelope
is not. The user supplies the hoop **stiffness/law** (`σ_hoop(ε_lat)`, optionally nonlinear with
yield), not a pressure.

**Caveats (a point material cannot see member geometry — record in [[LEDGER_quirks]]):**
- **Arching / effective confinement `ke` is a section/geometry effect**, not a material one. Hoop
  spacing and in-plane arching make `p` non-uniform; this must be baked into the hoop spring
  stiffness (`ke·ρ_hoop·…`). Physical/measurable input, but not free.
- **Rectangular ties are anisotropic** (`σ22≠σ33`): the symmetric `ε22=ε33` condensation is exact
  for **circular/spiral**, approximate for rectangular. A two-spring (`K_x,K_y`) condensation
  recovers the anisotropy and is the clean extension.
- **Cover vs core** handled the normal fiber way (confined view on core fibers, free-reduction
  view on cover fibers).
- **Robustness/state:** the inner condensation Newton inherits the softening-convergence
  difficulty (the high-`p` divergence Gate 2 saw); commit the **full** condensed lateral state
  (own plastic/damage history), like the `σ33=0` thickness stretch.

**Validation:** a confined-fiber-view leg reproducing the **Mander confined backbone** for
circular hoops across `p/fc` (the mechanism analog of the 3-D Gate 2), plus a reduce-to-plain-1D
check at zero hoop stiffness. Distinct from the PlaneStrain/AxiSym reductions (non-zero lateral
residual), so it lands as its own **P5 view**, after the 3-D triaxial physics is locked.

---

## 5. Validation & calibration battery

**Triple gate** (fork standard): numpy oracle `tests/_testbed/concrete3d_ref.py` (a from-the-
equations integrator — the **primary** accuracy reference, since no peer C++ model exists to
byte-match) + a standalone g++ build of the kernel + OpenSees pytest.

**[BLOCKING] Reduce-to-baseline gates — the achievable ones (drop "reduce-to-ASDConcrete3D"):**
1. all features off ⇒ **linear elastic** `σ=C₀:ε`, byte-identical.
2. `η=0` ⇒ inviscid, byte-identical to the non-viscous path.
3. damage-off, associated flow, MW→Drucker–Prager limit ⇒ match a **closed-form DP return-map
   oracle** to ~1e-10 (an independent leg that *does* hit a known surface).
4. Tier-3 explicit (`do_tangent=false`) committed stress == Tier-1 implicit on a single-element
   strain path, ~1e-12 (same kernel).
5. damage-only subsystem (plasticity off) ⇒ match a numpy re-impl of ASDConcrete3D's damage
   spine (tests the inherited damage code) — **only** the damage subsystem, never the surface.

**[BLOCKING] Triaxial battery — the heart, with TRUE σ₃ control.** The `_homog` driver is
all-DOF strain-prescribed and **cannot** impose constant lateral stress. Add a **lateral
stress-control Newton loop** (1 unknown by symmetry: iterate `ε_yy=ε_zz` so `σ_yy=σ_zz=σ₃,target`)
around the single-GP material call — trivial in the numpy oracle, replicated in Python around
the pyd for the OpenSees leg (or a meshed pressure-boundary Zone-B cell as confirmation).
Assert post-peak ductility **monotone in σ₃** against a *true* σ₃, not a held-strain surrogate.
Targets: **Kupfer 1969** biaxial envelope + **confined triaxial** at several σ₃ from Grassl's
CDPM2 papers.

**[MAJOR] Tiered Kupfer tolerance, parameters FROZEN at default:** ±8% on the two compression
peaks (equibiaxial ≈1.16 fc, asymmetric ≈1.27 fc at σ₂/σ₁≈0.5); ±15–20% on the tension–
compression and uniaxial-tension legs (larger experimental scatter). If the gate only passes by
retuning `e`, that is a **finding, not a pass**.

**[MAJOR/BLOCKING] Apex & Lode-corner gates (P1):** drive (a) pure hydrostatic tension to the
tensile apex — assert the return lands on the apex without NaN/non-convergence; (b) the Lode
corners θ=0/60° at high ρ — assert convergence and `f=0` to tolerance. Test **both** associated
and the shipping non-associated flow.

**[MAJOR] Gc mesh-objectivity on a CONFINED prism** (apply modest σ₃ so the crush band is
well-posed; unconfined compressive softening is snap-backy → make that a Tier-3 explicit demo).
Assert **energy** objectivity, not band-**width** objectivity (a single scalar `lch` cannot
deliver band-width objectivity — record this limit in [[LEDGER_quirks]]). Keep a
`CentralDifferenceLadruno` fallback ready for this gate.

**[MAJOR] Cross-platform tolerance floors** (Win→Linux drift history): g++-vs-oracle ~1e-7 abs /
1e-8 rel on smooth trajectories, ~1e-6 near the Lode vertices (θ→0/60°, where J₃ derivatives are
singular), 1e-9 **only** on the Lode-independent elastic/isotropic legs.

Plus: cyclic/unilateral stiffness-recovery, and rigid-rotation finite-strain objectivity (clean,
not xfail — isotropic plastic-damage co-rotates without J2's kinematic-backstress problem).

---

## 6. Parameters & defaults

| Param | Meaning | Default |
|---|---|---|
| `E, ν, ρ` | elasticity, density | required |
| `fc, ft` | uniaxial compressive / tensile strength | required |
| `e` | deviatoric eccentricity | derived from `fcc/fc` (default 1.16) per §4.1(b); ∈(0.5,1] |
| `m0` | friction | derived `3(fc²−ft²)/(fc·ft)·e/(e+1)` (coupled to `e`) |
| `Gf` | tensile fracture energy | fib correlation from `fc` (e.g. `Gf≈73·fc^0.18` N/m) |
| `Gc` | compressive fracture energy | **flag = weakest-calibrated knob**; default a published correlation, documented |
| `qh0` | initial yield (frac. of peak) | ~0.3 |
| `qh2,res` | residual friction level | pinned default (cite paper) |
| `Df` | dilatancy | default tied to a measured dilatancy angle |
| `Aₕ,Bₕ,Cₕ,Dₕ` | ductility-measure params | **CDPM2 published table, verbatim** |
| `η` | Duvaut–Lions viscosity | 0 (off ⇒ byte-identical to inviscid) |
| `-implex` | engage Tier-2 | off (Tier-1 default) |

Every literature default is tagged "recalibrate for HSC / fork data."

---

## 7. Phased roadmap (each phase has its exit gate)

- **P0 — ADR + scaffold.** This ADR; reserve `33017`; `LadrunoConcrete3DKernel.h` skeleton;
  `concrete3d_ref.py` numpy oracle; standalone g++ harness; the **lateral-σ₃ control loop**;
  pin the exact MW `f`/`m0`/`r(θ,e)` with equation numbers. **Gate:** un-hardened surface
  reproduces meridian ratio + Kupfer `fcc/fc` in the oracle.
- **P1 — MW plasticity, no damage.** Non-associated return map + ductility + **apex/corner
  handler**; analytic tangent. **Gate:** triaxial battery (Kupfer + confined σ₃, true control),
  DP-limit reduce-to-oracle 1e-10, apex/corner convergence, FD-tangent, drift check.
- **P2 — dual damage + crack-band.** `ω_t`/`ω_c`, `Gf`/`Gc` regularization, unilateral. **Gate:**
  actual-dissipated-energy = `Gf`/`Gc` after `regularize()` at 2–3 `lch`; gmsh SENB mesh-
  objectivity; confined-prism `Gc` objectivity at multiple σ₃; cyclic stiffness-recovery.
  (PlaneStrain/AxiSym views also land here, with their own reductions.)
- **P3 — robustness.** Tier-2 IMPL-EX (freeze plastic+damage), Duvaut–Lions, Tier-3 explicit.
  **Gate:** Tier-2 SPD-eigenvalue across snap-back; `η=0` & all-off byte gates; Tier-3==Tier-1
  stress 1e-12; 1-D overstress oracle.
- **P4 — finite-strain view.** `LogStrain` + `LadrunoBrick -geom finite` with the `bᵉ` out-
  contract fix. **Gate:** rigid-rotation objectivity; reduce-to-small-strain; Tier-2-under-
  finite `bᵉ` matches all-Tier-1.
- **P5 — confined-fiber view (§4.6).** Hoop-spring condensation (active + passive); 1-D fiber
  output for beam-column. **Gate:** reproduces the Mander confined backbone for circular hoops
  across `p/fc` (mechanism analog of Gate 2); reduce-to-plain-1D at zero hoop stiffness.
- **P6 (research) — automatic hybrid switch.** Only with the co-designed gate in §4.4.

---

## 8. Risk register (from the adversarial sweep — each with its resolution)

| # | Severity | Risk | Resolution |
|---|---|---|---|
| R1 | BLOCKING | "Port ASDConcrete3D's effective-stress plastic update / IMPL-EX wholesale" — **no such update exists** (damage-only model) | §2 reuse table: scalar machinery from ASD, return-map seam from J2Kernel, plasticity is new |
| R2 | BLOCKING | IMPL-EX freezing only damage ⇒ unregularized plastic snap-back, SPD promise false for non-associated flow | §4.4 Tier-2 freezes plastic multiplier + hardening + damage (true Oliver/Huespe); SPD-eigenvalue gate |
| R3 | BLOCKING | Finite-strain Tier-2 corrupts LogStrain `bᵉ` recovery (silent, compounding) | §4.4 kernel exposes implicit effective stress; LogStrain recovers `εᴱ` from it; `bᵉ` regression |
| R4 | BLOCKING | HardeningLaw abscissa as `κ_p` breaks `regularize()`/energy integral | §4.3 keep abscissa strain-like; `x(σ)` as a kernel-level rate multiplier |
| R5 | BLOCKING | "Same backbone container ⇒ free Gf/Gc objectivity" false (real plastic strain ≠ ASD's `q`) | §4.3 actual-dissipated-energy oracle = Gf/Gc after regularize at multiple `lch` |
| R6 | BLOCKING | Unilateral recovery dies in frozen-projector Tier-2; `s_rec` on wrong channel | §4.3 re-split converged effective stress in both tiers; `s_rec` via `g_close(σ)` |
| R7 | BLOCKING | Reduce-to-ASDConcrete3D / reduce-to-baseline byte gate impossible (different physics) | §5 oracle-based gates + DP-limit + elastic/η=0/Tier3 byte legs; damage-subsystem-only ASD check |
| R8 | BLOCKING | Confined-triaxial gate has no true σ₃ control (`_homog` is strain-prescribed) | §5 lateral stress-control Newton loop; P0 deliverable |
| R9 | BLOCKING | Auto hybrid switch: chatter / energy jumps / per-GP tier divergence, no gate | §4.4 **descoped from v1** → user-selected `-implex`; hybrid is gated Phase-5 research |
| R10 | BLOCKING | MW apex/Lode-corner return-map failure (Gf SENB drives to apex) | §4.2 dedicated apex/vertex sub-algorithm; §5 apex/corner gates in P1 |
| R11 | BLOCKING | Serialization: flat `Vector` can't hold variable-length borrowed sub-objects | §4.5 hybrid: per-sub-object scheme for var-length, flat block for scalars |
| R12 | BLOCKING | classTag: 33016 is reserved for `LogStrain2D` (ledger), not free | **33017** (verified against `classTags.h` + `LEDGER_implementations.md:75,89,90`) |
| R13 | MAJOR | Non-associated ⇒ non-symmetric tangent vs symmetric-solver assumption | §4.4 Tier-1 requires unsymmetric solver (documented); SPD is Tier-2-only |
| R14 | MAJOR | Demote-on-"loss of positive-definiteness" is a category error for non-assoc. flow | §4.4 demote on residual stagnation / local-Newton solvability / global stall |
| R15 | MAJOR | `x(σ)` ductility breaks `Gc` objectivity unless coupled to crack-band | §4.2 regularize first, then energy-preserving `x(σ)`; multi-confinement `Gc` gate |
| R16 | MAJOR | `m0`/`e` coupling; freezing both ⇒ surface inconsistent with strengths | §4.1 derive `e` from `fcc/fc` (option b); cite equation numbers |
| R17 | MAJOR | Calibration gaps (`Aₕ..Dₕ`, `Df`, `Gc`, `qh2`) — no defaults | §6 CDPM2 verbatim defaults; `Gc` flagged weakest |
| R18 | MAJOR | Cross-damage `cdf=0`/`alpha=1` hard-coded ⇒ coupling off by default | §4.3 fc/ft-derived `(cdf,alpha)` defaults; expose `alpha`; combined-T/C test |
| R19 | MAJOR | Two-scalar tangent ≠ J2 single-scalar cross-term; `∂PT/∂σ` dropped | §4.3 full dual-projector tangent; FD across reversal + eigenvalue crossing |
| R20 | MAJOR | v1 over-scoped (PlaneStrain/AxiSym; AxiSym ≠ vmap) | §4.5 v1 = 3D+finite only; reduced views deferred to P2 with own gates |
| R21 | MAJOR | Kupfer ±8% too tight across all legs; could silently retune `e` | §5 tiered tolerance, parameters frozen, Kupfer = validation not calibration target |
| R22 | MAJOR | g++-vs-oracle 1e-9 below cross-platform transcendental floor (Lode trig) | §5 1e-7/1e-8 smooth, 1e-6 near vertices, 1e-9 only Lode-independent legs |

---

## 9. Ledger / header / PR obligations

- **`LEDGER_implementations.md`** — add the `LadrunoConcrete3D` row with classTag **33017**,
  files, status, PR (when built).
- **`classTags.h`** — `#define ND_TAG_LadrunoConcrete3D 33017` with the standard `// N. Mora-Bowen
  (Ladruno) — …` comment.
- **`LEDGER_quirks.md`** — record: (a) single-scalar `lch` gives energy- not band-width-
  objectivity for compressive crush; (b) Tier-1 non-associated requires an unsymmetric solver;
  (c) the LogStrain `bᵉ` recovery needs the implicit effective stress under IMPL-EX.
- **`LEDGER_vanilla_files.md`** — only if the kernel out-contract forces a `LogStrainNDMaterial`
  edit (the `σ_eff,implicit` exposure) — mark with a `// Ladruno …` comment.
- **`stamp_headers.py`** — run on every new file (add to its GLOBS).
- **PRs** based on `ladruno`. One logical PR per phase; verify the ADR PR doesn't squash-strand
  the code (the [[feedback_stranded_commits_after_automerge]] trap).

---

## 10. Provenance

Scoped via a 15-agent adversarial workflow (`ladruno-concrete3d-adr-scoping`, run
`wf_2e5eda80-def`, 2026-06-16): 4 ground-truth recon readers (ASDConcrete3D reuse map,
classTags+sibling, kernel pattern, validation conventions) → 5 design panels (surface/return-map,
damage/regularization, robustness, architecture, validation) → 5 red-team critics (1.37M tokens
total). 15 blocking + 22 major findings; all blocking findings are resolved in the decision text
above (§4–§5) and itemized in §8. The synthesis draft's three worst errors — the false
"port IMPL-EX wholesale" premise, the impossible reduce-to-ASDConcrete3D byte gate, and a
phantom `33016` classTag — were caught by the recon/red-team split and corrected here.
