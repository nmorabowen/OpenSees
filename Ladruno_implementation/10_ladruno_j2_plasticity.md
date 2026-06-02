---
title: LadrunoJ2 — combined isotropic + Chaboche kinematic J2 plasticity
project: Ladruno
status: v1 implemented (3D)
priority: high
owner: nmora
tags:
  - implementation
  - material
---

# LadrunoJ2 — rate-independent J2 with combined hardening

## What

A single rate-independent von&nbsp;Mises (J2) `NDMaterial` that supersedes **both**
existing OpenSees J2 materials in one clean class:

- **Yield** on the shifted deviatoric stress **ξ = s − α** (relative stress).
- **Isotropic hardening**: Voce saturation + linear (analytic default), *or* a
  tabulated `σ_y(ε̄ᵖ)` curve (C¹ monotone spline).
- **Kinematic hardening**: **Chaboche** superposed Armstrong–Frederick backstress
  `α = Σₖ αₖ`, `k = 1…N` — architected for arbitrary `N`, **shipping `N = 3`**.
  Per-term limits recover single-AF (`N=1`) and linear Prager (`γₖ = 0`).
- **Integration**: implicit backward-Euler closest-point projection (CPP) with a
  **rigorous coupled** local Newton (Δγ + return direction) and an **analytic
  consistent tangent**. Serves implicit (uses tangent) and explicit (ignores it)
  global solvers with no compromise.
- **IMPL-EX**: *not* in v1, but the state layout is structured so an IMPL-EX path
  can be added later as a flag (see [Risks](#risks--open-questions)).

This is, in effect, **Abaqus `*PLASTIC, COMBINED` (Mises)** implemented as an
implicit-grade OpenSees `nDMaterial`. It is the flagship GREEN 3D material that the
finite-strain wrapper ([[09_finite_strain_material_wrapper]]) lifts to large strain,
and the natural constitutive partner for [[09_ladruno_brick]].

**Not in scope (follow-ups):** pressure-dependent yield (soil/concrete) — different
surface, different material; rate dependence (Perzyna / Johnson–Cook) — leave a hook
only; Ohno–Wang ratcheting threshold — documented future flag; a separate
explicit-only hypoelastic (Jaumann-rate) sibling — a different family entirely.

## Why

J2 is *the* workhorse plastic model in every production FEA code (LS-DYNA, Abaqus)
and the backbone of continuum **seismic steel** modeling. Today the fork has two
materials and neither is adequate:

- **`J2Plasticity`** (Ed Love, `SRC/material/nD/J2Plasticity.*` + dimensional
  specializations) — good Voce+linear isotropic law, robust local-Newton radial
  return, correct analytic consistent tangent, true `getInitialTangent`. **But zero
  kinematic hardening** ⇒ no Bauschinger, useless for cyclic.
- **`SimplifiedJ2`** (Gu & Qiu, `SRC/material/nD/SimplifiedJ2.*`) — has *linear*
  iso + kinematic, closed-form λ. **But carries real defects**: `getInitialTangent`
  returns the *current* (plastic) tangent (wrong for modal / initial-stiffness
  analyses); `revertToStart`/`revertToLastCommit` are empty stubs; several
  `exit(-1)` kernel-killers (the [[../.claude/.../MEMORY|known MPCO-style gotcha]]);
  `setParameter`/`updateParameter`/`Print` are stubs; 2D path is a plane-strain
  slice, not true plane stress; linear hardening only ⇒ over-predicts cyclic
  hardening, no ratcheting saturation, no mean-stress relaxation.

The missing capability is the **combination**: nonlinear isotropic **and**
*nonlinear* kinematic (Chaboche) hardening, done implicit-grade. That combination is
exactly what reproduces stable hysteresis loops and the correct ratcheting rate
under reversed cyclic loading — the response a BRB, a panel zone, or a continuum
steel coupon actually shows.

## Where

- **New code:** `SRC/material/nD/ladrunoJ2/`
  - `LadrunoJ2.{h,cpp}` — the material (3D canonical; 6-component Voigt internally).
  - `LadrunoJ2Kernel.h` — header-only return-map + tangent kernel (so the
    finite-strain wrapper and a future explicit sibling can reuse the *exact* inner
    map verbatim, the same pattern as `LogStrainKernel.h` / `EnergyBalanceKernel.h`).
  - `OPS_LadrunoJ2.cpp` — `OPS_*` parser (Tcl/Python registration).
  - `CMakeLists.txt`.
- **Modify (vanilla — log in [[LEDGER_vanilla_files]]):**
  - `SRC/classTags.h` — add `ND_TAG_LadrunoJ2 33011` (band ≥33000; materials at
    3301x after `LogStrainNDMaterial`=33010).
  - `SRC/api/elementAPI` material registry + `FEM_ObjectBroker` — register the
    `OPS_*` factory and broker entry.
  - `SRC/material/nD/CMakeLists.txt` — add the subdir.
- **Reference (copy patterns from):**
  - `SRC/material/nD/J2Plasticity.cpp` — local-Newton return-map skeleton,
    consistent-tangent assembly, `doInitialTangent`, the `IIdev`/`IbunI` machinery.
  - `SRC/material/nD/SimplifiedJ2.cpp` — Voigt backstress bookkeeping (what to keep,
    what to fix).
  - `SRC/material/nD/LogStrainNDMaterial.cpp` — modern Ladruno material idioms
    (kernel split, clean `sendSelf`/`recvSelf`, no `exit`).

## How

### Constitutive model (3 notations)

**Elastic / split.** Additive small-strain split `ε = εᵉ + εᵖ`. Volumetric–deviatoric
decoupling: `p = K·tr(ε)`, deviatoric trial `sᵗʳ = 2G·(e − eᵖₙ)`, `e = dev(ε)`.

**Yield surface** (shifted / relative stress **ξ = s − α**):

- Direct: `f = ‖ξ‖ − √(2/3)·σ_y(ε̄ᵖ)`
- Index: `f = √(ξ_{ij}ξ_{ij}) − √(2/3)·σ_y`,  `ξ_{ij} = s_{ij} − α_{ij}`
- Voigt: `f = ‖ξ‖ − √(2/3)·σ_y`, with `ξ = {ξ₁₁,ξ₂₂,ξ₃₃,ξ₁₂,ξ₂₃,ξ₃₁}` and the
  norm taken with the factor-2 weights on the off-diagonals (engineering-shear care).

**Flow (associative)** `ε̇ᵖ = γ̇·n`, `n = ξ/‖ξ‖`; iso rate `ε̄̇ᵖ = √(2/3)·γ̇`.

**Isotropic law** (selectable):
- `voce_linear`: `σ_y(ε̄ᵖ) = σ₀ + Q∞·(1 − e^{−b·ε̄ᵖ}) + H_iso·ε̄ᵖ`
  (Abaqus form; `Q∞=0,H_iso=0` ⇒ perfect; `Q∞=0` ⇒ linear; supersets Ed Love).
- `tabulated`: `σ_y` from `{ε̄ᵖᵢ, σᵢ}`, stored as a **monotone C¹ spline** so the
  tangent stays smooth (raw piecewise-linear kinks slow implicit Newton at
  breakpoints; LS-DYNA MAT_024 gets away with kinks only because it's explicit).

**Kinematic law** — three selectable modes on `α = Σₖ₌₁ᴺ αₖ`:

- **`af` (default)** — Chaboche superposed Armstrong–Frederick:
  `α̇ₖ = (2/3)·Cₖ·ε̇ᵖ − γₖ·αₖ·ε̄̇ᵖ` (dSNPO 6.223). `γₖ=0` ⇒ linear Prager term
  (`= SimplifiedJ2`); `N=1` ⇒ single AF; `N≥2` ⇒ full Chaboche. Each term
  self-saturates to `αₖ,sat = (2/3)Cₖ/γₖ`. **The recall term gives dynamic
  recovery / ratcheting / mean-stress relaxation** — the cyclic realism seismic
  steel needs.
- **`prager_nl` (oracle + cheap option)** — dSNPO §7.6 *nonlinear Prager*:
  `α̇ = (2/3)·Hᵏ(ε̄ᵖ)·ε̇ᵖ` via a sampled kinematic-hardening *stress* curve
  `ᾱ(ε̄ᵖ)` (dSNPO 6.224, 7.185). **Radial**, single scalar Δγ, matches HYPLAS
  `SUVMMX` / Box 7.5 / tangent eq 7.213 **bit-for-bit** ⇒ our gold validation
  oracle. Bounds the backstress but does **not** relax on reversal (no
  ratcheting) — offered as a robust cheap mode.
- **`prager_lin`** — classic linear Prager (`= af` with one term, `γ=0`).

> **AF ≠ nonlinear Prager.** dSNPO §7.6 (the model that otherwise *is* LadrunoJ2)
> implements `prager_nl`, **not** AF — it defers AF-coupled-with-isotropic to §12.3
> (Lemaitre damage) and §14.11 (finite strain). So we follow the book's *framework*
> (§7.2 general return map, §7.4.4 general consistent tangent, Box 7.5 structure)
> but the AF recall term comes from §6.6.4 theory + Kobayashi–Ohno (2002) + Abaqus.

### Backward-Euler return map (the kernel) — **scalar** Newton

Backward-Euler update of each backstress term and the deviatoric stress:

```
αₖⁿ⁺¹ = [ αₖⁿ + (2/3)·Cₖ·Δγ·n ] / (1 + γₖ·Δγ)        (closed per-term denominator)
sⁿ⁺¹  = sᵗʳ − 2G·Δγ·n ,    n = ξⁿ⁺¹/‖ξⁿ⁺¹‖ ,   ξⁿ⁺¹ = sⁿ⁺¹ − Σₖ αₖⁿ⁺¹
```

**The whole system reduces to ONE scalar equation in Δγ** (Kobayashi–Ohno 2002 — the
AF generalization of dSNPO's radial reduction 7.190–7.192). Substituting the
backstress update into `ξⁿ⁺¹`:

```
ξⁿ⁺¹ = M(Δγ) − θ(Δγ)·n ,   with
M(Δγ) = sᵗʳ − Σₖ αₖⁿ/(1+γₖΔγ)              (known tensor, fn of Δγ)
θ(Δγ) = 2G·Δγ + Σₖ (2/3)CₖΔγ/(1+γₖΔγ)       (known scalar)
```

Since `n ∥ ξⁿ⁺¹`, this forces **`n = M(Δγ)/‖M(Δγ)‖`** and `‖ξⁿ⁺¹‖ = ‖M(Δγ)‖ − θ(Δγ)`,
both closed-form in the scalar Δγ. The consistency residual is therefore scalar:

```
R(Δγ) = √(2/3)⁻¹·(‖M(Δγ)‖ − θ(Δγ)) − √(2/3)·σ_y(ε̄ᵖₙ + √(2/3)Δγ) = 0
```

solved by Newton, `dR/dΔγ` analytic via `dM/dΔγ = Σₖ αₖⁿ γₖ/(1+γₖΔγ)²`.

**Radiality is NOT required.** dSNPO's `ξⁿ⁺¹ ∥ ξᵗʳ` (7.192) is just the `γₖ=0`
special case of `n = M(Δγ)/‖M(Δγ)‖`. For AF the direction merely *shifts* with Δγ
(because each `αₖⁿ` is weighted by `1/(1+γₖΔγ)`) — but it stays a **closed-form
direction**, so this is a *scalar* Newton, **not** a coupled tensor solve. (This
corrects the earlier "≤7-dim coupled" framing — that was overly pessimistic.) Seed
Δγ from the radial estimate and bracket `Δγ ∈ [0, Δγ_max]` for reversal robustness.

**Consistent tangent.** dSNPO 7.213 structure (= `J2Plasticity`'s `c2·N⊗N + c3·(I_dev
− N⊗N) + K·I⊗I`):

```
D^alg = 2G(1 − Δγ·3G/q)·I_dev + 6G²(Δγ/q − 1/(3G+Hᵏ+Hⁱ))·N̄⊗N̄ + K·1⊗1
```

plus, for the AF modes, the `∂n/∂εᵗʳ` term arising because `n = M(Δγ)/‖M‖` depends
on Δγ (hence on ε) — analytic. In the `prager_nl` mode it collapses *exactly* to
7.213 (a bit-for-bit oracle check); in the `N=0` limit it must reduce to
`J2Plasticity` to **1e-12**.

### Public API

```tcl
# analytic Voce+linear iso, N kinematic backstresses
nDMaterial LadrunoJ2 $tag  $K $G  -iso voce $sig0 $Qinf $b $Hiso \
                                   -kin $N  $C1 $g1  $C2 $g2  $C3 $g3 \
                                   <-rho $rho>
# tabulated isotropic curve
nDMaterial LadrunoJ2 $tag  $K $G  -iso table {e0 s0  e1 s1 ...} \
                                   -kin 1  $C1 $g1
```

- Linear Prager term ⇒ pass `g_k = 0`. Single AF ⇒ `-kin 1 C1 g1`. Pure isotropic
  ⇒ `-kin 0`.
- `getCopy("ThreeDimensional"|"PlaneStrain"|"AxiSymmetric"|"PlaneStress"|"PlateFiber")`
  returns the dimensional view; **PlaneStress** needs the nested Δε₃₃ secant
  sub-iteration (LS-DYNA Theory eq 19.3.18–22 / Simo–Hughes) — flagged as the one
  non-trivial specialization.

### State (committed) and `sendSelf`/`recvSelf`

Per Gauss point: `εᵖ` (6), `ε̄ᵖ` (1), and `αₖ` (6·N). All serialized — `SimplifiedJ2`
correctly sends its single backstress; we extend to N. Hygiene fixes baked in: real
`getInitialTangent` (elastic), real `revertToStart`/`revertToLastCommit`, **no
`exit(-1)`** (return −1, opserr, Ladruno house style), working `setParameter`
(K,G,σ₀,Q∞,b,Cₖ,γₖ,ρ) for FE sensitivity, real `Print`.

### Testing / validation matrix

| # | Test | Reference / oracle |
|---|---|---|
| V0 | elastic round-trip (`f<0`) | `D^alg == D^e`, σ = D^e:ε |
| V1 | **reduce-to-J2Plasticity**: `N=0` (pure iso, Voce+linear) | stress+tangent reduce term-by-term; numeric ~**1e-7** (bounded by J2Plasticity's internal `γ*=(1−1e-8)` fudge, not 1e-12) |
| V2 | reduce-to-SimplifiedJ2: `N=1, γ=0` (linear Prager) + linear iso | analytic / `SimplifiedJ2` cross-check |
| V3 | uniaxial monotone | analytic Voce+linear backbone |
| V4 | uniaxial reversed cycle | Bauschinger offset, single-AF saturation `2·αₛₐₜ` |
| V5 | strain-controlled cyclic, N=3 | stable hysteresis loop shape (Chaboche calibration) |
| V6 | uniaxial ratcheting (stress-controlled) | AF over-prediction signature (documents Ohno–Wang need) |
| V7 | consistent-tangent FD check | `‖D^alg − D^FD‖ < tol` at several plastic states |
| V8 | objectivity / frame (via finite wrapper) | rigid rotation ⇒ no stress change |
| **V9** | **`prager_nl` mode vs dSNPO Box 7.5** (nonlin Prager + nonlin iso, sampled curves) | **bit-for-bit** vs hand-coded HYPLAS `SUVMMX`/`CTVMMX` (tangent 7.213) — the book oracle |

Smoke (L0): single-element LadrunoBrick + LadrunoJ2 push-pull, Zone-A pytest.

```cpp
// LadrunoJ2Kernel.h — reused verbatim by the finite-strain wrapper
struct LadrunoJ2State { double epsP[6]; double ebarP; double alpha[3][6]; };  // alpha[N][6]
// returnMap(eps, state_n) -> {stress[6], Dalg[6][6], state_np1}
//   SCALAR Newton on dGamma; direction n = M(dGamma)/||M(dGamma)|| (Kobayashi-Ohno);
//   analytic tangent (dSNPO 7.213 + dn/deps term)
```

## Risks / open questions

> [!question]
> **Reversal robustness of the scalar solve.** The Δγ Newton is scalar
> (Kobayashi–Ohno), but `M(Δγ)` can shrink near a reversal. Bracket
> `Δγ ∈ [0, Δγ_max]`, seed from the radial (`n=Mᵗʳ/‖Mᵗʳ‖`) estimate, fall back to
> bisection if `‖M(Δγ)‖ − θ(Δγ) ≤ 0`. *(No coupled tensor solve — superseded.)*

> [!question]
> **PlaneStress: nested vs projected.** dSNPO gives two routes — nested Δε₃₃
> iteration (§9.2.3) or the cleaner *plane-stress-projected* sub-space algorithm
> (§9.4, `SUVMPS`/`CTVMPS`). Prefer the projected route over the LS-DYNA secant loop.
> Build/validate 3D + PlaneStrain + AxiSymm first; PlaneStress is a second PR.

> [!question]
> **IMPL-EX hook (structure-only in v1).** Lay out `LadrunoJ2State` + the kernel so
> the *committed* `Δγⁿ` is stored and the kernel can expose an extrapolated
> `Δγ̃ⁿ⁺¹ = Δγⁿ·(Δtⁿ⁺¹/Δtⁿ)` with a **constant SPD** tangent — but ship no IMPL-EX
> code path. It earns its keep later for explicit dynamics and any future
> softening/damage extension (where implicit Newton stalls on indefinite tangents);
> it is *not* the default because it is O(Δt) and least accurate exactly at the
> cyclic reversals we care about.

- **Consistent tangent ↔ V1 1e-12.** The tangent must collapse to `J2Plasticity`'s
  closed form (= dSNPO 7.213 with `Hᵏ=0`) in the `N=0` limit; this regression is the
  cheap proof the iso path is untouched.
- **No `exit()` anywhere** — fork rule; `SimplifiedJ2`'s `exit(-1)` calls are the
  anti-pattern being corrected.
- **Build-control:** new `LEDGER_implementations` row (classTag 33011); `// Ladruno`
  comments on every `classTags.h` / broker / registry edit → `LEDGER_vanilla_files`;
  banner line in `Ladruno_scripts/banner_features.txt` → `patch_banner.py`.

## Implementation log

### v1 shipped (2026-06-01, PR #82) — ThreeDimensional
- `SRC/material/nD/LadrunoJ2.{h,cpp}` + `OPS_LadrunoJ2` parser; classTag 33011;
  wired into `classTags.h`, `material/nD/CMakeLists.txt`, `FEM_ObjectBrokerAllClasses`,
  `OpenSeesNDMaterialCommands`. Banner line + 3 ledgers updated.
- **Self-contained class** (does NOT inherit `J2Plasticity`); internal 3×3-symmetric
  tensors as 6 tensor-components `{00,11,22,01,12,02}`; tangent assembled in the
  `J2ThreeDimensional` rank-4→6×6 mapping (so the N=0 reduction is bit-faithful).
- Return map exactly as designed: **scalar Newton on Δγ**, `n=M(Δγ)/‖M(Δγ)‖`.
  AF backstress update `αₖ=(αₖⁿ+⅔CₖΔγ n)/(1+√⅔γₖΔγ)`.
- **Verified analytically** that residual + tangent reduce term-by-term to
  `J2Plasticity` at N=0 (the `2G·β1·IIdev + βNN·n⊗n` coeffs equal Ed Love's
  `2G+c3` / `c2−c3`). Confirmed numerically (V1 ~1e-7, bounded by J2Plasticity's
  internal `γ*=(1−1e-8)` fudge — NOT 1e-12 as the matrix optimistically said).
- Build: full from-scratch OpenSeesPy build green; `LadrunoJ2.cpp` compiled clean.
- Tests `tests/test_ladrunoJ2_material.py` (single stdBrick, 1/8-symmetry uniaxial,
  displacement-controlled): **5/5 pass** — elastic; reduce-to-J2Plasticity;
  monotonic linear-kin≡iso (pins the `(2/3)C` scaling); Bauschinger divergence;
  AF saturation → `σ_y0+C/γ`. The `(2/3)Cₖ` numerator ⇒ standard Chaboche `Cₖ,γₖ`.
- Re-landed after #82 stranding via **PR #87** (cherry-pick onto fresh `ladruno`);
  + adversarial-review hardening (‖M‖→0 guard restored from J2Plasticity,
  stress-scaled tolerance) + 3D mixed-shear test (battery → 6/6).

### Dimensional views shipped (2026-06-01, follow-up PR)
- **All five `getType()` views in one class** via a `dim` mode + `vmap[]` index
  table into the 6-comp tensor; the verified 3D `integrate()` is unchanged.
  PlaneStrain `{00,11,01}`, AxiSymmetric `{00,11,22,01}`, PlateFiber
  `{00,11,01,12,20}`, PlaneStress `{00,11,01}`.
- **PlaneStress / PlateFiber** enforce `σ₂₂=0` by a nested Newton on `eps₂₂`
  (`strain6[2] -= σ₂₂/Dtan[2][2]`, dSNPO §9.2.3 route) then **static condensation
  of the 33-dof done in `Dtan[6][6]`** (`Dtan[I][J] -= Dtan[I][2]Dtan[2][J]/Dtan[2][2]`),
  mirroring `J2PlaneStress`/`J2PlateFiber`. Committed `eps₂₂` carried in
  `sendSelf`/`revert*`. Member-sized return buffers replace the size-6 statics.
- Tests: **8/8 pass** — added `PlaneStrain` + `PlaneStress` quad reduce-to-J2Plasticity
  (single FourNodeQuad, mixed in-plane load incl. shear, into the plastic regime;
  matches `J2Plasticity` on disps + all GP stresses → validates the reduced mapping
  AND the condensation against the proven upstream specializations).
- **Still deferred**: tabulated isotropic; `prager_nl` oracle mode (dSNPO Box 7.5,
  V9); IMPL-EX code path; `setResponse` for backstress/`ε̄ᵖ`; AxiSymm/PlateFiber
  element-level tests (validated by construction — same machinery as PlaneStrain/Stress).

### Decisions locked (2026-06-01, design session)
- **Kinematic = Chaboche AF, design for arbitrary N, ship N=3** (`af` mode). Recovers
  AF (N=1) and linear Prager (γₖ=0). Plus a **`prager_nl` mode = dSNPO §7.6
  nonlinear-Prager** as the book-exact validation oracle + cheap radial option.
- **Return map = SCALAR Newton on Δγ** (Kobayashi–Ohno 2002): direction
  `n = M(Δγ)/‖M(Δγ)‖`, `M(Δγ)=sᵗʳ−Σαₖⁿ/(1+γₖΔγ)`. **Supersedes the earlier
  "coupled ≤7-dim" plan** — AF does *not* need a tensor solve; non-radiality is just
  a Δγ-dependent direction. dSNPO's radial reduction (7.192) is the `γ=0` case.
- **dSNPO alignment:** §7.6 (von Mises nonlinear mixed hardening) *is* this model but
  implements **nonlinear Prager, not AF** (AF deferred to §12.3 damage / §14.11 finite
  strain). We follow the book's framework (§7.2 return map, §7.4.4 + 7.213 tangent,
  Box 7.5 structure); the AF recall term is §6.6.4 theory + Kobayashi–Ohno. Plane
  stress = §9.4 projected (preferred) or §9.2.3 nested. Finite-strain kinematic lift
  = §14.11 via [[09_finite_strain_material_wrapper]].
- **IMPL-EX = structure-only hook** in v1; no code path yet.
- **Isotropic = Voce+linear (analytic default) + tabulated** (sampled curve handled
  dSNPO-style via 7.196 direct-curve trick; C¹ spline for tangent smoothness).
- classTag **33011**; files under `SRC/material/nD/ladrunoJ2/`; kernel split for
  finite-strain reuse.
