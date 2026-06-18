---
title: "LadrunoRCConcrete — nonlinear RC layer-shell constitutive material"
project: Ladruno
type: reference / user guide
status: shipped — Phase 1 / 2a / 2b.1 / 2b.2a / 2b.2b / 2b.2c.1-4a + IMPL-EX (4a) + Phase 3 (3a tension stiffening + 3b crack-band lch regularization); cyclic material physics COMPLETE. Remaining: quantitative Tran–Wallace experiment match, Phase 3b structural rotated-mesh objectivity gate (staged), 4b (finite view), 5 (solid-shell). MERGED PRs #155 #192 #239 #245 #246 #253 #263 #266 #273 (3a); Phase 3b in flight
classTag: 33015 (ND_TAG family)
material: "nDMaterial LadrunoRCConcrete — 3D / PlaneStress / PlateFiber views"
related:
  - "[[19_ladruno_rc_shell_adr]]"
  - "[[LadrunoConcrete3D_user_guide]]"
  - "[[Ladruno_materials_guide]]"
  - "[[LadrunoRebarBuckling_guide]]"
  - "[[LogStrain_guide]]"
tags:
  - material
  - rc-shell
  - mcft
  - plane-stress
  - compression-softening
  - aggregate-interlock
  - crack-band
  - reference
---

# LadrunoRCConcrete — nonlinear RC layer-shell concrete

> [!abstract] One-line summary
> `LadrunoRCConcrete` (nDMaterial, classTag **33015**) is a plastic-damage concrete
> material that keeps `ASDConcrete3D`'s verified spine (spectral effective-stress split,
> dual tension/compression damage, Lubliner–Lee–Fenves biaxial envelope, crack-band `lch`
> regularization) and adds the **four RC physics layers `ASDConcrete3D` lacks**: MCFT
> **compression (transverse-tension) softening** `β`, a **degrading aggregate-interlock
> shear-transfer** law (monotonic bound → cyclic friction-slip → X-cracking + wear),
> **tension stiffening** through the tabulated tension backbone, and a **consistent
> tangent** that honours the new strength couplings. It ships as **one multi-dim class** offering `ThreeDimensional`,
> `PlaneStress`, and `PlateFiber` views, so the **same material** that runs on a 3D brick
> drops into a `LayeredShell` section under `ASDShellQ4` and gives nonlinear RC walls and
> slabs **with zero element/section edit**. Every added layer is a **flag, default OFF** —
> with all flags off the material is trajectory-faithful to `ASDConcrete3D`.

This is the **descriptive reference** for the material: the gap it closes, the physics of
each layer, the command grammar, intended uses, recording, and the honest limits. For the
full decision log (D1–D6), the seam analysis, and the phased implementation plan, see the
source ADR [[19_ladruno_rc_shell_adr]]. For the spine it clones, see
[[LadrunoConcrete3D_user_guide]].

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why a material, not a new element|2. Why a material, not a new element]]
- [[#3 The constitutive stack|3. The constitutive stack]]
- [[#4 The four added physics layers|4. The four added physics layers]]
- [[#5 The multi-dim views and the shell seam|5. The views and the shell seam]]
- [[#6 Command syntax|6. Command syntax]]
- [[#7 Worked examples|7. Worked examples]]
- [[#8 State recording|8. State recording]]
- [[#9 Verification and validation|9. Verification]]
- [[#10 Limitations and boundaries|10. Limitations]]
- [[#11 References|11. References]]

---

## 1. Intended use cases

`LadrunoRCConcrete` is **plane-stress / layered-shell reinforced concrete** — the
constitutive engine for nonlinear analysis of RC walls, slabs, shells, and 2D membrane
panels. The canonical targets:

| Use case | Why this material |
|---|---|
| **Squat / shear-controlled RC walls (the headline)** | The structural shear `V` of an in-plane-loaded wall lives in the **membrane** block (`γxy`/`τxy`), so it is governed by the plane-stress *constitutive* law — not transverse shear, not element technology. The MCFT `β` softening + aggregate-interlock cap directly correct OpenSees' over-prediction of diagonal-strut capacity. |
| **Cyclic RC walls — pinching / hysteresis** | The fixed-crack friction-slip interlock (`-cyclic`) accumulates slip and degrades across reversals — the mechanism a rotating-coaxial model structurally cannot represent. X-cracking + cumulative wear (`-xcrack`) give bidirectional capping and gradual cyclic strength decay. |
| **RC slabs and shells in bending** | The same material in a `LayeredShell` (8–12 layers) under `ASDShellQ4` integrates cracked `σ(z)` through the thickness via the `PlateFiber` view + `σ33=0` condensation. |
| **2D membrane / continuum panels** | The `PlaneStress` view hosts a 2D continuum element directly and is the clean MCFT oracle target. |
| **3D continuum RC** | The `ThreeDimensional` view is directly comparable to `ASDConcrete3D` (the reduce-to-baseline gate) and composes with `LadrunoBrick`. |

> [!warning] When **not** to reach for it
> - **Punching / bearing / true 3D-stress crush** needs the through-thickness `σ33` state
>   a director shell cannot carry — that is the deferred `LadrunoSolidShell` host, not this
>   material on `ASDShellQ4`.
> - **Discrete boundary rebar where buckling matters** is a separate `PlateRebar`
>   ([[LadrunoRebarBuckling_guide]]) layer at the true bar depth — only **smeared web**
>   steel is homogenized inside this kernel for the MCFT composite `ε1`.
> - **Large-rotation cyclic walls:** the stored crack frame is a directional internal
>   variable subject to the dSNPO §14.11 objectivity boundary — use `ASDShellQ4
>   -corotational` (the supported objective route), not the finite-strain material view.

---

## 2. Why a material, not a new element

The deficiency that makes OpenSees over-predict squat-wall capacity and miss cyclic
pinching is **constitutive, not elemental**. `ASDShellQ4` (AGQ6-I membrane + 4-DOF EAS +
MITC4 transverse shear + Hughes–Brezzi drilling + optional corotational) already equals the
LS-DYNA `ELFORM=16` fully-integrated assumed-strain shell; `LayeredShellFiberSection`
already does standard director-stack integration. Rebuilding either re-treads solved
element technology for **zero new physics**.

So the keystone decision (ADR D1/D2): **keep the element and section at the frontier,
deliver the missing physics as a material** that drops into the 5-component `PlateFiber`
layer view the section already requests. The element computes the 8-component generalized
strain `E = [exx eyy gxy | kxx kyy kxy | gxz gyz]`, the section forms the per-layer plane
state `e(z) = e_m − z·κ` and hands each layer a `PlateFiber` strain — and
`LadrunoRCConcrete` answers on that seam, untouched element, untouched section.

> [!note] One class, many views (the `LadrunoJ2` / `ASDConcrete3D` pattern)
> Rather than three sibling classes, the material is **one `LadrunoRCConcrete`** whose
> view is selected by `getType()` / `getCopy(type)`. The in-plane return map is shared
> byte-for-byte between the `PlaneStress` and `PlateFiber` views; the 3D view is the
> reduce-to-`ASDConcrete3D` baseline. Physics lives in the header-only, OpenSees-free
> `LadrunoRCKernel.h` (namespace `ladruno_rc_kernel`), numpy-oracle-testable before any
> OpenSees link — cloned from the proven `LadrunoJ2Kernel.h`.

---

## 3. The constitutive stack

### 3.1 The inherited spine (cloned from `ASDConcrete3D`)

The base law is `ASDConcrete3D`'s plastic-damage spine, cloned **verbatim** (see
[[LadrunoConcrete3D_user_guide]] for the full theory):

- **Spectral effective-stress split** `σ̃ = σ̃⁺ + σ̃⁻` into tensile / compressive cones via
  the eigenprojectors of `σ̃`.
- **Dual damage** `dt` (tension) and `dc` (compression), each driven by its own
  equivalent-strain measure on a **tabulated backbone** (`-Te/-Ts/-Td`, `-Ce/-Cs/-Cd`).
- **Lubliner–Lee–Fenves biaxial envelope** (`fb = 1.16 fc`, `Kc`, `γ = 3(1−Kc)/(2Kc−1)`)
  raising the effective-compressive measure in the tension–compression quadrant.
- **Crack-band (Bažant–Oh) `lch` regularization** so dissipated energy is mesh-objective.
- **Crack-closure** spectral reassembly and an optional IMPL-EX path.

> [!important] Two spine-cloning gotchas (in [[LEDGER_quirks]])
> (1) The equivalent-strain measure needs the `/E` division — a `β`-ratio test cannot
> catch its omission; only an **absolute-stress** test can. (2) The effective-stress
> backbone `q` is **E-consistent by construction** (`buildBackbone` mirrors the
> `ASDConcrete3D` `HardeningLaw` constructor + `adjust()`), **not** `q = y/(1−d)` on raw
> points.

### 3.2 Backbones are tabulated (the `ASDConcrete3D` `HardeningLaw` format)

Compression and tension backbones are supplied as point lists — strain `e`, stress `s`,
and (optional) damage `d` — exactly like `ASDConcrete3D`:

```
-Ce e0 e1 e2 …   -Cs s0 s1 s2 …   [-Cd d0 d1 d2 …]     # compression (>=2 points)
-Te e0 e1 e2 …   -Ts s0 s1 s2 …   [-Td d0 d1 d2 …]     # tension (>=2 points)
```

`-Ce`/`-Cs` and `-Te`/`-Ts` are **required and must be equal length** (≥2). `-Cd`/`-Td`
are optional (pad to zero ⇒ pure plastic, no damage). The peak of `-Ts` is the tensile
strength `ft`; the peak of `-Cs` is `fc'`.

---

## 4. The four added physics layers

Each layer is a **flag, default OFF**. With every flag off, `LadrunoRCConcrete` reduces to
the cloned `ASDConcrete3D` spine (trajectory-faithful to ~1e-6).

### 4.1 MCFT compression softening — `-beta`

The verified gap: `ASDConcrete3D` never reduces compressive strength by transverse tensile
cracking. The Modified Compression Field Theory (Vecchio & Collins 1986) closes it with a
softening factor on the compressive **strength**:

$$\beta = \frac{1}{0.8 + 170\,\varepsilon_1} \le 1, \qquad \frac{\partial\beta}{\partial\varepsilon_1} = -170\,\beta^2$$

where `ε1` is the average principal **tensile** (cracking) strain transverse to the strut,
and the realised compressive peak becomes `|σc| = β·fc'`.

> [!important] β scales the STRENGTH axis — the blocking correctness gate (ADR D4)
> `β` must scale the **stress/strength axis** (the effective-compressive backbone value
> `q` fed into `dc`), **never** the *strain abscissa* fed to the equivalent-strain measure.
> Scaling the abscissa merely slides the lookup along a fixed non-linear curve, giving a
> realised peak `hc(β·xc)` ≠ `β·fc'`. This is proven **three ways** (numpy oracle,
> standalone g++, end-to-end OpenSees): hold confined compression, sweep `ε1`, assert
> realised `|σc| = β(ε1)·fc'` exactly. SFI-MVLEM=FSAM ≠ MCFT — do not mistake the FSAM
> material for an MCFT oracle.

> [!note] Avoiding double-count with the Lubliner envelope — `-lublinerReduced`
> The Lubliner envelope already raises the effective-stress measure in the
> tension–compression quadrant where struts live. Stacking `β` on top double-penalises
> that quadrant. `-lublinerReduced` reduces the t–c interaction so transverse tension is
> counted **once**; `β` and `Kc` are meant to be calibrated **jointly** against a
> tension–compression panel battery (Kupfer + Vecchio–Collins).

### 4.2 Aggregate-interlock shear transfer — `-interlock` (Phase 2a, monotonic)

A pure rotating-crack model is coaxial — it has no independent crack-plane shear stiffness.
To carry (and *bound*) cracked shear, the material freezes the in-plane crack **normal** at
first crossing of `crackStrain`, and thereafter **clips** the smeared (damage-reduced)
crack-plane shear `τ_sm = m_σ·σ_ip` to the MCFT crack-shear limit:

$$v_{ci,\max} = \frac{0.18\sqrt{f_c'}}{0.31 + 24w/(a_g+16)}, \qquad w = \langle\varepsilon_n\rangle\, s_\theta$$

`w` = crack width (Macauley of the strain normal to the **frozen** crack × crack spacing
`s_θ`), `a_g` = aggregate size (`-agg`, default 16). 2a is a **bound on the existing
shear**, not a replacement with bare-elastic `G·γ` — below the cap the stress is unchanged
(continuous, no discontinuity at cracking). The consistent tangent below the cap is the
baseline; at the cap it is the exact rank-1 removal
`Dtan -= (1−betaSrMin)·m_ε⊗(m_σ·Dtan)`, exact because `m_σ·m_ε = (c²+s²)² = 1`.

> [!tip] Interlock only engages under NON-proportional loading
> Under *proportional* shear the crack forms on the principal plane, so the crack-tangential
> shear `g_nt = 0` and interlock is inert. It activates when the loading path rotates the
> principal directions away from the frozen crack — i.e. in real panels / non-radial paths.

### 4.3 Cyclic friction-slip — `-cyclic` (Phase 2b.1)

`-cyclic` (a sub-flag of `-interlock`) promotes the crack-plane shear to an **independent
incremental friction-slip state**:

$$\tau_{cr} = \mathrm{clamp}\bigl(\tau_{cr}^{c} + G\,(\gamma_{nt} - \gamma_{nt}^{c}),\; \pm v_{ci,\max}(w)\bigr), \quad G = \frac{E}{2(1+\nu)}$$

with the slip origin frozen at crack capture. The crack width `w` is now **reversible**
(current opening, not the monotone max) so crack closure raises the cap (capacity
recovery). State = two committed scalars `{τcr, γcr}`. It diverges from FSAM in three
deliberate respects (kept for continuity with the 2a bound): the cap is MCFT `v_ci,max(w)`
not Coulomb `μ·f_n`; the slip stiffness is the full elastic `G` not `0.4·Ec`; an open crack
still transfers shear up to the width-degraded cap rather than collapsing to zero.

> [!note] Pinching is a PANEL effect
> At a material point with a fixed crack and constant normal strain the loop is *fat*; a
> pinched waist needs the crack to open/close during the cycle (principal rotation), which
> only happens in a meshed wall. The material gives the **ingredients** (reversal re-capping,
> friction unload at stiffness `G`, capacity recovery on closure); the waist emerges at the
> structural scale.

### 4.4 Second crack (X-cracking) + cyclic wear — `-xcrack` (Phase 2b.2b)

`-xcrack` (implies `-cyclic` ⇒ `-interlock`) adds an orthogonal crack 2 (normal `(−s,c)`)
that freezes when the strain normal to *that* plane reaches `crackStrain`. The interlock
cap then uses the **governing (most-open) crack's opening**
`v_ci,max(⟨max(εn1, εn2)⟩·s_θ)` — so the **reverse** shear direction is also capped (single
crack only capped one direction). There is **one** friction state and **one** governing cap
— no double-count.

Cyclic wear is an irreversible Archard-style knockdown driven by the **cumulative crack
sliding distance** `slipCum` (not peak slip, which saturates in cycle 1 at constant
amplitude):

$$k_n = \max\bigl(\text{degMin},\; 1 - \text{degKappa}\cdot\mathrm{clamp}(\text{slipCum}/\text{degSlipRef}, 0, 1)\bigr), \qquad v_{ci,\max} \leftarrow k_n\, v_{ci,\max}$$

Driven by the **committed** `slipCum` ⇒ tangent-neutral (the benign 2b.1 tangent is
unchanged). This produces gradual cyclic strength decay (verified caps `2.90 → 2.50 → 2.10`
over three cycles). Knobs: `-degKappa` (0.5), `-degSlipRef` (0.01), `-degMin` (0.1).

### 4.5 Crack-shear retention curve — `-shearRetention` (Phase 2b.2c)

The cyclic crack-plane **slip stiffness** `G_slip` in the friction predictor
`τ_cr = clamp(τ_cr^c + G_slip·Δγ_nt, ±v_ci,max(w))` is selectable. The `v_ci,max(w)` **cap is
unchanged** in every mode — only the slip stiffness changes. All modes reduce to `mcft`.

| `-shearRetention <mode>` | `G_slip` | meaning |
|---|---|---|
| `mcft` (default) | `G = E/2(1+ν)` | full elastic slip stiffness — the shipped DSFM-with-slip law |
| `const` | `μ·G` | classic constant shear retention (Rots/Červenka); `μ = -shearRetFactor` ∈ (0,1], default 0.4 |
| `dsfm` | `G·(0.31/denom)` | DSFM width-degraded: softens as the crack opens (`denom = 0.31 + 24w/(a_g+16)`); `→ G` at `w=0` |
| `rots` | — | rotating-coaxial: the fixed-crack shear is skipped entirely ⇒ identical to `-interlock` OFF |

`const`/`dsfm` imply `-interlock -cyclic`; `rots` implies `-interlock`. Omitting `-shearRetention`
(or `mcft`) is **2b.2b byte-identical**. Because `const`/`dsfm` parametrize the *slip* stiffness, they
only bite on the `-cyclic` path; on the 2a monotone-clip path they are inert (same `v_ci,max` plateau).

> [!example] Softer, width-aware cyclic crack shear
> ```tcl
> # constant 40%-of-G crack-shear retention (FSAM-like 0.4·Ec):
> nDMaterial LadrunoRCConcrete 1 $E $nu -Ce ... -Cs ... -Te ... -Ts ... \
>     -crackSpacing 50.0 -shearRetention const -shearRetFactor 0.4
> # DSFM width-degraded slip stiffness:
> nDMaterial LadrunoRCConcrete 2 $E $nu -Ce ... -Cs ... -Te ... -Ts ... \
>     -crackSpacing 50.0 -shearRetention dsfm
> ```

### 4.6 IMPL-EX robustness — `-implex`

Softening plastic-damage produces an **indefinite consistent tangent** on the softening
branch, so a vanilla Newton solve of a cyclic RC wall stalls (an `nDMaterial` with no
robustness escape cannot be driven through reversals of a meshed softening wall). IMPL-EX
(Oliver et al. 2008) is the escape: each step the damage thresholds `xt, xc` and the MCFT
`β` are **frozen** by an explicit extrapolation of the committed history,

$$x_{\text{ext}} = x_n + t_f\,(x_n - x_{n-1}), \qquad t_f = \frac{\Delta t_{n+1}}{\Delta t_n}\,\alpha,$$

so the returned tangent is the **secant** `W_B·C0` with no softening-rate term and no `β`
cross-term — it removes the dominant indefinite contribution and lets Newton converge. The
**true implicit** thresholds are recomputed at commit to advance the state and to measure
the IMPL-EX error (`implexError` response).

> [!important] What is and isn't frozen (honest scope)
> The **damage** (the softening source) is frozen — that is the robustness that matters.
> The spectral **projectors** PT/PC are **not** frozen here (they are recomputed from the
> live stress each call), consistent with this material's fixed-projector-secant philosophy
> (the implicit tangent also omits `∂P/∂ε`). So the tangent is not perfectly strain-constant
> under a rotating principal axis; a full `PT_commit` freeze for that case is a scoped
> follow-up. IMPL-EX is **exact on proportional paths** (the threshold grows linearly, so the
> extrapolation is exact — error ~1e-12) and lags `O(Δt)` at rate changes (e.g. cracking
> onset).

> [!warning] Static analysis: the time factor is guarded
> In a **static** analysis the "time" is the load factor, whose increment is erratic (and
> resets at `loadConst`), so the raw `Δt_{n+1}/Δt_n` would detonate the extrapolation. The
> material guards `t_f` (falls back to `α` for non-positive/non-finite increments, clamps to
> `2α`), which makes static IMPL-EX use a uniform `t_f≈α` — supply roughly uniform load steps.

`-implex` is **default OFF** ⇒ the material is the fully-implicit baseline, bit-identical.
Combine with `-numericalTangent` and the forward-difference path is automatically bypassed
(it is invalid under IMPL-EX). Pairs naturally with the cyclic interlock flags for cyclic
wall analysis.

### 4.7 Tension stiffening — `-tensStiff {vc|cm}` (Phase 3a)

Between cracks, bonded reinforcement carries part of the tension, so the **average** concrete
tensile stress stays **above** the bare fracture-energy softening curve. Phase 3a adds this as
an opt-in **stress floor** on the live in-plane principal tensile axis `p1`:

$$\sigma_{ts}(\varepsilon_1)=\begin{cases}\dfrac{f_t}{1+\sqrt{c\,\varepsilon_1}} & \texttt{vc}\ \text{(MCFT / Bentz)}\\[2ex]\dfrac{\alpha\,f_t}{1+\sqrt{500\,\varepsilon_1}} & \texttt{cm}\ \text{(Collins–Mitchell)}\end{cases}$$

with `f_t` = the tension-backbone peak and `ε1` the **composite** (reinforced) membrane
principal tensile strain — the *same* `ε1` the MCFT `β` uses (perfect-bond strain
compatibility), so it does not need a separate steel layer's strain. The floor injects
`Δ = σ_ts(ε1) − n^T σ n` along `p1` (rank-1, only when `Δ>0`) so the principal tensile stress
is pinned **up** to `σ_ts`; it never lowers the bare stress. It is active **only post-crack**
(`ε1 ≥ ε_cr`), so the elastic pre-crack branch is untouched. Equibiaxial (degenerate `p1`)
floors **both** in-plane normals to `σ_ts`. The consistent tangent adds the `dσ_ts/dε1` and
the `−d(n^Tσn)/dε` pinning terms (fixed-projector secant, omitting `dp1/dε`, like the `β`
tangent); under IMPL-EX the cross-term is dropped (frozen-`ε1` secant).

`-tensStiff` is **default OFF** ⇒ bit-identical baseline. Knobs: `-tensStiffC $c` (vc
coefficient, default 500, must be `>0`) and `-tensStiffAlpha $a` (cm scale, default 1); `vc`
with `c=500` and `cm` with `α=1` are the same curve.

> [!warning] Monotonic-scope (v1)
> `σ_ts` is a pure function of the **live** `ε1` (no `ε1max` memory). Because `σ_ts`
> *decreases* with `ε1`, on **unloading** the floor **re-inflates** (tracks `σ_ts(live ε1)`
> back up) — correct on the monotone loading branch (the slab / distributed-reinforcement
> pushover use case), but **not** a hysteretic cyclic-tension model. **Use `-tensStiff` for
> monotonic / pushover analyses.** When combined with the fixed-crack `-interlock`, TS uses
> the *live* `p1` while interlock uses the *frozen* crack normal, so under principal-axis
> rotation TS leaks a small shear onto the frozen plane that interlock then bounds —
> combined TS+interlock is validated for **proportional (non-rotating)** loading. An
> `ε1max`-envelope + secant-unload + frozen-plane TS is the documented cyclic upgrade.

### 4.8 Crack-band regularization — `-autoRegularization $lch_ref` (Phase 3b)

The softening backbones (`-Te/-Ts`, `-Ce/-Cs`) are authored as stress–strain curves, so without
regularization the dissipated energy depends on the element size (mesh-dependent softening). Phase
3b adds the opt-in Bažant–Oh **crack-band** scaling (a faithful clone of `ASDConcrete3D`'s
`-autoRegularization`): the post-peak branch is rescaled so the **specific** fracture energy
`g_reg = G_f0 · (lch_ref / lch)`, hence the **physical** dissipated energy `g_reg · lch` is
**mesh-objective**.

- `lch_ref` is the characteristic length the backbone was authored at (the strain axis is
  "calibrated" to a band of width `lch_ref`).
- `lch` is the element's `getCharacteristicLength()` (already EAS-aware on `ASDShellQ4`), latched
  **once** at the first `setTrialStrain` — unless an explicit `-lch $l` is supplied (which then
  also feeds the interlock `s_θ`). Each Gauss-point/layer copy regularizes with **its own**
  element's `lch`.
- The rescale runs the same iterative scaling + `adjust()` re-enforcement (E-cap, monotone
  plastic strain, non-decreasing damage) as the reference, so steep-softening backbones stay
  physically consistent. `lch == lch_ref` is an exact no-op.

`-autoRegularization` is **default OFF** ⇒ the backbones are used verbatim (bake `G_f/lch` in
yourself, the prior convention) and the material is bit-identical to baseline.

> [!warning] No silent fallback (ADR D5)
> If `-autoRegularization` is on but **no** `lch` can be resolved (no active element *and* no
> `-lch`), the material **refuses to run** (a `FATAL` message + a failed `setTrialStrain`, emitted
> once) rather than silently regularizing with a wrong/default length. Inside an element during
> analysis `ops_TheActiveElement->getCharacteristicLength()` is always available, so this only
> guards genuine misuse (a material exercised outside any element with no `-lch`).

> [!note] Scope (v1)
> This regularizes the **in-plane** softening energy via the element's scalar `lch`. A single
> scalar mis-regularizes inclined (~45°) struts by up to √2 (ADR D5), and through-thickness
> bending-crack energy is not size-objective on the director shell. The **structural** rotated-mesh
> (inclined-crack notched-panel) objectivity gate is staged (the material-point energy-objectivity
> `g_reg·lch = const` is proven across `lch`).

---

## 5. The views and the shell seam

`LadrunoRCConcrete` answers three `getType()` views from one class (classTag 33015):

| `getType()` | order | components | host / role |
|---|---|---|---|
| `ThreeDimensional` (`3D`) | 6 | `{σxx σyy σzz τxy τyz τzx}` | 3D continuum (`LadrunoBrick`/`stdBrick`); the reduce-to-`ASDConcrete3D` baseline |
| `PlaneStress` (`PlaneStress2D`) | 3 | `{σxx σyy τxy}` | 2D continuum / membrane quad; the clean MCFT oracle target — **not** a shell host |
| `PlateFiber` | 5 | `{σxx σyy τxy τyz τzx}` | **the shell path** — rides `LayeredShellFiberSection` `getCopy("PlateFiber")` |

The shell path:

```
ASDShellQ4 [-corotational]  →  section LayeredShell (unmodified)  →  per layer
   getCopy("PlateFiber")  →  LadrunoRCConcrete (PlateFiber view)  →  returnMapPS + transverse-shear block
```

The `PlateFiber` view condenses `σ33 = 0` with a **guarded** nested scalar Newton on `ε33`.
This is the one place where cloning the stock `PlateFiberMaterial` would be a bug: on the
softening branch `dd22 → (1−d)E → 0`, so the bare `1/dd22` step explodes and the stock loop
**silently returns 0 when unconverged**. The RC view guards the step and **propagates a
non-convergence code** instead.

> [!important] State purity under the σ33 Newton
> The condensation re-calls the in-plane update ~25× per converged point. All crack/slip
> state therefore depends **only** on the membrane strains and is written only to the
> output history — idempotent across the inner re-calls — so the Newton cannot corrupt the
> committed crack frame.

> [!note] `lch` resolution (ADR D5 Option A — wired in Phase 3b)
> The material accepts the element's **scalar in-plane** `lch` via
> `getCharacteristicLength()` (which already encodes `ASDShellQ4`'s EAS `/2` correction),
> or an explicit `-lch`. With `-autoRegularization` (§4.8) this `lch` drives Bažant–Oh
> crack-band scaling so the softening energy is mesh-objective; **no silent `lch` default** in a
> softening run (loud failure if unresolved). Through-thickness bending-crack energy is **not**
> size-objective on the director-shell host — that physics belongs to the future solid-shell.

---

## 6. Command syntax

```tcl
nDMaterial LadrunoRCConcrete $tag $E $nu \
    -Ce <e…> -Cs <s…> [-Cd <d…>] \
    -Te <e…> -Ts <s…> [-Td <d…>] \
    [-Kc $Kc] [-betaFloor $f] [-rho $rho] \
    [-beta] [-lublinerReduced] \
    [-secant | -numericalTangent] \
    [-interlock [-agg $ag] [-crackStrain $ec] [-crackSpacing $sθ] [-lch $l] [-betaSrMin $m]] \
    [-cyclic] \
    [-xcrack [-degKappa $k] [-degSlipRef $sr] [-degMin $dm]] \
    [-implex [-implexAlpha $a] [-implexControl $errTol $timeRedLim]] \
    [-tensStiff {vc|cm} [-tensStiffC $c] [-tensStiffAlpha $a]] \
    [-autoRegularization $lch_ref]
```

| Token | Meaning | Default |
|---|---|---|
| `$E $nu` | Young's modulus, Poisson's ratio | — (required) |
| `-Ce/-Cs/-Cd` | compression backbone strain / stress / damage points (≥2; `Ce`,`Cs` equal length) | — (`Ce`,`Cs` required) |
| `-Te/-Ts/-Td` | tension backbone strain / stress / damage points (≥2) | — (`Te`,`Ts` required) |
| `-Kc` | Lee–Fenves biaxial shape (`fb/fc` driver) | 2/3 |
| `-betaFloor` | floor on the MCFT `β` factor | 0.1 |
| `-rho` | smeared web-steel ratio (homogenized for the MCFT composite `ε1`) | 0.0 |
| `-beta` | **enable** MCFT compression softening (§4.1) | OFF |
| `-lublinerReduced` | reduce the Lubliner t–c interaction to avoid double-counting `β` | OFF |
| `-secant` / `-numericalTangent` | use secant / forward-difference tangent instead of the algorithmic one | algorithmic |
| `-interlock` | **enable** aggregate-interlock shear bound (§4.2) | OFF |
| `-agg` | aggregate size `a_g` (consistent length units) | 16.0 |
| `-crackStrain` | normal strain at which the crack freezes | 0.0 |
| `-crackSpacing` | crack spacing `s_θ` for `w = ⟨εn⟩·s_θ` (else `lch` → 1) | 0.0 |
| `-lch` | explicit characteristic length (else from the element) | 0.0 (auto) |
| `-betaSrMin` | residual interlock-shear floor (cap-utilization fraction) | 0.01 |
| `-cyclic` | **enable** cyclic friction-slip (implies `-interlock`) (§4.3) | OFF |
| `-xcrack` | **enable** second crack + cumulative-slip wear (implies `-cyclic`) (§4.4) | OFF |
| `-degKappa` / `-degSlipRef` / `-degMin` | Archard wear: max knockdown / slip reference / residual floor | 0.5 / 0.01 / 0.1 |
| `-implex` | **enable** IMPL-EX integration (robust constant secant for cyclic softening, §4.5) | OFF (fully implicit) |
| `-implexAlpha` | IMPL-EX extrapolation-factor multiplier | 1.0 |
| `-implexControl` | advisory error control `$errTol $timeRedLim` (warns when the implex error exceeds tol) | off |
| `-tensStiff` | **enable** tension stiffening: `vc` (Bentz) or `cm` (Collins–Mitchell) (§4.7) | OFF |
| `-tensStiffC` | `vc`-mode sqrt coefficient `c` (`> 0`; **ignored in `cm` mode** — `cm` hard-codes 500, a warning is emitted) | 500 |
| `-tensStiffAlpha` | `cm`-mode `α1·α2` scale | 1.0 |
| `-autoRegularization` | **enable** crack-band (Bažant–Oh) regularization at reference length `$lch_ref` (§4.8) | OFF |

> [!important] Flag implication chain
> `-xcrack` ⇒ `-cyclic` ⇒ `-interlock` (each law lives inside the previous block). All
> three OFF ⇒ trajectory-faithful to `ASDConcrete3D`. Units must be consistent — the
> `v_ci,max = 0.18√fc'/…` MCFT form is calibrated in **SI (MPa, mm)** (a [[LEDGER_quirks]]
> entry); supply `fc'`/`a_g`/spacing in a matching system.

---

## 7. Worked examples

### 7.1 RC shell wall — OpenSeesPy

```python
import openseespy.opensees as ops

E, nu, Kc = 30000.0, 0.2, 2.0/3.0          # MPa
# tabulated backbones (curved compression + softening tension)
CE = [0.0, 0.0007, 0.0020, 0.0100];  CS = [0.0, 24.0, 30.0, 5.0]
CD = [0.0, 0.0,    0.25,   1.0-5.0/45.0]
TE = [0.0, 0.0001, 0.0010];          TS = [0.0, 3.0, 0.5]
TD = [0.0, 0.0,    1.0-0.5/5.0]

ops.model("basic", "-ndm", 3, "-ndf", 6)
# … nodes …
ops.nDMaterial("LadrunoRCConcrete", 1, E, nu,
               "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
               "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", Kc,
               "-beta", "-lublinerReduced",          # MCFT compression softening
               "-interlock", "-agg", 16.0,           # aggregate-interlock bound
               "-cyclic", "-xcrack")                 # cyclic friction-slip + X-crack wear

# layered shell: 4 layers of a 0.1-thick wall, all the RC material
H, nlay = 0.1, 4
layers = []
for _ in range(nlay):
    layers += [1, H/nlay]
ops.section("LayeredShell", 1, nlay, *layers)
ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)          # corotational optional for large rotation
```

### 7.2 Membrane panel (PlaneStress view) — Tcl

```tcl
nDMaterial LadrunoRCConcrete 1 30000.0 0.2 \
    -Ce 0.0 0.0007 0.0020 0.0100  -Cs 0.0 24.0 30.0 5.0 \
    -Te 0.0 0.0001 0.0010         -Ts 0.0 3.0 0.5 \
    -beta -interlock -agg 16.0
# consumed by a 2D continuum element via its PlaneStress view
```

### 7.3 Plain baseline (all flags off ⇒ ≈ `ASDConcrete3D`)

```tcl
nDMaterial LadrunoRCConcrete 1 30000.0 0.2 \
    -Ce 0.0 0.0007 0.0020 0.0100 -Cs 0.0 24.0 30.0 5.0 -Cd 0.0 0.0 0.25 0.889 \
    -Te 0.0 0.0001 0.0010 -Ts 0.0 3.0 0.5 -Td 0.0 0.0 0.9
# no -beta / -interlock ⇒ trajectory-faithful to ASDConcrete3D (the reduce-to-baseline gate)
```

### 7.4 Cyclic RC wall — quasi-static EXPLICIT recipe (the cyclic solver)

A **cyclic** softening RC wall will **not** converge under implicit Newton — the
plastic-damage consistent tangent goes indefinite on the softening branch, and the
reversals stall (`-implex` helps but does not fully cure a meshed cyclic wall). The
monotonic solvers (`LadrunoArcLength`, `LadrunoIndirectControl`, `LadrunoDynamicRelaxation`,
`robust_drive`) trace **one** equilibrium path through a limit point — a load reversal is
not a single monotonic path, so they don't do cyclic. The right tool is **quasi-static
explicit** (`CentralDifferenceLadruno`): it forms no stiffness tangent, so the
indefinite-tangent stall simply doesn't exist and the reversals integrate through.

```python
import math, openseespy.opensees as ops
# ... build the meshed wall: ASDShellQ4 + section LayeredShell(LadrunoRCConcrete ... -beta
#     -interlock -cyclic -xcrack -rho 2.4e-9); fix the base; lock z + all rotations (planar) ...

# (1) element mass via the material -rho (NOT ops.mass — the dt_cr eigensolve needs element M)
# (2) ASDShellQ4 supplies no per-element dt_cr (criticalTimeStep() == -1): manual wave-speed bound
E, rho, h = 30000.0, 2.4e-9, 500.0           # MPa, tonne/mm^3, element size mm
dt = 0.2 * h / math.sqrt(E / rho)            # ~0.2 * transit time  (frac 0.1-0.2)

# (3) rigid top via per-node sp (NO equalDOF — stability/dt_cr ignore constraints)
#     held axial via a Constant pattern; cyclic drift via a slow cosine Path on a dt-based grid
ops.constraints("Transformation"); ops.numberer("RCM")
ops.system("Diagonal")                       # trivial M^-1
ops.algorithm("Linear")                      # exactly ONE solve/step
ops.integrator("CentralDifferenceLadruno", "-cflAbort", "-lump", "diagonal")
ops.analysis("Transient")
for _ in range(total_steps):                 # seg_time = steps_per_seg*dt MUST be >> structure period
    ops.analyze(1, dt)                       # (quasi-static: keep KE << strain energy)
```

> [!warning] Explicit gotchas (each is a real trap — see [[LEDGER_quirks]])
> - **Element mass, not nodal:** give the material `-rho`; `ops.mass(...)` alone leaves the
>   eigensolve with no element mass matrix.
> - **`ASDShellQ4` reports no `dt_cr`** (`criticalTimeStep()` → −1); using it blindly gives
>   `dt = 0.2·(−1) < 0` → instant blow-up. Use the manual `dt ≈ frac·h/√(E/ρ)`.
> - **No `equalDOF`/rigid ties** — `dt_cr` and the CD stability bound ignore constraints, so a
>   rigid top via `equalDOF` makes the bound a lie and the run detonates. Prescribe the rigid-top
>   drift by putting the **same `sp` value on every top node**.
> - **Mass-proportional damping only** — stiffness-proportional (`betaK`) Rayleigh collapses `dt_cr`.

This is the recipe validated in `tests/test_ladrunoRCConcrete_wall.py` (a 4×3 squat wall
completes the full ±drift schedule; the cyclic interlock degradation is load-bearing at the
panel scale — §9).

---

## 8. State recording

Through the element's `material` response (`eleResponse(ele, "material", gp, "<name>")` or
`recorder … material <gp> <name>`), or the section/shell `stresses` resultants:

| Response name(s) | Returns |
|---|---|
| `stress` / `strain` / `tangent` | the current stress, strain, and tangent of the active view |
| `damage` | `(dt, dc)` tension/compression damage scalars |
| `beta` / `compressionSoftening` | the MCFT softening factor `β` |
| `eps1` / `transverseTensileStrain` | the principal transverse tensile strain `ε1` |
| `betaShear` / `shearRetention` | shear-retention / cap-utilization factor |
| `crackState` / `crackNormal` | crack-1 frozen normal + state |
| `crackShear` | the crack-plane shear `τcr` and cap `v_ci,max` (reads the **trial** history) |
| `xcrackState` | crack-2 state + cumulative slip `slipCum` (X-crack mode) |

> [!note] `crackShear` reads trial state — exclude it from serialization round-trips
> `crackShear`/`tauCr` re-evaluate at the converged displacement (a ~1e-8 Penalty residual
> × stiff `G`), so they drift from the **committed** value by ~5e-4. A `save`/`restore`
> bit-exact test must compare the **committed** fields (crack frame, width, `{τcr,γcr}`,
> `{cracked2,slipCum}`) and drop `crackShear`. Serialization carries a hard-checked
> `RC_SCHEMA_VERSION` field that rejects mismatched vectors.

---

## 9. Verification and validation

Shipped across many phases, each with a Zone-A battery (material-point + shell-element) and
a numpy oracle (`tests/_testbed/rc_shell_ref.py`) the C++ matches step-by-step:

| Phase | What it pins | Tests |
|---|---|---|
| **1** (#155/#192) | `β` on the **strength** axis (`|σc| = β(ε1)·fc'` exact); reduce-to-`ASDConcrete3D` tension+compression byte-match with `β` off; biaxial β-softening ratio == `β(ε1)` | numpy oracle + standalone g++ + `tests/test_ladrunoRCConcrete_material.py` |
| **2a** (#239) | interlock-OFF reduce-to-baseline; `v_ci,max` cap vs closed-form; **off-axis** crack rotation; ON-vs-OFF ablation; standalone g++ FD on the cap + tangent-pinning identity | material-point battery |
| **2b.1** (#245) | reversal re-caps at ∓`v_ci,max` + energy dissipation; unload stiffness == `G`; crack-closure raises the cap; monotonic reduces to the 2a plateau; C++ ≡ numpy over a reversing path | material-point battery |
| **2b.2a** (#246) | **end-to-end in `ASDShellQ4` + `LayeredShell`**: membrane tension cracks+softens; cyclic `Nxy` saturates at `±v_ci,max·h` both signs; interlock load-bearing in the shell; serialization round-trip bit-exact | `tests/test_ladrunoRCConcrete_shell.py` |
| **2b.2b** (#253) | reduce-to-2b.1 (xcrack + `degKappa 0` = bit-identical); crack 2 captures under biaxial tension; cyclic strength degrades monotonically; C++ ≡ numpy with wear | material + shell batteries |
| **2b.2c.1** | `-shearRetention` modes: `const(μ=1)`≡`mcft` to 1e-9; `const` unload == `μ·G`; `dsfm` unload == `G·0.31/denom` (closed-form); cap still `±v_ci,max` every mode; `rots` ≡ interlock-off | material battery + numpy oracle |
| **2b.2c.2** | **rigid-rotation objectivity** — `σ(Q E Qᵀ) == Q σ(E) Qᵀ` over a cracked/cyclic/X-cracked history, frame rotated 30°/90°/127° (the corotational-element route is objective by construction) | `tests/test_ladrunoRCConcrete_objectivity.py` |
| **2b.2c.3** | **crack-closure on the normal** is already correct in the cloned spine (per-step spectral recompose = unilateral closure): compress-past-peak fully recovers; reopened tension follows the damaged envelope | material battery |
| **2b.2c.4a** | **meshed squat wall, quasi-static EXPLICIT** completes the full ±drift cyclic schedule (implicit walls at ~0.6 mm); cyclic interlock degradation load-bearing at panel scale (−28% energy, −11% peak vs monotone) | `tests/test_ladrunoRCConcrete_wall.py` (Zone-B) |
| **4a** (#263) | **IMPL-EX** (`-implex`): off-identical; tracks implicit on a smooth path; error active on rate change / 0 elastic; SPD secant under softening; save/restore continuation | `tests/test_ladrunoRCConcrete_implex.py` |
| **3a** | **tension stiffening** (`-tensStiff`): pre-crack untouched; post-crack `σ_xx == σ_ts(ε1)` closed-form (uniaxial + equibiaxial-both-normals); above-bare; `vc≡cm` defaults; `cm(α)`; PlateFiber-shell floor `Nxx==σ_ts·h`; TS+interlock proportional; FD tangent (g++); schema-v4 round-trip | numpy oracle T1 + standalone g++ + `tests/test_ladrunoRCConcrete_tensstiff.py` |
| **3b** | **crack-band regularization** (`-autoRegularization`): `g_reg·lch` mesh-objective across `lch` (50/25/12.5); `lch==lch_ref` no-op; element-`lch` path; steep-damage plastic-strain monotone (g++); parser guard; schema-v5 round-trip. *Structural rotated-mesh objectivity STAGED.* | numpy oracle R1 + standalone g++ + `tests/test_ladrunoRCConcrete_reg.py` |

> [!note] The load-bearing gate is β-on-the-strength-axis
> Proven three independent ways (numpy oracle, standalone g++ build of `LadrunoRCKernel.h`,
> end-to-end OpenSees). The forbidden abscissa-insertion **misses** the closed-form peak —
> that is the tripwire separating a correct MCFT from a knob that merely slides the lookup.

---

## 10. Limitations and boundaries

> [!caution] Known boundaries
> - **Pinching *shape* is a panel effect.** X-cracking + wear give bidirectional capping
>   and cyclic strength decay at the material point, but the pinched *waist* needs a meshed
>   wall where principal rotation opens/closes the crack within a cycle. Under proportional
>   homogeneous shear the crack sits on the principal plane and the interlock is inert. The
>   meshed wall + quasi-static explicit (§7.4, `tests/test_ladrunoRCConcrete_wall.py`) now
>   demonstrates the panel mechanism (the cyclic interlock degradation is load-bearing:
>   −28% hysteretic energy, −11% peak vs the monotone bound). What remains is a
>   **quantitative match to a named experiment** (Tran–Wallace RW-A20-P10 / a PEER squat
>   wall): specimen geometry + reinforcement layout + digitized measured loops, asserting
>   pinching shape + cumulative hysteretic energy. The solver is no longer the blocker.
> - **No through-thickness `σ33` crush.** The director-shell host cannot carry transverse
>   normal stress — punching / bearing / 3D crush is the deferred `LadrunoSolidShell`
>   (classTag 33020), not this material on `ASDShellQ4`.
> - **Large-rotation cyclic objectivity.** Scalar `dt`/`dc` are frame-indifferent, but the
>   stored crack normal + directional interlock state inherit the dSNPO §14.11 boundary.
>   The corotational **element** (`ASDShellQ4 -corotational`) is the supported objective
>   route; a finite-strain *material* view would be objective only for the isotropic spine.
> - **`lch` regularizes in-plane energy, not band width, and not through-thickness.** A
>   single scalar `lch` mis-regularizes inclined (~45°) struts by up to √2 unless an
>   explicit per-crack `-lch` is supplied; out-of-plane bending-crack energy is not
>   size-objective on the director-shell host (ADR D5).
> - **Smeared web steel only, inside the kernel.** Discrete boundary rebar where buckling
>   matters is a separate `PlateRebar(LadrunoRebarBuckling)` layer at the true depth, with
>   the overlapping concrete reduced (no ρ-weighted double-count).
> - **Tension stiffening** is available two ways: bake a Collins–Mitchell / Bentz
>   averaged-tension curve into the tabulated `-Te/-Ts` backbone, **or** use the dedicated
>   `-tensStiff {vc|cm}` knob (§4.7, Phase 3a) which floors the live principal tensile stress
>   to `σ_ts(ε1)` post-crack. The `-tensStiff` floor is **monotonic-scope** (re-inflates on
>   unload; use for pushover), and combined with `-interlock` is validated for proportional
>   loading. Default off ⇒ baseline.

> [!note] What's next (the deferred phases)
> The cyclic constitutive physics is complete and the `-tensStiff` tension-stiffening knob
> (Phase 3a) is shipped. The remaining work — directional `lch` (Phase 3b), the finite-strain view
> `LadrunoRCFiniteStrain` (Phase 4b), the through-thickness `LadrunoSolidShell` punching host
> (Phase 5), and the quantitative Tran–Wallace experiment calibration — is planned in the
> developer handout **[[rc_shell_phase3plus_handout]]** (goal, build, reuse, acceptance gates,
> gotchas, and effort per item).

---

## 11. References

**Theory**
- Vecchio & Collins (1986), *The Modified Compression-Field Theory for Reinforced Concrete
  Elements Subjected to Shear* — the `β` compression-softening law and `v_ci,max` crack-shear
  limit.
- Vecchio (2000), *Disturbed Stress Field Model* (DSFM) — the rotating + slip variant.
- Lubliner, Oliver, Oller & Oñate (1989); Lee & Fenves (1998) — the biaxial plastic-damage
  envelope inherited from `ASDConcrete3D`.
- Bažant & Oh (1983) — crack-band regularization.
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* — §9.4
  (plane-stress condensation hazard), §14.11 (large-rotation objectivity boundary).

**Within this repo**
- [[19_ladruno_rc_shell_adr]] — the source ADR: decisions D1–D6, the seam analysis, the
  full risk register, the phased plan (Phase 1 / 2a / 2b.1 / 2b.2a / 2b.2b).
- [[LadrunoConcrete3D_user_guide]] — the `ASDConcrete3D` spine this material clones.
- [[Ladruno_materials_guide]] — the materials catalog.
- [[LadrunoRebarBuckling_guide]] — the discrete boundary-rebar layer partner.
- [[LogStrain_guide]] — the finite-strain lift (isotropic-spine-only for RC).
- Source: `SRC/material/nD/LadrunoRCKernel.h` (header-only kernel),
  `LadrunoRCConcrete.{h,cpp}` (the multi-dim views, classTag **33015**). Oracle:
  `tests/_testbed/rc_shell_ref.py`. Tests: `tests/test_ladrunoRCConcrete_material.py`,
  `tests/test_ladrunoRCConcrete_shell.py`.
