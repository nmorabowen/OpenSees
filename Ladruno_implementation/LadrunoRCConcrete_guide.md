---
title: "LadrunoRCConcrete — nonlinear RC layer-shell constitutive material"
project: Ladruno
type: reference / user guide
status: shipped (Phase 1 + Phase 2a + 2b.1 + 2b.2a + 2b.2b; MERGED PRs #155 / #192 / #239 / #245 / #246 / #253)
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

### 4.5 IMPL-EX robustness — `-implex`

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

> [!note] `lch` resolution (ADR D5 Option A)
> The material accepts the element's **scalar in-plane** `lch` via
> `getCharacteristicLength()` (which already encodes `ASDShellQ4`'s EAS `/2` correction),
> or an explicit `-lch`. Through-thickness bending-crack energy is **not** size-objective on
> the director-shell host — that physics belongs to the future solid-shell. There is no
> silent `lch` default in a softening run.

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
    [-implex [-implexAlpha $a] [-implexControl $errTol $timeRedLim]]
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

Shipped across five phases, each with a Zone-A battery (material-point + shell-element) and
a numpy oracle (`tests/_testbed/rc_shell_ref.py`) the C++ matches step-by-step:

| Phase | What it pins | Tests |
|---|---|---|
| **1** (#155/#192) | `β` on the **strength** axis (`|σc| = β(ε1)·fc'` exact); reduce-to-`ASDConcrete3D` tension+compression byte-match with `β` off; biaxial β-softening ratio == `β(ε1)` | numpy oracle + standalone g++ + `tests/test_ladrunoRCConcrete_material.py` |
| **2a** (#239) | interlock-OFF reduce-to-baseline; `v_ci,max` cap vs closed-form; **off-axis** crack rotation (normal vs analytic principal dir); ON-vs-OFF ablation; standalone g++ FD on the cap + tangent-pinning identity | material-point battery |
| **2b.1** (#245) | reversal re-caps at ∓`v_ci,max` + energy dissipation; unload stiffness == `G`; crack-closure raises the cap; monotonic reduces to the 2a plateau; C++ ≡ numpy over a full reversing path | material-point battery |
| **2b.2a** (#246) | **end-to-end in `ASDShellQ4` + `LayeredShell`**: membrane tension cracks+softens; cyclic membrane shear saturates `Nxy` at `±v_ci,max·h` both signs; interlock load-bearing in the shell; serialization round-trip bit-exact | `tests/test_ladrunoRCConcrete_shell.py` |
| **2b.2b** (#253) | reduce-to-2b.1 (xcrack + `degKappa 0` = bit-identical); crack 2 captures under biaxial tension; cyclic strength degrades monotonically over cycles; C++ ≡ numpy with wear; shell-level cyclic decay | material + shell batteries |

> [!note] The load-bearing gate is β-on-the-strength-axis
> Proven three independent ways (numpy oracle, standalone g++ build of `LadrunoRCKernel.h`,
> end-to-end OpenSees). The forbidden abscissa-insertion **misses** the closed-form peak —
> that is the tripwire separating a correct MCFT from a knob that merely slides the lookup.

---

## 10. Limitations and boundaries

> [!caution] Known boundaries
> - **Pinching *shape* is a panel effect.** X-cracking + wear give bidirectional capping
>   and cyclic strength decay at the material point, but the pinched *waist* needs a meshed
>   wall where principal rotation opens/closes the crack within a cycle (deferred Phase
>   2b.2c). Under proportional homogeneous shear the crack sits on the principal plane and
>   interlock is inert.
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
> - **Tension stiffening** is supplied through the tabulated `-Te/-Ts` tension backbone
>   (bake a Collins–Mitchell / Bentz averaged-tension curve into the points) rather than a
>   separate flag; the default bare fracture-energy backbone reduces to baseline. A
>   dedicated `-tensStiff` knob is reserved in the ADR but not yet wired.

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
