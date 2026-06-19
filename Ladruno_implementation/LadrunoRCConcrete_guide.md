---
title: "Ladruno nonlinear layer-shell RC — implementation guide (theory · architecture · usage · quirks)"
project: Ladruno
type: reference / user + developer guide
status: shipped — Phase 1 / 2a / 2b.1 / 2b.2a / 2b.2b / 2b.2c.1-4a + IMPL-EX (4a) + Phase 3 (3a tension stiffening + 3b crack-band lch regularization + 3b-struct mesh-objectivity gate) + Phase 4b (finite-strain view LadrunoRCFiniteStrain, classTag 33018); cyclic material physics + regularization + finite-strain view COMPLETE. Remaining: quantitative Tran–Wallace experiment match, inclined/notched-panel structural objectivity (staged Zone-B), Phase 5 (LadrunoSolidShell through-thickness host). MERGED PRs #155 #192 #239 #245 #246 #253 #263 #266 #273 #277 #281 #282
classTag: 33015 (small-strain views) · 33018 (LadrunoRCFiniteStrain, finite-strain view)
material: "nDMaterial LadrunoRCConcrete — 3D / PlaneStress / PlateFiber views · nDMaterial LadrunoRCFiniteStrain — finite-strain (Hencky) view"
related:
  - "[[19_ladruno_rc_shell_adr]]"
  - "[[rc_shell_phase3plus_handout]]"
  - "[[LadrunoConcrete3D_user_guide]]"
  - "[[Ladruno_materials_guide]]"
  - "[[LadrunoRebarBuckling_guide]]"
  - "[[LogStrain_guide]]"
  - "[[LEDGER_quirks]]"
tags:
  - material
  - rc-shell
  - mcft
  - plane-stress
  - compression-softening
  - aggregate-interlock
  - crack-band
  - finite-strain
  - architecture
  - reference
---

# Ladruno nonlinear layer-shell RC — implementation guide

> [!abstract] What this stack is, in one paragraph
> The **nonlinear layer-shell RC** capability lets OpenSees model cracking, softening,
> cyclically-degrading reinforced-concrete **walls, slabs, shells and membrane panels** with
> **zero element or section edits** — the missing physics is delivered as a *material*. At its
> core is `nDMaterial LadrunoRCConcrete` (classTag **33015**): a plastic-damage concrete that
> keeps `ASDConcrete3D`'s verified spine (spectral effective-stress split, dual tension/compression
> damage, Lubliner–Lee–Fenves biaxial envelope, Bažant–Oh crack-band regularization) and adds the
> **physics `ASDConcrete3D` lacks** for RC — MCFT **compression (transverse-tension) softening**,
> a **degrading aggregate-interlock** shear-transfer law (monotonic bound → cyclic friction-slip →
> X-cracking + wear → retention curves), **tension stiffening**, **IMPL-EX** robustness, and a
> consistent tangent that honours the new couplings. It ships as **one multi-dim class** offering
> `ThreeDimensional` / `PlaneStress` / `PlateFiber` views, so the same material that runs on a 3D
> brick drops into a `LayeredShell` under `ASDShellQ4`. A **finite-strain companion**
> `LadrunoRCFiniteStrain` (classTag **33018**) lifts the same kernel to large strain for
> `LadrunoBrick -geom finite`. **Every added layer is a flag, default OFF** — with all flags off the
> material is trajectory-faithful to `ASDConcrete3D`.

This document is the **complete implementation reference** for the stack, in four parts:

- **Part I — Theory** — the constitutive gap and every physics layer, with equations.
- **Part II — Architecture** — how it is built: the header-only kernel, the "one kernel, many
  views" class map, the element/section seam, the finite-strain companion, state &
  serialization, registration, and the test harness.
- **Part III — Usage** — command grammar, worked examples, the cyclic explicit solver recipe,
  recording, units.
- **Part IV — Quirks** — the consolidated gotchas (the canonical copies live in [[LEDGER_quirks]]).

For the full decision log (D1–D6), the seam analysis, and the phased plan, see the source ADR
[[19_ladruno_rc_shell_adr]]. For the spine it clones, see [[LadrunoConcrete3D_user_guide]]. For
the forward plan (Phase 5 + remaining gates), see [[rc_shell_phase3plus_handout]].

---

## Contents

**Part I — Theory**
- [[#1 The constitutive gap|1. The constitutive gap]]
- [[#2 The inherited plastic-damage spine|2. The inherited plastic-damage spine]]
- [[#3 MCFT compression softening (`-beta`)|3. MCFT compression softening]]
- [[#4 Aggregate-interlock shear transfer|4. Aggregate-interlock shear transfer]]
- [[#5 Tension stiffening (`-tensStiff`)|5. Tension stiffening]]
- [[#6 Crack-band regularization (`-autoRegularization`)|6. Crack-band regularization]]
- [[#7 IMPL-EX robustness (`-implex`)|7. IMPL-EX robustness]]
- [[#8 The finite-strain (Hencky) lift|8. The finite-strain lift]]
- [[#9 Objectivity — the §14.11 split|9. Objectivity — the §14.11 split]]

**Part II — Architecture**
- [[#10 Material, not element — the keystone seam|10. Material, not element]]
- [[#11 One kernel, many views — the class map|11. One kernel, many views]]
- [[#12 The shell path and the σ33 condensation|12. The shell path]]
- [[#13 The finite-strain companion (`LadrunoRCFiniteStrain`)|13. The finite-strain companion]]
- [[#14 State, serialization and the IMPL-EX runtime|14. State & serialization]]
- [[#15 Registration and dispatch|15. Registration & dispatch]]
- [[#16 The verification harness|16. The verification harness]]

**Part III — Usage**
- [[#17 Command grammar|17. Command grammar]]
- [[#18 Worked examples|18. Worked examples]]
- [[#19 The cyclic quasi-static explicit recipe|19. The cyclic explicit recipe]]
- [[#20 State recording|20. State recording]]
- [[#21 Units|21. Units]]

**Part IV — Quirks**
- [[#22 Consolidated quirks|22. Consolidated quirks]]

**Appendices**
- [[#A. Phase and PR history + V&V matrix|A. Phase/PR history + V&V matrix]]
- [[#B. References|B. References]]

---
---

# Part I — Theory

## 1. The constitutive gap

The deficiency that makes OpenSees **over-predict squat-wall capacity** and **miss cyclic
pinching** is **constitutive, not elemental**. For an in-plane-loaded wall the structural shear
`V` lives in the **membrane** block (`γxy`/`τxy`), so it is governed by the plane-stress
*constitutive* law — not by transverse shear and not by element technology. Two physics are
missing from `ASDConcrete3D`:

1. **No compression softening from transverse cracking.** `ASDConcrete3D`'s compressive measure
   uses only the negative effective-stress principals; its biaxial Lubliner envelope is **not** the
   cracked-strut softening that MCFT captures. A diagonal compression strut that is simultaneously
   cracked transversely should lose strength — `ASDConcrete3D` does not model that, so it
   over-predicts strut capacity.
2. **No crack-plane shear-transfer / retention term.** A pure rotating-crack damage model is
   *coaxial* — it has no independent crack-plane shear stiffness, so it cannot carry (or *bound*)
   cracked shear, and it structurally cannot produce cyclic **pinching** (which needs a crack that
   opens/closes and slides during the cycle).

These gaps are closed by adding **four RC physics layers** (Part I §3–§6) on top of the cloned
spine (§2), each a flag that is **OFF by default** so the baseline is byte-faithful.

> [!tip] The headline use case
> **Squat / shear-controlled RC walls.** The MCFT `β` softening + the aggregate-interlock cap
> directly correct the over-predicted diagonal-strut capacity; the cyclic friction-slip + X-crack
> wear give the pinching and cyclic strength decay a rotating-coaxial model cannot. Slabs/shells in
> bending, 2D membrane panels, and 3D continuum RC are the other targets (Part III §18).

## 2. The inherited plastic-damage spine

The base law is `ASDConcrete3D`'s plastic-damage spine, cloned **verbatim** into the header-only
`LadrunoRCKernel.h` (see [[LadrunoConcrete3D_user_guide]] for the full theory):

- **Spectral effective-stress split** `σ̃ = σ̃⁺ + σ̃⁻` into tensile / compressive cones via the
  eigenprojectors `PT`/`PC` of `σ̃`.
- **Dual damage** `dt` (tension) and `dc` (compression), each driven by its own equivalent-strain
  measure on a **tabulated backbone** (`-Te/-Ts/-Td`, `-Ce/-Cs/-Cd`).
- **Lubliner–Lee–Fenves biaxial envelope** (`fb = 1.16 fc`, `Kc`, `γ = 3(1−Kc)/(2Kc−1)`) raising the
  effective-compressive measure in the tension–compression quadrant.
- **Crack-band (Bažant–Oh) `lch` regularization** so dissipated energy is mesh-objective (§6).
- **Crack-closure** spectral reassembly (per-step recompose = unilateral closure) and an optional
  IMPL-EX path (§7).

The assembled nominal stress is `σ = (1 − dt̄)·ST + (1 − dc̄)·SC` (with the MCFT `β` factor inserted
on the compressive cone, §3). Backbones are point lists (strain `e`, stress `s`, optional damage
`d`), exactly the `ASDConcrete3D` `HardeningLaw` format:

```
-Ce e0 e1 e2 …   -Cs s0 s1 s2 …   [-Cd d0 d1 d2 …]     # compression (>=2 points)
-Te e0 e1 e2 …   -Ts s0 s1 s2 …   [-Td d0 d1 d2 …]     # tension (>=2 points)
```

`-Ce`/`-Cs` and `-Te`/`-Ts` are **required and equal length** (≥2); `-Cd`/`-Td` are optional (pad to
zero ⇒ pure plastic, no damage). The peak of `-Ts` is `ft`; the peak of `-Cs` is `fc'`.

> [!important] Two spine-cloning gotchas (Quirk §22.1)
> (1) The equivalent-strain measure needs the `/E` division — a `β`-ratio test cannot catch its
> omission; only an **absolute-stress** test can. (2) The effective-stress backbone `q` is
> **E-consistent by construction** (`buildBackbone` mirrors the `ASDConcrete3D` `HardeningLaw`
> constructor + `adjust()`), **not** `q = y/(1−d)` on raw points.

## 3. MCFT compression softening (`-beta`)

The Modified Compression Field Theory (Vecchio & Collins 1986) softens the compressive **strength**
by the transverse cracking strain:

$$\beta = \frac{1}{0.8 + 170\,\varepsilon_1} \le 1, \qquad \frac{\partial\beta}{\partial\varepsilon_1} = -170\,\beta^2$$

where `ε1` is the average principal **tensile** (cracking) strain transverse to the strut, and the
realised compressive peak becomes `|σc| = β·fc'`. The assembled stress with `β` is
`σ = (1 − dt̄)·ST + β·(1 − dc̄)·SC` (Derivation B), and the consistent tangent carries the `−170 β²`
cross-term.

> [!important] β scales the STRENGTH axis — the blocking correctness gate (ADR D4)
> `β` must scale the **stress/strength axis** (the effective-compressive backbone value `q` fed into
> `dc`), **never** the *strain abscissa* fed to the equivalent-strain measure. Scaling the abscissa
> merely slides the lookup along a fixed non-linear curve, giving a realised peak `hc(β·xc) ≠ β·fc'`.
> Proven **three ways** (numpy oracle, standalone g++, end-to-end OpenSees): hold confined
> compression, sweep `ε1`, assert realised `|σc| = β(ε1)·fc'` exactly. **SFI-MVLEM = FSAM ≠ MCFT** —
> do not mistake the FSAM material for an MCFT oracle.

> [!note] Avoiding double-count with the Lubliner envelope — `-lublinerReduced`
> The Lubliner envelope already raises the effective-stress measure in the tension–compression
> quadrant where struts live. Stacking `β` on top double-penalises that quadrant. `-lublinerReduced`
> reduces the t–c interaction so transverse tension is counted **once**; `β` and `Kc` are meant to
> be calibrated **jointly** against a tension–compression panel battery (Kupfer + Vecchio–Collins).

## 4. Aggregate-interlock shear transfer

The crack-plane shear-transfer capability is built in four cumulative sub-laws; each is a flag, and
each lives **inside** the previous block (`-xcrack` ⇒ `-cyclic` ⇒ `-interlock`).

### 4.1 Monotonic shear bound — `-interlock` (Phase 2a)

The material freezes the in-plane crack **normal** `(c,s)` at the first crossing of `crackStrain`,
and thereafter **clips** the smeared (damage-reduced) crack-plane shear `τ_sm = m_σ·σ_ip` to the MCFT
crack-shear limit:

$$v_{ci,\max} = \frac{0.18\sqrt{f_c'}}{0.31 + 24w/(a_g+16)}, \qquad w = \langle\varepsilon_n\rangle\, s_\theta$$

`w` = crack width (Macauley of the strain normal to the **frozen** crack × spacing `s_θ`), `a_g` =
aggregate size (`-agg`, default 16). This is a **bound on the existing shear**, *not* a replacement
with bare-elastic `G·γ` (the bare-elastic-then-cap form was **rejected** in a 59-agent review for a
cracking-instant discontinuity + tangent inconsistency). Below the cap the stress is unchanged
(continuous). The consistent tangent below the cap is the baseline; at the cap it is the exact
rank-1 removal `Dtan -= (1−betaSrMin)·m_ε⊗(m_σ·Dtan)`, exact because `m_σ·m_ε = (c²+s²)² = 1`.

> [!tip] Interlock engages only under NON-proportional loading (Quirk §22.3)
> Under *proportional* shear the crack forms on the principal plane, so `g_nt = 0` and interlock is
> inert. It activates when the loading path rotates the principal directions away from the frozen
> crack — i.e. in real panels / non-radial paths.

### 4.2 Cyclic friction-slip — `-cyclic` (Phase 2b.1)

`-cyclic` promotes the crack-plane shear to an **independent incremental friction-slip state**:

$$\tau_{cr} = \mathrm{clamp}\bigl(\tau_{cr}^{c} + G_{\text{slip}}\,(\gamma_{nt} - \gamma_{nt}^{c}),\; \pm v_{ci,\max}(w)\bigr), \quad G = \frac{E}{2(1+\nu)}$$

with the slip origin frozen at crack capture (seeded to `clip(τ_sm)` for stress continuity). The
crack width `w` is now **reversible** (current opening, not the monotone max) so crack closure raises
the cap (capacity recovery). State = two committed scalars `{τcr, γcr}`. It diverges deliberately
from FSAM in three respects: the cap is MCFT `v_ci,max(w)` not Coulomb `μ·f_n`; the default slip
stiffness is the full elastic `G` not `0.4·Ec`; an open crack still transfers shear up to the
width-degraded cap rather than collapsing to zero.

> [!note] Pinching is a PANEL effect (Quirk §22.4)
> At a material point with a fixed crack and constant normal strain the loop is *fat*; a pinched
> waist needs the crack to open/close during the cycle (principal rotation), which only happens in a
> meshed wall. The material gives the **ingredients** (reversal re-capping, friction unload at
> stiffness `G`, capacity recovery on closure); the waist emerges at the structural scale.

### 4.3 Second crack (X-cracking) + cyclic wear — `-xcrack` (Phase 2b.2b)

`-xcrack` adds an orthogonal crack 2 (normal `(−s,c)`) that freezes when the strain normal to *that*
plane reaches `crackStrain`. The interlock cap then uses the **governing (most-open) crack's
opening** `v_ci,max(⟨max(εn1, εn2)⟩·s_θ)` — so the **reverse** shear direction is also capped (a
single crack only capped one direction). There is **one** friction state and **one** governing cap —
no double-count (the design-workflow's two-state SUM provably double-counted ≈2·v_ci,max).

Cyclic wear is an irreversible Archard-style knockdown driven by the **cumulative crack sliding
distance** `slipCum` (not peak slip, which saturates in cycle 1 at constant amplitude):

$$k_n = \max\bigl(\text{degMin},\; 1 - \text{degKappa}\cdot\mathrm{clamp}(\text{slipCum}/\text{degSlipRef}, 0, 1)\bigr), \qquad v_{ci,\max} \leftarrow k_n\, v_{ci,\max}$$

Driven by the **committed** `slipCum` ⇒ tangent-neutral. Verified caps `2.90 → 2.50 → 2.10` over
three cycles. Knobs: `-degKappa` (0.5), `-degSlipRef` (0.01), `-degMin` (0.1).

### 4.4 Crack-shear retention curve — `-shearRetention` (Phase 2b.2c.1)

The cyclic crack-plane **slip stiffness** `G_slip` in the friction predictor is selectable. The
`v_ci,max(w)` **cap is unchanged** in every mode; all modes reduce to `mcft`.

| `-shearRetention <mode>` | `G_slip` | meaning |
|---|---|---|
| `mcft` (default) | `G = E/2(1+ν)` | full elastic slip stiffness — the shipped DSFM-with-slip law |
| `const` | `μ·G` | classic constant retention (Rots/Červenka); `μ = -shearRetFactor` ∈ (0,1], default 0.4 |
| `dsfm` | `G·(0.31/denom)` | DSFM width-degraded (`denom = 0.31 + 24w/(a_g+16)`); `→ G` at `w=0` |
| `rots` | — | rotating-coaxial: the fixed-crack shear is skipped ⇒ identical to `-interlock` OFF |

`const`/`dsfm` imply `-interlock -cyclic`; `rots` implies `-interlock`. Omitting it (or `mcft`) is
2b.2b byte-identical. Because `const`/`dsfm` parametrize the *slip* stiffness, they bite only on the
`-cyclic` path; on the 2a monotone-clip path they are inert (same `v_ci,max` plateau).

## 5. Tension stiffening (`-tensStiff`)

Between cracks, bonded reinforcement carries part of the tension, so the **average** concrete tensile
stress stays **above** the bare fracture-energy softening curve. Phase 3a adds this as an opt-in
**stress floor** on the live in-plane principal tensile axis `p1`:

$$\sigma_{ts}(\varepsilon_1)=\begin{cases}\dfrac{f_t}{1+\sqrt{c\,\varepsilon_1}} & \texttt{vc}\ \text{(MCFT / Bentz)}\\[2ex]\dfrac{\alpha\,f_t}{1+\sqrt{500\,\varepsilon_1}} & \texttt{cm}\ \text{(Collins–Mitchell)}\end{cases}$$

with `f_t` the tension-backbone peak and `ε1` the **composite** (reinforced) membrane principal
tensile strain — the *same* `ε1` the MCFT `β` uses (perfect-bond compatibility), so it needs no
separate steel layer's strain. The floor injects `Δ = σ_ts(ε1) − n^T σ n` along `p1` (rank-1, only
when `Δ>0`); it is active **only post-crack** (`ε1 ≥ ε_cr`). Equibiaxial (degenerate `p1`) floors
**both** in-plane normals. `vc` with `c=500` ≡ `cm` with `α=1`.

> [!warning] Monotonic-scope (v1) (Quirk §22.13)
> `σ_ts` is a function of the **live** `ε1` (no `ε1max` memory). Because `σ_ts` *decreases* with `ε1`,
> on **unloading** the floor **re-inflates** — correct on the monotone loading branch (slab /
> distributed-reinforcement pushover) but **not** a hysteretic cyclic-tension model. **Use
> `-tensStiff` for monotonic / pushover analyses.** Combined with `-interlock`, TS uses the *live*
> `p1` while interlock uses the *frozen* crack normal — validated for **proportional (non-rotating)**
> loading; an `ε1max`-envelope + secant-unload + frozen-plane TS is the documented cyclic upgrade.

## 6. Crack-band regularization (`-autoRegularization`)

The softening backbones are authored as stress–strain curves, so without regularization the
dissipated energy depends on element size (mesh-dependent softening). Phase 3b adds the Bažant–Oh
**crack-band** scaling (a faithful clone of `ASDConcrete3D`'s `-autoRegularization`): the post-peak
branch is rescaled so the **specific** fracture energy `g_reg = G_f0 · (lch_ref / lch)`, hence the
**physical** dissipated energy `g_reg · lch` is **mesh-objective**.

- `lch_ref` is the length the backbone was authored at; `lch` is the element's
  `getCharacteristicLength()` (EAS-aware on `ASDShellQ4`), latched **once** at the first
  `setTrialStrain` — unless an explicit `-lch $l` is supplied (which then also feeds the interlock
  `s_θ`). Each Gauss-point/layer copy regularizes with **its own** `lch`.
- The rescale runs the same iterative scaling + `adjust()` re-enforcement (E-cap, monotone plastic
  strain, non-decreasing damage) as the reference, so steep backbones stay physically consistent.
  `lch == lch_ref` is an exact no-op.

> [!warning] No silent fallback (ADR D5) (Quirk §22.10)
> If `-autoRegularization` is on but no `lch` resolves (no active element *and* no `-lch`), the
> material **refuses to run** (a `FATAL` + a failed `setTrialStrain`, emitted once) rather than
> silently regularizing with a wrong length. The failure does **not** latch, so an algorithm
> step-retry cannot silently proceed un-regularized.

The **structural** mesh-objectivity is proven by the 3b-struct Bažant-bar gate (a row of `ASDShellQ4`
at three meshes dissipates a mesh-objective total energy; the un-regularized control halves and
snap-back-diverges). The residual: a single scalar `lch` mis-regularizes inclined (~45°) struts by up
to √2, and through-thickness bending-crack energy is not size-objective on the director shell — the
inclined/notched-panel gate and the solid-shell carry those.

## 7. IMPL-EX robustness (`-implex`)

Softening plastic-damage produces an **indefinite consistent tangent** on the softening branch, so a
vanilla Newton solve of a cyclic RC wall stalls. IMPL-EX (Oliver et al. 2008) is the escape: each
step the damage thresholds `xt, xc` and the MCFT `β` are **frozen** by an explicit extrapolation of
the committed history,

$$x_{\text{ext}} = x_n + t_f\,(x_n - x_{n-1}), \qquad t_f = \frac{\Delta t_{n+1}}{\Delta t_n}\,\alpha,$$

so the returned tangent is the **secant** `W_B·C0` with no softening-rate / no `β` cross-term — it
removes the dominant indefinite contribution and lets Newton converge. The **true implicit**
thresholds are recomputed at commit to advance the state and measure the IMPL-EX error
(`implexError` response).

> [!important] What is and isn't frozen (honest scope)
> The **damage** (the softening source) is frozen — the robustness that matters. The spectral
> **projectors** PT/PC are **not** frozen (recomputed from the live stress each call), consistent with
> the material's fixed-projector-secant philosophy. IMPL-EX is **exact on proportional paths**
> (error ~1e-12) and lags `O(Δt)` at rate changes.

> [!warning] Static analysis: the time factor is guarded (Quirk §22.7)
> In a static analysis the "time" is the load factor, whose increment is erratic (and resets at
> `loadConst`), so the raw `Δt_{n+1}/Δt_n` would detonate the extrapolation. The material guards `t_f`
> (falls back to `α`, clamps to `2α`) — supply roughly uniform load steps.

`-implex` is **default OFF** ⇒ fully-implicit, bit-identical. For a *meshed cyclic wall* `-implex`
helps but does not fully cure the reversals — the robust cyclic tool is quasi-static explicit
(Part III §19).

## 8. The finite-strain (Hencky) lift

For large-strain runs (`LadrunoBrick -geom finite`, the future solid-shell) the same physics is
available through `LadrunoRCFiniteStrain` (classTag **33018**, architecture in §13), via the Hencky
(logarithmic-strain) seam (de Souza Neto, Perić & Owen 2008, Box 14.3):

$$\mathbf B=\mathbf F\mathbf F^{\mathsf T},\quad \boldsymbol\varepsilon^e=\tfrac12\ln\mathbf B,\quad (\boldsymbol\tau,\mathsf D)=\texttt{returnMap3D}(\boldsymbol\varepsilon^e),\quad \boldsymbol\sigma=\boldsymbol\tau/J,\quad \mathsf c=\tfrac{1}{2J}[\mathsf D:\mathsf L:\mathsf B]$$

The **same** `LadrunoRCKernel` return map runs verbatim, so every RC flag + IMPL-EX behave
identically — only the strain measure differs. Because the RC spine carries **no tensorial plastic
strain**, the elastic left Cauchy–Green is simply the **total** `B = F Fᵀ` recomputed each step (no
`bᵉ` to track) — the multiplicative split is trivial. This is what makes the lift clean (and is why a
native subclass beats the generic `LogStrain` wrapper — §13).

## 9. Objectivity — the §14.11 split

The constitutive frame-indifference identity `σ(QEQᵀ) = Q σ(E) Qᵀ` splits in two:

- **Isotropic spine — objective.** Scalar `dt`/`dc` and the spectral split are frame-indifferent.
  Verified at the small-strain material point (rotated 30°/90°/127° over a cracked/cyclic/X-cracked
  history) **and** for the finite view (`σ(QF) = Q σ(F) Qᵀ` at two rotations).
- **Directional state — bounded by dSNPO §14.11.** The stored crack **normal** and the
  interlock/friction direction are *directional internal variables* in a fixed frame. The supported
  objective large-rotation route is the corotational **element** (`ASDShellQ4 -corotational`), which
  feeds the material the *de-rotated* strain `QᵀEQ` — objective by construction. The finite-strain
  *material* view (`setTrialF`) does **not** co-rotate the stored crack frame, so a *cracked* point is
  not frame-indifferent under large rotation (a documented xfail; the RC analog of
  `LadrunoJ2Finite`'s backstress co-rotation would flip it — deferred).

**Rule of thumb:** large-rotation **cyclic** walls → corotational element + small-strain
`LadrunoRCConcrete`. Monotonic / isotropic-damage-dominated large-strain solids → `LadrunoRCFiniteStrain`.

---
---

# Part II — Architecture

## 10. Material, not element — the keystone seam

`ASDShellQ4` (AGQ6-I membrane + 4-DOF EAS + MITC4 transverse shear + Hughes–Brezzi drilling + optional
corotational) already equals the LS-DYNA `ELFORM=16` fully-integrated assumed-strain shell;
`LayeredShellFiberSection` already does standard director-stack integration. Rebuilding either re-treads
solved element technology for **zero new physics**. So the keystone decision (ADR D1/D2): **keep the
element and section at the frontier, deliver the missing physics as a material** that drops into the
5-component `PlateFiber` layer view the section already requests.

```
element ASDShellQ4 (unmodified)                       # 8-comp generalized strain E = [exx eyy gxy | kxx kyy kxy | gxz gyz]
   └─ section LayeredShell (unmodified)               # per layer: plane state e(z) = e_m − z·κ
        └─ getCopy("PlateFiber")  →  LadrunoRCConcrete (PlateFiber view, order 5)
```

The element computes `E`, the section forms each layer's `PlateFiber` strain, and `LadrunoRCConcrete`
answers on that seam — **untouched element, untouched section** (ADR adversarial finding: a
`LayeredShellFiberSection` subclass that does not override `getCopy()` is *sliced* to vanilla on the
per-GP clone at `ASDShellQ4.cpp:740`, so v1 needs zero section edit).

## 11. One kernel, many views — the class map

Physics lives in the **header-only, OpenSees-free** `LadrunoRCKernel.h` (namespace
`ladruno_rc_kernel`), numpy-oracle-testable and g++-buildable **before any OpenSees link** — cloned
from the proven `LadrunoJ2Kernel.h`. The OpenSees-facing classes are thin wrappers over it.

| Class | classTag | base | role |
|---|---|---|---|
| `LadrunoRCConcrete` | **33015** | `NDMaterial` | the small-strain material — **one** class, three `getType()` views (`3D` / `PlaneStress` / `PlateFiber`) selected by `dim`/`vmap` |
| `LadrunoRCFiniteStrain` | **33018** | `FiniteStrainNDMaterial` | the finite-strain (Hencky) view — `ThreeDimensional`-only, driven by `setTrialF(F)` (§13) |

Kernel contents (all `inline`, plain doubles):

- `struct Params` — material constants + the two `HardeningLaw` backbones (`ht`, `hc`) + every flag.
- `struct RCHist` — the committed/trial history (the **entire** path-dependent state): effective
  stress, total strain, damage thresholds `xt/xc`, the frozen crack frame `{cracked, crackC, crackS,
  wmax}`, the cyclic friction `{tauCr, gammaCr}`, the X-crack/wear `{cracked2, slipCum}`, and the
  IMPL-EX `{xt_old, xc_old, eps1_old}`.
- `returnMap3D(P, eps6, in, …) → (sig6, Dtan6, out)` — the **pure** constitutive update (no global
  state; deterministic given the committed history `in`). Uses an *incremental* elastic predictor
  `SEFF = in.stress_eff + C0:(eps − in.strain)`, so the committed strain in `RCHist` is the de-facto
  committed reference.
- `buildBackbone` / `setupParams` — the E-consistent `HardeningLaw` constructor + `adjust()`.
- `regularize` / `fractureEnergy` / `adjustBackbone` — the Bažant–Oh crack-band scaling (§6).
- `elasticTangent`, `lublinerCriterion`, `computePjj`, the spectral helpers.

> [!note] Why one class, many views (the `LadrunoJ2` / `ASDConcrete3D` pattern)
> The in-plane return map is shared **byte-for-byte** between the `PlaneStress` and `PlateFiber`
> views; the `3D` view is the reduce-to-`ASDConcrete3D` baseline. A single `LadrunoRCConcrete` with a
> `dim` mode + `vmap` index map avoids three sibling classes and three drift surfaces.

## 12. The shell path and the σ33 condensation

| `getType()` | order | components | host / role |
|---|---|---|---|
| `ThreeDimensional` (`3D`) | 6 | `{σxx σyy σzz τxy τyz τzx}` | 3D continuum (`LadrunoBrick`/`stdBrick`); reduce-to-`ASDConcrete3D` baseline |
| `PlaneStress` | 3 | `{σxx σyy τxy}` | 2D continuum / membrane quad; the clean MCFT oracle target — **not** a shell host |
| `PlateFiber` | 5 | `{σxx σyy τxy τyz τzx}` | **the shell path** — rides `LayeredShellFiberSection getCopy("PlateFiber")` |

The `PlaneStress`/`PlateFiber` views enforce `σ33 = 0` with a **guarded** nested scalar Newton on
`ε33`. This is the one place cloning the stock `PlateFiberMaterial` would be a bug (Quirk §22.14): on
the softening branch `dd22 → (1−d)E → 0`, so the bare `1/dd22` step explodes and the stock loop
**silently returns 0 when unconverged**. The RC view guards the step and **propagates a
non-convergence code** instead.

> [!important] State purity under the σ33 Newton (Quirk §22.15)
> The condensation re-calls the in-plane update ~25× per converged point. All crack/slip state
> therefore depends **only** on the membrane strains and is written only to the *output* history —
> idempotent across the inner re-calls — so the Newton cannot corrupt the committed crack frame.

## 13. The finite-strain companion (`LadrunoRCFiniteStrain`)

`LadrunoRCFiniteStrain` (classTag **33018**) is a **native `FiniteStrainNDMaterial` subclass**, not a
use of the generic `nDMaterial LogStrain` wrapper. It is driven by `setTrialF(F)` from an element
(`LadrunoBrick -geom finite`), and implements the Hencky seam of §8 by reusing the RC kernel's
`returnMap3D` verbatim + `LogStrainKernel.h` for the kinematics (`hencky_voigt`, `assemble_material`,
`spatial_tangent_full`).

> [!important] Why native, NOT the `LogStrain` wrapper (Quirk §22.12)
> `LogStrainNDMaterial` commits the elastic left Cauchy–Green by recovering `εᵉ = C₀⁻¹:τ` through the
> inner's **initial** tangent — exact **only for a non-degrading (plasticity) inner** (its own
> comment: *"v1 assumes a linear-elastic inner law"*). The RC spine is a **damage** law
> (`τ = (1−d)·Dᵉ:εᵉ`), so that recovery is shrunk by `(1−d)` and the committed `bᵉ` drifts as damage
> grows. The native view sidesteps it: the RC kernel carries **no tensorial plastic strain**, so the
> elastic left Cauchy–Green is the **total** `B = F Fᵀ` recomputed each step — no `bᵉ` to track, and
> **no committed finite state beyond the same `RCHist`** as the small-strain class. Simpler than
> `LadrunoJ2Finite` (no `bᵉ`, no co-rotation channel-B).

The command grammar is **identical** to `LadrunoRCConcrete` (Part III §17) minus any dimensional-view
choice. The objectivity split (§9) applies: isotropic spine objective, directional crack state xfail.

## 14. State, serialization and the IMPL-EX runtime

**Committed state.** The entire path-dependence is the `RCHist` struct (§11) plus the regularization
latch `{regularizationDone, regLch}` and the IMPL-EX runtime `{dtime_n, dtime_n_commit, dtime_0,
commitDone, implexError}`. The finite view adds **nothing** to this (its committed elastic reference
is `RCHist.strain` = the committed Hencky strain).

**The commit cycle.** `commitState` (for `-implex`) re-integrates *implicitly* at the converged
strain (the trial held the explicit/extrapolated response), measures `implexError = max|dt_ex −
dt_im|`, rolls `n → n−1` (`xt_old ← xt`, etc.), then `histN = histTr`. `revertToLastCommit` restores
`histTr = histN`; `revertToStart` zeros the history and resets the elastic tangent.

**Serialization** (`sendSelf`/`recvSelf`) carries a **hard-checked schema version** that rejects
mismatched wire vectors. The two classes have independent schemas:

| material | schema | wire size | versions |
|---|---|---|---|
| `LadrunoRCConcrete` (33015) | `RC_SCHEMA_VERSION = 5` | `RC_DATA` (params + 2 backbones + RCHist + cEps33) | v1 base · v2 IMPL-EX (242→262) · v3 shearRetFactor · v4 tension stiffening · v5 regularization |
| `LadrunoRCFiniteStrain` (33018) | `RCF_SCHEMA_VERSION = 1` | `RCF_DATA = 269` (params + 2 backbones + RCHist; no `dim`/`cEps33`) | v1 |

> [!important] Serialization discipline (Quirks §22.5, §22.16)
> (1) `crackShear`/`tauCr` **read trial state** (a ~1e-8 Penalty residual × stiff `G` ⇒ ~5e-4 drift
> from the committed value), so a bit-exact `save`/`restore` test must compare the **committed**
> fields and drop `crackShear`. (2) Any new committed field bumps the schema version and must be
> written **and** read (the `RC_DATA` arithmetic must balance) — add a **non-default-value**
> round-trip test, since a dropped wire slot reverts to the ctor default and a default-valued
> round-trip can't catch it. (3) `recvSelf` must **not** call `setupParams` — the cached `ftPeak`
> (tension stiffening) is serialized, not recomputed.

## 15. Registration and dispatch

An `nDMaterial` is wired into **four** sites (fewer than an element — the modern
`OpenSeesNDMaterialCommands` registry serves **both** Tcl and Python). To add or move a material,
mirror an existing sibling across all four:

1. `SRC/classTags.h` — the `#define ND_TAG_… <tag>` (the **single source** for tag-collision
   avoidance; the same number is recorded in `LEDGER_implementations.md` and `testbed/manifest.yaml`).
2. `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — the `#include` + the `case ND_TAG_…:`
   returning `new …()` (for `recvSelf` reconstruction).
3. `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` — the `extern void* OPS_…();` + the
   `nDMaterialsMap.insert({"Name", &OPS_…})` (the parser dispatch).
4. `SRC/material/nD/CMakeLists.txt` — the `.cpp` and `.h` in the source list (a new `.cpp` requires a
   CMake **re-configure**, not just ninja).

Every fork-authored file also carries the **LADRUNO header stamp** (`Ladruno_scripts/stamp_headers.py`
— add the new file to its `GLOBS`) and every shipped feature gets a **splash-banner** line
(`Ladruno_scripts/banner_features.txt` → `patch_banner.py`). The build-control ledgers
(`LEDGER_implementations`, `LEDGER_quirks`, `LEDGER_vanilla_files`) and the `testbed/manifest.yaml`
row (the CI G9 gate keys off it) are updated **in the same PR**.

## 16. The verification harness

Each phase ships a layered gate set, designed so that a wrong number cannot pass silently:

- **numpy oracle** — `tests/_testbed/rc_shell_ref.py` reimplements the kernel in numpy; the C++ is
  matched step-by-step. Tension-stiffening/regularization add `run_T1`/`run_R1` oracles.
- **standalone g++** — `tests/_testbed/rc_{tensstiff,reg}_gpp.cpp` build *the kernel header alone*
  (no OpenSees) and FD-check the tangent. **Key insight:** a converged static stress is the residual
  *root* ⇒ tangent-independent, so **no OpenSees Zone-A test can catch a wrong tangent** — the only
  tangent gate is the standalone g++ (wired into CI by `tests/test_ladrunoRCConcrete_{tensstiff,reg}_cpp.py`,
  compile+run, skip if no g++).
- **Zone-A pytest** (`tests/test_ladrunoRCConcrete_{material,shell,implex,objectivity,reg,meshobj}.py`,
  `tests/test_ladrunoRCFiniteStrain.py`) — material-point + shell-element, runs in CI on Ubuntu.
- **Zone-B pytest** (`tests/test_ladrunoRCConcrete_wall.py`) — gmsh-meshed wall, quasi-static
  explicit, runs nightly.

The two-zone contract (`tests/conftest.py`): a test is `@zone_a` (upstreamable, no gmsh) or `@zone_b`
(needs gmsh, auto-skipped where absent). The full V&V matrix is Appendix A.

---
---

# Part III — Usage

## 17. Command grammar

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
    [-shearRetention {mcft|const|dsfm|rots}] [-shearRetFactor $μ] \
    [-implex [-implexAlpha $a] [-implexControl $errTol $timeRedLim]] \
    [-tensStiff {vc|cm} [-tensStiffC $c] [-tensStiffAlpha $a]] \
    [-autoRegularization $lch_ref]

# identical grammar (minus the view choice); ThreeDimensional-only, for -geom finite:
nDMaterial LadrunoRCFiniteStrain $tag $E $nu  <same backbones + flags as above>
```

| Token | Meaning | Default |
|---|---|---|
| `$E $nu` | Young's modulus, Poisson's ratio | — (required) |
| `-Ce/-Cs/-Cd` | compression backbone strain / stress / damage (≥2; `Ce`,`Cs` equal length) | — (`Ce`,`Cs` required) |
| `-Te/-Ts/-Td` | tension backbone strain / stress / damage (≥2) | — (`Te`,`Ts` required) |
| `-Kc` | Lee–Fenves biaxial shape (`fb/fc` driver) | 2/3 |
| `-betaFloor` | floor on the MCFT `β` factor | 0.1 |
| `-rho` | mass density (also the element mass for explicit, §19) | 0.0 |
| `-beta` / `-lublinerReduced` | enable MCFT compression softening (§3) / reduce the Lubliner t–c interaction | OFF / OFF |
| `-secant` / `-numericalTangent` | secant / forward-difference tangent | algorithmic |
| `-interlock` | enable aggregate-interlock shear bound (§4.1) | OFF |
| `-agg` / `-crackStrain` / `-crackSpacing` / `-lch` / `-betaSrMin` | `a_g` / crack-freeze strain / spacing `s_θ` / explicit `lch` / residual shear floor | 16 / 0 / 0 / auto / 0.01 |
| `-cyclic` | enable cyclic friction-slip (⇒ `-interlock`, §4.2) | OFF |
| `-xcrack` | enable 2nd crack + cumulative-slip wear (⇒ `-cyclic`, §4.3) | OFF |
| `-degKappa` / `-degSlipRef` / `-degMin` | Archard wear: max knockdown / slip reference / floor | 0.5 / 0.01 / 0.1 |
| `-shearRetention` / `-shearRetFactor` | crack-shear slip-stiffness mode (§4.4) / `const`-mode `μ` | mcft / 0.4 |
| `-implex` / `-implexAlpha` / `-implexControl` | enable IMPL-EX (§7) / extrapolation multiplier / advisory `$errTol $timeRedLim` | OFF / 1.0 / off |
| `-tensStiff` / `-tensStiffC` / `-tensStiffAlpha` | tension stiffening `vc`|`cm` (§5) / `vc` coeff `c>0` (ignored in `cm`) / `cm` scale | OFF / 500 / 1.0 |
| `-autoRegularization` | enable Bažant–Oh crack-band at reference length `$lch_ref` (§6) | OFF |

> [!important] Flag implication chain + defaults-off
> `-xcrack` ⇒ `-cyclic` ⇒ `-interlock` (each law lives inside the previous block); `const`/`dsfm`
> retention ⇒ `-interlock -cyclic`; `rots` ⇒ `-interlock`. **All flags OFF ⇒ trajectory-faithful to
> `ASDConcrete3D`** (the reduce-to-baseline gate). Units must be consistent — the `v_ci,max =
> 0.18√fc'/…` MCFT form is calibrated in **SI (MPa, mm)** (Quirk §22.2).

## 18. Worked examples

### 18.1 Cyclic RC shell wall — OpenSeesPy

```python
import openseespy.opensees as ops

E, nu, Kc = 30000.0, 0.2, 2.0/3.0          # MPa
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
               "-cyclic", "-xcrack",                 # cyclic friction-slip + X-crack wear
               "-autoRegularization", 50.0)          # mesh-objective softening (lch_ref = 50 mm)

H, nlay = 0.1, 4                                       # layered shell: 4 layers
layers = []
for _ in range(nlay):
    layers += [1, H/nlay]
ops.section("LayeredShell", 1, nlay, *layers)
ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)           # add -corotational for large rotation
```

### 18.2 Membrane panel (PlaneStress view) — Tcl

```tcl
nDMaterial LadrunoRCConcrete 1 30000.0 0.2 \
    -Ce 0.0 0.0007 0.0020 0.0100  -Cs 0.0 24.0 30.0 5.0 \
    -Te 0.0 0.0001 0.0010         -Ts 0.0 3.0 0.5 \
    -beta -interlock -agg 16.0
# consumed by a 2D continuum element via its PlaneStress view
```

### 18.3 Finite-strain RC solid (`-geom finite`)

```tcl
nDMaterial LadrunoRCFiniteStrain 1 30000.0 0.2 \
    -Ce 0.0 0.0007 0.0020 0.0100 -Cs 0.0 24.0 30.0 5.0 \
    -Te 0.0 0.0001 0.0010 -Ts 0.0 3.0 0.5 -beta
element LadrunoBrick 1 1 2 3 4 5 6 7 8 1 -formulation std -geom finite
```

### 18.4 Plain baseline (all flags off ⇒ ≈ `ASDConcrete3D`)

```tcl
nDMaterial LadrunoRCConcrete 1 30000.0 0.2 \
    -Ce 0.0 0.0007 0.0020 0.0100 -Cs 0.0 24.0 30.0 5.0 -Cd 0.0 0.0 0.25 0.889 \
    -Te 0.0 0.0001 0.0010 -Ts 0.0 3.0 0.5 -Td 0.0 0.0 0.9
# trajectory-faithful to ASDConcrete3D (the reduce-to-baseline gate)
```

## 19. The cyclic quasi-static explicit recipe

A **cyclic** softening RC wall will **not** converge under implicit Newton (the plastic-damage
consistent tangent goes indefinite on the softening branch; `-implex` helps but does not fully cure a
meshed cyclic wall). The monotonic solvers (`LadrunoArcLength`, `LadrunoIndirectControl`,
`LadrunoDynamicRelaxation`, `robust_drive`) trace **one** equilibrium path through a limit point — a
load reversal is not one monotonic path, so they don't do cyclic. The right tool is **quasi-static
explicit** (`CentralDifferenceLadruno`): it forms no stiffness tangent, so the indefinite-tangent
stall simply doesn't exist and the reversals integrate through.

```python
import math, openseespy.opensees as ops
# build the meshed wall: ASDShellQ4 + LayeredShell(LadrunoRCConcrete … -beta -interlock -cyclic
#   -xcrack -rho 2.4e-9); fix the base; lock z + all rotations (planar membrane).
E, rho, h = 30000.0, 2.4e-9, 500.0                     # MPa, tonne/mm^3, element size mm
dt = 0.2 * h / math.sqrt(E / rho)                      # manual wave-speed bound (frac 0.1-0.2)
ops.constraints("Transformation"); ops.numberer("RCM")
ops.system("Diagonal"); ops.algorithm("Linear")       # trivial M^-1, exactly ONE solve/step
ops.integrator("CentralDifferenceLadruno", "-cflAbort", "-lump", "diagonal")
ops.analysis("Transient")
for _ in range(total_steps):                           # seg_time = steps_per_seg*dt MUST be >> T_struct
    ops.analyze(1, dt)                                 # (quasi-static: keep KE << strain energy)
```

> [!warning] Explicit gotchas (Quirks §22.6) — each is a real trap
> - **Element mass, not nodal:** give the material `-rho`; `ops.mass(...)` alone leaves the `dt_cr`
>   eigensolve with no element mass.
> - **`ASDShellQ4` reports no `dt_cr`** (`criticalTimeStep()` → −1); blind use gives `dt = 0.2·(−1) <
>   0` → instant blow-up. Use the manual `dt ≈ frac·h/√(E/ρ)`.
> - **No `equalDOF`/rigid ties** — the CD stability bound ignores constraints, so a rigid top via
>   `equalDOF` makes the bound a lie. Prescribe the rigid-top drift by the **same `sp` on every top
>   node**.
> - **Mass-proportional damping only** — stiffness-proportional (`betaK`) Rayleigh collapses `dt_cr`.

Validated in `tests/test_ladrunoRCConcrete_wall.py` (a 4×3 squat wall completes the full ±drift
schedule; the cyclic interlock degradation is load-bearing at panel scale — −28% energy, −11% peak vs
the monotone bound).

## 20. State recording

Through the element's `material` response (`eleResponse(ele, "material", gp, "<name>")` or
`recorder … material <gp> <name>`), or the section/shell `stresses` resultants:

| Response name(s) | Returns |
|---|---|
| `stress` / `strain` / `tangent` | current stress, strain, tangent of the active view (finite view: Cauchy σ, Hencky εᵉ) |
| `damage` | `(dt, dc)` tension/compression damage scalars |
| `beta` / `compressionSoftening` | the MCFT softening factor `β` |
| `eps1` / `transverseTensileStrain` | the principal transverse tensile strain `ε1` |
| `betaShear` / `shearRetention` | shear-retention / cap-utilization factor |
| `crackState` / `crackNormal` | crack-1 frozen normal + state |
| `crackShear` | crack-plane shear `τcr` + cap (reads the **trial** history — exclude from bit-exact round-trips) |
| `xcrackState` | crack-2 state + cumulative slip `slipCum` |
| `implexError` | the last IMPL-EX extrapolation error |

## 21. Units

There is no built-in unit system; **be consistent**. The one hard constraint: the MCFT crack-shear
limit `v_ci,max = 0.18√fc'/(0.31 + 24w/(a_g+16))` is **calibrated in SI (MPa, mm)** — the `√fc'` and
the `24w/(a_g+16)` term are not dimensionless, so `fc'`, `a_g`, `crackSpacing`, and the implied crack
width must be supplied in a matching MPa/mm system (Quirk §22.2). The examples use MPa, mm,
tonne/mm³, N (the apeGmsh/STKO convention).

---
---

# Part IV — Quirks

> The canonical, dated copies live in [[LEDGER_quirks]]; this is the RC-stack-relevant subset,
> grouped. Each maps to a `### ` heading there.

## 22. Consolidated quirks

**Spine / constitutive**
1. **`/E` measure + E-consistent backbone `q`.** Cloning the `ASDConcrete3D` spine has two
   silent-but-fatal requirements: the equivalent-strain measure needs the `/E` division (a β-ratio
   test can't catch its omission — only an absolute-stress test), and `q` is E-consistent via
   `buildBackbone`/`adjust()`, not `q = y/(1−d)` on raw points.
2. **`v_ci,max` is SI-unit-bound** (MPa, mm) — supply `fc'`/`a_g`/spacing in a matching system.
3. **Interlock engages only NON-proportionally** — under proportional shear the crack sits on the
   principal plane (`g_nt = 0`) and the bound/friction is inert.
4. **Pinching is a PANEL phenomenon** — a fixed-crack material point at constant normal strain gives a
   *fat* loop; the pinched waist needs principal rotation in a meshed wall.

**Tension stiffening / regularization (Phase 3)**
5. **TS equibiaxial self-consistency coefficient is 0.5, not 1** — the naive rank-1 `(0.5,0.5)` floor
   reaches only halfway; use distinct inject/measure normals (`ts_inj=(1,1,0)`, `ts_meas=(0.5,0.5,0)`).
   A pure-uniaxial gate can't catch it — an **equibiaxial** test is required.
6. **TS pinning tangent must span all 6 columns** of the floored rows (the floor pins `nᵀσn`
   independent of bare stress ⇒ removes its `eps_zz` column too — FD-caught).
7. **TS is monotonic-scope** — the live-`ε1` floor re-inflates on unload; use for pushover. TS+interlock
   (live-`p1` vs frozen-crack) is validated proportional-only.
8. **`regularize()` must re-apply `adjust()` after scaling** — omitting it lets post-peak plastic strain
   go non-monotone (q/damage off ~6% on steep backbones), **invisible to every energy gate**.
9. **`getCopy` must PROPAGATE the regularization latch** (else a copy of a used instance double-scales);
   a **loud-fail must NOT latch** (else a step-retry silently proceeds un-regularized).
10. **Structural Bažant-bar recipe** — localize via a thinner **section** (not a weaker material, so the
    band factor cancels); the un-regularized **fine** mesh snap-back-diverges, so the mesh-dependence
    negative control uses the **coarser** meshes.

**Finite strain (Phase 4b)**
11. **A KINEMATIC-hardening material wrapped in `LogStrain` is NOT objective under large rotation** —
    the backstress doesn't co-rotate (the §14.11 boundary; the RC crack frame is the analog).
12. **The generic `LogStrain` wrapper is UNSOUND for a damage inner** — it recovers `bᵉ` via the
    inner's *initial* tangent, which a stiffness-degrading inner shrinks by `(1−d)`. Lift damage
    materials with a **native `FiniteStrainNDMaterial` subclass** (the RC kernel has no tensorial
    plastic strain ⇒ `B = F Fᵀ` total, no `bᵉ` to track).

**Shell / condensation**
13. *(see §5 / Quirk 7 above for the TS-unload note.)*
14. **Don't clone the stock `PlateFiberMaterial`** for the `σ33=0` Newton — on softening `dd22→0` the
    bare `1/dd22` explodes and the stock loop **silently returns 0** when unconverged; guard the step
    and propagate a non-convergence code.
15. **σ33-condensation state purity** — the inner update is re-called ~25×/point, so all crack/slip
    state must depend only on the membrane strains and write only to the output history (idempotent).
16. **`crackShear` reads trial state** — it drifts ~5e-4 from the committed value; drop it from
    bit-exact serialization round-trips and compare the committed crack fields.

**Solver / explicit (cyclic walls)**
17. **`timeSeries Path` returns 0 beyond its last node** — float-accumulated pseudo-time overshoots and
    collapses the prescribed strain; pad a hold node (or use `-useLast`).
18. **Cyclic softening shell walls wall on implicit Newton past first crack** — the consistent tangent
    goes indefinite; the cyclic tool is quasi-static **explicit** (`CentralDifferenceLadruno`).
19. **`ASDShellQ4` has no per-element `dt_cr`** (`criticalTimeStep()` = −1) → manual `dt ≈ 0.2·h/√(E/ρ)`;
    element mass via material `-rho`; no `equalDOF`; mass-proportional damping only.
20. **Prescribing ALL of a node's DOFs via `sp` under the `Transformation` handler ⇒ 0 free equations ⇒
    the process terminates** (segfault). Use `Penalty` (or `Lagrange`) for fully-prescribed
    homogeneous-strain / finite probes.

**Build / CI**
21. **Build against the main checkout** (the build tree lives there, not in the worktree): copy changed
    files in, **bump their mtime** (Copy-Item preserves the old mtime ⇒ ninja skips the recompile),
    `build.bat OpenSeesPy` via the PowerShell tool. If the main checkout is staler than `ladruno`,
    surgically edit *its* registration files rather than overwriting with worktree versions (which
    reference sibling classes the staler checkout lacks).
22. **A fork PR that sits goes DIRTY** (ladruno auto-merges hourly): merge `ladruno` in and **backfill
    any manifest row a sibling PR left behind** (the G9 gate fails on the inherited gap). Base fork PRs
    on `ladruno`, and after a squash-merge **re-branch** for the next phase (don't stack).

---
---

# Appendices

## A. Phase and PR history + V&V matrix

| Phase | Ships | Pins | Tests | PRs |
|---|---|---|---|---|
| **1** | the spine + MCFT `β` | `β` on the strength axis (`|σc|=β·fc'` exact); reduce-to-`ASDConcrete3D` byte-match; biaxial β ratio | numpy + g++ + `…_material.py` | #155 #192 |
| **2a** | `-interlock` bound | `v_ci,max` cap vs closed-form; off-axis crack; ON/OFF ablation; g++ FD on cap + tangent | material battery | #239 |
| **2b.1** | `-cyclic` friction-slip | reversal re-caps; unload stiffness `G`; closure cap recovery; C++≡oracle | material battery | #245 |
| **2b.2a** | shell integration | end-to-end `ASDShellQ4`+`LayeredShell`: tension cracks; cyclic `Nxy` caps both signs; round-trip | `…_shell.py` | #246 |
| **2b.2b** | `-xcrack` + wear | governing cap; cumulative-slip decay; C++≡oracle with wear | material + shell | #253 |
| **2b.2c.1** | `-shearRetention` | `const`/`dsfm`/`rots` modes vs closed-form; all reduce to `mcft` | material + oracle | (2c.1) |
| **2b.2c.2** | objectivity gate | `σ(QEQᵀ)=Qσ(E)Qᵀ` over cracked/cyclic/X-cracked, 30°/90°/127° | `…_objectivity.py` | (2c.2) |
| **2b.2c.3** | crack-closure verify | per-step recompose = unilateral closure; compress-past-peak recovers | material | (2c.3) |
| **2b.2c.4a** | meshed-wall explicit | quasi-static explicit completes ±drift; degradation load-bearing (−28% E, −11% peak) | `…_wall.py` (Zone-B) | #266 |
| **4a** | `-implex` | off-identical; tracks implicit; error on rate change; SPD secant; save/restore | `…_implex.py` | #263 |
| **3a** | `-tensStiff` | post-crack `σ_xx==σ_ts(ε1)`; equibiaxial both normals; PlateFiber `Nxx==σ_ts·h`; FD (g++) | numpy T1 + g++ + `…_tensstiff.py` | #273 |
| **3b** | `-autoRegularization` | `g_reg·lch` mesh-objective; no-op at `lch=lch_ref`; steep-damage monotone (g++) | numpy R1 + g++ + `…_reg.py` | #277 |
| **3b-struct** | structural objectivity | Bažant `ASDShellQ4` bar: energy mesh-objective with reg, halves without (negative control) | `…_meshobj.py` | #281 |
| **4b** | finite-strain view | `σ==(small-strain RC at ½lnB)/J`; elastic `K==FD`; isotropic objectivity ×2; det F guard | `test_ladrunoRCFiniteStrain.py` | #282 |

> [!note] The load-bearing gate is β-on-the-strength-axis
> Proven three independent ways (numpy oracle, standalone g++ build of `LadrunoRCKernel.h`, end-to-end
> OpenSees). The forbidden abscissa-insertion **misses** the closed-form peak — the tripwire
> separating a correct MCFT from a knob that merely slides the lookup.

**Remaining (forward plan, [[rc_shell_phase3plus_handout]]):** the inclined/notched-panel structural
objectivity gate (staged Zone-B); the quantitative **Tran–Wallace RW-A20-P10** experiment calibration
(specimen geometry + reinforcement layout + digitized loops); **Phase 5 `LadrunoSolidShell`** (classTag
33020) — the through-thickness `σ33` host for punching/bearing (genuine state-dependent EAS-on-`E33`,
ANS/MITC, multi-layer Lobatto; it wanted the Phase-4b 3D material view first, now in place).

## B. References

**Theory**
- Vecchio & Collins (1986), *The Modified Compression-Field Theory…* — the `β` compression-softening
  law and `v_ci,max` crack-shear limit.
- Vecchio (2000), *Disturbed Stress Field Model* (DSFM) — the rotating + slip variant.
- Lubliner, Oliver, Oller & Oñate (1989); Lee & Fenves (1998) — the biaxial plastic-damage envelope.
- Bažant & Oh (1983) — crack-band regularization.
- Bentz (2000); Collins & Mitchell (1991) — tension-stiffening averaged-tension curves.
- Oliver, Huespe & Cante (2008) — IMPL-EX.
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* — §9.4 (plane-stress
  condensation hazard), Box 14.3 (Hencky MATISU), §14.11 (large-rotation objectivity boundary).

**Within this repo**
- [[19_ladruno_rc_shell_adr]] — the source ADR: decisions D1–D6, the seam analysis, the full risk
  register, the phased plan.
- [[rc_shell_phase3plus_handout]] — the forward plan (Phase 5 + remaining gates): goal/build/reuse/
  acceptance/gotchas/effort per item.
- [[LadrunoConcrete3D_user_guide]] — the `ASDConcrete3D` spine this material clones.
- [[Ladruno_materials_guide]] — the materials catalog.
- [[LadrunoRebarBuckling_guide]] — the discrete boundary-rebar layer partner.
- [[LogStrain_guide]] — the generic Hencky wrapper (which RC does **not** use — §13).
- [[LEDGER_quirks]] — the canonical, dated quirk register.
- Source: `SRC/material/nD/LadrunoRCKernel.h` (header-only kernel), `LadrunoRCConcrete.{h,cpp}`
  (small-strain views, classTag **33015**), `LadrunoRCFiniteStrain.{h,cpp}` (finite view, classTag
  **33018**, reuses the kernel + `LogStrainKernel.h`). Oracle: `tests/_testbed/rc_shell_ref.py`.
