---
title: "Lemaitre Coupled Ductile Damage — LadrunoJ2 / LadrunoUniaxialJ2"
project: Ladruno
type: reference / user guide
status: shipped (v1, 3D + uniaxial, MERGED PR #115; validated A–E + perf, 26/26)
classTag: none (mode — no new class tag)
material: "nDMaterial LadrunoJ2 / uniaxialMaterial LadrunoUniaxialJ2 — -damage lemaitre"
related:
  - "[[15_lemaitre_ductile_damage_adr]]"
  - "[[LadrunoJ2_guide]]"
  - "[[LadrunoUniaxialJ2_guide]]"
  - "[[Ladruno_materials_guide]]"
  - "[[11_brick_asdconcrete_integration]]"
  - "[[LadrunoBrick_reference]]"
tags:
  - material
  - plasticity
  - damage
  - fracture
  - lemaitre
  - reference
---

# Lemaitre Coupled Ductile Damage — `LadrunoJ2` / `LadrunoUniaxialJ2`

> [!abstract] One-line summary
> Lemaitre damage is an **optional, OFF-by-default mode** — `-damage lemaitre $r $s $pD $Dc`
> (`+ -implex`) — added to **both** the 3D `nDMaterial LadrunoJ2` (classTag 33011) and
> its uniaxial twin `uniaxialMaterial LadrunoUniaxialJ2` (classTag 33000). It threads a
> single scalar damage field $D$ through the **effective stress** $\tilde{\boldsymbol\sigma}=\boldsymbol\sigma/(1-D)$
> so the J2–Chaboche return map degrades stiffness **inside** the flow and the damage
> rate becomes **triaxiality-driven** in 3D — turning the spine into a model of ductile
> fracture / void growth / necking-driven rupture. It is **not** a new class and has
> **no new classTag** (decision D1): with `-damage none` (the default) the material takes
> the exact undamaged code path, **bit-identical** to today.

This document is the **descriptive reference** for the damage mode: the continuum theory
(strain-equivalence, effective stress, the triaxiality function), the coupled return map
and tangent exactly as designed, the OpenSees wiring, the command syntax, intended uses,
and the honest limits. For the chronological decision log (D1–D4), the full risk register,
and the implementation plan see the source ADR [[15_lemaitre_ductile_damage_adr]]. For the
undamaged spine it builds on, see [[LadrunoJ2_guide]] (§10.2 there summarizes this mode) and
[[LadrunoUniaxialJ2_guide]].

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why a coupled mode, not a wrapper|2. Why a coupled mode, not a wrapper]]
- [[#3 Continuum theory|3. Continuum theory]]
- [[#4 The coupled return map|4. The coupled return map]]
- [[#5 The consistent tangent|5. The consistent tangent]]
- [[#6 Regularization and IMPL-EX|6. Regularization and IMPL-EX]]
- [[#7 OpenSees implementation|7. OpenSees implementation]]
- [[#8 Command syntax and usage|8. Command syntax]]
- [[#9 State recording and failure|9. State recording and failure]]
- [[#10 Verification and validation|10. Verification]]
- [[#11 Limitations and boundaries|11. Limitations]]
- [[#12 References|12. References]]

---

## 1. Intended use cases

Lemaitre damage is **continuum damage mechanics for ductile metals**: it predicts the
progressive loss of stiffness and load capacity as microvoids nucleate, grow, and coalesce
under plastic flow. The canonical targets:

| Use case | Why the Lemaitre mode |
|---|---|
| **Multiaxial ductile fracture (3D, the headline)** — notched bars, fillets, weld toes, plate tearing | The damage rate is driven by the **damage energy release rate $Y$**, which is **explicitly triaxiality-dependent** ($R_v(\sigma_H/\sigma_{eq})$). Confinement/notch acuity accelerates void growth exactly as real steel does — something a strain-watching wrapper cannot see. |
| **Necking-driven rupture / void growth** | The effective-stress coupling $\tilde{\boldsymbol\sigma}=\boldsymbol\sigma/(1-D)$ degrades stiffness *and* accelerates plastic flow ($\dot{\bar\varepsilon}^p=\dot\gamma/(1-D)$), producing the softening branch and localization that precede ductile failure. |
| **Low-cycle fatigue (LCF)** | Cyclic accumulation of $p$ past the threshold $p_D$ drives $\dot D=(Y/r)^s\dot p$ — validated against Coffin–Manson / shakedown behaviour (leg E). |
| **Energy-based fiber damage (1D, degenerate point)** | `LadrunoUniaxialJ2 -damage lemaitre …` gives a triaxiality-frozen ($\sigma_H/\sigma_{eq}\equiv\tfrac13$) energy-based damage law for fiber sections — the valid uniaxial *oracle*, not a multiaxial fracture model. |

> [!tip] Complementary to — not competing with — a fatigue wrapper
> Lemaitre is the **mechanistic, coupled** member of the damage family: continuum damage
> mechanics, triaxiality-driven, monotonic-ductile-tearing oriented. The parked
> `LadrunoFatigue` **rainflow** wrapper (phenomenological Coffin–Manson / Miner, cheap,
> symmetric-cyclic-seismic oriented) wraps the *undamaged* core. They are designed to
> **coexist** — Lemaitre is a core *mode*, the fatigue index stays a *wrapper*. v1 ships
> Lemaitre only.

> [!warning] When **not** to reach for it
> Lemaitre damage inherits all of LadrunoJ2's boundaries: it is a **pressure-independent
> (deviatoric J2)** plasticity model — not for concrete/soil/rock (use ASDConcrete3D /
> Drucker–Prager). It is **ductile** damage (void growth under plastic flow), not brittle
> cracking. And it does **not** erode elements — failure is capped and flagged, not deleted
> (§9).

---

## 2. Why a coupled mode, not a wrapper

Before this mode, the fork had a clean **undamaged** combined-hardening J2 spine
([[LadrunoJ2_guide]], [[LadrunoUniaxialJ2_guide]]): Voce + linear isotropic + Chaboche
Armstrong–Frederick kinematic hardening, a scalar-Newton return map
($\mathbf n=\mathbf M(\Delta\gamma)/\lVert\mathbf M(\Delta\gamma)\rVert$, Kobayashi–Ohno),
an analytic consistent tangent, and a committed IMPL-EX hook (`dGamma_n`). The J2 ADR
explicitly foresaw a softening/damage extension. Lemaitre is that continuation, and the ADR
chose a **coupled core mode** over a post-hoc wrapper for concrete reasons:

| | A `Fatigue`/`MinMax`-style wrapper | **Lemaitre coupled damage (this mode)** |
|---|---|---|
| coupling | post-hoc: watches $\Delta\varepsilon^p$, scales/zeros stress at $D=1$ | **in the return map**: $\tilde{\boldsymbol\sigma}=\boldsymbol\sigma/(1-D)$ degrades stiffness *and* accelerates plastic flow |
| driver | plastic-strain range / dissipated energy (rainflow) | **damage energy release rate $Y$** — explicitly **triaxiality-dependent** in 3D |
| 3D fracture | no — a wrapper cannot see $\sigma_H/\sigma_{eq}$ | **yes** — notch/triaxiality-driven ductile fracture, the headline |
| applies to | uniaxial fibers only | **both** `LadrunoJ2` (3D continuum) and `LadrunoUniaxialJ2` (fiber) — shared law |
| book home | phenomenological (Coffin–Manson / Miner) | **dSNPO §12.3** — AF + isotropic + Lemaitre, sitting right there |

The decisive property is the one [[LadrunoJ2_guide]] §10.2 highlights: the return map runs
on the **undamaged effective stress** (which carries no $D$), so the state update
**decouples exactly**. With damage OFF the material is byte-identical to the spine — the
**V0 tripwire** — and all the $(\Delta\gamma, D)$ coupling lives in the *tangent*. The
continuum element [[LadrunoBrick_reference]] gets ductile fracture **for free**: the damaged
material plugs into the exact same `setTrialStrain`/`getTangent` contract, with no element
change.

> [!note] Design decision D1 — a mode, not a sibling class
> The ADR locked damage as a **parameter mode** on the existing registered classes, not as
> new `LadrunoJ2Damage`/`LadrunoUniaxialJ2Damage` siblings. The payoff: **no `classTags.h`
> churn**, no duplicated kernel to drift, and every existing test stays bit-identical. The
> rejected sibling-class route would have needed two new class tags plus broker/registry
> entries — exactly the churn the mode decision avoids.

---

## 3. Continuum theory

Reference: dSNPO (de Souza Neto–Perić–Owen 2008) **§12.3** (the coupled model) and
**§12.4** (the integration algorithm); Lemaitre & Chaboche (1990); Lemaitre (1985). The
guiding idea is that this is the *existing* `LadrunoJ2` with **one scalar field $D$ threaded
through the effective stress** — the return-map skeleton is reused, not rebuilt.

> [!info] Notation
> Direct (bold), index, and Voigt forms are shown side by side, matching
> [[LadrunoJ2_guide]] §3. Damage $D$ is a scalar; $\boldsymbol\sigma$ is the **nominal**
> (actual) Cauchy stress carried by the cracked material, $\tilde{\boldsymbol\sigma}$ the
> **effective** stress carried by the still-intact ligament. $\mathbf s=\operatorname{dev}\boldsymbol\sigma$,
> $\boldsymbol\alpha$ the total backstress, $p\equiv\bar\varepsilon^p$ the accumulated
> equivalent plastic strain.

### 3.1 State, free energy, effective stress (strain equivalence)

The internal variables are those of the spine — elastic strain $\boldsymbol\varepsilon^e$,
accumulated plastic strain $p$, backstresses $\boldsymbol\alpha_k$ — **plus one scalar
damage** $D\in[0,D_c]$ ($D=0$ virgin, $D_c\approx0.2$–$0.5$ at rupture). The Helmholtz free
energy splits into a **damaged** elastic part and the (undamaged) plastic/hardening part:

$$\rho\psi = (1-D)\,\tfrac12\,\boldsymbol\varepsilon^e:\mathbf D^e:\boldsymbol\varepsilon^e \;+\; \rho\psi_p(p,\boldsymbol\alpha_k).$$

The **strain-equivalence** hypothesis (the strain of the damaged material under nominal
stress equals the strain of the virgin material under effective stress) means the $(1-D)$
factor multiplies the elastic stress:

- **Direct:** $\boldsymbol\sigma = (1-D)\,\mathbf D^e:\boldsymbol\varepsilon^e$, with the **effective stress** $\;\tilde{\boldsymbol\sigma}=\dfrac{\boldsymbol\sigma}{1-D}=\mathbf D^e:\boldsymbol\varepsilon^e$.
- **Voigt:** $\{\boldsymbol\sigma\}=(1-D)\bigl[\,K\,\mathbf 1\otimes\mathbf 1 + 2G\,\mathbb I_{\text{dev}}\,\bigr]\{\boldsymbol\varepsilon^e\}$.

> [!important] The effective stress is damage-independent at the trial state
> Because $\tilde{\boldsymbol\sigma}=\mathbf D^e:(\boldsymbol\varepsilon_{n+1}-\boldsymbol\varepsilon^p_n)$
> carries **no $D$**, the elastic trial predictor is unchanged from the spine. This is the
> single structural fact that lets the whole undamaged machinery survive — see §4.

### 3.2 Yield surface in *effective* stress space

Yield is evaluated on the **effective relative deviatoric stress**
$\tilde{\boldsymbol\xi}=\operatorname{dev}(\tilde{\boldsymbol\sigma})-\boldsymbol\alpha=\dfrac{\mathbf s}{1-D}-\boldsymbol\alpha$:

- **Direct:** $\Phi = \lVert\tilde{\boldsymbol\xi}\rVert - \sqrt{\tfrac23}\,\sigma_y(p)$, with $\sigma_y(p)=\sigma_0 + Q_\infty(1-e^{-bp}) + H_{\text{iso}}\,p$ (Voce + linear, shared `LadrunoHardening.h`).
- **Index:** $\Phi = \sqrt{\tilde\xi_{ij}\tilde\xi_{ij}} - \sqrt{\tfrac23}\,\sigma_y$, with $\tilde\xi_{ij}=s_{ij}/(1-D)-\alpha_{ij}$.
- **1D reduction:** $\Phi = \lvert\,\sigma/(1-D) - X\,\rvert - \sigma_y(p)$ (the $\sqrt{2/3}$ divides out, exactly as in [[LadrunoUniaxialJ2_guide]]).

Crucially, **damage scales the deviator magnitude, not its direction**. The
$\mathbf M(\Delta\gamma)$ radial machinery of the spine therefore carries over unchanged
(§4).

### 3.3 Flow, hardening, and damage evolution

| | form |
|---|---|
| plastic flow | $\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\dfrac{\partial\Phi}{\partial\boldsymbol\sigma} = \dfrac{\dot\gamma}{1-D}\,\mathbf N$, $\;\mathbf N=\tilde{\boldsymbol\xi}/\lVert\tilde{\boldsymbol\xi}\rVert$ — the $1/(1-D)$ **accelerates** plastic flow as $D$ grows |
| accumulated | $\dot p = \sqrt{\tfrac23}\lVert\dot{\boldsymbol\varepsilon}^p\rVert = \dot\gamma/(1-D)$ |
| isotropic | $\sigma_y(p)$ (Voce + linear, shared `LadrunoHardening.h`) |
| kinematic (AF) | $\dot{\boldsymbol\alpha}_k = \tfrac23 C_k\,\dot{\boldsymbol\varepsilon}^p_{\text{dir}} - \gamma_k\,\boldsymbol\alpha_k\,\dot p$ (unchanged form, written in effective space) |
| **damage** | $\dot D = (Y/r)^s\,\dot p$ for $p\ge p_D$, else $\dot D = 0$ (the damage threshold) |

### 3.4 The damage energy release rate $Y$ — the triaxiality carrier

The damage rate is driven by the **damage energy release rate** $Y$ (dSNPO §12.3):

- **Direct:** $Y = \tfrac12\,\boldsymbol\varepsilon^e:\mathbf D^e:\boldsymbol\varepsilon^e = \dfrac{\tilde\sigma_{eq}^2}{2E}\,R_v$.
- **Triaxiality function:** $\;R_v = \tfrac23(1+\nu) + 3(1-2\nu)\left(\dfrac{\sigma_H}{\sigma_{eq}}\right)^2$, with $\sigma_H=\tfrac13\operatorname{tr}\boldsymbol\sigma$ (hydrostatic) and $\sigma_{eq}=\sqrt{\tfrac32}\lVert\mathbf s\rVert$ (von Mises).

The ratio $\sigma_H/\sigma_{eq}$ is the **stress triaxiality** — high under confinement / at
notch roots, where $R_v$ (and hence the damage rate) is largest. This is *the*
differentiator: in 3D the damage rate genuinely depends on the stress state, not just on the
plastic-strain magnitude.

> [!note] The 1D degenerate point
> A `uniaxialMaterial` has **no triaxiality** — under pure uniaxial stress
> $\sigma_H/\sigma_{eq}\equiv\tfrac13$ *always*, so $R_v=\tfrac23(1+\nu)+\tfrac13(1-2\nu)=\text{const}$
> and $Y=\tilde\sigma^2/(2E)$ up to that constant. The 1D Lemaitre mode is therefore the
> valid **degenerate / energy-based** point and oracle — **do not market it as multiaxial
> fracture** (decision D3 caveat). The multiaxial fracture story lives only in 3D; the
> single-PR design just keeps both consistent at triaxiality $\tfrac13$.

### 3.5 The four damage constants

| Symbol | Name | Meaning |
|---|---|---|
| $r$ | energy denominator (sometimes $S$) | sets the damage *rate* — larger $r$ ⇒ slower damage growth |
| $s$ | exponent | shapes the $(Y/r)^s$ nonlinearity of the damage law |
| $p_D$ | plastic-strain threshold | no damage accrues until $p\ge p_D$ (a strain "incubation" gate) |
| $D_c$ | critical / rupture damage | damage is **capped** here and the point is flagged failed ($\approx0.2$–$0.5$) |

---

## 4. The coupled return map

> [!note] Decision D2 — fully-coupled backward-Euler, $D$ live
> The numerical scheme solves the $(\Delta\gamma, D)$ system **simultaneously** — $D$ is
> live inside every Newton iterate, *not* frozen. Because damage preserves the radial
> $\mathbf M(\Delta\gamma)$ direction (§3.2), this is the **same scalar-Newton-on-$\Delta\gamma$
> skeleton the spine already ships, promoted to a 2×2** — it is *not* a from-scratch tensor
> solve. The operator-split (frozen-$D$) scheme is kept only as the Newton **seed and a
> robustness fallback**, not the default.

Backward-Euler, elastic predictor / plastic corrector. The structural fact that makes the
coupling cheap: the **effective trial stress is damage-independent**
($\tilde{\boldsymbol\sigma}^{\text{tr}}=\mathbf D^e:(\boldsymbol\varepsilon_{n+1}-\boldsymbol\varepsilon^p_n)$,
no $D$), and damage scales the deviator *magnitude* not its *direction*. Damage therefore
adds **exactly one extra unknown** $D$ and **one extra equation** (the integrated evolution
law) to the shipped scalar Newton.

**1 — Trial.**
$$\tilde{\boldsymbol\sigma}^{\text{tr}}=\mathbf D^e:(\boldsymbol\varepsilon_{n+1}-\boldsymbol\varepsilon^p_n),\quad
\tilde{\boldsymbol\xi}^{\text{tr}}=\operatorname{dev}(\tilde{\boldsymbol\sigma}^{\text{tr}})-\boldsymbol\alpha_n,\quad
\Phi^{\text{tr}}=\lVert\tilde{\boldsymbol\xi}^{\text{tr}}\rVert-\sqrt{\tfrac23}\,\sigma_y(p_n).$$

**2 — Elastic branch.** If $\Phi^{\text{tr}}\le0$: $\;\boldsymbol\sigma=(1-D_n)\,\tilde{\boldsymbol\sigma}^{\text{tr}}$,
$\;\mathbf D^{\text{alg}}=(1-D_n)\,\mathbf D^e$. Done — no damage growth without plastic flow.

**3 — Plastic corrector, 2-unknown Newton in $(\Delta\gamma, D)$:**

- **Direction (unchanged from the spine):** $\;\mathbf n = \mathbf M(\Delta\gamma)/\lVert\mathbf M(\Delta\gamma)\rVert$, with the Kobayashi–Ohno tensor
  $$\mathbf M(\Delta\gamma) = \tilde{\mathbf s}^{\text{tr}} - \sum_k\frac{\boldsymbol\alpha_{k,n}}{1+\sqrt{\tfrac23}\gamma_k\Delta\gamma},\qquad \tilde{\mathbf s}^{\text{tr}}=\operatorname{dev}(\tilde{\boldsymbol\sigma}^{\text{tr}}).$$
- **Plastic-strain increment** carries the damage acceleration: $\;\Delta p = \sqrt{\tfrac23}\,\dfrac{\Delta\gamma}{1-D}.$
- **R₁ — consistency (effective space):** $\;\lVert\tilde{\boldsymbol\xi}(\Delta\gamma)\rVert - \sqrt{\tfrac23}\,\sigma_y(p_n+\Delta p) = 0.$
- **R₂ — integrated damage:** $\;D - D_n - \bigl(Y(\Delta\gamma)/r\bigr)^s\,\Delta p = 0$ for $p_n+\Delta p\ge p_D$, else $D-D_n=0$; with $Y(\Delta\gamma)=\dfrac{\tilde\sigma_{eq}(\Delta\gamma)^2}{2E}\,R_v(\sigma_H/\sigma_{eq})$.

Solve $[\,R_1, R_2\,]^{\mathsf T}=0$ for $(\Delta\gamma, D)$ by Newton with an analytic 2×2
Jacobian: every $\partial/\partial\Delta\gamma$ term is the spine's $\mathrm dM/\mathrm d\Delta\gamma$
extended by the $1/(1-D)$ factors, and the $\partial/\partial D$ terms couple through
$\Delta p$ and the degradation. **Seed** $\Delta\gamma$ from the radial / frozen-$D$
(operator-split) estimate and $D$ from $D_n$; **bracket** $D\le D_c$.

**4 — Degrade.** $\;\boldsymbol\sigma_{n+1}=(1-D_{n+1})\,\tilde{\boldsymbol\sigma}_{n+1}$, clamp $D_{n+1}$ at $D_c$.

> [!check] Why this is barely more than the spine
> It is the *same scalar-Newton-on-$\Delta\gamma$ skeleton*, promoted to a 2×2 by adding
> $D$ and the residual $R_2$ (dSNPO §12.4). The radial reduction the spine already exploits
> is exactly what keeps the coupled system small — there is **no return to a full tensor
> solve**. The operator-split (freeze $D=D_n$ over the return, then update $D$ explicitly)
> survives as the Newton seed and as the robustness fallback when the coupled 2×2 fails to
> converge near rupture ($D\to D_c$, $Y$ large).

---

## 5. The consistent tangent

For implicit global Newton the material must return the **algorithmic** tangent
$\mathbf D^{\text{alg}}=\partial\boldsymbol\sigma/\partial\boldsymbol\varepsilon_{n+1}$
consistent with the *fully-coupled* $(\Delta\gamma, D)$ solve. Differentiating the converged
coupled state yields the degraded elastoplastic operator **plus** a damage-sensitivity
rank-one term:

$$\boxed{\;\mathbf D^{\text{alg}} = (1-D_{n+1})\,\mathbf D^{\text{alg}}_{\text{plastic}} \;-\; \tilde{\boldsymbol\sigma}_{n+1}\otimes\frac{\partial D_{n+1}}{\partial\boldsymbol\varepsilon_{n+1}}\;}$$

where $\mathbf D^{\text{alg}}_{\text{plastic}}$ is the undamaged spine tangent (§5 of
[[LadrunoJ2_guide]]), and $\partial D_{n+1}/\partial\boldsymbol\varepsilon$ (together with
$\partial\Delta\gamma/\partial\boldsymbol\varepsilon$) comes from **inverting the 2×2 Newton
Jacobian** of §4 via the implicit-function theorem on $[R_1, R_2]=0$. The tangent is thus
**consistent with the coupled solve**, not a frozen-$D$ approximation.

- The $(1-D)\,\mathbf D^{\text{alg}}_{\text{plastic}}$ block is the elastoplastic stiffness, scaled down by the current damage.
- The $-\,\tilde{\boldsymbol\sigma}\otimes\partial D/\partial\boldsymbol\varepsilon$ rank-one term carries the **softening**: as damage rises with strain it subtracts stiffness along the effective-stress direction. It is what makes the global tangent *indefinite* on a softening branch — the reason the IMPL-EX escape (§6) exists.

> [!check] Exact reduction to the undamaged tangent (the cheapest regression, V0)
> With damage off ($r\to\infty$ or $p_D\to\infty$ ⇒ $R_2$ collapses to $D=D_n=0$), the
> rank-one term vanishes and $\mathbf D^{\text{alg}}$ reduces **bit-identically** to the
> shipped `LadrunoJ2` tangent. A finite-difference check $\lVert\mathbf D^{\text{alg}}-\mathbf D^{\text{FD}}\rVert$
> at several *damaged* states (test V5) confirms consistency away from $D=0$.

---

## 6. Regularization and IMPL-EX

> [!note] Decision D4 — characteristic-length + a real IMPL-EX path
> Softening localizes, and a *local* damage model is mesh-dependent. The mode addresses
> this with **two** complementary tools: characteristic-length regularization (for
> mesh-objective dissipated energy) and a genuine IMPL-EX code path (for tangent
> robustness).

### 6.1 Characteristic-length regularization

The damage evolution is scaled by the element's characteristic length, **reusing the
`getCharacteristicLength()` + [[LadrunoBrick_reference]] `lch` handshake already built for
ASDConcrete3D** ([[11_brick_asdconcrete_integration]]). The post-threshold softening is
calibrated so the **dissipated fracture energy** $G_f$ is mesh-objective: refine the mesh and
the global force–displacement curve and the energy to rupture stay (approximately) the same
(tests V6 / V7).

> [!warning] Honest disclaimer — energy, not width
> `lch` regularization makes the **dissipated energy** mesh-objective but **not the
> localization width**. The damage band still narrows to one element row as the mesh
> refines. True width-objectivity needs **nonlocal or gradient-enhanced damage** — a future
> ADR. This caveat is deliberately surfaced (in the banner and docs), never silently
> dropped.

### 6.2 IMPL-EX — `-implex`

Both J2 classes already commit the plastic multiplier $\Delta\gamma_n$ (the `dGamma_n` hook).
The **IMPL-EX** path is the first real consumer of that hook: it extrapolates
$\tilde{\Delta\gamma}=\Delta\gamma_n\cdot(\Delta t_{n+1}/\Delta t_n)$, freezes the flow
direction and $D$ from the extrapolation, and forms a **constant SPD tangent** — exactly the
robustness escape for the indefinite softening tangent of §5. It pairs naturally with explicit
dynamics (`CentralDifferenceLadruno` / `ExplicitBathe`), where softening + implicit Newton
would otherwise stall. Enable it with the optional `-implex` flag.

---

## 7. OpenSees implementation

### 7.1 Where the code lives

| Item | Value |
|---|---|
| Mode on (3D) | `nDMaterial LadrunoJ2` — classTag **33011** (unchanged; **no new tag**) |
| Mode on (1D) | `uniaxialMaterial LadrunoUniaxialJ2` — classTag **33000** (unchanged) |
| Shared damage law | `SRC/material/nD/LadrunoDamage.h` — header-only, OpenSees-free (mirrors `LadrunoHardening.h`) |
| 3D material | `SRC/material/nD/LadrunoJ2.{h,cpp}` (`damageModel` enum, default `NONE`) |
| 1D material | `SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}` ($R_v=\text{const}$) |
| Parsers | `OPS_LadrunoJ2` / `OPS_LadrunoUniaxialJ2` — parse `-damage lemaitre $r $s $pD $Dc <-implex>` |

> [!note] The shared `LadrunoDamage.h` contract
> The damage kinematics — `damageReleaseRate(sigEff_eq, sigH, nu, E) → Y` (with the $R_v$
> triaxiality function: 3D passes the real $\sigma_H/\sigma_{eq}$ ratio, 1D passes the frozen
> $\tfrac13$) and `damageIncrement(Y, dp, r, s, pD, p) → ΔD` plus its derivatives — live in a
> tiny header-only file that **both** `LadrunoJ2.cpp` and `LadrunoUniaxialJ2.cpp` `#include`.
> This makes the 3D↔1D damage reduction **byte-identical at triaxiality $\tfrac13$**, extending
> the spine's shared-hardening oracle from "isotropic backbone" to "isotropic backbone +
> damage" (test V4). The instant a 3D edit bypasses this header, the V4 1e-12 gate fails — the
> oracle-drift tripwire.

### 7.2 What the mode adds to the spine

Everything is **additive and OFF by default**. The mode adds: a `damageModel` enum (default
`NONE`), the four params $r, s, p_D, D_c$, one committed scalar $D_n$ (+ trial $D$), the §4
coupled branch inside `integrate()` (guarded so `NONE` takes the *exact* undamaged path), the
§5 tangent term, $D$ into `sendSelf`/`recvSelf`/`revert*`, and the `"damage"`/`"failed"`/`"Y"`
responses. State grows by **one committed scalar** $D_n$ plus the four params — no other layout
change.

> [!important] House hygiene (same as the spine)
> **No `exit()`** anywhere. Real `revert*`, `setParameter`, `Print`, `sendSelf`/`recvSelf`. And
> the **initial tangent for modal / eigen analysis is the virgin $\mathbf D^e$** (since $D=0$ at
> the start) — `getInitialTangent` must not return the degraded operator.

---

## 8. Command syntax and usage

### 8.1 Grammar

The damage mode is an **option block** appended to the normal `LadrunoJ2` /
`LadrunoUniaxialJ2` command (see [[LadrunoJ2_guide]] §7 for the full base grammar):

```tcl
-damage lemaitre  $r $s $pD $Dc   <-implex>
```

| Token | Meaning |
|---|---|
| `$r` | energy denominator (sometimes $S$) — larger ⇒ slower damage |
| `$s` | damage-law exponent |
| `$pD` | plastic-strain threshold below which no damage accrues |
| `$Dc` | critical/rupture damage (cap + failure flag, $\approx0.2$–$0.5$) |
| `-implex` | optional — use the IMPL-EX (constant-SPD-tangent) path for softening robustness |

> [!important] `-damage none` (the default, or simply omitting the block) ⇒ identical to today.
> $D\equiv0$, the exact undamaged code path, zero overhead, every existing test bit-identical.

### 8.2 Examples

```tcl
# --- 3D: combined hardening + Lemaitre ductile damage (triaxiality-driven) ---
nDMaterial LadrunoJ2 1  $K $G  -iso voce 350.0 60.0 10.0 0.0 \
                               -kin 3  60000.0 500.0  12000.0 120.0  2000.0 10.0 \
                               -damage lemaitre $r $s $pD $Dc

# --- 3D with IMPL-EX (softening robustness for explicit/quasi-static rupture) ---
nDMaterial LadrunoJ2 2  $K $G  -iso voce 350.0 60.0 10.0 0.0 \
                               -damage lemaitre $r $s $pD $Dc -implex

# --- uniaxial fiber: degenerate energy-based damage (triaxiality fixed at 1/3) ---
uniaxialMaterial LadrunoUniaxialJ2 3  $E  -iso voce 350.0 60.0 10.0 0.0 \
                                          -kin 1 20000.0 100.0 \
                                          -damage lemaitre $r $s $pD $Dc
```

```python
# OpenSeesPy — 3D ductile-damage material consumed by a LadrunoBrick element
import openseespy.opensees as ops
E, nu = 200000.0, 0.3
K, G = E/(3*(1-2*nu)), E/(2*(1+nu))
ops.nDMaterial("LadrunoJ2", 1, K, G,
               "-iso", "voce", 350.0, 60.0, 10.0, 0.0,
               "-kin", 3, 60000.0, 500.0, 12000.0, 120.0, 2000.0, 10.0,
               "-damage", "lemaitre", r, s, pD, Dc)
ops.element("LadrunoBrick", 1, *node_tags, 1)
```

> [!tip] Finite-strain ductile fracture
> Lemaitre damage composes with the Hencky log-strain lift. For large-strain ductile
> tearing use isotropic hardening (`-kin 0`) under the wrapper, since combined hardening is
> non-objective at large rotation (see [[LadrunoJ2_guide]] §10.1):
> ```tcl
> element → LadrunoBrick(-geom finite) → LogStrain( LadrunoJ2 -damage lemaitre … -kin 0 )
> ```

---

## 9. State recording and failure

### 9.1 Damage responses

Through the element's `material` response (e.g.
`eleResponse(ele, "material", gp, "<name>")` or `recorder … material <gp> <name>`), the mode
adds three responses on top of the spine's stress/strain/tangent/backStress/plasticStrain/
equivalentPlasticStrain:

| Response name(s) | Returns |
|---|---|
| `damage` / `D` | scalar Lemaitre damage $D$ |
| `failed` / `rupture` | $1$ if (damage on **and** $D\ge D_c$), else $0$ |
| `Y` / `energyReleaseRate` | the damage energy release rate $Y$ (effective) |

The uniaxial twin keeps its existing `plasticWork` (`Cwp`) accumulator as the energy-based
diagnostic.

### 9.2 Failure handling — cap and flag, no erosion

At rupture the damage is **clamped at $D_c$** and the point is flagged via the `"failed"`
response; the stress stays degraded. **Element erosion / deletion is explicitly out of scope**
— that is an element-level concern, not a material one. A failed integration point keeps
carrying a (small, capped) stress rather than vanishing.

---

## 10. Verification and validation

The mode shipped (PR #115) and was **validated A–E + perf, 26/26**, with the reproducibility
bundle in `Ladruno_implementation/lemaitre_validation/` (README + figures + code). The
Zone-A oracle matrix (V0–V9, source ADR §6) is the design contract; the highlights:

| # | Test | What it pins |
|---|---|---|
| **V0** | damage-off regression (`-damage none`) | **bit-identical** to the shipped `LadrunoJ2`/`LadrunoUniaxialJ2` batteries — the cheap tripwire |
| V1 | uniaxial monotone tension, damage on | semi-analytic Lemaitre coupon $\sigma(\varepsilon)$ with softening onset at $p_D$ |
| V2 | damage-off via $p_D=\infty$ mid-run | proves the threshold gate; pre-threshold response equals the undamaged one |
| **V3** | **triaxiality sweep (3D)**: uniaxial vs notched/confined at fixed $p$ | damage-rate ratio matches the $R_v(\sigma_H/\sigma_{eq})$ closed form — **the differentiator** |
| **V4** | **3D↔1D damage reduction** at triaxiality $\tfrac13$, identical $(r,s,p_D,D_c)$ | `LadrunoJ2` condensed-to-uniaxial vs `LadrunoUniaxialJ2` ~1e-12 (the shared `LadrunoDamage.h` contract) |
| V5 | consistent-tangent FD at several damaged states | $\lVert\mathbf D^{\text{alg}}-\mathbf D^{\text{FD}}\rVert<\text{tol}$ |
| **V6** | **mesh-objectivity**: notched bar, 2 mesh sizes, `LadrunoBrick`+`LadrunoJ2` | same global force–displacement with `lch` regularization (mirrors the ASDConcrete test) |
| V7 | dissipated energy at rupture ($D\to D_c$) | $=$ prescribed fracture energy $G_f$ (mesh-objective via `lch`) |
| V8 | IMPL-EX vs implicit on a softening step | bounded lag; SPD tangent; no Newton stall |
| V9 | `sendSelf`/`recvSelf` round-trip incl. $D$ | parallel / serialization hygiene |

The validation bundle covers: A — a uniaxial oracle (from-equations integrator vs closed
form); B — the triaxiality $R_v(\eta)$ law; C — a gmsh DEN-bar mesh-objectivity study (peak
and energy within a few percent across refinement); D — fiber-section backbone convergence;
E — cyclic / low-cycle-fatigue (Coffin–Manson + shakedown); plus F — a profiler perf pass.

> [!note] The load-bearing tests are V3, V4 and V6
> V3 (triaxiality sweep), V4 (3D↔1D reduction) and V6 (mesh-objectivity) are what prove this
> is **real triaxiality-driven fracture with mesh-objective energy**, not a tuning knob.

---

## 11. Limitations and boundaries

> [!caution] Known boundaries
> - **Localization width is not regularized.** `lch` makes the *dissipated energy*
>   mesh-objective, **not** the band width — true width-objectivity needs nonlocal/gradient
>   damage (a future ADR). This is the headline honest limit (§6.1).
> - **Pressure-independent ductile damage only.** It inherits LadrunoJ2's deviatoric-J2
>   boundary — not for concrete/soil/rock, and not a brittle-cracking model.
> - **The 1D mode is degenerate, not multiaxial.** Uniaxial Lemaitre is the
>   triaxiality-$\tfrac13$ energy-based oracle — do not use it as a fracture model; the
>   multiaxial fracture story lives only in 3D (decision D3 caveat).
> - **No element erosion.** Failure caps $D$ at $D_c$ and flags it; failed points are not
>   deleted (§9.2).
> - **Combined hardening + finite strain** remains the spine's §14.11 objectivity boundary
>   (the backstress does not co-rotate) — use isotropic hardening under the log-strain wrapper
>   for large-rotation ductile fracture.
> - **Near-rupture stiffening.** The coupled 2×2 can stiffen as $D\to D_c$ ($Y$ large);
>   mitigations are the operator-split seed, the $D\le D_c$ bracket, and the IMPL-EX escape.
> - **Deferred (separate ADRs):** nonlocal / gradient-enhanced damage; Gurson void-growth
>   (a different, pressure-dependent surface — dSNPO §12.5); anisotropic / tensorial damage;
>   the `LadrunoFatigue` rainflow wrapper (parked, complementary).

---

## 12. References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* — **§12.3**
  (Lemaitre ductile damage coupled to AF-kinematic + isotropic hardening) and **§12.4** (the
  integration algorithm). **The primary framework reference.**
- Lemaitre & Chaboche (1990), *Mechanics of Solid Materials* — continuum damage mechanics,
  the effective-stress / strain-equivalence formulation.
- Lemaitre (1985), *A Continuous Damage Mechanics Model for Ductile Fracture* — the ductile
  damage law $\dot D=(Y/r)^s\dot p$.
- Kobayashi & Ohno (2002) — the scalar reduction of the AF combined-hardening return map
  ($\mathbf n=\mathbf M/\lVert\mathbf M\rVert$) the coupled scheme inherits.

**Within this repo**
- [[15_lemaitre_ductile_damage_adr]] — the source ADR: decisions D1–D4, the full risk
  register, the V0–V9 matrix, the implementation plan.
- [[LadrunoJ2_guide]] — the undamaged 3D spine (§10.2 summarizes this mode).
- [[LadrunoUniaxialJ2_guide]] — the uniaxial twin (the 1D degenerate point).
- [[Ladruno_materials_guide]] — the materials catalog (§6 seeds this mode).
- [[11_brick_asdconcrete_integration]] — the `lch` characteristic-length handshake reused for
  regularization.
- [[LadrunoBrick_reference]] — the continuum element that consumes the damaged material for
  free.
- Source: `SRC/material/nD/LadrunoDamage.h`, `LadrunoJ2.{h,cpp}`,
  `SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}`. **No new classTag** — the mode rides on
  33011 (3D) and 33000 (1D). Validation bundle: `Ladruno_implementation/lemaitre_validation/`.
```