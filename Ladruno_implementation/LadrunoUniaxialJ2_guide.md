---
title: "LadrunoUniaxialJ2 — Uniaxial Combined-Hardening J2 Steel Material"
project: Ladruno
type: reference / user guide
status: shipped (v1 core)
classTag: 33000
material: uniaxialMaterial LadrunoUniaxialJ2
twin_of: "[[LadrunoJ2_guide|nDMaterial LadrunoJ2]]"
related:
  - "[[13_ladruno_uniaxial_j2_adr]]"
  - "[[LadrunoJ2_guide]]"
  - "[[10_ladruno_j2_plasticity]]"
  - "[[15_lemaitre_ductile_damage_adr]]"
tags:
  - material
  - uniaxial
  - plasticity
  - j2
  - steel
  - chaboche
  - reference
---

# LadrunoUniaxialJ2 — Uniaxial Combined-Hardening J2 Steel Material

> [!abstract] One-line summary
> `LadrunoUniaxialJ2` is the **uniaxial twin** of the 3D
> [[LadrunoJ2_guide|nDMaterial LadrunoJ2]]: a rate-independent 1-D von Mises (J2)
> `UniaxialMaterial` with **Voce + linear isotropic** and **Chaboche / Armstrong–
> Frederick kinematic** hardening, consumable by **fiber sections, trusses, and
> zeroLength** elements. It is deliberately **not** a "better Steel02" — its three
> jobs are (1) the **verification oracle** that pins the 3D Chaboche calibration to
> ~$10^{-12}$ under uniaxial stress, (2) a fiber material with **true
> multi-backstress ratcheting** (the one thing Menegotto–Pinto cannot do), and
> (3) the clean **core** onto which future steel features (fracture/LCF,
> local-buckling, weld/HAZ) bolt as separable layers.

This document is the **descriptive reference**: the theory, the reasoning that drove
the design, the algorithm exactly as coded, and the OpenSees wiring. For the full
decision record and roadmap see [[13_ladruno_uniaxial_j2_adr]]; for the 3D parent
see [[LadrunoJ2_guide]].

---

## Contents

- [[#1 Intended use cases — and what it deliberately is *not*|1. Intended use cases]]
- [[#2 Why it exists — what drove the design|2. Why it exists]]
- [[#3 Theory — the exact uniaxial reduction of 3D J2|3. Theory]]
- [[#4 Numerical formulation — the 1-D return map|4. Numerical formulation]]
- [[#5 The consistent algorithmic tangent|5. Consistent tangent]]
- [[#6 The oracle contract — shared hardening code|6. The oracle contract]]
- [[#7 OpenSees implementation|7. OpenSees implementation]]
- [[#8 Command syntax and usage|8. Command syntax]]
- [[#9 State recording and parameters|9. State recording]]
- [[#10 Optional extensions and roadmap|10. Extensions]]
- [[#11 Verification and validation|11. Verification]]
- [[#12 Limitations and boundaries|12. Limitations]]
- [[#13 References|13. References]]

---

## 1. Intended use cases — and what it deliberately is *not*

`LadrunoUniaxialJ2` is a **fiber/truss/spring** material for **uniaxial metal
plasticity with combined hardening**. It plugs into any OpenSees
`uniaxialMaterial` slot.

| Use case | Why LadrunoUniaxialJ2 |
|---|---|
| **Verification oracle** for the 3D `LadrunoJ2` Chaboche calibration | A uniaxial coupon has an *analytic* answer; the 3D return map does not. Driving the 3D material to a uniaxial stress state must reproduce this material's closed-form $\sigma(\varepsilon)$ to ~$10^{-12}$ with **identical** $(C_k,\gamma_k)$. This is the cheapest high-value regression the fork owns. |
| **Ratcheting / mean-stress relaxation** in a fiber | Multi-backstress AF with distinct recovery rates $\gamma_k$ is *the* classical ratcheting model. No OpenSees uniaxial material reproduces progressive strain accumulation under asymmetric cycling **mechanistically**. |
| **Cyclic metal fibers** (pressure vessels, piping, rails, rolling contact, asymmetric mechanical demand) | The features that survive a uniaxial reduction of J2 — ratcheting, mean-stress relaxation — govern exactly these problems. |
| **The substrate for future steel features** | Fracture/LCF index, local-buckling-aware degradation, and weld/HAZ zones are designed to bolt on as **separable layers** (wrappers or flags) without polluting the oracle core. The v1 core already exposes the `plasticStrain` / `plasticWork` / `backStress` hooks they need. |

> [!warning] What this material is **not**
> It is **not** a replacement for `Steel02` (Giuffré–Menegotto–Pinto, tag 26) or
> `Steel4` (tag 87). Earthquake demand on structural steel is broadly **symmetric
> low-cycle** — Bauschinger + low-cycle fatigue dominate, and those incumbents
> already cover it, calibrated to death. A standalone "uniaxial J2 beats Steel02"
> claim does **not** hold and is explicitly not the goal. Reach for `Steel02` /
> `Steel4` for ordinary seismic steel fibers; reach for `LadrunoUniaxialJ2` when you
> need the oracle property or genuine multi-AF ratcheting.

### How it relates to the incumbents

| material | model | cyclic mechanism | ratcheting |
|---|---|---|---|
| `Steel02` (tag 26) | Giuffré–Menegotto–Pinto | smooth $R(R_0,c_{R1},c_{R2})$ asymptote transition → **Bauschinger by geometry** | **no** |
| `Steel4` (tag 87) | M–P + separated iso/kin, asym, ult, fracture | iso + nonlinear-kin (M–P-style) | weak/approx |
| **`LadrunoUniaxialJ2`** | uniaxial von Mises + Voce iso + **Chaboche AF ($N$ backstresses)** | **mechanistic** backstress evolution | **yes (multi-AF)** |

---

## 2. Why it exists — what drove the design

Two facts drove the design, recorded in [[13_ladruno_uniaxial_j2_adr]]:

1. **The 3D model needs an analytic oracle.** Under a uniaxial stress state the 3D
   `LadrunoJ2` must reproduce a closed-form 1-D law *using the identical
   $(C_k,\gamma_k)$ and Voce parameters*. That cross-check pins the 3D Chaboche
   calibration and the $\tfrac23$ bookkeeping far more sharply than any 3D self-test
   — because the uniaxial coupon has an analytic answer and a 3D return map does not.
   The only way the cross-check is *real* is if the **hardening physics is literally
   the same code** (§6).

2. **No OpenSees uniaxial material ratchets mechanistically.** Menegotto–Pinto
   (`Steel02`) produces Bauschinger by curve geometry but cannot accumulate strain
   under asymmetric cycling. Multi-backstress AF can.

A third, architectural driver: the genuinely novel steel features on the roadmap
(fracture/LCF, local-buckling, weld/HAZ) all want a **clean, undamaged,
energy-consistent core** to compose over. So v1 ships the pure core, with the hooks
those layers need but none of their math.

And, as with the 3D class, all the `SimplifiedJ2` anti-patterns are corrected: a
**real `getInitialTangent()` $= E$** (elastic — correct for modal/eigen analysis),
real `revertToLastCommit` / `revertToStart`, real `setParameter` / `Print`, and
**no `exit()` anywhere** (a non-convergent solve warns and proceeds, never kills the
process).

---

## 3. Theory — the exact uniaxial reduction of 3D J2

The whole point of the theory is to show that the 3-D von Mises model of
[[LadrunoJ2_guide]] **collapses exactly** to the classical 1-D combined-hardening
model — which is *why* the oracle test can hit $10^{-12}$ and *why* the constants
$(C_k,\gamma_k,\sigma_y)$ carry over unchanged.

### 3.1 From 3D von Mises to 1D

Take a uniaxial stress state $\boldsymbol\sigma = \sigma\,(\mathbf e_1\otimes\mathbf e_1)$ (only $\sigma_{11}\ne0$). Its deviator is

- **Direct:** $\mathbf s = \operatorname{dev}\boldsymbol\sigma = \sigma\,(\mathbf e_1\otimes\mathbf e_1 - \tfrac13\mathbf 1)$
- **Index:** $s_{ij} = \sigma\,(\delta_{i1}\delta_{j1} - \tfrac13\delta_{ij})$
- **Voigt $\{11,22,33,12,23,13\}$:** $\{s\} = \sigma\,\{\tfrac23,-\tfrac13,-\tfrac13,0,0,0\}^{\mathsf T}$

so $\lVert\mathbf s\rVert = \sqrt{\mathbf s:\mathbf s} = |\sigma|\sqrt{\tfrac23}$.
In a uniaxial *process* the backstress inherits the same structure,
$\boldsymbol\alpha = X\,(\mathbf e_1\otimes\mathbf e_1 - \tfrac13\mathbf 1)$, with $X$
the **axial backstress** (units of stress). The relative stress is then collinear:

$$\boldsymbol\xi = \mathbf s - \boldsymbol\alpha = (\sigma - X)\,(\mathbf e_1\otimes\mathbf e_1 - \tfrac13\mathbf 1)
\quad\Longrightarrow\quad
\lVert\boldsymbol\xi\rVert = |\sigma - X|\sqrt{\tfrac23}.$$

Substituting into the 3-D yield $f = \lVert\boldsymbol\xi\rVert - \sqrt{\tfrac23}\,\sigma_y$:

$$f = \sqrt{\tfrac23}\bigl(|\sigma - X| - \sigma_y(\bar\varepsilon^p)\bigr) = 0
\quad\Longleftrightarrow\quad
\boxed{\,|\sigma - X| = \sigma_y(\bar\varepsilon^p)\,}$$

> [!important] The $\sqrt{2/3}$ **divides out**
> It appears on both sides of the 3-D yield equation and cancels. **The 1-D return
> map contains no $\sqrt{2/3}$ factor anywhere** — confirmed in the code. The
> uniaxial Young's modulus is recovered from the bulk/shear pair the 3-D class
> stores, $E = \dfrac{9KG}{3K+G}$. This is exactly *why* the same $(C_k,\gamma_k)$
> reproduce the 3-D response: the reduction is exact, not a re-calibration.

### 3.2 Governing equations (1-D)

| | Scalar form |
|---|---|
| elastic law | $\sigma = E\,(\varepsilon - \varepsilon^p)$ |
| yield | $f(\sigma,X,\bar\varepsilon^p) = \lvert\sigma - X\rvert - \sigma_y(\bar\varepsilon^p)$ |
| isotropic law | $\sigma_y(\bar\varepsilon^p) = \sigma_0 + Q_\infty(1-e^{-b\,\bar\varepsilon^p}) + H_\text{iso}\,\bar\varepsilon^p$ |
| flow | $\dot\varepsilon^p = \dot\gamma\,s,\quad s = \operatorname{sign}(\sigma - X)$ |
| accumulated | $\dot{\bar\varepsilon}^p = \lvert\dot\varepsilon^p\rvert = \dot\gamma$ |
| backstress (Chaboche AF) | $\dot X = \sum_k \dot X_k,\quad \dot X_k = C_k\,\dot\varepsilon^p - \gamma_k X_k\,\dot{\bar\varepsilon}^p = (C_k s - \gamma_k X_k)\,\dot\gamma$ |
| KKT | $\dot\gamma\ge0,\quad f\le0,\quad \dot\gamma f = 0$ |

Each AF term **self-saturates** to $X_{k,\text{sat}} = C_k/\gamma_k$ (monotonic) and
reverses on unloading through the $-\gamma_k X_k$ **recall** term — the dynamic
recovery that yields the Bauschinger effect, stable hysteresis, and ratcheting.

| Configuration | Recovers |
|---|---|
| `-kin 0` | pure isotropic (no Bauschinger) ⇒ uniaxial `J2Plasticity`/`Hardening` backbone |
| `-kin 1 C 0` ($\gamma_1=0$) | linear **Prager** kinematic ($\equiv$ `SimplifiedJ2`) |
| `-kin 1 C γ` | single **Armstrong–Frederick** |
| `-kin N …` ($N\ge2$) | full **Chaboche** (ship $N=3$, arch for $N\le8$) |

> [!note] $\sigma_y$ is **byte-identical** to the 3-D class
> Both materials evaluate the yield stress through the shared header
> `LadrunoHardening.h` (`yieldStressVoceLinear`, `yieldSlopeVoceLinear`). This is
> the "oracle contract" (§6): the isotropic backbone is the *same code* in 1-D and
> 3-D, so the $10^{-12}$ cross-check is meaningful.

### 3.3 The Bauschinger split

Under **monotonic** loading a linear kinematic modulus $C$ is identical to a linear
isotropic modulus $H_\text{iso}=C$ — both just steepen the backbone. The two only
**diverge on reversal**: the kinematic model yields earlier in reverse (narrower
hysteresis loop), because the shifted yield centre $X$ means $|\sigma - X|$ reaches
the reverse-side surface before $\sigma$ reaches the nominal $\sigma_y$. This is the
**Bauschinger effect** — captured by kinematic hardening, missed by isotropic. The
tests pin both halves (`test_linear_kinematic_equals_isotropic_monotonic` and
`test_bauschinger_kinematic_vs_isotropic`).

---

## 4. Numerical formulation — the 1-D return map

Implicit **backward-Euler** elastic-predictor / plastic-corrector. The crucial
simplification over the 3-D kernel: **in 1-D the flow direction is fixed by the trial
state** — the relative stress cannot rotate on a line — so there is **no $\mathbf M(\Delta\gamma)$ direction iteration**. Consistency reduces to **one scalar Newton on
$\Delta p$** (the accumulated-plastic-strain increment).

### 4.1 Elastic predictor

Freeze the committed plastic state and form the trial stress and relative stress:

$$\sigma^{\text{tr}} = E\,(\varepsilon_{n+1} - \varepsilon^p_n),
\qquad
\eta^{\text{tr}} = \sigma^{\text{tr}} - X_n,\quad X_n = \sum_k X_{k,n},
\qquad
f^{\text{tr}} = \lvert\eta^{\text{tr}}\rvert - \sigma_y(\bar\varepsilon^p_n).$$

The tolerance is **unit-aware** and never collapses to a fixed absolute:
$\text{stressScale} = \max(\lvert\eta^{\text{tr}}\rvert,\sigma_y)$,
$\text{tol}_Y = 10^{-12}\,\max(\text{stressScale},10^{-300})$.

If $f^{\text{tr}} \le \text{tol}_Y$ → **elastic**: $\sigma_{n+1} = \sigma^{\text{tr}}$,
$E_{\text{alg}} = E$, history unchanged. (If damage is on and $D_n>0$ the prior
damage degrades the stress/tangent but no *new* damage accrues on an elastic step.)

### 4.2 Plastic corrector — fixed direction, scalar Newton

The flow direction is **frozen** at $s = \operatorname{sign}(\eta^{\text{tr}})$. The
backward-Euler updates are

$$\sigma_{n+1} = \sigma^{\text{tr}} - E\,\Delta p\,s,
\qquad
X_{k,n+1} = \frac{X_{k,n} + C_k\,\Delta p\,s}{1 + \gamma_k\,\Delta p}.$$

Folding the per-term AF update into the consistency condition $|\sigma_{n+1} - X_{n+1}| = \sigma_y(\bar\varepsilon^p_n + \Delta p)$ gives a **single scalar
residual** in $\Delta p$:

$$\boxed{\;R(\Delta p) = \lvert\eta^{\text{tr}}\rvert - E\,\Delta p - \sigma_y(\bar\varepsilon^p_n + \Delta p) - \sum_k \frac{\Delta p\,a_k}{1 + \gamma_k\,\Delta p} = 0,
\qquad a_k = C_k - \gamma_k\,s\,X_{k,n}\;}$$

with the analytic derivative (note $a_k$ uses the **committed** $X_{k,n}$, so it is
constant in $\Delta p$):

$$\frac{\mathrm dR}{\mathrm d\Delta p} = -E - \sigma_y'(\bar\varepsilon^p_n + \Delta p) - \sum_k \frac{a_k}{(1 + \gamma_k\,\Delta p)^2}.$$

**Algorithm as coded** (`LadrunoUniaxialJ2::setTrialStrain`):

```text
s   = sign(eta_tr)
h0  = sig_y'(ebarP_n) + Σ_k C[k]          # seed hardening modulus
dp  = f_tr / (E + h0);  if dp < 0: dp = 0  # perfectly-plastic-style seed
for it in 0 .. 49:                         # maxIter = 50
    R  = |eta_tr| - E*dp - sig_y(ebarP_n + dp)
    dR = -E - sig_y'(ebarP_n + dp)
    for k:
        a_k = C[k] - gam[k]*s*X_n[k]
        den = 1 + gam[k]*dp
        R  -= dp*a_k/den
        dR -= a_k/den^2
    if |R| <= tolY: break                  # converged
    if dR == 0: warn "singular"; break      # NO exit
    dp -= R/dR
    if dp < 0: dp = 0
```

A non-convergent solve warns (`tag`, $|R|$) and **proceeds with the last $\Delta p$**
— never `exit()`.

### 4.3 State update

$$\bar\varepsilon^p_{n+1} = \bar\varepsilon^p_n + \Delta p,
\qquad
\varepsilon^p_{n+1} = \varepsilon^p_n + \Delta p\,s,
\qquad
X_{k,n+1} = \frac{X_{k,n} + C_k\,\Delta p\,s}{1 + \gamma_k\,\Delta p},$$

$$\sigma_{n+1} = \sigma^{\text{tr}} - E\,\Delta p\,s,
\qquad
W^p_{n+1} = W^p_n + \sigma_{n+1}\,(\Delta p\,s)\quad(\text{accumulated plastic work}).$$

The plastic work $W^p$ is tracked specifically as a hook for the future fracture/LCF
layer (§10).

---

## 5. The consistent algorithmic tangent

Differentiating the converged consistency condition gives the classical 1-D
elastoplastic tangent:

$$\boxed{\;E_{\text{alg}} = \frac{E\,h}{E + h},
\qquad
h = \sigma_y'(\bar\varepsilon^p_{n+1}) + \sum_k \frac{a_k}{(1 + \gamma_k\,\Delta p)^2}\;}$$

where $h$ is the **total plastic hardening modulus** (isotropic slope $+$ every
backstress's algorithmic kinematic modulus $H^{\text{kin}}_k = a_k/(1+\gamma_k\Delta p)^2$, with the same $a_k = C_k - \gamma_k s X_{k,n}$ as the Newton
residual). The isotropic slope is $\sigma_y'(\bar\varepsilon^p) = Q_\infty b\,e^{-b\bar\varepsilon^p} + H_\text{iso}$.

- **Elastic step** → $E_{\text{alg}} = E$ exactly.
- **Deep plastic** → $0 < E_{\text{alg}} < E$.
- **$\gamma_k\to0$ (linear)** → reduces to linear hardening $E_{\text{alg}} = E h/(E+h)$ with $h = H_\text{iso} + \sum_k C_k$.

> [!check] Why this is *cleaner* than the 3-D tangent and `Steel02`
> Because the direction $s$ is constant, there is **no** $\mathbf M_\perp\otimes\mathbf n$ cross-term (the 3-D AF complication vanishes) and **no** implicit
> $R$-curvature evaluation (`Steel02`). The result is the textbook scalar
> $E h/(E+h)$. The consistent-tangent FD test (`test_consistent_tangent_fd`)
> verifies it against a central difference of the one-step return at several
> plastic states, and confirms $E_{\text{alg}}=E$ elastically and $0<E_{\text{alg}}<E$ deep in plasticity.

---

## 6. The oracle contract — shared hardening code

The oracle property is only real if the hardening physics is *literally the same
code* as the 3-D class. That is achieved through a tiny header-only file:

- **`SRC/material/nD/LadrunoHardening.h`** — `Ladruno::yieldStressVoceLinear(p,σ0,Q∞,b,Hiso)` and `Ladruno::yieldSlopeVoceLinear(…)`. Consumed
  **verbatim** by both `LadrunoJ2` (3-D) and `LadrunoUniaxialJ2` (1-D).

The **return-map structure intentionally differs** and stays separate:

| | 3-D `LadrunoJ2` | 1-D `LadrunoUniaxialJ2` |
|---|---|---|
| flow direction | $\mathbf n = \mathbf M(\Delta\gamma)/\lVert\mathbf M\rVert$ (shifts with $\Delta\gamma$) | $s = \operatorname{sign}(\eta^{\text{tr}})$ (**fixed**) |
| consistency solve | scalar Newton on $\Delta\gamma$ (Kobayashi–Ohno) | scalar Newton on $\Delta p$ |
| $\sqrt{2/3}$ factors | present (tensor ↔ uniaxial measure) | **absent** (cancel out) |
| backstress denominator | $1 + \sqrt{2/3}\,\gamma_k\Delta\gamma$ | $1 + \gamma_k\Delta p$ |

The AF kinematic update is *provably equivalent* between the two but written in each
one's natural variable — the tensor multiplier $\Delta\gamma$ in 3-D, the accumulated
plastic-strain increment $\Delta p = \sqrt{2/3}\,\Delta\gamma$ in 1-D — so it is
**not** shared (forcing the 1-D law through the 6-component machinery would discard
the closed-form-direction win and *create* drift risk). The shared piece is exactly
and only the isotropic backbone, which is what the oracle test checks.

> [!example] The headline regression (V7)
> `test_oracle_vs_3D_LadrunoJ2` drives a 3-D `nDMaterial LadrunoJ2` (on a
> 1/8-symmetry unit-cube `stdBrick`) and this 1-D material along the **same** axial
> push-pull-cycle path with **identical** $(\sigma_0,Q_\infty,b,H_\text{iso})$ and a
> 2-term Chaboche $(C_1,\gamma_1,C_2,\gamma_2)$, and asserts the axial stresses
> agree to ~$10^{-7}$ (the residual is only the 3-D element's lateral-equilibrium
> iteration). If a future 3-D hardening edit ever breaks this, the shared header has
> been bypassed.

---

## 7. OpenSees implementation

### 7.1 Class and files

| Item | Value |
|---|---|
| Class | `LadrunoUniaxialJ2 : public UniaxialMaterial` |
| classTag | **`MAT_TAG_LadrunoUniaxialJ2 = 33000`** (`SRC/classTags.h`; Ladruno private *uniaxial* band — the per-registry 33000 does **not** collide with ELE/ND/INTEGRATOR 33000) |
| Source | `SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}` |
| Parser | `OPS_LadrunoUniaxialJ2` (in the `.cpp`) |
| Registration | `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`; Tcl via `TclModelBuilderUniaxialMaterialCommand.cpp` |
| Shared hardening | `SRC/material/nD/LadrunoHardening.h` |
| Damage helpers | `SRC/material/nD/LadrunoDamage.h` (optional Lemaitre) |
| Max backstresses | `MAXBACK = 8` |

### 7.2 The state-commit cycle

The material holds **committed** (`C*`) and **trial** (`T*`) history:

| Committed | Meaning |
|---|---|
| `CplasticStrain` | plastic strain $\varepsilon^p$ |
| `Cebarp` | accumulated plastic strain $\bar\varepsilon^p$ |
| `CX[MAXBACK]` | backstresses $X_k$ |
| `CdGamma` | last $\Delta p$ (IMPL-EX hook) |
| `CD` | Lemaitre damage $D$ |
| `Cwp` | accumulated plastic work $W^p$ (fracture/LCF hook) |

(`T*` mirrors these; plus trial state `Tstrain, Tstress, Ttangent`.)

- `setTrialStrain(strain)` → the return map of §4 (with a no-op early-out when the
  strain is unchanged to within `DBL_EPSILON`).
- `commitState` → copies every `T*` → `C*`.
- `revertToLastCommit` → copies `C*` → `T*` (real, not a stub).
- `revertToStart` → zeros all history, sets `Ttangent = E`.
- `getInitialTangent()` → returns **`E`** (elastic — correct for modal/eigen
  analysis; the bug `SimplifiedJ2` got wrong).
- `getCopy()` → a **full-state clone** (deep-copies all committed *and* trial
  history plus trial state — history is preserved, not reset to virgin; this matters
  for committed-clone backbone probes).

### 7.3 Serialization (`sendSelf` / `recvSelf`)

State packs into one `Vector` of size **42** ($= 1 + 5 + 1 + 5 + 2\cdot8 + 4 + 1 + 8 + 1$):

```
tag(1) | E,sig0,Qinf,bIso,Hiso (5) | nBack(1)
      | damage: dmgOn,r,s,pD,Dc (5)
      | Ckin[0..7] (8) | gKin[0..7] (8)
      | CplasticStrain,Cebarp,CdGamma,Cwp (4) | CD (1) | CX[0..7] (8)
      | Tstrain (1)                                       = 42
```

`recvSelf` unpacks the same order, then calls `revertToLastCommit()` (trial ← committed) and recovers `Tstress = E(Tstrain − CplasticStrain)`, `Ttangent = E`.
Trial history is reconstructed from committed, not separately serialized.

---

## 8. Command syntax and usage

### 8.1 Grammar

```tcl
uniaxialMaterial LadrunoUniaxialJ2 $tag $E  <options...>
```

The single positional double after the tag is the **Young's modulus** $E$ (a fiber
material is parameterized in $E$, not $K,G$; the oracle mapping
$E = 9KG/(3K+G)$ is documented for the cross-check, not exposed). Options in any
order:

| Option | Arguments | Meaning |
|---|---|---|
| `-iso voce` | `$sig0 $Qinf $b $Hiso` | Voce + linear isotropic law (only law in v1) |
| `-kin` | `$N  $C1 $g1  $C2 $g2 …` | $N$ Chaboche AF terms ($0\le N\le 8$), each a $(C_k,\gamma_k)$ pair |
| `-damage lemaitre` | `$r $s $pD $Dc` | optional Lemaitre ductile damage (default OFF); `-damage none` disables |

**Defaults** (everything omitted) → perfectly plastic, no kinematic, no damage:
`sig0=Qinf=b=Hiso=0`, `nBack=0`, damage off. Validation: `-iso` accepts only `voce`;
`-kin N` requires $0\le N\le8$ and exactly $N$ $(C,\gamma)$ pairs; `-damage lemaitre`
requires $r>0$ and $0<D_c\le1$. An unknown flag/law aborts parsing with a warning.

### 8.2 Examples

```tcl
# --- Pure isotropic (Voce + linear), no kinematic ---
uniaxialMaterial LadrunoUniaxialJ2 1  200000.0  -iso voce 350.0 120.0 8.0 500.0  -kin 0

# --- Single Armstrong–Frederick: Bauschinger + saturating backstress ---
#     saturated stress = sig0 + C/gamma = 350 + 20000/100 = 550
uniaxialMaterial LadrunoUniaxialJ2 2  200000.0  -iso voce 350.0 0.0 0.0 0.0  -kin 1  20000.0 100.0

# --- Full 3-term Chaboche + Voce iso (recommended cyclic / ratcheting setup) ---
uniaxialMaterial LadrunoUniaxialJ2 3  200000.0 \
        -iso voce 350.0 60.0 10.0 0.0 \
        -kin 3  60000.0 500.0  12000.0 120.0  2000.0 10.0

# --- Linear Prager kinematic (≡ SimplifiedJ2/Hardening, with correct hygiene) ---
uniaxialMaterial LadrunoUniaxialJ2 4  200000.0  -iso voce 350.0 0.0 0.0 0.0  -kin 1  2000.0 0.0
```

```python
# OpenSeesPy — three-term Chaboche in a fiber section
import openseespy.opensees as ops
ops.uniaxialMaterial("LadrunoUniaxialJ2", 1, 200000.0,
                     "-iso", "voce", 350.0, 60.0, 10.0, 0.0,
                     "-kin", 3, 60000.0, 500.0, 12000.0, 120.0, 2000.0, 10.0)
# ... ops.fiber(y, z, A, 1) inside a section, or ops.element("Truss", ..., 1)
```

> [!tip] Calibrating the Chaboche terms
> Pick the saturated cyclic stress $\sigma_\text{sat}$, set $\sum_k C_k/\gamma_k = \sigma_\text{sat}-\sigma_0$, and split it across a *fast* term (large $\gamma$, sharp
> elastic-plastic knee), a *medium* term, and a near-linear *slow* term (small
> $\gamma$) for the long-range hardening / ratcheting tail. For ratcheting studies
> the slow term's $\gamma$ governs the per-cycle accumulation rate.

---

## 9. State recording and parameters

### 9.1 Recordable responses

Through the element's `material` response (e.g.
`eleResponse(ele, "material", "<name>")` or `recorder Element … material <name>`):

| Response name(s) | ID | Returns |
|---|---|---|
| `backStress` / `alpha` | 100 | total backstress $\sum_k X_k$ (stress units) |
| `plasticStrain` | 101 | plastic strain $\varepsilon^p$ (signed) |
| `equivalentPlasticStrain` / `ebarP` | 102 | $\bar\varepsilon^p$ |
| `plasticWork` / `plasticEnergy` | 103 | accumulated plastic work $W^p$ |
| `damage` / `D` | 104 | Lemaitre damage $D$ |
| `failed` / `rupture` | 105 | $1$ if (damage on and $D\ge D_c$), else $0$ |
| `Y` / `energyReleaseRate` | 106 | Lemaitre energy release rate $\tilde\sigma^2/(2E)$ |

Stress / strain / tangent go through the `UniaxialMaterial` base class as usual.

```tcl
# accumulated equivalent plastic strain at a truss/fiber material
recorder Element -ele 1 -file ebar.out material equivalentPlasticStrain
recorder Element -ele 1 -file wp.out   material plasticWork
```

### 9.2 Settable parameters

`setParameter` / `updateParameter` expose `E`, `sigmaY` (`sig0`/`fy`/`Fy`), `Qinf`,
`b`, `Hiso`. (Kinematic $C_k,\gamma_k$ and damage parameters are not settable in v1 —
but the dispatch is deliberately shaped so a future buckling-length / slenderness key
can be added without reshaping it.)

---

## 10. Optional extensions and roadmap

### 10.1 Lemaitre ductile damage (optional, OFF by default)

`-damage lemaitre $r $s $pD $Dc` activates a **strain-equivalence / effective-stress**
ductile-damage mode — the exact 1-D reduction of the 3-D `returnMapDamaged`. In 1-D
the **triaxiality is fixed at $1/3$**, so the triaxiality factor $R_v\equiv1$ and the
energy release rate is simply $Y = \tilde\sigma^2/(2E)$ ($\tilde\sigma$ = effective /
undamaged stress). The return map runs on the undamaged effective stress (which
carries no damage), so with damage **OFF the behaviour is bit-identical** to the
undamaged material (no new classTag, no cost). When on, damage couples into the
stress via $\sigma = (1-D)\tilde\sigma$ and into the tangent via the extra
$-\tilde\sigma\,\partial D/\partial\varepsilon$ term. Full theory and validation:
[[15_lemaitre_ductile_damage_adr]].

### 10.2 IMPL-EX hook (structure-only)

The committed $\Delta p_n$ is stored (`CdGamma`) so an implicit–explicit path can be
added later. No IMPL-EX code path ships in v1 — only the state layout is prepared.

### 10.3 Future steel layers (separate ADRs — hooks only in v1)

These are **designed-for but not built**, kept as separable layers so the oracle
core stays pure:

1. **Fracture / low-cycle-fatigue index** (`LadrunoFatigue`-style **wrapper**) —
   accumulate damage $D\in[0,1]$ from plastic-strain half-cycles (Coffin–Manson /
   Miner) or dissipated energy; at $D=1$ rupture. The v1 core already exposes
   `plasticStrain`, `plasticWork`, and `backStress` for a wrapper to read.
2. **Local-buckling-aware degradation** — asymmetric compression-side knockdown
   (Dhakal–Maekawa / `ReinforcingSteel`-style). This must be a **core flag, not a
   wrapper**, because buckling rewrites the compression-side return asymmetrically;
   the buckling length enters as a **parameter** (slenderness $\lambda = L/D$), *not*
   by widening the scalar strain interface.
3. **Weld / HAZ** — largely a parameter-set + the fracture layer (reduced toughness);
   no new constitutive core.

See [[13_ladruno_uniaxial_j2_adr]] §6 for the full roadmap and design rationale.

---

## 11. Verification and validation

Zone-A pytest battery `tests/test_ladrunoUniaxialJ2_material.py` — a unit `Truss`
(L=1, A=1) uniaxial probe, plus a 3-D unit-cube driver for the oracle test:

| Test | What it pins |
|---|---|
| `test_elastic_uniaxial` | $\sigma = E\varepsilon$ below yield |
| `test_reduce_to_HardeningMaterial` | linear iso+kin reproduces upstream `Hardening` over a full reversal (~$10^{-9}$) |
| `test_linear_kinematic_equals_isotropic_monotonic` | linear $C\equiv H_\text{iso}$ monotonically (pins the $C_k$ slope scaling) |
| `test_bauschinger_kinematic_vs_isotropic` | the two **diverge on reversal**; kinematic loop is narrower |
| `test_af_monotonic_saturation` | single-AF monotonic stress → $\sigma_0 + C/\gamma$ |
| **`test_oracle_vs_3D_LadrunoJ2`** | **the headline**: 3-D `LadrunoJ2` driven uniaxially matches this material with identical $(C_k,\gamma_k,\text{iso})$, $N=2$ push-pull-cycle (~$10^{-7}$) |
| `test_state_recording` | `plasticStrain` / `equivalentPlasticStrain` / `backStress` ($\to C/\gamma$) / `plasticWork` recordable through the element |
| `test_consistent_tangent_fd` | $E_{\text{alg}}$ matches central-difference $\mathrm d\sigma/\mathrm d\varepsilon$ at several plastic states; $=E$ elastic, $0<E_{\text{alg}}<E$ plastic |
| `test_ratcheting_stress_controlled` | **the differentiator vs `Steel02`**: asymmetric stress cycles → monotone per-cycle plastic-strain accumulation; single-AF rate is ~constant (the documented over-prediction signature) |

---

## 12. Limitations and boundaries

> [!caution] Known boundaries
> - **Not a `Steel02`/`Steel4` replacement** — for ordinary symmetric seismic steel
>   fibers those calibrated incumbents are the right tool (§1).
> - **Rate-independent** — no viscoplastic rate effects.
> - **Ratcheting fidelity**: single-AF **over-predicts** uniaxial ratcheting (the
>   intended, documented signature); multi-term Chaboche with split $\gamma_k$
>   improves it, but **Ohno–Wang** is the real fix and is not implemented. Do not
>   claim quantitative ratcheting accuracy.
> - **Voce + linear isotropic only** in v1 — tabulated/Bézier isotropic curves are
>   gated on the 3-D class shipping them first (to keep the shared-header contract
>   honest).
> - **No fracture / buckling / weld layers** in v1 — those are separate ADRs; only
>   the read-out hooks (`plasticStrain`, `plasticWork`, `backStress`) exist.
> - **Softening tangent**: a user-supplied $H_\text{iso}<0$ or $Q_\infty<0$ can drive
>   the hardening modulus $h\to -E$ and the tangent $E_{\text{alg}}=Eh/(E+h)$
>   singular (inherited from the elastoplastic-tangent structure).

---

## 13. References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* —
  Ch. 3 §3.5 / Box 3.1 (1-D elastoplasticity), §6.6.4 (Armstrong–Frederick),
  §7 (return mapping). The reduction frame for §3.
- Simo & Hughes (1998), *Computational Inelasticity* — Box 1.5 (1-D combined
  hardening).
- Ibrahimbegovic (2009), *Nonlinear Solid Mechanics* — Ch. 3 (inelastic behaviour).
- Armstrong & Frederick (1966); Chaboche (1986, 1991) — nonlinear kinematic
  hardening.
- Ohno & Wang (1993) — the ratcheting-threshold model (the documented future fix).

**Within this repo**
- [[13_ladruno_uniaxial_j2_adr]] — the full ADR: context, decision, roadmap, test
  matrix.
- [[LadrunoJ2_guide]] — the 3-D parent / oracle target.
- [[10_ladruno_j2_plasticity]] — the 3-D design log.
- [[15_lemaitre_ductile_damage_adr]] — the optional damage mode.
- Source: `SRC/material/uniaxial/LadrunoUniaxialJ2.{h,cpp}`, shared
  `SRC/material/nD/LadrunoHardening.h`. classTag **33000**.
```
