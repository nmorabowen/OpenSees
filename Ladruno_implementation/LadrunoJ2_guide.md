---
title: "LadrunoJ2 — Combined-Hardening von Mises (J2) Plasticity"
project: Ladruno
type: reference / user guide
status: shipped (v1, 3D + all dimensional views, finite-strain lift)
classTag: 33011
material: nDMaterial LadrunoJ2
supersedes:
  - J2Plasticity
  - SimplifiedJ2
related:
  - "[[10_ladruno_j2_plasticity]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[15_lemaitre_ductile_damage_adr]]"
  - "[[13_ladruno_uniaxial_j2_adr]]"
tags:
  - material
  - plasticity
  - j2
  - vonMises
  - reference
---

# LadrunoJ2 — Combined-Hardening von Mises (J2) Plasticity

> [!abstract] One-line summary
> `LadrunoJ2` is a single, self-contained, **rate-independent von Mises (J2)**
> `nDMaterial` for OpenSees that unifies **nonlinear isotropic** (Voce + linear)
> **and nonlinear kinematic** (Chaboche / Armstrong–Frederick) hardening behind one
> clean class. It is the fork's flagship 3D plastic material — the OpenSees analogue
> of **Abaqus `*PLASTIC, COMBINED`** or **LS-DYNA `*MAT_PLASTIC_KINEMATIC` done
> right** — and the constitutive partner of [[10_ladruno_j2_plasticity|LadrunoBrick]]
> and the [[09_finite_strain_material_wrapper|finite-strain wrapper]].

This document is the **descriptive reference**: the continuum theory, the reasoning
that drove the implementation, the numerical algorithm exactly as coded, the OpenSees
wiring, and the intended use cases. For the chronological design/decision log and the
test matrix see [[10_ladruno_j2_plasticity]].

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why it exists — what drove the implementation|2. Why it exists]]
- [[#3 Continuum theory|3. Continuum theory]]
- [[#4 Numerical formulation — the return map|4. Numerical formulation]]
- [[#5 The consistent algorithmic tangent|5. Consistent tangent]]
- [[#6 OpenSees implementation|6. OpenSees implementation]]
- [[#7 Command syntax and usage|7. Command syntax]]
- [[#8 Dimensional views|8. Dimensional views]]
- [[#9 State recording and parameters|9. State recording]]
- [[#10 Optional extensions|10. Optional extensions]]
- [[#11 Verification and validation|11. Verification]]
- [[#12 Limitations and boundaries|12. Limitations]]
- [[#13 References|13. References]]

---

## 1. Intended use cases

`LadrunoJ2` is built for **metal plasticity under cyclic and reversed loading**,
where the *combination* of isotropic and kinematic hardening matters. The canonical
targets:

| Use case | Why LadrunoJ2 |
|---|---|
| **Cyclic / seismic steel** (continuum coupons, panel zones, plate buckling, weld regions) | Reproduces stable hysteresis loops, the **Bauschinger effect**, and the correct ratcheting rate under reversed straining — the behaviour real steel shows. |
| **Buckling-restrained braces (BRBs), shear links, fuses** | Nonlinear kinematic (Chaboche) hardening captures the smooth elastic-to-plastic transition and cyclic backstress saturation that linear-kinematic models over-predict. |
| **Monotonic ductile metal forming / limit analysis** | Voce + linear isotropic hardening gives a saturating monotonic backbone identical to `J2Plasticity` (and to Abaqus' isotropic form). |
| **Large-strain metal plasticity** | Lifts to finite strain through `nDMaterial LogStrain` (Hencky/log-strain) over LadrunoJ2 — see §10 and [[09_finite_strain_material_wrapper]]. |
| **Ductile fracture / damage studies** | Optional, OFF-by-default **Lemaitre coupled damage** mode (`-damage lemaitre …`) — see [[15_lemaitre_ductile_damage_adr]]. |
| **Implicit *and* explicit dynamics** | The material returns a true consistent tangent for implicit Newton solvers and a valid stress for explicit solvers that ignore the tangent — no compromise either way. |

> [!tip] When **not** to reach for LadrunoJ2
> It is a **pressure-independent** (deviatoric) yield model. For soils, concrete,
> rocks, or any pressure-sensitive material use a Drucker–Prager / Mohr–Coulomb /
> cap / ASDConcrete family instead. For rate-dependent (viscoplastic) metals you
> want a Perzyna or Johnson–Cook law — LadrunoJ2 leaves only a hook for that.

---

## 2. Why it exists — what drove the implementation

Before LadrunoJ2 the fork had **two** J2 materials and **neither was adequate** for
cyclic seismic work. The new class was driven by the need to fix both at once:

### 2.1 `J2Plasticity` (Ed Love) — good isotropic, **zero kinematic**

`SRC/material/nD/J2Plasticity.*` is a clean, correct implementation:

- Voce + linear isotropic hardening, $\sigma_y(\bar\varepsilon^p) = \sigma_\infty + (\sigma_0-\sigma_\infty)e^{-\delta\,\bar\varepsilon^p} + H\,\bar\varepsilon^p$.
- Robust local-Newton **radial return**.
- Correct analytic **consistent tangent**.
- A true `getInitialTangent` (elastic) — important for modal/initial-stiffness analyses.

**But it has no kinematic hardening at all.** That means no Bauschinger effect, no
mean-stress relaxation, no ratcheting — it is useless for cyclic loading. LadrunoJ2 is
designed to **reduce exactly to `J2Plasticity`** when kinematic hardening is switched
off (`-kin 0`), so the isotropic path is provably untouched (test V1).

### 2.2 `SimplifiedJ2` (Gu & Qiu) — has kinematic, but **defective**

`SRC/material/nD/SimplifiedJ2.*` has *linear* isotropic + *linear* kinematic
hardening and a closed-form multiplier, but it carries real defects:

- `getInitialTangent` returns the **current (plastic)** tangent — wrong for modal /
  initial-stiffness analyses.
- `revertToStart` / `revertToLastCommit` are **empty stubs** — state cannot be rolled
  back, breaking any adaptive step-cutting algorithm.
- Several `exit(-1)` **kernel-killers** — a failed Newton kills the whole process
  (and, in a Jupyter kernel, the session) instead of returning an error code.
- `setParameter` / `updateParameter` / `Print` are stubs.
- **Linear kinematic only** ⇒ over-predicts cyclic hardening, with no backstress
  saturation, no ratcheting saturation, no mean-stress relaxation.

### 2.3 The missing capability

The capability neither material had is the **combination**: *nonlinear* isotropic
**and** *nonlinear* kinematic hardening, done to implicit (analytic-tangent) grade.
That combination is exactly what reproduces stable hysteresis loops and the correct
ratcheting rate under reversed cyclic loading. LadrunoJ2 delivers it in one class,
with all the hygiene fixes baked in (real `getInitialTangent`, real `revert*`, **no
`exit()` anywhere**, working `setParameter`, real `Print`).

> [!note] Design philosophy
> *One class, one verified return-map kernel.* The inner stress-update lives in a
> header-only, OpenSees-free kernel (`LadrunoJ2Kernel.h`) so the **same** verified map
> serves the small-strain material, every dimensional view, and the finite-strain
> path — there is never a second implementation to drift out of sync.

---

## 3. Continuum theory

> [!info] Notation
> Direct (bold), index, and Voigt forms are shown side by side. Voigt ordering is
> the OpenSees/`J2ThreeDimensional` convention $\{11,22,33,12,23,31\}$. Internally
> the code stores symmetric tensors as **6 true-tensor components**
> $\{t_{00},t_{11},t_{22},t_{01},t_{12},t_{02}\}$ (off-diagonals are **true** tensor
> components, not engineering), and double contractions carry the factor-2 weight on
> the shear slots:
> $$\mathbf a : \mathbf b = a_0 b_0 + a_1 b_1 + a_2 b_2 + 2(a_3 b_3 + a_4 b_4 + a_5 b_5).$$
> The element-facing strain vector uses **engineering** shear $\gamma_{ij}=2\varepsilon_{ij}$; the conversion happens at the dimensional-view boundary (§8).

### 3.1 Additive split and elastic law

Small-strain rate-independent plasticity assumes an **additive decomposition** of the
total strain into elastic and plastic parts:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p,
\qquad \varepsilon_{ij} = \varepsilon^e_{ij} + \varepsilon^p_{ij}.$$

The elastic response is **isotropic linear** with a **volumetric–deviatoric split**.
Decompose stress into hydrostatic pressure $p$ and deviator $\mathbf s$, and strain
into volumetric and deviatoric parts $\mathbf e = \operatorname{dev}\boldsymbol\varepsilon$:

$$p = K\,\operatorname{tr}\boldsymbol\varepsilon, \qquad
\mathbf s = 2G\,(\mathbf e - \boldsymbol\varepsilon^p), \qquad
\boldsymbol\sigma = p\,\mathbf 1 + \mathbf s.$$

Here $K$ is the bulk modulus and $G$ the shear modulus. Plastic flow in J2 is
**isochoric** (volume-preserving), so $\boldsymbol\varepsilon^p$ is deviatoric and the
pressure is purely elastic — this is *why* the volumetric and deviatoric responses
decouple, and why the return map only ever touches the deviator.

### 3.2 The von Mises (J2) yield surface with backstress

Plasticity is **kinematic** when the yield surface *translates* in stress space. We
track that translation with a **backstress** (kinematic hardening) tensor
$\boldsymbol\alpha$ and define the **relative (shifted) deviatoric stress**:

$$\boldsymbol\xi = \mathbf s - \boldsymbol\alpha,
\qquad \xi_{ij} = s_{ij} - \alpha_{ij}.$$

The von Mises yield function measures the relative stress against the current yield
stress:

$$\boxed{\,f(\boldsymbol\xi,\bar\varepsilon^p) = \lVert\boldsymbol\xi\rVert - \sqrt{\tfrac23}\,\sigma_y(\bar\varepsilon^p)\,}$$

- **Direct:** $f = \lVert\boldsymbol\xi\rVert - \sqrt{2/3}\,\sigma_y$
- **Index:** $f = \sqrt{\xi_{ij}\xi_{ij}} - \sqrt{2/3}\,\sigma_y$
- **Voigt:** $f = \lVert\boldsymbol\xi\rVert - \sqrt{2/3}\,\sigma_y$, with the norm
  taken under the factor-2 shear weights above.

The $\sqrt{2/3}$ factor is the standard von Mises normalisation that makes
$\sigma_y$ equal to the **uniaxial** yield stress (so $\lVert\mathbf s\rVert = \sqrt{2/3}\,\sigma_y$ at uniaxial yield). $f<0$ is elastic, $f=0$ is the yield surface,
and $f>0$ is inadmissible.

### 3.3 Associative flow rule and consistency

The plastic strain evolves by an **associative** flow rule (normality — the plastic
strain rate is normal to the yield surface, which follows from the principle of
maximum plastic dissipation):

$$\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\mathbf n,
\qquad \mathbf n = \frac{\boldsymbol\xi}{\lVert\boldsymbol\xi\rVert},$$

where $\dot\gamma \ge 0$ is the **plastic multiplier rate** and $\mathbf n$ is the
unit flow direction. The accumulated **equivalent plastic strain** evolves as

$$\dot{\bar\varepsilon}^p = \sqrt{\tfrac23}\,\dot\gamma$$

(the $\sqrt{2/3}$ again ties the tensor measure to the uniaxial measure). Flow is
governed by the **Karush–Kuhn–Tucker** loading/unloading conditions and the
consistency condition:

$$\dot\gamma \ge 0,\qquad f \le 0,\qquad \dot\gamma\,f = 0,\qquad \dot\gamma\,\dot f = 0.$$

### 3.4 Isotropic hardening — Voce saturation + linear

The yield stress grows with accumulated plastic strain $\bar\varepsilon^p$. LadrunoJ2
uses the **Voce (exponential saturation) + linear** law (the Abaqus `*PLASTIC` form),
implemented in the shared header `LadrunoHardening.h`:

$$\boxed{\,\sigma_y(\bar\varepsilon^p) = \sigma_0 + Q_\infty\bigl(1 - e^{-b\,\bar\varepsilon^p}\bigr) + H_\text{iso}\,\bar\varepsilon^p\,}$$

with the isotropic plastic modulus

$$\frac{\mathrm d\sigma_y}{\mathrm d\bar\varepsilon^p} = Q_\infty\, b\, e^{-b\,\bar\varepsilon^p} + H_\text{iso}.$$

| Parameter | Meaning |
|---|---|
| $\sigma_0$ | initial yield stress |
| $Q_\infty$ | saturation increment (the surface grows by $Q_\infty$ as $\bar\varepsilon^p\to\infty$) |
| $b$ | saturation rate (how fast the exponential approaches $Q_\infty$) |
| $H_\text{iso}$ | linear isotropic modulus (unbounded linear growth) |

Special cases: $Q_\infty = H_\text{iso} = 0$ → **perfectly plastic**;
$Q_\infty = 0$ → **purely linear** isotropic; the law is a **superset** of the
`J2Plasticity` backbone (with the mapping $\sigma_\infty = \sigma_0 + Q_\infty$,
$\delta = b$, $H = H_\text{iso}$).

> Isotropic hardening alone produces a yield surface that **expands** symmetrically.
> It captures cyclic hardening/softening of the *envelope* but **cannot** reproduce
> the Bauschinger effect (early reverse yielding) — that needs kinematic hardening.

### 3.5 Kinematic hardening — Chaboche / Armstrong–Frederick

The backstress is a **superposition of $N$ Armstrong–Frederick (AF) terms** (the
Chaboche model):

$$\boldsymbol\alpha = \sum_{k=1}^{N}\boldsymbol\alpha_k,
\qquad
\dot{\boldsymbol\alpha}_k = \tfrac23 C_k\,\dot{\boldsymbol\varepsilon}^p - \gamma_k\,\boldsymbol\alpha_k\,\dot{\bar\varepsilon}^p.$$

Each AF term has two parameters:

- $C_k$ — the **initial kinematic modulus** (hardening rate of that term),
- $\gamma_k$ — the **recall / dynamic-recovery constant**.

The first term, $\tfrac23 C_k\dot{\boldsymbol\varepsilon}^p$, is the **hardening**
(Prager) term; the second, $-\gamma_k\boldsymbol\alpha_k\dot{\bar\varepsilon}^p$, is
the **dynamic-recovery (recall)** term that bends the backstress toward a saturation
value. Under **monotonic uniaxial** loading each term saturates to

$$\alpha_{k,\text{sat}} = \frac{2}{3}\frac{C_k}{\gamma_k}
\quad\Longleftrightarrow\quad
\sigma_{\text{back},\text{sat}} = \tfrac32\alpha_{k,\text{sat}} = \frac{C_k}{\gamma_k},$$

so the uniaxial saturated stress is $\sigma_0 + \sum_k C_k/\gamma_k$.

> [!important] The recall term is the whole point
> The $-\gamma_k\boldsymbol\alpha_k\dot{\bar\varepsilon}^p$ term is what gives
> **dynamic recovery, ratcheting, and mean-stress relaxation** — the cyclic realism
> seismic steel needs. Set $\gamma_k = 0$ and the term degenerates to **linear
> Prager** kinematic hardening ($\equiv$ `SimplifiedJ2`), which over-predicts cyclic
> hardening because the backstress grows without bound.

**Parameter limits recover every simpler model:**

| Configuration | Recovers |
|---|---|
| `-kin 0` | pure isotropic ⇒ **exactly `J2Plasticity`** |
| `-kin 1 C 0` (i.e. $\gamma_1=0$) | **linear Prager** kinematic ($\equiv$ `SimplifiedJ2`) |
| `-kin 1 C γ` | **single Armstrong–Frederick** |
| `-kin N …` ($N\ge2$) | **full Chaboche** superposition (ship default $N=3$) |

### 3.6 Combined hardening and the Bauschinger effect

With both mechanisms active the yield surface simultaneously **expands** (isotropic)
and **translates** (kinematic). Physically:

- On the **forward** branch the surface grows ($\sigma_y\uparrow$) and its centre
  shifts in the loading direction ($\boldsymbol\alpha$ grows).
- On **reversal** the shifted centre means the relative stress $\boldsymbol\xi = \mathbf s - \boldsymbol\alpha$ reaches the (smaller, reverse-side) surface **earlier**
  — yielding begins **before** the nominal $\sigma_y$ is reached. This is the
  **Bauschinger effect**, and it is exactly what a pure isotropic model misses.

> [!example] Monotonic ↔ reversal split (a key verification)
> Under *monotonic* loading a **linear kinematic** modulus $C$ is mathematically
> identical to a **linear isotropic** modulus $H_\text{iso}=C$ (dSNPO eq. 6.222) —
> both just steepen the backbone. The two only **diverge on reversal**: the kinematic
> model yields earlier in reverse (narrower hysteresis loop). LadrunoJ2's tests pin
> both halves: `test_linear_kinematic_equals_isotropic_monotonic` (they coincide
> forward) and `test_bauschinger_kinematic_vs_isotropic` (they diverge on reversal).
> This is the cleanest proof the $\tfrac23 C$ numerator scaling is correct.

---

## 4. Numerical formulation — the return map

The constitutive update is **implicit backward-Euler closest-point projection (CPP)**.
Given the committed state at step $n$ and a trial total strain at $n{+}1$, find the
new stress, plastic strain, backstresses, and consistent tangent. The classical CPP
for combined AF hardening looks like a coupled tensor system — but it **collapses to a
single scalar Newton iteration** on the plastic multiplier increment $\Delta\gamma$.

### 4.1 Elastic predictor (trial state)

Freeze plastic flow and compute the trial deviatoric stress and relative stress:

$$\mathbf s^{\text{tr}} = 2G\bigl(\mathbf e - \boldsymbol\varepsilon^p_n\bigr),
\qquad
\boldsymbol\xi^{\text{tr}} = \mathbf s^{\text{tr}} - \sum_k\boldsymbol\alpha_{k,n}.$$

Check the trial yield function

$$f^{\text{tr}} = \lVert\boldsymbol\xi^{\text{tr}}\rVert - \sqrt{\tfrac23}\,\sigma_y(\bar\varepsilon^p_n).$$

If $f^{\text{tr}} \le \text{tol}$ → **elastic step**: $\boldsymbol\sigma = K\operatorname{tr}\boldsymbol\varepsilon\,\mathbf 1 + \mathbf s^{\text{tr}}$, tangent $=$ elastic, done. (In code this is the *elastic fast-path*; tolerance is
stress-scaled, $\text{tol}_Y = 10^{-12}\max(\lVert\boldsymbol\xi^{\text{tr}}\rVert,\sqrt{2/3}\,\sigma_y)$.)

### 4.2 The plastic corrector reduced to one scalar equation

Backward-Euler updates for the deviatoric stress and **each** backstress term:

$$\mathbf s_{n+1} = \mathbf s^{\text{tr}} - 2G\,\Delta\gamma\,\mathbf n,
\qquad
\boldsymbol\alpha_{k,n+1} = \frac{\boldsymbol\alpha_{k,n} + \tfrac23 C_k\,\Delta\gamma\,\mathbf n}{1 + \sqrt{\tfrac23}\,\gamma_k\,\Delta\gamma}.$$

The per-term denominator $D_k = 1 + \sqrt{2/3}\,\gamma_k\Delta\gamma$ is the
backward-Euler discretisation of the recall term (using $\dot{\bar\varepsilon}^p = \sqrt{2/3}\,\dot\gamma$). Substituting both updates into the relative stress
$\boldsymbol\xi_{n+1} = \mathbf s_{n+1} - \sum_k\boldsymbol\alpha_{k,n+1}$ gives the
**Kobayashi–Ohno (2002)** reduction:

$$\boldsymbol\xi_{n+1} = \mathbf M(\Delta\gamma) - \theta(\Delta\gamma)\,\mathbf n,$$

$$\mathbf M(\Delta\gamma) = \mathbf s^{\text{tr}} - \sum_k\frac{\boldsymbol\alpha_{k,n}}{D_k}
\quad(\text{known tensor}),
\qquad
\theta(\Delta\gamma) = 2G\,\Delta\gamma + \sum_k\frac{\tfrac23 C_k\,\Delta\gamma}{D_k}
\quad(\text{known scalar}).$$

Because $\mathbf n \parallel \boldsymbol\xi_{n+1}$ and $\theta>0$, the flow direction
is forced to be **collinear with $\mathbf M$**:

$$\boxed{\;\mathbf n = \frac{\mathbf M(\Delta\gamma)}{\lVert\mathbf M(\Delta\gamma)\rVert},
\qquad
\lVert\boldsymbol\xi_{n+1}\rVert = \lVert\mathbf M(\Delta\gamma)\rVert - \theta(\Delta\gamma).\;}$$

> [!note] Radiality is **not** required (it's the special case)
> The classic *radial return* assumes $\boldsymbol\xi_{n+1} \parallel \boldsymbol\xi^{\text{tr}}$. That is just the $\gamma_k = 0$ special case of
> $\mathbf n = \mathbf M/\lVert\mathbf M\rVert$. For nonlinear AF the direction
> **shifts** with $\Delta\gamma$ (each $\boldsymbol\alpha_{k,n}$ is weighted by
> $1/D_k$), but it stays a **closed-form direction** — so this is a genuine *scalar*
> Newton, **not** a coupled tensor solve. (An earlier "≤7-dim coupled" framing was
> overly pessimistic and was superseded.)

### 4.3 The scalar residual and Newton iteration

The consistency condition $f=0$ becomes a single scalar equation in $\Delta\gamma$:

$$\boxed{\;R(\Delta\gamma) = \bigl(\lVert\mathbf M(\Delta\gamma)\rVert - \theta(\Delta\gamma)\bigr) - \sqrt{\tfrac23}\,\sigma_y\bigl(\bar\varepsilon^p_n + \sqrt{\tfrac23}\,\Delta\gamma\bigr) = 0\;}$$

solved by Newton with the analytic derivative

$$\frac{\mathrm dR}{\mathrm d\Delta\gamma} = \frac{\mathbf M : \mathbf M'}{\lVert\mathbf M\rVert} - \theta' - \tfrac23\,\sigma_y',$$

$$\mathbf M'(\Delta\gamma) = \sum_k\boldsymbol\alpha_{k,n}\,\frac{\sqrt{2/3}\,\gamma_k}{D_k^2},
\qquad
\theta'(\Delta\gamma) = 2G + \sum_k\frac{\tfrac23 C_k}{D_k^2}.$$

**Algorithm as coded** (`LadrunoJ2Kernel.h::returnMap`):

```text
dG = 0
for it in 0 .. 49:                      # maxIter = 50
    Dk[k]  = 1 + sqrt(2/3)*gam[k]*dG
    M[i]   = s_tr[i] - Σ_k alpha_n[k][i]/Dk[k]
    Mp[i]  =           Σ_k alpha_n[k][i]*sqrt(2/3)*gam[k]/Dk[k]^2
    normM  = sqrt(dotT(M,M))
    theta  = 2G*dG + Σ_k (2/3)*C[k]*dG/Dk[k]
    dtheta = 2G     + Σ_k (2/3)*C[k]/Dk[k]^2
    pbar   = ebarP_n + sqrt(2/3)*dG
    R      = (normM - theta) - sqrt(2/3)*sig_y(pbar)
    dR     = dotT(M,Mp)/normM - dtheta - (2/3)*sig_y'(pbar)
    if |R| <= tolY: break               # converged
    if dR == 0:     STATUS_SINGULAR; break
    dG -= R/dR
    if dG < 0: dG = 0                    # keep admissible
```

Convergence is on the **stress-scaled** residual $\lvert R\rvert \le \text{tol}_Y$.
Status codes: `STATUS_OK`, `STATUS_SINGULAR` ($\mathrm dR = 0$),
`STATUS_NO_CONVERGE` (hit 50 iterations).

### 4.4 State update at the converged multiplier

With converged $\Delta\gamma$ and $\mathbf n = \mathbf M/\lVert\mathbf M\rVert$:

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\mathbf n,
\qquad
\bar\varepsilon^p_{n+1} = \bar\varepsilon^p_n + \sqrt{\tfrac23}\,\Delta\gamma,$$

$$\boldsymbol\alpha_{k,n+1} = \frac{\boldsymbol\alpha_{k,n} + \tfrac23 C_k\,\Delta\gamma\,\mathbf n}{1 + \sqrt{2/3}\,\gamma_k\,\Delta\gamma},
\qquad
\boldsymbol\sigma_{n+1} = K\operatorname{tr}\boldsymbol\varepsilon\,\mathbf 1 + \mathbf s^{\text{tr}} - 2G\,\Delta\gamma\,\mathbf n.$$

### 4.5 Robustness guards

> [!warning] $\lVert\mathbf M\rVert \to 0$ fallback
> Near a sharp reversal $\mathbf M(\Delta\gamma)$ can shrink toward zero, making the
> flow direction $\mathbf n = \mathbf M/\lVert\mathbf M\rVert$ ill-defined. The
> kernel detects $\lVert\mathbf M\rVert \le \text{normFloor}$
> ($= 10^{-10}\,\text{stressScale}$) and **falls back to the elastic predictor at the
> committed state** rather than dividing by zero. The Newton derivative likewise
> guards $\lVert\mathbf M\rVert > \text{normFloor}$ before forming
> $\mathbf M:\mathbf M'/\lVert\mathbf M\rVert$. This guard was restored from
> `J2Plasticity` during adversarial review.

> [!warning] No `exit()` — ever
> Unlike `SimplifiedJ2`, a non-convergent or singular solve returns a **status code**
> that propagates an OpenSees warning (`opserr`) and a step-cut request, never a libc
> `exit(-1)`. This is a hard fork rule.

---

## 5. The consistent algorithmic tangent

For implicit (Newton) global solvers the material must return the **algorithmic**
(consistent) tangent $\mathbf D^{\text{alg}} = \partial\boldsymbol\sigma/\partial\boldsymbol\varepsilon$ — the *exact* linearisation of the
*discrete* return map, not the continuum tangent. Using the wrong tangent destroys the
quadratic convergence of global Newton. The structure follows dSNPO eq. 7.213,
extended with an Armstrong–Frederick cross-term.

After the state update, with converged $\mathbf n$, $\mathbf M'$,
$\lVert\mathbf M\rVert$, and $\Delta\gamma$, define:

$$h = \theta' + \tfrac23\,\sigma_y' - (\mathbf n : \mathbf M')
\qquad(\text{the hardening denominator, } = -\partial f/\partial\Delta\gamma > 0\ \text{for non-softening params}),$$

$$\mathbf M_\perp = \mathbf M' - (\mathbf n : \mathbf M')\,\mathbf n
\qquad(\text{the part of } \mathbf M' \text{ orthogonal to } \mathbf n).$$

The tangent is assembled from four building blocks:

$$\boxed{\;\mathbf D^{\text{alg}} = K\,(\mathbf 1\otimes\mathbf 1) + 2G\,\beta_1\,\mathbb I_{\text{dev}} + \beta_{NN}\,(\mathbf n\otimes\mathbf n) + \beta_{M_\perp N}\,(\mathbf M_\perp\otimes\mathbf n)\;}$$

with coefficients (writing $G^2 = G\cdot G$):

$$\beta_1 = 1 - \frac{2G\,\Delta\gamma}{\lVert\mathbf M\rVert},
\qquad
\beta_{NN} = \frac{4G^2\,\Delta\gamma}{\lVert\mathbf M\rVert} - \frac{4G^2}{h},
\qquad
\beta_{M_\perp N} = -\frac{4G^2\,\Delta\gamma}{h\,\lVert\mathbf M\rVert}.$$

- $K\,\mathbf 1\otimes\mathbf 1$ — the elastic **volumetric** block (plastic flow is
  isochoric, so the volumetric tangent is always elastic).
- $2G\,\beta_1\,\mathbb I_{\text{dev}}$ — the softened deviatoric block; $\beta_1$
  shrinks the deviatoric stiffness because the radial return relaxes the deviator.
- $\beta_{NN}\,\mathbf n\otimes\mathbf n$ — the rank-one **return-direction** term
  carrying the hardening through $h$.
- $\beta_{M_\perp N}\,\mathbf M_\perp\otimes\mathbf n$ — the **Armstrong–Frederick
  cross-term**: it is non-zero only when the flow direction *shifts* with
  $\Delta\gamma$ (i.e. when $\gamma_k\ne0$ and loading is **non-proportional**). In
  proportional/uniaxial paths $\mathbf M_\perp\approx\mathbf 0$ and the tangent
  reduces to the classical radial-return form.

> [!check] Exact reductions (the regression proofs)
> - **$N=0$ (pure isotropic)** → the tangent collapses **term-by-term to
>   `J2Plasticity`**: $\beta_{M_\perp N}\to0$, and $2G\beta_1$, $\beta_{NN}$ equal
>   Ed Love's $2G+c_3$, $c_2-c_3$ coefficients. The numerical match is $\sim10^{-7}$
>   (bounded by `J2Plasticity`'s internal $\gamma^\star=(1-10^{-8})$ fudge — **not**
>   $10^{-12}$).
> - The 6×6 matrix is assembled in the **engineering-shear `J2ThreeDimensional`
>   convention** (deviatoric projector $\mathbb I_{\text{dev}}$ has $0.5$ on the
>   shear-diagonal), so it maps engineering-shear strain increments to true-stress
>   increments — matching how OpenSees elements expect the material tangent.

The denominator $h>0$ holds for **non-softening** parameters. A user-supplied
$H_\text{iso}<0$ or $Q_\infty<0$ can drive $h\to0$ (an inherited limitation — see
[[#12 Limitations and boundaries|§12]] and `LEDGER_quirks`).

---

## 6. OpenSees implementation

### 6.1 Class and files

| Item | Value |
|---|---|
| Class | `LadrunoJ2 : public NDMaterial` (self-contained; does **not** inherit `J2Plasticity`) |
| classTag | **`ND_TAG_LadrunoJ2 = 33011`** (`SRC/classTags.h`) |
| Material source | `SRC/material/nD/LadrunoJ2.{h,cpp}` |
| Return-map kernel | `SRC/material/nD/LadrunoJ2Kernel.h` (header-only, OpenSees-free) |
| Shared hardening law | `SRC/material/nD/LadrunoHardening.h` (`yieldStressVoceLinear`, `yieldSlopeVoceLinear`) |
| Damage helpers | `SRC/material/nD/LadrunoDamage.h` (Lemaitre — optional) |
| Parser | `OPS_LadrunoJ2` in `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` |
| Broker | `FEM_ObjectBrokerAllClasses.cpp` (`case ND_TAG_LadrunoJ2`) |

> [!note] Why a header-only kernel?
> `LadrunoJ2Kernel.h` holds the return map + tangent as a **pure function** over plain
> `double[]` (`<math.h>` only — no OpenSees types). `LadrunoJ2::integrate()` is a thin
> *pack → call → warn* wrapper. This lets the **same verified map** drive the
> small-strain material **and** the finite-strain path (`nDMaterial LogStrain` over
> LadrunoJ2) with no second implementation to drift. It mirrors the
> `LogStrainKernel.h` / `EnergyBalanceKernel.h` pattern.

### 6.2 The state-commit cycle

OpenSees materials are driven by the element through the standard cycle. LadrunoJ2
holds **committed** history (subscript `_n`) separate from **trial** state:

| Committed (`*_n`) | Trial |
|---|---|
| `epsP_n[6]` plastic strain | `strain6[6]`, `epsP[6]`, `stress6[6]`, `Dtan[6][6]` |
| `ebarP_n` equiv. plastic strain | `ebarP` |
| `alpha_n[MAXBACK][6]` backstresses ($N\le 8$) | `alpha[MAXBACK][6]` |
| `dGamma_n` committed multiplier (IMPL-EX hook) | `dGammaTrial` |
| `D_n` Lemaitre damage | `Dtrial` |
| `cEps22` out-of-plane strain (condensed views) | — |

- `setTrialStrain` → maps the element strain into the 6-component tensor and calls
  `integrate()` (which calls the kernel).
- `commitState` → copies all trial → committed (locks in the step).
- `revertToLastCommit` → copies committed → trial, resets `dGammaTrial=0` (rolls back
  a rejected step — a **real** implementation, unlike `SimplifiedJ2`'s stub).
- `revertToStart` → zeros everything, rebuilds the elastic tangent.
- `getInitialTangent` → the **elastic** tangent $K\,\mathbf 1\otimes\mathbf 1 + 2G\,\mathbb I_{\text{dev}}$ (correct for modal / initial-stiffness analyses — the bug
  `SimplifiedJ2` got wrong).

### 6.3 Serialization (`sendSelf` / `recvSelf`)

State is serialized into a single `Vector` of size **91** for parallel/database use:

```
tag(1) | K,G,sig0,Qinf,bIso,Hiso,rho (7) | nBack(1) | dim(1)
      | damage: dmgOn,r,s,pD,Dc,regOn,lchRef (7)
      | Ckin[0..7] (8) | gKin[0..7] (8)
      | epsP_n[0..5] (6) | ebarP_n (1) | alpha_n[8][6] (48)
      | dGamma_n (1) | D_n (1) | cEps22 (1)            = 91
```

`recvSelf` reads the same layout, then re-runs `setupDim()`, `setStateToCommitted()`,
and rebuilds the elastic tangent.

---

## 7. Command syntax and usage

### 7.1 Grammar

```tcl
nDMaterial LadrunoJ2 $tag  $K $G  <options...>
```

The two positional doubles after the tag are the **bulk** and **shear** moduli
$K,G$. Options may appear in any order:

| Option | Arguments | Meaning |
|---|---|---|
| `-iso voce` | `$sig0 $Qinf $b $Hiso` | Voce + linear isotropic law (the only law in v1) |
| `-kin` | `$N  $C1 $g1  $C2 $g2 …` | $N$ Chaboche AF terms ($0\le N\le 8$), each a $(C_k,\gamma_k)$ pair |
| `-rho` | `$rho` | mass density (default 0) |
| `-damage lemaitre` | `$r $s $pD $Dc` | optional Lemaitre ductile damage (default OFF) — see [[15_lemaitre_ductile_damage_adr]] |
| `-autoRegularization` / `-lchRef` | `$lchRef` | crack-band length regularization for the damage mode |

**Defaults** (everything omitted) → perfectly plastic, no kinematic, no damage:
`sig0=Qinf=b=Hiso=0`, `nBack=0`, `rho=0`, damage off.

> The parser **always constructs a 3D (`DIM_3D`) material**. Dimensional views are
> chosen later when an element requests a typed copy via `getCopy(type)` — see §8.

### 7.2 Examples

```tcl
# --- Pure isotropic (≡ J2Plasticity): Voce + linear, no kinematic ---
nDMaterial LadrunoJ2 1  $K $G  -iso voce 350.0 120.0 8.0 500.0  -kin 0

# --- Single Armstrong–Frederick (Bauschinger, saturating backstress) ---
#     uniaxial saturated stress = sig0 + C/gamma = 350 + 20000/100 = 550
nDMaterial LadrunoJ2 2  $K $G  -iso voce 350.0 0.0 0.0 0.0  -kin 1  20000.0 100.0

# --- Full 3-term Chaboche + Voce isotropic (the recommended cyclic-steel setup) ---
nDMaterial LadrunoJ2 3  $K $G  -iso voce 350.0 60.0 10.0 0.0 \
                               -kin 3  60000.0 500.0  12000.0 120.0  2000.0 10.0 \
                               -rho 7.85e-9

# --- Linear Prager kinematic (≡ SimplifiedJ2, but with correct hygiene) ---
nDMaterial LadrunoJ2 4  $K $G  -iso voce 350.0 0.0 0.0 0.0  -kin 1  2000.0 0.0
```

```python
# OpenSeesPy — three-term Chaboche, consumed by a brick element
import openseespy.opensees as ops
E, nu = 200000.0, 0.3
K, G = E/(3*(1-2*nu)), E/(2*(1+nu))
ops.nDMaterial("LadrunoJ2", 1, K, G,
               "-iso", "voce", 350.0, 60.0, 10.0, 0.0,
               "-kin", 3, 60000.0, 500.0, 12000.0, 120.0, 2000.0, 10.0,
               "-rho", 7.85e-9)
ops.element("stdBrick", 1, *node_tags, 1)   # or LadrunoBrick
```

> [!tip] Calibrating Chaboche parameters
> A common starting point: pick the saturated cyclic stress $\sigma_\text{sat}$ and
> split $\sum_k C_k/\gamma_k = \sigma_\text{sat} - \sigma_0$ across 3 terms with
> decreasing $\gamma_k$ — a *fast* term (large $\gamma$, captures the sharp
> elastic-plastic knee), a *medium* term, and a near-linear *slow* term (small
> $\gamma$) for the long-range hardening / ratcheting tail.

---

## 8. Dimensional views

A **single class** drives every element type via a `dim` mode. The 3D return map
always runs on the full 6-component tensor; the strain/stress/tangent are mapped to
the element's reduced ordering through an index table `vmap[]`. The element selects
the view by calling `getCopy(type)`:

| `getType()` string | `dim` | `ncomp` | `vmap` (reduced → full) | $\sigma_{22}=0$? |
|---|---|---|---|---|
| `ThreeDimensional` / `3D` | `DIM_3D` | 6 | `{0,1,2,3,4,5}` | no |
| `PlaneStrain` | `DIM_PSTRAIN` | 3 | `{0,1,3}` | no |
| `AxiSymmetric` | `DIM_AXISYM` | 4 | `{0,1,2,3}` | no |
| `PlateFiber` | `DIM_PLATEFIBER` | 5 | `{0,1,3,4,5}` | **yes** |
| `PlaneStress` | `DIM_PSTRESS` | 3 | `{0,1,3}` | **yes** |

At the boundary the code converts **engineering shear → tensor shear** on input
($\gamma_{ij}\to\tfrac12\gamma_{ij}$ for the off-diagonal slots) and back on output.

### Plane-stress / plate-fiber condensation

For `PlaneStress` and `PlateFiber` the out-of-plane normal stress must vanish,
$\sigma_{22}=0$. LadrunoJ2 enforces this with a **nested Newton on
$\varepsilon_{22}$** (dSNPO §9.2.3 route):

```text
for it in 0 .. 24:
    integrate()                       # full 3D return map
    if |stress6[2]| <= tol22: break
    strain6[2] -= stress6[2] / Dtan[2][2]
```

then **statically condenses** the 33-dof out of the 6×6 tangent so the reported
reduced tangent is consistent with $\sigma_{22}=0$:

$$D_{IJ} \;\leftarrow\; D_{IJ} - \frac{D_{I2}\,D_{2J}}{D_{22}}.$$

The committed $\varepsilon_{22}$ is carried in `cEps22` so the nested solve restarts
from the converged guess each step. These reduced paths are validated against
`J2Plasticity` on a `quad` element for both `PlaneStrain` and `PlaneStress`
(`test_dimensional_view_reduce_to_J2Plasticity`).

---

## 9. State recording and parameters

### 9.1 Recordable responses

Through the element's `material` response (e.g.
`eleResponse(ele, "material", gp, "<name>")` or a `recorder … material <gp> <name>`):

| Response name(s) | ID | Returns |
|---|---|---|
| `stress` / `stresses` | 1 | stress vector (reduced view, true-stress shear) |
| `strain` / `strains` | 2 | strain vector (engineering shear) |
| `tangent` | 3 | algorithmic tangent (`ncomp×ncomp`) |
| `backStress` / `alpha` | 4 | **total** backstress $\sum_k\boldsymbol\alpha_k$ (reduced view) |
| `plasticStrain` | 5 | plastic strain $\boldsymbol\varepsilon^p$ (engineering shear) |
| `equivalentPlasticStrain` / `ebarP` | 6 | scalar $\bar\varepsilon^p$ |
| `damage` / `D` | 7 | scalar Lemaitre damage $D$ |
| `failed` / `rupture` | 8 | $1$ if (damage on and $D\ge D_c$), else $0$ |
| `Y` / `energyReleaseRate` | 9 | Lemaitre energy release rate (effective) |

```tcl
# record the accumulated equivalent plastic strain at Gauss point 1
recorder Element -ele 1 -file ebar.out material 1 equivalentPlasticStrain
```

### 9.2 Settable parameters (for FE sensitivity / updates)

`setParameter` / `updateParameter` expose the scalar elastic & isotropic parameters:
`K`, `G` (`mu`), `rho`, `sigmaY` (`sig0`), `Hiso`, `Qinf`, `b`. (The kinematic
$C_k,\gamma_k$ and damage parameters are **not** exposed as settable parameters in
v1.)

---

## 10. Optional extensions

### 10.1 Finite-strain lift (shipped)

LadrunoJ2 presents the **exact `J2ThreeDimensional` inner contract** the Hencky
log-strain wrapper expects (engineering-shear strain in / true stress out, linear
elasticity so $\mathbf C^e$ is constant, stateful $\boldsymbol\varepsilon^p$
subtraction, constant `getInitialTangent`). So **finite-strain combined-hardening J2
is simply**:

```tcl
nDMaterial LadrunoJ2 1 $K $G -iso voce ... -kin 3 ...
nDMaterial LogStrain 2 1        ;# Hencky/log-strain over LadrunoJ2
element LadrunoBrick ... 2 -geom finite
```

with **zero** change to `LogStrainNDMaterial`. See
[[09_finite_strain_material_wrapper]].

> [!warning] The §14.11 objectivity boundary
> The lift is **exact for the isotropic spine** (the elastic state co-rotates and
> isotropic yield depends only on $\lVert\mathbf s\rVert,\bar\varepsilon^p$, so a
> superposed rotation rotates the Cauchy stress rigidly — verified to ~$10^{-9}$).
> But **combined hardening is not exactly objective under large rotation**: the
> backstress $\boldsymbol\alpha$ is stored in the inner's *fixed* frame and does not
> co-rotate, so $\lVert\mathbf s-\boldsymbol\alpha\rVert$ is not rotation-invariant.
> This is a **framework** limit (dSNPO defers it to §14.11), pinned as a strict
> `xfail`. v1 ⇒ exact for no/small rotation; a finite-strain-native v2 that
> co-rotates $\boldsymbol\alpha$ each step is the planned fix (the kernel extraction
> is precisely the enabler).

### 10.2 Lemaitre coupled ductile damage (optional, OFF by default)

Adding `-damage lemaitre $r $s $pD $Dc` activates a **strain-equivalence (effective
stress)** ductile-damage mode (dSNPO §12.3/§12.4). Key property: the return map runs
on the **undamaged effective stress** (which carries no damage), so the state update
**decouples exactly** — with damage OFF the behaviour is **byte-identical** to the
undamaged material (no new classTag, no cost). All the $(\Delta\gamma, D)$ coupling
lives in the consistent tangent, $\mathbf D = (1-D)\mathbf D^{\text{alg}} - \tilde{\boldsymbol\sigma}\otimes\partial D/\partial\boldsymbol\varepsilon$. Full theory
and validation: [[15_lemaitre_ductile_damage_adr]].

### 10.3 IMPL-EX hook (structure-only)

The committed multiplier $\Delta\gamma_n$ is stored (`dGamma_n`) so an IMPL-EX
(implicit–explicit) path — extrapolated $\tilde{\Delta\gamma}_{n+1}$ with a constant
SPD tangent — can be added later for explicit dynamics or softening/damage robustness.
**No IMPL-EX code path ships in v1**; only the state layout is prepared.

### 10.4 Deferred backlog

Tabulated and Bézier/Bernstein isotropic curves; the `prager_nl` nonlinear-Prager
oracle mode (dSNPO Box 7.5); a native scalar `LadrunoUniaxialJ2`
([[13_ladruno_uniaxial_j2_adr]]); the plane-stress-projected route (dSNPO §9.4);
the finite-strain-native backstress-co-rotation v2.

---

## 11. Verification and validation

Zone-A pytest battery `tests/test_ladrunoJ2_material.py` (single `stdBrick` /
`quad`, displacement- and load-controlled):

| Test | What it pins |
|---|---|
| `test_elastic_uniaxial` | $\sigma_{xx} = E\,\varepsilon_{xx}$ below yield |
| `test_reduce_to_J2Plasticity` | $N=0$ reproduces upstream `J2Plasticity` (~$10^{-7}$) |
| `test_linear_kinematic_equals_isotropic_monotonic` | linear $C\equiv H_\text{iso}$ monotonically (pins the $\tfrac23 C$ scaling) |
| `test_bauschinger_kinematic_vs_isotropic` | the two **diverge on reversal**; kinematic loop is narrower |
| `test_af_monotonic_saturation` | single-AF monotonic stress → $\sigma_0 + C/\gamma$ |
| `test_state_recording` | `backStress` is deviatoric, axial → $\tfrac23(C/\gamma)$; `equivalentPlasticStrain` accumulates |
| `test_reduce_to_J2Plasticity_3D_mixed_shear` | full 3D-with-shear matches `J2Plasticity` (exercises `dotT` factor-2, `IIDEV6` shear block, $\mathbf n\otimes\mathbf n$ shear slots) |
| `test_dimensional_view_reduce_to_J2Plasticity` | `PlaneStrain` + `PlaneStress` `quad` match `J2Plasticity` (validates `vmap` + condensation) |

Additional kernel-level proofs (no OpenSees build, g++ + numpy): an independent
3×3-tensor oracle (`ladrunoj2_reference.py`) matches stress, full internal state
($\boldsymbol\varepsilon^p,\bar\varepsilon^p,\boldsymbol\alpha,\Delta\gamma$) and the
algorithmic tangent (vs independent finite-difference) over isotropic, single-AF,
three-term-Chaboche-with-shear, and linear-kinematic paths; plus closed-form
linear-kinematic and non-proportional (axial→shear) cross-term legs added during a
33-agent adversarial review that found **no correctness bugs**.

---

## 12. Limitations and boundaries

> [!caution] Known boundaries
> - **Pressure-independent** (deviatoric J2 only) — not for soils/concrete/rock.
> - **Rate-independent** — no viscoplastic rate effects (hook only).
> - **No ratcheting threshold** (Ohno–Wang) — single-AF over-predicts uniaxial
>   ratcheting; multi-term Chaboche mitigates but does not eliminate it. The
>   over-prediction signature is a documented, intended behaviour.
> - **Finite-strain combined hardening** is not exactly objective under large
>   rotation (§10.1, the §14.11 boundary). Exact for the isotropic spine.
> - **Softening tangent** $h\to0$: a user-supplied $H_\text{iso}<0$ or $Q_\infty<0$
>   can drive the tangent denominator $h = \theta' + \tfrac23\sigma_y' - \mathbf n:\mathbf M'$ to zero (inherited from the radial-return structure; logged in
>   `LEDGER_quirks`).
> - The `N=0` reduction matches `J2Plasticity` to ~$10^{-7}$, **not** $10^{-12}$,
>   bounded by `J2Plasticity`'s internal $\gamma^\star=(1-10^{-8})$ fudge.

---

## 13. References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* —
  Ch. 6 (mathematical theory of plasticity), §7.2 (general return map), §7.4.4 &
  eq. 7.213 (general consistent tangent), §9.2.3/§9.4 (plane stress), §12.3/§12.4
  (Lemaitre damage), §14.11 (finite-strain kinematic hardening). **The primary
  framework reference.**
- Belytschko, Liu, Moran & Elkhodary (2014), *Nonlinear Finite Elements* — Ch. 5
  (constitutive models, stress-update algorithms).
- Armstrong & Frederick (1966); Chaboche (1986, 1991) — the nonlinear kinematic
  hardening law.
- Kobayashi & Ohno (2002) — the scalar reduction of the AF combined-hardening return
  map ($\mathbf n = \mathbf M/\lVert\mathbf M\rVert$).
- Abaqus Theory Manual, `*PLASTIC, COMBINED` (Mises) — the functional equivalent.

**Within this repo**
- [[10_ladruno_j2_plasticity]] — design log, decisions, full test matrix.
- [[09_finite_strain_material_wrapper]] — the log-strain finite lift.
- [[15_lemaitre_ductile_damage_adr]] — the optional damage mode.
- [[13_ladruno_uniaxial_j2_adr]] — the planned native scalar sibling.
- Source: `SRC/material/nD/LadrunoJ2.{h,cpp}`, `LadrunoJ2Kernel.h`,
  `LadrunoHardening.h`. classTag **33011**.
```
