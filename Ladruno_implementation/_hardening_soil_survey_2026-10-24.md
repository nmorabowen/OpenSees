# Deep Research: Hardening Soil (`LadrunoHS` / `HSsmall`) Architecture, Dual `mc|mn` Lode Formulations & Validation Campaign

**Date:** 2026-10-24  
**Work Package:** WP-210 ([[210_ladruno_hs_plan]])  
**Scope:**
1. Comparative analysis of `LadrunoSANISAND` (`SRC/material/nD/LadrunoSANISAND.h`, ADR-86 / WP-138 / WP-151 / WP-152) and `LadrunoNorSand` (`Ladruno_implementation/144_ladruno_norsand_plan.md`, `Ladruno_implementation/_sand_model_survey_2026-09-27.md`) against the **Hardening Soil (`HS` / `HSsmall`)** family.
2. First-principles mathematical and numerical survey of Hardening Soil models across **PLAXIS** (`HS`, `HSsmall`, `GHS`), **FLAC/FLAC3D** (`Plastic-Hardening`, `PHSS`, `CYSoil`), **OptumG2/G3** (`HMC`), **ZSoil** (`HSS`), **numgeo**, and OpenSees's `ASDPlasticMaterial3D` (`StiffSoilShear` / `StiffSoilCap`).
3. **Architectural Decision:** Standalone unified `LadrunoHS` material (`ND_TAG_LadrunoHS = 33029`, with `-smallStrain` / `LadrunoHSSmall` alias for `HSsmall` and **`-lodeShape mc|mn`** for hexagonal Mohr–Coulomb vs. Panteghini–Lagioia convex Matsuoka–Nakai bounding comparisons) backed by an OpenSees-free `LadrunoHSKernel.h`.
4. **4-Tier Validation Campaign:** Complete benchmark suite, parameter tables, analytical oracles, Lode-bounding envelopes, and BVP targets linked to the Zotero **`Soil` (`LDCNC45G`)** and **`TIMs` (`HTBQAHHH`)** libraries.

---

## 1. Where Hardening Soil Sits Relative to `LadrunoSANISAND` and `LadrunoNorSand`

| Model | Lineage | Core State Variables | Primary Domain | Key Numerical / Physical Lesson in Fork |
|---|---|---|---|---|
| **`LadrunoSANISAND`** (`33019`–`33021`, shipped) | Dafalias & Manzari (2004) (`QYI8IQLH`) bounding-surface rotational kinematic hardening + fabric $\mathbf{z}$ | $(\boldsymbol{\sigma}, e, \boldsymbol{\alpha}, \boldsymbol{\alpha}_{\text{in}}, \mathbf{z})$ | Cyclic mobility, liquefaction, site response | Explicit SAS-ME (`IntScheme 129`, `LadrunoSANISANDSasME.cpp`) needed because discrete $\boldsymbol{\alpha}_{\text{in}}$ reseats create $b_0 / (\boldsymbol{\alpha} - \boldsymbol{\alpha}_{\text{in}}):\mathbf{n} \to \infty$ singularities and Zeno chatter under BVP Newton iteration (`151_sanisand_reseat_singularity.md`). Low-$p'$ tension cutoff added in WP-152 (`152_sanisand_tension_cutoff.md`). |
| **`LadrunoNorSand`** (`33023`–`33025`, planned WP-144) | Jefferies (1993) (`7TC285BN`) / Borja & Andrade (2006) (`KVS9MY4L`) | $(\boldsymbol{\varepsilon}^e, e, p_i)$ | Static liquefaction, tailings dams, CPT, dense-to-loose drained/undrained sand | Smooth single bullet-shaped yield surface; hyperelastic free energy $\psi(\varepsilon_v^e, \varepsilon_s^e)$; exact 3-invariant $(p, q, \theta)$ closest-point return map (`144_ladruno_norsand_plan.md`). |
| **`LadrunoHySAND`** (`33026`–`33028`, planned WP-145) | Simonin, Houlsby & Byrne (2026) | $(\boldsymbol{\sigma}, e, \boldsymbol{\alpha}_1 \dots \boldsymbol{\alpha}_N)$ | Cyclic ratcheting, offshore monopiles, thermodynamic hyperplasticity | $N$-surface kinematic hardening removes the discontinuous $\boldsymbol{\alpha}_{\text{in}}$ memory variable at the cost of $6N$ internal back-stress components (`145_ladruno_hysand_plan.md`). |
| **`LadrunoHS` / `HSsmall`** (`33029`, WP-210) | Schanz, Vermeer & Bonnier (1999) (`TU759V8E`) + Benz (2007) (`5PU9KPRT`) + Panteghini & Lagioia (2014) (`C93MGN9R`) | $(\boldsymbol{\sigma}, \gamma^p, p_p, [e], [\mathbf{H}_{6\text{D}}])$ | **Excavations, retaining walls, tunnels, foundations, embankments, SSI serviceability** in sands, silts, and stiff/soft clays | **Isotropic double-hardening** (hyperbolic shear cone $\gamma^p$ + elliptical volumetric cap $p_p$) + **MC/MN failure cutoff** + **Rankine tension cutoff** + **Benz (2007) paraelastic strain-history tensor $\mathbf{H}$** for $G/G_0(\gamma)$. |

---

## 2. Mathematical Formulation of `LadrunoHS` & Dual `-lodeShape mc|mn` Bounding

> [!NOTE]
> **Sign Convention:** Equations below use standard geotechnical **compression-positive** effective stresses $\sigma'_1 \ge \sigma'_2 \ge \sigma'_3$ ($p' = \tfrac{1}{3}I_1 > 0$ in compression). OpenSees C++ vectors (`getStress()`, `getStrain()`) use tension-positive and are sign-flipped at the kernel boundary.

### 2.1 Stress-Dependent Stiffness & Reference Moduli `[E]`
At reference pressure $p^{\text{ref}}$ (typically $100\text{ kPa}$), user inputs are $(E_{50}^{\text{ref}}, E_{\text{oed}}^{\text{ref}}, E_{ur}^{\text{ref}}, m, c', \phi', \psi, \nu_{ur}, K_0^{\text{nc}}, R_f)$. The dimensionless stress-dependency factor $s(\boldsymbol{\sigma}')$ is:

$$s(\boldsymbol{\sigma}') = \begin{cases} \left(\dfrac{c'\cos\phi' + \max(\sigma'_3, \sigma_{\min})\sin\phi'}{c'\cos\phi' + p^{\text{ref}}\sin\phi'}\right)^m & \text{if } \texttt{-stiffStress sigma3} \text{ (default, PLAXIS HS / FLAC PH)} \\[10pt] \left(\dfrac{c'\cos\phi' + \max(p', p_{\min})\sin\phi'}{c'\cos\phi' + p^{\text{ref}}\sin\phi'}\right)^m & \text{if } \texttt{-stiffStress p} \text{ (PLAXIS GHS / ZSoil / numgeo)} \end{cases}$$

$$E_{50} = E_{50}^{\text{ref}}\,s, \qquad E_{ur} = E_{ur}^{\text{ref}}\,s, \qquad E_i = \frac{2 E_{50}}{2 - R_f}, \qquad G_{ur} = \frac{E_{ur}}{2(1+\nu_{ur})}, \qquad K_{ur} = \frac{E_{ur}}{3(1-2\nu_{ur})}$$

*(Admissibility check enforced at construction per Obrzud & Truty 2018 [`D4IVV2XD`]: $E_{ur}^{\text{ref}} > E_i^{\text{ref}} = \frac{2 E_{50}^{\text{ref}}}{2 - R_f}$ so plastic shear strain $\gamma^p$ is strictly positive on primary loading).*

---

### 2.2 Dual Lode-Shape Formulations (`-lodeShape mc|mn`) as Rigorous Bounding Envelopes `[E, R]`

Why support **both** `-lodeShape mc` (Mohr–Coulomb hexagon) and `-lodeShape mn` (Panteghini & Lagioia 2014 [`C93MGN9R`] strictly convex Matsuoka–Nakai)?

1. **Exact Coincidence at Triaxial Compression ($\vartheta = +\pi/6$) and Triaxial Extension ($\vartheta = -\pi/6$) `[E]`:**
   Let $J_2 = \tfrac{1}{2}\mathbf{s}:\mathbf{s}$, $J_3 = \det\mathbf{s}$, and Lode angle $\vartheta = \frac{1}{3}\arcsin\left(-\frac{3\sqrt{3}}{2}\frac{J_3}{J_2^{3/2}}\right) \in [-\frac{\pi}{6}, +\frac{\pi}{6}]$ ($\vartheta = +\frac{\pi}{6}$ in TXC $\sigma'_1 > \sigma'_2 = \sigma'_3$; $\vartheta = -\frac{\pi}{6}$ in TXE $\sigma'_1 = \sigma'_2 > \sigma'_3$).
   - **In `-lodeShape mc` (Schanz et al. 1999 [`TU759V8E`] / PLAXIS `HS` / FLAC `PH`):**
     The shear yield and failure surfaces are written directly on principal stress pair $(1,3)$ (and $(1,2), (2,3)$ at corners):
     $$q_{\text{MC}} = \sigma'_1 - \sigma'_3, \qquad q_{f,\text{MC}}(\sigma'_3) = (c'\cot\phi' + \sigma'_3)\frac{2\sin\phi'}{1 - \sin\phi'}, \qquad q_{a,\text{MC}} = \frac{q_{f,\text{MC}}}{R_f}$$
     Or equivalently in $(p', q, \vartheta)$ invariant space with $q = \sqrt{3 J_2}$:
     $$g_{\text{MC}}(\vartheta, \phi') = \frac{3 - \sin\phi'}{2\sqrt{3}\sin\phi'}\left(\cos\vartheta - \frac{1}{\sqrt{3}}\sin\vartheta\sin\phi'\right), \qquad \Theta_{\text{MC}}(\vartheta, \phi') \equiv \frac{g_{\text{MC}}(\vartheta, \phi')}{g_{\text{MC}}(+\pi/6, \phi')}$$
   - **In `-lodeShape mn` (Panteghini & Lagioia 2014 [`C93MGN9R`] Eqs. 17–19 / `numgeo` [`2JN8QSRB`]):**
     With $m_{\text{MN}} = \frac{9 - \sin^2\phi'}{1 - \sin^2\phi'} = 9 + 8\tan^2\phi'$,
     $$\beta(\phi') = \sqrt{\frac{m_{\text{MN}} - 9}{m_{\text{MN}} - 3}} = \frac{2\sin\phi'}{\sqrt{3 + \sin^2\phi'}}, \qquad \gamma(\phi') = \frac{m_{\text{MN}}\sqrt{m_{\text{MN}} - 9}}{(m_{\text{MN}} - 3)^{3/2}} = \frac{(9 - \sin^2\phi')\sin\phi'}{(3 + \sin^2\phi')^{3/2}}$$
     $$g_{\text{MN}}(\vartheta, \phi') = \frac{2\sqrt{3}}{\beta(\phi')} \cos\left[\frac{1}{3}\arccos\big(-\gamma(\phi')\sin(3\vartheta)\big)\right], \qquad \Theta_{\text{MN}}(\vartheta, \phi') \equiv \frac{g_{\text{MN}}(\vartheta, \phi')}{g_{\text{MN}}(+\pi/6, \phi')}$$
   - **Exact Endpoint Parity `[R]`:**
     At TXC ($\vartheta = +\pi/6$, $\sin 3\vartheta = +1$) and TXE ($\vartheta = -\pi/6$, $\sin 3\vartheta = -1$):
     $$\Theta_{\text{MN}}(+\pi/6, \phi') = \Theta_{\text{MC}}(+\pi/6, \phi') = 1, \qquad \Theta_{\text{MN}}(-\pi/6, \phi') = \Theta_{\text{MC}}(-\pi/6, \phi') = \frac{3 + \sin\phi'}{3 - \sin\phi'} = \delta$$
     Thus, **every triaxial compression (`CIDC`/`CIUC`) and triaxial extension (`CIDE`/`CIUE`) test produces identical stress–strain curves under `-lodeShape mc` and `-lodeShape mn`!**
2. **Strict Bounding Property at Intermediate Lode Angles ($-\pi/6 < \vartheta < +\pi/6$, e.g. `PlaneStrain`) `[E, R]`:**
   - For all $-\pi/6 < \vartheta < +\pi/6$, the convex Matsuoka–Nakai locus strictly circumscribes the Mohr–Coulomb hexagon:
     $$\Theta_{\text{MN}}(\vartheta, \phi') < \Theta_{\text{MC}}(\vartheta, \phi') \implies q_{f,\text{MN}}(p', \vartheta) > q_{f,\text{MC}}(p', \vartheta)$$
   - Specifically, in **2D `PlaneStrain`** ($b = \frac{\sigma'_2 - \sigma'_3}{\sigma'_1 - \sigma'_3} \approx 0.28\text{--}0.38$, $\vartheta \approx +5^\circ\text{ to }+12^\circ$):
     - `-lodeShape mc` maintains $\phi_{\text{PS}}^{\text{MC}} = \phi_{\text{TXC}}$ (**lower bound** on plane-strain resistance, matching PLAXIS `HS` and Loukidis & Salgado 2009 [`JK8ZJP9D`]).
     - `-lodeShape mn` yields $\phi_{\text{PS}}^{\text{MN}} \approx 1.09\text{--}1.12\,\phi_{\text{TXC}}$ (**upper bound** on plane-strain resistance, matching Tatsuoka et al. 1986 [`2XI5UT3Y`] plane-strain compression experiments and explaining why $s_\gamma^{\text{circ}} < 1$ in Loukidis & Salgado 2011 [`9GMR7TZD`] when comparing circular vs. strip footings at the same relative density!).

---

### 2.3 Unified 4-Surface Yield & Potential System (`LadrunoHSKernel.h`) `[E]`

1. **Surface 1 — Hyperbolic Shear Hardening Cone ($f^s$):**
   $$f^s(\boldsymbol{\sigma}', \gamma^p) = \frac{2}{E_i(\sigma'_3)} \frac{q_{\text{eq}}}{1 - q_{\text{eq}} / q_a} - \frac{2 q_{\text{eq}}}{E_{ur}(\sigma'_3)} - \gamma^p = 0$$
   where $q_{\text{eq}} = \sigma'_1 - \sigma'_3$ under `mc` (or $\Theta_{\text{MN}}(\vartheta, \phi')\sqrt{3J_2}$ under `mn`).
   Non-associated plastic potential $g^s$ uses Rowe's stress–dilatancy $\sin\psi_m = \frac{\sin\phi_m - \sin\phi_{cv}}{1 - \sin\phi_m\sin\phi_{cv}}$ (with small-mobilization contraction cutoff and optional $e \ge e_{\max}$ dilatancy cutoff).
2. **Surface 2 — Ultimate Failure Truncation ($f^f$ at $q_{\text{eq}} = q_f$):**
   Truncates $f^s$ when $\gamma^p \ge \gamma^p_f = \frac{2q_f}{E_i(1-R_f)} - \frac{2q_f}{E_{ur}}$ so $q_{\text{eq}}$ never exceeds $q_f = R_f q_a$, transitioning to perfect plasticity with constant dilatancy $\psi$.
3. **Surface 3 — Elliptical Volumetric Cap ($f^c$, Associated Flow $g^c = f^c$):**
   $$f^c(\boldsymbol{\sigma}', p_p) = \frac{\tilde{q}^2}{\alpha^2} + (p' + a)^2 - (p_p + a)^2 = 0, \qquad \dot{p}_p = H_c(p_p)\,\dot{\varepsilon}_v^{p,c}$$
   where $a = c'\cot\phi'$, $\tilde{q} = \sigma'_1 + (\delta - 1)\sigma'_2 - \delta\sigma'_3$ (`mc`) or $\Theta_{\text{MN}}(\vartheta, \phi')\sqrt{3J_2}$ (`mn`), with $(\alpha, \beta)$ solved automatically at construction from $(E_{\text{oed}}^{\text{ref}}, K_0^{\text{nc}})$.
4. **Surface 4 — Rankine Tension Cutoff ($f^t$, Associated Flow $g^t = f^t$):**
   $$f^t_i(\boldsymbol{\sigma}') = -\sigma'_i - \sigma_t = 0 \quad (i = 1, 2, 3)$$
5. **Opt-In `HSsmall` Overlay (`-smallStrain $G0_ref $gamma07`, Benz 2007 [`5PU9KPRT`] / Benz, Vermeer & Schwab 2009 [`2U66T7JX`]):**
   Tracks the 6-component symmetric deviatoric strain-history tensor $\mathbf{H}$ via eigenvector projection in the principal frame of $\Delta\mathbf{e}$ to evaluate $\gamma_{\text{hist}}$ and modified Hardin–Drnevich ($a = 0.385$) tangent modulus:
   $$G_t^{\text{eff}}(\gamma_{\text{hist}}) = \max\left(\frac{G_0(\sigma'_3)}{\big(1 + 0.385\,\gamma_{\text{hist}}/\gamma_{0.7}\big)^2},\; G_{ur}(\sigma'_3)\right)$$

---

## 3. Why Standalone `LadrunoHS` (Not `ASDPlasticMaterial3D`) — Architectural Summary

| Criterion | `ASDPlasticMaterial3D` (`StiffSoilShear` / `StiffSoilCap`) | Standalone `LadrunoHS` (`LadrunoHSKernel.h`) |
|---|---|---|
| **Multi-Surface Coupling** | **Single-surface template** (`<EL, YF, PF>`); forced to split into 2 uncoupled materials (`StiffSoilShear` and `StiffSoilCap` in `ASD_material_definitions.cpp:433-461`). | **Coupled 4-surface Koiter active-set return map** ($f^s + f^f + f^c + f^t$) with simultaneous shear+cap yielding ($\Delta\lambda^s > 0, \Delta\lambda^c > 0$). |
| **Failure Surface & Cap Math** | Missing $f^f$ cutoff at $q_f$ (allows $q \to 0.99 q_f/R_f$ in `StiffSoilShear_YF.h:106`); conflates $\alpha$ with $\delta$ in `StiffSoilCap_YF.h:99`. | Exact $f^f$ truncation at $q = q_f$; distinct Lode $\delta$ and cap aspect ratio $\alpha$; constructor auto-calibration of $(\alpha, \beta)$ from $(E_{\text{oed}}^{\text{ref}}, K_0^{\text{nc}})$. |
| **Dual Lode Shapes (`mc` & `mn`)** | Fixed William–Warnke smooth approximation in `MC_YF` with 12-eval FD gradients. | Exact principal-stress hexagonal **`mc`** + Panteghini–Lagioia (2014) strictly convex **`mn`** with closed-form first & second derivatives. |
| **Dimensional Views** | **3D-only** (`"ThreeDimensional"`). Cannot be used in 2D `PlaneStrain` or `AxiSymmetric` elements. | **All views from 1 class** (`ThreeDimensional`, `PlaneStrain`, `AxiSymmetric`, `PlaneStress`) via `DimMode` (like `LadrunoJ2.h`). |
| **Compile Budget & MPI Wire** | 46-template `/bigobj` monster TU; `sendSelf`/`recvSelf` not implemented. | Fast standalone TU (`<2 s`); flat `Vector` wire serialization; WP-99 refusal latch; OpenSees-free C++ kernel testable standalone. |

---

## 4. Deep-Research Validation Campaign for `LadrunoHS`

Every benchmark below is tied to a reference and dataset in our Zotero **`Soil` (`LDCNC45G`)** and **`TIMs` (`HTBQAHHH`)** collections.

### Tier 1 — Analytical & Closed-Form Single-Point Oracles (`tests/verification/hs_oracle.py`)

| Gate ID | Test Path | Exact Analytical Reference | Target Tolerance | Zotero Key |
|---|---|---|---|---|
| **T1-A** | **Auto-Calibration of $(\alpha, \beta)$ & Virgin 1D Oedometer ($\varepsilon_2 = \varepsilon_3 = 0$)** | Along virgin $K_0^{\text{nc}}$ compression from $\sigma'_1 = p_0$ to $5 p^{\text{ref}}$ (where both $f^s=0$ and $f^c=0$ are simultaneously active), lateral stress must satisfy $\sigma'_3 / \sigma'_1 \to K_0^{\text{nc}}$ and axial tangent stiffness must satisfy $\mathrm{d}\sigma'_1 / \mathrm{d}\varepsilon_1 = E_{\text{oed}}^{\text{ref}}(\sigma'_1 / p^{\text{ref}})^m$. | $\|\sigma'_3/\sigma'_1 - K_0^{\text{nc}}\| < 10^{-3}$, $\|E_{\text{oed}}^{\text{num}}/E_{\text{oed}}^{\text{ana}} - 1\| < 10^{-3}$ | `TU759V8E` (Schanz et al. 1999), `D4IVV2XD` (Obrzud & Truty 2018) |
| **T1-B** | **Isotropic Drained Triaxial Compression (`CIDC`, $\sigma'_3 = \text{const}$)** | At constant $\sigma'_3$ (so $E_i, E_{ur}, q_f, q_a$ are constant and cap is inactive once calibrated for $OCR \gg 1$ or pure shear path), total axial strain is exact: $\varepsilon_1(q) = \frac{q}{E_{ur}} + \frac{\gamma^p}{1 + \sin\psi_m^*}$ with $\gamma^p(q) = \frac{2}{E_i}\frac{q}{1 - q/q_a} - \frac{2q}{E_{ur}}$ for $q < q_f$, and $q = q_f$ plateau for $\gamma^p \ge \gamma^p_f$. | Relative $L_\infty$ error on $q(\varepsilon_1) < 10^{-10}$ (`mc` and `mn` bit-identical at $\vartheta = +\pi/6$!) | `TU759V8E` (Schanz et al. 1999) |
| **T1-C** | **Rowe Stress–Dilatancy & $e_{\max}$ Dilatancy Cutoff** | In `CIDC` at constant $\sigma'_3$, plastic volumetric rate vs. plastic shear rate is $\frac{\dot{\varepsilon}_v^{p,s}}{\dot{\gamma}^p} = -\sin\psi_m = -\frac{\sin\phi_m - \sin\phi_{cv}}{1 - \sin\phi_m\sin\phi_{cv}}$, dropping to $0$ when $e \ge e_{\max}$. | $\|\dot{\varepsilon}_v^{p,s}/\dot{\gamma}^p + \sin\psi_m\| < 10^{-10}$; exact flat $\varepsilon_v$ plateau after $e = e_{\max}$. | `WBSTR5D2` (Schanz & Vermeer 1996) |
| **T1-D** | **`HSsmall` Monotonic & Cyclic Simple Shear Backbone ($G_s/G_0$ & Masing Damping $\xi$)** | Under pure elastic/paraelastic cyclic shear (high confinement / small strain $\gamma \in [10^{-6}, 10^{-2}]$), secant modulus follows $G_s/G_0 = \frac{1}{1 + 0.385\,\gamma/\gamma_{0.7}}$, tangent modulus floors at $G_{ur}$ for $\gamma \ge \gamma_{\text{cut}}$, and symmetric loop hysteretic damping ratio $\xi(\gamma) = \frac{\Delta W}{4\pi W_s}$ matches Benz (2007) Eq. (5.26). | Backbone error $< 10^{-10}$; Masing loop area error $< 10^{-3}$ | `5PU9KPRT` (Benz 2007), `2U66T7JX` (Benz et al. 2009) |
| **T1-E** | **`HSsmall` Multiaxial Strain-Path Rotation ($90^\circ$ and $180^\circ$ $\Delta\mathbf{e}$ Reversals)** | Strain path in $(e_{11}-e_{22}, e_{12})$ plane: load to $\gamma = 5\gamma_{0.7}$, then reverse $180^\circ$ (tangent modulus must jump back to $G_0$), vs. turn $90^\circ$ orthogonal (partial projection of $\mathbf{H}_{kl}$ per Benz 2007 Ch. 5). | Exact recovery $G_t \to G_0$ on $180^\circ$ reversal; orthogonal turn matches Benz (2007) analytical eigenvector projection. | `5PU9KPRT` (Benz 2007), `2U66T7JX` (Benz et al. 2009) |

---

### Tier 2 — Lode-Angle Bounding (`-lodeShape mc` vs. `mn`) & Multi-Surface Algorithmic Tangent Gates

| Gate ID | Verification Test | Expected Bounding / Numerical Behavior | Zotero Key |
|---|---|---|---|
| **T2-A** | **$\pi$-Plane Radial Failure & Yield Envelope Sweep ($\vartheta \in [-30^\circ, +30^\circ]$, $b \in [0, 1]$)** | Sweep 25 radial deviatoric paths at constant $p' = 100\text{ kPa}$ from TXC ($\vartheta = +30^\circ, b=0$) through Plane Strain ($b \approx 0.3$) to TXE ($\vartheta = -30^\circ, b=1$). Verify: (1) `mc` and `mn` coincide to $<10^{-12}$ at $\vartheta = \pm 30^\circ$; (2) for all $\vartheta \in (-30^\circ, +30^\circ)$, $q_{\text{MN}}(\gamma^p, \vartheta) > q_{\text{MC}}(\gamma^p, \vartheta)$ and $q_{f,\text{MN}}(\vartheta) > q_{f,\text{MC}}(\vartheta)$, with maximum ratio at $\vartheta \approx 0^\circ\text{--}10^\circ$ matching Matsuoka & Nakai (1974) and Panteghini & Lagioia (2014) Fig. 5b. | `C93MGN9R` (Panteghini & Lagioia 2014), `G95BCHF4` (Matsuoka & Nakai 1974), `RJRIJ5TU` (Mortara 2026) |
| **T2-B** | **Plane-Strain Compression (`PSC`, $\varepsilon_{zz} = 0$) Bounding Comparison** | Run single-element `PlaneStrain` compression on Loose & Dense Hostun Sand (`TU759V8E`) and Toyoura Sand (`2XI5UT3Y`): compare peak mobilized friction angle $\phi_{\text{PS}}$ and stiffness between `-lodeShape mc` (lower bound, $\phi_{\text{PS}}^{\text{MC}} = \phi_{\text{TXC}}$) and `-lodeShape mn` (upper bound, $\phi_{\text{PS}}^{\text{MN}} \approx \phi_{\text{TXC}} + 4^\circ$ for dense sand). | `2XI5UT3Y` (Tatsuoka et al. 1986), `C93MGN9R` (Panteghini & Lagioia 2014) |
| **T2-C** | **Central Finite-Difference Algorithmic Tangent Gate Across All 7 Active-Set Regimes** | Perturb all 6 Voigt strain components by $\pm \delta\varepsilon$ ($\delta\varepsilon = 10^{-7}$) and verify $\|\mathbb{C}^{\text{alg}} - \mathbb{C}^{\text{FD}}\|_F / \|\mathbb{C}^{\text{FD}}\|_F < 10^{-6}$ in: (1) Elastic / `HSsmall` paraelastic, (2) Smooth shear cone $f^s$, (3) Shear Lode apex (`mc` edge), (4) Ultimate failure $f^f$, (5) Elliptical cap $f^c$ only, (6) Simultaneous Shear + Cap corner ($f^s=0 \land f^c=0$), (7) Rankine tension cutoff $f^t$. | `C93MGN9R` (Panteghini & Lagioia 2014 §8), `9AV9TN8R` (Pérez-Foguet et al. 2000) |

---

### Tier 3 — Laboratory Element Calibration & Cross-Model Benchmark Suite

1. **Dataset 3.1 — Loose & Dense Hostun RF Sand** (Schanz, Vermeer & Bonnier 1999 [`TU759V8E`]; Benz 2007 [`5PU9KPRT`]).
2. **Dataset 3.2 — Brinkgreve, Engin & Engin (2010) [`ZRD4V6H8`] Empirical Sand Suite** ($RD = 30\%, 50\%, 80\%$).
3. **Dataset 3.3 — Karlsruhe Fine Sand Cross-Model Comparison** (Wichtmann & Triantafyllidis 2016 [`X7LKYZ53`]; Wichtmann et al. 2019 [`KYI7Q4TV`]; Machaček et al. 2021 [`2JN8QSRB`]).
4. **Dataset 3.4 — Soft & Stiff Bangkok Clay Undrained (`CIUC`) & Oedometer Suite** (Surarak et al. 2012 [`TTAH5N7C`]).

---

### Tier 4 — Boundary-Value Problem (BVP) Benchmarks (`2D PlaneStrain`, `AxiSymmetric`, and `$u\text{--}p$`)

1. **BVP-1 — Strip (`PlaneStrain`) vs. Circular (`AxiSymmetric`) Rough Footing Bearing Capacity & Load–Settlement Bounding (`mc` vs. `mn`)** (Loukidis & Salgado 2009 [`JK8ZJP9D`], 2011 [`9GMR7TZD`]; Panteghini & Lagioia 2014 [`C93MGN9R`]; `138_footing_sas_me_ab.md`).
2. **BVP-2 — DGGT Deep Tieback Excavation in Berlin Sand** (Schweiger 2002 [`4HC3USHC`]; Benz 2007 [`5PU9KPRT`] Ch. 7).

---

## 5. Summary of Items & Structured Notes in Zotero (`Soil` `LDCNC45G`)

| Zotero Item Key | Citation | DOI / Type | Validation Note Key |
|---|---|---|---|
| **`TU759V8E`** | **Schanz, T., Vermeer, P. A., & Bonnier, P. G. (1999).** *The hardening soil model: Formulation and verification.* Beyond 2000 in Computational Geotechnics, 281–296. | `10.1201/9781315138206-27` | `4QGISB79` |
| **`5PU9KPRT`** | **Benz, T. (2007).** *Small-Strain Stiffness of Soils and its Numerical Consequences.* PhD Thesis, Mitt. 55, Univ. Stuttgart. | `phdthesis` | `5Q4E8UZV` |
| **`2U66T7JX`** | **Benz, T., Vermeer, P. A., & Schwab, R. (2009).** *A small-strain overlay model.* IJNAMG, 33(1), 25–44. | `10.1002/nag.701` | `776T6FEC` |
| **`WBSTR5D2`** | **Schanz, T., & Vermeer, P. A. (1996).** *Angles of friction and dilatancy of sand.* Géotechnique, 46(1), 145–151. | `10.1680/geot.1996.46.1.145` | — |
| **`GKAKU3EH`** | **Schanz, T., & Vermeer, P. A. (1998).** *On the stiffness of sands.* Géotechnique, 48(3), 383–387. | `journalArticle` | — |
| **`ZRD4V6H8`** | **Brinkgreve, R. B. J., Engin, E., & Engin, H. K. (2010).** *Validation of empirical formulas to derive model parameters for sands.* NUMGE 2010, 137–142. | `inproceedings` | `TNNNS6FX` |
| **`TTAH5N7C`** | **Surarak, C., et al. (2012).** *Stiffness and strength parameters for hardening soil model of soft and stiff Bangkok clays.* Soils and Foundations, 52(4), 682–697. | `10.1016/j.sandf.2012.07.009` | `KADATH9D` |
| **`CIVJHKRF`** | **Calvello, M., & Finno, R. J. (2004).** *Selecting parameters to optimize in model calibration by inverse analysis.* Computers and Geotechnics, 31(5), 411–424. | `10.1016/j.compgeo.2004.03.004` | `49Q2HKRK` |
| **`D4IVV2XD`** | **Obrzud, R. F., & Truty, A. (2018).** *The Hardening Soil Model — A Practical Guidebook.* Z_Soil.PC 100701 report. | `techreport` | `6XFBNKZT` |
| **`4HC3USHC`** | **Schweiger, H. F. (2002).** *Results from numerical benchmark exercises in geotechnics.* NUMGE 2002, Paris, 305–314. | `inproceedings` | `UPJ39XZR` |
| **`2JN8QSRB`** | **Machaček, J., Staubach, P., Tafili, M., Zachert, H., & Wichtmann, T. (2021).** *Investigation of three sophisticated constitutive soil models...* Computers and Geotechnics, 138, 104276. | `10.1016/j.compgeo.2021.104276` | `XK96X6IF` |
| **`C93MGN9R`** | **Panteghini, A., & Lagioia, R. (2014).** *A fully convex reformulation of the original Matsuoka–Nakai failure criterion and its implicit numerically efficient integration algorithm.* IJNAMG, 38, 593–614. | `10.1002/nag.2228` | `EJSXQHSD` |
| **`X7LKYZ53`** | **Wichtmann, T., & Triantafyllidis, T. (2016).** *An experimental database for the development, calibration and verification of constitutive models for sand...* Acta Geotechnica, 11, 739–761. | `10.1007/s11440-015-0402-z` | `WZG92G8T` |
| **`JK8ZJP9D`, `9GMR7TZD`** | **Loukidis, D., & Salgado, R. (2009, 2011).** *Bearing capacity of strip and circular footings in sand using finite elements* / *Effect of relative density and stress level...* | Géotechnique | `HCY4BA8E` |
