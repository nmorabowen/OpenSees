---
title: "ADR 40 — Complex / state-space modal analysis for non-classically-damped systems (LadrunoComplexEigen): design spec"
project: Ladruno
type: ADR / design spec
status: draft — design only, NO code
priority: high
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"           # PARENT: §2a / ADR-A — the complex-modal quick win
  - "[[modal_gap_study/02_abaqus_theory]]"       # PRIMARY: TG §2.5.1 *COMPLEX FREQUENCY modal-projection QZ
  - "[[modal_gap_study/04_lsdyna_theory]]"       # LS-DYNA LCPACK=3 → ARPACK nonsymmetric quadratic (SMP-only)
  - "[[modal_gap_study/01_opensees_current_state]]" # ground truth: eigen / modalProperties / NO assembled-C
  - "[[44_ladruno_frequency_domain_adr]]"        # SIBLING (planned): FRF/SSD/random — shares the Φ + reduced-matrix idea
  - "[[42_ladruno_buckling_adr]]"                # SIBLING (planned): prestressed modal + linear buckling
  - "[[43_ladruno_feast_eigensolver_adr]]"       # SIBLING (planned): robust + parallel eigensolver; FEAST re-hosts THIS at scale
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
  - "[[LEDGER_vanilla_files]]"
tags: [adr, solver, dynamics, modal, complex-modes, non-classical-damping, state-space, qz, base-isolation, ssi, dampers, lapack, zggev]
updated: 2026-06-22
---

# ADR 40 — `LadrunoComplexEigen` (complex / state-space modal analysis)

> **Strategic role (load-bearing assessment — see [[modal_gap_study/00_SYNTHESIS]] §6).**
> **Domain-enabling / load-bearing for our research portfolio.** SSI, DRM, base isolation, and
> supplemental dampers are all *non-classically damped* — real modes give wrong damping ratios for
> every one of them, so this is the correct lens for the systems we actually model. It is also the
> prerequisite for reduced-order modeling of those systems (candidate ADR 45). **Recommended build
> first** — cheap serial proof (reuses existing `eigen`), low risk, directly serves the research.
> Re-hosted at scale by ADR 43's complex contours.

**Status:** draft. **Design only — no code has landed.** classTag **33019
(`LadrunoComplexEigen`, RESERVED, not yet built)** in the ND/solver band (33015 RCConcrete,
33016 LogStrain2D-reserved, 33017 LadrunoConcrete3D-reserved, 33018 LadrunoRCFiniteStrain —
33019 is the next free slot; verified against `LEDGER_implementations.md` 2026-06-22).

This is the **#1 target** of the modal-gap family ([[modal_gap_study/00_SYNTHESIS]] §5, ADR-A):
the lowest-core-risk, highest-value member. It adds **complex (damped) modal analysis** for
**non-classically-damped** systems by **projecting** the assembled $M$, $C$, $K$ onto a
**precomputed real undamped modal basis** $\Phi$ (from an ordinary `eigen` run) and solving the
small $2p\times2p$ **state-space generalized nonsymmetric (QZ)** eigenproblem with LAPACK.

> [!info] What this is, in one line
> After `eigen N`, project $\tilde M=\Phi^TM\Phi$, $\tilde C=\Phi^TC\Phi$ (full, non-classical),
> $\tilde K=\Phi^TK\Phi$; linearize to first-order state space; solve a $2p\times2p$ dense QZ
> with `dggev`/`zggev`; recover undamped frequency $\omega_0=|\lambda|$, damped frequency
> $\omega_d=\mathrm{Im}\,\lambda$, modal damping ratio $\zeta=-\mathrm{Re}\,\lambda/|\lambda|$,
> and **complex (phased) mode shapes** $\psi=\Phi z$. No assembly-core change, no new big
> eigensolver — just an **assembled-$C$ accessor**, a small **QZ wrapper**, and complex
> post-processing.

---

## 1. Driver & goal — why real modes are the wrong answer for damped structures

OpenSees `eigen` solves only the **undamped, symmetric** generalized eigenproblem
$(K-\omega^2 M)\phi=0$ (ground truth: `StaticAnalysis.cpp:285–330`,
`ArpackSolver.cpp:138–141` errors out for the non-generalized case; the assembled "A" is always
$K_t$ via `addKtToTang`, never $C$). The eigenvectors are **real**: every DOF reaches its peak
at the same instant, the mode has no phase. That is correct **only** when the damping is
**classical (proportional)** — i.e. $C$ is simultaneously diagonalizable with $M$ and $K$, so
$\Phi^TC\Phi$ is already diagonal and the modes decouple.

For the structures the fork most cares about, **the damping is *not* classical**:

- **Base isolation.** A nearly-rigid superstructure on a soft, highly-damped isolation layer
  (lead-rubber / friction-pendulum bearings idealized as concentrated dashpots) puts almost all
  the energy dissipation in *one localized interface*. $C$ is spatially lumped, not proportional
  to $M$ or $K$ — the textbook non-classical case. The first "isolation mode" has a damping
  ratio of 15–30%; treating it with a real mode and a hand-assigned $\zeta$ misstates both the
  **isolation period** and the **isolation-layer damping** that the whole design hinges on.
- **Supplemental dampers.** Fluid viscous dampers, viscoelastic dampers, and tuned mass dampers
  are *deliberately localized* energy sinks. Their entire purpose is to add damping where the
  structure has none; the resulting $C$ is the antithesis of mass/stiffness-proportional.
  Engineers need each mode's *actual* $\zeta$ (to verify the damper achieved the target supplemental
  damping) — which only a complex eigen delivers.
- **Soil–structure interaction (SSI).** Radiation damping into the half-space (dashpots /
  absorbing-boundary tractions on the foundation interface — cf. the fork's
  [[project_absorbing_pml_guide|Lysmer/ASD absorbing boundaries]]) is again a localized, frequency-
  rich $C$ concentrated at the soil interface. Foundation modes are heavily, non-classically damped.
- **Dissimilar-material / joint damping.** Composite assemblies and bolted/contact joints with
  per-material or interface damping produce $C$ that no single Rayleigh pair can represent.

**What the user actually needs (the deliverables this ADR exists to produce):**

1. **The true modal damping ratio $\zeta_k$ of each mode** — not a guessed global Rayleigh
   $\alpha,\beta$ pair, but the $\zeta$ that the localized dashpots / bearings / soil *actually*
   impose, mode by mode. This is the single most-requested number for an isolation or
   supplemental-damping design review.
2. **Damped natural frequency $\omega_{d,k}$** distinct from the undamped $\omega_{0,k}$ — for
   heavily-damped isolation modes the two differ enough to matter.
3. **Complex (phased) mode shapes** $\psi_k$ — the entries no longer peak simultaneously; the
   phase distribution is the *qualitative signature* of where the energy is being dissipated.
   Travelling-wave behavior across a damped interface is invisible to real modes.
4. A principled way to feed mode-by-mode $\zeta_k$ into downstream modal-superposition /
   response-spectrum work ([[44_ladruno_frequency_domain_adr|ADR 44]]) instead of a single
   blanket damping value.

Without this, the fork can *run* a damped transient (direct integration with element dashpots),
but it **cannot answer "what is the damping ratio of the isolation mode?"** — a question every
base-isolation and supplemental-damping project asks at the modeling-review stage.

---

## 2. Decision summary

**Adopt the Abaqus `*COMPLEX FREQUENCY` / LS-DYNA `LCPACK=3` strategy: modal-projection complex
eigen on a precomputed real basis, solved as a small dense QZ.** Concretely:

1. **Reuse the real undamped basis.** Require a prior `eigen N` (any backend) that has populated
   $\Phi=[\phi_1\dots\phi_p]$ on the nodes and $\omega_\alpha^2$ in the domain.
2. **Assemble $M$, $C$, $K$ as operands** and form the three small reduced matrices
   $\tilde M=\Phi^TM\Phi$, $\tilde C=\Phi^TC\Phi$ (**full**, because the damping is non-classical),
   $\tilde K=\Phi^TK\Phi$.
3. **Linearize the reduced quadratic problem** $(\lambda^2\tilde M+\lambda\tilde C+\tilde K)z=0$
   to a $2p\times2p$ first-order generalized **nonsymmetric** pencil and solve it with LAPACK
   **`dggev`** (real reduced matrices, the common case) or **`zggev`** (complex; future-proofs
   structural/complex-stiffness damping).
4. **Recover** $\omega_0=|\lambda|$, $\omega_d=\mathrm{Im}\,\lambda$,
   $\zeta=-\mathrm{Re}\,\lambda/|\lambda|$, complex modes $\psi=\Phi z$.

**Why NOT the full $N$-dimensional quadratic eigenproblem** $(\lambda^2M+\lambda C+K)\psi=0$:

- It is **nonsymmetric and complex at full size $N$** — none of OpenSees' eigensolvers handle it
  (ARPACK path is symmetric-generalized shift-invert only; `ArpackSolver.cpp:138–141`). Building a
  full-$N$ complex Krylov/QZ solver is a *major* core effort — exactly the [[43_ladruno_feast_eigensolver_adr|ADR 43]]
  (FEAST complex contours) scope, deferred.
- **Both mature codes deliberately avoid it.** Abaqus `*COMPLEX FREQUENCY` projects onto the real
  modal basis and solves the *reduced* problem (TG §2.5.1). LS-DYNA `LCPACK=3` folds damping into
  first-order state space and runs ARPACK on the **reduced/extracted** subspace — and even there it
  is **SMP-only** (`*CONTROL_IMPLICIT_SOLVER` Remark 8). The cost philosophy is identical to
  mode-based steady-state dynamics: trade an $N$-dim nonsymmetric solve for a $p$-dim one.
- The projection is **exact-enough** when the retained $p$ real modes span the dynamics of
  interest: damping that *perturbs but does not radically reshape* the low modes (the design regime
  for isolation / supplemental dampers) is captured faithfully. Completeness (enough $p$) is the
  user's responsibility and is checked (see §8/§9).
- It is a **low-core-risk quick win**: $p\ll N$, so the dense QZ is microseconds; the only genuinely
  new plumbing is the **assembled-$C$ accessor** (§4.6) and the small LAPACK wrapper.

**Deliberate v1 simplifications** (revisit later): serial only (the projection + dense QZ are cheap;
parallel re-hosting is [[43_ladruno_feast_eigensolver_adr|ADR 43]]'s job); Rayleigh + element/material
dashpots for $C$ (no structural / complex-stiffness damping in v1 — that needs `zggev` and belongs
with [[44_ladruno_frequency_domain_adr|ADR 44]]); the real undamped $\Phi$ as the projection basis
(no second-pass basis enrichment).

---

## 3. Scope-fence & classTag

**`#define LADRUNO_TAG_ComplexEigen 33019`** — analysis/solver object tag, **RESERVED, not yet
built**. Recorded here so a sibling never collides (next free after 33018 LadrunoRCFiniteStrain;
33016/33017 stay reserved for LogStrain2D / LadrunoConcrete3D per their ADRs).

**In scope (this ADR):**

- **Modal-projection complex eigen on a *precomputed* $\Phi$.** Consume a prior `eigen` result,
  project $\{M,C,K\}$, solve the $2p\times2p$ QZ, report complex frequencies / damping ratios /
  complex modes.
- An **assembled-$C$ accessor** path: trivial closed form for Rayleigh ($\alpha\tilde M+\beta\tilde K$
  + the $\beta_0,\beta_c$ variants), plus a general element-by-element `getDamp()` assembly for
  localized dampers (the key design decision, §4.6).
- Tcl + openseespy commands (`complexEigen` / `eigen -damped`) and complex output via an extended
  `modalProperties` query (§5).
- Serial, dense, small-$p$ reduced solve.

**NOT in scope (explicitly deferred):**

| Out of scope | Where it lives |
|---|---|
| **Full $N$-dim quadratic eigenproblem** $(\lambda^2M+\lambda C+K)\psi=0$ at scale | [[43_ladruno_feast_eigensolver_adr\|ADR 43]] (complex FEAST contours) |
| **Direct (un-projected) complex eigen** for large/nonsymmetric systems | [[43_ladruno_feast_eigensolver_adr\|ADR 43]] |
| **Parallel / distributed** complex eigen (SP/MP) | [[43_ladruno_feast_eigensolver_adr\|ADR 43]] — FEAST re-hosts THIS reduced case at scale |
| **Structural / complex-stiffness damping** ($F_s=i\,s\,F_{int}$) | needs `zggev` end-to-end; pairs with [[44_ladruno_frequency_domain_adr\|ADR 44]] SSD |
| **Gyroscopic / rotating-frame** ($C^R=C+M\Omega$) complex modes | future; LS-DYNA Theory §36.1 analogue, not a fork driver yet |
| **Frequency-domain response** (FRF/SSD/random) built on $\zeta_k$ | [[44_ladruno_frequency_domain_adr\|ADR 44]] |
| **Prestressed / geometric-stiffness** complex modes | compose with [[42_ladruno_buckling_adr\|ADR 42]] once both land |

**Sibling relationship.** ADR 40 (complex modes) and ADR 44 (frequency domain) both ride on the
real $\Phi$ + reduced matrices; ADR 44 *consumes* the $\zeta_k$ this ADR produces. ADR 43's complex
FEAST eventually offers the same complex eigenpairs at full scale, making this ADR the **serial,
projection-based** member that proves the physics cheaply first.

---

## 4. Formulation

### 4.1 The damped free-vibration problem

The damped, unforced semidiscrete system is

$$M\,\ddot u + C\,\dot u + K\,u = 0,\qquad u\in\mathbb R^N.$$

Seeking $u(t)=\psi\,e^{\lambda t}$ gives the **quadratic eigenproblem**

$$\big(\lambda^2 M + \lambda C + K\big)\,\psi = 0,$$

whose eigenvalues $\lambda$ and eigenvectors $\psi$ are **complex** when $C$ is non-classical.
Solving this at full size $N$ is nonsymmetric and expensive — avoided (§2).

### 4.2 State-space (first-order) form

Introduce the state $x=\begin{Bmatrix}u\\\dot u\end{Bmatrix}\in\mathbb R^{2N}$. The standard
symmetric-pencil linearization (preserving symmetry of $M$, $K$) is

$$\underbrace{\begin{bmatrix} C & M \\ M & 0 \end{bmatrix}}_{\mathcal A}\dot x
  + \underbrace{\begin{bmatrix} K & 0 \\ 0 & -M \end{bmatrix}}_{\mathcal B}x = 0
  \;\Longrightarrow\;
  \mathcal A\,x = -\lambda\,\mathcal B\,x,$$

i.e. a $2N\times2N$ generalized eigenproblem $\mathcal A\,\xi = -\lambda\,\mathcal B\,\xi$ with
$\xi=\begin{Bmatrix}\psi\\\lambda\psi\end{Bmatrix}$. (The $\mathcal A,\mathcal B$ blocks are
symmetric but the *pencil* is indefinite, so the eigenpairs are complex — a genuinely nonsymmetric
solve is still required. This is the same linearization Abaqus TG §2.5.1 uses.)

### 4.3 Projection onto the real undamped basis

Let $\Phi=[\phi_1,\dots,\phi_p]\in\mathbb R^{N\times p}$ collect the $p$ **real, mass-normalized**
undamped modes from a prior `eigen`, so

$$\Phi^T M\,\Phi = I_p,\qquad \Phi^T K\,\Phi = \mathrm{diag}(\omega_1^2,\dots,\omega_p^2).$$

(If `eigen` returned modes normalized to unit largest entry rather than unit modal mass, we
re-normalize internally: $\phi_\alpha \leftarrow \phi_\alpha/\sqrt{\phi_\alpha^TM\phi_\alpha}$;
`modalProperties` already forms $\phi^TM\phi$ — `DomainModalProperties.cpp:602`.)

Approximate the physical response in the reduced coordinates $z\in\mathbb C^p$ via $u\approx\Phi z$
and Galerkin-project the quadratic problem with $\Phi^T(\cdot)$:

$$\boxed{\;\tilde M=\Phi^T M\,\Phi=I_p,\qquad
        \tilde C=\Phi^T C\,\Phi,\qquad
        \tilde K=\Phi^T K\,\Phi=\mathrm{diag}(\omega_\alpha^2)\;}$$

The crucial object is the **reduced damping matrix** $\tilde C$. Its off-diagonal entries

$$\tilde C_{\alpha\beta}=\phi_\alpha^T C\,\phi_\beta\quad(\alpha\neq\beta)$$

are the **inter-modal damping coupling** that classical modal analysis throws away. If $C$ is
classical, $\tilde C$ is diagonal and the modes decouple (each gives a real mode with
$\zeta_\alpha=\tilde C_{\alpha\alpha}/2\omega_\alpha$); a **non-zero off-diagonal $\tilde C$ is the
mathematical fingerprint of non-classical damping** and the reason a complex solve is needed.

### 4.4 The reduced $2p\times2p$ QZ problem

Substituting the reduced matrices into the state-space form (§4.2) gives the small pencil

$$\underbrace{\begin{bmatrix} \tilde C & \tilde M \\ \tilde M & 0 \end{bmatrix}}_{\tilde{\mathcal A}\,(2p\times2p)}
  \begin{Bmatrix} z \\ \lambda z \end{Bmatrix}
  = -\lambda
  \underbrace{\begin{bmatrix} \tilde K & 0 \\ 0 & -\tilde M \end{bmatrix}}_{\tilde{\mathcal B}\,(2p\times2p)}
  \begin{Bmatrix} z \\ \lambda z \end{Bmatrix},$$

a **dense generalized nonsymmetric eigenproblem** $\tilde{\mathcal A}\,w = -\lambda\,\tilde{\mathcal B}\,w$
of size $2p$. Because $p\ll N$ (typically $p=\mathcal O(10\text{–}100)$), the QZ factorization is
negligible cost. Solve with LAPACK **`dggev`** (real $\tilde{\mathcal A},\tilde{\mathcal B}$ — the
Rayleigh + viscous-dashpot case, where everything is real) returning eigenvalues as
$(\alpha_k,\beta_k)$ pairs with $\lambda_k=\alpha_k/\beta_k$; or **`zggev`** when complex damping is
present (deferred, §3). `dggev` is *already linked* into the fork — it is the backend of
`FullGenEigenSOE`/`FullGenEigenSolver` (`fullGenLapack`) and of the per-element
`CriticalTimeStep` solve — so **no new external dependency**.

> [!note] Variant linearization
> The companion/$\mathcal L_1$ form
> $\begin{bmatrix}0 & I\\ -\tilde M^{-1}\tilde K & -\tilde M^{-1}\tilde C\end{bmatrix}$ turns the
> problem into a *standard* (non-generalized) eigenproblem solvable by `dgeev`, since $\tilde M=I_p$
> after mass-normalization makes $\tilde M^{-1}$ trivial. We prefer the **symmetric-block generalized
> $\tilde{\mathcal A},\tilde{\mathcal B}$ pencil + `dggev`** because (a) it never forms $\tilde M^{-1}$
> (robust if $\Phi$ is only approximately mass-orthonormal), and (b) the $(\alpha,\beta)$ return
> cleanly flags **infinite / rigid-body** eigenvalues ($\beta_k\to0$, §4.7). The `dgeev` route is a
> documented fallback.

### 4.5 Complex-eigenpair → frequency / damping / mode recovery

Each reduced eigenvector $w_k=\begin{Bmatrix}z_k\\\lambda_k z_k\end{Bmatrix}$ yields the top block
$z_k\in\mathbb C^p$, and the **physical complex mode** is the complex combination of real modes

$$\boxed{\;\psi_k=\Phi\,z_k\in\mathbb C^N.\;}$$

For underdamped modes the eigenvalues come in **complex-conjugate pairs**

$$\lambda_k = -\zeta_k\,\omega_{0,k}\;\pm\;i\,\omega_{d,k},$$

from which the engineering quantities read directly:

$$\omega_{0,k}=|\lambda_k|=\sqrt{(\mathrm{Re}\,\lambda_k)^2+(\mathrm{Im}\,\lambda_k)^2}
   \quad\text{(undamped natural circular frequency, rad/s)},$$
$$\omega_{d,k}=\mathrm{Im}\,\lambda_k \quad\text{(damped circular frequency, rad/s)},$$
$$\boxed{\;\zeta_k=-\dfrac{\mathrm{Re}\,\lambda_k}{|\lambda_k|}\;}
   \quad\text{(modal damping ratio)},\qquad
   T_{0,k}=\frac{2\pi}{\omega_{0,k}},\quad f_{d,k}=\frac{\omega_{d,k}}{2\pi}.$$

**Classification of $\lambda_k$:**
- **Underdamped** ($0<\zeta_k<1$): conjugate pair, $\omega_{d,k}\neq0$ — the oscillatory modes.
- **Overdamped** ($\zeta_k\ge1$): real, negative $\lambda_k$ (non-oscillatory), no conjugate
  partner; reported with $\omega_{d,k}=0$.
- We report one representative of each conjugate pair (the $+i\omega_d$ root), with a flag, so the
  user sees $p$ "physical" modes, not $2p$ duplicated roots.

A complex mode shape has a **phase distribution** — entries $\psi_{k,j}=|\psi_{k,j}|e^{i\theta_{k,j}}$
do not all reach their peak simultaneously. This phase spread is the qualitative output that real-mode
analysis cannot produce, and is exactly what reveals *where* the localized damping is acting.

### 4.6 Assembled-$C$ extraction — **the key design decision**

The projection needs $C$ as a usable operand, but **OpenSees never assembles a standalone $C$
matrix**. Confirmed in source: `EigenSOE` exposes only `addA`/`addM` (no `addC` —
`EigenSOE.h:48–54`); `EigenIntegrator` has only `formEleTangK`/`formEleTangM` (no `formEleTangC` —
`EigenIntegrator.cpp:191–204`); damping enters *only inside transient integrators* as a scaled
contribution to the effective tangent ($K_{\text{eff}}=c_1K+c_2C+c_3M$ via `addCtoTang(c2)` —
`Newmark.cpp:284–309`). So we must build $C$ ourselves. **Two complementary routes, chosen by what
the model contains:**

**Route A — Rayleigh damping (closed form, trivial).** When damping is purely Rayleigh, the global
$C=\alpha_M M+\beta_K K_t+\beta_{K0}K_0+\beta_{Kc}K_c$ — this is *exactly* `Element::getDamp()`
(`Element.cpp:210–231`, formula
$C_{\text{el}}=\alpha_M M+\beta_K K_t+\beta_{K0}K_0+\beta_{Kc}K_c$). Because projection is linear,
**we never assemble $C$ at all** — we project the *factors*:

$$\tilde C = \alpha_M\,\underbrace{\Phi^TM\Phi}_{\tilde M=I}
            + \beta_K\,\underbrace{\Phi^TK_t\Phi}_{\tilde K=\mathrm{diag}(\omega_\alpha^2)}
            + \beta_{K0}\,\Phi^TK_0\Phi + \beta_{Kc}\,\Phi^TK_c\Phi.$$

For the common pair $C=\alpha_M M+\beta_K K$ this collapses to the **diagonal**
$\tilde C_{\alpha\alpha}=\alpha_M+\beta_K\omega_\alpha^2$, i.e. $\zeta_\alpha=\alpha_M/2\omega_\alpha
+\beta_K\omega_\alpha/2$ — the classical result, recovered *for free* and used as the P1 unit-test
oracle (§7/§8). Rayleigh alone is classical, so this route is a **correctness anchor**, not the
motivating case.

**Route B — assembled $C$ from element/material dampers (the non-classical core).** Localized
dashpots (Viscous material in a `zeroLength`/`twoNodeLink`, `viscousDamper`, explicit dashpot
elements) contribute damping through their own `getDamp()` override, which returns a *real element
$C$ matrix* (e.g. `ZeroLength::getDamp` `ZeroLength.cpp:914–968`, `TwoNodeLink::getDamp`
`TwoNodeLink.cpp:724–753`; the underlying uniaxial law returns a scalar `getDampTangent()`,
`ViscousMaterial.cpp:154–162`, that the element transforms into its local $C$). To capture these we
**mirror the eigen-assembly pattern** (`StaticAnalysis::eigen` forms $K$ by looping FE_Elements and
calling `addKtToTang`→`addA`): a new **`ComplexEigenIntegrator::formEleTangC`** that calls
`theEle->addCtoTang(1.0)` (which internally invokes `myEle->getDamp()`,
`FE_Element.cpp:372–386`), assembled through a system matrix into a global $C$.

We do **not** widen the `EigenSOE` interface (no new `addC` virtual on the base). Instead the design
adds a small **`LadrunoDampingAssembler`** helper that:
1. takes the active `AnalysisModel` + `DOF_Numberer`,
2. loops `FE_Element`s, calls `addCtoTang(1.0)` to get each element $C_{\text{el}}$ + its DOF `ID`,
   and `addM`-style scatters it into a sparse/dense global $C$ in the same DOF ordering `eigen` used,
3. returns $C$ as a `Matrix`/sparse operand the projection driver consumes.

This is the **lighter-weight alternative** flagged in the recon ("a utility that loops FE_Elements
accumulating `getDamp()` into a matrix, bypassing the SOE"). It touches **no assembly core** — it
reuses `getDamp()`/`addCtoTang`, the very methods Newmark already calls every step.

**Route A + B compose:** the Rayleigh factors and the element-damper $C$ are simply summed (Newmark
does exactly this — `addCtoTang` already returns Rayleigh **plus** the element's material damping,
because `Element::getDamp()` is the Rayleigh part and overrides like `ZeroLength::getDamp` add the
material part on top). So **Route B's assembled $C$ already includes Route A** when the user has set
both — and Route A's closed form is just an *optimization / oracle* for the Rayleigh-only model. The
projection driver therefore always: (i) tries Route A factors if no element returns a non-Rayleigh
$C$; (ii) otherwise assembles the full $C$ via Route B and projects $\tilde C=\Phi^TC\Phi$ as three
dense mat-mats (or, more cheaply, $\tilde C=(C\Phi)^T\Phi$ exploiting symmetry).

> [!warning] Damping-flag quirks to honor (from `LEDGER_quirks` / [[project_damping_channels]])
> `getDamp()` only returns what the element is *configured* to return. **`zeroLength` `-doRayleigh`
> defaults to 0** ⇒ stiffness-Rayleigh gives that element **zero** damping; a `Viscous` material in a
> `zeroLength` with `-doRayleigh 0` contributes via the material branch only. The assembler must take
> `getDamp()` at face value (it already encodes these switches) and the docs must warn users that the
> complex-eigen $C$ is exactly the $C$ their transient analysis would feel — no more, no less. This is
> a *feature* (consistency with the actual dynamic model), but a documentation trap if unflagged.

### 4.7 Rigid-body / $\omega=0$ modes and normalization

- **Rigid-body / $\omega=0$ modes.** If $\Phi$ includes a zero-frequency mode ($\omega_\alpha=0$),
  the corresponding $\tilde K$ diagonal entry is 0; the state-space pencil handles this naturally —
  the QZ returns $\lambda=0$ (or an infinite eigenvalue with $\beta_k\to0$ from the trivial
  $-\tilde M$ block, which the $(\alpha,\beta)$ form of `dggev` flags as
  $|\beta_k|<\varepsilon_{\text{tol}}\|\tilde{\mathcal B}\|$). We **detect and tag** these (report
  $\omega_0=\omega_d=0$, $\zeta$ undefined) rather than dividing through. A mass-proportional
  Rayleigh term on a rigid-body mode gives a real, finite, non-oscillatory decay rate
  ($\lambda=-\alpha_M/2$) which we report as an overdamped/rigid mode.
- **Mode normalization (output convention).** Complex modes have an arbitrary complex scale
  ($\psi_k$ and $c\,\psi_k$, $c\in\mathbb C$, are the same mode). We normalize to a **canonical
  phase + magnitude**: scale so the largest-magnitude entry is real-positive
  ($\psi_k\leftarrow\psi_k\cdot e^{-i\arg\psi_{k,\text{max}}}$, then optionally
  $/|\psi_{k,\text{max}}|$). This makes the *relative* phases meaningful and reproducible. We also
  offer a **state-space mass normalization** ($w_k^{T}\tilde{\mathcal B}\,w_k=1$, the complex analogue
  of unit modal mass) as an option for downstream modal-superposition use.
- **Pairing.** Conjugate pairs are matched by $|\lambda_k-\bar\lambda_j|$ within tolerance and only
  the $+\mathrm{Im}$ representative is reported (with the count of physical vs. raw roots logged), so
  the user sees $p$ modes aligned 1:1 with the input real modes whenever damping is light.

---

## 5. Public API

### 5.1 Tcl

Two spellings, one engine (mirroring how `eigen` accepts both `genBandArpack` and `-genBandArpack`):

```tcl
# Prerequisite: a prior real eigen run that populated Phi on the nodes.
eigen  $numModes                ;# e.g.  eigen 12   (any backend)

# Complex modal solve on the existing Phi (uses all p modes from the last `eigen`):
complexEigen
# or, equivalently, fold into eigen with a flag:
eigen -damped $numModes         ;# runs the real eigen THEN the complex projection in one call

# Options:
complexEigen  -numModes $p              ;# use only the first p of the available real modes
              -normalize {phase|mass}   ;# complex-mode output normalization (default: phase)
              -solver    {dggev|zggev}  ;# real (default) vs complex QZ backend
              -tol       $eps           ;# infinite-eigenvalue / conjugate-pairing tolerance
```

**Returned (Tcl list, one entry per physical mode $k$, ascending $\omega_0$):**
the command returns a flat list the user can `lindex` / reshape, of
`{omega0_k  omega_d_k  zeta_k  Re(lambda_k)  Im(lambda_k)}` quintuples — analogous to how `eigen`
returns the eigenvalue list. Complex mode shapes are *not* dumped to the return list (too large);
they are pushed onto the nodes (§5.3) and queried via `modalProperties`/recorder.

### 5.2 openseespy

```python
# returns a list of dicts (or a structured numpy-friendly list), one per physical mode:
res = ops.complexEigen(p)                 # p optional → all available real modes
# or:
vals = ops.eigen('-damped', p)            # real eigen + complex projection

# res[k] = {
#   'omega0' : float,   # undamped natural circular freq (rad/s)
#   'omegaD' : float,   # damped circular freq (rad/s)
#   'fD'     : float,   # damped freq (Hz)
#   'T0'     : float,   # undamped period (s)
#   'zeta'   : float,   # modal damping ratio
#   'lambda' : complex, # the eigenvalue (Re, Im)
#   'kind'   : str,     # 'underdamped' | 'overdamped' | 'rigid'
# }
```

### 5.3 Complex output via extended `modalProperties` + complex mode shapes

- **Complex mode shapes** are stored on the nodes alongside the existing real eigenvectors: a new
  `Node::setComplexEigenvector(mode, Vector rePart, Vector imPart)` / `getComplexEigenvectors()`
  pair, mirroring `Node::setEigenvector`/`getEigenvectors` (`Node.cpp:1418–1440`). They are queryable
  per node/DOF and recordable (`recorder Node ... -dof ... complexEigen <mode> {real|imag|mag|phase}`).
- **`modalProperties` extension.** Add a `-complex` flag to the existing `modalProperties` command
  (`modal.cpp:32`): when present and a `complexEigen` has run, the printed/returned table gains
  columns **$\omega_{d}$, $\zeta$ (true, from the complex solve), $f_d$**, and the
  participation/effective-mass section additionally reports the **complex participation factors**
  $\tilde\Gamma_{k}=\psi_k^TMR/(\psi_k^TM\psi_k)$ (complex). The real-mode columns are unchanged so
  existing scripts keep working.

**Exact quantities returned by the complex query, per mode:** undamped natural frequency $\omega_0$
& period $T_0$; damped frequency $\omega_d$ & $f_d$; modal damping ratio $\zeta$; the complex
eigenvalue $\lambda$ (real & imaginary parts); mode kind (under/over/rigid); and the complex mode
shape $\psi_k$ (real + imaginary parts, or magnitude + phase) per node/DOF.

---

## 6. OpenSees integration points (file:line)

All citations against the worktree `SRC/` (audited 2026-06-22; see
[[modal_gap_study/01_opensees_current_state]] for the full trace).

**Where $\Phi$ lives (reuse, read-only):**
- Real eigen assembly + storage: `SRC/analysis/analysis/StaticAnalysis.cpp:244` (`eigen`), forms
  $K$/`M` at `:285–323`, solves `:330`, pushes eigenvectors to nodes via
  `theAnalysisModel->setEigenvector` `:339–345`. Twin: `BasicAnalysisBuilder::eigen`
  (`BasicAnalysisBuilder.cpp:916`).
- Per-node eigenvector storage: `Node::setEigenvector` `SRC/domain/node/Node.cpp:1418`,
  `Node::getEigenvectors` `:1436`; eigenvalues in `Domain::getEigenvalues()`.
- Mass-normalization helper already present: $\mathrm{diag}(\Phi^TM\Phi)$ at
  `SRC/runtime/commands/analysis/modal/DomainModalProperties.cpp:602`.

**How to get $M$/$K$ (reuse the eigen assembly pattern):**
- $K$: per-element `zeroTangent(); addKtToTang(1.0); addA(getTangent(),id)`
  (`StaticAnalysis.cpp:285–293`); integrator hook `EigenIntegrator::formEleTangK`
  (`EigenIntegrator.cpp:191–196`).
- $M$: per-element/DOF `addMtoTang(1.0); addM(...)` (`StaticAnalysis.cpp:300–323`);
  `EigenIntegrator::formEleTangM` (`EigenIntegrator.cpp:199–204`).
- We do **not** re-assemble $M$/$K$ at full size for the projection: $\tilde M,\tilde K$ come for
  free from the real eigen ($\tilde M=I$, $\tilde K=\mathrm{diag}(\omega_\alpha^2)$ from
  `Domain::getEigenvalues()`); only $C$ needs new assembly.

**How to get $C$ (NEW — the key plumbing, §4.6):**
- `Element::getDamp()` formula (the Rayleigh route): `SRC/element/Element.cpp:210–231`
  ($C=\alpha_M M+\beta_K K_t+\beta_{K0}K_0+\beta_{Kc}K_c$). Rayleigh factors set/stored:
  `Domain::setRayleighDampingFactors` `Domain.cpp:2119`, element storage
  `Element.cpp:116–121`, `Element.h:169`.
- Element-damper $C$ (the non-classical route): `FE_Element::addCtoTang` `SRC/analysis/fe_ele/FE_Element.cpp:372–386`
  (calls `myEle->getDamp()`); overrides `ZeroLength::getDamp` `ZeroLength.cpp:914–968`,
  `TwoNodeLink::getDamp` `TwoNodeLink.cpp:724–753`; scalar `ViscousMaterial::getDampTangent`
  `ViscousMaterial.cpp:154–162`.
- **NEW `ComplexEigenIntegrator::formEleTangC`** mirroring `EigenIntegrator::formEleTangK` but
  calling `addCtoTang(1.0)`; and **NEW `LadrunoDampingAssembler`** that scatters the per-element $C$
  into a global operand in the eigen DOF ordering (no `EigenSOE` interface change).

**Where to register the command:**
- xara runtime `eigen` parser: `SRC/runtime/commands/analysis/analysis.cpp:250` (add the `-damped`
  branch / a new `complexEigen` next to it); SOE/analysis construction
  `BasicAnalysisBuilder.cpp:871` (`newEigenAnalysis`).
- legacy/parallel + openseespy parser: `SRC/interpreter/OpenSeesCommands.cpp:270`
  (`OpenSeesCommands::eigen`) — add the Python `complexEigen` binding here.
- `modalProperties` driver to extend: `SRC/runtime/commands/analysis/modal/modal.cpp:32`;
  command registration `SRC/runtime/commands/domain/commands.cpp:121` and legacy table
  `SRC/tcl/commands.cpp:1033–1035`.

**New source (all fork-authored, stamped, ledgered):**
- `SRC/analysis/.../LadrunoComplexEigen.{h,cpp}` — the projection + QZ driver (classTag 33019).
- `SRC/analysis/.../LadrunoDampingAssembler.{h,cpp}` — Route-B global-$C$ assembly helper.
- `SRC/analysis/integrator/ComplexEigenIntegrator.{h,cpp}` — `formEleTangC` (mirrors `EigenIntegrator`).
- a thin LAPACK `dggev`/`zggev` wrapper (or reuse `FullGenEigenSolver`'s `dggev` linkage).

---

## 7. Phased roadmap with exit gates

Each phase is a **separate PR based on `ladruno`**, with an adversarial gate sized to risk
([[feedback_adversarial_gate_when]]): P1 introduces novel solver math → **full gate**; P2/P3 are
mechanical extensions → lighter gate carried by tests.

| Phase | Scope | Exit gate (must pass to merge) |
|---|---|---|
| **P0 — Skeleton & QZ wrapper** | `LadrunoComplexEigen` class shell (classTag 33019); LAPACK `dggev` wrapper on a hand-built $2p\times2p$ pencil; no OpenSees coupling yet. | Standalone unit: feed a known $\tilde M,\tilde C,\tilde K$ (numpy-generated), recover $\lambda$ to $1\text{e-}10$ vs `scipy.linalg.eig`. Stamp header, ledger row (RESERVED→active on P1). |
| **P1 — Rayleigh-only projection + QZ (the correctness anchor)** | Route A (project Rayleigh **factors**, no $C$ assembly): $\tilde C=\alpha_M I+\beta_K\mathrm{diag}(\omega^2)$; full `complexEigen` Tcl+py command consuming a prior `eigen`; recover $\omega_0,\omega_d,\zeta$. | **Classical oracle:** for $C=\alpha_M M+\beta_K K$ the complex solve must reproduce $\zeta_\alpha=\alpha_M/2\omega_\alpha+\beta_K\omega_\alpha/2$ and $\omega_{d,\alpha}=\omega_\alpha\sqrt{1-\zeta_\alpha^2}$ to $1\text{e-}8$ on a multi-DOF frame; complex modes collapse to real (zero phase spread). **Full adversarial gate.** |
| **P2 — Assembled $C$ for element dampers (the non-classical core)** | Route B: `ComplexEigenIntegrator::formEleTangC` + `LadrunoDampingAssembler` build global $C$ from `getDamp()`; $\tilde C=\Phi^TC\Phi$ (full); Route A+B compose. | **2-DOF non-classical analytic oracle** (§8) to $1\text{e-}8$; **damped SDOF**; cross-check $\tilde C$ vs `scipy` on the same projected matrices; honor `-doRayleigh` quirk (test that a `Viscous`-in-`zeroLength` model gives the dashpot's true $\zeta$, not zero). Lighter gate (tests carry it). |
| **P3 — Complex mode output + recorder + modalProperties** | `Node::setComplexEigenvector`/`getComplexEigenvectors`; `modalProperties -complex` columns; `recorder Node ... complexEigen <mode> {real|imag|mag|phase}`; complex participation factors. | **Base-isolated stick model** (§8): isolation-mode $\zeta$ and $\omega_d$ match the analytic 2-DOF + a **time-history log-decrement** cross-check to ~2–5%; recorder round-trips complex shapes; banner line + guide. |
| **P4 (optional/deferred)** | `zggov`/`zggev` complex path groundwork for structural damping; hand-off note to ADR 43 (FEAST) + ADR 44 (SSD consumes $\zeta_k$). | Documented hand-off; no functional gate. |

**Recommended stop point for v1:** end of **P3** — that delivers the full serial complex-modal
capability (the driver's actual need). P4 is a seam, not a feature.

---

## 8. Validation / oracle battery

Tiered like the fork's other material/solver ADRs (analytic closed form → projected-matrix
cross-check → physical decay). Tests live in `tests/test_ladrunoComplexEigen.py` (Zone-A) with a
numpy/scipy oracle in `tests/_testbed/complex_eigen_ref.py`.

**Tier 1 — analytic closed-form (the anchors):**

1. **Damped SDOF.** $m\ddot u+c\dot u+ku=0$: closed form $\lambda=-\zeta\omega_n\pm i\omega_n\sqrt{1-\zeta^2}$,
   $\zeta=c/2\sqrt{mk}$. The $p=1$ projection must return $\zeta$, $\omega_0=\omega_n$,
   $\omega_d=\omega_n\sqrt{1-\zeta^2}$ exactly. Also exercises the overdamped ($\zeta>1$, real
   $\lambda$) and rigid ($k=0$) branches.
2. **2-DOF non-classically-damped system (the headline oracle).** Two masses with a damper on
   *one* DOF only (so $C$ is **not** proportional) — the canonical non-classical textbook problem
   (e.g. Chopra / Veletsos–Ventura). The exact complex $\lambda_{1,2}$ and complex mode shapes are
   available in closed form (solve the $4\times4$ state-space by hand / symbolically). The
   projection (with $p=2$, i.e. $\Phi$ = the full real basis, so the projection is **exact, not
   approximate**) must reproduce both complex eigenpairs to $1\text{e-}8$ and exhibit the correct
   **non-zero phase spread** between the two DOFs. This is the test that *proves the non-classical
   path*, and the one a draft most easily gets wrong.
3. **Classical-damping degeneracy.** For any model with $C=\alpha_M M+\beta_K K$, every complex mode
   must collapse to a real mode (phase spread $<1\text{e-}8$) with
   $\zeta_\alpha=\alpha_M/2\omega_\alpha+\beta_K\omega_\alpha/2$ — guards against a non-classical
   bug masquerading as signal.

**Tier 2 — projected-matrix cross-check (same matrices, independent solver):**

4. Extract the very $\tilde M,\tilde C,\tilde K$ the C++ driver builds (dump via a debug response)
   and feed them to **`scipy.linalg.eig`** on the $2p$ state-space pencil; the C++ `dggev`
   eigenvalues/vectors must agree to $1\text{e-}10$. This isolates *projection correctness* from
   *QZ correctness*: if Tier-1 fails but Tier-2 passes, the bug is in $\tilde C$ assembly (Route
   A/B), not the eigensolver.

**Tier 3 — physical decay (end-to-end, no oracle matrices):**

5. **Base-isolated stick model.** A few-DOF superstructure stick on a soft, highly-damped
   isolation layer (lead-rubber bearing → linear spring + dashpot). Run `complexEigen` → get the
   isolation-mode $\zeta$ and $\omega_d$. **Cross-check by direct time-history:** release the
   structure from a modal initial condition, record a free-DOF response, fit the **logarithmic
   decrement** to extract $\zeta_{\text{TH}}=\delta/\sqrt{4\pi^2+\delta^2}$ and the period
   $\to\omega_{d,\text{TH}}$. Agreement to ~2–5% validates the whole chain against the actual
   integrated dynamics (the same `getDamp()`/$C$ the transient feels).
6. **Supplemental-damper frame.** A 2–3 story frame with a fluid-viscous damper on one story; check
   the added modal damping matches a hand calculation of the damper's contribution and that the
   damped story-drift mode has the expected phase lag across the damped story.

**Quirk regression (folds in [[project_damping_channels]] / `LEDGER_quirks`):**

7. `zeroLength` `-doRayleigh 0` with a `Viscous` material: the complex eigen must return the
   **dashpot's** $\zeta$ (material branch), and a stiffness-Rayleigh-only variant of the same model
   must return **zero** added damping from that element — proving the assembler honors the element's
   own damping switches rather than re-deriving $C$.

---

## 9. Risk register

> [!question] **R1 — Assembled-$C$ availability is the make-or-break dependency.**
> The whole non-classical case (Route B) hinges on every damping source faithfully returning its
> contribution through `getDamp()`/`addCtoTang`. Some elements *don't* (the `-doRayleigh` default-0
> family) or return damping only in `getResistingForceIncInertia` rather than `getDamp`. *Mitigation:*
> v1 explicitly scopes $C$ = "whatever `getDamp()` returns" (= exactly the transient's $C$), tests
> the known quirky paths (§8.7), and documents the contract loudly. Any element whose damping is not
> in `getDamp()` is out of scope and flagged in the guide.

> [!question] **R2 — Completeness: not enough real modes $p$.**
> The projection is only as good as the span of $\Phi$. If the damped response excites modes above
> the retained band, $\tilde C$ misses inter-modal coupling into the truncated space and $\zeta_k$
> is biased. *Mitigation:* reuse the **effective-mass participation check** `modalProperties` already
> computes (`DomainModalProperties.cpp:646–668`) to warn when cumulative participation $<90\%$;
> document that $p$ should over-cover the frequency band of interest; the 2-DOF oracle (§8.2) uses
> $p=N$ so the projection is *exact* there (decouples completeness error from algorithm error).

> [!question] **R3 — QZ conditioning of the $2p\times2p$ pencil.**
> Heavily-damped or near-defective (close-to-critically-damped) modes can make the nonsymmetric
> pencil ill-conditioned; `dggev` returns eigenvalues as $(\alpha,\beta)$ to expose
> infinite/indeterminate roots. *Mitigation:* prefer the no-$\tilde M^{-1}$ generalized form (§4.4);
> flag $|\beta_k|<\varepsilon\|\tilde{\mathcal B}\|$ as rigid/infinite; report `dggev` INFO and any
> failed-to-converge eigenvalues instead of silently returning garbage; conjugate-pair residual
> check $\|(\lambda^2\tilde M+\lambda\tilde C+\tilde K)z\|$ as a per-mode quality metric.

> [!question] **R4 — Non-orthogonality of complex modes.**
> Complex modes are **not** $M$-orthogonal (only the $2N$ state-space modes are
> $\tilde{\mathcal B}$-biorthogonal). Naively reusing real-mode `modalProperties` formulas (which
> assume $\Phi^TM\Phi=I$) on complex modes is wrong. *Mitigation:* the `-complex` participation uses
> the **state-space biorthogonality** (left/right eigenvectors of the pencil) or the explicit complex
> $\psi_k^TM\psi_k$ normalization — never the real-mode shortcut; documented in the guide and gated by
> the classical-degeneracy test (§8.3).

> [!question] **R5 — Mass-normalization assumption on $\Phi$.**
> If `eigen` returned modes normalized to unit largest entry, $\tilde M\neq I$ and $\tilde K$ is not
> $\mathrm{diag}(\omega^2)$. *Mitigation:* internally re-normalize via $\phi^TM\phi$ (the value
> `modalProperties` already forms) before projecting; never assume the incoming normalization.

Plus standard items: LAPACK `dggev` is already linked (no new dep, low ABI risk); serial-only by
design (parallel is ADR 43); backwards-compatible (purely additive command + node field + optional
`modalProperties` flag — existing models untouched).

---

## 10. Ledger / header / PR obligations

Per `CLAUDE.md` (build-control is part of the work, same PR as the change):

- **`LEDGER_implementations.md`** — add a row when P1 lands: feature **`LadrunoComplexEigen`**,
  kind *analysis/solver*, **classTag 33019**, files (`SRC/analysis/.../LadrunoComplexEigen.{h,cpp}`,
  `LadrunoDampingAssembler.{h,cpp}`, `ComplexEigenIntegrator.{h,cpp}`, the `dggev`/`zggev` wrapper,
  `tests/test_ladrunoComplexEigen.py`), status (RESERVED → active at P1), PR. Flip RESERVED→active
  only when code merges (the tag enters `SRC/classTags.h` then).
- **`LEDGER_vanilla_files.md`** — every touched upstream file gets a row + an in-source
  `// Ladruno ...` marker: `analysis.cpp` (parser branch), `OpenSeesCommands.cpp` (py binding +
  legacy parser), `modal.cpp` / `commands.cpp` (`-complex` flag), `Node.{h,cpp}`
  (`setComplexEigenvector`/`getComplexEigenvectors`), and the recorder registration. Keep the vanilla
  footprint minimal (most logic lives in the new fork-authored files).
- **`LEDGER_quirks.md`** — record findings as they surface: the `getDamp()`/`-doRayleigh` damping-
  contract subtlety (what the assembled $C$ does and does *not* include); the complex-mode
  non-orthogonality trap (R4); any `dggev` $(\alpha,\beta)$ infinite-eigenvalue handling lesson.
- **Banner** — add a line to `Ladruno_scripts/banner_features.txt` (e.g.
  `complex/state-space modal (non-classical damping) — eigen -damped`), then
  `python Ladruno_scripts/patch_banner.py` and rebuild
  `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`. Every `shipped` ledger row needs a
  matching banner line.
- **Header stamp** — run `Ladruno_scripts/stamp_headers.py` (add the new files to its GLOBS) so each
  new source file carries the LADRUNO ASCII-art four-author header ([[feedback_always_stamp_header]]);
  vanilla edits are exempt. `--check` in CI.
- **PRs** — base every PR on **`ladruno`** (the default branch), one logical PR per phase, never
  `--base` a stacked branch (the [[feedback_stranded_commits_after_automerge|stranded-commit]] /
  hand-stacking pitfalls). Verify `gh pr view <n> --json state` == OPEN before stacking a follow-up.

---

## 11. Cross-references

- **Parent / family overview:** [[modal_gap_study/00_SYNTHESIS]] — §2a + §5 ADR-A define this as the
  #1, low-core-risk quick win; this ADR is its realization.
- **Primary theory:** [[modal_gap_study/02_abaqus_theory]] GAP 1 — Abaqus `*COMPLEX FREQUENCY`
  modal-projection QZ (TG §2.5.1): the state-space linearization (§4.2), the projection (§4.3), the
  recovery formulas (§4.5), and the "must build: $C$ accessor + dense QZ + complex post + new command"
  porting note are taken directly from here.
- **Secondary theory:** [[modal_gap_study/04_lsdyna_theory]] §2 — LS-DYNA `LCPACK=3` → ARPACK
  nonsymmetric quadratic (first-order state space). Note their complex eigen is **SMP-only** even in
  LS-DYNA, reinforcing the §3 decision to keep v1 serial.
- **Ground truth:** [[modal_gap_study/01_opensees_current_state]] — the file:line audit behind §6,
  including the definitive "no assembled $C$, no `addC`, no `formEleTangC`" finding.
- **Sibling ADRs:**
  - [[42_ladruno_buckling_adr]] — prestressed modal + linear buckling; composes with this once both
    land (complex modes about a prestressed state).
  - [[43_ladruno_feast_eigensolver_adr]] — robust + parallel eigensolver. **ADR 43's complex FEAST
    contours will re-host THIS complex case at scale** (full-$N$ direct complex eigen, distributed),
    making ADR 40 the serial/projection proof-of-physics that FEAST later generalizes.
  - [[44_ladruno_frequency_domain_adr]] — FRF/SSD/random; **consumes the mode-by-mode $\zeta_k$ this
    ADR produces** instead of a blanket damping value, and shares the $\Phi$ + reduced-matrix machinery.
- **Fork context:** [[project_damping_channels]] (the six damping channels + the `-doRayleigh`
  default-0 quirk that the $C$-assembler must honor); [[project_absorbing_pml_guide]] (SSI / radiation
  damping = a motivating non-classical $C$ source).

> [!info] One open question flagged for review
> **Should the assembled-$C$ contract be "exactly `getDamp()`" (v1 decision) or should we attempt to
> capture damping sources that live outside `getDamp()` (e.g. integrator-level numerical damping,
> modal `modalDamping`, or elements that only fold damping into `getResistingForceIncInertia`)?**
> v1 deliberately scopes $C$ to `getDamp()` for a clean, transient-consistent definition — but a user
> who set damping via `modalDamping` or HHT $\alpha$ would get a complex eigen that *ignores* it,
> which could surprise. The proposed resolution (document loudly + warn if `modalDamping` is active
> while `complexEigen` runs) needs sign-off, since the alternative (synthesizing a $C$ from those
> channels) materially widens scope toward ADR 44 territory.
