# Abaqus Modal / Eigenvalue Algorithms — Formulation-Level Extraction

**Purpose.** Extract Abaqus's modal, eigenvalue, buckling, frequency-domain, and
large-scale/parallel modal algorithms at *formulation* level to inform a new
OpenSees ("Ladruno" fork) feature. This is a theory-mining dossier, not an
implementation.

**Source.** Abaqus 2016 Theory Guide, via the locally-installed
`abaqus-theory-procedures` skill (synthesized + section-cited reference files;
the manual PDFs are not on this machine). Primary reference files:
`modal-dynamics.md` (TG §2.5.1–§2.5.8), `buckling-riks.md` (TG §2.3.1–§2.3.2),
`dynamics.md` (TG §2.4.1–§2.4.5), `substructuring-submodeling-fracture.md`
(TG §2.14.1, §2.15.1). Citations carry through as `(TG §2.x.y)`. Where the skill
flags an externally-attached name (e.g. Craig–Bampton) that the manual does not
itself use, that caveat is preserved.

OpenSees correspondences are drawn from the skill's own Abaqus↔OpenSees crosswalk
(SKILL.md) plus the fork's machinery (serial ARPACK `eigen`, `EigenSOE`,
corotational/PDelta geometric stiffness, `modalProperties`).

---

## Architectural keystone: the linear-perturbation philosophy (TG §2.1.1)

Every Abaqus modal/buckling result rides on one architectural decision that
OpenSees does **not** share by default: Abaqus splits the load history into
**steps**, and **state (stress, strain, temperature) carries forward across
every step**. All linear analysis is treated as a **linear perturbation about
the current base state** (TG §2.1.1). Consequences that matter for every gap
area below:

- A FREQUENCY (modal) or BUCKLE step executed *after* a geometrically nonlinear
  static step **automatically inherits the preload** — i.e. **stress stiffening
  is free**. There is no separate "rebuild K with Kg" plumbing; the tangent at
  the base state already contains the initial-stress (geometric) term.
- The tangent stiffness assembled at the base state splits into three physically
  distinct contributions (TG §2.1.1, Eq. 2.1.1-4):

$$K^{NM} = \underbrace{K^{NM}_{\text{mat}}}_{\text{small-disp. material}}
        + \underbrace{K^{NM}_{\text{geom}}}_{\text{initial-stress / geometric}}
        + \underbrace{K^{NM}_{\text{load}}}_{\text{load stiffness}}$$

  - **Material** stiffness from $d\sigma = D{:}d\varepsilon$.
  - **Initial-stress (geometric)** stiffness from the current stress acting
    through the change of the strain-displacement operator $\beta^N$ — the
    source of stress stiffening.
  - **Load** stiffness, present only for configuration-dependent (follower,
    pressure, centrifugal) loads.

This is the conceptual frame that unifies prestressed modal, linear buckling, and
the substructure "linear perturbation about a self-equilibrated make-state."

---

## GAP 1 — COMPLEX / STATE-SPACE MODAL (non-classical damping) [#1 target]

### Problem statement

For systems with **non-classically-distributed damping** (damping that is *not*
simultaneously diagonalizable with $M$ and $K$ — e.g. localized dashpots,
contact/joint damping, dissimilar-material composite damping), the real undamped
modes do **not** decouple the equations. The damped free-vibration problem is the
quadratic eigenproblem

$$\big(\lambda^2 M + \lambda C + K\big)\,\psi = 0,$$

whose eigenvalues $\lambda$ and eigenvectors $\psi$ are **complex**. Solving this
$N$-dimensional quadratic problem directly is expensive and nonsymmetric.

### Abaqus's algorithm (`*COMPLEX FREQUENCY`, TG §2.5.1)

Abaqus does **not** solve the full quadratic problem. It uses a
**subspace-projection (reduced complex eigenproblem)** approach built on the
already-extracted **real undamped eigenvectors**:

**Step 0 — prerequisite real modal basis.** A prior FREQUENCY step has extracted
$p$ real, mass-normalized undamped eigenvectors $\phi_\alpha$ from the
symmetrized $(K - \omega^2 M)\phi = 0$ (damping neglected, $K$ symmetric — see
Gap 3). Collect them as columns of $\Phi = [\phi_1, \dots, \phi_p]$, with
$\Phi^T M \Phi = I$ (or $= \mathrm{diag}(m_\alpha)$) and
$\Phi^T K \Phi = \mathrm{diag}(\omega_\alpha^2 m_\alpha)$.

**Step 1 — project $M$, $C$, $K$ onto the real undamped eigenvectors.** Form the
small ($p \times p$) reduced matrices:

$$\tilde M = \Phi^T M \Phi, \qquad
  \tilde C = \Phi^T C \Phi, \qquad
  \tilde K = \Phi^T K \Phi.$$

With mass-normalized modes, $\tilde M = I$ and $\tilde K =
\mathrm{diag}(\omega_\alpha^2)$. The crucial term is the **reduced damping
matrix** $\tilde C$, which is in general **full (non-diagonal)** precisely
*because* the damping is non-classical — the off-diagonal entries
$\tilde C_{\alpha\beta} = \phi_\alpha^T C \phi_\beta$ ($\alpha \neq \beta$) are
the inter-modal damping coupling that classical modal analysis throws away.

**Step 2 — solve the small complex / QZ generalized nonsymmetric problem.** The
reduced quadratic eigenproblem

$$\big(\lambda^2 \tilde M + \lambda \tilde C + \tilde K\big)\, z = 0$$

is linearized to a $2p \times 2p$ **generalized nonsymmetric** (first-order /
state-space) form and solved with a dense QZ-type generalized nonsymmetric
eigensolver (TG §2.5.1). One standard state-space linearization is

$$\underbrace{\begin{bmatrix} \tilde C & \tilde M \\ \tilde M & 0 \end{bmatrix}}_{\mathcal A}
  \begin{Bmatrix} z \\ \lambda z \end{Bmatrix}
  = -\lambda
  \underbrace{\begin{bmatrix} \tilde K & 0 \\ 0 & -\tilde M \end{bmatrix}}_{\mathcal B}
  \begin{Bmatrix} z \\ \lambda z \end{Bmatrix},$$

i.e. a generalized eigenproblem $\mathcal A x = -\lambda \mathcal B x$ of size
$2p$. Because $p \ll N$, this dense nonsymmetric solve is cheap.

**Step 3 — recover complex modes, frequencies, damping ratios.** Each reduced
eigenvector $z_k$ maps back to a **physical complex mode**
$\psi_k = \Phi\, z_k$ (a complex combination of the real undamped modes). Each
eigenvalue is complex,

$$\lambda_k = -\zeta_k \omega_{0,k} \;\pm\; i\,\omega_{d,k},$$

from which the engineering quantities are read directly:

$$\omega_{0,k} = |\lambda_k| = \sqrt{(\mathrm{Re}\,\lambda_k)^2 + (\mathrm{Im}\,\lambda_k)^2}
  \quad\text{(undamped natural circular frequency)},$$
$$\omega_{d,k} = \mathrm{Im}\,\lambda_k \quad\text{(damped frequency)},$$
$$\zeta_k = -\frac{\mathrm{Re}\,\lambda_k}{|\lambda_k|}
  = -\frac{\mathrm{Re}\,\lambda_k}{\sqrt{(\mathrm{Re}\,\lambda_k)^2 + (\mathrm{Im}\,\lambda_k)^2}}
  \quad\text{(damping ratio)}.$$

Complex eigenvalues come in conjugate pairs for underdamped modes. A complex mode
shape has a **phase distribution** (entries do not all peak simultaneously) —
this is the qualitative signature of non-classical damping that real-modes
analysis cannot represent. Overdamped modes appear as real (non-oscillatory)
$\lambda_k$.

### Why the projection is exact-enough

The subspace spanned by the $p$ real undamped modes captures the damped response
to the extent the retained modes capture the dynamics; the method is accurate
when the damping perturbs but does not radically reshape the low modes. It trades
an $N$-dimensional nonsymmetric quadratic solve for a $p$-dimensional one — the
same cost philosophy as mode-based steady-state dynamics (Gap 4).

### OpenSees porting note (Gap 1)

- **Reusable.** The hard prerequisite — the real undamped modal basis $\Phi$ with
  $\Phi^T M\Phi$, $\Phi^T K\Phi$ — is *exactly* what serial ARPACK `eigen`
  (`-genBandArpack` / `-fullGenLapack`) already produces, and `modalProperties`
  already forms $\phi^T M \phi$. The projection $\Phi^T (\cdot) \Phi$ is three
  dense mat-mats once you can extract $\Phi$, $M$, $K$, $C$ as accessible
  operators.
- **Must build.** (a) Assembly/extraction of the **damping matrix $C$** as a
  first-class operand for projection — OpenSees applies damping (Rayleigh,
  element, modal) *inside* integrators rather than exposing an assembled $C$;
  this is the main plumbing gap. (b) A **dense nonsymmetric/QZ generalized
  eigensolver** for the $2p \times 2p$ reduced problem — LAPACK `dggev` (or
  `zggev` for complex) is the natural backend; OpenSees's eigen path is built for
  *symmetric* generalized problems (ARPACK shift-invert), so the small
  nonsymmetric solve is genuinely new. (c) Complex post-processing
  ($\omega_0, \omega_d, \zeta$, complex mode export). (d) A new analysis command,
  e.g. `complexEigen`, consuming a prior `eigen` result.
- **Strategy.** Because the heavy lifting is a *small* dense problem, this is a
  high-value/low-core-risk feature: no change to the assembly core, no new
  EigenSOE for the big problem — just add a $C$ accessor + a small QZ wrapper +
  a projection driver. This is the recommended #1 target.

---

## GAP 2 — PRESTRESSED MODAL + LINEAR (EIGENVALUE) BUCKLING

### Shared architecture: perturbation about a base state

Both prestressed modal and eigenvalue buckling exploit the same
linear-perturbation machinery (TG §2.1.1). Because state carries across steps,
the geometric (initial-stress) stiffness $K_{\text{geom}}$ produced by the base
state is *already in the tangent* when the perturbation step runs.

### Prestressed modal

After a (possibly nonlinear) static preload step, the FREQUENCY step solves the
symmetrized eigenproblem on the **base-state tangent**:

$$(K - \omega^2 M)\,\phi = 0, \qquad
  K = K_{\text{mat}} + K_{\text{geom}}(\sigma_{\text{base}}) + K_{\text{load}}.$$

Because $K$ now carries the initial-stress term, $K$ **need not be positive
definite** (a member in compression softens; a stressed cable stiffens) — the
eigensolver (Gap 3) is built to tolerate this. This is "stress stiffening for
free": the same Lanczos/subspace machinery, just on the preloaded tangent.

### Eigenvalue buckling (`BUCKLE` procedure, TG §2.3.1)

Estimates **elastic critical loads** of stiff, near-linear-prebuckling
structures. From an arbitrary preloaded **base state** carrying **dead/base
loads** $P$, a pattern of perturbation **"live" loads** $Q$ is applied, and
Abaqus solves the generalized eigenproblem for **load multipliers (buckling
factors)** $\lambda_i$ and mode shapes $\phi_i$:

$$\boxed{(K_0 + \lambda_i\, K_\Delta)\,\phi_i = 0}$$

**Assembly of the two matrices:**

- $K_0$ — the **base-state stiffness**: hypoelastic material tangent +
  initial-stress stiffness + load stiffness evaluated at the base state,

$$K_0 = K_{\text{hypoelastic}} + K_{\text{initial-stress}}(\sigma_{\text{base}}) + K_{\text{load}}.$$

- $K_\Delta$ — the **differential (geometric/incremental) stiffness**: the
  initial-stress stiffness produced by the **perturbation stresses** (the linear
  stress response to the live load $Q$) **plus the perturbation load stiffness**.
  Operationally: apply $Q$ as a linear perturbation, get the incremental stress
  field $\Delta\sigma(Q)$, and build $K_\Delta = K_{\text{geom}}(\Delta\sigma(Q))$
  (+ load-stiffness part).

**Critical load.** The estimated buckling load in mode $i$ is

$$P_{cr} = P_{\text{base}} + \lambda_i\, Q,$$

so $\lambda_i$ is the scalar by which the *live* load pattern must be multiplied
(on top of the dead load) to buckle. Several modes are extracted; the **lowest
$\lambda_i$** is usually the design quantity.

**Limitations / fidelity (TG §2.3.1).** Constitutive law may be
elastic/hypoelastic/hyperelastic, but **plasticity and rate effects are
ignored**. **Nonconservative (follower) load** contributions to the load
stiffness are nonsymmetric and are **symmetrized** because Abaqus's eigensolvers
require symmetric matrices. Eigenvalue buckling is **reliable only when
prebuckling deformation is small and the live-load response is essentially
linear up to bifurcation**; for significant prebuckling nonlinearity (large
deflection, plasticity, imperfection-sensitive shells) it over-/under-predicts
collapse — the recommended fallback is **nonlinear Riks** (`STATIC, RIKS`,
TG §2.3.2) seeded with imperfections from the buckling modes.

**Contrast: subspace iteration vs the buckling eigenproblem.** The *eigensolver*
used for $(K_0 + \lambda K_\Delta)\phi = 0$ is the same family as for FREQUENCY
(Lanczos default, or subspace iteration). Subspace iteration projects the pair
onto a fixed-dimension subspace $X$,
$K_0^* = X^T K_0 X$, $K_\Delta^* = X^T K_\Delta X$, solves the small problem by
Householder + shifted QR, maps back, and re-iterates (Gap 3) — convergence per
mode scales with the eigenvalue-ratio gap.

### OpenSees porting note (Gap 2)

- **Reusable.** OpenSees can assemble **geometric stiffness** via the
  **corotational** and **P-Delta** geometric transformations and via element
  `Kg` contributions; the fork already leans on corotational $K_g$. The `eigen`
  command already solves a symmetric generalized eigenproblem $(K - \omega^2
  M)\phi = 0$, which is structurally identical to $(K_0 + \lambda K_\Delta)\phi =
  0$ if you feed $K_\Delta$ where $M$ goes.
- **Must build.** There is **no single `BUCKLE` procedure** (no `LadrunoBuckle`)
  that packages base + differential stiffness in one step. To replicate Abaqus
  faithfully you need: (a) a driver that runs the dead-load base state, captures
  $K_0$; (b) a **linear perturbation** of the live load $Q$ to get
  $\Delta\sigma(Q)$ and assemble $K_\Delta = K_{\text{geom}}(\Delta\sigma)$
  *separately* from $K_0$ (today OpenSees workflows hand-assemble a single $K_g$
  and call `eigen`, conflating the two); (c) feed $(K_0, K_\Delta)$ to the
  symmetric generalized eigensolver and report $\lambda_i$ and $P_{cr} =
  P_{\text{base}} + \lambda_i Q$. (d) Optional follower-load symmetrization.
- **Prestressed modal** is closer to free: capture the preloaded tangent (with
  $K_{\text{geom}}$ already in it) and call `eigen` — the gap is mostly *workflow
  packaging* (a "perturbation step that inherits state") rather than new math.
- **Strategy.** The big missing piece is the **step/base-state architecture**
  (state carried into a perturbation that assembles a clean $K_\Delta$). A
  `LadrunoBuckle` analysis object that orchestrates dead → live-perturbation →
  generalized eigen is a natural, self-contained deliverable.

---

## GAP 3 — EIGENSOLVER INTERNALS + STURM GUARANTEE

The natural-frequency eigenproblem is solved after **symmetrizing**: damping $C$
is dropped and $K$ taken symmetric (TG §2.5.1):

$$(K - \omega^2 M)\,\phi = 0.$$

$K$ may include stress stiffening (Gap 2), so it **need not be positive
definite**. Abaqus offers two symmetric eigensolvers.

### Shifted block Lanczos (default; Grimes–Lewis–Simon), TG §2.5.1

**Spectral shift-invert transformation.** A shift $\sigma$ is chosen and the
problem is recast so the wanted (low) modes become the *dominant* (large) modes of
an operator amenable to Krylov iteration:

$$\boxed{(K - \sigma M)^{-1} M\,\phi = \theta\,\phi}, \qquad
  \boxed{\omega^2 = \sigma + \frac{1}{\theta}}.$$

Eigenvalues near the shift $\sigma$ map to large $\theta$, so Lanczos converges to
them first; sweeping the shift across the spectrum extracts bands of modes. Each
operator application requires a **factorization of the shifted matrix
$(K - \sigma M)$** (one $LDL^T$/Cholesky per shift) and back-substitutions.

**Block vectors for coincident / clustered modes.** Lanczos is run in **"runs"**
(each at a fixed shift) of **"steps"**, building a growing Krylov subspace. The
algorithm uses **block** (multiple simultaneous) Lanczos vectors so that
**coincident or closely-spaced eigenvalues** (common in symmetric structures) are
captured robustly — a single-vector Lanczos can miss multiplicity.

**Partial reorthogonalization.** To control loss of orthogonality among Lanczos
vectors (the classic Lanczos failure mode) without the cost of full
reorthogonalization, the method uses **partial reorthogonalization** — reorthog
only when monitored loss of orthogonality crosses a threshold.

**Sturm-sequence mode-count guarantee.** The decisive correctness feature: the
**number of negative pivots** in the $LDL^T$/Cholesky factorization of the
shifted matrix **equals the number of eigenvalues below the shift**:

$$\#\{\text{eigenvalues} < \sigma\} = \#\{\text{negative pivots of } (K - \sigma M)\}.$$

By comparing the negative-pivot count at two shifts bracketing a band against the
number of modes actually converged in that band, Abaqus **guarantees no modes are
missed** — a guarantee Arnoldi/Lanczos-without-Sturm (e.g. plain ARPACK) does not
provide.

### Subspace iteration (alternative), TG §2.5.1

A simultaneous inverse-power sweep on a **fixed-dimension subspace** (default
dimension $\min(2p,\, p+8)$ for $p$ requested modes). At each sweep $K$ and $M$
are projected onto the current subspace $X$,

$$K^* = X^T K X, \qquad M^* = X^T M X \;\Rightarrow\; (K^* - \omega^2 M^*)\,\psi = 0,$$

the small reduced problem is solved by **Householder tridiagonalization + QR with
shifts**, the Ritz vectors map back, and the sweep repeats. Convergence rate for
mode $i$ scales with the spectral gap $\lambda_i / \lambda_{m+1}$ (ratio of the
mode's eigenvalue to the first un-tracked eigenvalue). Cheaper to implement,
slower than Lanczos for many modes; no Sturm guarantee by itself.

### Modal variables (post-FREQUENCY), TG §2.5.2

Per mode $\alpha$, with $T_i$ the rigid-body vector for global direction $i$
(3 translations + 3 rotations about a **center of rotation**):

$$m_\alpha = \phi_\alpha^T M \phi_\alpha \qquad\text{(generalized / modal mass)}$$
$$\Gamma_{\alpha i} = \frac{\phi_\alpha^T M\, T_i}{m_\alpha} \qquad\text{(participation factor, dir. } i)$$
$$m_{\text{eff},\alpha i} = \frac{(\phi_\alpha^T M\, T_i)^2}{m_\alpha}, \qquad
  \sum_\alpha m_{\text{eff},\alpha i} \to m_{\text{total}}.$$

Eigenvectors may be normalized to **unit largest entry** or to **unit generalized
mass** — the choice does not affect physical results. The **effective mass**
measures each mode's contribution to total dynamic response in direction $i$;
summed over all extracted modes it should recover the total unrestrained mass. A
large **shortfall** warns that important modes were not extracted — the standard
**seismic mass-participation check** (target often ≥ 90%). **Composite modal
damping** assembles per-material critical-damping fractions into each mode,
mass-weighted by each material's mass contribution.

### OpenSees porting note (Gap 3)

- **Reusable.** OpenSees `eigen -genBandArpack` (shift-invert ARPACK,
  Arnoldi/Lanczos) **already implements the spectral transformation**
  $(K - \sigma M)^{-1} M \phi = \theta \phi$ — the core Lanczos machinery and the
  symmetric generalized solve are present. `EigenSOE`/`EigenSolver` is the
  pluggable abstraction. `modalProperties -print` **already reports** generalized
  mass, participation factors, and effective mass per the §2.5.2 formulas, and
  does the seismic mass check.
- **Must build / the real gap.** The **Sturm-sequence mode-count guarantee is
  not exposed** in OpenSees. To add it you need access to the **inertia
  (negative-pivot count) of the factorized shifted matrix $(K - \sigma M)$** from
  the linear-solver layer (e.g. the $LDL^T$/sparse factor's diagonal signs), then
  a driver that brackets a frequency band with two shifts and verifies the
  converged-mode count against the pivot-count difference. This is a genuinely
  valuable add: it converts "ARPACK returned $p$ modes" into "*these are
  provably all the modes in $[\omega_a, \omega_b]$*."
- **Optional.** Block-Lanczos handling of exact multiplicity is largely covered by
  ARPACK's implicitly-restarted Arnoldi in practice; subspace iteration is not
  worth porting given ARPACK. The center-of-rotation handling for rotational
  $T_i$ is already in `modalProperties`.

---

## GAP 4 — FREQUENCY DOMAIN

### Modal dynamic time history (exact piecewise-linear SDOF), TG §2.5.3, §2.5.5

Modal superposition projects onto a small set of modes; orthogonality
**decouples** the system into single-DOF oscillators:

$$\ddot q_\alpha + 2\xi_\alpha \omega_\alpha \dot q_\alpha + \omega_\alpha^2 q_\alpha = f_\alpha(t),
  \qquad
  \xi_\alpha = \frac{c_\alpha}{2 m_\alpha \omega_\alpha}, \quad
  f_\alpha = \frac{\phi_\alpha^T P}{m_\alpha}.$$

Assuming the load varies **piecewise-linearly within each increment**, each modal
ODE is integrated by an **exact closed-form solution** (particular + homogeneous),
with three damping cases (under-/critically-/over-damped) plus a special
**rigid-body** case ($\omega = 0$; only mass-proportional Rayleigh damping is
admissible). If the projected damping is **non-diagonal** (general coupling,
i.e. the Gap-1 situation in the time domain), a fully coupled form is solved with
a single one-time factorization. Physical responses recover by summing modal
contributions $u = \sum_\alpha \phi_\alpha q_\alpha$, $\sigma = \sum_\alpha
\sigma_\alpha q_\alpha$. **Base motion** is converted to an acceleration history
(relative-to-base formulation); total values add the base motion back. Structural
(complex-stiffness) damping is converted to equivalent viscous damping for the
time domain. Initial conditions are projected onto the modal basis (exact only if
\#modes = \#DOF).

### Steady-state dynamics (SSD), TG §2.5.7

Predicts harmonic response vs excitation frequency $\Omega$. **Mode-based SSD**
projects onto eigenmodes; per mode the complex harmonic balance gives a **complex
transfer function**:

$$-\Omega^2 m_\alpha q_\alpha + i\Omega c_\alpha q_\alpha + \omega_\alpha^2 m_\alpha q_\alpha
  = f_\alpha(\Omega),$$
$$\boxed{H_\alpha(\Omega) = \frac{1}{m_\alpha(\omega_\alpha^2 - \Omega^2) + i\,\Omega\, c_\alpha}}.$$

Physical peak amplitudes sum the modal amplitudes over a (biased) **frequency
sweep** clustered near resonances. Three flavors:

- **Mode-based SSD** — project onto extracted undamped modes (cheapest).
- **Subspace SSD** — project onto selected undamped modes but solve the *coupled*
  complex reduced system (handles modest frequency-dependence / non-diagonal
  damping).
- **Direct SSD** — solve the full complex harmonic equations
  $(-\Omega^2 M + i\Omega C + K)\,u = f(\Omega)$ **without modes** (needed for
  nonsymmetric systems or frequency-dependent material properties).

Modal, structural, and Rayleigh damping are all admissible; **structural damping**
(complex stiffness) fits naturally in the frequency domain.

### Random response, TG §2.5.8

For a statistically-defined (stationary, ergodic) excitation specified by an
input **cross-spectral density (CSD)** matrix $S_{in}(\Omega)$, the modal/transfer
functions propagate input → output:

$$\boxed{S_{out}(\Omega) = H^*(\Omega)\, S_{in}(\Omega)\, H(\Omega)^T}$$

(input CSD → output **power spectral density** PSD), and the **RMS** of any
response is the square-root of the integral of its PSD:

$$\sigma_x = \sqrt{\int S_{xx}(\Omega)\, d\Omega}.$$

The mean of each dynamic variable is zero (response about static equilibrium), so
variance = mean-square and outputs are reported as **RMS**. Related concepts:
autocorrelation; white / wide-band / narrow-band noise.

### Damping options (combinable additively), TG §2.5.4

1. **Direct critical-damping fractions** $\xi_\alpha$ per mode (1–10% typical) —
   a mathematical device tied to the eigenmodes; **not transferable** to
   nonlinear direct integration.
2. **Rayleigh** $C = \alpha_R M + \beta_R K$ — shares the undamped eigenvectors,
   so it maps exactly to modal fractions:
   $$\xi_\alpha = \frac{\alpha_R}{2\omega_\alpha} + \frac{\beta_R\,\omega_\alpha}{2}.$$
   Mass term damps low modes; stiffness term damps high modes.
   ($c_{cr} = 2\sqrt{mk} = 2m\omega$.)
3. **Composite modal damping** — per-material critical fraction mass-weighted into
   each mode.
4. **Structural (complex-stiffness) damping** — damping force
   $F_s = i\, s\, F_{\text{internal}}$, 90° out of phase with the elastic force;
   valid only for **sinusoidal** excitation, hence usable **only in SSD and random
   response**.

### Response spectrum (peak estimate), TG §2.5.6

Estimates **peak** response to base motion using the lowest modes. Spectrum
$S(\omega,\xi)$; for light damping $S_v \approx \omega S_d$, $S_a \approx \omega^2
S_d$. Modal peak under excitation along direction $k$:

$$q_{\alpha,\text{peak}} = a\, S(\omega_\alpha, \xi_\alpha)\, \sum_j n_{jk}\,\Gamma_{\alpha j},$$

with $a$ the excitation scale and $n_{jk}$ the direction cosines. Modal peaks are
non-simultaneous, so combine by **ABS** (most conservative), **SRSS** (well-
separated), **CQC** (Der Kiureghian, closely-spaced),
$R = \sqrt{\sum_\alpha \sum_\beta \rho_{\alpha\beta} R_\alpha R_\beta}$,
**Ten-Percent** (Reg. Guide 1.92), or **NRL** (regulatory).

### OpenSees porting note (Gap 4)

- **Reusable.** The modal basis + `modalProperties` (participation factors,
  $m_\alpha$) supply every ingredient for SSD/random/response-spectrum. OpenSees
  has `responseSpectrumAnalysis` (or scripted SRSS/CQC) and `rayleigh` /
  `modalDamping`. The Rayleigh↔modal-fraction conversion above is standard.
- **Must build.** OpenSees has **no native SSD, no random response, no native
  modal-dynamic superposition** (it is direct-integration-centric). All three are
  small dense post-processors on top of an existing `eigen` result:
  - **Mode-based SSD** = evaluate $H_\alpha(\Omega)$ over a frequency sweep and
    sum — a few hundred lines, no core change.
  - **Random response** = $S_{out} = H^* S_{in} H^T$ + numerical integration for
    RMS — likewise a post-processor.
  - **Modal transient** = exact piecewise-linear SDOF recurrence per mode — a
    clean, well-defined algorithm.
  - **Direct SSD** = needs the assembled complex operator
    $(-\Omega^2 M + i\Omega C + K)$ and a complex linear solve — this *does* touch
    the SOE/solver layer (complex system) and is the heaviest item; defer it.
  - **Structural (complex-stiffness) damping** is absent and only meaningful once
    a frequency-domain solver exists.
- **Strategy.** Mode-based SSD + random response + modal transient are a coherent,
  low-risk "frequency-domain toolkit" bundle riding entirely on `eigen` +
  `modalProperties`. Direct SSD and structural damping are a later, solver-touching
  phase.

---

## GAP 5 — LARGE-SCALE / PARALLEL (SUBSTRUCTURING, AMLS, SIM)

### Substructuring / superelements, TG §2.14.1

A **substructure (superelement)** is a group of elements whose **internal DOF
have been statically condensed out**, leaving only **retained DOF** that connect
to the rest of the model. Its response is a **linear perturbation about the
self-equilibrated state** at generation time (`SUBSTRUCTURE GENERATE`); reactions
at fixed retained DOF become a preload.

**Static (Guyan) condensation.** Partition stiffness into retained ($r$) and
internal ($i$) blocks and eliminate the internal block:

$$\boxed{K_{red} = K_{rr} - K_{ri}\, K_{ii}^{-1}\, K_{ir}}, \qquad
  P_{red} = P_r - K_{ri}\, K_{ii}^{-1}\, P_i.$$

**Dynamics → fixed-interface component modes (Craig–Bampton).** Static reduction
omits internal inertia, so for dynamics the basis is augmented with
**fixed-interface natural modes** $\Phi_i$ (generalized DOF $q$). The internal
displacement is the static (constraint-mode) part plus the modal part:

$$u_i = -K_{ii}^{-1} K_{ir}\, u_r + \Phi_i\, q,$$

the first term being the **constraint modes** $\Psi_c = -K_{ii}^{-1}K_{ir}$ and
$\Phi_i$ the **fixed-interface normal modes** (eigenvectors of $K_{ii} \phi =
\omega^2 M_{ii} \phi$). This is the standard **Craig–Bampton component-mode
synthesis** (the Theory Guide says "Guyan + generalized/normal modes" and does
**not** itself use the name "Craig–Bampton"). The reduction yields coupled reduced
mass/damping/stiffness over $a = (u_r, q)$:

$$M_{red}\,\ddot a + C_{red}\,\dot a + K_{red}\, a = P_{red}.$$

**Large-rotation substructures.** For geometric nonlinearity, a substructure
computes an equivalent **rigid-body rotation** from its retained nodes, subtracts
it to isolate the **strain-inducing displacements**, and rotates the reduced
stiffness/mass into the current frame; fixed-direction gravity is handled via
internally-generated unit-direction load cases scaled by the rotated gravity
direction (TG §2.14.1).

### AMLS (automated multi-level substructuring) and the SIM architecture

> Caveat: AMLS and the **SIM** architecture are Abaqus *eigenfrequency*
> capabilities for very large models; the skill's TG-cited reference files
> (§2.5.1, §2.14.1) document the underlying **Lanczos/subspace** eigensolvers and
> **Guyan + fixed-interface-mode** substructuring math but do not contain a
> dedicated AMLS section. The description below frames AMLS/SIM in terms of the
> formulations that *are* TG-documented (recursive substructuring + component
> modes); treat the AMLS-specific recursion as the standard published method, not
> a verbatim TG citation.

**AMLS** is the recursive, automated generalization of Craig–Bampton: the model's
DOF graph is **partitioned hierarchically** (a multi-level tree of substructures);
at each level the internal DOF are eliminated by Guyan condensation and the basis
augmented with a *truncated* set of fixed-interface modes, recursively merging
child substructures into parents. The result is a drastically reduced eigenproblem
solved once at the top, then back-transformed. AMLS trades a small, controlled
eigenvalue error for the ability to extract **thousands of modes on
multi-million-DOF models** far faster than global Lanczos — at the cost of
approximate (mode-truncated) accuracy, where shifted block Lanczos is exact.

**SIM** is Abaqus's internal high-performance **linear-algebra/eigen
architecture** (a sparse matrix + eigen infrastructure) that AMLS, Lanczos,
mode-based SSD, and complex-frequency are built on; it is what makes the
mode-based downstream procedures (Gaps 1, 4) cheap once the eigenbasis exists.

**Lanczos parallelism.** The shifted block Lanczos solver (Gap 3) parallelizes at
two levels: (a) the **factorization and back-substitution of the shifted matrix
$(K - \sigma M)$** use the parallel sparse direct solver (the dominant cost), and
(b) different **shifts / frequency bands** can be processed concurrently. Block
vectors also give BLAS-3 (matrix-matrix) density for cache/thread efficiency.
Domain decomposition feeds the parallel factorization.

### OpenSees porting note (Gap 5)

- **Reusable.** OpenSees has the static-condensation primitive in spirit
  (constraint handlers, `Subdomain`/`Actor` for parallel `OpenSeesSP`/`MP`), and
  the parallel sparse solvers (MUMPS) needed for the dominant factorization cost.
  The `eigen` path runs in the parallel build.
- **Must build (essentially all of it).** OpenSees has **no general
  superelement / Guyan / Craig–Bampton component-mode synthesis** — model
  reduction is done externally. A faithful port needs: (a) a **substructure
  generator** that partitions $r/i$ DOF, forms $K_{red} = K_{rr} -
  K_{ri}K_{ii}^{-1}K_{ir}$ and the constraint/fixed-interface modes, and emits a
  reusable reduced element; (b) the reduced $(M_{red}, C_{red}, K_{red})$ assembly
  for dynamics; (c) optional large-rotation corotational wrapper. **AMLS** (the
  recursive multilevel version) is a large research effort and should be scoped
  *after* a single-level Craig–Bampton superelement proves out. There is **no
  SIM-equivalent** unified eigen/linear-algebra layer; OpenSees relies on
  ARPACK + (MUMPS/UMFPACK) instead — adequate for single-level CMS, not for
  AMLS-scale recursion.
- **Strategy.** A **single-level Craig–Bampton superelement** (`LadrunoSuperElement`
  / substructure generate+use) is the realistic, high-value deliverable: it reuses
  the existing assembly + sparse factorization, requires no eigen-core change, and
  directly enables much larger modal models and component-level reuse. AMLS, SIM,
  and parallel multi-shift Lanczos are deliberately out of near-term scope.

---

## Consolidated OpenSees gap matrix

| Area | Abaqus algorithm (TG §) | OpenSees reuse | Must build | Risk/scope |
|---|---|---|---|---|
| **1. Complex modal** | Project $M,C,K$ on real modes → small QZ generalized nonsym. (§2.5.1) | `eigen` real basis + `modalProperties` ($\phi^TM\phi$) | $C$ accessor; $2p$ dense nonsym/QZ (`dggev`/`zggev`); complex post; `complexEigen` cmd | **Low core risk, high value — #1** |
| **2. Prestressed modal + buckling** | Perturbation about base state; $(K_0+\lambda K_\Delta)\phi=0$; $P_{cr}=P_{base}+\lambda Q$ (§2.1.1, §2.3.1) | corotational/PDelta $K_g$; symmetric generalized `eigen` | base-state→live-perturbation step arch.; clean $K_\Delta$ from $\Delta\sigma(Q)$; `LadrunoBuckle` | Medium (workflow > math) |
| **3. Eigensolver + Sturm** | Shift-invert block Lanczos $(K-\sigma M)^{-1}M\phi=\theta\phi$; Sturm neg-pivot count (§2.5.1) | `eigen -genBandArpack` (shift-invert ARPACK); `modalProperties` | **negative-pivot inertia from factor** → Sturm mode-count guarantee | Medium (needs solver-layer hook) |
| **4. Frequency domain** | $H_\alpha(\Omega)$ SSD; $S_{out}=H^*S_{in}H^T$ random; exact PWL modal transient (§2.5.3–8) | modal basis; `responseSpectrumAnalysis`; `rayleigh`/`modalDamping` | mode-based SSD + random + modal-transient post-processors; (direct SSD = complex solver, defer) | Low for mode-based; high for direct |
| **5. Large-scale / parallel** | Guyan + Craig–Bampton CMS; AMLS recursion; SIM; parallel multi-shift Lanczos (§2.14.1) | `Subdomain`/`SP`/`MP`; MUMPS; constraint handlers | single-level CMS superelement; (AMLS/SIM = out of scope) | High |

---

## Sources

All citations resolve to the Abaqus 2016 Theory Guide via the
`abaqus-theory-procedures` skill reference files:

- `references/modal-dynamics.md` — TG §2.5.1–§2.5.8 (eigensolvers, complex modes,
  modal variables, modal transient, damping, response spectrum, SSD, random).
- `references/buckling-riks.md` — TG §2.3.1–§2.3.2 (eigenvalue buckling, Riks).
- `references/dynamics.md` — TG §2.4.1–§2.4.5 (HHT, subspace dynamics, explicit).
- `references/substructuring-submodeling-fracture.md` — TG §2.14.1 (Guyan +
  fixed-interface CMS), §2.15.1.
- `SKILL.md` — TG §2.1.1 linear-perturbation architecture + Abaqus↔OpenSees /
  Ladruno crosswalk.

The manual PDFs are not present on this machine; equations above are restated from
the skill's synthesized, section-cited reference files. AMLS/SIM specifics (Gap 5)
are framed from the TG-documented substructuring + Lanczos formulations plus the
standard published method, and are flagged inline as not carrying a dedicated TG
section in the available references.
