---
title: "ADR 44 — Frequency-domain & modal-superposition response (FRF / SSD / random / LadrunoModalResponse): design spec"
project: Ladruno
type: ADR / design spec
status: draft — design only, NO code
priority: medium
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"          # §4 frequency domain — cheap post-processors on eigen
  - "[[modal_gap_study/02_abaqus_theory]]"      # PRIMARY: TG §2.5.3 PWL modal SDOF, §2.5.6 RS, §2.5.7-8 SSD/random
  - "[[modal_gap_study/04_lsdyna_theory]]"      # *FREQUENCY_DOMAIN_* family (FRF/SSD/random/RS), modal vs direct
  - "[[modal_gap_study/01_opensees_current_state]]" # ground truth: eigen + modalProperties + responseSpectrumAnalysis
  - "[[46_ladruno_complex_modal_adr]]"          # SIBLING (planned): non-classical damping → complex state-space modes
  - "[[42_ladruno_buckling_adr]]"               # SIBLING (planned): prestressed modal + linear buckling
  - "[[43_ladruno_feast_eigensolver_adr]]"      # SIBLING (planned): robust + parallel eigensolver (FEAST/complex-FEAST)
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
  - "[[LEDGER_vanilla_files]]"
tags: [adr, solver, dynamics, modal, frequency-domain, frf, ssd, random-vibration, psd, cqc, srss, nvh]
updated: 2026-06-22
---

# ADR 44 — Frequency-domain & modal-superposition response

> **Strategic role (load-bearing assessment — see [[modal_gap_study/00_SYNTHESIS]] §6).**
> **Deliverable layer — NOT load-bearing.** This ADR *consumes* the modal basis from ADRs 46/43;
> nothing downstream builds on it. It is the end-feature set (NVH, response spectrum, random
> vibration, fast modal transient) that monetizes the infrastructure. Genuinely useful, but **build
> it last — or only when a specific project needs it.**

**Status:** draft. **Design only — no code has landed.** classTag **33024 (RESERVED, not yet
built)** for the `LadrunoModalResponse` analysis/post-processor object. This ADR specifies the
modal-superposition *frequency-domain toolkit* that rides on the existing `eigen` +
`modalProperties` infrastructure. It is the fork's realization of [[modal_gap_study/00_SYNTHESIS]]
**ADR-D ("Frequency domain")** and is the cheap, low-core-risk member of the modal ADR family
(siblings: complex modes = ADR 46, buckling = ADR 42, parallel eigensolver = ADR 43).

> [!info] What this is, in one line
> Once `eigen` gives the mode shapes Φ and frequencies ω_α and `modalProperties` gives modal
> mass / participation, **FRF, steady-state dynamics (SSD), random response (PSD→RMS), exact
> modal-superposition transient, and native CQC/SRSS response-spectrum combination are all
> *small dense post-processors* on top of those results** — almost zero new core. They are
> dramatically cheaper than direct time integration for *linear* problems (NVH, equipment,
> wind, floor vibration, fast linear seismic) because the cost rides the one-time eigensolve.

---

## 1. Driver & goal

OpenSees today is **direct-integration-centric**. For a linear structure subjected to harmonic,
broadband, or earthquake excitation, the standard workflow is a full Newmark/HHT transient — even
though the response is fully determined by a handful of modes. Five large, real use-classes are
served far more cheaply by the frequency / modal-superposition domain:

1. **NVH / vibration serviceability** — find the response amplitude vs forcing frequency (FRF),
   identify resonances, check acceleration limits.
2. **Rotating equipment / machine foundations** — steady-state harmonic (SSD) response at and
   near operating speeds, including a frequency sweep biased toward resonances.
3. **Wind / turbulence / acoustic loading** — stationary random excitation specified as an input
   PSD/CSD; the engineering deliverable is **RMS response** and (optionally) peak / zero-crossing
   statistics.
4. **Floor-vibration / footfall** — combination of FRF + transient modal superposition.
5. **Fast linear seismic** — modal-superposition transient (exact piecewise-linear SDOF
   integration) replaces a costly nonlinear-grade direct run when the model is genuinely linear;
   and **native CQC/SRSS** finishes the response-spectrum workflow that
   `responseSpectrumAnalysis` deliberately leaves to the user.

The decisive property is **cost**: each of these is a dense operation on the `p`-dimensional modal
space (`p ≪ N`). The expensive step — the eigensolve — *already exists and is already run* for
`modalProperties`. So this whole toolkit is "free riding": ~zero assembly-core change, no new big
eigensolver, no SOE/solver-layer surgery. (Contrast the *direct* route in §3, which does need a
complex linear solver and is deferred.)

This is also a **parity gap**: Abaqus ships the entire stack (`*STEADY STATE DYNAMICS`,
`*RANDOM RESPONSE`, mode-based transient, response-spectrum with CQC/SRSS — TG §2.5.3–§2.5.8);
LS-DYNA ships `*FREQUENCY_DOMAIN_{FRF,SSD,RANDOM_VIBRATION,RESPONSE_SPECTRUM}`. OpenSees ships
*none* of FRF/SSD/random/native-CQC and has *no* native modal-superposition transient.

---

## 2. Decision summary

Build a family of **modal-superposition post-processors** that consume a prior `eigen` result and
the `DomainModalProperties` object, packaged behind a single fork-authored
`LadrunoModalResponse` driver (classTag **33024**, RESERVED) plus a few thin Tcl/Python command
wrappers. Concretely:

- **Reuse, don't duplicate.** `eigen` already produces Φ (nodal eigenvectors via
  `Node::getEigenvectors()`) and λ_α = ω_α² (`Domain::getEigenvalues()`,
  `Domain.h:211`). `DomainModalProperties` already computes generalized mass m_α
  (`diag(VᵀMV)`), participation factors Γ_αi, effective mass, eigenvector scale factors, and the
  rigid-body influence vector R. **All ingredients for FRF/SSD/random/RS already exist** —
  this ADR adds the *algebra layer* on top of them, not new infrastructure.
- **Extend `responseSpectrumAnalysis`, don't fork it.** The existing
  `ResponseSpectrumAnalysis` (Petracca/ASDEA) is explicitly *per-mode-only*: its source comment
  says "we just compute modal displacements without doing any (SRSS, CQC, etc.) modal combination
  … it's up to the user" (`ResponseSpectrumAnalysis.cpp:96–100`). We add an opt-in **combination
  stage** (CQC / SRSS / ABS / Ten-Percent) downstream of the per-mode loop, leaving the default
  byte-faithful to today's behavior.
- **Exact modal transient.** Each decoupled SDOF is integrated with the **exact piecewise-linear
  load recurrence** (Abaqus TG §2.5.3), not an approximate step scheme — it matches direct Newmark
  on a linear model to integration tolerance and is far cheaper.
- **Defer the direct (un-projected) route.** Direct SSD/harmonic needs the assembled *complex*
  operator `(K − Ω²M + iΩC)` and a complex linear solve — this touches the SOE/solver layer and is
  the heavy alternative. It is explicitly **out of scope here** and noted as future work that
  shares the complex-solver machinery with the complex-FEAST sibling (ADR 43).

**Recommended build order:** P1 (modal transient + CQC/SRSS) → P2 (FRF/SSD) → P3 (random/PSD→RMS).
P1 delivers the highest-demand items (fast linear seismic, native combination) and validates the
modal-recovery plumbing that P2/P3 reuse.

---

## 3. Scope-fence & classTag

**classTag: `LadrunoModalResponse` = 33024 — RESERVED, not yet built.** (33000–33018 are in use
per [[LEDGER_implementations]]; 33024 is clear of that band and of the HANDLER sub-band 33001/33002.)
A single analysis/post-processor object owns the modal-recovery, frequency-sweep, and combination
machinery; the public commands (§5) are thin wrappers around it. It is an *analysis driver*, not an
Element/Material/Integrator, so it carries a classTag only for object-broker registration symmetry
and recorder addressing — it does not enter the FE assembly.

### In scope

| Item | What | Phase |
|---|---|---|
| **Modal transient** | Exact piecewise-linear SDOF recurrence per mode; relative-to-base formulation for support motion; nodal/element response recovery `u=Σφ_α q_α` | P1 |
| **CQC / SRSS / ABS / Ten-Percent** | Modal-combination stage extending `responseSpectrumAnalysis` | P1 |
| **Modal FRF** | Complex transfer function H_α(Ω) per mode; nodal complex FRF, magnitude/phase | P2 |
| **Steady-state dynamics (SSD)** | Frequency sweep (linear/log/resonance-biased) of the harmonic response amplitude | P2 |
| **Random response** | Input PSD/CSD → modal response PSD → physical PSD → **RMS**; optional peak / zero-crossing-rate statistics | P3 |

### NOT in scope (noted as future)

- **Direct (full-system) SSD / harmonic** with a complex solver `(K − Ω²M + iΩC)u = f(Ω)` — needs
  a complex `LinearSOE` and complex sparse factorization; shares machinery with the complex-FEAST
  eigensolver. Deferred to a follow-up, links to **ADR 43** (robust+parallel/complex eigensolver).
- **Non-classical damping in the modal space** (full reduced Φᵀ C Φ with off-diagonal coupling).
  The *time-domain* coupled form is mentioned in §4.6 as a one-time-factorization extension; the
  *eigenvalue* treatment (complex modes ψ, complex ω_d, ζ) belongs to **ADR 46**. The default here
  is **classical (diagonal) modal damping**; non-classical escalation is flagged, not built.
- **Structural / complex-stiffness damping** (`F_s = i s F_int`) — only meaningful once a
  frequency-domain solver exists; admissible in SSD/random but deferred with direct SSD.
- **Parallel post-processing** — the toolkit is serial-first (matching today's
  `modalProperties`/`responseSpectrumAnalysis`, which are serial-only). Parallel reduction is an
  ADR-43 concern.

---

## 4. Formulation

Notation: `M, C, K` are the assembled system matrices; `Φ = [φ_1 … φ_p]` the real undamped mode
shapes from `eigen`; `λ_α = ω_α²` the eigenvalues; `m_α = φ_αᵀ M φ_α` the generalized (modal) mass;
`ξ_α` the modal damping ratio; `q_α(t)` the modal coordinate; `Ω` the forcing circular frequency.

### 4.1 Modal decoupling

For a linear system `M ü + C u̇ + K u = P(t)` with **classical** damping (C diagonalized by Φ),
substituting `u = Φ q` and pre-multiplying by `φ_αᵀ` decouples into `p` independent SDOF
oscillators:

$$
\ddot q_\alpha + 2\,\xi_\alpha \omega_\alpha\, \dot q_\alpha + \omega_\alpha^2\, q_\alpha
   = f_\alpha(t),
\qquad
f_\alpha(t) = \frac{\phi_\alpha^{\mathsf T} P(t)}{m_\alpha},
\qquad
\xi_\alpha = \frac{c_\alpha}{2 m_\alpha \omega_\alpha}.
$$

Physical recovery is the modal sum `u = \sum_\alpha \phi_\alpha q_\alpha`, and likewise for any
linear derived quantity (`σ = Σ σ_α q_α`, reaction, drift).

### 4.2 Exact piecewise-linear SDOF solution (modal transient)

Assume the modal load varies **linearly within each time step** `[t_n, t_{n+1}]`, `Δt = t_{n+1}-t_n`:

$$
f_\alpha(t) = f_n + \frac{f_{n+1}-f_n}{\Delta t}\,(t - t_n).
$$

Then each modal ODE has a closed-form (homogeneous + particular) solution, advanced by a constant
recurrence matrix per mode (Nigam–Jennings form). For the **under-damped** case
(`0 < ξ_α < 1`), with `ω_d = ω_α \sqrt{1-\xi_\alpha^2}`,
`e = \exp(-\xi_\alpha \omega_\alpha \Delta t)`, `s = \sin(\omega_d \Delta t)`,
`c = \cos(\omega_d \Delta t)`:

$$
\begin{bmatrix} q_{n+1} \\ \dot q_{n+1} \end{bmatrix}
=
\mathbf{A}_\alpha
\begin{bmatrix} q_{n} \\ \dot q_{n} \end{bmatrix}
+
\mathbf{B}_\alpha
\begin{bmatrix} f_{n} \\ f_{n+1} \end{bmatrix},
$$

with the propagation block

$$
\mathbf{A}_\alpha = e
\begin{bmatrix}
 c + \dfrac{\xi_\alpha\omega_\alpha}{\omega_d}\,s & \dfrac{1}{\omega_d}\,s \\[1.2em]
 -\dfrac{\omega_\alpha^2}{\omega_d}\,s & c - \dfrac{\xi_\alpha\omega_\alpha}{\omega_d}\,s
\end{bmatrix},
$$

and the load block `\mathbf{B}_\alpha` formed from the standard ramp-response integrals (closed
form in `e, s, c, ω_α, ξ_α, Δt`; computed once per mode for a fixed `Δt`). Because `A_α, B_α`
depend only on `(ω_α, ξ_α, Δt)`, the whole transient is a per-mode matrix–vector recurrence — no
iteration, no factorization, **exact** for piecewise-linear load.

**Critically-damped** (`ξ_α = 1`): the trigonometric block degenerates to its limit

$$
\mathbf{A}_\alpha = e\,
\begin{bmatrix}
 1 + \omega_\alpha \Delta t & \Delta t \\
 -\omega_\alpha^2 \Delta t & 1 - \omega_\alpha \Delta t
\end{bmatrix},
\qquad e = \exp(-\omega_\alpha \Delta t),
$$

with the corresponding ramp-load block (use this branch for `|ξ_α - 1| < \epsilon` to avoid the
`1/\sqrt{1-\xi^2}` blow-up).

**Over-damped** (`ξ_α > 1`): replace the circular functions with hyperbolic ones,
`ω_h = ω_α\sqrt{\xi_\alpha^2 - 1}`, `c→\cosh(ω_h Δt)`, `s→\sinh(ω_h Δt)`, `ω_d→ω_h`.

**Rigid-body / zero-frequency** (`ω_α = 0`, e.g. unrestrained mode): the ODE reduces to
`\ddot q_\alpha = f_\alpha`, integrated by exact double quadrature of the linear ramp:

$$
q_{n+1} = q_n + \dot q_n \Delta t + \tfrac{\Delta t^2}{6}\,(2 f_n + f_{n+1}),
\qquad
\dot q_{n+1} = \dot q_n + \tfrac{\Delta t}{2}\,(f_n + f_{n+1}).
$$

(Only **mass-proportional** damping is admissible on a rigid-body mode; stiffness-proportional
damping vanishes there — guard against `2ξ_αω_α·q̇` with `ω_α=0`.)

### 4.3 Base-motion (relative formulation)

For uniform support acceleration `ü_g(t)` along influence vector `R` (the rigid-body influence
already assembled by `modalProperties`), write `u = u_rel + R u_g`. The relative equation is

$$
M \ddot u_{\text{rel}} + C \dot u_{\text{rel}} + K u_{\text{rel}} = -\,M R\, \ddot u_g(t),
$$

so the modal load becomes

$$
f_\alpha(t) = -\,\frac{\phi_\alpha^{\mathsf T} M R}{m_\alpha}\,\ddot u_g(t)
            = -\,\Gamma_\alpha\, \ddot u_g(t),
$$

where `Γ_α = (φ_αᵀ M R)/m_α` is **exactly the participation factor `modalProperties` already
computes** (`DomainModalProperties::modalParticipationFactors()`). Total response adds the base
motion back: `u_total = u_rel + R u_g`.

### 4.4 Complex frequency-response function (FRF / SSD)

For harmonic excitation `f_α(t) = \hat f_α e^{i\Omega t}`, the steady modal amplitude is
`q_α = H_α(Ω)\,\hat f_α` with the **complex modal transfer function**

$$
\boxed{\;
H_\alpha(\Omega) = \frac{1}{\omega_\alpha^2 - \Omega^2 + 2 i\,\xi_\alpha \omega_\alpha \Omega}
\;}
\qquad\Longleftrightarrow\qquad
H_\alpha(\Omega) = \frac{1}{m_\alpha\big(\omega_\alpha^2-\Omega^2\big) + i\,\Omega\,c_\alpha}
$$

(the two forms agree under `f_α = φ_αᵀP/m_α`; `c_α = 2 m_α ξ_α ω_α`). The physical complex
response (displacement FRF) is the modal sum

$$
\hat u(\Omega) = \sum_{\alpha=1}^{p} \phi_\alpha\, H_\alpha(\Omega)\,
                 \frac{\phi_\alpha^{\mathsf T}\hat P(\Omega)}{m_\alpha},
$$

reported as magnitude `|\hat u|` and phase `\angle \hat u`. **SSD** evaluates this over a frequency
sweep `{Ω_k}` (linear, log, or resonance-biased clustering near each ω_α). Velocity/acceleration
FRFs follow by `iΩ` / `-Ω²` multipliers.

### 4.5 Random response — CSD → PSD → RMS

For a stationary, ergodic excitation defined by an input **cross-spectral-density** matrix
`S_{PP}(Ω)` (force/force, or for base input `S_{üü}(Ω)` scaled by `MR`), the output PSD propagates
through the (vector) frequency-response operator `H(Ω)` (which maps load DOFs to response DOFs via
the modal sum of §4.4):

$$
\boxed{\;
S_{uu}(\Omega) = H^{*}(\Omega)\; S_{PP}(\Omega)\; H(\Omega)^{\mathsf T}
\;}
$$

(`*` = complex conjugate). The mean of every dynamic variable is zero (response is about static
equilibrium), so variance = mean-square, and the **RMS** of any scalar response `x` is the
square-root of the integral of its auto-PSD over the analysis band:

$$
\sigma_x = \sqrt{\int_{\Omega_{\min}}^{\Omega_{\max}} S_{xx}(\Omega)\, d\Omega }
\qquad\text{(numerically: trapezoidal / Simpson over the sweep grid).}
$$

Optional second-moment statistics (Vanmarcke): with spectral moments
`λ_k = \int \Omega^k S_{xx}(\Omega)\,d\Omega`, the **mean zero-crossing rate** is
`ν_0 = \tfrac{1}{2\pi}\sqrt{\lambda_2/\lambda_0}` and the expected peak factor follows for a
specified duration — enabling a peak-response estimate `\hat x \approx p_f\,\sigma_x`.

### 4.6 CQC, SRSS, and combination rules (response spectrum)

Per-mode peak response (the quantity `responseSpectrumAnalysis` already produces, for excitation
direction `k`):

$$
R_\alpha = \big(\text{modal response of mode }\alpha\big)
         \;\propto\; \Gamma_{\alpha k}\, \frac{S_a(\omega_\alpha,\xi_\alpha)}{\omega_\alpha^2}.
$$

Because modal peaks are non-simultaneous, combine:

- **ABS** (most conservative): `R = \sum_\alpha |R_\alpha|`.
- **SRSS** (well-separated modes): `R = \sqrt{\sum_\alpha R_\alpha^2}`.
- **Ten-Percent** (Reg. Guide 1.92): SRSS plus `2\sum |R_i R_j|` over mode pairs whose periods are
  within 10%.
- **CQC** (Der Kiureghian; closely-spaced modes):

$$
R = \sqrt{\;\sum_{i}\sum_{j} \rho_{ij}\, R_i\, R_j\;},
\qquad
\rho_{ij} = \frac{8\sqrt{\xi_i \xi_j}\,(\xi_i + r\,\xi_j)\, r^{3/2}}
                 {\big(1-r^2\big)^2 + 4\,\xi_i \xi_j\, r\,(1+r^2) + 4\big(\xi_i^2+\xi_j^2\big) r^2},
\qquad r = \frac{\omega_j}{\omega_i}.
$$

(For equal damping `ξ_i=ξ_j=ξ`, `ρ_{ij} = \dfrac{8\xi^2(1+r)r^{3/2}}{(1-r^2)^2 + 4\xi^2 r (1+r)^2}`;
`ρ_{ii}=1`.) CQC requires modal damping ratios — supplied by the damping conversion in §4.7.

### 4.7 Damping conversions (how ξ_α is obtained)

The modal machinery needs `ξ_α` per mode. Supported channels (combinable additively):

- **Direct modal fractions** — user gives `ξ_α` per mode (or a single value applied to all).
- **Rayleigh** `C = a_0 M + a_1 K` → exact modal fraction (shares Φ):

$$
\boxed{\;\xi_\alpha = \frac{a_0}{2\,\omega_\alpha} + \frac{a_1\,\omega_\alpha}{2}\;}
$$

(mass term damps low modes, stiffness term damps high modes; `c_{cr}=2\sqrt{mk}=2m\omega`). This
reuses OpenSees's existing `rayleigh`/`modalDamping` inputs.

- **Composite (material/region) modal damping** — per-material critical fractions mass-weighted
  into each mode (as in `modalProperties` accounting); a future refinement, not P1-blocking.
- **Caplet — non-classical check.** If the assembled `C` is *not* diagonalized by Φ (off-diagonal
  `φ_αᵀ C φ_β ≠ 0` beyond tolerance), classical decoupling is approximate. The toolkit reports a
  **non-classicality index** and warns the user to escalate to **ADR 46** (complex modes). The
  *time-domain* coupled fallback (solve the small coupled reduced system with a one-time
  factorization, TG §2.5.3) is a documented extension but is **not** in P1–P3 scope.

---

## 5. Public API

All commands consume a prior `eigen` result; the combination/FRF/SSD/random commands additionally
consume a `modalProperties` object (the same `clientData` handoff that
`responseSpectrumAnalysis` already uses — `modal.cpp:92–93`). Both Tcl and openseespy bindings are
provided, mirroring the dual registration of the existing modal commands.

### 5.1 Modal-superposition transient — `modalResponseHistory`

```tcl
# Tcl
eigen $nModes
modalProperties -unorm
modalResponseHistory -dt $dt -nsteps $n \
    {-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 $xi2 ...} \
    {-load $patternTag | -baseAccel $tsTag -dir $dir} \
    -modes {1 2 3 ...}            ;# optional mode subset (default: all extracted)
```

```python
# openseespy
ops.modalResponseHistory('-dt', dt, '-nsteps', n,
                         '-rayleigh', a0, a1,
                         '-baseAccel', tsTag, '-dir', dir,
                         '-modes', 1, 2, 3)
```

**Behavior.** Builds per-mode `(A_α, B_α)` recurrences (§4.2), drives each modal SDOF over the
`nsteps`, recovers `u=Σφ_α q_α` and **commits one domain step per time station** (so existing
recorders capture displacement/stress/reaction histories exactly as in a direct run). Base motion
uses the relative formulation (§4.3); `-baseAccel` total output adds `R u_g` back.
**Output:** committed nodal/element response per step → ordinary recorders.

### 5.2 Frequency-response & steady-state — `frequencyResponse` / `steadyStateDynamics`

```tcl
# FRF: complex response per requested DOF over a frequency grid
frequencyResponse -freq {fmin fmax nf} {-lin | -log | -biased} \
    {-damp $xi | -rayleigh $a0 $a1 | -modalDamp ...} \
    {-load $patternTag | -baseAccel -dir $dir} \
    -node $nodeTag -dof $dof \
    -out $fileName {-magphase | -realimag}

# SSD: physical harmonic amplitude (peak) over the sweep; same parsing
steadyStateDynamics -freq {fmin fmax nf} -biased ... -out $fileName
```

**Behavior.** Evaluates `H_α(Ω_k)` (§4.4) at each sweep frequency, sums the modal contributions to
the requested response, and writes `(f, |·|, phase)` or `(f, Re, Im)` rows. `-biased` clusters
points near each ω_α to resolve resonant peaks. `steadyStateDynamics` is the amplitude-reporting
sibling (peak harmonic response vs Ω); `frequencyResponse` reports the full complex FRF.
**Output:** an ASCII/CSV table (one row per frequency) and/or a recorder-addressable complex field.

### 5.3 Random response — `randomResponse`

```tcl
randomResponse -freq {fmin fmax nf} {-lin | -log} \
    {-damp $xi | -rayleigh $a0 $a1 | -modalDamp ...} \
    -inputPSD $tsTag  {-baseAccel -dir $dir | -load $patternTag} \
    [-crossPSD $i $j $tsTagRe $tsTagIm] ... \
    -node $nodeTag -dof $dof \
    -out $fileName  [-stats]
```

**Behavior.** Forms `S_{uu}(Ω)=H^* S_{PP} H^{\mathsf T}` (§4.5) on the sweep grid, integrates the
auto-PSD to **RMS** for the requested response, and (with `-stats`) reports spectral moments,
zero-crossing rate `ν_0`, and an expected-peak estimate. Input PSD/CSD are supplied as
`timeSeries` (frequency-indexed); cross-PSD entries are optional real/imag pairs.
**Output:** RMS value(s) + optional response-PSD table + statistics block.

### 5.4 CQC/SRSS extension to `responseSpectrumAnalysis`

The existing per-mode command is preserved **verbatim** (default = today's behavior). Add an
optional combination stage:

```tcl
responseSpectrumAnalysis $tsTag $dir -combine {SRSS | CQC | ABS | TenPercent} \
    [-damp $xi | -modalDamp ...]   ;# CQC needs ξ per mode
# -> writes the COMBINED nodal response field (one committed state),
#    in addition to / instead of the per-mode states.
```

When `-combine` is absent, the loop is byte-identical to the current implementation
(`ResponseSpectrumAnalysis::analyze()`). When present, the per-mode `R_α` fields are combined by
§4.6 into a single committed result the user can record directly.

---

## 6. OpenSees integration points

Ground truth (worktree `SRC/`), to **extend** not duplicate:

- **Eigenvectors / eigenvalues (inputs).**
  - `Node::getEigenvectors()` — per-node Φ columns (read by `ResponseSpectrumAnalysis.cpp:257`).
  - `Domain::getEigenvalues()` — `Domain.h:211`; returns λ_α = ω_α² (read by
    `ResponseSpectrumAnalysis.cpp:86`, `:175`).
- **Modal properties (inputs).** `DomainModalProperties`
  (`SRC/runtime/commands/analysis/modal/DomainModalProperties.{h,cpp}`): exposes
  `eigenvalues()`, `generalizedMasses()` (`diag(VᵀMV)` → m_α), `modalParticipationFactors()`
  (Γ_αi, the §4.3 base-motion load), `eigenVectorScaleFactors()` (Vscale, applied at
  `ResponseSpectrumAnalysis.cpp:271`), `totalMass()`/`totalFreeMass()`. Header accessors at
  `DomainModalProperties.h:57–68`.
- **Per-mode recovery template (copy this).** `ResponseSpectrumAnalysis::solveMode()`
  (`ResponseSpectrumAnalysis.cpp:225–283`) is the exact pattern for "loop nodes, read scaled
  eigenvector, set trial displacement, commit step." The pressure-DOF skip
  (`:267`, `ndf==6 && node_ndf==4 && i==3`) and `min(node_ndf, ndf)` clamp (`:261`) are reused
  verbatim — these are non-obvious correctness details.
- **Step/commit cycle (copy this).** `beginMode()`/`endMode()`
  (`ResponseSpectrumAnalysis.cpp:192–222`): `analysisStep(0.0)` → `setTrialDisp` → `update()` →
  `commit()`. The modal-transient driver commits one such cycle per *time station*; FRF/SSD/random
  write tables rather than committing (no per-frequency domain state needed unless the user wants
  recorder capture of a single frequency).
- **Command registration (extend this).** `modal.cpp`: `modalProperties`
  (`modal.cpp:32`) installs `responseSpectrumAnalysis` via
  `Tcl_CreateCommand(... modal_props ...)` (`modal.cpp:92–93`). New commands register the same way,
  receiving the `DomainModalProperties*` as `clientData`. Python registration mirrors the existing
  modal Python bindings.
- **Damping inputs (reuse).** Rayleigh `a0,a1` and `modalDamping` already exist in the runtime; the
  §4.7 conversion `ξ_α = a0/2ω + a1·ω/2` is computed locally from ω_α.
- **New code.** `SRC/runtime/commands/analysis/modal/LadrunoModalResponse.{h,cpp}` (the driver +
  SDOF recurrence + FRF/PSD algebra) and small per-command wrappers alongside `modal.cpp`. The
  CQC/SRSS extension is added **inside** `ResponseSpectrumAnalysis.{h,cpp}` (an upstream Petracca
  file → triggers a `LEDGER_vanilla_files` row + `// Ladruno` markers; see §10).

No new EigenSOE, no assembly-core change, no SOE/solver change. (Direct SSD *would* need a complex
`LinearSOE` — out of scope, §3.)

---

## 7. Phased roadmap + gates

Each phase is gated; follow the fork's adversarial-gate policy ([[feedback_adversarial_gate_when]]):
P1 introduces novel SDOF math and an edit to an upstream file → **full gate**; P2/P3 are dense
post-processors mirroring P1 plumbing → lighter gate if P1 carried the recovery tests.

| Phase | Deliverable | Gate (must pass) |
|---|---|---|
| **P1** | `modalResponseHistory` (exact PWL SDOF, all damping cases incl. ω=0; base-motion relative form) **+** CQC/SRSS/ABS/TenPercent extension to `responseSpectrumAnalysis` | (a) SDOF recurrence vs analytic SDOF (under/crit/over-damped + rigid-body) to 1e-10; (b) modal transient vs **direct Newmark** on a linear frame matches to integration tol; (c) CQC reproduces a published closely-spaced-mode example; (d) `-combine` absent ⇒ byte-identical to current `responseSpectrumAnalysis` |
| **P2** | `frequencyResponse` + `steadyStateDynamics` (modal FRF, H_α(Ω), resonance-biased sweep) | (a) SDOF/2-DOF FRF magnitude+phase vs closed form; (b) FRF peak frequencies = ω_α (within sweep resolution); (c) SSD amplitude at Ω→0 equals static response |
| **P3** | `randomResponse` (PSD/CSD → response PSD → RMS, optional ν_0/peak) | (a) white-noise-through-SDOF RMS vs analytic `σ² = πS_0/(2ξω³)`; (b) RMS vs a long Monte-Carlo direct time-history (within MC sampling error); (c) `S_{uu}=H^*S_{PP}H^T` Hermitian/PSD-positive check |

P1 ships the highest-value features and is independently useful; P2/P3 are pure additions.

---

## 8. Validation / oracle battery

- **SDOF analytic (P1a).** Single mode, ramp/step/sine modal load; compare the §4.2 recurrence to
  the exact SDOF response in all four damping branches (under, critical, over, ω=0). Tolerance 1e-10
  (it is the *exact* PWL solution, so this is a coding check, not an accuracy bound).
- **Modal vs direct Newmark (P1b) — the headline gate.** A small linear 2-D/3-D frame under a
  recorded ground motion: run (i) full direct `Newmark` and (ii) `modalResponseHistory` with all
  modes. Roof-displacement and base-shear histories must match to **integration tolerance**
  (the only difference is exact-PWL vs Newmark's approximation; with enough modes and equal Δt they
  converge). This proves the participation/recovery plumbing end-to-end.
- **CQC closely-spaced (P1c).** Reproduce a textbook/Der-Kiureghian CQC example with two nearly
  equal frequencies; SRSS under-predicts, CQC matches the reference combined response.
- **FRF analytic (P2a).** SDOF and 2-DOF: `|H(Ω)|` and phase vs the closed-form transfer function,
  including the resonant peak height `1/(2ξ)` and the 90°-at-resonance phase.
- **Static limit (P2c).** SSD at Ω→0 reproduces the static displacement under the load amplitude.
- **White-noise RMS (P3a).** SDOF under flat input PSD `S_0`: numerical RMS vs the analytic
  `σ_x = \sqrt{\pi S_0 / (2 \xi \omega^3)}`.
- **Random vs Monte-Carlo (P3b).** Generate a long synthetic time history with the target input PSD,
  run direct integration, compute response RMS; compare to `randomResponse` RMS within Monte-Carlo
  sampling error.
- All oracles live in `tests/` as pytest (numpy), per the fork's Zone-A convention; the SDOF and
  FRF references are pure-numpy closed forms (self-contained, no external data).

---

## 9. Risk register

> [!question]
> **Modal truncation / missing mass.** Retaining only `p` modes truncates the high-frequency tail.
> The cumulative effective-mass ratio (already in `modalProperties`) flags shortfall; the standard
> remedy is a **static (missing-mass / residual-mode) correction** that adds the quasi-static
> response of the un-captured modes. *Sub-item:* expose the effective-mass shortfall in the command
> output and offer an optional residual-correction term (LS-DYNA `MISSING_MASS_CORRECTION`
> analogue). P1 reports the shortfall; the correction term is a P2/P3 refinement.

> [!question]
> **Accuracy vs direct integration.** Exact for *linear* systems with the retained modes; any
> material/geometric nonlinearity invalidates superposition. The commands must hard-check that the
> model is being used linearly (or at least warn) and document that nonlinearity ⇒ use direct
> integration. The P1b gate quantifies the linear-case agreement.

> [!question]
> **Non-classical damping.** Localized dashpots / dissimilar-material damping make `Φᵀ C Φ`
> non-diagonal; classical decoupling is then approximate. The toolkit computes a non-classicality
> index and **warns + points to ADR 46** (complex modes) rather than silently producing wrong
> phase. The coupled time-domain fallback (one-time factorization) is a documented, deferred
> extension.

> [!question]
> **Serial vs parallel.** Like `modalProperties`/`responseSpectrumAnalysis` today, this toolkit is
> **serial-first** and (under OpenSeesSP) sees only the master domain. Parallel reduction of the
> modal sums is an ADR-43 concern; until then the commands should refuse/warn under a partitioned
> domain rather than return partial results.

- **Frequency-sweep resolution.** Under-resolved sweeps miss sharp resonances (high-Q, low ξ) →
  underestimate FRF/SSD peaks and bias the random-response RMS integral. Mitigation: resonance-
  biased grid (`-biased`) clustering points around each ω_α; document a min points-per-peak rule.
- **PSD integration error.** RMS is an integral over the band; coarse grids or a band that omits a
  resonance corrupt `σ_x`. Mitigation: adaptive/biased grid + the spectral-moment consistency check
  in P3 gate (c).
- **Complex arithmetic / convention.** FRF sign and conjugate conventions (`e^{+iΩt}` vs `e^{-iΩt}`)
  must be fixed once and tested against the closed-form phase to avoid a silent sign flip — a
  classic frequency-domain bug; pin it in a `LEDGER_quirks` entry.

---

## 10. Ledger / header / PR obligations

- **`LEDGER_implementations.md`** — add a row: feature *Frequency-domain & modal-superposition
  response (FRF/SSD/random/modal-transient + CQC/SRSS)*, kind *analysis post-processor*, class tag
  **33024 (`LadrunoModalResponse`, RESERVED until P1 lands)**, files
  `SRC/runtime/commands/analysis/modal/LadrunoModalResponse.{h,cpp}` (+ command wrappers), status
  per-phase, PR link. Record 33024 here so it cannot collide.
- **`LEDGER_vanilla_files.md`** — the CQC/SRSS extension edits the **upstream** Petracca files
  `ResponseSpectrumAnalysis.{h,cpp}` (and possibly `modal.cpp` for the `-combine` parse). Add a row
  (file, why = "native CQC/SRSS/ABS/TenPercent combination stage", PR) and mark each edit in-source
  with a `// Ladruno …` comment so the table is grep-reconstructable.
- **`LEDGER_quirks.md`** — record any gotchas found during build (FRF conjugate/sign convention;
  pressure-DOF skip reuse; ω=0 rigid-body damping guard; sweep resolution vs resonance).
- **Banner** — add one line to `Ladruno_scripts/banner_features.txt` (e.g.
  `frequency domain: modal FRF / SSD / random / modal-transient + CQC/SRSS`), then
  `python Ladruno_scripts/patch_banner.py` and rebuild. Every `shipped` ledger row needs a matching
  banner line.
- **Header stamp** — run `Ladruno_scripts/stamp_headers.py` on the new `LadrunoModalResponse.*`
  files (add them to its GLOBS) so they carry the four-author LADRUNO header. The **edited upstream**
  files are exempt from the stamp (vanilla edits get `// Ladruno` markers only).
- **PR base** — base on `ladruno` (default branch), not `main`; one logical PR per phase; verify the
  prior PR is merged before stacking (see [[feedback_stranded_commits_after_automerge]]).

---

## 11. Cross-references

- [[modal_gap_study/00_SYNTHESIS]] **§4 (frequency domain)** and **ADR-D** in §5 — the parent
  recommendation that this ADR realizes; the "cheap post-processors on eigen" framing.
- [[modal_gap_study/02_abaqus_theory]] **GAP 4** — TG §2.5.3 (exact PWL modal transient), §2.5.6
  (response spectrum + CQC/SRSS/Ten-Percent), §2.5.7 (SSD + H_α(Ω)), §2.5.8 (random `S_out=H*S_in Hᵀ`),
  §2.5.4 (damping conversions).
- [[modal_gap_study/04_lsdyna_theory]] **§3** — `*FREQUENCY_DOMAIN_{FRF,SSD,RANDOM_VIBRATION,
  RESPONSE_SPECTRUM}` parity targets; modal-vs-direct split; `MISSING_MASS_CORRECTION`.
- [[modal_gap_study/01_opensees_current_state]] **§A/§B** — the `eigen` pipeline and the
  `modalProperties` / `responseSpectrumAnalysis` infrastructure this extends.
- **ADR 46** (`46_ladruno_complex_modal_adr`) — non-classical damping → complex state-space modes
  (escalation target when `Φᵀ C Φ` is non-diagonal).
- **ADR 42** (`42_ladruno_buckling_adr`) — prestressed modal + linear buckling (sibling; supplies
  preloaded Φ when a frequency-domain run rides a prestressed state).
- **ADR 43** (`43_ladruno_feast_eigensolver_adr`) — robust + parallel eigensolver (FEAST/complex-FEAST);
  owns the **direct-SSD complex linear solver** this ADR defers, and the parallel post-processing
  reduction.
