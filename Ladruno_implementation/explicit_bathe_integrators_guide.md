---
title: The Explicit Bathe (Noh–Bathe) Integrators — Theory, Implementation & Use
project: Ladruno
status: reference
priority: high
tags:
  - integrator
  - explicit
  - bathe
  - noh-bathe
  - time-integration
  - dynamics
  - theory
aliases:
  - ExplicitBathe
  - ExplicitBatheLNVD
  - Noh-Bathe explicit
---

# The Explicit Bathe Integrators

> [!abstract] What this document is
> A thorough, self-contained reference for the two Ladruno explicit time
> integrators built on the **Noh–Bathe** two-sub-step composite scheme:
> **`ExplicitBathe`** (classTag **33000**) and **`ExplicitBatheLNVD`** (classTag
> **33002**, the same scheme + FLAC local non-viscous damping for dynamic
> relaxation). Both live in the Ladruno private classTag band (≥33000). It covers the continuum/ODE **theory**, the **derivation** of the
> sub-step update, the **stability and dissipation** analysis that motivates the
> design, the actual **OpenSees implementation** (the `TransientIntegrator`
> contract, the `newStep → update → commit` life-cycle, the mass-only solve, the
> critical-time-step machinery), and the **intended use cases**.
>
> Companion document: [[central_difference_ladruno_guide]] (the non-dissipative
> explicit leap-frog sibling). Source-of-record ADRs:
> [[04_explicit_dynamics_and_energy_balance]] and [[Ladruno_explicit_roadmap]].

---

## 1. The problem these integrators solve

### 1.1 The semi-discrete equation of motion

After finite-element spatial discretization, the equilibrium of a structure in
motion is a system of coupled second-order ODEs in time. In the three notation
forms an analyst keeps in mind:

**Direct (matrix/operator) form**

$$
\mathbf{M}\,\ddot{\mathbf{u}}(t) \;+\; \mathbf{C}\,\dot{\mathbf{u}}(t) \;+\; \mathbf{F}^{\text{int}}(\mathbf{u}(t)) \;=\; \mathbf{P}(t)
$$

**Index (component) form** — with $a,b$ ranging over the $n_{\text{dof}}$ global degrees of freedom:

$$
M_{ab}\,\ddot{u}_b \;+\; C_{ab}\,\dot{u}_b \;+\; F^{\text{int}}_a(\mathbf u) \;=\; P_a(t), \qquad a = 1,\dots,n_{\text{dof}}
$$

**Linear (small-displacement) specialization** — when $\mathbf{F}^{\text{int}} = \mathbf{K}\mathbf{u}$:

$$
\mathbf{M}\ddot{\mathbf u} + \mathbf{C}\dot{\mathbf u} + \mathbf{K}\mathbf u = \mathbf P(t)
$$

Here $\mathbf{M}$ is the mass matrix, $\mathbf{C}$ the (viscous) damping matrix,
$\mathbf{F}^{\text{int}}$ the internal restoring force (= $\mathbf K\mathbf u$ in
the linear case), and $\mathbf{P}$ the applied external load. The physical
content is Newton's second law assembled element-by-element: *inertia +
dissipation + elastic restoring = applied force*, holding at every instant.

A **time integrator** turns this continuous-in-time system into a recursion that
marches the state $(\mathbf u_n, \dot{\mathbf u}_n, \ddot{\mathbf u}_n)$ from step
$n$ to $n+1$ over a step $\Delta t$.

### 1.2 Explicit vs. implicit — the fork in the road

> [!info] The defining question of a time integrator
> *Does advancing the state require solving a system involving the stiffness
> (tangent) matrix?*

- **Implicit** (Newmark-β with $\beta>0$, HHT-α, generalized-α, trapezoidal):
  the unknown acceleration/displacement at $t_{n+1}$ appears inside
  $\mathbf F^{\text{int}}$, so each step solves
  $\mathbf K_{\text{eff}}\,\Delta\mathbf u = \mathbf R$ — a global, often
  Newton-iterated, factorization of an effective stiffness. Unconditionally
  stable (for the right parameters), but every step is expensive and, in the
  nonlinear regime, *can fail to converge* near snap-through, contact, softening,
  or wave fronts.

- **Explicit** (central difference, Noh–Bathe, the schemes in this document): the
  state at $t_{n+1}$ is computed from **already-known** quantities at $t_n$ (and
  intermediate sub-steps). The only matrix that appears on the left-hand side is
  the **mass** $\mathbf M$. If $\mathbf M$ is **lumped (diagonal)**, the "solve"
  is a trivial element-wise division $\ddot{\mathbf u} = \mathbf M^{-1}\mathbf r$
  — no factorization, no iteration, no convergence failure. The price is
  **conditional stability**: $\Delta t$ must stay below a critical value
  $\Delta t_{\text{cr}}$ set by the highest natural frequency of the mesh.

Explicit integration is the roadmap's chosen path for the regimes where implicit
Newton is fragile: **collapse, impact/contact, wave propagation, soil–structure
interaction (SSI)**. The trade — many cheap, guaranteed-to-complete steps instead
of few expensive, possibly-non-converging ones — is exactly the right trade when
the physics is highly nonlinear or wave-dominated.

### 1.3 Why a *new* explicit scheme — the dissipation problem

Plain central difference (see [[central_difference_ladruno_guide]]) is **second
order accurate with zero algorithmic damping**. That is wonderful for energy
conservation but dangerous in practice: the **highest-frequency modes of an FE
mesh are spurious** — they are artifacts of the discretization, carry no physical
meaning, and (being the least accurately integrated) are exactly the modes most
likely to ring, pollute the solution, or seed instability. A non-dissipative
scheme keeps them forever.

What you want for robust explicit dynamics is a scheme that:

1. is still **explicit** (mass-only LHS, no tangent solve),
2. is still **second-order accurate** for the resolved low modes, **and**
3. provides **controllable, targeted high-frequency numerical dissipation** that
   quietly kills the spurious high modes while barely touching the physical low
   ones.

The **Noh–Bathe explicit scheme** (Gunwoo Noh & Klaus-Jürgen Bathe, *"An explicit
time integration scheme for the analysis of wave propagations,"* Computers &
Structures **129** (2013) 178–193, [doi:10.1016/j.compstruc.2013.06.007](https://doi.org/10.1016/j.compstruc.2013.06.007))
is precisely that scheme. It is the theoretical foundation of both
`ExplicitBathe` and `ExplicitBatheLNVD`.

---

## 2. Theory of the Noh–Bathe explicit scheme

### 2.1 The composite two-sub-step idea

The Noh–Bathe method splits each step $[t,\, t+\Delta t]$ into **two sub-steps**
separated by a parameter $p \in (0,1)$:

$$
t \;\xrightarrow{\;\text{sub-step 1}\;}\; t + p\,\Delta t \;\xrightarrow{\;\text{sub-step 2}\;}\; t + \Delta t
$$

- **Sub-step 1** ($t \to t+p\Delta t$) is an explicit **central-difference-like**
  advance: a Taylor predictor for displacement using the known acceleration
  $\mathbf a_t$, then a single mass-only solve for $\mathbf a_{t+p\Delta t}$, then
  a trapezoidal velocity correction.
- **Sub-step 2** ($t+p\Delta t \to t+\Delta t$) is a second explicit advance, but
  its **velocity update blends three accelerations** ($\mathbf a_t$,
  $\mathbf a_{t+p\Delta t}$, $\mathbf a_{t+\Delta t}$) with weights $q_0,q_1,q_2$
  chosen so the *composite* one-step map has the desired stability and dissipation.

The magic is in how the weights are chosen. By tuning $p$ (and the dependent
$q_i$) the composite map gains **frequency-selective damping**: at the parameter
$p\ne\tfrac12$ the scheme damps high-frequency content while remaining
second-order accurate and *more* stable (larger critical step) than plain central
difference.

### 2.2 The update equations (as implemented)

Define $\tau_1 = p\,\Delta t$ and $\tau_2 = (1-p)\,\Delta t$.

**Sub-step 1 — predictor + mass-only solve + corrector** ($t \to t+p\Delta t$):

$$
\begin{aligned}
\mathbf u_{t+p\Delta t} &= \mathbf u_t + \tau_1\,\dot{\mathbf u}_t + \tfrac{1}{2}\tau_1^{2}\,\ddot{\mathbf u}_t && \text{(displacement Taylor predictor)}\\[2pt]
\dot{\mathbf u}^{\,*}_{t+p\Delta t} &= \dot{\mathbf u}_t + \tau_1\,\ddot{\mathbf u}_t && \text{(trial velocity for the residual)}\\[2pt]
\mathbf M\,\ddot{\mathbf u}_{t+p\Delta t} &= \mathbf P_{t+p\Delta t} - \mathbf C\,\dot{\mathbf u}^{\,*}_{t+p\Delta t} - \mathbf F^{\text{int}}(\mathbf u_{t+p\Delta t}) && \text{(solve for acceleration)}\\[2pt]
\dot{\mathbf u}_{t+p\Delta t} &= \dot{\mathbf u}_t + \tfrac{1}{2}\tau_1\,\big(\ddot{\mathbf u}_t + \ddot{\mathbf u}_{t+p\Delta t}\big) && \text{(trapezoidal velocity correction)}
\end{aligned}
$$

**Sub-step 2 — predictor + mass-only solve + 3-point velocity update** ($t+p\Delta t \to t+\Delta t$):

$$
\begin{aligned}
\mathbf u_{t+\Delta t} &= \mathbf u_{t+p\Delta t} + \tau_2\,\dot{\mathbf u}_{t+p\Delta t} + \tfrac{1}{2}\tau_2^{2}\,\ddot{\mathbf u}_{t+p\Delta t}\\[2pt]
\dot{\mathbf u}^{\,*}_{t+\Delta t} &= \dot{\mathbf u}_{t+p\Delta t} + \tau_2\,\ddot{\mathbf u}_{t+p\Delta t}\\[2pt]
\mathbf M\,\ddot{\mathbf u}_{t+\Delta t} &= \mathbf P_{t+\Delta t} - \mathbf C\,\dot{\mathbf u}^{\,*}_{t+\Delta t} - \mathbf F^{\text{int}}(\mathbf u_{t+\Delta t})\\[2pt]
\dot{\mathbf u}_{t+\Delta t} &= \dot{\mathbf u}_{t+p\Delta t} + \tau_2\Big[\,q_0\,\ddot{\mathbf u}_t + \big(\tfrac12 + q_1\big)\ddot{\mathbf u}_{t+p\Delta t} + q_2\,\ddot{\mathbf u}_{t+\Delta t}\Big]
\end{aligned}
$$

The two left-hand-side solves are **identical in structure**: $\mathbf M$ times
the unknown acceleration equals a residual built from *known* displacements and a
*known* trial velocity. With lumped $\mathbf M$ these are diagonal divisions.

### 2.3 The coefficient family

The weights are fixed by $p$ through:

$$
q_1 = \frac{1 - 2p}{2\,p\,(1-p)}, \qquad
q_2 = \frac{1}{2} - p\,q_1, \qquad
q_0 = \frac{1}{2} - q_1 - q_2 .
$$

These three relations are exactly the lines in the constructor
(`ExplicitBathe.cpp:217–219`):

```cpp
q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
q2 = 0.5 - p * q1;
q0 = -q1 - q2 + 0.5;
```

> [!note] The $p=\tfrac12$ special case (no dissipation)
> At $p = \tfrac12$: $q_1 = 0$, $q_2 = \tfrac12$, $q_0 = 0$. The 3-point velocity
> update collapses to a trapezoidal rule and the whole method degenerates into
> **two equal central-difference half-steps** — i.e. zero algorithmic damping,
> identical in spirit to [[central_difference_ladruno_guide|CentralDifferenceLadruno]].
> Moving $p$ *away* from $\tfrac12$ is what switches the high-frequency
> dissipation on. This is why $p$ is the single "damping knob" of the scheme.

### 2.4 Why $p = 0.54$ is the default

`ExplicitBathe` and `ExplicitBatheLNVD` both default to **$p = 0.54$**. This is
the value Noh & Bathe identify as a near-optimal compromise:

- It preserves **second-order accuracy** in the low-frequency (resolved) range.
- It gives **strong, monotone high-frequency dissipation**: the spectral radius
  $\rho \to \rho_\infty < 1$ as $\omega\Delta t \to$ the stability limit, so the
  highest modes are annihilated within a few steps.
- It keeps the **stability limit large** — about twice that of central difference
  (see §2.5).

Values closer to $0.5$ damp less (and the scheme limits to non-dissipative);
values toward $1$ damp more aggressively but degrade low-frequency accuracy. The
admissible range is the open interval $p\in(0,1)$, enforced at parse time
(`p must be in (0,1)`); typical practice is $0.50\!-\!0.95$ with $0.54$
recommended.

### 2.5 Stability — the $\sim 2\times$ advantage

For the undamped linear problem, decompose into modal coordinates; each mode
$\omega_i$ behaves like a scalar SDOF oscillator. The composite one-step map of
the scheme has an **amplification matrix** $\mathbf A(\Omega)$ with
$\Omega \equiv \omega\,\Delta t$; the scheme is stable when its **spectral
radius** satisfies $\rho(\mathbf A) \le 1$ for every mode.

The conservative, always-safe reference is the **central-difference limit**

$$
\Delta t_{\text{cd}} \;=\; \frac{2}{\omega_{\max}}
$$

where $\omega_{\max}$ is the highest *element* natural frequency (see §4 — the
per-element eigenproblem). The key Noh–Bathe result, encoded as
`EB_NB_STABILITY_FACTOR = 2.0` in `ExplicitBathe.h`, is that **at the optimal $p$
the composite scheme is stable up to roughly twice that step**:

$$
\Delta t_{\text{cr}}^{\text{NB}} \;\approx\; 2\,\Delta t_{\text{cd}} \;=\; \frac{4}{\omega_{\max}} .
$$

Because each Noh–Bathe step costs *two* mass-only solves, a $2\times$ larger step
means **about the same number of solves per unit physical time as central
difference — but with built-in dissipation for free**. That is the headline
efficiency argument for the scheme.

> [!warning] The reported `dt_cr` is the conservative limit, not the NB limit
> The integrator's `criticalTimeStep()` / `getCriticalTimeStep()` always returns
> the **conservative** $\Delta t_{\text{cd}} = 2/\omega_{\max}$ (a guaranteed-safe
> lower bound), and the per-step report additionally prints the Noh–Bathe limit
> `EB_NB_STABILITY_FACTOR * dt_cd`. Sizing your step against the conservative
> value is always safe; the NB factor is the bonus headroom.

### 2.6 Damping and the critical step

With viscous damping active, the stable step shrinks. For modal damping ratio
$\xi$ the central-difference family obeys

$$
\Delta t_{\text{cr}} \;=\; \frac{2}{\omega_{\max}}\Big(\sqrt{1+\xi^{2}} - \xi\Big),
$$

which is the formula the shared eigensolver uses (`CriticalTimeStep.cpp:255`).
The damped value is reported as `damped_dt`, the undamped as `undamped_dt`. The
$\beta\mathbf K$ (stiffness-proportional) Rayleigh trap that collapses the
explicit step — $\xi = \beta\omega/2$ grows with frequency, so
$\Delta t_{\text{cr}} \to 2/(\beta\,\omega_{\max}^{2})$ — is discussed at length
for the sibling scheme in [[central_difference_ladruno_guide#The βK trap]]; the
same physics applies here, and mass-proportional $\alpha\mathbf M$ is the safe
choice.

---

## 3. What drove the implementation

The mathematics above is standard. The *engineering* decisions — why the OpenSees
classes look the way they do — are recorded in the ADR
[[04_explicit_dynamics_and_energy_balance]]. The load-bearing ones:

| Driver | Decision | Why |
|---|---|---|
| **Mass-only LHS** | `formEleTangent` adds **only** $\mathbf M$ to the tangent (`zeroTangent(); addMtoTang();`) | The whole point of explicit is a trivial $\mathbf M^{-1}$. Damping is moved to the **residual** at a *known* trial velocity, never onto the LHS (LHS-damping would make it implicit/coupled). |
| **Diagonal mass + `system Diagonal`** | Documented hard recipe | A consistent (non-diagonal) mass with a sparse solver defeats the purpose; the scheme is built for lumped mass where $\mathbf M^{-1}$ is element-wise. |
| **Exactly 2 solves/step ⇒ `algorithm Linear`** | `update()` guarded at `updateCount > 2` | The scheme is a fixed 2-sub-step recursion, *not* a Newton loop. A nonlinear algorithm would call `update()` repeatedly and corrupt the sub-step bookkeeping. The guard catches the misconfiguration with a clear error. |
| **Silent-wrongness of explicit** | Built-in **critical-time-step** estimate + **`EnergyBalanceRecorder`** | Explicit has no Newton residual to warn you when $\Delta t$ is too big — it just blows up (or worse, drifts). So the deliverable is *the scheme plus the diagnostics that prove a run stayed physical*: a sound $\Delta t_{\text{cr}}$ and an energy-closure monitor. |
| **Cold-start from rest** | One-time note if $\mathbf a_0 = \mathbf v_0 = \mathbf 0$ with load applied | An explicit scheme started with zero initial acceleration produces no motion on step 1 if a load is already present; the note tells the user to run an equilibrium step or supply a consistent $\mathbf a_0$. |
| **Circuit breakers** | NaN/Inf acceleration abort; optional kinetic-energy runaway abort (`-divergence f`) | Catch a blown-up step *immediately* instead of letting `NaN` silently propagate through the whole time history. |
| **Shared eigensolver** | `CriticalTimeStep.{h,cpp}` hoisted out and shared with the sibling | Was a hand-copied `extern` that could silently drift; now one definition (with the `DSYGVX`/`DGGEV` fallback, relative-β threshold, `-lump` option, MPI reduction). |
| **`ExplicitBatheLNVD` exists at all** | Wire up FLAC local damping (was dead commented-out code) | The LNVD variant's entire reason to exist — dynamic relaxation / quasi-static solving — was disabled in the inherited code; activating it is the differentiator (§6). |

---

## 4. The critical-time-step machinery (`CriticalTimeStep`)

Both Bathe integrators (and the sibling CD class) share one estimator:
`computeCriticalTimeStep()` in
[`SRC/analysis/integrator/CriticalTimeStep.cpp`](../SRC/analysis/integrator/CriticalTimeStep.cpp).

### 4.1 The per-element generalized eigenproblem

Rather than a global eigensolve (expensive), it exploits the theorem that the
**global maximum frequency is bounded by the maximum over elements**:
$\omega_{\max}^{\text{global}} \le \max_e \omega_{\max}^{(e)}$. So it scans every
element and solves the small **per-element generalized eigenproblem**

$$
\mathbf K^{(e)}\,\boldsymbol\phi \;=\; \lambda\,\mathbf M^{(e)}\,\boldsymbol\phi,
\qquad \omega_{\max}^{(e)} = \sqrt{\lambda_{\max}},
\qquad \Delta t_{\text{cd}}^{(e)} = \frac{2}{\omega_{\max}^{(e)}},
$$

and reports the **governing minimum** $\min_e \Delta t_{\text{cd}}^{(e)}$ together
with the tag of the binding element.

### 4.2 Numerical robustness

- **`DSYGVX` first, `DGGEV` fallback.** The lumped $\mathbf M$ is diagonal and
  (usually) positive-definite, so the symmetric-definite expert driver `DSYGVX`
  is tried first (requesting only the single largest eigenvalue via
  `RANGE='I', IL=IU=n`). If $\mathbf M$ is not positive-definite (zero-mass
  DOFs ⇒ indefinite pencil), it falls back to the general `DGGEV` with a
  **relative-β threshold** ($|\beta| > 10^{-12}\,\beta_{\max}$) to discard the
  spurious near-infinite eigenvalues of massless modes. *(This relative-β fix
  corrected a real silent bug where an exact `beta != 0` test let a near-massless
  DOF drive $\Delta t_{\text{cr}} \to 0$ and abort an otherwise-stable run.)*
- **Why `DSYGVX` not `DSYGV`** — a Ladruno-specific portability fix. The bundled
  `OTHER/LAPACK` reference library that the Linux/Ubuntu Zone-A CI build links
  statically ships `dsygvx.f` but **not** `dsygv.f`; on Windows MKL resolves
  either, so the difference only surfaces on Linux as an undefined-reference link
  error. `DSYGVX` is already proven to link (it is used by
  `SymmGeneralizedEigenSolver.cpp`).
- **Mass-before-stiffness aliasing (decision D8).** Many elements return a
  reference to the **same shared static matrix** from `getMass()` and
  `getInitialStiff()`/`getTangentStiff()`. Holding both references at once would
  alias — the second call clobbers the first. The estimator therefore **copies
  the mass diagonal first**, then fetches the stiffness. *(This caught a real bug
  where every distributed-mass model reported $\Delta t_{\text{cr}}=\infty$.)*
- **Damping reduction.** For each element it reads the Rayleigh factors, forms
  $\xi = \tfrac12(\alpha_M/\omega_{\max} + \beta_K\,\omega_{\max})$, and computes
  the damped step $\frac{2}{\omega_{\max}}(\sqrt{1+\xi^2}-\xi)$.
- **Parallel reduction.** Under a parallel build the scalar limits are
  `MPI_Allreduce(MPI_MIN)`-reduced across ranks so every rank sees the global
  governing step.

### 4.3 Lumping choice (`-lump`)

| Mode | Definition | Good for | Caveat |
|---|---|---|---|
| `rowsum` (Bathe default) | $M_{ii} \leftarrow \sum_j M_{ij}$ | Translational DOFs (bars, solids) — exact/conservative | Row sums of rotational DOFs (beams/shells) can be $\approx 0$ ⇒ non-conservative |
| `diagonal` | $M_{ii}$ taken directly (diagonal-of-consistent) | Rotational DOFs — strictly positive | **Not** mass-conserving ($\sum M_{ii} \ne$ element mass); can be non-conservative for translation. A true mass-conserving HRZ lump is a documented follow-up. |

> [!caution] Inherited `dt_cr` caveats
> The estimate **ignores constraints** (`equalDOF`, rigid diaphragms, MP
> constraints) and **pure nodal mass**, so a constrained or nodal-mass-only model
> may report a non-binding or `≤0` value. For rotational DOFs prefer
> `-lump diagonal`. Treat `dt_cr` as a strong guide, not an absolute guarantee,
> on beam/shell models.

---

## 5. OpenSees implementation

Both classes derive from `TransientIntegrator` and live in
[`SRC/analysis/integrator/`](../SRC/analysis/integrator/):

- `ExplicitBathe.{h,cpp}` — classTag `INTEGRATOR_TAGS_ExplicitBathe` (**33000**)
- `ExplicitBatheLNVD.{h,cpp}` — classTag `INTEGRATOR_TAGS_ExplicitBatheLNVD` (**33002**)
- `CriticalTimeStep.{h,cpp}` — the shared eigensolver (§4)

### 5.1 The `TransientIntegrator` life-cycle

OpenSees drives a transient analysis through a fixed sequence of virtual calls.
For an explicit scheme the mapping is:

```
analyze(nSteps, dt)
└── for each step:
      newStep(dt)        ── predictor for sub-step 1; sets trial response;
      │                     advances domain time by p·dt
      ├─ [algorithm Linear] formTangent → M only (formEleTangent)
      ├─ formUnbalance   ── residual R = P − C·v* − F_int(u)
      ├─ solve           ── M·a = R  (diagonal ⇒ a = M⁻¹R)
      ├─ update(a)  [1st] ── store a_{t+pΔt}; build sub-step-2 predictor;
      │                     advance time by (1−p)·dt; SECOND solve for a_{t+Δt};
      │                     final velocity update; push response to nodes
      └─ commit()        ── roll t+Δt state into t; commitDomain()
```

The unusual feature is that **`ExplicitBathe::update()` performs the *second*
solve internally** (`formUnbalance(); theLinSOE->solve();`), so the two sub-step
solves of one Noh–Bathe step map onto the framework's single
`formUnbalance/solve/update` cadence. This is why exactly **two solves per step**
are expected and `algorithm Linear` is mandatory.

### 5.2 Key methods, annotated

**`formEleTangent` / `formNodTangent` — the mass-only LHS:**

```cpp
int ExplicitBathe::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();   // ONLY the mass — no C, no K on the LHS
    return 0;
}
```

This is the line that makes the scheme explicit. Damping appears only inside the
residual, evaluated at the trial velocity `V_fake` set on the nodes before each
solve.

**`newStep(dt)` — sub-step 1 predictor.** Computes the per-step coefficients
`a0…a7` from $p$ and $\Delta t$, builds the displacement predictor
`U_tpdt = U_t + a0·V_t + a1·A_t`, sets the trial velocity
`V_fake = V_t + a0·A_t`, pushes them via `setResponse`, advances the domain clock
by $p\,\Delta t$. It also (optionally) refreshes $\Delta t_{\text{cr}}$, runs the
cold-start check, and applies the `-cflAbort` guard.

**`update(U)` — corrector + sub-step 2 + second solve.** `U` is the solved
$\ddot{\mathbf u}_{t+p\Delta t}$. It applies the trapezoidal velocity correction,
builds the sub-step-2 displacement predictor `U_tdt`, advances the clock by
$(1-p)\Delta t$, then **does the second mass-only solve itself**
(`formUnbalance(); solve();`) to get $\ddot{\mathbf u}_{t+\Delta t}$, runs the
NaN/Inf and kinetic-energy circuit breakers, applies the 3-point velocity update
with weights $q_0,(\tfrac12+q_1),q_2$, and pushes the final
$(\mathbf u,\dot{\mathbf u},\ddot{\mathbf u})_{t+\Delta t}$ to the nodes.

**`commit()` — roll forward.** $\mathbf u_t \leftarrow \mathbf u_{t+\Delta t}$,
likewise velocity and acceleration; resets `updateCount`; calls `commitDomain()`
(which commits path-dependent material state and fires recorders).

**`getVel()` — the modal-damping hook.** Returns the *committed* velocity
$\dot{\mathbf u}_t$. `IncrementalIntegrator::addModalDampingForce` assembles modal
damping into the residual via `setB`, so this must return a velocity known at
solve time.

**`getCriticalTimeStep()` — the public query.** Overrides the base (which returns
$-1$) to return the cached conservative `damped_minimum_critical_timestep`; the
`criticalTimeStep()` Tcl/Python command dispatches here.

### 5.3 State vectors

`ExplicitBathe` carries three time levels: committed $t$
(`U_t, V_t, A_t`), the intermediate $t+p\Delta t$
(`U_tpdt, V_tpdt, A_tpdt`), and the new $t+\Delta t$ (`U_tdt, V_tdt, A_tdt`), plus
`V_fake` (the trial velocity for `setResponse`). They are allocated/sized in
`domainChanged()` and seeded from the committed DOF disp/vel/accel.

### 5.4 Parsing, serialization, registration

- **Parsing** (`OPS_ExplicitBathe`): reads the leading numeric `p` with the typed
  getter `OPS_GetDoubleInput` (which works under both Tcl and OpenSeesPy, unlike
  `OPS_GetString` on a numeric Python arg), then the flags
  `-verbose -cfl -cflAbort -tangent -recompute N -lump rowsum|diagonal
  -divergence f`. `ExplicitBatheLNVD` additionally reads a second positional
  `alpha` (FLAC coefficient, default 0.8).
- **`sendSelf`/`recvSelf`**: `ExplicitBathe` ships just `p`; `ExplicitBatheLNVD`
  ships `(p, alpha)`. The dependent $q_i$ are recomputed on receipt.
- **Registration** (per the ADR): `classTags.h` (ExplicitBathe **33000**,
  ExplicitBatheLNVD **33002**, EnergyBalanceRecorder 26),
  `FEM_ObjectBrokerAllClasses.cpp` + `TclPackageClassBroker.cpp` (broker `case`s
  for `openseessp`/database restart), plus the `criticalTimeStep()` command
  wired into `OpenSeesCommands`, `PythonWrapper`, `TclWrapper`.

---

## 6. `ExplicitBatheLNVD` — FLAC local non-viscous damping

`ExplicitBatheLNVD` is the **same Noh–Bathe scheme** plus a **FLAC-style local
non-viscous damping** intended for **dynamic relaxation / quasi-static solving**.

### 6.1 The local damping law

At every solve the assembled residual $\mathbf r = \mathbf f^{\text{ext}} -
\mathbf f^{\text{int}}$ is modified component-wise:

$$
r_i \;\longleftarrow\; r_i \;-\; \alpha\,\lvert r_i\rvert\,\operatorname{sign}(v_i)
$$

where $v_i$ is the local trial velocity and $\alpha$ the FLAC local-damping
coefficient ($0 \le \alpha < 1$, classic default **0.8**). Physically this adds a
force **proportional to the magnitude of the local unbalanced force** and
**opposing the local velocity**. Its defining properties:

- It **removes kinetic energy** — it always opposes motion, damping oscillation
  toward equilibrium.
- It **vanishes at equilibrium**: as $\mathbf v \to \mathbf 0$ *and*
  $\mathbf r \to \mathbf 0$, the damping force $\to 0$. So it **does not bias the
  converged static solution** — unlike viscous damping, which leaves a residual
  $\mathbf C\dot{\mathbf u}$ term.
- It is **frequency-independent and mass/stiffness-independent** — the damping on
  each DOF self-scales to that DOF's own unbalanced force, which is what makes
  FLAC local damping so robust for multi-scale quasi-static problems.

The FLAC coefficient relates to an approximate fraction of critical damping by
$\alpha \approx \pi D$, so the default $\alpha = 0.8$ is roughly $D \approx 0.25$
of critical at the dominant response frequency.

### 6.2 Implementation — symmetric application to both sub-steps

The damping is injected through a **`formUnbalance()` override**:

```cpp
int ExplicitBatheLNVD::formUnbalance(void) {
    int res = this->IncrementalIntegrator::formUnbalance();   // standard residual
    if (res < 0) return res;
    lastUnbalanceNorm = theSOE->getB().pNorm(0);              // convergence monitor
    if (alpha_flac > 0.0) this->addLocalDamping();            // r_i -= α|r_i|sign(v_i)
    return res;
}
```

Because the solution algorithm calls `formUnbalance()` for **sub-step 1** and the
integrator's own `update()` calls it for **sub-step 2**, the damping is applied
**symmetrically to both Noh–Bathe sub-steps** — a subtlety the ADR (decision D3)
flags explicitly. `addLocalDamping()` gathers the global trial velocity from the
DOF groups, then forms $r_i - \alpha|r_i|\operatorname{sign}(v_i)$ and writes it
back with `setB`.

### 6.3 The convergence indicator

`getUnbalanceNorm()` exposes $\lVert\mathbf r\rVert_\infty$ at the most recent
solve. For dynamic relaxation you **watch this norm drive toward zero** to judge
that the model has reached static equilibrium — it is the LNVD analogue of a
Newton convergence test.

---

## 7. Intended use cases

### 7.1 `ExplicitBathe` — true explicit dynamics with HF dissipation

> [!success] Reach for `ExplicitBathe` when…
> - You are doing **wave propagation**, **impact/blast**, **collapse**, or **SSI**
>   where implicit Newton is fragile or the step count is enormous.
> - You want **controllable high-frequency numerical damping** to quietly kill
>   spurious mesh modes while keeping the physical response second-order accurate.
> - Your mass is (or can be) **lumped/diagonal** so $\mathbf M^{-1}$ is trivial.
> - You value the **$\sim2\times$ larger stable step** over plain central
>   difference at the cost of a second mass-only solve per step.

**Required recipe:**

```python
ops.mass(...)                                  # lumped/diagonal element or nodal mass
ops.system('Diagonal')                         # trivial M⁻¹ — not a sparse solver
ops.algorithm('Linear')                        # exactly 2 solves/step
ops.integrator('ExplicitBathe', 0.54, '-cfl')  # p=0.54; -cfl computes dt_cr
ops.analysis('Transient')

# driver-level adaptive sub-stepping against the queried dt_cr (ADR decision D5)
import math
dt_cr = ops.criticalTimeStep()
n = max(1, math.ceil(dt / (0.9 * dt_cr)))      # 0.9 = LS-DYNA-style TSSFAC
ops.analyze(n, dt/n)
```

Pair with `recorder EnergyBalance -file energy.txt -time` to watch the closure
residual `RES = ULW − (KE+IE+DW)` — spurious energy growth is the unmistakable
signature of an over-large $\Delta t$.

### 7.2 `ExplicitBatheLNVD` — quasi-static / dynamic relaxation

> [!success] Reach for `ExplicitBatheLNVD` when…
> - You need a **static equilibrium** but implicit static solvers struggle:
>   **snap-through**, **near-singular tangents**, **post-buckling**, strongly
>   **softening** materials, or initial **self-weight settling** of a large model.
> - You want to **drive the model to rest** using the explicit machinery, with the
>   FLAC local damping bleeding off kinetic energy without polluting the final
>   static field.

```python
ops.integrator('ExplicitBatheLNVD', 0.54, 0.8)   # p=0.54, FLAC alpha=0.8
# ... advance many steps; monitor convergence:
# the unbalanced-force norm → 0 signals static equilibrium reached
```

> [!tip] LNVD vs. true dynamics
> Use LNVD **only** to reach a static state. For genuine transient dynamics use
> the undamped `ExplicitBathe` (local damping would corrupt the real dynamic
> response).

### 7.3 When *not* to use these

- **Stiff/quasi-static problems where implicit is fine** — an unconditionally
  stable implicit scheme will take far larger steps.
- **Consistent-mass models with a sparse solver** — defeats the explicit premise;
  use lumped mass and `system Diagonal`.
- **Stiffness-proportional ($\beta\mathbf K$) Rayleigh damping** — it collapses
  the explicit critical step quadratically; prefer mass-proportional
  $\alpha\mathbf M$ (see [[central_difference_ladruno_guide#The βK trap]]).

---

## 8. Validation

The numerical battery (`Ladruno_scripts/_verify_explicit.py`, **9/9** at the time
of merge, ADR §How) checks:

| # | Check | Result |
|---|---|---|
| 1 | Order of accuracy (SDOF free-vibration log–log slope) | **1.99** ≈ 2 ✓ |
| 2 | Stability boundary vs. central difference | NB stable at $\omega\Delta t \ge 3.0$ vs. CD $1.95$ ✓ |
| 3 | Cross-check vs. Newmark/CD in the stable range | $<2.6\%$ ✓ |
| 4 | Spectral damping (light, well-behaved) | ✓ |
| 5 | 1-D wave speed | exact $=\sqrt{E/\rho}$ ✓ |
| 6 | LNVD → static field match | $0.0\%$ ✓ |
| 7 | Rigid-body momentum conservation | exact ✓ |
| 8 | `criticalTimeStep() = ℓ/c` on a 1-element bar | exact ✓ |
| 9 | Energy-recorder element-mass KE | correct ✓ |

PR #2 (recorder + integrators) merged to `ladruno` 2026-05-30; PR #3/#4 (queryable
`dt_cr` + tangent/periodic recompute + aliasing fix) followed. See the
implementation log in [[04_explicit_dynamics_and_energy_balance]].

---

## 9. Quick reference

### Command surface

```
integrator ExplicitBathe      <p=0.54> [-cfl] [-cflAbort] [-tangent]
                              [-recompute N] [-lump rowsum|diagonal]
                              [-verbose] [-divergence f]

integrator ExplicitBatheLNVD  <p=0.54> <alpha=0.8> [-cfl] [-cflAbort] [-tangent]
                              [-recompute N] [-lump rowsum|diagonal]
                              [-verbose] [-divergence f]

criticalTimeStep                         # returns conservative dt_cd = 2/ω_max
```

| Flag | Effect |
|---|---|
| `p` | sub-step / damping parameter, $0<p<1$, default 0.54 ($0.5$ = no dissipation) |
| `alpha` *(LNVD only)* | FLAC local-damping coefficient, $0\le\alpha<1$, default 0.8 |
| `-cfl` / `-criticalTimestep` | enable per-step $\Delta t_{\text{cr}}$ estimation + reporting |
| `-cflAbort` | hard-stop the run if $\Delta t$ exceeds the Noh–Bathe limit |
| `-tangent` | size $\Delta t_{\text{cr}}$ from the current **tangent** stiffness (default: initial) |
| `-recompute N` | refresh $\Delta t_{\text{cr}}$ every $N$ steps (implies `-tangent`) — tracks softening/stiffening |
| `-lump rowsum\|diagonal` | element-mass lumping for the $\Delta t_{\text{cr}}$ pencil (default rowsum) |
| `-verbose` | per-step $\Delta t$ / energy reporting |
| `-divergence f` | abort if the kinetic-energy proxy grows by factor $f$ in one step |

### Key constants & tags

| Symbol | Value | Where |
|---|---|---|
| `EB_NB_STABILITY_FACTOR` | `2.0` | `ExplicitBathe.h` — NB stability advantage over CD |
| default $p$ | `0.54` | both parsers |
| default $\alpha$ (LNVD) | `0.8` | `ExplicitBatheLNVD` parser |
| classTag `ExplicitBathe` | **33000** | `classTags.h` (Ladruno band ≥33000; was 61) |
| classTag `ExplicitBatheLNVD` | **33002** | `classTags.h` (Ladruno band ≥33000; was 63) |
| classTag `EnergyBalanceRecorder` | 26 | `classTags.h` |

---

## 10. References & cross-links

- **Primary:** Gunwoo Noh & Klaus-Jürgen Bathe, *"An explicit time integration
  scheme for the analysis of wave propagations,"* Computers & Structures **129**
  (2013) 178–193, [doi:10.1016/j.compstruc.2013.06.007](https://doi.org/10.1016/j.compstruc.2013.06.007).
- **Stability / spectral analysis:** Hughes, *The Finite Element Method: Linear
  Static and Dynamic Finite Element Analysis* (1987), Ch. 9 (amplification matrix,
  period elongation/amplitude decay). Belytschko et al., *Nonlinear Finite
  Elements for Continua and Structures*, Ch. 6 (explicit stability, damped central
  difference). Bathe, *Finite Element Procedures*, Ch. 9.
- **FLAC local damping:** Itasca FLAC theory (local non-viscous damping for
  dynamic relaxation), $\alpha \approx \pi D$.
- **In-repo ADRs:** [[04_explicit_dynamics_and_energy_balance]] (the source-of-
  record decision log), [[Ladruno_explicit_roadmap]] (the broader explicit-
  dynamics program), [[central_difference_ladruno_guide]] (the non-dissipative
  sibling), [[12_damping_channels]] (how damping enters OpenSees).
- **Source:** [`ExplicitBathe.{h,cpp}`](../SRC/analysis/integrator/),
  [`ExplicitBatheLNVD.{h,cpp}`](../SRC/analysis/integrator/),
  [`CriticalTimeStep.{h,cpp}`](../SRC/analysis/integrator/).
