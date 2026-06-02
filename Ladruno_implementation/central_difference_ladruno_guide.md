---
title: CentralDifferenceLadruno — Theory, Implementation & Use
project: Ladruno
status: reference
priority: high
tags:
  - integrator
  - explicit
  - central-difference
  - leap-frog
  - time-integration
  - dynamics
  - theory
aliases:
  - CentralDifferenceLadruno
  - robust central difference
  - explicit leap-frog
---

# CentralDifferenceLadruno

> [!abstract] What this document is
> A thorough, self-contained reference for **`CentralDifferenceLadruno`**
> (classTag **33003**, Ladruno private band ≥33000): a single, clean **explicit
> leap-frog central-difference** integrator. It is the "central difference done
> right" — with a **correct first step**, a **built-in critical-time-step guard**,
> **clean full-step velocity output**, and **energy-balance discipline** — without
> modifying any upstream OpenSees class. This document covers the continuum/ODE
> **theory**, the **leap-frog derivation**, the **stability and dissipation**
> analysis, the **starter** that the legacy classes get wrong, the **OpenSees
> implementation**, and the **intended use cases**.
>
> Companion: [[explicit_bathe_integrators_guide]] (the dissipative Noh–Bathe
> sibling). Source-of-record ADR: [[05_robust_central_difference]]. Broader
> program: [[Ladruno_explicit_roadmap]] and [[04_explicit_dynamics_and_energy_balance]].

---

## 1. The problem and why a new class

### 1.1 The semi-discrete equation of motion

After FE spatial discretization the structure obeys the coupled second-order ODE
system (the same starting point as every transient integrator):

**Direct form**

$$
\mathbf{M}\,\ddot{\mathbf{u}}(t) + \mathbf{C}\,\dot{\mathbf{u}}(t) + \mathbf{F}^{\text{int}}(\mathbf u(t)) = \mathbf{P}(t)
$$

**Index form**

$$
M_{ab}\,\ddot u_b + C_{ab}\,\dot u_b + F^{\text{int}}_a(\mathbf u) = P_a(t), \qquad a = 1,\dots,n_{\text{dof}}
$$

with $\mathbf M$ the mass, $\mathbf C$ the viscous damping, $\mathbf F^{\text{int}}$
the internal restoring force ($=\mathbf K\mathbf u$ in the linear case), and
$\mathbf P$ the external load. The physics is Newton's second law assembled over
the mesh: *inertia + dissipation + restoring = applied*.

### 1.2 The explicit central-difference idea

Central difference is the oldest and simplest explicit scheme. The classical
"two-level" form approximates the time derivatives of $\mathbf u$ with centered
finite differences about $t_n$:

$$
\dot{\mathbf u}_n \approx \frac{\mathbf u_{n+1} - \mathbf u_{n-1}}{2\,\Delta t},
\qquad
\ddot{\mathbf u}_n \approx \frac{\mathbf u_{n+1} - 2\mathbf u_n + \mathbf u_{n-1}}{\Delta t^{2}} .
$$

Substituting into the equation of motion at $t_n$ and solving for
$\mathbf u_{n+1}$ gives an update that needs only **already-known** quantities at
$t_n$ and $t_{n-1}$, with the mass (and, in the classical form, the damping) on
the left. If $\mathbf M$ is lumped/diagonal the solve is a trivial $\mathbf M^{-1}$
division — the hallmark of an explicit method.

### 1.3 The leap-frog (velocity-staggered) form

`CentralDifferenceLadruno` implements the equivalent but more robust **leap-frog**
("velocity Verlet" / staggered) form, in which the velocity lives at the
**half-steps** $t_{n\pm1/2}$ and the displacement at the integer steps:

$$
\boxed{\;
\begin{aligned}
\mathbf a_n &= \mathbf M^{-1}\big(\mathbf P_n - \mathbf C\,\mathbf v_{n-1/2} - \mathbf F^{\text{int}}(\mathbf u_n)\big) &&\text{(mass-only solve; damping lagged at the known half-step velocity)}\\[2pt]
\mathbf v_{n+1/2} &= \mathbf v_{n-1/2} + \Delta t\,\mathbf a_n &&\text{(half-step velocity kick)}\\[2pt]
\mathbf u_{n+1} &= \mathbf u_n + \Delta t\,\mathbf v_{n+1/2} &&\text{(full-step drift)}
\end{aligned}
\;}
$$

This is exactly **Newmark with $\beta = 0$, $\gamma = \tfrac12$** — the
non-dissipative explicit member of the Newmark family. Compared with the
two-level form, the leap-frog form:

- keeps the velocity as a **first-class primary state** (cleaner damping and
  momentum handling),
- needs only **one** mass-only solve per step,
- treats viscous damping by **lagging** it at the *known* half-step velocity
  $\mathbf v_{n-1/2}$ — which keeps the scheme genuinely explicit (no $\mathbf C$
  on the LHS), the LS-DYNA-style treatment.

### 1.4 Why this is a *new* class, not a patch

OpenSees already ships **six** explicit CD-family classes, and **every one of
them is missing something** (ADR [[05_robust_central_difference]] §Why):

| Class | tag | First step | Damping | `dt_cr` | Velocity out | Mass |
|---|---|---|---|---|---|---|
| `CentralDifference` | 5 | warns `Ut-1 = Ut` ⚠️ | Rayleigh flags stored but **unused** | none | full-step | lumped |
| `CentralDifferenceAlternative` | 17 | V from committed | none | none | half-step | lumped |
| `CentralDifferenceNoDamping` | 18 | V from committed | none | none | half-step | lumped |
| `ExplicitDifference` | 55 | warns `Ut-1 = Ut` ⚠️ | Rayleigh via residual (lagged) ✓ | none | half-step → crude reconstruct | lumped |
| `ExplicitDifferenceStatic` | 62 | — | FLAC local | none | — | lumped |
| `NewmarkExplicit` | 19 | clean (2-level) ✓ | **C on LHS** (coupled/implicit-damped) | none | full-step | any |

The genuinely-missing combination is **(1) a correct explicit starter + (2) a
`dt_cr` guard + (3) clean full-step velocity output + (4) energy-balance
discipline** in one *truly explicit leap-frog* scheme. The Noh–Bathe work already
built the infrastructure for (2) and (4) (`CriticalTimeStep`,
`EnergyBalanceRecorder`); this class brings the leap-frog scheme up to that
standard. Because the **sibling-fork policy forbids patching upstream classes in
place**, it lands as a new fork class (classTag 33003) rather than a fix to
`ExplicitDifference`.

> [!note] What it is deliberately NOT
> The **coupled / implicit-damped** central-difference scheme (damping
> $\mathbf C$ on the LHS, consistent-mass with a real solve) **already exists** as
> `NewmarkExplicit(0.5)` — an adversarial sweep confirmed a coupled mode would
> have been bit-identical to it (default $\gamma=0.5$, `addCtoTang(c2); addMtoTang(c3)`,
> no starter defect). So that mode was **dropped**: for implicit-damped CD, the
> answer is documentation — `integrator NewmarkExplicit 0.5` — not new code. This
> class is the *explicit leap-frog* scheme **only** (decisions C2/C7).

---

## 2. Theory

### 2.1 Accuracy and dissipation ($\gamma = \tfrac12$)

The leap-frog scheme is **second-order accurate** and, with $\gamma=\tfrac12$, has
**exactly zero algorithmic (numerical) dissipation** — the amplitude of every
resolved mode is preserved. Its only frequency error is a slight **period
*shortening*** of order

$$
\frac{\Delta T}{T} \;\approx\; -\frac{\Omega^{2}}{24}, \qquad \Omega \equiv \omega\,\Delta t,
$$

i.e. the discrete frequency is slightly **higher** than the true one (the
*opposite* sign to the trapezoidal rule's period *elongation*). This is verified
in the ADR (theory item T6) and is the canonical Hughes Ch. 9 result.

> [!info] Zero dissipation is a double-edged sword
> **Good for** wave propagation and energy-conservation studies — no artificial
> energy bleed. **Bad for** un-damped spurious high-frequency mesh modes — they
> ring forever. If you need controllable high-frequency dissipation, reach for the
> dissipative sibling [[explicit_bathe_integrators_guide|ExplicitBathe]] instead.
> This trade-off is stated right in the class header.

### 2.2 Stability — the exact $2/\omega_{\max}$ limit

Modal decomposition reduces the linear system to independent scalar SDOF
oscillators. For the undamped leap-frog scheme the amplification factor has unit
modulus (no dissipation) precisely while

$$
\boxed{\;\Delta t \;<\; \Delta t_{\text{cr}} \;=\; \frac{2}{\omega_{\max}}\;}
$$

and the scheme goes unstable for $\Delta t > 2/\omega_{\max}$. **There is no
Noh–Bathe $2\times$ bonus here** — the central-difference stability limit is
*exactly* $2/\omega_{\max}$, so the integrator hard-codes **stability factor 1.0**
(`CentralDifferenceLadruno.cpp`, comment at line 61). $\omega_{\max}$ is the
highest *element* natural frequency from the per-element eigenproblem (§4).

For a 1-D rod of element length $\ell$ and wave speed $c=\sqrt{E/\rho}$ the limit
is the **CFL condition** $\Delta t < \ell/c$ — the time for a wave to cross one
element. Lumped mass gives $\omega_{\max} = 2c/\ell$ (limit $\ell/c$) whereas
consistent mass gives $\omega_{\max} = 2\sqrt3\,c/\ell$ (limit $\ell/\sqrt3\,c$) —
another reason lumped mass is preferred: it **raises** the stable step (ADR T5).

### 2.3 Damping reduces the critical step

With modal damping ratio $\xi$ the stable step shrinks to

$$
\Delta t_{\text{cr}} \;=\; \frac{2}{\omega_{\max}}\Big(\sqrt{1+\xi^{2}} - \xi\Big)
$$

(LS-DYNA Eq. 24.29; standard in Belytschko Ch. 6 and Hughes Ch. 9). This is the
formula the shared eigensolver applies per element (`CriticalTimeStep.cpp:255`),
and the integrator reports the **damped** value whenever damping is active. In an
explicit scheme, **damping costs you step size** — the opposite of the implicit
intuition.

### 2.4 The βK trap

For Rayleigh damping $\mathbf C = \alpha\mathbf M + \beta\mathbf K$ the modal ratio
is

$$
\xi_i \;=\; \frac{\alpha}{2\,\omega_i} \;+\; \frac{\beta\,\omega_i}{2}.
$$

- The **mass-proportional** term $\alpha/2\omega$ **decreases** with frequency, so
  at $\omega_{\max}$ it is benign — it barely moves $\Delta t_{\text{cr}}$.
- The **stiffness-proportional** term $\beta\omega/2$ **grows** with frequency, so
  at $\omega_{\max}$ the factor $(\sqrt{1+\xi^2}-\xi) \to 1/(2\xi)$ and

$$
\Delta t_{\text{cr}} \;\longrightarrow\; \frac{2}{\beta\,\omega_{\max}^{2}} \quad\text{— a *quadratic collapse* of the stable step.}
$$

LS-DYNA declares classical $\beta\mathbf K$ damping "impractical" for explicit
analysis for exactly this reason. The integrator **detects this data-driven** —
when the damped $\Delta t_{\text{cr}}$ falls below `CDL_BETAK_WARN_RATIO = 0.90`
times the undamped one — and emits a one-time warning recommending
mass-proportional $\alpha\mathbf M$. It does **not** refuse; `-cflAbort` is the
hard guard (decision C6).

```cpp
// CentralDifferenceLadruno.cpp — βK trap detection (paraphrased)
if (damped_dt < 0.90 * undamped_dt) {
    opserr << "NOTE damping is significantly reducing the critical time step ...\n"
              "  stiffness-proportional (betaK) Rayleigh collapses dt_cr ~ 2/(betaK*omega_max^2);\n"
              "  prefer mass-proportional alphaM ...\n";
}
```

### 2.5 The first-step starter — the correctness fix

This is the single biggest reason the class exists. The leap-frog recursion needs
a velocity at the **back half-step** $t_{-1/2}$ to begin, but the user supplies
$\mathbf v_0$ at the integer step $t_0$. The correct starter is

$$
\mathbf a_0 = \mathbf M^{-1}\big(\mathbf P_0 - \mathbf C\,\mathbf v_0 - \mathbf F^{\text{int}}(\mathbf u_0)\big),
\qquad
\boxed{\;\mathbf v_{-1/2} = \mathbf v_0 - \tfrac12\,\Delta t\,\mathbf a_0\;}
$$

The legacy classes punt on this — `CentralDifference::domainChanged()` literally
prints `"WARNING: assuming Ut-1 = Ut"`, which is equivalent to seeding
$\mathbf v_{-1/2} = \mathbf v_0$ (or $\mathbf u_{-1} = \mathbf u_0$). That silently
**loses an order of accuracy on the first step and double-counts initial
conditions** whenever $\mathbf v_0$ or $\mathbf a_0$ is nonzero (e.g. a load
already applied at rest, or a prescribed initial velocity).

> [!important] Why the starter is on the *first step*, not in `domainChanged()`
> Computing $\mathbf a_0$ requires (a) a time step $\Delta t$ — which is **not
> known** in `domainChanged()` — and (b) a **formed and factored SOE**, which does
> not exist there either. That is *exactly* why the legacy classes punt. So
> `CentralDifferenceLadruno` defers the starter to the **first `newStep()`**
> (behind a `firstStep` flag), where it forms its own tangent, does **one extra
> mass-only solve** at the committed configuration to get $\mathbf a_0$, then seeds
> $\mathbf v_{-1/2} = \mathbf v_0 - \tfrac12\Delta t\,\mathbf a_0$ before the
> normal advance (decisions C3/B1). This mirrors `ExplicitBathe`'s
> deferred-first-solve pattern and is the fix the legacy classes never made.

### 2.6 Two velocities, kept separate

A subtle but important design point (decision C5). The half-step staggering means
there are **two distinct velocities**, and conflating them is a classic bug:

1. **The lagged half-step velocity** $\mathbf v_{n-1/2}$ (`getVel()`). At solve
   time this is the *known* velocity, and it is what the residual's damping force
   is evaluated against. `IncrementalIntegrator::addModalDampingForce` assembles
   modal damping into the residual via `setB`, and the element Rayleigh force is
   likewise built at the velocity set on the nodes — so `getVel()` **must** return
   the lagged half-step value. `EnergyBalanceRecorder`'s damping-work term `DW`
   integrates against this same lagged velocity for energy closure.

2. **The clean full-step velocity** $\mathbf v_n$ (node/recorder output). The
   physically meaningful velocity *at the integer step* is the average of the two
   surrounding half-steps:

$$
\mathbf v_{n} \;=\; \tfrac12\big(\mathbf v_{n-1/2} + \mathbf v_{n+1/2}\big) \;=\; \mathbf v_{n-1/2} + \tfrac12\,\Delta t\,\mathbf a_n
\;\equiv\; \frac{\mathbf u_{n+1} - \mathbf u_{n-1}}{2\,\Delta t}.
$$

   This centered, second-order-accurate value is pushed to the nodes via
   `setResponse`, **fixing the legacy half-step / crude-reconstruct output** of
   `ExplicitDifference`. The beauty is that it is available **immediately after the
   solve** (no future value needed) because it equals
   $\mathbf v_{n+1/2} + \tfrac12\Delta t\,\mathbf a_{n+1}$ in the implementation's
   advance-then-solve bookkeeping.

---

## 3. What drove the implementation

| Driver | Decision | Why |
|---|---|---|
| **Genuinely missing combination** | New fork class `CentralDifferenceLadruno`, classTag 33003 | The four desiderata (correct starter, `dt_cr`, clean velocity, energy discipline) exist in *no* single upstream class; sibling-fork policy forbids patching them in place (C1). |
| **Avoid redundancy** | **No** coupled/`-damping` mode | The sweep proved a coupled mode == `NewmarkExplicit(0.5)`; one classTag, one behavior (C2/C7). |
| **Correct first step** | Starter on first `newStep()`, not `domainChanged()` | No `dt`, no factored SOE in `domainChanged()`; the legacy `assuming Ut-1=Ut` bug is the thing being fixed (C3/B1). |
| **Mass-only LHS** | `formEleTangent`: `zeroTangent(); addMtoTang();` | Trivial diagonal $\mathbf M^{-1}$; damping enters the residual at the lagged half-step velocity, never the LHS. |
| **Single solve/step** | `update()` guarded at `updateCount > 1` ⇒ `algorithm Linear` | CD is one solve/step (unlike ExplicitBathe's two). A nonlinear algorithm would re-enter `update()` and corrupt the leap-frog state. |
| **`dt_cr` valid before `analyze`** | Computed **unconditionally in `domainChanged()`** | `CriticalTimeStep` runs its **own** LAPACK eigensolve (it needs neither `dt` nor the global SOE), so `criticalTimeStep()` can be queried before the first `analyze` — dropping the `analyze(1,1e-9)` priming hack the Bathe classes needed (C4). |
| **βK is a trap** | Data-driven warn (`damped_dt < 0.9·undamped_dt`) + report damped value | Domain exposes no Rayleigh getter, so the collapse signature itself is the trigger; warn, don't refuse (C6). |
| **Silent-wrongness of explicit** | Built-in `dt_cr` + NaN/Inf abort + `-divergence` KE guard + energy-recorder ready | No Newton residual to warn you; the diagnostics *are* the deliverable. |
| **Clean full-step output** | Separate `Vfull` buffer pushed via `setResponse`; `getVel()` returns lagged `Vhalf` | Fixes the legacy half-step output bug; keeps the modal-damping hook consistent (C5). |

---

## 4. The critical-time-step machinery (shared)

`CentralDifferenceLadruno` reuses the **same** `computeCriticalTimeStep()`
eigensolver as the Bathe integrators
([`CriticalTimeStep.cpp`](../SRC/analysis/integrator/CriticalTimeStep.cpp)) — see
[[explicit_bathe_integrators_guide#4. The critical-time-step machinery (CriticalTimeStep)]]
for the full treatment. In brief:

- Per-element generalized eigenproblem
  $\mathbf K^{(e)}\boldsymbol\phi = \lambda\,\mathbf M^{(e)}\boldsymbol\phi$;
  governing step $= \min_e 2/\sqrt{\lambda_{\max}^{(e)}}$ (global $\omega_{\max}$
  is bounded by the per-element max).
- **`DSYGVX`** (symmetric-definite, largest eigenvalue only) with a **`DGGEV`**
  fallback and a relative-β threshold for indefinite/singular pencils. `DSYGVX`
  (not `DSYGV`) because the bundled Linux reference LAPACK ships `dsygvx.f` only.
- **Mass copied before stiffness** (shared-static-matrix aliasing fix, D8).
- Per-element Rayleigh factors → damped step
  $\frac{2}{\omega_{\max}}(\sqrt{1+\xi^2}-\xi)$.
- `MPI_Allreduce(MPI_MIN)` across ranks under a parallel build.

**Differences specific to this class:**

- **Default lumping is `diagonal`** (diagonal-of-consistent), not `rowsum` — the
  explicit default is chosen to be robust for the rotational DOFs of beams/shells,
  where row-sum lumping can give near-zero (non-conservative) diagonal mass.
- **Stability factor 1.0** — the reported `getCriticalTimeStep()` is the exact CD
  limit (the damped value when damping reduces it, else $2/\omega_{\max}$), with
  **no** Noh–Bathe $2\times$ headroom.

```cpp
// cdl_reported(): pick the binding limit
//   undamped == 2/omega_max ; damped applies the (sqrt(1+xi^2)-xi) reduction.
//   report damped when damping actually reduces the step, else undamped.
static inline double cdl_reported(double damped, double undamped) {
    if (undamped == inf || undamped <= 0.0) return -1.0;       // not applicable
    if (damped != inf && damped > 0.0 && damped < undamped) return damped;
    return undamped;
}
```

> [!caution] Inherited caveats
> `dt_cr` ignores constraints (`equalDOF`/rigid diaphragms/MP) and pure nodal
> mass; a pure-nodal-mass model reports a non-applicable (`≤0`) value. For
> rotational DOFs use `-lump diagonal` (already the default here).

---

## 5. OpenSees implementation

`CentralDifferenceLadruno : public TransientIntegrator` lives in
[`SRC/analysis/integrator/CentralDifferenceLadruno.{h,cpp}`](../SRC/analysis/integrator/),
classTag `INTEGRATOR_TAGS_CentralDifferenceLadruno = 33003`.

### 5.1 The life-cycle (one solve per step)

```
analyze(nSteps, dt)
└── for each step:
      newStep(dt)
      │   ├─ [first step only] starter:
      │   │     setResponse(u0, v0, 0); formTangent; formUnbalance; solve
      │   │     → a0 = M⁻¹(P0 − C v0 − F_int(u0));  v_{-1/2} = v0 − ½dt·a0
      │   ├─ leap-frog advance:  v_{n+1/2} = v_{n-1/2} + dt·a_n
      │   │                      u_{n+1}   = u_n + dt·v_{n+1/2}
      │   └─ setResponse(u_{n+1}, v_{n+1/2}, 0); updateDomain(t+dt, dt)
      ├─ [algorithm Linear] formTangent → M only (formEleTangent)
      ├─ formUnbalance   ── R = P − C·v_{n+1/2} − F_int(u_{n+1})   (accel set to 0 ⇒ no inertia term)
      ├─ solve           ── M·a_{n+1} = R   (diagonal ⇒ a = M⁻¹R)
      ├─ update(a)  [1×] ── NaN/Inf + KE circuit breakers; clean full-step
      │                     v_{n+1} = v_{n+1/2} + ½dt·a_{n+1}; setResponse;
      │                     carry a_{n+1} as a_n for the next step
      └─ commit()        ── commitDomain()
```

The defining structural choices:

- **`Azero`** — a persistent zero vector. Before each solve the integrator sets
  the nodal acceleration to zero via `setResponse(Ut, Vhalf, Azero)` so the
  residual carries **no inertia term** ($\mathbf M\mathbf a$ is *not* subtracted)
  — the mass appears only on the LHS, and the solve directly yields
  $\mathbf a_{n+1} = \mathbf M^{-1}\mathbf r$.
- **Advance-then-solve** ordering: `newStep()` advances $\mathbf u$ to
  $\mathbf u_{n+1}$ and sets the trial state; the framework then solves for
  $\mathbf a_{n+1}$; `update()` finalizes. This records $\mathbf u_{n+1}$ at
  $t_{n+1}$, so $N$ analyze steps → $N$ records (no off-by-one).

### 5.2 Key methods, annotated

**`formEleTangent` / `formNodTangent` — mass-only LHS:**

```cpp
int CentralDifferenceLadruno::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();   // ONLY the mass — damping is NOT on the LHS
    return 0;
}
```

**`domainChanged()`** — allocates the five state vectors
(`Ut, Vhalf, Aprev, Vfull, Azero`), seeds $\mathbf u_n,\mathbf v_{n-1/2}(=\mathbf v_0),
\mathbf a_n$ from the committed DOF state, sets `firstStep = true`, and **runs the
`dt_cr` eigensolve unconditionally** so `criticalTimeStep()` is valid before any
`analyze`. It does **not** solve for $\mathbf a_0$ (no `dt`, no factored SOE).

**`newStep(dt)`** — on the first call, runs the starter (§2.5); then does the
leap-frog advance `Vhalf += dt·Aprev; Ut += dt·Vhalf`, sets the trial response
with zero accel, advances the clock. Optionally refreshes/aborts on `dt_cr`.

**`update(U)`** — `U` is the solved $\mathbf a_{n+1}$. Runs the NaN/Inf abort,
forms the clean full-step velocity `Vfull = Vhalf + 0.5·dt·U`, pushes
$(\mathbf u_{n+1}, \mathbf v_{n+1}, \mathbf a_{n+1})$ to the nodes, runs the
optional kinetic-energy runaway guard, and **carries `Aprev = U`** for the next
step's leap-frog kick. Guarded at `updateCount > 1` (exactly one solve/step).

**`getVel()`** — returns `Vhalf` (the lagged half-step velocity), the
modal-damping / residual-evaluation hook (§2.6).

**`getCriticalTimeStep()`** — returns `cdl_reported(damped, undamped)` (§4).

**`commit()`** — time was already advanced inside `newStep`/`update`; just calls
`commitDomain()`.

**`sendSelf`/`recvSelf`** — serialize the seven run-control scalars
(`compute_critical_timestep, verbose, cflAbort, divergenceFactor, cflUseTangent,
cflRecomputeEvery, lumping`) for `openseessp`/database restart.

### 5.3 State vectors

| Vector | Meaning |
|---|---|
| `Ut` | displacement $\mathbf u_n$ (advanced to $\mathbf u_{n+1}$ in `newStep`) |
| `Vhalf` | **half-step** velocity $\mathbf v_{n-1/2}$ — primary state; what `getVel()` returns |
| `Aprev` | acceleration $\mathbf a_n$ carried between steps (drives the velocity kick) |
| `Vfull` | **full-step** output velocity $\mathbf v_n = \tfrac12(\mathbf v_{-}+\mathbf v_{+})$ |
| `Azero` | persistent zero vector → sets nodal accel to 0 so the residual has no inertia term |

### 5.4 Registration (8 sites, no upstream class touched)

`classTags.h` (33003, Ladruno band ≥33000); `FEM_ObjectBrokerAllClasses.cpp`
(+include +case); `runtime/runtime/TclPackageClassBroker.cpp` (+include +case);
`interpreter/OpenSeesCommands.{h,cpp}` (forward-decl + string dispatch to
`OPS_CentralDifferenceLadruno()`); `tcl/commands.cpp` (extern + `integrator`
branch); integrator `CMakeLists.txt` + `Makefile`. *(Python `integrator` dispatch
is string-based, so `PythonWrapper.cpp` needs no edit.)*

---

## 6. Intended use cases

> [!success] Reach for `CentralDifferenceLadruno` when…
> - You need a **second-order, energy-conserving, zero-dissipation** explicit
>   scheme — **wave propagation**, **momentum-preserving** dynamics, or any study
>   where you must *not* artificially bleed energy.
> - You have **nonzero initial conditions** (prescribed $\mathbf v_0$, or a load
>   applied at rest needing a consistent $\mathbf a_0$) — the **correct starter**
>   reproduces the analytical trajectory **from step 1**, where every legacy CD
>   class is wrong.
> - You want a **built-in `dt_cr` guard** queryable *before* the first `analyze`,
>   with optional `-cflAbort` and tangent/periodic refresh for nonlinear runs.
> - You want **clean full-step velocity output** at the recorder, not the legacy
>   half-step / crude-reconstruct value.

**Required recipe:**

```python
ops.mass(...)                                   # lumped/diagonal element or nodal mass
ops.system('Diagonal')                          # trivial M⁻¹
ops.algorithm('Linear')                         # exactly ONE solve/step
ops.integrator('CentralDifferenceLadruno', '-cflAbort', '-lump', 'diagonal')
ops.analysis('Transient')

import math
dt_cr = ops.criticalTimeStep()                  # valid even before the first analyze (C4)
n = max(1, math.ceil(dt / (0.9 * dt_cr)))       # 0.9 ≈ LS-DYNA TSSFAC
ops.analyze(n, dt/n)                             # driver-level sub-stepping
```

> [!tip] Decision guide — which explicit integrator?
> | Want | Use |
> |---|---|
> | Zero dissipation, energy/momentum conservation, waves | **`CentralDifferenceLadruno`** |
> | Controllable high-frequency dissipation (kill spurious modes) | [[explicit_bathe_integrators_guide\|`ExplicitBathe 0.54`]] |
> | Quasi-static / dynamic relaxation (snap-through, settling) | [[explicit_bathe_integrators_guide\|`ExplicitBatheLNVD 0.54 0.8`]] |
> | Coupled / implicit-damped CD (C on LHS, consistent mass) | `NewmarkExplicit 0.5` |

### 6.1 When *not* to use it

- **Spurious high modes ring and pollute the result** → you want dissipation; use
  `ExplicitBathe`.
- **Stiffness-proportional ($\beta\mathbf K$) Rayleigh damping** → it collapses the
  step quadratically (§2.4); switch to mass-proportional $\alpha\mathbf M$.
- **Stiff/quasi-static where implicit converges fine** → an unconditionally stable
  implicit scheme takes far larger steps.
- **Consistent mass + sparse solver** → defeats the explicit premise; use lumped
  mass + `system Diagonal`.

---

## 7. Validation

The acceptance battery (`Ladruno_scripts/_verify_explicit.py`, **CDL-1…10**,
ported to the Zone-A pytest `tests/test_centralDifferenceLadruno_integrator.py`):

| # | Check | Expected |
|---|---|---|
| 1 | Order of accuracy (SDOF free-vibration log–log slope) | ≈ **2.0** |
| 2 | Stability boundary | stable for $\omega\Delta t < 2.0$, unstable just above |
| 3 | Damped SDOF decay vs. $e^{-\xi\omega t}$ + damped period, largest stable step $=\tfrac{2}{\omega}(\sqrt{1+\xi^2}-\xi)$ | match at $\xi=0.02,0.1,0.5$ |
| 4 | Cross-check vs. Newmark/CD in the stable range | $<$ a few % |
| 5 | 1-D wave speed | exact $=\sqrt{E/\rho}$ |
| 6 | Rigid-body momentum (zero stiffness) | conserved exactly |
| 7 | `criticalTimeStep() = ℓ/c` on a 1-element bar | exact |
| 8 | **First-step correctness** (nonzero $\mathbf v_0/\mathbf a_0$ reproduces the analytic trajectory from step 1) | ✓ — *the test every legacy class fails* |
| 9 | **βK trap + guard**: `getCriticalTimeStep()` drops $\sim$quadratically as $\beta_K$ rises; the C6 warn fires; $\alpha_M$ stays benign | ✓ |
| 10 | Energy closure ≈ 1% via `EnergyBalanceRecorder` (multi-DOF) | ✓ |

Merged to `ladruno` via **PR #22** (the coupled mode was dropped during the
adversarial sweep). The CDL-1…10 battery was later ported into the Zone-A CI
(numpy-added, CI green). Implementation log in [[05_robust_central_difference]].

---

## 8. Quick reference

### Command surface

```
integrator CentralDifferenceLadruno  [-cfl] [-cflAbort] [-tangent]
                                     [-recompute N] [-lump rowsum|diagonal]
                                     [-verbose] [-divergence f]

criticalTimeStep                      # exact CD limit 2/ω_max (or damped value)
```

**No positional argument and no `-damping` flag** — this is one explicit scheme.
(For coupled/implicit-damped CD use `integrator NewmarkExplicit 0.5`.)

| Flag | Effect |
|---|---|
| `-cfl` / `-criticalTimestep` | enable per-step $\Delta t_{\text{cr}}$ reporting (the value is *always* computed in `domainChanged`) |
| `-cflAbort` | hard-stop if $\Delta t$ exceeds the CD limit (stability factor 1.0) |
| `-tangent` | size $\Delta t_{\text{cr}}$ from the current **tangent** stiffness |
| `-recompute N` | refresh $\Delta t_{\text{cr}}$ every $N$ steps (implies `-tangent`) |
| `-lump rowsum\|diagonal` | element-mass lumping (**default diagonal** for rotational-DOF robustness) |
| `-verbose` | per-step $\Delta t$ / max-acceleration reporting |
| `-divergence f` | abort if the kinetic-energy proxy grows by factor $f$ in one step |

### Key facts

| Item | Value |
|---|---|
| classTag | **33003** (Ladruno private band ≥33000) |
| Scheme | Newmark $\beta=0,\ \gamma=\tfrac12$ — explicit leap-frog |
| Order / dissipation | 2nd order, **zero** algorithmic dissipation (period shortening $\approx -\Omega^2/24$) |
| Stability | $\Delta t < 2/\omega_{\max}$ (factor **1.0**, no NB $2\times$) |
| Damped stability | $\tfrac{2}{\omega_{\max}}(\sqrt{1+\xi^2}-\xi)$ |
| Solves / step | **1** (⇒ `algorithm Linear`) |
| Starter | $\mathbf v_{-1/2} = \mathbf v_0 - \tfrac12\Delta t\,\mathbf a_0$ on first `newStep` |
| βK warn ratio | `CDL_BETAK_WARN_RATIO = 0.90` |
| Default lumping | `diagonal` |

---

## 9. References & cross-links

- **Stability / spectral theory:** Hughes, *The Finite Element Method: Linear
  Static and Dynamic Finite Element Analysis* (1987), Ch. 9 (amplification matrix,
  period elongation/shortening, damped stability). Belytschko et al., *Nonlinear
  Finite Elements for Continua and Structures*, Ch. 6 (central-difference
  leap-frog box, undamped $2/\omega_{\max}$, damped $(2/\omega)(\sqrt{1+\xi^2}-\xi)$,
  consistent-vs-lumped $2c/\ell$ vs $2\sqrt3 c/\ell$). Bathe, *Finite Element
  Procedures*, Ch. 9 (the starter $^{-1}\!\mathbf U$, $\Delta t_{\text{cr}}=T_n/\pi$).
  Ibrahimbegovic (2009), Ch. 6 (CD ≡ Newmark $\beta=0,\gamma=\tfrac12$).
- **Explicit / CFL / damping reduction:** LS-DYNA Theory §24.2 (leap-frog update),
  §24.3 Eq. 24.29 (damped stability), §22.1 (solid `dt_cr`), `*CONTROL_TIMESTEP`
  TSSFAC, `*DAMPING_PART_STIFFNESS` (βK impractical for explicit).
- **In-repo ADRs:** [[05_robust_central_difference]] (the source-of-record
  decision log + implementation log), [[04_explicit_dynamics_and_energy_balance]]
  (the shared `dt_cr` + energy machinery), [[Ladruno_explicit_roadmap]],
  [[explicit_bathe_integrators_guide]] (the dissipative sibling),
  [[12_damping_channels]] (how damping enters OpenSees).
- **Source:**
  [`CentralDifferenceLadruno.{h,cpp}`](../SRC/analysis/integrator/),
  [`CriticalTimeStep.{h,cpp}`](../SRC/analysis/integrator/); contrast
  [`CentralDifference.cpp`](../SRC/analysis/integrator/CentralDifference.cpp)
  (the legacy `assuming Ut-1=Ut` starter) and `NewmarkExplicit.cpp` (the
  coupled-CD equivalent users are pointed to).
