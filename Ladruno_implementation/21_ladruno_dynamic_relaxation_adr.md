---
title: ADR — Quasi-static dynamic relaxation (LadrunoDynamicRelaxation)
project: Ladruno
status: proposed — scoped, no code
priority: medium
owner: nmora
tags:
  - implementation
  - integrator
  - transient
  - quasi-static
  - dynamic-relaxation
  - kinetic-damping
  - softening
  - adr
---

# ADR — Quasi-static dynamic relaxation (`LadrunoDynamicRelaxation`)

**Status:** proposed — design scoped, **no code yet** (sibling of ADR-20
`LadrunoArcLength`) · **Registry:** `TransientIntegrator` · **classTag:**
`INTEGRATOR_TAGS_LadrunoDynamicRelaxation = 33005` (next free in the ladruno
integrator band; siblings ExplicitBathe=33000 … CentralDifferenceLadruno=33003,
LadrunoArcLength=33004) · **Supersedes nothing** — adds a fork integrator
alongside stock `Newmark`/`CentralDifference` · **Oracle:** the transient
dynamic-relaxation run already standing in `tests/_proto_arch_snapthrough.py`
(`run_dynamic_relaxation`, AL-6 in ADR-20; far-branch `uy ≈ −0.217` under a
1.15×crit load) — `LadrunoDynamicRelaxation` *promotes that oracle to a
first-class integrator*.

> [!note] As-designed scope (v1 target)
> One self-contained **`TransientIntegrator` leaf** that drives a model to **static
> equilibrium by explicit, matrix-free dynamic relaxation**: a fictitious-mass
> leap-frog march with **kinetic (Cundall) damping** that bleeds off kinetic energy
> at each energy peak, terminating at static rest. It runs under the **stock
> `DirectIntegrationAnalysis` driver** with `algorithm Linear` + `system Diagonal`,
> exactly as `CentralDifferenceLadruno` already does — **no analysis-core surgery**.
> Per pseudo-step the SOE solve is a trivial diagonal `M*⁻¹` apply; **K is never
> assembled or factorized in the stepping loop**. It is the **implicit-free**
> counterpart to `LadrunoArcLength`: where the arc-length integrator clears a limit
> point by regularizing `K_T`, this integrator sidesteps the singular tangent
> entirely by never forming one.
>
> **Deliverable = integrator (33005) + a thin `Ladruno_scripts` relax-to-static
> helper.** The stock transient driver has **no convergence early-exit**
> (verified — see decision 5), so static-rest detection is **script-owned** (a
> Python/Tcl loop polling the integrator's residual/KE query), mirroring
> `LadrunoArcLength`'s Layer-B script-owned retry. This is a leaf, zero-core-touch
> property, not a limitation we can engineer away without core surgery.
>
> **Hard constraint honored:** `DirectIntegrationAnalysis::analyze` and
> `StaticAnalysis::analyze` are **not touched, not cloned, not subclassed.** The
> only vanilla edit is the standard integrator-registration row in the broker +
> Tcl/Python dispatch — the same 8-site seam every ladruno integrator already uses
> (ledgered in `LEDGER_vanilla_files.md`).

---

## 1. Context

### 1.1 The phenomenon, and what is / isn't in scope

Dynamic relaxation (DR) reaches a **static** equilibrium by integrating a
**fictitious damped transient** to rest:

$$\mathbf M^{*}\ddot{\mathbf u} + \mathbf C^{*}\dot{\mathbf u} + \mathbf f_\text{int}(\mathbf u) = \mathbf f_\text{ext},$$

with `M*`, `C*`, and the pseudo-time step chosen for fastest decay, **not**
physical accuracy. When the transient comes to rest `ü = u̇ = 0` and the surviving
equation is the static equilibrium `f_int(u) = f_ext`. DR is the classic
matrix-free engine for problems where forming/factorizing `K_T` is undesirable or
impossible: severe softening, contact, form-finding, very large models, and
**limit points / snap-through** — the same payload ADR-20 targets, approached from
the explicit side.

OpenSees ships the transient machinery (`Newmark`, `CentralDifference`) and the
fork already ships two explicit relatives — `CentralDifferenceLadruno` (33003, a
true leap-frog explicit integrator) and `ExplicitDifferenceStatic` (33001, a
quasi-static explicit scheme that applies **user-supplied Rayleigh damping**
`alphaM`/`betaK` via `setRayleighDampingFactors`, verified at
`ExplicitDifference.cpp:154-156`). Neither is **Cundall kinetic-damping** DR with a
**scale-free fictitious mass** and a **clean static-rest termination**. This ADR
adds exactly that.

**Why a new integrator and not just `ExplicitDifferenceStatic`?** That one makes
the user *guess* the Rayleigh constants — the classic DR failure mode (too little →
endless ringing, too much → crawling convergence, and the right value drifts as the
structure softens). This ADR's whole differentiator is **auto-tuning**: a
**Gershgorin scale-free fictitious mass** + **parameter-free Cundall kinetic
damping** + **automatic static-rest termination**, so it is a one-line `integrator
LadrunoDynamicRelaxation` with *nothing to tune*. The auto-tuning *is* the
contribution; the explicit march itself is shared (and ~70% reused, §4).

**Industry siblings.** DR is the one place where the cleaner sibling is
**LS-DYNA**, not Abaqus. Recall the fork's mapping: `LadrunoArcLength ↔ Abaqus
*STATIC, RIKS`; `-stabilize ↔ *STATIC, STABILIZE`. A true *auto-tuned* DR solver,
though, is **`*CONTROL_DYNAMIC_RELAXATION`** (LS-DYNA's named DR — gravity/bolt
preload, form-finding) — **Abaqus exposes no named Cundall-DR feature.** Abaqus
meets the need two other ways: **Abaqus/Explicit quasi-static** (`*DYNAMIC,
EXPLICIT` + mass scaling, policed by `ALLKE/ALLIE`) — same physics, but the user
hand-tunes mass scaling and watches the energy ratio; and **Abaqus/Standard
quasi-static implicit dynamics** (`*DYNAMIC, application=QUASI-STATIC`, HHT with
numerical dissipation) — an *implicit* cousin of `-stabilize`, not of DR. So
`LadrunoDynamicRelaxation` is an **LS-DYNA sibling and an Abaqus gap-filler**: it
packages the auto-tuned DR algorithm neither Abaqus path exposes, in the
form-finding / preload / softening niche where DR dominates.

**In scope (v1):**
- A `TransientIntegrator` leaf running matrix-free explicit DR under the stock
  transient driver.
- A **Gershgorin row-sum fictitious lumped mass** `M*` built once in
  `domainChanged()` (scale-free / optimal-convergence), with `lumped` and `unity`
  fallbacks.
- **Kinetic (Cundall) damping** — track `KE = ½ vᵀM*v`, zero all velocities at each
  KE peak. Parameter-free, robust through limit points/softening.
- A **dual static-rest termination** (true static residual `‖f_ext − f_int‖`
  *and* KE decay) surfaced via a `getResponse`/integrator query for a
  script-owned relax loop.

**Out of scope (recorded as follow-ups in §8):**
- **Viscous-critical damping** (`C* = c·M*` auto-tuned to ~critical for the lowest
  fictitious mode). Faster than kinetic when the frequency estimate is good, but
  fragile when it is off and adds a knob to mistune. Deferred to v2; kinetic is the
  robust default and ships first.
- **Per-region / per-DOF KE zeroing** (global KE can mask one region ringing while
  another is quiet → over-aggressive global zeroing). Reuses `EnergyBalanceRecorder`
  `MeshRegion` machinery; a follow-up, not v1.
- **In-engine self-termination** (a new `analysis` type that stops at rest). The
  stock transient driver has no convergence concept; an in-engine stop needs a new
  analysis type — **rejected** on the leaf-only invariant (§7), same call ADR-20
  made for its in-engine retry driver.
- **True snap-back path following** — DR traces a *dynamics-regularized* path; the
  rest states approximate the quasi-static branch but DR is **not** a path-follower.
  Same honest caveat as ADR-20 §3.4. The dissipation-controlled arc-length follower
  (ADR-20 §8) remains the tool for that.

### 1.2 Why a `TransientIntegrator`, not a `StaticIntegrator` — settled by source

This is the load-bearing architecture decision and it is **decided by verified
evidence, not preference.** A `StaticIntegrator` DR was proposed (the
"static-internal-loop / parasite-on-Linear" stance) and is **rejected**:

| Claim under review | Verdict | Decisive evidence |
|---|---|---|
| A `StaticIntegrator` can own the outer solve loop / veto the SOE solve | **refuted** | `StaticAnalysis::analyze` (StaticAnalysis.cpp:166-214) hardcodes `newStep()`→`theAlgorithm->solveCurrentStep()`→`commit()`; the algorithm (`Linear.cpp:131`) calls `theSOE->solve()` **unconditionally**. The integrator never owns the loop. (FC1, confirmed) |
| A `StaticIntegrator` can get its mass into the LHS via the nodal path | **refuted** | `StaticIntegrator::formNodTangent` is a hard-error stub returning −1 (StaticIntegrator.cpp:96-102); only a `TransientIntegrator`-style `formTangent` runs the DOF_Group loop that calls `addMtoTang(getMass())`. (FC4 / F2, confirmed) |
| `StaticIntegrator` DR can run a velocity-zeroing kinetic damper | **refuted** | `StaticIntegrator` never references `setVel`/`setResponse`; `StaticAnalysis` carries no velocity concept, so a velocity write is read by nobody. The `setResponse(U,V,A)` plumbing kinetic damping needs is `TransientIntegrator`-only (CDL.cpp:443,497). (F6, confirmed) |
| The "parasite on Linear" smuggle (identity SOE + march hidden in `update()`) | **buildable but rejected** | Verified *possible* (FC2): override `formEleTangent`→identity, hide the march in `update()`. But it forces an un-enforceable "use `algorithm Linear` only" constraint, the outer `ConvergenceTest` measures the **identity** residual (~0) not the true unbalance, and it reuses far less. The reviewers' own honest verdict: "the worse of the two viable architectures … recommend AGAINST." |

By contrast, the `TransientIntegrator` path is a **live, shipped existence proof**:
`CentralDifferenceLadruno` is a `TransientIntegrator` (CDL.h:117) that runs the
exact matrix-free M-only scheme through the **stock** `DirectIntegrationAnalysis`
driver with zero core edits — `formEleTangent = zeroTangent(); addMtoTang()` puts
**mass only** on the LHS, `system Diagonal` makes the per-step solve a literal
`X = B/aᵢᵢ` (DiagonalDirectSolver.cpp:124-126), and `algorithm Linear` guarantees
exactly one solve/step (F1, confirmed). `ExplicitDifferenceStatic` (33001) is a
*second* shipped existence proof of **quasi-static DR specifically** done as a
`TransientIntegrator` (F2 confirmed: a Rayleigh-damped quasi-static explicit
scheme, registered + dispatched transient — what `LadrunoDynamicRelaxation`
auto-tunes rather than asks the user to set).

**Honest framing of "matrix-free."** Neither integrator *bypasses* the SOE — the
algorithm owns the solve. "Matrix-free" means the global LHS is the **diagonal
fictitious mass** `M*`, so `theSOE->solve()` degenerates to an element-wise
`M*⁻¹·R` apply; **K is never assembled or factorized in the stepping loop.** This
is the correct restatement, verified against every relevant `Linear`/Diagonal/CDL
path. We do not claim SOE bypass.

---

## 2. Decision

1. **Ship one `TransientIntegrator` leaf, `LadrunoDynamicRelaxation` (classTag
   33005)**, that performs quasi-static dynamic relaxation under the stock
   `DirectIntegrationAnalysis` driver with `algorithm Linear` + `system Diagonal`.
   Mirror `CentralDifferenceLadruno`'s leap-frog skeleton verbatim and
   `ExplicitDifferenceStatic`'s quasi-static intent. **No `StaticIntegrator`, no
   analysis-core surgery, no cloned driver** (§1.2, §7).

2. **`M*` = the integrator's own fictitious lumped diagonal, built once in
   `domainChanged()`, default Gershgorin row-sum.**
   `m_i = (Δt²/4)·Σ_j |K_ij|` per DOF, harvested by a **one-time element-stiffness
   probe** that reuses the `CriticalTimeStep` per-element scan
   (`computeCriticalTimeStep`, the same path CDL calls in `domainChanged`,
   CDL.cpp:292 — SOE-free, no global factorization, F3 confirmed). This gives every
   DOF the same fictitious frequency `≈ 2/Δt` ⇒ scale-free critically-fast
   convergence and automatic stability (`ω·Δt ≈ 2` by construction). **`M*` is the
   integrator's cached `Vector`, applied by the integrator — it does NOT trust
   `Element::getMass()`** (consistent-by-default, density-scaled, zero on
   zero-density research models — the exact ADR-20 §2.5 BLOCKER, FC4/F4 confirmed
   the `getMass` hooks cannot carry an artificial mass). Modes:
   `-mass {gershgorin|lumped|unity}`, default `gershgorin`.

3. **Damping = kinetic (Cundall), the v1 default and only v1 mode. NO viscous
   constant.** This is the whole point of the stance: nothing to mistune,
   maximally robust through limit points/softening, and the convergence test is the
   **true** static residual `‖f_ext − f_int‖` with **no `f_v` pollution** (contrast
   ADR-20's viscous road, whose `ConvergenceTest` sees the polluted residual). The
   KE-peak detector lives **inside** the integrator (its own `KE = ½ vᵀM*v` from
   `M*` and the full-step velocity), because a `Recorder` is a passive post-commit
   observer with no integrator handle and cannot drive control flow (F5 confirmed:
   `EnergyBalanceRecorder::record(int,double)` has no `Integrator`/SOE/velocity
   handle). `EnergyBalanceRecorder` is the **external watchdog** that logs the KE
   sawtooth and `RES→0`, not the driver.

4. **Termination = a dual test, surfaced to a script-owned relax loop.** The
   integrator computes each commit (i) the true static residual norm
   `‖f_ext − f_int‖/‖f_ext‖_ref` and (ii) the KE level; convergence = both below
   tol. Because the stock transient driver loops a **fixed N** with **no
   convergence early-exit** (F5/F6 confirmed: `DirectIntegrationAnalysis::analyze`
   only bails on a **negative** result), rest-detection is owned by a thin
   `Ladruno_scripts` Python/Tcl helper that calls `analyze(chunk, dt)` and polls a
   `getResponse`/integrator query. This is the **same script-owned pattern**
   ADR-20 landed on for Layer-B retry — accepted, zero-touch.

5. **Reuse the `CriticalTimeStep`/`CTSLumping` probe for `M*` AND for the dt
   stability advisory.** The probe is SOE-free and already wired into the CDL
   `domainChanged` path. **Honest scope cut (F5/FC5):** `CriticalTimeStep`'s loop
   does not *return* an assembled global diagonal mass today — it deletes `mdiag`
   each iteration and returns a scalar dt. So `M*` harvesting is **reuse-of-pattern,
   not drop-in**: v1 adds a small global-diagonal accumulation pass alongside the
   existing per-element scan (integrator-local code in `domainChanged`, **not** a
   `CriticalTimeStep` API change — keeps it a leaf). Optionally refresh `M*` every
   `-recompute N` steps for strong softening (mirrors CDL's `-recompute`), **re-zero
   velocities on refresh** or it injects energy (R6).

6. **No new analysis type, no new SOE, no constraint-handler change.** Works with
   the existing `Linear` algorithm, the `DiagonalSOE`, and the stock
   `DirectIntegrationAnalysis` loop exactly as `CentralDifferenceLadruno` does
   today. The single-solve `updateCount>1` guard is reused verbatim.

7. **Police artificial dissipation / certify rest with the existing
   `EnergyBalanceRecorder`.** It already computes `KE = ½ vᵀMv` from `node->getVel()`
   (EnergyBalanceKernel.h:96,133; F4 confirmed) — the same vector DR pushes to the
   nodes via `setResponse(...,Vfull,...)` before commit. A converging relaxation
   shows a decaying KE sawtooth and `RES→0`. **Mass-weighting caveat (F4):** the
   recorder uses `½ vᵀMv` while CDL's internal proxy is an unweighted `½ v·v`; DR's
   internal detector must use the **mass-weighted** form to agree with the recorder
   and to locate the true KE peak on non-uniform `M*`.

---

## 3. The model

### 3.1 The DR march (leap-frog, kinetic-damped)

Per pseudo-step `n` (verbatim CDL leap-frog skeleton, no damping force term):

$$
\mathbf a_n = (\mathbf M^{*})^{-1}\big(\mathbf f_\text{ext} - \mathbf f_\text{int}(\mathbf u_n)\big),
\quad
\mathbf v_{n+1/2} = \mathbf v_{n-1/2} + \Delta t\,\mathbf a_n,
\quad
\mathbf u_{n+1} = \mathbf u_n + \Delta t\,\mathbf v_{n+1/2}.
$$

`a_n` is assembled matrix-free: `formUnbalance`→`formElementResidual` gives
`−f_int` via `addRtoResidual` (FE_Element.cpp:486-505: `R += −getResistingForce()`,
no tangent touched, F6 confirmed); external load enters via `addPtoUnbalance`; **no
`C` term anywhere** (kinetic damping is velocity-zeroing, not viscous). The
`(M*)⁻¹` apply is realized through the **diagonal SOE solve** (`addMtoTang(M*)`
LHS + `system Diagonal`), not a hand-rolled divide — the verified CDL mechanism.

### 3.2 Kinetic (Cundall) damping

Track `KE_n = ½ vₙᵀ M* vₙ` using the **centered full-step** velocity
`vₙ = ½(v_{n−1/2} + v_{n+1/2})` — the same `Vfull` CDL outputs via
`setResponse(Ut, Vfull, U)` before commit (F7 confirmed: `Vfull = Vhalf + ½Δt·a`,
pushed to nodes as the last `setResponse` before `commit`, so the recorder's KE and
the detector's KE derive from the **identical** vector). Keep `KE_{n-1}, KE_n`.

A **KE peak** is detected when `KE_n < KE_{n-1}` (energy started falling ⇒ the
march just passed an equilibrium-energy maximum). On a peak:

1. (optional `-interp`) **parabolic 3-point KE fit** → back-step `u` to the true
   KE-max configuration before zeroing (reduces overshoot ringing — Cundall
   standard);
2. **zero all velocities** (`Vhalf→Zero()`, `Vfull→Zero()`, keep `Aprev`) and
   continue from rest.

This converges to static equilibrium in `O(√(cond))` peaks — parameter-free,
provably robust through limit points and softening.

### 3.3 Termination

Static rest = **both** hold over a short window:
- **true static residual** `‖f_ext − f_int(u)‖_∞ / ‖f_ext‖_ref < tolR` — the
  physically meaningful test, immune to the `f_v`-pollution that the ADR-20 viscous
  road suffers (here there is **no** `f_v`);
- **KE decay** — successive KE-peak magnitudes ratio `< tolKE` (the dynamics have
  rung down).

The integrator exposes both via `getResponse('residualNorm')` /
`getResponse('kineticEnergy')` and a boolean `converged` query. The relax helper
(§5) polls these; the integrator does **not** try to stop the driver in-engine.

### 3.4 Honest limits (carry into §6 tests)

- **Pseudo-`Δt` is fictitious** and folds into the `M*` scale; with `-mass
  gershgorin` the user sets `dt=1` and stability is guaranteed by construction
  (`ω·Δt ≈ 2`). With `-mass unity` the real `dt_cr` matters again → keep CDL's
  `dt_cr` advisory + `-cflAbort` available.
- **Gershgorin on softening/indefinite `K`:** take `|K_ij|` (absolute) so row sums
  stay positive, and **floor** `m_i`; optionally `-recompute` to track softening.
- **Global KE zeroing can be over-aggressive** on heterogeneous models (one region
  rings while another is quiet). Accepted for v1 (Cundall is global by design),
  documented; per-region KE is a §8 follow-up.
- **DR is not a path-follower.** Rest states approximate the quasi-static branch;
  true snap-back is *regularized but not traced*. Same caveat as ADR-20 §3.4.

---

## 4. Architecture

### 4.1 Class skeleton

```
class LadrunoDynamicRelaxation : public TransientIntegrator {  // classTag 33005
  // --- leap-frog state (verbatim from CentralDifferenceLadruno) ---
  Vector *Ut, *Vhalf, *Aprev, *Vfull, *Azero;
  int    updateCount;                       // single-solve guard (>1 ⇒ error)
  double dtPseudo;                          // fictitious step (folds into M*)
  // --- fictitious mass (the integrator's OWN, NOT getMass) ---
  Vector *Mstar;                            // cached lumped diagonal, built in domainChanged
  int    massMode;                          // gershgorin | lumped | unity
  int    recomputeN;                        // 0 = build once
  // --- kinetic damping + termination ---
  double KEprev, KEnow, KEpeakPrev;         // KE = 1/2 v^T M* v (mass-weighted!)
  bool   interpPeak;                        // -interp parabolic back-step
  double tolR, tolKE; int maxPeaks, peakCount;
  double resNorm;                           // ||fext - fint|| cached at commit
  // overrides
  int newStep(double dt);                   // leap-frog advance + setResponse(Ut,Vhalf,Azero)
  int update(const Vector &a);              // accept solved accel; KE peak detect + zero
  int commit();                             // setResponse(Ut,Vfull,U); cache resNorm/KE
  int formEleTangent(FE_Element*);          // zeroTangent(); add M* (NOT addMtoTang/getMass)
  int formNodTangent(DOF_Group*);           // add nodal slice of M*
  int domainChanged();                      // allocate state; build M* via CTS-pattern probe
  int getResponse(...) / setResponse(...);  // residualNorm, kineticEnergy, converged
  double getCriticalTimeStep();             // dt_cr advisory (reuse CTS)
  Vector &getVel();                         // returns Vhalf (modal-damping hook)
  int sendSelf(...); int recvSelf(...); int Print(...);
};
```

### 4.2 The fictitious-mass seam (`domainChanged` / `formEleTangent`)

`M*` is **not** `Element::getMass()`. In `domainChanged()` run the SOE-free
per-element probe (reusing the `CriticalTimeStep` scan, F3) and accumulate a global
diagonal:

```
// for each element e: K_e = e->getTangentStiff()  (or getInitialStiff())
//   for each local row i (global DOF g): Mstar[g] += (dt*dt/4) * sum_j |K_e(i,j)|
// floor each Mstar[g] >= eps;  (lumped: use e->getMass() diag * scale;  unity: 1)
```

Then put `M*` on the LHS so the diagonal SOE solve is `M*⁻¹·R`. Because the verified
`addMtoTang`/`addM_Force` hooks **always pull `getMass()`** (FC4/F4: no overload
accepts a substitute matrix), DR **cannot** route `M*` through them. Instead the
integrator pokes `M*` into the assembled diagonal itself (the
`IncrementalIntegrator` base already carries the `isDiagonal/diagMass/mV`
cached-diagonal pattern — `doMv`, IncrementalIntegrator.cpp:652-678 — exactly
"own M* as a cached diagonal and apply it yourself"). This sidesteps the `getMass`
trap entirely; it is the cleanest matrix-free path and the reason DR is robust on
zero-density research models where a `getMass`-based `M*` would be 0.

> **Build-plan detail:** `formEleTangent` may still `zeroTangent()` and contribute
> a placeholder so the diagonal SOE has a well-formed LHS; the authoritative `M*`
> diagonal is the integrator's cached `Vector`. Prototype this seam first — it is
> the one genuinely new mechanic over the CDL skeleton.

### 4.3 Reuse map (CDL donor ≈ 70%)

REUSE verbatim/near-verbatim from `CentralDifferenceLadruno`:
- state vectors `Ut/Vhalf/Aprev/Vfull/Azero` + `domainChanged` allocate-and-seed
  from committed DOF (CDL.cpp:233-282);
- the leap-frog advance + `setResponse(Ut,Vhalf,Azero)`→`updateDomain`→solve→
  `update(a)` cycle (CDL.cpp:438-522);
- the first-step starter `a0 = M*⁻¹(P0 − f_int(u0))`, `v_{−1/2} = v0 − ½Δt·a0`
  (CDL.cpp:410-434) — **minus the `C·v0` term** (no damping);
- the `updateCount>1` single-solve guard (CDL.cpp:455-463);
- `getVel()` returning `Vhalf` (CDL.cpp:539-542);
- NaN/Inf circuit-breaker + `divergenceFactor` KE-runaway guard (CDL.cpp:482-512);
- `sendSelf/recvSelf` Vector-pack pattern;
- the `CriticalTimeStep`/`CTSLumping`/`CTSResult` `dt_cr` probe (CDL.cpp:292),
  **extended** to also emit the Gershgorin `M*` and the `dt_cr` advisory.

REUSE from `ExplicitDifferenceStatic`: the proof that a `TransientIntegrator` +
`algorithm Linear` quasi-static DR works end-to-end with zero core edits, and its
`formNodalUnbalance`-override precedent for legally modifying residual assembly.

NEW (DR-specific): the Gershgorin `M*` builder; the mass-weighted KE-peak detector
+ velocity-zeroing + optional parabolic back-step; the integrator-owned `M*`
diagonal poke (bypassing `getMass`); the dual residual+KE termination query.

### 4.4 Build-control obligations (REQUIRED, same PR)

- `SRC/classTags.h`:
  `#define INTEGRATOR_TAGS_LadrunoDynamicRelaxation 33005 // N. Mora-Bowen (Ladruno) — quasi-static dynamic relaxation, kinetic damping; integrator band >=33000 (siblings 33000 ExplicitBathe … 33004 LadrunoArcLength). Tag bands are PER-REGISTRY; 33005 in the integrator registry does not collide with any ELE_TAG.`
  (verified next free integrator tag — F5/F7: 33005 absent from the whole SRC tree.)
- **Broker:** `FEM_ObjectBrokerAllClasses::getNewTransientIntegrator` (NOT
  `getNewStaticIntegrator`, NOT `getNewIncrementalIntegrator`, NOT the
  `FEM_ObjectBroker.cpp` stub which returns 0) — add
  `case INTEGRATOR_TAGS_LadrunoDynamicRelaxation: return new LadrunoDynamicRelaxation();`
  (F7 confirmed: CDL lives in `getNewTransientIntegrator`; placement is by base
  class, so a default/`recvSelf`-able ctor is required).
- **Tcl + Python dispatch:** `SRC/interpreter/OpenSeesCommands.{h,cpp}`
  (fwd-decl `void* OPS_LadrunoDynamicRelaxation();` + string dispatch into the
  **transient** branch, next to the CDL case, `setIntegrator` on the transient
  analysis); `SRC/runtime/runtime/TclPackageClassBroker.cpp` (+include +case);
  `SRC/tcl/commands.cpp` (extern + integrator branch). Each vanilla touch → a
  `LEDGER_vanilla_files.md` row, marked `// Ladruno`.
- `SRC/analysis/integrator/CMakeLists.txt` + `Makefile`: add the new `.cpp/.h/.o`.
- `LEDGER_implementations.md`: new row (Integrator / 33005 / files / status / PR).
- Banner: add a line to `Ladruno_scripts/banner_features.txt` →
  `python Ladruno_scripts/patch_banner.py` (regenerates `tclMain.cpp` +
  `PythonModule.cpp`).
- `Ladruno_scripts/stamp_headers.py`: add the new files to GLOBS + rerun (LADRUNO
  header stamp is non-optional for fork-authored files).

This is **exactly** CDL's 8-site seam (F8 confirmed present for
`CentralDifferenceLadruno`), plus the ledger/banner/stamp obligations.

---

## 5. Public API (proposed)

```tcl
# v1 default: Gershgorin M*, kinetic (Cundall) damping, dual residual+KE termination
integrator LadrunoDynamicRelaxation

# full surface (all optional):
integrator LadrunoDynamicRelaxation \
    -mass     gershgorin|lumped $scale|unity \
    -tolR     1e-5  -tolKE 1e-4 \
    -maxPeaks $N    -recompute $N \
    -dt       1.0   -interp  -cflAbort  -verbose

# required driver recipe (mirrors CDL):
system Diagonal ; numberer Plain ; constraints Transformation
test NormUnbalance 1e-6 1 ; algorithm Linear
analysis Transient
```

```python
ops.integrator('LadrunoDynamicRelaxation', '-mass', 'gershgorin',
               '-tolR', 1e-5, '-tolKE', 1e-4)
```

**The real UX — `Ladruno_scripts/relax_to_static.py` helper (script owns
termination, decision 4):**

```python
def relax_to_static(chunk=50, dt=1.0, maxIter=200000, tolR=1e-5, tolKE=1e-4):
    ops.analysis('Transient')
    done = 0
    while done < maxIter:
        ops.analyze(chunk, dt)                       # stock transient driver
        R  = ops.getResponse('residualNorm')         # DR query: ||fext - fint||
        KE = ops.getResponse('kineticEnergy')        # DR query: 1/2 v^T M* v
        if R < tolR and KE < tolKE:
            return True
        done += chunk
    return False
```

So: `analysis Transient` + a relax loop, **not** a new analysis type. The helper is
the documented entry point; the integrator + `getResponse` query are the in-engine
half.

---

## 6. Testing / oracle matrix (Zone-A)

> **Prototype status:** the snap-through fixture is **already standing and green**
> in `tests/_proto_arch_snapthrough.py` (shallow von Mises `corotTruss`, `E·A`
> tuned so the limit load ≈ 3.80). `run_dynamic_relaxation(P_over_crit=1.15×crit)`
> already snaps through to the far stable branch (`uy ≈ −0.217`, currently via
> stock `Newmark` + Rayleigh) — this is the DR oracle. The DR-1/DR-6 legs become
> pytest fixtures by **swapping the `Newmark` runner for `LadrunoDynamicRelaxation`**
> when code lands, asserting the **same far-branch `uy`**.

| ID | Check | Oracle / pass |
|----|-------|---------------|
| DR-1 | `LadrunoDynamicRelaxation` (defaults) relaxes the von Mises truss under a **sub-critical** load to static rest | converges to the stock-`LoadControl` equilibrium `uy` (same point both reach), `‖f_ext − f_int‖ < tolR`, `KE → 0` |
| DR-2 | Gershgorin `M*` builder unit test | `m_i = (Δt²/4)Σ_j|K_ij|`, floored, on a 2-element patch — closed form |
| DR-3 | KE-peak detector + velocity-zeroing | synthetic KE series → peak located, `Vhalf/Vfull` zeroed at peak; mass-weighted KE matches `EnergyBalanceRecorder` KE (F4) bit-for-bit |
| DR-4 | One pseudo-step is exactly one diagonal solve | `updateCount>1` guard fires; SOE LHS is `M*` (diagonal), K never factorized (assert via a probe) |
| DR-5 | True static residual termination (no `f_v` pollution) | at "converged" `‖f_ext − f_int‖` is genuinely `< tolR` (contrast ADR-20 viscous road) |
| **DR-6** | **Snap-through** — `LoadControl` **fails** at λ=3.80; DR under 1.15×crit constant load **snaps through** to the far branch | far-branch `uy ≈ −0.217` (matches the AL-6 prototype already in `_proto_arch_snapthrough.py`) |
| DR-7 | DR rest point == `DisplacementControl` reference branch point | within tol of the stable-branch equilibrium the prototype's reference path traces |
| DR-8 | **`M*` independence (BLOCKER, FC4/F4)** — a **zero-density / nodal-mass-only** model still relaxes | `-mass gershgorin` engages (artificial `M*` ≠ `getMass`); a `getMass`-based `M*` would give `M*=0` and never move |
| DR-9 | `EnergyBalanceRecorder` watchdog | KE column shows the decaying sawtooth (rise-zero-rise), `RES → 0` over the run |
| DR-10 | `-recompute N` refresh re-zeros velocities | KE does **not** jump at a refresh step (R6 — no injected energy) |
| DR-11 | `sendSelf`/`recvSelf` round-trip (all flags) | state restored byte-faithful |
| DR-12 | Gershgorin convergence rate vs `-mass unity` | gershgorin reaches rest in **fewer total pseudo-steps** (scale-free claim) |

**Canonical model:** the shallow von Mises truss
(`tests/_proto_arch_snapthrough.py`) for DR-1…DR-9/DR-12 — **already built**; a
single zero-density / nodal-mass-only variant for DR-8. Both tiny, deterministic,
Zone-A-portable.

---

## 7. Rejected alternatives

### 7.1 `StaticIntegrator`-internal-loop ("parasite on Linear") — REJECTED

Smuggle a full DR march inside `update()`, feed the `Linear` algorithm an identity
`formEleTangent` so `theSOE->solve()` is inert, return negative on stall so
`StaticAnalysis` reverts. **Verified buildable (FC2) but rejected**:
- The `StaticIntegrator` interface gives DR **nothing it needs**: no mass access
  (`formNodTangent` is a hard-error stub, FC4), no solve ownership (the algorithm
  owns it, FC1), no velocity plumbing for kinetic damping (`setVel`/`setResponse`
  are `TransientIntegrator`-only, F6).
- It forces an **un-enforceable** "use `algorithm Linear` + `system Diagonal` only"
  user constraint, and the outer `ConvergenceTest` measures the **identity**
  residual (~0) not the true unbalance — mis-wiring silently "converges" at
  non-equilibrium (the one genuinely fragile trick).
- The reviewers' own honest verdict: "the worse of the two viable architectures …
  recommend AGAINST unless a strict no-dt static command API is a non-negotiable
  product requirement." It is not.

### 7.2 Viscous-critical damping in v1 — DEFERRED (not rejected)

`C* = c·M*` auto-tuned to ~critical for the lowest fictitious mode converges in
fewer pseudo-steps when `ω₁` is well-estimated, but is **fragile** when the
estimate is off and adds a knob to mistune. Kinetic damping is parameter-free and
limit-point-proof, so it ships first; viscous-critical is a v2 opt-in
(`-damping viscous`) reusing CDL's `betaK`-trap guard (forbid stiffness-proportional
damping, which would re-introduce K-assembly and break matrix-free-ness).

### 7.3 In-engine self-termination / cloned driver — REJECTED

A new `analysis LadrunoRelax` type that stops at rest in-engine. The stock
transient driver has **no convergence concept** (F5/F6: fixed-N loop, bail only on
negative result), so an in-engine stop needs a new analysis type + parser dispatch
+ ~200 lines of cloned `analyze()`. **Rejected** on the leaf-only invariant — the
same call ADR-20 §7 made. Script-owned termination via the relax helper costs one
small Python file and **zero** core touch.

### 7.4 `EnergyBalanceRecorder`-driven KE-peak detection — REJECTED (impossible)

A `Recorder` is a passive post-commit observer: `record(int commitTag, double
timeStamp)` has **no** `Integrator`/`AnalysisModel`/SOE handle and its return value
does not steer the driver (F5 confirmed). It **cannot** zero velocities at a KE
peak. The detector therefore lives inside the integrator; the recorder is the
external watchdog only.

---

## 8. Follow-ups (deferred)

1. **Viscous-critical damping mode** (`-damping viscous`, §7.2) — faster on
   well-conditioned models; v2.
2. **Per-region / per-DOF KE zeroing** — reuse `EnergyBalanceRecorder` `MeshRegion`
   machinery so a quiet region is not zeroed by a loud one's peak (R/3.4).
3. **Native finite-strain / large-rotation DR** consuming `LadrunoBrick -geom
   finite` softening payloads (the intended consumers — ASDConcrete / Lemaitre).
4. **True-equilibrium residual already native here** — note for ADR-20: DR's
   convergence test *is* the true `‖f_ext − f_int‖` (no `f_v`), so DR is the clean
   oracle for the ADR-20 §8.4 follow-up.
5. **Lift `M*`-build into a shared `CriticalTimeStep` global-diagonal accessor**
   (one vanilla row) if a second consumer wants the assembled lumped mass — only if
   demand materializes (today it is integrator-local, F5).

---

## 9. Relationship to other ladruno work

- **`LadrunoArcLength`** ([[20_ladruno_arclength_stabilized_adr]]): the **implicit**
  answer to the same limit-point payload — regularize `K_T` and let Newton pass.
  This ADR is the **explicit / matrix-free** counterpart — never form `K_T` at all.
  Its AL-6 transient dynamic-relaxation oracle **is** this integrator; the two ADRs
  share `tests/_proto_arch_snapthrough.py`.
- **`CentralDifferenceLadruno`** ([[05_robust_central_difference]]) and
  **`ExplicitDifferenceStatic`**: the donor leaf skeletons — ~70% reuse — and the
  two shipped existence proofs that a `TransientIntegrator` matrix-free quasi-static
  scheme runs on the stock driver with zero core surgery.
- **`EnergyBalanceRecorder`** ([[04_explicit_dynamics_and_energy_balance]]): the
  mandatory KE/dissipation watchdog (decision 7) — the same role it plays for
  ADR-20's stabilized mode.
- **`LadrunoBrick` / Béziers**: consumers — softening solid/concrete
  ([[11_brick_asdconcrete_integration]]) and Lemaitre-damaged
  ([[15_lemaitre_ductile_damage_adr]]) payloads whose limit points DR clears by
  relaxation rather than by factorizing an indefinite tangent.
