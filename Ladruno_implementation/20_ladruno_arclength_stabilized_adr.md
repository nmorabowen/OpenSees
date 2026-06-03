---
title: ADR — Adaptive + stabilized static path-following (LadrunoArcLength)
project: Ladruno
status: proposed — scoped, no code
priority: medium
owner: nmora
tags:
  - implementation
  - integrator
  - static
  - arc-length
  - softening
  - stabilization
  - adr
---

# ADR — Adaptive + stabilized static path-following (`LadrunoArcLength`)

**Status:** proposed — design scoped, **no code yet** (sibling of the contact /
E-FEM scoping ADRs) · **Registry:** `StaticIntegrator` · **classTag:**
`INTEGRATOR_TAGS_LadrunoArcLength = 33004` (next free in the ladruno integrator
band; siblings ExplicitBathe=33000 … CentralDifferenceLadruno=33003) ·
**Supersedes nothing** — adds a fork integrator alongside stock `ArcLength`
(33004 leaves vanilla `INTEGRATOR_TAGS_ArcLength` untouched) · **Oracle:** stock
`ArcLength` (bit-identical when all extra flags are off) + a transient
dynamic-relaxation run (for the stabilized mode).

> [!note] As-designed scope (v1 target) — revised after the agent-wise review
> One self-contained `StaticIntegrator` leaf class that is a **strict superset of
> stock `ArcLength`**, with **two mutually-exclusive robustness modes** plus a
> script-driven retry hook. Default (no flags) ⇒ **bit-identical to `ArcLength`**:
> 1. **Arc-length mode (`-adapt Jd ℓmin ℓmax`)** — stock cylindrical constraint +
>    Ramm desired-iteration radius adaptation ("Layer A"): grow/shrink the arc
>    radius between steps from the converged corrector-iteration count of the last
>    step. This is the rule stock `LoadControl` / `DisplacementControl` /
>    `MinUnbalDispNorm` already carry but `ArcLength` never received.
> 2. **Stabilized mode (`-stabilize f [-adaptStab]`)** — Abaqus-`STABILIZE`-style
>    **artificial viscous stabilization** ("the viscous road"): inject a dashpot
>    force `c·M*·Δu/Δt` into residual + tangent so the tangent stays
>    positive-definite through a **limit point**. `c` auto-calibrated to a
>    dissipated-energy fraction `f`; `-adaptStab` lets `c` track the run.
>    **DECISION (review §2.8):** `-stabilize` and the arc-length quadratic are
>    **mutually exclusive** — exactly as Abaqus attaches `STABILIZE` to plain
>    `*STATIC` and **never** to Riks. With `-stabilize` on, the class runs
>    **viscous-regularized incremental load control** (predictor `Δλ` from the
>    adaptive rule, **no** quadratic corrector constraint), so the stabilized
>    `K_T`/residual never pollute the arc-length geometry. Combining them is a
>    non-goal.
> 3. **`reduceStep()` / `increaseStep()` mutators + a `revertToLastStep()`
>    override** — let a *script* drive cut-and-retry across a failed step
>    **without reconstructing the integrator**. The override restores
>    `deltaLambdaStep` / `deltaUstep` / `currentLambda` / `signLastDeltaLambdaStep`
>    from a **committed snapshot** (review §2.10: stock `ArcLength` mutates these
>    *before* its own `b24ac<0` bail, so "state is preserved for free" is false —
>    we must snapshot at `commit()` and restore on revert). "Layer B" delivered
>    **without cloning the analysis driver** — see §7.
>
> **`M*` = unity-density artificial mass (review BLOCKER §2.5).** A
> `StaticIntegrator` assembles **element** mass only — `formNodTangent()` is a
> hard-error stub, so **nodal `mass` and zero-density elements contribute NOTHING**;
> and `Element::getMass()` is **consistent-by-default and density-scaled**. So the
> real mass matrix is unusable as `M*`. v1 builds a **unity-density artificial
> lumped `M*`** (geometry-based, density-independent — Abaqus-faithful) in
> `domainChanged()`. This was demoted from a "documented upgrade" to the v1 default.
>
> **Hard constraint honored:** the analysis core (`StaticAnalysis::analyze`) is
> **not touched, not cloned, not subclassed.** The only vanilla edit is the
> normal integrator-registration row in the broker + Tcl/Python dispatch — the
> same seam every ladruno integrator already uses (ledgered in
> `LEDGER_vanilla_files.md`).
>
> > **Reviewer-confirmed correctness caveats (folded into §3.4):** stabilized
> > results are **step-size / Δt-dependent** (not a step-converged equilibrium
> > path); the convergence test sees the **stabilized** residual `λp − f_int − f_v`,
> > so `‖true unbalance‖` is nonzero by `‖f_v‖` at "convergence"; the
> > positive-definiteness claim holds for **snap-through / limit points**, **not**
> > true snap-back (which still needs a path-following constraint). DDM sensitivity
> > is **out of scope** in stabilized mode (guarded), preserving `ArcLength`'s
> > `formEleResidual` sensitivity override only on the unstabilized path.

---

## 1. Context

### 1.1 The phenomenon, and what is / isn't in scope

A static push past a **limit point** (snap-through) or **snap-back** drives the
tangent stiffness `K_T` singular and then indefinite. Standard Newton + load
control diverges there; the classic remedy is **path-following** — augment the
`n` equilibrium equations with one constraint equation that lets the load factor
`λ` become an unknown, so the solver tracks arc length along the equilibrium
path rather than along the load axis.

OpenSees ships the constraint-equation family — `ArcLength`, `ArcLength1`,
`HSConstraint`, `MinUnbalDispNorm` — but with **two gaps** this ADR closes:

- **No radius adaptation on the true arc-length constraint.** Stock `ArcLength`
  freezes `arcLength2` at construction
  ([`ArcLength::newStep`](../SRC/analysis/integrator/ArcLength.cpp), the
  `dLambda = sqrt(arcLength2/(dUhat·dUhat + alpha2))` line). The desired-
  iteration cut/grow rule lives only in `LoadControl` /
  `DisplacementControl` / `MinUnbalDispNorm`, never grafted onto `ArcLength`.
- **No in-engine robustness for *local* instability.** A *global* arc-length
  constraint cannot parameterize localized softening (one band snaps back while
  the bulk loads); radius cutting does not fix the constraint geometry. The
  industry answer for that case is **not** more arc-length — it is **artificial
  viscous stabilization** (Abaqus `*STATIC, STABILIZE`), which regularizes `K_T`
  back to positive-definite and lets **ordinary Newton** pass the instability.

**In scope:**
- Adaptive arc radius (Layer A) on the stock cylindrical arc-length constraint.
- Integrator-level artificial viscous stabilization (the viscous road), with an
  auto-calibrated, optionally-adaptive damping factor.
- Script-drivable step-size mutators so a failed step can be re-attempted from
  the last converged state without integrator reconstruction (Layer B, no
  analysis-core surgery).

**Out of scope (recorded as follow-ups in §8):**
- **Dissipation-/energy-release-controlled path following** (Gutiérrez 2004;
  Verhoosel–Remmers–Gutiérrez 2009) — constrain incremental *dissipated energy*
  instead of a displacement norm. This is the genuinely state-of-the-art tool for
  localized/fracture softening and would be a *second* integrator; it slots into
  the same `newStep`/`update` seam. Deferred — bigger, and orthogonal.
- **Indirect / CMOD control** (control a monotone relative DOF) — needs a
  weighted/relative constraint; deferred.
- **Material-side viscous regularization** (Duvaut–Lions) — fixes the
  *ill-posedness* of the softening band (mesh objectivity), a different goal from
  numerical convergence; lives in the material, not here. Complementary, not
  competing.
- **A cloned / forked `StaticAnalysis` driver** — explicitly **rejected**, see
  §7.

### 1.2 Why integrator-level, not element-level, not analysis-level

Three distinct things are all called "stabilization"; only one belongs here:

| "Stabilization" | Lives in | Cures | In this fork |
|---|---|---|---|
| **Hourglass control** (`Kstab`) | element | zero-energy modes from reduced integration | ✅ already shipped on `LadrunoBrick -hourglass`, Bézier char-length |
| **Material viscous regularization** (Duvaut–Lions) | material | ill-posedness / mesh-dependence of softening band | future, material-side (§8) |
| **Artificial viscous stabilization** (`STABILIZE`) | **integrator** | unstable *equilibrium* (singular/indefinite `K_T`) | **this ADR** |

The artificial-viscous form is a **residual + tangent modification**, and the
`StaticIntegrator` is exactly the object that forms residual and tangent
(`formEleResidual` / `formEleTangent`). So it fits in a leaf integrator class
with **no analysis-core change** — the decisive property given the project's
"leaf additions, never core surgery" track record (every prior ladruno
integrator is a self-contained leaf). It is *not* element-side: like Abaqus
`STABILIZE`, it is blind to which element it hits.

---

## 2. Decision

1. **Ship one `StaticIntegrator` leaf, `LadrunoArcLength` (classTag 33004)**, a
   strict superset of stock `ArcLength`. All additions are default-off and gated;
   with no extra flags the class is **bit-identical to `ArcLength`** (the fork's
   standard identity gate). Oracle = stock `ArcLength`.

2. **Layer A (adaptive radius) is the floor.** Port the
   `factor = Jd / J_last` rule (clamped to `[ℓmin, ℓmax]`) from `LoadControl`
   into `newStep()`, acting on the arc radius `√arcLength2`. This brings OpenSees
   arc-length to **Abaqus modified-Riks parity** (Riks already auto-adapts the
   increment from iteration count) — useful, but understood as parity, not a
   differentiator.

3. **The viscous road is the robustness payload.** Add Abaqus-`STABILIZE`-style
   artificial viscous stabilization **inside the integrator** (§3.3, §4.2). It
   regularizes `K_T` so ordinary Newton clears a **snap-through / limit point** in
   the stock loop, **dissolving the need to clone the analysis driver** for that
   case. **Scoped (review §2.9):** this is sound for snap-through and *local* limit
   points; it is **not** a cure for true **snap-back** (the path turns back in both
   load and displacement — load-incrementing Newton still cannot follow it), which
   remains the job of an arc-length / dissipation-controlled follower.

3b. **Stabilized mode and the arc-length quadratic are mutually exclusive**
   (review §2.8). `-stabilize` ⇒ viscous-regularized **load control** (no
   constraint quadratic); `-stabilize` off ⇒ stock arc-length (+`-adapt`). This
   keeps the stabilized `K_T`/residual out of the constraint geometry — otherwise
   `deltaUhat = K_T⁻¹p̂` and `deltaUbar = K_T⁻¹R` would be built from a *distorted*
   tangent and the quadratic would no longer measure true arc length.

4. **Layer B without core surgery.** Provide `reduceStep(f)` / `increaseStep(f)`
   mutators **and an explicit `revertToLastStep()` override** so a **script** loop
   can cut-and-retry a failed step from the last *committed* state, with **no
   integrator reconstruction**. The override restores the per-step state from a
   `commit()`-time snapshot (review §2.10 — stock `ArcLength` overwrites
   `deltaLambdaStep`/`deltaUstep` *before* its own `b24ac<0` bail, so the naïve
   "state preserved for free" claim is false). The cloned `StaticAnalysisAdaptive`
   driver is **rejected** (§7).

5. **`M*` = unity-density artificial lumped mass (review BLOCKER §2.5 + §2.6).**
   A `StaticIntegrator` sees **element** mass only (`formNodTangent` is a
   hard-error stub) and `getMass()` is **consistent-by-default + density-scaled**,
   so the real mass matrix is unusable (zero-density elements / nodal `mass` →
   `M*=0`). Build a geometry-based unity-density lumped `M*` in `domainChanged()`.
   This was the original "documented upgrade" (old D2), now promoted to v1 default.

6. **DDM sensitivity is out of scope in stabilized mode** (review §2.11). Stock
   `ArcLength` overrides `formEleResidual` for sensitivity; the unstabilized path
   reproduces that verbatim, the stabilized path **guards** (`if (stabilize &&
   sensitivityFlag) error`) — deriving `∂f_v/∂h` is a non-goal.

7. **Police the artificial dissipation with the existing
   `EnergyBalanceRecorder`** — track viscous-dissipation / strain-energy and warn
   when stabilization is doing too much (the watchdog already exists; real
   synergy). No new recorder.

8. **No new analysis type, no new SOE, no constraint-handler change.** Works with
   the existing `Newton`/`NewtonLineSearch` algorithms, any `LinearSOE`, and the
   stock `StaticAnalysis` loop exactly as `ArcLength` does today. `M*` flows
   through `addMtoTang`/`addM_Force` (never a direct SOE diagonal poke), so it is
   transparent to the Transformation/Penalty constraint handlers (review §2.20).

---

## 3. The model

### 3.1 Stock arc-length recap (the identity baseline)

Cylindrical constraint enforced per step:

$$\Delta\mathbf{u}^{\mathsf T}\Delta\mathbf{u} + \alpha^2\,\Delta\lambda^2 = \ell^2$$

Predictor (`newStep`): solve `K_T\,\hat{\mathbf u} = \hat{\mathbf p}` (reference
load), then

$$\Delta\lambda^{(1)} = \pm\,\frac{\ell}{\sqrt{\hat{\mathbf u}^{\mathsf T}\hat{\mathbf u} + \alpha^2}},$$

sign from `signLastDeltaLambdaStep`. Corrector (`update`): the quadratic in
`δλ` ([`ArcLength::update`](../SRC/analysis/integrator/ArcLength.cpp), the
`a,b,c,b24ac` block) with the positive-`θ` root selection. The notorious failure
mode `b24ac < 0` ("imaginary roots … initial load increment was too large") is
exactly the limit-point breakdown the additions below attack.

### 3.2 Layer A — Ramm desired-iteration adaptation

At the start of a new step, before computing `Δλ^{(1)}`:

$$\ell \leftarrow \mathrm{clamp}\!\Big(\ell\cdot\big(\tfrac{J_d}{J_\text{last}}\big)^{p},\ \ell_\text{min},\ \ell_\text{max}\Big),\qquad p\in\{1,\tfrac12\}$$

`J_d` = desired iterations/step, `J_last` = corrector iterations the last
converged step took. `J_last` is obtained the same way `LoadControl` does it — a
counter bumped once per `update()` call (reset to 0 in `newStep()`,
[`LoadControl.cpp:132,157`](../SRC/analysis/integrator/LoadControl.cpp)), so **no
`ConvergenceTest` coupling is required**. `p = 1` reproduces the stock-LoadControl
behavior; `p = ½` is the gentler Crisfield variant (offer both, default `p = 1`).

> **Guard (review §2.1):** the division by `J_last` must be protected —
> `J_last ← max(1, J_last)`. The new exposure vs stock `LoadControl` is a
> **predictor-only convergence** (the `ConvergenceTest` passing before any
> `update()`), which leaves `J_last = 0` and would make the *next* `newStep()`
> divide by zero / produce an infinite radius. (Stock `LoadControl` only guards
> `numIncr == 0` in its constructor; that does not cover this case.)

> **Naming note (review §2.7, uncertain-but-accurate):** stock `LoadControl`
> applies `factor` to `deltaLambda` (a load increment); `ArcLength` stores
> `arcLength2` (the *square* radius). The port therefore scales `arcLength2 ←
> arcLength2 · (Jd/Jlast)^{2p}` (or scales `√arcLength2` by `(Jd/Jlast)^p`) — write
> the implementation against `arcLength2`, not a notional `ℓ`.

### 3.3 The viscous road — artificial stabilization

Append a vanishing artificial dashpot to the static equilibrium:

$$\mathbf K_T\,\delta\mathbf u \;=\; \lambda\hat{\mathbf p} - \mathbf f_\text{int} \;-\; \underbrace{c\,\mathbf M^{*}\,\frac{\Delta\mathbf u}{\Delta t}}_{\mathbf f_v},
\qquad
\mathbf K_T \leftarrow \mathbf K_T + \frac{c}{\Delta t}\,\mathbf M^{*}.$$

- `M*` = **unity-density artificial lumped mass** (review BLOCKER §2.5/§2.6 —
  *not* `getMass()`, which is consistent-by-default, density-scaled, and zero for
  zero-density elements; and a `StaticIntegrator` never assembles nodal mass).
  Δu is the global per-step accumulated increment, reusing `ArcLength`'s own
  `deltaUstep` member (review §2.3 — `FE_Element` exposes **no** displacement
  accessor; `addM_Force(deltaUstep, −c/Δt)` slices the element DOFs internally).
- `Δt` is **fictitious** (static analysis has no real time) — a pseudo-velocity
  scale. It is folded into `c`, so the user sees **one** knob, not two.
- `f_v` is **small on well-converged stable increments** (small Δu) and
  **concentrates at the instability** — it does **not** vanish there (review §2.2:
  through the snap the per-step Δu is largest, so the injected dissipation is
  largest precisely where it acts — see the §3.4 limits, with which this must stay
  consistent).

**Calibration of `c` (energy-fraction target).** Following Abaqus: choose `c` so
the viscous dissipation over a reference increment is a small fraction `f` of the
model's strain energy:

$$\sum_\text{steps}\mathbf f_v^{\mathsf T}\,\Delta\mathbf u \;\approx\; f\cdot E_\text{strain},\qquad f_\text{default}=2\times10^{-4}.$$

In practice: from the first increment, estimate `E_strain ≈ ½ Δu^T f_int` and back
out `c` from the dashpot work of that same increment. **Assumption (review §2.4):**
this requires the first increment to be **stable and near-linear from the
reference state** — it is invalid if the model is pre-loaded into a nonlinear
regime, so the first-increment calibration must be taken from a genuinely elastic
start (or `c` supplied directly). `-adaptStab` then rescales `c` each step to hold
the running dissipation ratio near `f`. Default (non-adaptive) holds `c` fixed.

**Why this clears a limit point in the stock loop.** Adding `+ (c/Δt)·M*` makes
`K_T` positive-definite through a **snap-through / limit point**, so **ordinary
Newton converges without the arc-length retry geometry**. This is precisely why
Abaqus reaches for `STABILIZE` instead of Riks for *local* instabilities — and why
this case needs no cloned retry-driver (§7). **It does not follow a true
snap-back** (§3.4, review §2.9).

### 3.4 Honest limits (reviewer-confirmed)

- **Step-size / Δt dependence (review §2.13).** Because `f_v = (c/Δt)M*·Δu` scales
  with the per-step increment and the fictitious `Δt`, the stabilized solution is
  **not a step-converged equilibrium path** — refining the load stepping changes
  the answer. Standard practice (carry into §6): verify a stabilized result by a
  `c`-reduction or step-refinement convergence check; AL-5's peak/energy tolerance
  is defined **relative to `c`/step size**, not absolutely.
- **The convergence test sees the *stabilized* residual (review §2.12).**
  `CTestNormUnbalance`/`CTestEnergyIncr` norm `theSOE->getB() = λp − f_int − f_v`,
  so Newton is satisfied while the **true** static unbalance `‖λp − f_int‖` is
  still nonzero by `‖f_v‖`. v1 documents this and recommends tightening tol
  relative to `f`; exposing a true-equilibrium residual for the test is a §8
  follow-up.
- Stabilization **adds artificial energy** — too-large `f` smears the response
  (rounds the snap, raises the peak). The `EnergyBalanceRecorder` watchdog
  (decision 7) is mandatory practice, not optional.
- **Snap-through ≠ snap-back (review §2.9).** Positive-definiteness clears a limit
  point under load/displacement increments; a true snap-back (path turns back in
  *both* λ and u) still requires a path-following constraint (arc-length or the
  dissipation-controlled follower, §8). `-stabilize` there gives a *regularized*
  but **wrong** monotone path.

---

## 4. Architecture

### 4.1 Class skeleton

```
class LadrunoArcLength : public StaticIntegrator {        // classTag 33004
  // --- stock ArcLength state (verbatim) ---
  double arcLength2, alpha2, deltaLambdaStep, currentLambda;
  int    signLastDeltaLambdaStep;
  Vector *deltaUhat,*deltaUbar,*deltaU,*deltaUstep,*phat;
  double a,b,c_,b24ac;                                    // quadratic coeffs
  // --- Layer A ---
  bool   adapt;  double Jd, ellMin, ellMax, pExp;  int numIncrLastStep;
  // --- viscous stabilization (MUTUALLY EXCLUSIVE with the arc-length quadratic) ---
  bool   stabilize, adaptStab;  double fTarget, cVisc, dtPseudo;
  double Estrain0, dissipVisc;                            // calibration + watchdog
  Vector *Mstar;                 // unity-density artificial lumped mass (§2.5 BLOCKER)
  // --- Layer B committed snapshot (§2.10) ---
  double cDeltaLambdaStep, cCurrentLambda;  int cSign;  Vector *cDeltaUstep;
  int    sensitivityFlag, gradNumber;                    // reproduce ArcLength's path
  // overrides
  int newStep();                 // §3.2 radius adapt + predictor (load-ctrl if stabilize)
  int update(const Vector&);     // §3.1 quadratic if !stabilize; else load-ctrl update
  int formEleTangent(FE_Element*);   // K (+ (c/dt)M* if stabilize)   §4.2
  int formEleResidual(FE_Element*);  // -fint (- c M* deltaUstep if stabilize) §4.2
  int domainChanged();           // ArcLength vectors + phat probe + M* build
  int commit();                  // snapshot committed per-step state (§4.3)
  int revertToLastStep();        // restore snapshot (§4.3) — NOT a no-op
  // mutators (Layer B, script-driven) — §4.3
  int setArcLength(double); int reduceStep(double f); int increaseStep(double f);
  int sendSelf(...); int recvSelf(...); int Print(...);
};
```

### 4.2 The stabilization seam (`formEleTangent` / `formEleResidual`)

This is the whole trick and the one piece worth prototyping first. The
`StaticIntegrator` base returns `K` from `formEleTangent`. **Correction (review
§2.16):** `Newmark::formEleTangent` (CURRENT_TANGENT) forms **three** terms —
`addKtToTang(c1) + addCtoTang(c2) + addMtoTang(c3)` — not `c1·K + c3·M`. We reuse
just the `addKtToTang(1.0) + addMtoTang(cVisc/Δt)` **subset** (no damping term):

```
int LadrunoArcLength::formEleTangent(FE_Element *ele) {
  ele->zeroTangent();
  ele->addKtToTang(1.0);                            // stock static tangent
  if (stabilize) ele->addMtoTang(cVisc/dtPseudo);   // + (c/dt) M*   (M* = artificial)
  return 0;
}
int LadrunoArcLength::formEleResidual(FE_Element *ele) {
  this->StaticIntegrator::formEleResidual(ele);     // -fint   (+ sensitivity path, unstabilized only)
  if (stabilize)                                     // f_v = (c/dt) M* * deltaUstep
    ele->addM_Force(*deltaUstep, -cVisc/dtPseudo);  // GLOBAL vector; slices myID + applies getMass internally
  return 0;
}
```

- `addMtoTang(double)` / `addM_Force(const Vector&, double)` already exist on
  `FE_Element` (lines 70/85; used by the transient integrators) → **no element API
  change**. **`addM_Force` takes a GLOBAL-sized vector** and slices the element
  DOFs internally via `myID` ([`FE_Element.cpp:766-799`](../SRC/analysis/fe_ele/FE_Element.cpp)),
  so we pass the whole `deltaUstep` (review §2.3/§2.15/§2.18).
- **`M*` (review BLOCKER §2.5/§2.6):** `addM_Force`/`addMtoTang` apply
  `myEle->getMass()`, which is **consistent-by-default and density-scaled** — and
  a `StaticIntegrator` never assembles **nodal** mass (`formNodTangent` is a
  hard-error stub). So the *real* mass is unusable as `M*`. v1 instead builds a
  **unity-density artificial lumped `M*`** in `domainChanged()` and applies it
  through these same calls (i.e. it cannot rely on `getMass()`; the artificial
  mass is the integrator's own — design detail for the build plan: either a
  per-element unity-density lumped diagonal cached at `domainChanged`, or a global
  identity-scaled `M*`). **The `addM_Force`-with-`getMass` shortcut above is
  illustrative; the real path uses the artificial `M*`.**
- **Δu source (review §2.3):** `FE_Element` exposes **no** displacement accessor —
  the integrator owns the global per-step increment `deltaUstep` (stock
  `ArcLength` already accumulates it, `ArcLength.cpp:294`) and passes it in.
- **Sensitivity (review §2.11):** the stabilized branch must guard
  `if (stabilize && sensitivityFlag) error` and the unstabilized branch must
  reproduce stock `ArcLength::formEleResidual`'s `sensitivityFlag` override
  verbatim.

### 4.3 Layer B mutators + `revertToLastStep` override (no analysis clone)

`reduceStep(f)`: `arcLength2 *= f*f;` (with an `ℓ ≥ floor` guard) — pure scalar
mutation, no reallocation. **But the "state is preserved for free" claim is false
(review §2.10):** stock `ArcLength::newStep` overwrites `deltaLambdaStep`
(`ArcLength.cpp:167`) and `deltaUstep` (175-176) **before** `update()`'s `b24ac<0`
bail, and the corrector mutates them incrementally (294-296) before any failure.
So a failed step leaves these **polluted**. v1 therefore:

1. **snapshots** `{deltaLambdaStep, deltaUstep, currentLambda,
   signLastDeltaLambdaStep, arcLength2}` at `commit()`, and
2. **overrides `revertToLastStep()`** to restore them from that snapshot.

`reduceStep` then operates on the *restored* committed radius. A script drives:

```python
ok = analyze(1)
tries = 0
while ok != 0 and tries < maxTries:
    integ.reduceStep(0.5)     # operates on committed state restored by revertToLastStep
    ok = analyze(1)           # StaticAnalysis reverts (revertToLastCommit + revertToLastStep) before retry
    tries += 1
```

`StaticAnalysis::analyze` does call `revertToLastCommit()` + `revertToLastStep()`
on failure (verified, review §2 staticanalysis-revert: accurate) — but it returns
**immediately** after, so the *script* owns the retry loop (it is **not** in-engine
retry; review §2.21). The mutator is exposed to Tcl/Python via a small
`integrator`-object command (or `setParameter`), TBD in the build plan.

### 4.4 Build-control obligations (REQUIRED, same PR)

- `SRC/classTags.h`:
  `#define INTEGRATOR_TAGS_LadrunoArcLength 33004 // N. Mora-Bowen (Ladruno) — adaptive + viscous-stabilized arc-length; integrator band >=33000 (siblings 33000 ExplicitBathe, 33001 ExplicitDifferenceStatic, 33002 ExplicitBatheLNVD, 33003 CentralDifferenceLadruno). NB 33004 is independently reused in the ELE_TAG space (ELE_TAG_LadrunoIMKBeam2d=33004) — tag bands are PER-REGISTRY, no collision.`
  (review §2.17/§2.22: 33004 is the next free *integrator* tag but also live as an *element* tag — document the per-registry independence explicitly.)
- **Broker (review §2.14 — corrected):** the real dispatch is
  **`FEM_ObjectBrokerAllClasses::getNewStaticIntegrator`**
  (`SRC/.../FEM_ObjectBrokerAllClasses.cpp`), **not** `FEM_ObjectBroker.cpp` (whose
  `getNew*Integrator` are stubs returning 0) and **not** `getNewIncrementalIntegrator`.
  Add `case INTEGRATOR_TAGS_LadrunoArcLength: return new LadrunoArcLength(1.0);` —
  so the class **must provide a default/`recvSelf`-able ctor**.
- **Tcl + Python dispatch (review §2.19):** primary seam is
  `SRC/interpreter/OpenSeesCommands.cpp` (the `OPS_*` factory dispatch, shared by
  Python and modern Tcl); legacy `SRC/tcl/commands.cpp` only if the old path is
  wanted. One `integrator LadrunoArcLength …` branch → **`LEDGER_vanilla_files.md`
  row** (the only vanilla touch; mark with `// Ladruno` comment).
- `LEDGER_implementations.md`: new row (Integrator / 33004 / files / status / PR).
- Banner: add a line to `Ladruno_scripts/banner_features.txt` →
  `python Ladruno_scripts/patch_banner.py`.
- `Ladruno_scripts/stamp_headers.py`: add the new files to GLOBS + rerun (LADRUNO
  header stamp is non-optional for fork-authored files).

---

## 5. Public API (proposed)

```tcl
# identity — bit-identical to stock ArcLength:
integrator LadrunoArcLength $arc $alpha

# Layer A — Ramm adaptive radius:
integrator LadrunoArcLength $arc $alpha -adapt $Jd $ellMin $ellMax  ;# -p 1|0.5

# the viscous road — auto-calibrated stabilization:
integrator LadrunoArcLength $arc $alpha -stabilize 2.0e-4
integrator LadrunoArcLength $arc $alpha -stabilize 2.0e-4 -adaptStab

# everything:
integrator LadrunoArcLength $arc $alpha -adapt 5 1e-3 1e-1 -stabilize 2e-4
```

```python
integrator('LadrunoArcLength', arc, alpha,
           '-adapt', Jd, ellMin, ellMax,
           '-stabilize', 2.0e-4, '-adaptStab')
```

Knobs (all optional, default-off):
- `-adapt Jd ℓmin ℓmax` (+ `-p {1|0.5}`) — arc-length mode radius adaptation.
- `-stabilize f` (+ `-adaptStab`) — viscous road; `f` = dissipated-energy
  fraction (Abaqus default `2e-4`). **Mutually exclusive with `-adapt`'s
  constraint** (§decision 3b): `-stabilize` ⇒ regularized load control.
- mutators `setArcLength` / `reduceStep` / `increaseStep` via the integrator
  object command — Layer B.

---

## 6. Testing / oracle matrix (Zone-A)

> **Prototype status:** the snap-through fixtures (AL-1/AL-5/AL-6) are **already
> standing and green** in `tests/_proto_arch_snapthrough.py` (shallow von Mises
> truss, `corotTruss`, `E·A` tuned so the limit load ≈ 3.80). Verified pre-code:
> LoadControl **fails** at λ=3.80; DisplacementControl traces the **full path**
> incl. the unstable branch (limit λ=3.80119) and stock `ArcLength` cross-checks
> to 3.77 (<1%); a transient dynamic-relaxation run **snaps through** to the far
> stable branch (uy=−0.217) under a 1.15×crit load — the AL-6 stabilization
> oracle. The prototype becomes the pytest fixture when code lands.

| ID | Check | Oracle / pass |
|----|-------|---------------|
| AL-1 | `LadrunoArcLength arc alpha` (no flags) vs stock `ArcLength` on the von Mises truss | node-for-node, λ-for-λ **bit-identical**. **Qualified (review §2.23):** identity holds only when `adapt==false && stabilize==false`, with every added counter/snapshot guarded so it cannot perturb the disabled path. |
| AL-2 | Predictor `Δλ^{(1)}` formula | closed form `√(arcLength2/(û·û+α²))` |
| AL-3 | Layer A radius update | `arcLength2·(Jd/Jlast)^{2p}` clamped + **`Jlast←max(1,·)` guard** (§3.2) — unit test on the scalar rule |
| AL-4 | Adaptive run vs fixed-radius on the truss | adaptive reaches the same path in **fewer total factorizations** |
| AL-5 | **Snap-through through the limit point** — fixed `ArcLength` stalls/`b24ac<0`; `-stabilize` (regularized load control) converges with plain Newton | converges past λ=3.80; peak within tol **defined relative to `c`/step** (not absolute — §3.4 step-dependence) |
| AL-6 | … vs the **transient dynamic-relaxation** run (mass + Rayleigh) | far-branch displacement agrees within tol (the stabilized-static oracle; prototype: uy≈−0.217) |
| AL-5b | **Step-refinement / `c`-reduction convergence (review §2.13)** | halving the step (or `c`) moves the stabilized peak monotonically toward the AL-5-reference limit load — proves the result is `c`-controlled, not arbitrary |
| AL-7 | Energy pollution bound | `EnergyBalanceRecorder` viscous-diss/strain-energy ≤ a few × `f` |
| AL-8 | `-adaptStab` holds the dissipation ratio near `f` | running ratio within band |
| AL-9 | Layer-B `reduceStep` + `revertToLastStep` restores committed state through a failed step | `deltaLambdaStep`/`deltaUstep`/sign restored to the **committed** values (not the polluted in-flight ones — §4.3) |
| AL-10 | `sendSelf`/`recvSelf` round-trip (all flags) | state restored byte-faithful |
| AL-11 | **`M*` independence (review BLOCKER §2.5)** — a zero-density / nodal-`mass`-only model still stabilizes | `-stabilize` engages (artificial `M*` ≠ `getMass`); contrast: a `getMass`-based `M*` would give `M*=0` and no effect |

**Canonical model:** the shallow von Mises truss (smallest limit-point problem)
for AL-1…AL-8/AL-11 — **already built**; a single-element softening column
(ASDConcrete or Lemaitre-damaged brick) for the snap-**back** legs. Both tiny,
deterministic, Zone-A-portable.

---

## 7. Rejected: the cloned `StaticAnalysisAdaptive` driver

The original Layer-B idea was to clone `StaticAnalysis` into a fork driver owning
an in-engine cut-and-retry loop. **Rejected** — re-based after the review on the
cost argument as primary (review §2.21):

1. **(Primary) Cost / vanilla parity.** A new `analysis LadrunoStatic` type needs
   a dispatch branch in the interpreter parser regardless, and the payoff (copying
   ~200 lines of `analyze()` — a single non-virtual method over private members)
   is poor. The fork's leaf-only invariant is worth more than the convenience of
   an in-engine retry loop.
2. **(Covers the common case, with exceptions) The viscous road.** `-stabilize`
   makes `K_T` positive-definite through a **snap-through / limit point**, so
   ordinary Newton in the *stock* loop converges — Abaqus's `STABILIZE`-over-Riks
   logic. **Documented exceptions:** it does **not** help (a) zero-density /
   nodal-mass models if `M*` were `getMass`-based — fixed by the artificial `M*`
   (§2.5); (b) true **snap-back** (§3.4/§2.9). So it does not *universally* remove
   the need for cut-and-retry.
3. **(Covers the residual, with a caveat) The mutators.** `reduceStep()` +
   `revertToLastStep()` + the existing `revert`-on-failure in `StaticAnalysis`
   give script-driven retry with **zero analysis-core change**. Caveat (review
   §2.21): this is **script-owned**, not in-engine — `analyze(1)` returns to the
   script on the first failure; the loop lives in Python/Tcl, not C++.

If a future need for *fully* in-engine retry survives all three, revisit — the
`reduceStep`/`increaseStep`/`revertToLastStep` contract is deliberately
integrator-agnostic so it could later be lifted into a `StaticAnalysis` flag (one
vanilla row, reused by LoadControl/DisplacementControl too).

---

## 8. Follow-ups (deferred, ranked by softening value)

1. **Dissipation-/energy-release-controlled path following** (Gutiérrez 2004;
   Verhoosel–Remmers–Gutiérrez 2009) — a *second* `StaticIntegrator` that
   constrains incremental dissipated energy instead of a displacement norm. The
   one path-follower that genuinely beats geometric arc-length on localized /
   multi-crack softening. Same `newStep`/`update` seam; the high-value next step.
2. **Material viscous regularization** (Duvaut–Lions rate-dependent softening) —
   fixes the *physics* (mesh objectivity of the band), complementary to this
   ADR's *numerical* stabilization. Material-side.
3. **Indirect / CMOD control** — monotone relative-DOF constraint for clean
   snap-back without artificial dissipation.
4. **True-equilibrium residual for the convergence test (review §2.12)** — expose
   `‖λp − f_int‖` (without `f_v`) so the `ConvergenceTest` measures real static
   unbalance in stabilized mode, instead of the `f_v`-polluted SOE residual.
5. **Lift `reduceStep`/`revertToLastStep` contract into a stock `StaticAnalysis`
   flag** if true *in-engine* retry is ever wanted for all integrators (one
   vanilla row, reused by LoadControl/DisplacementControl too).

---

## 9. Relationship to other ladruno work

- **Explicit stack** ([[05_robust_central_difference]], explicit Bathe family,
  `Ladruno_explicit_roadmap`): the *other* answer to softening — brute-force
  quasi-static explicit, which sidesteps singular tangents entirely. This ADR is
  the **implicit** counterpart; the transient dynamic-relaxation oracle (AL-6) is
  the bridge between them.
- **`EnergyBalanceRecorder`** ([[04_explicit_dynamics_and_energy_balance]]): the
  mandatory watchdog for artificial-dissipation pollution (decision 6).
- **`LadrunoBrick` / Béziers**: consumers — softening solid/concrete
  ([[11_brick_asdconcrete_integration]]) and Lemaitre-damaged
  ([[15_lemaitre_ductile_damage_adr]]) elements are the intended payloads whose
  limit points this integrator is meant to clear. NB their hourglass `Kstab` is a
  *different* stabilization (§1.2) — do not conflate.
```
