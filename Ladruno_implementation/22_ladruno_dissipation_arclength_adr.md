---
title: ADR — Dissipation / energy-release controlled path-following (`LadrunoDissipationArcLength`)
project: Ladruno
status: proposed — scoped, no code
priority: medium
owner: nmora
tags:
  - implementation
  - integrator
  - static
  - arc-length
  - dissipation
  - energy-release
  - softening
  - localization
  - snap-back
  - adr
---

# ADR — Dissipation / energy-release controlled path-following (`LadrunoDissipationArcLength`)

**Status:** proposed — design scoped, **no code yet** (sibling of ADR-20
`LadrunoArcLength` and ADR-21 `LadrunoDynamicRelaxation`) ·
**Registry:** `StaticIntegrator` ·
**classTag:** `INTEGRATOR_TAGS_LadrunoDissipationArcLength = 33007` (next free in
the ladruno integrator band; siblings ExplicitBathe=33000,
ExplicitDifferenceStatic=33001, ExplicitBatheLNVD=33002,
CentralDifferenceLadruno=33003, LadrunoArcLength=33004,
LadrunoDynamicRelaxation=33005, **LadrunoIndirectControl=33006** — 33006 was
taken by the indirect/CMOD-control integrator *after* the first drafts of this
ADR were sketched, so the canonical free numeral is now **33007**, re-verified
against `../SRC/classTags.h` line 1131) ·
**Supersedes nothing** — adds a **standalone fork integrator** alongside stock
`ArcLength`, stock `LoadControl`, and the fork's own `LadrunoArcLength` (33004);
33007 leaves vanilla `INTEGRATOR_TAGS_ArcLength`, vanilla
`INTEGRATOR_TAGS_LoadControl`, and the 33004 / 33006 fork siblings *untouched* ·
**Oracle:** stock `LoadControl` (the mandatory phase-1 elastic predictor is
**bit-identical** to `LoadControl` when no localization has yet occurred) +
`LadrunoArcLength` (donor of ~70 % of the seam, verbatim) + a standing
`tests/_proto_dissipation_arclength.py` single-DOF closed-form softening oracle.

> [!note] As-designed scope (v1 target)
> One self-contained `LadrunoDissipationArcLength : public StaticIntegrator`
> **leaf class** — a *sibling*, not a subclass, of `LadrunoArcLength` (the two
> share copied plumbing, not an inheritance edge; see §1.2 and §4.3). It follows
> equilibrium past **limit points and severe snap-back** by stepping a fixed
> quantum of *dissipated energy* `Δτ` per increment (Gutiérrez 2004 /
> Verhoosel-Remmers-Gutiérrez 2009), so the constraint *sees only the failing
> material* and is insensitive to elastic unloading of the bulk — the regime
> where geometric (norm-based) arc-length picks the wrong root. The control
> quantity `Δτ` is computed **purely integrator-local** (no element / material /
> velocity query, no domain sweep). A small state machine runs a mandatory
> **phase-1 force-control** elastic start, then auto-switches to dissipation
> control once measurable energy is being released. Default invocation with no
> extra robustness flags runs the simplest path: phase-1 load control + a
> one-way latch into fixed-`Δτ` dissipation control.
>
> **Hard constraint honored:** `StaticAnalysis::analyze` is **not touched, not
> cloned, not subclassed**. The only vanilla edits are the standard
> integrator-registration rows in the broker (`getNewStaticIntegrator`) plus the
> Tcl/Python dispatch sites — the same N-site seam every ladruno integrator
> already uses (ExplicitBathe … LadrunoIndirectControl), ledgered in
> `LEDGER_vanilla_files.md`. The constraint, the bordered corrector, and the
> switch state machine live **entirely inside the leaf** `newStep`/`update`/
> `commit`/`revertToLastStep` overrides.
> > **Reviewer-confirmed caveats:** (a) the *cold start is singular* — at `u_0=0`
> > the bordered scalar `w = −½ f̂ᵀu_n = 0`, so dissipation control **cannot** be
> > the first step; the phase-1 force-control start is *load-bearing*, not a
> > nicety (BLOCKER, §3.4). (b) The Gutiérrez `Δτ` form is exact for the **damage**
> > constitutive class with a *constant* reference load `f̂` (no follower loads);
> > plasticity/cohesive dissipation uses a different expression and is a §8
> > follow-up. (c) The bordered tangent is **non-symmetric** even when `K` is
> > symmetric; the two-solve corrector sidesteps this by only ever factoring `K`
> > (§3.3, §7.3). (d) The switch state machine — not the algebra — is the
> > highest-scrutiny code; `mode` and the running scalar `phatDotU` **must** be in
> > the commit/revert/sendSelf payload or a `reduceStep` retry resumes in the
> > wrong phase with a stale `f̂ᵀu_n` (BLOCKER, §3.4).

---

## 1. Context

### 1.1 The phenomenon, and what is / isn't in scope

Quasi-static simulation of **softening / localizing** structures (cracking
concrete, decohesion, phase-field fracture, fiber-section strength loss) produces
load–displacement paths with **limit points** (load control fails: the load
would have to decrease) and **snap-back** (displacement control fails: the
controlled DOF is non-monotone). Geometric arc-length (spherical / cylindrical
Riks-Crisfield, `‖Δu‖² + α²Δλ² = ℓ²`) handles isolated limit points but
**breaks under sharp localization**: when one crack band grows while the elastic
bulk *unloads*, the displacement norm `‖Δu‖²` is dominated by the large
elastic-unloading motion of the bulk, *not* by the localizing DOFs. The
constraint surface becomes nearly tangent to the equilibrium path, the predictor
picks the wrong root, and Newton diverges or jumps branches. With **multiple
competing cracks** the global norm cannot distinguish which crack is active, and
the path has several sharp snap-backs the norm cannot separate.

**Dissipation control** (Gutiérrez 2004) replaces the geometric norm with a
constraint on *energy actually dissipated this step*. Elastic unloading
contributes **zero** dissipation, so the constraint is **insensitive to the
elastic bulk**, follows the localization regardless of how it redistributes
spatially, and **needs no a-priori knowledge of where the crack is** — the
original motivation versus indirect/CMOD control (which needs the crack location
pre-chosen; that is the *separate* `LadrunoIndirectControl`, 33006). Because
dissipation is monotone (`ḋ ≥ 0`), stepping by a fixed energy quantum gives
uniform "progress" along the failure path even through severe snap-back.

**In scope (v1):**
- A standalone `StaticIntegrator` leaf, `LadrunoDissipationArcLength`, classTag 33007.
- The **incremental dissipation constraint** `φ^D = ½(λ_n u_{n+1} − λ_{n+1} u_n)·f̂ − Δτ = 0` (Gutiérrez/Verhoosel; matches the brief's stated `g`).
- **Two-solve bordered corrector** (reuse the existing `K⁻¹f̂` predictor solve; bordered matrix *never* assembled; no `LinearSOE`/`Solver` change).
- **Integrator-local `Δτ`** sourcing from `phat`, `deltaUstep`, `deltaLambdaStep`, `currentLambda`, plus one new running scalar `phatDotU = f̂ᵀu`.
- The **phase-1 force-control → dissipation-control** switch state machine (cold-start cure) with the `Δτ^D > a·Δτ^U` switch criterion (default `a = 0.1`).
- Script-driven cut-and-retry via `reduceStep` (reuse the Layer-B snapshot/revert pattern verbatim).
- **Damage** constitutive class as the validated target.

**Out of scope (recorded as follow-ups in §8):**
1. **Switch-*back*** to force control on elastic reload / crack arrest (re-stiffening). v1 uses a *one-way latch* + a clean-fail denominator floor (§3.4, §8.1).
2. **Plasticity / cohesive** dissipation expression (plastic-work form, Verhoosel 2009) — v1 targets damage only.
3. **Geometrically-nonlinear** dissipation constraint variant (Verhoosel 2009 large-strain coefficients).
4. The May-Vignollet-De Borst (2016) **internal-energy `φ^U` dual control** as a *richer* elastic phase (v1 uses plain force control for phase 1).
5. **Adaptive `Δτ`** (iteration-count scaling, Crisfield-style) — v1 uses fixed `Δτ` + script `reduceStep`.
6. **Sherman-Morrison** rank-one bordered inverse (v1 uses the two-solve; same cost when `K` is factored).
7. **Energy-balance cross-check** wiring to `EnergyBalanceRecorder` as a watchdog.

### 1.2 Why a standalone sibling class, not a mode folded into `LadrunoArcLength`

The competing angle was to *fold* a `-dissipation Δτ` mode into the existing
`LadrunoArcLength` (33004), saving a classTag and all dispatch edits. **Rejected
(see §7.1), settled by source:**

| Claim under review | Verdict | Decisive evidence |
|---|---|---|
| The arc-length corrector and the dissipation corrector are the same algebra | **FALSE** | `LadrunoArcLength::update` solves a *quadratic* constraint — `a,b,c` coefficients, discriminant `b24ac<0` imaginary-root bail, two-root `dlambda1/2`, and a `theta1` angle pick (`../SRC/analysis/integrator/LadrunoArcLength.cpp:382–419`). The dissipation constraint is **linear in `δλ`** (a single explicit root, no discriminant, no angle test). Folding carries dead `a/b/c/b24ac/alpha2` members on the dissipation path. |
| Folding is "zero new vanilla edits" so it's strictly cheaper | **MISLEADING** | The send/recv payload is `Vector(27)` with exactly **one** free reserved slot `data(26)` (`LadrunoArcLength.cpp:675, 704`). Dissipation control needs ~8 new fields (`Δτ`, `Δτ^U` seed, `a`, `mode`, `phatDotU`, snapshots). Folding forces a `Vector(35)` widening *with* a legacy-read size guard on the shared 33004 class — a migration risk to a *shipped* integrator. A new leaf gets a clean payload for free. |
| Mode-flags keep the class readable | **FALSE at this size** | `LadrunoArcLength` already carries two mutually-exclusive modes (adaptive arc + `-stabilize` viscous, ADR-20 §). A third, *algebraically disjoint* mode (different constraint, different corrector, different predictor, different commit-time state machine) tips it from "superset of `ArcLength`" into a multiplexer. The house pattern is **leaf-per-method** (six integrator siblings already). |
| There is precedent for a *separate* class per control method | **TRUE** | `LadrunoIndirectControl` (33006) is itself a *separate* `StaticIntegrator` for a *different control quantity* (`c·U`), shipped as its own leaf rather than a `LadrunoArcLength` mode (`../SRC/classTags.h:1131`). Dissipation control is one more control quantity → one more leaf. |

**Decision: a new leaf, `LadrunoDissipationArcLength` (33007).** Leaf additions,
never core surgery; the *constraint* is what differs, and a different constraint
earns a different leaf.

---

## 2. Decision

1. **Build a standalone `LadrunoDissipationArcLength : public StaticIntegrator`,
   classTag `INTEGRATOR_TAGS_LadrunoDissipationArcLength = 33007`** — sibling of
   `LadrunoArcLength` (33004) and `LadrunoIndirectControl` (33006), *not* a
   subclass or a folded mode (review §1.2; rejected alternative §7.1). Registered
   via `getNewStaticIntegrator` (it is a `StaticIntegrator`, **not** the
   `getNewIncrementalIntegrator` path, **not** the broker stub that returns 0).

2. **The constraint is the incremental-dissipation form**
   `φ^D = ½(λ_n u_{n+1} − λ_{n+1} u_n)·f̂ − Δτ = 0` (Gutiérrez 2004 Eq. 19/21;
   matches the brief's stated `g` up to the trivial transpose `f̂ᵀu = uᵀf̂` and an
   overall sign convention). It is **bilinear** in the unknowns (committed
   `u_n, λ_n` against trial `u_{n+1}, λ_{n+1}`), hence **linear in `δλ`** at
   corrector level — the load-bearing simplification versus the quadratic
   arc-length constraint (review §3.2).

3. **Corrector = two-solve bordered Newton, bordered matrix never assembled.**
   Reuse the existing `phat`-solve idiom (`δu_f = K⁻¹f̂`, computed once per step,
   already in `LadrunoArcLength::newStep`/`update`). The dissipation update adds
   only `h = −½λ_n f̂` (vector), `w = −½ f̂ᵀu_n` (scalar), and the single explicit
   `δλ = −(φ^D + hᵀδu_r)/(w + hᵀδu_f)`. **No `LinearSOE`/`Solver` change**, and
   the non-symmetric border never reaches a symmetric solver (review §3.3, §7.3).

3b. **`Δτ` is sourced purely integrator-local — no element / material / velocity
   query.** `f̂ = phat`, `Δu = deltaUstep`, `Δλ = deltaLambdaStep`,
   `λ_n = currentLambda − deltaLambdaStep`, and the one missing scalar `f̂ᵀu_n` is
   carried as a *running accumulator* `phatDotU += phat ^ deltaUstep` updated at
   `commit()`. `EnergyBalanceKernel.h` is **explicitly rejected** as the source —
   it is velocity- and element-query-based, intended for transient analysis, and
   would re-introduce the per-element sweep this design avoids (review §3.2,
   ENERGY-report). The precedent already exists: `LadrunoArcLength::commit`
   forms `Estrain0 = ½|deltaUstepᵀ phat|·|currentLambda|` from exactly this
   triple (`LadrunoArcLength.cpp:550`).

4. **A mandatory two-phase predictor with a switch.** Phase 1 = **force
   control** (`φ^F = λ − Δτ^F`, `∂φ^F/∂λ = 1`, nonsingular from `u_0 = 0`).
   Auto-switch to dissipation control when the measured incremental dissipation
   crosses the threshold `Δτ^D > a·Δτ^U` (May-De Borst 2016 Eq. 29; `a ∈ (0,1)`,
   default `0.1`, smaller for more brittle response). The first dissipative `Δτ`
   is seeded from the last elastic step's internal-energy rate (review §3.5).
   This switch is what makes the bordered system **well-posed**; the cold start is
   singular without it (BLOCKER, review §3.4).

4b. **v1 ships a one-way latch (force → dissipation), not the full switch-back.**
   On elastic reload the measured dissipation rate `→ 0`, driving the bordered
   denominator `w + hᵀδu_f → 0`. v1 does **not** bounce back to force control
   (that is §8.1); instead it **fails cleanly** — when `|denom| < tol` the
   corrector returns `−1`, handing control to `StaticAnalysis::analyze`'s
   `revertToLastStep` + the script `reduceStep` retry, rather than emitting a
   silent NaN. Failing safe is the v1 contract (review §3.4).

5. **`formTangent` / `formUnbalance` are pure passthroughs to the base** — no
   `-stabilize` viscous road is carried over from `LadrunoArcLength` (the
   dissipation constraint *is* the regularizer; the viscous dashpot would be
   redundant and is omitted, review §7.2 / SEAM-report §4 "orthogonal/optional").

6. **Default-off / backward-compat = pure leaf addition.** There is no "reduces
   to `ArcLength`" identity (it's a different method invoked by a different
   command). The compatibility guarantee is: (a) the leaf touches no stock class;
   (b) the **phase-1 force-control predictor is bit-identical to stock
   `LoadControl`** until the first measurable dissipation, giving a clean identity
   gate (DA-1). Stock `ArcLength`, `LoadControl`, and the 33004 / 33006 fork
   siblings are byte-untouched (review §6).

7. **Same-PR build-control footprint** is the standard 9-item integrator seam
   (§4.4), with three new `LEDGER_quirks` rows (non-symmetric bordered tangent;
   cold-start singularity; one-way-latch ill-conditioning on reload).

---

## 3. The model

Notation (fixed): `f̂ = phat` = reference (unit) load vector; external load
`f_ext = λ f̂`; residual `r(u,λ) = f_int(u) − λ f̂`. Within increment `n+1`,
`u_n, λ_n` are the **committed** (start-of-step) values; `Δu = deltaUstep`,
`Δλ = deltaLambdaStep` are the step accumulators; iterate `i` gives corrections
`δu, δλ`. The residual displacement solve is `δu_r = K⁻¹(−r)` (the `dU` handed to
`update`); the reference-load solve is `δu_f = K⁻¹f̂` (the `deltaUhat` re-solve).

### 3.1 Stock `LoadControl` recap (the identity baseline)

In the elastic regime there is no dissipation, so v1's phase-1 predictor is plain
force control: prescribe `Δλ = Δτ^F`, set `Δu = Δλ · δu_f`, correct with
`δλ = 0` (residual-only corrector). This is *exactly* what stock `LoadControl`
does, so DA-1 can assert bit-identity on an elastic problem.

### 3.2 The dissipation constraint (Gutiérrez 2004)

For a damage law `σ = g(d) C:ε`, the global dissipation rate reduces to a form in
**global vectors only**:

$$
\dot{\mathcal{E}}^D \;=\; \tfrac{1}{2}\bigl(\lambda\,\dot{u}^{\mathsf T} - \dot\lambda\,u^{\mathsf T}\bigr)\hat f .
$$

Discretizing across the step (the midpoint-consistent bilinear form, exact for
the damage class) and equating to a prescribed energy quantum `Δτ`:

$$
\boxed{\;\varphi^D(u,\lambda) \;=\; \tfrac{1}{2}\bigl(\lambda_n\,u_{n+1}^{\mathsf T} - \lambda_{n+1}\,u_n^{\mathsf T}\bigr)\hat f \;-\; \Delta\tau \;=\; 0\;}
$$

This is the user's stated `g` (with `f = f̂`, `du ≡ u_{n+1}`, `dλ ≡ λ_{n+1}`).
**Sign conventions locked:** `Δτ > 0` always (dissipation is monotone; Gutiérrez
Eqs. 14, 16); the bilinear form pairs **committed** `u_n, λ_n` against **trial**
`u_{n+1}, λ_{n+1}`, which is precisely why the Jacobian row depends on committed
state (§3.3). The overall sign of `φ^D` is convention; what matters is `φ^D = 0`.

**Integrator-local evaluation** (the elegant part — *confirmed* against source).
Every term is already live in the seam (ENERGY-report §2):

| Gutiérrez symbol | seam variable |
|---|---|
| `f̂` reference load | `*phat` (built once in `domainChanged`) |
| `Δu` step disp increment | `*deltaUstep` |
| `Δλ` step load increment | `deltaLambdaStep` |
| `λ_n` step-start load factor | `currentLambda − deltaLambdaStep` (or committed `cCurrentLambda`) |
| `f̂ᵀu_n` (scalar only) | running accumulator `phatDotU`, updated `phatDotU += phat ^ deltaUstep` at `commit()` |

The only quantity *not* held as a literal vector is `u_n`, but `φ^D` needs only
the **scalar** `f̂ᵀu_n`, captured by `phatDotU`. Then `f̂ᵀΔu = phat ^ deltaUstep`,
giving:

$$
\Delta\tau_{\text{meas}} \;=\; \tfrac12\!\left[(\lambda_n)\,(\hat f^{\mathsf T}\!\Delta u) \;-\; (\Delta\lambda)\,(\hat f^{\mathsf T} u_n)\right]
\;=\; \tfrac12\!\left[(\texttt{currentLambda}-\texttt{deltaLambdaStep})(\texttt{phat}\,{}^\wedge\,\texttt{deltaUstep}) - \texttt{deltaLambdaStep}\cdot\texttt{phatDotU}\right].
$$

No `Element`/`NDMaterial`/`UniaxialMaterial` query, no domain sweep, no
velocities. This is exact only when `f̂` is a **constant** reference load (no
follower / displacement-dependent loads) — which `phat` is, computed once in
`domainChanged` (§3.4 caveat).

### 3.3 The bordered corrector (two-solve, linear in `δλ`)

Linearize the augmented system `H(u,λ) = [r ; φ^D] = 0` about iterate `i`:

$$
\underbrace{\begin{bmatrix} K & -\hat f \\ h^{\mathsf T} & w \end{bmatrix}}_{\text{bordered }K_T}
\begin{bmatrix}\delta u \\ \delta\lambda\end{bmatrix}
= -\begin{bmatrix} r \\ \varphi^D \end{bmatrix},
\qquad
h = \dfrac{\partial\varphi^D}{\partial u} = -\tfrac12\lambda_n\,\hat f,
\quad
w = \dfrac{\partial\varphi^D}{\partial\lambda} = -\tfrac12\,u_n^{\mathsf T}\hat f .
$$

(The committed `u_n, λ_n` are *constants* under this derivative — that is why `h`,
`w` use the start-of-step state.) **Two-solve (Riks/Batoz-Dhatt) decomposition**,
reusing the factored `K`:

$$
K\,\delta u_r = -r,\qquad K\,\delta u_f = \hat f,
\qquad
\boxed{\;\delta\lambda = -\dfrac{\varphi^D + h^{\mathsf T}\delta u_r}{\,w + h^{\mathsf T}\delta u_f\,}\;},
\qquad
\delta u = \delta u_r + \delta\lambda\,\delta u_f .
$$

`δu_f = K⁻¹f̂` is exactly the existing `deltaUhat` re-solve. **The bordered matrix
is never assembled** — only the symmetric `K` is ever factored, and the scalar
`δλ` is a closed form. **Single root, no discriminant, no two-root angle test, no
imaginary-root bail** — the entire `a/b/c/b24ac` + `theta1` block of
`LadrunoArcLength::update` (lines 382–419) is *deleted* for this constraint.

After computing `δλ`/`δu`, the seam continues verbatim: accumulate
`deltaUstep += δu`, `deltaLambdaStep += δλ`, `currentLambda += δλ`, apply
`incrDisp / applyLoadDomain / updateDomain`, and — *critically* — write the
augmented increment back with `theLinSOE->setX(*deltaU)` so the `ConvergenceTest`
measures the full path-following step (SEAM-report §1, line 439).

### 3.4 Honest limits (reviewer-confirmed)

- **Cold start is singular (BLOCKER).** At `u_0 = 0`, `w = −½ f̂ᵀu_n = 0`, and the
  predictor dissipation `∂φ^D/∂λ|_{n=1} = ½ u_1ᵀf̂ = 0`. The bordered denominator
  is `0` and dissipation control is *undefined* for the first step. **Cure:**
  the mandatory phase-1 force-control start (§3.5). This is not optional polish;
  it is what makes the method well-posed.
- **Switch state machine is the highest-risk code.** `mode` (LOAD vs
  DISSIPATION) and `phatDotU` **must** be in the commit / revert / sendSelf
  payload. If they are not, a `reduceStep` cut-and-retry (which calls
  `revertToLastStep`) resumes in the wrong phase with a stale `f̂ᵀu_n`, silently
  corrupting the constraint. (BLOCKER — the #1 implementation pitfall across all
  three source proposals.)
- **One-way latch on reload (v1 limitation, fails safe).** v1 does not switch
  *back* to force control when a crack arrests and the structure re-stiffens
  (`Δτ_meas → 0`, denominator ill-conditioned). v1 detects `|denom| < tol` and
  returns `−1` (clean fail → revert/`reduceStep`), rather than producing garbage.
  Full switch-back is §8.1.
- **Constitutive scope.** The `½(λΔuᵀf̂ − Δλ f̂ᵀu_n)` form is exact for **damage**
  with a *constant* `f̂`. **Plasticity / cohesive** dissipation is plastic work, a
  *different* expression (Verhoosel 2009) — §8.2. **Follower / geometrically
  nonlinear** loads break the `f̂`-constant assumption — §8.3.
- **Non-symmetric bordered tangent.** `K_T` is non-symmetric (`h ≠ −f̂`, `w ≠ 0`)
  even when `K` is symmetric. The two-solve sidesteps this (only `K` is factored);
  any future *assembled* bordered system would need a non-symmetric solver
  (logged as a `LEDGER_quirks` row, §7.3).

### 3.5 The predictor and the load→dissipation switch

1. **Phase 1 — force control** from the first step: `φ^F = λ − Δτ^F`, `δλ = 0`
   corrector (residual-only), bit-identical to `LoadControl` (§3.1). Nonsingular
   from `u_0 = 0`.
2. **Switch criterion.** After each elastic step, measure `Δτ^D` (§3.2) and the
   internal-energy increment `Δτ^U = ½ λ_{n+1} u_{n+1}ᵀf̂` (the precedent
   `Estrain0`-style scalar). Switch LOAD→DISSIPATION when

   $$
   \Delta\tau^D \;>\; a\,\Delta\tau^U, \qquad a \in (0,1)\ \text{(default }0.1\text{; smaller }\Rightarrow\text{ more brittle / earlier switch)} .
   $$

3. **Seed `Δτ`.** The first dissipative step uses `Δτ = a·Δτ^U` from the last
   elastic step (May-De Borst Eq. 28/29).
4. **v1 latch:** once in DISSIPATION, stay (one-way), with the `|denom| < tol`
   clean-fail guard. Switch-*back* (`Δτ^D ≤ 0` ⇒ return to LOAD with hysteresis)
   is §8.1.

The commit payload therefore carries: `u_n` (via `phatDotU`), `λ_n`
(`cCurrentLambda`), the active `mode` flag, `Δτ`, the `Δτ^U` seed, and `a`.

---

## 4. Architecture

### 4.1 Class skeleton

```cpp
// SRC/analysis/integrator/LadrunoDissipationArcLength.h
class LadrunoDissipationArcLength : public StaticIntegrator {   // classTag 33007
 public:
  LadrunoDissipationArcLength();                                 // recvSelf-able default ctor
  LadrunoDissipationArcLength(double dTau, double aRatio,
                              double dLambda0F);                 // Δτ, switch ratio a, phase-1 force incr
  ~LadrunoDissipationArcLength();

  // --- StaticIntegrator / IncrementalIntegrator seam (overrides) ---
  int newStep(void);              // predictor: phase-1 force OR dissipation predictor
  int update(const Vector &dU);   // corrector: δλ=0 (LOAD) | bilinear δλ (DISSIPATION)
  int domainChanged(void);        // REUSE verbatim: vector sizing + phat probe + zero-load guard
  int commit(void);               // update phatDotU; evaluate switch; snapshot
  int revertToLastStep(void);     // restore snapshot (mode + phatDotU + step accumulators)
  int formTangent(int flag = CURRENT_TANGENT); // pure passthrough to base (no -stabilize road)
  int formUnbalance(void);                     // pure passthrough to base

  // --- object/parallel ---
  int sendSelf(int cTag, Channel &);   // flat Vector pack (widened payload, legacy size guard)
  int recvSelf(int cTag, Channel &, FEM_ObjectBroker &);
  void Print(OPS_Stream &, int flag = 0);

  // --- script mutators (mirror LadrunoArcLength accessor/mutator pattern) ---
  int reduceStep(double f);            // dTau *= f (floor-guarded); the cut-and-retry hook
  int increaseStep(double f);          // dTau *= f (cap-guarded)
  int setDissipationIncrement(double); // dTau = v
  double getDissipationIncrement(void) const;   // Δτ
  double getDissipation(void) const;            // accumulated Σ Δτ
  double getCurrentLambda(void) const;

 protected:
  // --- reference-load + path-following vectors (REUSE from LadrunoArcLength) ---
  Vector *phat;          // reference load, built in domainChanged (unit-λ unbalance probe)
  Vector *deltaUhat;     // K^{-1} f̂  (δu_f)
  Vector *deltaUbar;     // K^{-1}(−r) (δu_r), = dU handed to update()
  Vector *deltaU;        // augmented increment of this iterate
  Vector *deltaUstep;    // per-step displacement accumulator (= Δu)
  double  deltaLambdaStep;   // per-step load-factor accumulator (= Δλ)
  double  currentLambda;     // end-of-iterate load factor
  int     nEqn;

  // --- dissipation-control state (NEW) ---
  double  dTau;          // Δτ  target energy quantum (the control parameter)
  double  aRatio;        // a   switch ratio (default 0.1)
  double  dLambda0F;     // phase-1 force-control increment
  double  phatDotU;      // running f̂ᵀu  (the only new "vector-like" state, as a scalar)
  double  dissipTotal;   // Σ Δτ  (diagnostic)
  int     mode;          // MODE_LOAD | MODE_DISSIPATION

  // --- Layer-B committed snapshot (REUSE pattern; carries mode + phatDotU) ---
  double  cCurrentLambda, cDeltaLambdaStep, cPhatDotU, cDTau;
  int     cMode;
  Vector *cDeltaUstep;
  bool    haveSnapshot;
};
```

### 4.2 The `update()` constraint seam

```cpp
int LadrunoDissipationArcLength::update(const Vector &dU) {
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE     *theSOE   = this->getLinearSOE();

  (*deltaUbar) = dU;                                  // δu_r = K^{-1}(−r)         [REUSE]
  theSOE->setB(*phat); theSOE->solve();              // δu_f = K^{-1} f̂           [REUSE]
  (*deltaUhat) = theSOE->getX();

  double dLambda;
  if (mode == MODE_LOAD) {                            // phase-1: == LoadControl
    dLambda = 0.0;                                    // residual-only corrector
  } else {                                            // MODE_DISSIPATION
    double lambda_n = currentLambda - deltaLambdaStep;          // λ_n
    // h^T x  =  (−½ λ_n f̂)^T x  =  −½ λ_n (phat ^ x)
    double hTdur = -0.5 * lambda_n * ((*phat) ^ (*deltaUbar));  // hᵀδu_r
    double hTduf = -0.5 * lambda_n * ((*phat) ^ (*deltaUhat));  // hᵀδu_f
    double w     = -0.5 * phatDotU;                             // w = −½ f̂ᵀu_n
    double phiD  = 0.5*( lambda_n*((*phat)^(*deltaUstep))       // φ^D = Δτ_meas − Δτ
                         - deltaLambdaStep*phatDotU ) - dTau;
    double denom = w + hTduf;
    if (fabs(denom) < DTAU_DENOM_TOL) return -1;      // clean-fail → revert/reduceStep
    dLambda = -(phiD + hTdur) / denom;                // single explicit root (linear in δλ)
  }

  (*deltaU) = (*deltaUbar);
  deltaU->addVector(1.0, *deltaUhat, dLambda);        // δu = δu_r + δλ δu_f        [REUSE]
  (*deltaUstep) += (*deltaU);
  deltaLambdaStep += dLambda;
  currentLambda   += dLambda;

  theModel->incrDisp(*deltaU);                        // apply to domain           [REUSE]
  theModel->applyLoadDomain(currentLambda);
  theModel->updateDomain();
  theSOE->setX(*deltaU);                              // augmented increment → ConvTest [REUSE]
  return 0;
}
```

### 4.3 Reuse map (`LadrunoArcLength` ≈ 70 %)

**REUSE verbatim** (copy unchanged from `LadrunoArcLength`, confirmed
SEAM-report §4):
- `domainChanged()` *in full* — vector sizing from `getNumEqn()`, the unit-λ
  `phat` unbalance probe (`LadrunoArcLength.cpp:492–513`), the **zero-load
  guard**, `nEqn` caching, snapshot seeding.
- The predictor `phat`-solve idiom (`formTangent(); setB(*phat); solve();
  deltaUhat = getX();`).
- The corrector residual + `phat` re-solve idiom and `setX(*deltaU)` write-back.
- All domain-application calls (`incrDisp / applyLoadDomain / updateDomain`,
  `getCurrentDomainTime / setCurrentDomainTime`).
- The Layer-B **commit-snapshot / revert** pattern and `reduceStep /
  increaseStep / setX-of-control-parameter` mutators — `dTau` plays the exact
  role `arcLength2` plays in 33004 (and inherits the rule: **do not restore the
  control parameter on revert**, so `reduceStep` actually shrinks).
- The `sendSelf / recvSelf` flat-`Vector` pack/unpack skeleton (resized + remapped
  for the new fields; a legacy size guard since the 33004 payload was `Vector(27)`).

**NEW (dissipation-specific):**
- `phatDotU` running accumulator + its update at `commit()`.
- The `φ^D / h / w / δλ` block in `update()` (§4.2) — *replaces* the quadratic
  `a/b/c/b24ac/theta1` block.
- The phase-1 → dissipation `mode` state machine + switch criterion in `commit()`.
- `dTau`, `aRatio`, `dLambda0F`, `mode` members and their script mutators.

**OMITTED (carried by 33004, not wanted here):** the entire `-stabilize` viscous
road (`buildArtificialMass`, `Mstar/cVisc/cOverDt`, the `formTangent/formUnbalance`
diagonal pokes) — the dissipation constraint *is* the regularizer; the viscous
dashpot is redundant (§7.2). `formTangent/formUnbalance` are therefore pure base
passthroughs.

### 4.4 Build-control obligations (REQUIRED, same PR)

The implementation PR — based on **`ladruno`**, not `main` — must, in the *same
PR*:

1. **`../SRC/classTags.h`** — add
   `#define INTEGRATOR_TAGS_LadrunoDissipationArcLength 33007 // N. Mora-Bowen (Ladruno) — dissipation / energy-release controlled arc-length; ladruno integrator band >=33000 (siblings 33000 ExplicitBathe … 33005 LadrunoDynamicRelaxation, 33006 LadrunoIndirectControl); StaticIntegrator. Tag bands are PER-REGISTRY — 33007 in the integrator registry does not collide with any ELE_TAG/MAT_TAG/ND_TAG/RECORDER_TAGS use of 33007.`
2. **Broker** — `FEM_ObjectBrokerAllClasses::getNewStaticIntegrator` (this is a
   `StaticIntegrator`; **NOT** `getNewTransientIntegrator`, **NOT**
   `getNewIncrementalIntegrator`, **NOT** the `FEM_ObjectBroker.cpp` stub that
   returns 0). The default ctor must be `recvSelf`-able.
3. **Tcl + Python dispatch** — `../SRC/interpreter/OpenSeesCommands.{h,cpp}` (the
   `OPS_*` factory, shared by Python + modern Tcl),
   `../SRC/runtime/runtime/TclPackageClassBroker.cpp`, and legacy
   `../SRC/tcl/commands.cpp` if the old path is wanted — the same N-site seam
   every ladruno integrator uses ("CDL's 8-site seam plus ledger/banner/stamp").
4. **`../SRC/analysis/integrator/CMakeLists.txt` + `Makefile`** — add the
   `.cpp/.h/.o`.
5. **`LEDGER_vanilla_files.md`** — one row per vanilla file touched (broker + each
   dispatch site), each in-source edit marked with a `// Ladruno …` comment.
6. **`LEDGER_implementations.md`** — new row: feature
   `LadrunoDissipationArcLength` / kind=Integrator / classTag 33007 / files /
   status / PR.
7. **`LEDGER_quirks.md`** — three new rows: (a) **non-symmetric bordered tangent**
   (a symmetric solver breaks an *assembled* bordered system; two-solve avoids
   it); (b) **cold-start singularity** (`w = −½ f̂ᵀu_n = 0` at `u_0=0` ⇒ phase-1
   force control mandatory); (c) **one-way-latch ill-conditioning** on elastic
   reload (`denom → 0`, v1 clean-fails).
8. **Banner** — add a line to `Ladruno_scripts/banner_features.txt`, then
   `python Ladruno_scripts/patch_banner.py` (regenerates `FEATURES-START/END` in
   `../SRC/tcl/tclMain.cpp` + `../SRC/interpreter/PythonModule.cpp` — **do not
   hand-edit**). Every `shipped` impl-ledger row needs a matching banner line.
9. **`Ladruno_scripts/stamp_headers.py`** — add the new fork-authored files to
   GLOBS + rerun (the LADRUNO ASCII + four-author header is non-optional;
   `--check` for CI).

---

## 5. Public API (proposed)

**Tcl:**
```tcl
# Dissipation-controlled arc-length, fixed energy quantum dTau, switch ratio a:
integrator LadrunoDissipationArcLength $dTau                ;# minimal (a defaults 0.1)
integrator LadrunoDissipationArcLength $dTau -a 0.05        ;# more brittle, earlier switch
integrator LadrunoDissipationArcLength $dTau -forceIncr $dLambda0F   ;# phase-1 force step

# script-driven cut-and-retry (relies on StaticAnalysis revertToLastStep on failure):
set ok [analyze 1]
while {$ok != 0} { LadrunoDissipationArcLength reduceStep 0.5; set ok [analyze 1] }
```

**Python (OpenSeesPy):**
```python
ops.integrator('LadrunoDissipationArcLength', dTau)                 # minimal
ops.integrator('LadrunoDissipationArcLength', dTau, '-a', 0.05)
ops.integrator('LadrunoDissipationArcLength', dTau, '-forceIncr', dLambda0F)

ok = ops.analyze(1)
while ok != 0:
    ops.LadrunoDissipationArcLength('reduceStep', 0.5)
    ok = ops.analyze(1)
```

**Knobs (all optional, default-off / safe defaults):**
- `dTau` *(required)* — target dissipation increment per step (energy units; physical, so bounds are interpretable).
- `-a <ratio>` — LOAD→DISSIPATION switch ratio `a ∈ (0,1)`, default `0.1`; smaller ⇒ switch sooner / more brittle.
- `-forceIncr <dLambda0F>` — phase-1 force-control load increment (default: a small fraction of unit load).
- `reduceStep f` / `increaseStep f` / `setDissipationIncrement v` — runtime `Δτ` mutators (the cut-and-retry hook).

---

## 6. Testing / oracle matrix (Zone-A)

> **Prototype status:** `tests/_proto_dissipation_arclength.py` stands and runs
> green *today* against a hand-integrated single-DOF closed-form softening oracle
> (the DA-3 reference), and becomes the pytest fixture the moment 33007 lands.
> The identity-gate DA-1 reuses a stock `LoadControl` elastic run as its oracle.

| ID | Check | Oracle / pass |
|---|---|---|
| **DA-1** | **Identity gate** — elastic problem, phase-1 force control only | byte-identical `u(λ)` vs stock `LoadControl` (no localization ⇒ no switch); proves the leaf is a no-op until dissipation starts |
| **DA-2** | `Δτ_meas` accumulator correctness | `phatDotU`-based `½(λ_n f̂ᵀΔu − Δλ f̂ᵀu_n)` matches a brute-force `½(λ_n f̂ᵀu_{n+1} − λ_{n+1} f̂ᵀu_n)` recomputed from full state, < 1e-10 |
| **DA-3** | **Single-DOF softening bar (closed form)** — canonical oracle | traced `(u,λ)` incl. the snap-back loop matches the analytic damage backbone < 1e-6; `Σ Δτ` matches analytic dissipated energy |
| **DA-4** | Cold-start well-posedness | analysis advances from `u_0=0` with no singular-solve warning; phase-1 active for step 1 (proves the BLOCKER cure) |
| **DA-5** | Switch firing | `mode` flips LOAD→DISSIPATION at `Δτ^D > a·Δτ^U`; sweeping `a` smaller fires earlier |
| **DA-6** | **L-shaped panel** (single dominant crack, snap-back) | mesh-objective peak load (< ~2 %) and total dissipated energy (< ~5 %) across one 2× refine; stable softening branch traced |
| **DA-7** | **SEN beam / Nooru-Mohamed** (curved mixed-mode crack) | the case geometric arc-length *fails* — DA-7 traces the full softening curve; cross-check: stock `ArcLength` diverges on the same model (the discriminating test) |
| **DA-8** | One-way-latch clean-fail on reload | on a re-stiffening model, `|denom|<tol` ⇒ `update` returns `−1` ⇒ `revertToLastStep`/`reduceStep` engage (no NaN, no silent garbage) |
| **DA-9** | `reduceStep` retry / snapshot round-trip | after a forced non-converge, `revertToLastStep` restores `mode`, `phatDotU`, `deltaUstep`, `currentLambda`; the next `reduceStep`-shrunk step converges (proves the snapshot carries the state machine) |
| **DA-10** | `sendSelf/recvSelf` round-trip | pack→unpack reproduces all dissipation + Layer-B state bit-identically (widened payload + legacy size guard) |

**Canonical model:** the **single-DOF closed-form softening bar (DA-3)** is the
canonical oracle — it is the smallest model with an analytic snap-back and an
analytic dissipated-energy budget, so it pins *both* the path and the `Δτ`
accounting in one check.

---

## 7. Rejected alternatives

### 7.1 Fold a `-dissipation` mode into `LadrunoArcLength` (33004) — REJECTED

Tempting (no new classTag, no dispatch edits), but the constraint algebra is
*disjoint*: dissipation is **linear in `δλ`** (single explicit root) while
arc-length is **quadratic** (discriminant + two roots + `theta1` angle pick,
`LadrunoArcLength.cpp:382–419`). Folding strands dead `a/b/c/b24ac/alpha2`
members on the dissipation path and forces a payload widening (`Vector(27)`→`(35)`)
*on a shipped class* with a legacy-read migration guard. The house pattern is
**leaf-per-control-method** — `LadrunoIndirectControl` (33006) is itself a
separate leaf for a separate control quantity. *Leaf-only invariant: a different
constraint earns a different leaf.* **REJECTED** (review §1.2).

### 7.2 Carry the `-stabilize` viscous road into the dissipation leaf — DEFERRED (not rejected)

`LadrunoArcLength` injects an artificial dashpot `c·M*·Δu/Δt` to hold the tangent
positive-definite through limit points. For a dissipation integrator this is
**redundant** — the dissipation constraint already follows the negative-tangent
branch directly. Carrying it would mean two regularizers fighting. Omitted from
v1; could be revisited only if a pathological problem needs *both*. The
SEAM-report classifies it as "orthogonal/optional" and recommends omission.
**DEFERRED (not rejected).**

### 7.3 Assemble the full bordered system + use a non-symmetric solver — REJECTED (unnecessary)

`K_T` is genuinely non-symmetric (`h ≠ −f̂`, `w ≠ 0`), so an assembled bordered
system would require swapping in a non-symmetric solver — a `LinearSOE`/`Solver`
change touching core. The **two-solve** decomposition reuses the *symmetric* `K`
factorization and computes `δλ` as a closed-form scalar, so no solver change and
no non-symmetry trap. (Sherman-Morrison is the equivalent-cost alternative but
maps less cleanly onto the `SOE`/`Solver` split — §8.6.) **REJECTED** as
unnecessary core surgery; the quirk is logged so nobody re-discovers it (§4.4-7).

### 7.4 Source `Δτ` from `EnergyBalanceKernel.h` — REJECTED (wrong machinery)

`EnergyBalanceKernel.h` is a domain-sweeping, **velocity- and element-query-based**
transient energy-balance recorder (KE/IE/DW/ULW from per-element `getMass`,
`getDamp`, `getResistingForce` and nodal velocities). Quasi-static continuation
has no velocities, and pulling it in re-introduces the per-element sweep the whole
design avoids. The Gutiérrez `Δτ` needs *only* `{phat, deltaUstep,
deltaLambdaStep, currentLambda, phatDotU}` — all integrator-local. **REJECTED.**

---

## 8. Follow-ups (deferred, ranked by value)

1. **Switch-*back* to force control on reload** — the full hysteretic state
   machine (May-De Borst `InternalEnergyNegative` guard) so re-stiffening / crack
   arrest does not clean-fail. Seam: the `commit()` switch logic. *(Highest value
   — removes the one v1 robustness gap.)*
2. **Plasticity / cohesive dissipation constraint** — the plastic-work expression
   (Verhoosel 2009) so the method works on `LadrunoJ2` / Lemaitre / cohesive
   interfaces, not only damage. Seam: the `φ^D / h / w` block in `update()`.
3. **Geometrically-nonlinear constraint variant** — Verhoosel 2009 large-strain
   coefficients (the `f̂`-constant assumption relaxed). Seam: `φ^D` evaluation +
   `domainChanged` `phat` probe.
4. **May-De Borst internal-energy `φ^U` dual control** as a richer elastic phase
   (smoother switch than plain force control). Seam: phase-1 predictor in
   `newStep`.
5. **Adaptive `Δτ`** — `Δτ_new = Δτ_old·(N_desired/N_last)^p` (Crisfield), with
   `[Δτ_min, Δτ_max]` clamps to prevent stall/overshoot. Seam: `commit()` (the
   Layer-A Ramm hook ports directly to drive `dTau`).
6. **Sherman-Morrison rank-one bordered inverse** — single-solve alternative to
   the two-solve (equivalent cost when `K` is factored). Seam: `update()`.
7. **Energy-balance cross-check** — wire `EnergyBalanceRecorder` as a *watchdog*
   that `Σ Δτ` matches the recorder's global dissipation. Seam: a `getResponse`
   on the integrator + a test, no core change.

---

## 9. Relationship to other ladruno work

- **`LadrunoArcLength`** ([[20_ladruno_arclength_stabilized_adr]]): the donor of
  ~70 % of the seam (verbatim `domainChanged`/`phat` probe, predictor & residual
  solve idioms, `setX` write-back, Layer-B snapshot/`reduceStep`/revert,
  `sendSelf/recvSelf`). 33004 stays a *quadratic*-constraint sibling; 33007 is its
  *linear*-constraint sibling. They share code by copy, not inheritance (§1.2,
  §4.3).
- **`LadrunoDynamicRelaxation`** ([[21_ladruno_dynamic_relaxation_adr]]): the
  *other* way past limit points — pseudo-transient kinetic damping rather than a
  static constraint. Complementary: DR is matrix-free / explicit-friendly and
  needs no tangent; dissipation arc-length is implicit and traces the *exact*
  static path with energy as the abscissa. A model that defeats one often yields
  to the other.
- **`LadrunoIndirectControl`** (33006): the *sibling control quantity* — indirect
  / CMOD control follows snap-back via a *pre-chosen* weighted DOF set `c·U`;
  dissipation control needs **no a-priori crack location** and is insensitive to
  the elastic bulk. Different control quantity ⇒ separate leaf (§1.2 precedent).
- **`EnergyBalanceRecorder`** ([[04_explicit_dynamics_and_energy_balance]]):
  intentionally *not* the `Δτ` source (§7.4), but the natural §8.7 **watchdog**
  for cross-checking `Σ Δτ` against the global energy budget.
- **Softening consumers** — `LadrunoBrick` + ASDConcrete3D
  ([[11_brick_asdconcrete_integration]]), Lemaitre ductile damage
  ([[15_lemaitre_ductile_damage_adr]]), and the Bézier elements: these are the
  materials/elements whose snap-back this integrator is *for*. The v1 damage
  scope aligns directly with the ASDConcrete3D / crack-band consumers; the §8.2
  plasticity follow-up aligns with Lemaitre and `LadrunoJ2`.
