---
title: Robust explicit central-difference integrator (CentralDifferenceLadruno)
project: Ladruno
status: ready-to-implement
priority: high
tags:
  - implementation
  - integrator
  - explicit
  - central-difference
---

# Robust explicit central-difference integrator (`CentralDifferenceLadruno`)

> **Design / ADR (pre-implementation).** A single, clean **explicit leap-frog
> central-difference** integrator — a sibling-fork class (classTag 33003, Ladruno
> private band ≥33000) that delivers
> the one combination *no* existing OpenSees class has: a **correct first step + a
> built-in critical-timestep guard + clean full-step velocity output + energy-balance
> discipline**, without modifying any upstream file.
>
> **Status 2026-05-30: design approved, then NARROWED after an adversarial sweep.**
> The sweep (4 independent reviewers + targeted code verification) established that a
> *coupled* (C-on-LHS) mode would have duplicated the existing `NewmarkExplicit(γ=½)`
> — same effective tangent (`addCtoTang(c2); addMtoTang(c3)`), no starter defect,
> consistent-mass capable. So the **coupled mode was dropped**: for the
> implicit-damped / consistent-mass central-difference case, **use the existing
> `integrator NewmarkExplicit 0.5`**. This class is the *explicit leap-frog* CD only —
> the genuinely-missing, non-redundant piece. The numerical core was independently
> re-derived and verified sound; the surviving fixes (starter location, leap-frog
> state, registration list) are folded in below.

## What

A new `TransientIntegrator` subclass `CentralDifferenceLadruno` (classTag **33003**,
Ladruno private band ≥33000),
the explicit **leap-frog** central-difference scheme done right:

- **Single explicit scheme** (no mode switch): `M` alone on the LHS → a trivial
  diagonal `M⁻¹` solve every step. Viscous damping enters the residual at the known
  half-step velocity (the LS-DYNA-style genuinely-explicit form, same family as the
  existing `ExplicitDifference`).
- **Correct first step.** Half-step starter `v₋½ = v₀ − ½Δt·a₀` with
  `a₀ = M⁻¹(P₀ − C v₀ − Fᵢₙₜ(u₀))`, computed on the **first step** (where Δt is set
  and the SOE is formed — *not* in `domainChanged()`; see C3/B1). Kills the legacy
  `assuming Ut-1 = Ut` error.
- **`dt_cr` built in.** Reuse `CriticalTimeStep::computeCriticalTimeStep` verbatim;
  report the CD limit `2/ω_max` (stability factor **1.0**, no Noh–Bathe 2× bonus).
  Implement the `getCriticalTimeStep()` override; the `criticalTimeStep()` Py/Tcl
  command already exists and dispatches to it. Same `-cfl / -cflAbort / -tangent /
  -recompute N / -lump` surface as `ExplicitBathe`.
- **Two velocities, kept separate.** (a) **Node/recorder velocity** = full-step
  `vₙ = ½(v₋½ + v₊½)` (≡ `(uₙ₊₁−uₙ₋₁)/2Δt`), pushed via `setVel()`/`setResponse()` —
  fixes the legacy half-step/offset output. (b) **`getVel()`** = the modal-damping
  hook (`IncrementalIntegrator::addModalDampingForce` → residual via `setB`), which
  returns the half-step `v₋½` known at solve time.
- **βK guard.** Warn (with the dt_cr-collapse explanation) when stiffness-proportional
  Rayleigh (`betaK`/`betaKi`/`betaKc`≠0) is used, and report the damping-reduced
  `damped_dt`. Mass-proportional `αM` is safe.
- **Energy-balance ready.** Works with `EnergyBalanceRecorder` (classTag 26); the
  damping-work term `DW` integrates against the same lagged `getVel()` velocity.

**Explicitly NOT in scope:**
- **Coupled / implicit-damped CD** (C on the LHS, consistent-mass with a real solve):
  that scheme already ships as **`NewmarkExplicit(0.5)`** — use it; we do not rebuild it.
- Mass scaling (roadmap §5.1); sub-cycling inside the integrator (Bathe-ADR D5 —
  stays driver-level); batch/SoA dispatch (§5.2); **touching any upstream class**
  (`ExplicitDifference`, `NewmarkExplicit`, the 4 other CD classes stay frozen — the
  sibling-fork policy is the reason this is a new class and not an in-place patch).

## Why

The explicit CD landscape in OpenSees, *including* the one the original plan missed:

| Class | tag | First step | Damping | dt_cr | Velocity out | Mass |
|---|---|---|---|---|---|---|
| `CentralDifference` | 5 | warns `Ut-1=Ut` | Rayleigh flags stored, unused | none | full-step | lumped |
| `CentralDifferenceAlternative` | 17 | V from committed | none | none | half-step | lumped |
| `CentralDifferenceNoDamping` | 18 | V from committed | none | none | half-step | lumped |
| `ExplicitDifference` | 55 | warns `Ut-1=Ut` | Rayleigh via residual (lagged) | none | half-step→crude reconstruct | lumped |
| `ExplicitDifferenceStatic` | 62 | — | FLAC local | none | — | lumped |
| **`NewmarkExplicit`** | **19** | **clean (2-level, no defect)** | **C on LHS (proper)** | **none** | **full-step** | **any** |

Honest delta: the *explicit leap-frog scheme itself* already exists (`ExplicitDifference`),
and the *coupled/implicit-damped* scheme already exists and is robust on most axes
(`NewmarkExplicit`). What is **missing from every one of them** is the combination:
**(1) a correct explicit starter, (2) a `dt_cr` guard, (3) clean full-step velocity
output, (4) energy-balance discipline.** The Noh–Bathe work already built the
infrastructure for (2) and (4) (`CriticalTimeStep`, `EnergyBalanceRecorder`); this
class brings the explicit leap-frog scheme up to that standard. Because the
sibling-fork policy forbids patching the upstream `ExplicitDifference` in place (the
reason "harden ExplicitDifference" was rejected), it lands as a new fork class.

For the coupled/implicit-damped case, the answer is documentation, not code:
`integrator NewmarkExplicit 0.5`.

## Decisions (settled at review — 2026-05-30, post-sweep)

| # | Decision | Rationale |
|---|----------|-----------|
| C1 | **New clean fork class**, `CentralDifferenceLadruno`, classTag 33003 (Ladruno private band ≥33000); all upstream classes frozen | Sibling-fork policy (no upstream diffs); the reason this is new code, not an in-place patch |
| C2 | **Single explicit leap-frog scheme** — no `-damping` mode switch | The sweep showed a coupled mode == `NewmarkExplicit(0.5)`; bundling two schemes under one classTag was rejected as a split-personality tool. One scheme, one behavior |
| C3 | **Starter on the FIRST STEP, not `domainChanged()`** (B1 fix): `a₀ = M⁻¹(P₀−Cv₀−Fᵢₙₜ(u₀))` then `v₋½ = v₀ − ½Δt·a₀` | `domainChanged()` has no Δt set and no factorized SOE (this is exactly why legacy CD punts). `ExplicitBathe` precedent: defer the first acceleration solve to the first `update()` behind a `firstStep` flag. `dt_cr` *can* stay in `domainChanged()` (it does its own LAPACK eigensolve, not the global SOE) |
| C4 | **Reuse `CriticalTimeStep`**, stability factor 1.0; **implement `getCriticalTimeStep()` override** | Base `TransientIntegrator::getCriticalTimeStep()` returns `-1.0` (verified `TransientIntegrator.h:71`); the `criticalTimeStep()` command (`OPS_criticalTimeStep`, wired in OpenSeesCommands/PythonWrapper/TclWrapper) already dispatches to it. CD limit is exactly `2/ω_max` |
| C5 | **Two velocities**: node/recorder = full-step `vₙ` (via `setVel`); `getVel()` = half-step `v₋½` (modal-damping hook, lagged) | `addModalDampingForce` (IncrementalIntegrator.cpp:525→556) assembles into the residual via `setB`, so `getVel()` must return a solve-time-known (lagged) velocity. `EnergyBalanceRecorder` `DW` integrates against that same lagged velocity for closure |
| C6 | **Diagonal/lumped mass required** (error otherwise); **βK = warn-and-proceed** + report `damped_dt` | Explicit's whole point is trivial `M⁻¹`. βK Rayleigh (`ξ=βω/2` grows with frequency) collapses `dt_cr` ~`2/(βω²)` (T4); warn, auto-pick `damped_dt`, let `-cflAbort` protect. `αM` is safe |
| C7 | **Coupled/implicit-damped CD is OUT** — document `NewmarkExplicit 0.5` | It already exists and is robust; rebuilding it was the redundancy the sweep caught |

## Where

(Footprint verified against the actual `ExplicitBatheLNVD` wiring — 12 files.)

- **New code**: `SRC/analysis/integrator/CentralDifferenceLadruno.{h,cpp}`
- **Modify (registration)**:
  - `SRC/classTags.h` — `#define INTEGRATOR_TAGS_CentralDifferenceLadruno 33003`
  - `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — `#include` + `case`
  - `SRC/runtime/runtime/TclPackageClassBroker.cpp` — `#include` + `case`  *(corrected path — NOT `SRC/tcl/`)*
  - `SRC/interpreter/OpenSeesCommands.cpp` — string dispatch → `OPS_CentralDifferenceLadruno()`
  - `SRC/interpreter/OpenSeesCommands.h` — forward decl `void* OPS_CentralDifferenceLadruno();`  *(was missed)*
  - `SRC/tcl/commands.cpp` — `extern` decl + Tcl `integrator` dispatch branch  *(was missed — without it the integrator is unreachable from Tcl)*
  - `SRC/analysis/integrator/CMakeLists.txt` + `Makefile`
  - *(PythonWrapper.cpp NOT needed — Python `integrator` dispatch is string-based)*
- **Reference (copy patterns from)**:
  - `ExplicitDifference.cpp` — the leap-frog M-only scheme, lagged-Rayleigh residual, `getVel()` returning the half-step velocity. **The closest sibling; this class is essentially it + correct starter + dt_cr + clean velocity output + energy discipline.**
  - `ExplicitBathe.{h,cpp}` — `dt_cr` wiring, `getCriticalTimeStep()`, the `firstStep`-gated first solve, options parsing, broker, the explicit recipe.
  - `CentralDifference.cpp:284` — the single-`update()`-per-step guard (CD allows exactly **one** solve/step; use this guard, *not* ExplicitDifference's looser `>2`).
- **Build**: no new target/dep. Full installer link still blocked by the Ladruno
  link error ([[04_explicit_dynamics_and_energy_balance]]); per-TU `cl.exe` compile-verify
  is the interim gate.

## How

### Algorithm (explicit leap-frog)

Primary state carried: `uₙ`, `v₋½` (half-step velocity), `aₙ`; plus `uₙ₋₁` only for the
full-step output velocity. Coefficient `c3 = 1` (solve for acceleration; `M a = R`).

```
aₙ      = M⁻¹ ( Pₙ − C·v_{n−½} − Fᵢₙₜ(uₙ) )      // damping lagged at known v_{n−½}
v_{n+½} = v_{n−½} + Δt·aₙ
uₙ₊₁    = uₙ + Δt·v_{n+½}
```

- `formEleTangent`: `zeroTangent(); addMtoTang();` → pure diagonal `M⁻¹` (matches
  `ExplicitDifference.cpp:124-131`). Damping is **not** on the LHS (that would be the
  coupled scheme = `NewmarkExplicit`); it is the residual term `−C·v_{n−½}`, formed
  from the committed Rayleigh factors at the lagged half-step velocity.
- Stability: `Δt_cr = (2/ω_max)(√(1+ξ²)−ξ)` — damping reduces the step (T3/T4).

### Initialization (B1 fix)

`domainChanged()` does **only**: allocate state, seed `uₙ, v₋½(=v₀), aₙ` from committed
disp/vel/accel, and run `CriticalTimeStep::computeCriticalTimeStep(model, useTangent,
lumping)` (cache `dt_cr` so `getCriticalTimeStep()` is valid before the first `analyze`).
It does **not** solve for `a₀` — Δt is unset and the SOE is unfactored there.

On the **first `newStep`/`update`** (Δt set, SOE formed/factored), behind a `firstStep`
flag:
```
a₀  = M⁻¹ ( P₀ − C·v₀ − Fᵢₙₜ(u₀) )       // the first ordinary explicit solve
v₋½ = v₀ − ½Δt·a₀                         // back half-step starter (replaces v₀ seed)
```
Then proceed with the standard step. This mirrors `ExplicitBathe`'s deferred-first-solve
pattern and is the correctness fix the legacy classes never made.

### Public API

```python
ops.integrator('CentralDifferenceLadruno')                  # plain explicit CD
ops.integrator('CentralDifferenceLadruno', '-cflAbort', '-lump', 'diagonal')
dt_cr = ops.criticalTimeStep()                              # existing command (C4)
n = max(1, ceil(dt / (0.9 * dt_cr)))                        # 0.9 = LS-DYNA TSSFAC
ops.analyze(n, dt/n)                                        # driver-level sub-stepping
# coupled / implicit-damped CD instead?  ->  ops.integrator('NewmarkExplicit', 0.5)
```

Options (mirror `ExplicitBathe`): `-cfl`, `-cflAbort`, `-tangent`, `-recompute N`,
`-lump <rowsum|diagonal>`, `-verbose`, `-divergence f`. **No `-damping` flag.**
Required recipe: lumped/diagonal mass, `system Diagonal`, `algorithm Linear`
(exactly one solve/step — guard at `updateCount>1`), `dt < dt_cr`.

### Testing / acceptance (extend `Ladruno_scripts/_verify_explicit.py`)

1. **Order of accuracy ≈ 2.0** (SDOF free vibration, log–log error slope).
2. **Stability boundary (undamped)**: stable for `ωΔt < 2.0`, unstable just above
   (strict `<`; at exactly 2.0 the undamped scheme is only marginally stable).
3. **Damped SDOF decay**: matches `e^{−ξωt}` envelope + damped period; largest stable
   step matches `(2/ω)(√(1+ξ²)−ξ)` across ξ = 0.02, 0.1, 0.5.
4. **vs Newmark/CD cross-check** < a few % in the stable range.
5. **1-D wave speed** exact `= √(E/ρ)`.
6. **Rigid-body momentum** conserved exactly (zero stiffness).
7. **`criticalTimeStep() = ℓ/c`** exact on a 1-element bar.
8. **First-step correctness (the key test)**: nonzero `v₀`/`a₀` SDOF reproduces the
   analytical trajectory from step 1 — the test every legacy class fails.
9. **βK trap + guard (T4)**: `getCriticalTimeStep()` drops ~quadratically as `betaK`
   rises (→ `2/(βω²)`); the C6 warn fires for `betaK≠0`; `alphaM` stays benign.
10. **Energy closure ≈ 1%** via `EnergyBalanceRecorder` (multi-DOF). *Gated by the
    Ladruno link blocker — runs once the full link is unblocked.*

## Theory & cross-code review (hardening + adversarial sweep)

> Reviewed against CD theory (Belytschko *Nonlinear FE*, Ch. 6 — stability of the
> damped central difference; Bathe *FE Procedures*, Ch. 9 — starter & amplification;
> Ibrahimbegovic Ch. 6; Hughes *The FEM*, Ch. 9 — spectral analysis, period error)
> and LS-DYNA (Theory §§21–24; `*CONTROL_TIMESTEP`/`*DAMPING_*`). An adversarial sweep
> (numerical / C++ feasibility / consistency / red-team) then stress-tested the plan;
> its outcomes are folded into the Decisions above and noted per-item here.

**T0 — Adversarial-sweep outcomes (2026-05-30).** *Scope:* coupled mode dropped
(== `NewmarkExplicit(0.5)`, verified: default γ=0.5, `addCtoTang(c2);addMtoTang(c3)`,
no `Utm1`/starter defect) → narrowed to explicit-only (C2/C7). *Feasibility:* the
starter cannot live in `domainChanged()` (no Δt, no factorized SOE) → moved to the
first step (C3). *False alarm:* `getCriticalTimeStep()` *is* a base virtual and
`criticalTimeStep()` *is* a live command — the plan's reuse claim was correct.
*Numerical core:* independently re-derived and **verified sound** (boundary, damped
limit, starter, period error, βK collapse all confirmed).

**T1 — This is the leap-frog explicit CD** (LS-DYNA Theory §24.2, Eqs. 24.7–24.12):
`aₙ = M⁻¹(Pₙ − F − H); v_{n+½}=v_{n−½}+Δt·aₙ; uₙ₊₁=uₙ+Δt·v_{n+½}`. LS-DYNA stores
velocity at half-steps by design — our `getVel()` returns that half-step value. The
*coupled* (C-on-LHS, central-velocity) variant is `NewmarkExplicit(0.5)`, out of scope.

**T2 — Cheap `dt_cr` fallback = LS-DYNA solid formula** (Theory §22.1):
`Δt_e = L_e/√(Q+√(Q²+c²))`, `L_e(hex)=V_e/A_max^face`, dilatational
`c=√(E(1−ν)/((1+ν)(1−2ν)ρ))`; global `Δt = 0.9·minₑΔt_e` (TSSFAC=0.9, matches our
`0.9·dt_cr`). The shipped path uses the `K v=λMv` pencil; the `ℓ_e/c_e` estimate
(dilatational speed) is the deferred cheap alternative.

**T3 — Damping reduces `dt_cr` in this (explicit) scheme.** `Δt_cr=(2/ω_max)(√(1+ξ²)−ξ)`
— LS-DYNA Eq. 24.29 ("damping reduces the critical time step size"); standard result
in Belytschko Ch. 6 (damped central difference) and Hughes Ch. 9. Verified numerically
to machine precision at ξ=0.02/0.1/0.5. Report `damped_dt` whenever damping is active.

**T4 — βK is a trap; αM is safe.** For `C=αM+βK`, `ξ_i = α/2ω_i + βω_i/2`. The βK term
grows with frequency → at ω_max, `(√(1+ξ²)−ξ)→1/(2ξ)`, so `Δt_cr→2/(βω_max²)` —
quadratic collapse. LS-DYNA declares classical βK "impractical" in explicit (Vol I
`*DAMPING_PART_STIFFNESS`). αM (`ξ=α/2ω`) damps low frequencies and is safe. → C6 warn.

**T5 — Lumped mass is doubly preferred.** Diagonal `M⁻¹` is trivial, and lumped *lowers*
ω_max vs consistent (rod: `2c/ℓ` vs `2√3·c/ℓ`) → larger stable step. C6 requires it.

**T6 — Zero algorithmic dissipation (γ=½).** 2nd-order accurate, **no** numerical
damping, slight *period shortening* (≈ −Ω²/24, sign verified — opposite to trapezoidal
elongation). Good for waves/energy conservation; bad for un-damped spurious high modes
→ for controllable HF dissipation, reach for `ExplicitBathe`. State this in the header.

**T7 — Bulk/artificial viscosity (out of scope).** LS-DYNA's von Neumann–Richtmyer +
Landshoff `q` (Theory §21) feeds `Q` in the timestep and further shrinks `Δt_e` in
compression — an element/material concern; flagged so wave/shock users know the pencil
`dt_cr` is optimistic when bulk viscosity is active.

## Risks / open questions

- **Implementation traps from the sweep (must honor when coding):**
  - The explicit-mode state must carry the half-step velocity `v_{n−½}` as primary
    (like `ExplicitDifference`'s `Utdot`), **and** the first `v_{n−½}` must be the
    back-half-step starter `v₀ − ½Δt·a₀`, *not* `v₀`. Mis-seeding it (or advancing it
    with a stale accel) silently loses an order / double-counts damping.
  - Single solve/step: guard `updateCount > 1` (CD), not ExplicitDifference's `>2`.
  - `getVel()` is pure-virtual on the base → must be implemented or it won't compile.
- **`dt_cr` caveats inherited** from `CriticalTimeStep`: ignores constraints/MP;
  row-sum lumping non-conservative for rotational DOFs (use `-lump diagonal`); pure
  nodal-mass models report `NOT_APPLICABLE`.
- **D8 aliasing**: copy element mass before fetching stiffness/damping from the same
  element (shared static `theMatrix`) — handled inside `CriticalTimeStep`; replicate in
  any starter/energy code.
- **Energy test (#10) is link-blocked** by the Ladruno error until the full link
  is restored.
- **Backwards compatibility**: none affected — new class, new tag, no upstream edits.

> [!note] Resolved at review (no longer open)
> Default damping mode, βK warn-vs-refuse, `getVel()` convention, consistent-mass — all
> settled (see Decisions). The biggest resolution was *scope*: coupled mode dropped in
> favor of the existing `NewmarkExplicit(0.5)`.

## References

- **LS-DYNA Theory** (`Dropbox/UANDES EC/Lit Review/Explicit/LS-DYNA/LS-DYNA Theory.pdf`):
  §24.2 leap-frog update (p. 501); §24.3 damped stability Eq. 24.29 (pp. 502–503); §22.1
  solid `dt_cr` (pp. 489–490); §21 bulk viscosity (pp. 483–488). **Vol I R16**:
  `*CONTROL_TIMESTEP` TSSFAC (pp. 1890–1897); `*DAMPING_PART_STIFFNESS` βK warning (p. 1935).
- **FEM texts** (`C:\Users\nmora\Desktop\FEM expert\books\`):
  - **Belytschko et al., *Nonlinear FE*** — Ch. 6: CD algorithm (leap-frog box), undamped
    `Δt_cr=2/ω_max`, **damped `(2/ω)(√(1+ξ²)−ξ)`**, consistent-vs-lumped (`2c/ℓ` vs
    `2√3c/ℓ`), 2nd-order. *NOTE: the local PDF is the 1st ed. (2001); reconcile exact
    §/Eq. numbers against it — the result, not the page number, is what's verified.*
  - **Bathe, *FE Procedures* (1982)** — Ch. 9 / Table 9.1: starter `⁻¹U` (`a₃=Δt²/2`),
    `Δt_cr=Tₙ/π=2/ω_max` (Ex. 9.12); damping doesn't change stability character for small ξ
    (p. 540). *Scanned image PDF — render pages to read.*
  - **Ibrahimbegovic (2009)** — Ch. 6: CD ≡ Newmark(β=0,γ=½); lumped `h≤ℓ/(√2c)` vs
    consistent `h≤ℓ/(√6c)`.
  - **Hughes, *The FEM: Linear Static & Dynamic* (1987), Ch. 9** — amplification-matrix
    spectral analysis, damped stability, **period elongation/amplitude decay** (backs T6).
- **In-repo**: `NewmarkExplicit.cpp` (the coupled-CD equivalent we point users to);
  `ExplicitDifference.cpp` (closest leap-frog sibling); `CriticalTimeStep.{h,cpp}` +
  `EnergyBalanceRecorder` ([[04_explicit_dynamics_and_energy_balance]]).

## Implementation log

**2026-05-30 — implemented as designed (explicit leap-frog only).** New sibling-fork
class `CentralDifferenceLadruno` (classTag **33003**, Ladruno private band ≥33000),
no upstream class touched.

- **New files**: `SRC/analysis/integrator/CentralDifferenceLadruno.{h,cpp}`.
- **Registration (8 sites, verified against `ExplicitBatheLNVD`)**: `classTags.h` (33003, Ladruno band ≥33000);
  `FEM_ObjectBrokerAllClasses.cpp` (+include +case); `runtime/runtime/TclPackageClassBroker.cpp`
  (+include +case); `interpreter/OpenSeesCommands.{h,cpp}` (fwd-decl + string dispatch);
  `tcl/commands.cpp` (extern + `integrator` branch); integrator `CMakeLists.txt` + `Makefile`.
- **Scheme bookkeeping**: advance-then-solve (records `u_{n+1}` at `t_{n+1}`, like
  `ExplicitDifference`, so N analyze steps → N records, no off-by-one). The clean
  full-step output velocity is `v_{n+1} = v_{n+1/2} + ½ Δt a_{n+1}` (the centered
  `(u_{n+2}-u_n)/2Δt` identity — available right after the solve, no future value
  needed). `getVel()` returns the lagged half-step `Vhalf` used by the current solve.
- **Starter (C3/B1)**: on the first `newStep()` (behind `firstStep`), the integrator
  forms+factors its own tangent and does ONE extra solve at the committed config to get
  `a₀ = M⁻¹(P₀−Cv₀−Fᵢₙₜ(u₀))`, then seeds `v₋½ = v₀ − ½Δt a₀` before the normal advance.
  domainChanged() only allocates/seeds/runs the dt_cr eigensolve (no Δt, no factored SOE).
- **dt_cr (C4)**: computed in `domainChanged()` UNCONDITIONALLY (so `criticalTimeStep()`
  is valid before `analyze`); stability factor 1.0 (no Noh–Bathe 2×). Reports `damped_dt`
  when damping reduces the step, else `undamped_dt`. `-recompute N`/`-tangent` refresh in
  `newStep()`.
- **βK guard (C6)**: data-driven — Domain exposes no Rayleigh getter, so we flag the trap
  when `damped_dt < 0.9·undamped_dt` (βK collapses it at ω_max; αM barely moves it).
- **Single-solve guard**: `updateCount > 1` (CD), not ExplicitDifference's `>2`.
- **Options**: `-cfl -cflAbort -tangent -recompute N -lump rowsum|diagonal -verbose
  -divergence f`. NO `-damping`, NO positional arg. Default lump = `diagonal`.
- **Diagonal-mass requirement**: documented in the header / `Print` (recipe: `system
  Diagonal`); no cheap runtime introspection of mass diagonality exists, so it is a
  documented requirement rather than a hard error (noted for the PR).
- **Build gate**: `CentralDifferenceLadruno.cpp` per-TU compile-verified standalone with
  `cl.exe` (flags/includes lifted from the Ninja build, repointed at the worktree). Full
  link still blocked by the Ladruno link error — interim gate only, as planned.
- **Tests**: `Ladruno_scripts/_verify_explicit.py` extended with the CDL-1..10 battery
  (CDL-8 first-step / CDL-9 βK collapse are the differentiators; CDL-10 energy closure is
  the MPCO-link-gated skip). They need a rebuilt `.pyd` that registers the class to run.

*(Move to `Ladruno_internal/implemented_robust_central_difference.md` once the full link
is unblocked and the runtime battery passes.)*
