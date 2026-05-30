---
title: Robust central-difference integrator (CentralDifferenceLadruno)
project: Ladruno
status: ready-to-implement
priority: high
owner: nmora
tags:
  - implementation
  - integrator
  - explicit
  - central-difference
---

# Robust central-difference integrator (`CentralDifferenceLadruno`)

> **Design / ADR (pre-implementation).** Plan for a new, clean explicit
> central-difference integrator that fixes the structural defects shared by the
> five existing CD-family classes, reuses the `dt_cr` machinery built for the
> Noh–Bathe work ([[04_explicit_dynamics_and_energy_balance]]), and supports
> **two selectable viscous-damping treatments** (coupled-LHS and explicit-lagged).
> Sits beside `ExplicitBathe` in the explicit core; the Noh–Bathe scheme stays the
> recommended default, this is the textbook-CD reference/peer it never had.
> **Status 2026-05-30: design approved** — all four review questions settled (see
> Decisions); ready to implement on `ladruno` (via this worktree branch → PR).

## What

A new `TransientIntegrator` subclass `CentralDifferenceLadruno` (classTag **64**)
implementing the standard second-order explicit **central-difference** scheme,
done correctly:

- **Correct first step.** Proper half-step starter `u₋₁ = u₀ − Δt·v₀ + ½Δt²·a₀`
  with `a₀ = M⁻¹(P₀ − C v₀ − Fᵢₙₜ(u₀))`. No more `assuming Ut-1 = Ut` warning.
- **First-class, selectable damping** (the decision below): `-damping coupled`
  (C on the effective LHS — textbook robust CD) or `-damping explicit` (C in the
  residual at the known half-step velocity — leap-frog/LS-DYNA style, keeps a pure
  diagonal solve). Rayleigh via `setRayleighDampingFactors`; modal damping via the
  `getVel()` hook.
- **Lumped-mass discipline.** Detect a non-diagonal mass / non-`Diagonal` system
  and warn (coupled mode still works with a consistent solver; explicit mode
  requires diagonal M and says so).
- **`dt_cr` built in.** Reuse `CriticalTimeStep::computeCriticalTimeStep` verbatim;
  report the CD limit `2/ω_max` with stability factor **1.0** (no Noh–Bathe 2×
  bonus). Same `-cfl / -cflAbort / -tangent / -recompute N / -lump` surface as
  `ExplicitBathe`, and the queryable `getCriticalTimeStep()`.
- **Consistent full-step velocity output** `vₙ = (uₙ₊₁ − uₙ₋₁)/(2Δt)` so recorders
  and the energy balance see a physical velocity (not an offset half-step value).
- **Energy-balance ready.** Works with `EnergyBalanceRecorder` (classTag 26)
  unchanged — the closure residual is the silent-wrong detector.

**Not in scope** (deferred / roadmap-gated): mass scaling (roadmap §5.1, a separate
decorator/integrator); sub-cycling inside the integrator (D5 — stays driver-level);
batch/SoA dispatch (§5.2); touching any of the five existing CD classes (left frozen).

## Why

OpenSees today has **five** partial central-difference implementations and **none
is robust**:

| Class | tag | First step | Damping | dt_cr | Velocity out |
|---|---|---|---|---|---|
| `CentralDifference` | 5 | warns `Ut-1=Ut` (wrong) | Rayleigh flags **stored, never used** in stepping | none | full-step |
| `CentralDifferenceAlternative` | 17 | V from committed vel | none | none | half-step (t−½) |
| `CentralDifferenceNoDamping` | 18 | V from committed vel | none | none | half-step |
| `ExplicitDifference` | 55 | warns `Ut-1=Ut` | Rayleigh via residual (lagged) | none | half-step→corrected |
| `ExplicitDifferenceStatic` | 62 | — | FLAC local non-viscous | none | — |

Shared defects this integrator closes: (1) **wrong initial conditions** for any
nonzero `v₀`/`a₀`; (2) damping that is either ignored, silently dead, or only
lagged-explicit with no coupled option; (3) **no critical-timestep guard** — you
learn `Δt` was too big only when it blows up; (4) inconsistent half-step velocity
output that corrupts post-processing. The Noh–Bathe work already built the cure for
(3) (`CriticalTimeStep`) and the detector for instability (`EnergyBalanceRecorder`);
this integrator is the standard-CD peer that should have had them from the start —
useful as the canonical verification baseline (CD is the textbook reference every
explicit scheme is measured against) and for users who want plain CD with discipline.

## Decisions (settled at review — 2026-05-30)

> Four review questions resolved: **default = `coupled`**; **βK = warn-and-proceed**;
> **`coupled` accepts consistent mass (slow-path warning)**; **`getVel()` is the
> modal-damping hook → modal damping is lagged in both modes** (code-confirmed). See
> the Risks section for the resolved-question detail. Numbering note: C1–C7 are grouped
> by topic (C2/C7 are the damping pair), not strictly sequential.

| # | Decision | Rationale | Consequence |
|---|----------|-----------|-------------|
| C1 | **New clean class** `CentralDifferenceLadruno`, classTag 64; leave the 5 legacy CD classes untouched | Sibling-fork policy (same as MPCO_Ladruno); zero regression risk to existing models; smallest review surface | New entry in `classTags.h`, both brokers, command parsing. No diff to upstream-named files |
| C2 | **Selectable damping**: `-damping coupled` (default) vs `-damping explicit` | Two textbook-distinct schemes (T1): **coupled** = C on the effective LHS (3-level u-form) → stability *independent of damping* (T3, Belytschko §6.6.7) but a real solve if C has K-structure; **explicit** = C lagged at the half-step → always diagonal/trivial but damping *reduces* `dt_cr` by `(√(1+ξ²)−ξ)` (T3, Eq. 6.6.43). Pick per the model | One enum + branch in `formEleTangent` (coupled adds C; explicit moves C·v to RHS) and the residual assembly; mode-specific `dt_cr` reporting (C7) |
| C7 | **Mode-aware `dt_cr` + βK guard (from review, sharpened by T3).** In `explicit` mode report `CTSResult.damped_dt`; in `coupled` mode report the undamped `2/ω_max` (stability is damping-independent there). **Hard-warn when stiffness-proportional Rayleigh (`betaK`/`betaKi`/`betaKc`≠0) is used**, with mode-specific text: in `explicit` it *collapses* `dt_cr` (~`2/(βω²)`, T4); in `coupled` it *destroys diagonality* of `M_eff` (real solve, no longer truly explicit, T3). Recommend mass-proportional `alphaM` or modal damping | `αM` is safe in both modes: it keeps C ∝ diagonal M (coupled stays diagonal) and `ξ=α/2ω` *decreases* with frequency. `βK` is the trap LS-DYNA forbids in explicit (Vol I p.1935) | `CriticalTimeStep` already computes both `damped_dt`/`undamped_dt`; select by mode + add the warning. Document `αM`-safe / `βK`-dangerous in the help string |
| C3 | **Correct half-step starter** in `domainChanged()`/first `newStep()` | The defining defect of the legacy CD classes; required for any nonzero IC | One extra mass solve at t=0 to get `a₀`; needs `P₀`, `Fᵢₙₜ(u₀)`, `C v₀` |
| C4 | **Reuse `CriticalTimeStep` as-is**, stability factor **1.0** | Already hardened (DSYGV→DGGEV, lumping, MPI reduce, D8 aliasing-safe); CD limit is exactly `2/ω_max` | `#include <CriticalTimeStep.h>`; no `EB_NB_STABILITY_FACTOR` |
| C5 | **Two velocities, kept separate** (resolved by code review). (a) **Node/recorder velocity** = full-step `vₙ=(uₙ₊₁−uₙ₋₁)/2Δt`, pushed via `theModel->setVel()`/`setResponse()` — fixes the legacy half-step output bug. (b) **`getVel()`** = the *modal-damping* hook (`IncrementalIntegrator::addModalDampingForce`, line 525 → residual via `setB`), so it returns the velocity known at solve time (lagged/half-step) | Modal damping is thus applied **explicitly/lagged in BOTH C-modes** — intrinsic, matches Bathe-ADR D2 | Store `uₙ₋₁`; set node vel to `vₙ`; `getVel()` returns lagged vel; `EnergyBalanceRecorder` `DW` must integrate against the *same* lagged velocity |
| C6 | **Lumped-mass guard**: `-damping explicit` **requires** diagonal M (error if not). `-damping coupled` **accepts** a consistent (non-diagonal) M — does a real factor-once-per-Δt `M_eff` solve — but **warns it is not the fast path** (resolved at review) | Explicit mode's whole point is the trivial `M⁻¹`; coupled mode is mathematically sound with consistent mass and that capability distinguishes this class from the diagonal-only siblings | Diagonal-dominance / `system` check at `domainChanged()`: hard-error (explicit) vs warn (coupled) |

## Where

- **New code**:
  - `SRC/analysis/integrator/CentralDifferenceLadruno.{h,cpp}`
- **Modify**:
  - `SRC/classTags.h` — add `#define INTEGRATOR_TAGS_CentralDifferenceLadruno 64`
  - `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — include + `case`
  - `SRC/tcl/TclPackageClassBroker.cpp` — same (parallel/restart)
  - `SRC/interpreter/OpenSeesCommands.cpp` — string match `"CentralDifferenceLadruno"` → `OPS_CentralDifferenceLadruno()`
  - Dir `CMakeLists.txt` + `Makefile` in `SRC/analysis/integrator/`
- **Reference (copy patterns from)**:
  - `SRC/analysis/integrator/ExplicitBathe.cpp` — `dt_cr` wiring, options parsing, `domainChanged()` state alloc, `getCriticalTimeStep()`, broker, the explicit recipe
  - `SRC/analysis/integrator/CentralDifference.cpp` — the coupled-LHS `addMtoTang(c3)+addCtoTang(c2)` pattern (mode `coupled`)
  - `SRC/analysis/integrator/ExplicitDifference.cpp` — the lagged-residual / M-only pattern + `getVel()` modal-damping hook (mode `explicit`)
- **Build**: no new target or external dep. Full installer link still **blocked by
  the [[04_explicit_dynamics_and_energy_balance]] / MPCO_Ladruno link error** —
  per-TU `cl.exe` compile-verify is the interim gate, same as the Bathe work.

## How

### Algorithm

State stored: `uₙ₋₁, uₙ, vₙ, aₙ` (+ work vectors). Coefficients `c2 = 1/(2Δt)`,
`c3 = 1/Δt²`.

**Central-difference kinematics** (3-level, full step):

```
vₙ = (uₙ₊₁ − uₙ₋₁) / (2Δt)         (= c2 (uₙ₊₁ − uₙ₋₁))
aₙ = (uₙ₊₁ − 2uₙ + uₙ₋₁) / Δt²      (= c3 (uₙ₊₁ − 2uₙ + uₙ₋₁))
```

EOM at step n:  `M aₙ + C vₙ + Fᵢₙₜ(uₙ) = Pₙ`.

**Mode `coupled` (default — proper LHS treatment).** Substitute both relations and
collect `uₙ₊₁`:

```
(M·c3 + C·c2) uₙ₊₁ = Pₙ − Fᵢₙₜ(uₙ) + M·c3·(2uₙ − uₙ₋₁) + C·c2·uₙ₋₁
└──── M_eff ────┘
```

- `formEleTangent`: `zeroTangent(); addMtoTang(c3); addCtoTang(c2);` — exactly the
  standard `CentralDifference` effective matrix, but now with a correct starter and
  `dt_cr`. With lumped M **and** diagonal/Rayleigh-mass C, `M_eff` is diagonal →
  trivial solve. With a consistent C it stays a real (cheap, constant-per-Δt) solve.
- Damping is treated **implicitly** (C straddles step n via the central velocity, on
  the LHS), so per Belytschko §6.6.7 (p. 403) the **stability limit is independent of
  damping** — `dt_cr` stays the undamped `2/ω_max` (T3). The price: with C on the LHS
  the scheme is "not truly explicit for a damped system" — `M_eff` is diagonal **only**
  when C is diagonal-compatible (lumped M + mass-proportional `αM`). With `βK` Rayleigh,
  C carries stiffness structure → `M_eff` is non-diagonal → a real (factor-once-per-Δt)
  solve, losing the explicit advantage. So: choose `αM` to stay fast *and* unconditional.

**Mode `explicit` (leap-frog / LS-DYNA style).** Keep `M` alone on the LHS, move the
damping force to the residual using the known half-step velocity (lagged):

```
M·c3 · uₙ₊₁ = Pₙ − Fᵢₙₜ(uₙ) − C v_{n−½} + M·c3·(2uₙ − uₙ₋₁)
```

- `formEleTangent`: `zeroTangent(); addMtoTang(c3);` → pure diagonal `M⁻¹` (matches
  `ExplicitDifference`/`CentralDifferenceNoDamping`).
- Damping is explicit/conditionally-stable → **reduces `dt_cr`**; documented, and the
  `-cfl` report should annotate the reduction. Requires diagonal M (C6).

### Initialization (the robustness fix, C3)

In `domainChanged()` (after state alloc, seeding `uₙ, vₙ, aₙ` from committed
disp/vel/accel) compute the starter:

```
a₀ = M⁻¹ ( P₀ − C v₀ − Fᵢₙₜ(u₀) )          // one mass solve at t=0
u₋₁ = u₀ − Δt v₀ + ½ Δt² a₀
```

Then run `CriticalTimeStep::computeCriticalTimeStep(model, useTangent, lumping)` and
cache `dt_cr` so `getCriticalTimeStep()` is valid **before** the first `analyze`
(same contract as `ExplicitBathe`, ADR item D6).

### Public API

```python
# default: coupled damping, no dt_cr report
ops.integrator('CentralDifferenceLadruno')

# explicit-lagged damping (pure diagonal solve), with CFL guard
ops.integrator('CentralDifferenceLadruno', '-damping', 'explicit',
               '-cflAbort', '-lump', 'diagonal')

# query the stable step (valid after domainChanged, i.e. after model build)
dt_cr = ops.criticalTimeStep()          # reuses the existing command (D6)
n = max(1, ceil(dt / (0.9 * dt_cr)))
ops.analyze(n, dt/n)                     # driver-level adaptive sub-stepping (D5)
```

Options (mirror `ExplicitBathe` where they overlap): `-damping <coupled|explicit>`,
`-cfl`, `-cflAbort`, `-tangent`, `-recompute N`, `-lump <rowsum|diagonal>`,
`-verbose`, `-divergence f`.

Required recipe (explicit mode): lumped/element mass, `system Diagonal`,
`algorithm Linear` (exactly one solve/step for CD), `dt < dt_cr`.

### Testing / acceptance (extend `Ladruno_scripts/_verify_explicit.py`)

Analytical references — CD must match these to pass:

1. **Order of accuracy ≈ 2.0** (log–log slope of error vs Δt, SDOF free vibration).
2. **Stability boundary**: stable for `ωΔt < 2.0`, unstable just above (vs Noh–Bathe ≈ 3.0).
3. **Damped SDOF decay** matches the closed-form `e^{−ξωt}` envelope and damped
   period — run in **both** damping modes; coupled should hold accuracy at larger
   Δt than explicit, explicit should show the expected `dt_cr` reduction.
4. **Coupled vs explicit vs Newmark** cross-check < a few % in the stable range.
5. **1-D wave speed** exact `= √(E/ρ)`.
6. **Rigid-body momentum** conserved exactly (zero stiffness).
7. **`criticalTimeStep() = ℓ/c`** exact on a 1-element bar.
8. **First-step correctness**: nonzero `v₀`/`a₀` SDOF reproduces the analytical
   trajectory from step 1 (the test the legacy classes fail).
9. **Energy closure ≈ 1%** via `EnergyBalanceRecorder` over a multi-DOF run.
10. **Mode-dependent damped stability (T3)**: on a damped SDOF, sweep Δt to find the
    stability boundary. `explicit` mode must match `(2/ω)(√(1+ξ²)−ξ)`; `coupled` mode
    must stay at `≈2/ω` *independent of ξ* (verify across ξ = 0.02, 0.1, 0.5).
11. **βK trap, both modes (T4)**: in `explicit` mode confirm the stable step drops
    ~quadratically as `betaK` rises (tracks `2/(βω²)` at large βK); in `coupled` mode
    confirm stability holds but `M_eff` becomes non-diagonal (solve required). Confirm
    the C7 hard-warn fires for `betaK≠0`, and that `alphaM` is benign in both modes.
12. **Mode equivalence (T1)**: with C=0, `coupled` and `explicit` produce bitwise-
    close trajectories (they are the same scheme); they diverge only once damping ≠ 0.

## Theory & cross-code review (hardening pass)

> Reviewed against first-principles CD theory (Newmark β=0, γ=½ spectral analysis;
> Belytschko §6, Ibrahimbegovic §6, Hughes Ch. 9) and LS-DYNA's explicit solver
> (Theory §§21–24, Vol I `*CONTROL_TIMESTEP` / `*DAMPING_*`). The two findings that
> changed the plan are **T3** (damping always shrinks `dt_cr`) and **T4** (βK is a
> trap). Citations are PDF pages.

**T1 — CD is Newmark(β=0, γ=½); the two modes are two algebraically-distinct schemes that coincide only when C=0.**
The "coupled" mode is the **3-level displacement form**

```
vₙ = (uₙ₊₁−uₙ₋₁)/2Δt ,  aₙ = (uₙ₊₁−2uₙ+uₙ₋₁)/Δt²   (C straddles step n)
```

The "explicit" mode is the **leap-frog half-step form** — exactly LS-DYNA's native
scheme (Theory §24.2, Eqs. 24.7–24.12, PDF p. 501):

```
aₙ = M⁻¹(Pₙ − C v_{n−½} − Fᵢₙₜ(uₙ));  v_{n+½} = v_{n−½} + Δt aₙ;  uₙ₊₁ = uₙ + Δt v_{n+½}
```

With **C = 0 both are identical**. With C ≠ 0 they differ: coupled uses the *central*
velocity (C on LHS); leap-frog uses the *trailing half-step* velocity (C on RHS,
fully explicit). This is precisely the `coupled`/`explicit` switch (C2). LS-DYNA
stores velocity at half-steps by design — confirming our `explicit` mode is the
"native" explicit scheme and motivating the C5 velocity-reconstruction care.

**T2 — `dt_cr` cheap estimate (the deferred `ℓ_e/c_e`) is the LS-DYNA solid formula.**
LS-DYNA Theory §22.1 (Eqs. 22.1–22.8, PDF pp. 489–490):

```
Δt_e = L_e / √(Q + √(Q²+c²)) ,   global Δt = TSSFAC · minₑ Δt_e ,   TSSFAC = 0.9 default
L_e(hex) = V_e / A_max^face       (volume / largest-face area)
c = √( E(1−ν) / ((1+ν)(1−2ν)ρ) )  (DILATATIONAL wave speed — note the (1−ν)/… factor)
```

Two takeaways: (a) our example's `0.9 * dt_cr` safety factor matches LS-DYNA's
default TSSFAC exactly — keep it; (b) the deferred `ℓ_e/c_e` estimate (Bathe ADR
item 2) should use the **dilatational** speed `√((λ+2μ)/ρ)`, not the shear speed,
and `L_e = V/A_max` for hexes. This is the cheap fallback when the eigensolve is
overkill.

**T3 — damping's effect on `dt_cr` is MODE-DEPENDENT (the sharp, textbook-correct result).**
This is the load-bearing finding. Belytschko (2nd ed.) treats the two cases separately:

- **`explicit` mode (half-step-lagged C — genuinely explicit, LS-DYNA's scheme):**
  damping *reduces* the step. Spectral/`z`-transform analysis gives
  `Δt_cr = (2/ω_max)(√(1+ξ²) − ξ) ≤ 2/ω_max` — **Belytschko §6.6.6, Eq. (6.6.43),
  p. 401** (and Box 6.2, p. 339); identical to **LS-DYNA Theory §24.3, Eq. 24.29,
  pp. 502–503** ("damping reduces the critical time step size"). Belytschko states
  verbatim that *the velocity lag* is what decreases the stable step, "both for
  explicit damping by a linear law such as `C·v` and for any damping in a material law."

- **`coupled` mode (C on the LHS, i.e. damping implicit — the 3-level u-form):**
  the step is **INDEPENDENT of damping** — **Belytschko §6.6.7, p. 403**. The price:
  with C on the LHS the scheme is "not truly explicit for a damped system" — `M_eff`
  needs a real solve unless C is diagonal-compatible. (Bathe agrees qualitatively:
  for small ξ, damping "does not change the overall stability characteristics,"
  §9.4.2, p. 540.)

**This corrects my own first correction.** The original draft's "coupled doesn't erode
`dt_cr`" was actually *right for coupled mode* — what was wrong was generalising it.
The accurate trade-off:

| | `dt_cr` vs damping | solve cost | βK behaviour |
|---|---|---|---|
| `coupled` (C on LHS) | **independent of damping** | diagonal *only if* C diagonal (αM); real factor-once solve if C has K-structure (βK) | safe for `dt_cr`, but βK ⇒ non-diagonal `M_eff` ⇒ loses explicitness |
| `explicit` (lagged C) | **reduced** by `(√(1+ξ²)−ξ)` | always diagonal/trivial | βK ⇒ ξ=βω/2 huge at ω_max ⇒ step **collapses** (T4) |

**Design consequence**: report `CriticalTimeStep.damped_dt` in `explicit` mode (it
already computes the damped pencil); in `coupled` mode the undamped `2/ω_max` is the
honest bound. → this is what C7 now encodes.

**T4 — stiffness-proportional Rayleigh (βK) is a trap in explicit; mass-proportional (αM) is safe.**
For Rayleigh `C = αM + βK`, the modal damping is `ξ_i = α/(2ω_i) + βω_i/2`:

- **αM term**: `ξ = α/2ω` *decreases* with frequency → damps rigid-body/low modes,
  negligible effect on `ω_max` → **safe** (LS-DYNA `*DAMPING_GLOBAL`/`_PART_MASS`,
  p. 1930–1932).
- **βK term**: `ξ = βω/2` *grows* with frequency. At `ω_max`, plug into T3: for large
  ξ, `√(1+ξ²)−ξ → 1/(2ξ)`, so `Δt_cr → 1/(ξω_max) = 2/(βω_max²)` — a **quadratic
  collapse** of the stable step. LS-DYNA therefore declares true classical βK
  *"impractical"* in explicit and substitutes a bounded per-element high-frequency
  damping, requiring a manual TSSFAC cut otherwise (Vol I `*DAMPING_PART_STIFFNESS`,
  p. 1934–1935). → drove the C7 hard-warn.

**T5 — lumped mass is doubly preferred (confirms C6).** Lumped mass both (a) makes
`M`/`M_eff` diagonal → trivial `M⁻¹`, and (b) *lowers* `ω_max` versus consistent mass
(e.g. a 1-D bar: `ω_max^lumped = 2c/L` vs `ω_max^consistent = 2√3·c/L`), so lumped
gives a **larger** stable step *and* a cheaper solve. Consistent mass is legal in
`coupled` mode (real `M_eff` solve) but is the slow path and shrinks `dt_cr` by √3 —
worth the C6 warning.

**T6 — CD has zero algorithmic dissipation (γ=½); this is a feature and a footgun.**
β=0, γ=½ gives second-order accuracy with **no** numerical damping and a slight
*period shortening* (≈ −Ω²/24). Good: energy-conserving, clean for wave propagation.
Bad: spurious high-frequency modes (from stiff/odd elements) never decay and can
pollute the solution — which is exactly why Noh–Bathe (`p`) exists. **So the honest
positioning is: `CentralDifferenceLadruno` = the accurate, dissipation-free CD
baseline; reach for `ExplicitBathe` when you need controllable high-frequency
dissipation.** State this in the class header.

**T7 — bulk/artificial viscosity (out of scope, noted for completeness).** LS-DYNA
adds a von Neumann–Richtmyer + Landshoff viscous pressure `q` (Theory §21, Q1=1.5,
Q2=0.06) that smears shocks and *also* feeds `Q` in the timestep (Eq. 22.3), further
shrinking `Δt_e` in compression. This is an element/material concern, not an
integrator one — flagged only so wave/shock users know `dt_cr` from the pure
`K v = λ M v` pencil is optimistic when bulk viscosity is active.

## Risks / open questions

> [!done] Default damping mode — **RESOLVED: `coupled`** (2026-05-30).
> Default is `coupled`; paired with `αM` it is *both* unconditional-in-damping (T3) *and*
> diagonal — the genuine best case. `explicit` stays available via `-damping explicit`
> for users who want the LS-DYNA-style pure-leapfrog scheme.

> [!done] βK handling (C7) — **RESOLVED: warn-and-proceed** (2026-05-30).
> Hard-warn (mode-specific text), auto-select `damped_dt` in `explicit` mode so
> `-cflAbort` still protects the run; do **not** refuse. A future bounded per-element
> high-frequency damping (LS-DYNA's substitute for βK) is a separate roadmap item.

> [!done] Modal-damping velocity — **RESOLVED by code review** (2026-05-30).
> `getVel()` is the modal-damping hook (`IncrementalIntegrator::addModalDampingForce`,
> line 525): the force is assembled straight into the SOE residual via `setB`, so the
> velocity must be one known at solve time → **lagged/half-step**. Modal damping is
> therefore applied explicitly in *both* C-modes (intrinsic, not a defect; matches
> Bathe-ADR D2). The full-step `vₙ` is the *node/recorder* velocity (set via
> `setVel()`), kept separate — see C5. `EnergyBalanceRecorder` `DW` must integrate
> against the same lagged velocity for closure.

> [!done] Consistent mass in `coupled` mode — **RESOLVED: allow + warn** (2026-05-30).
> `coupled` accepts a consistent (non-diagonal) M and does a real factor-once-per-Δt
> `M_eff` solve, with a "not the fast path" warning. `explicit` still requires diagonal
> M. See C6.

- **`dt_cr` caveats inherited** from `CriticalTimeStep`: ignores constraints/MP;
  row-sum lumping non-conservative for rotational DOFs (use `-lump diagonal` for
  beams/shells); pure nodal-mass models report `NOT_APPLICABLE`.
- **D8 aliasing**: copy the element mass to a local before fetching stiffness/damping
  from the same element (shared static `theMatrix`). Already handled inside
  `CriticalTimeStep`; replicate the pattern in any starter/energy code here.
- **Backwards compatibility**: none affected — new class, new tag, no edits to
  existing integrators.
- **Build**: link still gated by the MPCO_Ladruno error; compile-verify per TU.

## References

- **LS-DYNA Theory Manual** (`Dropbox/UANDES EC/Lit Review/Explicit/LS-DYNA/LS-DYNA Theory.pdf`):
  §24.2 central-difference update / leap-frog (PDF p. 501); §24.3 stability incl. the
  damped bound Eq. 24.29 (pp. 502–503); §22.1 solid-element `dt_cr` Eqs. 22.1–22.8
  (pp. 489–490); §21 artificial bulk viscosity (pp. 483–488); §24.4 nodal/subcycling
  timestep (pp. 504–505).
- **LS-DYNA Keyword Vol I R16** (`…/LS-DYNA_Manual_Vol_I_R16.pdf`): `*CONTROL_TIMESTEP`
  TSSFAC/DT2MS/IMSCL (pp. 1890–1897); `*DAMPING_*` — βK-vs-timestep warning on p. 1935,
  mass-proportional on pp. 1930–1932.
- **CD theory** (textbooks at `C:\Users\nmora\Desktop\FEM expert\books\`):
  - **Belytschko et al., 2nd ed. (2014)** — the primary source. Box 6.1 p. 333 (CD
    algorithm, leap-frog); Eq. (6.2.13) p. 335 (undamped `Δt_cr=2/ω_max`); §6.6.6
    **Eq. (6.6.43) p. 401** (damped `(2/ω)(√(1+ξ²)−ξ)`, the half-step-lagged case);
    **§6.6.7 p. 403** (damping-implicit ⇒ step independent of damping, "not truly
    explicit"); §6.6.9 p. 405 + Table 6.1 (consistent mass shrinks `Δt_cr`: rod
    `ω_max` = `2c/ℓ` lumped vs `2√3·c/ℓ` consistent); §6.2.4 p. 336 (2nd order).
  - **Bathe, *FE Procedures* (1982)** — Table 9.1 **p. 502** (starter `u₋₁`, constants;
    `a₃=Δt²/2`); Eq. (9.13) p. 503 / Ex. 9.12 p. 539 (`Δt_cr=Tₙ/π=2/ω_max`); p. 540
    (damping doesn't change stability character for small ξ). *Note: scanned image PDF,
    no text layer — render pages to read.*
  - **Ibrahimbegovic (2009)** — §6.2.1 p. 405, Table 6.1 (CD ≡ Newmark β=0,γ=½;
    lumped `h≤ℓ/(√2·c)` vs consistent `h≤ℓ/(√6·c)`).
  - **Hughes, *The FEM: Linear Static & Dynamic* (1987), Ch. 9** — the canonical
    time-integration reference; the definitive treatment of the amplification-matrix
    spectral analysis, the damped stability bound, and **period elongation / amplitude
    decay** (backs T6: CD has γ=½ ⇒ zero amplitude decay and slight period error).
    Use to deepen T3/T6 once direct PDF reading is available (post-restart).
- **In-repo**: [[04_explicit_dynamics_and_energy_balance]] (the `CriticalTimeStep` /
  `EnergyBalanceRecorder` infrastructure reused here); `SRC/analysis/integrator/CriticalTimeStep.h`.

## Implementation log

*(empty — fill in once the plan is approved and execution starts; move to
`Ladruno_internal/implemented_robust_central_difference.md` when done)*
