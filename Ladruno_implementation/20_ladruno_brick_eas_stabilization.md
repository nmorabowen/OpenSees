---
title: LadrunoBrick EAS — stabilization for inelastic localization
project: Ladruno
status: rejected
priority: low
owner: nmora
tags:
  - implementation
  - element
  - eas
  - stabilization
---

# LadrunoBrick EAS — stabilization for inelastic localization

> [!failure] **REJECTED (2026-06-03).** This ADR scoped, then implemented, then
> *refuted* a scalar `β·Kαα⁰` tangent regularization for `-formulation eas`. The
> DEN-bar gate showed there is no reproducible element-tangent stall to cure (bare
> `eas` traverses with a reasonable solver) and that usable β degrades convergence.
> The code was removed. See the **Implementation log** for the data, the surviving
> **β-independence theorem**, and the **ADR 19 re-diagnosis**. The scoping/analysis
> below is kept verbatim as the record of what was tried and why it seemed plausible.

> Scoping doc. Followed the EAS-vs-ssp/bbar damage comparison in
> [[19_ladruno_brick_eas_simo_rifai]] (PR-2 follow-up), which *reported* the true-EAS
> element (`-formulation eas`) "not robust on notched inelastic problems" — a claim
> the DEN-bar gate here **overturned** (it was a solver/tolerance artifact, not an
> element instability; see the log).

## What

Add an optional **stabilization** to the true Simo–Rifai EAS element so it
traverses notched / localization-dominated **inelastic** problems (where the bare
E9 currently stalls) without sacrificing its smooth-problem accuracy or the patch
test. Default **off** (β=0 ⇒ bit-identical to today's `eas`); opt-in via a knob
(e.g. `-formulation eas -stab <β>`).

**In scope:** small-strain EAS, the stabilization term + its knob + validation
against the ADR 19 DEN-bar comparison. **Not in scope:** finite-strain (enhanced-`F`)
EAS and its separate compressive-hourglassing problem (ADR 19 §3); richer mode
sets E12/E21/E30 (ADR 19 §7).

## Why

The ADR 19 comparison (steel DEN bar, J2 + Lemaitre damage) established, with
evidence, that the bare E9 element:

- is **rank-sufficient elastically** — a free element has exactly 6 zero-energy
  modes and an eigen-spectrum identical to `std`/`bbar` (so this is **not** elastic
  hourglassing);
- handles **homogeneous** plasticity+damage fine (single element → full load,
  damage 0.48, = `bbar`);
- but **stalls just past yield onset on the notch** — a genuine instability of the
  enhanced modes under **non-homogeneous inelastic** tangents. An inner-Newton line
  search only partially mitigated it (reach 0.169→0.231, warnings 89k→28k; not a
  cure) at ~13× cost, so it was reverted.

`bbar`/`ssp` are robust there, so EAS is not *needed* for these problems — but EAS
is the most accurate member of the family on bending/incompressibility, and a
stabilized EAS that is *also* robust on inelastic localization would make it the
default solid element rather than a special-case one. That is the prize.

## Mechanism (the thing the stabilization must fix)

EAS's stability guarantee (Reddy & Simo 1995) is for the **elastic** operator. The
enhanced sub-problem stiffness (live in `formEAStrue` as
`Kaa = Σ_gp Mᵀ dd M · dvol`, [LadrunoBrick.cpp:2680](../SRC/element/ladrunoBrick/LadrunoBrick.cpp))

```
K_αα = ∫ Mᵀ C M dV
```

is positive-definite for the elastic `C`, but as the material tangent `C → C_ep`
(plastic) and especially `C_softening` (damage) degrades in a **high-gradient**
region, `K_αα` can lose positive-definiteness. `K_αα` enters the element at **two**
distinct sites, and they fail in distinct ways:

**(a) Inner-Newton conditioning.** The per-element Newton update
`dalpha = −K_αα⁻¹ residE` ([:2695](../SRC/element/ladrunoBrick/LadrunoBrick.cpp))
loses descent when `K_αα` goes indefinite → the inner solve for `α` flails (the
ADR 19 28k–89k inner warnings). **This is what the reverted line search targeted —
and it "helped but did not cure"** (reach 0.169→0.231). So (a) is *not* the dominant
killer.

**(b) Global-tangent indefiniteness.** The static condensation
([:2734](../SRC/element/ladrunoBrick/LadrunoBrick.cpp))

```
K* = K_dd − K_dα K_αα⁻¹ K_αd
```

subtracts `K_dα K_αα⁻¹ K_αd` from `K_dd`. As `K_αα` shrinks (degrading tangent),
`K_αα⁻¹` *grows*, so the subtracted block grows and can drive `K*` indefinite —
a spurious soft/negative mode of the **global** tangent → the global Newton cannot
converge → stall. The line search left (b) untouched, which is why it failed to
cure. **(b) is the stall driver.**

> [!important] This is **not** a rank deficiency. ADR 19's rank test is explicit: a
> free `eas` element has **exactly 6 zero-energy modes** and a first deformation
> eigenvalue identical to `std`/`bbar` — full rank, **no spurious zero-energy
> (hourglass) mode**. The failure is loss of *positive-definiteness of a full-rank
> tangent under inelastic localization*, not a rank-deficient kinematic mode. This
> distinction drives the choice of cure (§How): a **tangent regularization**, not an
> hourglass-mode stabilization (Reese–Wriggers cures the latter — the wrong family
> here; see §How).

A stabilization must floor `K*` back to soundness (and optionally recondition the
inner solve) **without** perturbing the converged solution on smooth problems —
which, as §How shows, it does *exactly*, because `K_αα` never enters the element
**residual**.

## How — the lever, and why it is a modified-Newton not a physics change

### The lever: regularize `K_αα` *at the solve sites only*

The try-first cure is a constant **Tikhonov floor** on `K_αα`:

```
K_αα_stab = K_αα + β · K_αα⁰ ,    K_αα⁰ = ∫ Mᵀ C₀ M dV          (elastic, constant)
```

`C₀` = initial tangent; `K_αα⁰` is **geometry-only in small strain** (`M` and `dvol`
are cached), so it is a single 9×9 matrix computed **once in `buildEAStrue`**
alongside `easJ0inv`/`easJ0det` ([LadrunoBrick.cpp:2508](../SRC/element/ladrunoBrick/LadrunoBrick.cpp))
— the same "cache the elastic stabilizer at setDomain" pattern as the SSP `Kstab`.
`β` is dimensionless, default **0** (bit-identical to today).

**Critical implementation detail — apply `β·K_αα⁰` at the `.Solve()` sites, never to
the accumulated `Kaa`:** keep `Kaa` as the *true* material Hessian `Σ Mᵀ dd M·dvol`
and form `Kaa + β·Kaa0` only as the linear-solve operator at the two sites that use
it — the inner-Newton update ([:2695](../SRC/element/ladrunoBrick/LadrunoBrick.cpp))
and the condensation `Kaa.Solve(Kad, KaaInvKad)`
([:2733](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)). This keeps the **residual**
`residE = −∫Mᵀσ` and the element internal force `f* = ∫Bᵀσ` **untouched** — which is
the whole reason the converged state is `β`-independent (next).

### Why this is a modified-Newton, not a change to the element physics

`K_αα` is a **Jacobian**, never a **residual**. It enters (i) the inner-Newton
*search direction* and (ii) the element *tangent* `K*` — but the converged element
state is defined by the **residual** condition `h(d,α) = ∫Mᵀσ = 0`, and the global
state by `R(d) = ∫Bᵀσ − f_ext = 0`. Neither residual contains `K_αα`. Therefore:

> [!important] **β-independence theorem (convex regime).** Wherever the enhanced
> sub-problem `h(d,α)=0` has a **unique** root (i.e. `K_αα` stays positive-definite —
> elastic, hardening, and *every accuracy gate*: patch test, bending,
> near-incompressible, reduce-to-`std`), the converged `α(d)`, hence `f*(d)`, hence
> the global root `d*`, are **exactly independent of `β`**. `β` changes only the
> *iteration matrix*, exactly like a modified-Newton / Levenberg–Marquardt damping.
> **The patch test and all smooth-accuracy results are preserved to machine
> precision for any `β` — proven, not "believed safe."** (This retires open
> question #1.)

What `β` *does* do is interpolate the **global tangent** between the two robust
endpoints:

```
β = 0    →  K* = K_dd − K_dα K_αα⁻¹ K_αd        (exact EAS; can go indefinite)
β → ∞    →  K* → K_dd                           (the std displacement tangent; locks but is sound)
```

so `β>0` *floors* the subtracted block `K_dα (K_αα+βK_αα⁰)⁻¹ K_αd`, pulling `K*`
back toward the positive-definite `K_dd` and curing the §Mechanism(b) stall — while
the converged answer stays put.

**The one honest caveat — the localization regime is non-convex.** Under softening
`K_αα` *does* go indefinite; there `h=0` can be multi-valued and the regularized
Newton may converge to a **different (the stable) branch** than the bare element.
That is exactly the behaviour we *want* — but it means in the notch the converged
answer is **mildly `β`-dependent**, so it must be validated against `bbar` and shown
to **plateau** over a `β`-window (the DEN-bar gate, §Testing). Net: `β`-exactness
where we need *accuracy*, `β`-as-branch-selector where we need *robustness*.

> [!warning] Holds only for a **converged** global Newton. Under `algorithm Linear`
> or a loose residual tolerance, an inexact tangent biases the result like any
> modified-Newton (cf. ADR 19 PR-2: `update()` was made non-trivial precisely
> because `algorithm Linear` never re-forms at the final `u`). Document that `-stab`
> assumes an iterating algorithm with a tight test.

### Why **not** Reese–Wriggers / hourglass stabilization

The EAS literature's stabilizers cure a **different** disease:

- **Reese–Wriggers (2000) / Reese (2005)** — stabilization built from the
  **hourglass modes**, scaled by a parameter; cures a *spurious zero-energy mode*.
- **Korelc–Šolinc–Wriggers (2010)** — a *modified enhanced field* `M` that stays
  rank-stable in **finite** deformation.
- **Glaser–Armero (1997)** — stabilized EAS for the **finite compressive** instability.
- **Wriggers–Reese (1996)** — first *documented* the instability (diagnostic).

All target a **rank deficiency / kinematic zero-energy mode** — which ADR 19's rank
test **rules out** for our small-strain element (6 ZEMs, no spurious mode,
§Mechanism). Adding hourglass-mode energy here would inject spurious stiffness
(changing the converged answer) to fix a mode that does not exist. **So the correct
family is tangent regularization (the `β·K_αα⁰` floor above), not hourglass
control.** If the constant scalar `β` proves insufficient, the right escalation is
**within the regularization family**, in order:

1. **State-dependent `β`** — scale `β` by a local loss-of-PD indicator (smallest
   eigenvalue of `K_αα`, or the cached damage `max(dₜ,d_c)` already available via
   `setResponse("damage")`), so the floor is ~0 in the convex bulk and rises only
   where `K_αα` degrades. Keeps β-exactness in the bulk sharper than a global scalar.
2. **A dissipation-/arc-length-controlled global solver** for the softening branch
   (orthogonal to the element; a solver-side fix for the genuinely non-convex
   problem). Korelc-style modified-`M` belongs to the **finite-strain** follow-up
   (ADR 19 §3), not here.

### Knob / API

```
element('LadrunoBrick', tag, n1..n8, matTag, '-formulation', 'eas' [, '-stab', beta])
```

`beta` default 0 (current behaviour, bit-identical). Parsed in
`OPS_LadrunoBrick.cpp` next to the other formulation options; serialized with the
formulation ordinal so DB/parallel streams round-trip. Document a recommended value
once the DEN-bar `β`-sweep finds the plateau.

## Testing

The merge gate is the **ADR 19 DEN-bar comparison, re-run with `-stab`**:

1. **Robustness (the headline):** stabilized `eas` traverses the notched bar to
   the same elongation as `bbar` (currently stalls at u≈0.17/1.5), with the
   converged peak load / dissipated energy **within bbar's mesh-objectivity band**
   (~2–5%, the ADR 19 numbers).
2. **β=0 is bit-identical** to the current element — the `-stab` code path with
   `β=0` must reproduce today's `eas` to the last bit (regression guard; cheapest
   test, run first).
3. **β-independence in the convex regime — to MACHINE PRECISION (the theorem's
   teeth):** for the *accuracy* gates that live in the convex regime — patch test
   (distorted 2×2×2, `∫M dV=0`), bending-beats-`std`, near-incompressible
   (ν=0.45 bending), reduce-to-`std` (α→0) — assert the converged `u`, `α`, and
   reactions are **identical for β ∈ {0, 1e-3, 1e-1, 1}** to ~Newton-tolerance
   (`≤ a few × 1e-8`), **not** merely "within engineering tolerance." This is the
   direct empirical check of the β-independence theorem; if it fails, the
   stabilization is leaking into the residual (a bug — `β·K_αα⁰` was added to the
   accumulated `Kaa` instead of only at the `.Solve()` operator). Reuses the
   existing convex-regime tests in `test_ladrunoBrick_eas.py`, parametrized over β.
4. **β-sweep plateau (the non-convex gate):** on the DEN bar, sweep β and show
   (a) the **minimum β** that clears the stall, and (b) that the converged
   load–displacement response **plateaus** (converges to a β-stable, bbar-matching
   branch) as β grows past that — i.e. a usable window, not a moving target. A
   *cliff* or a monotonic drift with β means scalar β is selecting different
   softening branches → escalate to state-dependent β (§How).
5. **Convergence-rate sanity:** record global iteration counts vs β — confirm β
   trades quadratic convergence for robustness gracefully (counts rise modestly,
   no oscillation), consistent with the modified-Newton framing.

## Risks / open questions

> [!success] **RETIRED — patch-test / smooth-accuracy bias.** The original open
> question ("does `β·K_αα⁰` bias the converged `α`?") is resolved by the
> β-independence theorem (§How): `K_αα` is a Jacobian, never a residual, so the
> converged state in the convex regime — where every accuracy gate lives — is
> **exactly** β-independent. Demoted from "must prove" to a **machine-precision
> regression assertion** (§Testing #3). *Contingent on the implementation applying
> `β·K_αα⁰` only at the `.Solve()` operator, never to the accumulated `Kaa`.*

> [!question] **Open — is scalar β enough on the notch?** The genuinely open
> question is the **non-convex** one: does a single global scalar β clear the stall
> *and* land on the bbar-matching branch with a stable plateau (§Testing #4), or
> does β drift the softening branch (→ state-dependent β, §How escalation 1)?
> Note this is **not** the "Reese–Wriggers deformation-dependent vs scalar" question
> from the original draft — that was the wrong family (hourglass control for a
> non-existent zero-energy mode, §How). Escalation, if needed, stays **inside the
> regularization family**.

- **Over-stiffening is NOT an accuracy risk (correction to the original draft).** A
  large β does *not* "lose EAS accuracy" — by the theorem the converged answer is
  β-independent in the convex regime; it only **slows global convergence** (pushes
  `K*` toward `K_dd`). So the β-sweep window is bounded by *convergence cost*, not
  *accuracy loss*. The only place a large β changes the *answer* is the softening
  branch — captured by gate #4, not by a bending-accuracy gate.
- **Tuning vs automatic:** a hand-tuned global β is a mild footgun *for convergence
  speed* (not accuracy); the principled automatic version is the state-dependent β
  of §How (eigenvalue- or damage-scaled), which keeps β≈0 in the bulk. Decide after
  the scalar DEN-bar sweep.
- **Interaction with `-autoRegularization`:** orthogonal by construction — `β·K_αα⁰`
  regularizes the *element tangent* (a Jacobian), the crack-band `lch` scaling
  reshapes the *material* softening law (a residual). Still assert it in gate #1
  (the DEN bar runs with `-autoRegularization` on).
- **Inner-solve coverage:** applying `β·K_αα⁰` at the inner-Newton site too (free —
  doesn't move the fixed point) reconditions §Mechanism(a) at no cost; verify it
  does not *slow* the inner Newton enough to matter (it shouldn't — the inner loop
  is ≤12 iters and the floor only helps when `K_αα` is near-singular).

## Implementation log

### 2026-06-03 — refinement pass (still scoping, no code)

Sharpened the mechanism and the cure after re-reading the live `formEAStrue`
against ADR 19's rank evidence. Net changes to the plan:

- **Reframed the cure as a tangent regularization / modified-Newton**, not an
  "enhanced sub-problem stabilization." `K_αα` is a Jacobian, never a residual ⇒
  the converged state is **exactly β-independent in the convex regime** (the
  β-independence theorem) ⇒ **open question #1 (patch test) is retired**, demoted to
  a machine-precision regression assertion.
- **Split the failure into (a) inner-Newton conditioning vs (b) global-`K*`
  indefiniteness**; identified **(b) as the stall driver** (the reverted line search
  targeted (a) and "helped but didn't cure" — consistent).
- **Reclassified Reese–Wriggers as the wrong family** (hourglass control cures a
  zero-energy mode ADR 19's rank test rules out); escalation now stays inside the
  regularization family (state-dependent β → dissipation/arc-length solver).
- **Corrected the over-stiffening risk:** large β costs *convergence speed*, not
  *accuracy* (theorem) — the β-window is bounded by cost, and answer-change only
  happens on the softening branch.

**Concrete implementation checklist (when picked up), all in `LadrunoBrick`:**

1. `LadrunoBrick.h` — add `double easStabBeta;` (default 0) and `Matrix easKaa0(9,9);`.
2. `buildEAStrue` ([:2508](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)) — after
   caching `easJ0inv`/`easJ0det`, accumulate `easKaa0 = Σ_gp Mᵀ C₀ M · dvol` using
   `getInitialTangent()` (geometry+elastic, constant in small strain).
3. `formEAStrue` ([:2577](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)) — keep `Kaa`
   as the true Hessian; build a local `KaaStab = Kaa; KaaStab.addMatrix(1.0,
   easKaa0, easStabBeta);` and use **`KaaStab.Solve(...)`** at **both** the inner
   update ([:2695](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)) and the
   condensation ([:2733](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)) — **and the
   `getInitialStiff` condensation** ([:2650](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)).
   **Do not** touch `residE` or the `f* = ∫Bᵀσ` assembly (that is what guarantees
   β-independence). The inner-Newton **convergence test stays on the true `residE`**
   (`r ≤ tolRel·r0 + tolAbs`) — unchanged.
4. `OPS_LadrunoBrick.cpp` — parse `-stab <beta>` (EAS-only; warn/ignore for other
   formulations), store into `easStabBeta`.
5. `sendSelf`/`recvSelf` — ship `easStabBeta` (one `dData` slot) so DB/parallel
   round-trip; `easKaa0` is rebuilt on the receive side in `buildEAStrue` (already
   called in `setDomain`/recv path).
6. Ledger + banner: this is a *new knob on an existing shipped feature*, not a new
   classTag — note it on the `LadrunoBrick eas` ledger row; no banner line needed
   unless we want to advertise `-stab`.

**Start:** the scalar `β·K_αα⁰` experiment on the DEN bar, gated by the §Testing
battery — run gate #2 (β=0 bit-identical) and gate #3 (β-independence to machine
precision) *first* as a correctness firewall, then gate #1/#4 for the actual cure.

### 2026-06-03 — scalar-β implemented, tested against the DEN bar, **REFUTED and REMOVED**

> [!failure] **OUTCOME: scalar `β·Kαα⁰` does not work and was removed.** The lever
> was fully implemented, unit-passed (β=0 bit-identical; β-independence to machine
> precision), then tested on the actual ADR 19 DEN-bar gate — where it **failed**:
> there is no reproducible element-tangent stall to cure, and usable β values
> *degrade* convergence. The code (`-stab`, `easStabBeta`, `easKaa0`, the test) was
> stripped back out. The β-independence theorem and the ADR 19 re-diagnosis below are
> the surviving results. **Status flipped to `rejected`.**

#### What was built (now reverted)

Implemented the try-first scalar lever in `LadrunoBrick` (branch
`guppi/eas-stab-impl`). What landed, vs the checklist above:

- **State:** `easStabBeta` (default 0) + `easKaa0(9,9)` members; `easKaa0 =
  Σ MᵀC₀M·dvol` accumulated once in `buildEAStrue` (getInitialTangent + the cached
  centroid map), rebuilt on the receive side. `-stab β` parsed in
  `OPS_LadrunoBrick.cpp` (eas-only — warned+ignored for other formulations; negative
  β clamped). β shipped in the EAS-guarded extra vector, widened `Vector(9)→(10)`
  (`alphaCommit` + β); non-eas streams still byte-identical.

- **DEVIATION from the checklist — applied at the condensation site ONLY, not the
  inner-Newton, not `getInitialStiff`.** The checklist (and the original §How "free
  at the inner site too") said to use `Kaa+βKaa0` at *both* the inner-Newton update
  and the condensation. Implementing it exposed a flaw in the "free" claim:

  - **Inner Newton (reverted):** the convex inner sub-problem is *linear*, so the
    true-`Kaa` Newton converges in **one** step. Damping the direction with
    `βKaa0` turns that into a geometric iteration with ratio ≈ `β/(1+β)` — for β=1
    that needs ~26 iters, exceeding `maxIters=12`, so α would **not** fully
    converge and gate #3 (β-independence to machine precision) would *fail by
    construction*. So the inner Newton stays on the **true `Kaa`** → α converges
    exactly and fast → α (hence `f*`) is rigorously β-independent. The inner site is
    **not** "free"; it is actively harmful in the convex regime. (Inner-solve
    robustness for mechanism (a) — if ever needed — must come from a method that
    doesn't slow convex convergence, e.g. an eigenvalue-gated floor, not a constant
    `βKaa0`.)
  - **`getInitialStiff` (left pure):** the elastic tangent at α=0 is PD by
    construction (never the indefinite one), so regularizing it buys no robustness
    and would needlessly perturb modal/eigen and initial-stiffness analyses with an
    opt-in localization knob. `-stab` touches only the nonlinear condensed tangent.

  Net: `K* = Kdd − Kda (Kaa+βKaa0)⁻¹ Kad` at the **final condensation only**
  (`formEAStrue` tang path). This is the minimal lever that (i) keeps the converged
  state exactly β-independent in the convex regime — the inner Newton owns α, the
  condensation owns only the *tangent* handed to the global solver — and (ii) floors
  K* toward Kdd to cure mechanism (b), the stall driver. `β·Kaa0` is formed as a
  local `KaaStab` copy at the `.Solve()` so the accumulated `Kaa` is never mutated.

- **β=0 fast path:** every site branches `if (easStabBeta != 0.0)`, so β=0 runs the
  original code verbatim (gate #2 bit-identical).

- **Unit tests (passed):** `tests/test_ladrunoBrick_eas_stab.py` — parser/guard, β=0
  bit-identical, β-independence to machine precision (distorted patch + bending +
  hardening-J2 over β∈{0,1e-3,1e-1,1}), all under `algorithm Newton` + tight residual
  test. These confirmed the *correctness* of the lever (the β-independence theorem),
  not its usefulness.

#### The DEN-bar gate — the refutation

Ran the ADR 19 DEN-bar comparison (`test_lemaitre_notched_bar` geometry; steel J2 +
Lemaitre + `-autoRegularization`) with `-formulation eas -stab β` over a β sweep, two
solvers, and four difficulty configs. The headline tables (reached elongation of 1.5):

```
adaptive solver (step-cut + KrylovNewton/LineSearch), coarse mesh:
  β:        0      1e-3    1e-2    1e-1    1       10
  reached:  1.50   1.50    0.78    1.50    1.50    0.57     ← bbar ref = 1.50
  peak[N]:  23916  23916   23916   23916   23916   23916    ← IDENTICAL (β-independent peak)

plain Newton (fixed step), reached elongation by β:
  config                  β=0    1e-2   1e-1   0.3    1      3
  fine,   reg, small-step 1.50   1.50   1.50   0.23   0.06   0.06
  fine, NOreg, small-step 1.50   1.50   1.46   0.20   0.06   0.06
  coarse,NOreg, BIG-step  0.69   0.24   0.45   1.50   0.15   0.09   ← β=0.3 rescues...
  fine,  NOreg, BIG-step  0.49   0.11   0.11   0.06   0.06   0.06   ← ...but here nothing helps
```

**Findings:**

1. **No reproducible element-tangent stall.** Bare `eas` (β=0) traverses to full
   elongation in *every small-step* config, any mesh, with or without regularization.
   A stall appears only under **big steps** — a Newton basin-of-attraction / solver
   issue, not an enhanced-mode instability. ADR 19's "stall at 0.17, damage≈0 (at
   yield onset)" was the **inner-Newton absolute-tolerance bug** that ADR 19 PR-2
   itself fixed (the relative criterion), plus coarse stepping — **not** the
   "instability of the enhanced modes under non-homogeneous inelasticity" that ADR 19
   §PR-2-follow-up (and this ADR's §Why) inferred. **That inference is now retracted.**

2. **β is erratic and usually harmful.** The β-independence theorem holds beautifully
   (peak load bit-identical across all β — the peak is convex-regime). But in the
   *non-convex* softening regime β selects a branch, non-monotonically: the very same
   β=0.3 that rescues `coarse-BIG-step` (0.69→1.50) **destroys** the two small-step
   configs (1.50→0.23/0.20) and helps nothing in `fine-BIG-step`. There is **no β
   that is universally safe**, and usable magnitudes (β≳0.3) typically convert a
   converging run into a stall.

3. **Why, by theory (the lesson):** Newton needs a tangent *consistent* with the
   residual, **not** a positive-definite one. Regularizing `K* → Kdd` makes it
   inconsistent with the exact (β-independent) `f*`, which *degrades* Newton — exactly
   the observed β≳1 collapse. And `Kαα⁰` is *elastic*-scale, so flooring a genuinely
   degraded `Kαα` needs β~O(1), at which point `β·Kαα⁰` dominates and the tangent is
   fully inconsistent. There is **no β window** between "inert" and "dominant-and-
   inconsistent." The "floor K* toward Kdd to help" intuition is backwards for a
   Newton solver. The real cure for these big-step stalls is solver-side (line search,
   arc-length/dissipation control, adaptive step-cutting) — already in OpenSees.

#### Decision

**Removed.** Scalar `β·Kαα⁰` tangent regularization is not a viable cure for EAS on
inelastic localization. The `-stab` knob, `easStabBeta`/`easKaa0` state, the
serialization widening, and `tests/test_ladrunoBrick_eas_stab.py` were reverted;
`LadrunoBrick` is back to the bare-`eas` state. **Surviving results:** (a) the
**β-independence theorem** (a clean statement of why `Kαα`-side modifications never
perturb the converged answer in the convex regime); (b) the **ADR 19 re-diagnosis**
(the notched-bar stall = the already-fixed inner-Newton tolerance bug + coarse
stepping, **not** an enhanced-mode instability — `eas` is robust on this DEN bar with
a reasonable solver). **Do not re-attempt scalar β.** The only EAS instability that
genuinely needs a stabilization is **finite-strain compressive hourglassing** (a real
zero-energy mode, de Souza Neto §15.2.5) — a *different* mechanism cured by a
*deformation-dependent* Reese–Wriggers term, out of scope for small-strain `eas`.
