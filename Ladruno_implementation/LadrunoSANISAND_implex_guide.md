---
title: LadrunoSANISAND — IMPL-EX stress update (-implex)
project: Ladruno
status: draft — PR #798 open, not yet merged
priority: high
adr: ADR 92
tags:
  - material
  - nd-material
  - soil
  - sanisand
  - manzari-dafalias
  - implex
  - integration
---

# LadrunoSANISAND — IMPL-EX stress update (`-implex`)

**What it is:** an alternative stress-integration path for `LadrunoSANISAND`, on by a flag. In
the normal (implicit) path, every global Newton iteration re-runs the substepped return-mapping
integrator at the trial strain — and at low confinement that integrator can seize, taking up to
125 state-determination passes and 2000+ seconds in a single `analyze(1)` (ADR 90 GATE U, #791).
Under `-implex` the response handed to the element on every global iteration is instead a
**linear extrapolation** built once from the last committed step:

```
sigma~ = sigma_n + Ce(p_n) : ( d_eps - f * d_eps_p(n) ),   f = (dt_{n+1}/dt_n) * implexAlpha
```

`Ce(p_n)` is the pressure-dependent elastic operator frozen at the committed mean stress —
symmetric, positive definite, constant within the step. Nothing on this path touches `mAlpha`,
`mFabric` or `mAlpha_in`, so a global Newton step against it converges in one iteration and the
assembled tangent is symmetric (unlocks `system Pardiso -matrixType sym`, ADR-75 P1d). The true
substepped return still runs — exactly once, at `commitState` — to obtain the actual `sigma(n+1)`
and advance history; **the committed state is always the implicit one.** IMPL-EX never plasticizes
history on the extrapolated path.

**What it does not do.** IMPL-EX (Oliver, Huespe & Cante 2008) is a **regularizer**, not a cure —
it is a first-order-in-`dt` perturbation with the structure of an artificial viscosity, and *that
is why* it robustifies softening. **The step size is therefore a regularization parameter, and
that parameter has no length in it.** Any width, band or post-peak branch read off an `-implex`
leg is regularized by `dt`, and mesh dependence is not removed. It does not remove the low-`p`
wall either — it relocates the cost: on the fork's own footing-corner deck the controlled arm
traded reach for a bounded refusal count rather than an unbounded ladder (§7 below). Background:
[[92_ladruno_sanisand_implex_adr]].

---

> [!warning] **Push idiom is a precondition for `-implexControl` (measured 2026-09-07, Esmeralda).**
> Drive a prescribed-settlement push as the fork's campaigns do — `LoadControl(-ds)` on an `sp`
> pattern under a `Linear` series with the `Transformation` handler — **not** with
> `DisplacementControl`. Under `DisplacementControl` the load-factor prediction from the frozen
> elastic tangent puts an O(1) trial strain on a near-zero-stiffness free-surface ring whatever
> the step size, so the control refuses at iteration 1 and the leg walls (four legs walled at
> s/B 0.0009–0.002 with error rising 30× while ds shrank 16×). Under `LoadControl(-ds)` the same
> deck, engine and material walk at the full step with the error scaling with ds and the curve
> overlaying the implicit twin to line width. The D2 sign-change guard is not the reason
> (zero firings either way); the pseudo clock is exact on this idiom.
>
> **P2 update (2026-09-07).** `-implexFloor implicit` (the default since P2, §11) closes the
> *other* half of this failure mode — the floor no longer commits an O(1)-error state that the
> next linear solve had to close with a runaway, non-scaling strain, which is what fed the
> self-sustaining loop this warning's own investigation traced. That is a fix to what the
> material **commits** at the floor; it does nothing to what `DisplacementControl` **predicts**
> as a trial strain, which is a property of the integrator, not the floor policy. The idiom
> warning above stands unchanged: still drive a push with `LoadControl(-ds)`, not
> `DisplacementControl`.

## 1. The command

```tcl
nDMaterial LadrunoSANISAND $tag  <23 constants>  \
    -Presidual $pr -Pmin $pmin -honorTolR $h -maxSubsteps $N \
    -implex  <-implexControl $tol $reductionLimit>  <-implexAlpha $a> \
    <-implexDt pseudo|strain|user <$dt>>  \
    <-implexFloor implicit|accept|refuse>  <-implexGuard on|off>  <-implexTrialGuard on|off>  \
    <-flipAlphaIn init|vanilla>  <-implexFlipAbsorb on|off>  \
    <-implexFactor fixed|control>
```

```python
ops.nDMaterial("LadrunoSANISAND", tag, *constants,
               "-maxSubsteps", 20000,
               "-implex", "-implexControl", 0.1, 0.01)
```

`-implex` off (the default) is **byte-identical** to today's `LadrunoSANISAND` — the wrappers'
`setTrialStrain` calls `ladrunoTrialUpdate()`, which with the flag off is exactly
`integrate(); return ladrunoUpdateStatus();`, in that order. So you can add the flags to a model
generator unconditionally and only turn `-implex` on where you mean it.

| flag | meaning | default | notes |
|---|---|---|---|
| `-implex` | turn extrapolation on | off | everything below is refused if given without it |
| `-implexControl $tol $reductionLimit` | refuse a step whose extrapolation error exceeds `$tol` | off; `tol=0.1`, `reductionLimit=0.01` if given bare `-implexControl` values must still be supplied | `$tol > 0`; `$reductionLimit` in `(0, 1]` |
| `-implexAlpha $a` | scales the extrapolated plastic-strain increment | `1.0` | `1.0` = standard IMPL-EX, `0.0` = purely elastic predictor; must be `>= 0` |
| `-implexDt pseudo\|strain\|user <$dt>` | source for `dt_{n+1}` in `f` | `pseudo` | see §5 |
| `-implexFloor implicit\|accept\|refuse` | what a Gauss point commits when `-implexControl` hits the reduction floor with nothing left to cut | `implicit` | ADR-92 P2-1; see §11 |
| `-implexGuard on\|off` | force `f = 0` (elastic predictor) on a step whose committed predecessor showed a loading reversal or `Kp <= 0` | `on` | ADR-92 P2-2; see §11 |
| `-implexTrialGuard on\|off` | on a trial whose `-implexControl` error exceeds `tol` (floor not reached), retry that Gauss point with `f = 0` before refusing | `on` | ADR-92 P2-6; see §11 |
| `-reversalTol $tol` / `-reversalRel $rel` | magnitude guard on the loading-reversal reset (`α_in := α_n`): skip the reset when `‖Δε‖ < max($tol, $rel·‖Δε_lastCommitted‖)` | `tol=1e-10`, `rel=0.05` | ADR-92 P2-5/P2-5b; relative because a hold's per-point strain increment is Newton-tolerance-scale noise (measured median 4e-9, max 1.4e-6) that no fixed absolute threshold clears — see §11 |
| `-flipAlphaIn init\|vanilla` | at the `updateMaterialStage 0 -> 1` flip, leave initialisation to the sign test (`vanilla`, deterministic on a real deck) or force `α_in := α` unconditionally at every point (`init`, a declared modelling variant) | `vanilla` | ADR-92 P2-7; see §11 |
| `-implexFlipAbsorb on\|off` | under `-implex`, whether the flip's first plastic trial also runs a zero-increment companion return to absorb the drift-correction jump (`implexGuards[5]` counts it when `on`) | `off` | ADR-92 P2-7c; opt-in — `on` unconditionally changes the committed state at the flip and fails ADR-92 gate 5 (zero-free-DOF ON/OFF identity); see §11 |
| `-implexFactor fixed\|control\|controlIter` | how `f` is CHOSEN: `fixed` = the clock ratio `alpha*dt_{n+1}/dt_n` (the pre-P2-9 operator, and the only mode gate-passed); `control` = the closed-form minimiser of `\|\|sigma~(f) - sigma_impl\|\|`, computed ONCE at the first trial of the step and frozen — **R3 REFUTED** (biased by the elastic-predictor first iterate); `controlIter` = the same minimiser recomputed at EVERY trial from that iterate's own `d_eps` — **R3 PASSES**, at a wall-time/Newton-churn cost; both control modes keep the clock ratio as the upper bound `f_max` | `fixed` | ADR-92 P2-9; **requires `-implexControl`** (refused without it, not silently downgraded); see §12 |

## 2. What the nine words mean

| flag | reads/writes | scales |
|---|---|---|
| `-implex` | `mImplexOpt.enabled` | switches `ladrunoTrialUpdate()` onto the extrapolated path |
| `-implexControl` | `mImplexOpt.control`, `.errorTol`, `.reductionLimit` | governs the in-step refusal (§7) |
| `-implexAlpha` | `mImplexOpt.alpha` | the extrapolation factor `f`'s scale; not a substep-size knob |
| `-implexDt` | `mImplexOpt.dtSource` (+ `.dtUser`) | what `dt_{n+1}/dt_n` is computed from |

Giving `-implexControl`, `-implexAlpha`, or `-implexDt` **without** `-implex` is refused at parse
time (a flag whose value nothing would read is exactly the "claims to have done something it did
not" defect this fork's parsers exist to make impossible).

## 3. Hard requirements

**`-maxSubsteps` is mandatory** for both companion schemes IMPL-EX qualifies (`ADR 92 D3`):

- **Scheme 1 (`ModifiedEuler`, the deck default).** `-implex` on IntScheme 1 with
  `-maxSubsteps <= 0` (or omitted, whose default is `0` = uncapped) is a **hard parse-time
  refusal**. The companion runs at `commitState`, where no global Newton is left to react if it
  seizes — it must be able to *fail* rather than force-accept at `dT_min = 1e-6`, which is
  precisely what `-maxSubsteps` (ADR-86b / #792 T1) buys.
- **Scheme 2 (`BackwardEuler_CPPM`).** *Permitted* but not the default — P0 measured 58–74 % of
  its calls on the low-confinement corner path taking the low-`p` branch, whose Newton is
  disabled by a literal `errFlag = 0` (`ManzariDafalias.cpp:2264`), so it silently falls through
  to `explicit_integrator` (i.e. `ModifiedEuler` again) and costs a 19-unknown Newton everywhere
  else it doesn't. Same `-maxSubsteps > 0` requirement applies, refused the same way.
- **Every other scheme (0/3/4/5/6/7/8/9/45) is refused with a sentence.** They carry no
  error-controlled substepping, so the companion could not report a failed return and
  `-implexControl` would have nothing to refuse.

**Only `LadrunoBrick` propagates a refusal.** `-implexControl` (and the D2 sign-change guard, and
a companion cap-hit) return the sentinel `LADRUNO_MATERIAL_REFUSED` (`-33086`,
`SRC/material/LadrunoMaterialStatus.h`). Whether that sentinel does anything depends entirely on
the *element*:

- **Propagates it (subdivision engages):** `LadrunoBrick`.
- **Silently accepts it (Newton converges on a refused state, nothing in any log):** `SSPbrick`
  (`SSPbrick.cpp:445`), `Brick` (`Brick.cpp:1069`), `BbarBrick` — and everything else that does
  not specifically check the material's return code.

On a non-propagating element, `-implexControl` still *measures* and *records* the error
(`implexError`), it just cannot cut the step. If your element is not `LadrunoBrick`, read
`implexError` yourself rather than trusting the analysis to stop.

## 4. The stage rule

`-implex` is **inert at stage 0.** `mElastFlag` (a static, flipped for every SANISAND instance at
once by `updateMaterialStage`) gates `integrate()`'s elastic branch; while it is `0` the
extrapolated path is simply unreachable, so gravity and a `LoadControl 0.0` re-equilibration are
bit-identical with the flag on or off. IMPL-EX's own history (`d_eps_p(n)`, `dt` bookkeeping)
initializes at the stage flip to `updateMaterialStage 1`, not before it. The flip handling is per
instance and lazy, so it does not depend on how `updateMaterialStage` is dispatched.

**The first plastic step after the stage flip is exempt from `-implexControl` refusal.**
`d_eps_p(n) = eps_p(n) - eps_p(n-1)` is exactly `0` on that one step (there is no committed
plastic history yet), so `sigma~ = sigma_n + Ce:d_eps` — a pure elastic predictor, with no
extrapolation to be wrong about. The `implexError` measured there is not an extrapolation error;
it is the companion's own drift-correction jump from wherever the elastic stage left the stress to
the first plastic return, and it does not shrink with `d_eps`. Refusing on it means refusing
forever: this was measured directly (`_adr92_p1_bvp_gate_rerun.md`) — before the fix, the
registered arm refused step 1 of every stage at `implexError` 0.13–0.30 against `tol = 0.05`
before any real extrapolation had happened. The error is still computed and reported through
`implexError` / `implexDetail` / `avgImplexError`; only the *refusal* is suppressed, and only on
that one step. Every later step in the stage is primed and refused normally.

## 5. `-implexDt` and the sign-change refusal

`f = (dt_{n+1}/dt_n) * implexAlpha` needs a `dt`, and three sources are offered:

- **`pseudo` (default).** `ops_Dt`, the domain's pseudo-time increment — the `ASDConcrete3D`
  convention. Correct under a settlement-controlled `LoadControl` pattern, including under
  ladder subdivision, because pseudo-time is then proportional to the settlement increment.
- **`strain`.** The contravariant norm of the strain increment.
- **`user $dt`.** A fixed value the deck supplies (also settable at run time via
  `setParameter "implexDt"`).

**Guards, both computed from the frozen-once-per-step `dt`:**

- `dt_{n+1} == 0` (a hold): `f = 0` — no strain advanced, no plastic flow predicted.
- `dt_n == 0` (first step, or first after a hold): falls back to `f = implexAlpha`.
- **A monotone negative clock is legal, and is not refused.** A settlement deck driven by
  `LoadControl(-ds)` has `dt < 0` on every step; two negative increments give the same *positive*
  ratio a monotone-positive clock would. The refusal fires only on a **sign change** between
  consecutive steps' `dt` — a load factor that has turned round (a limit point under
  `DisplacementControl` or arc length), where `dt_{n+1}/dt_n` stops being the extrapolation the
  operator assumes. `DisplacementControl` and arc-length integrators pass a `dt = d(lambda)` that
  is not proportional to the applied increment and can change sign at a limit point — refused,
  with a sentence, rather than silently extrapolating garbage; such a deck should pass
  `-implexDt user` or `-implexDt strain` instead.

(An earlier build gated on `dt > 0.0`, which silently froze `f == 1.0` for the life of any
`LoadControl(-ds)` leg — the exact deck this campaign uses. Fixed; see `LEDGER_quirks.md`.)

## 6. Reading the responses

Five material responses. `implexError` is the per-integration-point value at the last commit;
`avgImplexError` is a process-wide running mean (non-destructive read — every Gauss point a
recorder touches reports the same number). `implexDetail` splits the error and reports the clamp
and `f`; `implexRefusals` is the process-wide refusal ledger — it, not the throttled `opserr`
lines (10 per process, plus one per new subdivision rung), is the only reliable count once a run
generates thousands of refusals. `implexGuards` (ADR-92 P2, §11) is the same kind of ledger for
the three P2 events — none of them prints anything per occurrence (they are designed behaviour,
not warnings), so this response is the only record any of them fired at all.

| response | slots | meaning |
|---|---|---|
| `implexError` | 1 | total error, this material's last commit |
| `avgImplexError` | 1 | process-wide running mean over all commits |
| `implexDetail` | 6 | `[0]` total error · `[1]` deviatoric leg · `[2]` volumetric leg (`sqrt(3)\|dp\|`) · `[3]` `p_min` clamp fired on the last pass (0/1) · `[4]` clamp fire count, ever · `[5]` the `f` **actually used** for the last extrapolation, frozen for this step (reads `0` on a guarded step, §11; under `-implexFactor control\|controlIter` this is `f*`, not the clock ratio — §12) |
| `implexRefusals` | 4 | `[0]` total refusals · `[1]` D2 sign-change · `[2]` `-implexControl` past tolerance · `[3]` companion hit `-maxSubsteps` |
| `implexGuards` | 7 | `[0]` floor fallbacks (P2-1, `-implexFloor implicit`) · `[1]` guard firings (P2-2, `f = 0` after a reversal/softening commit) · `[2]` holds preserved (P2-3, zero-`dt` commits left alone) · `[3]` reversal resets restored (P2-5, `-reversalTol`) · `[4]` trial-time `f = 0` fallbacks (P2-6, `-implexTrialGuard`) · `[5]` hold-skip commits (P2-5c, once per point per hold) · `[6]` control-factor back-offs (P2-9, steps where `f* < 0.5·f_max`) |

Python:

```python
r = ops.eleResponse(eleTag, "material", intPtNum, "implexDetail")
total, dev, vol, clampFired, clampCount, f = r

refusals = ops.eleResponse(eleTag, "material", intPtNum, "implexRefusals")
n_total, n_signchange, n_control, n_companion = refusals
```

Tcl:

```tcl
set r [eleResponse $eleTag material $intPtNum implexDetail]
lassign $r total dev vol clampFired clampCount f

set refusals [eleResponse $eleTag material $intPtNum implexRefusals]
lassign $refusals nTotal nSignChange nControl nCompanion
```

A recorder over the process-wide counter, once per step, is the practical way to track
`implexRefusals` on a long run: `recorder Element -ele $ele -file refusals.out -material $ip
implexRefusals`.


### 6.1 State diagnostics: `psi` and `yieldDistance` (not IMPL-EX specific)

Two more scalar responses, added for the TIMs proposed-model request (2026-09-07, F4). They are
read-only and answer from the **committed** state, so a recorder sees the values that fed the
last committed update.

| response | slots | meaning |
|---|---|---|
| `psi` (alias `stateParameter`) | 1 | the state parameter `psi = e - e_c(p')`, the model's own `GetPSI` with `p' = p + p_residual` floored at 1e-10 — the psi behind `M^b` and `M^d`. With the fork's default `p_r = 0` this is plain `e - e_c(p)` from `state[24]` and the mean stress. |
| `yieldDistance` (alias `yieldFunction`) | 1 | the yield-function value `f = |s - p' alpha| - sqrt(2/3) m p'` on the committed pair: negative inside the cone, `~mTolF` (1e-7 default) on it, never positive after a converged return. |

```python
psi = ops.eleResponse(eleTag, "material", intPtNum, "psi")[0]
f   = ops.eleResponse(eleTag, "material", intPtNum, "yieldDistance")[0]
```

Vanilla `ManzariDafalias` does not answer either name (empty response). Both are inherited by the
3D and plane-strain wrappers. Test: `tests/test_ladruno_sanisand_responses.py`.

## 7. Choosing the tolerance

The registered `-implexControl` operating point (`tol = 0.05`, `reductionLimit = 0.01`) was swept
against three looser tolerances on the same footing-corner deck (`h1.0_e0.6944`, build
`afb95c40c9`), reference `control` (no `-implex`) reaching `s/B = 0.0678` (`WALL` termination):

| `tol` / `reductionLimit` | mode | s/B (depth) | mean overlay dev % (excl. step 1) |
|---|---|---|---|
| 0.05 / 0.01 (registered) | BUDGET | 0.028 | 2.10 |
| 0.05 / 0.1 | BUDGET | 0.028 | 2.10 |
| **0.1** | BUDGET | **0.076** | 1.87 |
| 0.2 | BUDGET | 0.150 | 1.91 |
| 0.5 | TARGET | 0.250 | 2.29 |

**The registered `0.05` fails on reach, not accuracy.** It never gets past `s/B = 0.028` against
`control`'s own `0.0678` — hitting the subdivision *budget*, not a bad extrapolation — while every
tolerance in the sweep, `0.05` included, tracks `control` to a 1.9–2.3 % mean overlay deviation
once the shared step-1 elastic-predictor outlier is excluded. **`0.1` is the tightest tolerance
tested that beats `control`'s own depth while staying under a 5 % mean deviation.**
**Decided (WP-92d):** the C++ default is now `0.1`; use it as the deck
default unless you have a specific reason to run tighter. `reductionLimit` (a floor relative to
the deck's **first** increment) measured **inert at `tol = 0.05`** on this deck — the floor sits
two orders below the working step size at depth and never gets a chance to bind before `tol`
already refuses; loosening it 10x (`0.01 -> 0.1`) with `tol` held at `0.05` produced a
byte-identical run. Do not expect `reductionLimit` alone to buy you reach. Full numbers:
`_adr92_p1_bvp_gate_rerun.md`, "Operating-point sweep" section.

## 8. The reading hazard — read this before quoting a limit point

**An `-implex` curve satisfies equilibrium with the *extrapolated* stress, not the implicit one.**
This is the serious risk in this feature, stated plainly by the ADR (§8): a limit point, plateau,
or capacity read off an `-implex` leg is not evidence of anything on its own. **Confirm it on the
implicit material** (up to the last settlement the implicit solver itself reaches), and **print
`implexError` beside every verdict** — the same discipline the fork's dynamic-relaxation lesson
already forced on this campaign (a solver that finds an exact equilibrium on a wrong path is the
most dangerous instrument the campaign owns). Two arms disagreeing at depths beyond where the
implicit run reaches proves nothing either way; the only honest comparison is over the overlap.

**Reporting condition (decided with the TIMs footing act, 2026-09-07).** Every reported IMPL-EX
curve names its floor-fallback and `f = 0`-guard counts (`implexGuards[0]` and `[1]`, §6/§11)
beside the verdict, and any limit point is confirmed against the implicit twin **over the overlap
only** — the same overlap-only rule stated two paragraphs up, now extended to cover what P2 added:
a curve whose depth outruns its own guard counts, or whose counts are not reported at all, is not
a verdict yet.

## 9. Known limits

- **No plateau measured.** On the fork's own footing-corner deck, no arm — `control`, the
  uncontrolled `-implex` leg, or the registered controlled leg — reaches a plateau on the
  matched-window `t_init` tail (`PLATEAU_FRAC = 2 %`; all three run far above it). `-implex` is
  a solver-robustness device on this evidence, not (yet) a way to see a capacity this deck could
  not otherwise measure.
- **Three mutation survivors are owed tests, not waived** (`_adr92_p1_mutation_gate.md`, score
  0.750 against a 0.60 floor, 9 of 12 hand-mutants killed): **M4** — no deck in the battery arms
  the `-implexControl` reduction floor (`mImplexDt0`), so a subdivision ladder driven by
  `reductionLimit` is untested; **M5** — a refused trial returning a bare `0` instead of
  `LADRUNO_MATERIAL_REFUSED` survives, i.e. the battery pins the refusal's *symptoms* but not its
  return-code *contract*; **M10** — the re-arm-after-refusal line is redundant only because every
  test's failed step goes through `Domain::revertToLastCommit()` (which re-arms anyway), so a
  caller that retries without reverting is untested.
- **Parallel `sendSelf`/`recvSelf` is untested.** The wire grew from 5 to 22 slots to carry the
  flags and `d_eps_p`; the roundtrip test is skipped on every run and, on a zero-free-DOF deck,
  blind by construction even when it isn't (`_adr92_p1_redblue_review.md`). Do not assume a
  parallel (`OpenSeesMP`/`OpenSeesSP`) IMPL-EX run reproduces a serial one until this is measured.
- **`LadrunoSANISANDPlaneStrain` is routed but untested.** The flags reach the 2D wrapper; no
  battery exercises it yet.
- **Cyclic response lags by one step.** `alpha`, `z` and `alpha_in` advance only on the committed
  path, so a reversal is not detected on the extrapolated state. Monotonic pushover is the
  measured target; cyclic use needs its own reversal test before being trusted.

## 10. Verification

`tests/test_ladruno_sanisand_implex.py`. Mutation-gated as part of ADR-87 D2 — PASSED at score
0.750 against the 0.60 floor (`_adr92_p1_mutation_gate.md`); the three survivors are listed in §9
above. BVP-level evidence (not a unit test, a full boundary-value gate on the fork's own
footing-corner deck) is `_adr92_p1_bvp_gate_rerun.md`: the ladder-removal claim confirmed on both
the uncontrolled arm (142/142 converged steps on rung 1, 0 subdivisions) and the registered
controlled arm (504/504 converged steps on rung 1, 0 rung-2/3 — every subdivision attempt was a
material refusal, none a `CTestNormUnbalance` failure). P0's numpy oracle
(`adr92_p0_oracle/sanisand_implex_oracle.py`, `_adr92_p0_oracle_results.md`) is the C++'s
reference on the deck-default paths, matched to `1e-8`; note it has **no `p_min` clamp**, so
parity is meaningful only where the C++ clamp is idle (`LEDGER_quirks.md`).

## 11. ADR-92 P2 — floor policy, the guard, the hold rule, `stressCorrection`

P2 closes the four defects/limits the Esmeralda census and the fork-side probes found on 2026-09-07
(`_adr93_seat_replay.md`, ADR 93 Log 2026-09-06/07; see `92_ladruno_sanisand_implex_adr.md` §"P2
(owed)" for the full evidence table). Shipped in PR #807 (`87b9cf846`).

### Floor policy — `-implexFloor`

At the `-implexControl` reduction floor — error still past `tol`, but `|dt|` already cut below
`reductionLimit * |dt0|`, so there is nothing left to cut — three policies decide what that Gauss
point commits:

- **`implicit` (the default).** The companion return is already computed at this point — it is
  `sigImplicit`, the state the error was just measured against — so the Gauss point delivers it
  (stress, elastic strain, plastic bookkeeping) for that step, under the unchanged frozen `Ce`.
  **Why the default:** the implicit return passes the state at which the control is refusing, i.e.
  it is admissible by construction, so `refuse` would stop IMPL-EX exactly where the implicit
  material itself is still walking forward. `implicit` keeps the committed history internally
  consistent — no O(1) equilibrium-gap commit, and therefore no equilibrium-gap loop for the next
  linear solve to feed on — at the cost of the operator seeing a nonlinear residual at that one
  point for +1–2 Newton iterations, counted (`implexGuards[0]`).
- **`accept`** — the pre-P2 behaviour: commit the extrapolation whatever its error. ADR 93's
  2026-09-07 census measured this feeding a self-sustaining loop (the committed state sits out of
  equilibrium by O(1); the next linear solve closes the gap with a strain that does not scale with
  `ds`, 2–9× per step at the ring; that strain re-triggers the floor; repeat). Kept for reproducing
  a pre-P2 run or isolating the loop itself — diagnostic, not a recommended operating point.
- **`refuse`** — return `LADRUNO_MATERIAL_REFUSED` at the floor too: an honest wall instead of a
  creeping curve. Counted in the `implexRefusals` control bucket, not `implexGuards`.

### The softening / reversal guard — `-implexGuard`

At every commit the material checks its own just-committed state for two conditions and, if
either holds, arms a flag for the **next** step: a loading reversal (`mAlpha_in_n` moved —
`ManzariDafalias::commitState()` reassigns it exactly on the commits that detected one) or
softening (`Kp <= 0`, from the base's own `Kp = (2/3) p h (b:n)`, reproduced verbatim from
`ManzariDafalias.cpp:1373`/`:4954` — the source-true expression, not the cheaper
`(alpha - alpha_in):n` sign proxy, because the proxy misses a softening state reached through
`b:n < 0` past the bounding surface, which is exactly the kind of point this guard exists for).

With `-implexGuard on` (the default), an armed step extrapolates with `f = 0` — a pure elastic
predictor — instead of the previous plastic increment, which is an increment of a branch the
material has already left. The tangent identity is unaffected (`f` is still a constant within the
step, and it is still exactly `Ce`); only the prediction's accuracy is traded, not the step or the
operator's symmetry. A guarded step reads `implexDetail[5] == 0` (that slot is `f`, frozen for the
step). Measured on the seat replay (element 4095, GP 8): error 0.4625 as shipped, 0.029 — under
`tol` — with the guard firing. `-implexGuard` only fires on a *committed* predecessor's state;
`-implexTrialGuard` (P2-6, `on` by default, `implexGuards[4]`) covers the trial that first
*reaches* a softening/reversing point mid-step by retrying that Gauss point with `f = 0` before
the control refuses it.

### The hold rule — P2-3

A zero-strain-increment commit (`analyze` at `dt = 0`, a `LoadControl 0.0` re-equilibration, or
any step `ladrunoImplexTrial()` measures as `mImplexDt == 0`) no longer overwrites
`mImplexDtCommit` or `mImplexDEpsP`. Before P2, storing a zero there made the **next** step's `f`
fall back to `alpha` against a zero history — an extrapolation with the plastic increment silently
switched off on a step nobody asked to change; ADR 93's census measured a curve running
4/21/29 % above the plain leg after a hold. The history and the clock now describe the last step
that actually moved, which is the only step either can honestly describe. A hold is clock-safe
under P2: it increments `implexGuards[2]` and changes nothing else.

The IMPL-EX side of P2-3 was the smaller half. Widening the ADR 93 census (2026-09-06/07) to the
**implicit** column found it was worse: vanilla `ManzariDafalias::integrate()` (`:1005-1013`) resets
`α_in := α_n` on loading reversal with no magnitude guard on the strain increment that decides the
sign, so a hold's round-off noise fires the reset directly — 28–54 % of 34 560 points on Esmeralda
146458 — with no P2-3-style fallback to catch it, sending `h → ∞` and stiffening the implicit column
2.5x for tens of steps. This is P2-5, tracked in `LEDGER_quirks.md` and the ADR 92 P2 table; the fix
is a subclass magnitude guard, `-reversalTol` (default 1e-10 on `‖Δε‖`), counted in `implexGuards[3]`.
Built in #807, pending acceptance.

**P2-5b supersedes the threshold, not the mechanism.** Measured on the fork's R3 footing (1600
GPs, `708152eac`): a hold's per-point strain increment is Newton-tolerance-scale noise (median
4e-9, max 6.4e-8 IMPL-EX / 1.4e-6 implicit), so a fixed `-reversalTol` cannot clear it — at 1e-10
`alpha_in` still reset at 42 % / 9.5 % of points on a hold, and even 1e-7 leaves 2.2 % resetting on
the implicit arm. The guard is now relative to the last committed increment, `‖Δε‖ <
max(reversalTol, reversalRel·‖Δε_lastCommitted‖)`, with `-reversalRel` defaulting to `0.05` — a
hold's increment is `<= 1e-2` of the previous step, a genuine reversal is `~1×`, a halved retry is
`0.5×`, so `0.05` separates a hold from real motion with margin on both sides. The pre-hold
reference is kept across a zero-increment commit so a run of holds does not drift the baseline.
Still building; no holds inside a reported push on either material until the hold acceptance
passes (hold probe `alpha_in` changed = 0 on both arms).

### `alpha_in` at the stage flip is decided by the sign test — deterministic, not noise (`-flipAlphaIn`, P2-7)

Vanilla `ManzariDafalias` never explicitly initialises `α_in` at the `updateMaterialStage 0 -> 1`
flip; it relies on the loading-reversal sign test inside `integrate()` firing on the first plastic
step. That decision is decided by the sign test on the first plastic increment: deterministic on
a real deck, and noise only on an exactly-zero re-equilibration where the fork's synthetic return
uses an exactly zero increment and is neutral. The earlier, wider P2-5/5b/5c reversal-noise guard
had suppressed that reset unconditionally at every state (not just primed ones), which is what
left `α_in = 0` and the implicit path 23-33 % soft from step 1 — a guard-scope defect, not a
defect in the sign test. With the guard confined to PRIMED states, Esmeralda 887fea475 (a real
deck) shows the sign test setting `α_in := α` at 28 629/34 560 points on step 1, identical every
run, and the implicit twin's first-step stiffness returns to the pre-P2 number to the digit
(6.511 / 11.539 / 16.117 / 20.528). There is no defect at the flip left to fix by default, so
`-flipAlphaIn vanilla` (the default) leaves the sign test in control and reproduces real
`ManzariDafalias` exactly; `-flipAlphaIn init` (opt-in) forces `α_in := α` unconditionally at
every point at the flip on both the implicit and IMPL-EX paths — a declared modelling variant,
not a defect fix. Every P2-7 curve names which flag it used.

**The zero-increment companion return at the flip is opt-in, default off (`-implexFlipAbsorb`,
P2-7c).** The first cut of P2-7 had this absorption run unconditionally under `-implex`: at the
flip, the Gauss point's first plastic trial also ran a zero-increment companion return to absorb
the drift-correction jump immediately rather than carry it into the first real step. That
unconditionally changed the *committed* state at the flip, and only when `-implex` was on — which
fails ADR-92 gate 5 (on a zero-free-DOF deck, `-implex` ON and OFF must commit the same state).
Absorbing at the flip is a modelling choice, not a defect fix (the guard-scope fix above is the
defect fix), so it does not ship as default. `-implexFlipAbsorb off` (default) leaves the flip
byte-identical to pre-P2-7 behaviour, and the un-primed first step's committed error (0.24 on the
R3 probe) remains, exempt, as the price of not absorbing. `-implexFlipAbsorb on` runs the
zero-increment companion return as before, counted in `implexGuards[5]`, and brings that error to
~0.05 at the cost of an ON/OFF difference at the flip. **The flip-handled marker itself is now
serialized** (carried through `sendSelf`/`recvSelf` and both `getCopy` forms, per the ADR-86
six-override rule), so a database-restored or MPI instance sees the marker already set and does
not redo the flip's init/absorb work on a redundant `updateMaterialStage` re-assert.

### `stressCorrection` now works — P2-4

`setParameter -val 0 -ele $eleTag stressCorrection` (the idiomatic element-forwarded route) and
the tag-guarded direct route both now reach the flag. Two independent defects made it a silent
no-op before P2: (1) `ManzariDafalias::updateParameter` reads `info.theInt` for this id, but every
interpreter path (`OPS_updateParameter` → `Parameter::update(double)`) writes only
`info.theDouble`, so the flag never actually changed regardless of what the deck asked for; and
(2) the base's own `setParameter` registration requires `argv[1]` to carry the material tag, which
an element-forwarded call never supplies (`Brick::setParameter` hands the material
`argv = {"stressCorrection"}`, `argc = 1`), so the idiomatic route could not even reach id 9 to
begin with. `LadrunoSANISAND::setParameter` now claims the id without the tag guard, and
`LadrunoSANISAND::updateParameter` reads `theDouble` (`ON <=> theDouble != 0.0`; ids `1`
(`updateMaterialStage`) and `5` (`materialState`) are untouched — they already read the field the
interpreter writes). Both routes land on the same base flag now.

## 12. ADR-92 P2-9 — the control-informed extrapolation factor (`-implexFactor`)

**Status: implemented and measured, not shipped as default.** `-implexFactor fixed` is
the default and is byte-identical to every build before P2-9. Two opt-in control modes
exist; both require `-implexControl`. `control` (f* frozen at the first trial of the
step) was run through the plan's Fork R3 registered arm and **REFUTED**: depth 0.052 <
the 0.076 bar and overlay 11.1 % mean deviation, both worse than `fixed` on the same
deck (`_adr92_p2_9_r3_results.md`, Leg 1). `controlIter` (f* recomputed at every Newton
trial from that trial's own `d_eps`) **PASSES** the same bars — depth 0.115, overlay
1.70 % — but costs materially more wall time and Newton churn (see below) and its
Esmeralda dense-refuse arm is still owed before it can be considered for shipping. Both
modes are documented here as measured findings, not as recommended settings; `fixed`
remains the thing to reach for.

### The operator

Under `-implexControl` the companion stress `σ_impl` is already computed at every trial, and the
extrapolated stress is **affine in `f`** on the frozen `Ce`:

```
σ~(f)            = σ_n + Ce:(Δε − f·Δε_p(n))
σ~(f) − σ_impl   = A − f·B,     A = σ_n + Ce:Δε − σ_impl,   B = Ce:Δε_p(n)
f*               = clamp( (A:B) / (B:B), 0, f_max ),   f_max = alpha·dt_{n+1}/dt_n
```

`f_max` is exactly today's `f` — the clock ratio — so `control` never *raises* the factor; it only
decides how much of it to spend. The inner product is `DoubleDot2_2_Contr`, the one `GetNorm_Contr`
and therefore `implexError` itself are built on, so the quantity being minimised is the numerator
of the control's own error measure.

- `f* → 0` exactly where the committed plastic increment points the **wrong way** (the P2-2 case,
  reached with **no model-specific trigger**) or is stale;
- `f* → f_max` where the history is right — today's default behaviour;
- in between it is a graded factor rather than a threshold.

`B:B == 0` (an un-primed history, or a purely elastic one) means there is nothing to choose:
`f_max` stands untouched. That is what keeps the P0 oracle's elastic rows byte-identical.

### Frozen per step (`control`) — and what "the step" means

Under `-implexFactor control`, `f*` is computed **once**, at the first trial of the step, and held
for the rest of it. Later Newton iterates reuse it, so `f` is a constant within the step and the
delivered operator is still exactly `Ce` — the property that removed the subdivision ladder. That
first trial is the elastic predictor, though, and its companion plastic increment is biased small
(the oracle's GD.4 finding, §"R3 verdict" below): freezing `f*` there is what R3 measured as a
REFUTED mode, not a shipped one.

"First trial of a step" reuses the **existing** arm (`mImplexStepArmed`), not a second notion of
freshness: it is the first `setTrialStrain` carrying a non-zero strain increment after a
`commitState()` **or** a `revertToLastCommit()`. A driver that halves its increment and retries
reverts first, so the retry re-arms and recomputes `f*` against its own, new `f_max`. Every refusal
site in `ladrunoImplexTrial()` already re-arms the step, so a refused trial also recomputes.

### Recomputed per iterate (`controlIter`)

`-implexFactor controlIter` is the plan's priced alternative, built and measured rather than left
open: `f*` is recomputed at **every** Newton trial of the step, from that trial's own `d_eps`, not
just the first. `f_max` (the clock ratio, after the P2-2 guard) is still stored once per step in
`mImplexCtlFMax` — the upper bound does not change trial to trial, only the anchor `d_eps` does —
and the back-off census `implexGuards[6]` still fires at most once per step (on the first pass
only), so it stays comparable across modes. Recomputing per iterate spends the step-linearity
property `control` preserved: the delivered *operator* is still `Ce` (the tangent identity is
untouched), but the extrapolated *stress* is no longer affine in `d_eps` within the step, since `f`
itself now varies trial to trial. R3's Leg 2 measured this trade: it removes the first-iterate bias
(overlay 11.1 % → 1.70 %, depth 0.052 → 0.115, both clearing the plan's bars) at the cost of ~13x
the wall time and a jump in Newton max-iteration stalls (1 → 89 over a comparable step count) and in
`n_material_refused`/converged-step (2.31 → 7.39) — see the R3 results doc for the full table.

An abandoned companion never reaches the `f*` arithmetic under either control mode: if
`ladrunoImplexTrial()`'s companion probe hits the substep cap (`mSubstepCapHitInME`), the trial is
refused (`LADRUNO_MATERIAL_REFUSED`, counted in `noteRefusalCompanion()`) and the step is re-armed
**before** the W1b block that builds `A`/`B`/`f*` is ever reached — so `f*` is never built from a
failed companion return, in `control` or `controlIter`.

### Precedence — the P2-2 guard is NOT bypassed

Ordering inside `ladrunoImplexArmStep()` / `ladrunoImplexTrial()`, decided here and recorded so a
reader does not have to infer it:

1. the clock ratio is formed (`f = alpha·dt_{n+1}/dt_n`, `0` on a hold);
2. **the P2-2 guard runs, unconditionally and first.** If it fires it sets `f = 0` and bumps
   `implexGuards[1]`, exactly as before;
3. whatever survives is `f_max`. Under `control` or `controlIter`, `f*` is chosen inside
   `[0, f_max]`.

So **when the guard fires it wins outright** — `f_max = 0` and the clamp can only return `0`. `f*`
replaces the guard's *degree* only where the guard had nothing to say, and it acts one step
**earlier**: the guard reads the committed *predecessor*, so on the reversal step itself it has not
fired yet, while `f*` sees the wrong-way history at the trial. Neither control mode weakens the
no-control fallback, and with `-implexControl` off the guard remains the only mechanism (which is
why `-implexFactor control\|controlIter` is *refused* without `-implexControl` rather than
downgraded).

Everything downstream — the W7 refusal, the `-implexControl` reduction floor, P2-6's trial-time
`f = 0` fallback, and `implexDetail[5]` — reads the **`f` actually used**, so none of that
machinery can be short-circuited by this flag. In particular a P2-6 fallback that fires after `f*`
was chosen still overwrites `f` with `0` and `implexDetail[5]` still reports `0`.

### Reporting

| where | what |
|---|---|
| `implexDetail[5]` | the `f` **actually used** for the last extrapolation — `f*` under `control`/`controlIter`, the clock ratio under `fixed`, `0` if P2-2 or P2-6 acted |
| `implexGuards[6]` | count of steps where the operator "backed off", i.e. `f* < 0.5·f_max`, in EITHER control mode (counted once per step, on the first pass). Not counted when `f_max == 0` (there was no choice to make — P2-2's own slot `[1]` records that) |
| `Print` / construction echo | `-implexFactor = fixed\|control\|controlIter`, beside the other IMPL-EX flags |

### Refusals

| deck says | result |
|---|---|
| `-implexFactor control\|controlIter` with no `-implexControl` | **refused** at construction (and on every `getCopy`/`recvSelf` clone — the check is not gated on `verbose`) |
| `-implexFactor …` with no `-implex` | refused, on the same list as `-implexGuard` / `-implexFlipAbsorb` |
| `-implexFactor <anything else>` | refused with the fixed/control/controlIter explanation |
| the companion probe hits the substep cap (`mSubstepCapHitInME`), under either control mode | the trial is refused (`LADRUNO_MATERIAL_REFUSED`, `noteRefusalCompanion()`) and the step re-armed **before** `f*` is computed — an abandoned companion never seeds `f*` |

### Design forks the plan left open, and how they were settled

- **Recompute `f*` per iterate instead of per step?** Done, as `-implexFactor controlIter` (see
  above) — priced and measured rather than left open. It costs the linearity of the global step
  (the extrapolated stress is no longer affine in `d_eps` within a step), which is the property
  `control` preserved and IMPL-EX was adopted for, but it is a separate flag value, not a change to
  `control`'s own semantics, exactly as this section originally proposed.
- **What if `f_max` is negative?** It cannot be — a sign-changed clock is refused by D2 before this
  point and a zero `dt` gives `f = 0` — but the clamp defensively takes `max(f_max, 0)` so the
  upper bound can never sit below the lower one.
- **Wire format.** `factorMode` crosses `sendSelf`/`recvSelf` in a new slot (`data(32)`, the vector
  widened 32 → 33) on the same rule as the rest of `mImplexOpt`. The per-step arm for the `f*`
  computation is transient and is **not** sent, like `mImplexStepArmed` itself.

### R3 verdict — `control` REFUTED, `controlIter` PASSES, Esmeralda owed

The plan's Fork R3 registered arm (`_adr92_p2_9_r3_results.md`) is the decisive measurement, run on
the same hypoplastic-bearing BVP deck as the P1/P2 `tol0.1` legs:

| | bar | `control` (Leg 1, `a6a53948e`) | `controlIter` (Leg 2, `9e73060a9`) |
|---|---|---|---|
| depth `s/B` | `>= 0.076` | 0.0521 — **fails** | 0.1149 — **passes** |
| overlay mean \|dev\| | `<= 2 %` (PASS) / `> 5 %` (REFUTE) | 11.12 % — **REFUTED** | 1.70 % — **PASSES** |
| `n_material_refused`/converged step | (informational) | 2.31 (down from `fixed`'s 17.97) | 7.39 |
| wall time | (informational) | 141.6 s | 1858.4 s (~13x) |
| Newton max-iter stall markers | (informational) | 1 | 89 |

**Why `control` fails:** `f*` is frozen on the first trial of the step, which is the elastic
predictor — the companion's plastic increment there is biased small (near-zero `A`), so the
closed-form minimiser collapses toward `f* ≈ 0` on exactly the steps where the real history is
non-trivial. This is the same bias the numpy oracle found independently (`GD.4`, see
`_adr92_p2_9_oracle_results.md`): frozen on a bad first iterate, `f*` made the path error up to
100–450x *worse* than today's `f`; recomputed on the converged `d_eps` it was 1.0–2.1x *better*.
The fork measurement and the oracle measurement agree on the mechanism.

**Why `controlIter` passes, and what it costs:** recomputing `f*` at every trial removes the
first-iterate bias (the anchor `d_eps` is no longer the elastic predictor's by the time the step
converges), which is why its overlay and depth both clear the bars — in fact its overlay (1.70 %)
beats every arm measured in this campaign, including plain `fixed`. The cost is real: ~13x the wall
time and two orders of magnitude more Newton non-convergence stalls per comparable step count, plus
more (not fewer) material refusals per converged step than `control`. `controlIter` is therefore the
first `-implexFactor` mode to independently clear both R3 bars, but it is not shipped as default —
the Esmeralda dense-refuse arm (the TIMs-owed arm at build hash with `-implexFactor controlIter`) is
the outstanding gate before any shipping decision.
