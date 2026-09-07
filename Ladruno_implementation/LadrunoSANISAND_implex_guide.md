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
    <-implexFloor implicit|accept|refuse>  <-implexGuard on|off>  <-implexTrialGuard on|off>
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
initializes at the stage flip to `updateMaterialStage 1`, not before it.

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
| `implexDetail` | 6 | `[0]` total error · `[1]` deviatoric leg · `[2]` volumetric leg (`sqrt(3)\|dp\|`) · `[3]` `p_min` clamp fired on the last pass (0/1) · `[4]` clamp fire count, ever · `[5]` `f`, frozen for this step (reads `0` on a guarded step, §11) |
| `implexRefusals` | 4 | `[0]` total refusals · `[1]` D2 sign-change · `[2]` `-implexControl` past tolerance · `[3]` companion hit `-maxSubsteps` |
| `implexGuards` | 5 | `[0]` floor fallbacks (P2-1, `-implexFloor implicit`) · `[1]` guard firings (P2-2, `f = 0` after a reversal/softening commit) · `[2]` holds preserved (P2-3, zero-`dt` commits left alone) · `[3]` reversal resets restored (P2-5, `-reversalTol`) · `[4]` trial-time `f = 0` fallbacks (P2-6, `-implexTrialGuard`) |

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
