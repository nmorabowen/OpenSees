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
    -Presidual $pr -pRe $pre -Pmin $pmin -honorTolR $h -maxSubsteps $N \
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
| `-flipAlphaIn init\|vanilla` | at the `updateMaterialStage 0 -> 1` flip, force `α_in := α` at every point (`init`) or leave initialisation to vanilla's loading-reversal sign test (`vanilla`, which reproduces real `ManzariDafalias` and reads the sign of round-off after elastic holds) | **`init`** (since WP-112; was `vanilla`) | ADR-92 P2-7, WP-112 (TIMs F14); `vanilla` warns once per Gauss point when it meets a round-off `α − α_in`; see §11 |
| `-implexFlipAbsorb on\|off` | under `-implex`, whether the flip's first plastic trial also runs a zero-increment companion return to absorb the drift-correction jump (`implexGuards[5]` counts it when `on`) | `off` | ADR-92 P2-7c; opt-in — `on` unconditionally changes the committed state at the flip and fails ADR-92 gate 5 (zero-free-DOF ON/OFF identity); see §11 |
| `-pRe $p` | **elastic-only** confinement floor (ADR-93 II.1): the three `GetElasticModuli` overloads read `G, K ~ sqrt(max(p + $p, p_min)/P_atm)` and **nothing else in the model changes** | `0` = OFF | not IMPL-EX-specific and not gated on `-implex`; `-Pelastic` is accepted as a synonym. Refused below `0` and refused if given twice (it is a constitutive constant, and the echo can report only one value); **warns** above `0.1*P_atm`; prints a NOTE when `pRe <= p_min`, where the clamp already dominates as `p -> 0`. **Inert at stage 0** — see §3.1 and §4 |
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
- **Scheme 2 (`BackwardEuler_CPPM`).** *Permitted* but not the default. Same `-maxSubsteps > 0`
  requirement applies, refused the same way (`LadrunoSANISAND.cpp:2039-2048`, current line numbers -- see the source, not this note, if they move again) — and, as of WP-108,
  the constructor's own "`-maxSubsteps` has NO EFFECT with IntScheme 2" warning no longer prints,
  because it was **false** (WP-105 / F12; `schemeReachesModifiedEuler()` used to answer `false` for
  `s == 2`, see `LEDGER_quirks.md`). `BackwardEuler_CPPM`'s own retry ladder
  (`ManzariDafalias.cpp` ~2472-2588) falls back to `explicit_integrator` on non-convergence or
  ladder exhaustion, and that switch does not enumerate `INT_BackwardEuler` — it hits `default:` ->
  `ModifiedEuler`, the same seam `-maxSubsteps`/`-honorTolR` read. WP-105 measured what scheme 2
  buys and what it costs at a material point on a GIVEN strain increment (the commit-time companion's
  own regime): **3.7–4.3x less discretisation error and 4.2–7.6x less wall time** than scheme 1 at
  `dEz = 1e-4`, rising to 7–30x / 10–13x at `4.6e-4`; the low-`p` explicit fallback fires on **0 %**
  of steps until the point is pinned at `p_min`, where it fires on ~53 % and both schemes become the
  same operator. What it costs UNDER A GLOBAL NEWTON (an increment merely proposed, not given): a
  CPPM step that cannot return recurses through up to 512 half-increments before falling back,
  **silently** — `integrate()` discards the return value, `debugFlag` is compiled off, and nothing
  refuses. Measured at 12–134 s for a single failing step and at outright stalls (8 of 8 arms) where
  scheme 1 completes (1 of 8); on the CP1/ADR-95 bearing leg it reached `s/B = 4e-5` in 1347 s
  against the baseline's `0.019` in 1267 s, with `ds` pinned at 25x the subdivision floor. **Use it
  where the increment is given (a companion return, a prescribed-strain probe); do not make it the
  primary integrator of a load- or displacement-controlled BVP without a cap on the ladder.** Full
  tables: `Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md`.
- **Every other scheme (0/3/4/5/6/7/8/9/45) is refused with a sentence.** They carry no
  error-controlled substepping, so the companion could not report a failed return and
  `-implexControl` would have nothing to refuse.

**A TRIAL-time refusal only works on an element that forwards `setTrialStrain`'s return code.**
`-implexControl` (and the D2 sign-change guard, and a companion cap-hit caught at the trial)
return the sentinel `LADRUNO_MATERIAL_REFUSED` (`-33086`,
`SRC/material/LadrunoMaterialStatus.h`). Whether that sentinel does anything depends entirely on
the *element*. Audited at source on `9c2f964ea` — this list used to say "only `LadrunoBrick`",
which was wrong in both directions (WP-99 / F7):

**The full audited roster — all 52 NDMaterial-hosting elements, 26 forward / 1 sentinel-only / 25
discard — is the table "Element refusal roster" in `LEDGER_quirks.md`, and that is the only
authoritative copy.** Examples, so this section reads on its own:

- **Propagates ANY nonzero code (subdivision engages):** the fork continuum elements
  `LadrunoBrick20`, `LadrunoQuad`/`CST`/`LST` (`update()` and the EAS path), `LadrunoUP`,
  `BezierTet10`, `BezierTri6` — and, the group that matters most here because u-p is SANISAND's
  canonical host, **every vanilla u-p element that has its own `update()`**: `FourNodeQuadUP`
  (`:419`), `BBarFourNodeQuadUP` (`:369`), `Nine_Four_Node_QuadUP` (`:572`),
  `Twenty_Eight_Node_BrickUP` (`:983`).
- **Propagates ONLY the sentinel:** `LadrunoBrick` — deliberately, per ADR-33/34, so
  `ASDConcrete3D`'s negative "best-state" codes do not fail a step. A material that returns some
  *other* nonzero value is silently swallowed here. It is the only element in this class.
- **Silently accepts it (Newton converges on a refused state, nothing in any log):** `Brick`
  (= `stdBrick`, `Brick.cpp:1069` → `return 0` at `:1073`), `BbarBrick`, `BrickUP`, and the u-p
  elements that drop the code — `BBarBrickUP` (no `update()` override at all),
  `SSPquadUP` and `SSPbrickUP` (an `update()` that calls `setTrialStrain` and returns 0
  regardless) — plus
  `SSPquad`, `SSPbrick`, `FourNodeTetrahedron`, `EnhancedQuad`, `NineNodeMixedQuad`,
  `LadrunoSolidShell`.

On a non-propagating element, `-implexControl` still *measures* and *records* the error
(`implexError`), it just cannot cut the step — read `implexError` yourself rather than trusting
the analysis to stop.

**At COMMIT time no element propagates anything**, because `Domain::commit()` is a bare
`elePtr->commitState();`. That is why a commit-time companion failure LATCHES the material
instead — see §9.

### 3.1 `-pRe` — a stiffness floor, and the three things it is not (ADR-93 II.1)

A cohesionless sand has no stiffness at zero confinement, so at the free-surface ring beside a
footing `G, K -> 0` and the substep controller spends thousands of substeps proving a stress
that carries nothing to relative tolerance. `-pRe $p` puts a declared floor under the ELASTIC
moduli only: `pn = tr(sigma)/3 + pRe` inside the three `GetElasticModuli` overloads, clamped at
`m_Pmin` after the addition, so the effective floor is `sqrt(max(p + pRe, p_min)/P_atm)`.

**It is not `-Presidual`**, which floors the STRENGTH side (`GetF`, `psi`, `M^b`, `M^d`, `D`,
the `D_factor` sigmoid, the low-`p` integrator guards) and never reached the moduli at all —
the asymmetry ADR 93 §1 row 1 measured. **It is not a cohesion**: it adds no strength and does
not move the critical-state line. **It is not `-Pmin`**, which is a clamp on the STRESS; `-pRe`
changes no stress, only the tangent the stress is integrated with. Default `0` is
byte-identical: the parameter enters as `+ 0.0`, and the one derived quantity recomputed for it
(the initial `mCe` at `p = P_atm`, which the base fixes before the fork's last write lands) sits
behind an early return.

**"Elastic-only" scopes the VARIABLE, not the EFFECT.** `m_PreElastic` is read at the three
`GetElasticModuli` overloads and in no other expression — but `K` and `G` leave those functions
by reference and are then used by the plastic machinery: `Stress_Correction` (including its
low-`p` rescue `dLambda = (p_min - p)/K`), `IntersectionFactor` /
`IntersectionFactor_Unloading`, `GetElastoPlasticTangent`, and the plastic multiplier itself,

```
NextDGamma = (2G n:de_dev - K de_v (n:r)) / (Kp + 2G(B - C tr(n^3)) - K D (n:r)),   Kp UNfloored.
```

So the floor changes `L` mid-path — `sqrt((p + pRe)/p)` is `1.41` at `p' = 1` kPa (+41 % on `G`),
`1.095` at 5 kPa, `1.077` at the BVP ring's 6.25 kPa — and the committed curve moves with it.
What it does **not** change is the destination: `eta = M^b` at the bounding state is a strength
identity in which the moduli do not appear, which is why §7.4's capacity gate reads `< 1e-5`
while §7.9's pre-failure curve reads `+6–8 %`. **Those two numbers are the same fact, not a
trade** — a stiffness floor is capacity-neutral and path-changing by construction, and that is
exactly why ADR-93 §7.6 replaces the "committed curve inside the certificate tolerance"
requirement with capacity neutrality rather than negotiating the tolerance.

**It does nothing in the ELASTIC stage — by vanilla's design, not by omission.** See §4: while
`mElastFlag == 0` the three `GetElasticModuli` overloads take the branch
`G = G0*P_atm*(2.97-e)^2/(1+e)` **without** the `sqrt(pn/P_atm)` factor, so the gravity / K0 leg
is pressure-INDEPENDENT and `pn` — hence `-pRe`, and `-Pmin` with it — is computed and unused. A
stage-0 leg is bit-identical at `pRe = 0` and `pRe = 1e6` alike, and `-Pmin 0.0101` vs `10.0`
changes nothing there either. There is no confinement dependence for a confinement floor to
floor. Where the factor **is** live (`mElastFlag == 1`) the initial elastic operator is re-derived
with the floor and scales by exactly `sqrt((P_atm + pRe)/P_atm)`: measured `2.000000000` on all
six modes of an unstrained `stdBrick` at `pRe = 3*P_atm`, after `updateMaterialStage ... 1` +
`reset()`.

**Explicit dynamics.** The floor raises the assembled stiffness, so it shortens
`CentralDifferenceLadruno`'s critical time step in the same proportion — `dt_cr` scales as
`1/sqrt((p + pRe)/p)`, up to 41 % shorter at `p' = 1` kPa. Budget for it before adopting a value
on an explicit deck.

**What WP-106 measured, so you know what to expect — this is not free speed.** On the ADR-93
ring path the floor cuts NOTHING: that dumped point is confined (min `p` 6.4–6.5 kPa), the
substep count moves under 1 % and in the WRONG direction, and committed stress moves ~3.5–4 % at
`pRe = 1` kPa — because `sqrt((6.5 + 1)/6.5) = 1.074` is a 7 % move on `G` wherever `p` is
small-ish, not only where the model has no answer. **On a real strip footing it is worse:** at
the campaign's own 7.65 kPa surcharge the live free-surface ring sits at `p ~ 6.25` kPa, and
`pRe = 1` kPa there costs **1.93× the substeps per step**, takes the worst single Gauss point
from 1387 to **18 831**, reaches LESS settlement in the same wall budget, and moves the
load–settlement curve **+5.9 %**. The floor pays only where a point genuinely reaches `p -> 0`
— at `p0 = 0.5` kPa it is 324× cheaper — and an embedment or surcharge that has already removed
that state removes the reason for the floor with it. Adopting a value is a declared modelling
statement about small-strain stiffness at low confinement, and it has to be paid for on the deck
that uses it: report the substep census (the `substeps` response) and the peak resultant on both
arms first. The WP-106 numbers are in `93_ladruno_sanisand_zero_confinement_adr.md` §7, with
the instruments in `Ladruno_implementation/wp106_pre_floor/`.

**Parallel / restore.** The request crosses the fork wire (`Vector(35)`, slot 34) and is carried
by `getCopy(const char*)` to every Gauss point, so an MP rank and a database-restored material
run the same elastic law as the process that wrote them — the ADR-86 §3 defect, which is what
this family of flags exists to make impossible.

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

**"Process-wide" means *within one model*, since WP-104.** `wipe` zeroes every process-wide
counter above — `implexRefusals` slots 0–3 and 5, all of `implexGuards`, `avgImplexError` and its
commit-round marker — so a fresh material in a new model reads `[0,0,0,0,0,0]` whatever ran
before it in the process (before WP-104 it inherited the previous model's totals, measured as
`[9,0,0,9,0,9]` on a fresh tag by the apeGmsh live test). `reset()` / `revertToStart()` do **not**
zero them: they rewind the *same* model, and a multi-leg campaign reads its legs as deltas over a
running total (`LEDGER_quirks`, "read it as DELTAS"). The 10-per-process `opserr` throttles are
untouched by either. Regression: `tests/test_wp104_implex_refusals_wipe_reset.py` and the
classic-Tcl twin `tests/tcl/wp104_implex_refusals_wipe.tcl`.

Since WP-86d, `implexError`, `avgImplexError`, `implexDetail` and `implexRefusals` also emit
`output.tag("ResponseType", ...)` in `setResponse` (the `FSAM`/`ASDConcrete3DMaterial` idiom), so a
`recorder Element -xml`/`-file` or the fork's own `recorder ladruno` names each column instead of
falling back to the generic `C1..Cn`. `implexGuards` is unchanged (out of WP-86d's scope).
| response | slots | meaning | ResponseType name(s) |
|---|---|---|---|
| `implexError` | 1 | total error, this material's last commit | `implexError` |
| `avgImplexError` | 1 | process-wide running mean over all commits | `avgImplexError` |
| `implexDetail` | 6 | `[0]` total error · `[1]` deviatoric leg · `[2]` volumetric leg (`sqrt(3)\|dp\|`) · `[3]` `p_min` clamp fired on the last pass (0/1) · `[4]` clamp fire count, ever · `[5]` the `f` **actually used** for the last extrapolation, frozen for this step (reads `0` on a guarded step, §11; under `-implexFactor control\|controlIter` this is `f*`, not the clock ratio — §12) | `implexDetail_total`, `implexDetail_dev`, `implexDetail_vol`, `implexDetail_clampFired`, `implexDetail_clampCount`, `implexDetail_f` |
| `implexRefusals` | 6 | `[0]` total GENUINE refusals (`[1]+[2]+[3]`) · `[1]` D2 sign-change · `[2]` `-implexControl` past tolerance · `[3]` companion hit `-maxSubsteps` · `[4]` **this integration point's** commit-refusal latch, 0/1 — the only per-instance slot (WP-99) · `[5]` POST-latch refusals, process-wide; fires once per Newton iteration per point and is deliberately **not** folded into `[0]` (WP-99) | `implexRefusals_total`, `implexRefusals_signChange`, `implexRefusals_control`, `implexRefusals_companion`, `implexRefusals_commitLatched`, `implexRefusals_latched` |
| `implexGuards` | 7 | `[0]` floor fallbacks (P2-1, `-implexFloor implicit`) · `[1]` guard firings (P2-2, `f = 0` after a reversal/softening commit) · `[2]` holds preserved (P2-3, zero-`dt` commits left alone) · `[3]` reversal resets restored (P2-5, `-reversalTol`) · `[4]` trial-time `f = 0` fallbacks (P2-6, `-implexTrialGuard`) · `[5]` hold-skip commits (P2-5c, once per point per hold) · `[6]` control-factor back-offs (P2-9, steps where `f* < 0.5·f_max`) | (none yet — out of WP-86d's scope) |

Python:

```python
r = ops.eleResponse(eleTag, "material", intPtNum, "implexDetail")
total, dev, vol, clampFired, clampCount, f = r

refusals = ops.eleResponse(eleTag, "material", intPtNum, "implexRefusals")
n_total, n_signchange, n_control, n_companion, commit_latched, n_latched = refusals
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

Since WP-86d both also emit a `ResponseType` tag, so a recorder names the column `psi` /
`yieldDistance` (the response's canonical name, not its alias) instead of the generic `C1`.

| response | slots | meaning | ResponseType name |
|---|---|---|---|
| `psi` (alias `stateParameter`) | 1 | the state parameter `psi = e - e_c(p')`, the model's own `GetPSI` with `p' = p + p_residual` floored at 1e-10 — the psi behind `M^b` and `M^d`. With the fork's default `p_r = 0` this is plain `e - e_c(p)` from `state[24]` and the mean stress. | `psi` |
| `yieldDistance` (alias `yieldFunction`) | 1 | the yield-function value `f = |s - p' alpha| - sqrt(2/3) m p'` on the committed pair: negative inside the cone, `~mTolF` (1e-7 default) on it, never positive after a converged return. | `yieldDistance` |

```python
psi = ops.eleResponse(eleTag, "material", intPtNum, "psi")[0]
f   = ops.eleResponse(eleTag, "material", intPtNum, "yieldDistance")[0]
```

Vanilla `ManzariDafalias` does not answer either name (empty response). Both are inherited by the
3D and plane-strain wrappers. Test: `tests/test_ladruno_sanisand_responses.py`.

### 6.2 `substepStats` — what the integrator cost, per point, since `revertToStart` (WP-127, F20(a))

`substeps` (ADR-86b) reports the LAST update only, and is zeroed by the next one — including the
zero-increment settle pass `analyze` pushes through every material after a failed step. So right
after a failed `analyze` it reads 0 everywhere (the TIMs ring dump). `substepStats` is the
post-mortem counter: **every column is per integration point** (one material instance), none is
process-wide; the cumulative ones count since `revertToStart` (`reset()`), are **not** reset by
`revertToLastCommit`, cross `getCopy` and the MP/database wire. Reading them never changes a number
(byte-identity pinned, `tests/test_ladruno_sanisand_replay_counters.py`). Columns 0-16 instrument
`ModifiedEuler` (IntScheme 1, and 0 through `MaxEnergyInc`); columns 17-27 (WP-130) instrument
`BackwardEuler_CPPM` (IntScheme 2, and the `-meFallback cppm` retry). 28 columns in all.

| slot | name (`substepStats_*`) | meaning |
|---|---|---|
| 0 | `updates` | material updates (`integrate()` calls, any stage/scheme) |
| 1 | `meCalls` | `ModifiedEuler()` calls |
| 2 | `substeps` | substep ATTEMPTS; closes as `[3]+[4]+[5]+[7]+[8]+[9]` |
| 3 | `accepted` | passed the error test |
| 4 | `rejectedErr` | failed it at `dT > dT_min`, retried smaller |
| 5 | `forcedAtDTmin` | **failed it at `dT == dT_min` (1e-6) and was ACCEPTED anyway** — elastic tangent, `alpha` re-derived, see `LEDGER_quirks` "ACCEPTS a substep that FAILED" |
| 6 | `forcedClampMc` | of `[5]`, those where the radial `eta -> Mc` stress clamp fired |
| 7 | `rejectedLowP` | `p < p_r` inside a substep: `dT` cut by 10 |
| 8 | `abandonedLowP` | `p < p_r` at `dT == dT_min`: ModifiedEuler **returned at `T < 1`**, silently |
| 9 | `capHits` | `-maxSubsteps` fired (update refused) |
| 10 | `entryPminClamps` | stress below `p_min + p_r` on entry: rebuilt at `p_min` |
| 11 | `pnResets` | `explicit_integrator`'s `p_n < p_r` reset (`sigma := p_min I`, `alpha := 0`) |
| 12 | `maxSubstepsOneUpdate` | most substeps any single update took |
| 13 | `lastSubsteps` | substeps of the last update that entered ModifiedEuler |
| 14 | `lastForcedAtDTmin` | `[5]` for that update |
| 15 | `lastAbandonedLowP` | `[8]` for that update |
| 16 | `lastCapHit` | 0/1 for that update |
| 16' | `lastCapHit` = **2** | the cap hit was RESCUED by `-meFallback cppm` (refused cap hits = `capHits` - `meFallbackOk`) |
| 17 | `cppmCalls` | top-level `BackwardEuler_CPPM` calls (IntScheme 2 plastic-branch updates + ME fallbacks) |
| 18 | `cppmNewtonFail` | local Newton (+ `Check`) did not return a valid state, at any halving level |
| 19 | `cppmHalvings` | recursive half-increment calls that did work |
| 20 | `cppmExplicitFail` | **vanilla's SILENT explicit fallback** after a failed Newton / exhausted ladder (F12 §5.2) |
| 21 | `cppmExplicitLowP` | the trial-`p < p_min` branch's explicit integration (by design) |
| 22 | `cppmRefusals` | the CPPM REFUSED the update (`-cppmOnFail refuse`, or inside the ME fallback) |
| 23 | `meFallbacks` | ModifiedEuler hit `-maxSubsteps` and the increment went to the CPPM |
| 24 | `meFallbackOk` | ... and the CPPM returned it (the update stands) |
| 25 | `lastCppmRefused` | 0/1 for the last update whose top-level CPPM call left the elastic branch |
| 26 | `cppmGuessTries` | `-cppmStart explicit`: local Newton restarted from the explicit guess |
| 27 | `cppmGuessOk` | ... and the gate ACCEPTED the root (admissible, within 2 % of the explicit walk) |
| 28 | `cppmLineSearchCuts` | `-cppmLineSearch on`: step halvings the search took |

`cppmOptions` (response id 33100, per instance, 7 values): onFail, halvings, lineSearch,
meFallback, start, tangentFixed, latchCause (0 -implex companion, 1 CPPM refusal).

```python
s = ops.eleResponse(ele, "material", ip, "substepStats")
substeps, forced, abandoned, caps = s[2], s[5], s[8], s[9]
```

### 6.3 Replaying one material point: `ladrunoSANISANDReplay` (WP-127, F21)

Puts a **private copy** of the `nDMaterial LadrunoSANISAND` prototype into a given state and drives
one strain increment through the same `setTrialStrain` an element uses. No element, domain or
analysis; nothing survives the call (the stage flag is forced to 1 and `ops_Dt` set for the call,
then both restored).

```
ladrunoSANISANDReplay $matTag -convention compressionPositive|tensionPositive
    -sigma s11 s22 s33 s12 s23 s31   -alpha a..6   -alphaIn ai..6   -fabric z..6
    -voidRatio $e   -dStrain d11 d22 d33 g12 g23 g31
    <-type 3D|PlaneStrain> <-trace $maxRecords (10000)> <-dt $dt (1.0)>
    <-primed 0|1 (1)> <-prevIncrNorm $norm (0)>
```

- **`-convention` is required.** `compressionPositive` = the model's internal `mSigma`
  (and compression-positive strain); `tensionPositive` = what `eleResponse ... stress` returns and
  an element strain. `alpha`, `alpha_in`, `z` are ratios and are never flipped. Shear strain is
  engineering (gamma). **The TIMs ring CSVs are compression-positive** despite their README
  (`LEDGER_quirks`, finding A).
- `alpha`, `alpha_in`, `z` are projected to their deviatoric parts (warning above 1e-6 relative);
  the given traces are returned.
- The state is committed through the base `commitState()`, so `K`, `G` and `e` are exactly what a
  converged step leaves for the next one: replaying step k+1 from the committed state of step k
  reproduces the analysis' stress (pinned to 1e-9, `test_replay_reproduces_an_analysis_step`).
- `-primed`/`-prevIncrNorm` feed the ADR-92 P2-5 reversal-noise guard; the defaults (armed, 0) are
  a plastic point with no history. `-dt 0` makes the call a hold.
- `-type PlaneStrain` requires `d33 = g23 = g31 = 0`.

Returns one flat list (format 1): `[1, rc, nStats, nRec, 5, dropped]` (`nStats` = 28 since WP-130), the `substepStats` columns of
this one update, 34 state values (`sigma` in the request convention, `alpha`, `alpha_in`, `z`, `e`,
`p` (compression-positive), `q`, `f` before, `f` after, path code, elastic ratio, given `tr(alpha)`,
`tr(alpha_in)`, `tr(z)`), then `nRec` records `T, dT, err, code, atDTmin`. Path codes: -1 not the
explicit path, 0 elastic, 1 start outside the yield surface, 2 elastic->plastic, 3 plastic,
4 unload-then-plastic, 5 `p_n < p_r` reset. Trace codes: 0 accept, 1 reject (error), 2 forced at
`dT_min`, 3 forced + `Mc` clamp, 4/5 low-p cut, 6 low-p abandon, 7 cap. The trace buffer lives only
for the call and is capped (`dropped` counts the rest). Classic Tcl prints the list with `%35.20f`,
so read tiny `err`/`dT` from Python.

Python helper (reads the attached CSVs, runs the documented probes):
`Ladruno_scripts/sanisand_replay.py` — `replay(...)`, `read_ring_csv(...)`, `probes(delta)`;
`python -S <bootstrap> Ladruno_scripts/sanisand_replay.py --delta 1e-5`.

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

### `-implex` without `-implexControl` used to walk past its own failed commits — WP-99 (F7)

Until WP-99 a commit-time companion failure propagated **nowhere**. `Domain::commit()`
(`SRC/domain/domain/Domain.cpp`) is `elePtr->commitState();` with the return value dropped, so no
element — fork or vanilla — can refuse a step at commit. When `ladrunoImplexCommit()` found the
`-maxSubsteps` cap hit, it committed the partially integrated state `ModifiedEuler` had left at
`T < 1`, returned `LADRUNO_MATERIAL_REFUSED` into that discarding caller, and the analysis
continued reporting every step converged. The 10-warning-per-process budget meant a long run said
so ten times and then went quiet. What that buys, measured: the TIMs plane-strain strip
(`LadrunoQuad` bbar `PlaneStrain`, 9 720 Gauss points, `-maxSubsteps 1000`, `-Pmin 0.0101`) reached
**25.9 million** commit-time companion cap hits and drew a **straight-line load–settlement curve to
2 674 kPa** with **every step reported converged** — a number with no mechanics behind it at all.
Since WP-99 such a commit **fails the step, on every element type**. Two mechanisms, and review
round 1 measured why both are needed:

1. **The commit is aborted.** The material declares the refusal out of band —
   `ladrunoNoteCommitRefusal()` in `SRC/material/LadrunoMaterialStatus.h` — and `Domain::commit()`
   checks that counter after its element loop and returns a failure.
   `AnalysisModel::commitDomain()` turns that into `-2` and `StaticAnalysis` /
   `DirectIntegrationAnalysis` / `VariableTimeStepDirectIntegrationAnalysis` all revert and return
   `-4`. This does not go through the element at all, which is the point: it works under a
   DISCARD element exactly as under a FORWARD one.
2. **The instance latches.** The trial is restored from the committed state,
   `ManzariDafalias::commitState()` is skipped, and every later `setTrialStrain` on that point
   returns `LADRUNO_MATERIAL_REFUSED`. That is the second line of defence, for a driver that
   ignores what `analyze()` returns.

**Mechanism 2 alone was not enough, and this is the measurement that decided it.** Two stacked
`stdBrick` — a DISCARD element, see the roster in `LEDGER_quirks.md` — with the lower element
starved and `algorithm Linear` ran **20 further accepted steps** with `analyze() == 0` throughout,
the refusing element frozen as a rigid inclusion, and *more quietly than before the WP*: a latched
`commitState()` returns early, so the old 10-per-process cap warnings stopped firing too.

**What "commits nothing" means, precisely:** *that Gauss point* commits nothing. `Domain::commit()`
walks the nodes first and then the elements, and every integration point is its own material
object, so by the time one refuses, the nodes and its sibling points have already committed
(measured: `[1,0,0,1]` latched across the four points of one `LadrunoQuad`). The model state at
that commit is **inconsistent**, not merely un-advanced — which is why it is aborted rather than
repaired.

The latch is sticky and is cleared only by `revertToStart()`: a driver that subdivides keeps being
refused and gives up, which is the intended outcome, because the step the latch is about was
already *accepted* and there is nothing to revert to. **`-implexControl` is the way to get a
*recoverable* refusal**: it catches the same cap one phase earlier, at the trial, so the step it
belongs to fails and a retry at a smaller increment is meaningful. Read the latch through slot 4
(`commitLatched`) of the `implexRefusals` response; genuine companion cap hits stay in slot 3, and
post-latch refusals — which fire once per Newton iteration per point — get their own slot 5 rather
than polluting the counters (review round 1 measured 244 of them against 4 real cap hits before
that split).

**Why not simply propagate the element's `commitState()` return code?** Because ADR-33/34 forbids
it. `ASDConcrete3D` and friends return *negative* "best-state" codes from commits that are
perfectly valid, and failing the step on those was measured to break mesh-objectivity gates that
had been green for months. The fork's rule is that only a **declared** refusal fails a step, never
any nonzero code — and at commit there is no sentinel-filtering element in the path to tell the two
apart, so the declaration has to arrive out of band. That is exactly what the counter is.

### Self-weight bearing decks — ADR-92 F10

- **On a SELF-WEIGHT bearing deck, `-implexControl` can be what STOPS the run — and control-off
  buys reach, not a confirmed curve.** Measured on a `B = 1.5` m self-weight strip
  (`gamma' = 9.81`, K0 = 0.455, 2 280 Gauss points, `-maxSubsteps 1000`, `-Pmin 0.0101`,
  ADR-92 F10): with `-implexControl 0.05 0.01` and the fork's usual ADR-63 D16 halve/double
  controller the leg refuses 724 times and dies on the harness step FLOOR at `s/B = 0.0085`;
  with `-implexControl` simply **removed** — bare `-implex`, the same doubling controller,
  nothing else changed — the same deck reaches `s/B = 0.0500` in **104 steps, 0 subdivisions,
  0 failed attempts, 58 s**, refusal ledger `0/0/0/0`. That is how the ADR-95 campaign which
  reached `s/B = 0.15` on this material was run (`sanisand_path_diag.py` passes only `-implex`).
  **Two things follow, and they are different things.** (a) *Termination*, and it is measured:
  the commit-time companion — the thing §3's hard requirement is about, which `-maxSubsteps`
  buys — integrated every increment it was handed on the control-off arms
  (`implexRefusals[3] = 0` on the campaign's B, C, D, E, K, L, N, N1 legs; `<= 42` anywhere,
  M 42 / F1 29 / I 12 / H 6 / F3 3 / J 2). So the discipline for a control-off leg is
  **read the companion bucket `implexRefusals[3]` at the end of every run**, which since
  WP-99 / PR #838 (merged as `c75edc95c`) is **belt-and-braces**: a capped companion commit now
  ABORTS the run — `Domain::commit()` fails and `analyze()` returns `-4`, the subsection above —
  so on any build from that merge on the watch confirms what the engine already enforces. On an
  older build it is the ONLY thing that would catch it: the F10 campaign itself ran on
  `9c2f964`, which predates #838, and there a capped companion commit is silent.
  (b) *Accuracy is untouched by any of this.* §8 and the constructor echo this
  material prints on every control-off run (`LadrunoSANISAND.cpp:2092`-`:2095`) say IMPL-EX was
  measured unusable from `d_eps = 5e-4` at `p0 = 5 kPa`, and **that deck sits inside that
  range**: its minimum `p'` is 6.374 kPa (1.27x the corner) and the control-off leg's strain
  increment crosses `5e-4` at `s/B = 0.0012` and runs at 2.6-4x the corner to the target, with
  no implicit anchor past `s/B = 0.00227`. **Do not generalise "the control is not needed at low
  confinement" from that campaign** — it measured what terminates, not what is right.
- **If you keep `-implexControl` on a self-weight deck, its refusal count is set by your
  CONTROLLER's growth rule, and its FLOOR is set by the `implexPrimed` test.** `implexError` is
  first order in the step (measured on a controlled refinement at a fixed committed state: max
  error 1.13e-2 / 5.78e-3 / 3.04e-3 / 1.67e-3 at `ds` = 8e-5 / 4e-5 / 2e-5 / 1e-5 m, onto a
  `dt`-independent floor of ~4e-4), while `-implexControl` bounds it **absolutely** — so a
  halve-on-failure / double-after-N controller can only find the bound by crossing it, with the
  clock ratio `f = dt_{n+1}/dt_n` sitting at 2 on exactly the crossing step. Measured at
  `tol 0.05`: growth ×2 → 724 refusals and `FLOOR` at 0.0085; ×1.25 → 253 and 0.0113; **×1.0 → 6
  refusals and 0.0265 with zero subdivisions**. (With the control OFF, ×2 and ×1.0 agree to
  0.383 % and ×2 is 6× faster, so the growth rule is innocent on its own.) The `FLOOR` itself is
  a separate mechanism: Gauss points whose committed plastic history is 1e-12…1e-21 pass
  `implexPrimed`'s bare `> 0.0` test (`LadrunoSANISAND.cpp:3021`), lose the un-primed exemption,
  and are refused on an error that does **not** decay with `dt` (0.2243 at `|dt| = 4e-5` →
  0.2143 at `2e-5`) — the asymptote `:3005`-`:3014` documents. No subdivision clears that.
  Secondary levers, measured: `tol = 0.1` (the shipped default) +37 % reach; `tol = 0.5` reaches
  the target with 18 refusals and agrees with the control-off arm to **0.395 % mean / 1.625 %
  max** over `0.002 ≤ s/B ≤ 0.05`, i.e. it makes the control nearly inert; `reductionLimit 0.5`
  +71 %, but understand it as "switch the control off after one halving" (at the shipped `0.01`
  the material floor `reductionLimit·|dt0| = 2e-7` m is the harness `DS_MIN` exactly, which is
  why it never fires). NOT levers: `-implexFactor controlIter` (+5 % for 2.9× the wall time),
  and more confinement (a 100 kPa surcharge leg has **zero** over-tolerance Gauss points at every
  step size tested and still refuses 2 754 times).
- **Never quote a load from a leg that refused.** Equilibrium is found on the extrapolated
  stress, so a leg that refused, halved and re-grew reaches a given settlement on a different
  strain path and reads a **stiffer** curve: +1.72 % at `s/B = 0.0030`, +2.56 % at 0.0050,
  +5.99 % at 0.0070 and **+19.40 % at 0.0085** against the refusal-free arm on the same deck.
  Compare reach and refusal counts freely; compare `q` only against an arm that ran refusal-free.
  Full tables and the three-candidate verdict: [[92b_implex_selfweight_wall_note]].
- **F10's bare-`-implex` recipe did not transfer to a finer strip (TIMs, 2026-09-18).** On the
  TIMs act's own self-weight strip, bare `-implex` (control off, as above) aborts at
  `s/B = 0.0004` on the footing-edge Gauss point on every build: the fork's F10 deck was too
  coarse to resolve the edge, so its 0.0500 reach says nothing about a mesh that does.

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

### `IntScheme 2` (BackwardEuler_CPPM) -- qualified per increment, refuted as the BVP integrator -- WP-105 (F12)

Scheme 2 is **not** a drop-in replacement for the deck default, and it is not simply worse either --
it depends entirely on whether the strain increment is *given* or *proposed*. On a replayed strain
path (zero free DOF, the increment supplied rather than found by a global Newton) scheme 2
integrates the **same model to the same limit** as scheme 1 -- `1.3e-3` / `2.9e-3` maximum relative
stress deviation over the whole path at `dEz = 1e-5` (`p0 = 100` / `20 kPa`), terminal `eta` within
`4.2e-4` / `2.0e-4` -- and at the campaign's own increment (`dEz = 1e-4`) it is **3.7-4.3x more
accurate and 4.2-7.6x cheaper** than scheme 1; at `4.6e-4`, **7-30x more accurate and 10-13x
cheaper**, with scheme 1 itself the one leaving its own bounding surface at `p0 = 20 kPa`
(`eta/M^b = 1.056`). **It is refuted as the primary integrator of a load-controlled BVP.** Under a
global Newton at `dEz >= 1e-4` it stalls in **8 of 8** free-standing drained-triaxial arms (scheme 1
stalls in 1 of 8), and each failing step burns **12-134 s** grinding `BackwardEuler_CPPM`'s
recursive-halving ladder against a 30 ms normal step -- up to **4400x**. On the real CP1/ADR-95
bearing leg (`x10z8`, `h1.0_e0.6944`, 1200 s budget) it committed 11 steps to `s/B = 4e-5` against
the scheme-1 baseline's 51 steps to `s/B = 0.019` in the same wall clock -- **475x shallower for the
same wall clock** -- with `ds` pinned at 25x the subdivision floor and every one of its committed
steps on the relaxed rung 3. **Use it where the increment is already given** -- a prescribed-strain
material-point study, or (open question) as the commit-time `-implex` companion itself, which runs
at `commitState` on an increment nothing proposes off-path and therefore sits in the regime where
scheme 2 wins; `-implex` was OFF in every WP-105 arm, so that companion question is unresolved, not
answered favourably. Do not reach for scheme 2 as the strip's primary integrator on this evidence.
Full numbers, the replay/free-standing/floor/bearing tables, and the "could not verify" list:
[[Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md]] (also see
`LEDGER_quirks.md` for the two related defects this same study found: the `-maxSubsteps` inertness
warning is wrong for scheme 2, and a CPPM non-convergence is invisible end to end).

### `IntScheme 2` under a global Newton -- the WP-130 flags (TIMs F18(c)/(d))

```
nDMaterial LadrunoSANISAND ... 2 2 ...                  (IntScheme 2, TanType 2)
    <-cppmOnFail explicit|refuse>   default explicit (vanilla)
    <-cppmHalvings n>               0..9, default 9 (vanilla: up to 2^9 half-increments)
    <-cppmLineSearch on|off>        default off
    <-cppmStart trial|explicit>     default trial (vanilla)
    <-cppmTangent fixed|vanilla>    default FIXED (owner decision); vanilla = ManzariDafalias' WRONG SIGN, reproduction only
nDMaterial LadrunoSANISAND ... 1 ... -maxSubsteps N
    <-meFallback cppm|off>          default off; needs IntScheme 1 and -maxSubsteps > 0
```

All defaults are vanilla's control flow, **byte-identical** (seven IntScheme-2 decks incl. a
free-DOF Newton deck, `tests/wp130_sanisand_byteid.py`) -- **except the tangent sign**: on
LadrunoSANISAND `-cppmTangent fixed` is the DEFAULT (owner decision, WP-130), so an IntScheme 2 +
TanType 2 deck hands its elements a different (correct) tangent than before. Every deck NOT on
IntScheme 2 + TanType 2, and vanilla `nDMaterial ManzariDafalias` everywhere, is bit-identical;
`-cppmTangent vanilla` reproduces the old binary bit for bit. On zero-free-DOF decks only the
`tangent` response changes (its sign); with free DOF the global Newton path changes. None is qualified with `-implex` (the
parser refuses the combination). A flag that could not act on the deck is refused.

- **`-cppmTangent fixed` -- the DEFAULT on LadrunoSANISAND.** Vanilla's CPPM hands the element MINUS its
  algorithmic tangent (`NewtonSol`: `Cep = -1.0 * CSigma`): a negative-definite stiffness, so the
  global Newton diverges from its first iteration and only a Krylov/relaxed rung ever commits a
  step. `fixed` hands out `+CSigma` -- the right SIGN -- and also corrects the low-p D_factor
  derivative in the local Jacobian (vanilla: wrong sign, and it drops the dilative D < 0 branch).
  **It is still not a fully consistent tangent** (review r1): it is one local iterate stale (and the
  local convergence norm mixes strain and stress units, so at the default TolR 1e-7 the staleness
  reaches 0.27-0.53 relative on a shear column); the void-ratio dependence (eps -> e -> psi) is
  missing from dR/deps (1e-4..1e-3 on the volumetric column); and after a SUCCESSFUL halving the
  tangent handed out is the second half-increment's (O(1) errors). The plane-strain wrapper hands
  out the same object (FD-checked). `-cppmTangent vanilla` is kept for reproduction only.
- **A CPPM refusal under a DISCARDING element** (SSPquad, stdBrick, BbarBrick, SSPbrick, BrickUP,
  LadrunoSolidShell) is caught at `commitState`: the refusal is declared to `Domain::commit()`
  (WP-99's channel), the commit aborts and the point latches -- analyze < 0, nothing drifts. **The
  latch is sticky until `reset()` (revertToStart)**: a smaller step does NOT clear it (measured:
  a 1000x smaller step still returns -3). Under a discarding element the only recovery is a
  restart; use a forwarding element (quad, LadrunoQuad/Brick, u-p family) so the step is cut.
- **`-cppmOnFail refuse`**: where vanilla, after a failed local Newton and the halving ladder,
  integrates the increment explicitly and reports success, the material REFUSES
  (`LADRUNO_MATERIAL_REFUSED`), so a forwarding element fails `Domain::update` and the step is
  cut. With `-cppmHalvings 0` that happens on the first try: measured 8-22 ms per refused step on a
  one-quad deck against 5.6 s at the defaults. The trial-`p < p_min` explicit branch is kept (it is
  the designed low-p route, and ModifiedEuler's own `-maxSubsteps` guards it).
- **`-cppmStart explicit`** -- NOT in the recommended recipe (review r1): when the local Newton from
  the elastic trial fails, retry it once from a 50-substep ForwardEuler guess before halving. A root
  found that way is ONE backward-Euler step over an increment on which the ladder would have
  halved, so it is less accurate. The acceptance gate (dGamma >= 0, p > 0, and agreement with the
  explicit walk to 2 %) rejects most of the bad ones, but on the review's 300-increment oracle set
  the gated guess's error is more than twice the default ladder's on 105 of the 171 increments
  where a guess was accepted (74 if the excess must also exceed 0.01 absolute); median relative
  error 0.021 vs 0.007; the worst oracle-converged increment is 0.54 with the guess vs 0.51
  without; the largest single-increment degradation is 0.09 -> 0.40 (confirmation review round 2
  recounted these from `wp130_f18c/review_r1/p2_guess_vs_oracle_after_gate.txt`). Use it only where
  speed is worth that.
- **`-cppmLineSearch on`**: backtracking (halving, at most 8 cuts) on the residual norm the local
  convergence test reads; a full step is taken if no cut helps.
- **`-meFallback cppm`** (IntScheme 1, F10b(b)): when ModifiedEuler hits `-maxSubsteps`, the SAME
  increment goes to `BackwardEuler_CPPM` (halving allowed, NO explicit exit); the update is refused
  only if the CPPM fails too. One-element test: a leg that `-maxSubsteps 20` refuses at step 1 runs
  all 10 steps with the fallback, stress within 1.3 % of the uncapped integration.
- **Recommended recipe** (IntScheme 2 under a global Newton): `2 2 ... -cppmOnFail refuse
  -cppmHalvings 3 -cppmLineSearch on` (`-cppmTangent fixed` is the default). With SAS-ME
  (IntScheme 129, WP-129) on the same model, the three refusal sources -- SAS-ME, the ModifiedEuler
  cap, the CPPM -- all reach the element with the same code and are all caught at commit under a
  discarding element; the latch warning names which one. `refuse` without
  `-cppmHalvings` now bounds the ladder at 3 by itself (<= 15 local Newtons per refused update).
- **What it buys on a BVP** (`Ladruno_files/testbed/hypo_bearing/wp130_f18c/`): F12's bearing deck (x10z8, `h1.0_e0.6944`, 1200 s budget, TanType 2, driver unchanged), the RECOMMENDED recipe (`fixed` default + `-cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on`, no `-cppmStart`; build 6726f5e24, `wp130_f18c/tables_recipe.md`), measured back to back with an IntScheme-1 control on the same box, which was at 100 % CPU (so compare these two with each other only): the recipe reaches s/B 0.00293 / 0.00421 / 0.00523 / 0.00626 at 300 / 600 / 900 / 1200 s against IntScheme 1's 0.00138 / 0.00250 / 0.00442 / 0.00698 -- ahead at 300, 600 and 900 s, BEHIND at 1200 s -- with 3.2 global iterations per committed step against 16.8 (481 of 527 steps on the plain Newton rung), load-settlement within 0.2-1.4 % of IntScheme 1, and its refusals had spent 76 of the driver's 80 pinned subdivisions when the wall stopped it. Vanilla IntScheme 2 reaches 0.00002 and `-cppmTangent fixed` alone 0.00378 (earlier, unloaded runs). The global Newton is NOT quadratic even with the fixed tangent: median observed order 1.14 on the last three residuals (9 % of committed calls >= 1.8); see the four tangent error sources. (Pre-round-1 note, superseded: an arm WITH `-cppmStart explicit` -- whose accuracy review round 1 measured and rejected -- reached 0.00876 in 1081 s on an unloaded box.).
  Recipe measured there: `2 2 ... -cppmTangent fixed -cppmOnFail refuse -cppmHalvings 3
  -cppmStart explicit -cppmLineSearch on`. Without `fixed`, no combination of the other flags got
  past s/B 0.0002.
- **WP-128's smallest reproducer** (`sigma = 0.0101 I`, `alpha = alpha_in = z = 0`, plane-strain
  `d eps_yy = 1e-4`): ModifiedEuler returns `alpha/alpha^b` 5.10 in one accepted substep; the CPPM
  returns 0.18 with rc 0 in one local Newton (every variant), against ~0.27 from WP-128's alpha-aware
  references (`wp130_f18c/q128_reproducer.txt`). The implicit return does not escape the bounding
  surface there; its one-step error is its own.

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

### `alpha_in` at the stage flip — `init` by default since WP-112 (`-flipAlphaIn`, P2-7 / F14)

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

**Default moved to `init` (WP-112, TIMs F14, 2026-09-18) — the paragraph above is the P2-7c
record, and it missed one case.** The sign test reads only the SIGN of
`(α_n − α_in_n) : Ce : Δε`, with no magnitude guard on either factor, and vanilla runs it in the
elastic stage too. A `LoadControl(0)` hold's `Δε` is solver noise, so each hold sets
`α_in := α_n` at a coin-flip of points, and afterwards `α_n − α_in_n` is only the round-off by
which `α` has moved since. On the first plastic step the direction of a round-off perturbation
then picks the branch — and the thread count of MKL's solve is such a perturbation. The TIMs
self-weight strip (9 720 Gauss points) read the **first push step at 1.511 / 1.824 / 1.824 /
1.489 kPa at 1 / 2 / 4 / 8 MKL threads under `vanilla`** on Windows (1.597 / 1.824 / 1.824 on
Linux), and the branches it opened were 30 % apart by `s/B = 0.035`; under **`init` the same leg
reads 1.824 / 14.339 / 36.586 kPa at rows 1 / 8 / 15 on every thread count and both builds**. The
fork reproduces the mechanism on a 12 × 6 `LadrunoQuad` self-weight strip
(`tests/test_ladruno_sanisand_flip_determinism.py`): two elastic holds leave 229 of 288 Gauss
points with `‖α − α_in‖` below `1e-12` of `max(‖α‖, m)`, and vanilla's first push step then
reads 4.107 / FAIL / 8.332 / 9.483 kN/m after 0 / 1 / 2 / 3 holds, while `init` reads 9.659 kN/m
after every one of them (to 3e-14) and is bit-identical at 1 / 2 / 4 / 8 threads for ten steps.
So `init` is the default now. `-flipAlphaIn vanilla` stays, for reproducing real
`ManzariDafalias` (A/B against vanilla decks and golden files), and prints a warning once per
Gauss point (10 per process) when a plastic-stage trial meets `0 < ‖α_n − α_in_n‖ ≤ 1e-8 ·
max(‖α_n‖, m)`. **The R3 numbers do not move:** P2-7c measured the Esmeralda implicit twin under
`init` identical to `vanilla` to the digit (6.511 / 11.539 / 16.117 / 20.528 kN, "the RC14 price
on this column is zero"). What does move is the IMPL-EX dense refuse arm's `fixed` reference wall,
0.01689 under `vanilla` vs 0.01754 under `init` (also P2-7c) — P2-9's ship/refute bars were set
against the former, and `controlIter` has not been re-measured under `init`. One limit stays: on
the fork's deck the `init` curve still moves with the hold count from step 4 on (step 10: 35.5 –
36.1 kN/m after 0 – 4 holds), because the holds change the committed state by round-off and later
branch decisions amplify it. `init` removes the flip's sign lottery, not every sensitivity of the
model to its state.

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

**Status: CLOSED 2026-09-08 — measured, not shipped as default.** `-implexFactor fixed` is
the default and is byte-identical to every build before P2-9; that does not change here.
Two opt-in control modes exist; both require `-implexControl`. `control` (f* frozen at the
first trial of the step) was run through the plan's Fork R3 registered arm and
**REFUTED**: depth 0.052 < the 0.076 bar and overlay 11.1 % mean deviation, both worse
than `fixed` on the same deck (`_adr92_p2_9_r3_results.md`, Leg 1). `controlIter` (f*
recomputed at every Newton trial from that trial's own `d_eps`) **PASSES** the same
bars — depth 0.115, overlay 1.70 % — but costs materially more wall time and Newton churn
(see below). Its Esmeralda dense-refuse arm has now been measured
(`_adr92_p2_9_esmeralda_results.md`): honest wall **0.01755**, above the plan's
refutation bar (0.0169) but 0.85 % short of its ship bar (0.0177). Per the pre-registered
decision rule's otherwise branch, **`controlIter` does not ship as default** — it is
recorded here as a graded guard with its measured gain (reach +3.9 %, refusal churn down
~400x) against its measured cost (~2.9x wall time on this deck; the R3 registered arm's
~13x figure does not generalise). Both control modes are documented as measured findings,
not as recommended settings; `fixed` remains the default and the thing to reach for. See
"When to reach for `controlIter`" below for the practical guidance.

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

### R3 verdict — `control` REFUTED, `controlIter` PASSES R3

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
first `-implexFactor` mode to independently clear both R3 bars, but R3 was not the final gate —
the Esmeralda dense-refuse arm was.

### Esmeralda verdict — CLOSED, otherwise branch, not shipped

TIMs' Esmeralda dense/loose-refuse arm (`_adr92_p2_9_esmeralda_results.md`, engine `179da6ffb`,
PR #822) is the decision rule's actual gate (`_adr92_p2_9_control_informed_f_plan.md` §2: "ships
if the dense refuse wall reaches >= 0.0177 ... otherwise the ADR records the factor as a graded
guard"). Measured dense honest wall on `controlIter` (leg 146607, control 0.1/0.01): **0.01755**.
That clears the refutation bar (`< 0.0169`) but misses the ship bar (`>= 0.0177`) by 0.85 %, so
neither branch of the rule's `if` fires and the **otherwise branch is decisive**:

- **P2-9 does not ship as a default.** `-implexFactor fixed` remains the default — already the
  code state, so this closeout is documentation only, no code change owed.
- **`controlIter` is recorded as a graded guard with its measured gain**, not promoted to
  default: reach +3.9 % over the P2-7c fixed-f wall (0.01689 -> 0.01755), overlay
  comparable-to-better (+0.28 % mean vs +0.31 %), and a collapse in refusal churn (refusals
  42 545 -> 102, failed attempts 248 -> 11) at a measured cost of ~2.9x wall time on this deck
  (1 933 s -> 5 620 s at 17.2 it/step) — the R3 registered arm's ~13x figure does **not**
  generalise to this deck.
- `control` (frozen f*) remains **REFUTED** and is not a candidate in any form.
- The loose arm ends at 0.03921, unchanged at 0.039 to the third figure exactly as
  pre-registered — a PASS, not a finding. A companion tol-0.01 dense variant (146608) reaches
  further (0.01932) by spending its full subdivision budget rather than refusing honestly, but
  sits +1.44 % HIGH against the twin (stiffer, not closer) — a finding, not a candidate
  configuration.
- **P2-8's fixed threshold (`-implexGuardKp`, listed, not built) remains the documented
  fallback** if a future measurement wants a graded guard without `controlIter`'s cost.

### When to reach for `controlIter`

`controlIter` is not a default and not a general recommendation, but it is a real, measured
tool for one specific situation: a **deep, dense push toward a softening seat** where the P2-2
guard's committed-predecessor threshold is under-firing and refusals are dominating wall time.
On Esmeralda it bought a few percent more reach at comparable-or-better overlay accuracy while
cutting refusals by roughly two orders of magnitude, for a wall-time premium of roughly 3x on a
production-scale deck (not the smaller R3 deck's ~13x — that figure does not generalise, and
should not be quoted outside the R3 deck). Reach for it when:

- the deck is a dense, quasi-static push through a shear zone approaching a softening or
  reversal point (the seat P2-2/P2-9 both target), and
- refusal churn (not overlay accuracy) is the binding cost, and
- ~3x wall time is affordable for the run.

**Know what you are spending.** Measured on Esmeralda's dense column, iterations per committed
step (cumulative Newton count, failed attempts included) run **2.2 for `fixed`, 17.2 for
`controlIter`, 52.4 for the implicit twin** — the same ~8x also holds on the loose column
(2.1 -> 17.0). IMPL-EX exists to buy cheap steps from a frozen operator, and recomputing `f*`
per trial spends most of that: `controlIter` keeps roughly a 3x per-step edge over implicit
where `fixed` keeps ~24x. Wall time rises only 2.9x rather than 8x because the refusals it
removes were wasted work. If your run is already iteration-bound rather than refusal-bound,
this trade is against you — measure before switching.

It is **not** a fix for the p = 0 confinement ring (ADR 93, `93_ladruno_sanisand_zero_confinement_adr.md`)
— that wall is the material's, not the extrapolation factor's, and `controlIter` does not touch
it. It requires `-implexControl` (the companion computation this factor is built on) and is
refused without it. `-implexFactor fixed` (the clock-ratio default, unchanged since before
P2-9) remains what every deck should reach for unless the situation above applies.


## 13. Choosing an IntScheme — and SAS-ME (`IntScheme 129`, WP-129)

The scheme is the 20th positional argument (`IntScheme`, after the 18 model parameters and the
tag). What each one is, measured against WP-134's independent oracle (`uw_model`: the DM04 rate
equations with the UW constitutive additions, integrated exactly):

| IntScheme | what | use it? |
|---|---|---|
| **1** ModifiedEuler (the fork's default) | explicit Heun, stress-only error at a hardcoded `1e-4` (unless `-honorTolR 1`), moduli frozen at the committed state (U9), a loading stage with a negative denominator read as elastic + uncapped step growth (the "err = 0 path"), force-accept at `dT_min`, a drift correction that can give up with `f > 0` | the calibrated default; know its quirks rows. Ring states: α can leave the bounding surface (WP-128). Benign states: up to 15–100 % stress error on 1e-4 increments against the oracle (WP-129 §13.3) |
| **2** BackwardEuler_CPPM | implicit; under TanType 2 a SIGN-FIXED (WP-130 `-cppmTangent fixed`, the LadrunoSANISAND default) but NOT fully consistent tangent -- one local iterate stale, no void-ratio term, the second half's tangent after a halving (§9) | accurate per increment; under a global Newton use §9's recipe (WP-105; WP-130) |
| **45** RungeKutta45 | explicit Sloan RK45 | **not a reference**: dT_min 1e-3 hard-coded, Mc-clamp force-accept, no drift correction, and `dAlpha3/dAlpha4` never computed (α weights sum to 301/336) |
| 3, 5 | RK4 / Forward Euler, no error control | no |
| 0, 4, 6–9 | MaxEnergy / MaxStrain wrappers | no; IntScheme 4 is even non-deterministic (uninitialised moduli) |
| **129 SAS-ME** | this section | when the answer at low confinement / after reversals matters more than the cost |

### 13.1 Syntax

```tcl
nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m $h0 $ch $nb \
    $A0 $nd $z_max $cz $Rho  129 $TanType $JacoType $TolF $TolR \
    <-errFloor $sigRef> <-alphaBoundTol $kappa> <-alphaProject 0|1> \
    <-sasAlphaIn reseat|bracket|stale> <-sasErrorVars full|stress> \
    <-sasHFloor $cA> <-sasReseatHyst $cRev> <-sasSoftCap $kappa> \
    <-maxSubsteps $n> <-Pmin ...> <-Presidual ...> ...
```

The three `-sas{HFloor,ReseatHyst,SoftCap}` flags are WP-151's opt-in DM04 variant (§13.4). All are
OFF by default, and with them off SAS-ME is byte-identical to WP-129.

- `TolR` IS the substep tolerance of the PLASTIC part (`-honorTolR` is inert and warned); the
  elastic part is exact (closed form, below), so it has no tolerance to honour. Recommended range
  **1e-4 to 1e-7**: `1e-4` is ModifiedEuler's scale, `1e-7` lands on the oracle to ~1e-7. Below
  ~1e-8 the first-substep error of a large low-p increment (~6e2·dT² at 20 kPa, 1e-3 shear) cannot
  meet the tolerance above `dT_min = 1e-6`, so the update refuses (`errorAtDTmin`) or hits
  `-maxSubsteps`: a global cut then handles it, at a cost. The default `1e-7` is inside the range.
- `-errFloor` σ_ref of the stress error `‖dσ₂−dσ₁‖ / max(2‖σ‖, σ_ref)`; default `P_atm/101`
  (1 kPa at P_atm 101 — exactly ModifiedEuler's implicit floor, WP-128 §5.1). α and z use the unit
  reference (`max(2‖α‖, 1)`): the α error is a stress error in units of p. **The floor is not the
  lever** — at low p the cost is stability-limited (WP-128 §5.3).
- `-alphaBoundTol κ` (default 0.1): ρ_α = √(3/2)‖α‖ / α^b(θ_α, ψ) with α's OWN Lode angle. An
  accepted substep is rejected (refused at dT_min) only when PLASTIC FLOW carried α outward past
  1 + κ (ρ_α with the substep's end surface is larger than with its start α) — never because ψ moved
  the surface: the continuum itself carries α outside when elastic compression raises ψ (review of
  #871: proportional compression from ρ_α 0.999 at 20 kPa gives 1.13 at 430 kPa, 2.04 at 12.5 MPa).
- `-alphaEntryTol κ_e` (default 2): a START with ρ_α > 1 + κ_e is refused
  (`startAlphaOutsideBounding`); 1 + κ < ρ_α ≤ 1 + κ_e is counted (`entryOverKappa`), not refused.
  Why 2: on that compression path ρ_α reaches 3 only past ~30 MPa, far outside the model's range,
  while the dumped TIMs states b8 1950/2-3 sit at 6.8/7.3 (b:n = −8.2).
- `-alphaProject 1`: instead, project α radially onto the bounding surface (the deviatoric stress
  follows by p·Δα, so f, n, p and ψ are unchanged); counted. OFF by default: it rewrites history,
  and the stress jump is ≥ 9 % of p·‖α‖ whenever it fires (review of #871).
- `-sasAlphaIn`: `reseat` (DEFAULT) = the paper's rule — wherever (α − α_in):n < 0 a new loading
  process starts, α_in := α there; integrate()'s once-per-increment trial test is undone.
  `bracket` = keep UW's trial test and only use h = 1e10 where (α − α_in):n ≤ 0. `stale` and
  `-sasErrorVars stress` reproduce ModifiedEuler's defects G and E — attribution only.
- `-implex` is refused with 129 (not qualified as a companion). `-maxSubsteps` binds (refusal).
- TanType 1 and 2 both return the continuum tangent at the end state; 0 = Ce at the end state.

### 13.2 What it does, per substep

1. Elastic predictor, EXACT: with G = g·√max(p + pRe, p_min) and K = cG, √p is linear in the
   volumetric strain (√x = √x₀ + c·g·t·dε_v/2 above p_min, linear below) and ∫G dt is exact for the
   deviatoric part (review of #871: one Heun step was 4–32 % off, independent of TolR). Elastic if
   f_trial ≤ TolF. Otherwise: on the surface with (α − α_in):n < 0, α_in := α (the oracle's t = 0
   rule); the loading test on the TRUE gradient ∂f/∂σ = n − ⅓(n:α + √(2/3)m)I; the intersection by
   Pegasus on the SAME exact path (unload-then-reload: 64 samples to bracket the exit), so the
   plastic part starts on the surface.
2. Two Heun stages, each evaluated entirely at its own state (K, G, n, b, d, h, D, B, C — U9).
   Stage classification from N = ∂f/∂σ : C : dε: N ≤ 0 elastic (α, z unchanged); N > 0, H > 0
   plastic; N > 0, H ≤ 0 has no plastic solution — REFUSED at stage 1, cut at stage 2.
3. Error on σ, α, z; accept iff err ≤ TolR; q = clamp(0.9√(TolR/err), 0.1, 1.1), no growth after a
   rejection.
4. Drift correction (consistent with σ, α, z; then normal); if neither direction reduces |f| the
   substep is cut, refused at dT_min — never returned with f > TolF. Both-sided (|f| ≤ TolF) only
   for an all-plastic substep that STARTED on the surface; otherwise only f > TolF is corrected.
5. ρ_α check (above). α_in, the paper's rule, decided only ON the surface: a stage whose start is
   on the surface with (α − α_in):n < 0 re-seats there; a reversal detected at stage 2 of a plastic
   substep cuts the substep (so it is located at a substep start), re-seating at dT_min only; an
   accepted substep ending on the surface with it negative re-seats at its end.

Refusal codes (`sasStats` column `lastRefuseCode`, and the warning text): 1 startOutsideYield,
2 startAlphaOutsideBounding, 3 startInadmissible (trace / tension / non-finite), 4 errorAtDTmin,
5 loadingNonPosH, 6 tensionAtDTmin, 7 driftFailed, 8 alphaOutsideAtDTmin, 9 maxSubsteps. A refusal
leaves the trial on the committed state and returns `LADRUNO_MATERIAL_REFUSED` (element roster:
LEDGER_quirks "element refusal roster").

**Discarding elements** (SSPquad, stdBrick, BbarBrick, the SSP/brick u-p variants, LadrunoSolidShell,
...: the roster) drop that code, so their Newton "converges" on the refused state. WP-129 (review of
#871) makes the COMMIT refuse instead: `commitState` sees the refused update, declares it to
`Domain::commit()` (the WP-99 channel), the analysis step fails, and the point latches (cleared by
`revertToStart`). The same now holds for the ModifiedEuler `-maxSubsteps` cap, which used to commit
the strain without the stress. Use a forwarding element (quad, LadrunoQuad/CST/LST, LadrunoBrick,
the u-p family) to get a recoverable, cuttable refusal.

### 13.3 Measured (WP-129, `Ladruno_files/testbed/wp129_sasme/`)

- **Oracle, benign** (K0 states 20/50/100 kPa × active/passive/shear × 1e-5/1e-4): SAS-ME at TolR
  1e-7 within 2e-7 relative of `uw_model`; at TolR 1e-4 within 5e-5. ModifiedEuler (campaign
  options): 6–15 % on 1e-4 increments, 100 % on the reproducer.
- **WP-128 reproducer** (σ = 0.0101 I, α = α_in = z = 0, dε_yy = 1e-4): ρ 0.252, η 0.531 in 251
  substeps (oracle 0.252 / 0.531; ModifiedEuler 5.14 in 1 substep).
- **Ring** (80 rows × ± iso, ± shear at 1e-6, 1e-5): 624/640 integrated, the 16 of b8 1950/2-3
  refused `startAlphaOutsideBounding`; max f at exit 1e-7, max ρ_α 0.983 (= the start value); median
  substeps 17 (ModifiedEuler 4), p95 153 (82), max 520 (1336). Against the oracle (`uw_model`, 624
  admissible cases): median 8e-6 / p95 5e-5 / max 2.3e-4 relative in σ at TolR 1e-4; median 4e-9 /
  p95 3e-8 at TolR 1e-7.
- **Reversal chains** (WP-128 vertUnload / extShear from p0 2 kPa at 1e-4): ρ ≤ 0.61 / 0.46
  (ModifiedEuler 5.2 / 7.1); the increments that drive p to the floor are REFUSED
  (38–39 of 90: `tensionAtDTmin` / `errorAtDTmin`) where ModifiedEuler resets the stress to p_min·I — use smaller
  increments there, or accept the global cutback.
- **Cost**, smooth monotonic chains: 4–6 substeps per 1e-5 increment (ModifiedEuler 1–3); the
  profile split at a ring state is ~60 % stages (half state-dependent quantities), 10 % drift,
  6 % α check; at a deep state the tangent and drift are ~11 % each.


### 13.4 Re-seat regularization — WP-151 R1, an OPT-IN DM04 variant

> [!important] The TIMs setting under SAS-ME (`IntScheme 129`): the FULL set
> The owner said "if R1 makes sense, let's use it", held the merge to see the full set past the old wall, and
> then authorized it ("ok, when ready merge"). All three were relayed by the TIMs orchestrator on 2026-09-28.
> κ stays an owner/TIMs choice.
> - **The configuration is the full set: `-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5`.** A partial set is
>   an ablation, not a lighter fix. On the footing (B/8, E_B settings; E_B walls at s/B 0.0508):
>   - the floor alone walls at 0.0453 and the hysteresis alone at 0.0461, both EARLIER than DM04;
>   - floor + hysteresis without the cap walls at 0.0525, on a b:n < 0 post-peak point near a reversal;
>   - the full set had 0 `loadingNonPosH` at s/B 0.054 and was still hardening.
>   - **The post-wall gate is pending:** the full set to its end, and κ 0.25 / 0.75 beside 0.5.
> - The flags stay opt-in in the code (default OFF, byte-identical to DM04).
> - **The cap is a constitutive choice**, for the owner/TIMs to decide.
>   - It enforces H ≥ κX: post-peak softening per unit plastic strain is bounded at (1−κ) of the elastic
>     projection X.
>   - It binds only within about a cone radius of a reversal, where DM04 gives K_p → −∞. It never binds in
>     DM04's regular softening.
> - **More footing evidence** ([[151_sanisand_reseat_singularity]] §9):
>   - q–s within ±0.22 % of E_B below s/B 0.03;
>   - 0 NonPosH past E_B16's wall on B/16;
>   - cost within 1.4 % of E_B.
>   - A B/4 leg reached s/B 0.127. That is a coarser mesh than E_B's B/8, so it does not compare with E_B's
>     wall like for like.
> - **Do not substitute a recalibration of the Lode parameter c.** A c = 0.80 footing walls too (s/B 0.048),
>   on compression-side states where DM04 runs the same re-seat sequence (memo §2.5.1).
> - It is not a cure for mesh-dependent localization (WP-150). It also does not help low-confinement surface
>   points: the B/4 leg's eventual limiter is `errorAtDTmin` at p′ → 0, a p′-floor question.

**What it is for.** Near the peak, with the campaign set's thin yield cone (m = 0.005) and
near-neutral or rotating loading, the exact DM04 rate equations re-seat α_in again and again in
finite pseudo-time (a Zeno accumulation). Meanwhile b:n → 0⁺ and |dα/dt| → ∞. The discrete image of
this is the `loadingNonPosH` refusal wall of the WP-138 footing: the refusals match the oracle's
failures one-to-one on the wall states. Full study:
[[151_sanisand_reseat_singularity]].

```tcl
... 129 $TanType $JacoType $TolF $TolR -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 ...
```

| flag | equation (ρ_c = √(2/3)·m, the yield-cone radius) | DM04 |
|---|---|---|
| `-sasHFloor c_A` | h = b0 / max((α−α_in):n, c_A·ρ_c), bounded everywhere | h = b0/((α−α_in):n), ∞ at a re-seat |
| `-sasReseatHyst c_rev` | α_in := α only when (α−α_in):n < −c_rev·ρ_c (a FINITE reversal) | … when < 0 |
| `-sasSoftCap κ` | where b:n < 0: h ≤ (1−κ)X/(⅔p\|b:n\|), so H = K_p + X ≥ κX | no cap |

- **Use all three together.**
  - On the wall fan (5 states × 64 trials), floor alone and hysteresis alone each leave 97–102 of 320
    trials singular; together they leave 0.
  - Past the old wall, though, the footing meets genuine b:n < 0 states near a reversal, where the floored
    h still drives H ≤ 0. At fh's final wall point, the oracle fails 32 of 64 trials without the cap and 0
    with it (κ 0.25–0.75).
  - On the footing, the floor alone turns the refusals into `maxSubsteps` (572 log mentions vs 91) and
    walls at 0.0453.
  - The hysteresis alone keeps `loadingNonPosH`: its first comes at s/B 0.0362, where E_B has its first at
    0.0363, and it walls at 0.0461.
- **Recommended: c_A = 1, c_rev = 1, κ = 0.5.** No parameter is fitted: c_A and c_rev are in units of the
  calibrated m.
  - **c_rev ≤ 1 for cyclic work.** At c_rev = 2, ten CVSS cycles at γa 1e-5 never re-seat (0, against DM04's
    20), and τ differs by 6.5 % of τ_max.
- **What changes in calibrated behaviour** (oracle, campaign set, p0 100 kPa; memo §6):
  - The elastic range (γ ≲ 3e-6) is identical.
  - Just past it R1 is slightly SOFTER, because it starts plastic flow at a re-seat where DM04 is still
    elastic (h = ∞ there). At the same strain:
    - CVSS τ is −1.1 % at γ 1e-5, −0.33 % at 3e-5 and −0.09 % at 1e-4;
    - undrained TC q is −1.3 % at ε_a 3e-6 and −0.5 % at 1e-5.
  - Cyclic CVSS at γa 1e-5 (c_rev 1) stays within 1.2 % of τ_max. At γ ±0.1 %, drained cycles are identical
    to 4 digits.
  - Monotonic tests to large strain: |Δq| ≤ 2.7e-4·q_max. Undrained cycles to liquefaction are identical.
  - The cap never binds in an element test.
  - **Take G0 from the elastic range** (γ ≲ 3e-6), where R1 = DM04.
- **Census** (`sasStats`, appended columns): `sas_hFloored` counts stages where the floor bound,
  `sas_hSoftCapped` stages where the cap bound, and `sas_reseatHeld` the re-seat DECISIONS a sub-threshold
  reversal held (one reversal can be counted at several stage and substep tests, so this is not a count of
  distinct reversals). **`sasStats` is now 36 long** (columns 0–32 unchanged). A consumer that hard-codes 33
  breaks: the WP-138 deck driver's census did ("broadcast (36,) into (33,)"). Read the length from the
  response, or zip against `sanisand_replay.SAS_NAMES`.
- It removes the singular set and the re-seat chatter. It does **not** regularize strain
  localization (mesh dependence): that is WP-150 R2/R3.

### 13.5 Tension cutoff (separation) — WP-152, an OPT-IN constitutive choice for near-surface sand

**What it is for.** After R1 (§13.4), the binding limiter of a dilating sand under the TIMs footing is the free
surface. A few surface Gauss points outside the footing edge go to p′ → 0, and SAS-ME refuses them:
- code 6 `tensionAtDTmin`;
- code 3 when the committed p ≤ 0;
- codes 4 and 9, when the α and fabric error, or the substep count, blows up as p → 0.

Vanilla hides this: `Stress_Correction` silently resets such a point to σ = (p_min + p_r)·I with α = 0. The cutoff
does the same thing openly, reversibly and counted. It is FLAC's tension cutoff done at the material level. Plan and
evidence: [[152_sanisand_tension_cutoff]].

```tcl
... 129 $TanType $JacoType $TolF $TolR -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 -sasTensionCutoff $pSep $pContact ...
```

- **E2 is the operative trigger, not tension.** At p → 0 a free-surface point fails SAS-ME's accuracy/cost limit
  (codes 4/9) BEFORE its trajectory crosses p = 0; the α and fabric error terms do not scale with p. With p_sep = 0 every
  Toyoura footing leg stopped at s/B ≈ 0.0044 on two top-row points just outside the edge, codes 4/9 only, zero code 6
  (WP-152, 7f1562c81). So the cutoff is a LOW-CONFINEMENT SEPARATION: a point whose update fails at committed
  p < p_sep under a non-compressing increment is treated as separated. E1 (tension) is a backstop. Do not use p_sep = 0.
- **SAS-ME's only low-p test is p + p_r > 0.** `-Pmin` is not an admissibility threshold under IntScheme 129: it only
  floors the elastic moduli.
- **Entry masks ONLY low-p/tension refusals** (gated after the 2026-09-29 review):
  - E1: code 6, or code 3 whose cause is p0 ≤ 0, **at committed p0 ≤ p0max** (`-sasSepMaxP0`, default and minimum
    p_contact: the entry jump stays at the re-contact scale) **and under a non-compressing increment** (a zero-stress
    start under gravity is a deck error and refuses). Above the bound one increment carried a well-confined point
    through p = 0: a step to cut, so it REFUSES.
  - E2: code 4 or 9 while the committed p0 < p_sep **and the increment does not compress** (tr Δε ≤ 1e-10·‖Δε‖,
    compression positive; the tolerance keeps isochoric shear, whose trace from B·u is round-off, non-compressing). An accuracy or cost failure of a low-p point being compressed (the B/8 top row sits at p′ 0.2–0.35 in
    situ, below p_sep 0.5) REFUSES. `p_sep = 0` disables E2, which leaves a pure tension cutoff.
  - A qualifying refusal that a gate holds back is counted (`sepHeldHighP`, `sepHeldCompressing`, per update call).
  - Code 5 (`loadingNonPosH`, the α_in singularity), code 2, a non-finite start, codes 7 and 8, and codes 4/9 at
    p0 ≥ p_sep still refuse, at any p.
- **While separated** (g = tr ε − tr ε_entry, compression positive):
  - no shear; model p = p_c(g) = p_min + K(p_contact)·max(g, 0). OPEN (g ≤ 0) the point sits at p_min and absorbs the
    strain; CLOSING (g > 0) it reloads isotropically and elastically. The point keeps its weight (body forces act on the
    nodes) and its place in the mesh.
  - The closing branch is CONTINUOUS on purpose. The first build jumped p_min → p_contact at re-contact, and wherever
    the point has compliance around it that leaves a band of load with NO equilibrium (for σ = p_contact the
    neighbours must yield, which reopens the gap). A free-node Newton column failed there.
  - α = α_in = 0 and the fabric is kept.
  - The tangent is C_e at p_min while open (the model's own moduli floor, a declared Newton regularisation) and,
    while closing and on the exit iterate, the consistent bulk K(p_contact) with shear G(p) (a regularisation: the
    stress carries no shear). `tangentEP` returns the same from the committed state.
  - The entry itself is a jump (from the SAS-ME trial at p0 to p_min). With p0 ≤ p_sep (E2) or ≤ p_contact (E1) it
    is small; on the footings the entry steps took more Newton iterations (median 17–28 vs 7–15) but no more step
    cuts.
  - Known limitation: a stage flip back to the elastic stage (`updateMaterialStage 0`) while separated leaves the
    separation flag set; the elastic stage then evolves the stress from p_min I and stage 1 overwrites it. Do not flip
    back while points are separated.
- **Re-contact** happens at g ≥ g_c = (p_contact − p_min)/K(p_contact): SAS-ME restarts from σ = p_re·I,
  p_re = p_c(g) ≥ p_contact, α = α_in = 0.
  - p_contact > p_sep, and entry needs a qualifying refusal, so a point cannot chatter.
  - Re-contact is VOLUMETRIC only: isochoric shear of a separated point never closes the gap (tested).
- **Parameters:**
  - The parser requires 0 ≤ p_sep < p_contact, p_contact > p_min, p0max ≥ p_contact, and SAS-ME. The cutoff is
    refused with `-implex`, with `-Presidual ≠ 0` (the separated state would carry tension), with `-pRe ≠ 0` (the
    moduli read tr σ/3 + p_Re, so K(p_contact) would be K(p_contact + p_Re)), and without `-sasHFloor > 0` (re-contact
    lands on α = α_in = 0, where h is otherwise the 1e10 sentinel).
  - Starting values: p_sep 0.5 kPa (TIMs D1's p-floor bound) and p_contact 1.0 kPa.
  - Report the limit load at p_sep and at p_sep/2 (the D1 rule, < 2 %).
- **Census** (`sasStats`, appended). Each transition is counted ONCE, when it commits; an element may call the
  update several times per step, and a cut step counts nothing:
  - `sas_sepEntriesTension` (E1) and `sas_sepEntriesLowP` (E2);
  - `sas_sepExits`;
  - `sas_sepActive`, the committed 0/1. Its sum over points is the number of points separated now;
  - `sas_sepLastCode` (the refusal code the last committed entry masked: 3, 4, 6 or 9) and `sas_sepMaxP0` (the
    largest committed p0 at an entry);
  - `sas_sepHeldHighP`, `sas_sepHeldCompressing` (per update CALL: qualifying refusals the gates refused. A held
    refusal fails its step, so it never commits; counting at commit would always read 0).
  - **`sasStats` is now 44 long.**
- **`sasOptions`** is 12 values long: indices 9–11 are p_sep, p_contact and p0max.
- **InitialStateAnalysis:** `revertToStart` under ISA keeps the separation state with the stress it produced; a plain
  `reset` starts the point NORMAL.
- **Solution control:** use a force test (`NormUnbalance`, as the footing driver does). A separated cluster carries
  only p_min, so a displacement test cannot tell a converged separated region from a stalled one.
- **It is a constitutive choice about near-surface sand.** The owner and TIMs decide p_sep and p_contact, and how
  separated points are reported in the capacity. The p_r bracket (§13.4's cross-check) stays available.
