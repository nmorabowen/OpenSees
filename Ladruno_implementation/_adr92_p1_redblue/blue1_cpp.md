# BLUE TEAM 1 — verification of RED 1 against the C++ (PR #798, `wp/92c-implex-p1`)

Worktree `C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr92-p1-implex`, HEAD
`70f6bf0e8`. All line numbers re-verified in that tree; where red's number is off by a few
lines the corrected number is given. Read-only: no edit, no build, no analysis run. Every
number below comes from `git`/`sed`/`grep` on the tree or from `python3.12` arithmetic on
committed CSV/JSON artefacts already in the repo.

Scoreboard: **F1 CONFIRMED, F2 CONFIRMED, F3 CONFIRMED (empirically, with arithmetic),
F4 CONFIRMED, F5 CONFIRMED (and worse), F6 CONFIRMED, F7 CONFIRMED (and the PR's own memo
asserts the opposite), F8 CONFIRMED (empirically), F9 PARTIAL, F10 PARTIAL — mechanism real,
claimed consequence REFUTED.**

---

## F1 — commit-time refusal is dead code. **CONFIRMED.**

Full call chain traced, all four links:

`SRC/analysis/analysis/StaticAnalysis.cpp:214`
```cpp
	result = theIntegrator->commit();
```
`SRC/analysis/integrator/IncrementalIntegrator.cpp:254`
```cpp
    return theAnalysisModel->commitDomain();
```
`SRC/analysis/model/AnalysisModel.cpp:655-659`
```cpp
    if (myDomain->commit() < 0) {
	opserr << "WARNING: AnalysisModel::commitDomain - Domain::commit() failed\n";
	return -2;
    }
```
`SRC/domain/domain/Domain.cpp:2241-2245`
```cpp
    Element *elePtr;
    ElementIter &theElemIter = this->getElements();
    while ((elePtr = theElemIter()) != 0) {
      elePtr->commitState();
    }
```
— return not assigned, not tested. `Domain::commit()` ends `return 0;` at `Domain.cpp:2309`.
The **only** `< 0` return `Domain::commit()` can produce is the ADR-73 porous-overlay branch
at `Domain.cpp:2277-2283` (`overlayRes < 0`), which requires a `PATTERN_TAG_LadrunoPorousOverlay`
load pattern in the model. With no such pattern the function is a total function returning 0,
so `AnalysisModel.cpp:655`'s test is unreachable and `LadrunoSANISAND::commitState()`'s
`LADRUNO_MATERIAL_REFUSED` is discarded twice over (once by the element loop, once by the
constant return).

Staleness on the early return at `LadrunoSANISAND.cpp:1641-1660` is exactly as claimed. The
early `return LADRUNO_MATERIAL_REFUSED` at `:1660` precedes every one of:

```cpp
:1670    int res = ManzariDafalias::commitState();     // <- skipped: mEpsilon_n, mSigma_n,
:1673    mImplexDEpsP = mEpsilon_n;                    //    mEpsilonE_n, mAlpha_n, mFabric_n,
:1674    mImplexDEpsP.addVector(1.0, mEpsilonE_n, -1.0);//   mAlpha_in_n, mDGamma_n all frozen at n
:1675    mImplexDEpsP.addVector(1.0, epsPOld,     -1.0);
:1677    mImplexDtCommit  = mImplexDt;
:1678    mImplexStepArmed = true;
:1679    mImplexTrialDone = false;
```
so after a refused commit `mImplexStepArmed` is **false** (cleared at `:1337`) and
`mImplexTrialDone` is **true**. Consequences, exactly as red states: the next step's
`ladrunoImplexTrial()` does not re-arm, so `mImplexFactor` stays frozen at the failed step's
value; `dEps = mEpsilon - mEpsilon_n` at `:1563-1565` spans two steps; `mImplexDEpsP` is the
increment of step *n−1*; and `commitState()` at `:1702` still routes to `ladrunoImplexCommit()`.
The comment at `:1642-1645` describes the hazard the code then commits.

Not refuted anywhere. No element propagates `commitState()` upward in a way `Domain` reads.

---

## F2 — `f` degrades to `alpha` on a negative clock. **CONFIRMED**, with the BVP leg quantified.

`SRC/material/nD/LadrunoSANISAND.cpp:1325-1335`:
```cpp
    mImplexDt = dt;
    if (mImplexDt0 <= 0.0 && dt > 0.0)
        mImplexDt0 = dt;

    if (mImplexDtCommit > 0.0)
        mImplexFactor = (dt / mImplexDtCommit) * mImplexOpt.alpha;
    else
        mImplexFactor = mImplexOpt.alpha;

    if (dt == 0.0)
        mImplexFactor = 0.0;
```
Sign-blind as charged. `ops_Dt` is confirmed to be the signed pseudo-time increment:
`Domain.cpp:2081-2082` `currentTime = timeStep; dT = currentTime - committedTime;`, published
at `Domain.cpp:2125` `ops_Dt = dT;` and again at `Domain.cpp:2392` in `Domain::update()`. The
deck drives `ops.integrator("LoadControl", -ds)` (`sanisand_tau0_band.py:722`), so `dt < 0` on
every push step and the `> 0.0` gate is never taken.

Gravity cannot rescue `mImplexDtCommit`: gravity runs at stage 0, where
`ladrunoTrialUpdate()` takes `if (!this->ladrunoImplexActive())` (`:1274`) and sets
`mImplexTrialDone = false`, so `commitState()` at `:1702` takes the base path and
`ladrunoImplexCommit()` — the only writer of `mImplexDtCommit` — never runs. Hence on the BVP
deck `mImplexDtCommit <= 0` and `mImplexDt0 == 0.0` for the entire analysis.

**What the noctl leg actually ran.** Reconstructed from
`adr92_bvp/implex_noctl/a2_h1.0_e0.6944_curve.csv` (142 rows). Note the CSV's `ds_mm` column
is written *after* the growth update (`sanisand_tau0_band.py:749-755`), so the applied
increment is recovered from the settlement column instead — `d_s(i) = s(i) − s(i−1)`:

| applied ds | steps |
|---|---|
| 0.020 mm | 1–6 |
| 0.040 / 0.080 / 0.160 / 0.320 / 0.640 / 1.280 / 2.560 mm | 6 steps each, 7–48 |
| 5.000 mm | 49–141 (93 steps) |
| 4.400 mm | 142 (truncated by `ds = min(ds, smax - s_now)`, `:721`) |

Steps at which `dt_{n+1}/dt_n ≠ 1`: **9 of 142** — steps 7, 13, 19, 25, 31, 37, 43 (ratio
exactly 2.0), step 49 (2.56→5.00, ratio 1.953125), step 142 (5.00→4.40, ratio 0.88).
`nfail = 0`, `nsub = 0`, `nrelax = 0` in `a2_h1.0_e0.6944.json`, so there were **no** halvings
and no reverts to complicate the sequence.

Spec `f = (dt_{n+1}/dt_n)·alpha` with `alpha = 1` (json `implex: true`, alpha echoed as 1 in
the engine log). Delivered `f ≡ 1.0` on all 142 steps. So the delivered operator differs from
the specified one on those 9 steps: the extrapolated plastic strain `f·d_eps_p(n)` was **0.5×**
specified on steps 7–43, **0.512×** on step 49 and **1.136×** on step 142.

Red's stated failure scenario (a `ds *= 0.5` halving giving delivered `f = 1.0` where 0.5 is
correct) did **not** occur in this leg — there were no halvings. The defect is real and fired,
but in the *growth* direction, and 8 of the 9 errors are **under**-extrapolation (delivered
closer to implicit than specified), not over-. That is the one place red's severity narration
is off; the mechanism, the gate and the "measured on a different operator" conclusion stand.

Refutation criterion checked and not met: `tests/test_ladruno_sanisand_implex.py:746` is
`ops.integrator('LoadControl', -1.0)` after a positive history, i.e. a **sign change** that
takes the `:1511` refusal and returns before `armStep()` computes any ratio. No test in the
file asserts `implexDetail[5]`.

---

## F3 — sign-blind reduction floor. **CONFIRMED**, and it is *the* cause of the ctl arm's death.

```cpp
:1326    if (mImplexDt0 <= 0.0 && dt > 0.0)
:1327        mImplexDt0 = dt;
```
```cpp
:1610    if (mImplexOpt.control) {
:1611        this->ladrunoImplexMeasureError(mSigma, sigImplicit, mEpsilon);
:1612        if (mImplexError > mImplexOpt.errorTol) {
:1615            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
:1616                return LADRUNO_MATERIAL_REFUSED;
```
Parser doc, `:311-313`, is a magnitude statement (`"as a fraction of the FIRST dt seen"`)
implemented as a signed comparison. Both branches red predicts are reachable:

* **Negative clock throughout** (the campaign deck): `mImplexDt0` never leaves `0.0` (traced in
  F2), so `mImplexDt0 <= 0.0` is permanently true and the floor **never disengages** — refuse
  forever.
* **Positive `dt0`, then negative `dt`**: `mImplexDt < 0 < reductionLimit·mImplexDt0` makes the
  `>=` false, so the control **never refuses**. (The `:1511` sign-change guard fires on the
  first sign flip only; once `mImplexDtCommit` is itself negative the guard is quiet and W7 is
  permanently disabled.)

**Arithmetic proof from the PR's own committed artefact.** `adr92_bvp/implex_ctl/a2_h1.0_e0.6944.json`:
`implex_control = [0.05, 0.01]`, `ds_min = 2e-07`, `mode = "FLOOR"`,
`verdict = "step collapsed to the DS_MIN floor at s/B = 0.00504 (every ladder rung failed at
ds = 0.000156 mm)"`, `nsub = 25`, `nfail = 75`. The leg was still being **refused** at
`ds = 1.5625e-7 m`. A magnitude-correct floor would be
`reductionLimit · |dt_0| = 0.01 × 2.0e-5 = 2.0e-7 m` (`ds_base = 2e-05`, json), so at
`1.5625e-7 m < 2.0e-7 m` the material was supposed to **stop refusing and accept**. It did not,
because `mImplexDt0 == 0.0` short-circuits the test. The last halving that killed the leg is
precisely the one the reduction limit exists to prevent.

Consequence for the results memo: `_adr92_p1_bvp_gate_results.md`'s finding 3 — *"With the
control ON at `tol = 0.05` the leg seizes at `s/B = 0.00504` … the tolerance as set does not
[work]"* — is **confounded**. The seizure is at least partly a `reductionLimit` that was
inoperative, not a tolerance that was too tight. The proposed remedy in §4.1 (a tolerance
sweep 0.2/0.5/1.0) is aimed at the wrong knob until F3 is fixed.

---

## F4 — refusal is dropped by non-propagating elements. **CONFIRMED**, with one scope correction.

Sentinel consumers, whole tree:
```
$ grep -rn "LADRUNO_MATERIAL_REFUSED" SRC/ --include=*.cpp --include=*.h
SRC/element/ladrunoBrick/LadrunoBrick.cpp:1034,1080,1111,1182,1827,3319
SRC/material/LadrunoMaterialStatus.h:74
SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h:2198,2350
SRC/material/nD/LadrunoSANISAND.cpp (the producer)
```
No other element. Confirmed at both named call sites:

`SRC/element/UWelements/SSPbrick.cpp:443-447` — return not even captured:
```cpp
	Vector strain(6);
	strain = Bnot*u;
	theMaterial->setTrialStrain(strain);

	return 0;
```
`SRC/element/brick/Brick.cpp:1069-1073` — captured and discarded:
```cpp
    success = materialPointers[i]->setTrialStrain( strain ) ;

  } //end for i gauss loop

  return 0;
```
`mImplexTrialDone = true` is set at `:1453`, first statement of `ladrunoImplexTrial()`, before
any refusal, so `commitState()` at `:1702` routes to `ladrunoImplexCommit()` for a step the
material declined. `ladrunoRestoreTrialFromCommitted()` at `:1344-1357` begins
`mSigma = mSigma_n;`. No construction-time guard on element class exists (the material has no
handle on its element).

**Scope correction, in the code's favour on two points.** (a) The two *loud* refusal sites
(D2 `:1523-1526`, companion `:1657-1660`) leave `mSigma` **independent of the trial strain**,
so on a non-propagating element the internal force stops responding to `dU` at those points
and Newton will usually stall rather than converge — a loud failure more often than a silent
one. The genuinely silent case is the **`-implexControl` W7 site** at `:1616`, which returns
the sentinel while leaving `mSigma = sigma~` (strain-dependent, post-clamp): there Newton
converges normally and a non-propagating element commits an answer the material refused, with
nothing anywhere in the log. Red identifies this inconsistency in F8 but files it under the
weaker site here. (b) The BVP gate itself is unaffected — `sanisand_tau0_band.py:597` builds
`ops.element("LadrunoBrick", ...)`, the one element that does propagate; the ctl arm's engine
log carries 10 `WARNING LadrunoBrick::update - element … the material REFUSED the trial strain`
lines and then its budget-10 suppression notice.

Verdict CONFIRMED: the return *is* dropped, no guard exists, and the W7 path is a true silent
wrong answer on `SSPbrick`/`stdBrick`/`BbarBrick`.

---

## F5 — `-maxSubsteps` mandatory on scheme 1, advisory on scheme 2. **CONFIRMED, and worse than filed.**

```cpp
:1182    if (s == 1 && mMaxSubsteps <= 0) {
:1189        return -1;
:1190    }
```
```cpp
:1192    if (s == 2 && verbose) {
...
:1205        if (mMaxSubsteps <= 0)
:1206            opserr << ... " -maxSubsteps is strongly recommended." << endln;
:1211    }
```
No `return -1`. **The aggravation red missed:** the scheme-2 no-cap warning is nested inside
`if (s == 2 && verbose)`. `verbose` defaults true only on the deck-command path;
`LadrunoSANISAND.h:199-200` declares `setLadrunoImplexOptions(const LadrunoImplexOptions &opt,
bool verbose = true)` and it is called with `false` at `:826`, `:834` (both `getCopy`) and
`:980` (`recvSelf`). So a parallel/`recvSelf` rank or any restored model gets **no warning at
all**, while the scheme-1 hard refusal at `:1182` still fires on those paths.

Red's refutation criterion is met against the code, not for it. `ManzariDafalias.cpp:1536`:
```cpp
        if (mMaxSubstepsInME > 0 && mSubstepsTakenInME > mMaxSubstepsInME) {
            mSubstepCapHitInME = true;
```
With `mMaxSubsteps == 0` the flag is unreachable, so `ladrunoImplexTrial()`'s only
companion-failure detector (`:1546 if (mSubstepCapHitInME)`) is dead and `-implexControl`
degrades to the error-magnitude test at `:1612` alone. Confirmed.

---

## F6 — `takeAverageError()` destroys the accumulator on read. **CONFIRMED.**

`LadrunoSANISAND.cpp:1074-1080`:
```cpp
    double takeAverageError(void)
    {
        double a = (count > 0) ? (sumError / (double)count) : 0.0;
        sumError = 0.0;
        count    = 0;
        return a;
    }
```
Call sites, both per-integration-point:
* `:1916` inside `getResponse(LadrunoSanisandAvgImplexErrorResponseID, …)` — `setResponse`
  (`:1899-`) builds one `MaterialResponse` per requested material point, so `-ele all -material 1
  avgImplexError` on a 200-element × 8-GP mesh performs 1600 destructive reads per record; the
  first returns the true mean, the other 1599 return exactly `0.0`.
* `:1767` inside `setParameter("avgImplexError")` — one extra destructive read at
  parameter-creation time.

`accumulate()` (`:1063-1071`) is a non-atomic read-modify-write on a function-local `static`
singleton (`:1057-1061`), i.e. a data race the moment ADR-75b Lane 3 threads the element loop.
`maxError` has no reset path at all: `ladrunoImplexInitState()` (`:1100-1116`, called from
`revertToStart`) touches only per-instance members, and `LadrunoImplexGlobals` exposes no
reset — so `implexError` (`:1762`) is a since-process-start maximum surviving `ops.reset()`.
Confirmed as filed, including the `maxError` rider.

---

## F7 — refusal warnings are unthrottled. **CONFIRMED**, and the PR's own memo states the opposite.

`grep -n "static int|WarnCount|warnCount|budget" SRC/material/nD/LadrunoSANISAND.cpp` returns
exactly one hit — `:100 static int numLadrunoSANISANDMaterials = 0;` — which is the ordinary
class-echo counter. Neither refusal block has a counter:

* D2 sign change, `:1512-1522`: a single `opserr` statement spanning 11 source lines,
  terminated by one `endln` at `:1522`. It emits **one line of 527 characters** per event (tag
  and the two `dt` values are the only variable parts). At the 4896 refusals commit `3c788778f`
  reports, that is **2.58 MB** of synchronous console output for one push step.
* Commit-time companion failure, `:1647-1656`: same shape, one `opserr`, no counter.

The precedent is in the direct parent, `ManzariDafalias.cpp:1451-1466`:
```cpp
        static int ladrunoClampWarnCount = 0;                               // Ladruno
        if (ladrunoClampWarnCount < 10) {                                   // Ladruno
```
with the rationale at `:1439-1450` (`"a 50k-element mesh holds ~400k instances, and a
per-instance echo of exactly this shape measured 83 MB of stderr during ADR-86 PR-1"`).

Two further facts red did not have. (1) `LadrunoBrick.cpp:952-960`, the element-side reporter,
justifies its own throttle by asserting *"An unthrottled line here would bury the material's own
(already throttled) diagnostic under its own noise"* — a premise that is false for these two
new sites. (2) `Ladruno_implementation/_adr92_p1_bvp_gate_results.md:84-85` states *"the
material's refusal message is throttled to 10 per process so the engine log cannot supply it
either."* No such throttle exists; the 10-line budget in the ctl engine log is
`LadrunoBrick`'s (`LadrunoBrick.cpp:970-980`), not the material's. The memo's stated reason for
disclaiming arm A's refusal count is wrong even though the disclaimer itself is prudent.

---

## F8 — the `-implexControl` refusal is silent. **CONFIRMED, empirically.**

`:1610-1617`, quoted whole — there is nothing between the tolerance test and the return:
```cpp
    if (mImplexOpt.control) {
        this->ladrunoImplexMeasureError(mSigma, sigImplicit, mEpsilon);
        if (mImplexError > mImplexOpt.errorTol) {
            // Below the reduction limit there is nothing left to cut, so refusing
            // again would only turn a bounded inaccuracy into a dead analysis.
            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
                return LADRUNO_MATERIAL_REFUSED;
        }
    }
```
No `opserr`, no `ladrunoRestoreTrialFromCommitted()`, no tangent refreeze — against
`:1523-1526` and `:1657-1660`, which do all three.

Empirical confirmation from the committed artefact: `implex_ctl/a2_h1.0_e0.6944_engine.log`
is 394 lines covering 25 refused/subdivided steps and 75 failed rungs, and contains **zero**
`LadrunoSANISAND` refusal lines (`grep -c "REFUSES this step"` → 0; `grep -c "COMPANION hit"`
→ 0). Every diagnostic naming a material refusal in that log comes from the element
(`LadrunoBrick::update`, 10 lines + suppression notice) and says only *that* a refusal
happened, never which of the three sites produced it.

Severity nuance worth stating: the missing state restore is defensible — on a propagating
element the step will be reverted anyway, and leaving `sigma~` is what makes the W7 path
usable — and, given F7, an *unthrottled* warning here would be the worst site in the file to
add one (W7 is the highest-frequency refusal by construction). The defect is real; the fix is
a throttled one-line notice, not a copy of the D2 block.

---

## F9 — a zero-`dt` step wipes the extrapolation history. **PARTIAL.**

Mechanism confirmed. `:1478-1487`:
```cpp
    if (mImplexStepArmed) {
        ...
        if (this->GetNorm_Cov(dArm) > 0.0) {
            this->ladrunoImplexArmStep();       // clears mImplexStepArmed
        } else {
            mImplexDt     = 0.0;
            mImplexFactor = 0.0;
        }
    }
```
and `:1677 mImplexDtCommit = mImplexDt;` with no zero guard, so a hold does commit
`mImplexDtCommit = 0` and (via `integrate()` on a zero increment) `mImplexDEpsP ≈ 0`. Red is
right about the state transitions.

Where it overreaches, on two counts:

1. **`f = alpha` after `dt_n == 0` is the specified fallback, not a defect.** The code's own
   D2 note at `:1445-1449` — *"`dt_n == 0` the ratio is undefined (first step after a hold, or
   the first step of all). Fall back to `f = alpha`, which is what the P0 oracle's default
   does"* — and the ADR §2 formula `f = (dt_{n+1}/dt_n)·alpha`
   (`92_ladruno_sanisand_implex_adr.md:86`) is simply undefined at `dt_n = 0`. Falling back to
   `alpha` is a documented choice matching the P0 oracle.
2. **`d_eps_p(n) = 0` after a genuine zero-strain step is the correct datum, not a corruption.**
   IMPL-EX extrapolates the *previous converged step's* plastic increment. If that step moved
   no strain, its plastic increment really is zero, and predicting zero plastic flow is the
   operator behaving as specified. Red's phrase *"a pure elastic predictor — indistinguishable
   from `-implexAlpha 0`"* describes correct behaviour on that one step, and standard IMPL-EX
   does the same after a null step.

What survives: a `LoadControl 0.0` hold placed **after** `updateMaterialStage 1` does silently
reset the extrapolation history for every Gauss point, so the first genuinely plastic step
after such a hold runs with no extrapolation and its reported `implexError` is not
representative of the steps around it. That is a documentation/diagnostic issue (nothing marks
the step), not the wrong-answer red claims. Downgrade to MINOR.

---

## F10 — `revertToLastCommit` gated on `enabled`, not `Active()`. **PARTIAL — mechanism real, consequence REFUTED.**

The two facts are correct. `:1708-1721`:
```cpp
LadrunoSANISAND::revertToLastCommit(void)
{
    int res = ManzariDafalias::revertToLastCommit();   // an empty `return 0` today
    if (!mImplexOpt.enabled)
        return res;
    ...
    mEpsilon = mEpsilon_n;
    this->ladrunoRestoreTrialFromCommitted();
```
and the parent is indeed inert — `ManzariDafalias.cpp:542-546`:
```cpp
int ManzariDafalias::revertToLastCommit (void)
{
    // need to be added
    return 0;
}
```
So the ADR's *"There is no branch to get wrong: the extrapolated path is simply not reachable"*
(`:1252-1257`) is literally false — this entry point runs under `-implex` at stage 0 while
`ladrunoImplexActive()` (`:1259-1263`, `enabled && mElastFlag != 0`) is false.

**But the divergence red claims does not occur, and the decisive line is
`ManzariDafalias.cpp:4864-4868`:**
```cpp
    const double& eG = mUseCurrentVoidRatioInG ? en : m_e_init;
    if (mElastFlag == 0)
        G = m_G0 * m_P_atm * pow((2.97 - eG),2) / (1 + eG);
    else
        G = m_G0 * m_P_atm * pow((2.97 - eG),2) / (1 + eG) * sqrt(pn / m_P_atm);
    K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
```
At `mElastFlag == 0` the moduli are **independent of the stress argument** (no
`sqrt(pn/P_atm)`) and, since `mUseCurrentVoidRatioInG` is false in all four constructors,
independent of the void-ratio argument as well. `ladrunoImplexFreezeTangent()` (`:1378-1384`)
calls `GetElasticModuli(mSigma_n, mVoidRatio, K, G)` and writes
`mCe = GetStiffness(K,G); mCep = mCe; mCep_Consistent = mCe;` — bit-for-bit the same three
writes `elastic_integrator`/`initialize` perform (`ManzariDafalias.cpp:938-941`). So the tangent
rewrite at stage 0 is a **no-op in value**, not a change of answer. And `mEpsilon = mEpsilon_n`
is immediately overwritten: `Domain::revertToLastCommit()` ends `return this->update();`
(`Domain.cpp:2385-2386`), which re-runs every element's `setTrialStrain` from the reverted
(committed) displacements, i.e. it writes back the same `B·u` that produced `mEpsilon_n`. There
is no recorder sampling point in between — recorders fire from `Domain::commit()`
(`Domain.cpp:2296-2299`), not from the revert.

Residual, genuine: `ladrunoRestoreTrialFromCommitted()` also re-derives `mVoidRatio` at `:1356`,
which vanilla's revert does not touch, so `getVoidRatio()` differs in the window between the
revert and the trailing `update()`. Nothing in the tree reads it there. Verdict: the ADR
sentence must be corrected and the gate should be `ladrunoImplexActive()`, but the stage-0
byte-identity claim is **not** broken by this path.

---

## The question red did not ask: was the noctl leg's `f` the spec value on every step?

**No — on 133 of 142 steps, yes; on 9, no.** Delivered `f ≡ 1.0` throughout (F2). Spec
`f = (dt_{n+1}/dt_n)·alpha` equals 1.0 on step 1 (`dt_n = 0`, documented `alpha` fallback,
`:1445-1449`) and on every constant-`ds` step. It differs only at the 9 increment changes
enumerated in F2 — the seven ramp-up doublings (steps 7, 13, 19, 25, 31, 37, 43; spec `f = 2.0`),
the capped doubling at step 49 (spec `f = 1.953125`) and the truncated final step 142
(spec `f = 0.88`). All nine sit in the ramp-up or the last step; **all 93 tail steps at
`ds = 5 mm` (steps 49–141) ran at the spec value.**

What this does and does not do to `_adr92_p1_bvp_gate_results.md`:

* **Survives.** Every claim in its "MAY say" list. The ladder claim ("0 of 252 steps needed a
  fallback rung") is a property of the step being *linear in the trial strain*, and a constant
  `f` is still constant within the step — the tangent identity `d sigma~/d eps = Ce(p_n)` holds
  regardless of `f`'s numeric value, so the frozen-SPD argument is untouched. Likewise
  `TARGET @ s/B 0.25`, `142 steps`, `0/80 subdivisions`, `tail = 95.9 %`.
* **Must be qualified.** Any statement that the leg exercised "the ADR §2 operator". It
  exercised §2's operator with `f` pinned to `alpha` instead of `(dt_{n+1}/dt_n)·alpha`; on 8
  of the 9 divergent steps that is a **halved** extrapolation, i.e. the leg was *more implicit*
  than specified at exactly the moments the step size grew. Any `implexError` figure quoted
  from this leg is the error of `f ≡ alpha`, not of the specified operator, and cannot be used
  to calibrate `-implexControl`'s tolerance (memo §4.1).
* **Must be withdrawn as stated.** The memo's finding 3, *"the tolerance as set does not
  [work]"*, per F3: the ctl arm died still refusing at `ds = 1.5625e-7 m`, below its own
  intended floor of `0.01 × 2.0e-5 = 2.0e-7 m`, because `mImplexDt0` was never set on a
  negative clock. The two arms do not bracket the truth about `tol = 0.05`; the ctl arm
  measured a broken floor.
* **Also wrong on a matter of fact.** Memo lines 84-85, "the material's refusal message is
  throttled to 10 per process" — see F7. The W7 refusal that produced those 25 subdivisions
  emits no message at all.
