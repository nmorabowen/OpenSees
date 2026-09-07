# RED TEAM 1 — C++ IMPL-EX on LadrunoSANISAND (PR #798, `wp/92c-implex-p1`)

Range reviewed: `235934592..70f6bf0e8 -- SRC/`. All line numbers are in the worktree
`C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr92-p1-implex`.

Vanilla footprint is clean (only the three `LadrunoSANISAND*` files changed; no classTag
touched; the 3308x ids are response/parameter ids, not class tags). That is the only thing
below that is not a defect.

---

## 1. BLOCKER — the commit-time companion refusal is DEAD CODE: `Domain::commit()` discards every element's `commitState()` return and unconditionally returns 0. The material silently stops advancing while the model does not.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1640-1662`, `SRC/domain/domain/Domain.cpp:2244`, `SRC/domain/domain/Domain.cpp:2309`

**CLAIM** When the companion hits `-maxSubsteps` inside `commitState()` — the *default* configuration, since `-implexControl` is off by default and `-maxSubsteps > 0` is *mandatory* for scheme 1 — the material refuses the commit, but nothing in OpenSees can see that refusal, so the analysis advances `committedTime`, records, and leaves this Gauss point's committed state frozen one step behind the rest of the model, permanently and silently.

**EVIDENCE**

`LadrunoSANISAND.cpp:1640-1662` (the refusal):
```cpp
    if (mSubstepCapHitInME) {
        // Reachable only with -implexControl OFF ...
        opserr << ... " The committed state is NOT advanced and this commit is refused ("
               << LADRUNO_MATERIAL_REFUSED << ")." ...
        this->ladrunoRestoreTrialFromCommitted();
        mSigma = sigTilde;          // leave the element holding what it equilibrated with
        ...
        return LADRUNO_MATERIAL_REFUSED;
    }
```
`Domain.cpp:2241-2245` (the return is not even assigned):
```cpp
    Element *elePtr;
    ElementIter &theElemIter = this->getElements();
    while ((elePtr = theElemIter()) != 0) {
      elePtr->commitState();
    }
```
`Domain.cpp:2309`: `Domain::commit()` ends `return 0;` — unconditionally. So
`AnalysisModel::commitDomain()`'s `if (myDomain->commit() < 0)` (`AnalysisModel.cpp:656`) can
never fire.

Note also what is *not* reset on that path: `ManzariDafalias::commitState()` is skipped, so
`mEpsilon_n` stays at step *n* while the element keeps feeding step *n+1* strains; and
`mImplexDtCommit`, `mImplexStepArmed = true`, `mImplexTrialDone = false` and `mImplexDEpsP`
(all updated only at `:1673-1682`, after the early return) are left stale.

**FAILURE SCENARIO** `nDMaterial LadrunoSANISAND 1 ... -implex -maxSubsteps 1000` on
`SSPbrick`/`stdBrick`, IntScheme 1 (the deck default and the only fully-qualified companion).
At some step the extrapolated global solve hands the commit-time return an increment that
needs > 1000 ModifiedEuler substeps. Output: 1 warning line per offending Gauss point, then
the analysis reports the step as converged and committed. From then on, that point's
`mEpsilon_n` is one step stale, so the *next* step's `dEps = mEpsilon - mEpsilon_n` spans two
steps, `f` is the stale frozen factor from the failed step (`mImplexStepArmed` was never
re-armed), and `sigma~` jumps by an increment the model never applied. The error compounds
every subsequent step. The final stress/settlement curve is wrong and nothing in the return
codes, the recorders or `analyze()`'s status says so. The comment at `:1642` ("refusing at
commit leaves the analysis believing a step advanced that did not") describes exactly what the
code then does.

**WHAT WOULD REFUTE IT** A code path showing an element or the Domain propagating
`Element::commitState()`'s return to `StaticAnalysis`/`TransientAnalysis`. (`Domain::commit()`
returning `0` on `Domain.cpp:2309` and the bare `elePtr->commitState();` at `:2244` would both
have to be wrong.) Or a demonstration that `mSubstepCapHitInME` cannot be set at commit under
`-implex -maxSubsteps N` without `-implexControl` — which contradicts the code's own comment
at `:1641`.

---

## 2. BLOCKER — `f` is NEVER the specified `dt_{n+1}/dt_n` on a monotonically-negative pseudo-clock; it silently degrades to `alpha` on every step. This is the fork's own campaign deck, and it is the deck the D2 "fix" (3c788778f) was measured on.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1329-1332`

**CLAIM** `ladrunoImplexArmStep()` gates the ratio on `mImplexDtCommit > 0.0`, not on
`mImplexDtCommit != 0.0`. Commit 3c788778f fixed the *refusal* to accept a consistently
negative clock but left the *factor* computation sign-blind, so every prescribed-settlement
deck (`LoadControl(-ds)`) runs with `f == implexAlpha == 1.0` forever, including through the
ladder's `ds *= 0.5` subdivisions where the correct `f` is 0.5.

**EVIDENCE** `LadrunoSANISAND.cpp:1325-1334`:
```cpp
    mImplexDt = dt;
    if (mImplexDt0 <= 0.0 && dt > 0.0)
        mImplexDt0 = dt;

    if (mImplexDtCommit > 0.0)
        mImplexFactor = (dt / mImplexDtCommit) * mImplexOpt.alpha;
    else
        mImplexFactor = mImplexOpt.alpha;
```
The D2 guard three hundred lines later explicitly *blesses* the negative clock
(`LadrunoSANISAND.cpp:1518-1521`): *"A consistently negative clock is NOT refused: a deck that
drives settlement downward has dt < 0 on every step and a positive, correct ratio."* The ratio
is never computed on that branch.

The deck: `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py:722`
`ops.integrator("LoadControl", -ds)`, with `ds *= 0.5` on failure at `:735` and a non-uniform
final step `ds = min(ds, smax - s_now)` at `:721`. `ops_Dt = currentTime - committedTime`
(`Domain.cpp:2081`), so `ops_Dt = -ds < 0` on every step.

The spec (`Ladruno_implementation/92_ladruno_sanisand_implex_adr.md:86`) is
`f = (dt_{n+1} / dt_n) * implexAlpha`.

**FAILURE SCENARIO** `sanisand_tau0_band.py` with `-implex`. Step *n* runs at `ds = 1 mm`
(`dt_n = -1e-3`); it fails, the ladder halves to `ds = 0.5 mm` (`dt_{n+1} = -5e-4`). Correct
`f = 0.5`; delivered `f = 1.0`. The extrapolated plastic strain `f*d_eps_p(n)` is **2x** what
the operator specifies, so `sigma~` is under-predicted by `Ce:(0.5*d_eps_p(n))` — a first-order
error in the step, exactly the quantity IMPL-EX's accuracy claim rests on. The BVP-gate
measurement in commit 70f6bf0e8 ("the ladder is GONE, and the leg reached its TARGET") was
therefore taken on `f ≡ 1.0`, i.e. on a different operator from the one the ADR specifies.

**WHAT WOULD REFUTE IT** Evidence that `ops_Dt` is positive on a `LoadControl(-ds)` deck, or a
test that asserts `implexDetail(5)` (`mImplexFactor`) equals 0.5 after a halving on a
negative-clock deck. The existing battery does not: `tests/test_ladruno_sanisand_implex.py:746`
only exercises a *sign change* (positive history, then `LoadControl(-1.0)`), which takes the
refusal branch and never reaches the ratio.

---

## 3. BLOCKER — `-implexControl`'s reduction floor is sign-blind in the same way, so on the same decks the control is either always-refuse (dead analysis, no message) or never-refuse (no control at all).

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1326-1327` and `:1615`

**CLAIM** `mImplexDt0` is only ever set from a *positive* `dt`, and the floor test compares the
*signed* `mImplexDt` against `reductionLimit * mImplexDt0`. On a negative-`ops_Dt` deck
`mImplexDt0` stays `0.0` forever, so the floor's `mImplexDt0 <= 0.0` short-circuit makes the
material refuse unboundedly; on a mixed-sign history where `mImplexDt0` did get set positive,
`mImplexDt < 0 < reductionLimit*mImplexDt0` makes the floor fire on the first refusal, so
`-implexControl` never refuses anything at all.

**EVIDENCE**
```cpp
:1326    if (mImplexDt0 <= 0.0 && dt > 0.0)
:1327        mImplexDt0 = dt;
```
```cpp
:1612    if (mImplexOpt.control) {
:1613        this->ladrunoImplexMeasureError(mSigma, sigImplicit, mEpsilon);
:1614        if (mImplexError > mImplexOpt.errorTol) {
:1615            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
:1616                return LADRUNO_MATERIAL_REFUSED;
```
The parser documents `reductionLimit` as *"the floor below which the material stops refusing,
as a fraction of the FIRST dt seen"* (`:311-313`) — a statement about magnitude, implemented as
a signed comparison.

**FAILURE SCENARIO** `sanisand_tau0_band.py` + `-implex -implexControl 0.05 0.01`. Every step
has `dt < 0`, so `mImplexDt0 == 0.0`. Any step whose error exceeds 0.05 is refused; the ladder
halves; the error is still > 0.05 at the low-confinement corner (P0 measured 1.26 at
p0 = 5 kPa); the floor that is supposed to stop the refusals at `0.01*dt0` never engages, so
the leg dies at `DS_MIN` having advanced nothing — the identical symptom 3c788778f was written
to cure (4896 refusals of the first push step), reproduced through a second sign-blind guard
that the fix did not touch.

**WHAT WOULD REFUTE IT** A path that sets `mImplexDt0` from `|dt|`, or evidence that
`-implexControl` is only ever used with `-implexDt user`/`strain` (both of which are
non-negative by parser check, `:352`, `:1311`) — but the default source is `pseudo` and nothing
refuses the combination.

---

## 4. MAJOR — every trial-time refusal is a *silent wrong answer* on any element that does not propagate `-33086` (SSPbrick, stdBrick, BbarBrick, the quads): the material hands the element the COMMITTED stress and still runs the companion at commit.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1453`, `:1523-1526`, `:1547-1552`;
`SRC/element/UWelements/SSPbrick.cpp:445-447`

**CLAIM** `mImplexTrialDone = true` is set at the top of `ladrunoImplexTrial()`, before any
refusal, and the D2 and companion-failure refusals restore the trial state to the committed
one. An element that ignores the return code therefore equilibrates on `sigma_n` (a stress
independent of the trial strain, with tangent `Ce`) and then commits a *companion* return for
that step as if nothing happened. The echo's caveat ("only an element that PROPAGATES it can
act on -- today LadrunoBrick", `:1226-1229`) understates this: the refusal is not inert, it
substitutes a wrong stress.

**EVIDENCE**
```cpp
:1450 LadrunoSANISAND::ladrunoImplexTrial(void)
:1452     // commitState() owes a companion return after this.
:1453     mImplexTrialDone = true;
```
```cpp
:1523        this->ladrunoRestoreTrialFromCommitted();
:1524        double Kr = 0.0, Gr = 0.0;
:1525        this->ladrunoImplexFreezeTangent(Kr, Gr);
:1526        return LADRUNO_MATERIAL_REFUSED;
```
`ladrunoRestoreTrialFromCommitted` (`:1347`) begins `mSigma = mSigma_n;`.
`SSPbrick.cpp:443-447` — the element the UW SANISAND family ships with:
```cpp
	Vector strain(6);
	strain = Bnot*u;
	theMaterial->setTrialStrain(strain);

	return 0;
```
Only `LadrunoBrick.cpp` tests for the sentinel (`grep -rn LADRUNO_MATERIAL_REFUSED SRC/`).

**FAILURE SCENARIO** `element SSPbrick ... 1` with `-implex -implexControl 0.05 0.01` on a
DisplacementControl leg (D2's stated target, which the guard *cannot* catch pre-limit-point,
`:1502-1505`). A step reverses sign → the material refuses → SSPbrick returns 0 → Newton sees a
constant internal force `∫B^T sigma_n` with a non-zero stiffness `Ce`, converges or stalls on a
meaningless residual → `commitState()` still sees `mImplexTrialDone == true` and commits the
companion. The reported equilibrium is an equilibrium of the *previous step's* stress field.
No non-zero return code anywhere.

**WHAT WOULD REFUTE IT** A refusal on `setTrialStrain` that leaves `mSigma` on a defensible
value (the `-implexControl` W7 site at `:1615` does — it leaves `sigma~` — which is precisely
the inconsistency), or a guard that refuses `-implex` at construction on elements that do not
propagate the sentinel. Neither exists.

---

## 5. MAJOR — D3's "the companion MUST be able to fail" rule is enforced for IntScheme 1 and merely *warned* for IntScheme 2, even though the code's own text says scheme 2's ladder ends in an uncapped ModifiedEuler. `-implex IntScheme 2` with no `-maxSubsteps` is accepted.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1182-1190` vs `:1192-1211`

**CLAIM** The same argument that makes `-maxSubsteps > 0` a *hard refusal* on scheme 1 applies
verbatim to scheme 2 by the code's own admission, but scheme 2 gets a warning and proceeds.
The consequence is worse than slowness: `mSubstepCapHitInME` is the *only* failure detector
`-implexControl` has (`:1543`), and it can never be set when the cap is 0, so
`-implex -implexControl` on scheme 2 with no cap silently cannot detect a companion failure at
all — it degrades to an error-magnitude check only.

**EVIDENCE**
```cpp
:1182    if (s == 1 && mMaxSubsteps <= 0) {
:1183        opserr << ... ": -implex on IntScheme 1 REQUIRES -maxSubsteps > 0 ..."
:1189        return -1;
:1190    }
```
```cpp
:1205        if (mMaxSubsteps <= 0)
:1206            opserr << ... ": -implex with IntScheme 2 and -maxSubsteps 0. Scheme 2's own"
:1207                      " retry ladder ends in explicit_integrator, so the companion can"
:1208                      " still spend an unbounded number of ModifiedEuler substeps in one"
:1209                      " commit. -maxSubsteps is strongly recommended." << endln;
```
(no `return -1`). The setLadrunoImplexOptions echo itself notes P0 measured "58-74 % of
BackwardEuler_CPPM's calls on the low-confinement corner path taking its low-p branch ...
integrates by explicit_integrator (ModifiedEuler)" (`:1195-1199`) — i.e. on the campaign's own
corner, scheme 2 *is* ModifiedEuler most of the time.

**FAILURE SCENARIO** `-implex -implexControl 0.05 0.01` on `IntScheme 2` with no
`-maxSubsteps`: accepted with a "strongly recommended" note; the commit-time companion can
run ~1e6 substeps (ADR-90 GATE U measured 34-minute single updates) with no cap and no
detector; the user believes `-implexControl` is guarding them.

**WHAT WOULD REFUTE IT** Evidence that `mSubstepCapHitInME` can be set with `mMaxSubsteps == 0`
(`ladrunoUpdateStatus`'s own comment at `:1838` says the opposite: *"at the default cap 0
nothing branches on either value"*), or a scheme-2 path that reports failure some other way.

---

## 6. MAJOR (diagnostic) — `avgImplexError` destroys the process accumulator on every read, and both read sites are per-integration-point, so a recorder over more than one Gauss point reports ~0 for all but the first.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1074-1082`, `:1767`, `:1916`

**CLAIM** `takeAverageError()` zeroes `sumError`/`count`. `getResponse` calls it once per
Gauss point per record, and `setParameter` calls it once at parameter-creation time. The number
a deck reads is therefore a function of how many points asked and in what order, not of the
analysis.

**EVIDENCE**
```cpp
:1074    double takeAverageError(void)
:1075    {
:1076        double a = (count > 0) ? (sumError / (double)count) : 0.0;
:1077        sumError = 0.0;
:1078        count    = 0;
:1079        return a;
:1080    }
```
```cpp
:1916        out1(0) = LadrunoImplexGlobals::instance().takeAverageError();
```
```cpp
:1767            param.setValue(LadrunoImplexGlobals::instance().takeAverageError());
```
`accumulate()` is also non-atomic on a process-wide singleton (`:1066-1072`), which the fork's
own ADR-75b Lane-3 threaded-assembly work makes a live concern.

**FAILURE SCENARIO** `recorder Element -ele all -material 1 avgImplexError` on an 8-point brick
mesh: point 1 of element 1 returns the true average since the last record; the remaining
~400k reads return exactly `0.0`, because `count` was reset. The recorded column is a field of
zeros with one non-zero entry, and it looks like "IMPL-EX error is negligible". `maxError` is
never reset at all, so `implexError` (the `Parameter` at `:1762`) is a since-process-start
maximum that no `revertToStart` clears.

**WHAT WOULD REFUTE IT** A caller that guarantees exactly one read per record. There is none —
`setResponse` (`:1908`) is reached once per requested integration point.

---

## 7. MAJOR — the two refusal warnings are UNTHROTTLED per Gauss point per iteration, in the one file that documents an 83 MB stderr incident from exactly this shape and installs a process-wide budget for it.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1512-1522` and `:1647-1656`; contrast
`SRC/material/nD/UWmaterials/ManzariDafalias.cpp:1444-1466`

**CLAIM** Both `opserr` blocks sit inside a Gauss-point update inside a Newton iteration
inside a load step, with one material instance per integration point — the exact three-deep
nesting the ADR-86 PR-1 throttle was written for — and neither has a counter.

**EVIDENCE** The D2 refusal (`:1512`) prints ~11 lines of prose with no guard. The precedent in
the same family (`ManzariDafalias.cpp:1451-1466`):
```cpp
        static int ladrunoClampWarnCount = 0;                               // Ladruno
        if (ladrunoClampWarnCount < 10) {                                   // Ladruno
```
with the rationale spelled out at `:1444-1453`: *"a 50k-element mesh holds ~400k instances, and
a per-instance echo of exactly this shape measured 83 MB of stderr during ADR-86 PR-1"*.

**FAILURE SCENARIO** Measured, by this PR's own history: commit 3c788778f reports *"The first
push step was refused 4896 times"*. Each of those refusals executed the `:1512` block. On the
50k-element mesh the same warning is a multi-GB stderr stream, and on Windows the console write
is synchronous — the "seizure" the whole ADR-86b cap exists to prevent, reintroduced as I/O.

**WHAT WOULD REFUTE IT** A throttle elsewhere on the path, or evidence that `opserr` itself
rate-limits. Neither exists.

---

## 8. MAJOR — a `-implexControl` refusal is completely silent (no warning, no state restore, no re-arm), the opposite of the other two refusal sites.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1612-1618`

**CLAIM** The W7 refusal returns the sentinel with no `opserr` at all, no
`ladrunoRestoreTrialFromCommitted()`, and no `mImplexStepArmed` reset — while the D2 site
(`:1512-1526`) and the commit site (`:1647-1660`) all three print, restore and refreeze. So the
one refusal a user is most likely to hit produces zero diagnostic output, and combined with
finding 3 a deck can die at `DS_MIN` with nothing in the log naming the cause.

**EVIDENCE**
```cpp
:1614        if (mImplexError > mImplexOpt.errorTol) {
:1615            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
:1616                return LADRUNO_MATERIAL_REFUSED;
:1617        }
```
Compare `:1523-1526` and `:1657-1660`, both of which restore and refreeze before returning the
same value.

**FAILURE SCENARIO** `-implex -implexControl 0.05 0.01`, a leg that fails at `DS_MIN`. The log
contains no line explaining why any step was rejected. The only evidence is
`implexDetail`/`implexError`, which the deck has to have asked for in advance. A user
debugging this cannot distinguish an `-implexControl` refusal from a solver failure.

**WHAT WOULD REFUTE IT** A warning emitted on that path (there is none between `:1612` and
`:1618`), or a documented decision that this site is deliberately quiet — the surrounding
comment (`:1615-1616`) discusses only the floor, not the silence.

---

## 9. MINOR/MAJOR — any zero-`dt` step silently destroys the extrapolation history: `mImplexDtCommit` is set to 0 and `d_eps_p(n)` collapses to ~0, so the *next* step reverts to `f = alpha` and a purely elastic predictor.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1479-1489`, `:1673-1682`, `:1329-1332`

**CLAIM** A `LoadControl 0.0` hold, or any step in which a given Gauss point's strain increment
happens to be exactly zero, still sets `mImplexTrialDone = true` (`:1453`) and therefore still
commits: `mImplexDtCommit = mImplexDt = 0.0` and `mImplexDEpsP` = the companion's plastic
increment over a zero increment ≈ 0. Both the ratio and the extrapolated plastic strain for the
*following* step are then wrong, and nothing reports it.

**EVIDENCE**
```cpp
:1483            this->ladrunoImplexArmStep();       // clears mImplexStepArmed
:1484        } else {
:1485            mImplexDt     = 0.0;
:1486            mImplexFactor = 0.0;                // no strain advanced, no plastic flow predicted
:1487        }
```
```cpp
:1679    mImplexDtCommit  = mImplexDt;
```
then `:1329` `if (mImplexDtCommit > 0.0)` is false → `f = alpha`.

**FAILURE SCENARIO** The staged-gravity idiom this family uses everywhere ends with a
`LoadControl 0.0` hold *after* `updateMaterialStage 1` (so IMPL-EX is active). That hold
commits `dt_n = 0` and `d_eps_p(n) = 0`. The first real plastic step after the hold then runs
with `f = 1.0` on `d_eps_p = 0`, i.e. a pure elastic predictor — indistinguishable from
`-implexAlpha 0`. The IMPL-EX operator the deck asked for does not run on that step and
`implexError` reports the (small) error of an operator nobody requested.

**WHAT WOULD REFUTE IT** A guard that skips the `mImplexDtCommit`/`mImplexDEpsP` update when
`mImplexDt == 0`. The commit path at `:1673-1682` has none.

---

## 10. MINOR — `revertToLastCommit()` is gated on `mImplexOpt.enabled`, not on `ladrunoImplexActive()`, so the ADR's stage-0 inertness claim ("the extrapolated path is simply not reachable") is false for the revert path; it also writes `mEpsilon`, which vanilla never does.

**file:line** `SRC/material/nD/LadrunoSANISAND.cpp:1713-1720`; claim at `SRC/material/nD/LadrunoSANISAND.cpp:1249-1257`

**CLAIM** During elastic stage 0 (gravity, the `LoadControl 0.0` hold), a `-implex` deck takes a
materially different `revertToLastCommit` from the same deck without `-implex`: it overwrites
the element's trial strain with the committed one and rewrites all three tangent slots with
`Ce(mSigma_n, mVoidRatio)` computed under `mElastFlag == 0`. The stage-0 byte-identity claim is
therefore established by inspection of `ladrunoImplexActive()` only, and one of the three
overridden entry points does not use it.

**EVIDENCE**
```cpp
:1713    if (!mImplexOpt.enabled)
:1714        return res;
...
:1720    mEpsilon = mEpsilon_n;
:1721    this->ladrunoRestoreTrialFromCommitted();
```
against the claim at `:1252-1257`: *"While it is 0 the base takes its elastic_integrator branch
and -implex is inert, so gravity and the `LoadControl 0.0` re-equilibration are bit-identical
with the flag on and off. There is no branch to get wrong: the extrapolated path is simply not
reachable."*

`mEpsilon = mEpsilon_n` has no vanilla counterpart — `ManzariDafalias::revertToLastCommit()` is
an empty `return 0`, and no ManzariDafalias path ever writes the trial strain from the
committed one. Any caller that invokes `material->revertToLastCommit()` and then reads
`getStrain()` (or a recorder that samples between the revert and the Domain's trailing
`update()`) sees a different value than on the `-implex`-off build.

**WHAT WOULD REFUTE IT** A measurement showing stage-0 gravity + a failed step + revert is
bit-identical with and without `-implex`. The battery does not contain one (its gravity checks
do not force a revert during stage 0).

---

### Two things I checked and could NOT break (recorded so the blue team does not have to)

- The `p_min` clamp threshold. The comment at `:1568-1571` claims parity with the base's
  ModifiedEuler pre-loop guard; that guard is `p = one3*trace + m_Presidual; if (p < m_Pmin +
  m_Presidual)` (`ManzariDafalias.cpp:1438-1439`), which is algebraically `one3*trace < m_Pmin`
  — exactly the implex test at `:1584`. The claim holds.
- `mVoidRatio` re-derivation in `ladrunoRestoreTrialFromCommitted` (`:1357`) matches
  `ManzariDafalias::commitState()`'s own expression (`ManzariDafalias.cpp:510`) with
  `mEpsilon == mEpsilon_n` after commit. Exact.
