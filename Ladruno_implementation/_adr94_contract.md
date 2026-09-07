# ADR-94 R5 — host-element contract for ASDPlasticMaterial3D failure returns

Build `52314165a`. Plan: `94_asdplastic_review_plan.md` sec. 5, R5. Tests:
`tests/test_adr94_contract.py` (Zone-A, `t0m`, 8/8 pass, ~1.6 s). Refs:
`SRC/material/LadrunoMaterialStatus.h`, `_adr94_hlist_R1B.md` (H4),
`_adr94_hlist_R1A.md` (H7, H8).

## 1. Integrator × failure mode → return code

`ASDPlasticMaterial3D.h`, lines 1386-3950. `LADRUNO_MATERIAL_REFUSED = -33086`
(`SRC/material/LadrunoMaterialStatus.h`) is returned from exactly **two**
sites, both inside `Backward_Euler`; every other failure path returns a bare
`-1`.

| Integrator | Failure mode | Site(s) | Code |
|---|---|---|---|
| Forward_Euler | singular local tangent / NaN guard | 1467, 1545 | bare -1 |
| Forward_Euler_Subincrement | singular local tangent / NaN guard (mirrors FE) | 1649, 1714 | bare -1 |
| Backward_Euler | `special_return` fallback-to-vertex under `strict_convergence` | 2198 | **sentinel** |
| Backward_Euler | singular local tangent / NaN guard (trial stress) | 2284, 2326 | bare -1 |
| Backward_Euler | scalar-Newton exhaustion under `strict_convergence` | 2350 | **sentinel** |
| Backward_Euler_LineSearch | split-loop exhaustion (all halvings failed) | 2588 | bare -1 |
| Runge_Kutta_45_Error_Control_old | substep niter exceeded | 2902 | bare -1 |
| Modified_Euler_Error_Control | NaN+`dT<dT_min` / niter exceeded / final NaN | 3223, 3258, 3377 | bare -1 |
| Runge_Kutta_45_Error_Control | NaN+`dT<dT_min` / final NaN | 3744, 3931 | bare -1 |

`strict_convergence` only reaches `Backward_Euler` (`be_strict`, 2081/2086/
2184/2338) — every other integrator's failure paths above fire regardless of
the flag (this is H5/H7/H8's own finding, reconfirmed here structurally).

## 2. Host propagation (the R5 question)

| Host | Mechanism | Sentinel (-33086) | Bare -1 (everything else) |
|---|---|---|---|
| `TenNodeTetrahedron` | `success += ...setTrialStrain(...)`, returns `success` | propagated | **propagated** |
| `stdBrick` (`Brick.cpp`) | `success = ...setTrialStrain(...)`; `update()` unconditionally `return 0` | swallowed | swallowed |
| `LadrunoBrick` | 6 call sites (`updateHypo` SSP-centroid + per-GP loop, `formEAStrue` condensed + full loops), each `if (... == LADRUNO_MATERIAL_REFUSED)` | propagated (element prints its own "...REFUSED the trial strain..." warning, returns -1) | **swallowed at the element level** — no warning, call treated as success by `LadrunoBrick` itself |

Pinned structurally (regex over the three hosts' source, all 3 pass).

**Runtime pin (`Backward_Euler_LineSearch` exhaustion, bare -1, site 2588 —
the one mode this review could force reliably; see §3):** on the ADR-84 MC
tet leg, `TenNodeTetrahedron` shows H8's 2-good/18-fail pattern. On a
`LadrunoBrick` cube driven the same way, `analyze()` STILL fails every step —
but `LadrunoBrick`'s own refusal message never appears (child-process
capture). **The failure is only visible because the uncorrected stress
happens to unbalance the global residual enough to blow the Newton budget —
an accident of this rig, not a guarantee**: H5/H7 show the opposite accident
(a bad local state satisfying the global residual, rc=0). The structural
swallow is real regardless of which way any rig's global Newton falls.

## 3. Method note — FE-NaN / ME-max-iter / RK45-dT_min not runtime-isolated

Runtime pins for the FE NaN-guard, ME max-iter, and RK45 `dT_min` paths were
attempted and **not reproduced** in budget. VonMises converges in ~1-2
iterations at any step size (never exhausts/NaNs). MohrCoulomb needs far
more (peak ~91 at default niter 100) — even niter 20-30 already refuses on
the FIRST step, so there is no stable small/big-step window for a multi-step
reproducer (a niter=10, 1e-4-vs-1e-3 window exists for step 1 only, closing
once any plastic strain commits). Forcing a real NaN needs an ill-
conditioned Newton (`H_iso` steeper than `2G`, H7) at overflow scale, but
risks the H7 elastic-predictor commit (rc=0) instead — observed on VonMises,
`H_iso=-120000`, utop to -500. These sites stay confirmed-by-reading only
(§1); H15's aborted `StiffSoilShear_YF` attempt (an unrelated NaN on BE's
first step) is the best lead for a future reproducer.

## 4. Revert semantics — measured consequences of H4

**`ops.reset()` leaves the material inconsistent with geometry (CONFIRMED,
new).** `Domain::revertToStart()` zeros nodal displacements correctly;
`ASDPlasticMaterial3D::revertToStart()` is a no-op ("not implemented", -1,
ignored by `OPS_resetModel()`). Querying stress right after `reset()`
exercises `TenNodeTetrahedron`'s self-heal (recomputes strain from the now-
zero nodal displacement) against the STALE material Commit state: the result
is neither ~zero (a real reset) nor the pre-reset value (untouched) — a
third, inconsistent number.

**Cutback after a forced global-Newton failure — measured, inconclusive.** A
step failed via an impossible `NormDispIncr` budget leaves a dirty trial
stress `revertToLastCommit()`'s no-op body never clears. Retrying the
identical step plus the remaining steps reaches a final stress differing
from a never-failed reference by ~6e-9 relative — the level of the
`NormDispIncr 1e-8` tolerance itself, so this probe cannot separate H4's
broken revert from ordinary Newton-truncation noise (test asserts
`1e-9 < diff/scale < 1e-6` so either extreme flags for re-triage).

## 5. `LadrunoBeginAugment` / `LadrunoEndAugment`

Not contact-gated: both unconditionally call
`Domain::setContactAugmenting(true/false)` — no contact-handler check,
idempotent, always `return 0`. Legal (and a functional no-op for the
material's own integrator/revert path) on a plain, contact-free model; the
flag only changes whether `Domain::commit()` fires recorders / bumps
`commitTag`. Not measured (budget): whether repeated held-load re-commits
multiply H12's per-tag `cout` diagnostic — flagged as a follow-up.

## 6. Summary for the R6 verdict

`LadrunoBrick` is the ONLY host that ever distinguishes a material failure,
and only at the two `Backward_Euler` sentinel sites; every other
integrator's bare -1 is treated as success on all three hosts (`stdBrick`
always; `LadrunoBrick` for anything but the sentinel). Whether a swallowed
failure is ever OBSERVED depends on whether the bad state happens to
unbalance the global residual (H5/H7: silent rc=0; here: accidental rc!=0) —
"my model would have failed anyway" is not a safe inference. Fix direction:
widen `LADRUNO_MATERIAL_REFUSED` to every site in §1, rather than loosening
`LadrunoBrick`'s checks to `< 0` (loses the "commit guaranteed unchanged"
guarantee unless every bare-`-1` site is re-audited).
