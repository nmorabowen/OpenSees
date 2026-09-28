# WP-110 / F15a -- GetElastoPlasticTangent measurement probe: results (revised)

> **Superseded as a test (WP-110 phase 1):** the probe file `tests/test_manzari_ep_tangent_probe.py` became the gate
> `tests/test_manzari_ep_tangent_gate.py`, which reads the ENGINE's tangent via `eleResponse(ele, 'tangent')`
> instead of a numpy transcription. The `get_ep_tangent_workbench` oracle cited below lives on there as
> `workbench_tangent`. This note stays as the measurement record.

**Status: CONFIRMED, and a SECOND, independent defect found.**
`ManzariDafalias::GetElastoPlasticTangent` (ManzariDafalias.cpp:5110) has two
separate Voigt covariant/contravariant mistakes, both in the plastic
correction term `aCep = aC - (Macauley/temp3) * Dyadic2_2(temp1, temp2)`:

1. **Numerator (originally reported, still confirmed):** `R` is converted to
   covariant TWICE — `R = ToCovariant(temp0)` (:5131) then
   `temp1 = DoubleDot4_2(aC, ToCovariant(R))` (:5133, wraps the *already*
   covariant `R` again). This doubles `temp1`'s shear (Voigt 3,4,5) entries a
   second time, so it can only ever corrupt aCep's **shear rows** — measured
   at an exact **2.000x** ratio on every shear-row/normal-column entry.

2. **Denominator (NEW, found only after the coordinator's Workbench
   cross-check flagged a residual):**
   `temp3 = DoubleDot2_2_Contr(temp2, R) + Kp;` (ManzariDafalias.cpp:5139)
   uses the CONTRAVARIANT double-dot (`DoubleDot2_2_Contr`, which re-doubles
   Voigt indices 3-5) on a mismatched pair: `temp2` (a stress-like vector from
   contracting a covariant strain direction through `aC`) against `R`
   (covariant, already shear-doubled). This is a different Voigt-convention
   slip from defect 1, and it survives even after defect 1 is fixed. Measured:
   `temp3` (engine, as-written -- unaffected by defect 1) is a **constant
   ~1.0757x** too large relative to the Workbench's properly-contracted
   `denom = Kp + Q:Ce:R` (`ddot(QCe, R)`, a true Frobenius double-dot), at
   BOTH states tested (1.07577 at p'~20, 1.07568 at p'~214 -- consistent to 4
   significant figures, i.e. state-independent in this sample). Since `temp3`
   is a scalar dividing the WHOLE rank-1 correction, this rescales *every*
   entry of the correction term uniformly by ~7.6%, regardless of defect 1.

**Numerator-only fix ("fixed" in the original note) is INSUFFICIENT.** It
still uses the wrong `temp3`, so it is still ~6-8% off the true continuum
tangent. Both defects must be fixed to match FD.

## What changed since the first note, and why

The first note's ~5-7% "fixed vs FD" residual was reported at IntScheme=2 with
a *central* finite difference and, for the "HIGH" case, at an unrealistic
p'=13218 kPa. Two problems were found on review:

- **Harness artifact (confirmed, real):** `integrate()`'s reversal check
  (ManzariDafalias.cpp ~1005-1013: `if (alpha_n - alpha_in_n):Ce:d_eps < 0,
  reset alpha_in := alpha_n`) fires or not depending on the SIGN of the probed
  strain increment relative to the ongoing loading direction. A symmetric
  central difference `(F(+h)-F(-h))/(2h)` therefore differences two states
  with DIFFERENT `alpha_in` (confirmed by direct read-back: at every one of
  the 6 probe directions, `alpha_in(+h) != alpha_in(-h)`), straddling a real
  kink in the response rather than sampling one smooth branch. This
  contaminates the off-diagonal FD entries badly (LOW p', central FD at
  h=1e-8: row3 = `[2327, -1165, 2372, 16417, 0, 1.5]` -- compare to any of the
  three analytic candidates below; none matches within 40%).
- **Unrealistic test state:** p'=13218 kPa is far outside where these Toyoura
  parameters were tuned; the original note's per-block error table was partly
  measuring that, not just the bug.

**Fix applied here:** (a) realistic states, p'~20 and ~214 kPa, eta~1.3, WITH
genuine shear demand (`n`'s Voigt-shear component ~0.21, so the numerator bug
is actually visible -- a state with zero shear, e.g. pure coaxial triaxial
compression, makes `R`'s shear component zero and the numerator bug
disappears identically, which is what the *first* realistic-state attempt
without shear accidentally measured); (b) a *one-sided* finite difference
whose sign is chosen, per probe direction, to keep `alpha_in` identical to the
committed base state (no reversal reset) -- see `check_alphain.py` logic
folded into `final_compare.py`.

## Method addition: the Workbench formula, transcribed as a third candidate

Per the coordinator's instruction, `C4-sanisand/sanisand_cep.py` (read-only,
copied to scratch, not run) was diffed term by term against
`GetElastoPlasticTangent`. It computes `Dep6` with proper 3x3 tensor algebra
throughout (no Voigt ToCovariant/DoubleDot2_2_Contr detours):

```
CeR  = 2G*dev(R) + K*tr(R)*I      # Ce:R
QCe  = 2G*dev(Q) + K*tr(Q)*I      # Q:Ce,  Q = n - (n:r)/3 * I
denom = Kp + ddot(QCe, R)         # Q:Ce:R, a true Frobenius double-dot
Dep6  = Ce6 - outer(CeR_vec, QCe_vec) / denom
```

`Kp`, `R`, `B`, `C`, `D`, `n`, `d`, `b`, `h`, the state evaluated at
(mSigma, mAlpha, mFabric, mAlpha_in, e) and the loading test are otherwise
IDENTICAL to `GetStateDependent`/`GetElastoPlasticTangent` -- only the Voigt
bookkeeping in the final assembly differs. This "workbench" formula was added
as a third numpy candidate (`get_ep_tangent_workbench` in
`tests/test_manzari_ep_tangent_probe.py`) alongside "buggy" (as-written) and
"fixed" (numerator-only correction).

## Results table (row 3 = dSigma_xy / dEps_*, the shear-shear diagnostic row)

Two states, p'~20 and p'~214 kPa, eta~1.3, WITH shear demand
(n_shear = dep['n'][3] ~ 0.214 at both):

### LOW: p'=20.03 kPa, q=25.90, eta=1.293, dGamma=1.60e-6

| source | col0 (xx) | col1 (yy) | col2 (zz) | col3 (xy, diag) |
|---|---:|---:|---:|---:|
| buggy (as-written) | 8643.2 | -4339.5 | 8833.8 | 14663.9 |
| fixed (numerator only) | 4321.6 | -2169.8 | 4416.9 | 15866.0 |
| **workbench (both fixed)** | **4649.0** | **-2334.2** | **4751.6** | **15774.9** |
| FD, central diff, h=1e-8 (contaminated) | 2327.2 | -1165.3 | 2372.0 | 16417.3 |
| **FD, one-sided reversal-safe, h=1e-8** | **4654.4** | **-2330.7** | **4744.0** | **15766.5** |

One-sided FD vs workbench: **0.12%, 0.11%, 0.16%, 0.05%** relative error.
One-sided FD vs fixed (numerator-only): **7.1%, 0.9%, 7.4%, 0.63%**.
One-sided FD vs buggy: **85.7%, 86.1%, 86.2%, 6.99%**.

### HIGH: p'=213.56 kPa, q=292.9, eta=1.372, dGamma=7.41e-5

| source | col0 (xx) | col1 (yy) | col2 (zz) | col3 (xy, diag) |
|---|---:|---:|---:|---:|
| buggy (as-written) | 29755.6 | -12585.8 | 30336.4 | 47887.5 |
| fixed (numerator only) | 14877.8 | -6292.9 | 15168.2 | 51808.1 |
| **workbench (both fixed)** | **16003.7** | **-6769.1** | **16316.1** | **51511.4** |
| FD, central diff, h=1e-8 (contaminated) | 8004.7 | -3383.6 | 8153.9 | 53615.2 |
| **FD, one-sided reversal-safe, h=1e-8** | **16009.4** | **-6767.1** | **16307.7** | **51501.8** |

One-sided FD vs workbench: **0.04%, 0.03%, 0.05%, 0.02%**.
One-sided FD vs fixed: **7.1%, 0.02%, 7.0%, 0.02%**.
One-sided FD vs buggy: **85.9%, 85.9%, 86.0%, 7.06%**.

**IntScheme cross-check:** repeating the one-sided-FD-vs-workbench comparison
with IntScheme=1 (ModifiedEuler) instead of 2 (BackwardEuler_CPPM) gives the
same result to within 0.1% at both states -- the continuum-tangent limit does
not depend on which scheme produced the probed states, as expected (FD with
h->0 samples the local rate response, not the discretization).

## Verdict

**CONFIRMED, with a correction to the earlier verdict's scope: two defects,
not one.**

1. Numerator double `ToCovariant` on `R` (:5131, :5133) -- exact 2x error on
   aCep's shear rows. (Original finding, still holds.)
2. Denominator Voigt-mismatch (:5139, `DoubleDot2_2_Contr(temp2, R)`) --
   ~7.6% systematic overestimate of `temp3` relative to the true `Q:Ce:R`,
   uniformly weakening the entire plastic correction. (New.)

The apparent "fixed is still 5-7% off FD" was **both** real (defect 2, ~6-8%,
confirmed present and independent of defect 1) **and** partly a harness
artifact (a central-difference FD straddling `integrate()`'s alpha_in
reversal branch, which by itself produced 40-100%+ noise in the off-diagonal
FD columns and had nothing to do with either defect). With the harness fixed
(one-sided, reversal-consistent FD) and both defects corrected in the numpy
oracle (the "workbench" formula), FD matches to **0.02%-0.17%** at both
p'~20 and p'~214 kPa, eta~1.3, with genuine shear content -- as good a
confirmation as the two independent tangent derivations (this probe's
transcription and the Workbench's `sanisand_cep.py`) can give without a
third, independently-authored implementation.

## F15(d) -- BVP replay: SANISAND self-weight strip footing, `-flipAlphaIn init`

**Question:** the Workbench ran a SANISAND self-weight strip footing (B/8
quad, implicit, `-flipAlphaIn init`) and reported TanType 1/2 parting from
TanType 0 near s/B ~ 0.0075 (+7 % at 0.009, +13 % at 0.010, +22 % at 0.0115),
not brought back by a tighter tolerance, and not reproducible run-to-run at
TanType 1/2 while TanType 0 was. Does that survive the F15 tangent fix
(commit `dee04dbe3`, this WP's own P0/P1)?

**Deck:** no fork BVP deck matched "B/8 quad" or "-flipAlphaIn" by name, so
this replay reuses the fork's own SANISAND self-weight strip-footing driver,
`Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py` (the ADR-90/WP-A2
tau=0 collapse-band study, also the base of WP-106's `-pRe` BVP gate) --
`LadrunoBrick -formulation bbar`, plane strain, rough rigid footing, K0
self-weight staged before the push, the R3-graded mesh, `NormUnbalance`
ladder Newton -> NewtonLineSearch -> KrylovNewton. **The driver had no
`-flipAlphaIn` passthrough**, so one was added: `run_leg(..., flip_alpha_in=)`
and a `--flipAlphaIn init|vanilla` CLI flag, emitted only when non-default
`vanilla` so every pre-existing leg's material command stays byte-identical
(same rule as `-pRe`/`-implexFactor`). Shrunk for wall-clock: single leg
`h1.0_e0.6944` (h0 = 1.0 m, Gorini's calibrated e_init = 0.6944, the coarsest
of the deck's 3x2 grid), `--sfrac 0.012` (vs the deck's own 0.25 default), no
surcharge, `IntScheme 1` (ModifiedEuler, the deck default). No C++ changes,
no rebuild.

**Binaries**, both confirmed via `ops.ladrunoBuild()`:
post-fix = this worktree's `dist\bin` (`dee04dbe39117d16835f7caa52a6b71aab72f1ba`,
F15 fix committed); pre-fix = the installed
`C:\Program Files\Ladruno\OpenSees\bin` (`48c0e99bc8e28bbb4fdf965015f285f01e16e90d`).

**A genuine seizure, not a flake:** the first TanType 2 attempt (post-fix,
uncapped substeps, the deck's own default) hung completely -- zero progress
past its s/B = 0.005 checkpoint for over 20 minutes of active CPU (a
threaded Pardiso solve, confirmed still `Responding` and burning CPU, not
deadlocked) before being killed. This is the deck's own documented GATE-U
failure mode ("every uncapped leg of this deck seizing inside ModifiedEuler,
so uncapped ... would simply spend their wall budget and measure the
budget" -- the module docstring). Both TanType 2 legs below therefore use
`--maxsubsteps 20000`, WP-106's own precedent for this exact deck; TanType 0
never needed it (uncapped throughout, no seizure, `nsub 0` in the JSON either
way).

**Legs run** (all `--tantype`/`--flipAlphaIn init`, `h1.0_e0.6944`, one
foreground `python3.12 -u` call each, ~3.5-5 min wall):

| leg | binary | TanType | maxsubsteps | wall_s | steps | failed rungs | relaxed |
|---|---|---:|---:|---:|---:|---:|---:|
| postfix_tan0 | dee04dbe3 (post-fix) | 0 | 0 (uncapped) | 214 | 46 | 10 | 0 |
| postfix_tan2_run1 | dee04dbe3 (post-fix) | 2 | 20000 | 252 | 46 | 14 | 3 |
| postfix_tan2_run2 | dee04dbe3 (post-fix) | 2 | 20000 | 252 | 46 | 14 | 3 |
| prefix_tan0 | 48c0e99bc (pre-fix) | 0 | 0 (uncapped) | 216 | 46 | 10 | 0 |
| prefix_tan2 | 48c0e99bc (pre-fix) | 2 | 20000 | 297 | 46 | 18 | 2 |

Artifacts under `Ladruno_files/testbed/hypo_bearing/wp110_f15d/<leg>/`
(`a2_h1.0_e0.6944_curve.csv`, `_field.csv`, `.json`, engine logs) plus the
five `*_stdout.log` driver transcripts.

### q (kPa) at s/B, linearly interpolated from each leg's own step sequence

| leg | 0.005 | 0.0075 | 0.009 | 0.010 | 0.0115 | end (0.012) |
|---|---:|---:|---:|---:|---:|---:|
| postfix_tan0 | 205.061 | 294.085 | 346.670 | 381.095 | 432.124 | 449.454 |
| postfix_tan2_run1 | 204.819 | 295.173 | 346.177 | 380.619 | 431.843 | 448.887 |
| postfix_tan2_run2 | 204.819 | 295.173 | 346.177 | 380.619 | 431.843 | 448.887 |
| prefix_tan0 | 205.061 | 294.085 | 346.670 | 381.095 | 432.124 | 449.454 |
| prefix_tan2 | 204.824 | 294.130 | 346.392 | 380.905 | 433.842 | 450.462 |

TanType 2 vs TanType 0 at the same binary (% difference):

| s/B | post-fix (2 vs 0) | pre-fix (2 vs 0) |
|---|---:|---:|
| 0.005 | -0.12 % | -0.12 % |
| 0.0075 | +0.37 % | +0.02 % |
| 0.009 | -0.14 % | -0.08 % |
| 0.010 | -0.12 % | -0.05 % |
| 0.0115 | -0.06 % | +0.40 % |
| end | -0.13 % | +0.22 % |

### Verdict

**1. The families do NOT part here, before or after the fix.** TanType 0 and
TanType 2 track within +/-0.4 % at every checkpoint on both binaries -- nowhere
near the Workbench's reported +7 / +13 / +22 %. This replay does **not**
reproduce the reported divergence at this scale/configuration. Candidate
reasons the gap is real but this deck doesn't show it: the Workbench's
"B/8 quad" mesh, footing width and boundary conditions are unknown and may
differ from this deck's rough-footing R3-graded strip; the Workbench may have
run the denser `e_init = 0.60` leg (this deck's own "ill-posed case", not
tried here) rather than Gorini's calibrated 0.6944; a surcharge or a larger
push reach past s/B = 0.012 (this leg is still climbing steeply, nowhere
near a peak); or an uncapped substep budget that seizes (see below) rather
than a clean converged answer at each checkpoint.

**2. TanType 2 IS reproducible here, bit-for-bit, once substeps are capped.**
`postfix_tan2_run1` and `postfix_tan2_run2` are identical to every reported
significant figure in `q_foot_kPa`, `q_base_kPa`, and every substep-census
column across all 46 steps (`diff`'d directly); only the cumulative `wall_s`
timing column differs, as expected from run-to-run scheduling noise. The
Workbench's reported *non*-reproducibility is more consistent with the
uncapped-substep GATE-U seizure hitting a wall-clock cutoff at a different
point each run (a machine-timing artifact) than with genuine numerical
nondeterminism in the fixed tangent's converged path -- this replay found no
nondeterminism once the seizure is avoided.

**3. The fix's footprint shows up in Newton cost, not in the converged q.**
Pre-fix TanType 2 needed 18 failed ladder rungs vs post-fix's 14 (+29 %) to
commit the same 46 steps, echoing (at a much smaller scale) the WP-110 F15
integrator finding of 283 -> 103 Newton iterations over 40 steps on the
drained-triaxial gate. The two binaries' TanType 0 legs (elastic tangent,
untouched by the fix) are byte-identical, as expected.
