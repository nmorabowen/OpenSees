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
