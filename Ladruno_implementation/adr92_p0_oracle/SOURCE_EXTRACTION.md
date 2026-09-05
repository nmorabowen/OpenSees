# ADR-92 P0 — `ManzariDafalias` source extraction

The oracle must reproduce **the code**, not Dafalias & Manzari (2004). This file is the
transcription the oracle is built from. All line numbers are
`SRC/material/nD/UWmaterials/ManzariDafalias.cpp` at `ladruno` `3f003d110`.

> **Trust level.** The items marked **[V]** were re-read and verified by the orchestrator.
> Everything else is a delegated transcription and may contain errors — which is exactly
> what **gate G0** exists to catch: if the oracle's `implicit` path reproduces the binary to
> 1e-8 on the recorded strain path, the transcription was right. **Re-read the source for
> anything that does not reconcile; do not "fix" a formula to make G0 pass.**

---

## 1. Elastic moduli — `GetElasticModuli`, three identical overloads (`:4716`, `:4759`, `:4779`)

```
pn = tr(sigma)/3                      # NO m_Presidual here
pn = (pn <= m_Pmin) ? m_Pmin : pn     # :4721 / :4763 / :4783
eG = mUseCurrentVoidRatioInG ? e : m_e_init
mElastFlag == 0:  G = G0*P_atm*(2.97-eG)^2/(1+eG)
mElastFlag == 1:  G = G0*P_atm*(2.97-eG)^2/(1+eG) * sqrt(pn/P_atm)
K  = (2/3)*(1+nu)/(1-2nu) * G
```

**`mUseCurrentVoidRatioInG` is `false` in all four constructors and never set** ⇒ `eG` is
**always `m_e_init`**, never the current void ratio, even though the current one is passed
in. This is ADR-86 §7.3's open defect; it is preserved deliberately. **The oracle must
reproduce it.**

`GetStiffness(K,G)`: isotropic, `a = K + 4G/3` on the normal diagonal, `G` on the shear
diagonal, `b = K - 2G/3` on the normal off-diagonals.

## 2. State parameter — `GetPSI(e,p)` (`:4697`)

`psi = e - (e0 - lambda_c * (p/P_atm)^ksi)`

**One call site only** — `:4946`, inside `GetStateDependent`, with
`p = tr(sigma)/3 + m_Presidual` floored at `small = 1e-10` (**not** at `m_Pmin`), `:4936-4937`.

## 3. Yield surface — `GetF(stress, alpha)` (`:4686`)

```
s = dev(sigma);  p = tr(sigma)/3 + m_Presidual      # no clamp
f = ||s - p*alpha||_contr - sqrt(2/3)*m*p
```

## 4. `GetStateDependent` (`:4927-4996`) — the constitutive core

```
p          = tr(sigma)/3 + m_Presidual ; p = max(p, small)
n          = GetNormalToYield(sigma, alpha)             # :4870
AlphaAlphaInDotN = (alpha - alpha_in) : n
psi        = GetPSI(e, p)
cos3Theta  = GetLodeAngle(n) = clamp(sqrt(6)*tr(n.n.n), -1, 1)     # :4704
g(c3t,c)   = 2c / ((1+c) - (1-c)*c3t)                              # :4679
alphaBtheta = g*Mc*exp(-nb*psi) - m
alphaDtheta = g*Mc*exp( nd*psi) - m
b0         = G0*h0*(1 - ch*e)/sqrt(p/P_atm)            # uses the CURRENT e, unlike G
d          = sqrt(2/3)*alphaDtheta*n - alpha
b          = sqrt(2/3)*alphaBtheta*n - alpha
h          = (|AlphaAlphaInDotN| < small) ? 1.0e10 : b0/AlphaAlphaInDotN
A          = A0 * (1 + <fabric:n>)                     # <> = Macauley
D          = A * (d:n)
if p < 0.05*P_atm:  D *= 1/(1 + exp(7.6349 - 7.2713*(101.0/P_atm)*p))   # :4971-4985
B          = 1 + 1.5*(1-c)/c * g * cos3Theta
C          = 3*sqrt(1.5)*(1-c)/c * g
R          = B*n - C*(n.n - I1/3) + (D/3)*I1
```

`GetNormalToYield` (`:4870`): `n = [dev(sigma) - p*alpha] / ||.||`, `p` including
`m_Presidual`; if `|p| < small`, `n = 0`; the norm is replaced by 1.0 if below `small`.

The `D_factor` sigmoid is the dimensional-dilatancy defect of ADR-86 §7.2 — a no-op at
`P_atm = 101`. Reproduce it as written.

## 5. `ForwardEuler` (`:1313-1362`) — the single-step explicit update

```
CurVoidRatio      = e_init - (1+e_init)*tr(CurStrain)
NextVoidRatio     = e_init - (1+e_init)*tr(NextStrain)
NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
aC = GetStiffness(K,G)
GetStateDependent(CurStress, CurAlpha, CurFabric, CurVoidRatio, alpha_in, ...)
dVolStrain = tr(dStrain) ; dDevStrain = dev(dStrain)
p  = tr(CurStress)/3 + m_Presidual
r  = 0                                              # SEE THE BUG BELOW
Kp = (2/3)*p*h*(b:n)
temp4 = Kp + 2G*(B - C*tr(n.n.n)) - K*D*(n:r)
if |temp4| < small: temp4 = small
NextDGamma = (2G*(n:dDevStrain)_mixed - K*dVolStrain*(n:r)) / temp4
dSigma  = 2G*dDevStrain^contr + K*dVolStrain*I1
          - <NextDGamma> * ( 2G*(B*n - C*(n.n - I1/3)) + K*D*I1 )
dAlpha  = <NextDGamma> * (2/3)*h*b
dFabric = -<NextDGamma> * cz * <-D> * (z_max*n + CurFabric)
dPStrain = NextDGamma * ToCovariant(R)
NextStress = CurStress + dSigma ; NextAlpha += dAlpha ; NextFabric += dFabric
NextElasticStrain = CurElasticStrain + dStrain - dPStrain
```

> [!bug] **[V] `r` is identically ZERO in `ForwardEuler`, and this is a real vanilla defect**
> `:1330-1332` reads
> ```cpp
> Vector r(6);
> if (p > small)
>     Vector r = GetDevPart(CurStress) / p;   // <-- declares a SECOND, block-scoped r
> ```
> The inner declaration shadows the outer and dies at the closing brace, so the `r` used at
> `:1334`, `:1337`, `:1342` and `:1350` is the zero vector. Both `(n:r)` terms — the
> volumetric coupling in `temp4` and in `NextDGamma` — are silently dropped.
>
> **Scope, measured, not assumed:** `ModifiedEuler` (`:1483`) and `GetElastoPlasticTangent`
> (`:4839`) both use the corrected two-statement form `r = GetDevPart(...); r /= p;` with the
> single-statement version left commented out directly above — i.e. **this trap was hit and
> fixed twice and missed once**. `explicit_integrator`'s `default:` branch is `ModifiedEuler`
> (`:1060-1061`), so scheme 2's silent fallback does **not** land on the buggy path either.
> **The defect reaches only the schemes that select `ForwardEuler`: 5, and the FE inner
> variants of `MaxStrainInc` / `MaxEnergyInc` (9, 4). It does NOT affect the deck default
> (scheme 1) and it does NOT affect the TIMs campaign.** Say it that way and no wider.
> Owed a `LEDGER_quirks` row and, separately from ADR-92, a two-line vanilla fix.
>
> **The oracle reproduces it verbatim** in its `forward_euler` path — an oracle that "fixes"
> it will fail G0 on scheme 5 and, worse, will look right.

## 6. `ModifiedEuler` (`:1365-1692`) — the substepping controller, the deck default

- `T = 0, dT = 1, dT_min = 1e-6`, `TolE = mHonorTolRInME ? mTolR : 1e-4` (`:1380`).
  `-honorTolR` defaults 0, so **`TolE = 1e-4` and `mTolR` is ignored** — a deck asking for
  `1e-10` silently gets six orders looser (ADR-86 PR-2 D7).
- Loop `while (T < 1.0)` (`:1453`).
- Pre-loop low-p guard (`:1399-1449`): if `tr/3 + m_Presidual < m_Pmin + m_Presidual`, then
  `NextStress = dev(NextStress) + m_Pmin*I1`, throttled warning (first 10).
- Substep void ratio interpolates from **`CurStrain + T*dStrain`** (`:1469`) — the Ladruno
  fix; upstream had `NextStrain + T*dStrain`.
- Two Heun stages, each the `ForwardEuler` formula but with `r` **correctly** assigned
  (`:1483`). `NextDGamma < -small` ⇒ set 0, elastic `dSigma`, `mUseElasticTan = true`.
- Error: `stressNorm = ||NextStress||`;
  `err = ||dSigma2 - dSigma1||` if `stressNorm < 0.5`, else `/(2*stressNorm)`.
- **[V] Reject with `dT == dT_min` FORCE-ACCEPTS** (`:1649-1662`): sets
  `mUseElasticTan = true`, applies the averaged plastic strain, takes `NextStress = nStress`,
  then caps the stress ratio — `eta = sqrt(13.5)*||dev||/tr`, and if `eta > Mc`,
  `NextStress = tr/3*I1 + (Mc/eta)*dev` — recomputes `NextAlpha` from the capped stress, and
  **advances `T += dT` anyway**. This is #792 T1's target: the update always "succeeds", so
  nothing upstream can react. Reproduce it, and count how often it fires.
- Accept: `Stress_Correction` (drift correction), `T += dT`,
  > **[V] CORRECTED 2026-09-05.** A first version of this line said `Stress_Correction` was
  > "inert unless `mStressCorrectionInUse`". **It is ON by default** — `mStressCorrectionInUse
  > = true` in all four constructors (`:217`, `:303`, `:368`, `:429`), switchable only via
  > `setParameter` (`:868`). The P0 oracle found this when its G0 residual would not close
  > without it. The gate did its job.
  running consistent-tangent recursion `aCep_Consistent = aCep_step*(aD*aCep_Consistent + T*IImix)`.
- Growth/shrink `q = max(0.8*sqrt(TolE/err), floor)`, floor **0.1 on reject, 0.5 on accept**;
  `dT = max(q*dT, dT_min)`, then `dT = min(dT, 1-T)` on accept.
- **[V] No step-total `dGamma` is accumulated.** Every substep overwrites `NextDGamma`; the
  caller sees the last substep's value. **This is ADR-92 §C2 and it is why D1 extrapolates
  plastic strain.**

## 7. `BackwardEuler_CPPM` (`:2189-2461`) — the only implicit return

- Unknowns `x = [stress(6), alpha(6), fabric(6), dGamma]` (19).
- Residual (`NewtonRes`, `:2972-3035`):
  ```
  R1 = eStrain - TrialElasticStrain + dGamma*ToCovariant(R)     (6)
  R2 = alpha  - curAlpha  - dGamma*(2/3)*h*b                    (6)
  R3 = fabric - curFabric - dGamma*(-cz*<-D>*(z_max*n + fabric))(6)
  R4 = GetF(stress, alpha)                                      (1)
  ```
- `NewtonIter2` (`:2876`): `MaxIter = 30`; tolerance `tolR_loc = ||R_0||*mTolR + mTolR`
  (relative **plus** absolute); the line-search block is entirely commented out.
- Elastic trial moduli come from **`CurStress`** (`:2223-2226`), i.e. the committed stress.
- Low-p branch at `:2234` tests `tr/3 < m_Pmin` **without** `m_Presidual` — the one place that
  omits it. **[V] Inside that branch the Newton is disabled outright**: the call to
  `NewtonIter2_negP` is commented out and replaced by a literal `errFlag = 0;` (`:2264-2266`,
  with the author's note *"tension-cutoff surface … not working properly. Using explicit
  integrator for the time being"*), so **every low-p step under scheme 2 is integrated by
  `explicit_integrator`, i.e. `ModifiedEuler`**, and flagged as success. Scheme 2 is
  therefore NOT an implicit return where the campaign's problem lives.
- **[V] Non-convergence NEVER propagates.** The retry ladder is: (1) a 50-step forward-Euler
  warm start then re-Newton; (2) recursive bisection into two `BackwardEuler_CPPM` calls at
  `implicitLevel+1` (guard `mMaxSubStep = 10`); (3) **`explicit_integrator(...); errFlag = 1;`**
  (`:2434-2437`) — it gives up and silently returns an explicit answer flagged as success.
  **This is load-bearing for D3:** "the implicit companion cannot fail" is true only in the
  sense that it cannot *report* failure. P0 must instrument how often the fallback fires.

## 8. `integrate()` (`:957-994`)

```
trialDirection = mCe * (mEpsilon - mEpsilon_n)          # mCe = the COMMITTED elastic tangent
if (mAlpha_n - mAlpha_in_n) : trialDirection < 0:  mAlpha_in = mAlpha_n     # reversal
else:                                              mAlpha_in = mAlpha_in_n
if mElastFlag == 0:  elastic_integrator(...)
elif mScheme == 2:   BackwardEuler_CPPM(...)
else:                explicit_integrator(...)
```

`mElastFlag` is a **static** member shared by every instance — `updateMaterialStage` flips it
for all Gauss points at once.

`elastic_integrator` (`:1005-1020`): `NextStress = CurStress + aC:dStrain` with the moduli
from `GetElasticModuli(CurStress, NextVoidRatio, ...)`; `NextAlpha = dev(NextStress)/p` only
if `p > small`; **fabric and `dGamma` are not arguments and are not touched**.

## 9. Committed vs trial, and `commitState` (`:467-511`)

| trial | committed |
|---|---|
| `mEpsilon` `mSigma` `mEpsilonE` `mAlpha` `mAlpha_in` `mDGamma` `mFabric` | `mEpsilon_n` `mSigma_n` `mEpsilonE_n` `mAlpha_n` `mAlpha_in_n` `mDGamma_n` `mFabric_n` |

`mVoidRatio`, `mCe`, `mCep`, `mCep_Consistent`, `mK`, `mG` have **no trial/committed pair** —
single copies overwritten in place. `commitState` recomputes
`mVoidRatio = e_init - (1+e_init)*tr(mEpsilon)` and rebuilds `mK`, `mG` from `mSigma`
(`:481-486`), and clears `mUseElasticTan` if `tr(mSigma) > 0.01*P_atm`.

**For ADR-92 §2 this is the good news:** `eps_p(n) = mEpsilon_n - mEpsilonE_n` is committed
and exact, so D1 needs exactly one new committed vector (`d_eps_p`) and no vanilla hook.

## 10. Low-stress device map

`m_Pmin = 1e-4*P_atm`, `m_Presidual = 1e-2*P_atm` at `:889-890` — both overridden by
`LadrunoSANISAND` (`-Pmin 0.0101`, `-Presidual 0.0` for the campaign).

- Clamps the *pressure used for the moduli*: `:4721`, `:4763`, `:4783`.
- Rebuilds the *stress tensor* at the floor: `:1100` (`explicit_integrator`), `:1429`
  (`ModifiedEuler`), `:1945` (`RK45`), `:2251` (`BackwardEuler`), `:2630`
  (`Stress_Correction` bailout), `:5001-5003` (`Elastic2Plastic`).
- Perturbs `p` (and therefore `psi`, `M^b`, `M^d`, `D`): every `+ m_Presidual` site, of which
  the consequential one is `:4936` feeding `GetPSI`. At the campaign's `-Presidual 0.0` all of
  these are inert — **which is why the oracle must be run at `Presidual = 0` AND at the
  vanilla `1.01` to show the difference is the campaign's, not ours.**
