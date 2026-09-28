# SANISAND outside OpenSees: implementations, integration practice, oracle plan, and WP-129 recommendations

Survey for the Ladruno fork's TIMs soil work (F18, WP-128/129). Written 2026-09-26. Read-only research: web
search/fetch plus the fork's local docs and `ManzariDafalias.cpp`. No external code was downloaded, built or
run. Where a paper was paywalled or the host returned HTTP 403, this is stated and the claim is marked
**[unverified]**.

Tags used throughout:
- **[E]** established practice, with a source I read (full text or an authoritative abstract).
- **[E-sec]** established, but confirmed only through a secondary source that cites the original.
- **[R]** from my recollection of a paper I could not re-open in this session. Treat these as unverified until
  someone reads the original.
- **[I]** my own inference, including inferences from reading the fork's C++.

Line numbers are for this worktree's HEAD (`3dd487633`, WP-128 branch, which carries WP-127's instrumentation).
They are about 50 lines later than the `ladruno` numbers the orchestrator quoted (for example, the error block
the orchestrator gives as `:1866-1873` is at `:1972-1981` here). The function names are stable.

---

## 0. Executive summary

1. **Inventory.** I found only one public SANISAND-2004 implementation with source code and a separate
   lineage from the UW C++: the **Martinelli–Miriano–Tamagnini (2013) Abaqus UMAT, updated by Mašín**, on
   soilmodels.com (Fortran, free registration, released "as-is"). Everything else falls into one of four
   groups:
   - binary-only or licence-gated: the Itasca FLAC3D DM04 UDM, PM4Sand for FLAC, and IC MAGE M11
     (restricted, PLAXIS);
   - closed: numgeo;
   - a port of the OpenSees code: the `cgl-sd/SANISand04` MATLAB repo;
   - OpenSees itself: SAniSandMS.

   numgeo's SANISAND-MSf defaults match `ManzariDafalias`'s constants exactly (stress tolerance 1e-4,
   yield-function tolerance 1e-7, minimum substep 1e-6, 50 drift-correction iterations). That suggests a
   shared lineage, so numgeo is not an independent check **[I]**.
2. **Sloan–Abbo–Sheng (SAS 2001) against our `ModifiedEuler`.** Our integrator has SAS's skeleton:
   Heun plus an embedded Euler error, Pegasus-style intersection, an unloading test and consistent drift
   correction. It differs from SAS in six places that matter at the ring (§3):
   1. **The error omits the internal variables (α, z).** SAS measures stress **and** hardening variables.
   2. **The growth factor has no upper bound, and the error has no floor.** A substep whose two stages both
      take the elastic branch has err = 0, so q = +∞ and the next substep takes the whole rest of the
      increment. SAS caps q at 1.1 [E-sec] and, I believe, floors R at the machine epsilon EPS [R].
   3. **Failures at `dT_min` are force-accepted.** They are not reported as failures.
   4. **Drift correction may exit with f > FTOL silently.**
   5. **A start state with f > TolF goes straight into substepping.** No correction and no refusal comes first.
   6. **The loading/unloading decision is taken from the sign of Λ.** When the plastic-modulus denominator
      turns negative (Kp < 0 with η beyond M^b), a genuinely loading substep is classified as elastic, and
      α is re-derived to follow r with no bounding check. This is a strong candidate for how
      η/M^b ≈ 6 arises (finding B) **[I, code-read; WP-128 should confirm]**.
3. **Oracle.** The best independent oracle is (O1) a **paper-derived stiff-ODE reference integrator**, written
   from the Dafalias & Manzari (2004) equations and not from the C++ (2–3 days, no licences). The best
   cross-implementation check is (O2) the **Tamagnini/Mašín UMAT driven by Mašín's `triax` element driver or
   Niemunis' Incremental Driver** (2–4 days including the owner's registration, download and licence
   decision). O2 has convention traps, listed in §4.3. The biggest is that UW's low-p `D_factor` sigmoid
   is not part of DM2004 and is active below 0.05·P_atm = 5.05 kPa, which covers most of the ring.
4. **What WP-129 should adopt** (§5):
   1. Include α (and z) in the error, with a mixed absolute/relative scale. **[E]** SAS / Carow & Rackwitz;
      Fellin et al. 2009; the fork's own `RungeKutta45`.
   2. Cap step growth and floor the error. **[E-sec]** / **[R]** SAS.
   3. Classify loading from the elastic-trial numerator, and treat a non-positive denominator explicitly
      (PM4Sand sets Kp = 0 outside the bounding surface). **[E]** PM4Sand manual §2.7.
   4. Add an α-admissibility check or projection after each accepted substep. **[E]** in PM4Sand's FLAC
      implementation; its application to SANISAND is **[I]**.
   5. Refuse, don't force-accept, at `dT_min` and on a failed drift correction, with a counted and forwarded
      refusal. Abaqus's PNEWDT cutback is the established idiom.

---

## 1. Q1: implementations outside OpenSees

"Single-point replay" means: can it be started from an arbitrary (σ, α, α_in, z, e) and driven by a
prescribed dε?

| # | Implementation | Model | Language / host | Licence / access | Integration | Maintenance | Single-point replay from a given state? | Independent of UW `ManzariDafalias`? |
|---|---|---|---|---|---|---|---|---|
| 1 | **Martinelli, Miriano, Tamagnini (2013); Mašín update (2015); PLAXIS interface.** soilmodels.com "SANISAND Abaqus UMAT and Plaxis implementations" | DM2004 (19 parameters, 36 state variables, "Ptmult" artificial cohesion, 0–5 kPa) | Fortran: `umat.for`, `mod_sas.for`, `stvarlini.for`, `usr_add.for`, `usrmod.for` | Free registration. No OSS licence stated; Mašín calls it "as-is, for experienced users" | **[unverified]** Not documented on the page. The Tamagnini/Mašín UMATs on the same site use explicit adaptive RK (RKF-23) with error control; it is plausible that this one does too. Read `umat.for` to confirm. (`mod_sas` may mean SAniSand, not Sloan–Abbo–Sheng.) | Package last updated 2018-08-29; forum help only | **Yes, through any UMAT driver.** `triax` (Mašín: `straincontrol_wholetensor`, state-variable initialisation, Linux recompilation with your `umat.f`) or Niemunis' Incremental Driver (Fortran source, free). An "expert" initialisation mode sets the state variables manually | **Yes.** Separate code lineage (Perugia) |
| 2 | Itasca UDM "Dafalais-Manzari" (DM04), Cheng, Dafalias & Manzari (2013) | DM2004 | C++ DLL for FLAC3D 7/9 | Download from Itasca's UDM library; no source seen; needs a FLAC3D licence | Not documented; FLAC is explicit dynamic relaxation | UDM v4.0, 2023-04-28 | Only inside FLAC3D (a one-zone model). Impractical | Probably yes (Z. Cheng), but it cannot be inspected |
| 3 | PM4Sand (Boulanger & Ziotopoulou), FLAC/FLAC2D DLL | DM2004-derived, plane strain, no Lode angle | C++ DLL | Free DLL, no source | Explicit, **no substepping**. Zone-averaged drift correction plus a bounding-surface projection (§2c) | v3.3, June 2023, active | FLAC only | Different model. Useful for **practice**, not as an oracle |
| 4 | IC MAGE Model 11 (Tantivangphaisal, Imperial) in IC-UMIP | DM2004 "with control parameters added for BVPs" | Fortran UDSM for PLAXIS | **Restricted files** on Zenodo (10.5281/zenodo.7765065); UMIP v3.6 (10.5281/zenodo.13890543) | UMIP states it uses "Modified Euler sub-stepping with automatic error control", citing Sloan | 2023–2024 | Only through PLAXIS, unless Imperial shares a driver | Yes. But gated, and needs PLAXIS |
| 5 | numgeo (Machaček & Staubach) | SANISAND ("Sanisand-2"), SANISAND-F, SANISAND-MSf; also hypoplasticity + ISA | Own FE code; element-test inputs provided | Documentation public; no licence or source statement found | Modified Euler with error control (default), Forward Euler, constant-substep explicit. Defaults: stress tolerance 1e-4, F tolerance 1e-7, minimum substep 1e-6, up to 50 drift iterations, p_min 0.5 kPa, compression positive, optional M_b ≥ M_c fix | Active (docs dated 2026-09) | Element tests, yes. Loading an arbitrary state is unclear ("user-defined initial states" exists) | **Doubtful [I].** Its defaults equal UW's constants, and its benchmark page compares against OpenSees |
| 6 | `github.com/cgl-sd/SANISand04` | DM2004 | MATLAB | No licence file seen | Modified Euler / RK45 with error control (README), nSub = 10 | 6 commits | Yes (`inputData.m`) | **No.** The README says it is ManzariDafalias.cpp switched to MATLAB |
| 7 | Carow & Rackwitz (2021) research code | Li (2002) two-surface model, zero elastic range | Fortran, run in Niemunis' Incremental Driver and ANSYS | Not published | Explicit SAS without the yield-surface steps; implicit BE with numerical Jacobian and halving substeps | Paper only | n/a | Different model. Useful as the **open-access description of SAS error control** (§3) |
| 8 | SANISAND-Z implicit DIRK UMAT (C&G 2023, S0266352X23006560) | SANISAND-Z | Abaqus UMAT | Paper only; code not found | Two-stage DIRK with an embedded error estimate and substepping | 2023 | n/a | Different model |
| 9 | Petalas & Dafalias (2019), C&G 112 | SANISAND-Z | Research code | Paper only | Backward Euler with damped Newton (per the SANISAND-Z DIRK paper's summary) | — | n/a | Different model |
| 10 | Zhang et al. (2026), IJNAG 10.1002/nag.70354 | "Highly nonlinear sand models" including SANISAND-04 (per the search snippet) | — | Paywalled (403); code availability unknown | Fully implicit, line-search Newton, complex-step Jacobian | 2026 | Unknown | Probably yes. **[unverified]** |
| 11 | NTUA-SAND (Andrianopoulos, Papadimitriou & Bouckovalas 2010), IJNAG 34:1586, 10.1002/nag.875 | Bounding surface, vanished elastic region | FLAC UDM | Not public | Adaptive explicit substepping: second-order modified Euler vs first-order Euler error | Paper | FLAC only | Different model |
| 12 | Kratos, MOOSE, Code_Aster / MFront | — | — | — | **No SANISAND found.** Kratos GeoMechanics can host a UMAT/UDSM DLL, so it could run #1 | — | — | — |
| 13 | ISA-hypoplasticity (Fuentes) and Sand/Clay hypoplasticity UMATs (soilmodels.com) | Close relatives, not bounding surface | Fortran UMAT | Free registration | Explicit adaptive RK with error control (RKF-23 for the Mašín UMATs, per Acta Geotech. 10.1007/s11440-018-0684-z) | Active | Yes, same drivers as #1 | n/a |

Also in OpenSees but relevant: SAniSandMS (Liu et al. 2019) uses RK4 with error control per the OpenSees
documentation. Our own `RungeKutta45` in `ManzariDafalias.cpp` has the comment "after Sloan (J. Abell @ UANDES
added)".

**Bottom line for Q1.** Candidate #1 is the only one with source code and an independent lineage. Candidates
#4 and #2 are independent but gated behind PLAXIS or FLAC and closed. #5 and #6 share UW's lineage.

---

## 2. Q2: how mature codes handle our failure modes

### (a) The substep error norm near σ → 0

- **SAS / Sloan 1987.** The relative error is the maximum over the stress error and the hardening-variable
  error, each normalised by the magnitude at the end of the substep **[E-sec: Carow & Rackwitz 2021, eq. 48,
  citing Sloan 1987 and Sloan et al. 2001]**. I recall SAS also takes `max(…, EPS)`, with EPS the machine
  constant, so that R never reaches 0 in the step-size formula **[R]**.
- **The singularity is acknowledged.** Carow & Rackwitz (2021, §3.3.1) note that the relative normalisation is
  singular when a norm goes to zero, "for example in boundary value problems with a free surface". For that
  case they point to the **mixed relative plus absolute error of Fellin, Mittendorfer & Ostermann (2009)**,
  C&G 36:698, 10.1016/j.compgeo.2008.10.002 **[E]**. That is our ring exactly.
- **General ODE practice.** Hairer, Nørsett & Wanner, *Solving ODEs I*, §II.4 scales each component as
  `sc_i = ATOL_i + RTOL·max(|y_i|, |ŷ_i|)` and uses the RMS or max of `err_i/sc_i` **[E, textbook]**. That scale
  is continuous. Ours switches abruptly at ‖σ‖ = 0.5 (a value that carries units) and ignores α.
- **Our own `RungeKutta45`** already takes the maximum of the stress error and the α error (the orchestrator's
  `:2359-2362` on `ladruno`). `ModifiedEuler` does not (**finding E**).

### (b) Drift correction back to the yield surface

- **SAS.** After each accepted substep, SAS applies Potts & Gens' consistent correction, which updates stress
  **and** the hardening variables. If that increases |F|, SAS falls back to a normal correction, and it iterates
  until |F| ≤ FTOL. FTOL is small; about 1e-9 is the value "recommended by Sloan et al." per the 2026 Zenodo
  archive (§2 bullet on intersection). My recollection is that exceeding the iteration limit is an error exit,
  not an acceptance **[R]**.
- **Ours.** `Stress_Correction` has the same consistent-then-normal structure (`:2924` ff., up to 50
  iterations, a bisection fallback). **If neither correction reduces |f| it `return`s silently with f > TolF**
  (the "Couldn't decrease the yield function" branch). TolF is 1e-7 in stress units.
- **PM4Sand (FLAC).** The drift correction projects α along the zone-averaged stress ratio **[E, manual §3.2]**.
- **Numerical alternatives.** Sołowski, Sheng & Sloan: "NICE" (next increment corrects error), a
  drift-reduction scheme **[E, title and abstract]**.

### (c) Keeping α inside the bounding surface

- **DM2004 itself.** It relies on Kp ∝ (α^b − α):n going negative outside the bounding surface, which pulls α
  back. Boulanger & Ziotopoulou describe this as one of two ways to bound the stress state: allow it to sit
  "marginally outside" with a negative modulus. The other way is to enforce consistency on the bounding surface.
- **PM4Sand v3.3.** It does both of the following:
  - **Kp is set to zero, never negative, when the stress ratio is outside the bounding surface**, for numerical
    stability (manual §2.7, below eq. 38) **[E]**.
  - After each FLAC step, if the zone-averaged stress ratio lies outside the bounding surface, **it is projected
    back along the normal to the bounding surface** (manual §3.2) **[E]**.

  PM4Sand also scales α_in at initialisation so that its stress ratio is ≤ 0.9·M^b (§2.5), and it warns that
  initial states outside the bounding and dilatancy lines cause problems (§3.6) **[E]**.
- **Kan & Taiebat (2014), C&G 55:103.** They address bounding-surface **overshooting** (spurious stiffness
  after small reversals) with "clouds" of loading surfaces **[E, abstract]**.
- **Chen, Ghorbani, Zhang & Kodikara (2022), C&G 152:105008.** They address the same overshooting inside
  explicit substepping **[E, abstract]**.

  Both are the related reversal problem, not α escaping α^b.
- **I found no published SAS adaptation that adds an explicit α ∈ bounding-surface admissibility check for
  SANISAND.** PM4Sand's projection is the nearest established practice.

### (d) h → ∞ when (α − α_in):n → 0 after a reversal

- **DM2004 / UW.** h = b0/((α−α_in):n). UW caps it at h = 1e10 when |(α−α_in):n| < 1e-10
  (`GetStateDependent`, about `:5385`). **h is negative when (α−α_in):n < 0**, which can happen inside a
  substep because α_in is re-seated only once per global increment (`integrate()`, `:1027-1033`) **[I, code]**.
- **PM4Sand.** It regularises the denominator as `exp((α−α_in^app):n) − 1 + C_γ1` with **C_γ1 = h0/200**,
  "to avoid division by zero". It also uses an apparent α_in to avoid over-stiffening after small cycles
  (§2.5, eqs. 37–38) **[E]**. This changes the model.
- **Li (2002), as used by Carow & Rackwitz.** A reversal makes Kp infinite only for an instant: the projection
  centre jumps to r, and the next increment is plastic again **[E, Carow §2.7.2]**.
- **Mathematics [I].** As h → ∞, Λ ~ 1/h and hΛ stays finite (α tracks r through consistency), so the rate
  form itself is not singular. The discrete danger is the product (2/3)hΔΛ per substep. The update
  dα = (2/3)ΔΛ·h·(α^b − α) is explicit Euler on a relaxation towards α^b. It overshoots α^b once
  (2/3)hΔΛ > 1 and is unstable once (2/3)hΔΛ > 2.

### (e) p' → 0 and tension

- **PM4Sand.** Mean stress floored at **0.5 kPa or 0.005 × the initial consolidation stress**. Nominal shear
  resistance in surface zones through c_hg **[E, §3.2, §3.6]**.
- **numgeo SANISAND-MSf.** Minimum mean effective stress 0.5 kPa by default **[E, docs]**.
- **Tamagnini/Mašín UMAT.** `Ptmult` artificial cohesion, "0 if the analysis runs smoothly, up to 5 kPa to
  stabilise" **[E, package readme as summarised on the soilmodels page]**.
- **The fork** already has `-Pmin`, `-Presidual` and `-pRe` (ADR-86/93). No external code offers anything the
  fork lacks here. The difference is that the others **floor early (0.5 kPa or more)** and **say so**.

### (f) Implicit CPPM for bounding-surface models under a global Newton

- **Carow & Rackwitz 2021, Algorithm 1.** Backward Euler with a numerical Jacobian (Pérez-Foguet et al. 2000).
  If the local Newton fails, ΔT is halved down to ΔT_min, and **at ΔT_min the result is accepted with a
  warning** **[E]**. Their verdict: for a given accuracy, the explicit update is significantly more efficient
  than the implicit one for zero-elastic-range bounding-surface models **[E]**.
- **Implicit substepping for critical-state models** (arXiv 2504.17476, Eng. Comput. 2026): when the local
  Newton fails, double the number of substeps and retry. The consistent tangent is linearised through the
  substeps so that the global Newton stays quadratic **[E]**.
- **SANISAND-Z:** Petalas & Dafalias (2019) use damped Newton; the 2023 DIRK paper uses an embedded error with
  substepping **[E, abstract]**.
- **Zhang et al. 2026:** line-search Newton plus a complex-step Jacobian for SANISAND-type models
  **[E, abstract only]**.
- **Manzari & Prachathananukit (2001)**, IJNAG 25:525: for a two-surface cyclic sand model, closest-point
  projection stayed stable at large increments where explicit substepping struggled **[E, abstract]**.
- **Refusal in the host code.** Abaqus UMATs set **PNEWDT < 1** to request a global cutback. That is the
  standard idiom for "refuse instead of force". Whether the Tamagnini/Mašín UMAT uses it is **[unverified]**.

### (g) Tangents

- SAS-style codes and Carow & Rackwitz hand back the **continuum** tangent at the end of the increment
  **[E, Carow §4.1]**.
- Consistent tangents for adaptive explicit RK:
  - Fellin & Ostermann (2002);
  - Monforte et al. (2025), IJNAG 10.1002/nag.70016;
  - a C&G 2026 paper, S0266352X26004155.

  All are **[E, abstracts]**.
- Numerical differentiation: Pérez-Foguet et al. (2000) **[E-sec]**.
- The fork's `TanType 2` chain (`aCep_Consistent = ½(aCep1+aCep2)·(…)`) is a homegrown variant (finding D).

### Intersection and the "silent elastic" failure (context for (b) and (c))

A 2026 reproducibility archive (Kaewhanam et al., Zenodo 10.5281/zenodo.21296348 and 21379836, MIT and CC-BY,
Python material-point scripts) reports two failures of SAS's intersection machinery on transformed-stress
bounding-surface models:
1. a **loud** failure: Pegasus exhausts its iterations at FTOL = 1e-9 near the double-precision floor;
2. a **silent** failure: plastic increments are classified as elastic without warning.

Its fixes are a safeguarded Pegasus with bracket-collapse acceptance, and a C0-continuous invariant **[E,
archive description; the paper itself is unverified, possibly a preprint]**. Its "silent" class is the same
failure class as our f > 0 point integrated elastically with err 0.

---

## 3. SAS vs our `ModifiedEuler`, item by item

Sources for SAS:
- Sloan, Abbo & Sheng (2001), Eng. Comput. 18(1/2):121–154, 10.1108/02644400110365842.
- Sloan (1987), IJNME 24:893.

**The SAS full text was not reachable in this session** (Emerald paywall; the Newcastle author copy returned
403). SAS items below are therefore [E-sec] where Carow & Rackwitz (2021, open access, 10.14279/depositonce-15079)
restate them, and [R] otherwise.

"Ours" means `ManzariDafalias::ModifiedEuler` and `explicit_integrator` in this worktree.

| # | Item | SAS (2001) | Ours | Consequence at the ring |
|---|---|---|---|---|
| 1 | Variables in the error | Stress **and** hardening/internal variables, max of the relative norms [E-sec, Carow eq. 48] | **Stress only.** dα and dz are computed per stage but their difference is ignored (finding E). `RungeKutta45` does include α | Large α drift passes the test. Stage disagreement in α at h → ∞ is invisible |
| 2 | Normalisation | Magnitude at the end of the substep (the higher-order solution) [E-sec] | `GetNorm_Contr(NextStress)`, i.e. the **start** of the substep (`:1972`) | Minor |
| 3 | Guard near zero | R = max(…, **EPS**) [R]. The σ → 0 problem is acknowledged; fix via mixed abs/rel (Fellin 2009) [E-sec] | Hard switch: absolute below ‖σ‖ = 0.5, relative above (`:1972-1981`). Unit-bearing, discontinuous | Tolerance jumps by 2‖σ‖ across 0.5 kPa |
| 4 | Step factor on success | q = 0.9√(STOL/R), clamped to **0.1 ≤ q ≤ 1.1** [E-sec, Carow eqs. 50–52 citing SAS] | q = max(0.8√(TolE/err), 0.5), **no upper cap** (`:2036`), then `dT = min(dT, 1−T)` | **With err = 0 (both stages elastic), q = +∞, and the next substep swallows the rest of the increment [I, code-verified]** |
| 5 | Step factor on failure | q = max(0.9√(STOL/R), 0.1) [E-sec]. No growth allowed in the step that follows a failure [R] | q = max(0.8√(TolE/err), 0.1). No post-failure growth limit | Oscillating accept/reject |
| 6 | Minimum step | ΔT = max(qΔT, ΔT_min). I believe reaching ΔT_min is reported as failure [R] | `dT_min = 1e-6`. **A failed substep at dT_min is accepted**: elastic tangent, η clamped to **Mc** (not M^b, compression side only), α re-derived as `CurAlpha + 3(dev σ/tr σ − dev σ₀/tr σ₀)`, no bounding check (`:1983-2011`) | Can manufacture an inadmissible committed α. Counted since WP-127 (`LMS_FORCED_DTMIN`) |
| 7 | Tolerance source | STOL is user input | `TolE = mHonorTolRInME ? mTolR : 1e-4` | Known (F18) |
| 8 | Drift correction | After each accepted substep; consistent (Potts–Gens) with normal fallback; iterate to FTOL; failing is an error [R] | Same structure (`Stress_Correction`); **returns silently with f > TolF** when neither correction reduces f; absolute TolF 1e-7 | **f > 0 exits as success** |
| 9 | Elastic→plastic intersection | Pegasus to FTOL [E, Zenodo archive citing SAS] | `IntersectionFactor` (`explicit_integrator`, `:1238-1247`) | — |
| 10 | Elastoplastic unloading | Test cosθ < −LTOL, subdivide NSUB, then Pegasus [R] | Test `n:dσ/‖dσ‖ > −√TolF` → plastic, else `IntersectionFactor_Unloading` (`:1252-1276`) | — |
| 11 | Start state outside the yield surface (f₀ > FTOL) | Assumed impossible after drift correction [R] | Debug print only, then substepping straight from the illegal state (path 1, `:1229-1236`) | Illegal states propagate |
| 12 | Loading/unloading inside a substep | Classical hardening: the denominator is positive, and Λ comes from ⟨·⟩ of a positive-denominator formula | **Sign of Λ** (`NextDGamma < 0 → elastic`, `:1826-1834`, `:1901-1907`). The elastic branch sets dα = 3(dev σ/tr σ − …): **α follows r with no bounding check** | **If the denominator Kp + 2G(B − C·tr n³) − K·D·n:r < 0 (Kp < 0 with η beyond M^b, magnified by large h), a loading substep with a positive numerator gets Λ < 0 and is treated as elastic. Both stages agree, err = 0, q = ∞, and α is dragged outward with r. This is a plausible generator of η/M^b ≈ 6 [I, WP-128 to confirm]** |
| 13 | Tangent returned | Continuum at the end [E-sec] | TanType 0/1/2 (WP-110 end-state continuum; product chain for 2) | Finding D |

### 3.1 "SAS-ME": a spec for a new `IntScheme` in `ManzariDafalias`

Proposal: a new scheme number, byte-identical defaults elsewhere, all constants as flags with defaults. Each
line is tagged [E] established, [E-sec], [R] or [I].

**Inputs.** STOL (= TolR, always honoured), ATOL_σ (kPa, default = STOL·P_atm) [I], ATOL_α (dimensionless,
default STOL) [I], FTOL (relative, default 1e-8) [I], LTOL 0.01 [R], dT_min 1e-4 [R], NSUB 10 [R],
MAXITS_drift 10 [R], safety 0.9 and q ∈ [0.1, 1.1] [E-sec], c_hΛ = 0.5 [I].

**0. Entry check.**
- Compute f₀ = F(σ₀, α₀). If f₀ > FTOL_abs, run the drift correction on the start state first. If that fails,
  refuse with code `START_INADMISSIBLE` [I; SAS assumes admissible starts].
- Check α₀ admissibility: `(α^b(θ, ψ) − α₀):n ≥ −tol_b`. If it fails, refuse, or project if `-alphaProject on`
  [I].

**1. Elastic predictor and intersection.** Keep the present Pegasus and unloading logic, with a safeguarded
bracket-collapse acceptance (Zenodo 2026) [E].

**2. Substep loop, T ∈ [0, 1].** For each stage k = 1, 2:
- **Loading criterion from the elastic trial** (numerator) N_k = 2G·n:de − K·dε_v·n:r [I, standard];
  denominator H_k = Kp + 2G(B − C·tr n³) − K·D·n:r.
- If N_k ≤ 0: elastic stage. α and z are **unchanged** (SAS), not re-derived from r [I]. If the elastic stage
  would leave f > 0, treat it as an intersection (cut the substep).
- If N_k > 0 and H_k > 0: plastic, with Λ = N_k/H_k.
- If N_k > 0 and H_k ≤ 0 (softening beyond α^b, or h < 0): either use a PM4Sand-style Kp := max(Kp, 0) under
  flag `-kpFloor zero` [E for PM4Sand, I for SANISAND], or reject the substep (reduce dT). Never classify it
  as elastic [I].
- Evaluate h with (α − α_in):n re-checked per stage. If (α − α_in):n < 0 inside the increment, reseat
  α_in := α at that stage (per-substep reversal detection) behind a flag [I].
- **Stiffness step limit:** if (2/3)·h·Λ_k > c_hΛ, reject with q = c_hΛ/((2/3)hΛ_k) [I; explicit-Euler
  stability of the α relaxation].

**3. Error.**

    R = max( max_i |Δσ₂−Δσ₁|_i / (ATOL_σ + STOL·max(|σ̂_i|, |σ̃_i|)),
             max_i |Δα₂−Δα₁|_i / (ATOL_α + STOL·max(|α̂_i|, |α̃_i|)),
             same for z,
             EPS )

The ½ factor is optional; keep SAS's. [E-sec for stress + internal variables; E for mixed abs/rel scaling
(Fellin 2009; Hairer); I for the exact form and the defaults.]

**4. Accept if R ≤ 1** (the mixed form folds STOL into the scale). Then:
- Update σ, α, z, ε^e.
- Run the drift correction (consistent, then normal, then refuse).
- Run the **α admissibility check**. If it fails, reject the substep (or project under the flag).
- Set q = clamp(0.9/√R, 0.1, 1.1). If the previous substep failed, q = min(q, 1) [R].

**5. Reject if R > 1.** q = max(0.9/√R, 0.1). If dT is already at dT_min: **refuse** with code
`DTMIN_FAIL`. Hand the increment to the element and global cutback (F7 refusal roster; like Abaqus PNEWDT)
[E idiom; SAS's own behaviour R]. Keep today's force-accept only as `-onDtMin accept`, counted and warned.

**6. Exit.**
- Assert |f| ≤ FTOL and α admissible. Otherwise refuse [I].
- Tangent: continuum at the end-state (TanType 0/1). TanType 2 is out of scope, or use numerical
  differentiation (Pérez-Foguet) later.

**Counters.** One per refusal code, per path, plus the existing WP-127 census.

**What is new math versus bookkeeping.** Steps 2 (the H ≤ 0 branch and the hΛ limit) and 4 (the α
check/projection) change behaviour on points like 1950/2–3. Steps 3, 5 and 6 are standard SAS hygiene.

---

## 4. Q3: oracle feasibility for the 80 ring states

### 4.1 What the data allow

- b8 has 40 rows. **Only two (element 1950, Gauss points 2 and 3) have η/M^b > 1** (6.13 and 5.75), and one of
  those has a tensile normal component in compression-positive terms (σ_xx = −0.71 kPa).
- b16's worst is 0.92.
- In every row the shear components yz and zx are zero (plane strain).
- **The CSVs carry no strain increment.** Every oracle run needs a prescribed dε: the F18(a) protocol, a
  fan of directions, or the dε from WP-128's chain driver.
- For the two inadmissible rows, an oracle cannot "validate" the start state. It can only show what each code
  does from it: refuse, project or drift. The valuable oracle work is on the 78 admissible rows, plus a
  re-driven path from an admissible ancestor state.

### 4.2 Candidates, ranked by value per effort

**O1: a paper-derived reference integrator.** Recommended first. 2–3 days, no licences.
- **What.** A Python implementation of the DM2004 rate equations (Dafalias & Manzari 2004, JEM 130(6):622,
  10.1061/(ASCE)0733-9399(2004)130:6(622)), written **from the paper, not from the C++ or WP-128's port**.
  Integrate them with SciPy `solve_ivp(method="Radau")` at rtol 1e-10 and scaled atol, using event functions
  for the loading-index sign change and the α_in reversal.
- **Toggles.** Switch on UW deviations explicitly: `D_factor`, `Presidual`, `Pmin`.
- **Why.** It is independent of both the UW lineage and any licence, and it answers the physics question:
  what the rate equations do from state X under dε.
- **Risk.** Discontinuous right-hand side at loading/unloading; the events handle it. Cost is irrelevant for
  80 points.

**O2: the Tamagnini/Mašín UMAT under `triax` (or Niemunis' Incremental Driver).** Recommended second.
2–4 days including owner steps.

Owner steps (approval needed, because this downloads and executes external code):
1. Register on soilmodels.com and accept its terms.
2. Download "SANISAND Abaqus UMAT and Plaxis implementations" (1.09 MB, Fortran source and PDFs, updated
   2018-08-29) and `triax` (Mašín; source included; free).
3. Decide on the licence status. None is stated, so use it internally only and do not vendor it into the fork.

Agent steps:
1. Read `umat.for` and its PDFs to extract the integration scheme, the 19-parameter order and the 36
   state-variable layout. This settles the [unverified] cells in §1.
2. Build `triax` on Linux (Esmeralda, gfortran) with `umat.for` in place of `umat.f`.
3. Write a converter from the CSV to triax input (expert state-variable initialisation, then
   `straincontrol_wholetensor`).
4. Run the 78 admissible rows under the same dε as the C++ replay.
5. Compare σ, α and z at the end of the increment.

Alternative driver: Niemunis' Incremental Driver (Fortran source, free; used by Carow & Rackwitz 2021).

**Not recommended as an oracle:**
- numgeo: closed, and apparently the same lineage.
- cgl-sd MATLAB: a port of our code.
- Itasca DM04 UDM: needs FLAC3D, DLL only.
- IC MAGE M11: restricted, needs PLAXIS. Worth an email to Imperial (Taborda) only if O2 fails.
- PM4Sand: a different model.

### 4.3 Convention traps for O2 (and partly O1)

| Trap | Ours (CSV / `ManzariDafalias`) | Abaqus UMAT convention | Action |
|---|---|---|---|
| Sign | Compression positive internally. The CSV is **compression positive**, despite what its README says (finding A) | Tension positive | Negate σ. **α, α_in and z are stress-ratio-like, so negate them too** (r = s/p flips with s, and z aligns with n). Keep e. Check how the UMAT defines p |
| Voigt order | xx, yy, zz, xy, yz, zx | 11, 22, 33, 12, 13, 23 | Swap components 4 and 5 (zero in all 80 rows, but map them anyway) |
| Shear strain | Covariant (engineering γ) strains; contravariant tensor stresses, α and z | Engineering γ in STRAN/DSTRAN; stress tensor components | Confirm the UMAT's storage of α and z shear components (tensor or engineering) from its docs |
| Lode function | g = 2c/((1+c) − (1−c)·cos3θ), cos3θ = √6·tr(n³), clamped | Unknown | cos3θ changes sign with n under the sign flip. Verify triaxial compression gives g = 1 in the UMAT |
| Low-p dilatancy | **`D_factor` sigmoid for p < 0.05·P_atm** (UW addition, not in DM2004); most ring points are below 5.05 kPa | Probably absent | Expect D to differ. In O1, make it a toggle. In O2, compare against a fork build with D_factor off, or accept the gap |
| p_r / Pmin | Fork: Presidual 0, Pmin 0.0101 kPa (as run) | Ptmult cohesion | Set Ptmult = 0, and document any internal UMAT floor |
| ν | 0.312885 on the ring (the Jaky substitution) | — | Use the CSV README's ν, not 0.3 |
| Elastic stage / flip | `-flipAlphaIn init` | "Standard" initialisation mode may overwrite state variables | Use the "expert" mode |
| Units | kPa, P_atm = 101 | Parameter p_atm | Same units |

---

## 5. Q4: recommendations for WP-129

Ordered by expected effect on the ring. Established practice is kept apart from inference.

1. **Close the "silent elastic" hole first.** This is inference from the code, backed by established
   analogues.
   - Classify stages from the elastic-trial numerator. Never map a negative denominator to "elastic". Never
     let an elastic stage re-derive α from r without a bounding check. (Items 12 and 6 in §3.)
   - Established analogue: PM4Sand sets Kp = 0 outside the bounding surface instead of letting it go negative
     (manual §2.7) [E]. The 2026 intersection-failure archive shows this failure class is real and silent
     elsewhere [E].
   - Test: WP-128's replay of 1950/2–3 and its ancestors, with the path counter for "Λ < 0 while N > 0".
2. **Error norm: add α and z, and use a mixed absolute/relative scale.**
   - Established: SAS measures stress and internal variables [E-sec]. The fork's own `RungeKutta45` does so
     too. Fellin et al. (2009) mixed abs/rel error for vanishing norms [E-sec]. Hairer's scaled norm [E].
   - The 0.5 kPa switch becomes a continuous ATOL (F18(a)'s `-errFloor` is exactly the ATOL_σ of §3.1).
   - Inference: ATOL_σ should be tied to the global Newton tolerance per point, which F18(a) already asks.
     Without α in the norm, a floor alone may make the ring cheaper *and* more wrong.
3. **Step control hygiene.** Cap q at 1.1 and floor R at EPS [E-sec / R]. Allow no growth right after a
   failure [R]. Add a stiffness step limit on (2/3)hΔΛ [I]. This kills the q = ∞ jump (item 4 of §3), which
   is code-verified.
4. **Refuse, don't force-accept.** Apply this at `dT_min`, when drift correction cannot reach FTOL, and when a
   start state is outside the yield surface or inadmissible.
   - Established idiom: Abaqus PNEWDT cutback [E]. Counter-practice exists: Carow & Rackwitz's implicit scheme
     accepts at ΔT_min *with a warning* [E], and ours accepts silently. SAS's own rule is [R].
   - The fork already forwards refusals (F7) and counts `LMS_FORCED_DTMIN` (WP-127). Make refusal the new
     scheme's default and keep acceptance behind a flag.
5. **α-admissibility after each accepted substep.** Reject the substep; optionally project α onto α^b along n.
   - Established only in PM4Sand's FLAC zone-level correction [E]. Its use in SANISAND and the
     reject-before-project ordering are [I].
   - A projection changes the constitutive answer. Keep it opt-in and counted, as ADR-86/93 did for the p
     clamps.
6. **h regularisation.** Keep UW's cap of 1e10 as default. Offer per-stage re-evaluation of the sign of
   (α − α_in):n [I]. Do **not** import PM4Sand's C_γ1: it changes the model and its calibration [E, PM4Sand
   §2.7].
7. **Low-p floors.** Nothing new is needed. External codes floor at 0.5 kPa (PM4Sand, numgeo) or add cohesion
   up to 5 kPa (Ptmult) [E], which matches the fork's existing `-Pmin`/`-Presidual`/`-pRe`. At most, document
   the comparison in the guide.
8. **Tangent.** Continuum end-state tangent for the new scheme (SAS practice) [E-sec]. The consistent RK
   tangent (Monforte 2025, C&G 2026) is a later work package [E, abstracts].
9. **CPPM (WP-130), for context.** Established remedies are halving or doubling substeps on local-Newton
   failure (Carow Algorithm 1; arXiv 2504.17476) and damped or line-search Newton (Petalas & Dafalias 2019;
   Zhang 2026). Refusing early through the host is consistent with them.

---

## 6. What I could not verify

- **The SAS (2001) full text:** the exact error formula (the ½ factor and the EPS term), the ΔT_min rule
  (fail or accept), the post-failure growth limit, the default constants (STOL, FTOL, LTOL, DTmin, NSUB,
  MAXITS), and the drift-correction exit rule. The Emerald paywall and a 403 on the Newcastle author copy
  blocked it. The step factor, q bounds and error variables are confirmed only through Carow & Rackwitz 2021.
- **The Tamagnini/Mašín UMAT's integration scheme** and state-variable layout. It needs registration to read
  `umat.for`. Whether `mod_sas` means Sloan–Abbo–Sheng is unknown.
- **The Itasca DM04 UDM's** scheme and source availability. Whether IC MAGE M11 has a driver outside PLAXIS.
  numgeo's licence and whether it derives from the UW code (inferred from identical defaults only).
- **Andrianopoulos et al. 2010** details beyond the abstract: whether drift correction is used, and their
  tolerances.
- **Zhang 2026 and the SANISAND-Z DIRK paper** (403): whether code exists, and which SANISAND version.
- **The 2026 Zenodo "intersection failures" paper's** publication status.
- Whether any published code enforces α ∈ bounding surface for SANISAND specifically. I found none.
- The §3 item-12 mechanism (a negative denominator taken as "elastic", with α dragged along with r). It is
  code-read inference and has not been executed. WP-128 owns that test.

---

## Sources

**SAS family and integration**
- Sloan, Abbo & Sheng (2001). *Refined explicit integration of elastoplastic models with automatic error
  control*, Eng. Comput. 18(1/2):121–154. https://doi.org/10.1108/02644400110365842 (abstract only here)
- Sloan (1987). *Substepping schemes…*, IJNME 24:893. https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240505
- Carow & Rackwitz (2021). *Comparison of implicit and explicit numerical integration schemes for a bounding
  surface soil model without elastic range*, C&G 140:104206. Open access (read in full, §3):
  https://doi.org/10.14279/depositonce-15079
- Fellin, Mittendorfer & Ostermann (2009). *Adaptive integration of constitutive rate equations*,
  C&G 36:698. https://www.sciencedirect.com/science/article/abs/pii/S0266352X08001547
- Sołowski & Gallipoli (2010). *Explicit stress integration with error control for the BBM*, Parts I and II,
  C&G 37. https://www.sciencedirect.com/science/article/abs/pii/S0266352X09001219
- Zhao, Sheng, Rouainia & Sloan (2005). *Explicit stress integration of complex soil models*, IJNAG 29:1209.
  https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.456
- Andrianopoulos, Papadimitriou & Bouckovalas (2010). IJNAG 34:1586. https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.875
- Manzari & Prachathananukit (2001). IJNAG 25:525. https://onlinelibrary.wiley.com/doi/10.1002/nag.140
- Kan & Taiebat (2014). C&G 55:103. https://www.sciencedirect.com/science/article/abs/pii/S0266352X13001213
- Chen, Ghorbani, Zhang & Kodikara (2022). C&G 152:105008. https://www.sciencedirect.com/science/article/abs/pii/S0266352X22003457
- Implicit substepping for critical state models (2025/2026). https://arxiv.org/html/2504.17476
- SANISAND-Z implicit DIRK with substepping (2023). https://www.sciencedirect.com/science/article/abs/pii/S0266352X23006560
- Petalas & Dafalias (2019). C&G 112. https://www.sciencedirect.com/science/article/abs/pii/S0266352X19301120
- Zhang et al. (2026). IJNAG. https://onlinelibrary.wiley.com/doi/10.1002/nag.70354
- Monforte et al. (2025). IJNAG. https://onlinelibrary.wiley.com/doi/10.1002/nag.70016
- Kaewhanam et al. (2026). Reproducibility archives: https://zenodo.org/records/21296348 and https://zenodo.org/records/21379836

**Models and implementations**
- Dafalias & Manzari (2004). JEM 130(6):622. https://doi.org/10.1061/(ASCE)0733-9399(2004)130:6(622)
- Boulanger & Ziotopoulou (2023). *PM4Sand v3.3*, UCD/CGM-23/01. Read §2.5, §2.7, §3.2, §3.6:
  https://itasca-software.s3.amazonaws.com/udm-library/Boulanger_Ziotopoulou_PM4Sand_v3.3_CGM-23-01.pdf
- SoilModels SANISAND page: https://soilmodels.com/sanisand/
- SoilModels UMAT/PLAXIS download: https://soilmodels.com/download/plaxis-umat-sanisand/
- SoilModels update note: https://soilmodels.com/upated-version-of-sanisand-umat-including-plaxis-interface/
- SoilModels forum thread: https://soilmodels.com/subroutine-of-sanisand-model/
- SoilModels FLAC3D download: https://soilmodels.com/download/sanisand-flac3d-download/
- Itasca DM04 UDM: https://www.itascainternational.com/software/udm-library/dafalais-manzari
- IC MAGE M11: https://zenodo.org/records/7765065
- IC MAGE UMIP: https://zenodo.org/records/13890543
- numgeo SANISAND benchmark: https://j-machacek.github.io/numgeo/2026-09/benchmark/constitutive-models/sanisand.html
- numgeo SANISAND-MSf reference: https://j-machacek.github.io/numgeo/2026-09/reference/material/mechanical/sanisand_msf.html
- cgl-sd MATLAB port: https://github.com/cgl-sd/SANISand04
- OpenSees SAniSandMS documentation: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/SAniSandMS.html
- Niemunis Incremental Driver: https://soilmodels.com/idriver/
- Mašín `triax`: https://soilmodels.com/triax/
- Hypoplastic UMAT integration (RKF-23): https://link.springer.com/article/10.1007/s11440-018-0684-z

**Local (fork, read-only)**
- `Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md` and its `README.md` and CSVs
- `Ladruno_implementation/127_tims_2d_requests_plan.md`
- `Ladruno_implementation/86_ladruno_sanisand_tims_report.md`
- `Ladruno_implementation/92_ladruno_sanisand_implex_adr.md`
- `Ladruno_implementation/93_ladruno_sanisand_zero_confinement_adr.md`
- `SRC/material/nD/UWmaterials/ManzariDafalias.cpp`: `integrate`, `explicit_integrator`, `ModifiedEuler`,
  `Stress_Correction`, `GetStateDependent`, `GetF`, `g`
