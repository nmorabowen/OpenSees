---
title: "ADR 92 — IMPL-EX integration for LadrunoSANISAND: a stress update that cannot seize, and a tangent that cannot lose ellipticity"
project: Ladruno
type: ADR
status: "ACCEPTED — P0 complete, CP1 measured, D0 DISCHARGED 2026-09-05; P1 (C++) OPEN on the CP1 ladder warrant"
priority: high
owner: nmora
orchestrator: "Opus 5 session `wp/92-sanisand-implex` (owns the context end to end)"
requested_by: "TIMs Workbench, act `work/ape/response-curve-matrix` (gupi claude-code), 2026-09-05"
related:
  - "[[_adr92_tims_request_2026-09-05]]"
  - "[[_adr92_sanisand_implex_scoping]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[_adr90_tau0_qu_band]]"
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[31_ladruno_concrete3d_adr]]"
  - "[[75_ladruno_sparse_direct_strategy_adr]]"
  - "[[LEDGER_implementations]]"
tags: [adr, material, nd-material, soil, sanisand, manzari-dafalias, implex, integration, tims]
aliases: [ADR-92, "LadrunoSANISAND -implex"]
updated: 2026-09-05
---

# ADR 92 — IMPL-EX for `LadrunoSANISAND`

> [!warning] Status — **PROPOSED. No C++ is open.**
> Numbered **92, not 91**: ADR 91 is the shell stiffness modifiers (`wp/91-shell-modifiers`,
> C++ already written). The consumer's request was written as "ADR 91" and is preserved
> unedited at `[[_adr92_tims_request_2026-09-05]]`; the fork-side review that corrects six
> of its statements about this code base is `[[_adr92_sanisand_implex_scoping]]`.
>
> **D0 gates the C++ on PR #792's T8.** WP-86b is attacking the same seizure from the cheap
> end (a substep-count cap that lets `ModifiedEuler` fail instead of force-accepting at
> `dT_min`, plus the consistent-tangent default). Until its GATE U re-run reports, the size
> of the problem IMPL-EX is being asked to solve is unmeasured.
>
> **P0 is COMPLETE** (`[[_adr92_p0_oracle_results]]`, 2026-09-05): G0 reproduces the binary
> to round-off on the deck-default paths; **D1 = A by 18–22x**; **D3 is reversed** — scheme
> 2's low-p Newton is disabled at `:2264` and 58–74 % of its calls at the corner are
> `ModifiedEuler` in disguise; and **the extrapolated stress crosses the `p_min` floor into
> tension** (`min p = −1.37 kPa` on a `+0.0101` state), so P1 gains a clamp. The seizure
> mechanism reproduces at one Gauss point only at low confinement x 10x the nominal
> increment (`p0 = 5`, `1e-3`: 163 force-accepts, implicit `q/p -> 0`) — and IMPL-EX-A is
> already broken one row earlier there (`5e-4`), so **`-implexControl` is not optional at the
> corner** and the cost claim in §2 stays unmeasured until #792 T8 / P3. **CP1 is next.**

---

## 1. Driver

`_adr90_tau0_qu_band` (GATE U, #791) could not measure a collapse load on softening
SANISAND because **every leg seized inside the constitutive integrator**: 0 of 80 pinned
subdivisions used, terminal step 6400-25000x above the floor, longest single `analyze(1)`
**2056 s = 59 % of that leg's entire budget in one step**, with up to **125
state-determination passes** per step and a substepped return collapsing toward
`dT_min = 1e-6` (`ManzariDafalias.cpp:1320`). The controller had every resource it was
given and never used any of it.

The TIMs act adds a second, independent finding from the same week: the matrix-free routes
do not substitute. `LadrunoDynamicRelaxation` relaxes each hold to a residual below
`1e-6 kN` and still returns equilibria that depend on the damping law, the mass safety
factor and the hold size by factors of 2.6 and 27 — because every sweep commits, and
SANISAND's `alpha`, `z` and `e` are path variables that read a fictitious oscillation as a
load cycle. **A solver that finds an exact equilibrium on a wrong path is the most
dangerous instrument the campaign owns**, and the cure is C++ in the integrator, not a
flag (`DirectIntegrationAnalysis.cpp:259` commits unconditionally).

So the requirement is a stress update that (a) terminates in bounded work at every Gauss
point regardless of confinement, (b) hands the global solver an operator that stays
positive definite where the material has lost ellipticity, and (c) does not manufacture
history. IMPL-EX (Oliver, Huespe & Cante 2008) is the standard answer to all three, and
**the fork already ships a working implementation of it** in `ASDConcrete3DMaterial`
(`-implex`, `-implexControl`, `-implexAlpha`, `implexError` / `avgImplexError`), written
for softening concrete for the same reason.

## 2. The operator

Let `n` be the last committed step. Committed and available on the base class:
`mEpsilon_n`, `mEpsilonE_n` (so `eps_p(n) = mEpsilon_n - mEpsilonE_n`), `mSigma_n`,
`mAlpha_n`, `mFabric_n`, `mVoidRatio`, `mAlpha_in_n`.

**Extrapolated (the response handed to the element, every global iteration):**

```
f          = (dt_{n+1} / dt_n) * implexAlpha                     # 1.0 at a fixed increment
d_eps_p~   = f * d_eps_p(n)                                      # d_eps_p(n) = eps_p(n) - eps_p(n-1)
sigma~     = sigma_n + Ce(p_n) : ( (eps_{n+1} - eps_n) - d_eps_p~ )
d sigma~ / d eps_{n+1} = Ce(p_n)                                  # constant in the step
```

**Incremental, not total.** SANISAND's elasticity is hypoelastic — the code integrates
`dsigma = Ce(p) : deps_e` with the moduli at the committed stress (`elastic_integrator`
`:1008-1011`, `BackwardEuler_CPPM` `:2223-2226`) — so the stress must be advanced from
`sigma_n`, never rebuilt from a total elastic strain. *(A first draft of this section wrote
the total form `Ce(p_n):(eps_{n+1} - eps_p~)`; the Fable review of 2026-09-05 caught it.
That form discards the accumulated pressure-dependent history and its error does not
vanish as `dt -> 0`.)*

`Ce(p_n)` is SANISAND's pressure-dependent elastic operator **frozen at the committed mean
stress**. It is symmetric and positive definite wherever `p_n >= p_min`, which is
everywhere by construction of the clamp. Nothing on this path touches `mAlpha`, `mFabric`
or `mAlpha_in`, so no fabric is accumulated and no reversal is detected on an extrapolated
state.

**Companion, at `commitState` only:** run the true return from state `n` with the actual
strain increment, obtain `sigma(n+1)` and the full internal state, set
`d_eps_p(n+1) = eps_p(n+1) - eps_p(n)`, and record

```
implexError = || sigma~ - sigma_implicit || / ( || sigma_implicit || + P_atm * eps_norm )
```

**The un-primed step is exempt from `-implexControl` refusal.** The first plastic step
after a stage flip (`updateMaterialStage 1`) has no committed plastic history yet —
`d_eps_p(n) = eps_p(n) - eps_p(n-1)` is `0` by construction — so `implexError` measured
there is not an extrapolation error at all, it is the companion's own drift-correction
jump from wherever the elastic stage left the stress to the first plastic return. Refusing
on that number subdivides a step that has nothing wrong with its extrapolation (`f`
was never even exercised) purely because the *companion* moved. `-implexControl` therefore
does not enforce its tolerance on this one step; the error is still computed and recorded
(so it remains visible in `implexError`/`avgImplexError`), only the refusal is suppressed.
Every later step in the same stage is primed and refused normally. (Fixed by `afb95c40c`
after the BVP re-run gate showed the registered arm refusing on step 1 of every stage —
`implexError` 0.13-0.30 against `tol 0.05` — before any real extrapolation had occurred.)

**Why the plastic strain and not `dGamma`** (the request's choice, corrected — scoping §C2):
`mDGamma` is the step total only under `BackwardEuler_CPPM` (zeroed at `:2220`, solved as
`Delta(18)` at `:2274`). Under every
substepped explicit scheme it is *the last substep's* multiplier — `ForwardEuler`
reassigns it fresh at `:1342`, `ModifiedEuler` never accumulates it — and at the corner the
substep count swings by orders of magnitude between steps, so extrapolating it extrapolates
noise. The plastic-strain form is integrator-agnostic, needs **no virtual hook into
vanilla**, and inherits the guards for free.

**What this buys, in order of confidence:**

1. **Bounded work per step.** The expensive return runs **once per committed step** instead
   of once per state-determination pass — up to 125 of them at the seizure. That is a ~100x
   cut in constitutive work before any change to the cost of a single return, and it is the
   same lever #792's cap pulls from the other end. **Still unmeasured.** P0 shows the
   mechanism (163 force-accepts at `p0 = 5`, `d eps = 1e-3`) but not IMPL-EX surviving it —
   A breaks at `5e-4` at that confinement while the implicit is still 7 % off; the BVP re-run
   (`_adr92_p1_bvp_gate_rerun.md`) counts solver-ladder work, not constitutive-integrator
   work per step, so this claim is not touched by that data either — it is still a P3 / T8
   number and the campaign must not read the item-2 confirmation below as covering it.
2. **A global step that is linear.** Newton converges in one iteration on a frozen operator;
   the ladder never fires, and the "rung 3 commits states nothing afterwards converges from"
   pathology (ESMERALDA §30-31) cannot arise. **Confirmed at BVP level on the fixed binary,
   both arms** (`_adr92_p1_bvp_gate_rerun.md`): control-OFF, `-implex` alone, 142/142
   converged steps on rung 1, 0 subdivisions; the registered `-implex -implexControl 0.05
   0.01` arm, 504/504 converged steps on rung 1, 0 rung-2/3 — the 81 subdivided-and-abandoned
   attempts (`nfail = 243 = 3 x 81`) were all material refusals, none a
   `CTestNormUnbalance` failure, so no converged step ever left rung 1 on either arm.
3. **A symmetric global matrix** — the non-associated consistent tangent disappears. This
   unlocks `system Pardiso -matrixType sym` (ADR-75 P1d: 1.94-1.96x vs UmfPack, -42 % peak
   memory, exact) on the 21 058-DOF coarse and 175 290-DOF fine cells. **Unclaimed by the
   request; measured at P3 — on the drained `LadrunoBrick` legs only.** `LadrunoUP`'s u-p
   tangent is unsymmetric regardless of the material (ADR 71), so the `U-L` row keeps its
   general solver.
4. **The corner stops being a solver event — with one guard P0 found.** At Gauss points
   pinned at `p_min` the operator is the floor's operator — small, positive definite. But the
   extrapolated correction is **not** bounded by the floor: nothing clamps `sigma~`, and P0
   measured it crossing into tension (`min p = −1.37 / −0.16 / −0.09 kPa` at 40 / 80 / 160
   steps on a `+0.0101 kPa` state — first order, O(1)–O(10) relative). **P1 applies the
   code's own device to the extrapolated stress: `sigma~ = dev(sigma~) + p_min*I1` whenever
   `tr(sigma~)/3 < p_min`.** *(Review, 2026-09-05: that clamp repairs the isotropic part
   only; the runaway quantity is the unbounded `f·d_eps_p(n)`, which distorts the deviator
   by the same order. P0 decides between the clamp and **bounding `f`** so the whole
   `sigma~` stays admissible — the second is the deeper fix and acceptance 1b must then
   check the deviator too, not just `tr`.)* The ring still flows; what disappears is
   Newton's missing descent direction, not the mechanics.

### P2 addendum: trial-direction correction (variant B) — REJECTED

P2 asked whether extrapolating the flow **direction** at the elastic-trial stress
(`R_tr`, variant B), rather than freezing the whole committed `d_eps_p(n)` tensor (A,
as specified above), fixes the one place A is provably wrong: A's extrapolated
volumetric increment carries the wrong **sign** at every resolved phase-transformation
crossing (`_adr92_p2_direction_oracle.md`). The oracle (`--gate GE`, drained TX
compression, companion scheme 1) measured both variants against the implicit
companion at the same crossings:

| `p0` | `d eps_z` | A, `implexError` at crossing / path max | B, at crossing / path max | A/B (path max) |
|---|---|---|---|---|
| 100 kPa | 1e-5 | 2.80e-6 / 1.39e-3 | 4.00e-5 / 3.09e-2 | **0.045 (B is 22x worse)** |
| 5 kPa | 1e-5 | 9.22e-5 / 1.18e-2 | 1.04e-4 / 7.06e-2 | **0.17 (B is 6x worse)** |

**Kept A.** The sign error is real (§4 of the oracle memo: A, B and a third variant C
all get the crossing's volumetric sign wrong) but inconsequential — phase
transformation is *defined* by `D -> 0`, so the wrong-signed term is `O(1e-10)` against
a deviator `O(1e-5)`, and the crossing is the **quietest** step on the whole path (500x
below A's own path-max error). B does not repair the sign either (`alpha` is committed
in every variant, and the sign lives in `alpha`, not in where `R` is evaluated), misses
the terminal `q/p` by 9-54 % where A matches the implicit companion to four figures, and
forfeits ADR §2 benefit #2 above: B's true tangent differs from the frozen `Ce(p_n)` the
element is actually handed by 44-229 % of `max|Ce|`, so a "linear" B step is 2-78 %
non-linear in truth. Reopens only on BVP evidence a triaxial-ramp oracle cannot give:
a Gauss point whose *committed* `D` oscillates in sign step to step, not a single
monotone crossing.

**And the cost, stated first because the campaign must print it (scoping §C5):**
IMPL-EX is a first-order-in-`dt` perturbation of the constitutive response with the
structure of an artificial viscosity. That is *why* it robustifies softening. **An IMPL-EX
leg is therefore a regularized leg with the step size as the regularization parameter, and
that parameter has no length in it.** Every width, band and post-peak branch read off one
is regularized by `dt`. The request's line disclaiming ADR 90 is wrong in this direction;
this ADR carries the disclosure instead, and ADR 90's reopened D2 gains a second candidate
regularizer that is already half-built.

## 3. Where it goes, corrected

`LadrunoSANISAND` subclasses `ManzariDafalias` (ADR 86; ND tags 33019 / 33020 / 33021).
The base's scheme map, verified at `ManzariDafalias.cpp:40-49` and the dispatch at
`:984-993` / `:1031-1057` — **the request had 1 and 2 swapped and mislabelled 0 and 3**:

| `mScheme` | integrator | kind |
|---|---|---|
| 0 / 4 / 6 | `MaxEnergyInc` | explicit, substepped |
| **1** | **`ModifiedEuler`**, error-controlled substepping, `dT_min = 1e-6` | **explicit — the deck default (`:93`) and the one that seizes** |
| **2** | **`BackwardEuler_CPPM`** | **implicit — the only one** |
| 3 / 5 | `RungeKutta4` / `ForwardEuler` | explicit, no error control |
| 7 / 8 / 9 | `MaxStrainInc` | explicit, substepped |
| 45 | `RungeKutta45` (Abell) | explicit, error control |

Three seams the request does not mention and the implementation must honour:

- **Stage 0 is elastic.** `updateMaterialStage -stage 0` sets `mElastFlag = 0` and
  `integrate()` takes the `elastic_integrator` branch (`:978`). `-implex` must be **inert**
  during gravity and during the `LoadControl 0.0` re-equilibration, and must initialise its
  history at the stage switch, not before it.
- **Elastic moduli must be frozen at `p_n`.** The base evaluates `G`, `K` on the current
  stress; left alone, the delivered tangent is not the operator the stress was built with.
  Test it as an identity to machine precision.
- **The refusal contract already exists.** ASD returns a bare
  `EC_IMPLEX_Error_Control = -10` (`ASDConcrete3DMaterial.cpp:59-61`, `:1679-1684`). The
  fork settled that question the other way three commits ago in #792: the sentinel
  `LADRUNO_MATERIAL_REFUSED (-33086)` in `SRC/material/LadrunoMaterialStatus.h`, propagated
  **only by exact value**, with a process-budgeted report. `-implexControl` uses it, and
  then `-maxSubsteps` and `-implexControl` share one refusal path.

## 4. Deck syntax

```
nDMaterial LadrunoSANISAND $tag  <23 constants>  \
    -Presidual $pr -Pmin $pmin -honorTolR $h -maxSubsteps $N \
    -implex  <-implexControl $tol $reductionLimit>  <-implexAlpha $a> \
    <-implexDt pseudo|strain|user>
```

`-implex` off (default) is **byte-identical** to today. `implexError` and `avgImplexError`
join the material responses on the `ASDConcrete3DMaterial.cpp:2073-2077` template.

## 5. Decisions

| | decision | resolution |
|---|---|---|
| **D0** | Sequence against WP-86b (#792) | **P0 opens now; C++ (P1+) is GATED on #792 T8**, the GATE U re-run with the substep cap and the consistent-tangent default. Owner checkpoint **CP1**. |
| **D1** | Extrapolated history variable | **Plastic strain** (§2). `dGamma` + frozen flow direction stays on the table as the textbook alternative and is measured against it at P0; it would require the companion to be scheme 2 and a vanilla hook. |
| **D2** | Time source | `-implexDt {pseudo\|strain\|user}`, default `pseudo` = `ops_Dt` (ASD behaviour). Correct for the TIMs deck **including under subdivision**, because it drives settlement by `LoadControl` on a prescribed-settlement SP pattern so pseudo-time is proportional to settlement. Guarded at `dt = 0` (holds, stage switch); **refused** on integrators that solve for the load factor (`DisplacementControl`, arc length), where `dt = d(lambda)` is not proportional to the increment and changes sign past a limit point. **Why refuse rather than clamp:** the June `LadrunoRCConcrete` entry in `LEDGER_quirks` (§ "IMPL-EX in a STATIC analysis") proved a clamp-and-degrade fix (`tf` falls back to `alpha`, capped at `2·alpha`) for a *cyclic wall* whose static steps are meant to be uniform. Here the ratio is the extrapolation itself — a wrong `dt` is a wrong answer that passes every gate — so SANISAND refuses where the ratio cannot be trusted and degrades only where it can (`dt = 0` holds).
**The guard is on a SIGN CHANGE, not on `dt > 0` — a monotone negative clock is legal.**
The original C++ (`3c788778f`'s predecessor) gated the extrapolation factor
`f = dt_{n+1}/dt_n * implexAlpha` on `mImplexDtCommit > 0.0`, so on any `LoadControl(-ds)`
leg (settlement driven by a negative pseudo-time, which the campaign's decks all use) the
ratio was never computed and `f` silently froze at `1.0` for the leg's entire life — found
by the red/blue review (B1) reconstructing 9/142 steps with a true ratio != 1, all run at
`f = 1`. Fixed by `2473ce46c`: the gate at `LadrunoSANISAND.cpp:1329` is now
`mImplexDtCommit != 0.0`, and the factor is the **sign-consistent** ratio
`dt_{n+1}/dt_n` — two increments of the same sign give the same positive ratio a
monotone-positive clock would, and a sign **change** (a limit point, or a step that
crosses back through a hold) is what gets refused, not negativity itself. |
| **D3** | Companion integrator | **REVERSED by P0.** Default **scheme 1 (`ModifiedEuler`) with `-maxSubsteps` required** (#792 T1) so the companion cannot seize. Scheme 2 is *permitted* but its low-p Newton is disabled (`:2264`, literal `errFlag = 0`) and on the corner path 58–74 % of its calls integrate by `explicit_integrator` — it is not an implicit return where the campaign's problem lives, and it costs a 19-unknown Newton everywhere else. Schemes 3 / 5 / 7 / 8 / 9 refused with a sentence (5 additionally carries the zero-`r` defect). |
| **D4** | Class tags | **None new.** Flags on 33019 / 33020 / 33021. The wire format grows: `sendSelf` / `recvSelf` / both `getCopy` forms carry the flags and `d_eps_p`, per the ADR-86 six-override rule and the FSPM `getCopy` lesson. |
| **D5** | Vanilla footprint | **Zero** under D1. If P0 overturns D1, one flag seam in `ManzariDafalias.h` on the `mHonorTolRInME` / `mMaxSubstepsInME` pattern plus a `LEDGER_vanilla_files` row. |
| **D6** | Relationship to ADR 90 | Cross-reference, not firewall. §2's disclosure is a **P0 deliverable**, not a P3 one, because the act will need it the first time it quotes an IMPL-EX curve. |
| **D7** | Number / branch | **ADR 92**, branch `wp/92-sanisand-implex`. The stray untracked `91_ladruno_sanisand_implex_adr.md` in the main checkout is superseded by this file and should be deleted by the owner. |

## 6. Plan

Each phase ends in a written artefact; the owner decides at CP1 and CP2. **The orchestrator
holds the context across all phases** and briefs every delegated agent from it; no agent is
given the campaign to re-derive.

| phase | deliverable | gate to leave it |
|---|---|---|
| **P0** *(COMPLETE)* | Single-Gauss-point IMPL-EX oracle in numpy on the D-L cell's 23 constants, driven by the act's Level-0 drained triaxial path. Measures: first-order convergence of the IMPL-EX/implicit difference under increment halving; the error at the act's `1e-4 m` increment; behaviour as `p' -> p_min`; **D1 vs the `dGamma` form, head to head**. Plus the §2 disclosure text. | The convergence exponent is measured (not asserted), D1 is decided on numbers, and the error at the campaign's increment is known. Artefact `_adr92_p0_oracle_results.md`. |
| **CP1** | Owner checkpoint. Reads P0 **and** #792 T8 together. | Is IMPL-EX unblocking, or an optimisation? Phasing and risk budget are set here. |
| **P1** | `-implex` on `LadrunoSANISAND3D`: extrapolated stress **(incremental, clamped at `p_min`)**, frozen `Ce(p_n)`, companion at `commitState` **(scheme 1 + `-maxSubsteps`)**, `implexError` / `avgImplexError`, `getCopy` / `sendSelf` / `recvSelf`, stage-0 inertness, **and `-implexControl`** — moved up from P2 because P0 measured A unusable from `d eps = 5e-4` at `p0 = 5 kPa`, and the corner is where this campaign lives, so the control is not optional there. **The P0 oracle is the reference: `implex_A` on the binary's own paths to 1e-8.** | `-implex` unset is **byte-identical** on every existing SANISAND deck; tangent identity to machine precision; Zone-A green. |
| **P2** | `LadrunoSANISANDPlaneStrain` (the 2D act needs it on the strip); the cyclic/reversal test. **`-implexControl` moved up into P1** (P0 measured A unusable from `d eps = 5e-4` at `p0 = 5 kPa`, and the corner is where this campaign lives). | Reversal test green; the plane-strain lane carries the flags. |
| **P3** | Esmeralda: the corner patch; the coarse bare `D-L` leg against job 146299's wall at `s/B = 0.0206`; the `U-L` coupled row inside `LadrunoUP`; **the symmetric-solver measurement of §2.3**. | A WP1 plateau, or a named reason there is none — plus `implexError` reported beside every verdict. |
| **CP2** | Close-out. | Ledgers, banner row if shipped, ADR status. |

**Not in scope:** IMPL-EX on vanilla `ManzariDafalias`; any change to `LadrunoUP`; ADR 90's
rate regularization, which addresses localization width and is a different instrument.

## 7. Acceptance

The request's §5 list survives, remapped: its tests 1, 3, 6, 7 -> P0/P1; test 2 -> P2;
tests 4, 5 -> P3. Added by this ADR:

1. **Tangent identity.** Returned tangent == numerical `d sigma~ / d eps`, machine precision
   (P0 measured `3.5e-11` on the oracle; the C++ must match).
1b. **Floor clamp.** On the P0 G3 path, `tr(sigma~)/3 >= p_min` at every iteration of every
   step; the oracle without the clamp reaches `−1.37 kPa` and is the negative control.
1c. **Oracle parity.** `implex_A` in the C++ reproduces the P0 oracle's `implex_A` on the
   recorded binary paths to `1e-8` — the same G0 discipline, one level up.
7. **The BVP gate — the one that can refute the cost case.** On the CP1 deck
   (`h1.0_e0.6944`, `Q = 10`, cap 1000), with `-implex -implexControl`, measure the ladder
   decomposition CP1 measured without it. **Prediction: steps past rung 1 fall from 61 % to
   near zero and the failed-rung share of iterations falls from 89 % to single digits.** If
   the ladder still fires at CP1's rate, IMPL-EX's cost case is refuted at BVP level and P1
   closes as a correctness-only feature. Report the same table either way.
2. **Stage-0 inertness.** Gravity and the `LoadControl 0.0` re-equilibration are bit-identical
   with `-implex` on and off.
3. **`dt` guards.** A `dt = 0` step and a halved step both extrapolate by the right factor;
   `DisplacementControl` is refused with a sentence.
4. **Refusal.** `-implexControl` past tolerance returns `LADRUNO_MATERIAL_REFUSED`, the
   element propagates only that value, subdivision engages, the committed state is unchanged.
5. **Symmetry.** Under `-implex` the assembled tangent is symmetric to round-off, and
   `system Pardiso -matrixType sym` reproduces the general solver's answer.
6. **Vanilla untouched.** `ManzariDafalias` decks bit-identical (the ADR-86 gate).
7. **BVP ladder-removal gate, registered arm.** Not in the original list; added once the
   fixed binary (`2473ce46c` + `afb95c40c`) made a same-binary registered-arm run possible.
   `_adr92_p1_bvp_gate_rerun.md`'s status line, quoted verbatim: "All three arms COMPLETE.
   Registered arm (-implex -implexControl 0.05 0.01, build afb95c40c): gate
   --registered-arm VERDICT = PARTIAL -- 0.0% past rung 1 (converged-only) but 13.8% on
   attempts (81/585), 99.0% failed-rung iterations; terminated BUDGET at s/B 0.02754.
   Control-OFF arm (build 2473ce46c): VERDICT = PREDICTION MET (0.0%/0.0%), as before."
   PARTIAL, not PASS: read alongside item 2 above, not as a substitute for it — the
   converged-only 0.0 % is real (no converged step left rung 1 on the registered arm
   either) but the 13.8 %/99.0 % pair records that the ladder still fires and burns
   through all three rungs before every one of the 81 abandoned attempts is subdivided
   away, exactly the §8 risk below stated it would.
8. **Mutation gate (ADR-87 D2).** `_adr92_p1_mutation_gate.md`: **PASSED, score 0.750
   (9 of 12 hand-mutants killed)** against the 0.60 floor, run on
   `tests/test_ladruno_sanisand_implex.py` (20 baseline-passing detectors). Three
   mutants survived and are owed tests, not waived: **M4** (the `-implexControl`
   reduction floor, `mImplexDt0`, is never armed by any deck in the battery — no test
   drives a subdivision ladder), **M5** (a refused trial returns a bare `0` instead of
   the `LADRUNO_MATERIAL_REFUSED` sentinel — the battery pins the refusal's symptoms
   but not its return-code contract), and **M10** (the re-arm-after-refusal line is
   redundant only because every test's failed step goes through
   `Domain::revertToLastCommit()`, which re-arms anyway — a caller that retries
   without reverting is untested). See `_adr92_p1_mutation_gate.md` §4 for the full
   audit of each.

### `-implexControl` operating point (measured 2026-09-06)

The registered arm's tolerance (`tol = 0.05`, `reductionLimit = 0.01`) was swept
against three looser tolerances and a 10x-looser floor, same deck and build
(`afb95c40c9`), reference `control` at `s/B = 0.0678` (WALL) —
`_adr92_p1_bvp_gate_rerun.md`'s "Operating-point sweep" section:

| tol / rLimit | mode | s/B (depth) | mean overlay dev % (excl. step 1) |
|---|---|---|---|
| 0.05 / 0.01 (registered) | BUDGET | 0.028 | 2.10 |
| 0.05 / 0.1 | BUDGET | 0.028 | 2.10 |
| **0.1** / 0.01 | BUDGET | **0.076** | 1.87 |
| 0.2 / 0.01 | BUDGET | 0.150 | 1.91 |
| 0.5 / 0.01 | TARGET | 0.250 | 2.29 |

**The registered `0.05` fails on reach, not accuracy** — it never gets past
`s/B = 0.028` (BUDGET, the substep cap, not a bad extrapolation) against `control`'s
own `0.0678`, while every tolerance in the sweep, `0.05` included, tracks `control`
to a 1.9–2.3 % mean overlay deviation once the shared step-1 elastic-predictor
outlier is excluded. **`0.1` is the tightest tolerance tested that beats `control`'s
depth while staying under a 5 % mean deviation** (`0.076` vs `0.0678`, at `1.87 %`).
**Decided:** default `0.1` since WP-92d (the owner's ready-flip on this
recommendation) -- the C++ struct default (`LadrunoImplexOptions::errorTol`)
and the documented deck default in the P1/P2 guide and any campaign driver now
both read `0.1`. Separately: `reductionLimit`
as defined (a floor relative to the **first** increment) is inert at `tol = 0.05` on
this deck — the floor sits two orders below the working step size at depth and
never gets a chance to bind before the `tol` criterion already refuses — and should
be re-based on the **current** step's nominal increment in P2, not the first one.

## 8. Risks

- **The companion sees a different strain path.** Under `-implex` the global step is solved
  on the elastic operator, so the strain increments handed to the commit-time return can be
  larger and differently directed than the ones implicit Newton would have found. The
  companion may therefore be *harder* per step even though it runs far less often; P0
  cannot see this (single Gauss point, prescribed strain) — it is a P3 measurement.
- **`-implexControl` needs an in-step implicit or a one-step lag.** ASD computes the implicit
  solution on every `setTrialStrain` when control is on (`:1665-1684`), which is the cost
  IMPL-EX was meant to remove. The P2 design must choose: pay it on the (now 1-2)
  iterations, or check the error at commit and shrink the *next* step (a-posteriori,
  one-step lag, cannot refuse the step it measured).
- **First-order lag — priced.** At the nominal campaign increment (`1e-4`) the extrapolation
  costs `5e-5` in stress and `1e-3` in `eta`, below the substepper's own error; at `5e-4` it is
  `4e-3 / 3e-2`, at `1e-3` `3e-2 / 0.13` (P0 §4). The corner Gauss point sees 10–100x the
  nominal, so `-implexControl` at `0.05` will halve the step exactly where the wall was.
  Correct behaviour, bounded by the reduction limit — IMPL-EX trades the wall for a cost
  there rather than removing it. **At `p0 = 5 kPa` A is unusable from `5e-4` up** (P0 §4,
  `q/p` 0.09 vs 2.07), so at the corner `-implexControl` is a requirement, not an option.
  **Confirmed at BVP level:** the registered arm (`-implexControl 0.05 0.01`) reaches
  only `s/B = 0.02754` before hitting its subdivision budget (80, one leg over at 81),
  with `n_material_refused = 10270` over 504 converged steps — against the same deck's
  unpoliced `-implex` arm, which reaches the full `s/B = 0.25` target on zero refusals.
  The control does exactly what this paragraph predicted: it does not remove the wall,
  it relocates it, trading depth for a bounded, counted refusal cost instead of an
  unbounded ladder (`_adr92_p1_bvp_gate_rerun.md`).
- **Cyclic response lags by one step.** `alpha`, `z` and `alpha_in` advance only on the
  committed path. Monotonic pushover is the target; cyclic use needs the P2 reversal test
  before it is claimed.
- **Two tangents in one mesh.** A Drucker-Prager crust beside `-implex` sand: the global
  matrix is positive definite only if every material's is, and the crust's consistent tangent
  is not always so. Run the crust elastic-tangent, or the crust implicit with rungs 0-1 kept.
- **The floor overshoot** (P0 §5): without the P1 clamp, the free-surface ring receives
  tensile mean stress every iteration. With it, the ring is still the least accurate part of
  an IMPL-EX field (first order, O(1)–O(10) relative at the floor) and the disclosure says so.
- **A reading hazard, and it is the serious one.** An IMPL-EX curve satisfies equilibrium
  with the *extrapolated* stress. A limit point called on an IMPL-EX leg must be confirmed by
  the implicit material up to the last settlement the implicit solver reaches, and
  `implexError` must be printed beside every WP1 verdict. This is the same failure mode as
  the DR leg that passed every gate on a wrong path — the fork must not hand the campaign a
  second one.

## 9. Staffing — the orchestrator owns the context

**This session (Opus 5, `wp/92-sanisand-implex`) is the orchestrator and holds the context
end to end**: the GATE U evidence, the TIMs act's five days, #792's state, and the six
corrections. Delegated agents receive a self-contained brief written from that context and
return a measurement or a diff; **none is given the campaign to re-derive**, and none owns a
decision.

| work item | agent | model / effort | why this level |
|---|---|---|---|
| ADR, scoping, decisions, all owner-facing reading | **orchestrator (this session)** | Opus 5, high | The context is the deliverable; it does not survive a handoff. |
| **P0** oracle — IMPL-EX vs implicit at one Gauss point, D1 head to head | `general-purpose` | **Opus, high** | Constitutive algebra where a plausible-looking wrong answer is the failure mode, and it decides D1. The ADR-90 WP-A pattern. |
| P0 source read-back — confirm the oracle's `Ce`, flow and clamp against `ManzariDafalias.cpp` | `Explore` | **Sonnet, medium** | Bounded read-only lookup with exact file:line answers; cheap, and the orchestrator checks it against what it already knows. |
| **P1** C++ — the flag, the operator, the six overrides | `general-purpose` in this worktree | **Opus, high** | Touches the wire and `getCopy`, where ADR-86 has already been bitten once. |
| P1 / P2 tests — `tests/test_ladruno_sanisand_implex.py` | `general-purpose` | **Sonnet, medium** | The fork's rule: the author does not write the tests. Acceptance list is fully specified in §7, so judgment is not the binding resource. |
| Adversarial review of §2's guards and §8's reading hazard | `general-purpose` | **Fable, one pass** | One cheap independent pass on the two places a silent wrong answer would hide. |
| `/code-review high` on the P1 head | slash command | — | The standing gate. |
| Builds (`build.bat OpenSees OpenSeesPy`), Esmeralda submission, merge | **the human** | — | ADR-87: the owner merges; agents do not. Builds are launched from a real terminal (see the build-launch traps). |

## 10. Ledger obligations

- `LEDGER_implementations` — one row: `LadrunoSANISAND -implex`, ADR 92, files
  `LadrunoSANISAND{,3D,PlaneStrain}.{h,cpp}`, **no new classTag** (D4).
- `LEDGER_vanilla_files` — **no row expected** under D1 (D5). If P0 overturns D1, one row for
  the `ManzariDafalias.h` seam.
- `LEDGER_quirks` — the `mDGamma` finding of §2 (it is not the step total under substepped
  schemes) is a fork-wide gotcha and is owed a quirks entry regardless of whether this ADR
  ships.
- Banner: a `shipped` row only at CP2.

## Log

- **2026-09-05 (later)** — P0 complete (`[[_adr92_p0_oracle_results]]`): G0 PASS to
  round-off, G1 order 1.7–2.1, G5 `5.7e-11`; **D1 = A** (18–22x over `dGamma` on scheme 1);
  **D3 reversed** (scheme 2 is explicit at low `p`, `:2264`); **floor clamp added to P1**;
  seizure not reproducible at one Gauss point, cost claim stays for T8/P3. The Fable review
  caught the total-vs-incremental stress form before the oracle measured it (78 % at
  `p0 = 5`). Builder terminated by its session limit at G0; Fable carried the gates.
- **2026-09-05** — Requested by the TIMs act after the strict-ladder diagnostic (ESMERALDA
  §31) located the wall at the footing's corner and the relaxation route was found to load
  SANISAND's memory (§33). Fork-side scoping found six corrections, a live overlapping work
  package (#792) and an ADR-number collision; renumbered 91 -> 92, D0-D7 taken, P0 opened,
  C++ gated on #792 T8. Written on `wp/92-sanisand-implex` cut from `ladruno` at `3f003d110`.

## P2 (owed) — what the Esmeralda census and the fork-side probes found, 2026-09-07

Evidence: `_adr93_seat_replay.md` (gate met, error reproduced), ADR 93 Log 2026-09-06/07,
`LEDGER_quirks.md` rows of 2026-09-07. Four defects/limits, each with acceptance data.

| # | item | mechanism | fix shape | acceptance | status |
|---|---|---|---|---|---|
| P2-1 | **The floor branch commits an O(1) error** | `-implexControl`'s "nothing left to cut ⇒ accept" commits whatever error remains; the committed state is then out of equilibrium by O(1); the next linear solve closes the gap with a ds-independent strain (ring: 2–9× per ds, relaxing over ~15 rows, growing on repeats); the control refuses to the floor again ⇒ self-sustaining loop | at the floor, that Gauss point delivers the **implicit** stress for the step under the frozen `Ce` (SPD kept, +1–2 iterations, no O(1) commit possible); alternative: refuse at the floor (honest wall) | Esmeralda loose 146456 rows 1101–1135: the strain-per-ds excursions must vanish; dense 146457 must walk past 0.0178 or stop honestly | built in #807 (87b9cf846), acceptance pending Esmeralda |
| P2-2 | **Extrapolation at a softening / reversal point** | at `Kp ≤ 0` (post-peak dense) with an `α_in` reset between steps, the previous plastic increment is the wrong thing to extrapolate: seat 4095/8, error 0.46 reproduced; alpha 0.5 → 0.22, direction variant → no change, **f = 0 → 0.029 (under tol)** | `f = 0` on a step whose committed predecessor showed `Kp ≤ 0` or an `α_in` reset (elastic predictor there); report the count as a new `implexRefusals`-style census | the seat row's error under tol without refusal; the P0 oracle rows unchanged elsewhere (byte-identity where `Kp > 0`) | built in #807 (87b9cf846), acceptance pending Esmeralda |
| P2-3 | **A zero-dt hold corrupts the next step** | a zero-increment commit stores `dt_n = 0` ⇒ next `f = alpha` with a zero history; RED-1 F9 was wrongly downgraded; NOTE the 4 / 21 / 29 % curve divergence first attributed to this was later shown to be harness-level on BOTH materials (reads inert, holds perturb the implicit twin too, Esmeralda 146459/146460), so P2-3 stands on the fallback mechanism alone | on a zero-increment commit keep the previous `dt_n` and `Δε_p(n)` | a hold inside a push is byte-inert on the following steps; guide warning until then | built in #807 (87b9cf846), acceptance pending Esmeralda |
| P2-4 | **`setParameter stressCorrection` is a no-op** | `ManzariDafalias::updateParameter` reads `theInt` (`:897`, `:864`), interpreters set `theDouble` | `LadrunoSANISAND::updateParameter` override accepting `theDouble`; zero vanilla footprint | a positive control: the first-step `Q` moves when the flag is set | built in #807 (87b9cf846), acceptance pending Esmeralda |
| P2-5 | **No loading-reversal reset on a round-off strain increment** | vanilla `ManzariDafalias::integrate()` `:1005-1013` dots `(α_n − α_in_n)` with `Ce·Δε` and resets `α_in := α_n` on a negative sign with no magnitude guard, so on a zero-increment step (hold) the sign is round-off noise and the reset fired at 28–54 % of 34 560 points on Esmeralda 146458, sending `h → ∞` and stiffening the implicit column 2.5× for tens of steps | subclass guard `‖Δε‖ < -reversalTol` (default 1e-10) restores `α_in_n`, counted in `implexGuards[3]` Extended (Esmeralda 146580): the same relative gate applies to the P2-2 guard's reversal/softening flags — with `alpha_in` resets down sixfold the post-hold jump stayed 1.9× because the guard's flag was set from the hold's noise and the next step ran the elastic predictor at ~1000 extra points (`f=0` count 340 → 986, committed error 0.015 → 0.071); a hold-preserved commit sets no guard flag; at `-reversalTol 1e-7` the guard fires ~10 000 point-calls per ordinary 1e-4 m step on a 34k-point deck (stated, harmless). | a hold inside an implicit push must leave every point's `α_in` unchanged and the load increment per step unchanged across the hold | built in #807 (8bfdfbc17); hold acceptance pending Esmeralda (alpha_in-equality count after a 0.005 hold must not exceed a monotonic push's 9642/34560; per-step load increment unchanged across the hold; census curve overlays the plain implicit twin); superseded by P2-5b for the threshold; the mechanism stands |
|P2-5c|**a zero pseudo-time increment is a hold: no loading-reversal test and no guard flags at any point, both paths.** Fork probe on d5bd259f6 (R3 footing, 1600 GPs): P2-5b is the active ingredient (the `-reversalRel 0` mutant reproduces the absolute-only rate) and cleans the guard-flag channel (post-hold `f=0` delta 0, error at baseline), but `alpha_in` still reset at 136/1600 (IMPL-EX) and 88/1600 (implicit) at points whose own last increment was tiny. A hold is a global fact the material can see: `ops_Dt == 0.0` on a `LoadControl(0.0)` step and on `DisplacementControl`'s dλ = 0 (P2-3 already keys on it)|skip the reversal test and set no guard flag on `ops_Dt == 0` at every point, both paths; no new flag (echo states the rule); per-call skips counted in `implexGuards[3]`, and **hold-skip commits once per point per hold in a new `implexGuards[5]`** (Vector(6)) so a harness can assert every point took the skip. Implicit-side band edge (Esmeralda 146574, absolute default 1e-10): 7802/34 560 resets after one hold and a subdivision cascade to 1.6e-6 m with a 1.7× tangent — the absolute default did almost nothing on the implicit path|built in #807 (d30c66582); **fork probe PASSED** — `alpha_in` changed 0/1600 on both arms (and on the `-reversalRel 0` mutant: the hold rule is unconditional), `implexGuards[5]` +1600 per hold, post-hold guard channel at baseline; Esmeralda census rerun pending. Acceptance was: fork probe `alpha_in` changed = 0/1600 on both arms at default flags; Esmeralda census pair same-leg before/after count unchanged across each hold and dQ/ds unchanged||| **Esmeralda hold acceptance on 887fea475 (IMPL-EX census 146593, holds at step 1 / 0.0001 / 0.005): PASS** — `implexGuards[5]` +34 560 at every hold, `alpha_in` equalities unchanged across the holds (28 629 → 28 267 at step 1 is round-off in `alpha`), dQ/ds 71.2k before / 71.2k after the 0.005 hold, and the census leg sits on the plain leg's own curve (71.4–72.1k): holds are free. Implicit census pending. |
| P2-5b | **Reversal-noise guard relative to the last committed increment** | P2-5's absolute threshold (1e-10) cannot work because a hold's per-point strain increment is Newton-tolerance-scale noise, measured on the fork's R3 footing (1600 GPs, `708152eac`): median 4e-9, max 6.4e-8 (IMPL-EX) / 1.4e-6 (implicit); at 1e-10 `alpha_in` still reset at 42 % / 9.5 % of points on a hold, at 1e-7 still 2.2 % on the implicit arm; Esmeralda's hold acceptance on `8bfdfbc17` failed the same way (1702/34 560, tangent jump 3.7× → 2.2×) | `‖Δε‖ < max(reversalTol, reversalRel·‖Δε_lastCommitted‖)`, `-reversalRel` default 0.05 (a hold ≤ 1e-2 of the previous step; a genuine reversal ~1×; a halved retry 0.5×), pre-hold reference kept across a zero-increment commit | hold probe `alpha_in` changed = 0 on both arms, Esmeralda census legs' same-leg before/after count unchanged across each hold and dQ/ds unchanged | built in #807 (d5bd259f6). **Esmeralda 146585 (IMPL-EX census, default flags): PASSES at the 0.005 hold** — `alpha_in == alpha` 0 → 0 of 34 560, dQ/ds 69.5k before / 69.4k after, committed error 1e-3 unchanged, `implexGuards[5]`-equivalent hold-preserved +34 560; the two EARLY holds (step 1: 0 → 23 796; s/B 0.0001: 16 498 → 30 219) still reset because the points' own previous increments are tiny, costing 2.5 % softer through 0.005 — the P2-5c case |
|P2-7|**the stage flip absorbs the drift correction immediately.** the IMPLICIT path is 23–33 % soft from step 1 on every engine since P2-5 (Esmeralda 146574/146586: first rows 1.296 / 2.6 kN vs 6.511 on c162833ed; step 1 subdivides). Source: the elastic stage sets `α = dev(σ)/p` every step (`ManzariDafalias.cpp:1055`) so at the flip `α = r_gravity ≠ 0`, `α_in = 0`; `Elastic2Plastic()` never touches `α_in` (init lines `:5139-5141` commented out); the only initialisation is the sign test in `integrate()` `:1005-1013` on the first plastic evaluation, which on a zero-increment re-equilibration is round-off noise — vanilla initialised `α_in := α` at roughly half the points by chance; P2-5/5b/5c turned that into never, so `(α − α_in):n` starts large and `Kp` small. The first P2-7 attempt (hold-style skip at the flip, 691f4064d) would have done the same and was not built.|at the flip `α_in := α` deterministically at every point on both paths (`-flipAlphaIn init|vanilla`, default `vanilla`; `init` is an opt-in modelling variant), the reversal-noise guard applies only to PRIMED states (after the first plastic commit since the flip), and under `-implex` an **opt-in** zero-increment companion return at the flip absorbs the drift (`-implexFlipAbsorb on|off`, default `off`; history zero, counted in `implexGuards[5]` only when `on`)|built in #807 (887fea475); Esmeralda 887fea475 with the guard confined to primed states: default `vanilla` returns the implicit twin to the pre-P2 number to the digit (6.511 / 11.539 / 16.117 / 20.528) — no modelling change on the default path; the flip's sign test is deterministic on a real deck (28 629/34 560 points set `α_in := α` at step 1, identical every run), so vanilla is a real loading-direction decision, not noise; `init` remains available and diverges as documented under RC14; Esmeralda dense P2-6 legs give ~9.9 kN at step 2 with no refusal and overlay the twin; hold census passes at every hold|||| **TIMs RC14 — resolved on evidence: default `vanilla`; `init` remains available as a declared modelling option.** Esmeralda 887fea475 (real deck, guard confined to primed states) shows the flip's sign test is DETERMINISTIC, not round-off noise: 28 629/34 560 points set `α_in := α` at step 1, identical every run, because it is deciding a genuine continuing-loading direction — and the implicit twin's first-step stiffness then matches the pre-P2 number to the digit (6.511 / 11.539 / 16.117 / 20.528). There is no defect at the flip to fix by default, so the fork does not ship a modelling change as default: `-flipAlphaIn vanilla` is default and reproduces real `ManzariDafalias` exactly; `init` (deterministic `α_in := α` at every point, unconditionally) stays available as a declared, opt-in modelling choice, and every P2-7 curve names the flag in its header. **P2-7c:** measured on the R3 footing (200 elements): `ops.updateMaterialStage('-material', 1, '-stage', 1)` reached ONE element's 8 Gauss points (`SRC/domain/component/MaterialStageParameter.cpp:76` — it breaks the domain-element loop as soon as `theEle->setParameter(...)` returns other than `-1`, i.e. after the first element that accepts the parameter); vanilla never noticed because `mElastFlag` is a class-wide static, so one dispatch flips every instance, but P2-7's per-instance work ran at 8/1600 points and `init` vs `vanilla` were indistinguishable on 887fea475. Fix: each instance detects the flip itself at its first plastic trial (`mElastFlag == 1` and not yet seen), does the init and, under `-implex` with `-implexFlipAbsorb on`, the zero-increment companion return, then marks itself — dispatch-independent, so database-restored and MPI instances are covered too; the flip-handled marker is itself now serialized (`sendSelf`/`recvSelf`/both `getCopy` forms), so a database-restored or MPI instance does not re-run the flip work on a redundant `updateMaterialStage` re-assert. **`-implexFlipAbsorb` made opt-in, default `off` (P2-7c):** ON unconditionally changed the committed state at the flip, which broke ADR-92 gate 5 (a zero-free-DOF deck must commit the same state with `-implex` ON and OFF) — the absorption is a modelling choice, not a defect fix, so it does not ship as default. With the default `off` and `-flipAlphaIn vanilla`, the flip is byte-identical to pre-P2-7 behaviour except that the reversal-noise guard is now confined to primed states (the actual fix); the un-primed first step's committed error (0.24 on the R3 probe) remains and is documented as the price of not absorbing. `-implexFlipAbsorb on` recovers it to ~0.05 at the cost of an ON/OFF difference at the flip. Status: shipped in #807 (P2-7c). |
| P2-6 | **Trial-time `f = 0` fallback before a control refusal** | with P2-1..5 the loop is gone but the dense Esmeralda leg 146569 crawled at 0.8 µm/step to its budget near s/B 0.019 because the trial error switches from < 1e-7 to > 0.1 between 0.8 and 1.6 µm at states throughout the shear zone: the P2-2 guard acts only on a COMMITTED predecessor, so a trial that first reaches a softening/reversing point is extrapolated with full `f` and refused | when the control's error exceeds `tol` (floor not reached) recompute `σ~` with `f = 0` for that Gauss point, re-measure against the same companion result, deliver it if under `tol` (`implexDetail[5] = 0`, counted in `implexGuards[4]`), else the existing refusal/floor logic; flag `-implexTrialGuard on\|off` (default on) | the 146569 leg must return toward the 1e-4 step past 0.018 with committed errors under tol and the fallback counts reported | built in #807 (708152eac); Esmeralda dense rerun with the trial guard pending |

Acceptance on 87b9cf846, dense q10 (Esmeralda 146569/146570): PASSES on correctness — no
committed row above tol, state consistent, strain O(ds), loop gone (13 floor fallbacks at
0.01773–0.01774 then none; refuse arm honest wall at 0.01774, within 0.1 % of the pre-P2
engine); FAILS on performance — crawl at 0.8 µm/step to the budget near 0.019; P2-6 is the
answer to the second half.

Not in P2: ADR 93's own question (the free-surface ring at the pressure floor, no plateau) —
untouched by all four; the campaign returns to it once P2-1/2 land.

**Default decision:** `-implexFloor` defaults to `implicit`, not `refuse` or the pre-P2 `accept`,
because the implicit return at the floor passes the state the control is refusing at — it is
admissible by construction — so `refuse` would stop IMPL-EX exactly where the implicit material
itself is still walking forward, and `implicit` closes the self-sustaining gap-closing loop P2-1
measured without paying for an extra return map (the companion is already computed for the error
comparison). `accept` remains available only to reproduce pre-P2 behaviour or isolate the loop.
