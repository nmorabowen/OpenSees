---
title: "ADR 92 — IMPL-EX integration for LadrunoSANISAND: a stress update that cannot seize, and a tangent that cannot lose ellipticity"
project: Ladruno
type: ADR
status: "PROPOSED — P0 COMPLETE (D1 = A confirmed, D3 REVERSED, floor-clamp item added); C++ GATED on #792 T8 (CP1)"
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
| **P1** | `-implex` on `LadrunoSANISAND3D`: extrapolated stress **(incremental, clamped at `p_min`)**, frozen `Ce(p_n)`, companion at `commitState` **(scheme 1 + `-maxSubsteps`)**, `implexError` / `avgImplexError`, `getCopy` / `sendSelf` / `recvSelf`, stage-0 inertness. **The P0 oracle is the reference: `implex_A` on the binary's own paths to 1e-8.** | `-implex` unset is **byte-identical** on every existing SANISAND deck; tangent identity to machine precision; Zone-A green. |
| **P2** | `-implexControl` through `LADRUNO_MATERIAL_REFUSED`; `LadrunoSANISANDPlaneStrain` (the 2D act needs it on the strip); the cyclic/reversal test. | Error-driven refusal provably triggers subdivision and the committed state is intact across it. |
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
