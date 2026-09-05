---
title: "Scoping review — IMPL-EX for LadrunoSANISAND (TIMs request of 2026-09-05)"
project: Ladruno
type: scoping-review
status: "REVIEW COMPLETE — request is SOUND in direction, needs 6 corrections and a renumber before it becomes an ADR; sequencing depends on #792"
priority: high
owner: nmora
requested_by: "TIMs Workbench, act `work/ape/response-curve-matrix`"
related:
  - "[[86_ladruno_sanisand_adr]]"
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[_adr90_tau0_qu_band]]"
  - "[[31_ladruno_concrete3d_adr]]"
tags: [scoping, adr, material, sanisand, implex, tims]
updated: 2026-09-05
---

# Scoping review — IMPL-EX for `LadrunoSANISAND`

> [!success] Owner decisions taken 2026-09-05
> **D7 = 92** (shell modifiers keeps 91) and **D0 = ADR + P0 oracle now, C++ after #792 T8**.
> The ADR that carries all of D0-D7 is `[[92_ladruno_sanisand_implex_adr]]`; the request as
> received is `[[_adr92_tims_request_2026-09-05]]`. This memo is kept as the evidence behind
> the corrections, not as the live plan.

The request as received is `Ladruno_implementation/91_ladruno_sanisand_implex_adr.md`,
**untracked on `ladruno` at `3f003d110`** in the main checkout. It is already written in
ADR shape by the consumer. This memo is the fork-side review: what survives, what must
change, and what the work package looks like.

**Verdict.** The direction is right and the driver is the best-evidenced one the fork has
(`_adr90_tau0_qu_band` §6.3). The document must not be adopted as written: **six of its
statements about this code base are wrong or unsafe**, one of them load-bearing on the
design, and it does not know about the work package that landed the same afternoon.

---

## 1. Sequencing — the request does not know about #792

**WP-86b / PR #792 (draft, opened 2026-09-05 18:17, branch `wp/86b-sanisand-integrator`)
is attacking the same seizure from the cheap end** and is four commits in:

| #792 item | state | what it does to this request |
|---|---|---|
| T1 `-maxSubsteps N` on `LadrunoSANISAND` (vanilla seam `mMaxSubstepsInME`, default 0 = uncapped) | committed `ee841bead` | Makes `ModifiedEuler` **able to fail** instead of force-accepting at `dT_min`; the D16 subdivision controller finally gets to act. This is the first thing to try, and it is nearly free. |
| `LADRUNO_MATERIAL_REFUSED = -33086` + `LadrunoBrick::update()` propagation, throttled reporting | committed `49937326e` | **The refusal convention `-implexControl` needs already exists.** ADR-92 must reuse it, not invent one. It also carries the ADR-33/34 lesson the hard way: a blanket `< 0` propagation is wrong and broke two green tests. |
| T2 `TanType` default 0 → 2 (consistent) | committed | Every fork SANISAND deck to date ran Newton as **modified** Newton. The GATE U seizure numbers were measured on the elastic tangent. |
| T8 GATE U re-run with cap + new default | not run | **This is ADR-92's baseline.** Until T8 reports, nobody knows how much of the wall survives the cheap fix. |

**Recommendation D0: do not open C++ on IMPL-EX until #792's T8 has reported.** If a capped,
consistent-tangent leg reaches a plateau, IMPL-EX is a performance and robustness option
rather than an unblocking one, and its phasing and its risk budget both change. The
scoping, the ADR and the oracle work below are not blocked by T8 and should proceed now.

---

## 2. Number collision — take 92, not 91

`91_ladruno_shell_stiffness_modifiers_adr.md` **also exists, also untracked**, in the
`wp/91-shell-modifiers` worktree, with `LadrunoShellModifierSection.{h,cpp}` already
written. Neither 91 is on `ladruno`, so this is cheap to fix now and expensive later —
see what ADR-78 cost by being three ADRs. The shell lane is further along and has the
branch name. **IMPL-EX takes 92.** (88 is cited throughout `SRC/` by #778; 89 is proposed
for Track T; 90 is the overstress wrapper.)

Branch for this work: **`wp/92-sanisand-implex`** — `wp/91-*` is taken.

---

## 3. The six corrections, each verified at source

### C1 — the integration-scheme map in §2 is wrong in every entry

The request says: *"0 forward Euler, 1 backward Euler implicit (`BackwardEuler_CPPM`),
2 backward Euler with stability considerations, 3 constrained-strain-increment explicit …
the modified-Euler substepping behind scheme 0/3."*

What `ManzariDafalias.cpp:40-49` and the dispatch at `:984-993`, `:1031-1057` actually say:

| `mScheme` | integrator | kind |
|---|---|---|
| 0 / 4 / 6 | `MaxEnergyInc` (MFE / FE / RK inner) | explicit, substepped |
| **1** | **`ModifiedEuler`** — error-controlled substepping, `dT_min = 1e-6` (`:1320`) | **explicit** — and the deck default (`:93`, `oData[0] = 1`) |
| **2** | **`BackwardEuler_CPPM`** | **implicit — the only one** |
| 3 | `RungeKutta4` | explicit, no error control |
| 5 | `ForwardEuler` | explicit, no error control |
| 7 / 8 / 9 | `MaxStrainInc` | explicit, substepped |
| 45 | `RungeKutta45` (Abell) | explicit, error control |

Consequences: the request's *"IMPL-EX is meaningful only on schemes 1 and 2"* names one
explicit and one implicit scheme; its *"refuse `-implex` on `mScheme` 0 and 3"* leaves
1, 5, 7, 8, 9 and 45 open. **The scheme the TIMs deck actually runs — and the one that
seizes — is 1, which is explicit substepping, not an implicit return.** Everything §2
says about "the implicit return at commit" is therefore a change of integrator as well
as a change of scheme, and that must be stated and defaulted deliberately.

### C2 — `mDGamma` is not the step-total plastic multiplier (load-bearing)

The request extrapolates `dGamma`, "position [25] of the recorded `gp_state`". Under
`BackwardEuler_CPPM` that is the step total (`:1171` `NextDGamma = 0;`, then solved).
Under every substepped explicit scheme it is **the last substep's multiplier**:
`ForwardEuler` assigns `NextDGamma = (…)/temp4` (`:1342`) fresh on each substep and
`ModifiedEuler` never accumulates it. Extrapolating it across a step whose substep count
changed by orders of magnitude — exactly what happens at the corner — extrapolates noise.

**D1 (recommended): extrapolate the plastic strain, not the multiplier.**
`eps_p = mEpsilon − mEpsilonE` is committed and available; store one extra committed
vector `d_eps_p_commit` and use

    eps_p~(n+1) = eps_p(n) + f · d_eps_p(n),      sigma~ = Ce(p_n) : (eps(n+1) − eps_p~(n+1))

This is integrator-agnostic — it works with the companion being scheme 1, 2 or 45 — needs
no access to the flow direction, no new virtual hook into vanilla, and carries the fabric
and reversal guards for free because nothing on the extrapolated path touches `mAlpha_in`
or `mFabric`. It is the same operator as the request's, with a better-defined history
variable. The alternative (extrapolate `dGamma` and freeze `n`, `R`, `D` at state n) is
the textbook form and stays on the table, but it **requires** the companion to be scheme 2.

### C3 — the time source works for this deck by luck, and must be declared

`ASDConcrete3DMaterial` reads `ops_Dt` (`:1624-1631`) with a `dtime_is_user_defined`
escape hatch via `setParameter`. `ops_Dt` is `Domain`'s `dT = currentTime − committedTime`
(`Domain.cpp:2054`, published at `:2125` and `:2392`), i.e. **the pseudo-time increment**.

- The TIMs harness drives settlement by `LoadControl` on a prescribed-settlement SP
  pattern, so pseudo-time is proportional to settlement and `dtime_n/dtime_n_commit` is
  exactly the right extrapolation ratio — **including under the ladder's subdivision**,
  which halves the load factor and therefore halves `ops_Dt`. Good news, but it is a
  property of that deck, not of the material.
- `LoadControl 0.0` holds and the `updateMaterialStage` re-equilibration inject `dT = 0`
  steps. `time_factor` then goes to 0 (an elastic step), and the next step divides by a
  committed 0 and silently falls back to 1.0. Both behaviours need a guard and a test.
- True `DisplacementControl` sets `dT = Δλ`, which is not proportional to the displacement
  increment and **changes sign past a limit point**. Refuse, or offer a strain-increment
  time source.

**D2:** `-implexDt {pseudo|strain|user}`, default `pseudo` (the ASD behaviour) with an
explicit refusal on integrators whose load factor is solved for rather than prescribed.

### C4 — the step-request mechanism exists, but not the one the request describes

ASD signals over-tolerance by returning `EC_IMPLEX_Error_Control = -10` from
`setTrialStrain` (`ASDConcrete3DMaterial.cpp:59-61`, `:1679-1684`) and lets the element
and the analysis interpret it. The fork has just decided that a bare negative code is
the wrong contract (#792 `49937326e`, and `LEDGER_quirks` on ADR-33/34): it must be the
sentinel `LADRUNO_MATERIAL_REFUSED (-33086)` in `SRC/material/LadrunoMaterialStatus.h`,
propagated only by exact value, with a process-budgeted report. **ADR-92 inherits that
convention and must not re-open it.** This also makes `-implexControl` and `-maxSubsteps`
share one refusal path, which is the right outcome.

### C5 — "ADR 90 … is not touched by it" is false in the direction that matters

IMPL-EX's extrapolation is a first-order-in-Δt perturbation of the constitutive response
with the same structure as an artificial viscosity; that is why it robustifies softening
in the first place (Oliver, Huespe & Cante 2008). **An IMPL-EX leg is a regularized leg
with the step size as the regularization parameter.** Everything the campaign reads off
such a leg — the localization width that GATE U measured as non-convergent
(`w(h/2)/w(h) = 0.675–0.824`), the matched-settlement band, any post-peak branch — is
regularized by Δt, and by a parameter with no length in it.

This does not block the work. It changes two things: the ADR must carry the disclosure in
the same words the act will have to print beside its verdicts, and **ADR 90's D2 gains a
second candidate regularizer that is already half-built.** The request's §6 "not in scope"
line disclaiming ADR 90 should become a cross-reference, not a firewall.

### C6 — the tangent claim is right, with one condition and one unclaimed prize

*Condition.* `dsigma~/deps = Ce` holds **only if `G` and `K` are frozen at the committed
`p_n`**. SANISAND's elasticity is pressure-dependent (`G ∝ sqrt(p')`), so evaluating the
moduli at the trial pressure — which the base does on the current stress — reintroduces a
strain dependence and the tangent is no longer the operator that was handed out. Freeze
them; test it as an identity (numerical vs returned tangent, to machine precision).

*Prize the request does not claim.* Under `-implex` the non-associated consistent tangent
disappears, so **the global matrix becomes symmetric**. That unlocks `system Pardiso
-matrixType sym` (ADR-75 P1d: 1.94–1.96× vs UmfPack, −42 % peak memory, exact). On a
21 058-DOF coarse cell and a 175 290-DOF fine cell that is a second, independent speedup
on top of the iteration-count collapse, and it should be measured in the same gate.

---

## 4. What the cost argument actually is

GATE U §6.3, restated as a budget: the worst step cost **2056 s** across up to **125
state-determination passes**, each pass running a substepped return that collapses toward
`dT_min = 1e-6`. IMPL-EX moves the expensive companion return to **one call per committed
step** and makes the global step linear (Newton converges in one iteration on a frozen
operator). The arithmetic is a ~100× cut in constitutive work per step **before** any
change to the cost of a single return, and it is the same lever #792's cap pulls from the
other end. Both should be measured against the same GATE U legs, which is why D0 puts T8
first.

---

## 5. Decisions the ADR must record

| | decision | recommendation |
|---|---|---|
| **D0** | Sequence against #792 | Scope and ADR now; **no C++ until T8 reports**. |
| **D1** | Extrapolate `dGamma` or plastic strain | **Plastic strain** (C2) — integrator-agnostic, no vanilla hook. |
| **D2** | Time source | `-implexDt {pseudo\|strain\|user}`, default `pseudo`; refuse solved-load-factor integrators. |
| **D3** | Companion integrator under `-implex` | Default **scheme 2** (`BackwardEuler_CPPM`) — the only implicit return; allow 1 with `-maxSubsteps` set, refuse 3/5/7/8/9. |
| **D4** | Class tags | **None new.** Flags on the existing 33019/33020/33021 (`LadrunoSANISAND`, 3D, PlaneStrain); the wire format grows, so `sendSelf`/`recvSelf`/`getCopy` gates per ADR-86. |
| **D5** | Vanilla footprint | Target **zero** under D1. If the `dGamma` route is chosen instead, one flag seam in `ManzariDafalias.h` on the `mHonorTolRInME` / `mMaxSubstepsInME` pattern, plus a `LEDGER_vanilla_files` row. |
| **D6** | Relationship to ADR 90 | Cross-reference, not firewall (C5). Disclosure text is a deliverable of P0, not P2. |
| **D7** | ADR number / branch | **92**, `wp/92-sanisand-implex` (§2). |

---

## 6. Proposed work package

**P0 — oracle and design, no C++.** Single-Gauss-point IMPL-EX of the D-L cell's SANISAND
constants in numpy (`python3.12`), driven by the act's Level-0 drained triaxial path:
first-order convergence of the IMPL-EX/implicit difference under increment halving, the
size of the error at the act's `1e-4 m` increment, and the behaviour at `p' → p_min`.
Decides D1 and D3 on measurement rather than argument. This is the ADR-90 WP-A pattern and
it is the cheapest place to be wrong.

**P1 — `-implex` on `LadrunoSANISAND3D`.** Extrapolated stress, frozen elastic tangent,
companion return at `commitState`, `implexError` / `avgImplexError` responses
(`ASDConcrete3DMaterial.cpp:2073-2077` is the template), `getCopy`/`sendSelf`/`recvSelf`.
Byte-identity with `-implex` unset is the gate.

**P2 — `-implexControl` and plane strain.** Error-driven refusal through
`LADRUNO_MATERIAL_REFUSED`; `LadrunoSANISANDPlaneStrain` for the 2D act.

**P3 — the campaign gates on Esmeralda.** The corner patch, the coarse bare `D-L` leg
against job 146299's wall at `s/B = 0.0206`, the `U-L` coupled row inside `LadrunoUP`, and
the symmetric-solver measurement of C6.

The request's own §5 acceptance list is good and mostly survives; its tests 1, 3, 6, 7 map
onto P0/P1, test 2 onto P2, tests 4 and 5 onto P3.

---

## 7. What the consumer is still owed a request for

- **OQ2 from TIMs** — the tolerance band. GATE U could not adjudicate without it and
  neither can this: "the IMPL-EX curve is within X % of the implicit one" needs an X that
  the campaign, not the fork, sets.
- **Which leg is the acceptance leg.** The request names job 146299 (`s/B = 0.0206`) as
  the wall to pass. After #792's T8 that wall may have moved; the ADR should name the
  criterion (a WP1 plateau) rather than the number.
