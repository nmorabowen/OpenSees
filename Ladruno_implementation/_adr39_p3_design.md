# ADR-39 P3 — implementation plan: IMPL-EX Coulomb friction (the v1 SHIP)

> Adds tangential Coulomb friction to the shipped frictionless NTS penalty contact
> (P2b/P2.5). Grounded in the friction oracle `contact_prototypes/proto_p3_friction.py`
> (6/6: stick, slip-caps-μN, dissipation, non-symmetric tangent vs FD, IMPL-EX→implicit,
> sliding-block incline) + the ADR §friction grounding. Parent ADR `39_..._adr.md`;
> driver `_adr39_loop_state.md`. P3 = a CRITICAL JUNCTION → adversarial design gate before C++.

> **⇒ DESIGN GATE RUN (2 source-grounded reviewers: state/lifecycle + mechanics/signs) →
> verdict SALVAGEABLE; folded below.** Core architecture (engine-owned keyed state, commit-at-
> Domain-hook, clean return map) CONFIRMED sound. Caught BEFORE any C++: a silent energy-
> injection SIGN bug, a wrong slip ORIGIN, and IMPL-EX being useless-but-harmful in the explicit
> ship. The 7 folds are marked `[GATE]` inline; the verdict + fix list is at the end.

## Scope (v1 = explicit ship; the de-risk split the design gate mandated)

- **v1 P3 SHIPS the EXPLICIT leg: friction FORCE only.** Under `CentralDifferenceLadruno`
  the FE tangent routes to mass-only (`formEleTangent` = `addMtoTang`), so the friction
  tangent is NEVER assembled in explicit. ⇒ v1 needs ONLY the tangential traction in the
  residual — exactly the P2a/P2b pattern (residual exact, tangent main-term-only). This is
  the whole reason friction was always an explicit-first ship (ADR Q-FRIC / B1).
- **Implicit frictional Newton with the consistent NON-SYMMETRIC tangent = P3.5 (DEFERRED).**
  The oracle has the tangent (`∂tT/∂gT=(μN·kt/‖tT*‖)(I−n̂⊗n̂)`, `∂tT/∂gN=−μ·kn·n̂` — no
  symmetric partner ⇒ rigorously non-symmetric, needs `FullGeneral`/`UmfPack`). Not a v1 gate.
- **v1 uses the DIRECT implicit return map per step, NOT IMPL-EX `[GATE MAJOR-2]`.** IMPL-EX
  exists to give a constant SECANT tangent so implicit Newton converges cheap — but the explicit
  ship NEVER assembles the tangent (mass-only), so IMPL-EX buys nothing here while COSTING the
  documented one-step onset overshoot, which in explicit (no Newton to re-sync) is a real spurious
  force impulse advected forward by the leap-frog. The direct return map is exact at the cone every
  step, has no overshoot, and DROPS the `dlam/dlam_old/dlam_trial` state entirely. Keep the IMPL-EX
  functions in the kernel header for the P3.5 implicit leg only; do not wire them into the residual.
  (This also dissolves the firstStep-double-getResidual BLOCKER: with no extrapolation seed, the
  trial write `gpT_trial = gpT + Δ` is a pure idempotent function of committed state.)

## The constitutive update — CLEAN Coulomb radial return (NOT the ASDimplex damage form)

Author fresh in the kernel (the reference `ZeroLengthContactASDimplex::updateInternal`
uses an "apparent damage" formulation carrying the quarantined bugs — shadowed dE1/dE2,
ddata(39) OOB; we don't lift it). Per active pair, given `N = kn·⟨−gN⟩₊` (penetration),
tangential penalty `kt`, friction `μ`, committed plastic slip `gpT` (tangent-plane vector):

```
tT*  = kt·(gT − gpT)                 # trial elastic tangential traction
f*   = ‖tT*‖ − μN                    # cone
STICK (f*≤0): tT = tT*, gpT unchanged, dlam = 0
SLIP  (f*>0): n̂ = tT*/‖tT*‖ ; dlam = f*/kt ; tT = μN·n̂ ; gpT += dlam·n̂
```

`[GATE]` **v1 uses the direct map above (no IMPL-EX).** The near-zero slip guard MUST be a
PHYSICALLY-SCALED tolerance, not a denormal floor: treat `‖tT*‖ < tol·μN` (or `tol·kt·Lchar`)
as stick (`tT=0`, no normalize) so a numerically-tiny trial can't normalize a meaningless `n̂`
or hit `0/0`. `[GATE MINOR-2]` `N` for the cone `μN`: v1 uses the CURRENT-step penetration
`N=kn⟨−gN⟩₊` (matches the oracle, zero extra state; explicit has no within-step normal-shear
solve to destabilize). Using committed `N` for the cap is a smoothness refinement (decouples
the cap from normal-penalty ringing) — noted for P3.5, not v1. `[GATE MAJOR-4]`

The IMPL-EX variant (`dlam~ = max(0, dlam_c + (Δt/Δt_c)(dlam_c − dlam_old))`,
`tT = tT* − kt·dlam~·n̂*`) is retained in the kernel header for the **P3.5 implicit leg only**.

## Tangential kinematics — ENGAGEMENT-referenced slip `[GATE MAJOR-1]`

Tangential slip is measured from the **ENGAGEMENT** config (the step the pair first becomes
active), NOT the global undeformed reference. The gate caught the bug: a slave that engages
LATE (separated, drifts tangentially, then penetrates) has a global `(u_s−ū)` that already
contains the entire pre-contact tangential approach → a huge spurious stick traction at first
contact (physically the traction at engagement must be 0). Small-slide scope bounds the tangent-
plane ROTATION, not the pre-contact drift magnitude, so it does NOT excuse a wrong origin.

```
d    = (u_s − ū)                     # relative displacement, current − reference (global)
gT   = d − (d·n) n                   # tangential component on the CURRENT normal n (global 3-vec)
gTeff= gT − gT0                      # slip measured from the ENGAGEMENT config
tT*  = kt·(gTeff − gpT)              # trial elastic tangential traction
```

`gT0` = the tangential measure captured ONCE at first activation, persisted per pair (one
`double[3]`). `gpT` (committed) accumulates plastic slip relative to engagement. `ū` = the
projected master point built each call from REFERENCE coords + the current shape weights `N`
(the SAME `N` the B-operator uses, so normal + friction act at one consistent point) — NOT
stored. Convective caveat: as the slave slides, the projection point migrates within the facet,
so `(u_s−ū)` conflates physical slip with projection migration; bounded under the v1 small-slide
scope (slave stays within one facet, `‖gTeff‖ ≪ facet size` — assert + document). Finite sliding
(Laursen §6) is deferred. Plastic slip `gpT` lives in the tangent plane; large frame rotation →
epoch re-emit + re-projection, a documented v1 limit.

## Per-pair STATE on the engine (the first path-dependent contact state — gate crux)

Adapters are stateless VIEWS rebuilt every `handle()` (P1 design), so friction state lives
on the Domain-owned `LadrunoContactDomain`, keyed STABLY so it survives adapter rebuilds:

```
key   = (contactTag, slaveNodeTag, segIndex)          # segIndex = GLOBAL segment ordinal
state = { double gpT[3];        // committed plastic slip from engagement (global, tangent plane)
          double gpT_trial[3];  // trial (written by the adapter each getResidual)
          double gT0[3];        // ENGAGEMENT-config tangential origin [GATE MAJOR-1]
          bool   engaged; }     // has gT0 been captured (persisted) — set at first activation
```
`[GATE MAJOR-2]` v1 carries NO `dlam` fields (direct return map, not IMPL-EX). `active`/engaged-
this-step is PER-EVALUATION SCRATCH (recomputed every getResidual), never persisted `[GATE MINOR-1]`;
only `engaged`+`gT0`+`gpT` persist.

- `std::map<key,state>` on `LadrunoContactDomain`; lazily created per pair.
- **`segIndex` = the GLOBAL segment ordinal `seg`** (the `mTags(seg*nps..)` block index), NOT
  the bucket-sort candidate position `ci` — `ci` ordering isn't rebuild-stable; `seg` is (the
  surface is immutable post-construction). `[GATE BLOCKER-2]`
- **Adapter (getResidual):** look up the slot; if first activation capture `gT0=gT`, set `engaged`;
  read committed `gpT`; compute `gTeff`,`N`; run the DIRECT return map; WRITE `gpT_trial` as a PURE
  function of committed state (`gpT_trial = gpT + Δ`, never `+=` — idempotent across the firstStep
  double-eval) `[GATE BLOCKER-1/state]`; assemble the residual.
- **`mu ≤ 0` short-circuits at the TOP of the friction block, BEFORE any slot touch** `[GATE MAJOR-3]`:
  guarantees byte-identical-to-P2b frictionless result + zero map growth when friction is off +
  dodges the `n̂=tT*/‖tT*‖` 0/0 (the ‖M‖→0 class of bug from LadrunoJ2).
- **`Domain::commit()` → `LadrunoContactDomain::commit()`** (hook wired in P1b, today a counter):
  for every slot `gpT = gpT_trial`. **`revertToLastCommit()`:** `gpT_trial = gpT` (drop trial).
  NOTE explicit CDL has no convergence retry — a Python adaptive-step retry wrapper MUST call
  `revertToLastCommit()` (document; add an abort-then-retry test). `[GATE MAJOR-2-state]`
- **Dead-slot GC every `handle()`** `[GATE MAJOR-1-state]`: the engine survives `domainChanged`
  AND across `analyze()` calls, so a re-meshed second analysis leaks old slots (the ADR-30 `theEQs`
  leak class). The handler hands the engine the live key-set each rebuild; the engine erases slots
  not in it. (For v1 small-slide the candidate set is reference-coord-static, so within one analysis
  the pair set is stable; GC matters across analyses.)
- **Cross-contact shared-slave warning** `[GATE BLOCKER-2]`: a slave node in two contacts' slave
  sets gets two slots → double traction (a modeling error the engine can't auto-resolve); warn at
  handle() (cheap multiset check, mirrors the ndf!=3 warning).
- **Adapter→engine lifetime** `[GATE BLOCKER-2]`: `ops.wipe()`→`Domain::clearAll()` DELETES the
  engine; adapters (owned by AnalysisModel) must not outlive it holding a raw ptr. Adapter holds a
  `Domain*` and LAZILY re-fetches `theDomain->getLadrunoContactDomain()` in getResidual (null-check)
  rather than caching the engine ptr. Add a wipe-then-reanalyze test.
- **`segIndex`-flip slip-loss** `[GATE MINOR-3]`: if the broad phase flips a slave's nearest segment
  step-to-step, the key changes and slip history resets to stick on the new facet. Bounded under
  small-slide (slave shouldn't cross a facet seam); warn if an active pair's segIndex changes mid-contact.

## Residual assembly (SIGN — the gate's BLOCKER) `[GATE BLOCKER-1]`

`IncrementalIntegrator::formElementResidual` does `theSOE->addB(eleResidual, getID())` —
the element residual is ADDED to the RHS force with NO sign flip; CDL then does
`a = M⁻¹(P − Cv − F_int)` with the residual ON the RHS as an applied nodal force.

The return map gives `tT` along the trial slip direction `n̂* = tT*/‖tT*‖`, which points
ALONG the slave's tangential MOTION. So assembling `+tT` on the slave (mirroring the normal)
would push the slave ALONG its motion → `a = g(sinθ + μcosθ)` on an incline = silent ENERGY
INJECTION (the frictionless `mu=0` test passes regardless → false confidence). The earlier
"mirror the normal's `+`" reasoning was a non-sequitur: the normal's `n` opposes the motion by
construction, but `n̂*` FOLLOWS it.

FIX: friction must OPPOSE `n̂*`. Carry the negation in ONE auditable place — the kernel
returns the FORCE the contact APPLIES to the slave, already negated: `t_fric = −tT` (with
`‖tT‖=μN` at slip). Then the FE mirrors the normal block exactly:
`resid(slave) += tn·n + t_fric ; resid(master_i) += −N_i·(tn·n + t_fric)`. On the incline
this yields `M a = mg sinθ·t + (−μN)·t ⇒ a = g(sinθ − μcosθ)` ✓.

**The incline test MUST exercise the REAL OpenSees assembly** (CDL + LadrunoContactFE + `addB`,
reading the node acceleration from the solved SOE) — NOT re-implement the oracle's hand-written
`a=(drive−fric)/m`. Mirroring the oracle's subtraction in the test would let the sign bug ship
green. This is the single most important guardrail of the rung.

## FE adapter / parser

- `LadrunoContactFE` SEGMENT mode gains `kt,mu` + the engine-slot key; getResidual adds the
  friction traction. getTangent UNCHANGED for v1 (normal `kn·BᵀB` only; friction tangent =
  P3.5). The adapter needs a back-pointer to `LadrunoContactDomain` to reach its state slot.
- Parser: `contact tag m s kn kt mu …` ALREADY parses `kt,mu` (positional, P1b) — today the
  handler ignores them. P3 threads `ct.kt, ct.mu` into the segment-adapter ctor. `-kn auto`
  composes (auto kn, explicit kt/mu). No new parser surface.

## Tests — `tests/test_adr39_contact_p3_friction.py`

1. **stick** — a tangential load below `μN`: zero slip, slave held (tT = kt·gT < μN).
2. **slip caps at μN** — tangential load above `μN`: steady slide, friction force == μN.
3. **sliding block on an incline** — `a = g(sinθ−μcosθ)` (the analytic gate; explicit CDL).
4. **energy dissipation** — slip dissipates ≥0 (EnergyBalance / work check), stick conserves.
5. **frictionless regression** — `mu=0` ⇒ byte-identical to the shipped P2b result.
6. **commit/revert** — friction state commits on `Domain::commit`, reverts on a failed step.

## Files (P3)
- `SRC/domain/contact/LadrunoContactKernel.h` — add the friction return map (pure fns;
  `frictionReturnMap`, `frictionImplex`). Header-only, oracle-faithful.
- `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — per-pair friction state store + the
  real `commit()/revertToLastCommit()` (replacing the P1b counters) + a slot accessor.
- `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — kt/mu + engine back-pointer + friction
  traction in getResidual.
- `SRC/analysis/handler/LadrunoContactHandler.cpp` — pass kt/mu + the ContactDomain ptr + key.
- `tests/test_adr39_contact_p3_friction.py`.

## DESIGN GATE verdict + folds (2 reviewers: state/lifecycle + mechanics/signs → SALVAGEABLE)

Core architecture CONFIRMED sound: engine-owned per-pair state keyed `(contactTag,slaveTag,
segIndex)` surviving stateless-adapter rebuilds, committed at the `Domain::commit` choke point,
clean Coulomb radial return. The 9 findings, all FOLDED above (inline `[GATE ...]`):

1. **BLOCKER (sign)** — `addB` sums the residual; `+tT` on the slave ACCELERATES it (energy
   injection, `a=g(sinθ+μcosθ)`). FIX: friction OPPOSES `n̂*`; kernel returns the already-negated
   applied force; incline test through REAL CDL+addB (not the oracle's hand-subtraction).
2. **BLOCKER→dissolved (firstStep double-eval)** — CDL's starter calls getResidual at u_0 then
   u_1 (twice, no commit between). FIX: trial = pure fn of committed (`gpT_trial=gpT+Δ`); dropping
   IMPL-EX removes the only non-idempotent (multiplier-seed) hazard.
3. **MAJOR (engagement reference)** — global-reference `gT` corrupts late-engaging contact; store
   `gT0` at first activation.
4. **MAJOR (drop IMPL-EX from v1)** — useless (tangent discarded in explicit) + harmful (onset
   overshoot = real impulse); use the direct return map; removes `dlam` state.
5. **MAJOR (GC leak)** — prune dead slots each handle() (live key-set), ADR-30 `theEQs` leak class.
6. **MAJOR (mu≤0 byte-identity)** — short-circuit before any slot touch; also dodges 0/0.
7. **MAJOR (revert path)** — explicit has no auto-retry; document + test abort-then-retry revert.
8. **MINOR (current N)** — v1 uses current penetration N; committed-N cap = P3.5 smoothness.
9. **MINORs** — `segIndex`=global ordinal (not `ci`); cross-contact shared-slave warning;
   adapter lazy-refetches the engine via `Domain*` (wipe deletes the engine); near-zero guard
   physically scaled; `segIndex`-flip slip-loss + `kt≤kn` explicit-stability note.

**Status: design HARDENED, ready for C++.** NEXT = code (kernel direct return map + engagement
`gT0` + engine state/GC + mu≤0 short-circuit + negated friction force) → build → test (incline
through real assembly = the key gate, + stick/slip/energy/frictionless-regression/commit-revert/
wipe-reanalyze) → code gate → PR base ladruno.
