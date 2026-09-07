---
title: "ADR 93 — LadrunoSANISAND at zero confinement: what the model should do once p reaches the floor"
project: Ladruno
type: adr
status: "BRAINSTORM — problem stated, candidates listed, no decision; owner + TIMs to choose the P0 experiments"
priority: high
owner: nmora
related:
  - "[[86_ladruno_sanisand_adr]]"
  - "[[86_ladruno_sanisand_handoff]]"
  - "[[86_ladruno_sanisand_pr3_tripwire_memo]]"
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p1_bvp_gate_rerun]]"
  - "[[_adr92_cp1_surcharge_results]]"
  - "[[84_ladruno_mc_tension_cutoff_adr]]"
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
tags: [adr, sanisand, low-confinement, free-surface, integrator, brainstorm]
updated: 2026-09-06
---

# ADR 93 — LadrunoSANISAND at zero confinement

> [!abstract] **The question this ADR exists to answer**
> A cohesionless sand has no strength at zero confinement. Every strip-footing leg in the
> ADR-90/92 campaign is decided at the free-surface ring beside the footing edge, where the
> model reaches that state and the fork's integrator either **seizes** (thousands of substeps
> then a forced accept), **clamps** (pressure pasted back to `p_min`, deviator kept — no longer
> the model's answer), or, under IMPL-EX, **refuses** until the step is gone. IMPL-EX (ADR 92)
> removed the solver ladder and reaches the target settlement, but it cannot make the ring
> carry load, and the surcharge (ADR-86b §5b option B) bought separability, not confinement.
> Three implicit legs with the substep cap in the trial are the only clean runs, at s/B 0.0025.
> **What should the material do there?** This ADR lists the candidates with the experiment
> that decides each. It takes no decision.

## 1. What happens today at the ring (source-verified facts, nothing new)

| fact | where | consequence at p → 0 |
|---|---|---|
| Elastic moduli `G, K ∝ sqrt(p / P_atm)` | `ManzariDafalias::GetElasticModuli` | stiffness → 0; the element carries nothing and its tangent is singular-ish |
| Plastic modulus and dilatancy carry `p` via `psi`, `M^b`, `M^d`, `D` | `GetStateDependent` | the cone collapses to a point; the return has nothing to return to |
| `m_Presidual` enters ~30 mean-stress sites as `p = tr(σ)/3 + p_r` | header note, `LadrunoSANISAND.h:28-62` | vanilla's 1.01 kPa is an apparent cohesion `c = p_r tan φ ≈ 0.95 kPa` **and** bounds the `D_factor` dilatancy sigmoid from below (floor 0.4278); the fork's default `p_r = 0` (NTUASand02) drops that floor by 886× (**D5a, still open**) |
| `m_Pmin` clamp: deviator preserved, pressure pasted in | `ModifiedEuler`, `RungeKutta45` | a clamped Gauss point is a projection, not a constitutive answer; `p_min = 1e-3 P_atm = 0.101 kPa` in the fork |
| Substep error norm is **already mixed, with a hardcoded switch**: `err = ‖Δσ₂−Δσ₁‖` if `‖σ‖ < 0.5` else `‖Δσ₂−Δσ₁‖ / (2‖σ‖)`, against `TolE = 1e-4`; `dT_min = 1e-6`, `q = max(0.8·sqrt(TolE/err), 0.1)` | `ManzariDafalias.cpp:1746-1753`, `:1419`, `:1655` | the switch `0.5` is **unit-bearing** (0.5 kPa on a kPa deck = 5e-3 P_atm) and not declared; at the ring the absolute branch demands `1e-4` stress units per substep while the ring's strain increment per step is large (it has no stiffness to resist), so `dT` collapses → 1000–5000 substeps on the campaign deck, 64 000 unregularised; `-maxSubsteps` (ADR-86b) caps the count and force-accepts. P0 (ADR 92) showed the seizure does **not** reproduce at a prescribed-strain Gauss point with deck-sized increments: it is the *BVP path* (ring strain per step) that drives it |
| `BackwardEuler_CPPM` low-p branch returns `errFlag = 0` without solving | `ManzariDafalias.cpp:2264` (P0, ADR 92) | the "implicit" companion is the same explicit return in this regime |
| A material cannot refuse at commit | `Domain::commit()` discards `commitState()`'s return (ADR-92 review B3, quirks row) | any "rewind the step at commit" design needs a pre-commit gate or a vanilla change |
| IMPL-EX freezes `Ce(p_n)` and extrapolates `d_eps_p` | ADR 92 | the global step is linear regardless of the ring; the companion still seizes or is capped; `-implexControl` refuses the ring's error and spends the subdivision budget (registered arm: s/B 0.0275, 10 270 refusals) |
| Surcharge outside the footing | ADR-86b §5b B, CP1 | separable at collapse (`Q·N_q`), leaves the material untouched, but 10 kPa is 100× the floor and the ring beside the edge still sits at `p ≈ 0` |

**The seizure is a numerics symptom of a physics fact.** The model is *correct* to have no
stiffness and no strength at `p = 0`; the integrator is *wrong* to spend 5 000 substeps proving
it to relative tolerance on a stress that carries nothing. The two must be separated in any
answer: what the ring *is* (physics) and how much it may *cost* (numerics).

## 2. Constraints any answer must respect

1. **The soil is calibrated.** `G0 = 264.32` etc. were fitted with vanilla's `p_r = 1.01 kPa`
   and its `D_factor` floor. Any floor change is a modelling move on a calibrated material
   (ADR-86b §5b), so it must be declared, wired, echoed and *measured* against the calibration
   rows, not inherited.
2. **Zero vanilla footprint** on `ManzariDafalias.cpp` unless a ledger row is opened for it;
   `LadrunoSANISAND` wins by last write on protected data and by overriding virtuals.
3. **Parallel and database parity**: whatever is added goes on the fork wire (`Vector(97+)`),
   or an MP worker runs a different material than the serial one (the `p_r` wire gap already
   bit once).
4. **Falsifiable at three levels**, cheapest first: the P0 numpy oracle at a Gauss point
   (`adr92_p0_oracle/`, reproduces the binary to 1e-13), the G0 binary probe
   (`probe_binary_triaxial.py`), then the GATE U strip-footing deck (`sanisand_tau0_band.py`).
   No candidate goes to the deck before it has a Gauss-point number.
5. **Report what the ring does to the answer**, not only whether the leg runs: the ring's
   share of the bearing resultant, the clamp / cap / refusal census per step, and the overlay
   against an implicit control where one exists (the ADR-92 gate's discipline).

## 3. Candidates

Grouped by *where* the change lives. For each: mechanism, what it changes physically, cost,
the decisive experiment, and a first verdict to be argued with.

### I. Integrator-side (bound the cost, change no physics)

**I.1 Declare and scale the substep norm's absolute floor.** The norm is already mixed
(§1): absolute below `‖σ‖ = 0.5`, relative above, i.e. an effective `σ_ref = 0.5` stress
units that is hardcoded and unit-bearing. Make it a declared, `P_atm`-scaled parameter
(`-substepStressRef`, echoed, on the wire) and measure the trade-off. *Physics:* none — a
Gauss point below `σ_ref` is integrated to an absolute accuracy `TolE·σ_ref` that is
negligible for the BVP because that point carries nothing. *Cost:* substeps scale roughly
as `1/sqrt(σ_ref)` under the `q` rule, so a 10× cut needs ~100× `σ_ref` (≈ 50 kPa on this
deck), which is no longer "negligible". Expect a modest gain, not a cure. *Experiment:*
oracle — replay the corner Gauss point's *actual BVP strain path* (P0 showed a
prescribed-strain path does not seize) under `σ_ref ∈ {0.5, 5, 50}` stress units: substeps
per step, committed stress, and the ring's stress error against the tightest run. *Verdict:*
**run first because it is an hour and it calibrates every other candidate**: it tells us how
much of the seizure is the norm (fixable at zero physics cost) and how much is the ring's
strain per step (which only a stiffness floor, II.1, or an elastic skin, I.3, can reduce).

**I.2 Substep-count cap in the trial (exists: `-maxSubsteps`, ADR-86b).** Bounds cost,
force-accepts the last substep, reports `mSubstepCapHitInME`. *Verdict:* keep as the safety
net behind I.1, not as the answer — a capped point is an unconverged point, and the three
clean implicit legs at s/B 0.0025 are the existence proof of what it can and cannot buy.

**I.3 Elastic-below-p_switch.** Below a declared `p_switch`, skip the plastic return and
integrate elastically with floored moduli. *Physics:* the ring becomes a weak elastic skin —
a fiction, but a bounded and declared one. *Cost:* trivial. *Experiment:* oracle + GATE U leg,
report the ring's share of the resultant. *Verdict:* second-tier; it answers "how much may it
cost" by declaring the ring irrelevant, which must be *shown* (ring share), not assumed.

**I.4 Commit-time rewind.** Refuse a step whose companion failed and rewind. Needs a
pre-commit gate (the integrator asks each material `canCommit()` before `Domain::commit()`)
or a vanilla change to propagate the return. *Verdict:* infrastructure, not an answer to the
ring; park under ADR 92 P2.

### II. Material-side (declare what the ring is)

**II.1 Decoupled floors: stiffness floor ≠ strength floor.** Today one number (`p_r`) both
stiffens (`G(p + p_r)`) and strengthens (`M^b, M^d, D` through `psi(p + p_r)`). Split them:
`p_r,e` in the elastic moduli only (keeps the tangent regular, adds no strength, no cohesion),
`p_r,p = 0` in the plastic state functions. *Physics:* a sand with a small non-zero small-strain
stiffness at zero confinement — defensible (suction, interlock) and *calibratable separately*.
*Cost:* one new parameter on the wire; ~6 of the ~30 `p_r` sites move. *Experiment:* oracle
sensitivity of the G0 rows to `p_r,e ∈ {0, 0.1, 1} kPa` with `p_r,p = 0`; then the deck.
*Verdict:* **the most promising material-side option** because it separates the two things
the vanilla `p_r` conflates and leaves the calibration of strength untouched.

**II.2 Settle D5a — the `D_factor` sigmoid at `p_r = 0`.** The tripwire memo measured the
factor of 886 and listed four experiments. Whatever this ADR chooses, the sigmoid's floor is
part of the ring's physics and must be decided *with* it, not after. *Verdict:* fold the four
D5a experiments into this ADR's P0.

**II.3 Consistent tension/low-p cutoff (ADR-84 pattern).** Replace the `p_min` clamp
(pressure pasted in) with a composite return: Rankine-style cap at `p_min` returned
*consistently* (stress and internal variables together, tangent consistent), the way ADR 84
did for Mohr–Coulomb. *Physics:* same floor, but the clamped point is now a constitutive
answer with a tangent, not a projection. *Cost:* a `special_return` in the fork subclass;
medium. *Experiment:* oracle — clamp census and tangent symmetry/consistency on the corner path.
*Verdict:* worth doing *after* I.1 shows whether clamping is still frequent once the norm is
fixed; if clamps go to zero this is moot.

**II.4 Declared residual pressure (ADR-86b §5b option A).** Already refused as the *default*;
remains available per deck (`-Presidual`). *Verdict:* only as a control arm in experiments,
never as the campaign's answer, for the reasons in §5b (moves `psi`, unmasks/masks `D_factor`).

### III. BVP-side (give the ring something to stand on)

**III.1 Crust.** A thin surface layer of a different, cohesive material (MC with small `c'`,
ADR 84) or of the same sand at declared `p_r`. *Physics:* real soils have one; the campaign's
idealised half-space does not. *Cost:* deck only. *Experiment:* GATE U leg with a crust of
thickness `t ∈ {0.1, 0.25} B`, report `q(s)` and the crust's share of the resultant.
*Verdict:* the honest BVP answer if the *question* is a real footing; the wrong answer if the
question is the idealised bearing-capacity benchmark, because the crust's `c'` enters `N_c`
and cannot be subtracted the way a surcharge can.

**III.2 Surcharge (option B, in use).** Separable, insufficient at the ring at 10 kPa; a
larger `Q` moves the answer toward `Q·N_q` and away from `N_γ`, which is the term the campaign
wants. *Verdict:* keep as the control arm; do not raise it to buy confinement.

**III.3 Ravelling / element erosion.** Elements whose Gauss points sit at the floor lose their
stiffness (or are deactivated) — sand at a free surface ravels rather than carrying load.
*Physics:* plausible for the ring, but it is a *failure* model layered on a constitutive one.
*Cost:* element-level flag; medium. *Verdict:* defer; only if I.1 + II.1 leave a ring that
still dominates the step.

### IV. Time-integration side (change how the ring is *paid for*)

**IV.1 IMPL-EX hybrid: extrapolate only where `p < p_switch`.** Implicit everywhere except
the ring. *Verdict:* the companion still seizes at the ring, so this buys nothing I.1 does not;
and the ADR-92 gate already shows the whole-domain IMPL-EX is accurate to 2–5 % where it can be
checked. Park.

**IV.2 Viscoplastic regularisation (ADR 90, Duvaut–Lions).** Bounds the plastic strain *rate*;
the inviscid backbone return still runs at the ring. *Verdict:* orthogonal; does not address
this ADR's question.

## 4. Proposed P0 (Gauss-point only, one session, no C++)

All in `adr92_p0_oracle/sanisand_implex_oracle.py`, which already reproduces the binary:

1. **Extract the corner path.** From a CP1 leg's field dumps, the strain history of the
   worst ring Gauss point (`x = B/2`, top row). If the dumps lack it, run one short deck leg
   with a strain recorder on that element.
2. **Replay with I.1.** `σ_ref ∈ {0.5 (today), 5, 50}` stress units, i.e. `{5e-3, 5e-2, 0.5} P_atm`:
   substeps per step, committed `σ`, `e`, `α`; the difference in committed stress is the
   *cost of the floor in accuracy*, the substep ratio is the *saving*.
3. **Replay with II.1.** `p_r,e ∈ {0, 0.1, 1} kPa`, `p_r,p = 0`: same outputs, plus the
   tangent's smallest eigenvalue along the path (does the floor keep it regular?).
4. **D5a (II.2).** The tripwire memo's four experiments on the same path.
5. **Control arms.** Vanilla `p_r = 1.01` (option A) and `p_r = 0` unmodified.

Decision rule, written before the run: if I.1 at a *defensible* `σ_ref` (≤ 1e-2 P_atm) cuts
the ring's substeps by ≥ 10× with committed-stress change below the deck's push tolerance,
P1 = I.1 in C++ and II.1 is deferred. If it does not — the a-priori estimate says it will
not, because the driver is the ring's strain per step, not the norm — P1 = II.1 + II.2
together (a stiffness floor reduces the ring's strain per step at its source), with I.1's
declared `σ_ref` kept as the cost bound. III.1 is decided by the *question*, not by P0: ask TIMs whether the
benchmark is a real footing or the idealised half-space.

## 5. Open questions for the owner and TIMs

- Is the campaign's question `N_γ` on an idealised half-space (then no crust, no cohesion, the
  ring must be *paid for* not *stiffened*: I + II.1) or a real footing (then III.1 is legitimate
  and the whole ring problem is a modelling choice)?
- Who owns D5a? It is a calibration question on Gorini's material; the fork can measure, the
  decision to change the sigmoid's floor is his.
- Is a declared elastic floor `p_r,e` acceptable to the material's author as a fork parameter
  (default 0, echoed) — or must it stay a deck-level flag only?

## Log

- 2026-09-06 — Opened after ADR-92 P1 shipped (#798–#800): IMPL-EX reaches the target but the
  ring cannot be made to carry load; surcharge bought separability not confinement. Candidates
  listed, P0 proposed, no decision.
- 2026-09-06 (later) — **External evidence, TIMs Esmeralda footing (session
  `nonlinear-response-curve-planning-2b60bc-95`, engine c162833ed, B = 1.5 m, 0.5 m H8 bbar
  elements, 21 058 DOF, DisplacementControl + strain clock, -Presidual 0 and the 1.01 vanilla
  twin):** with the cap raised to 20000 (their original 1000 was the misuse: zero cap hits
  after) the registered `-implexControl 0.1` arm still walls at s/B 0.0009 on **4082 control
  refusals on the tolerance alone**, throttled to ds 4e-7–1.6e-6 m; the vanilla twin throttles
  to 3e-6 m; **a 10 kPa surcharge changes nothing**. Their ring is DENSE (e = 0.60). On this
  ADR's gate (B = 2 m, 1 m elements, loose e = 0.6944) the worst ring element's plastic
  deviatoric strain is 2.5–3.1e-4 per mm of settlement on every arm; the controlled arm
  accepted 2.3–3.6e-5 per step at ds 0.09–0.15 mm. Scaled to their elements, their step
  size sits in the same regime, so **density, not step size, separates the two rings** —
  consistent with CP1's dense deep legs dying at half the loose depth and GATE U's clamp
  firing on dense only. This puts II.2 (the `D_factor` floor at `p_r = 0`, D5a) ahead of I.1
  as the first thing P0 must separate.
