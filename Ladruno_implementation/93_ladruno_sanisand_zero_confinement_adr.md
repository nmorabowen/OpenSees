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
| Elastic moduli `G, K ∝ sqrt(p / P_atm)`, with `p = tr(σ)/3` floored at **`m_Pmin`** and **`m_Presidual` absent** | `ManzariDafalias::GetElasticModuli`, all three overloads (`:4834-4835`, `:4876-4877`, `:4896-4897`) — **source-verified 2026-09-06, `_adr93_p0_ring_replay`** | stiffness → 0; the element carries nothing and its tangent is singular-ish. **There is NO stiffness floor today**: `p_r` does not enter here, so `G, K → 0` with `p` whatever `p_r` is, and the only thing standing between the moduli and zero is `p_min = 1e-4 P_atm` |
| Plastic modulus and dilatancy carry `p` via `psi`, `M^b`, `M^d`, `D` | `GetStateDependent` | the cone collapses to a point; the return has nothing to return to |
| `m_Presidual` enters ~30 mean-stress sites as `p = tr(σ)/3 + p_r` — **all of them plastic-side**: the yield function `GetF`, `psi`, `M^b`, `M^d`, `D`, the `D_factor` sigmoid's argument, and the `ModifiedEuler` / `BackwardEuler` low-`p` guards. It does **NOT** enter `GetElasticModuli` (row 1), so it floors STRENGTH and never STIFFNESS | header note, `LadrunoSANISAND.h:28-62`; 38 `m_Presidual` sites in `ManzariDafalias.cpp`, none in `GetElasticModuli` | vanilla's 1.01 kPa is an apparent cohesion `c = p_r tan φ ≈ 0.95 kPa` **and** bounds the `D_factor` dilatancy sigmoid from below (floor 0.4278); the fork's default `p_r = 0` (NTUASand02) drops that floor by 886× (**D5a, still open**) |
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

**I.5 Floor-aware `-implexControl`: refuse on error only where the stress is above the
floor; count the rest.** The Esmeralda measurement (Log, 2026-09-07) shows the ring's
extrapolation error is **O(1) in ds** at the wall: the committed error *rises* 30× while the
step shrinks 16×, so the control refuses every trial down to 2e-7 m and the leg walls. A
fixed discrepancy however small the increment is the signature of a **non-smooth event** —
the `p_min` projection in the companion (deviator kept, pressure pasted in) or a dilatancy
sign flip — that a linear extrapolation from the last committed step cannot see. Refusing on
it buys nothing: no step size cures it, and the uncontrolled arm has already shown the global
answer is insensitive to that point (overlay 2–5 % where checkable) because it carries
nothing. *Mechanism:* at the trial, if the companion clamped or `p_impl ≤ k·p_min`, do not
refuse; increment a fourth census (`implexRefusals` grows a "floor-exempt" slot) and keep
measuring the error. Elsewhere the control is unchanged. *Physics:* none — it is ADR-86b's
option C (accept-and-count) **localised to the floor points only**, where the model has no
answer to protect. *Cost:* trivial; wire slot. *Experiment (pre-registered):* the registered
arm with I.5 must reach the control-OFF arm's depth with the same overlay (≤ 5 % mean vs the
implicit control over its reach) and report the exempt count per step; if the overlay
degrades, the exempted points were carrying load and I.5 is refuted. *Verdict (revised 2026-09-07):* **PARKED — its motivating evidence was withdrawn.** The O(1)
event was DisplacementControl's trial states, not the ring's material state (Log); under the
fork's push idiom the control's error is O(ds) and the registered arm walks. Keep the census
idea only if a floor-bound ring ever shows up under LoadControl(-ds).

### II. Material-side (declare what the ring is)

**II.1 A stiffness floor, which today does not exist.** *(Reworded 2026-09-06 against the source;
the earlier text claimed `p_r` "both stiffens (`G(p + p_r)`) and strengthens" and that II.1 was a
**split** of that coupling. It is not: `p_r` never reaches `GetElasticModuli` (§1 row 1), so today
`G, K → 0` with `p` whatever `p_r` is and the only floor is `p_min`.)* So II.1 is not a split — it
is a **NEW, elastic-only parameter `p_r,e`** entering the moduli as `G, K ∝ sqrt((p + p_r,e)/P_atm)`
and nowhere else, with the plastic floor `p_r,p` left free to be 0. *Physics:* a sand with a small
non-zero small-strain stiffness at zero confinement — defensible (suction, interlock), adds no
strength and no cohesion, and *calibratable separately* precisely because it touches nothing the
strength calibration used. *Cost:* one new parameter on the wire; **3 sites** (the three
`GetElasticModuli` overloads), not "~6 of the ~30". *Experiment:* oracle sensitivity of the G0 rows
to `p_r,e ∈ {0, 0.1, 1} kPa` with `p_r,p = 0`, **on a Gauss point that reaches the floor** — the
P0 replay ran the arms but on a confined point, where `λ_min(Ce) ≥ 3.97e4 kPa` and the question
cannot be posed (`_adr93_p0_ring_replay` §3). *Verdict:* **still the most promising material-side
option**, because it is the only candidate that puts a floor under the stiffness at all.

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
6. **The O(1) event (I.5).** On the Esmeralda worst-Gauss-point path (their dump, `implexDetail`
   slots 3–5 + `implexRefusals` per step), replay in the oracle and separate the absolute
   discrepancy `‖σ~ − σ_impl‖` from the denominator `‖σ_impl‖ + P_atm‖ε‖`: if the absolute
   part is O(ds) and only the ratio is O(1), the *measure* is the problem and I.1's `σ_ref`
   logic applies to the control too; if the absolute part is O(1), find the event (clamp
   fired? `D` sign?) and I.5 is the answer. Lane E's oracle already showed the dilatancy
   sign flip is quiet at `p0 = 5` and `100 kPa` — it has not been run at `0.1 kPa`.

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
  firing on dense only. **Corrected 2026-09-07:** their loose twin (e 0.6944, job 146445) walls at s/B 0.00196
  vs dense 0.0009–0.0011 — density moves the wall by 2×, it does not remove it; step growth
  walls at the same places, and with the surcharge EARLIER (0.00044) because growth reaches
  larger ds sooner. Their six legs ran the PSEUDO clock under DisplacementControl with zero D2
  firings (the 889 600 were on the cap-1000 legs).
- 2026-09-07 — **The measurement that reorders this ADR (Esmeralda, loose job 146445 and
  dense 146450):** `implex_err_max` at commit does **not** scale with ds at the wall. Last 40
  rows: ds 1.25e-5 → 4–9e-4; 2.5e-5 → 1.5–2.2e-3; 5e-5 → 6e-3 then 2e-2; then at 2.5e-5 it
  RISES to 2.6e-2 and 8.7e-2, and at 3.1e-6 it is still 4.5e-2 while Q turns over. Error ×30
  while ds ÷16. The extrapolation error at the ring is an **O(1) event**, not O(ds): the
  ring Gauss point crosses a non-smooth state (the `p_min` projection or the dilatancy sign)
  where a linear extrapolation is wrong by a fixed amount however small the increment. Added
  candidate **I.5 floor-aware control** (accept-and-count localised to floor points) as the
  cheapest lift and P0 item 6 to separate absolute discrepancy from denominator. Their
  worst-Gauss-point path dump (`implexDetail` 3–5, `implexRefusals`, per step) is P0's input.
- 2026-09-07 (03:00) — **Reversal, and the reading that stands (Esmeralda jobs 146451 dense,
  146453 loose, 146438 implicit twin; same engine/mesh/material, `p_r` 1.01, cap 20000,
  control 0.1/0.01, pseudo clock, surcharge 10):** under **`LoadControl(-ds)` on the
  prescribed-settlement sp** — the fork's own push idiom — the registered IMPL-EX arm walks
  straight past every DisplacementControl wall: s/B 0.0079 (dense) and 0.0082 (loose) at
  25 min, at the full 1e-4 m step, committed ring error steady at 3e-4–1e-3, clamp never
  fired at the ring element, zero D2 refusals, and the only control refusals are the early
  growth attempts (refused at iteration 1, accepted six steps later). **The error IS O(ds)
  here**, and the dense IMPL-EX leg overlays the implicit leg to line width over the whole
  overlap. So the O(1) event of the previous entry was **not the material's state at the
  ring**: it was the trial states DisplacementControl puts the ring through — a load-factor
  prediction from the frozen elastic tangent lands an O(1) trial strain on a near-zero-
  stiffness ring whatever ds is, and the control refuses at iteration 1 (their inference from
  the two refusal patterns and the two error scalings; not measured per iterate). Candidate
  I.5 is parked. For the ledger and the guide: **the fork's push idiom is a precondition for
  `-implexControl`; under DisplacementControl the control sees trial strains that do not
  scale with the step.** Both entries written. The ring question of this ADR — no plateau,
  what the material does at the floor — is unchanged; what changed is that IMPL-EX's control
  is no longer in the way of measuring it. Paths: `labs/ape/response-curve-matrix/level3/
  D-L-dl-vt-{dense,gorini}-q10-sp-{implex,implicit}/coarse/out/` on the TIMs worktree,
  `ring_point.csv` per step at Gauss point 5 of element 4047; ESMERALDA.md §48–49.
- 2026-09-06 — **P0 ran on those dumps and is INCONCLUSIVE for this ADR's question**
  (`_adr93_p0_ring_replay`): the instrumented point is chosen by GEOMETRY (top row, first element
  outboard of the footprint) and is **confined throughout** — `p` rises 6 → 51 kPa (dense) / 35
  (gorini), `clamp_fired` = 0 on all 362 steps, the `D_factor` sigmoid never fires. So I.1's
  substep cut is 1.00× *because the point never seizes*, II.1's floor is a ≤ 8 % perturbation of
  `G` on a regular tangent, and D5a is unmeasurable; the §4 decision rule **stays unapplied** and
  the tables are kept only as a confined-point control arm. The reproduction gate reached 3e-5
  (not 1e-6): `ring_point.csv` carries no `getState` slots, so `alpha`/`z`/`alpha_in` had to be
  fitted. **What P0 needs: the point with MINIMUM committed `p` at the wall — chosen by state, not
  geometry — from a leg that walls, with `alpha`(6), `z`(6), `alpha_in`(6) and `e` dumped every
  step.** One thing the run did settle: §1 row 1 above, and the II.1 rewording that follows from
  it. Open for ADR 92: the IMPL-EX legs replay ~100× worse than their implicit twins (2.4e-3 /
  6.8e-3 against 2.8e-5 / 2.5e-5) on fitted states, unexplained.
- 2026-09-07 (endpoints + the first NAMED O(1) event under the fork idiom; Esmeralda jobs
  146451/146455 dense, 146453/146456 loose, implicit twin 146438): the dense fork-idiom arm
  walked at 1e-4 m to **s/B 0.0176** (Q 2019 kN half-model, tangent still rising) and there
  the control refused at iteration 1–2 on every halving; with the floor at 1.95e-7 and budget
  200 it crept to 0.0185 on 900 accepted steps of ~1e-6 m. **Step 315 (s/B 0.017577,
  ds 3.9e-7 m) is a COMMITTED row with `implex_err_max = 1.111`** over every SANISAND point,
  then 1.7e-2, 1.5e-2 on rows 327–328, 8e-8 elsewhere: at a 0.4 µm increment the companion
  moved one Gauss point O(1) away from `σ_n` (the extrapolation term is negligible at
  f ≈ 4e-3). Not a step-size effect: under a frozen `Ce` the global step is linear, iterate 1
  *is* the converged trial, and as ds → 0 both `σ~` and the companion's input tend to `σ_n` —
  so an O(1) discrepancy means the companion moves the **committed state itself** on a
  near-zero increment. The instrumented ring point had 3e-5 on that row; the seat is elsewhere
  and unnamed. Three fork-side candidates, none a step-size question: (a) the committed state
  is off the surface as the companion sees it and the default-ON `Stress_Correction` fires on a
  zero increment (also fits the replay's 100× reproducibility gap and the un-primed stage-flip
  mechanism), (b) the `p_min` projection pasting a pressure into a point far below the floor,
  (c) a companion force-accept at `dT_min`. **Decisive census, no fork change:** the argmax
  point at commit on rows with error > tol (element, GP, `implexDetail` 0–5, stress, and the
  `alpha` / `fabric` / `alpha_in` / `state` responses), a zero-increment hold at the wall read
  over every point ((a) fires on a zero increment, (b)/(c) do not), and a third arm with
  `setParameter stressCorrection` off, labelled a model change. Queued on Esmeralda. If (a):
  the fix is ADR 92 P2 (commit the companion's corrected state consistently), not this ADR.
  The loose leg's stop at 0.0297 was a harness artefact (growth tries charged to the budget),
  rerun without growth in flight. Whether the implicit twin passes 0.0176 (at 0.0174,
  subdividing) decides whose wall it is.
- 2026-09-07 05:30 — **Whose wall: IMPL-EX's.** The implicit DisplacementControl twin of the dense
  q10 cell (146438, same everything, no `-implex`) walked through 0.0176 without noticing, at
  s/B 0.01823 and 2107 kN, tangent still rising, steady at 5e-5 m per step. The state where the
  control refuses every step and commits the 1.11 error at 0.4 µm is one the implicit return
  handles routinely. So the step-315 event is a property of the **committed state under
  IMPL-EX**, not of the material at that state — which moves it out of this ADR and into ADR 92
  P2, and promotes candidate (a) (off-surface committed state + default-ON stress correction).
  Loose IMPL-EX leg 146456 shows the same signature deeper (s/B 0.0396, 2497 kN). A fork-side
  zero-increment probe on the P1 test deck (hold after a plastic history; `implexError` on the
  hold must be ≈ 0 if the committed state is consistent) is running now; Esmeralda's census,
  zero-increment hold and stress-correction-off arm follow.
- 2026-09-07 (fork-side probes; scratchpad only) — **Hypothesis (a) REFUTED on a clean history.**
  Zero-free-DOF confine-first deck, `stdBrick`, ds 1e-3, 30 plastic steps, zero companion cap
  hits (max 17 572 substeps under the 20 000 cap): a tiny follow-on increment of 1e-6 / 1e-8 /
  1e-10 gives `|Δσ|/|σ| ÷ d_eps` = 145 / 152 / 152 for `-implex`, identical for the implicit
  twin (bit-identical committed state, expected with no free DOFs), and 145 / 146 / 146 on an arm that set `stressCorrection` to 0 via `setParameter` — **which is a
  silent no-op on this build** (quirks row 2026-09-07: `updateParameter` reads `theInt`), so that
  arm is NOT a stress-correction-off control and its small differences from arm A are unexplained
  (a fresh material instance, not the parameter). The linear, on-surface signature over three
  decades stands on the `-implex` and implicit arms; the correction's role is UNTESTED. (A first attempt on the free-DOF settlement column was not evidence: a
  `dt = 0` hold measures no `implexError` by construction, the hold moves nodes to close the
  committed state's equilibrium gap, and that column hits the companion cap on 19/30 commits
  even at 20 000 — three quirks rows written.) **Loose Esmeralda leg 146456** ended at s/B
  0.0397 on budget: clean to 0.0394, then O(1) committed errors by the dozen (0.2–5.6 on ~250
  of the last 660 rows), Q wandering 2527 → 2462 → 2554 kN over 0.0004 s/B, and the first two
  companion cap hits at 20 000.
  **What the arithmetic now says.** On a committed tiny step after a full one, the error is
  `≈ f·‖Ce:Δε_p(n)‖ / ‖σ‖` with `f = dt_{n+1}/dt_n`. At 3.9e-7 after a 1e-4 step, `f ≈ 4e-3`
  and a seat-sized `Δε_p(n) ~ 1e-3` give ~2e-3, not 1.11. An O(1) committed error at that
  step needs either `f ≈ 1` (the committed `dt_n` lost — the `mImplexDtCommit == 0 → f = alpha`
  fallback fired inside a refusal chain) or `‖Δε_p(n)‖ ~ 0.1–0.5` (a stale `epsPOld`, i.e. the
  history increment became a history *total*). Both are **bookkeeping under long refusal
  chains ending in the floor-accept branch** (`|dt| < reductionLimit·|dt0|` accepts whatever
  the error is) — a path the clean probe never exercised (zero refusals) and the P1 battery
  covers only for a single refusal (M10). The material at the seat is not implicated; ADR 92
  P2 is. **Discriminator in Esmeralda's argmax census:** `implexDetail[5]` (= `f`) on rows
  315 / 1130 / 1142 / 1290: `f ≈ 1` names the clock bookkeeping, `f ≈ 4e-3` names the
  history term (then `pstrains` differences across the chain give `‖Δε_p(n)‖`).

