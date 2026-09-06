---
title: "The Band You Cannot Declare — regularization: what we measured, and what we shipped instead"
project: Ladruno
status: issued
date: 2026-09-06
audience: TIMs project team
tags:
  - report
  - evidence
  - regularization
  - localization
  - sanisand
  - integrator
---

# The Band You Cannot Declare

**Answer to fork item 10 (vault 64 §6) and to the regularization decision of vault 65 D5.**

This is written for the people who will use the result, not for the people who wrote the C++.
It leads with the only question you actually have — *do I need to change anything, and is my
P1 affected?* — and answers it before explaining anything.

---

## The short version

1. **Item 10 as framed is not achievable.** A rate-regularized, τ-declared *band width* does not
   exist in quasi-statics. We measured why, and the reason is not an implementation detail — it
   survives every variant of the method, including the one you asked for.
2. **You very probably do not need it.** Under your own ultimate criterion (ADR 65 D6: `q` at a
   fixed `s/B` when no peak forms), the quantity that feeds P1 **already converges under mesh
   refinement with no regularization at all**.
3. **A regularizer would have cost you accuracy, not bought it.** At the only Deborah numbers
   where it bites, it inflates the collapse load by **+10 to +17 %** — the same order, and a
   fixed sign, as the 11.75 % Rayleigh artifact you spent notes 13/14 identifying and removing.
4. **What actually blocks your campaign is the constitutive integrator, and we fixed that.**
   Not ill-posedness. Your legs were dying inside `ModifiedEuler`, invisibly.

Nothing in your existing decks is invalidated by this report. Three things in §5 change how you
should *write* decks from here on, and one of them (§5.1) means your past runs were solving with
a different tangent than you thought.

---

## 1. Your deliverable already converges

We ran the τ = 0 (unregularized) strip footing on the element ADR 65 D4 freezes
(`LadrunoBrick -formulation bbar`), at three resolutions (`h0` = 1.0 / 0.5 / 0.25) and two
densities, and compared the load at **matched settlement** — which is what D6's criterion asks
for when no peak forms.

| density | matched-settlement load band, coarse → fine |
|---|---|
| `e_init = 0.6944` (Gorini's) | **7.80 % → 5.28 % → 2.86 %** |
| `e_init = 0.60` (dense) | **5.42 % → 3.55 %** |

Monotone, contracting, from above — the same qualitative behaviour as the b-bar hexahedron on
the exact Prandtl–Reissner problem (note 71). **This is the P1 input, and it is convergent
without a regularizer.**

Two caveats we will not hide: the `s/B = 0.002` row is **not resolved** (replicate scatter on
relaxed-tolerance steps is 5.93 %, against 0.000 % on clean steps), while `s/B = 0.005` and
`0.010` are; and the solver-configuration floor is 0.8–1.4 %. And **we cannot tell you whether
this is "inside tolerance"** — that number is your OQ2 and it has not been supplied. We did not
substitute one.

## 2. The band width does not converge, and no viscosity fixes that

Same runs, measuring the localization width instead of the load: the refinement ratio
`w(h/2)/w(h)` sits at **0.675–0.824** at every settlement and both densities (1.0 would be
converged, 0.5 is "one element"). Yielding *volume* is far closer to objective (3.95 % spread)
than width — which is precisely why the energy-based `lch` gates the fork already ships cannot
serve this consumer.

The decisive measurement is the next one. In a controlled 1-D study, **at fixed De**, the
converged width moves:

- **×4.3** with the imperfection *amplitude*
- **×3.5** with the imperfection *zone length*

So viscosity does stop the band collapsing onto a single element — but **what it collapses to is
chosen by the imperfection field, not by τ**. A "band width of τ" is not a thing that can be
delivered. This result is independent of which viscoplastic variant is used, so it applies to
the true Duvaut–Lions formulation as much as to any wrapper.

This is de Borst, Sluys, Mühlhaus & Pamin (1993) — rate regularization is a *transient* remedy
whose effect "rapidly decreases for slow loading rates" — now with numbers on your material.

## 3. The regularizer would have biased P1 upward

From the same study, the peak load relative to the unregularized run:

| De | 3e-4 | 1e-3 | 3e-3 |
|---|---|---|---|
| bias in the load | **+9.7 %** | **+12.5 %** | **+17 %** |

3e-4 is the *lowest* De at which the width converges at all — below it, nothing converges. So
there is no operating point that regularizes without inflating. Compare note 14: the artifact
you removed was 11.75 %, and your own conclusion there was that targeting it "is not defensible
if the macroelement is meant to reproduce the *soil*." The same sentence applies here, with the
sign fixed and the parameter chosen by the analyst.

One more property worth knowing: the limits **do not commute**. For any fixed mesh there is a De
below which the answer is mesh-locked again, so "refine the mesh and lower De" is not a path to
a limit — it is a two-parameter surface you would have to report.

## 4. What was actually stopping you — and what we shipped

All six of our τ = 0 legs **seized far short of any peak** (deepest `s/B` = 0.0228 against a
0.25 target, still hardening at 10.7 kPa/mm). The controller was not the problem: every leg used
**0 of its 80** available subdivisions, ending with its step 6400–25000× above the floor. The
binding resource was inside the material — `ManzariDafalias`'s substepped `ModifiedEuler`
collapses toward `dT_min = 1e-6` with **no limit on substep count**, and on reaching the floor it
*force-accepts* and reports success. One `analyze(1)` took **34.3 minutes** — 59 % of that leg's
whole budget in a single step — and the analysis was never told anything was wrong.

Shipped (PR #792, on `ladruno`):

- **`-maxSubsteps N`** on `LadrunoSANISAND` — bounds the substep count for one material update.
  Default **0 = uncapped = exactly today's behaviour**, so nothing you have changes unless you
  ask for it. Past the cap the material **refuses the increment** (it does not force-accept),
  leaves its committed state untouched, and says so — which lets your step controller cut the
  load step, the tool that was being denied the information it needed.
- **The refusal is now carried** by `LadrunoBrick` on every path, and `LoadPath` / `ArcLength`
  no longer discard it (`LoadControl` / `DisplacementControl` / `Newmark` already checked).
- **Tolerance guidance** (§5.2) and the **`-Presidual` decision brief** (§5.3).

Measured effect on your own deck, coarse legs, both densities:

| | before | after |
|---|---|---|
| subdivisions used | 0 / 80 | 2, 16, 18, 81 |
| depth reached | — | **2.1–2.6× deeper** |
| worst single step | 759 s | **94 s** |
| matched-settlement `q` | — | **bit-identical** where both runs reach, 0.28–0.94 % beyond (inside the solver floor) |

**And still no peak, on any leg.** The integrator fix buys reach and honesty, not a mechanism.
If your ultimate criterion ever requires a genuine limit point on this problem, that is still
not reachable and is a separate question; D6's settlement criterion is the one that works today.

## 5. Three things that change how you write decks

### 5.1 The tangent you were using was not the one you thought

`LadrunoSANISAND`'s parser defaulted `TanType` to **0 — the elastic tangent**. Every fork
SANISAND deck that omitted the positional argument has been running `algorithm Newton` as
*modified* Newton. It is invisible on a single-element, fully-prescribed deck (no global solve),
and expensive on a boundary-value problem: we measured **800 vs 283 Newton iterations** on a
drained triaxial (2.83×), and at a tighter tolerance the elastic-tangent leg cannot finish the
push at all while the consistent one can.

The default is now **2 (consistent)** and the material echoes which tangent it chose. Vanilla
`ManzariDafalias` keeps its own default of 0 — deliberately, so the archived record stays
citable — so the two now differ. **Emit `TanType` explicitly in every deck.** And note that the
consistent tangent of a non-associated model is genuinely unsymmetric: pair it with
`FullGeneral`, `UmfPack` or `Pardiso -matrixType 0`, never a symmetric solver. The converged
answer is unchanged by the switch (we gated that); only the path to it is.

### 5.2 `NormDispIncr` is unreachable on SANISAND

The displacement-increment residual stalls at a median **6.6e-7 m** and never reaches a tight
tolerance — and it is not mesh-neutral, so the same tolerance means different things at
different `h`. Use `NormUnbalance` scaled to the model's own weight (`γ'V`), as the R3 gate does.

### 5.3 `-Presidual 0` has a visible onset, and it is your decision

At the fork's `-Presidual 0`, a free-surface Gauss point gets clamped at `s/B ≈ 0.0153` on the
dense coarse leg — past every number we report, but a faster deck reaches it. Three options
(declare a non-zero `p_r`; add a surcharge; accept the clamp and disclose it) are tabulated in
the engineering handoff with a fork-side recommendation. We did **not** change the default.

### 5.4 One boundary on the new flag

`-maxSubsteps` only does what it says under an element that propagates a material refusal —
today that is `LadrunoBrick`. Under `stdBrick`, `BrickUP` or the `QuadUP` family the refusal is
silently swallowed and the analysis would accept a partially-integrated state. The material
prints this precondition when the cap fires. **Do not combine `-maxSubsteps` with those
elements.**

---

## 6. What we need from you

- **OQ2 — the tolerance.** Band width relative to `B`, and the band on `q` you are willing to
  accept. Without it, §1's contracting sequence cannot be declared "converged enough", and no
  one else can declare it for you.
- **OQ1 — the ultimate criterion** (still Prof. Gorini's). It scales P1 by up to 3.3× and is
  unchanged by anything here.
- **A decision on ADR 65 D5.** On this evidence, "Duvaut–Lions wrapper with declared τ" is not
  the right instrument. Our recommendation, for your side to accept or reject: amend D5 to
  **disclosure** — report the three-mesh band at matched settlement, with provenance, and drop
  the regularization requirement — and retire the acceptance case that was owed to us, since
  there is no longer a wrapper for it to accept. If you would rather keep a conditional path
  open, the honest form is: *revisit only if a genuine limit point becomes reachable and the
  matched-settlement band then falls outside OQ2's tolerance.*

## 7. Where the numbers live

Everything quoted here is on `ladruno`, with the driver, the per-step curves and the field dumps
committed beside the prose:

| what | where |
|---|---|
| the decision, the theory, the alternatives | `Ladruno_implementation/90_ladruno_viscoplastic_regularization_adr.md` |
| the 1-D study (width vs imperfection, bias, De window) | `_adr90_a0_results.md`, `_adr90_p0b_results.md` |
| the τ = 0 strip-footing campaign, and the capped re-run | `_adr90_tau0_qu_band.md` (§7b) |
| the driver and its run records | `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py`, `adr86b_t8/` |
| the integrator change | PR #792; `LEDGER_implementations.md` (ADR-86b row) |

Every number names the engine build it was taken on. The campaign engine was pinned by git hash
and the driver aborts on a mismatch — the provenance lesson from the T1 incident is wired in, not
remembered.

---

> **One sentence.** Rate regularization cannot give you a band width you declare, your collapse
> load did not need it, and the thing that was actually costing you was an integrator that could
> spend 34 minutes on one step and call it success — which is now capped, honest, and 2.5× faster
> to the same answer.
