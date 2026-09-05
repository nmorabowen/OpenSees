---
title: "ADR 90 / WP-A2 — the tau = 0 collapse-load band on softening SANISAND"
project: Ladruno
type: measurement report (ADR-90 un-descope gate)
status: "RESULTS 2026-09-05 — verdict in §1"
owner: nmora
adr: ADR 90
engine_build: 548fe911427e90a2edfead05cb3672a738d25b6d
driver: Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py
related:
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[_adr90_regularization_planning_brief]]"
  - "[[_adr90_a0_results]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[86_ladruno_sanisand_apegmsh_emitter_guide]]"
  - "[[testbed/00_canonical_testbed]]"
  - "[[LEDGER_quirks]]"
tags: [adr90, regularization, sanisand, localization, collapse-load, measurement]
updated: 2026-09-05
---

# ADR 90 / WP-A2 — the tau = 0 collapse-load band on softening SANISAND

> [!abstract] The question, and why it is the un-descope gate
> ADR 90's warrant rests on a premise nobody had measured: that the
> **unregularized** collapse load of the campaign's own problem is *not*
> already mesh-converged inside the campaign's own tolerance. The fork has that
> measurement for `DruckerPrager` — `tests/test_r3_prandtl_collapse_gate.py`
> reads 1.0849 / 0.9977 / 0.9514 of exact Prandtl–Reissner at
> h0 = 1.0 / 0.5 / 0.25, every leg a genuine plateau — and did **not** have it
> for `LadrunoSANISAND`, whose non-associated **softening** is the reason
> ADR 90 exists.
>
> **Nothing in `SRC/` was changed and no build was made.** This measures ADR 90
> §3.1 option (d) — *"do nothing to the lane; report the band"* — which that
> section itself calls *"the negative control of the whole ADR, not an option"*.

---

## 1. VERDICT

> **τ = 0 `q_u` does not converge — because on this deck it cannot be MEASURED
> AT ALL. Every leg SEIZES, on every mesh, at both densities, far short of its
> peak.** The deepest leg reached `s/B = 0.0203` of a `0.25` target, still
> hardening at ~10 kPa/mm. **No leg is a capacity**; nothing below may be quoted
> as a collapse load.

The seizure is **not** the stepping controller. Every leg used **0 of its 80
pinned subdivisions** and ended with its step still 800–12800× above the floor;
the one leg that reached its own wall budget cleanly did so with `free-advance =
yes`. The binding resource is the **constitutive integrator**: single
`analyze(1)` calls cost **11–28 minutes**, because SANISAND's substepped
`ModifiedEuler` return collapses to its `dT_min = 1e-6` floor
(`ManzariDafalias.cpp:1380`, `while (T < 1.0)`) at Gauss points inside the
developing plastic zone, and the Newton ladder then repeats that up to 125
state-determination passes per step. That is a **seizure** by
`00_canonical_testbed` §1c, and it is the boundary-value-problem form of exactly
what ADR 90 §1.1 already reports for the A0 bar: *"the number of uniform load
steps the rate-independent problem needs in order to converge at all grows
4000 → 64000 across the same four meshes, while **every** regularized run
completes at 250 steps on **every** mesh."*

**What IS measured — at matched settlement, pre-peak:**

| quantity | e_init = 0.6944 (Gorini) | e_init = 0.60 (dense) |
|---|---|---|
| three-mesh band in `q` at s/B = 0.002 | **7.80 %** | **5.42 %** |
| … at s/B = 0.005 | **5.28 %** | **3.55 %** (2 meshes) |
| … at s/B = 0.010 | **2.86 %** (2 meshes) | — |
| width ratio `w2(h/2)/w2(h)` | **0.675 / 0.802** (0.002), **0.692 / 0.811** (0.005), **0.722** (0.010) | **0.702 / 0.824** (0.002), **0.722** (0.005) |

So: **the pre-peak LOAD contracts under refinement (7.8 % → 5.3 % → 2.9 %),
while the localization WIDTH does not converge at any settlement or either
density** — `w2(h/2)/w2(h)` sits at **0.68–0.82** throughout, against **1.0** for
a converged width and **0.5** for a strictly one-element band. Width convergence
and load convergence are different gates, which is precisely ADR 90 §3.2's own
wording, now confirmed on the consumer's own material and element.

**And a second, independent obstruction.** The dense coarse leg reaches the
material's low-p floor at `s/B ≈ 0.0153` and starts logging *"the result at this
integration point is set by the clamp, not by the model"* (§3.1) — the
free-surface low-p pathology the WP brief named in advance as a candidate
finding, appearing exactly where predicted. It is past every number reported
here, but it means a faster deck would hit it next.

**Against the R3 comparison basis, stated exactly.** R3's per-resolution band is
**± 3 %**, but that is a *tolerance on a known centre*, not a convergence
criterion: R3's own three-mesh spread is 1.0849 → 0.9514, i.e. **13.1 % of their
mean**, and is accepted because it is monotone and every leg is a genuine
plateau. The bands above (2.9–7.8 %) sit inside that spread — **but they are
pre-peak load values, not collapse loads, so the comparison is context and not a
verdict.** **The campaign's own tolerance is OQ2 and has NOT been supplied on
either side** (ADR 90 §7.5: TIMs must set the target band width relative to B,
the ramp duration, the De and the tolerance bands). Any statement of the form
*"inside the campaign's tolerance"* is therefore **unavailable by construction**,
and this report does not substitute a number of its own.

**Consequence for ADR 90.** The wrapper's warrant is **strengthened, and on a
different axis than expected**:

* **C5** (*collapse load converges in h at fixed De*) cannot be tested against a
  τ = 0 control, **because the τ = 0 control does not terminate**. The negative
  control of ADR 90 §3.1(d) is not merely worse than the regularized lane — on
  this deck it is unrunnable. A regularizer that lets the problem *finish* is a
  deliverable before it is an accuracy claim.
* **C4** (*width converges at fixed De, and does not at τ = 0*) has its negative
  half **MEASURED here**: 0.68–0.82 at every settlement and both densities, on
  `LadrunoBrick -formulation bbar` — the element ADR 90 S3/D7 freezes. The
  fork's existing `lch` energy gates cannot see this (brief F6), and now there
  is a number instead of an argument.
* The claim ADR 90 may **not** make from this WP is that τ = 0 `q_u` is far from
  converged *in value*. It is not measured. What is measured is that it is not
  **reachable**.

---

## 2. What was run

**Engine.** `ladrunoBuild()` = `548fe911427e90a2edfead05cb3672a738d25b6d`, whose
only C++ difference from `origin/ladruno` is a comment. The binary is bound
explicitly (`sys.path.insert` on that worktree's `dist/bin`), the hash is
**asserted** at the top of every leg, and it is stamped into every leg JSON and
every CSV header. Nothing was built for this WP.

**Driver / summariser.** `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py`
and `sanisand_tau0_summary.py`. Raw output — curve CSVs, per-checkpoint
plastic-strain field CSVs, leg JSONs, engine logs — in
`Ladruno_files/testbed/hypo_bearing/adr90_a2/`, with the first launch kept as an
independent replicate in `adr90_a2_run1/` and the tolerance control in
`adr90_a2_tolref/`.

**Deck.** Plane-strain strip footing, `LadrunoBrick -formulation bbar` (ADR 90
S3), one element thick, `u_y = 0` on every node. The mesh, the domain (B = 2 m,
14.5 B of clearance, 10 B of depth, grading 1.35) and the ADR-63 D16 stepping
guards — *including the pinned `SUBDIV_BUDGET = 80` and its calibration history*
— are **imported** from `tests/test_r3_prandtl_collapse_gate.py`, not copied.
`_PARAMS` is Gorini's calibrated set, **imported** from
`tests/test_ladruno_sanisand.py`.

| | h0 = 1.0 | h0 = 0.5 | h0 = 0.25 |
|---|---|---|---|
| hexes | 200 | 390 | 782 |
| DOF | 1386 | 2592 | 5040 |

**Three deliberate divergences from R3**, each because R3's choice is right for
Prandtl–Reissner and wrong here:

1. **Weighted soil, not weightless.** R3 runs `gamma = 0` under a surcharge
   because Prandtl–Reissner *is* a weightless problem with an exact answer.
   SANISAND cannot be run that way: `G ∝ sqrt(p)` and its whole state (`psi`,
   `M^b`, `M^d`) is pressure-dependent, so `p → 0` is a pathology, not a
   simplification. Self weight `gamma = rho·g = 2.0 × 9.81 = 19.62` kN/m³ via
   `eleLoad -type -selfWeight` with `-b 0 0 -gamma`.
2. **Rough footing** (`u_x = 0` on the footing nodes). R3's is smooth because
   the exact Prandtl answer is the smooth one; there is no exact answer here, so
   the deck uses the rough footing the campaign's problem actually has.
   **No number in this report may be compared to R3's.**
3. **`DROP_TRUNCATE` disabled.** R3 cuts the curve where `q` first falls 2 %
   below its max — a bottomed-out-run detector on a plateau problem, and a
   *signal deleter* on a softening one.

**Staging** follows all three ADR-86 emitter-guide deck rules
(`86_ladruno_sanisand_apegmsh_emitter_guide` §2): gravity ramped over 10 steps
at material stage 0; then `updateMaterialStage -material 1 -stage 1` explicitly
for the tag; then the push. `-Presidual 0.0` and `-Pmin 1e-3·P_atm` are the
ADR-86 fork defaults, emitted explicitly.

**Two densities.** Gorini's `e_init = 0.6944`, and a denser `e_init = 0.60`. At
`p ≈ 20` kPa the state parameter is `psi = −0.123` vs `−0.217`, so `M^b/M_c` is
1.54 vs 2.14 and `M^d/M_c` is 0.49 vs 0.29 — the dense leg dilates hard and
softens hard.

**Push.** R3's pattern exactly: `sp` inside a load pattern under one-argument
`LoadControl`, pseudo-time *is* the imposed footing settlement, adaptive
subdivision against the pinned budget. Target `s/B = 0.25`, with an early stop
the moment a peak is passed (§4). `DS_MAX = 5 mm` is pinned here, not inherited
(R3's 1 mm is sized for a 0.30 m target on a two-orders-cheaper material).

---

## 3. Controls

Every control passes on every leg, by wide margins.

| leg | resultant identity | 1-D geostatic patch | η/M_c at flip | `Outside Bounding` | low-p `CLAMPING` | p_min (kPa) |
|---|---|---|---|---|---|---|
| h1.0_e0.6944 | 4.44e-16 | 1.15e-13 | 0.6425 | 0 | 0 | 6.248 |
| h1.0_e0.60 | 4.44e-16 | 1.77e-13 | 0.6425 | 0 | **10** (§3.1) | 6.248 |
| h0.5_e0.6944 | 4.44e-16 | 5.36e-13 | 0.6425 | 0 | 0 | 3.124 |
| h0.5_e0.60 | 4.44e-16 | 4.49e-13 | 0.6425 | 0 | 0 | 3.124 |
| h0.25_e0.6944 | 4.44e-16 | 1.68e-12 | 0.6425 | 0 | 0 | 1.562 |
| h0.25_e0.60 | 4.44e-16 | 2.28e-12 | 0.6425 | 0 | 0 | 1.562 |

*(the summariser reads these counts from the engine log at leg end; for the five
legs killed mid-leg it prints `-1` = "not read", so the counts above were
grepped from the engine logs directly.)*

* **Resultant identity** `ΣR_z(base) = γ·V` — 4.44e-16 on every mesh. It is a
  cheap catch for the body-force convention (`LadrunoBrick -b` is per unit
  *volume*, never multiplied by rho — [[LEDGER_quirks]]) and nothing more:
  §1d of the testbed doc is explicit that it is **structurally blind** to
  load-distribution error.
* **The 1-D geostatic patch** is the field check that identity cannot make. A
  laterally-rollered, base-fixed box under self weight has the closed form
  `sigma_zz = gamma·z`, `sigma_xx = sigma_yy = K0·sigma_zz`, `K0 = nu/(1−nu)` —
  exact even for a pressure-dependent material, because `nu` is constant so `K0`
  is depth-independent. Under `bbar` the field comes back **piecewise constant
  at the element centroid value**, matched to 1.15e-13 … 2.28e-12 over every
  Gauss point of every element.
* **This is also the gravity-state MESH-SENSITIVITY measurement**, and it is the
  strongest form available: it compares every Gauss point against a
  *mesh-independent closed form* rather than comparing meshes to each other, so
  the whole three-mesh spread is bounded by **2.3e-12**, twelve orders inside
  the 0.1 % requirement.
* **η/M_c = 0.6425 at the stage flip on every mesh and both densities**, and
  **zero** `Outside Bounding` events — so `Elastic2Plastic` never inflated the
  calibrated `M_c` and every leg ran Gorini's soil. At the gravity state the
  shallowest Gauss point sits at 1.56 kPa on the finest mesh, 15× the
  `-Pmin + -Presidual` floor of 0.101 kPa; five of the six legs then fire **zero**
  low-p clamps for the whole push, and the sixth is §3.1.
* **Footing-vs-base reaction cross-check.** `q_foot ≡ −q_base` by global
  equilibrium during the push. Measured max mismatch 6.2e-3 on the leg that
  terminated cleanly. It is not decoration: the mismatch spikes to ~3.5e-2 on
  steps carried by the relaxed ladder rung, and those are exactly the steps that
  fail to reproduce (§5.4).

### 3.1 The low-p clamp DOES fire — on the dense coarse leg, past s/B = 0.0153

One control is **not** clean, and it is the one the WP brief named in advance as
a candidate finding ("SANISAND low-p at the free surface"). The
`h1.0_e0.60` leg logged **10** copies of

> `ManzariDafalias::ModifiedEuler() - material tag 1: mean stress p = 0.100973
> is below the floor m_Pmin + m_Presidual = 0.101; CLAMPING the stress to
> p = 0.101 (deviator preserved). **The result at this integration point is set
> by the clamp, not by the model.**`

**Where, exactly.** All ten fire at load factor `−0.0420013`, i.e. at
`s = 0.0306 m`, `s/B = 0.0153` — on the step *after* the leg's last converged
row, during the grinding step that was in progress when the process was reaped.
So they are past the deepest converged point of the deepest dense leg, and
**every checkpoint reported in §6 for that leg (s/B ≤ 0.010) is upstream of
them.** No number in this report is set by the clamp.

**Why it matters anyway.** It is a *second, independent* obstruction, and it is
not a solver artefact: at `-Presidual 0.0` (the ADR-86 fork default) the floor
is `-Pmin = 1e-3·P_atm = 0.101` kPa, and the dilating dense sand drives a
free-surface Gauss point to it. Beyond `s/B ≈ 0.0153` the dense coarse leg's
answer is partly the clamp's, not the model's. **So even a deck fast enough to
reach a peak would need this settled first** — by a declared `-Presidual`, a
small surcharge, or an accepted and disclosed clamp — and that is a prerequisite
for any τ = 0 SANISAND collapse load, regularized or not.

That the loose legs and both finer meshes fire **zero** clamps is the
discriminator: this is a *dense-dilatant, coarse-element, large-settlement*
effect, not a property of the deck as such.

---

## 4. Capacity: the three-clause rule, adapted for softening

`00_canonical_testbed` §1c requires three clauses, and clause 1 — *a flat tail* —
is written for a **perfectly plastic plateau**. A softening material does not
plateau; it **peaks and descends**. Clause 1 is therefore satisfied here by
*either* a plateau (|tail dq/ds| < 2 % of the initial tangent) *or* a **passed
peak**: `q` has fallen ≥ 5 % below `q_max` and stayed there for ≥ 15 further
steps. One dip is a solver artefact; a sustained descent is a post-peak branch.
`PEAK` joins `TARGET` and `BUDGET` as an admissible termination mode. `FLOOR`,
`WALL` and `KILLED` remain **seizure** and are never a capacity, whatever the
tail says.

**No leg satisfied clause 1 by either route.** All six are pre-peak.

---

## 5. Two defects in the deck that were measuring the solver, not the soil

Both were found by instrumenting the run rather than by reading its curve, and
both are on [[LEDGER_quirks]]. The wall times in §6 are the times **after** these
were fixed, and the first one silently affects every existing fork SANISAND deck.

### 5.1 `TanType` defaults to the ELASTIC tangent

`ManzariDafalias3D::getTangent()` returns `mCe` for `TanType == 0`, `mCep` for
`1`, and `mCep_Consistent` otherwise (`ManzariDafalias3D.cpp:134-142`). The
**parser default is 0** (`LadrunoSANISAND.cpp:117`, deliberately copied from
`OPS_ManzariDafaliasMaterial` so a renamed deck behaves identically), while the
**null and parallel constructors default to 2** (`ManzariDafalias.cpp:365`,
`:426`). A deck that emits only the 18 positional parameters therefore hands
`algorithm Newton` an elastic tangent — modified Newton, with linear
convergence, dressed as full Newton.

It hides because **on a zero-free-DOF material-point deck it costs nothing**, and
every existing fork SANISAND deck is of that shape
(`tests/test_ladruno_sanisand.py` passes `*_PARAMS` plus flags and no positional
optionals). On a boundary-value problem it was the difference between a
measurement and a non-measurement: at the parser default the h0 = 0.25 leg
advanced **11 steps to s/B = 1.6e-4 in 350 s**. This deck now emits the full
positional block with `TanType 2`; `mCep_Consistent` is unsymmetric under
non-associated flow, which is why the solver is `Pardiso -matrixType 0`.

### 5.2 The convergence test could not be met

The WP brief asks for `NormDispIncr` "per the R3 gate". **The R3 gate does not
use `NormDispIncr`** — it uses `NormUnbalance` at `1e-5 × the applied load`. And
on SANISAND the displacement form is not merely tight but *unreachable*: the
substepped `ModifiedEuler` return (hardcoded substep tolerance `TolE = 1e-4`)
makes the discrete stress–strain map only piecewise smooth, so Newton stalls
rather than converging quadratically. Measured on the h0 = 0.5 leg over 47 failed
attempts, the residual displacement norm stalls at a **median of 6.6e-7 m** (min
3.4e-8, max 1.3e-4). The run that nominally used `1e-8` was carried by the
relaxed third ladder rung on **18 of its 26 steps** — 65 of every 125
state-determination passes spent failing two rungs that could not succeed.

A displacement-norm tolerance is also **not mesh-neutral**, which is
disqualifying inside a mesh-convergence study: the norm runs over the free DOFs,
and there are 3.6× more of them at h0 = 0.25 than at h0 = 1.0, so the same
nominal number is a different physical requirement on each mesh of the sequence.
The push therefore uses **`NormUnbalance` at a tolerance relative to the model's
own weight `gamma·V`**, identical on all three meshes by construction.

### 5.3 The A/B that pinned the settings

h0 = 0.5, 400 s of wall each, same box, same `DS_MAX`:

| test / tolerance | reach s/B | steps | subdivisions | relaxed |
|---|---|---|---|---|
| `NormDispIncr` 1e-8 m (TanType 0) | 0.00106 | 26 | 1 | 18/26 |
| `NormUnbalance` 1e-6 ·γV (TanType 2) | 0.00218 | 31 | 0 | 11/31 |
| `NormUnbalance` 1e-5 ·γV (TanType 2) | 0.00442 | 37 | 0 | — |

`1e-5 ·γV` is pinned because it is the only setting under which the step reaches
`DS_MAX` at all, and *reach* is the binding constraint on whether this study can
answer its own question. One leg was re-run at `1e-6 ·γV` as the tolerance
control, and because both runs snapshot at the same checkpoint list the
comparison is **interpolation-free, at exactly matched settlement**:

| s/B (actual) | q at 1e-5 ·γV | q at 1e-6 ·γV | difference |
|---|---|---|---|
| 0.002180 | 53.0409 | 52.6062 | **0.82 %** |
| 0.005060 | 115.6139 | 114.7235 | **0.77 %** |
| 0.010180 | 225.5879 | 222.4336 | **1.42 %** |

There is also a reading across a *complete* change of solver configuration
(TanType 0 + `NormDispIncr 1e-8` → TanType 2 + `NormUnbalance 1e-5·γV`): `q` at
the matched `s/B = 0.002` checkpoint, h0 = 1.0, read **56.045** vs **55.756**
kPa — **0.52 %**.

> [!warning] Provenance caveat on that one number
> The `56.045` comes from the **pre-fix campaign attempt**, whose output
> directory was **deleted** when the settings were re-pinned and the campaign
> relaunched. It is a measurement that was made and logged, but it is **not
> reproducible from the data committed with this report**, and by this fork's
> own rule (ADR 90 S4, *provenance is output*) it must not carry the same weight
> as the rest. **The reproducible tolerance evidence is the three rows above**,
> which live in `adr90_a2_tolref/` and can be re-derived from the committed
> checkpoint CSVs. Deleting a superseded output directory instead of renaming it
> was the mistake; `adr90_a2_run1/` is kept for exactly this reason.

**So the solver-configuration floor on any band quoted here is ~0.8–1.4 %**
(from the on-disk tolerance control), and no band below 1.4 % should be read as
a physical difference.

### 5.4 Reproducibility, and where it fails

The campaign was launched twice (different thread counts: 4 processes ×
`OMP_NUM_THREADS=4`, then 6 × 3). Comparing every checkpoint the two launches
share:

| leg | s/B | run 1 q | run 2 q | Δ | note |
|---|---|---|---|---|---|
| h1.0_e0.6944 | 0.002 / 0.005 / 0.010 / 0.020 | 55.756 / 120.602 / 232.128 / 450.392 | identical | **0.000 %** | clean steps |
| h0.5_e0.6944 | 0.002 / 0.005 | 53.041 / 115.614 | identical | **0.000 %** | clean steps |
| h0.25_e0.6944 | 0.005 | 113.056 | 114.395 | **1.18 %** | relaxed rung |
| h0.25_e0.6944 | 0.002 | 54.819 | 51.570 | **5.93 %** | relaxed rung; run 1 also carried a 3.5 % footing-vs-base mismatch on that step |

`w2` behaves the same way (0.000 % on clean steps, 0.58 % and 1.45 % on the two
relaxed ones). **The deck is deterministic; the scatter enters only through the
relaxed ladder rung**, where a threaded PARDISO factorization that is not
bit-reproducible across thread counts is amplified by a marginally-converged
step. The footing-vs-base mismatch is the diagnostic that flags those steps.

**This bounds the resolution of the whole study**: a band of 5.4 % at
`s/B = 0.002` is *not* clearly resolved above a 5.9 % run-to-run scatter on the
finest mesh's relaxed step. The bands at `s/B = 0.005` and `0.010`, where the
shared checkpoints reproduce to 0.000–1.18 %, **are** resolved.

---

## 6. Results

### 6.1 The 3 × 2 table — every leg, and why none is a capacity

| leg | DOF | q_max (kPa) | s_end/B | of target | mode | plateau | peak | ds/floor | subdiv | steps | wall (s) | CAPACITY |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h1.0_e0.6944 | 1386 | 503.90 | 0.0228 | 0.25 | KILLED | NO | NO | — | **0**/80 | 50 | 2401 | **NO** |
| h1.0_e0.60 | 1386 | 621.91 | 0.0153 | 0.25 | KILLED | NO | NO | — | **0**/80 | 44 | 951 | **NO** |
| h0.5_e0.6944 | 2592 | 198.88 | 0.0089 | 0.25 | **WALL** | NO | NO | **12800** | **0**/80 | 43 | 3503 | **NO** |
| h0.5_e0.60 | 2592 | 286.43 | 0.0076 | 0.25 | KILLED | NO | NO | — | **0**/80 | 38 | 1414 | **NO** |
| h0.25_e0.6944 | 5040 | 114.40 | 0.0051 | 0.25 | KILLED | NO | NO | — | **0**/80 | 38 | 2240 | **NO** |
| h0.25_e0.60 | 5040 | 138.03 | 0.0038 | 0.25 | KILLED | NO | NO | — | **0**/80 | 31 | 2295 | **NO** |

`q_max` is **where the run stopped, not a capacity** — every leg is still
hardening steeply. `KILLED` = the harness reaped the background process at
~1 hour; `WALL` = the leg's own 3300 s budget bound. Both are seizure modes and
neither is admissible. **Note the `subdiv` column: 0/80 on every single leg.**
The pinned controller budget — the guard R3's docstring exists to protect — was
never touched. Nothing about this result is a stepping artefact.

### 6.2 The matched-settlement tables — the deliverable that survives

Legs on different meshes stop at different settlements, so every leg snapshots
`q` and the width at the **same** `s/B` list. Values are the union of both
launches (run 2 preferred where both exist).

**Load, `q_foot` (kPa), and the three-mesh band:**

| s/B | e0.6944: h1.0 / h0.5 / h0.25 | band | e0.60: h1.0 / h0.5 / h0.25 | band |
|---|---|---|---|---|
| 0.002 | 55.76 / 53.04 / 51.57 | **7.80 %** | 84.12 / 79.68 / 81.78 | **5.42 %** |
| 0.005 | 120.60 / 115.61 / 114.40 | **5.28 %** | 194.39 / 187.61 / — | **3.55 %** |
| 0.010 | 232.13 / 225.59 / — | **2.86 %** | 401.98 / — / — | — |
| 0.020 | 450.39 / — / — | — | — | — |

The `e0.6944` sequence is **monotone decreasing under refinement at every
settlement it is measured at** — the same "converging from above" shape R3 finds
for DruckerPrager — and the band **contracts** as the settlement grows.
`e0.60` at `s/B = 0.002` is **not** monotone (84.12 / 79.68 / 81.78), and that
non-monotonicity is inside the ±5.9 % relaxed-step scatter of §5.4, so it is not
interpretable either way.

**Localization width, `w2` (m) at z = −0.625 m, and its refinement ratios:**

| s/B | e0.6944: h1.0 / h0.5 / h0.25 | w(h/2)/w(h) | e0.60: h1.0 / h0.5 / h0.25 | w(h/2)/w(h) |
|---|---|---|---|---|
| 0.002 | 6.3801 / 4.3070 / 3.4563 | **0.675, 0.802** | 6.0533 / 4.2489 / 3.5027 | **0.702, 0.824** |
| 0.005 | 8.1332 / 5.6254 / 4.5649 | **0.692, 0.811** | 8.2060 / 5.9236 / — | **0.722** |
| 0.010 | 9.7337 / 7.0264 / — | **0.722** | 10.4772 / — / — | — |
| 0.020 | 11.3830 / — / — | — | — | — |

**A converged width would give ratios of 1.0; a strictly one-element band would
give 0.5. Every measured ratio is 0.675–0.824**, at every settlement, on both
densities — the width is still tracking the mesh. It is *not* one element wide
(`w2/h0` runs 6.4 → 8.6 → 13.8, i.e. the plastic zone is 6–14 elements across at
these settlements), because at `s/B ≤ 0.02` the zone is still a diffuse
pre-peak plastic bulb rather than a formed shear band. **This is the negative
half of ADR 90's C4, measured on the campaign's own element** — but it is a
*pre-peak* measurement, and a post-peak one is owed once the deck can reach a
peak at all.

**Yielding extent** (`eps_q^p > 1e-5`) — the count is not a mesh-convergent
quantity, the volume is:

| s/B | e0.6944 volume (m³): h1.0 / h0.5 / h0.25 | spread | element count |
|---|---|---|---|
| 0.002 | 126.68 / 111.73 / 112.45 | 12.5 % | 97 / 214 / 504 |
| 0.005 | 309.03 / 297.07 / 302.20 | **3.95 %** | 162 / 335 / 716 |
| 0.010 | 516.95 / 473.29 / — | 8.8 % | 194 / 379 / — |

The yielding **volume** is far closer to mesh-converged than the **width** is
(3.95 % at s/B = 0.005 against width ratios of 0.69/0.81 at the same
settlement). That is the WP-A2 form of the distinction ADR 90 §3.3 draws between
the fork's shipped `lch` gates and what this consumer needs: *how much soil
yields* is nearly mesh-objective already; *how thick the mechanism is* is not.

### 6.3 The seizure, quantified

| observation | measurement |
|---|---|
| longest single `analyze(1)` call | **~11 min** (h1.0_e0.6944, 83 s → 765 s of leg wall for one step at ds = 0.64 mm) and **20–28 min** on several legs late in the run |
| state-determination passes per step near the stall | up to **125** (25 Newton + 40 line-search + 60 Krylov) |
| substep floor the material collapses to | `dT_min = 1e-6` (`ManzariDafalias.cpp:1380`), i.e. up to **10⁶ return maps per Gauss point per pass** |
| subdivisions used, every leg | **0 / 80** |
| terminal step above the floor, cleanest leg | **12800×** |
| deepest settlement reached, any leg | **s/B = 0.0228** of a 0.25 target (**9 %**) |

**The two facts together are the finding**: the controller had every resource it
was given and never used any of it, and the run still could not advance. The
cost is inside the constitutive integrator, and it grows as the plastic zone
develops — which is the rate-independent ill-posedness making itself felt as
*non-termination* rather than as a wrong number.

---

## 7. What this does and does not license ADR 90 to say

**MAY say, on this evidence:**

* The τ = 0 negative control of §3.1(d) **does not terminate** on the consumer's
  own material, element and mesh sequence, inside a 55–90 min per-leg budget on
  a 24-core box, at either density. *(§1, §6.1, §6.3)*
* The τ = 0 localization **width does not converge** under h-refinement:
  `w(h/2)/w(h) = 0.675–0.824` at every settlement and both densities.
  *(§6.2)* — the negative half of **C4**, on `LadrunoBrick -formulation bbar`.
* The τ = 0 **pre-peak load** band contracts under refinement (7.8 % → 5.3 % →
  2.9 % for Gorini's `e_init`), and the sequence is monotone from above.
  *(§6.2)*
* The **yielding volume** is much closer to mesh-objective than the width is
  (3.95 % vs ratios of 0.69/0.81 at the same settlement), which is why energy-
  regularizing `lch` gates cannot serve this consumer. *(§6.2, brief F6)*

**MUST NOT say:**

* Anything about the **value** of the τ = 0 collapse load, or its convergence in
  value. `q_u` is not measured; no leg is a capacity.
* That any band here is "inside" or "outside" the campaign's tolerance. **OQ2 is
  unsupplied**; there is no tolerance to be inside of.
* That the 5.42 % non-monotone `e0.60` band at `s/B = 0.002` means anything —
  it is inside the 5.93 % relaxed-step scatter. *(§5.4)*
* That the contraction 7.80 → 5.28 → 2.86 % is a *rate* of convergence. Each
  band sits on a different number of meshes (3, 3, 2) and the last is only ~2×
  the 1.4 % solver-configuration floor of §5.3. The **direction** is evidence;
  the **exponent** is not.
* Anything comparing these numbers to R3's. Different problem (weighted vs
  weightless), different footing (rough vs smooth), different material.

**Owed next, in order:**

1. **A cheaper deck that can reach a peak.** The domain is R3's — 60 m × 20 m
   for a 2 m footing, sized for a *weightless* half-space. A weighted footing
   mechanism is far more local; a purpose-sized domain (and a `dT_min` /
   substep-cap study on `ManzariDafalias`) is the obvious lever, and neither
   touches the physics under test. Until then **every** τ = 0 SANISAND collapse
   number in this campaign is unreachable, not merely expensive.
2. **Settle the low-p floor before any deeper push** (§3.1). At the ADR-86
   default `-Presidual 0.0` the dense coarse leg is already being clamped at
   `s/B ≈ 0.0153`, and a clamped integration point is the model's answer no
   longer. A declared `-Presidual`, a small surcharge, or an explicit accepted
   clamp — but a choice, made and disclosed.
3. **A post-peak width measurement.** §6.2's ratios are pre-peak; the C4
   negative control ADR 90 actually wants is on a formed band.
4. **OQ2 from TIMs** — without it no band, regularized or not, can be adjudicated.

---

## 8. Reproducing

```bash
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_summary.py --selfcheck
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py \
    --out <dir> --legs h1.0_e0.6944 --wall 3300
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_summary.py <dir>
```

`--selfcheck` pins the width metric's calibration: a k-element top hat must read
exactly `k·hx` (measured exact to 1e-12 for k = 1, 2, 3, 5, 8 at two element
sizes). The six legs were run as six concurrent processes
(`OMP_NUM_THREADS=3` each, 24-core box), so **the wall times in §6.1 are upper
bounds under mild contention and are not a performance measurement**; the
*ordering* of the seizure evidence in §6.3 does not depend on them.

> [!warning] Infrastructure note
> Both campaign launches were **killed from outside at ~1 hour of background
> runtime**. The first lost everything it had computed, because `run_leg` wrote
> its JSON only on return; the driver now writes the leg record at **every**
> checkpoint, atomically, marked `partial: true`, and the summariser maps a
> partial leg to mode `KILLED` — a non-result, printed and excluded from any
> band, exactly as a wall-clock stop is. Budget any re-run in ≤ 50-minute
> chunks, or run in the foreground.
