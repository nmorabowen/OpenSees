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

## 1. VERDICT

RESULTS-PENDING

## 2. What was run

**Engine.** `ladrunoBuild()` = `548fe911427e90a2edfead05cb3672a738d25b6d`, whose only C++
difference from `origin/ladruno` is a comment. The binary is bound explicitly
(`sys.path.insert` on that worktree's `dist/bin`), the hash is **asserted** at the top of every
leg, and it is stamped into every leg JSON and every CSV header. Nothing was built for this WP.

**Driver / summariser.** `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py` and
`sanisand_tau0_summary.py`. Raw output — curve CSVs, per-checkpoint plastic-strain field CSVs, leg
JSONs, engine logs — in `Ladruno_files/testbed/hypo_bearing/adr90_a2/`, with the tolerance control
in `adr90_a2_tolref/`.

**Deck.** Plane-strain strip footing, `LadrunoBrick -formulation bbar` (ADR 90 S3 freezes the
element), one element thick, `u_y = 0` on every node. The mesh, the domain (B = 2 m, 14.5 B of
clearance, 10 B of depth, grading 1.35) and the ADR-63 D16 stepping guards — *including the pinned
`SUBDIV_BUDGET = 80` and its calibration history* — are **imported** from
`tests/test_r3_prandtl_collapse_gate.py`, not copied. `_PARAMS` is Gorini's calibrated set,
**imported** from `tests/test_ladruno_sanisand.py`.

| | h0 = 1.0 | h0 = 0.5 | h0 = 0.25 |
|---|---|---|---|
| hexes | 200 | 390 | 782 |
| DOF | 1386 | 2592 | 5040 |

**Three deliberate divergences from R3**, each because R3's choice is right for Prandtl–Reissner
and wrong here:

1. **Weighted soil, not weightless.** R3 runs `gamma = 0` under a surcharge because
   Prandtl–Reissner *is* a weightless problem with an exact answer. SANISAND cannot be run that
   way: `G ∝ sqrt(p)` and its whole state (`psi`, `M^b`, `M^d`) is pressure-dependent, so
   `p → 0` is a pathology, not a simplification. Self weight
   `gamma = rho·g = 2.0 × 9.81 = 19.62` kN/m³ via `eleLoad -type -selfWeight` with `-b 0 0 -gamma`.
2. **Rough footing** (`u_x = 0` on the footing nodes). R3's is smooth because the exact Prandtl
   answer is the smooth one; there is no exact answer here, so the deck uses the rough footing the
   campaign's problem actually has. **No number in this report may be compared to R3's.**
3. **`DROP_TRUNCATE` disabled.** R3 cuts the curve where `q` first falls 2 % below its max — a
   bottomed-out-run detector on a plateau problem, and a *signal deleter* on a softening one. It
   would delete exactly the branch this study exists to see.

**Staging** follows all three ADR-86 emitter-guide deck rules
(`86_ladruno_sanisand_apegmsh_emitter_guide` §2): gravity ramped over 10 steps at material stage
0; then `updateMaterialStage -material 1 -stage 1` explicitly for the tag; then the push.
`-Presidual 0.0` and `-Pmin 1e-3·P_atm` are the ADR-86 fork defaults, emitted explicitly.

**Two densities.** Gorini's `e_init = 0.6944`, and a denser `e_init = 0.60`. At `p ≈ 20` kPa the
state parameter is `psi = −0.123` vs `−0.217`, so `M^b/M_c` is 1.54 vs 2.14 and `M^d/M_c` is 0.49
vs 0.29. The dense leg dilates hard and softens hard: it is the ill-posed case, and the one that
matters.

**Push.** R3's pattern exactly: `sp` inside a load pattern under one-argument `LoadControl`,
pseudo-time *is* the imposed footing settlement, adaptive subdivision against the pinned budget.
Target `s/B = 0.25`, with an early stop the moment a peak is passed (§4).

## 3. Controls

RESULTS-PENDING-CONTROLS

## 4. Capacity: the three-clause rule, adapted for softening

`00_canonical_testbed` §1c requires three clauses, and clause 1 — *a flat tail* — is written for a
**perfectly plastic plateau**. A softening material does not plateau; it **peaks and descends**.
So clause 1 is satisfied here by *either* a plateau (|tail dq/ds| < 2 % of the initial tangent)
*or* a **passed peak**: `q` has fallen ≥ 5 % below `q_max` and stayed there for ≥ 15 further
steps. One dip is a solver artefact; a sustained descent is a post-peak branch. `PEAK` joins
`TARGET` and `BUDGET` as an admissible termination mode — the run stopped because the mechanics
finished, not because a guard bound. `FLOOR` and `WALL` remain **seizure** and are never a
capacity, whatever the tail says.

## 5. Two defects in the deck that were measuring the solver, not the soil

Both were found by instrumenting the run rather than by reading its curve, and both are on
`[[LEDGER_quirks]]`. They are recorded here because a reader who sees the wall times below needs
to know they are the times **after** these were fixed, and because the first one silently affects
every existing fork SANISAND deck.

### 5.1 `TanType` defaults to the ELASTIC tangent

`ManzariDafalias3D::getTangent()` returns `mCe` for `TanType == 0`, `mCep` for `1`, and
`mCep_Consistent` otherwise (`ManzariDafalias3D.cpp:134-142`). The **parser default is 0**
(`LadrunoSANISAND.cpp:117`, deliberately copied from `OPS_ManzariDafaliasMaterial` so a renamed
deck behaves identically), while the **null and parallel constructors default to 2**
(`ManzariDafalias.cpp:365`, `:426`). A deck that emits only the 18 positional parameters therefore
hands `algorithm Newton` an elastic tangent — modified Newton, with linear convergence, dressed as
full Newton.

It hides because **on a zero-free-DOF material-point deck it costs nothing**, and every existing
fork SANISAND deck is of that shape (`tests/test_ladruno_sanisand.py` passes `*_PARAMS` plus flags
and no positional optionals). On a boundary-value problem it was the difference between a
measurement and a non-measurement: at the parser default the h0 = 0.25 leg advanced 11 steps to
`s/B = 1.6e-4` in 350 s. This deck now emits the full positional block with `TanType 2`.
`mCep_Consistent` is unsymmetric under non-associated flow, which is why the solver is
`Pardiso -matrixType 0`.

### 5.2 The convergence test could not be met

The WP brief asks for `NormDispIncr` "per the R3 gate". **The R3 gate does not use
`NormDispIncr`** — it uses `NormUnbalance` at `1e-5 × the applied load`. And on SANISAND the
displacement form is not merely tight but *unreachable*: the substepped `ModifiedEuler` return
(hardcoded substep tolerance `TolE = 1e-4`) makes the discrete stress–strain map only piecewise
smooth, so Newton stalls rather than converging quadratically. Measured on the h0 = 0.5 leg, over
47 failed attempts, the residual displacement norm stalls at a **median of 6.6e-7 m** (min 3.4e-8,
max 1.3e-4). The run that nominally used `1e-8` was carried by the relaxed third ladder rung on
**18 of its 26 steps** — 65 of every 125 state-determination passes spent failing two rungs that
could not succeed.

A displacement-norm tolerance is also **not mesh-neutral**, which is disqualifying inside a
mesh-convergence study: the norm runs over the free DOFs, and there are 3.6× more of them at
h0 = 0.25 than at h0 = 1.0, so the same nominal number is a different physical requirement on each
mesh of the sequence. The push therefore uses **`NormUnbalance` at a tolerance relative to the
model's own weight `gamma·V`**, identical on all three meshes by construction.

### 5.3 The A/B that pinned the settings

h0 = 0.5, 400 s of wall each, same box, same `DS_MAX`:

| test / tolerance | reach s/B | steps | subdivisions | relaxed |
|---|---|---|---|---|
| `NormDispIncr` 1e-8 m (TanType 0) | 0.00106 | 26 | 1 | 18/26 |
| `NormUnbalance` 1e-6 ·γV (TanType 2) | 0.00218 | 31 | 0 | 11/31 |
| `NormUnbalance` 1e-5 ·γV (TanType 2) | 0.00442 | 37 | 0 | — |

and **the answer moves by a median of 0.3–0.75 %** between all three at matched settlement. The
max excursions (2.5–8 %) are interpolation artefacts of comparing two adaptive step sequences on a
curve still rising at ~25 kPa per 0.001 of s/B — an s/B mismatch of 4e-5 is worth 3 % there — and
they vanish at the peak, where `dq/ds → 0`. An independent confirmation at a *matched checkpoint*
rather than an interpolated point: `q` at `s/B = 0.002`, h0 = 1.0, reads **56.045** kPa under
`TanType 0` + `NormDispIncr 1e-8` and **55.756** kPa under `TanType 2` + `NormUnbalance 1e-5·γV` —
**0.52 %** apart across a complete change of solver configuration.

`1e-5 ·γV` is pinned because it is the only setting under which the step reaches `DS_MAX` at all,
and *reach* is the binding constraint on whether this study can answer its own question. One leg
is re-run at `1e-6 ·γV` as the tolerance control.

RESULTS-PENDING-TOLREF

## 6. Results

RESULTS-PENDING-TABLES

## 7. What this does and does not license ADR 90 to say

RESULTS-PENDING-DISCUSSION

## 8. Reproducing

```bash
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_summary.py --selfcheck
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py \
    --out <dir> --legs h1.0_e0.6944,h1.0_e0.60 --wall 5400
python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_summary.py <dir>
```

The three meshes were run as three concurrent processes (`OMP_NUM_THREADS=4` each, 24-core box),
one per resolution, coarse first so a partial answer existed early. **Wall times below are
therefore upper bounds under mild contention and are not a performance measurement.**
