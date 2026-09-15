---
title: "ADR-92 F10 — the IMPL-EX self-weight wall: `-implexControl` is what stops this deck, and the wall is the control's own primed test"
project: Ladruno
type: measurement note
status: "MEASURED 2026-09-14; REWRITTEN after adversarial review round 1 (leg N added, verdict re-ranked)"
priority: high
owner: nmora
amends: 92_ladruno_sanisand_implex_adr
requested_by: "TIMs Workbench, F10 (self-weight strip footing act)"
engine: "ladruno tip 9c2f964eae3bbd1a055c3ede81381a6c601a982b (pinned release build, nothing built for this WP; PREDATES PR #838)"
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[93_ladruno_sanisand_zero_confinement_adr]]"
  - "[[95_prandtl_reissner_campaign_report]]"
  - "[[_adr90_tau0_qu_band]]"
  - "[[LEDGER_quirks]]"
tags: [adr, material, sanisand, implex, soil, stepping, tims, measurement]
aliases: [ADR-92 F10, "implex self-weight wall"]
updated: 2026-09-14
---

# ADR-92 F10 — the IMPL-EX self-weight wall

> [!info] Numbering
> **`92b`, not a new ADR number.** This amends ADR 92; the folder's ADR high-water
> mark is **97** (`97_ladruno_asdp_closest_point_adr.md`) and is untouched. The
> `<n><letter>` suffix follows `40a`/`40b`, `75a`/`75b`/`75c`, `80a`–`80d`, `49a`.
> `92b` rather than `92a` because `_adr92_p1_*` / `_adr92_p2_9_*` already occupy
> the underscore-prefixed working-note space for P1/P2.

> [!warning] Review round 1 changed the verdict, not the craft
> The first cut of this note ranked the *stepping controller* as the cause and
> recommended pinning the growth factor. That was measured on a leg set in which
> **every arm had `-implexControl` on**, because the driver hard-wired it onto
> `-implex`. Adding the missing arm (**leg N**: bare `-implex`, doubling
> controller, everything else identical) reaches the target settlement in 58 s
> with **zero refusals and zero subdivisions**. The controller is a *co-factor of
> the control*, not the primary cause, and the first question is whether this
> deck needs the control at all. §§0, 3–9 are rewritten accordingly.

**A DIAGNOSIS. No C++, no banner, no flag default changed.** Test bed:
`Ladruno_files/testbed/hypo_bearing/adr92_f10/` (driver, reducer, artefacts,
README with the deck table and the declared differences from the TIMs deck).

---

## 0. The report, and the answer in three sentences

TIMs measured, on a plane-strain strip footing on SELF-WEIGHT SANISAND under
`-implex -implexControl`, that the control's refusal counter climbs from the
first push step (204 120 at step 3 on 9 720 Gauss points), that every refusal
cuts the step, and that the harness step falls below its floor at `s/B = 0.0125`
— while the same material class and build carried the weightless ADR-95 slab
campaign to `s/B = 0.15`. Minimum `p'` in the block is 3.8 kPa, so this is not
the ADR-93 apex wall. **Reproduced** on the fork's own deck: leg B walls on the
step FLOOR at `s/B = 0.0085` with 724 refusals, every one in the `control`
bucket.

1. **Removing `-implexControl` makes this deck run — and the ADR-95 reference
   campaign never used it either.** Bare `-implex` with the *same doubling controller* (leg N)
   reaches the target `s/B = 0.0500` in **104 steps, 0 subdivisions, 0 failed
   attempts, 58 s** — and the refusal ledger reads **0 / 0 / 0 / 0**, the
   **companion bucket explicitly zero**, which is the number that matters on
   this build (§3). That is a TERMINATION result, not an accuracy one: ADR-92 §8
   and the engine's own control-off echo (`:2029`-`:2031`) make a separate
   accuracy claim that this deck sits inside rather than outside (§3, "what leg N
   does not establish"). The ADR-95 leg D1 that reached `s/B = 0.15`
   (`95_prandtl_reissner_campaign_report.md:176`, 0 failed, 296 s) was itself
   control-OFF under a doubling controller
   (`sanisand_path_diag.py:91`-`:104` passes only `-implex`;
   `deformed_snapshot.py:294` is `ds = min(2*ds, dmax)`).
2. **With the control on, the doubling controller is a real co-factor — but only
   then.** Control on: growth ×2 walls at 0.0085 with 724 refusals; growth ×1.0
   reaches 0.0265 with 6. Control **off**: growth ×2 and growth ×1.0 (leg N1)
   both reach the target and agree to **0.383 % mean / 0.779 % max**, with ×2
   **6× faster** (58 s vs 364 s). The growth rule is harmless until an absolute
   per-step error bound is placed on a quantity that grows with the step.
3. **The FLOOR seizure itself is the control's own `implexPrimed` test.** Leg B's
   throttled warning lines carry refusals whose error does **not** decay with
   `dt` — 0.2243 at `|dt| = 4e-5` and 0.2143 at `2e-5` (step halved, error moved
   4.5 %) at `|d_eps_p(n)| = 9.09e-13`, and 0.063 at `f = 0` with
   `|d_eps_p(n)| = 9.58e-21`. That is exactly the signature
   `LadrunoSANISAND.cpp:2889`-`:2898` documents as the reason the un-primed step
   is exempt ("it ASYMPTOTED at ~0.076 instead of decaying … a companion jump
   that is independent of the increment … a dead analysis"). `implexPrimed` at
   `:2905` is a bare `> 0.0` sign test, so a Gauss point carrying a
   numerically-zero plastic history is "primed", loses the exemption, and is
   refused on a drift no subdivision can shrink.

---

## 1. The deck

`B = 1.5` m strip, `15B x 12B` box, plane strain, 285 `LadrunoBrick
-formulation bbar`, **2 280 Gauss points**; `gamma' = 9.81` kN/m3 as a body
force, 7.65 kPa surcharge outside the footprint, 18.4 kN/m of footing dead load;
`LadrunoSANISAND` on the ADR-92 CP1 / ADR-86 §5 Gorini set
(`tests/test_ladruno_sanisand.py::_PARAMS`, `e_init = 0.6944`), `IntScheme 1
TanType 2`, `-Presidual 0 -Pmin 0.0101 -maxSubsteps 1000`; rough rigid footing;
`confine -> updateMaterialStage 1 -> push`; `LoadControl(-ds)` on an `sp`
pattern under a `Linear` series with the `Transformation` handler (the guide's
push-idiom precondition — **never `DisplacementControl`**); R3's adaptive
controller with `DS_BASE = 2e-5`, `DS_MIN = 2e-7`, `DS_MAX = 1e-3`,
`GROW_AFTER = 6`, `SUBDIV_BUDGET = 80`.

Controls, per leg: gravity resultant identity, the **1-D geostatic patch**
(`3.1e-13` measured — the check the resultant is structurally blind to), and
`eta_max/M_c = 0.7275 < 1` at the stage flip.

**The K0 device.** `nu* = K0/(1+K0)` is passed as the material's Poisson ratio,
so the stage-0 elastic gravity state is exactly `K0` at every depth (exact
because `nu` is constant, even though SANISAND's moduli go as `sqrt(p)`).
**Declared difference:** the fork's `LadrunoSANISAND` takes `nu` positionally and
offers no way to put it back after the K0 stage, so `nu*` is held for the whole
leg rather than being temporary. Remaining declared differences (calibration,
mesh density, footing kinematics) are in the test bed's README; **measured `psi`
is `-0.128 … -0.109`** against TIMs' quoted `~ -0.13`.

**The at-rest census matches the report.** At `K0 = 0.455`, **2 280/2 280
= 100.00 %** of Gauss points sit on the DILATANT side of `M^d` (median
`eta = 0.857`, median `M^d = 0.666`) — TIMs' 99.96 %. `p'` at rest spans
**6.4 – 106 kPa**, and the minimum over the committed fixed-`ds` censuses is
**6.374 kPa** (`out/f10_census_B_ds4e-05.csv`), falling to 6.06 kPa at
`ds = 2e-4`. **That minimum is 1.27× the P0 low-confinement corner's 5 kPa**, and
§3 is where that matters.

---

## 2. Which criterion fires

`-implexControl` reaches `noteRefusalControl()` from three sites in
`SRC/material/nD/LadrunoSANISAND.cpp`:

| line | test | counted in |
|---|---|---|
| `:2905` | `implexPrimed = GetNorm_Cov(mImplexDEpsP) > 0.0` — the exemption test | (gates everything below) |
| `:2907` | `implexPrimed && mImplexError > errorTol` | (gates the three below) |
| `:2919` → `:3019` | `... && dtAbs >= reductionLimit * mImplexDt0` — **above** the reduction floor | `implexRefusals[2]` (`control`) |
| `:2919` else → `:3154` | **at** the floor, `-implexFloor implicit` (the default) delivers the companion and does **not** refuse | `implexGuards[0]` |
| `:2919` else → `:3123` | at the floor under `-implexFloor refuse` — **unreachable at the default floor mode**, listed for completeness | `implexRefusals[2]` (`control`) |
| `:2716` / `:3235` | the companion hit `-maxSubsteps` | `implexRefusals[3]` (`companion`) |

Measured across all seventeen legs:

* **every refusal is the tolerance branch above the floor.** Leg B: 724 total,
  724 `control`, **0 `companion`, 0 `signChange`**;
* **the reduction limit cannot bind on this deck, and that is arithmetic, not a
  measurement.** `reductionLimit * |dt0| = 0.01 x 2e-5 = 2e-7 m`, which is
  **exactly the harness's own `DS_MIN`** (the driver's line 150). The controller
  declares FLOOR at the same step size at which the material's floor branch
  would first become reachable, so `implexGuards[0] = 0` on every leg except
  **F2**, where `reductionLimit` was raised to `0.5` — it then fired **151**
  times and bought +71 % reach. §9 recommendation 3 says what that setting
  actually is;
* **the companion never failed** at `-maxSubsteps 1000`: `implexRefusals[3] = 0`
  on B, C, D, E, K, L, N and N1, and ≤ 42 anywhere (M).

The throttled warning line (`:3040`-`:3046`) names both quantities. The full set
for the five legs that produced them is committed as
`out/refusal_warnings_<leg>.txt` (49 lines for leg B, 11.6 kB) — the first cut of
this note gitignored them with the raw logs, and they are the evidence for §4.

---

## 3. Leg N — the arm the first cut was missing

`build()` used to append `-implexControl` unconditionally whenever `-implex` was
on, so **no leg could run IMPL-EX the way the fork's own reference campaign runs
it.** `-implexControl` is now a leg knob.

| | leg B (control on) | **leg N (control OFF)** | leg N1 (control OFF, growth 1.0) |
|---|---|---|---|
| `s/B` | 0.0085 | **0.0500 (TARGET)** | 0.0500 (TARGET) |
| steps / subdivisions / failed | 364 / 53 / 169 | **104 / 0 / 0** | 3 750 / 0 / 0 |
| refusals total / ctl / **comp** / sign | 724 / 724 / **0** / 0 | **0 / 0 / 0 / 0** | 0 / 0 / 0 / 0 |
| wall | 142 s (to `FLOOR`) | **58 s** | 364 s |
| `q_end` | 195.90 kPa @ 0.0085 | 924.91 kPa @ 0.05 | 920.55 kPa @ 0.05 |

`N1` vs `N` over `0.002 <= s/B <= 0.05`: **0.383 % mean, 0.779 % max**. So with
the control off the doubling controller is not merely survivable, it is
*accurate* and 6× cheaper.

### What leg N does and does not establish

**What it establishes — TERMINATION, by measurement.** The guide §3 concern is
that the commit-time companion must be able to *fail* rather than force-accept at
`dT_min`, because it runs where no global Newton is left to react; that concern is
answered by `-maxSubsteps` (this deck sets 1000) and by reading the result.
Measured: `implexRefusals[3]` (the companion bucket) is **0 on B, C, D, E, K, L,
N and N1**, and `<= 42` anywhere in the campaign (M 42, F1 29, I 12, H 6, F3 3,
J 2). So on *this deck* the companion integrated every increment it was handed,
control-off terminates cleanly, and it is **6× cheaper** than the same arm with
the growth pinned (58 s vs 364 s) and reaches 5.9× further than the controlled
arm. The discipline that makes that safe to repeat:

> **Run control-off, and read `implexRefusals[3]` at the end of every leg.** The
> pinned build `9c2f964` predates PR #838 (`wp/99-refusal-propagation-audit`,
> still OPEN at the time of writing), so a capped companion commit is otherwise
> **silent**; once #838 lands the run aborts instead and the read becomes a
> belt-and-braces check rather than the only one.

**What it does NOT establish — ACCURACY.** `-implexControl` is not only a
termination device; ADR-92 §8 and the engine's own constructor echo make an
accuracy claim about the *extrapolation*, printed on every control-off run
including leg N's (`LadrunoSANISAND.cpp:2029`-`:2031`):

> `-implexControl` OFF … P0 measured IMPL-EX unusable from `d_eps = 5e-4` at
> `p0 = 5 kPa`, so at a low-confinement corner the control is a requirement, not
> an option

**Nothing measured here touches that claim, and this deck sits inside its
range**, not decades outside it:

* the deck's **minimum `p'` is 6.374 kPa** (`out/f10_census_B_ds4e-05.csv`) —
  **1.27× the P0 corner's 5 kPa**, not a comfortable margin. (`-Pmin = 0.0101` kPa
  is the wrong yardstick and the round-1 text used it; the corner is the yardstick);
* leg N's strain increment **crosses `5e-4` at step 24** (`s/B = 0.0012`, where
  its `ds` first reaches 0.32 mm — verified in `out/f10_N.csv`) and, per the
  re-verification, runs at **1.28e-3 to 2.0e-3, i.e. 2.6–4× the P0 corner**, from
  `s/B = 0.00248` to the target, where `ds` is pinned at the 1.0 mm cap;
* and there is **no implicit anchor beyond `s/B = 0.00227`** (leg G), so nothing
  on this deck checks the extrapolation over the range where it exceeds the
  corner.

So the honest statement is narrow: **control-off is the configuration that
terminates and is far cheaper on this deck, and its accuracy over most of its
range is unanchored.** A deck whose companion *does* cap is a different case, and
the control (or #838) is how you find that out — but a control-off curve that
nobody has confirmed against an implicit twin is a reach result, not an accurate
one, exactly as ADR-92 §8 already requires.

---

## 4. The FLOOR seizure — a refusal that does not decay with `dt`

Leg B's 49 throttled warning lines, reduced (`out/refusal_warnings_B.txt`):

| `|d_eps_p(n)|` | `|dt|` | `f` | `implexError` |
|---|---|---|---|
| 2.40e-21 | 4e-5 | 0.5 | 0.0516 |
| 9.58e-21 | 4e-5 | 0 | 0.0631 |
| 9.06e-13 | 4e-5 | 2 | 0.1413 |
| **9.09e-13** | **4e-5** | 2 | **0.2243** |
| **9.09e-13** | **2e-5** | 1 | **0.2143** |
| 6.66e-12 | 2e-5 | 1 | 0.2111 |
| 2.86e-07 | 4e-5 | 0 | 0.2987 |
| 2.86e-07 | 2e-5 | 0 | 0.1491 |
| 2.86e-07 | 1e-5 | 0 | 0.0835 |
| **2.86e-07** | **5e-6** | 0 | **0.1724** |
| … 39 more rows, all at `|d_eps_p(n)| >= 5e-8` | 4e-5 / 8e-5 | 0/1/2 | 0.050 – 0.101 |

Two families, and they behave completely differently:

* **the ordinary family** (`|d_eps_p(n)| >= 5e-8`, 39 of 49 lines) refuses at
  `|dt|` 4e-5 and 8e-5 with errors just over `tol`, and its error decays with the
  step — the 2.86e-07 point runs 0.299 / 0.149 / 0.084 as `|dt|` goes
  4e-5 → 2e-5 → 1e-5, i.e. first order. Halving works on these;
* **the seizure family** carries a plastic history of **1e-12 to 1e-21** — i.e.
  numerically zero — and its error does **not** decay: at `9.09e-13` the step
  halves 4e-5 → 2e-5 and the error moves from 0.2243 to 0.2143, **4.5 %**. Even
  the asymptote the source itself documents at `:2889`-`:2898` as the reason the
  un-primed step is exempt: *"A pure elastic predictor's error must scale with
  `d_eps`; this one does not, which is the signature of a companion jump that is
  independent of the increment."*
* **and one point that spans both families.** `|d_eps_p(n)| = 2.86e-07` decays
  first order over three rungs — 0.2987 / 0.1491 / 0.0835 at `|dt|` 4e-5 / 2e-5 /
  1e-5 — and then **rebounds to 0.1724 at 5e-6**, i.e. it is *non-monotone in
  `dt`*. It refuses at every one of those four rungs. That single point is the
  clearest thing in the table: the error is first order until it is not, and a
  controller that only knows how to halve has no way to tell which regime it is
  in. It is also why a leg can burn rungs all the way to `DS_MIN` while its
  subdivision budget still has room.

`implexPrimed` (`:2905`) is `GetNorm_Cov(mImplexDEpsP) > 0.0`. `9.09e-13` passes
it. So a Gauss point that took essentially no plastic strain in the previous
committed step, and then yields in this one, is treated as primed, loses the
exemption, and is refused on the very quantity the exemption exists to tolerate.
**The controller then halves to no purpose until it reaches `DS_MIN`: that is the
`FLOOR` mode of legs B, E, F1, I and M.** The subdivision budget is *not* spent
(53 of 80 on leg B) — the *floor* is reached, which is the signature of an error
that stopped responding to `dt`.

**Why the refinement probe of §5 cannot see this.** The probe reports
`n_over_005 = 0` at every step size, and its worst point (`ele 135 gp 5`,
`p' = 93.6` kPa) never refuses in any leg. The probe measures the *bulk* error
field at a clean state; the seizure is a handful of points with a
numerically-zero plastic history, reached only along a path that has already been
refusing. Two further caveats on that probe, recorded so it is not over-read:
its maximum-error point runs at **`f = 0` at every `ds`** (so it measures the
guarded elastic predictor's drift, not the extrapolation), and `f` itself scales
with `ds` through the probe, so `O(ds)` and `O(f)` are **not separated** by it.

### The bulk field IS first order in the step

Walk to `s/B = 0.00851` on a refusal-free constant-`ds` path (638 steps of 2e-5,
`q = 164.686` kPa, reproduced to the digit by all seven runs), then take **one**
step of size `ds`, from a separate process per `ds`:

| `ds` (m) | median err | max err | GPs > 0.05 | max err / previous |
|---|---|---|---|---|
| 8e-5 | 2.14e-05 | 1.126e-02 | 0 | — |
| 4e-5 | 5.74e-06 | 5.778e-03 | 0 | 1.95 |
| 2e-5 | 1.35e-06 | 3.039e-03 | 0 | 1.90 |
| 1e-5 | 1.90e-06 | 1.672e-03 | 0 | 1.82 |
| 5e-6 | 2.54e-06 | 9.902e-04 | 0 | 1.69 |
| 2e-6 | 2.95e-06 | 5.851e-04 | 0 | 1.69 |
| 5e-7 | 1.66e-05 | 4.309e-04 | 0 | 1.36 |

First order down to ~5e-6, then a `dt`-independent floor of ~4e-4 — three orders
below `tol`. **That is the field the control is bounding, and it is well behaved;
the seizure lives outside it.** Note also that at this settlement the maximum
error at a step four times the base is 0.011, i.e. 4.4× under `tol` — so leg B's
`FLOOR` at this very settlement is a property of its *path*, not of the state.

---

## 5. The error field, and what the census can and cannot support

Per-Gauss-point census at a fixed `ds`, control tolerance set so high nothing can
refuse. Counts are Gauss points with `implexError > 0.05`, out of 2 280, at
**push step 1 and push step 3** (step 1 is the un-primed step `:2907` exempts;
step 3 sits at `s/B ≈ 8e-5` for `ds = 4e-5`):

| deck | `p'` min/med/max kPa | dilatant at rest | `2e-5` s1/s3 | `4e-5` s1/s3 | `8e-5` s1/s3 | `2e-4` s1/s3 |
|---|---|---|---|---|---|---|
| **B** self-weight `K0 = 0.455` | 6.4 / 29.2 / 106 | **100 %** | 32 / 0 | **80 / 4** | 136 / 4 | 328 / 81 |
| **A** weightless, uniform 10 kPa | 6.3 / 6.5 / 7.0 | 100 % | 64 / 0 | 204 / 4 | 596 / 8 | 1764 / 199 |
| **C** self-weight `K0 = 0.818` | 8.9 / 40.4 / 146 | **0 %** | 40 / 0 | 80 / 0 | 144 / 0 | 340 / 20 |
| **H** self-weight + 100 kPa | 70 / 93 / 170 | 100 % | 0 / 0 | 0 / 0 | 0 / 0 | 8 / 0 |

**Read the "4 of 2 280" correctly.** It is *push step 3* at `ds = 4e-5`,
`s/B ≈ 8e-5`. The same step size at **step 1** has **80** over tolerance (32 over
0.1). Excluding step 1 is defensible — `:2907` exempts it by construction — but
it is load-bearing for anything said about candidate (1), so the step is named
every time the number appears.

**Where those four points are (corrected).** All four are one state replicated
four times by symmetry: elements 120 and 180, two Gauss points each, at
`x = ±1.00` m, `z = -0.25` m. The footing half-width is **0.75 m**, so they sit
`dx = 0.25` m **OUTSIDE** the footing edge, not under it (the first cut of this
note said `dx = 0`, which was wrong). `p' = 9.897` kPa, `eta/M^d = 1.101`,
`psi = -0.126`, error 0.0542 of which **0.0542 is deviatoric and 0.0006
volumetric**, and **`f = 0` on all four** — the P2-2 guard had already switched
the extrapolation off. They are evidence for the guard/control conflict (§10,
defect 2), not for the dilatant-extrapolation candidate.

**The census metric predicts nothing about reach, and this note does not pretend
it does.** Max error at `ds = 4e-5` step 3 against the leg's own refusal count:
A 0.050 / 1 545, B 0.054 / 724, C 0.022 / 1 727, H 0.006 / 2 754. The slope `c`
in `err ≈ c·ds` is the **same** for A and B (~1.3e3 m⁻¹) and A reaches 2.2×
further; H has the smallest error field of all four and the **largest** refusal
count. **The first cut's claim that "self-weight sets where the maximum
admissible step sits" is WITHDRAWN: untested, and the H row contradicts it.**
What the census does support is narrow, and is stated as such: the error field is
smaller at higher confinement and smaller in the contractant-at-rest state, and
it is monotone in `ds`.

---

## 6. The leg table

`h0 = 0.5` m, 2 280 Gauss points, wall budget 400 s/leg, target `s/B = 0.05`.
`TARGET`/`BUDGET` are admissible modes; `FLOOR`/`WALL` are seizure/stop.
`relaxed` = converged steps carried by the third ladder rung (a 10× looser test).

> [!warning] The `q` columns are NOT comparable across legs
> A leg that refused, halved and re-grew arrives at a given settlement on a
> **different strain path** and reads a **stiffer** curve: at `s/B = 0.0085`,
> leg B reads 195.90 kPa against leg N's 164.07 — **+19.4 %**. See §7.

| leg | control | change vs B | s/B | mode | steps | subdiv | relaxed | wall s | refusals tot / ctl / comp | floor fb | f=0 guard |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **N** | **OFF** | **`-implexControl` removed; growth still 2.0** | **0.0500** | **TARGET** | **104** | **0** | 0 | **58** | **0 / 0 / 0** | 0 | 2 624 |
| **N1** | OFF | N + growth 1.0 | 0.0500 | TARGET | 3 750 | 0 | 0 | 364 | 0 / 0 / 0 | 0 | 82 300 |
| A | on | weightless + uniform 10 kPa | 0.0189 | BUDGET | 656 | 81 | 2 | 381 | 1545 / 1545 / 0 | 0 | 123 055 |
| **B** | on | **the reported configuration**, `tol 0.05` | **0.0085** | **FLOOR** | 364 | 53 | 5 | 142 | **724 / 724 / 0** | **0** | 31 491 |
| C | on | `K0 = 0.818` via `nu* = 0.45` | 0.0209 | BUDGET | 511 | 81 | 0 | 188 | 1727 / 1727 / 0 | 0 | 13 484 |
| D | on | `-implexFactor controlIter` | 0.0089 | WALL | 102 | 13 | 6 | 408 | 725 / 725 / 0 | 0 | 5 289 |
| E | on | hold at the flip + 2e-6 first step, growth 1.5 | 0.0139 | FLOOR | 431 | 40 | 4 | 122 | 506 / 506 / 0 | 0 | 17 697 |
| F1 | on | `tol = 0.1` (the shipped default) | 0.0117 | FLOOR | 236 | 37 | 0 | 110 | 342 / 313 / 29 | 0 | 24 393 |
| F2 | on | `reductionLimit = 0.5` | 0.0145 | BUDGET | 615 | 81 | 5 | 201 | 1246 / 1246 / 0 | **151** | 45 016 |
| F3 | on | `tol = 0.5` | 0.0500 | TARGET | 114 | 2 | 0 | 169 | 18 / 15 / 3 | 0 | 2 991 |
| G | on | **implicit** (no `-implex`) | 0.0023 | **WALL** | 29 | **0** | **10** | 413 | 0 / 0 / 0 | 0 | 0 |
| H | on | + 100 kPa uniform surcharge | 0.0135 | BUDGET | 496 | 81 | 0 | 129 | 2754 / 2748 / 6 | 0 | 24 254 |
| I | on | `-implexGuard off` | 0.0157 | FLOOR | 522 | 78 | 5 | 166 | 1032 / 1020 / 12 | 0 | 0 |
| J | on | `-implexGuard off` **and** `-implexTrialGuard off` | 0.0185 | BUDGET | 499 | 81 | 0 | 183 | 12758 / 12756 / 2 | 0 | 0 |
| K | on | growth factor pinned at 1.0 | 0.0265 | WALL | 1 984 | 0 | 0 | 400 | 6 / 6 / 0 | 0 | 57 109 |
| L | on | constant `ds0 = 8e-5`, growth 1.0 | 0.0187 | WALL | 2 794 | 3 | 0 | 400 | 144 / 144 / 0 | 0 | 51 716 |
| M | on | growth factor 1.25 | 0.0113 | FLOOR | 371 | 25 | 2 | 169 | 253 / 211 / 42 | 0 | 21 281 |

**Leg G was never walled.** Its `WALL` is the wall *clock*: `nsub = 0`, no failed
ladder, and its step was **still doubling at termination** (0.02 → 0.32 mm). It
is simply slow, and **10 of its 29 converged steps ran on the relaxed rung**
(10× looser test) after burning two full-Newton rungs first — so its 413 s
carries a large ladder overhead, and its reach is not a statement about what the
implicit material can do.

**Reach per wall second, honestly:** leg N `8.6e-4 /s`, leg G `5.5e-6 /s` — a
factor of **~156**, not the "11×" the first cut quoted against leg K. The right
comparison is against the control-off arm, because that is the recommended
configuration.

---

## 7. Accuracy — what is and is not established

**Two arms reached the target, so one long-window comparison exists.** Over
`0.002 <= s/B <= 0.05`:

* `F3` (`tol 0.5`) vs `N` (control OFF): **0.395 % mean, 1.625 % max**;
  `q_end` 921.03 vs 924.91 kPa; 18 refusals vs 0. **`tol = 0.5` is very nearly an
  inert control on this deck** — it costs 111 s of extra wall time and buys a
  0.4 % difference;
* `N1` (control OFF, growth 1.0) vs `N`: **0.383 % mean, 0.779 % max**.

**The refuse/halve/regrow path STIFFENS the curve**, monotonically:

| s/B | leg B `q` | leg N `q` | B/N − 1 |
|---|---|---|---|
| 0.0030 | 64.22 | 63.14 | **+1.72 %** |
| 0.0050 | 103.30 | 100.72 | **+2.56 %** |
| 0.0070 | 145.53 | 137.31 | **+5.99 %** |
| 0.0085 | 195.90 | 164.07 | **+19.40 %** |

and at leg B's own terminal settlement the arms split by whether they refused:
B +19.4 %, L +20.2 %, I +18.5 %, F2 +12.1 % against N; F1 +2.4 %, M +1.5 %,
J +0.7 %, K +0.3 %, F3 −0.7 %. The refusal-free constant-`ds` walk of §4 reads
164.69 kPa at `s/B = 0.008507`, i.e. with N to 0.4 %. **The first cut's "every
arm agrees to 1–3 %" holds only on `0.00045 <= s/B <= 0.00227`** and is corrected
here.

**Not established.** Against the *implicit* material: the implicit leg reached
only `s/B = 0.00227`, so no `-implex` curve on this deck is confirmed against an
implicit twin over any useful window, and **no capacity, peak or plateau is
claimed on any leg**. Also not established: mesh convergence (one mesh,
`h0 = 0.5`), behaviour past `s/B = 0.05`, and parallel behaviour (ADR-92's
`sendSelf` limit stands). Leg E's load column is offset by a harness artefact
(its reference reaction is captured before its hold) — its reach is usable, its
`q` is not.

---

## 8. Verdict against the three candidates, re-ranked

| candidate | verdict | evidence |
|---|---|---|
| **(1)** the DILATANT-AT-REST extrapolation error, "wrong at every point at once" | **Aggravator of the error field, not the cause of the wall.** 100 % dilatant at rest confirmed (TIMs: 99.96 %), and the contractant twin (C) has a strictly smaller error field at every step size. But the wall is not an error-field phenomenon at all: the *same* error field with the control OFF (leg N) produces zero refusals and the target settlement. And at the step that refuses, **push step 3** has 4 over-tolerance points of 2 280 — one state replicated 4× by symmetry, **outside** the footing edge, all at `f = 0` (step 1 has 80, and is exempt). | §3, §5 |
| **(2)** the substepper's error control vs the control's tolerance | **PARTLY — and it is the CONTROL's error control, not the SUBSTEPPER's.** The substepper is exonerated on the arms that matter: `implexRefusals[3] = 0` on B, C, D, E, K, L, N and N1, and `<= 42` anywhere in the campaign (M 42, F1 29, I 12, H 6, F3 3, J 2) — so the companion is not what refuses leg B. But the control's own machinery **is** the wall: its `implexPrimed` gate (`:2905`) is a bare sign test, so points with a 1e-12…1e-21 plastic history lose the un-primed exemption and are refused on an error that does not decay with `dt` (0.2243 → 0.2143 for a halved step). That is the `FLOOR` of legs B, E, F1, I and M. **The first cut called this candidate REFUTED; that was wrong.** | §2, §4 |
| **(3)** the `nu*` K0 device | **REFUTED.** At `K0 = 0.455` the device's `nu* = 0.31271` **is** the material's own calibrated `nu = 0.3129` to three decimals, so leg B is simultaneously the "K0 reached by the material's own elasticity" control and there is nothing anomalous left behind. (The leg-C half of the first cut's argument is **withdrawn**: leg C changes `K0` *and* holds `nu = 0.45` for the whole push, so it is not a clean one-variable test of the device.) | §1 |
| **(4)** the stepping controller | **A CO-FACTOR OF THE CONTROL, not an independent cause.** Control on: growth ×2 → 724 refusals and `FLOOR` at 0.0085; ×1.25 → 253 and 0.0113; ×1.0 → 6 and 0.0265. Control **off**: ×2 and ×1.0 both reach the target and agree to 0.383 %, with ×2 six times faster. The growth rule is harmless until an absolute per-step bound is placed on an error that grows with the step. | §3, §6 |
| **(5) — the primary answer** | **`-implexControl` is what stops this deck.** Bare `-implex` reaches `s/B = 0.05` in 104 steps, 0 subdivisions, 0 failed attempts, 58 s, with the companion bucket verified at zero — which is how the ADR-95 campaign that reached 0.15 was run. Scope: that is a TERMINATION result; §8's accuracy claim is untested here and this deck sits inside its range (§3). | §3 |

**A caveat we cannot close.** TIMs report `K0 = 0.82` walling *earlier* (0.0082)
than `K0 = 0.455` (0.0125); this deck measures the opposite ordering by 2.5×
(C 0.0209 vs B 0.0085). The likely reason is the one difference we cannot
reproduce: TIMs' `nu*` is **temporary** — restored to the calibrated value after
the K0 stage — so their `K0 = 0.82` leg pushes on the calibrated `nu` from a
`K0 = 0.82` initial stress, while ours pushes on `nu = 0.45` throughout. The fork
takes `nu` positionally with no runtime setter, so this is not testable here. The
two orderings should not be quoted as agreement or disagreement.

---

## 9. Recommendation

**For a self-weight SANISAND push, in order.**

1. **Run bare `-implex`, read the companion bucket, and do not call the result a
   capacity.** Measured: target settlement, 104 steps, 0 subdivisions, 58 s;
   `implexRefusals[3] = 0`. Check that bucket at the end of every leg — on a build
   predating PR #838 a capped companion commit is otherwise silent; once #838
   lands the run aborts instead. This is the configuration the ADR-95 reference
   campaign used, and on this deck it is 156× the reach-per-wall-second of the
   implicit leg. **But see §3:** the control is also an accuracy device, this
   deck's minimum `p'` (6.374 kPa) is 1.27× the P0 corner and leg N's strain
   increment runs up to 4× the corner's `5e-4` with no implicit anchor past
   `s/B = 0.00227`. Control-off buys reach and termination; it does not buy a
   confirmed curve.
2. **If you want the control, either pin the growth factor at 1.0 or use
   `tol = 0.5`.** Growth 1.0 at `tol 0.05`: 0.0265, 6 refusals, 400 s. `tol 0.5`
   at growth ×2: the target, 18 refusals, 169 s, and within **0.395 % mean /
   1.625 % max** of the control-off arm — i.e. **`tol 0.5` makes the control
   nearly inert on this deck**, which is the honest description of what it buys.
   A gentler growth factor (×1.25) is **not** enough.
3. **`reductionLimit = 0.5` is not a tuning knob, it is "switch the control off
   after one halving."** It buys +71 % with `tol` untouched, but what it does is
   move the floor above the working step so the P2-1 floor branch — which
   *delivers the companion and does not refuse* — takes over almost immediately.
   Say that plainly rather than presenting it as a tolerance-preserving win. (At
   the shipped `0.01` the floor sits at `2e-7` m = the harness `DS_MIN` exactly,
   which is why it never fires.)
4. **`-implexGuard off` (+85 %) is a real lever and is NOT recommended**: the
   guard exists for ADR-93's softening seat, which this deck does not test.
5. **`-implexFactor controlIter` is not the tool here** (0.0089 vs 0.0085 at 2.9×
   the wall time), as the guide §12's own "does not generalise" caveat predicts.
6. **The implicit lane is the right answer when the curve must be a capacity**,
   and this campaign cannot say more than that: its implicit leg reached
   `s/B = 0.00227` in 413 s with 10 of 29 steps on the relaxed rung and its step
   still doubling — slow and ladder-heavy, but **not walled**.

---

## 10. Two defects, one promoted

1. **`implexPrimed` is a bare sign test — PROMOTED to the mechanism of the FLOOR
   seizure.** `LadrunoSANISAND.cpp:2905`:
   `const bool implexPrimed = (this->GetNorm_Cov(mImplexDEpsP) > 0.0);`. The
   exemption at `:2907` exists because the companion's drift-correction jump does
   not scale with `d_eps` (`:2889`-`:2898`), but any non-zero plastic history,
   however small, forfeits it. Measured on leg B at `|d_eps_p(n)|` of 9.09e-13,
   6.66e-12, 9.58e-21 and 2.40e-21, with errors 0.21–0.22 that move 4.5 % when
   the step halves. The matching fix is the shape P2-5b already used for the
   reversal reset — a **relative** test, `||d_eps_p(n)|| > rel * ||d_eps||` — and
   it needs its own gate and mutation score. **Not fixed in this WP.**
2. **`-implexGuard`'s `f = 0` and `-implexControl` work against each other.** The
   guard forces `f = 0` on a reversing/softening predecessor; `sigma~` is then a
   pure elastic predictor and the control refuses the step *for being*
   inaccurate — and P2-6's trial-time fallback cannot help, because it only runs
   when `mImplexFactor != 0.0` (`:2939`). **30 of leg B's 49 throttled refusal
   lines are at `f = 0`**, and so are all four of §5's over-tolerance census
   points. Either the guard should exempt the point it just guarded from the
   tolerance for that step (the same shape as the un-primed exemption), or the
   control should measure a guarded point against the elastic predictor's own
   expected error. **Design question for ADR 92, not a fix here.**

Both are in `LEDGER_quirks.md`.

---

## 11. Campaign hygiene

* **Wall times sum to ≈ 84 minutes** across 17 legs, 16 censuses and 7 probes,
  not the 75 the test bed's README first said. Four legs (D, G, K, L) terminated
  on the wall clock, and **the box was not attested idle** — treat every wall
  number as an upper bound on a busy desktop, and do not compare wall times
  across legs that ran at different times.
* **Legs L and M were re-run.** Their first run shared `out/f10_L.csv` with a
  second, accidentally-launched process: the CSVs carried interleaved rows and a
  torn line, exactly the failure `hypo_bearing/README.md` already documents for
  the ADR-79 runner. The driver now carries that runner's guard (refuse a CSV
  another process touched in the last 180 s; `F10_FORCE=1` overrides) and the
  reducer drops torn lines. The re-run changed leg M not at all (`0.011286`,
  identical to the digit — the physics was deterministic, only the file was
  damaged) and leg L from 0.0185 to 0.0187; both wall times rose once they had
  the box to themselves (M **109.3 → 169.4 s**, the 109.3 read from the
  pre-re-run `f10_M.json` at `26d5c607f`, not from the interleaved batch log).
