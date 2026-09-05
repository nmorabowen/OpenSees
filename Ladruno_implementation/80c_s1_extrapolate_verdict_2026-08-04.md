---
title: ADR-80 S1 — LadrunoLoadControl -extrapolate is BUILT, SAFE, and DOES NOT MEET ITS ACCEPTANCE GATE
project: Ladruno
status: complete (negative result)
priority: high
owner: nmora
tags:
  - implementation
  - solver
  - integrator
  - constraints
  - validation-gate
aliases:
  - ADR-80 S1 verdict
updated: 2026-08-04
parent: "[[80_ladruno_sp_imposition_strengthening_adr]]"
---

# ADR-80 S1 verdict — the predictor works, and it is not enough

> ADR-80's acceptance gate says, in its own words: *"Anything less than a large
> improvement means the mechanism is wrong or incomplete — **report that, do not
> tune toward the preferred conclusion**."* This is that report.

## What was built

`LadrunoLoadControl` (classTag **33015**), a fork-owned `StaticIntegrator` that
replicates stock `LoadControl` and adds `-extrapolate <frac>`, the Abaqus
`EXTRAPOLATION=LINEAR` analogue. Predictor ordering copied from
`DisplacementControl::newStep` (`incrDisp → applyLoadDomain → updateDomain`).

**Not** a subclass of `LoadControl` — its members are private, and
`LadrunoArcLength` sets the fork precedent of deriving from `StaticIntegrator`
directly. **Vanilla footprint is therefore ZERO** apart from the unavoidable
registration points (`classTags.h`, the two interpreter ladders, CMake).

Plus `ladrunoLoadControl` — a runtime command (`setDeltaLambda`, `deltaLambda`,
`extrapolate`, `armed`), for the reason in §2.

## 1. The safety property HOLDS

`-extrapolate 0.0` is **bit-identical to stock `LoadControl`** in every idiom
measured — 43 / 23 / 224 on the adaptive march, 60 iterations on the fixed one,
matching stock exactly. The default is off; existing decks cannot be affected.

## 2. A design defect found by measurement: re-issuing the integrator kills the predictor

`integrator LadrunoLoadControl ...` **constructs a new object**. A caller-driven
adaptive ladder that re-issues it every step — which is exactly what the G5
harness and the fork's robust-solve driver do — therefore destroys the stored
increment every step and **the predictor never fires at all**:

| idiom | `frac` | increments | cutbacks | iterations | predictor armed |
|---|---|---|---|---|---|
| re-issue | 0.0 | 43 | 23 | 224 | 0 |
| re-issue | **1.0** | 43 | 23 | **224** | **0** |

Byte-identical to stock. ADR-80 §6 asked "does the predictor interact badly with
a caller-driven adaptive march?" — the answer is worse than *badly*: under the
fork's own driver shape it is **inert**.

**Fix shipped:** the `ladrunoLoadControl setDeltaLambda` runtime command, so an
adaptive caller resizes the step without reconstructing the object. Modelled on
`ladrunoArcLength`'s existing Layer-B reduceStep command. With it the predictor
demonstrably fires — `armed` reports **63** armed steps.

⚠ Registered in `PythonWrapper.cpp` / `TclWrapper.cpp` but **not** in
`SRC/tcl/commands.cpp`, so `OpenSees.exe` does not have it. That gap is
**pre-existing and not specific to this work** — `ladrunoArcLength` is missing
from the `.exe` for the same reason. Idiom C is therefore measured from Python.

## 3. The result, with the predictor genuinely firing

| idiom | `frac` | armed | increments | cutbacks | iterations |
|---|---|---|---|---|---|
| runtime-cmd | 0.0 | 0 | 43 | 23 | 224 |
| runtime-cmd | **1.0** | **63** | 42 | **23** | **201** |
| fixed march | 0.0 | — | 10 | 0 | 60 |
| fixed march | **1.0** | — | 10 | 0 | **35** |

- **Fixed march:** excess over the elastic control falls 40 → 20, a **50 %**
  reduction. Real.
- **Adaptive march with cutbacks — the gate that matters:** iterations fall
  224 → 201 (**−10 %**) and **cutbacks do not move at all: 23 → 23.**

The converged answer is 0.150000 everywhere; nothing about accuracy changes.

## 4. Verdict: the acceptance gate FAILS

The stated criterion was that `-extrapolate 1.0` should drive
**43 / 23 / 224 toward 6 / 0 / 12**, with **cutbacks going to zero as the
headline number**. Measured: **42 / 23 / 201**. Cutbacks are untouched.

This is not a wiring failure — `armed` = 63 proves the predictor ran on 63
steps. It is a **capability** failure: linear extrapolation of the previous
increment *reduces* the driven layer's overstrain but does not *remove* it, and
what survives is still enough to trip the spurious yield that drives the
cutbacks.

### This also partly walks back my own G3 correction

[[80b_sp_gates_g1g2g3_2026-08-04]] argued that because the collapsed tangent is
downstream of the overstrain, a predictor removes **both** halves and its ceiling
is the full excess. That argument is sound for a predictor that *eliminates* the
overstrain. **A linear predictor does not eliminate it**, so the argument
describes an unreachable ideal rather than this remedy. The reasoning was right
about what bounds the ceiling; it was wrong to imply `-extrapolate` would reach
it. Both statements now stand corrected in the record rather than quietly
re-scoped.

## 5. Consequence: ADR-80 D6 is triggered

D6 says `TransformationFE::getTangForce()` — the Kratos-style `b −= K·Δu_D`
route — is *"deferred unless S1's acceptance gate fails."* **It has failed.**
That route is principled and exact: it forms the tangent at the **committed**
state and multiplies by the prescribed increment, so the prescribed motion
reaches the RHS with **no constitutive evaluation in the lagging state at all** —
the overstrain is eliminated, not merely reduced. Its cost is the blast radius
D6 records (behaviour changes for every model with a non-homogeneous `sp`).

**Recommendation:** keep `LadrunoLoadControl` — it is safe, off by default, and
buys a genuine 50 % on fixed marches — but do **not** present it as the fix for
the ADR-80 mechanism, and open D6/Candidate C as the next step rather than
tuning `frac`.

> ✅ **FOLLOWED, AND IT WORKED** — [[80d_p3_tangent_predictor_verdict_2026-08-04]].
> Candidate C shipped as `-tangentPredictor` on this same class: cutbacks
> **23 → 0**, iterations **224 → 12 = the elastic control exactly**. This
> section's diagnosis is confirmed in full: a *linear* predictor reduces the
> overstrain, a *tangent* predictor eliminates it, and only the second one
> reaches the ceiling §3 described.
>
> Two numbers in the table above are also re-measured there and reproduce
> exactly, including `-extrapolate 1.0` being **completely inert (60, not 35)**
> under the re-issue idiom.

## 6. What NOT to do

- **Do not tune `frac` upward** hunting for the gate. `frac` = 0.5 gave 56
  iterations against `frac` = 1.0's 35 on the fixed march; the trend is
  monotone and 1.0 is already the Abaqus default. There is no setting that
  turns −10 % into cutbacks-to-zero.
- **Do not quote the fixed-march 50 %** as S1's headline. The field failure mode
  is cutbacks, and on cutbacks this measures **zero**.
