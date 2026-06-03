---
title: "LadrunoRebarBuckling — Rebar-Buckling Wrapper UniaxialMaterial"
project: Ladruno
type: reference / user guide
status: shipped (DM + GA monotonic + cyclic re-straightening)
classTag: 33001
material: uniaxialMaterial LadrunoRebarBuckling
related:
  - "[[14_ladruno_rebar_buckling_adr]]"
  - "[[LadrunoUniaxialJ2_guide]]"
tags:
  - material
  - uniaxial
  - buckling
  - steel
  - rebar
---

# LadrunoRebarBuckling — Rebar-Buckling Wrapper UniaxialMaterial

`uniaxialMaterial LadrunoRebarBuckling` (classTag **33001**) is a
**stress-modifying wrapper** around *any* tension–compression
`UniaxialMaterial`. In compression past a slenderness-dependent onset it applies
a reinforcing-bar **buckling-average** degradation — the section-average
compressive stress of a bar buckling between ties drops well below the bare-bar
stress — while leaving the wrapped material byte-untouched. It is the clean,
opt-in, geometric overlay companion to the bare steel material (designed for
[[LadrunoUniaxialJ2_guide|LadrunoUniaxialJ2]], composes over `Steel02`/`Steel4`).

See the ADR [[14_ladruno_rebar_buckling_adr]] for the design record and the full
test matrix.

## Contents

1. Intended use — and what it is *not*
2. The wrapper idea
3. Theory — Dhakal–Maekawa (`-model dm`)
4. Theory — Gomes–Appleton (`-model ga`)
5. Cyclic re-straightening (v2)
6. Consistent tangent
7. Command syntax and examples
8. Recording and parameters
9. Composition and stacking
10. Verification & validation
11. Limitations and boundaries
12. References

---

## 1. Intended use — and what it is *not*

**Use it for** the average compressive stress–strain degradation of a single
reinforcing-bar fiber buckling between restraints (ties), with cyclic
re-straightening — i.e. rebar / RC fiber sections where bar buckling matters.

**It is NOT:**
- **Steel-profile local (plate) buckling** — a 2-D plate phenomenon with
  post-buckling reserve and boundary-condition dependence; use element-level IMK
  ([[project_ladruno_imk_beam|LadrunoIMKBeam]] / `Bilin`) or a shell/continuum
  route, not this fiber wrapper.
- **Global member buckling** (`KL/r`) — that is geometry, the corotational
  element's job (e.g. `LadrunoBrick -geom corot` or corotational beam transforms
  with an imperfection). **Do not double-count** member buckling with this fiber.

## 2. The wrapper idea

Dhakal–Maekawa computes the buckled *average* stress as a modification of the
bare-bar curve at the same strain — the bar's material still follows its own
σ–ε; the section average drops because the bar bows. So the natural decomposition
is **core = the bar material (a pure oracle), wrapper = the geometric buckling
average on top**:

```
σ_buckled = f(σ_bare, ε, λ)            λ = s/d   (tie spacing / bar diameter)
```

Benefits: the wrapped material stays bit-for-bit unchanged (its own verification
is preserved); the layer is **stackable and opt-in** (`Fatigue ∘ RebarBuckling ∘
J2`); it is **reusable** over any bare-bar material; and `-lsr 0` is an exact
identity gate (zero overhead, pure pass-through).

The formulae are ported from `ReinforcingSteel::Buckled_stress_Dhakal` /
`_Gomes`, but applied in **engineering strain on a generic wrapped bar** (where
`ReinforcingSteel` works in log strain on its own monolithic backbone), so the
match to `ReinforcingSteel` is a tight-tolerance formula port, **not** bit-for-bit.

## 3. Theory — Dhakal–Maekawa (`-model dm`, default)

`λ = s/d`; `ε_y = f_y/E`; `e = ε − e_cross`, where `e_cross = ε_max − σ_max/E` is
the tensile-unload anchor (the max tensile strain excursion). Define
`k = √(f_y/E · 2000)`.

**Onset (deep-branch) strain:**
```
ε* = −max(7, 55 − 2.3·k·λ)·ε_y          (more slender ⇒ shallower, earlier onset)
```

**Buckled average** (`f*L = σ_bare(ε*)`, the bare stress at the onset strain,
read by a committed-clone probe of the wrapped bar):
```
f* = f*L · α · (1.1 − 0.016·k·λ),   floored if f* > −0.2 f_y ⇒ f* = −0.2 f_y
```
- **Pre-onset / pass-through:** `e ≥ −ε_y` ⇒ `σ = σ_bare` (no change until the
  bar yields in compression).
- **Intermediate branch** `−ε_y > e ≥ ε*`:
  `σ = σ_bare·[1 − (1 − f*/f*L)·(e+ε_y)/(ε*+ε_y)]` — a multiplicative reduction
  `r ∈ [resid, 1]` that ramps from `r=1` at `e=−ε_y` down toward `f*/f*L` at `ε*`.
- **Deep branch** `e < ε*`: ramps toward the `−0.2 f_y` residual floor.

> ⚠️ The buckling reduction **begins at the yield strain `−ε_y`** (the
> intermediate branch), *not* at `ε*`. `ε*` only marks where the deep branch
> starts. So "no buckling effect" means the extreme strain is shallower than
> `ε_y`, not shallower than `ε*`.

`α` is the single DM residual-shape factor (`-alpha`, = `ReinforcingSteel`'s
`beta`, default 1.0).

## 4. Theory — Gomes–Appleton (`-model ga`)

Alternative rebar law (ported from `Buckled_stress_Gomes`). Gate: pass-through
unless `ε < e_cross` (buckles for *any* compression past the anchor, unlike DM's
`e < −ε_y`). With `fs_buck = √(32/(e_cross − ε)) / (3π·λ)`,
`factor = min(1, fs_buck)·β·(1 − r) + r`:
```
σ = f_sup·g − (factor + g)·(f_sup·g − σ_bare)/(1 + g)
```
where `f_sup = σ_max` (anchor stress), `g = fsu_fraction` (`-fsufrac`), `r =
reduction` (`-reduction`). `r = 1 ⇒ no buckling`; `r = 0 ⇒ full GA`. GA needs
**no `-fy`**. (Faithful quirk: upstream shadows its `beta`/`gama` with hardcoded
locals, so the GABuck `beta` argument is dead — reproduced.)

## 5. Cyclic re-straightening (v2)

When a buckled bar is unloaded and reloaded in tension it **straightens** before
it can carry full load — the buckling reduction recovers back to the bare curve
over a finite tension-side strain span. The wrapper models this with a
`{PASS, BUCKLING, RESTRAIGHTEN}` committed-latch state machine over the buckling
**raise** `g = σ − σ_bare ≥ 0` (buckling reduces the compressive magnitude, i.e.
raises σ toward 0):

```
σ = σ_bare + g ,        g ∈ [0, g_law] ,    g_law = σ_v1 − σ_bare
```

- **BUCKLING** = exactly the monotonic law above (`g = g_law`), so a monotonic
  history is bit-identical to the no-cyclic behaviour.
- On reversal the raise `g` is **held** on the compression side, then **decayed
  to 0** over the span `L_rs` with a smoothstep `B(q) = 1 − (3q²−2q³)`
  (C1-clean), tracking the live `σ_bare`, and rejoins the bare curve.
- A partial reload followed by re-compression re-enters BUCKLING
  **C0-continuously** (a strain-based smoothstep carry of the residual raise back
  onto the deepest-bow envelope) — no stress jump.

**The span `L_rs`** is set by `-restraighten`:
- `-restraighten c <value>` (default, `c=1.0`): `L_rs = max(c·|e_rev − e_cross_rs|, ε_y)`.
- `-restraighten lambda`: `L_rs = f(λ)·ε_y`, `f(λ) = 0.5 + 0.10·max(0, λ−5)`,
  clamped to `[0.5, 6]·ε_y`.

> Practical note: with the default `c=1` and a virgin bar (`e_cross_rs≈0`) a deep
> buckle gives a very long `L_rs` (≈ the whole compressive excursion). For most
> cyclic work `-restraighten lambda` is the better-scaled default. The exact
> shape/`f(λ)` are physically motivated but not yet fit to the Dhakal–Maekawa
> 2002 cyclic companion — both forms are exposed precisely so they can be
> calibrated.

## 6. Consistent tangent

The tangent is **analytic** (improving on `ReinforcingSteel`, which
finite-differences it):
- **DM monotonic:** `dσ/dε = r·k_bare + σ_bare·(∂r/∂ε)` (closed-form `∂r/∂ε` per
  branch; the deep-branch floor gives a flat residual ⇒ near-zero tangent).
- **GA monotonic:** hybrid `ρ·k_bare + (∂σ/∂factor)·factor'`, with `factor'` from
  a cheap closed-form central FD of the scalar `factor` (no material probe).
- **RESTRAIGHTEN Phase 2:** `dσ/dε = k_bare + h_rev·B'(q)/L_rs` (C1-clean at both
  seams).

## 7. Command syntax and examples

```tcl
uniaxialMaterial LadrunoRebarBuckling $tag $matTag -lsr $s_over_d \
    <-model dm|ga> <-alpha $a> <-reduction $r> <-fsufrac $g> \
    <-fy $fy> <-E $E> <-restraighten c $c | -restraighten lambda>
```
- `$matTag` — the wrapped bare-bar `uniaxialMaterial` (already defined).
- `-lsr` — slenderness `s/d` (the primary input); **`-lsr 0` ⇒ exact pass-through.**
- **DM requires `-fy`** (a generic wrapped material exposes no yield-stress
  getter); GA does not. `-E` should equal the wrapped bar's elastic modulus (it
  sets `ε_y` and the deep-branch slope); if omitted it is taken from the bar's
  initial tangent.

```tcl
# RC bar fiber: J2 steel -> rebar buckling -> (optional) fatigue
uniaxialMaterial LadrunoUniaxialJ2    10  $E -iso voce $fy $Qinf $b -kin 3 ...
uniaxialMaterial LadrunoRebarBuckling 11  10  -lsr 8.0 -model dm -fy $fy -E $E \
                                              -restraighten lambda
uniaxialMaterial Fatigue              12  11  -min -0.05 -max 0.05
# the fiber section consumes tag 12 (or 11 without fatigue)

# Gomes-Appleton, full knock-down:
uniaxialMaterial LadrunoRebarBuckling 21  20  -lsr 10.0 -model ga -reduction 0.0 -E $E
```

```python
# OpenSeesPy — wrapped bar fiber in a section
ops.uniaxialMaterial("Steel02", 1, fy, E, 0.01, 18.0, 0.925, 0.15)
ops.uniaxialMaterial("LadrunoRebarBuckling", 2, 1,
                     "-lsr", 8.0, "-fy", fy, "-E", E, "-restraighten", "lambda")
ops.fiber(y, z, A_bar, 2)        # inside a fiber section
```

## 8. Recording and parameters

**Recordable responses** (`-G`/`eleResponse … material … <resp>`):
`stress`, `strain`, `tangent`, `stressStrain`, and `reduction` (alias
`buckling`) — the effective reduction `σ/σ_bare`. Unrecognised responses are
forwarded to the wrapped bar, so core quantities (`plasticStrain`, `backStress`,
…) remain reachable through the wrapper.

**Settable parameters** (`setParameter`): `lsr`, `alpha`, `restraightenC`;
anything else is forwarded to the wrapped bar.

## 9. Composition and stacking

Because it is a wrapper, it composes: `Fatigue ∘ RebarBuckling ∘ {J2|Steel02}`.
The buckling reduction propagates through an outer `Fatigue` layer, and Fatigue's
`-min`/`-max` rupture triggers on the *buckled* response (verified, test B6).
Serialization (`sendSelf`/`recvSelf`) carries the wrapper state **and** the
nested wrapped material via the broker (test B7 / B4h).

## 10. Verification & validation

- **Material battery** (`tests/test_ladrunoRebarBuckling_material.py`, 29):
  identity gate (B0), pre-onset pass-through (B1), DM/GA formula-port oracles
  (B2/B2b/GA2), slenderness & reduction trends (B3/GA3), consistent-tangent FD
  (B5/GA5/B4e), composition (B6), serialization round-trip (B7/B4h), and the v2
  cyclic battery (B4a reduce-to-v1, B4b rejoin, B4c hold-raise, B4d Phase-2
  bracket, B4f C0 re-entry, B4-GA law-agnostic, B4g structural Newton, B4i `L_rs`
  floor).
- **Fiber-section integration** (`tests/test_ladrunoRebarBuckling_section.py`, 3):
  B4j `zeroLengthSection` moment–curvature (buckled peaks shed capacity vs the
  identity gate, smaller loop area); B4k `forceBeamColumn` RC cantilever cyclic;
  B4l multi-element `dispBeamColumn` RC column (hinge localizes at the base).

## 11. Limitations and boundaries

- **Softening ⇒ localization.** The post-onset negative tangent localizes; for
  fiber-section elements this is usually acceptable (it is the section average),
  but the same regularization caveat as other softening materials applies; use
  adaptive step-cutting (and IMPL-EX where available) for robustness.
- **Calibration deferred.** The cyclic span/shape (`L_rs`, smoothstep) are
  physically scaled but not yet fit to the DM 2002 cyclic companion;
  `-restraighten` exposes both span forms for calibration.
- **`-E` must match the wrapped bar's modulus** and **DM requires `-fy`.**
- **Not** plate buckling, **not** global member buckling (see §1).

## 12. References

- Dhakal, R.P. & Maekawa, K. (2002). *Reinforcing bar stress–strain relations
  including buckling.* J. Struct. Eng. 128(9), 1139–1147 (+ the cyclic companion).
- Gomes, A. & Appleton, J. (1997). *Nonlinear cyclic stress–strain relationship
  of reinforcing bars including buckling.* Eng. Struct. 19(10).
- OpenSees `ReinforcingSteel` (`Buckled_stress_Dhakal` / `_Gomes`) — the ported
  formulae.
- ADR: [[14_ladruno_rebar_buckling_adr]].
