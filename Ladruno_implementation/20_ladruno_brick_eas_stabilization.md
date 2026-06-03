---
title: LadrunoBrick EAS — stabilization for inelastic localization
project: Ladruno
status: scoping
priority: low
owner: nmora
tags:
  - implementation
  - element
  - eas
  - stabilization
---

# LadrunoBrick EAS — stabilization for inelastic localization

> Scoping doc (no code yet). Follows the EAS-vs-ssp/bbar damage comparison in
> [[19_ladruno_brick_eas_simo_rifai]] (PR-2 follow-up), which found the true-EAS
> element (`-formulation eas`) **not robust on notched / high-strain-gradient
> inelastic problems**. This ADR scopes the fix.

## What

Add an optional **stabilization** to the true Simo–Rifai EAS element so it
traverses notched / localization-dominated **inelastic** problems (where the bare
E9 currently stalls) without sacrificing its smooth-problem accuracy or the patch
test. Default **off** (β=0 ⇒ bit-identical to today's `eas`); opt-in via a knob
(e.g. `-formulation eas -stab <β>`).

**In scope:** small-strain EAS, the stabilization term + its knob + validation
against the ADR 19 DEN-bar comparison. **Not in scope:** finite-strain (enhanced-`F`)
EAS and its separate compressive-hourglassing problem (ADR 19 §3); richer mode
sets E12/E21/E30 (ADR 19 §7).

## Why

The ADR 19 comparison (steel DEN bar, J2 + Lemaitre damage) established, with
evidence, that the bare E9 element:

- is **rank-sufficient elastically** — a free element has exactly 6 zero-energy
  modes and an eigen-spectrum identical to `std`/`bbar` (so this is **not** elastic
  hourglassing);
- handles **homogeneous** plasticity+damage fine (single element → full load,
  damage 0.48, = `bbar`);
- but **stalls just past yield onset on the notch** — a genuine instability of the
  enhanced modes under **non-homogeneous inelastic** tangents. An inner-Newton line
  search only partially mitigated it (reach 0.169→0.231, warnings 89k→28k; not a
  cure) at ~13× cost, so it was reverted.

`bbar`/`ssp` are robust there, so EAS is not *needed* for these problems — but EAS
is the most accurate member of the family on bending/incompressibility, and a
stabilized EAS that is *also* robust on inelastic localization would make it the
default solid element rather than a special-case one. That is the prize.

## Mechanism (the thing the stabilization must fix)

EAS's stability guarantee (Reddy & Simo 1995) is for the **elastic** operator. The
enhanced sub-problem stiffness

```
K_αα = ∫ Mᵀ C M dV
```

is positive-definite for the elastic `C`, but as the material tangent `C → C_ep`
(plastic) and especially `C_softening` (damage) degrades in a **high-gradient**
region, `K_αα` can lose positive-definiteness. The static condensation

```
K* = K_dd − K_dα K_αα⁻¹ K_αd
```

then amplifies that into a spurious soft/negative mode of the **global** tangent →
the global Newton cannot converge → stall. A stabilization must keep the enhanced
sub-problem well-posed (and `K*` sound) **without** perturbing the converged
solution on smooth problems.

## How — candidate approaches (literature + recommended first attempt)

The EAS hourglass/stability literature (mostly finite-strain, but the lever is the
same):

- **Reese–Wriggers (2000) / Reese (2005)** — a stabilization built from the
  hourglass modes, scaled by a stabilization parameter; the canonical EAS-hourglass
  cure.
- **Korelc–Šolinc–Wriggers (2010)** — an EAS brick with a *modified, automatically
  stabilized* enhanced field that stays stable in finite deformation.
- **Glaser–Armero (1997)** — early stabilized EAS for the compressive instability.
- **Wriggers–Reese (1996)** — first documented the instability (diagnostic, not a
  cure).

For our **small-strain inelastic** case the simplest effective lever is a
**`K_αα` stabilization** that keeps the enhanced sub-problem positive-definite:

```
K_αα_stab = K_αα + β · K_αα⁰
```

where `K_αα⁰ = ∫ Mᵀ C₀ M dV` is the **elastic** enhanced stiffness (constant,
cacheable at setDomain like the SSP `Kstab`), `C₀` the initial tangent, and `β` a
small dimensionless parameter (default 0). Properties to verify:

- **Patch test preserved:** on a constant-strain state `α → 0` regardless of `β`
  (the enhanced residual is already 0), so the stabilization does not pollute the
  patch test. *(Must confirm — the stabilization changes the `α` *path*, not the
  converged `α=0`, so it should hold.)*
- **Smooth-problem accuracy:** small `β` should leave bending/incompressibility
  results within tolerance of `β=0` (the enhanced modes still do their job; the
  stabilization only floors the sub-problem stiffness).
- **Robustness:** `β>0` keeps `K_αα` invertible through plastic/softening tangents
  ⇒ `K*` stays sound ⇒ the global Newton traverses the notch.

This is the **try-first** option (cheap, one cached matrix, one knob). If a constant
`β·K_αα⁰` blend over-stiffens the enhanced response (loses EAS accuracy) or still
fails, escalate to the Reese–Wriggers hourglass-mode stabilization (deformation-
dependent, more faithful but more code), then Korelc-style modified-`M`.

### Knob / API

```
element('LadrunoBrick', tag, n1..n8, matTag, '-formulation', 'eas' [, '-stab', beta])
```

`beta` default 0 (current behaviour, bit-identical). Document a recommended value
once tuned on the DEN bar.

## Testing

The merge gate is the **ADR 19 DEN-bar comparison, re-run with `-stab`**:

1. **Robustness (the headline):** stabilized `eas` traverses the notched bar to
   the same elongation as `bbar` (currently stalls at u≈0.17/1.5).
2. **Patch test still exact** (distorted 2×2×2, `∫M dV=0`) — `test_ladrunoBrick_eas.py`.
3. **Smooth accuracy preserved:** bending-beats-`std`, near-incompressible, and
   reduce-to-`std` results unchanged within tolerance vs `β=0`.
4. **β=0 is bit-identical** to the current element (regression guard).
5. **β-sweep:** show the minimum `β` that gives robustness, and that accuracy
   degrades gracefully (not a cliff) as `β` grows.

## Risks / open questions

> [!question]
> Does a constant `β·K_αα⁰` blend preserve the patch test exactly, or does it bias
> the converged `α` away from 0 on a constant-strain state? (Believed safe — the
> enhanced residual is already 0 there — but must be proven, not assumed.)

> [!question]
> Is a single scalar `β` enough, or does robustness need a *deformation-dependent*
> (Reese–Wriggers) stabilization to avoid over-stiffening bending while still
> flooring the notch sub-problem? Start scalar; escalate only if it fails the
> accuracy gate.

- **Over-stiffening:** too much stabilization turns `eas` back toward a
  fully-integrated/`std`-like response, losing the reason to use EAS. The β-sweep
  must show a usable window.
- **Tuning vs automatic:** a hand-tuned `β` is a footgun; the Korelc approach is
  parameter-free but much more code. Decide after the scalar experiment.
- **Interaction with `-autoRegularization`:** the stabilization must not interfere
  with the crack-band lch scaling (it acts on the enhanced *kinematics*, the
  regularization on the *material* — believed orthogonal).

## Implementation log

*(empty — scoping only. Pick up after deciding to invest; start with the scalar
`β·K_αα⁰` experiment on the DEN bar, gated by the §Testing battery.)*
