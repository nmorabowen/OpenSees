---
title: LadrunoShellModifier — ETABS-style shell stiffness modifiers
project: Ladruno
status: shipped — PR #793
priority: medium
adr: ADR 91
tags:
  - section
  - shell
  - wall
  - cracked-section
  - etabs
  - aci318
  - asce41
---

# LadrunoShellModifier — ETABS-style shell stiffness modifiers

**What it is:** a section *decorator*. You wrap an existing order-8 plate/shell section and
scale its eight stiffness terms independently, plus its mass. It is the OpenSees equivalent of
ETABS's area section property modifiers, and it exists so elastic shell models of RC walls can
carry code-mandated **cracked-section** stiffness (ACI 318-25 §6.6.3.1.1, ASCE 41 Table 10-5)
without you having to fake it by editing `E` or the thickness.

Background: [[91_ladruno_shell_stiffness_modifiers_adr]].

---

## 1. The command

```tcl
section LadrunoShellModifier $tag $innerSecTag \
    <-f11 v> <-f22 v> <-f12 v> <-m11 v> <-m22 v> <-m12 v> <-v13 v> <-v23 v> <-mass v>
```

```python
ops.section("LadrunoShellModifier", tag, innerSecTag, "-f11", 0.35, "-f22", 0.35, ...)
```

All nine flags are optional, order-independent, and default to `1.0`. **All defaults is
byte-identical to the inner section** — that is the first thing the test suite checks, so you can
wrap unconditionally in a model generator and only pass the flags you actually mean.

Then use the wrapper's tag wherever you would have used the inner one:

```tcl
section ElasticMembranePlateSection 1  $E $nu $h $rho
section LadrunoShellModifier       10 1  -f11 0.35 -f22 0.35 -f12 0.35
element ShellMITC4  1  $n1 $n2 $n3 $n4  10
```

## 2. What the nine names mean

They are the ETABS names, and they map one-to-one onto the OpenSees plate resultant ordering:

| flag | ETABS meaning | scales the resultant |
|---|---|---|
| `-f11`, `-f22` | membrane direct | `FXX`, `FYY` |
| `-f12` | membrane **shear** (not the Poisson term) | `FXY` |
| `-m11`, `-m22` | bending | `MXX`, `MYY` |
| `-m12` | twisting | `MXY` |
| `-v13`, `-v23` | transverse (out-of-plane) shear | `VXZ`, `VYZ` |
| `-mass` | mass only | section `rho` |

`-f12` is the in-plane **shear** modifier. This trips people up: `f12` looks like it should be
the `1-2` coupling (the Poisson term), and it is not. That matches ETABS.

### The one that will actually bite you: `m` does nothing to a wall

**A wall or deep beam loaded in its own plane is cracked with `f11`, not `m11`.**

In-plane bending of a shell is carried by **membrane** action — `sigma_xx` varying over the
depth of the member. The `m` modifiers are *out-of-plane plate* bending, and a wall bending in
its own plane has none. So if you crack a shear wall with `m11 = m22 = m12 = 0.35` because those
are "the bending ones", **you have applied no cracking at all, and nothing warns you.**

This is measured, not asserted. A flexure-controlled cantilever (L/d = 10) modelled as a frame
and as a ShellMITC4 mesh of the same member:

| case | tip deflection | ratio vs its own gross |
|---|---|---|
| FRAME gross | 0.053333 m | 1.0000 |
| FRAME `A,I × 0.35` | 0.152381 m | 2.8571 |
| SHELL gross | 0.051911 m | 1.0000 |
| SHELL `f11,f22,f12 × 0.35` | 0.148318 m | 2.8571 |
| SHELL `m11,m22,m12 × 0.35` | 0.051911 m | **1.0000 — no change** |

`1/0.35 = 2.8571`. The frame route (ETABS frame modifiers scale `A` and `I33` directly) and the
shell membrane route land on the *same* softening to 8+ significant figures, and the `m` route
does exactly nothing. `v13`/`v23` likewise do nothing in plane.

The residual −2.7% between frame and shell is discretisation, not the modifier: a bilinear
membrane locks in bending, so a coarse shell mesh reads stiff. It converges away under
refinement (−12.6% at 20×2 → −0.5% at 120×24) while the softening ratio stays 2.8571 at *every*
mesh density. The modifier scales straight through the discretisation gap rather than distorting
it.

Pinned by `tests/test_ladrunoShellModifier_frame_equivalence.py` (G9/G10).

Out-of-plane modifiers are for shells actually loaded out of plane — slabs, and walls carrying
face loads or spanning between floors. A model with both actions wants both sets.

## 3. Recipes

**Cracked shear wall, in-plane (the common case).** ACI 318-25 §6.6.3.1.1 gives walls
0.70·Ig uncracked / 0.35·Ig cracked. For a wall modelled with shells, the in-plane stiffness is
the *membrane* block:

```tcl
section LadrunoShellModifier 10 1  -f11 0.35 -f22 0.35 -f12 0.35
```

**Cracked in-plane, near-gross out-of-plane.** Walls often crack in their own plane long before
they crack out of plane:

```tcl
section LadrunoShellModifier 10 1  -f11 0.35 -f22 0.35 -f12 0.35 \
                                   -m11 0.70 -m22 0.70 -m12 0.70
```

**Wall that should carry gravity but not in-plane shear.** A common ETABS idiom for isolating
load paths — keep the vertical membrane stiffness, remove the shear:

```tcl
section LadrunoShellModifier 10 1  -f12 0.01
```

Prefer a small number to a hard `0.0` unless you genuinely want the singularity — see §5.

**Mass without stiffness change** (e.g. tributary-mass tuning):

```tcl
section LadrunoShellModifier 10 1  -mass 0.85
```

## 4. It works on nonlinear sections too

The decorator is not restricted to `ElasticMembranePlateSection`. It wraps **any** order-8
plate section — `LayeredShell`, `PlateFiber`, `RCLayeredMembraneSection` — and it is exact for
them: the modifier is applied as a change of variables, so the tangent it reports is the true
consistent tangent of the modified law, not an approximation. There is no convergence penalty.

That said, modifiers are a *linear-elastic modelling device*. Applying them to a nonlinear
section rescales the strain that section sees, so it will yield at a different apparent
deformation. That is a legitimate thing to want (transverse-shear-only modification, say), but
it is a modelling decision you should make deliberately, not a free knob.

## 5. Zero is allowed, and it does what it says

A modifier of exactly `0.0` is accepted — ETABS permits it — and warns once per `section`
command naming the mode it kills. The result really is a singular section in that mode, and your
solve may fail with `factorization failed, matrix singular` depending on what else restrains the
model.

Note that zeroing `f11` removes the whole `FXX` row **and column**, including the Poisson
coupling to `FYY`. A shell patch with `f11 = 0` therefore loses more than one deformation mode.
If you want "very soft" rather than "absent", use `0.01`.

## 6. What it will refuse

Four named refusals, all of which fail at parse time with a message naming the cause:

| | condition |
|---|---|
| **R1** | the inner section is not order 8 (e.g. you passed a beam section) |
| **R2** | the inner section is order 8 but its response codes are not exactly `{FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ}` — it is not a plate section |
| **R3** | any modifier is negative |
| **R4** | the inner section tag does not exist |

## 7. There is no weight modifier — read this before looking for it

ETABS has ten modifiers; this has nine. The missing one is `weight`, and it is missing on
purpose.

In OpenSees the shell self-weight body force is derived from **the same `getRho()`** that builds
the mass matrix (`ShellMITC4.cpp:1725-1755`: `momentum = appliedB; momentum *= rhoH`). A
`-weight` flag could therefore not be independent of `-mass`; it would silently alias it. Shipping
an argument that quietly does something other than its name is worse than not shipping it.

**Scale self-weight at the load level instead** — the `eleLoad -type -selfWeight` factors, or the
load factor on the gravity pattern.

## 8. Reporting

No special handling needed. `section force` / `section deformation` recorders on a wrapped
section report the **modified** resultants, which is the ETABS reporting convention — the forces
you read back are consistent with the reduced stiffness you asked for.

## 9. One caveat worth knowing (ADR 91 OQ-1)

The modifiers are applied as a **congruence**: `D' = S·D·S` with `S = diag(√f11 … √v23)`. So the
diagonal terms scale exactly as you would expect (`f11` scales `D(0,0)` by `f11`), but the
off-diagonal Poisson coupling `D(0,1)` scales by `√(f11·f22)`.

When `f11 == f22` — which is the overwhelmingly common case, and every recipe in §3 above — this
is indistinguishable from scaling the whole membrane block by that factor, and the question does
not arise. It only matters if you set `f11` and `f22` to very different values.

The reason it is done this way rather than scaling only the diagonal: diagonal-only scaling
destroys positive definiteness. With `f11 = f22 = 0.1` and `ν = 0.2`, leaving the coupling at
full value makes the membrane block indefinite — and that is exactly the cracked-wall regime this
feature exists to serve. The congruence form is provably definiteness-preserving (Schur product
theorem). CSI's documentation does not state which convention ETABS uses; if you have an ETABS
model with strongly unequal `f11`/`f22` to compare against, that comparison is the open question
in ADR 91 and worth reporting back.

## 10. Verification

`tests/test_ladrunoShellModifier_section.py` and `..._structural.py`, 15 tests, ADR 91 §9:
identity at all-1.0 (vs the bare section, `ShellMITC4` + `ASDShellQ4`), the congruence against an
independent numpy `S·D·S` oracle with all nine modifiers distinct, exact analytic doubling at
uniform 0.5, positive-definiteness under `f11=f22=0.01`, parity with the upstream `Ep_mod`
scalar, mass scaling by √2 in eigenvalues, the four refusals, and nonlinear passthrough on
`LayeredShell`.

Mutation-gated as family `SHELLMOD` (ADR-87 D2) — see `Ladruno_scripts/mutation_build.bat` for
building a mutant on Windows.
