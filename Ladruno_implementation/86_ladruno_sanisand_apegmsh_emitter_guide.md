---
title: LadrunoSANISAND — implementation guide for apeGmsh
project: Ladruno
status: ready to implement — the fork dependency is MERGED (PR #767, #768)
priority: medium
adr: ADR 86
tags:
  - apegmsh
  - emitter
  - material
  - sanisand
  - manzari-dafalias
  - soil
---

# LadrunoSANISAND — implementation guide for apeGmsh

**What we are asking for:** `nDMaterial LadrunoSANISAND` as a first-class material in the
`apeSees(fem)` bridge, emitted to both the Tcl and openseespy backends.

Most of this note is not about the emitter, though. It is about **three deck-construction rules
that any Manzari-family emitter must follow** — one of which silently changes the material's
bounding surface, and one of which can make the sand behave elastically for an entire analysis
while reporting success. Those apply to `ManzariDafalias` too, and apeGmsh has almost certainly
been emitting decks that hit at least the first.

Background: [[86_ladruno_sanisand_adr]] (spec), [[86_ladruno_sanisand_tims_report]] (findings, in
consumer language).

---

## 1. The material

A thin subclass of `ManzariDafalias` whose only difference is that the two low-stress constants
`m_Presidual` and `m_Pmin` — hardcoded upstream, not settable, not on the wire — become optional
deck arguments, carried across the wire, and echoed at construction. Defaults `0.0` and
`1.0e-3 * P_atm`.

```tcl
nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m \
    $h0 $ch $nb $A0 $nd $z_max $cz $Rho \
    <$IntScheme $TanType $JacoType $TolF $TolR> \
    <-Presidual $pr> <-Pmin $pmin>
```

The **first 18 positional arguments and the 5 positional optionals are identical to
`ManzariDafalias`**, so if apeGmsh already emits Manzari, this is the same parameter object with a
different command name plus two flags. Class tags 33019 / 33020 / 33021 (base / 3D / PlaneStrain).

Registered in **both** interpreters, so both emitters work with no further wiring.

### Argument-ordering rule the parser enforces

Positional optionals must come **before** any `-flag`. A positional after a flag is a **hard parse
error**, by design:

```tcl
# OK
... $Rho 1 0 1 1e-7 1e-7 -Presidual 0.0 -Pmin 0.101
# OK  (flags only, positionals defaulted)
... $Rho -Presidual 0.0 -Pmin 0.101
# ERROR — positional after a flag
... $Rho -Presidual 0.0 1 0 1 1e-7 1e-7
```

This matters for an emitter that assembles arguments from a dict or a dataclass: **emit the
positional block first, then the flags.** Do not let key ordering decide it.

### Always emit both flags explicitly, even for defaults

Two reasons, and the second is the one that bites:

- The material echoes what it is running at construction, so an explicit deck makes the log
  self-documenting.
- **The `-Pmin` default is `1.0e-3 * P_atm`, ten times vanilla `ManzariDafalias`'s `1.0e-4 *
  P_atm`.** Any A/B comparison between the two materials that does not pin `-Pmin` in *both* legs
  is moving two variables, not one.

### `-honorTolR`

Accepts **`0` or `1`** as of ADR-86 PR-3, which wires the base flag seam PR-2 opened.
(Until PR-3 it accepted only `0`, and `-honorTolR 1` was a hard parse error rather than a silent
no-op. If you are reading an older copy of this guide, that is why.)

`0` keeps vanilla behaviour: `ModifiedEuler` runs its hardcoded substep tolerance `TolE = 1e-4`
and ignores the deck TolR. `1` makes it honour the deck TolR instead.

**Emit `0` unless you have a reason.** Two cautions if you emit `1`:

- It is read at **exactly one site**, inside `ModifiedEuler`, so it does nothing on any IntScheme
  that does not route there -- 2, 3, 4, 5, 6, 7, 8, 9 and 45 (45 already honours TolR anyway).
  The material warns at construction when you ask for it on one of those.
- The parser default is `TolR = 1e-7`, so `-honorTolR 1` on a default deck tightens the substep
  tolerance by 1000x. Measured on our confine-first deck: the answer moves 2.641e-04 at
  `TolR = 1e-6`, and is bit-identical at `TolR = 1e-4` (both operands of the seam are then the
  same double). Emit an explicit TolR if you emit the flag.

---

## 2. Three deck rules — these apply to `ManzariDafalias` as well

### 2.1 Do not shear during the elastic stage

**What goes wrong:** if the stress ratio `eta = q/p` is already above `M_c` at the moment the deck
calls `updateMaterialStage ... -stage 1`, `Elastic2Plastic` **silently raises `M_c`** to
accommodate it — measured on our own test deck at **1.3309 -> 1.99878, a 50 % inflation, once per
Gauss point.** From then on the material runs against a bounding surface that is not the one that
was calibrated, and nothing in the output says so.

**Rule:** reach the target stress state under stage 0 with a path that keeps `eta` well under
`M_c`, then flip, then apply the deviatoric history.

Ordinary `K0` gravity is safe — at `K0 = 0.5`, `eta = 0.75` against `M_c = 1.331`. What is not safe
is applying a deviatoric *strain* ramp during stage 0, which is exactly what a naive
single-element calibration deck does.

**Recommended emitted ordering:**

1. stage 0, isotropic or `K0` consolidation to the target confinement, run to equilibrium
2. `updateMaterialStage -material <tag> -stage 1`
3. the deviatoric / cyclic history

Our own test decks were rebuilt around exactly this ordering, and it produces **zero** inflation
events where the previous ordering produced 8 per leg.

### 2.2 `updateMaterialStage` must be explicit, per tag, at every stage boundary

`mElastFlag` is a **`static`** member of `ManzariDafalias`, so constructing *any* material of the
family resets the stage flag for *every* instance in the process.

Consequences for an emitter:

- Never rely on stage state surviving a later material construction.
- If materials are created lazily, or in a different order than the stage calls, the flag can be
  clobbered after it was set.
- Emit `updateMaterialStage -material <tag> -stage <n>` **explicitly for each tag**, at each
  boundary — not once for the model.

This is an upstream property of the whole UW soil family (`CycLiqCP`, `SAniSandMS`,
`BoundingCamClay`, `DruckerPragerThermal` all share it), not something the fork can fix in one
material.

### 2.3 A proportional strain ramp with `p_residual = 0` never yields

**What goes wrong:** every *elastic* step re-sets the back-stress ratio to `s/(p + p_r)`. On a
perfectly proportional strain path the stress path is exactly radial, so with `p_r = 0` that ratio
is constant along the ray — the yield surface is welded to the stress and the yield function stays
negative **forever**. Zero plastic strain, on every integration scheme, with the analysis
converging and reporting success throughout.

**Where an emitter hits this:** single-element calibration or verification decks — staged gravity,
then a single proportional strain ramp. That is a very natural thing for a mesh tool to generate.

**Mitigation:** the ordering in §2.1 fixes it too. Confining hydrostatically first leaves
`alpha = 0` at the flip, so the material yields normally. If apeGmsh ships single-element
verification decks for this family, they should confine first.

Note this is a property of *proportional* decks, not of `p_r = 0`. `NTUASand02` ships `p_r = 0`
with identical code and is used for staged foundation work, where the stress direction changes
constantly and the freeze never occurs.

---

## 3. Two things that affect apeGmsh regardless of whether you adopt the new material

### 3.1 Golden-file regression tests on `IntScheme 45` will move

PR-2 repaired a void-ratio interpolant in `ManzariDafalias`'s substepping loops. Measured:

| scheme | shift |
|---|---|
| `IntScheme 1` (ModifiedEuler, the default) | 2.8e-6 relative |
| `IntScheme 3` (RungeKutta4) | none |
| **`IntScheme 45`** (RungeKutta45) | **6.33e-3 — about 0.6 %** |

If apeGmsh holds golden files produced by RK45 Manzari decks, they need regenerating. This is a
correction, not a regression.

### 3.2 Vanilla `ManzariDafalias` MP decks are not equivalent to their serial twins

`ManzariDafalias`'s wire format has **no slot for `m_Presidual`**, and the broker path leaves it at
zero. So an `OpenSeesMP` worker — or any model restored from a `FileDatastore` — runs at
`p_r = 0` while the serial process beside it runs at 1.01 kPa. Measured: a restored vanilla
material changes its stress **5.5 %** the instant it is restored.

`LadrunoSANISAND` carries the constant on the wire and does not have this problem. **If apeGmsh
emits MP or database-restore workflows for this material family, prefer the new class.**

---

## 4. Suggested shape in the bridge

A typed primitive mirroring the existing Manzari one, plus the two constants:

```python
LadrunoSANISAND(
    tag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m,
    h0, ch, nb, A0, nd, z_max, cz, rho,
    int_scheme=1, tan_type=0, jaco_type=1, tol_f=1e-7, tol_r=1e-7,
    p_residual=0.0,          # emitted as -Presidual, always
    p_min=None,              # None -> 1e-3 * P_atm; emitted as -Pmin, always
)
```

Emission order: `tag`, the 18 positional doubles, the 5 positional optionals, then the flags.

Reference parameter set (Gorini's, calibrated, used throughout our verification):

```
G0 264.32   nu 0.3129   e_init 0.6944   Mc 1.33090   c 0.71    lambda_c 0.027
e0 0.83     ksi 0.45    P_atm 101.0     m 0.005      h0 1.3    ch 0.968
nb 3.5      A0 0.05     nd 5.75         z_max 12.5   cz 1100.0 Rho 2.0
```

---

## 5. What we are NOT asking for

- No changes to `g.loads`, `g.constraints`, `g.masses` or the FEMData broker — this is a material,
  it lands in the material registry like any other.
- No new ndf handling. `LadrunoSANISAND3D` (33020) and `LadrunoSANISANDPlaneStrain` (33021) are
  requested through the ordinary `getCopy("ThreeDimensional")` / `"PlaneStrain"` type strings, the
  same as vanilla Manzari.
- No parallel-specific wiring. The class is registered in the object broker, so MP and database
  restore work through the normal path.

---

## Log

- 2026-08-27 — Written after PR #767 and PR #768 merged.
