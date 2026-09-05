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
    <-Presidual $pr> <-Pmin $pmin> <-honorTolR 0|1> <-maxSubsteps $n>
```

The **first 18 positional arguments and the 5 positional optionals occupy the same slots as
`ManzariDafalias`**, so if apeGmsh already emits Manzari, this is the same parameter object with a
different command name plus the flags. Class tags 33019 / 33020 / 33021 (base / 3D / PlaneStrain).

> [!warning] One positional DEFAULT differs, as of ADR-86b: **`TanType`**
> `OPS_ManzariDafaliasMaterial` defaults `TanType = 0`, which is the **ELASTIC** tangent
> (`ManzariDafalias3D::getTangent()` returns `mCe`), so a deck that emits only the 18 positional
> parameters runs `algorithm Newton` as de-facto **modified Newton**. `LadrunoSANISAND` now defaults
> it to **2** (the consistent elastoplastic tangent). Vanilla is deliberately unchanged — its golden
> files depend on 0, so the two parsers now disagree on this one slot on purpose.
>
> **Emit `TanType` explicitly anyway.** An emitter that names all five positionals is immune to this
> and to any future default move, and the material echoes the tangent it will run at construction.
> **`TanType 2` REQUIRES an unsymmetric system.** `mCep_Consistent` is unsymmetric under a
> non-associated flow rule, so emit `system FullGeneral`, `UmfPack`, or
> `Pardiso -matrixType 0` — a symmetric solver would be factorising a different matrix than the
> tangent describes. There is **no runtime guard** for this and there will not be: the material
> cannot see the SOE. It is an emitter rule.
>
> It changes the **iteration path only** — the substepped stress update is untouched — so it is free
> accuracy-wise. Two measurements, and they are different claims: the answer moved 0.52 % at a
> matched-settlement checkpoint on the WP-A2 footing deck against a ~7x wall-time saving (a
> *comparison across a whole configuration change*), and on a single-element free-DOF drained
> triaxial pushed through 40 accumulated `LoadControl` steps the two tangents converge to the
> **same** answer to within 4.5e-3 relative in displacement (7.0e-4 in stress) — MEASURED, not the
> two-orders-tighter-than-tolerance figure an earlier draft of this gate assumed. Over 40 steps of
> SANISAND's path-dependent fabric/backstress evolution the answer discrepancy tracks the deck's own
> `1e-3`-relative `NormUnbalance` tolerance (roughly 4-5x it, not a small fraction of it), so the
> gate's floor is set at `10x` that tolerance rather than as an independent constant
> (`tests/test_ladruno_sanisand_integrator.py::test_tantype_does_not_change_the_converged_answer`).
> The second measurement is the one that warrants moving a default.

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

### `-maxSubsteps` (ADR-86b)

A **substep-COUNT cap** on `ManzariDafalias::ModifiedEuler`. Default **`0` = uncapped**, which is
vanilla, so an emitter that omits the flag changes nothing at all.

Vanilla bounds the adaptive substep only from *below* (`dT_min = 1e-6`), never by a count, so one
Gauss-point update can legally cost ~10^6 return maps. On a softening boundary-value problem it
does: ADR-90 GATE U measured single `analyze(1)` calls of **11 to 34 minutes** on a strip footing,
with the stepping controller using **0 of its 80** subdivisions — the integrator never *failed*, it
merely did not terminate in useful time, so nothing upstream could react.

With a positive cap, an update that exceeds it **returns failure** rather than force-accepting: it
stops at `T < 1`, leaves the **committed** state untouched, and returns a refusal so the element can
fail the step and the analysis' own step-cut / subdivision logic gets its chance. One throttled
`opserr` line per process (budget 10) names the material tag, `T` and `dT`.

> [!warning] PRECONDITION — do not emit a cap without also emitting the right element
> **A cap is only safe under an element that PROPAGATES a material refusal.** Today that is **`LadrunoBrick` only**. Under vanilla `Brick` (its `update()` discards the return code), `BrickUP` / `QuadUP` (`setTrialStrain` is called inside a *void* `formResidAndTangent`, `BrickUP.cpp:1069`) and `stdBrick`, a capped update returns early at `T < 1` and the element assembles a **partially-integrated** stress / `alpha` / `fabric` and a partial `aCep_Consistent` as if converged. That is strictly **worse** than the un-capped force-accept it replaces, which at least always drove `T` to 1. **On those elements a capped run is INVALID, not merely un-cut.**
>
> The material cannot see which element holds it, so **nothing checks this at run time**. It is an
> emitter rule. The construction `Print()` record states it; the `opserr` line when the cap fires
> states it; neither can enforce it.

- **Only emit it on IntScheme 0 or 1.** Like `-honorTolR`, the seam is read at exactly one site,
  inside `ModifiedEuler`. The material warns at construction if you ask for it anywhere else.
  (IntScheme 7 is *called* `INT_MAXSTR_MFE` and does **not** route there — `MaxStrainInc` has no
  case for it and falls through to `ForwardEuler`. Read the switch, not the name.)
- **Emit `LadrunoBrick`** if you emit the cap — see the precondition above.
- `eleResponse <ele> material <gp> substeps` returns `[substeps_taken, cap_hit]` for the last update
  at that integration point, so a driver can size the cap from a measurement instead of a guess.

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

---

## 6. Convergence-test guidance (ADR-86b, from ADR-90 GATE U) — measured, not advisory

**Do not emit `NormDispIncr` for a SANISAND deck.** On this material it is not merely tight, it is
**unreachable**, and it is **not mesh-neutral** — the second of which is disqualifying inside a
mesh-convergence study.

**Why unreachable.** The default integration scheme is a *substepped* `ModifiedEuler` return with a
hardcoded substep tolerance `TolE = 1e-4` (exactly what `-honorTolR` above exists to expose). The
stress a Gauss point returns is therefore only a **piecewise-smooth** function of the strain
increment; the assembled residual inherits that, and Newton stops converging quadratically and
stalls. Measured on a strip-footing leg (h0 = 0.5, 390 hexes), over 47 failed convergence attempts,
the residual displacement norm **stalls at a median of 6.6e-7 m** (min 3.4e-8, max 1.3e-4). The run
nominally set to `NormDispIncr 1e-8` was in fact carried by the relaxed third rung of its algorithm
ladder on 18 of its 26 steps — 65 of every 125 state-determination passes spent failing two rungs
that could not succeed. A study in that state is measuring its own convergence test.

**Why not mesh-neutral.** The displacement norm runs over the **free DOFs**, and there are 3.6x more
of them at `h0 = 0.25` than at `h0 = 1.0`, so the same nominal number is a different physical
requirement on every mesh of a refinement sequence.

**Emit `NormUnbalance`, scaled to a deck-intrinsic force** — the model's own weight `gamma*V`, or
the applied load. That is what the fork's own R3 collapse gate actually uses
(`tests/test_r3_prandtl_collapse_gate.py`: `tol = 1.0e-5 * want`), despite "NormDispIncr per the R3
gate" appearing in more than one downstream brief. A force tolerance scaled by `gamma*V` is
identical on all three meshes by construction.

Measured on one deck at a fixed wall budget, same binary:

| test / tolerance | reach `s/B` | steps | subdivisions | relaxed rung |
|---|---|---|---|---|
| `NormDispIncr 1e-8` m (TanType 0) | 0.00106 | 26 | 1 | 18/26 |
| `NormUnbalance 1e-6 · gamma·V` (TanType 2) | 0.00218 | 31 | 0 | 11/31 |
| `NormUnbalance 1e-5 · gamma·V` (TanType 2) | 0.00442 | 37 | 0 | — |

and the answer moved by a **median of 0.3–0.75 %** at matched settlement between all three: the
looser setting buys reach, not error.

**There is no runtime warning for this, deliberately.** It is a deck rule, and a material cannot see
what convergence test the deck installed.

---

## 7. `-Presidual` on a deep push — the decision the consumer owes (ADR-86 D5a tripwire)

**Nothing changed in the code.** The fork default stays **`-Presidual 0.0`** (a cohesionless sand
has no cohesion), and ADR-86b does not move it. This section exists so the choice is made
consciously rather than by drift.

At `-Presidual 0` the low-`p` floor is `-Pmin` alone, `1e-3 * P_atm = 0.101` kPa. On a **dense**
(`e_init = 0.60`, `psi = -0.22`) sand the dilating soil under a footing edge drives a shallow Gauss
point onto that floor and the material logs

```
WARNING ManzariDafalias::ModifiedEuler() - material tag 1: mean stress p = 0.100973 is below
the floor m_Pmin + m_Presidual = 0.101; CLAMPING the stress to p = 0.101 (deviator preserved).
The result at this integration point is set by the clamp, not by the model.
```

**Measured onset: `s/B ≈ 0.0153`, on the coarsest mesh (`h0 = 1.0`) and the dense density only.**
The loose (`e_init = 0.6944`) legs and both finer meshes fired **zero** clamps over the same and
deeper settlements — so it is a dense-dilatant / coarse-element / large-settlement effect, not a
property of the deck. It sits past every number ADR-90 GATE U reported. A **faster** deck reaches it
next, which is precisely what the ADR-86b substep cap and `TanType` default are for.

Three legitimate answers; **pick one and disclose it**, do not leave it implicit:

1. **A declared non-zero `-Presidual`.** It restores a numerical regulariser as well as an apparent
   cohesion, and shifts `M^b` at low `p`. Tripwire: `p_residual` also bounds the `D_factor`
   dilatancy sigmoid **from below** (`86_ladruno_sanisand_pr3_tripwire_memo.md` measures vanilla's
   1.01 kPa holding it at 0.4278 against 4.83e-4 at `p_r = 0`, a factor of 886), so this is not a
   one-parameter move.
2. **A small surcharge** on the free surface, keeping the shallow points confined. Changes the
   boundary-value problem, not the material — usually the honest choice for a footing study.
3. **An explicitly accepted clamp**, with the event count reported alongside the result.

**It is silent in Python.** The message goes to `opserr`. Capture and COUNT it —
`ops.logFile(path, '-noEcho')`, run, redirect away, then read the file — rather than hoping a human
sees the console. The gravity state is no warning at all: the shallowest Gauss point sat at
1.56–6.25 kPa, 15–60x the floor, on every mesh, before the push ever started.

## Log

- 2026-08-27 — Written after PR #767 and PR #768 merged.
- 2026-09-05 — ADR-86b: `-maxSubsteps`, the `TanType` default move (0 → 2), section 6
  (convergence-test guidance) and section 7 (the `-Presidual` decision). All four come from
  ADR-90 GATE U (`_adr90_tau0_qu_band.md`).
