---
title: LadrunoPorousOverlay — user guide (persistent-fluid staggered u-p)
project: Ladruno
status: reference
priority: high
adr: ADR-73
pattern_tag: 33022
tags:
  - porous-media
  - biot
  - u-p
  - consolidation
  - liquefaction
  - fixed-stress
  - staggered
  - explicit
  - load-pattern
---

# LadrunoPorousOverlay — user guide

`LadrunoPorousOverlay` (PATTERN tag **33022**, ADR-73) is a **pore-pressure field
with its own life-cycle**. It is *not* an element: it is a `LoadPattern` that owns
a separately-discretized pressure field over a snapshot of your existing solid
mesh, advances that field with a fixed-stress staggered split, and pushes the
pore-pressure force `+Q·p` back onto the skeleton nodes each step. The solid
elements never gain a pressure DOF — any 3/4-node 2D or 8-node 3D solid cell
(`LadrunoQuad`/`CST`/`Brick`, vanilla `quad`/`stdBrick`/`SSPquad`, all at
`ndf = ndm`) qualifies, with **zero changes to the element**.

Two things this buys you that a monolithic u-p element cannot: (1) the fluid
survives `remove element`, so gap flow / sand boils / staged excavation are honest
(the monolithic element deletes the fluid continuum with the skeleton — measured
14× trapped-pressure error in the unsafe direction, ADR §1.1); (2) an **explicit
lane** — the solid marches under `CentralDifferenceLadruno` matrix-free-diagonal
while the fluid solves implicitly at commit (ADR §1.1/§3.4). For standard implicit
statics and the validated benchmark battery, the monolithic `LadrunoUP` element
remains the tool — see the division of labor in §1.

---

## Contents

1. [What it is / when to use it](#1-what-it-is--when-to-use-it)
2. [Quick start](#2-quick-start)
3. [Command reference](#3-command-reference)
4. [Initialization recipes](#4-initialization-recipes)
5. [Layered profiles](#5-layered-profiles)
6. [The iterated driver (P2)](#6-the-iterated-driver-p2)
7. [The explicit lane (P3)](#7-the-explicit-lane-p3)
8. [Why not CentralDifference + u-p elements](#8-why-not-centraldifference--u-p-elements)
9. [Recording p](#9-recording-p)
10. [Element removal](#10-element-removal)
11. [Quirks / troubleshooting](#11-quirks--troubleshooting)

---

## 1. What it is / when to use it

The overlay and the monolithic `LadrunoUP` element (ADR-71) are siblings that solve
the same Biot u-p physics through the same shared kernel (`LadrunoUPKernel.h`), but
with opposite trade-offs. `LadrunoUP` ties the fluid to the element and solves one
coupled unsymmetric system; the overlay decouples the fluid life-cycle and solves
a small **SPD** fluid system separately (symmetric solvers return; factorization
measured ×3.3 cheaper than the coupled unsymmetric solve at 40×40, ADR §1.2/§3.3).

### Division of labor vs LadrunoUP

| Job | Use | Why |
|---|---|---|
| Implicit statics, steady seepage, consolidation on an intact mesh | **LadrunoUP** | THE reference tool; Taylor–Hood, the validated B1–B5 battery live here (ADR §2.6/§6) |
| Element removal below the water table (gap flow, sand boils, staged excavation, progressive collapse) | **overlay** | fluid survives `remove element` — the monolithic element seals the flow path (ADR §1.1) |
| Explicit dynamics of a saturated body (large models, `CentralDifferenceLadruno`) | **overlay** (`-fsL zero`) | solid stays matrix-free diagonal; fluid solves implicitly at commit — LadrunoUP cannot ride the fork explicit stack (ADR §1.1/§3.4) |
| Companion pore field to a solid model you do not want to re-mesh with u-p elements | **overlay** | element-type-agnostic snapshot; zero element changes |

The overlay does **not** deprecate `LadrunoUP`. The overlay's gates cross-check
*against* it as the reference operator (spike E6 methodology; ADR §6).

### The one hard modeling rule

**One overlay per hydraulically connected water body.** The p-field lives
per-overlay (its own CSR system, its own drained set); there is no flow between two
overlays by construction. A layered deposit under one water table is therefore
**ONE** overlay spanning all layers, with per-layer permeability and porosity via
`-layer` blocks (§5) — not one overlay per layer. Two overlays claiming the same
element is a snapshot-time fatal (double-counted fluid). See §11.

---

## 2. Quick start

Minimal 1-D Terzaghi column: a strip of plain `quad` cells at `ndf = 2`, a strip
load on top, drained at the top nodes, overlay over every cell. The overlay applies
`+Q·p`; the solid analysis is ordinary `Transient`/`Newmark`.

**openseespy**

```python
import opensees.openseespy as ops

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
# ... build nodes, plain 'quad' elements (ndf=2), fix the base, strip load on top ...

# gather region element tags and the top (drained) node tags
region_eles = [...]        # every solid cell the fluid lives in
top_nodes   = [...]        # p-fixed boundary (>= 1 per connected region for statics)

ops.pattern("Plain", 1, 1)         # the mechanical load pattern (strip load), separate
# ... ops.load(...) on the loaded nodes ...

ops.pattern("LadrunoPorousOverlay", 77,
            "-region",  *region_eles,
            "-Kf",      2.2e9,      # fluid bulk modulus
            "-rhoF",    1000.0,     # fluid density
            "-perm",    1e-9, 1e-9, # k-bar = k_hydraulic / gamma_w, ndm values
            "-poro",    0.4,        # porosity n
            "-moduli",  1.0e7, 0.25,# REQUIRED: skeleton E, nu for L / stab moduli
            "-drained", *top_nodes) # p = 0 boundary

ops.system("ProfileSPD")           # symmetric solver is fine — fluid is SPD-separate
ops.numberer("Plain")
ops.constraints("Transformation")
ops.algorithm("Linear")
ops.integrator("Newmark", 0.6, 0.3025)
ops.analysis("Transient")

for _ in range(nsteps):
    assert ops.analyze(1, dt) == 0
```

**Tcl**

```tcl
pattern LadrunoPorousOverlay 77 \
    -region {*}$regionEles \
    -Kf 2.2e9 -rhoF 1000.0 -perm 1e-9 1e-9 -poro 0.4 -moduli 1.0e7 0.25 \
    -drained {*}$topNodes
system ProfileSPD
numberer Plain
constraints Transformation
algorithm Linear
integrator Newmark 0.6 0.3025
analysis Transient
analyze $nsteps $dt
```

Plain `analyze` gives the **fs1** lane (one fixed-stress sweep per step, O(Δt)
splitting error). No driver needed. The P1 physics battery measures the overlay
reproducing the frozen fixed-stress toy to **3.46e-7** under a massless solid, with
splitting error 1.00× the frozen reference (ADR §12 P1).

---

## 3. Command reference

The overlay is a pattern command; the type token is consumed before your `$tag`.
**Unknown flags are FATAL** (a mistyped porous flag silently changes physics —
ADR-71/73 house rule). Verified against
`SRC/domain/pattern/ladrunoPorousOverlay/OPS_LadrunoPorousOverlay.cpp`.

```
pattern LadrunoPorousOverlay $tag \
    -region $e1 $e2 ...                       ;# REQUIRED: solid cells (flat list / {*}$list)
    -Kf $Kf -rhoF $rhoF -poro $n \            ;# REQUIRED: fluid bulk, fluid density, porosity
    -perm $k1 $k2 <$k3> \                     ;# REQUIRED: k-bar = k_hydraulic/gamma_w, ndm values
    -moduli $E $nu \                          ;# REQUIRED (v1): skeleton moduli for L / stab
    -drained $nd1 $nd2 ... \                  ;# REQUIRED: p-fixed nodes (>=1 per connected region for statics)
    <-alpha $biotAlpha>                       ;# default 1.0
    <-Ks $Ks>                                 ;# grain bulk; default <=0 => incompressible grains
    <-thick $t>                               ;# 2D only (fatal in 3D); default 1.0 — MUST match the solid thickness
    <-layer $e1 ... -perm $k1 $k2 <$k3> <-poro $n> <-moduli $E $nu>> ...   ;# per-layer override, repeatable (§5)
    <-stab auto <$a0> | off | $alpha>         ;# equal-order stabilization; default auto (a0 = 0.25)
    <-fsL classic | oedometric | zero | $scale>  ;# split modulus; default classic (§3, §6, §7)
    <-staticMode hold | steady>               ;# static-analysis commit policy (§4); DEFAULT MARCHES
    <-onRemoval keep | drain $kFactor>        ;# fluid life-cycle at element death; default keep (§10)
    <-fluidUpdate implicit|explicit>          ;# default implicit; explicit = P3b matrix-free lane (§7)
    <-subcycle auto | $N>                     ;# fluid advance every N commits; default 1 (§7, §10)
    <-pInit steady | hydrostatic $gw <$zw> | $nd $val ...>   ;# initial p (§4)
    <-fluidBody $b1 $b2 <$b3>>                ;# seepage body drive (acceleration form)
    <-dynSeepage on | off>                    ;# default off (ADR-71 §12 amendment)
    <-record $file <$everyN>>                 ;# CSV gate output (§9)
```

Key defaults and gotchas, all verified in the parser:

- **Required:** `-region`, `-Kf`, `-rhoF`, `-poro`, `-perm`, `-moduli`, `-drained`.
  A missing one names itself in the fatal.
- **`-moduli $E $nu` is required in v1.** The overlay does not own the solids'
  materials, so the L and stabilization moduli come from explicit user input
  (deliberate v1 simplification; ADR plan WP1). Validated `E > 0`, `-1 < nu < 0.5`.
- **`-perm` takes exactly `ndm` values** and uses the ADR-71 convention
  `k-bar = k_hydraulic / gamma_w` (velocity/pressure-gradient units). Values `>= 0`.
- **`-thick` is 2D-only** and *must match your solid elements' thickness* — the
  overlay cannot read thickness from an arbitrary solid element. It is fatal in 3D.
- **`-staticMode` default MARCHES the fluid every commit.** This is correct for
  transient consolidation. It is **wrong under a static analysis** (§4) — there is
  no auto-detect; the overlay warns loudly on the first march instead.
- **`-poro` must satisfy `0 < n <= alpha`** (storage `1/Qbar` domain); a violation
  is fatal. Same check per-layer.
- **`-substep` was REMOVED** (its accuracy claim was refuted at P0 E7.3b — error
  flat over M = 1…10). The parser rejects it with a pointer to the ADR §12 P0 entry.
  Use `-subcycle` for the explicit lane.
- **`-fluidUpdate explicit`** — the fully-matrix-free Xu-2021-class lane, shipped
  at P3b (lumped-S* forward p-step at load-application time; §7 owns the recipe
  and advisories). `-fluidUpdate implicit` is the default.

**The pattern's load factor / TimeSeries is ignored by design.** The overlay owns
its own force amplitudes — the applied nodal force is always `+Q·p_committed` at
full amplitude, the pore-pressure field itself, not a factored load. Attaching a
`timeSeries` or expecting `loadConst`-style ramping does nothing; a one-time notice
prints. See §11.

---

## 4. Initialization recipes

### The static-time trap (read this before gravity staging)

Under a **static** analysis the domain "time" is the **load factor**, not seconds:
`LoadControl` hands the commit hook `dt = Δλ`. The overlay's default
`-staticMode` marches the fluid, so a static gravity stage would integrate the mass
balance over `Δt = Δλ` — silently wrong physics (measured: node p ran 24.8 → 49.9
over λ 0.1 → 1.0 with no error, ADR §11 (P1 panel) finding P-1). A `LoadPattern` cannot see
the integrator, so this cannot be auto-detected; the overlay prints a **loud
one-time advisory on the first real march** naming the `dt`, and you fix it by
choosing a `-staticMode`:

- **`-staticMode hold`** — freeze p during static commits (forces still applied from
  the committed p; no fluid marching).
- **`-staticMode steady`** — re-solve steady seepage `H·p = f_seep` per static commit.

Transient fluid marching happens only under transient analyses. The recommended
gravity-init sequence:

```tcl
# stage 1: static gravity with the fluid held/steady (NOT marching on load factor)
pattern LadrunoPorousOverlay 77 ... -staticMode steady -pInit steady ...
integrator LoadControl 0.1
analysis Static
analyze 10                         ;# gravity brought on over lambda 0 -> 1

# stage 2: switch to transient consolidation / shaking — fluid now marches on real dt
wipeAnalysis
integrator Newmark 0.6 0.3025
analysis Transient
analyze $nsteps $dt
```

### Setting the initial pressure

The overlay owns H, so it solves its own steady state — no `InitialStateAnalysis`
interplay and no staged-`sp` Transformation trap, because the solid analysis never
sees a p DOF (ADR §3.5):

- `-pInit steady` — solve `H·p = f_seep` with the drained set.
- `-pInit hydrostatic $gw <$zw>` — hydrostatic column, unit weight `$gw`, optional
  phreatic level `$zw` (default 0).
- `-pInit $nd1 $val1 $nd2 $val2 ...` — explicit nodal values (region nodes).

### The mixture-gravity / staggered-twin body-force trap

Solid self-weight stays **user-side with mixture density** (unchanged family
convention). The overlay only adds seepage body terms via `-fluidBody`
(acceleration form, scaled by `rhoF` internally — matching `LadrunoUP -fluidBody`).

If you build the staggered twin of a monolithic `LadrunoUP` model (plain `quad` +
overlay), **do not copy the monolithic `-body` value into the quad's trailing
`b1 b2`**: the plain quad's `b` is body FORCE per unit VOLUME, but `LadrunoUP -body`
is an ACCELERATION (multiplied by `rho_mix` inside the element). Copying the number
gives the solid weight per unit volume instead of `rho_mix·accel`, and with a
hydrostatic overlay `+Q·p` the net solid load can cancel to ~zero. Measured before
the fix: settlement **1e-9 vs 3.5e-4** (u-trace rel diff exactly 1.0), nothing
errored (ADR §12 P2 item 4; LEDGER_quirks 2026-07-18). Recipe: the quad gets
`b2 = rho_mix * bY_accel` (full mixture weight as force/volume); the overlay's
`-fluidBody` keeps the acceleration form.

---

## 5. Layered profiles

A layered deposit under one water table is **ONE overlay**. Assembly is per-cell, so
heterogeneous k, porosity, and moduli cost nothing — the storage `n/Kf`, the L
modulus, and the α-stab modulus are already per-element. `Kf` and `rhoF` stay
overlay-global (same water). Each `-layer` block names its element set and its own
`-perm` (required per layer) plus optional `-poro`/`-moduli`; cells not in any
layer use the global values.

```tcl
pattern LadrunoPorousOverlay 77 \
    -region {*}$allEles \
    -Kf 2.2e9 -rhoF 1000.0 \
    -perm 1e-9 1e-9 -poro 0.40 -moduli 1.0e7 0.25 \   ;# global (fallback) props
    -layer {*}$sandEles  -perm 1e-4 1e-4 -poro 0.45 -moduli 3.0e7 0.30 \
    -layer {*}$clayEles  -perm 1e-10 1e-10 -poro 0.55 -moduli 5.0e6 0.35 \
    -drained {*}$surfaceNodes
```

Every `-layer` must set its own `-perm` (fatal otherwise). A layer element outside
`-region`, or claimed by two layers, is a snapshot-time fatal (ADR §11 (P1 panel)
finding P-2). Separate overlays are for **deliberately-separate** water bodies only
(perched zones, isolated aquifers) — never for layers of one connected deposit
(§11).

---

## 6. The iterated driver (P2)

Plain `analyze` gives fs1 (one sweep, O(Δt)). For the **exact monolithic
backward-Euler** answer, drive the fixed-stress iteration to convergence with
`LadrunoStaggeredAnalyze` (ADR §12 P2). At convergence the L modulus cancels and
the fixed point IS single-step backward Euler.

**Tcl**

```tcl
LadrunoStaggeredAnalyze $nSteps $dt -tol 1e-6 -maxIter 500 -pScale $q -verbose
LadrunoStaggeredAnalyze -stats
```

**openseespy**

```python
rc = ops.LadrunoStaggeredAnalyze(nSteps, dt, "-tol", 1e-6, "-maxIter", 500, "-pScale", q)
stats = ops.LadrunoStaggeredAnalyze("-stats")   # [nSteps, totalFluidSolves, meanIters, maxIters, lastResidual, maxResidual]
```

Defaults: `tol = 1e-6`, `maxIter = 500`, `pScale = 1.0`. `-stats` returns the last
run's telemetry (zeros before any run).

- **TRANSIENT-only.** The driver requires an active transient analysis
  (`Newmark`/etc.); a static analysis is a loud fatal (the load-factor trap again,
  §4). Empty driven set (no qualifying overlay) is a loud fatal.
- **Integrator WHITELIST — Newmark / HHT / GeneralizedAlpha only.** The driver
  repeats a same-step `newStep`/`revertToLastStep` pair per iteration; TRBDF2's
  two-step history is not restored (silently degrades to trapezoidal per iteration)
  and the explicit family's `revertToLastStep` is a no-op — both are now a **loud
  fatal** (ADR §12 P2 item 5). Only whitelisted integrators drive.
- Measured iteration counts (classic L, no divergence, no maxIter hit anywhere):
  column **11.47 mean / 19 max**, footing **5.05 mean**; `-fsL oedometric` column
  **4.25 mean** (ADR §12 P2 item 2). Against monolithic BE (the frozen toy
  `qs_mono`), the driver at `-tol 1e-10` matches to **3.44e-7 (p) / 1.1e-8 (u)**
  compressible and **3.39e-7 (p)** near-undrained; against the *real* `LadrunoUP`
  it is a mutual-Δt-convergence leg (observed orders **0.99/0.99**) because no
  stock Newmark parameterization reproduces BE rates on the fluid row (ADR §12 P2
  item 1).

### Stage-flip moduli recipe (PDMY / PM4Sand staging)

`updateMaterialStage` **cannot reach the overlay** — `MaterialStageParameter`
registers only the first accepting element and never scans load patterns (the
ADR-71 sibling-broadcast trap). After you flip a soil constitutive stage, the
overlay keeps its stage-0 moduli, so its L is stale (stable, but
convergence-degraded — L is a preconditioner-like term, not a physics error). The
transport contract: **re-set the overlay moduli explicitly via the parameter
route** after the flip.

```tcl
updateMaterialStage -material $matTag -stage 1        ;# flip the solid material
parameter 9001 loadPattern 77 E                       ;# target the overlay's E
updateParameter 9001 $E_stage1                        ;# re-set -> marks moduli dirty, rebuilds L/stab lazily
```

Parameter targets: `E`, `nu` (overlay-global); `layerE $i`, `layerNu $i` (1-based
`-layer` index). `updateParameter` re-validates (`E > 0`, `0 <= nu < 0.5`). A
BIT-TWIN gate confirms this path is exactly equivalent to constructing with the new
E (ADR §12 P2 item 4). The staged PDMY column tracks the monolithic ADR-71 P4
reference to **rel_p 6.1e-3 / rel_u 5.5e-3** at the same Δt (ADR §12 P2 item 4).

Note `layerNu 0` is **unreachable** via the parameter route — the `Layer` struct
uses `nu <= 0` as an "inherit global" sentinel, so `layerNu 0` is rejected loudly
(the global `nu` accepts 0). See §11.

### openseespy failure surface (split)

`ops.LadrunoStaggeredAnalyze(...)` **raises** `OpenSeesError` for setup misuse
(static analysis active, no transient analysis) but **returns a negative int** (no
exception) for run-time aborts (empty driven set, solve/fluid failure, maxIter:
−1/−2/−3/−6/−7). `assert rc == 0` catches only the second class — guard both. See
§11.

---

## 7. The explicit lane (P3)

For explicit dynamics of a saturated body, run the solid under
`CentralDifferenceLadruno` (matrix-free diagonal — the fork explicit stack stays
untouched, p never enters the integrator) and set the overlay to **`-fsL zero`**.
This is the Zienkiewicz–Paul–Chan 1988 implicit-p/explicit-u split: L ≡ 0 makes the
committed fluid advance degenerate to plain backward-Euler SPD at each commit
(ADR §12 P3, §3.1).

```tcl
pattern LadrunoPorousOverlay 77 -region {*}$eles -Kf 2.2e9 -rhoF 1000.0 \
    -perm 1e-4 1e-4 -poro 0.4 -moduli 2.5e7 0.30 -drained {*}$topNodes \
    -fsL zero -stab auto 0.25
system Diagonal
integrator CentralDifferenceLadruno
analysis Transient
set dtcr [criticalTimeStep]          ;# OVERLAY-AWARE undrained pencil
analyze $nsteps [expr {0.5*$dtcr}]   ;# run at <= 0.5x the pencil
```

```python
ops.integrator("CentralDifferenceLadruno")
ops.analysis("Transient")
dtcr = ops.criticalTimeStep()        # overlay-aware
ops.analyze(nsteps, 0.5 * dtcr)
```

### The overlay-aware critical time step

`-fsL zero` **bypasses the oedometric floor** and prints a one-time loud advisory:
this is the explicit-lane setting; a quasi-static implicit fs1 march with L = 0 is
the naive drained split and **diverges in ~4 steps** at soil coupling (ADR §3.2;
§11). `LadrunoStaggeredAnalyze` **refuses** `-fsL zero` overlays with a loud fatal
(iterating with L = 0 is that same drained split).

`criticalTimeStep()` is **overlay-aware**: it prices the **discrete undrained
pencil**, not the drained one, by folding a per-cell undrained stiffness
augmentation `ΔK_e = Q_e·S*_e⁻¹·Q_eᵀ` into the shared element scan. The report
prints drained-vs-undrained and the implied factor. Run at **≤ 0.5×** the reported
value (the L=0 stability boundary measured exactly 1.000× the undrained pencil at
P0; the 0.5× is the safety margin — ADR §12 P0 E7.2).

The per-element factor **legitimately exceeds** the closed-form material bound
`√(1 + K_f/(n·M_oed))`: measured drained/undrained pencil ratio **26.24 = 21.43 ×
√1.5** on the e72 soft-soil column — the material formula (21.43) is a *lower bound*
on the discrete factor, and the √1.5 excess is a mode-shape effect for the unit Q4
at ν = 0.25 (ADR §12 P3 item 1). The material formula is **documentation-only**
(measured ~1.85× conservative as a Δt predictor, ADR §12 P0 E7.4); the advisory uses
the discrete pencil. The augmentation is **stab-invariant** (measured 1.3e-14): the
per-cell stab Laplacian annihilates the constant-p mode that carries the dominant
coupling, so α-stab cannot throttle the advisory (ADR §12 P3 item 2).

### Zero-drainage secular-pumping advisory (important)

On a **zero-drainage, undamped** column (k-bar ≈ 1e-11, no Rayleigh) there is **no
asymptotically stable Δt at all**. The frozen `+Q·p` force injects O(Δt) splitting
energy per step, and with no dissipation channel it accumulates secularly:
steps-to-blowup measured **48k @ 0.4× pencil, 30k @ 0.5×, 7.8k @ 1.0×, 855 @ 3.0×**
(ADR §12 P3 item 3). A short march looks stable; a long enough one always diverges.
The "L = 0 stable at ≤ 0.5× pencil" pin is **horizon-relative** on this pathology
(0.5× buys ~4× the blow-up horizon of 1.0×, ~35× of 3.0×) — it is not
unconditional stability there. **Any physical drainage or damping absorbs the
input**: the drained control stays bounded, the ZS84 two-leg gate marches clean and
Δt-converges (orders 1.05/1.01), and real soils are fine. Mitigation for a
near-zero-k undamped model: add any small `alphaM` Rayleigh, or use the implicit
lane. See §11.

### SMS + overlay: LUMPED SMS composes (P3b); consistent Olovsson is warned

**Fixed at P3b (ADR §12 P3b item 5).** Both `LadrunoMassScaling` builders now pass
the per-element undrained augmentation through the `elementCriticalDt` Kadd seam,
so SMS sizing prices the **undrained** pencil for overlay-owned cells (a one-time
INFO line — `SMS sizing priced the UNDRAINED pencil for N overlay-owned elements
(ADR-73 P3b)` — replaces the old blanket UNSUPPORTED warning).

- **Lumped `CentralDifferenceSMS` is the supported combination**: the certified
  `dtTarget` march is battery-hard-gated stable (gate (d): 8/8 elements scaled,
  4000-step march bounded; the unscaled control at the same `dtTarget` diverges
  at step 31).
- **Consistent `CentralDifferenceSMSConsistent` under-delivers on overlay cells**
  (measured ~×1.83/step divergence at the certified `dtTarget`): the Olovsson
  centroid-preserving blocks scale only non-rigid element modes, and the undrained
  coupling mode carries a rigid-translation component that stays unscaled. It now
  prints a loud one-time warning when it scales overlay-owned elements — use
  lumped SMS instead, or size dt from the `criticalTimeStep()` report. See §11.

### Subcycling

`-subcycle N` advances the fluid every N solid commits with the accumulated Δu
(the fluid is the slow field — its diffusion CFL is slack by orders on real soils).
`-subcycle auto` picks `N = max(1, floor(0.09·h²/(c_v·Δt)))` from the measured
θ = 0.089 (error doubles at that θ; ADR §12 P0 E7.3a). Error grows ~N^1.2, all N ≤
50 stable on the L=0 lane. Removal events force a sync regardless of N (§10).

### `-fluidUpdate explicit` — the fully matrix-free fluid step (P3b)

Shipped at P3b: the fluid advances by a **lumped-S\* forward step** — no CG, no
factorization anywhere in the transient march (the fluid step is a local axpy).
Recipe:

```tcl
pattern LadrunoPorousOverlay 77 -region {...} -Kf $Kf -rhoF $rf -perm $k $k \
    -poro $n -moduli $E $nu -drained {...} -stab auto 0.25 \
    -fluidUpdate explicit
integrator CentralDifferenceLadruno -cfl
set dtcr [criticalTimeStep]          ;# undrained pencil, min-folded with the
analyze $nsteps [expr {0.5*$dtcr}]   ;#   fluid diffusion limit (dual CFL)
```

What to know (every claim measured — ADR §12 P3b):

- **Pairing/sequencing**: the explicit fluid advances at *load-application*
  time on the trial window Δu (the commit hook only syncs). This is invisible
  to you but load-bearing: a commit-time advance is one force lag behind the
  CD pairing and unconditionally unstable. The shipped lane is the frozen-toy
  E7.4 scheme to machine precision (twin gate 6.7e-14). **CD-family
  integrators are the supported lane** (under an implicit integrator the
  advance pairs with the predictor state — not gated, not recommended).
- **Dual CFL**: the report prints the explicit-fluid diffusion limit
  `h²/(2·ndm·c_v)` and its slack factor (measured 1.5e4×–1.5e8× across
  realistic k̄ — it never governs on real soils); `-subcycle N` must keep
  N·Δt below it. The governing value returned by `criticalTimeStep()` already
  min-folds it.
- **Δt margin**: run at ≤ 0.5× the advisory as usual. The advisory is
  *conservative* for this lane: the measured boundary sits above the pencil
  (toy constant 1.32×; ≥ 2.5× on the battery's all-x-fixed column — the
  per-element pencil is BC-blind, so constrained meshes gain margin).
- **`-fsL` is inert** (no L, no iteration — one-time notice if you pass a
  non-default value); `-stab` is march-inert (row-sum lumping annihilates the
  stab Laplacian) but still shapes the advisory (stab-invariant anyway);
  `-staticMode steady` keeps its small CG solve at static commits (gravity
  init recipes unchanged); `-subcycle`, `-onRemoval`, `-pInit`, `-record`,
  and the P4 recorder channels (`recorder ladruno -overlay`, `Monitor
  -overlay`) work unchanged (CSV-twin gate exact under the lane).
- **Zero-drainage pathology largely dissolves**: the §11 secular-pumping
  legs that kill the implicit `-fsL zero` lane at 48k/30k/7.8k steps are
  bounded to 60k/45k/15k caps under this lane — near-zero-k undamped models
  are better served by `-fluidUpdate explicit`.
- `LadrunoStaggeredAnalyze` **refuses** FU_EXPLICIT overlays (nothing to
  iterate); MP remains per-process serial v1 (the halo design is documented
  in ADR §8, not built).

---

## 8. Why not CentralDifference + u-p elements

You might try running a monolithic u-p element under an explicit integrator instead.
Don't, for the honest-p element:

**Honest-p `LadrunoUP` + upstream `CentralDifference` = Richardson-unstable pore
pressure.** Honest-p carries no p-mass (the fluid row is first-order,
`S·p' + H·p = ...`), so central difference's second-order leapfrog on that row is
the Richardson / DuFort-Frankel explicit scheme for a pure diffusion operator —
**unconditionally unstable**. No Δt, damping, or mass-scaling choice fixes it; it is
structural (ADR §12 P3 item 5; §11).

**The trap:** it **tracks the reference to ~1e-4 first, then dies.** Measured
non-finite at step **~1299** (Δt = 2e-4) / **~209** (Δt = 1e-3). A 400-step
validation march "passes" and lies — the overlay's own P3 run-1 did exactly that
before the incubation was found (ADR §12 P3 item 5). If you must validate this
combination, march thousands of steps.

The overlay explicit lane (§7) is the honest replacement — solid stays diagonal,
fluid stays implicit SPD at commit. **Head-to-head** on the ZS84 column at matched
Δt (ADR §12 P3 item 4, cited — not re-measured here):

| Lane | rel-vs-implicit | cost/step |
|---|---|---|
| Overlay (`-fsL zero`, `CentralDifferenceLadruno`) | 6.1e-3 | 0.074 ms |
| Upstream `CentralDifference` + `FourNodeQuadUP` | 1.33e-2 | 0.181 ms |

The overlay lane measured **2.2× less error at 2.4× less cost per step**. (The
predicted S→0 failure of the incumbent did *not* manifest at K_f×1e3 — reported
honestly; the demo pins the overlay's robustness, not an incumbent failure.)

---

## 9. Recording p

The overlay's p is **not a nodal DOF** — ordinary nodal recorders cannot see it.
Three ways to capture it.

### `-record` CSV (always available)

The `-record $file <$everyN>` flag writes a CSV: a header row of region node order,
then `time, p(node1..nodeN)` per commit. The file opens truncate (a restart
overwrites). This is the gate-facing output and needs nothing else.

### Dedicated recorder channels (P4)

> The `-overlay` recorder channels below are the **pinned P4 surface** (ADR plan
> §4.2/§4.4). They are built in the P4 recorder work; if your build predates it,
> use `-record` CSV above.

`LadrunoRecorder` gains an `-overlay` channel that writes the p-field into the
standard self-describing HDF5 layout
(`RESULTS/ON_NODES/overlayPressure_<tag>/{ID,DATA,STEP,TIME}`), honoring
`-envelope`/`-kind` exactly like other channels (`|p|max` envelope is the
liquefaction use case):

```tcl
recorder ladruno $file -overlay $tag1 $tag2   ;# explicit overlay tags
recorder ladruno $file -overlay               ;# bare = every 33022 overlay (fatal if none)
```

`LadrunoMonitorRecorder` gains a live-tail `-overlay` mode (columns
`overlay<tag>.p.node<n>`):

```tcl
recorder Monitor -overlay $tag <-nodes {$n1 ...}> -sink $file -every $N
;# EXCLUSIVE of -node/-region/-dof/-resp. Optional -nodes: restrict the
;# columns to a subset of REGION nodes (default = all region nodes; a
;# non-region tag in the subset is fatal).
```

An unknown/non-33022 tag after `-overlay`, a bare `-overlay` with no overlay in the
domain, mixing Monitor `-overlay` with `-node`, or a Monitor `-nodes` containing a
non-region node are all parser-time fatals (unknown-flag-FATAL house rule). A
consequence: **create the overlay pattern BEFORE the recorder command** — both
`-overlay` forms validate against the domain at parse time. Recorders
read the freshly-committed pⁿ⁺¹; under `-subcycle N>1` they read the last synced p
between windows (those are the forces the skeleton actually feels).

Further semantics, all by design:

- **Serial-only v1.** The overlay is per-process; under OpenSeesMP a partition
  without the 33022 pattern cannot serve `-overlay` channels (the recorder
  prints a loud error and writes no overlay data on that partition — the
  analysis itself is unaffected).
- **Static stages:** recorder rows appear per static commit and their `TIME`
  axis carries the domain "time" = the **load factor λ** there (the ⟨A-3⟩
  convention, same as every static recorder channel).
- **Pattern removed mid-run** (`remove loadPattern`): the `-overlay` channel
  and the Monitor both keep their frozen column set and stream **0.0** for the
  vanished overlay (advisory live-tail semantics — the run is never stalled).

### Energy accounting note

If you read an ADR-69 `EnergyBalanceRecorder` on an overlay run, the `+Q·p`
coupling work rides the **ULW (external-load-work) channel** — the overlay injects
its forces through `Node::addUnbalancedLoad`, which the ADR-69 kernel's
`ULW = ∫ vᵀP_ext dt` reads. **Closure holds** (measured `|ERR%| = 0.000 ≤ 2` on the
ZS84 overlay run, ADR §12 P3 item 7); only the *attribution* is merged — you cannot
separate coupling work from genuine external-load work in the per-channel numbers.
See §11.

---

## 10. Element removal

The overlay's fluid measure is a **geometry snapshot** — it survives
`remove element`. At each commit the region is rescanned against the domain; newly
dead cells drop their Q coupling automatically (no skeleton to push on), and their
storage/conductance is handled per `-onRemoval`:

- **`-onRemoval keep`** (default) — the dead cell stays in the fluid measure at soil
  k (a water-filled gap). The dead-cell GP set never shrinks the fluid continuum.
- **`-onRemoval drain $kFactor`** — the dead cell's permeability is scaled by
  `$kFactor` (crack-as-drain). Higher kFactor drains faster: measured `Tv90 =
  1.004 / 0.938 / 0.933` for kFactor `{1, 10, 100}`, strictly monotone (ADR §12 P0
  E7.5). The `drain` variant re-factors the (region-sized, cheap) fluid operator.

Under `-subcycle`, a removal event **forces a fluid sync** regardless of N — the
window is caught up so the discontinuous Q change is paired with the right Δu.

The removal transient carries no spurious numerical artifact: the single-commit
jump vanishes ~Δt^2.6 (no 1/Δt impulse); the physical stress redistribution
dominates (ADR §12 P0 E7.5).

---

## 11. Quirks / troubleshooting

Every row traces to a `LEDGER_quirks.md` entry or an ADR §12 log entry.

| Symptom | Cause | Fix / status |
|---|---|---|
| Attaching a `timeSeries` / `loadConst` to the overlay does nothing | The overlay owns its force amplitudes; `applyLoad` ignores `time` and any series (`+Q·p` at full amplitude). The surface cannot even structurally attach a series. | By design; a one-time notice prints. |
| Layered deposit consolidates wrong; no cross-layer flow | Modeled as **one overlay per layer** — each has an independent p-field / drained set. | Use ONE overlay spanning all layers with `-layer` blocks (§5). Element overlap across overlays is a snapshot fatal. |
| Fluid "consolidates" during a static gravity stage; loud march advisory | `-staticMode` default MARCHES; under a static analysis domain "time" is the load factor, so it integrates over Δλ. | Set `-staticMode hold` or `-staticMode steady` for static stages (§4). |
| `-fsL zero` under a quasi-static implicit `analyze` diverges in ~4 steps | L = 0 is the naive drained split (diverges ~10 orders at soil coupling τ ≈ 10³). | Explicit-lane-only setting. Use `-fsL classic\|oedometric` for implicit/driver lanes; reserve `zero` for `CentralDifferenceLadruno` (§7). Parser warns once. |
| `LadrunoStaggeredAnalyze` ignores my `-subcycle N` | The driver syncs the fluid every step by construction. | Intended; one-time advisory prints. A pending window is caught up at driver entry. |
| Overlay L stale after `updateMaterialStage` (slow convergence) | `updateMaterialStage` cannot reach a `LoadPattern` (registers only the first accepting element). | Re-set moduli via `parameter $p loadPattern $tag E $newE` after the flip (§6). |
| `parameter ... layerNu $i` + `updateParameter ... 0` rejected loudly | `Layer` uses `nu <= 0` as an "inherit global" sentinel; a stored 0 would silently mean "unset". | By design; global `nu` accepts 0, per-layer `layerNu` cannot be exactly 0. |
| Staggered twin of a `LadrunoUP` model: solid barely settles (or double-loads) | Copied monolithic `-body` (ACCELERATION) into plain quad `b1 b2` (FORCE/VOLUME); with hydrostatic `+Q·p` the net can cancel to ~zero. | Quad gets `b2 = rho_mix * bY_accel`; overlay `-fluidBody` keeps the acceleration form (§4). |
| Honest-p `LadrunoUP` under `CentralDifference` blows up — but only after ~1300 steps | CD leapfrogs the mass-less (first-order) fluid row = Richardson/DuFort-Frankel diffusion = unconditionally unstable. Tracks the reference (~1e-4) first, then dies. | Structural — no Δt/damping fixes it. Use the overlay explicit lane, or keep `LadrunoUP` implicit. A 400-step validation lies; march thousands (§8). |
| No "pore-coupling work" channel in the energy breakdown; ULW looks inflated | The `+Q·p` forces ride the ADR-69 ULW (external-load-work) channel by construction. | Closure ERR is trustworthy; per-channel attribution of coupling work is not separable (P4 may split it) (§9). |
| SMS certifies a `dtTarget` that then blows up (CONSISTENT builder) | Fixed at P3b for lumped `CentralDifferenceSMS` (Kadd seam wired — sizing prices the undrained pencil; certified march hard-gated stable). `CentralDifferenceSMSConsistent` still under-delivers: the Olovsson centroid-preserving blocks don't scale the coupling mode's rigid-translation component (measured ×1.83/step at the certified dtTarget). | Use lumped SMS with overlays, or size from the `criticalTimeStep()` report × 0.5. The consistent builder warns loudly when it scales overlay cells (§7, §11). |
| `-fsL zero` undamped near-zero-k column diverges no matter how small Δt | No dissipation channel to absorb the frozen-force O(Δt) splitting energy → secular pumping; no asymptotically stable Δt (blow-up horizon 48k @ 0.4× … 855 @ 3.0× pencil). | Add any small `alphaM` Rayleigh, or use the implicit lane. Real soils drain/damp and are fine (§7). |
| `ops.LadrunoStaggeredAnalyze(...)` sometimes raises, sometimes returns a negative int | Setup misuse (static active / no transient analysis) raises `OpenSeesError`; run-time aborts (empty driven set, solve/fluid fail, maxIter) return −1/−2/−3/−6/−7. | Check for both: `assert rc == 0` catches only the run-time class (§6). |
| Deformable transient DB (FileDatastore) restart used to diverge non-deterministically | Was an upstream quad/tri element-`rho` serialization bug (garbage restored mass), amplified by the overlay's large `+Q·p`. | FIXED in #577 (rho serialized + zero-init, all six quad/tri elements); deformable transient DB restart is supported again fork-wide. |

---

### Provenance

All numbers and behaviors above trace to: the shipped parser
(`OPS_LadrunoPorousOverlay.cpp`) and driver (`LadrunoStaggeredDriver.cpp`); ADR-73
§12 implementation-log entries (P0 #575, P1 #576, P2 #580, P3 #581, P3b); the frozen
measured constants block (`_adr73_implementation_plan.md` E7.1–E7.6); and the
`LEDGER_quirks.md` rows dated 2026-07-14 through 2026-07-19. Where a behavior is
unsupported (consistent-SMS+overlay, implicit-integrator explicit-fluid pairing,
deformable ArcLength) or a
pathology has no stable Δt (undamped zero-drainage), the guide says so plainly
rather than softening it.
