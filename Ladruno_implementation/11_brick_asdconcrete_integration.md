---
title: LadrunoBrick × ASDConcrete3D — integration findings & plan
project: Ladruno
tags:
  - plan
  - ladruno-brick
  - asdconcrete
  - element-material-contract
---

# LadrunoBrick × ASDConcrete3D

Driving doc for using [[09_ladruno_brick|LadrunoBrick]] with **ASDConcrete3D**
(Petracca's plastic-damage concrete/masonry model,
`SRC/material/nD/ASDConcrete3DMaterial.{h,cpp}`). ASDConcrete3D **softens**
(tension cracking + compression crushing, stiffness degradation), which changes
what the element owes the material. This records what we verified in the source
and the resulting work list. Companion of [[10_ladruno_j2_plasticity]] (the
non-softening flagship material).

> Status: **(a)+(b)+(c) SHIPPED.** (b) single-point material-eval fix = PR #94;
> (a) Tier-A damage-scaled `Kstab` + (c) Zone-A validation
> (`tests/test_ladrunoBrick_asdconcrete.py`) = this branch. Remaining: the meshed
> Zone-B cases (notched 3-pt bend convergence + hourglass-band monitor). See §6.

## TL;DR

1. **Characteristic-length handshake WORKS for LadrunoBrick out of the box** —
   fracture-energy regularization will function (lch = min edge). Still must be
   *validated* for mesh-objectivity; it is wired, not proven.
2. **Higher-order Ladruno elements (BezierTri6 / BezierTet10) do NOT** — the base
   `getCharacteristicLength()` returns ~half the element size for them (mid-edge
   nodes). They need an override before any fracture-regularized material. See
   [[04_bezier_elements]] / [[06_bezier_tet10]].
3. **Tier-A `Kstab` (damage-scaled hourglass stabilization) is warranted** for the
   softening regime, and ASDConcrete3D makes it *safe* (secant default + IMPLEX +
   exposed damage). Not urgent until reduced/EAS is actually pointed at concrete.
4. **Cost gotcha:** the `eas` path runs 8 redundant ASDConcrete3D return-maps per
   element at the *same* centroid strain — expensive for this material. Fixing it
   is a clear win independent of everything else.

## 1. The characteristic-length (lch) handshake — VERIFIED

ASDConcrete3D regularizes the softening branch so the dissipation per unit crack
area is mesh-objective: it scales the fracture energy by the parent element's
characteristic length, `g_f = G_f / lch`. It obtains `lch` through a global
active-element pointer, **once**, on the first strain push (then frozen):

```cpp
// ASDConcrete3DMaterial.cpp:1613-1620 (inside setTrialStrain)
if (!regularization_done) {           // one-time guard (.h:386)
    if (ops_TheActiveElement)
        lch = ops_TheActiveElement->getCharacteristicLength();
    regularization_done = true;       // FROZEN for the instance's life
    if (auto_regularize) { ht.regularize(lch, lch_ref); hc.regularize(lch, lch_ref); }
}
```

The **one-time-frozen** semantics raise the stakes: whatever element is active at
the *first* `setTrialStrain` sets `lch` permanently. Verified there is no
wrong-element window (see below) — but it means a future refactor that pushes
strain from outside `update()` could silently freeze a stale `lch`.

Traced end-to-end for LadrunoBrick — **all three links hold** (independently
re-verified by an adversarial agent sweep, see §7):

| Link | Where | Result |
|---|---|---|
| Element pushes strain ONLY in `update()` | `LadrunoBrick.cpp:656/691/713/772` (`setTrialStrain`); finite via `setTrialF` 1023. `formResidAndTangent`/`formUri`/`formEAS`/`formPhysical` are READ-ONLY (`getStress`/`getTangent`) | ✅ first material touch can only come from `update()` |
| Domain sets the active pointer before `update()` | `Domain.cpp:2263` `ops_TheActiveElement = theEle;` then `theEle->update()` (adjacent lines, nothing between) | ✅ pointer correct at regularization time |
| `getCharacteristicLength()` returns a sane value | base `Element.cpp:682` → **min node-to-node distance**; **bug-fixed by Petracca 2020** | ✅ ≈ min edge for the 8-node hex |

So `lch = min edge length` for LadrunoBrick, and there is **no window where `lch`
could be frozen from the wrong element** — the dangerous candidate (`formResid…`
called by the analysis FE_Element *outside* `Domain::update()`) does not push
strain, so it cannot trip the one-time guard. Two caveats:
- It is **wired, not validated** — must confirm mesh-objectivity numerically (see §5).
- For a badly **distorted** hex, "min edge" can under-estimate the band width in
  the crack-normal direction. Acceptable for v1; revisit if distorted meshes
  mis-regularize.

## 2. Higher-order elements break the base lch (cross-element finding — ✅ FIXED in ladruno)

> **RESOLVED on `ladruno` (independently, 2026-06-01):** BezierTri6 / BezierTet10
> now override `getCharacteristicLength()` — BezierTri6 returns `sqrt(2·A)`
> (right-isosceles edge), BezierTet10 returns `cbrt(6·V)` (right-tet leg); both
> recover the true edge on right-simplex meshes and err safe (under-estimating)
> on curved Bézier edges. Logged in LEDGER_quirks. The analysis below is retained
> as the rationale; the gap it identified is closed. LadrunoBrick remains exempt.


The base `getCharacteristicLength()` = **minimum** node-to-node distance. That is
correct for **corner-only** elements (LadrunoBrick hex, the 4-node tet), where it
equals the min edge. It is **wrong for elements with mid-edge / interior nodes**:

- **BezierTet10** — nodes 0–3 corners, 4–9 mid-edge (`BezierTet10.cpp:24-26`,
  `edgeV` table 101-108). Min node distance = corner→adjacent-mid-edge.
- **BezierTri6** — nodes 0–2 corners, 3–5 mid-edge (`BezierTri6.cpp:25-30`). Same.

For a **straight-sided, centered** mid-edge node that's ≈ **½ the shortest edge**;
the exact factor is geometry-dependent (curved Bézier edges / off-center control
points shift it), but the *direction* — `lch` under-estimated — is robust.
Neither overrides `getCharacteristicLength()` today, and both consume
`NDMaterial`: **BezierTet10** = `"ThreeDimensional"` → ASDConcrete3D directly;
**BezierTri6** = `"PlaneStress"/"PlaneStrain"` → ASDConcrete3D still reachable via
the `NDMaterial` base wrapper around a 3D copy. A ~2× under-estimate of `lch`
mis-regularizes the fracture energy (`g_f = G_f/lch` doubled) → wrong/too-brittle
softening, mesh-dependent results.

**Fix is element-level** (`getCharacteristicLength()` is virtual on `Element`):
override in each higher-order element to measure over the **corner nodes only**
(e.g. min/representative corner-to-corner edge), not all nodes. LadrunoBrick needs
**no** change. Tracked as a Bezier follow-up, not part of the brick work.

## 3. Tier-A damage-scaled `Kstab` (the softening question)

Background: [[09_ladruno_brick]] `eas`/`uri` use a hourglass/EAS stabilization
stiffness `Kstab` condensed from the **initial elastic** tangent (constant, no
per-step α state — the SSPbrick trick). The 1-pt constant-strain response always
uses the true current material; `Kstab` only stabilizes the hourglass (bending)
modes.

**Problem with softening:** a cracked element keeps *elastic* hourglass-bending
stiffness while its bulk carries ~0 stress → the crack can't localize, stress
locks in the band. The fix (LS-DYNA type-6 / Puso, Theory §3.2 eq 3.28 `2μQ_hg`
with the *current* μ) is to soften the stabilization with the material.

**Two tiers** (from the comparison work):
- **Tier A — scale the elastic `Kstab` by the current degradation.** Cheap, no new
  state. `Kstab ← max(floor, s)·Kstab_elastic`.
- **Tier B — full Simo–Rifai EAS** (rebuild condensation from the consistent
  tangent, commit α per element). Reintroduces the heavy α state AND is *less
  robust* at finite strain (Wriggers–Reese 1996 EAS large-strain hourglassing).
  **Not recommended** speculatively.

**ASDConcrete3D makes Tier A both more-needed and safe:**
- Returns the **secant** operator by **default** (`ASDConcrete3DMaterial.h`:
  "False (default) = use the secant matrix") → positive-semidefinite, monotone
  degradation, **never negative** → the scale degrades cleanly to ~0, no sign flip.
- Supports **IMPLEX** (frozen extrapolated damage over the step → positive operator).
- **Exposes damage directly**: `getAvgDamage()` / `getMaxDamage()` each return a
  **`Vector(2)` = [tension damage dₜ, compression damage d_c]** (`crackingDamage()`
  per branch). So a *scalar* degradation must combine the two — e.g.
  `d = max(dₜ, d_c)`, or a stress-state-weighted blend — rather than a single
  `d_avg`. Cleanest robust choice: `s = 1 − max(dₜ, d_c)`.

**Recommended Tier-A form:** `Kstab ← max(floor, 1 − max(dₜ,d_c)) · Kstab_elastic`,
`floor ≈ 1–5%` so cracked elements keep *some* hourglass control (never fully
unstabilized → no spurious hourglassing in the band). Tune `floor` with the
hourglass-energy report (§5). Design the damage query **generically** (any
material that reports a scalar/Vector damage via `getResponse`) so LadrunoBrick
isn't hard-coupled to ASDConcrete3D — if a material doesn't report damage, fall
back to the current constant elastic `Kstab`.

> **SHIPPED (a), 2026-06-02.** Tier A is implemented exactly as recommended:
> `Kstab ← max(floor, 1−max(dₜ,d_c))·Kstab_elastic`, `floor = 1%`
> (`HG_DAMAGE_FLOOR` in `LadrunoBrick.cpp`), in `formEAS` + `formUri` STIFFNESS.
> Damage is read generically via a cached `materialPointers[0]->setResponse(
> "damage", …)` Response (built in `setDomain`); a material with no `"damage"`
> channel falls back to the constant elastic `Kstab` (`s=1`) — so elastic/J2 are
> bit-unchanged. NB for `uri+stiffness` the elastic base is the *initial* shear
> (not the current secant), otherwise the secant would drive `κ→0` and the floor
> would be meaningless. For *implicit* quasi-static concrete, `std`/`bbar` (full
> integration) still sidesteps the whole `Kstab` question — 8 independent damage
> points, no hourglass — at 8× material cost.

## 4. Cost gotcha — 8 redundant return-maps per element (eas) — ✅ FIXED

> **FIXED (this branch, 2026-06-02).** `LadrunoBrick::update()` now evaluates the
> constitutive model **once** (material slot 0) for the single-point formulations;
> slots 1-7 are output mirrors. Added `isSinglePoint()`; per-GP output paths
> (`getResponse` 3/4/6/7, the `material`/`integrPoint` `setResponse`+`setParameter`
> channels, `displaySelf`) map GP queries to slot 0. Test
> `tests/test_ladrunoBrick_singlepoint_output.py` (6 gates incl. a guard that the
> mirror does NOT leak into std/bbar). Full LadrunoBrick suite 90/90.
> **Consequence:** `eas`/`uri-stiffness`/`uri-viscous` now genuinely cost **1**
> material eval/element — the cost premise in the handout below is updated
> accordingly. (Memory note: 7 phantom material instances per element are still
> *allocated* and committed — cheap state-rolls, no return maps — a separate
> memory-only optimization if it ever matters.)

The original (now-fixed) behavior — single-point paths fanned the centroid strain
to all 8 material copies; the multi-point paths do genuine per-GP work:

| Formulation | # `setTrialStrain` | strain | redundant? |
|---|---|---|---|
| `eas` | 8 | centroid (identical) | **yes** (`LadrunoBrick.cpp:651-656`) |
| `uri` + `stiffness` | 8 | centroid (identical) | **yes** (`:698-713`) |
| `uri` + `viscous` | 8 | centroid (identical) | **yes** (same branch) |
| `std` | 8 | per-GP distinct | no (real 2×2×2) |
| `bbar` | 8 | per-GP distinct | no |
| `uri` + `physical` | 8 | per-GP distinct | no (assumed strain at 2×2×2) |
| `-geom finite` | 8 × `setTrialF` | per-GP distinct F | no |

With a cheap elastic material the 8 identical evals are free; with **ASDConcrete3D
+ IMPLEX** (one of the most expensive materials in OpenSees) it means **8 identical
return-maps per element** — throwing away most of reduced integration's cost
advantage, which was the reason to pick `eas`/`uri` for explicit concrete.

**Fix (done):** the single-point paths now evaluate slot 0 once and mirror it in
the output paths (lazy populate in `getResponse`), as above.

## 5. Validation plan (and the hourglass report is the instrument)

The hourglass-energy / viscous-dissipation report (PR #86, `"hourglassEnergy"` /
`"hgDissipation"`) is exactly the dial for tuning the Tier-A `floor`:

1. **Single-element uniaxial tension** → stress-strain follows the *regularized*
   softening backbone; dissipated energy `= G_f/lch × V = G_f × crack area`.
   Proves the lch handshake end-to-end.
2. **Notched 3-point bend (mode-I)**, 2–3 mesh refinements → load–CMOD curves must
   **converge** (fracture-energy objectivity). The acceptance test for §1/§2.
3. **Watch `hourglassEnergy` in the cracked band**: spikes ⇒ `floor` too low
   (spurious hourglassing); crack won't localize / over-stiff ⇒ `floor` too high.
   Empirical floor-setting loop.

## 6. Scoped next steps (decide a/b/c)

- ~~**(b) Fix the 8× redundant ASDConcrete3D eval**~~ **DONE (2026-06-02)** —
  single-point paths do 1 eval; output mirrored from slot 0;
  `tests/test_ladrunoBrick_singlepoint_output.py`; suite 90/90.
- ~~**(a) Tier-A `max(floor, 1−max(dₜ,d_c))`-scaled `Kstab`**~~ **DONE (this branch,
  2026-06-02).** Implemented in `formEAS` (scales `easKstab` in both stiffness and
  the internal stabilization force) **and** `formUri` STIFFNESS branch (for a
  softening material rebases `κ` on the *initial* elastic shear, then ×`s`, so the
  floor actually protects — using the current secant shear would let `κ→0` at full
  damage and defeat the floor). `floor = 1%` (`HG_DAMAGE_FLOOR`). Damage read
  **generically** via a cached `materialPointers[0]->setResponse("damage", …)`
  Response built in `setDomain` (rebuilt on the receive side; **no NDMaterial base
  change**); `damageScale()` returns `max(floor, 1−max(d_i))`, or **1.0** when the
  material has no `"damage"` channel ⇒ elastic/J2 behave **exactly** as before
  (192 Zone-A tests, incl. the full LadrunoBrick suite, unchanged). uri+viscous is
  damping (not stiffness) so it is deliberately **not** degraded. The
  `"hourglassEnergy"` report carries the same scale (the §5 tuning instrument).
- ~~**(c) Concrete validation tests**~~ **DONE (Zone-A, `tests/test_ladrunoBrick_asdconcrete.py`).**
  (1) **lch handshake / mesh-objectivity** — a cube pulled to failure at L and 2L
  dissipates energy in the ratio **≈4** (∝ crack AREA = `G_f·A`), decisively not 8
  (∝ volume); at L=lch_ref the absolute dissipation equals the backbone specific
  energy × volume. Proves the regularization consumes the element `lch`.
  (2) **Tier-A** — for an identical prescribed deformation the stabilization-energy
  ratio between an ASDConcrete3D element and an elastic one (same E,ν) equals
  **exactly** `max(floor, 1−max(dₜ,d_c))` at every step, from intact (≈1) down to
  the floor (0.01) — never zero. **Remaining for Zone-B (needs apeGmsh meshing):**
  the **notched 3-point bend** with 2–3 mesh refinements (load–CMOD convergence)
  and the **hourglass-energy-fraction-in-a-cracked-band** monitor (`<5–10%` of
  internal energy) — the multi-element localization study was too solver-fragile
  to ship as a deterministic Zone-A test (snap-back / non-convergence at fine
  meshes), so it is deferred to a meshed Zone-B case.

Cross-element follow-up: ~~override `getCharacteristicLength()` in BezierTri6 /
BezierTet10~~ **DONE on ladruno** (§2) — both now supply a corner-node element
size, so they can drive ASDConcrete3D too.

## 7. Verification record (adversarial agent sweep)

Every source claim in §§1–4 was independently re-verified by three skeptical
read-only agents (instructed to *refute*), each reading the actual source. All
**CONFIRMED**; refinements already folded above:

- **lch handshake (§1) — CONFIRMED.** Regularization is **one-time, guarded by
  `regularization_done`** (`ASDConcrete3DMaterial.h:386`, set in `setTrialStrain`
  `.cpp:1616`). LadrunoBrick's *only* strain-push site is `update()`
  (`setTrialStrain` at 656/691/713/772; `setTrialF` 1023); all assembly routines
  are read-only. `Domain::update()` sets `ops_TheActiveElement` on the line above
  `theEle->update()`. **No wrong-element window.**
- **Bezier override needed (§2) — CONFIRMED, and since FIXED on ladruno.** No
  override or alternate lch path existed; node ordering and the `½-edge` direction
  confirmed; LadrunoBrick correctly exempt. ladruno has since added the overrides
  (BezierTri6 `sqrt(2·A)`, BezierTet10 `cbrt(6·V)`).
- **ASDConcrete3D secant-default + IMPLEX + damage exposure (§3) — CONFIRMED.**
  `tangent=false` default in `.h`, parser (`.cpp:409`, only `-tangent` flips it),
  and constructor (`.cpp:1559`); `getTangent()` returns the damaged secant
  `C = W·D0` by default (`.cpp:2461-2469`); IMPLEX present; `getAvgDamage()` /
  `getMaxDamage()` return `Vector(2)=[dₜ,d_c]`.
- **8× redundant eval (§4) — CONFIRMED**, and it covers `eas` **and** `uri`
  stiffness/viscous (table above).

---

# USER HANDOUT — using LadrunoBrick for concrete (ASDConcrete3D)

> Practical, consume-later reference. Synthesized from a 3-agent source audit
> (element technology, ASDConcrete3D settings, command-level usage). For the
> *why*, see §§1–4; this section is the *how*.

## H.0 Material-eval cost per formulation (fix (b) landed)

Per element, per state determination (the constitutive return-map count — the
dominant cost for ASDConcrete3D):

| Formulation | Material evals | Damage points |
|---|---|---|
| `eas`, `uri+stiffness`, `uri+viscous` | **1** (centroid) | 1 (smeared over the element) |
| `std`, `bbar`, `uri+physical` | 8 (per-GP) | 8 (independent) |

The single-point formulations **now genuinely save ~8× material cost** (fix (b),
§4) — relevant because ASDConcrete3D + IMPLEX is expensive. The tradeoff is
**fidelity**: a single-point element resolves damage at ONE point, and — until
Tier-A `Kstab` (a) lands — stabilizes its hourglass modes with a *constant
elastic* `Kstab` that over-stiffens cracked zones. So the choice is now a real
**cost vs. fidelity** decision, not cost-parity.

## H.1 Formulation decision guide

| Analysis | Use | Why |
|---|---|---|
| **Implicit quasi-static** (push-over, monotonic/cyclic capacity) | **`bbar`** (default), `std` if ν low & no bending | Full 2×2×2 → 8 independent damage points, **no hourglass modes**, best-conditioned tangent for Newton+arc-length, sidesteps the un-implemented damage-scaled `Kstab` entirely. `bbar` also cures volumetric locking as ν→0.5 / dilatant crushing. Costs 8 material evals — buy the robustness; drop to `eas` (1 eval) only if that cost bites *and* Tier-A `Kstab` has landed. |
| **Implicit, near-incompressible (ν→0.5) or coarse bending pre-peak** | **`eas`** | Cures shear+volumetric across all ν; constant initial-tangent `Kstab` keeps the tangent well-conditioned. Verify on a coarse crushing case (constant `Kstab` can mildly over-stiffen a fully-cracked element). |
| **Explicit dynamic** (impact/blast/fast crushing, or violent-snapback fallback) | **`eas`** (cheapest: 1 eval, cures locking) if you can accept 1 damage point + monitor hourglass; **`uri`+`physical`** (8 evals) when you need genuine per-GP damage + material-degrading stabilization | Real cost/fidelity split now (§H.0): `eas` is 1 material eval but smears damage to one point and uses a constant elastic `Kstab` (Tier-A pending); `physical` pays 8 evals for 8 damage points and stabilization that degrades with the material. Both need `-lumped`. |
| **Avoid for softening** | `uri` + `stiffness` | Elastic `Kstab` does NOT degrade → over-stiffens cracked elements, props up the load, corrupts the crack-band energy. Only if the element never fully damages, and *monitor*. |
| **Explicit-only, never implicit/eigen** | `uri` + `viscous` | Rank-deficient tangent (no hourglass stiffness) — singular under any static/eigen step. |

**One-liner:** *Implicit → `bbar` (or `eas` for ν→0.5). Explicit → `uri+physical` (or `eas`). Never `uri+stiffness`/`viscous` under softening-implicit.*

## H.2 Command syntax + worked examples

```python
element('LadrunoBrick', tag, n1..n8, matTag,
        '-formulation', <std|bbar|uri|eas>,            # default std
        '-hourglass',   <stiffness|physical|viscous> [, coeff],  # uri only; coeff default 0.05
        '-lumped',                                      # diagonal mass (required for explicit)
        '-b', bx, by, bz,                               # body force / unit volume
        '-damp', dampTag,                               # std/bbar ONLY (dropped+warned for uri/eas)
        '-geom', <linear|corot|finite>)                 # default linear; corot is std/bbar only
```

Implicit quasi-static concrete (recommended):
```python
ops.nDMaterial('ASDConcrete3D', 1, *concreteArgs)
ops.element('LadrunoBrick', 101, 1,2,3,4,5,6,7,8, 1, '-formulation', 'bbar')
```

Explicit dynamic concrete:
```python
ops.nDMaterial('ASDConcrete3D', 1, *concreteArgs)
ops.element('LadrunoBrick', 201, 1,2,3,4,5,6,7,8, 1,
            '-formulation', 'uri', '-hourglass', 'physical', '-lumped')
ops.integrator('CentralDifferenceLadruno')      # or 'ExplicitBathe'
dt = 0.9 * ops.criticalTimeStep()               # global op, NOT an element method
```

> `-geom corot` (large rotation / small strain) is **std/bbar only** (uri/eas
> guarded off) — usable with a small-strain concrete material **if** strains stay
> small. `-geom finite` needs a `FiniteStrainNDMaterial`; ASDConcrete3D is
> small-strain, so for **finite strain** concrete is not yet supported here.

## H.3 ASDConcrete3D robust recipe (material side)

```python
ops.nDMaterial('ASDConcrete3D', tag, E, v, '-rho', rho,
    '-Te', *Te, '-Ts', *Ts, '-Td', *Td,           # tension backbone (strain, stress, damage)
    '-Ce', *Ce, '-Cs', *Cs, '-Cd', *Cd,           # compression backbone
    '-autoRegularization', lch_ref,               # MUST pass for mesh-objective softening
    '-implex', '-implexControl', 0.05, 0.01,      # robust integration + error-driven step cuts
    '-Kc', 0.6667)
```

| Setting | Code default | Recommended | Why |
|---|---|---|---|
| secant vs `-tangent` | **secant** | keep secant | Positive, monotone-degrading operator — robust in softening. (`-tangent` is ignored under IMPLEX.) |
| `-autoRegularization $lch_ref` | **OFF** (factory) | **ON** | Regularizes `G_f` by the element `lch` (= LadrunoBrick min edge, §1). Without it softening is **mesh-dependent**. `lch_ref` = the size your `G_f` was calibrated for. |
| `-implex` `-implexControl tol tred` | OFF | ON, `0.05 0.01` | IMPLEX gives a positive operator through localization; control drives step cuts (returns `-10`). Pair with an adaptive step. Without control a too-large step drifts silently. |
| `-eta` | `0.0` | small (∼ fraction of `dt`) for dynamic | Duvaut–Lions viscous regularization of softening; default rate-independent. |
| `-Kc` | `0.6667` (2/3) | keep | Lubliner shape factor; valid range [2/3, 1]. |
| `-cdf` | `0.0` | raise→1 only if you want tension/compression damage coupling | cross-damage factor. |

**Record at minimum:** `damage` (`[d⁺ d⁻]`), `crackWidth`, and (under IMPLEX) `implexError`.
Per-GP from the element: `eleResponse(tag, 'material', gp, 'damage')`.

## H.4 Solver guidance

- **Implicit post-peak: arc-length or local displacement/CMOD control is mandatory.** Steep concrete softening produces snap-back that plain load/displacement control **cannot traverse**. Control the crack-opening DOF when global arc-length struggles on a single sharp band.
- **Explicit central-difference** sails through negative-definite tangents — use it for impact/blast, *and* as a fallback when softening snap-back simply won't converge implicitly (load slowly / mass-scale; keep kinetic ≪ internal energy).

## H.5 Mesh guidance

- **Crack-band upper bound on element size:** `lch < E·G_f / f_t²` (Hillerborg length). Too-large elements make the regularized softening slope exceed the elastic stiffness → **spurious GP-level snap-back** (numerical, not physical). ASDConcrete3D silently floors the energy if `lch` is too large (`gmin`), quietly flattening softening — so size elements sanely. The single-point formulations use the *whole-element* `lch` (larger → hits the bound sooner) — another reason coarse `uri` meshes are touchy in softening.
- **Keep cracking-zone hexes compact and near-cubic** (aspect ratio ≲ 3). Distortion degrades every formulation and the `uri`/`eas` stabilization assumes near-parallelepiped geometry; push distortion into the elastic far-field.
- Resolve the band with **a few elements across** where affordable — `lch` keeps the *energy* objective even at one-element bands, but path/orientation are far better resolved.

## H.6 The hourglass-in-softening hazard + monitoring (for any `uri`, and `eas` under heavy crushing)

The danger is structural: reduced integration's zero-energy modes are normally resisted by the deviatoric/stabilization stress — **the very stress damage drives to zero**. So a cracked element either over-stiffens (elastic `stiffness` `Kstab`) or lets the hourglass mode grow into the cracked bulk as a fake, mesh-dependent mechanism. **This is the whole reason to prefer full integration for softening.**

Monitor every `uri` run:
```python
ops.eleResponse(tag, 'hourglassEnergy')   # Vector(1)
```
| Formulation | What it reports |
|---|---|
| `uri+stiffness` / `eas` | **stored** stabilization energy (instantaneous) |
| `uri+viscous` | **cumulative dissipated** energy (committed accumulator) |
| `std` / `bbar` / `uri+physical` | `0` |

**Acceptance:** hourglass energy **< ~5%** of internal energy (hard ceiling **~10%**); watch the ratio *over time* (a late spike = hourglass feeding on the cracking event). Also confirm dissipation localizes into **one physically-sensible crack band**, not a **checkerboard/keystone** pattern. Above ceiling or checkerboard → **fall back to `bbar`/`std`** for that region.

## H.7 Limitations today (read before trusting `uri`/`eas` for concrete)

1. ~~**`Kstab` is constant elastic, not damage-scaled**~~ **FIXED (a)** — `eas` and `uri+stiffness` now degrade `Kstab` by `max(floor, 1−max(dₜ,d_c))` (floor 1%), so cracked zones soften and can localize while keeping residual hourglass control; the `"hourglassEnergy"` report reflects the degraded value. `uri+viscous` is damping (unchanged). Materials without a `"damage"` channel are unaffected. *Full integration (`bbar`/`std`) remains the most robust choice for implicit softening (8 independent damage points, no hourglass at all).*
2. ~~**8× redundant material eval**~~ **FIXED (b)** — single-point formulations now do 1 material eval (§4); slots 1–7 are output mirrors. (7 phantom material *instances* still allocated — memory-only optimization remains.)
3. **Geometry for concrete:** `-geom corot` (large-rotation/small-strain, std/bbar only, shipped on ladruno PR #88) works with a small-strain concrete material as long as strains stay small; **finite strain** concrete is unsupported (`-geom finite` needs a `FiniteStrainNDMaterial`, which ASDConcrete3D is not).

## H.8 Starting recipes (copy-paste)

```python
# ---- IMPLICIT quasi-static concrete (recommended baseline) ----
ops.nDMaterial('ASDConcrete3D', 1, E, v, '-rho', rho,
               '-Te',*Te,'-Ts',*Ts,'-Td',*Td, '-Ce',*Ce,'-Cs',*Cs,'-Cd',*Cd,
               '-autoRegularization', lch_ref, '-implex', '-implexControl', 0.05, 0.01)
ops.element('LadrunoBrick', eTag, *nodes, 1, '-formulation', 'bbar')
# ... numberer/system/constraints, then arc-length or displacement control for post-peak

# ---- EXPLICIT dynamic concrete ----
ops.element('LadrunoBrick', eTag, *nodes, 1,
            '-formulation', 'uri', '-hourglass', 'physical', '-lumped')
ops.integrator('CentralDifferenceLadruno')
dt = 0.9 * ops.criticalTimeStep()
# ... monitor eleResponse(eTag,'hourglassEnergy') < 5-10% of internal energy
```
