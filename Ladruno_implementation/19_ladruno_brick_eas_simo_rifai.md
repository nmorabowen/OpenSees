---
title: LadrunoBrick — true Simo–Rifai EAS (and the ssp/eas rename)
project: Ladruno
status: implemented
priority: medium
owner: nmora
tags:
  - implementation
  - element
  - eas
  - locking
---

# LadrunoBrick — true Simo–Rifai EAS (and the `ssp`/`eas` rename)

> Authored from a 10-agent planning + adversarial sweep (5 design dimensions ×
> {plan, refute}), grounded in de Souza Neto §15.2, Belytschko Ch 8,
> Ibrahimbegovic Ch 2, the in-tree `EnhancedQuad` (2-D 4-mode Simo–Rifai), and
> the live `LadrunoBrick` source. The adversarial pass found one **blocker in the
> theory** (a wrong natural→physical map) and several citation/robustness fixes;
> all are folded in below and called out in **[FIX]** notes so the corrections are
> auditable.

## What

Two coupled changes to `LadrunoBrick` (classTag 33002, one element, formulation
axis):

1. **Rename** the existing `-formulation eas` → `-formulation ssp`. Today's "eas"
   is *not* EAS: it is a verbatim port of UWelements `SSPbrick` — a **Stabilized
   Single-Point** element (one centroid Gauss point; nine "enhanced" modes
   condensed **once** at `setDomain` against the *initial* elastic tangent `C(0)`
   into a **constant** `easKstab`/`easBnot`; no per-step `α` state, no
   per-iteration condensation). See `buildEAS`/`formEAS`
   ([LadrunoBrick.cpp:1896](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
   [:2374](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)). The name is borrowed;
   `ssp` is the honest name.

2. **Add** a *true* small-strain **Simo–Rifai Enhanced Assumed Strain** element
   under the freed `-formulation eas`: full 2×2×2 integration, 9 enhanced
   parameters `α` as committed element-internal state, the internal condition
   `∫ Mᵀσ dV = 0` driven to zero by a **per-element Newton inside every global
   iteration**, and static condensation of `α`.

**Scope (v1):** small strain, geometry = linear, **E9** enhanced-mode set.
**Not in scope (designed-for, deferred):** richer mode sets (E12/E21/E30) and the
finite-strain (enhanced deformation gradient) EAS. Finite EAS is a separate
follow-up because Q1/E9 is documented to **hourglass under finite compressive
strain** (de Souza Neto §15.2.5, Fig 15.10(a)).

## Why

- `bbar` cures volumetric locking but not bending (transverse-normal/Poisson)
  locking; `std` locks both ways on coarse meshes; `ssp` is a stabilization, not
  a variationally consistent mixed element, and freezes its enhancement against
  `C(0)` so it cannot track a nonlinear tangent. A true EAS hex is the
  industry-standard coarse-mesh bending/incompressibility cure (Abaqus C3D8I
  lineage) and completes the LadrunoBrick formulation matrix with a
  *consistent* mixed element.
- It is the natural consumer of the finite-strain seam already built for
  `[[project_ladruno_brick]]` (`-geom finite`, `FiniteStrainNDMaterial::setTrialF`),
  giving a future enhanced-`F` large-strain element with no analysis-core change.

## Where

- **Modify (rename `eas`→`ssp`):**
  - `SRC/element/ladrunoBrick/LadrunoBrick.h` — enum `Formulation`
    ([:72](../SRC/element/ladrunoBrick/LadrunoBrick.h)); class doc strings
    ([:32](../SRC/element/ladrunoBrick/LadrunoBrick.h),
    [:38](../SRC/element/ladrunoBrick/LadrunoBrick.h)).
  - `SRC/element/ladrunoBrick/LadrunoBrick.cpp` — `formulationName`
    ([:114](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)); the `// eas …` doc at
    [:29](../SRC/element/ladrunoBrick/LadrunoBrick.cpp); every dispatch on
    `Formulation::EAS` ([:249](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
    [:465](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
    [:748](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
    [:760](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
    [:909](../SRC/element/ladrunoBrick/LadrunoBrick.cpp),
    [:2792](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)).
  - `SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp` — parser
    ([:109](../SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp)), usage strings
    ([:33](../SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp),
    [:66](../SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp)), and the
    finite/corot reserved-combination guards
    ([:210](../SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp),
    [:252](../SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp)).
- **New code (true EAS):** new methods on `LadrunoBrick` — `formEAStrue`,
  `computeBenh` (the 6×9 `M` operator), the shared `condenseEAS` helper, `α`
  state. No new file, no new classTag.
- **Reference (copy patterns from):**
  `SRC/element/fourNodeQuad/EnhancedQuad.{h,cpp}` — the proven 2-D Simo–Rifai:
  `formResidAndTangent` inner-`α` Newton + condensation
  ([EnhancedQuad.cpp:1028](../SRC/element/fourNodeQuad/EnhancedQuad.cpp)),
  `computeBenhanced` ([:1392](../SRC/element/fourNodeQuad/EnhancedQuad.cpp)),
  `getInitialStiff` ([:1276](../SRC/element/fourNodeQuad/EnhancedQuad.cpp)).
- **Ledgers / banner / tests:** `Ladruno_implementation/LEDGER_implementations.md`,
  `Ladruno_scripts/banner_features.txt` (+ `patch_banner.py`), the 9 test files
  under `tests/` that hard-code `-formulation eas` (see §6).

---

## How — the design (5 dimensions)

### 1. Theory — variational basis and the enhanced-strain interpolation

Voigt convention throughout = the element's existing **{xx, yy, zz, xy, yz, zx},
engineering shear** (`LadrunoBrick::computeB`,
[:1383](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)).

**[FIX — citation faithfulness]** de Souza Neto §15.2 (eqs 15.39–15.44, 15.57–15.58,
15.82, 15.85) is the **finite-strain enhanced-*gradient*** principle: independent
fields `(u, H̃, P)`, energy `Π̃ = ∫[ρψ(C̃) + P:H̃]dV`, `F̃ = F + H̃`, `H̃`
generally **unsymmetric**. What v1 implements is its **small-strain symmetric
reduction**, whose primary citation is **Simo & Rifai (1990)**; the book equation
numbers are kept only where they transfer exactly (the E9 mode *set*, eq 15.85).

#### 1.1 Functional and additive split (Simo–Rifai 1990)

Three-field Hu–Washizu with independent `(u, ε, σ)`:

```
Π(u, ε, σ) = ∫_Ω [ W(ε) − σ : (ε − ∇ˢu) ] dV − Π_ext(u)
```

EAS re-parametrises the assumed strain by an **additive enhanced split**:

```
ε = ∇ˢu + ε̃ = B d + M(ξ) α                                        (E.1)
```

`ε̃ = M α` is the enhanced (incompatible) strain; `α ∈ ℝ⁹` are element-internal.

#### 1.2 Stationarity (the three conditions)

```
δσ :  ∫_Ω δσ : ε̃ dV = 0                          (L2-orthogonality)   (E.2)
δu :  ∫_Ω Bᵀ σ dV = f_ext                         (global equilibrium) (E.3)
δα :  ∫_Ω Mᵀ σ dV = 0   per element               (internal/enhanced)  (E.4)
```

with `σ = ∂W/∂ε` at the **total** strain `ε = Bd + Mα`. Because `M` is
interelement-discontinuous, (E.4) is element-local and `α` is condensed out.

#### 1.3 Patch test = mean-zero condition

**[FIX — citation]** Patch test (constant-stress reproduction) holds iff each
enhanced mode is L2-orthogonal to the assumed-stress space; for the Simo–Rifai
constant-per-element stress choice this reduces to

```
∫_Ω M(ξ) dV = 0      (each of the 9 columns integrates to zero)        (E.5)
```

Primary citation **Simo & Rifai (1990)**; the general orthogonality is de Souza
Neto eq (15.58), and (E.5) is its small-strain, constant-stress specialization.

#### 1.4 Discrete system + condensation

```
K_dd = ∫ Bᵀ C B dV,  K_dα = ∫ Bᵀ C M dV,  K_αα = ∫ Mᵀ C M dV,  h = ∫ Mᵀ σ dV
K̄ = K_dd − K_dα K_αα⁻¹ K_αd,   r̄ = f_int − K_dα K_αα⁻¹ h               (E.6)
```

For a symmetric material tangent `C = Cᵀ` (hyperelastic / associative /
incremental-potential), `K_αα` and `K̄` are **symmetric** (de Souza Neto
Remark 15.7). **The variational distinction from `ssp`:** there `K_αα` is
condensed once against `C(0)` and `α` is never solved; **true EAS drives (E.4) to
zero by an inner Newton at every global iteration** — exactly `EnhancedQuad`'s
`do … while(‖residE‖>tol)` loop ([:1028](../SRC/element/fourNodeQuad/EnhancedQuad.cpp)).

#### 1.5 The E9 enhanced modes (de Souza Neto eq 15.85)

The nine modes are the single-entry **gradient** tensors `E_i` (eq 15.85); v1
takes their symmetric part. Each of the three **normal** strains is enriched by
each natural coordinate — **no enhanced shear in v1**:

```
            α1     α2     α3     α4     α5     α6     α7     α8     α9
          ┌                                                            ┐
   ε̃_xx  │  ξ      0      0      η      0      0      ζ      0      0   │
   ε̃_yy  │  0      ξ      0      0      η      0      0      ζ      0   │
M_ξ =     │  0      0      ξ      0      0      η      0      0      ζ   │   (E.7)
   ε̃_zz  │  ...(see legend)...                                        │
   ε̃_xy  │  0      0      0      0      0      0      0      0      0   │
   ε̃_yz  │  0      0      0      0      0      0      0      0      0   │
   ε̃_zx  │  0      0      0      0      0      0      0      0      0   │
          └                                                            ┘
```

**[FIX — unambiguous legend]** column→tensor map (eq 15.85), order-locked (not
cosmetic): `α1=ξ·e_xx, α2=ξ·e_yy, α3=ξ·e_zz, α4=η·e_xx, α5=η·e_yy, α6=η·e_zz,
α7=ζ·e_xx, α8=ζ·e_yy, α9=ζ·e_zz`. (The partial sketch that earlier preceded this
matrix is deleted — it had a stray ambiguity.)

**Why E9:** it is the minimal set that gives each normal strain a *linear*
variation in the natural coordinates — precisely the field the trilinear `Q1`
cannot represent and the source of volumetric + bending locking. It has an exact
textbook anchor (eq 15.85) and its 2-D parent Q1/E4 *is* the installed,
validated `EnhancedQuad`. It deliberately omits the enhanced-shear modes (the
three zero rows) — those are E12/E21/E30 (§7) and carry the finite-strain
hourglass risk.

**[FIX — softened stability claim]** Q1/E9 passes the constant-stress patch test
in small strain and is proven stable/convergent in the *infinitesimal* theory
(Reddy & Simo 1995); it is documented to **hourglass in finite strain** (de Souza
Neto §15.2.5, Fig 15.10(a), the Q1/E9 3-D brick) — hence small-strain v1 and the
E12+ deferral gate.

#### 1.6 Natural→physical map — **[BLOCKER FIX, authoritative]**

The map is the **one-sided enhanced-gradient transform** (de Souza Neto eq 15.82,
small-strain symmetric reduction), which is **exactly what `EnhancedQuad`
implements**:

```
M_i(ξ) = sym[ (j0 / j(ξ)) · J0⁻ᵀ · E_i(ξ) ]   in Voigt                 (E.8)
```

- `j0 = det J0` at the **centroid** ξ=0;  `j(ξ) = det J(ξ)` at the Gauss point;
- `J0⁻ᵀ` = inverse-transpose of the **centroid** Jacobian, applied **one-sidedly**
  to the single-entry gradient tensor `E_i` (e.g. `E_1 = ξ·e_xx`), then `sym(·)`.

> **Do NOT use the two-leg strain (Voigt) transform `(j0/j)·T0⁻ᵀ·M_ξ`.** The
> original plan proposed it and claimed equivalence with the gradient transform;
> the adversarial pass proved this **false and self-contradictory**: for a sheared
> `J0`, `sym(J0⁻ᵀ E_1) ≠ T0⁻ᵀ{1,0,0,0,0,0}` (cross terms differ), and
> `EnhancedQuad::computeBenhanced` ([:1409–1444](../SRC/element/fourNodeQuad/EnhancedQuad.cpp))
> demonstrably builds the **one-sided** map (a *column* of `J0⁻ᵀ` scaled by `L/j`,
> laid into the B-slots) — the Wilson incompatible-mode lineage. Both maps pass
> the patch test (both are centroid-constant, so the odd-function cancellation in
> (E.5) works regardless), so the **patch test cannot distinguish them** — only one
> is Simo–Armero Q1/E9, and that is (E.8). A 6×6 `T0` is *not* shipped in v1; if it
> is ever revived for enhanced-shear modes it must be derived from first
> principles and unit-tested against a brute-force `A εₙ Aᵀ` before being trusted.

**Why (E.8) gives `∫M dV=0` exactly:** change variables `dV = j(ξ)dξ`,

```
∫_V M_i dV = ∫_□ (j0/j) sym(J0⁻ᵀ E_i) · j dξ = j0 · sym(J0⁻ᵀ ∫_□ E_i dξ) = 0   (E.9)
```

because (i) `1/j` cancels the variable Jacobian, leaving the **constant**
centroid operator `j0 J0⁻ᵀ` outside the integral, and (ii) every `E_i` is a single
**odd** natural coordinate, so `∫_□ E_i dξ = 0`. The **centroid-only** evaluation
of `j0`, `J0` is exactly what lets them factor out; per-GP evaluation would
destroy the cancellation and fail the patch test on distorted meshes.

### 2. Element architecture (mirror `EnhancedQuad`)

**Enum & factory.** `enum class Formulation { STD, BBAR, URI, SSP, EAS }` —
ordinals `{0,1,2,3,4}`. **`SSP` deliberately takes ordinal 3** (the old `EAS`
slot) so serialized streams stay compatible (§6.2). Parser: `"ssp"` → `SSP`;
`"eas"` → the new `EAS`.

**Integration points.**
- `isSinglePoint()` ([:741](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)):
  `if (formulation == Formulation::SSP) return true;` — true EAS **falls through to
  `return false`** (8 live material points; this single flip routes `update()`,
  the per-GP output mirrors, and damage averaging to real per-GP data). **[FIX]**
  add a unit assertion that `-formulation eas` reports 8 live GPs, so a future
  mechanical rename cannot silently flip it back.
- `damageScale()`/`easKstab`/`damageResponse` are **SSP-only**. True EAS has **no
  artificial stabilization** — localization comes from the 8 live points; do not
  wire the Tier-A damage lever into the EAS path.

**State (the only new persistent data):** `Vector alpha(9)`, `Vector alphaCommit(9)`.
- `commitState`: `alphaCommit = alpha`.
- `revertToLastCommit`: `alpha = alphaCommit` (**[FIX]** — EnhancedQuad has a
  latent missing-revert here; we fix it).
- `revertToStart`: zero both.
- `update()`/`getTangentStiff()`/`getResistingForce()` for the EAS path call
  `formEAStrue`.

**`formEAStrue(tang_flag, useInitialTangent)`** — faithful 3-D port of
`EnhancedQuad::formResidAndTangent`:
1. Loop 8 GPs, build `M_i(ξ)` per (E.8), set strain `ε = Bd + Mα`, pull `σ, C`.
2. **Inner Newton on `α`** (`d` fixed) until `‖residE‖ < tol`. **[FIX — pin sign
   convention to match the reference verbatim]:** `residE = −∫Mᵀσ dV`;
   `K_αα.Solve(residE, dalpha)`; `alpha += dalpha`. (Drop the alternative
   `+h / −dalpha` variant entirely so a line-by-line diff vs `EnhancedQuad` stays
   trivial.)
3. Condense per (E.6): `K̄ = K_dd − K_dα K_αα⁻¹ K_αd`,
   `r̄ = ∫Bᵀσ − K_dα K_αα⁻¹ h`.
4. Route `K̄, r̄` through the seam-3 `globalizeStiff`/`globalizeForce` exactly like
   the other formulations.

**`getInitialStiff` (EAS):** condense `getInitialTangent` at `α=0` with **no inner
Newton**, mirroring `EnhancedQuad::getInitialStiff`
([:1276–1283](../SRC/element/fourNodeQuad/EnhancedQuad.cpp)). **[FIX]** do not
assert `h=0`; either document the stress-free-reference assumption (`σ(0)=0`,
consistent with how `std` builds `BᵀD_init B`) or run one condensation pass at
`α=0` without the claim.

### 3. Finite-strain seam (designed-for, **deferred**)

Reuse the existing finite path: `isFinite()`, `theGeom` (`SolidTransformation`),
`deformationGradient()`, `FiniteStrainNDMaterial::setTrialF(F)`
([:1093](../SRC/element/ladrunoBrick/LadrunoBrick.cpp), .h
[:224](../SRC/element/ladrunoBrick/LadrunoBrick.h)).

- **Enhanced deformation gradient** (Simo–Armero–Taylor; de Souza Neto eq
  15.89/15.90): `F̃ = F + H̃`, with the modes mapped as
  `H_i(ξ) = (j0/j(ξ)) · F0 · J0⁻ᵀ E_i(ξ) J0⁻¹`. At small strain
  `F̃ → I + ∇u + H̃` and `sym H̃` collapses to **exactly** the additive enhanced
  strain `ε̃ = Σ αᵢ sym(J0⁻ᵀ E_i)` of §1.6 — the seam is consistent by
  construction.
- **Shared, geometry-agnostic `condenseEAS()`:** the block reduction
  `K̄ = K_uu − K_uα K_αα⁻¹ K_αu`, `r̄ = r_u − K_uα K_αα⁻¹ h` operates on
  already-assembled 24×24 / 9×9 blocks, so small- and finite-EAS call it
  verbatim. **What must be abstracted now:** (a) `α` state + commit cycle +
  serialization, (b) `condenseEAS()`, (c) a `computeBenh`/`computeGenh` that the
  finite path overrides for the `H_i` gradient map. **What waits:** the enhanced-`F`
  kinematics and the material driver call.
- **[FIX — solver symmetry]** finite-EAS `K̄_T` is **symmetric** for
  hyperelastic/incremental-potential materials (Remark 15.7) and unsymmetric only
  for non-potential dissipative ones. This is **independent** of the multiplicative
  F-bar unsymmetry (eq 15.10). The element should query the material/formulation
  to decide solver symmetry, not assume "already unsymmetric."
- **[FIX — no false hook]** the SSP `Kstab` is a constant-`C(0)` condensation lever
  and `damageScale()` is a softening-concrete degrade; **neither is the
  Glaser–Armero per-step artificial hourglass control** that finite EAS needs to
  tame the **compressive hourglassing** (Fig 15.10(a)). The code has a *place* to
  add such a term, not the term itself — do not claim otherwise.
- **Non-goals v1:** any `-geom finite`/`corot` + `eas` combo (the parser keeps
  rejecting them, §6); enhanced-shear modes; hourglass stabilization.

### 4. Validation (the merge gate)

Built `opensees.pyd` on CPython 3.12; gmsh 4.15.2 available; Zone-A/Zone-B testbed.

1. **Constant-stress patch test** — distorted (non-parallelepiped) hex patch,
   constant-strain BCs. Must reproduce constant stress to **machine precision**
   (the direct check of `∫M dV=0`). This is the test that catches a wrong `M(ξ)`.
2. **EnhancedQuad oracle** — **[BLOCKER NOTE]** `EnhancedQuad` is the *only* in-tree
   oracle for the **per-iteration condensation algorithm** (commit/iterate/condense).
   Match a thin-extruded 3-D EAS against 2-D plane-strain `EnhancedQuad` on the
   `{xx,yy,xy}` sub-block to **1e-6**.
3. **Bending beats `std`/`bbar`** — slender single-layer cantilever / pure bending;
   reference = Euler–Bernoulli tip deflection. EAS must approach it where `std`
   locks.
4. **Cook's membrane** (extruded 1-hex-thick, `E=1, ν=1/3`): tip `u_C ≈ 23.9`,
   accept `u_C(fine) ∈ [23.6, 24.2]`.
5. **Near-incompressible** `ν → 0.499, 0.4999`: no volumetric locking vs reference.
6. **Reduce-to-`std`** on a regular cube under uniform strain: converged `α → 0`.
   **[FIX — realistic tolerances]** `‖α‖` to the **inner-Newton tolerance**
   (~1e-8…1e-10, one number reused across tests 1/6/7 and tied to the C++ inner
   tol), and `|u_eas − u_std| ≤ a few × Newton_tol` (**not** 1e-12). Reserve 1e-12
   for a *white-box* check that the condensation correction `K_dα K_αα⁻¹ K_αd` is
   `< 1e-10·‖K‖` on the regular cube.
7. **State cycle** — **[FIX]** use a **path-dependent** material (`LadrunoJ2`) on a
   **distorted** hex under bending/near-incompressible load so `α` is genuinely
   large; assert `commit`/`revert`/`revertToStart` of `α` to ~1e-9, **and** assert
   `‖α‖ ≫ tol` and `α₁ ≠ α₂` (else the test is vacuous). `sendSelf`/`recvSelf`
   round-trip of `alphaCommit`.
8. **White-box `eleResponse`** — new C++ `setResponse` branches `"alpha"` (the 9
   parameters) and `"intMdV"` (Σ `M(ξ_gp) w_gp detJ_gp` over 8 GPs, asserted
   `< 1e-12·V`).

### 5. (covered above — finite seam)

### 6. The rename migration (`eas` → `ssp`)

#### 6.1 Touch points (exhaustive — a miss is a silent regression)
- **Enum/dispatch:** `LadrunoBrick.h:72`; all `Formulation::EAS` dispatch sites in
  `.cpp` (listed in *Where*); `formulationName` `.cpp:114`.
- **Doc strings that name the selector** (will *lie* after the rename): `.h:32`,
  `.h:38`; `.cpp:29`; **[FIX]** also the `.cpp` header-banner usage copies at
  `.cpp:32` and `OPS_LadrunoBrick.cpp:33`/`:66` — separate copies the first pass
  missed.
- **Parser + guards:** `OPS_LadrunoBrick.cpp:109` (`"eas"`); the reserved-combo
  guards at `:210`/`:215` and `:252`/`:257` retext `uri/eas` → `uri/ssp`.
- **Tests (9 files) flip `-formulation eas` → `ssp`** and *keep* asserting
  single-point semantics: `test_ladrunoBrick_singlepoint_output.py` (mirrors slot
  0), `_hourglass_energy.py` (Kstab energy), `_asdconcrete.py` + `_asdconcrete_bend.py`
  (damageScale/Kstab), `_finite.py:286` (rejects `-geom finite`), `_bending.py`,
  `_element.py`, `_nonlinear.py`, `_recorder_loads.py`. **New** `eas` tests (§4) are
  separate, fresh files for true EAS.
- **Banner:** edit `Ladruno_scripts/banner_features.txt` (do **not** hand-edit the C
  strings), then `python Ladruno_scripts/patch_banner.py`, then rebuild. Two lines:
  one for `ssp` (stabilized single-point hex) and one for `eas` (true Simo–Rifai
  EAS). Each `shipped` ledger row needs a matching banner line.
- **Ledger:** update the `LadrunoBrick eas` row in `LEDGER_implementations.md` to
  `ssp` + correct its description ("enhanced assumed strain" → "stabilized
  single-point"); add a new row for the true-EAS `eas`. Honesty pass on every
  "enhanced assumed strain" comment that actually describes the SSP port.

#### 6.2 Back-compat — **[the elegant part]**
**[FIX]** `idData(28)` packs `formulation + 10·massType + 100·hgType +
1000·geomMethodID` ([LadrunoBrick.cpp:2508](../SRC/element/ladrunoBrick/LadrunoBrick.cpp)),
so the formulation **ordinal** is what travels in DB/parallel streams. Ordering
`{STD=0,BBAR=1,URI=2,SSP=3,EAS=4}` makes an **old `eas` stream (ordinal 3)
deserialize as `SSP`** — i.e. the single-point element it actually was. New `eas`
models write ordinal 4. **No reorder of STD/BBAR/URI; classTag stays 33002**
(formulation is an internal axis — no `classTags.h`/broker change).

**`α` serialization — [FIX, the architecture blocker]:** do **not** reuse
`idData(25)` (it is the existing Rayleigh-damping flag, 0/1). Instead send a
**guarded extra `Vector(9)` carrying `alphaCommit`** (the *committed* state — DB
sends happen post-commit, so this is the right field) **only when
`formulation == EAS`**; `recvSelf` decodes the formulation from `idData(28)`
first, conditionally receives the vector, then sets `alpha = alphaCommit =
received`. Non-EAS streams are byte-identical to today.

#### 6.3 Deprecation policy for the reused token

> [!question]
> `-formulation eas` now selects a **different element** than before. Two options:
> **(A) immediate semantic switch** — `eas` = true EAS now; emit a **one-time**
> info notice the first time it is constructed ("`eas` now selects true
> Simo–Rifai EAS; the legacy stabilized single-point element is `-formulation
> ssp`"). **(B) hard-error window** — `-formulation eas` produces **no element** +
> a migration warning for one release, forcing an explicit `ssp`-or-`eas` choice.
> **Recommendation: (A)** for a fast-moving research fork (true EAS is a strict
> upgrade; serialized models already reload correctly as `ssp` via §6.2). If (B)
> is chosen, **[FIX]** add a regression test asserting `-formulation eas` yields no
> element + the warning, and note the deprecation branch short-circuits *before*
> the finite/corot guards.

## Risks / open questions

> [!question]
> **E9 vs richer sets for incompressibility.** E9 passes the patch test but its
> near-incompressible performance vs `ssp`/`bbar` is unproven. If E9 underperforms
> on volumetric locking, hold the E12/E21/E30 deferral gate — do **not** jump to
> E21/E30 and inherit hourglass instability. Decide E12 first (cheapest, adds the 3
> enhanced-shear modes) only if a shear-locking benchmark demands it.

> [!question]
> **Inner-Newton robustness.** Iterating `α` per global step changes convergence
> vs the frozen `ssp` condensation. Softening/non-associative materials make `K_αα`
> unsymmetric (Remark 15.7) and may need an unsymmetric inner solve; pick and
> document the inner tolerance (the single number reused in validation §4).

- **Map authority:** (E.8) one-sided gradient transform is authoritative; the
  two-leg `T0` form is **not shipped** and is unvalidated for shear. Any future
  enhanced-shear work must re-derive and unit-test the operator first.
- **Finite seam latency:** keep `H̃ = M_grad α` (gradient, possibly unsymmetric)
  re-exposable at the seam; do **not** bake `sym(M)α` so deep that the finite
  enhanced-gradient map cannot be restored.
- **Numerical:** symmetric solver is correct for the symmetric small-strain `K̄`;
  the unsymmetric concern is a finite-strain / dissipative-material follow-up.

## Implementation log

- Sequencing: **PR-1 = rename `eas`→`ssp`** (mechanical + enum-ordinal back-compat
  §6.2 + retarget the 9 tests + banner/ledger) — lands green with zero behavior
  change. **PR-2 = true EAS under `eas`** (state + `formEAStrue` + condensation +
  `computeMenh` per (E.8) + §4 validation). Splitting keeps the rename bisectable
  and the new physics isolated.

### PR-1 — DONE (2026-06-02, merged PR #150)

Renamed `eas`→`ssp` across source/tests/banner/ledger/guide; enum `SSP=3` (old EAS
ordinal) for serialized back-compat; `-formulation eas` errored with a migration
hint (interim). Verified: build clean + 118 brick tests green. (Merge `ad19d65f3`.)

### PR-2 — DONE (2026-06-03)

True Simo–Rifai EAS shipped under `-formulation eas` (enum `EAS=4`). Built and
green: `tests/test_ladrunoBrick_eas.py` **6/6** + full brick regression **124/124**.

What landed, vs the plan above:
- **Mode set / map:** implemented the authoritative **one-sided** map
  `M_i = sym[(j0/j)·J0⁻ᵀ·E_i]` (§1.6 / E.8) as the **incompatible-modes (Wilson)
  form** the 2-D `EnhancedQuad` actually uses — **E9 = 3 natural-direction bubbles
  × 3 dofs** in `computeMenh` (6×9). The two-leg `T0` strain transform is NOT
  shipped (the blocker the adversarial sweep caught). `easJ0inv`/`easJ0det` cached
  in `buildEAStrue` (setDomain).
- **Solve + condense:** `formEAStrue` runs the inner Newton on α
  (`residE = −∫Mᵀσ`, `Kaa.Solve`, `alpha += dalpha`, ≥2 iters, tol 1e-10), then
  `K* = Kdd − Kda·Kaa⁻¹·Kad`; the inner Newton drives `h→0` so `f* = ∫Bᵀσ` with no
  explicit residual-condensation term (the `EnhancedQuad` contract — the plan's
  `r̄ = r_u − K_uα K_αα⁻¹ h` term is exactly the dropped `h≈0`). Condensation reuses
  the inner-loop `Kaa` (Solve preserves the matrix, as `EnhancedQuad` relies on).
  `getInitialStiff` condenses the initial elastic tangent at α=0, no Newton.
- **Deviation from plan — `update()` is NOT a no-op.** The plan (mirroring
  `EnhancedQuad`) had `update()` do nothing and let the form pass strain the
  materials. That commits **stale (zero) material state under `algorithm Linear`**,
  where no force/tangent form runs at the final `u` before `commitState` (the
  8-GP stress test caught all-zero stresses). Fix: `update()` for EAS calls
  `formEAStrue(0,…)` to solve α + strain the 8 GPs — the standard `update()`
  contract — so commit is correct under **any** algorithm. (`EnhancedQuad` has the
  latent Linear-staleness; we don't.)
- **State / serialization:** `alpha`/`alphaCommit(9)` exactly per §2/§6.2 — guarded
  extra `Vector(9)` of `alphaCommit` sent only for `formulation==EAS`, gated on the
  decoded ordinal in `recvSelf`; `idData(25)` (Rayleigh flag) untouched.
- **Deprecation policy:** PR-1's interim hard-error on `-formulation eas` is now
  replaced by EAS itself — `eas` builds true EAS (the §6.3 option-A end state). The
  PR-1 deprecation test was flipped to assert `eas` builds.
- **Validation (§4) — what was built vs deferred:** patch test (distorted 2×2×2,
  free interior, `∫M dV=0`), reduce-to-`std` (α→0), bending-beats-`std`→Euler,
  near-incompressible non-locking (nu=0.45 bending), 8-live-GP, α commit/revert with
  J2 — all green. **Deferred (noted in the test file):** Cook's membrane band, the
  `EnhancedQuad` 2-D plane-strain sub-block oracle (§4 item 2), and the white-box
  `intMdV`/`alpha` `eleResponse` branches (no C++ response added — the distorted
  patch test is the operative `∫M dV=0` gate).
- **Test-harness lessons (not element bugs):** `ops.fix`+`ops.sp` on the same dof =
  two competing Penalty terms → half the imposed value (use `sp` only); a free
  uniaxial single cube is statically determinate so volumetric locking does NOT
  manifest there (use bending at high ν).

### PR-2 follow-up — EAS vs ssp/bbar/SSPbrick on Lemaitre ductile damage (2026-06-03)

Compared `eas` against `ssp`, `bbar`, and upstream `SSPbrick` on the Lemaitre
double-edge-notched (DEN) bar (`test_lemaitre_notched_bar.py` geometry; steel
E=2e5, J2 + `-damage lemaitre` + `-autoRegularization`). Findings:

- **`ssp` ≡ upstream `SSPbrick` under damage too.** Peak load 32055 N vs 32081 N
  (0.08%), dissipated energy 44615 vs 44635 N·mm — our port reproduces SSPbrick on
  a path-dependent + softening problem, not just elastically.
- **`bbar` is the well-resolved reference:** 8 damage points capture the ligament
  localization (damage 0.13–0.31 across the element); `ssp`'s **single centroid
  damage point under-resolves it** (0.05–0.07) and over-stiffens (peak 32 kN vs
  bbar's 24 kN). Both are mesh-objective (W within 2–5% coarse↔fine).
- **`eas` STALLS on the notched bar** (reaches only u≈0.17/1.5 coarse, 0.09 fine —
  just past yield onset at the notch, damage≈0).

  **Mechanism — investigated, NOT elastic hourglassing (an earlier draft of this
  note said "small-strain EAS hourglassing"; that is wrong and is corrected here):**
  - *Rank test:* a free `eas` element has **exactly 6 zero-energy modes** (rigid
    body) and its first deformation eigenvalue (1.2821e4) **matches `std`/`bbar`
    exactly** — `eas` is rank-sufficient with **no spurious elastic mode**, as
    small-strain EAS theory guarantees (Reddy & Simo 1995). So there is no elastic
    hourglass.
  - *Homogeneous nonlinearity is fine:* a single steel element (uniaxial, same
    J2+damage) reaches full load with damage 0.48 — identical to `ssp`/`bbar`
    (`test_eas_drives_j2_lemaitre_damage`). So plasticity+damage per se is not the
    problem.
  - *The failure is specific to NON-HOMOGENEOUS (high strain-gradient) PLASTIC
    states* (the notch). The inner-Newton on `α` is fragile there (the spurious
    warnings were concentrated at the notch). A **backtracking line search** on the
    inner Newton was tested: it helped (reach 0.169→0.231, inner warnings 89k→28k)
    but did **NOT** clear the stall, at ~13× unit-suite cost — so it was reverted.
  - **Conclusion:** a genuine **instability of the enhanced modes under
    non-homogeneous inelasticity** — EAS's stability guarantee is for the *elastic*
    operator; under plastic/softening tangents in high-gradient regions the enhanced
    sub-problem destabilizes. Inner-solve robustification only partially mitigates
    it; the real cure is a dedicated **EAS stabilization scheme** — scoped in
    **ADR 20** ([[20_ladruno_brick_eas_stabilization]]).

  > [!warning] **RETRACTED (2026-06-03 — see [[20_ladruno_brick_eas_stabilization]]).**
  > This "genuine instability of the enhanced modes" conclusion is **wrong**. The
  > ADR 20 DEN-bar sweep showed bare `eas` traverses the notched bar to full
  > elongation under any reasonable solver (small enough step), at every mesh/
  > regularization. The observed stall was the **inner-Newton absolute-tolerance bug**
  > (next bullet — fixed *this same PR*) plus **coarse stepping** (a Newton
  > basin/solver issue), **not** an enhanced-mode instability. The rank/eigen evidence
  > above already showed the element is rank-sufficient with no spurious mode; the
  > "destabilizes under inelastic tangents" leap was unsupported. A scalar `β·Kαα⁰`
  > stabilization was implemented and **refuted** (it has no stall to cure and
  > degrades convergence). Use `eas` freely on these problems with a normal adaptive
  > solver.

- **Bug the comparison caught + fixed (this PR):** the inner-Newton convergence test
  used an **absolute** tol `1e-10` on `‖∫Mᵀσ‖` — scale-dependent (units force×length²),
  unreachable at steel scale ⇒ 300k+ spurious "did not converge" warnings. Replaced
  with a **relative** criterion `‖∫Mᵀσ‖ ≤ tolRel·r0 + tolAbs` (tolRel 1e-8, r0 = the
  first-iteration residual; condense against the converged `Kaa`). The PR-2 unit
  tests only passed before because they used E≈1e3.

**Guidance (UPDATED 2026-06-03 — supersedes the earlier "eas not robust on notched"
claim):** `eas` is the most accurate member of the family on bending /
near-incompressible response and is **also robust on notched / localization-dominated
inelasticity** with a normal adaptive solver (step-cut + line search) — the DEN-bar
sweep in ADR 20 confirmed it reaches full elongation where the earlier note claimed it
stalled. `ssp`/`bbar` remain fine choices (and `ssp` is cheaper, single-point); but
there is **no robustness reason to avoid `eas`** here. The previously-"scoped" EAS
stabilization (ADR 20) was implemented and **rejected** — scalar `β·Kαα⁰` has no
stall to cure and degrades convergence; see [[20_ladruno_brick_eas_stabilization]].

**Still deferred (future PRs):** enhanced-`F` finite EAS (the §3 seam; Q1/E9
compressive hourglassing — the *one* place EAS genuinely needs a deformation-dependent
Reese–Wriggers stabilization); richer mode sets E12/E21/E30 (§7); the §4 oracle/Cook
refinements and white-box responses. *(EAS inelastic-localization stabilization is no
longer a deferred item — it was tried and rejected, ADR 20.)*
