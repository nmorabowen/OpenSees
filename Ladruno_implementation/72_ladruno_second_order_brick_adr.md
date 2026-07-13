---
title: "ADR 72 — Second-order bricks: LadrunoBrick20 (20-node serendipity hex) + the quadratic-hex strategy"
status: draft — planning, NO code
---

# ADR 72 — Second-order bricks on the Ladruno solid family

**Status:** draft. Planning only — no code yet. Fills the oldest open slot in the
element roadmap: the [[README]] "Ladruno brick element(s) — higher-order hex"
placeholder row (there since the Bézier days) and the matching
[[ladruno_apegmsh_contract]] heads-up row. 3D sibling of the ADR-70 quadratic
triangle (`LadrunoLST`, T6): 2D now has Q4/T3/T6; 3D has H8
([[09_ladruno_brick|LadrunoBrick]]) and the quadratic Bézier simplex
([[06_bezier_tet10|BezierTet10]]) — the quadratic hex is the missing member.

Companion theory loaded for this ADR: `/fem-mechanics-expert` (Felippa AFEM
Ch. 18, Cook Ch. 6/9/13, Hughes Ch. 4/7/9, Belytschko Ch. 8, Oñate Ch. 8,
Bathe Ch. 5), `/abaqus-theory` element technology (TG §3.2.1–3.2.6; AUG
§28.1.1/§28.1.4/§39.1.2 — C3D20/C3D20R/C3D27 family), `/ls-dyna` (Theory §1.1,
§3.7, §22.1/§34; Vol I `*SECTION_SOLID` ELFORM 23/24, `*ELEMENT_SOLID_H20`).
All manual/book citations below were re-verified against the sources in this
ADR's research pass (2026-07-10); items a source does not print are flagged
**[computed]**.

---

## 1. Context & problem

The fork's 3D continuum story is first-order: `LadrunoBrick` (33002) carries
every anti-locking formulation (`std|bbar|uri|ssp|eas`), the `-geom finite`
axis, explicit hardening (bulk viscosity, hourglass energy, HRZ/SMS), and the
coupling/host APIs — but it is trilinear. The only second-order hex anywhere in
the build is upstream `Twenty_Node_Brick` (Lu & Yang 2004, `ELE_TAG 49`), a
2004-era element with real gaps:

- **Consistent mass only** (`getMass` = 27-pt consistent; no lumped path), and
  a dense 60×60 Rayleigh `getDamp` assembled per call.
- **27-pt full integration only.** No reduced-integration option — 3.4× the
  material evaluations of the 2×2×2 rule that Abaqus calls "the most
  cost-effective element in Abaqus" for smooth problems (TG §3.2.4; AUG
  §28.1.1), and full integration **volumetric-locks**: C3D20 is on Abaqus'
  explicit lock list for near-incompressible Mises/Hill response (AUG §28.1.1
  "Fully incompressible material behavior").
- No `-lumped`/HRZ, no `-damp` (Damping objects), no `getCharacteristicLength`
  (crack-band lch), no `-geom` axis, no element-contract awareness, legacy
  static shared buffers, and the old u-p-era code style.

Why the fork wants a second-order hex at all (the demand, ranked):

1. **Stress accuracy per DOF on smooth problems.** Quadratic solids converge
   O(h³) in displacement / O(h²) in energy vs h²/h¹ for H8 (Belytschko §8.2.3),
   capture stress concentrations and curved geometry with far fewer elements
   (TG §3.2.1), and a single-layer H20 mesh reproduces thick-beam bending
   exactly where H8 needs stacks (Oñate §8.11, Fig. 8.23). Foundations,
   D-regions, thick members, machine-part-like details — the implicit statics
   workhorse.
2. **Family completeness.** ADR-70 gives the 2D family its quadratic member
   (T6). gmsh/apeGmsh produce second-order hex meshes natively; the contract
   doc has carried the "higher-order hex" heads-up row since 2026-05.
3. **Eigen / implicit dynamics accuracy** (modal SSI, frequency-domain work):
   quadratic elements + consistent mass give 4th-order-class eigen accuracy
   (Hughes §7.3.2).

And the honest counter-demand this ADR must respect (all three sources agree):

- **Explicit dynamics is NOT the target.** Abaqus/Explicit excludes
  second-order hexes entirely (only modified C3D10M tets; AUG §28.1.1);
  LS-DYNA dropped 16/20-node solids in **1979** for "5× cost, the zero-energy
  modes of the reduced 8-point rule, the high frequency content which drove the
  time step size down, and numerical noise from ad hoc mass lumping" (LS-DYNA
  Theory §1.1; Vol I Intro), and still tags its re-added quadratic solids
  (ELFORM 24/25/26) "under development" in R16. The fork's explicit-quadratic
  slot is already owned by the Bézier line (positive Bernstein lumping) — this
  element only needs a *correct* (positive-mass) explicit path, not a fast one.
- **Softening/localization is NOT the target.** For C⁰ constitutive response
  (elastoplastic), convergence caps at h²/h¹ — "there appears to be no benefit
  in going to higher-order elements for nonsmooth materials" (Belytschko
  §8.2.4); Abaqus recommends first-order elements for large-strain plasticity /
  limit loads (TG §3.2.4). Crack-band softening stays on the H8 family.

---

## 2. Decision (summary)

1. **New element `LadrunoBrick20`** (`ELE_TAG_LadrunoBrick20 = 33018` —
   33017 is claimed by the same-day ADR-71 `LadrunoUP` draft; reserve at P0 per
   the ADR-70/#540 precedent). 20-node
   serendipity hex, 60 DOF, one `NDMaterial` per GP, in
   `SRC/element/ladrunoBrick/` beside `LadrunoBrick`. A **separate class**, not
   a node-count switch inside `LadrunoBrick` — the 8-node element is shipped,
   heavily validated, serialization-frozen; the house pattern is siblings
   (Quad/CST/LST), not variable-topology classes.

2. **`-formulation {std|uri}`** (single selector, family vocabulary):
   - `std` — full **3×3×3 (27-pt)** Gauss displacement element. **Reduce-to
     anchor:** reproduces upstream `Twenty_Node_Brick` stiffness/resisting
     force/consistent mass to ~1e-12 on a distorted hex (the same headline gate
     `LadrunoBrick std` ran against `Brick`).
   - `uri` — **uniform 2×2×2 (8-pt)** reduced integration, the C3D20R analog
     and the production configuration for smooth problems (~30% of the
     std assembly cost, and stresses sampled at the Barlow points). **No
     hourglass control, by design** — Abaqus applies none to C3D20R (TG §3.2.4)
     — but the 6 spurious modes are *gated and documented*, not hand-waved
     (§3.2, P2 gates). `-hourglass` on this element is a **clear parser error**
     (asymmetry with `LadrunoBrick -formulation uri` is deliberate: H8@1-pt has
     12 communicable hourglass modes that MUST be stabilized; H20@8-pt has 6
     modes that are non-communicable in solid meshes).
   - Reserved keywords (parser error "not implemented", API stable): `bbar`
     (Hughes §4.5.2 quadratic B-bar / selective 2×2×2 volumetric — see §3.4 for
     why H20+bbar is a half-measure and H27+bbar is the real vehicle) and `i14`
     (the Irons 14-pt rank-sufficient middle rule, LS-DYNA Theory §3.7's
     choice — see §7).

3. **Mass = formulation-independent full integration.** Consistent mass always
   integrates ∫ρNᵀN with the 27-pt rule (matching Abaqus: "the mass matrix and
   distributed loads use full integration" even on reduced elements, AUG
   §28.1.1). **`-lumped` = HRZ diagonal scaling via the shared
   [[35_ladruno_hrz_lumped_mass_adr|LadrunoMassLumping.h]]** — NEVER row-sum:
   row-sum lumping of the H20 gives **negative corner masses (−M/8 each)**
   [computed; the 2D Q8 analog −M/12 is book-printed in Hughes §7.3.2 /
   Cook Fig. 13.3-4]. HRZ on the cube gives corners 7/248, mid-edges 2/31 of
   total mass [computed, method validated against Cook's printed Q8 values] —
   all positive, and the P0 oracle pins those exact fractions.

4. **Geometry: v1 `-geom linear` only.** `-geom finite` is a later phase that
   should ride a **shared 3D finite-strain kernel** extracted from
   `LadrunoBrick::formResidAndTangentFinite` (the ADR-70 lesson: don't inline
   the hardest mechanics twice). No corot.

5. **Node & GP ordering = upstream `Twenty_Node_Brick` order** (corners 1–8,
   then bottom / top / vertical mid-edges — the same layout as Abaqus C3D20 and
   LS-DYNA `*ELEMENT_SOLID_H20` Fig. 19-31; exact edge sequence re-verified
   against `brcshl` at P0). This makes the reduce-to gate per-GP aligned and
   reuses the recorder's already-tabulated rules: `Hexahedron_GaussLegendre_3`
   (27-pt, brcshl order — tabulated for `Twenty_Node_Brick` in
   `Ladruno_ElementResults.h`) for `std`, and the 8-node Brick's
   `Hexahedron_GaussLegendre_2` for `uri`. gmsh's hex20 ordering differs →
   apeGmsh carries the permutation (contract doc row).

6. **27-node H27 sibling: deliberately NOT in v1, slot kept warm** (33019 is
   the natural tag). The theory scorecard actually *favors* H27 for three
   futures — explicit (row-sum = HRZ = Lobatto lumping, all positive:
   1/216 corners, 4/216 mid-edges, 16/216 mid-faces, 64/216 center — "the only
   [lumping] … recommended" combination, Cook pp. 374–376), distorted meshes
   (Oñate §8.11: H27 "generally preferable … on distorted geometries"), and
   near-incompressibility (selective-2×2×2 constraint ratio r = 3 = optimal vs
   H20's 1.5) — but each of those demands is hypothetical today, H27 costs 35%
   more DOF, and its own reduced form is the one with genuinely dangerous
   modes (C3D27R: "three unconstrained, propagating hourglass modes", AUG
   §28.1.1). Decide by demand; see §7.

7. **Out of scope v1** (each with a §7 note): quadratic contact/tie faces (the
   documented corner-force sign-flip pathology), embedded-rebar/node host APIs,
   u-p (upstream `Twenty_Eight_Node_BrickUP` exists), bulk-viscosity flags,
   cubic elements (never: "little is gained by going beyond the second-order
   elements", TG §3.2.1).

---

## 3. Formulation (grounded in theory)

### 3.1 Shape functions & completeness

Serendipity H20 (natural coords ξ,η,ζ; corners at ±1, mid-edges at 0):

- corners: `N_i = ⅛(1+ξξ_i)(1+ηη_i)(1+ζζ_i)(ξξ_i+ηη_i+ζζ_i−2)`
- mid-edges (ξ_i = 0 class): `N_i = ¼(1−ξ²)(1+ηη_i)(1+ζζ_i)` + cyclic.

(Felippa AFEM §18.3 Eqs. 18.4–18.7; Oñate Eqs. 8.56a,b.) Complete to
**quadratic** (all 10 cubic-space monomials of degree ≤ 2) plus exactly 10
incomplete higher terms; the full triquadratic (27-monomial) basis needs H27
(Oñate p. 274 — which is why H20 "saves 21 nodal variables per element" and is
the industry serendipity default). Both pass the linear patch test on arbitrary
(straight-edge) geometry; **quadratic completeness in physical coordinates
degrades with edge curvature and off-center midside nodes**, and H20 degrades
faster than H27 (Cook pp. 180, 187–188; Oñate §8.11 p. 296). Gates use
straight-edged meshes; distortion sensitivity is documented, not fought.

### 3.2 Quadrature, rank, and the spurious-mode contract

| rule | pts | K rank (≤ 6/pt) | deficiency vs 54 | status |
|---|---|---|---|---|
| 3×3×3 | 27 | 54 (full) | 0 | `std` — Bathe's "reliable" order |
| Irons 14-pt | 14 | 54 (84 ≥ 54) | 0 **[computed]** | reserved `i14` (§7) |
| 2×2×2 | 8 | 48 | **6 spurious modes** **[computed]** | `uri` |

(Felippa AFEM §18.7: "minimum rank-sufficient rules for the 8-node and 20-node
hexahedra are p = 2 and p = 3"; Belytschko Eqs. 8.3.18–8.3.26.)

The `uri` spurious-mode contract, stated honestly:

- Abaqus ships C3D20R with **no hourglass control**: the 3D second-order
  reduced modes "can propagate in a single stack of elements … [but] rarely
  cause trouble … no special techniques are used in Abaqus to control them"
  (TG §3.2.4). The 2D analog (Q8@2×2, 1 mode) is "non-communicable" in any
  2+-element assembly (Hughes §4.6 Ex. 2).
- The textbooks document when "rarely" bites: a stiff footing on soft soil can
  arouse the Q8 mode (Hughes Fig. 4.6.4); a minimally-restrained single element
  is outright unstable at the reduced rule (Bathe Fig. 5.46); near-mechanisms
  can overshoot 500× (Cook pp. 192–193); hourglassing is "particularly
  troublesome in eigenvalue extraction" (TG §3.2.4); point loads excite the
  modes far more than distributed loads (TG §3.2.4).
- **Our contract:** P2 pins (a) exactly 6 zero-energy modes on the single
  element, (b) zero spurious modes in a restrained multi-element block, and (c)
  a *demonstrated* single-stack pathology case in the test battery as living
  documentation. Usage guidance (factory header + guide): `uri` for smooth
  production meshes ≥ 2 elements thick; `std` for eigen extraction, single-
  element-thick members, point-loaded corners, and any soft-support contrast.

### 3.3 Stress accuracy — the Barlow-point argument for `uri`

The optimal (Barlow) stress points of quadratic hexes are **the 2×2×2 Gauss
points** — ξ,η,ζ = ±1/√3, "rigorously true for rectangular elements … very
good choices" when distorted (Cook §6.13 p. 195; Barlow 1976 via Belytschko
§9.7.4). Stresses there are superconvergent (same order as displacements —
O(h³)-class vs O(h²) generic). So `uri` doesn't merely tolerate 8 points — it
evaluates the material **exactly where the element is most accurate**, which is
the theoretical basis of Abaqus' "most cost-effective element" claim
(TG §3.2.4) and of `std`'s own practice of reporting/extrapolating from GPs.
Per-GP recorder output at the 8 `uri` points is therefore the *premium* stress
product, not a degraded one. (Nodal extrapolation, where wanted, uses the
trilinear r = √3ξ map over the 8 GPs — Cook Eqs. 6.13-2…-5.)

### 3.4 Locking

- **`std` (27-pt) volumetric-locks near incompressibility** — constraint
  ratio r = 12 added dof / 27 constraints ≈ 0.44 **[computed** by Cook §9.5's
  method**]**; Abaqus lists C3D20 among the elements that "will lock … as the
  material becomes more incompressible" (AUG §28.1.1) — note the first-order
  B-bar fix upstream applies to first-order elements ONLY ("full integration
  … in first-order elements … is more accurately called selectively reduced
  integration", TG §3.2.4); C3D20 gets no such rescue, and neither does
  upstream `Twenty_Node_Brick`.
- **`uri` mostly relieves it** but is not a mixed element: C3D20R "develop[s]
  volumetric locking … after significant straining", often accompanied by an
  hourglass-*looking* mode; remedy = refine where plastic strains are large
  (AUG §28.1.1). Checkerboard pressure at neighboring GPs is the diagnostic.
- **Why no `bbar` in v1:** the quadratic B-bar exists (Hughes §4.5.2 Ex. 3 —
  project the dilatational rows onto a linear pressure field; the Q9/P1-class
  element) and generalizes to nonlinear form, but the constraint count says
  selective 2×2×2 on H20 lands at r = 1.5 — the 3D analog of Q8's known
  sub-optimal 6/4, historically "viewed by many as more a 'trick' … some bad
  experiences … for the serendipity element" (Hughes p. 226). The ratio-optimal
  quadratic brick is **H27** + selective/B-bar (r = 3), i.e. incompressibility
  demand should pull the H27 sibling (or the fork's H8 `bbar`/`eas`), not a
  half-measure on H20. Keyword reserved; decision recorded.
- **Shear/bending: no cure needed.** Quadratic interpolation natively carries
  the linear strain field — one element through the thickness suffices in
  bending (TG §3.2.1); parasitic shear is a first-order disease. This is the
  positive half of the element's reason to exist: coarse-mesh bending accuracy
  that H8 needs `eas` + mesh stacks to approach, without internal-mode state.

### 3.5 Mass & lumping (the explicit-correctness core)

- **Consistent** (default): 27-pt ∫ρNᵀN, formulation-independent. Abaqus uses
  consistent mass for ALL second-order solids in Standard (lumping is
  first-order-only there, TG §2.4.1); consistent-mass eigenfrequencies converge
  from above at 4th-order-class rates (Hughes §7.3.2).
- **Row-sum is forbidden** — this is the load-bearing quadratic-mass fact:
  m_a = ∫ρN_a dΩ goes **negative at H20 corners (−M/8 each; mid-edges +M/6)**
  [computed; 2D analog −M/12 / +M/3 book-printed: Hughes Ex. 8 §7.3.2, Cook
  Fig. 13.3-4]. Negative lumped masses give exponentially growing SDOF modes —
  "vitiating the true solution" (Hughes Remark 1 §7.3.2) — and LS-DYNA's 1979
  "numerical noise resulting from the ad hoc mass lumping" verdict is the same
  story lived in production (Theory §1.1).
- **`-lumped` = HRZ** (diagonal-of-consistent scaled to conserve directional
  mass): "always produces positive lumped masses … the only lumping method
  that can be recommended for arbitrary elements" (Hughes §7.3.2 Eqs.
  7.3.47–48). The fork already ships the exact algorithm as the shared
  header-only `Ladruno::hrzLump` ([[35_ladruno_hrz_lumped_mass_adr|ADR 35]],
  `LadrunoMassLumping.h`) — the element's `-lumped` path calls it on the 27-pt
  consistent mass (positivity guard passes by construction: consistent
  diagonals are Gram entries ∫ρN_a² > 0). Cube fractions to pin in the oracle:
  **corners 7/248 ≈ 2.82%, mid-edges 2/31 ≈ 6.45%** [computed; method
  reproduces Cook's printed Q8 values 3/76 & 16/76 exactly]. Known cost: HRZ on
  serendipity elements sacrifices the optimal-lumping convergence rate (no
  positive optimal lumping *exists* for serendipity — Cook pp. 373–374); H27's
  Lobatto lumping is the book-sanctioned positive+optimal combination (§2.6,
  §7). The domain-level integrator `-hrz` flag (ADR 35/52) composes identically
  and needs zero element work.

### 3.6 Explicit dynamics — supported, honestly discouraged

- **Δt penalty:** the quadratic rod's lumped ω_max = 2√6·c/h ⇒ Δt_crit factor
  **1/√6 per element length, 2/√6 ≈ 0.82 per equal nodal spacing** vs linear
  elements (Hughes §9.2 Eqs. 9.2.5–6, p. 514; Belytschko Table 6.1) — "still
  the advantage is with linear elements." No closed-form ω_max exists for
  quadratic continua (Hughes p. 517); LS-DYNA falls back to geometric
  characteristic lengths (tet-10: "characteristic length … half that of a
  4-noded tetrahedron ⇒ time step smaller by a factor of 2", Vol I Remark 11).
  **The fork needs no such heuristic:** `CriticalTimeStep` (ADR 65) already
  solves the per-element algebraic pencil (60-DOF K,M eigen), which is exact
  for any element — an architectural advantage worth stating. SMS composes the
  same way (pencil-based).
- **Wavefronts:** "lower-order displacement elements are more adept at modeling
  [propagating strain discontinuities] … higher-order elements … produce more
  numerical noise" (Cook §13.14 p. 418); higher-order pays off in hyperbolic
  problems only with smooth data + consistent mass (Belytschko p. 486) — which
  explicit can't use.
- **Position:** explicit runs are *permitted* (HRZ mass is positive; CDL/SMS
  pencils are exact) and P3 proves one, but the guide steers explicit users to
  `LadrunoBrick uri` / Bézier. This mirrors the industry line (Abaqus: no
  second-order hexes in Explicit at all; LS-DYNA: supported but "under
  development", default answer stays linear + hourglass control). No bulk
  viscosity in v1.

### 3.7 Nonlinear materials — where quadratic pays and where it doesn't

- **Pays:** smooth elastic/hardening problems (J2 with hardening, confined
  concrete pre-softening, foundations, contactless machine-part stress) — the
  h³/h² rates and Barlow sampling hold.
- **Doesn't:** C⁰ constitutive response caps rates at h²/h¹ (Belytschko
  §8.2.4); localization/softening wants the most element-edge directions per
  DOF ("first-order elements are usually recommended" for large-strain
  plasticity/limit loads, TG §3.2.4). `getCharacteristicLength` returns
  lch = ∛V (the `LadrunoBrick` convention) so crack-band materials *regularize
  consistently*, but the guide states plainly: softening concrete / shear
  bands → H8 family (`LadrunoBrick` + ASDConcrete3D/LadrunoConcrete3D), not
  H20. (Crack-band theory is calibrated on constant-strain fields; a band
  localizing across a sub-cell of a quadratic element makes ∛V the wrong
  band width.)

### 3.8 Consistent loads & the contact pathology (why contact faces are out)

Consistent (variational) nodal loads on quadratic faces are sign-mixed:
uniform pressure p·A on a serendipity face pulls the corners **backwards**
(corner −pA/12, mid-side +pA/3 — textbook values; the sign flip itself is
manual-documented, AUG §39.1.2 Fig. 39.1.2-6). Consequences the ADR bakes in:

- **SelfWeight / `-b` / face tractions must go through the consistent N-path**
  (they do — same code shape as `LadrunoBrick`); users hand-lumping gravity
  onto nodes of quadratic meshes get the wrong answer by construction.
  Document loudly in the guide.
- **NTS contact against quadratic faces is a known failure mode** — the
  sign-mixed nodal forces make per-node contact decisions "ambiguous"; Abaqus'
  own fix is converting C3D20→C3D27 slave faces (midface node ⇒ all-positive
  distribution) or going surface-to-surface/mortar (AUG §39.1.2); LS-DYNA's is
  part-ID-defined contact surfaces + mortar segments with true quadratic
  geometry (Vol I Remark 11; Remark 14). The fork's NTS engine
  ([[39_ladruno_contact_domain_adr|ADR 39]]ff) builds linear facets — **v1
  simply does not offer LadrunoBrick20 faces to the contact/tie engines**
  (factory/docs guard). The clean future options, in order: tie the quadratic
  region to a linear contact skin ([[62_ladruno_kinematic_mesh_tie_adr|ADR
  62]]), mortar with quadratic segments (ADR 41 extension), or the H27 route.

### 3.9 Cross-code positioning

| Concern | Abaqus | LS-DYNA | upstream OpenSees | this ADR |
|---|---|---|---|---|
| Quadratic hex | C3D20 (27-pt) / C3D20R (8-pt, the recommended default) / C3D27 variable-node | ELFORM 23 (H20, R8+) / ELFORM 24 (H27 S/R, "under development") | `Twenty_Node_Brick` 27-pt only | `LadrunoBrick20 -formulation {std|uri}` |
| Hourglass control on reduced | none on C3D20R (modes benign) | 14-pt rule chosen to *eliminate* modes (Theory §3.7) | n/a | none; modes gated + documented (§3.2) |
| Volumetric locking | C3D20 locks (listed); C3D20R late-strain lock; hybrids H for incompressible | SRI on ELFORM 24 | locks (27-pt, no B-bar) | `uri` relief; `bbar` reserved → H27/H8 for hard incompressibility |
| Mass | consistent (2nd-order never lumped in Standard) | row-sum lumped everywhere (noise verdict, 1979) | consistent only | consistent 27-pt; `-lumped` = HRZ (positive by construction) |
| Explicit | **excluded** (C3D10M only) | supported, niche | n/a | permitted, discouraged; pencils exact via `CriticalTimeStep` |
| Contact on quadratic faces | auto-convert →C3D27 or S2S | part-ID surfaces / mortar quadratic segments | n/a | **out of scope v1** (linear-facet engine) |
| Stress output | GPs; extrapolation caveats | 8-value d3plot cap (NINTSLD), SOLSIG for 23 | per-GP `setResponse` | per-GP tree + `stress3D6`; recorder rules already tabulated |

The deliberate divergences from both codes: (i) HRZ element-level lumping
(neither ships it for these elements — Abaqus refuses to lump, LS-DYNA row-sums
and eats the noise); (ii) algebraic per-element Δt pencils instead of
characteristic-length heuristics; (iii) an explicit reduce-to anchor against a
living upstream element.

---

## 4. Architecture

```
SRC/element/ladrunoBrick/
  LadrunoHex20Shape.h            # NEW, PURE header (no OpenSees deps): serendipity
                                 #   N/∂N∂ξ, Jacobian, B; 27-pt (brcshl order) +
                                 #   8-pt (Brick order) tables; g++/numpy-testable
                                 #   (LadrunoMassLumping.h / LogStrainKernel.h idiom)
  LadrunoBrick20.{h,cpp}         # NEW element (33018)
  OPS_LadrunoBrick20.cpp         # NEW Tcl/Python factory
  LadrunoBrick.{h,cpp}           # untouched
```

**Public API**

```
element('LadrunoBrick20', tag, n1..n20, matTag,
        ['-formulation', <std|uri>]        # default std
        ['-lumped']                        # HRZ lumped mass (default consistent)
        ['-b', bx, by, bz]                 # body force
        ['-damp', dampTag])                # Damping objects, per GP
```

**Element responsibilities** (thin, mirrors `LadrunoBrick`): node/GP plumbing
sized 20/27; state cycle over `materialPointers[27]` (std) or `[8]` (uri);
`formResidAndTangent` = plain B·u Gauss loop (no seams needed until `-geom
finite` — the kinematics ledger stays trivially linear in v1); consistent-mass
+ HRZ paths; `zeroLoad`/`addLoad(SelfWeight)`/`appliedB`; `setResponse` tree =
`material`/`stresses`/`strains` per-GP + `stress3D6`/`strain3D6` averages;
`getCharacteristicLength` = ∛V; `Print` JSON; `sendSelf`/`recvSelf` (nodes,
formulation ordinal, massType, b[], per-GP materials); `displaySelf`. The
upstream `setParameter` `i<4` damping-loop bug class is not inherited (loop
bounds from one `NGP` constant).

**Kernel responsibilities** (`LadrunoHex20Shape.h`, written once, pure):
shape values/derivatives at arbitrary ξ, both GP tables + weights, Jacobian/
inverse/det, the 6×60 B fill, and the natural-coordinate node table. Consumed
by the element, the P0 oracle (standalone g++), and later by the embedded-host
APIs / basisInfo if demanded.

**Recorder path:** class-tag route, zero new math — add `33018` beside
`ELE_TAG_Twenty_Node_Brick` in the fork-owned `Ladruno_ElementResults.h` tag
lists (rule `Hexahedron_GaussLegendre_3` for std; `Hexahedron_GaussLegendre_2`
for uri — both already tabulated) and to the four viz recorders'
`vtktypes[...] = *_QUADRATIC_HEXAHEDRON` maps (VTK/VTKHDF/PVD/Gmsh, additive
one-liners). The frozen `MPCORecorder.cpp` is NOT touched
([[feedback_vanilla_footprint]] — LadrunoRecorder is the blessed path).

**apeGmsh:** gmsh emits hex20 via `Mesh.SecondOrderIncomplete=1` (its default
second-order hex is the 27-node); node-order permutation gmsh→OpenSees lives on
the apeGmsh side; contract doc row updates from "TBD" to the command + ordering
note when P1 lands.

---

## 5. Class tags & registration

- `ELE_TAG_LadrunoBrick20 = 33018` (ELE registry; classTags.h comment says
  "33017-33019 remain free", and the same-day ADR-71 `LadrunoUP` draft claims
  33017 → this ADR takes 33018, H27 pencils 33019). Reserve **at P0** with the
  standard annotated `#define` (the ADR-70 precedent: reservation landed with
  #540, not the ADR). Whichever of ADR-71/ADR-72 reaches P0 second must
  re-check the tag frontier.
- `LadrunoHex20Shape.h`: no classTag (pure helper).
- Vanilla touches (strictly additive, `// Ladruno`, ledgered):
  `SRC/classTags.h`; `FEM_ObjectBrokerAllClasses.cpp`; Tcl + Python element
  command registries; `SRC/element/ladrunoBrick/CMakeLists.txt`;
  `VTK_Recorder.cpp` / `VTKHDF_Recorder.cpp` / `PVDRecorder.cpp` /
  `GmshRecorder.cpp` vtktype maps. Fork-owned: `Ladruno_ElementResults.h`.
- `Ladruno_scripts/banner_features.txt` + `patch_banner.py` at ship;
  [[LEDGER_implementations]] row at P1; stamp headers on all new files
  ([[feedback_always_stamp_header]]).

---

## 6. Phases & exit gates

| Phase | Scope | Gate |
|---|---|---|
| **P0** | `LadrunoHex20Shape.h` + numpy oracle (standalone g++, no OpenSees build); tag reservation | N: partition of unity, Kronecker-delta at nodes, Σ∂N = 0, degree-2 completeness vs symbolic; J/detJ on a distorted hex vs numpy; K assembled from the kernel: **rank 54 @ 27-pt, rank 48 @ 8-pt (6 spurious modes, eigenvectors dumped & catalogued)**; 14-pt note verified or dropped; consistent-mass diagonals > 0; **HRZ cube fractions == 7/248 & 2/31; row-sum corner negativity −M/8 demonstrated** (the pin that justifies HRZ-only); GP tables match brcshl (27) / Brick (8) orders |
| **P1** | element `-formulation std` (27-pt) + registration + factory + recorder wiring | **reduce-to `Twenty_Node_Brick`**: K, resisting force, consistent M to ~1e-12 on a distorted hex under mixed loads; constant-strain patch (distorted 2-element mesh); rank/6-RBM; SelfWeight/`-b` vs closed form (corner shares near-zero/negative — asserted, documented); lch = ∛V; per-GP recorder round-trip via LadrunoRecorder (rule reuse, GP coords vs `GLOBAL_GP_COORDS` oracle); sendSelf/recvSelf round-trip; ledger + banner |
| **P2** | `-formulation uri` (8-pt) | single element: rank 48, exactly 6 zero-energy modes; restrained ≥2×2×2 block: zero spurious global modes (non-communicability pin); **single-stack pathology case demonstrated & documented** (the honest counter-example); coarse bending: 1-layer cantilever tip disp ≥ 0.98·analytic where `LadrunoBrick std` locks and `eas` needs refinement (Oñate Fig. 8.23 replication); ν = 0.4999 confined block: `uri` relieves, `std` locks (both documented); Barlow check: `uri` GP stresses beat `std` GP stresses vs the analytic bending field; cost: assembly time ratio ≈ 0.3–0.4× std. **R3-deferred debts due here**: (a) recorder metadata seam — the legacy `getGeometryAndIntRuleByClassTag` arm hardwires 33018→`Hexahedron_GaussLegendre_3` (27 GPs); `uri` (8 material points, same classTag) MUST answer the `basisInfo`/CustomIntegrationRule probe instead, and the P1 battery's unconditional `NUM_GP==27` pin must become formulation-aware; (b) blocked/sparse BᵀDB assembly (per-node 3×3 blocks from `dNdx`, symmetry exploited — upstream uses `addMatrixTripleProduct`) — do it when the assembly is touched for `uri` anyway, re-anchored by the reduce-to gate; (c) battery table dedup — `test_ladrunoBrick20_element.py` re-transcribes `NODE_XI`/GP27/edges from the oracle; import `tests/hex20_reference.py` (sympy now in Zone-A) so ordering corrections can't drift. |
| **P3** | dynamics: `-lumped` (HRZ), explicit proof-of-life | eigen: consistent brackets from above, HRZ below, both → analytic bar/beam frequencies; HRZ mass vector matches P0 fractions through the element path; wave bar under `CentralDifferenceLadruno` runs stably at the pencil Δt (measured Δt ratio vs equal-node H8 mesh reported, ≈ 0.8 rod-theory ballpark); `criticalTimeStep()` returns the 60-DOF pencil value; energy-balance channels close (no hourglass channel — asserted zero); **betaK-Rayleigh clobber gate** (the ADR-66/#562 family bug: `getResistingForceIncInertia` under `rayleigh 0 βk 0 0` must keep M·a and −Q — Brick20's R3 F12/F13 static-buffer split should make it structurally immune; assert it like the plane family's ×4 gate) |
| **P4** *(deferred, demand-gated — each its own mini-ADR/PR)* | `-geom finite` on a shared 3D kernel extracted from `LadrunoBrick`; `bbar`/H27 sibling (33019) if incompressibility or contact demand materializes; `i14` rule; embedded-host APIs (`getInterpolationWeights/Gradients` sized 20 + `LadrunoEmbeddedKernel` dispatch); quadratic tie/mortar faces | (scoped when opened) |

Adversarial-gate policy ([[feedback_adversarial_gate_when]]): **no full Opus
gate planned** — P0/P1 are anchored by textbook oracles + the upstream
reduce-to (the strongest coverage class we have), and the `uri` risk surface is
pinned by explicit P2 gates; escalate to a full gate only if a P2 gate
surprises (e.g., the mode census or the stack pathology disagrees with the
§3.2 predictions). Vanilla touches are additive one-liners.

Each phase is one PR off `ladruno`.

---

## 7. Risks / open questions

> [!question]
> **`i14` (Irons 14-pt) — worth a third formulation?** Rank-sufficient
> [computed], ~half the 27-pt cost, and it is what LS-DYNA's H20-based element
> actually uses ("necessary to eliminate the zero energy modes … nearly the
> same accuracy as the twenty-seven point Gauss rule", Theory §3.7; Cook p. 189
> recommends it "especially if the element is made very thin in one
> direction"). Costs: a nonstandard GP layout (new recorder rule tabulation +
> basisInfo/GLOBAL_GP_COORDS entries) and stress points that are neither Barlow
> nor familiar. Lean: **reserve the keyword, revisit only if `std` assembly
> cost actually shows up in a profile** (ADR-40 data says element `update()`
> dominates fiber/solid lanes — 27 vs 14 material evaluations is the real
> difference, and for smooth problems `uri`'s 8 already win).

> [!question]
> **H27 sibling trigger.** Three futures pull toward it (explicit Lobatto
> lumping, distorted-mesh robustness, r = 3 incompressibility; §2.6) — but
> which demand arrives first decides the *shape* of the element (explicit →
> Lobatto-lumped `std`-only H27, possibly Bézier/Bernstein per the Kadapa
> line; incompressibility → H27 + `bbar`; contact → H27 faces for the
> all-positive load distribution). Don't build speculatively; record the
> decision tree here. NB the one variant to never build: reduced H27
> (C3D27R's "three unconstrained, propagating hourglass modes", AUG §28.1.1;
> ≥ 27-mode deficiency at 2×2×2 [computed]).

> [!question]
> **Softening materials on H20 — warn or allow silently?** lch = ∛V
> regularizes consistently, but crack-band theory on quadratic fields is
> shaky (§3.7). Options: (a) docs-only guidance, or (b) a one-time `opserr`
> advisory when a material with a "damage" response is attached.
> **RESOLVED 2026-07-10 (user, U1): BOTH** — the guide documents the caveat
> AND the element emits a one-time `opserr` advisory at `setDomain` when the
> attached material exposes a "damage" response channel (the cached-Response
> probe pattern from `LadrunoBrick::damageResponse`). Advisory, not a
> rejection — the run proceeds. Lands with P1.
> *P1 implementation note*: the probe requires a **dual-scalar** damage
> channel (response Size ≥ 2, the {d⁺,d⁻} / {ω_t,ω_c} crack-band concrete
> family: ASDConcrete3D, LadrunoConcrete3D). Rationale: `LadrunoJ2` answers
> the `"damage"` probe with a size-1 Lemaitre channel even when its damage
> law is OFF, so a mere-presence probe would false-positive on every J2
> metal. Consequence accepted: a Lemaitre-configured LadrunoJ2 (scalar
> ductile damage, not crack-band) does not trigger the advisory — its
> regularization story is the (F)-family's, documented there, not here.
> Known limitations (R3): the probe runs on the virgin material at `setDomain`
> — a material that lazily allocates its damage vector on first strain would
> silently skip the advisory (none of the current family does); and the
> probe+`Size≥2` discriminator is inline in `LadrunoBrick20::setDomain` — when
> the H27 sibling (33019) or another consumer lands, hoist it to a shared
> `isCrackBandMaterial(NDMaterial*)` predicate beside the H8
> `damageResponse` pattern instead of a third copy.

- **Node-ordering verification is a P0 gate, not an assumption.** The corner/
  mid-edge pattern (bottom ring 9–12, top 13–16, vertical 17–20) is common to
  upstream/Abaqus/LS-DYNA, but the exact edge sequence within each ring must be
  read out of `brcshl`/`shp3dv.cpp` once and pinned; a silent mismatch would
  poison the reduce-to gate *and* the recorder GP pairing. gmsh's differing
  hex20 order is apeGmsh's job (permutation + `SecondOrderIncomplete`).
- **`uri` + eigen extraction:** spurious low-frequency modes are the documented
  failure (TG §3.2.4). The guide says eigen → `std`; P2's mode census makes the
  failure legible. Domain `-hrz` flag composes with either.
- **Static workspace:** stiff/resid/mass statics sized 60 follow the house
  norm; the ADR-40 de-static audit sweeps this element with the rest of the
  library when it happens (no special debt added).
- **Upstream coexistence:** `Twenty_Node_Brick` stays untouched (reduce-to
  anchor + back-compat); we add, never modify.
- **Backwards compatibility:** new class + new tag ⇒ zero impact on existing
  models.

---

## 8. References

- Felippa, *Advanced FEM* (2000): Ch. 18 hexahedra (§18.3 H20 shape functions,
  §18.7 rank-sufficient rules), Ch. 15 Remark 15.1 (27-pt H20 ≈ 22× H8 K-cost).
- Cook, Malkus & Plesha 3rd ed.: §6.11 quadrature/Irons 14-pt, §6.12 spurious
  modes (Q8/Q9 census, near-mechanisms), §6.13 Barlow points & extrapolation,
  §9.5 constraint counting, §13.3 mass lumping (HRZ recipe, Fig. 13.3-3/-4),
  §13.10/13.14 explicit Δt & wave guidance.
- Hughes, *The FEM* (1987): §4.3.7/4.4/4.5.2 constraint ratios, Malkus–Hughes
  equivalence, quadratic B-bar (Ex. 3); §4.6 rank/stability (footing example);
  §7.3.2 lumping (row-sum negativity Ex. 8, HRZ 7.3.47–48, Lobatto Ex. 7);
  §9.2 explicit Δt (quadratic rod 9.2.5–6, no closed forms p. 517).
- Belytschko, Liu, Moran & Elkhodary (2014): §8.2.3–8.2.4 convergence rates &
  the C⁰-material cap, §8.3 rank/patch/completeness (Ciarlet–Raviart), Table
  6.1 element frequencies, §9.7.4 Barlow.
- Oñate, *Structural Analysis with the FEM* Vol. 1 (2009): §8.5 H20/H27 shape
  functions, §8.10–8.11 quadrature + the H20-vs-H27 distortion verdict
  (p. 296), Fig. 8.23 single-layer bending.
- Bathe (1982): Table 5.5 "reliable" orders, Fig. 5.46 reduced-rule
  instability, §5.8.2 GP stress accuracy.
- Abaqus 2016: TG §3.2.1/3.2.4 (integration, Barlow, spurious modes, element
  selection), §2.4.1 (mass); AUG §28.1.1/§28.1.4 (C3D20/R/H, C3D27 family,
  locking lists, Explicit exclusion, auto-conversion), §39.1.2 (quadratic-face
  contact pathology, Fig. 39.1.2-6).
- LS-DYNA: Theory §1.1 (the 1979 higher-order verdict), §3.7 (14-pt rule),
  §22.1 (Δt), §34 (row-sum); Vol I R16 `*SECTION_SOLID` ELFORM 23–29 + Remarks
  11/16, `*ELEMENT_SOLID_H20` Fig. 19-31, `*CONTROL_IMPLICIT_CONSISTENT_MASS`,
  mortar Remarks (11-131).
- Fork: [[09_ladruno_brick]] (the element pattern + seams),
  [[70_ladruno_plane_finite_triangles_adr]] (quadratic sibling + shared-kernel
  lesson), [[35_ladruno_hrz_lumped_mass_adr]] (`LadrunoMassLumping.h`),
  [[65_ladruno_explicit_dt_strategies_adr]] (`CriticalTimeStep` pencils),
  [[04_bezier_elements]]/[[06_bezier_tet10]] (the explicit-quadratic lane),
  [[ladruno_element_contract]] + `Ladruno_ElementResults.h` (recorder rules),
  [[ladruno_apegmsh_contract]] (the reserved hex row), upstream
  `SRC/element/brick/Twenty_Node_Brick.{h,cpp}` (reduce-to anchor).

## 9. Implementation log

*(filled in as phases land; move to `Ladruno_internal/` when complete)*

- 2026-07-10 — ADR drafted. Research pass: Abaqus TG/AUG + LS-DYNA Theory/Vol I
  + Felippa/Cook/Hughes/Belytschko/Oñate/Bathe (three parallel research
  agents; all citations re-verified against sources; [computed] flags mark
  derived-not-printed numbers). Key calls: separate `LadrunoBrick20` class
  (33018), `std|uri` only in v1, HRZ-only lumping (row-sum −M/8 corners),
  explicit permitted-but-discouraged, contact faces excluded, H27 kept as a
  demand-gated sibling (33019 candidate).
- 2026-07-10 — Renumbered 71→72 same day: a parallel session drafted ADR-71
  (`71_ladruno_up_family_adr.md`, LadrunoUP Biot u-p, uncommitted on its own
  worktree) first and claimed ELE 33017 → this element moved to 33018, the H27
  pencil to 33019.
