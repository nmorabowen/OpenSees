---
title: "ADR 66 — LadrunoSolidShell (ELE_TAG 33020): 8-node ANS/EAS solid-shell — the through-thickness σ33 host for punching/bearing"
project: Ladruno
type: ADR / scoping
status: "G9 EXPLICIT PATH LANDED (2026-07-06, SRC + tests): probing exposed the real D7 gap — the P5.1 element was CONSISTENT-mass only, and under the explicit recipe (system Diagonal) the SOE keeps just the raw diagonal = 8/27 of the element mass, so explicit dynamics ran 3.375x light and the '-lump hrz' pencil over-reported the stable step 1.84x (measured blow-up at 0.9x, true boundary 0.54x = the consistent-eig value; the CDL '-lump' flag models the PENCIL only — new LEDGER_quirk). FIX: element '-lumped' flag (LadrunoBrick donor: row-sum diagonal rho*int N_a dV == HRZ on the trilinear shape; matching lumped inertia force; massType in the spare idData(6) wire slot, old streams byte-compatible). tests/test_ladrunoSolidShell_explicit.py: mass-correctness (free axial period == 4L/c to 5%; the deficient path rings 1.84x fast), pencil == t/c + thickness-halving + bisection (stable 0.9x, unstable <= 1.5x), massType roundtrip (physics-discriminating: dropped flag = blow-up in ~25 steps), SMS sliver-mesh speedup (dtTarget = bulk step, injection only at the 0.01-wide sliver column; the UNIFORM-mesh 20x-target arm was measured under-delivering 42x — the documented scales-the-bulk trap, reframed), Concrete3D softening band to full localization under CentralDifferenceSMS at 10x the pencil (Tier-3 + D7 cure end-to-end, 0.7*ft seed vs the dynamic localization overshoot).\n\nG3 BENCHMARKS LANDED (2026-07-06, test-only): tests/test_ladrunoSolidShell_benchmarks.py — Scordelis-Lo quarter roof ans 0.948/0.974/0.989/0.992 (4/8/16/24; ~1% at 16x16, inside the ADR band; std control 0.123 = 8x lock) + pinched cylinder w/ diaphragms eighth model ans 0.378/0.745/0.927/0.970 (16x16 = 0.927 == the canonical MITC4 ~0.93 — the expected class result for a shear/thickness-tied element WITHOUT in-plane membrane enhancement; slow-from-below tail = the O3 in-plane-EAS backlog signature, documented not defect; std control 0.069 = 11x lock). First curved-mesh coverage (radial zeta-fibers, ANS off its exact flat-face class). Monotone-convergence + std-discrimination gates, bands pinned to measurement with margin.\n\nP5.2 GATES LANDED (2026-07-06, TEST-ONLY — no source change needed): tests/test_ladrunoSolidShell_softening.py — Concrete3D hosting parity (solid-shell GP stress == NDTest replay of the same strain path, all 6 comps; alpha stays exactly 0 on uniform fields), LadrunoRCConcrete hosting smoke, G5 softening-band energy mesh-objectivity (weakened-column Bazant band, in-plane-SQUARE elements dx=2/1/0.5, -autoRegularization + -implex: W = 0.01752/0.01739/0.01737 -> 0.9% spread across the 4x range, vs the fixed-lch control at 0.25x = the band-volume ratio; band omega_t = 1.0, zero neighbor leakage), G6 EAS alpha-stability BOTH arms (fully-prescribed bending ramped ~4.5x eps_f deep post-peak on the steep -lch 50 law: alpha finite/bounded/settling; free-DOF uniaxial-stress arm through the indefinite-Kaa freeze path: alpha pinned at 0, no corruption) — the P5.1 freeze-alpha guard SUFFICES, the fallback ladder (freeze-post-peak/mode-removal/Kstab) stays dormant. Three new LEDGER_quirks: sqrt(A)-lch needs in-plane-square band meshes (O4), hairline band seeds (0.95 ft) drive bulk GPs into the return-map apex (use 0.8 ft), multi-element implicit softening cut-crawls without -implex (13 min -> 16 s).\n\nP5.1 BUILT + ADVERSARIALLY REVIEWED (2026-07-02/03): element + registration + gates 16/16 green (Concrete3D 28/28 non-regressed). A 5-lens / 89-agent adversarial review found ONE critical defect — getResistingForceIncInertia clobbered the shared static resid during the Rayleigh re-entry, silently dropping element inertia under stiffness-proportional Rayleigh (no dynamic test existed to catch it). FIXED via the LadrunoBrick snapshot pattern + a dynamic-Rayleigh regression gate PROVEN discriminating (fails `analyze -3` on the reverted-buggy binary, passes on the fixed). Also fixed: null-material parse-time probe, recvSelf stream validation + Ki invalidation, Print/broker null guards, one-time center-Jacobian geometry check, constructor nz guard. Test coverage strengthened: rigorous no-re-solve alpha-survival roundtrip, gamma_23 (y-axis) transverse-shear gate, getInitialStiff==getTangentStiff elastic gate, setResponse forwarding. LEDGER_vanilla_files backfilled. Left with skeptic-refutation backing: setDomain-continues-on-bad-node (OpenSees norm, donor same), per-GP-detJ-in-hot-loop (center check added; concave=user error), unknown-option warn-continue (fork idiom).\n\nP5.1 BUILT (2026-07-02): element + registration + gates green on a fresh worktree build (Concrete3D battery 28/28 non-regressed). Corrects one stale ADR-19 risk-register claim (LadrunoBrick now HAS true state-dependent Simo-Rifai EAS — the machinery is an adapt, not net-new). P5.1 deviations: SolidTransformation seam not plumbed (v1 assembles directly in the global frame; parser rejects -geom corot/finite — seam wiring moves to P5.4); state-dependent EAS (committed alpha) shipped in P5.1 rather than P5.2 (the formEAStrue adapt made it cheaper than the initial-tangent stub); G3 pinched-cylinder/Scordelis-Lo deferred to a follow-up test PR. Observed: ~15% stiffening residual on the 27-degree trapezoidal-fiber mesh (in-plane trapezoid coupling; Betsch-Stein removes only the curvature-thickness part) — gate band set accordingly."
related:
  - "[[19_ladruno_rc_shell_adr]]"        # the parent — Phase 5 mandate, D1 host decision, 33020 reservation
  - "[[31_ladruno_concrete3d_adr]]"      # LadrunoConcrete3D (33017) — the headline triaxial material this element unlocks in shell geometry
  - "[[09_ladruno_brick]]"               # LadrunoBrick (33002) — formEAStrue Simo-Rifai machinery + damage-scaled Kstab precedent
  - "[[rc_shell_phase3plus_handout]]"    # the forward handout whose Phase-5 section this ADR supersedes
  - "[[35_ladruno_hrz_lumped_mass_adr]]" # HRZ lumping — explicit mass path
  - "[[36_ladruno_selective_mass_scaling_adr]]" # SMS — the cure for the thickness-bound explicit dt_cr
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
tags: [adr, element, solid-shell, ans, mitc, eas, locking, punching, rc-shell, through-thickness, proposed]
updated: 2026-07-01
---

# ADR 66 — `LadrunoSolidShell` (ELE_TAG 33020)

**One sentence.** An 8-node, 3-translational-DOF-per-node solid-shell element — ANS transverse-shear
+ ANS thickness-strain + EAS-on-`E33` — that carries the full triaxial stress state (`σ33 ≠ 0`)
through shell-like geometry, as the **narrow punching/bearing/3D-crush specialist** the director
shell (`ASDShellQ4` + `LayeredShellFiberSection`) structurally cannot be.

**Status.** PROPOSED, scoping only, no code. This is the ADR that ADR 19 Phase 5 deferred; it
supersedes the Phase-5 section of `rc_shell_phase3plus_handout.md`. The ELE_TAG **33020**
reservation from ADR 19 is honored (ELE slots 33016–33019 remain free for other elements; tag bands
are per-registry).

---

## 1. Driver

ADR 19 built the RC shell stack as *material-on-frontier-host* and identified exactly **one
elemental blind spot** it could not discharge: the director-shell seam condenses `σ33 = 0` per layer
(the `PlateFiber` view's guarded σ33-Newton), so it cannot represent

- **punching shear** at slab–column connections (a triaxial compression-shear state around the
  loaded area — the ADR 19 documented xfail this element discharges);
- **bearing / anchorage-zone crush** (concentrated loads, post-tensioning anchorages);
- **through-thickness confinement** — and this is where the fork's materials are *ahead of the
  element inventory*: `LadrunoConcrete3D` (33017) carries a Menétrey-Willam 3-invariant surface with
  **confinement-dependent ductility by mechanism** ("Mander by mechanism", proven to ≤2.9% against
  the Mander formula in the P5 confined-fiber gates), but no shell-geometry host can feed it a
  triaxial state today. The solid shell is the missing consumer.

A plain brick mesh (`LadrunoBrick`) can of course do all of this — at the cost of either
catastrophic locking at shell slenderness (full integration, Belytschko Ch.8) or hourglass-managed
reduced integration with 3+ elements through the thickness. The solid shell exists to make
**one-to-few elements through the thickness at shell aspect ratios** a converged, locking-free
proposition.

**What this element is NOT** (locked in ADR 19 D1, restated so it does not creep): not a co-equal
flexural host — walls and slabs in bending stay on `ASDShellQ4` + `LayeredShellFiberSection`; not a
new director shell; not a 6-DOF drilling element (that would defeat brick/contact parity and the
`SolidTransformation` seam). Any flexural claim must first pass the moment-curvature benchmark
against `LayeredShellFiberSection` (§7, gate G7).

## 2. Prior art (why this exact recipe)

The 8-node/3-DOF ANS+EAS solid-shell is the settled industrial and academic answer:

- **Abaqus SC8R "continuum shell"** — 8-node, displacement-only, assumed-strain treatment +
  reduced integration; the direct commercial precedent for "shell response from a brick topology".
- **LS-DYNA thick-shell (`*ELEMENT_TSHELL`) family** — assumed-strain solid-shells with full 3D
  constitutive input; same "thin element, thick material" division ADR 19 adopted for the director
  host. (Per-ELFORM details cited loosely by design — verify against the on-disk Theory Manual
  before using as durable citations, per the ADR 19 D6 honesty rule.)
- **Academic spine** (the formulation this ADR pins to):
  - Dvorkin & Bathe (1984) — ANS transverse shear (the MITC4 sampling, reused here on the
    solid-shell midsurface);
  - Betsch & Stein (1995) — ANS thickness strain (bilinear interpolation of `E33` from the four
    corner sampling points; cures curvature-thickness/trapezoidal locking);
  - Simo & Rifai (1990) — EAS; here a minimal `E33`-directed mode set curing Poisson-thickness
    locking;
  - Hauptmann & Schweizerhof (1998), Klinkel–Gruttmann–Wagner (1999, 2006), Vu-Quoc & Tan (2003) —
    the assembled ANS+EAS solid-shell, including materially nonlinear use.

The locking inventory this recipe must clear (and the cure assignment):

| Locking mode | Trigger | Cure in this element |
|---|---|---|
| Transverse shear | thin bending, `t/L → 0` | ANS (Dvorkin-Bathe edge-midpoint sampling of `E13`,`E23`) |
| Curvature-thickness (trapezoidal) | curved/warped mesh, sloped lateral edges | ANS (Betsch-Stein corner sampling of `E33`) |
| Poisson thickness | bending with `ν ≠ 0` (constant `E33` cannot follow linear-in-ζ in-plane strain) | EAS on `E33` (linear-in-ζ enhanced mode) |
| Volumetric | incompressible limit | out of v1 scope (concrete ν≈0.2; document; B-bar interaction with ANS is a research topic) |
| In-plane shear/membrane | element distortion | v1 = full 2×2 in-plane integration, accept mild stiffness; optional in-plane EAS modes deferred (§9 O3) |

## 3. The stale-claim correction (load-bearing for the budget)

ADR 19's D1 re-scope note and the Phase-5 handout both carry:

> "`LadrunoBrick`'s 'EAS' is NOT a reusable template. It is the SSPbrick stabilized-single-point
> scheme … constant `Kstab` … Genuine EAS-on-`E33` … is **net-new Simo-Rifai code**, not a port."

**This is now half-stale.** It was true of the *original* `-formulation eas` (renamed `ssp`,
2026-06-02). Since 2026-06-03, `LadrunoBrick -formulation eas` is **true Simo-Rifai EAS**
(`LadrunoBrick.cpp: formEAStrue`, ~line 2614): full 2×2×2 integration, per-GP `ε = B·u + M·α`, an
**inner Newton solving the 9 enhanced parameters against the live material stress**
(`∫MᵀσdV = 0`, relative tolerance), static condensation `K* = K_dd − K_da K_aa⁻¹ K_ad` **rebuilt
every form pass from the current (damaged) consistent tangent**, committed/reverted/serialized
`alphaCommit`, an initial-tangent path for `getInitialStiff`, an `update()` hook so `commitState`
captures the converged α under any algorithm, and the one-sided mode map
`M_i = sym[(j0/j)J0⁻ᵀE_i]` with cached centroid `easJ0inv/easJ0det`.

**Consequence:** the solid-shell EAS is an **adapt of a proven in-fork pattern, not net-new code.**
What is genuinely new (and stays in the budget):

1. the **mode set** — `E33`-directed thickness modes in the shell-natural frame (not the Wilson E9
   volumetric set), starting at 1 mode (`M33 ∝ ζ`, mapped through `(j0/j)J0⁻ᵀ`), extensible to the
   3–5-mode variants with bilinear in-plane variation if the plate-bending gates demand it;
2. the **ANS↔EAS interaction** — `B` is ANS-modified *before* the enhanced field is added; the
   orthogonality condition `∫M dV = 0` must hold against the ANS-modified strain measure (patch-test
   requirement);
3. the **softening-stability gate** (G6) — unchanged from the handout: enhanced modes can grow
   unstably when the material tangent goes indefinite post-peak; `LadrunoBrick`'s EAS has only been
   exercised through `LadrunoJ2`-class hardening and moderate `ASDConcrete3D` damage, not deep RC
   softening at shell slenderness.

The risk register keeps the *softening* clause of the old trap and retires the *net-new* clause.

## 4. Decisions

**D1 — Topology and DOFs (locked, from ADR 19).** 8 nodes, 3 translational DOF/node, nodes 1–4 =
bottom face, 5–8 = top face; the thickness direction ζ is defined by the node ordering (the 1-4 →
5-8 fibers). No director, no rotational DOFs, no drilling. Brick/contact/`SolidTransformation`
parity is the point.

**D2 — Constitutive seam: `NDMaterial` (3D), NOT a Section.** The element consumes any
`ThreeDimensional` `NDMaterial` per Gauss point — `LadrunoConcrete3D` (33017, the headline),
`LadrunoRCConcrete` 3D view (33015), `LadrunoJ2` (33011), `ASDConcrete3D`, elastic. This is the
brick contract, deliberately *not* the `LayeredShellFiberSection` contract: the entire value of the
element is the un-condensed triaxial state. Smeared reinforcement rides *inside* the RC kernels
(`rho_x/rho_y` params, per-element); discrete bars are `LadrunoEmbeddedRebar` (33005) in the solid
mesh — both already shipped.

**D3 — Integration: 2×2 in-plane × selectable `n_z` through thickness.**
`-nz $n` (default **2** = the elastic-exact minimum; RC guidance **5–9**) with
`-quadz {gauss|lobatto}` (default **lobatto** for `n_z ≥ 3` — face points capture surface crack
initiation; plain Gauss for `n_z = 2`). One material state per GP, `4·n_z` states total.
*Per-GP-layer material assignment (`-layers`, LayeredShell-style) is deliberately deferred* (§9 O2):
the honest solid answer to resolution is stacking elements through the thickness, and the
moment-curvature benchmark (G7) decides whether per-layer assignment ever earns its serialization
complexity.

**D4 — Locking cures (locked): ANS shear + ANS thickness + EAS-on-`E33`,** per the §2 table. All
three are *mandatory and always-on* — a solid-shell without them is just a bad brick, and "flags to
turn cures off" only multiplies the test matrix. The single debug escape is `-formulation std`
(plain 2×2×`n_z` displacement brick behavior, for A/B-ing locking in the gates), mirroring the
`LadrunoBrick` selector idiom.

**D5 — Geometry seam: `SolidTransformation`, v1 = `linear`.** The element routes through the
seam-3 globalize exactly like `LadrunoBrick` (the EAS assembly is already seam-routed there).
`-geom corot` ships **guarded**: the ADR 19 finding stands — a single nodal-cloud polar `R` is
ill-conditioned for `t ≪ L` (`cond(S) ~ (L/t)²`), so `update()` computes the eigenvalue ratio of
`S` and **loud-fails past a threshold** (no silent wrong answer; LEDGER_quirk on land). Thin
large-rotation use routes to `-geom finite` (+ `LogStrain`/native finite materials) or stays on the
director shell. Punching/bearing — the design load case — is a small-to-moderate-rotation regime;
`linear` covers the headline.

**D6 — Characteristic length: direction-split, not scalar.** The element's
`getCharacteristicLength()` must not hand a crack-band material the thickness dimension. v1
contract: **`lch = √(in-plane area at the GP)`** (in-plane projected), which is correct for the
dominant in-plane crack bands, plus a documented quirk that through-thickness softening bands (the
punching cone's vertical branch) regularize over `t/n_elem_z` — mesh guidance, not code, in v1.
A true directional `lch(n)` API is a cross-cutting fork upgrade (also wanted by ADR 19's √2-strut
residual) and stays out of scope here (§9 O4).

**D7 — Explicit path: real `criticalTimeStep()` + SMS pairing.** Unlike `ASDShellQ4` (returns −1;
LEDGER_quirk), this element **implements a per-element estimate**: `dt_cr ≈ min-altitude / c_d`
with the thickness fiber dominating at shell aspect ratios. That bound is brutal (`t` small) —
which is exactly the load case **selective mass scaling** (ADR 36, `CentralDifferenceSMS`/
`ExplicitBathe -sms`) was built for; the explicit demo gate (G9) pairs them. Lumped mass via the
HRZ path (ADR 35) like every fork solid.

**D8 — Serialization and state.** Committed per-GP material states (`4·n_z` copies) +
`alphaCommit[nEAS]` + the transformation state, schema-versioned from day one
(hard-checked version field — the RC-material lesson), with a **non-default-value round-trip
test** in the first PR that adds state (a default-valued round-trip cannot catch a dropped wire
slot).

**Rejected alternatives** (inherited from ADR 19 D1, kept for the record): new director shell
(zero physics over `ASDShellQ4`); 6-DOF drilling solid-shell; bolting σ33 onto `ASDShellQ4`;
full-integration unenhanced solid-shell (locks); reduced-integration 1-point solid-shell with
hourglass control only (the `ssp`/`Kstab` route — cannot represent the linear-in-ζ bending strain
that IS the shell's job; `Kstab` stays available as a *stabilization idea*, re-tuned, only if the
EAS route hits an unforeseen wall).

## 5. Formulation specification (what P5.1/P5.2 implement)

Natural coordinates `(ξ, η, ζ)`, ζ through thickness. Compatible Green-Lagrange (v1: small-strain)
covariant strain from the standard trilinear map, then three modifications **in this order**:

1. **ANS transverse shear** (Dvorkin-Bathe): replace covariant `E13`, `E23` with the midsurface
   edge-midpoint sampled interpolation —
   `E13(ξ,η) = ½(1−η)·E13|_B + ½(1+η)·E13|_D`, `E23(ξ,η) = ½(1−ξ)·E23|_A + ½(1+ξ)·E23|_C`,
   sampling points `A(−1,0,0), B(0,−1,0), C(1,0,0), D(0,1,0)` at ζ per GP.
2. **ANS thickness strain** (Betsch-Stein): replace covariant `E33` with the bilinear corner
   interpolation `E33(ξ,η) = Σ_i ¼(1+ξ_i ξ)(1+η_i η)·E33|_(ξ_i,η_i,0)`, corners `(±1,±1,0)`.
3. **EAS on `E33`** (Simo-Rifai, adapted `formEAStrue`): enhanced covariant strain
   `Ẽ33 = ζ·α₁` (v1 single mode; slots for the bilinear-in-plane extension
   `ζ·{ξ,η,ξη}·α₂..₄` reserved), pushed to the physical frame via the proven one-sided map
   `M = sym[(j0/j) J0⁻ᵀ E ζ]`, inner Newton `∫MᵀσdV = 0`, condensation each form pass, committed
   `alphaCommit`.

Covariant→physical transformation at each GP via the Jacobian (the standard solid-shell `T`
matrix); material called with the full physical 6-strain; engineering-shear convention at the
`NDMaterial` boundary (the fork-standard ×2, watch the LadrunoConcrete3D tensor-column quirk if
that material's kernel path is ever bypassed — the wrapper already handles it).

Mass: consistent trilinear-brick mass → HRZ lump on demand. Loads: body + surface pressure on the
two faces (punching needs face pressure).

## 6. Reuse map (verified on this branch)

| Asset | State | Role here |
|---|---|---|
| `LadrunoBrick::formEAStrue` + `alphaCommit` machinery | shipped, reviewed | the EAS pattern: inner Newton, condensation, commit/revert/serialize, initial-tangent path, `update()` hook. **Adapt** (new mode map + ANS interaction). |
| `SolidTransformation{,Corot,Finite}` | shipped (`SRC/element/solidTransformation/`) | seam-3 geometry; corot gains the `cond(S)` guard |
| `LadrunoConcrete3D` (33017) incl. `-implex`/`-eta` | shipped, byte-verified | the headline GP material; its robustness tiers are exactly what post-peak punching needs |
| `LadrunoRCConcrete` 3D view (33015) / `LadrunoRCFiniteStrain` (33018) | shipped | RC-with-smeared-steel GP material; finite view for `-geom finite` |
| `LadrunoEmbeddedRebar` (33005) | shipped | discrete bars through the solid-shell mesh |
| HRZ lumping (ADR 35), SMS (ADR 36), `CentralDifferenceLadruno` | shipped | the explicit path + the thickness-dt cure |
| `ASDShellQ4::getCharacteristicLength` idiom | shipped | the lch-latch pattern the D6 contract mirrors |
| Registration surface (4 core files + banner + manifest + ledger) | documented in [[project memory / LEDGER]] | ELE-flavored repeat of the LadrunoConcrete3D checklist |

## 7. Acceptance gates

| # | Gate | Phase | Pass criterion |
|---|---|---|---|
| G1 | **Patch tests** — membrane, bending, and **thickness** (uniform σ33 through a distorted mesh) | P5.1 | exact to solver tol; the thickness patch is the one that catches a broken `∫M dV = 0` |
| G2 | **Slenderness locking sweep** — clamped/SS plate, `t/L ∈ {10⁻¹, 10⁻², 10⁻³}`, 1 element through thickness | P5.1 | center deflection → analytic thin-plate value, no shear/Poisson-locking degradation; `-formulation std` A/B shows the disease |
| G3 | **Pinched cylinder + Scordelis-Lo** (the classic obstacle course) | P5.1 | reference values within the standard ~1–2% at the canonical meshes |
| G4 | **Curved/trapezoidal-mesh bending** (Betsch-Stein test) | P5.1 | no curvature-thickness locking vs the flat-mesh result |
| G5 | **Softening snap-back objectivity** — notched tension band of solid-shells, `LadrunoConcrete3D` `-autoRegularization`, 3 mesh densities | P5.2 | dissipated energy mesh-objective (the Bažant-bar pattern from RC Phase 3b-struct, adapted); converges under arc-length or quasi-static explicit |
| G6 | **EAS mode stability under softening** — single element driven deep post-peak; monitor α | P5.2 | α bounded and convergent; a failure here triggers the documented fallback ladder (freeze-α-post-peak à la IMPL-EX → mode removal → Kstab-stabilized variant) |
| G7 | **Moment-curvature benchmark vs `LayeredShellFiberSection`** — same RC wall section, solid-shell stack vs director shell | P5.3 | quantified agreement/deviation published in the guide; **gates any flexural claim** (D1) |
| G8 | **Punching headline** — slab-column benchmark (literature specimen TBD at P5.3, e.g. a Guandalini/Muttoni-series slab), `LadrunoConcrete3D` | P5.3 | punching capacity within the scatter band of the reference; the director-shell xfail documented and discharged |
| G9 | **Explicit path** — G8 or G5 geometry under `CentralDifferenceLadruno` + SMS; real `criticalTimeStep()` | P5.3 | completes; dt_cr estimate validated against a bisected stable dt; SMS speedup quantified |
| G10 | **6-DOF↔3-DOF connection recipe** — `ASDShellQ4` panel butted to a solid-shell patch, moment across the seam | P5.3 | documented, *validated* recipe (rigid-link/`equalDOF`/`LadrunoTie -hermite` candidates); moment continuity error quantified |

Zone-A pytest throughout; oracle-first where a numpy reference is buildable (G1/G2 closed forms,
G5's energy bookkeeping), element-battery otherwise. Every PR: ledgers + manifest row + header
stamps + schema round-trip (non-default values).

## 8. Phasing

- **P5.1 — skeleton + anti-locking (elastic). [BUILT 2026-07-02]** Element class, parser
  (`element LadrunoSolidShell $tag $n1..$n8 $matTag <-nz $n> <-quadz gauss|lobatto>
  <-formulation ans|std> <-geom linear>`), ANS shear + ANS thickness, **state-dependent EAS**
  (the formEAStrue adapt shipped whole — live-σ inner Newton + committed serialized α — since it
  was cheaper than a throwaway initial-tangent stub; the *softening-stability* gate stays P5.2),
  registration surface, G1/G2/G4 (11 gates green; G3 pinched-cylinder/Scordelis-Lo landed
  2026-07-06 in `tests/test_ladrunoSolidShell_benchmarks.py` — SL 0.989 @16x16, PC 0.927 @16x16
  == canonical MITC4, std controls lock 8-11x). Gate discoveries recorded in LEDGER_quirks: the 1-element-thick interior patch must
  use a traction-consistent field (`σ·e_z = 0`), and the site `.pth` boot module pins one
  worktree's pyd (PMI_RANK bypass). Rho rides the material (`getRho`), no `-rho` flag; `-lch`
  is the material's own override (the element supplies √A_midsurface via
  `getCharacteristicLength`, verified quantitatively against the LadrunoConcrete3D crack-band
  latch). *Defines ELE_TAG 33020; ledger row in progress.*
- **P5.2 — state-dependent EAS + nonlinear materials. [GATES LANDED 2026-07-06, test-only]**
  The `formEAStrue` adaptation itself shipped with P5.1; P5.2 landed as
  `tests/test_ladrunoSolidShell_softening.py`: hosting parity (NDTest replay) +
  `LadrunoRCConcrete` smoke, G5 (0.9% energy spread across a 4× in-plane-square mesh range with
  `-autoRegularization`+`-implex`; fixed-lch control at the 0.25× band-volume ratio), G6 both arms
  (bending deep post-peak + indefinite-Kaa freeze path) — α bounded/settling, the P5.1 freeze-α
  guard suffices, the fallback ladder stays dormant. The heavy P5.1 adversarial review already
  covered the ANS↔EAS/formANS machinery; P5.2 added no source, so the review burden was the gate
  battery itself (three new LEDGER_quirks recorded).
- **P5.3 — validation + integration.** G7–G10: benchmarks, connection recipe, explicit demo,
  user guide (`LadrunoSolidShell_guide.md`), banner line.
  **[G9 LANDED 2026-07-06, #504]** `-lumped` + explicit/SMS battery (see §7 and LEDGER_quirks).
  **[G7 LANDED 2026-07-07, test-only]** `tests/test_ladrunoSolidShell_flexure.py` — the D1
  flexural gate PASSED with the deviation quantified: same RC wall section (t=200,
  ρ_tot=0.6% at z=±80) as ASDShellQ4+LayeredShell/PlateRebar vs a solid-shell
  through-thickness stack (3D RCConcrete at the GPs, Steel01 interface trusses); elastic
  anchor exact both arms (layered midpoint-rule Σh³/12 deficit predicted, ≈2% at 5 layers);
  classic M–κ (dummy-node `equationConstraint` section driver, N held exactly): dM ≤ 2.2% of
  Mmax through deep cracking + steel yield at N=0 and at N=−0.1·fc·Ag, ε₀-migration ≤ 5%;
  restrained (ε₀=0) resultants on a 6-element stack: dM 2.1%, dN 16% of ft·t·b. The IMPL-EX
  reporting quirk (ASDShellQ4 = post-commit state ⇒ implex invisible to prescribed rigs;
  trial-state hosts show it) is pinned discriminatingly — parity gates run implicit-both.
  Remaining: G8 (blocked on the O5 specimen decision), G10, guide, banner.
- **P5.4 (optional, demand-gated) — kinematic upgrades.** `-geom corot` beyond the guard
  (shell-aware corotation), `-geom finite` validation with `LadrunoRCFiniteStrain`; per-layer
  material assignment is OFF the list (G7/O2 decided against it — see §9).

Effort: **L overall** (unchanged from the handout) but P5.2 shrinks materially via the
`formEAStrue` adapt (§3). Multi-PR, one phase ≈ one PR, all `--base ladruno`.

## 9. Open questions (decide at the flagged phase)

- **O1 (P5.2):** EAS mode count — does the single `ζ·α₁` mode clear G2/G4 at `ν = 0.2`, or do the
  bilinear in-plane thickness modes earn their 3 extra α? Start at 1, let the gates decide.
- **O2 (P5.3, G7): DECIDED 2026-07-07 — element stacking; per-GP-layer assignment is not
  needed.** Measured on the G7 restrained-section rig: with ONE core element through
  [−80, 80] the resultant-N parity error vs the layered reference is ≈62% of ft·t·b deep
  post-crack (a single EAS ζ-mode + the nodal-uz mean cannot relax the kinked cracked-state
  σ₃₃ profile), while a 6-element stack (nz=3 lobatto each) converges it to 16% (≈4%
  relative); M(κ) is forgiving either way (≤3% vs ≤2%). Stacking both resolves the
  through-thickness state AND is where interface rebar naturally attaches — per-layer
  assignment would add serialization complexity without addressing the σ₃₃ kinematics.
- **O3 (backlog):** in-plane EAS/incompatible modes for distortion insensitivity — only if a G3
  distorted-mesh variant shows unacceptable in-plane stiffness.
- **O4 (cross-cutting, own ADR if pursued):** directional `lch(n)` API for crack-band materials —
  wanted by this element's through-thickness bands *and* ADR 19's √2-strut residual.
- **O5 (P5.3, G8):** the punching specimen — pick one with published load-rotation curves and a
  well-documented failure cone; digitized data availability decides (same constraint as the
  Tran–Wallace item V).

## 10. Risk register

| Risk | Mitigation |
|---|---|
| EAS α-growth instability in deep softening (the surviving half of the ADR 19 trap) | G6 gate + fallback ladder (freeze-α, mode removal, re-tuned damage-scaled Kstab as last resort) |
| ANS↔EAS orthogonality broken → thickness patch-test failure | G1 thickness patch is first-class, runs from P5.1 |
| Thin-shell corot ill-conditioning | D5 `cond(S)` loud-fail guard; punching regime doesn't need corot |
| Thickness-bound explicit dt_cr makes explicit useless | D7 SMS pairing, G9 quantifies |
| Scope creep toward "co-equal flexural host" | D1 lock + G7 must pass *before* any flexural claim appears in docs |
| Indefinite/non-symmetric softening tangents stall implicit punching runs | inherit the material toolbox: `-implex`, `-eta`, quasi-static explicit (all shipped on 33017) |
| Fixed-projector secant + eigenprojector degeneracy at equibiaxial slab states (ADR 19 cross-cutting risk) | states are *material*-side and already handled there; element adds no new spectral machinery |

---

*Build-control: this ADR adds no code — no ledger/manifest/banner change. The 33020 manifest row,
LEDGER_implementations row, and header-stamped sources land with P5.1. File map (from ADR 19):
`SRC/element/ladrunoSolidShell/LadrunoSolidShell.{h,cpp}` + `OPS_LadrunoSolidShell.cpp` +
`CMakeLists.txt`; tests `tests/test_ladrunoSolidShell_{element,locking,softening}.py`; oracle
additions in `tests/_testbed/` as needed.*
