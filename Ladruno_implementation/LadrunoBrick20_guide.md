# LadrunoBrick20 — 20-node serendipity quadratic hex (user guide)

Second-order displacement brick (H20, ELE_TAG 33018, ADR 72). Quadratic
accuracy for smooth implicit problems — stress concentrations, curved
geometry, bending with one element through the thickness — at the standard
serendipity cost (60 DOF, 27 Gauss points).

```tcl
element LadrunoBrick20 $tag $n1 ... $n20 $matTag \
    <-formulation std|uri> <-b $bx $by $bz> <-damp $dampTag> <-lumped>
```

```python
ops.element('LadrunoBrick20', tag, *nodes20, matTag,
            '-formulation', 'uri', '-b', bx, by, bz)
```

- **Node order** = upstream `20NodeBrick` / Abaqus C3D20 / gmsh-via-apeGmsh:
  corners 1-8 (bottom ring CCW, then top ring CCW), mid-edges 9-12 (bottom
  ring), 13-16 (top ring), 17-20 (vertical). Wrong ordering dies at
  `setDomain` with a negative-Jacobian error.
- **Material**: any 3D `NDMaterial`, cloned per GP (27 copies under `std`,
  8 under `uri`), `getCopy("ThreeDimensional")`.

## Scope — what works and what is reserved

| Surface | Status |
|---|---|
| `-formulation std` (27-pt full integration) | **live (P1)** — reduces to upstream `20NodeBrick` at ~1e-15 on the same mesh |
| `-formulation uri` (2×2×2 reduced, Barlow points; alias `reduced`) | **live (P2)** — the C3D20R analog; anchored to the P0 sympy oracle at ≤1e-12; mass/body-force/volume stay 27-pt |
| `-lumped` (HRZ mass) | **live (P3)** — HRZ diagonal lump of the 27-pt consistent mass, positive by construction; cube fractions pinned at 7/248 (corners) / 2/31 (mid-edges) through the element path |
| `-hourglass` | **never** — hard error. H20@2×2×2 spurious modes are non-communicable in meshes (Abaqus applies no control to C3D20R either); a control would mask, not fix, single-element-stack pathologies |
| Explicit dynamics | **permitted since P3, honestly discouraged** — see the dynamics section below |
| `-geom finite` | not offered on this element (ADR 72 §2: quadratic + finite strain is a documented anti-pattern; use the H8 family) |
| Gmsh recorder (`recorder gmsh`) | **wired since PR #564** — MSH type 17 + write permutation, for this element and upstream `Twenty_Node_Brick` |
| Embedded hosts (`LadrunoEmbeddedNode`/`Rebar` inside H20) | deferred — `getInterpolationWeights/Gradients` land in P4 (ADR §6); embed in the H8 family meanwhile |

## std vs uri — the selection table (ADR 72 §3.2/§3.3)

`uri` integrates stiffness and stress at the 8 **Barlow points** (2×2×2),
where quadratic-hex stresses are superconvergent — it is both ~3× cheaper per
tangent AND the *premium* stress product. The price: a lone `uri` element has
6 zero-energy modes. They cannot communicate in any assembly ≥ 2 elements
thick (pinned by the P2 battery), but they CAN propagate along a
single-element-thick stack, and they poison eigen extraction.

| Situation | Use |
|---|---|
| Smooth production meshes ≥ 2 elements thick (statics, implicit dynamics) | **`uri`** — cheaper, superconvergent GP stresses |
| Eigenvalue / modal extraction | **`std`** (spurious modes appear as garbage low modes under `uri`) |
| Single-element-thick members (walls, slabs modeled 1-hex thick), single stacks | **`std`** — the uri modes propagate there (demonstrated in the battery) |
| Point loads / point-restrained corners, stiff-on-soft (footing-on-soil) contrasts | **`std`** (point actions excite the modes) |
| ν → 0.5 (near-incompressible, confined plastic flow) | **`uri`** relieves the volumetric lock that `std` suffers (ADR §3.4); for truly incompressible demand use the H8 `bbar` family |
| Minimally-restrained single elements (unit tests, patch rigs) | **`std`** — a determinately-restrained lone `uri` element is a mechanism |

Measured numbers from the P2 battery (this build): coarse bending (5×1×1
cantilever, ν=0, tip/beam-theory) `uri` 0.999 / `std` 0.996 where the H8
`LadrunoBrick std` shear-locks at 0.333 (H8 `eas`: 0.989); near-incompressible
(4×1×1, ν=0.4999) `uri` 0.891 vs `std` 0.749 — uri *halves* the locking error
but is not a mixed element ("relieves … after significant straining it
re-locks", ADR §3.4): refine, or use the H8 `bbar`/`eas` family for true
incompressibility (at ν=0.3 the two formulations agree within 0.012 — no
formulation gap at moderate ν, both →~0.99 on refinement). Barlow
superconvergence: GP σ_xx RMS error 3.7e-10 (`uri` — exact for a linear
bending field) vs 3.9e-2 (`std`). Tangent-formation cost ratio uri/std ≈ 0.63
measured through the full analyze round-trip (pure assembly ≈ 0.3–0.4; the
round-trip's fixed overhead dilutes it). See
`tests/test_ladrunoBrick20_uri.py` S4/S5/S6/S9 logs.

## Dynamics, `-lumped`, and explicit runs (ADR 72 §3.5/§3.6 — P3)

**Mass models.** Default = the 27-pt **consistent** mass (formulation-
independent; eigenfrequencies converge from above at 4th-order-class rates,
Hughes §7.3.2). `-lumped` = the **HRZ** diagonal lump of that consistent block
via the shared `Ladruno::hrzLump` (ADR 35). Why HRZ is the *only* lumping
offered: plain row-sum lumping of an H20 gives **negative corner masses,
−M/8 each** (mid-edges +M/6) — negative lumped masses produce exponentially
growing modes, which is exactly LS-DYNA's 1979 "numerical noise from ad hoc
mass lumping" verdict on quadratic solids. HRZ scales the strictly-positive
consistent diagonal to conserve directional mass, so every entry is positive
by construction: on the unit cube, corners carry 7/248 and mid-edges 2/31 of
the element mass (pinned ≤1e-14 through the element path in the P3 battery).

**HRZ accuracy caveat** (Cook p. 373): no positive *optimal* lumping exists
for serendipity elements — HRZ trades the optimal-lumping convergence rate
for guaranteed positivity. Measured (fixed-free axial rod, exact analytic):
consistent f₁/analytic = 1.00026 (nx=2) → 1.000003 (nx=6) from **above**;
HRZ = 0.9536 → 0.9947 from **below**. Both bracket the analytic and both
converge; the HRZ error is simply larger at a given mesh. Eigen/modal work
should stay on the consistent mass (and `std` — see the selection table).

**Explicit runs are permitted — and discouraged** (ADR 72 §3.6). The HRZ mass
is positive, `criticalTimeStep()` returns the exact per-element 60-DOF
algebraic pencil (ADR 65 — verified ≤1e-10 vs an independent numpy
generalized-eig oracle, measured ~3e-16, cube and distorted, std and uri),
and the P3 battery runs a 10-element wave bar 2500 steps at 0.9×pencil under
`CentralDifferenceLadruno` with ADR-69 energy closure (drift 0.28% of peak,
residual ≤3.8%, and NO hourglass channel — this element has no hourglass
machinery at all). But the Δt economics are poor: the measured **dt_cr ratio
vs an equal-node-spacing `LadrunoBrick` H8 mesh is ≈0.50** — *below* the 1-D
rod-theory 2/√6 ≈ 0.82 ballpark, because the tiny 3-D HRZ corner masses
(7/248 vs the rod's order-1/6 shares) push ω_max higher still. Half the
stable step per node, ~3× the per-element cost, and higher-order elements
add wavefront noise (Cook §13.14). **For production explicit work use
`LadrunoBrick -formulation uri -lumped` or the Bézier line** (positive
Bernstein lumping) — this element's explicit path exists to be *correct*,
not fast.

## Response surface

`forces`, `stiff`, `material $gp ...`, `stresses`/`strains` (per-GP: 27×6
brcshl order under `std` — pairs 1:1 with the recorder rule
`Hexahedron_GaussLegendre_3` and upstream `20NodeBrick` GP numbering; 8×6
z-fastest lexicographic under `uri` — the `Hexahedron_GaussLegendre_2` /
Barlow set), `stress3D6`/`strain3D6` (volume-averaged), `charLength` (= ∛V,
27-pt volume in both formulations). The element self-describes to the
LadrunoRecorder via `basisInfo`/`integrationPoints`, so `.ladruno` GP
metadata (NUM_GP, GP_PARAM, GLOBAL_GP_COORDS) follows the formulation
automatically — mixed std/uri meshes are fine.

## Caveats (read before production use)

1. **Consistent nodal loads look wrong — they aren't.** Uniform pressure /
   self-weight on a quadratic face puts **near-zero or sign-flipped forces on
   corner nodes** (H20 volume shares: corners −V/8·ρb, mid-edges +V/6·ρb).
   Reactions and totals are exact; per-node values are not intuitive. This is
   also why NTS contact on quadratic faces is deferred (ADR 72 §5).
2. **Softening / crack-band materials** (ASDConcrete3D, LadrunoConcrete3D):
   the element emits a **one-time advisory** at `setDomain`. lch = ∛V is
   handed to the material consistently, but crack-band regularization theory
   on quadratic strain fields is shaky and localization wants first-order
   elements (Abaqus and the textbooks agree). Prefer `LadrunoBrick` (H8) for
   softening-dominated runs; the advisory does not stop the analysis.
3. **Large plastic strain / severe distortion**: quadratic elements lose
   their accuracy edge on C⁰ material response (rates cap at the linear-
   element order) and distort worse. They pay off in the smooth/elastic
   regime.
4. **Cost**: ~22× the H8 stiffness-formation cost per element (27 GPs × 60
   DOF); fewer elements needed for the same accuracy on smooth problems —
   that trade is the whole point.

Full decision record: `72_ladruno_second_order_brick_adr.md`. Tests:
`tests/test_ladrunoBrick20_element.py` (P1 std battery) +
`tests/test_ladrunoBrick20_uri.py` (P2 uri contract battery, specs
`_adr72_p2_test_specs.md`) + `tests/test_ladrunoBrick20_dynamics.py` (P3
HRZ / explicit / energy / betaK-clobber battery) +
`tests/test_hex20_kernel_cpp.py` (kernel-vs-sympy-oracle gate).
