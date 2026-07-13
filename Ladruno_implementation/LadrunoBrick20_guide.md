# LadrunoBrick20 — 20-node serendipity quadratic hex (user guide)

Second-order displacement brick (H20, ELE_TAG 33018, ADR 72). Quadratic
accuracy for smooth implicit problems — stress concentrations, curved
geometry, bending with one element through the thickness — at the standard
serendipity cost (60 DOF, 27 Gauss points).

```tcl
element LadrunoBrick20 $tag $n1 ... $n20 $matTag \
    <-formulation std> <-b $bx $by $bz> <-damp $dampTag> <-lumped>
```

```python
ops.element('LadrunoBrick20', tag, *nodes20, matTag,
            '-formulation', 'std', '-b', bx, by, bz)
```

- **Node order** = upstream `20NodeBrick` / Abaqus C3D20 / gmsh-via-apeGmsh:
  corners 1-8 (bottom ring CCW, then top ring CCW), mid-edges 9-12 (bottom
  ring), 13-16 (top ring), 17-20 (vertical). Wrong ordering dies at
  `setDomain` with a negative-Jacobian error.
- **Material**: any 3D `NDMaterial`, cloned per GP (27 copies),
  `getCopy("ThreeDimensional")`.

## v1 scope (P1) — what works and what is reserved

| Surface | Status |
|---|---|
| `-formulation std` (27-pt full integration) | **live** — reduces to upstream `20NodeBrick` at ~1e-15 on the same mesh |
| `-formulation uri` (2×2×2 reduced, Barlow points) | reserved — **P2** (parser says so) |
| `-lumped` (HRZ mass) | reserved — **P3**; accepted at parse, `getMass` errors once and falls back to consistent |
| `-hourglass` | **never** — hard error. H20@2×2×2 spurious modes are non-communicable in meshes (Abaqus applies no control to C3D20R either); a control would mask, not fix, single-element-stack pathologies |
| Explicit dynamics | unsupported until P3 (needs the HRZ lump; row-sum lumping of H20 gives **negative corner masses −M/8** and is forbidden) |
| `-geom finite` | not offered on this element (ADR 72 §2: quadratic + finite strain is a documented anti-pattern; use the H8 family) |
| Gmsh recorder (`recorder gmsh`) | **not wired** — the recorder's 20-node MSH mapping is defective for ALL 20-node bricks (wrong element type 12-vs-17 + node-order permutation missing); use VTK/PVD/LadrunoRecorder (all verified) until the GmshRecorder fix lands |
| Embedded hosts (`LadrunoEmbeddedNode`/`Rebar` inside H20) | deferred — `getInterpolationWeights/Gradients` land in P4 (ADR §6); embed in the H8 family meanwhile |

## Response surface

`forces`, `stiff`, `material $gp ...`, `stresses`/`strains` (per-GP, 27×6,
brcshl GP order — pairs 1:1 with the recorder rule `Hexahedron_GaussLegendre_3`
and with upstream `20NodeBrick` GP numbering), `stress3D6`/`strain3D6`
(volume-averaged), `charLength` (= ∛V).

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
`tests/test_ladrunoBrick20_element.py` (element battery, 18 gates) +
`tests/test_hex20_kernel_cpp.py` (kernel-vs-sympy-oracle gate).
