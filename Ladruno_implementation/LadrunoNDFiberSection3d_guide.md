---
title: LadrunoNDFiberSection3d and LadrunoUniaxialFiber — shear-profile NDFiber section and bar fibres
project: Ladruno
tags:
  - section
  - concrete
  - shear
---

# `section LadrunoNDFiber` and `nDMaterial LadrunoUniaxialFiber` (WP-245)

MCFT shear-critical fibre-beam program (#1131), rcbf-3. Classes `LadrunoNDFiberSection3d`
(SEC 33001) and `LadrunoUniaxialFiber` (nD 33030). Both are in openseespy and in the classic
Tcl exe.

## 1. What it is for

A shear-critical RC column modelled with fibre sections needs three things that the stock
`NDFiber` section does not give together:

1. **A shear-strain profile.** Stock `NDFiberSection3d` imposes the same shear strain on every
   fibre (`γxy = √α γy`), so the cover fibres at the faces carry full shear. Model C (the program's
   oracle, `research/06_modelC_oracle_spec.md` §1.2) imposes a parabola that vanishes at the
   faces.
2. **Model C's shear resultant**, `V = Σ τ A`, which Model C's calibration (k_β, c_vci, Gf) was
   done with.
3. **Steel bars in the same section.** `NDFiber` sections take only nD materials, so a
   `UniaxialMaterial` bar cannot sit next to the concrete.

The concrete fibre is `LadrunoRCConcrete` with its BeamFiber view (WP-243,
[[LadrunoRCConcrete_guide]] §12.1).

## 2. Syntax

```tcl
section LadrunoNDFiber $tag -shearProfile uniform   [-alpha $a]                      [-noCentroid] { fibres }
section LadrunoNDFiber $tag -shearProfile parabolic -shearDepth $hy $hz
                            [-shearResultant modelC|conjugate]                       [-noCentroid] { fibres }

nDMaterial LadrunoUniaxialFiber $tag $uniTag [-rho $rho]
```

```python
ops.section("LadrunoNDFiber", 1, "-shearProfile", "parabolic", "-shearDepth", h, b)
ops.fiber(y, z, A, concreteTag)        # stock fibre / patch / layer commands fill it
ops.nDMaterial("LadrunoUniaxialFiber", 6, steelTag)
ops.fiber(yb, zb, Abar, 6)
```

| Token | Values | Rule |
|---|---|---|
| `-shearProfile` | `uniform` \| `parabolic` | **Required.** There is no silent default. |
| `-shearDepth` | `hy hz` (> 0) | The section's extent along y and z. Required for `parabolic`, refused for `uniform`. |
| `-shearResultant` | `modelC` (default) \| `conjugate` | `parabolic` only. |
| `-alpha` | `a` (> 0) | The stock shear shape factor. `uniform` only. |
| `-noCentroid` | — | As stock: fibre coordinates are used as given. |

The parser fails closed: an unknown token, a missing value or a wrong combination makes the
command fail with the token named (WP-167). There is no `-GJ`/`-torsion`: an ND section carries
torsion through its fibres. 3D models only (`ndm 3`); `ndm 2` is refused.

## 3. Kinematics and resultants

A fibre at `(y, z)` from the centroid takes the section deformation
`e = (ε0, κz, κy, γy, γz, θ)` (stock order P, MZ, MY, VY, VZ, T):

```
εxx = ε0 − y κz + z κy
γxy = φy(y) γy − z θ
γxz = φz(z) γz + y θ
```

- `uniform`: `φ = √α` (stock).
- `parabolic`: `φy = 1.5 (1 − (2y/hy)²)`, `φz = 1.5 (1 − (2z/hz)²)`. It is zero at the faces and 1.5 at
  mid-depth, evaluated at the fibre centroid with no renormalisation, as Model C does: the mean
  over 40 equal strips is `1 + 1/3200`.

The resultant uses a force map `b`, with `s = Σ A bᵀσ` and `K = Σ A bᵀ D a`:

| Rule | Shear rows of `b` | Elastic rectangle, `K_VV` | Tangent |
|---|---|---|---|
| uniform | `√α` (= `a`) | `α G A` | symmetric if D is |
| parabolic `modelC` | `1` | `G A (1 + 1/(2n²))` → `G A` (k = 1) | **unsymmetric** |
| parabolic `conjugate` | `φ` (= `a`) | `2.25 G A (1 − 2μ₂ + μ₄)` → `6/5 G A` | symmetric if D is |

(`n` is the number of equal strips; `μ₂, μ₄` are the midpoint means of `u²` and `u⁴`.)

> [!warning] Neither rule gives Timoshenko's 5/6
> k = 5/6 comes from a parabolic shear *stress* with the energy of the exact solution. A fibre
> section imposes the *strain*. The work-conjugate force for a parabolic strain is
> `G A mean(φ²) = 6/5 G A`: the kinematic constraint stiffens the section. Model C's `V = Σ τ A`
> is not work-conjugate and gives `G A`. The rcbf-3 plan's "k = 5/6" gate was therefore replaced
> by the closed forms in the table above (decided 2026-11-22).

`modelC` is the default because the oracle and its calibration use it. Its tangent is
unsymmetric (`∂V/∂γ = Σ D11 φ A`, while `∂(work)/∂γ` would carry `φ²`), so use an unsymmetric
solver. On the desktop that is plain `system Pardiso`: with no flag it stores the full matrix and
factors it as real unsymmetric (MKL mtype 11). `UmfPack`, `SparseGEN`, `Mumps` and `FullGeneral` also
work. **Do not use `system Pardiso -symmetric` or `-spd` with `modelC`.** They store only the upper
triangle, so the solve uses a symmetrised tangent. Newton then still converges, but on a wrong
Jacobian. The fork's Pardiso warns ("an assembled element matrix is UNSYMMETRIC") at the first
assembly where that happens. On a cracking RC cantilever (forceBeamColumn, 40 strips, axial load,
displacement control), plain `Pardiso` and `UmfPack` give identical histories (89 vs 90 Newton
iterations over 15 steps). `-symmetric` warned with an asymmetry of 1e-3 of the matrix scale. The
`conjugate` resultant is symmetric whenever the fibre tangent is.

Fail-closed geometry: a fibre centroid outside `|y| ≤ hy/2`, `|z| ≤ hz/2` (measured from the
centroid) would get a negative `φ`. Such a section refuses its first state determination, so
`analyze` fails with the fibre named. DDM sensitivities (`getStressResultantSensitivity` etc.)
use the stock √α kinematics, so they are refused under the parabolic profile. `setParameter alpha`
is refused.

## 4. The bar fibre

`LadrunoUniaxialFiber` is the bar. The fibre's axial strain drives the uniaxial law, and the
fibre has **zero shear stress and zero shear stiffness** (no dowel action, as in Model C §1.4):

| View | Strain | Stress | Tangent |
|---|---|---|---|
| `BeamFiber` (order 3) | `(ε, γxy, γxz)` | `(σ(ε), 0, 0)` | `diag(Eₜ, 0, 0)` |
| `BeamFiber2d` (order 2) | `(ε, γ)` | `(σ(ε), 0)` | `diag(Eₜ, 0)` |

It works in the stock `NDFiber` sections as well. Every other view (`ThreeDimensional`,
`PlaneStress`, `PlateFiber`, …) is refused. Responses: `stress`/`strain`/`tangent` as for any nD
material; any other response name is forwarded to the uniaxial material. `-rho` overrides the
uniaxial material's density.

Model C subtracts each bar's area from the nearest concrete strip. Do the same when you mesh:
the concrete fibre's area should be net of the bars it contains.

## 5. Model C layout (the rcbf-0 oracle's specimens)

```python
h = b = 457.2; n = 40; dy = h / n; cover_tie = ...
ops.section("LadrunoNDFiber", 1, "-shearProfile", "parabolic", "-shearDepth", h, b)
for i in range(n):
    y = -h/2 + dy*(i + 0.5)
    core = abs(y) <= h/2 - cover_tie
    ops.fiber(y, 0.0, b*dy - bars_in_strip(i), coreTag if core else coverTag)   # ties: core only
for (yb, Ab) in bars:
    ops.fiber(yb, 0.0, Ab, barTag)                                             # LadrunoUniaxialFiber
```

Here `coreTag` is a `LadrunoRCConcrete` with `-ties ρt 0 fyt Es` (the transverse ties, y only),
and `coverTag` is the same concrete without ties.

## 6. Verification

- `tests/_testbed/rc_section_gpp.cpp` (g++, header-only; wrapper
  `test_ladrunoNDFiberSection3d_cpp.py`):

  | Gate | Check | Result |
  |---|---|---|
  | S1 | elastic closed forms | 4e-15 |
  | S2 | conjugate k → 6/5 | — |
  | S3 | conjugate symmetric; modelC unsymmetric by design | 5e-17 |
  | S4 | uniform = the stock unrolled formulas | 2e-16 |
  | S5 | chain rule by central FD, both resultants | 1e-12 |
  | S6 | Model C assembly vs strips driven directly on a cracking RC path with commits (40/40) | 2e-16 |
  | S7 | Vy/Vz invariance | 1e-9 |
  | S8 | band check | — |

- `tests/test_ladrunoNDFiberSection3d.py` (Zone-A, `zeroLengthSection`):
  - elastic closed forms;
  - Model C assembly against 40 independently driven stock NDFiber strips;
  - uniform = stock `NDFiber`;
  - bar fibre: axial force, no shear, and Steel01 past yield;
  - parser, `ndm 2` and out-of-band refusals;
  - a `forceBeamColumn` cantilever whose tip deflection is `PL³/3EI + PL/(k G A)` with the section's own k, to 1e-12 (the unsymmetric modelC tangent on UmfPack);
  - a two-round database round trip (elastic concrete in shear under the conjugate resultant, plus yielded and unloaded Steel01 bars).
  - the classic Tcl exe gives the same section.

## 7. Scope

- 3D only. A 2D section (`P, MZ, VY`) can follow if rcbf-4 needs it; the 3D section reproduces
  2D in-plane loading exactly.
- Not wired into `SRC/runtime` (OpenSeesRT), which this fork does not build.
- The comparison with the rcbf-0 oracle's *strip law* (rotating MCFT with Popovics) is a band,
  not a gate: the two concrete laws differ (WP-243 note). The section-level gate is S6: same
  fibres, Model C's assembly.
