---
title: LadrunoMCFTStrip and LadrunoMCFTSteel — the Model C MCFT strip fibre and its bar law
project: Ladruno
tags:
  - material
  - concrete
  - shear
---

# `nDMaterial LadrunoMCFTStrip` and `uniaxialMaterial LadrunoMCFTSteel` (WP-254)

MCFT shear-critical fibre-beam program (#1131). Classes `LadrunoMCFTStrip` (nD 33031) and
`LadrunoMCFTSteel` (uniaxial 33004). Both are in openseespy and in the classic Tcl exe.

## 1. What it is for

A fibre column whose cyclic shear-strength loss comes from the section mechanics, after the
"Model C" fibre beam (a force-based Timoshenko element with MCFT strips, calibrated on the
Sezen & Moehle 2002 and Lynn et al. 1996 columns). The program's first route ran the fork's
general 3D concrete (`LadrunoRCConcrete`, BeamFiber view: three condensed lateral strains, a
damaged-secant tangent) in the same section; it was accurate where it ran but rarely got past
the shear failure (three lateral unknowns, discontinuous root choices, negative axial stiffness
after the failure). This material is the simpler strip of Model C itself: **one** lateral
unknown, the strip law and root rule exactly as Model C has them.

```
LadrunoForceBeamColumn3d (HingeRadau Lp Lp, Linear transf)
 └─ section LadrunoNDFiber -shearProfile parabolic -shearResultant modelC -noCentroid   (WP-245)
     ├─ strips:  nDMaterial LadrunoMCFTStrip   (core: -ties; cover: none)
     └─ bars:    nDMaterial LadrunoUniaxialFiber -> uniaxialMaterial LadrunoMCFTSteel
```

## 2. Syntax

```
nDMaterial LadrunoMCFTStrip $tag $fc $fp $ep $e20 -crackSpacing $smx $smy
    [-ties $rho $fyt $fut] [-tieSteel $Es $esh $esu]
    [-agg $a] [-betaC $C] [-cvci $c] [-Gxz $G] [-tangentStep $h]

uniaxialMaterial LadrunoMCFTSteel $tag $fy $fu [-Es $Es] [-esh $esh] [-esu $esu]
```

| Input | Meaning | Default |
|---|---|---|
| `$fc` | f'c (> 0); `Ec = 4700 √fc`, `fcr = 0.33 √fc`, `ecr = fcr/Ec` (MPa) | — |
| `$fp $ep` | compression peak stress and strain, magnitudes (Mander core / f'c, 0.002 cover) | — |
| `$e20` | strain at 0.2 fp (> ep). Model C regularises it: `e20 = Gf/(0.6 fp Lp) − 0.8 fp/Ec + ep`, `Gf = 0.5·2fc` (cover), `0.5·1.7·2fc` (core) N/mm, `Lp` the hinge length | — |
| `-crackSpacing $smx $smy` | crack-spacing parameters (Model C: bar-layer spacing, tie spacing) | required |
| `-ties $rho $fyt $fut` | smeared tie ratio across the strip and its steel; omit for cover strips | none |
| `-tieSteel $Es $esh $esu` | tie steel backbone | 200000 0.008 0.12 |
| `-agg` | maximum aggregate size in the `v_ci` formula | 19 |
| `-betaC` | `C` in `β = min(1, 1/(0.8 + C e1))` | 170 |
| `-cvci` | `v_ci,max = c √fc / (0.31 + 24 w/(a + 16))` | 0.18 |
| `-Gxz` | elastic out-of-plane shear modulus (3D view) | `Ec/2` |
| `-tangentStep` | forward-difference step of the tangent | 1e-7 |

Units: the formulas carry Model C's MPa/mm constants (`4700 √fc`, `0.33 √fc`, `24 w/(a + 16)`).
Use N, mm, MPa.

## 3. The law

Tension positive; `c` is a compression magnitude. Per call, given the fibre strains
`(eps_x, gamma_xy)`:

1. Principal strains of `(eps_x, eps_y, gamma_xy)`, coaxial and rotating.
2. `f2 = −β f_env(c)`, Popovics to `(ep, fp)`, then linear to `0.2 fp` at `e20` and flat;
   `β` on the CURRENT `e1`. `f2 = 0` if `e2 > 0`.
3. `f1 = Ec e1` to `ecr`, then `fcr/(1 + √(500 (e1 − ecr)))`; `f1 = −f_env(−e1)` if `e1 < 0`.
4. History (committed only): `em1` max principal tension, `em2` max principal compression;
   compression un/reloads linearly to the Palermo–Vecchio plastic offset, tension on the secant
   to the origin.
5. Crack check `sum`: `f1 ≤ v_ci,max tan θ + ρ max(fyt − f_st, 0)`, `w = e1 s_θ`,
   `s_θ = 1/(|sin θ|/smx + |cos θ|/smy)`, θ the strut angle from the member axis.
6. `eps_y` from `σ_y + ρ f_st(eps_y) = 0` (§4); stresses by the coaxial transformation.

The tie and bar steel: elastic, plateau, quadratic hardening to `fu` at `esu`, `1e-3 fy` past
`esu` (fracture by total strain); elastic unloading from a committed plastic strain, clipped to
the backbone of each sign at the current total strain; no Bauschinger effect, no buckling.

## 4. The root rule

Light ties (`ρ fyt < fcr`, every column of the calibration set) give the transverse residual
a **fold** at cracking: the uncracked root disappears and `eps_y` jumps to the cracked one.
Which root is taken decides the response, so the rule is Model C's, exactly:

- before the first commit, bisect `[−0.004, 0.15]`;
- after it, take the 81 candidates `eps_y,committed + {0, ±geomspace(1e−7, 0.15, 40)}` clipped to
  `[−0.01, 0.15]`, and the **increasing** sign change (r ≤ 0, then r ≥ 0) whose midpoint is
  nearest the committed root; none: bisect `[−0.004, 0.15]`;
- 60 halvings, mid-point returned. There is no "no root" failure: the bisection end point is
  the answer, as in Model C.

The C++ scan evaluates the residual lazily, interval by interval in order of distance, and
stops at the first increasing one; it is the same interval as the full scan (g++ gate K1, 4000
random states, bit for bit). A call costs about 60 residual evaluations when a root sits next
to the committed one, at most 141.

## 5. Tangent and the out-of-plane shear

The condensed tangent `∂(σx, τxy)/∂(εx, γxy)` is a forward difference, each perturbed state
re-solving `eps_y` from the committed memory: Model C's own section tangent, done per fibre.
At a fold it is large (the root jumps); that is the model, not a defect. `getInitialTangent`
is the uncracked one, `diag(Ec, Ec/2, Gxz)` (Model C's membrane has no Poisson effect).

The strip is a membrane in x–y. In the 3D (BeamFiber) view `τxz = Gxz γxz`, elastic and
uncoupled: it gives a 3D section weak-axis shear and torsion stiffness and carries no damage.
Load the member in its x–y plane (the section's y is the strip's transverse direction).

## 6. Building the Model C section

- 40 strips through the depth `h` at `y_i = −h/2 + (i + ½) h/40`, area `b h/40` minus the area of
  any bar whose nearest strip it is (ties go to the lower index).
- Two strip materials: core strips (`|y| ≤ h/2 − cover to the tie outside`: Mander `fcc, ecc`,
  core `e20`, `-ties ρ'' fyt fut`) and cover strips (`fc, 0.002`, cover `e20`, no ties).
- Bars: `LadrunoMCFTSteel` through `LadrunoUniaxialFiber`, one fibre per bar layer (3-2-3).
- `section LadrunoNDFiber … -shearProfile parabolic -shearDepth h b -shearResultant modelC
  -noCentroid`: Model C takes moments about the geometric centre; with the stock area centroid,
  the middle bar layer (subtracted from the lower of the two central strips) moves it by
  `Ab dy/(b h)`.
- Model C strengthens the far-end section by 1.001 to break the symmetric bifurcation: scale
  that section's fibre areas by 1.001.

`V = Σ τ A` makes the section tangent unsymmetric: use an unsymmetric system (`UmfPack`,
`SparseGEN`, `FullGeneral`).

## 7. Responses

`stress`, `strain`, `tangent`; `mcft` = `(eps_y, e1, e2, w, f_st, f1, f2, β, θ, bracketed)`
(`bracketed` 1 when the scan found the root near the committed one, 0 for the full range);
`eps_y`; `crackWidth`; `tieStress`.

## 8. Verification

| Gate | What | Result |
|---|---|---|
| `mcft_strip_gpp.cpp` K1–K6 | lazy = full scan; steel; path independence; equilibrium; tangent; root rule | K1 0/4000 differ; K4 2e-16; K5 uncracked 7e-15, cracked forward vs central 5e-5 |
| `test_ladrunoMCFTStrip_cpp.py` strips | core and cover strips, 5 paths each (monotonic, crushing, cyclic, random, with uncommitted probes) vs `tests/_testbed/mcft_strip_ref.py` | bit-identical |
| `test_ladrunoMCFTStrip_cpp.py` section | Model C column section, 264 steps through shear softening | 1e-15 relative (summation order) |
| `test_ladrunoMCFTStrip.py` | the section through OpenSees vs the reference; bar law; responses; views; parsers; round trip; classic Tcl; `ndTest` | see the file |

Breaking the crack-width constant `24 → 24.1` fails 9 of the 12 g++-path tests.
