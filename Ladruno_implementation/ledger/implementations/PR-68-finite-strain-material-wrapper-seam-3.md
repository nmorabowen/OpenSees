---
wp: PR-68
title: "Finite-strain material wrapper (seam 3)"
pr: "#68, #70, #73, #76"
status: "LogStrain shipped (#70); full-tensor channel + finite wiring (PR #76)"
class_tags: ["ND_TAG_LogStrainNDMaterial"]
section: "table"
legacy_seq: 67
---
| **Finite-strain material wrapper (seam 3)** — logarithmic (Hencky) strain-space adaptor lifting any GREEN 3D small-strain `NDMaterial` to finite strain (**spatial** multiplicative plasticity, dSNPO 2008 Box 14.3 *MATISU*), reusing the inner return map *verbatim*: `bᵉᵗʳ=F_Δ bᵉ_n F_Δᵀ → εᵉ=½ln bᵉ → unchanged return map → τ → σ=J⁻¹τ`. Material side only: concrete `LogStrainNDMaterial` + fork-local base `FiniteStrainNDMaterial` (adds `setTrialF(F)`; **interface-only header shipped in #68**); adaptor owns committed `bᵉ` per GP + the repeated-eigenvalue degeneracy branch (D3, ours). **Element side = the [[09_ladruno_brick]] `finite` geometry method** (`LadrunoBrick`, classTag 33002, [[solid_transformation_wrapper]]); the two meet at seam 3 — now WIRED end-to-end (LadrunoBrick `-geom finite` drives LogStrain via `setTrialF`; FD-tangent green). `LogStrainNDMaterial` MERGED (#70), broker + contract-lock (#73). **Seam-3 tangent channel REVISED: added `FiniteStrainNDMaterial::getSpatialTangentTensor(double c[3][3][3][3])` (+ kernel `spatial_tangent_full`) because `c=(1/2J)[D:L:B]` is non-minor-symmetric in (k,l) and the 6×6 `getTangent` is lossy** — the element consumes the full 4th-order `c`. Plan + GREEN/YELLOW/RED matrix in `09_finite_strain_material_wrapper.md`. | Material | `ND_TAG_LogStrainNDMaterial` **33010** (element tag 33002 = companion) | `SRC/material/nD/{FiniteStrainNDMaterial.h,LogStrainNDMaterial.{cpp,h},LogStrainKernel.h}` | LogStrain shipped (#70); full-tensor channel + finite wiring (PR #76) | [#68](https://github.com/nmorabowen/OpenSees/pull/68), [#70](https://github.com/nmorabowen/OpenSees/pull/70), [#73](https://github.com/nmorabowen/OpenSees/pull/73), [#76](https://github.com/nmorabowen/OpenSees/pull/76) |
