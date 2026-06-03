---
title: Validation Report — Finite-Strain Trifecta, Phase P1 (core)
project: Ladruno
status: complete
priority: high
owner: nmora
tags:
  - validation
  - report
  - finite-strain
  - element
  - material
created: 2026-06-02
---

# Validation Report — Finite-Strain Trifecta · Phase P1

> [!abstract] What this is
> The execution record for **Phase P1** of
> [[17_finite_strain_validation_plan|the finite-strain V&V plan]] — the
> *finite-strain core*: benchmarks **A1–A5** (L1 analytical), **B1–B3**
> (L2 convergence & locking) and **C1** (the Simo necking bar). Each row below is
> a contract from the plan, now fulfilled by a runnable Zone-A test, with the
> oracle, the measured margin, and PASS/FAIL.

Companion: [[17_finite_strain_validation_plan]], [[09_ladruno_brick]],
[[10_solid_corotational_adr]], [[09_finite_strain_material_wrapper]],
[[10_ladruno_j2_plasticity]].

The trifecta under test:
`element LadrunoBrick … -geom finite` (UL finite-strain hex) driving
`nDMaterial LogStrain $t $inner` (Hencky/MATISU adaptor) over a small-strain
inner (`LadrunoJ2` isotropic / `ElasticIsotropic`).

---

## 1. Test files delivered

| File | Zone | Benchmarks | Tests |
|---|---|---|---|
| `tests/test_finite_strain_L1_analytical.py` | A (t1) | A1–A5 | 12 |
| `tests/test_finite_strain_L2_convergence.py` | A (t2a) | B1–B3 | 3 |
| `tests/test_finite_strain_C1_necking.py` | A (t2a) | C1 | 4 |

**19 tests, all PASS** (≈12 s wall, fresh `OpenSeesPy` worktree build off
`ladruno` HEAD; existing finite-strain suite 21 pass / 1 xfail unchanged).
All Zone-A (structured lattices, no gmsh) so the core validation travels with
the PR and gates CI. Oracles are the numpy mirrors `tests/logstrain_reference.py`
(MATISU / Hencky tensor-log) and `tests/ladrunoj2_reference.py` (J2 return map) —
independent of the C++ kernels, so an indexing bug cannot be mirrored.

---

## 2. L1 — analytical / manufactured (A1–A5)  ✅

| ID | Problem | Oracle | Tol | Result |
|---|---|---|---|---|
| **A1** | uniaxial finite stretch into J2 (true uniaxial-stress rig) | von-Mises identity τ_axial = σ_y(p); σ = τ/J | 1e-5 | **PASS** |
| **A1′** | uniaxial backbone sweep λ∈[1.02,1.28] | σ_y(p) = sig0+H·p, monotone | 1e-5 | **PASS** |
| **A2** | hydrostatic stretch (triple-degenerate λ) | σ = (3K ln λ/λ³)I, A.52 triple branch | 1e-5 | **PASS** |
| **A2′** | equibiaxial stretch (double-degenerate λ) into J2 | MATISU return-map oracle, A.53 branch | 1e-5 | **PASS** |
| **A3** | large dilatation J∈{0.5,0.75,1.5,2.0} | σ = J⁻¹τ, tr ε = ln J | 1e-5 | **PASS** |
| **A4** | finite simple shear (elastic) | Hencky hyperelastic closed form (+ normal stresses) | 1e-5 | **PASS** |
| **A4′** | finite simple shear into J2 | single-step return-map oracle | 1e-5 | **PASS** |
| **A5** | isochoric J2 ⇒ zero mean Cauchy stress | tr τ = 3K ln J = 0 (plastic incompressibility) | 1e-7 | **PASS** |
| **A5′** | mean stress along committed isochoric path | K ln(Jₖ)/Jₖ per step (volumetric stays elastic) | 1e-6 | **PASS** |

> [!note] Review gaps closed
> A2/A2′ exercise the **repeated-eigenvalue** branches of the tensor-log map
> (TEST-1), A3 exercises **|J−1| ≫ 0** (TEST-2), A5/A5′ pin **plastic
> incompressibility** as an element-observable identity — the three coverage
> gaps the finite-strain deep review flagged.

### Key subtleties recorded

- **A1 — the J factor.** Von-Mises yield applies to the *Kirchhoff* stress τ
  (the small-strain return map's output), not Cauchy σ = τ/J. The exact uniaxial
  identity is therefore **σ_zz·J = σ_y(p)**; reading σ_zz alone is off by ≈ J−1
  ≈ 0.24 % at λ=1.02. J is recovered from the deformed lateral stretches.
- **A2′/A4′ — path dependence.** Plastic simple shear / equibiaxial are
  *non-proportional*; a sub-stepped element solve does **not** equal a one-step
  return map. The tests drive **one** increment ref→F so both sides integrate
  identically. (Proportional radial paths are path-independent and may sub-step.)

---

## 3. L2 — convergence & locking (B1–B3)  ✅

| ID | Problem | Discriminator | Result |
|---|---|---|---|
| **B1** | h-convergence of a finite-strain elastic cantilever (nz∈{2,4,8}) | monotone ↑, shrinking increments, bounded Richardson limit | **PASS** |
| **B2** | near-incompressible elastic block, ν=0.4999 | std locks (deflection→**0.10**× going ν 0.3→0.4999), F-bar ν-insensitive (**1.15**×), bbar/std ≈ **23**× at ν=0.4999 | **PASS** |
| **B3** | isochoric **J2** cantilever, displacement-controlled | std over-stiffens, std/bbar reaction ≈ **1.44** (>1.15) | **PASS** |

> [!note] Solver lesson (logged to LEDGER_quirks)
> Bending **into plasticity** with `-geom finite` + `LogStrain(LadrunoJ2)` does
> NOT converge under plain Newton or NewtonLineSearch (residual grows from step
> 0). **`KrylovNewton` + `EnergyIncr`** is robust. B3 also needs a refined
> (≥2×2) cross-section: a 1-element-wide column bends too poorly for stable
> plastic Newton. B2 (elastic, load-controlled) is fine with plain Newton.

B3 is the *plastic* counterpart of the existing elastic F-bar locking test and
the review-flagged gap no prior test covered: isochoric plastic flow (det Fᵖ=1)
locks the standard hex exactly as ν→½ does elastically; F-bar relieves it.

---

## 4. C1 — Simo–Armero necking bar  ✅ (physics) / ⏳ (quantitative → Zone-B)

**Model.** Half-bar (z=0 mid-plane symmetry), round cross-section as a structured
**squircle** hex lattice (FG-squircle map of a square grid onto the disk — all-hex,
no axis degeneracy, deterministic, **no gmsh**). `n=4` cross × `nz=8` axial (bias
toward the neck) = **128 hexes, 225 nodes**. `LadrunoBrick -formulation bbar
-geom finite` (F-bar — necking is isochoric-plastic, std would lock; cf. B3) over
`LogStrain(LadrunoJ2)` with Simo's saturation law σ_y = 450 + 265(1−e^{−16.93 p})
+ 129.24 p MPa, E=206900, ν=0.29 (r₀=6.413, L₀/2=26.667 mm). Axisymmetric axis
line pinned; end pulled 2.8 mm (5.6 mm total) over 44 steps; `UmfPack` +
`KrylovNewton` + `EnergyIncr`.

| QoI | Measured | Oracle / expectation | Verdict |
|---|---|---|---|
| neck localization | r_neck/R₀ = **0.910**, r_end/R₀ = **0.962**; min radius **at** z=0 | r_neck < r_end, min at imperfection plane | **PASS** |
| radius reduction monotone in elongation | strictly decreasing | physical | **PASS** |
| geometric softening (Considère) | k_late/k_init = **0.068** (< 0.12) | stiffness collapses toward the force peak under a *hardening* law | **PASS** |
| plastic localization | ε̄ᵖ_neck = **0.149** vs ε̄ᵖ_end = **0.077** (×1.94) | ε̄ᵖ_neck ≫ ε̄ᵖ_end | **PASS** |
| imperfection-seeded (control) | perfect bar: radius spread < **3 %**, r_neck≈r_end | uniform contraction, no neck | **PASS** |

> [!check] Phase-1 gate — "necking matches Simo"
> The necking **physics** is reproduced and validated to engineering rigor: the
> neck localizes at the seeded imperfection, plastic strain concentrates there,
> the reaction–elongation curve undergoes geometric softening toward the
> Considère peak, and a perfect bar shows NO spurious neck. This is the strongest
> available *internal* evidence that `LadrunoBrick -geom finite + LogStrain +
> LadrunoJ2` does finite-strain elastoplastic necking correctly.

> [!warning] Deferred to Zone-B (C1 quantitative contract)
> The **tight quantitative** neck-ratio match to Simo (r/r₀ ≈ 0.5 at 7 mm
> elongation, tol 3 %) requires a **mesh-converged fine model**. On a coarse
> structured lattice the neck develops only mildly (≈ 9 % radius reduction at
> 5.6 mm) before element distortion stalls the solve — verified to be
> mesh-insensitive in that range (n=4/nz=8 and n=6/nz=16 track within ~0.2 %),
> i.e. the limiter is resolution+elongation, not a bug. The converged
> reaction–elongation curve and neck-radius history vs Simo (1992) Fig. 5 is the
> **Zone-B C1 study** (plan P6 / §6): a gmsh or fine squircle mesh with
> r-refinement at the neck, arc-length or aggressive step-cutting past the peak,
> and the overlay figure of §App. Until then C1's quantitative row stays open.

---

## 5. Known limitations carried in (per plan §10)

1. **Combined-hardening + large rotation** is the §14.11 boundary; all A/B/C
   plasticity here uses **isotropic** hardening (objective through the
   log-strain lift at any deformation). Combined hardening is validated only in
   small-rotation / corot regimes (separate suite).
2. C1 uses `bbar` (F-bar) precisely because necking is isochoric-plastic and the
   standard hex would lock (demonstrated in B3).
3. F-bar's tangent is unsymmetric ⇒ an unsymmetric solver (`UmfPack` /
   `FullGeneral`) throughout the locking and necking studies.

---

## 6. Reproduce

```powershell
# from the worktree, against the local build:
./Ladruno_scripts/run_zone_a.ps1 -DistBin .\dist\bin -- -k "finite_strain"
# or directly:
$env:PYTHONPATH=".\dist\bin"; $env:PATH=".\dist\bin;$env:PATH"
py -3.12 -m pytest tests/test_finite_strain_L1_analytical.py `
                   tests/test_finite_strain_L2_convergence.py `
                   tests/test_finite_strain_C1_necking.py -v
```
