# WP-210 — `LadrunoHS` / `LadrunoHSSmall`: Unified Hardening Soil & Small-Strain Overlay Material (`-lodeShape mc|mn`)

**Status: IN PROGRESS.** Recorded 2026-10-24 on the owner's decision after the deep-research survey ([[_hardening_soil_survey_2026-10-24]]). Reserves **`ND_TAG_LadrunoHS = 33029`** in `SRC/classTags.h` (single multi-dimensional `NDMaterial` serving `ThreeDimensional`, `PlaneStrain`, `AxiSymmetric`, `PlaneStress`, and `PlateFiber` via `dim_`).

## Why Standalone `LadrunoHS` (Not `ASDPlasticMaterial3D`)

1. `ASDPlasticMaterial3D<EL, YF, PF>` is a single-surface template with one yield function and one scalar `dLambda`. Because Hardening Soil requires coupled shear-cone + failure-truncation + elliptical-cap + Rankine-tension yielding, `ASD_material_definitions.cpp:433-461` had to split `StiffSoil` into two uncoupled materials (`StiffSoilShear` and `StiffSoilCap`).
2. `StiffSoilShear_YF.h` omits the ultimate Mohr–Coulomb truncation surface at $q = q_f$ (allowing $q$ to climb toward $0.99\,q_f/R_f$), and `StiffSoilCap_YF.h` conflates the Lode parameter $\delta = \frac{3+\sin\phi'}{3-\sin\phi'}$ with the cap aspect ratio $\alpha$.
3. `ASDPlasticMaterial3D` is 3D-only (`"ThreeDimensional"`), relies on 12-evaluation finite-difference gradients, and requires the user to hand-supply calibrated $(\alpha, \beta)$ instead of solving for them automatically from $(E_{\text{oed}}^{\text{ref}}, K_0^{\text{nc}})$.
4. Standalone `LadrunoHS` (backed by OpenSees-free `SRC/material/nD/LadrunoHSKernel.h`) couples all 4 yield surfaces, auto-calibrates $(\alpha, \beta)$ at construction, supports both `-lodeShape mc` (Mohr–Coulomb hexagon, Schanz et al. 1999) and `-lodeShape mn` (Panteghini & Lagioia 2014 strictly convex Matsuoka–Nakai) for rigorous lower/upper bounding comparisons, includes Benz's (2007) `HSsmall` 6D strain-history tensor $\mathbf{H}$ behind `-smallStrain $G0_ref $gamma07` (and convenience alias `nDMaterial LadrunoHSSmall`), and serves all dimensional views (`ThreeDimensional`, `PlaneStrain`, `AxiSymmetric`, `PlaneStress`) from a single class.

## Command Syntax

```tcl
nDMaterial LadrunoHS $tag \
    $E50_ref $Eoed_ref $Eur_ref $m \
    $c $phi_deg $psi_deg \
    <-nu $nu_ur>           ;# default 0.20
    <-pref $p_ref>         ;# default 100.0 (in model stress units)
    <-Rf $Rf>              ;# default 0.90
    <-K0nc $K0nc>          ;# default 1 - sin(phi)
    <-sigmat $sigma_t>     ;# default 0.0 (Rankine tension cutoff; alias -sigma_t)
    <-OCR $ocr>            ;# default 1.0
    <-pop $pop>            ;# default 0.0 (pre-overburden pressure)
    <-e0 $e0 -emax $emax>  ;# optional void ratio + dilatancy cutoff
    <-smallStrain $G0_ref $gamma_07> ;# activates Benz (2007) HSsmall overlay
    <-lodeShape mc|mn>     ;# default mc (Schanz 1999); mn = Panteghini-Lagioia (2014) convex Matsuoka-Nakai
    <-stiffStress sigma3|p>;# default sigma3 (PLAXIS HS); p = mean-stress (GHS)
    <-alphaCap $alpha -betaCap $beta> ;# optional override of auto-calibrated cap parameters
    <-rho $rho>            ;# default 0.0
    <-implex>              ;# Oliver et al. (2008) IMPL-EX
```

*(Convenience alias `nDMaterial LadrunoHSSmall $tag $E50_ref $Eoed_ref $Eur_ref $m $c $phi_deg $psi_deg $G0_ref $gamma_07 <optional flags...>` dispatches to the same `LadrunoHS` class).*

## Phased Implementation Plan

| Phase | Content | Gate |
|---|---|---|
| **P0** | Python reference oracle & $(\alpha, \beta)$ calibrator (`tests/verification/hs_oracle.py` + `tests/test_ladruno_hs_oracle.py`) covering Tier 1 (T1-A..T1-E) and Tier 2 (T2-A $\pi$-plane bounding `mc` vs. `mn`, T2-B `PlaneStrain` compression) | All Tier 1 & Tier 2 oracle gates green in pytest |
| **P1** | OpenSees-free C++17 header-only `SRC/material/nD/LadrunoHSKernel.h` (4-surface return map, auto-calibrated $(\alpha, \beta)$, dual `-lodeShape mc\|mn`, Benz 2007 `HSsmall` 6D strain-history tensor $\mathbf{H}$, consistent algorithmic tangent $\mathbb{C}^{\text{alg}}$) + standalone C++ check `tests/ladruno_hs_kernel_check.cpp` | Kernel parity vs. Python oracle $\le 10^{-10}$; central FD tangent Gate T2-C across all 7 active-set regimes $\le 10^{-6}$ |
| **P2** | OpenSees wrapper `SRC/material/nD/LadrunoHS.{h,cpp}`, `DimMode` (`3D`, `PlaneStrain`, `AxiSymmetric`, `PlaneStress`), `OPS_LadrunoHS` / `OPS_LadrunoHSSmall`, `sendSelf`/`recvSelf`, WP-99 refusal latch, `-implex`, `updateParameter` (`stage`/`firstStage`/`updateMaterialStage`) | Build all 5 targets via `Ladruno_scripts\build.bat`; material checklist & quirk pattern gates green |
| **P3** | Tier 3 laboratory element suite (Hostun sand, Brinkgreve RD sweep, Karlsruhe fine sand vs. `LadrunoSANISAND`, Bangkok clay `CIUC`) & Tier 4 BVP benchmarks (Strip/Circular footing `mc` vs. `mn` bounding, Berlin Sand excavation) | Warrant package + user guide + `banner_features.txt` + Zone-A green |
