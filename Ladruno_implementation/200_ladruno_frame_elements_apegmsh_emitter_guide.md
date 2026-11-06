---
title: "Ladruno Frame & Wall Elements — apeGmsh Emitter & Post-Processor Guide"
project: Ladruno
status: ready-for-apegmsh
priority: high
created: 2026-04-20
adr: "WP 200"
tags: [apegmsh, emitter, postprocessor, frame, beam, column, wall, mvlem, imk, cohesive-hinge, ladruno-hdf5]
---

# WP-200: Ladruno Frame & Wall Elements — `apeGmsh` Emitter & Post-Processor Guide

**Target Audience:** `apeGmsh` preprocessor/emitter and `.ladruno` HDF5 post-processor developers.  
**Scope:** Unified reference for the Ladruno frame and macro-wall element ecosystem (`LadrunoDispBeamColumn`, `LadrunoForceBeamColumn`, `LadrunoIMKBeam`, `ForceBeamColumn`, `MVLEM` / `SFI_MVLEM` / `MVLEM_3D` / `SFI_MVLEM_3D`), supporting hinge/fiber materials (`LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `LadrunoUniaxialJ2`, `LadrunoRebarBuckling`, `LadrunoBondSlip`, `ASDConcrete1D`, `ASDSteel1D`), elemental loads (including partial trapezoidal and 3D `LadrunoMembrane` tributary slab loads), `setParameter`/`updateParameter`, `eleResponse` queries, `.ladruno` HDF5 recorder schemas, and OpenMP/IMPLEX integration.

---

## 1. Overview & Decision Matrix: When to Emit Which Element

### 1.1 Frame & Wall Element Selection Matrix

| Modeling Goal / Structural Component | Recommended Element Command | Dim / DOFs | Formulation & Key Features |
| :--- | :--- | :--- | :--- |
| **General nonlinear RC / steel frame (Displacement-based)** with distributed plasticity, axial-flexural ($P\text{--}M$) fiber coupling, fracture-energy ($G_f$) regularization, exact span-load quartic bubble enhancement (`-spanLoad`), and **either** $P\text{--}\delta$ axial bowing (`-nl`) **or** embedded cohesive hinges | `LadrunoDispBeamColumn` *(or `LadrunoDispBeamColumn2d` / `3d`)* | 2D (3 DOF)<br>3D (6 DOF) | Displacement-based fiber beam-column (`-lch ip|element|$val` regularization for `ASDConcrete1D`/`ASDSteel1D`, `-spanLoad` / `-noSpanLoad` quartic bubble + initial fixed-end strain, `-nl` nonlinear bowing strain, or mutually exclusive `-hinge` / `-hingeI` / `-hingeJ` / `-hingeIJ` / `-hingeY*` / `-hingeBiaxial*` with `-pm` $P\text{--}M$ interaction and `-implex`). |
| **Force-based (flexibility) fiber frame with Ladruno hinges, $P\text{--}M$ scaling, & shear** | `LadrunoForceBeamColumn` *(or `LadrunoForceBeamColumn2d` / `3d`)* | 2D (3 DOF)<br>3D (6 DOF) | Standalone Ladruno flexibility-based fiber beam-column (`ELE_TAG_LadrunoForceBeamColumn2d = 33024`, `ELE_TAG_LadrunoForceBeamColumn3d = 33025`) with exact intra-span equilibrium, `-lch ip|element|$val` / `-totalElementLch`, Tier-2 embedded cohesive hinges (`-hinge*`, `-hingeY*`, `-hingeBiaxial*`) in series with the flexibility core, `-pm $P0 $alpha`, `-shear $GAy <$GAz>`, `-implex`, and `-cMass`. |
| **Concentrated-plasticity frame** (Ibarra–Medina–Krawinkler peak-oriented / bilinear / pinched, or cohesive end hinges) with elastic interior + optional Timoshenko shear / inelastic series shear hinges + rigid end zones | `LadrunoIMKBeam` *(or `LadrunoIMKBeam2d`)* | 2D (3 DOF)<br>3D (6 DOF) | Elastic Euler-Bernoulli / Timoshenko interior in series with 2 (2D) or 4 (3D) rotational hinge slots (`Zi, Zj, Yi, Yj` or coupled 3D `-hingeBiaxial*`), optional series inelastic shear hinges (`-shearHinge`, `-shearHingeY`, `-shearHingeZ`), optional rigid end offsets (`-endZones`), exact $2\times 2$ / $4\times 4$ static condensation, exact fixed-end forces (`q0`/`p0`), optional `-pm` $P\text{--}M$ interaction, and `-implex`. |
| **Upstream force-based (flexibility) fiber frame** with fracture-energy ($G_f$) auto-regularization | `forceBeamColumn` *(or `ForceBeamColumn2d` / `3d`)* | 2D (3 DOF)<br>3D (6 DOF) | Upstream flexibility-based fiber beam-column (`-iter maxIters tol`, `-mass`, `-cMass`, `-lch ip|element|$val` or `-totalElementLch` characteristic length for `ASDConcrete1D`/`ASDSteel1D`). |
| **2D planar RC shear wall** with uniaxial vertical macro-fibers + decoupled horizontal shear spring | `MVLEM` | 2D (3 DOF, 2 nodes) | $m$ vertical macro-fibers (concrete + steel uniaxial materials) + 1 shear spring at height $c \cdot h$. |
| **2D planar RC shear wall with shear–flexure interaction (SFI)** | `SFI_MVLEM` | 2D (3 DOF at 2 ext. nodes + $m$ internal 1-DOF nodes) | $m$ 2D RC panel `NDMaterial` (`FSAM`) macro-fibers capturing biaxial $\sigma_x\text{--}\sigma_y\text{--}\tau_{xy}$ interaction. |
| **3D RC shear wall / core wall** with uniaxial vertical fibers + elastic out-of-plane plate bending | `MVLEM_3D` | 3D (6 DOF, 4 quad nodes) | 4-node 24-DOF quad wall element; in-plane `MVLEM` kinematics + out-of-plane Kirchhoff plate bending (`-ThickMod`, `-Poisson`). |
| **3D RC shear wall / core wall with SFI** | `SFI_MVLEM_3D` | 3D (6 DOF at 4 ext. nodes + $m$ internal 1-DOF nodes) | 4-node external quad + $m$ internal panel DOFs; in-plane `SFI_MVLEM` (`FSAM`) + out-of-plane Kirchhoff plate bending (`-ThickMod`, `-Poisson`). |

---

## 2. Element Command Signatures (Tcl & Python `opensees.pyd`)

> **Parser Contract (WP-200):** All Ladruno beam parsers (`OPS_LadrunoDispBeamColumn2d`, `OPS_LadrunoDispBeamColumn3d`, `OPS_LadrunoForceBeamColumn`, `OPS_LadrunoIMKBeam2d`, `OPS_LadrunoIMKBeam`) and `OPS_ForceBeamColumn2d/3d` use fail-closed option parsing. Unknown `-flag` tokens or truncated/invalid option arguments raise an immediate error and return `nullptr`.

### 2.1 `LadrunoDispBeamColumn` (`LadrunoDispBeamColumn2d` / `LadrunoDispBeamColumn3d`)

* **Class Tags (`SRC/classTags.h`):**
  * `ELE_TAG_LadrunoDispBeamColumn2d = 33013`
  * `ELE_TAG_LadrunoDispBeamColumn3d = 33014`
* **Dispatch Rule:** `element LadrunoDispBeamColumn` automatically dispatches to `LadrunoDispBeamColumn2d` when `(ndm == 2 && ndf == 3)` and to `LadrunoDispBeamColumn3d` when `(ndm == 3 && ndf == 6)`. Explicit `LadrunoDispBeamColumn2d` and `LadrunoDispBeamColumn3d` commands are also registered.

#### 2D Signature (`ndm = 2, ndf = 3`)

```tcl
element LadrunoDispBeamColumn $eleTag $iNode $jNode $transfTag $beamIntTag \
    <-mass $rho> <-cMass> <-damp $dampingTag> <-lch ip|element|$val> <-nl> \
    <-spanLoad | -noSpanLoad> \
    <-hinge $matTag> \
    <-hingeI $matTagI> <-hingeJ $matTagJ> <-hingeIJ $matTagI <$matTagJ>> \
    <-hingeLoc $xiI <$xiJ>> <-implex> <-pm $P0 <$alpha>>
```

**Python (`opensees.pyd`) — Example with Dual End Hinges:**
```python
ops.element(
    "LadrunoDispBeamColumn", eleTag, iNode, jNode, transfTag, beamIntTag,
    "-mass", rho, "-cMass", "-lch", "ip",
    "-hingeI", matTagI, "-hingeJ", matTagJ, "-hingeLoc", xiI, xiJ,
    "-implex", "-pm", P0, alpha,
)
```

**Python (`opensees.pyd`) — Example with Geometric Axial Bowing (`-nl`, no hinges):**
```python
ops.element(
    "LadrunoDispBeamColumn", eleTag, iNode, jNode, transfTag, beamIntTag,
    "-mass", rho, "-cMass", "-lch", "ip", "-nl",
)
```

#### 3D Signature (`ndm = 3, ndf = 6`)
Accepts all 2D flags above plus minor-axis (`-hingeY`, `-hingeYI`, `-hingeYJ`, `-hingeYIJ`) and coupled 2-DOF biaxial (`-hingeBiaxial`, `-hingeBiaxialI`, `-hingeBiaxialJ`, `-hingeBiaxialIJ` taking an `NDMaterial` tag of class `LadrunoCohesiveHingeBiaxial`):

```tcl
element LadrunoDispBeamColumn $eleTag $iNode $jNode $transfTag $beamIntTag \
    <-mass $rho> <-cMass> <-damp $dampingTag> <-lch ip|element|$val> <-nl> \
    <-spanLoad | -noSpanLoad> \
    <-hinge $matTagZ> <-hingeI $matTagZI> <-hingeJ $matTagZJ> <-hingeIJ $matTagZI <$matTagZJ>> \
    <-hingeY $matTagY> <-hingeYI $matTagYI> <-hingeYJ $matTagYJ> <-hingeYIJ $matTagYI <$matTagYJ>> \
    <-hingeBiaxial $ndMatTag> <-hingeBiaxialI $ndMatTagI> <-hingeBiaxialJ $ndMatTagJ> \
    <-hingeBiaxialIJ $ndMatTagI <$ndMatTagJ>> \
    <-hingeLoc $xiI <$xiJ>> <-implex> <-pm $P0 <$alpha>>
```

> **Mutual-Exclusion Rules (`LadrunoDispBeamColumn2d` / `3d`):**
> 1. **`-nl` vs. `-hinge*`:** `-nl` (nonlinear $\frac{1}{2}\theta^2$ axial bowing strain) and any embedded hinge flag (`-hinge*`, `-hingeY*`, `-hingeBiaxial*`) are **mutually exclusive** (`WARNING LadrunoDispBeamColumn2d/3d: -hinge and -nl are mutually exclusive`). Do not emit both on the same element.
> 2. **`-hingeBiaxial*` vs. `-hinge*` / `-hingeY*` (3D):** Coupled 2-DOF biaxial hinges (`-hingeBiaxial*` / `-hingeC*`) are **mutually exclusive** with scalar uniaxial hinges (`-hinge*` / `-hingeY*`).
> 3. **`-hingeY*` requires `-hinge*` (3D):** Weak-axis uniaxial hinges (`-hingeY`, `-hingeYI`, `-hingeYJ`, `-hingeYIJ`) are an add-on to the strong-axis hinge at the same end (`-hingeY` requires `-hinge`/`-hingeI`, and `-hingeYJ` requires `-hingeJ`/`-hingeIJ`).

#### Option Flag Reference (`LadrunoDispBeamColumn2d` / `3d`)

| Flag | Arguments | Default | Description |
| :--- | :--- | :--- | :--- |
| `-mass` | `$rho` (float) | `0.0` | Mass per unit length. |
| `-cMass` *(or `-cmass`)* | *(none)* | `0` (lumped) | Consistent mass matrix flag. |
| `-damp` | `$dampingTag` (int) | `0` (none) | Elemental `Damping` object tag. |
| `-lch` | `'ip'`, `'element'`, or `$val` ($>0$) | `'ip'` (`lchMode = 0`) | Characteristic length $L_{ch}$ for `ASDConcrete1D` / `ASDSteel1D` fracture-energy regularization (`ip`: $L_{ch,i} = w_i \cdot L$; `element`: $L_{ch,i} = L$; `$val`: user-specified positive length — accepts either numeric float or string in Python/Tcl). |
| `-spanLoad` / `-noSpanLoad` | *(none)* | `-spanLoad` (`true`) | Enables (`-spanLoad`, default) or disables (`-noSpanLoad`) exact span-load internal strain $\mathbf{e}_{0,i} = \mathbf{k}_{s0,i}^{-1}\mathbf{s}^{\text{fixed}}(\xi_i)$ and discrete-orthogonalized coupled axial + quartic flexural bubble enhancement $\phi_b(\xi) = 1 - 6\xi + 6\xi^2 - (c_0 + c_1\xi)$ on $(\beta_P, \beta_Z)$ in 2D and $(\beta_P, \beta_Z, \beta_Y)$ in 3D, eliminating both 1-element midspan moment blindness and neutral-axis axial locking on cracked/yielding RC/T-beam fiber sections under distributed/point span loads. |
| `-nl` | *(none)* | `false` | Enables nonlinear geometric axial bowing strain ($\frac{1}{2}\theta_z^2$ in 2D; $\frac{1}{2}(\theta_z^2 + \theta_y^2)$ in 3D). Mutually exclusive with all `-hinge*` flags. |
| `-hinge` | `$matTag` (int) | disabled | Single-hinge mode about local $z$: clones `UniaxialMaterial` `$matTag` for strong-axis bending. |
| `-hingeI` *(aliases: `-hinge1`, and `-hingeZI` in 3D)* | `$matTagI` (int) | disabled | Dual-end hinge mode about local $z$ at end $I$ ($2\times 2$ Hermite static condensation). |
| `-hingeJ` *(aliases: `-hinge2`, and `-hingeZJ` in 3D)* | `$matTagJ` (int) | disabled | Dual-end hinge mode about local $z$ at end $J$. |
| `-hingeIJ` *(alias `-hingeZIJ`)* | `$matTagI` `<$matTagJ>` | disabled | Dual-end hinges about local $z$ at both ends $I$ and $J$ (if `$matTagJ` is omitted, uses `$matTagI` at both ends). |
| `-hingeLoc` | `$xiI` `<$xiJ>` (floats $\in [0, 1]$) | `0.0`, `1.0` | Normalized beam coordinates $\xi_I, \xi_J \in [0, 1]$ for dual-end hinges (`Hermite` moment interpolation $b_1(\xi) = 6\xi - 4$, $b_2(\xi) = 6\xi - 2$). If `$xiJ` is omitted, defaults to $\xi_J = 1 - \xi_I$. |
| `-hingeY` | `$matTagY` (int, **3D only**) | disabled | Single-hinge mode about local $y$ (minor axis). Requires `-hinge`. |
| `-hingeYI` *(alias `-hingeY1`)* / `-hingeYJ` *(alias `-hingeY2`)* | `$matTagY` (int, **3D only**) | disabled | Dual-end hinge mode about local $y$ at end $I$ and/or end $J$. Requires corresponding $z$-axis hinge (`-hingeI`/`-hingeJ`). |
| `-hingeYIJ` | `$matTagYI` `<$matTagYJ>` (**3D only**) | disabled | Dual-end hinges about local $y$ at both ends $I$ and $J$ (defaults `$matTagYJ = $matTagYI` if omitted). Requires `-hingeIJ` (or `-hingeI` + `-hingeJ`). |
| `-hingeBiaxial` *(alias `-hingeC`)* | `$ndMatTag` (int, **3D only**) | disabled | Coupled 2-DOF biaxial hinge (`LadrunoCohesiveHingeBiaxial`) in single-hinge mode. Mutually exclusive with `-hinge*`/`-hingeY*`. |
| `-hingeBiaxialI` *(aliases `-hingeCI`, `-hingeBiaxial1`)* / `-hingeBiaxialJ` *(aliases `-hingeCJ`, `-hingeBiaxial2`)* | `$ndMatTag` (int, **3D only**) | disabled | Coupled 2-DOF biaxial hinge (`LadrunoCohesiveHingeBiaxial`) at end $I$ or end $J$ in dual-hinge mode. |
| `-hingeBiaxialIJ` *(alias `-hingeCIJ`)* | `$ndMatTagI` `<$ndMatTagJ>` (**3D only**) | disabled | Coupled 2-DOF biaxial hinges at both ends $I$ and $J$ (defaults `$ndMatTagJ = $ndMatTagI` if omitted). |
| `-implex` *(or `-IMPLEX`)* | *(none)* | `false` | Enables element-level IMPLEX explicit-implicit extrapolation on embedded hinges (`LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `Bilin`, `ModIMKPeakOriented`). |
| `-pm` *(or `-PM`)* | `$P0` `<$alpha>` | `P0 = 0.0` (off), `alpha = 0.5` | Axial-moment ($P\text{--}M$) interaction on embedded hinges: scales hinge capacity by the parabolic factor $f_{pm}(N) = \operatorname{clamp}\!\left(1 - \alpha \frac{N}{P_0} - (1 - \alpha)\left(\frac{N}{P_0}\right)^2,\; 0.05,\; 2.5\right)$, where basic axial force $N < 0$ in compression ($N > 0$ in tension) and $P_0 > 0$. |

---

### 2.2 `LadrunoIMKBeam` (`LadrunoIMKBeam2d` / `LadrunoIMKBeam`)

* **Class Tags (`SRC/classTags.h`):**
  * `ELE_TAG_LadrunoIMKBeam = 33003` (3D, 2 nodes, 6 DOFs/node)
  * `ELE_TAG_LadrunoIMKBeam2d = 33004` (2D, 2 nodes, 3 DOFs/node)
* **Dispatch Rule:** `element LadrunoIMKBeam` dispatches to `LadrunoIMKBeam2d` when `ndm == 2` and to 3D `LadrunoIMKBeam` when `ndm == 3`. Explicit `element LadrunoIMKBeam2d` is also registered for 2D.

#### 2D Signature (`ndm = 2, ndf = 3`)
```tcl
element LadrunoIMKBeam $eleTag $iNode $jNode $A $E $Iz $transfTag \
    <-matZ $tag> <-matZi $tag> <-matZj $tag> \
    <-hinge $loc> <-shear $GAvy | -shear $G $Avy> \
    <-shearHinge $matTagVy> \
    <-endZones $dI $dJ> <-mass $m> <-cMass> \
    <-pm $P0 <$alpha>> <-implex>
```

**Python (`opensees.pyd`):**
```python
ops.element(
    "LadrunoIMKBeam", eleTag, iNode, jNode, A, E, Iz, transfTag,
    "-hinge", "both", "-matZi", matZiTag, "-matZj", matZjTag,
    "-shear", GAvy, "-shearHinge", matVyTag, "-endZones", dI, dJ,
    "-mass", m, "-cMass", "-pm", P0, alpha, "-implex",
)
```

#### 3D Signature (`ndm = 3, ndf = 6`)
```tcl
element LadrunoIMKBeam $eleTag $iNode $jNode $A $E $G $Jx $Iy $Iz $transfTag \
    <-matZ $tag> <-matZi $tag> <-matZj $tag> \
    <-matY $tag> <-matYi $tag> <-matYj $tag> \
    <-hingeBiaxial $ndTag> <-hingeBiaxialI $ndTagI> <-hingeBiaxialJ $ndTagJ> \
    <-hingeBiaxialIJ $ndTagI <$ndTagJ>> \
    <-hinge $loc> <-shear $GAvy $GAvz | -shear $G $Avy $Avz> \
    <-shearHingeY $matTagVy> <-shearHingeZ $matTagVz> <-shearHinge $matTagVy <$matTagVz>> \
    <-endZones $dI $dJ> <-mass $m> <-cMass> \
    <-pm $P0 <$alpha>> <-implex>
```

**Python (`opensees.pyd`):**
```python
ops.element(
    "LadrunoIMKBeam", eleTag, iNode, jNode, A, E, G, Jx, Iy, Iz, transfTag,
    "-hinge", "both", "-matZi", matZiTag, "-matZj", matZjTag,
    "-matYi", matYiTag, "-matYj", matYjTag,
    "-shear", GAvy, GAvz, "-shearHingeY", matVyTag, "-shearHingeZ", matVzTag,
    "-endZones", dI, dJ,
    "-mass", m, "-cMass", "-pm", P0, alpha, "-implex",
)
```

#### Option Flag Reference (`LadrunoIMKBeam2d` / `LadrunoIMKBeam`)

| Flag | Arguments | Default | Description |
| :--- | :--- | :--- | :--- |
| `-matZ` | `$matTag` (int) | `0` (elastic) | Uniaxial rotational hinge material tag for both ends (`Zi` and `Zj`) about local $z$. |
| `-matZi` / `-matZj` | `$matTag` (int) | overrides `-matZ` | Individual end-$I$ or end-$J$ hinge material tag about local $z$. |
| `-matY` | `$matTag` (int, **3D only**) | `0` (elastic) | Uniaxial rotational hinge material tag for both ends (`Yi` and `Yj`) about local $y$. |
| `-matYi` / `-matYj` | `$matTag` (int, **3D only**) | overrides `-matY` | Individual end-$I$ or end-$J$ hinge material tag about local $y$. |
| `-hingeBiaxial` *(alias `-hingeC`)* | `$ndMatTag` (int, **3D only**) | disabled | Coupled 2-DOF biaxial `NDMaterial` hinge (`LadrunoCohesiveHingeBiaxial`) at both ends (or masked by `-hinge i\|j`). Mutually exclusive with `-matZ*`/`-matY*`. |
| `-hingeBiaxialI` / `-hingeBiaxialJ` / `-hingeBiaxialIJ` | `$ndMatTagI <$ndMatTagJ>` (**3D only**) | disabled | Coupled 2-DOF biaxial `NDMaterial` hinge at end $I$, end $J$, or both ends. |
| `-hinge` | `'both'`, `'i'`, `'j'`, `'none'` (case-insensitive) | `'both'` | Active hinge location mask for symmetric `-matZ`/`-matY`/`-hingeBiaxial` (`'i'` disables end-$J$ hinges; `'j'` disables end-$I$ hinges; `'none'` makes symmetric slots purely elastic). |
| `-shear` | 2D: `$GAvy` *(1 arg)* or `$G $Avy` *(2 args)*<br>3D: `$GAvy $GAvz` *(2 args)* or `$G $Avy $Avz` *(3 args)* | `0.0` (Euler-Bernoulli) | Effective shear rigidity $GA_{vy}$ (and $GA_{vz}$ in 3D). When $>0$, activates Timoshenko shear flexibility $\Phi_z = \frac{12 E I_z}{G A_{vy} L_c^2}$ (and $\Phi_y = \frac{12 E I_y}{G A_{vz} L_c^2}$). |
| `-shearHinge` / `-shearHingeY` / `-shearHingeZ` | 2D: `-shearHinge $matTagVy`<br>3D: `-shearHingeY $matTagVy`, `-shearHingeZ $matTagVz`, or `-shearHinge $matTagVy <$matTagVz>` | disabled | Inelastic series shear hinge (`UniaxialMaterial` relating shear displacement jump $\Delta_v$ to clear-span shear force $V = -(M_i + M_j)/L_c$), adding compliance $L_c^2 / K_{s,\text{hinge}}$ in series with the flexural + Timoshenko element core. |
| `-endZones` *(alias `-endZone`)* | `$dI $dJ` (floats $\ge 0$) | `0.0 0.0` | Rigid end-zone offsets $d_I, d_J$ at ends $I$ and $J$ ($d_I + d_J < L$); clear flexible span is $L_c = L - d_I - d_J$ with kinematic transformation $\mathbf{T}$ mapping chord rotations to clear-span hinge rotations. |
| `-mass` | `$m` (float) | `0.0` | Mass per unit length. |
| `-cMass` *(or `-cmass`)* | *(none)* | `0` (lumped) | Consistent 2-node Hermite beam mass matrix flag. |
| `-pm` | `$P0` `<$alpha>` | `P0 = 0.0` (off), `alpha = 0.5` | Axial-moment ($P\text{--}M$) interaction: scales hinge capacity by $f_{pm}(N) = \operatorname{clamp}\!\left(1 - \alpha \frac{N}{P_0} - (1 - \alpha)\left(\frac{N}{P_0}\right)^2,\; 0.05,\; 2.5\right)$ ($N < 0$ in compression, default $\alpha = 0.5$). |
| `-implex` | *(none)* | `false` | Enables element-level IMPLEX explicit plastic hinge rotation extrapolation during global Newton step. |

---

### 2.3 `LadrunoForceBeamColumn` (`LadrunoForceBeamColumn2d` / `LadrunoForceBeamColumn3d`) & `forceBeamColumn`

* **Class Tags (`SRC/classTags.h`):**
  * `ELE_TAG_LadrunoForceBeamColumn2d = 33024`
  * `ELE_TAG_LadrunoForceBeamColumn3d = 33025`
  * `ELE_TAG_ForceBeamColumn2d = 73`
  * `ELE_TAG_ForceBeamColumn3d = 74`
* **Command Signatures (`LadrunoForceBeamColumn`):**

```tcl
# 2D (ndm == 2, ndf == 3):
element LadrunoForceBeamColumn $eleTag $iNode $jNode $transfTag $beamIntTag \
    <-mass $rho> <-cMass> <-damp $dampingTag> <-iter $maxIters $tol> \
    <-subdivide $numSub <$subFac>> <-lch ip|element|auto|$val> <-totalElementLch> \
    <-hinge $matTag> <-hingeI $matTagI> <-hingeJ $matTagJ> <-hingeIJ $matTagI <$matTagJ>> \
    <-hingeLoc $xiI <$xiJ>> <-implex> <-pm $P0 <$alpha>> <-shear $GAy>

# 3D (ndm == 3, ndf == 6):
element LadrunoForceBeamColumn $eleTag $iNode $jNode $transfTag $beamIntTag \
    <-mass $rho> <-cMass> <-damp $dampingTag> <-iter $maxIters $tol> \
    <-subdivide $numSub <$subFac>> <-lch ip|element|auto|$val> <-totalElementLch> \
    <-hinge $matTagZ> <-hingeI $matTagZI> <-hingeJ $matTagZJ> <-hingeIJ $matTagZI <$matTagZJ>> \
    <-hingeY $matTagY> <-hingeYI $matTagYI> <-hingeYJ $matTagYJ> <-hingeYIJ $matTagYI <$matTagYJ>> \
    <-hingeBiaxial $ndMatTag> <-hingeBiaxialI $ndMatTagI> <-hingeBiaxialJ $ndMatTagJ> \
    <-hingeBiaxialIJ $ndMatTagI <$ndMatTagJ>> \
    <-hingeLoc $xiI <$xiJ>> <-implex> <-pm $P0 <$alpha>> <-shear $GAy <$GAz>>
```
*(Explicit `LadrunoForceBeamColumn2d` and `LadrunoForceBeamColumn3d` commands, as well as the legacy inline `$numSecs $secTag1 ... $secTagN $transfTag` form, are also supported.)*

---

### 2.4 Macro-Fiber Wall Elements (`MVLEM`, `SFI_MVLEM`, `MVLEM_3D`, `SFI_MVLEM_3D`)

* **Class Tags (`SRC/classTags.h`):**
  * `ELE_TAG_MVLEM = 162` (2D, 2 external nodes, 3 DOFs/node)
  * `ELE_TAG_SFI_MVLEM = 163` (2D, 2 external 3-DOF nodes + $m$ internal 1-DOF nodes created automatically)
  * `ELE_TAG_MVLEM_3D = 212` (3D, 4 external quad nodes `iNode jNode kNode lNode`, 6 DOFs/node)
  * `ELE_TAG_SFI_MVLEM_3D = 213` (3D, 4 external 6-DOF quad nodes + $m$ internal 1-DOF nodes)
* **WP-200 Enhancements:**
  * All four elements set `ops_TheActiveElement = this` in `setDomain()` and `update()` and implement `getCharacteristicLength() = h` (wall element height from bottom to top node), enabling automatic fracture-energy regularization ($L_{ch} = h$) for `ASDConcrete1D` and `ASDSteel1D` fibers inside wall macro-fibers.

#### Signatures

**`MVLEM` (2D):**
```tcl
element MVLEM $eleTag $Dens $iNode $jNode $m $c \
    -thick $t1 ... $tm \
    -width $b1 ... $bm \
    -rho $rho1 ... $rhom \
    -matConcrete $matC1 ... $matCm \
    -matSteel $matS1 ... $matSm \
    -matShear $matShearTag
```

**`SFI_MVLEM` (2D):**
```tcl
element SFI_MVLEM $eleTag $iNode $jNode $m $c \
    -thick $t1 ... $tm \
    -width $b1 ... $bm \
    -mat $fsamTag1 ... $fsamTagm \
    <-Coupling $couplingFlag>
```
*(`-Coupling` or `-coupling` is optional: `0` = uncoupled diagonal panel tangent `D(0,0), D(1,1), D(2,2)` (default), `1` = axial $\sigma_x\text{--}\sigma_y$ coupling `D(0,1), D(1,0)`, `2` = full $3\times 3$ axial+shear tangent coupling.)*

**`MVLEM_3D` (3D 4-Node Wall):**
```tcl
element MVLEM_3D $eleTag $iNode $jNode $kNode $lNode $m \
    -thick $t1 ... $tm \
    -width $b1 ... $bm \
    -rho $rho1 ... $rhom \
    -matConcrete $matC1 ... $matCm \
    -matSteel $matS1 ... $matSm \
    -matShear $matShearTag \
    <-CoR $c> <-ThickMod $tMod> <-Poisson $nu> <-Density $Dens>
```
*(Defaults: `c = 0.4`, `tMod = 0.63`, `nu = 0.25`, `Dens = 0.0`. Note exact capitalization `-ThickMod` — lowercase `-thickMod` is rejected. Node ordering: `iNode, jNode` at bottom edge; `kNode, lNode` at top edge counter-clockwise.)*

**`SFI_MVLEM_3D` (3D 4-Node SFI Wall):**
```tcl
element SFI_MVLEM_3D $eleTag $iNode $jNode $kNode $lNode $m \
    -thick $t1 ... $tm \
    -width $b1 ... $bm \
    -mat $fsamTag1 ... $fsamTagm \
    <-CoR $c> <-ThickMod $tMod> <-Poisson $nu> <-Density $Dens> <-Eave $Eave>
```
*(Defaults: `c = 0.4`, `tMod = 0.63`, `nu = 0.25`, `Dens = 0.0`, `Eave = 0.0` [computed from panels]. Note exact capitalization `-ThickMod`.)*

---

## 3. Supporting Hinge & Fiber Material Signatures

### 3.1 `uniaxialMaterial LadrunoCohesiveHinge`
* **Class Tag:** `MAT_TAG_LadrunoCohesiveHinge = 33003`
* **Signature:**
```tcl
uniaxialMaterial LadrunoCohesiveHinge $matTag $Mc $Gf \
    <-penalty $K | -Kpen $K | -penaltyRatio $r> \
    <-exp | -exponential | -linear | -lin> \
    <-pm $P0 <$alpha>> <-implex>
```
* **Defaults & Options:**
  * `$Mc` ($>0$), `$Gf` ($>0$): Cracking/yield moment capacity $M_c$ and rotational fracture energy $\mathcal{G}_f$ (exact area under the monotonic $M\text{--}[\![\theta]\!]$ envelope).
  * `-penalty $K` *(alias `-Kpen $K`)* or `-penaltyRatio $r`: Initial pre-peak penalty rotational stiffness $K_{\text{pen}}$. By default (`-penaltyRatio 1000.0`), $K_{\text{pen}} = r \cdot K_{\text{pen,floor}}$ where $K_{\text{pen,floor}} = M_c^2 / (2\mathcal{G}_f)$ guarantees positive softening energy $E_{\text{soft}} = \mathcal{G}_f - M_c^2 / (2 K_{\text{pen}}) > 0$.
  * `-exp` *(alias `-exponential`, default)* or `-linear` *(alias `-lin`)*: Selects exponential softening $M(\kappa) = M_c \exp(-a_{\text{exp}}(\kappa - \kappa_0))$ or linear softening $M(\kappa) = M_c (1 - (\kappa - \kappa_0)/\kappa_f)$ with origin-directed secant unloading/reloading.
  * `-pm $P0 <$alpha>`: Parabolic axial-moment interaction scaling $f_{pm}(N) = \operatorname{clamp}\!\left(1 - \alpha \frac{N}{P_0} - (1 - \alpha)\left(\frac{N}{P_0}\right)^2,\; 0.05,\; 2.5\right)$ where $N < 0$ in compression (`P0` default `0.0` = off, `alpha` default `0.5`).
  * `-implex`: Enables explicit-implicit IMPLEX stress extrapolation ($\tilde{M}_{n+1}$) and initial-tangent substitution during global Newton iterations.

### 3.2 `nDMaterial LadrunoCohesiveHingeBiaxial`
* **Class Tag:** `ND_TAG_LadrunoCohesiveHingeBiaxial = 33004`
* **Signature:**
```tcl
nDMaterial LadrunoCohesiveHingeBiaxial $matTag $Mcz $Gfz $Mcy $Gfy \
    <-exp | -exponential | -linear | -lin> \
    <-penaltyRatio $r> <-bk $eta | -eta $eta | -modemix $eta> \
    <-pm $P0 <$alpha>> <-implex>
```
* **Defaults & Options:**
  * `$Mcz, $Gfz, $Mcy, $Gfy` ($>0$): Major-axis ($z$) and minor-axis ($y$) yield moments and rotational fracture energies.
  * `-exp` *(alias `-exponential`, default)* or `-linear` *(alias `-lin`)*: Exponential or linear softening on the elliptical effective jump norm $r = \sqrt{(\alpha_z / a_{0,z})^2 + (\alpha_y / a_{0,y})^2}$.
  * `-penaltyRatio $r`: Ratio above the per-axis energy floor $K_{\text{pen},z} = r \cdot M_{cz}^2 / (2\mathcal{G}_{fz})$, $K_{\text{pen},y} = r \cdot M_{cy}^2 / (2\mathcal{G}_{fy})$ (default `1000.0`).
  * `-bk $eta` *(aliases `-eta $eta`, `-modemix $eta`)*: Benzeggagh–Kenane mode-mix exponent $\eta > 0$ for mixed-mode fracture energy $\mathcal{G}_{f,\text{mix}} = \mathcal{G}_{fz} + (\mathcal{G}_{fy} - \mathcal{G}_{fz}) w_y^\eta$ (default `1.0` = linear interpolation).
  * `-pm $P0 <$alpha>`: Biaxial $P\text{--}M_z\text{--}M_y$ capacity scaling $f_{pm}(N) = \operatorname{clamp}\!\left(1 - \alpha \frac{N}{P_0} - (1 - \alpha)\left(\frac{N}{P_0}\right)^2,\; 0.05,\; 2.5\right)$ (`P0` default `0.0`, `alpha` default `0.5`).
  * `-implex`: Enables IMPLEX biaxial extrapolation.

### 3.3 `uniaxialMaterial ModIMKPeakOriented` & `Bilin`
* **Class Tags:** `MAT_TAG_ModIMKPeakOriented = 67`, `MAT_TAG_Bilin = 2007`
* **Signatures (with WP-200 `-implex` support):**
```tcl
uniaxialMaterial Bilin $matTag $Ke $AsPos $AsNeg $My_pos $My_neg \
    $LamdaS $LamdaD $LamdaA $LamdaK $Cs $Cd $Ca $Ck \
    $Thetap_pos $Thetap_neg $Thetapc_pos $Thetapc_neg \
    $KPos $KNeg $Thetau_pos $Thetau_neg $PDPlus $PDNeg <$nFactor> <-implex>

uniaxialMaterial ModIMKPeakOriented $matTag $Ke $alfaPos $alfaNeg $My_pos $My_neg \
    $Ls $Ld $La $Lk $Cs $Cd $Ca $Ck \
    $thetaPpos $thetaPneg $thetaPCpos $thetaPCneg \
    $ResfacPos $ResfacNeg $fracDispPos $fracDispNeg $DPos $DNeg <$nFactor> <-implex>
```
* **Usage with `LadrunoIMKBeam` & `LadrunoDispBeamColumn`:** Any OpenSees `UniaxialMaterial` representing a moment–rotation ($M\text{--}\theta$) spring can be assigned to `-matZ`, `-matZi`, `-matZj`, `-matY`, `-matYi`, `-matYj` in `LadrunoIMKBeam2d` / `LadrunoIMKBeam` (or `-hinge*` in `LadrunoDispBeamColumn`). Both `Bilin` and `ModIMKPeakOriented` accept a trailing `-implex` flag and respond to `'implexError'` (and `'energy'` / `'dissipatedEnergy'`).

### 3.4 Fiber / Reinforcing Bar / Bond-Slip Uniaxial Materials

| Material Command | Class Tag | Tcl / Python Signature |
| :--- | :--- | :--- |
| `LadrunoUniaxialJ2` | `33000` | `uniaxialMaterial LadrunoUniaxialJ2 $tag $E <-iso voce $s0 $Qinf $b $Hiso> <-kin $N $C1 $g1 ...> <-damage lemaitre $r $s $pD $Dc \| -damage none> <-implex>` |
| `LadrunoRebarBuckling` | `33001` | `uniaxialMaterial LadrunoRebarBuckling $tag $matTag <-lsr $lsr> <-model dm\|ga> <-alpha $a> <-reduction $r> <-fsufrac $g> <-fy $fy> <-E $E> <-restraighten lambda \| -restraighten c $val>` |
| `LadrunoBondSlip` | `33002` | `uniaxialMaterial LadrunoBondSlip $tag $tau_max $s1 $s2 $s3 $tau_f $alpha <-Gf $Gf> <-s0 $s0>` |
| `ASDConcrete1D` | `231` | `uniaxialMaterial ASDConcrete1D $tag $E <-fc $fc> <-ft $ft> <-Te $Te -Ts $Ts <-Td $Td>> <-Ce $Ce -Cs $Cs <-Cd $Cd>> <-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha> <-eta $eta> <-tangent> <-autoRegularization $lch_ref>` |
| `ASDSteel1D` | `235` | `uniaxialMaterial ASDSteel1D $tag $E $sy $su $eu <-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-auto_regularization> <-buckling $lch <$r>> <-fracture <$r>> <-slip $matTag <$r>> <-K_alpha $K_alpha> <-max_iter $max_iter> <-tolU $tolU> <-tolR $tolR>` |

---

## 4. Elemental Loads & Tributary Slab Transfer (`eleLoad`)

All frame families (`LadrunoDispBeamColumn2d/3d`, `LadrunoForceBeamColumn2d/3d`, `LadrunoIMKBeam2d/3d`, and `ForceBeamColumn2d/3d`) support full elemental span loads with exact fixed-end forces (`q0`) and basic reactions (`p0`), plus 3D `LadrunoMembrane` tributary slab load transfer via `"localAxes"`:

| Load Type | OpenSees `eleLoad` Syntax | Supported Elements | Notes |
| :--- | :--- | :--- | :--- |
| **2D Uniform Load** (`LOAD_TAG_Beam2dUniformLoad`) | `eleLoad -ele $eleTag -type -beamUniform $Wy <$Wx>` | `LadrunoDispBeamColumn2d`, `LadrunoForceBeamColumn2d`, `LadrunoIMKBeam2d`, `ForceBeamColumn2d` | `$Wy` transverse (+ve along local $+y$), `$Wx` axial (+ve $I \to J$). |
| **2D Partial / Trapezoidal Load** (`LOAD_TAG_Beam2dPartialUniformLoad`) | `eleLoad -ele $eleTag -type -beamUniform $Wya $Wxa $aOverL $bOverL <$Wyb $Wxb>` | `LadrunoDispBeamColumn2d`, `LadrunoForceBeamColumn2d`, `LadrunoIMKBeam2d`, `ForceBeamColumn2d` | Cancellation-free shifted-origin moment integrals over $[a, b] \subset [0, L]$. |
| **2D Point Load** (`LOAD_TAG_Beam2dPointLoad`) | `eleLoad -ele $eleTag -type -beamPoint $Py $aOverL <$Px>` | `LadrunoDispBeamColumn2d`, `LadrunoForceBeamColumn2d`, `LadrunoIMKBeam2d`, `ForceBeamColumn2d` | Point load at $x = a/L \in [0, 1]$. |
| **3D Uniform Load** (`LOAD_TAG_Beam3dUniformLoad`) | `eleLoad -ele $eleTag -type -beamUniform $Wy $Wz <$Wx>` | `LadrunoDispBeamColumn3d`, `LadrunoForceBeamColumn3d`, `LadrunoIMKBeam`, `ForceBeamColumn3d` | `$Wy, $Wz` transverse along local $y, z$; `$Wx` axial. |
| **3D Partial / Trapezoidal Load** (`LOAD_TAG_Beam3dPartialUniformLoad`) | `eleLoad -ele $eleTag -type -beamUniform $Wya $Wza $Wxa $aOverL $bOverL <$Wyb $Wzb $Wxb>` | `LadrunoDispBeamColumn3d`, `LadrunoForceBeamColumn3d`, `LadrunoIMKBeam`, `ForceBeamColumn3d` | Trapezoidal 3D distributed load over $[a/L, b/L]$; used directly by `LadrunoMembrane` two-way slab tributary load distribution! |
| **3D Point Load** (`LOAD_TAG_Beam3dPointLoad`) | `eleLoad -ele $eleTag -type -beamPoint $Py $Pz $aOverL <$Px>` | `LadrunoDispBeamColumn3d`, `LadrunoForceBeamColumn3d`, `LadrunoIMKBeam`, `ForceBeamColumn3d` | 3D concentrated load at $x = a/L \in [0, 1]$. |
| **Uniform Distributed Moment** (`LOAD_TAG_BeamUniformMoment`) | *(Internal `BeamUniformMoment` elemental load)* | `LadrunoDispBeamColumn2d/3d`, `LadrunoIMKBeam2d/3d` | Adds distributed couple reactions `mz` (and `my` in 3D) to `p0`. |

> **Important Inelastic Behavior (`LadrunoIMKBeam2d` / `LadrunoIMKBeam`):** Fixed-end moments `q0` enter inside the $2\times 2$ (2D) or $4\times 4$ (3D) hinge static-condensation loop:
> $$M_{\text{end}} = \mathbf{K}_{\text{el}}\,(\mathbf{v} - \boldsymbol{\theta}_h) + \mathbf{q}_0 = \mathbf{M}_h(\boldsymbol{\theta}_h)$$
> Consequently, a fixed-fixed `LadrunoIMKBeam` subjected to a large gravity span load $w L^2 / 12 > M_y$ will properly yield its end hinges even when the end nodes do not rotate ($\theta_I = \theta_J = 0$).

---

## 5. `setParameter` / `updateParameter` Reference

| Element / Material | Direct Parameter Names | Forwarded Sub-Component Syntax |
| :--- | :--- | :--- |
| `LadrunoDispBeamColumn2d`<br>`LadrunoDispBeamColumn3d`<br>`LadrunoForceBeamColumn2d`<br>`LadrunoForceBeamColumn3d` | `'rho'`, `'lch'`, `'GAy'` / `'GAvy'` (`LadrunoForceBeamColumn`), `'GAz'` / `'GAvz'` (3D `LadrunoForceBeamColumn`) | • `'section' $secNum $arg...` (1-based section index)<br>• `'sectionX' $x $arg...` (physical coordinate $x \in [0,L]$)<br>• `'hinge'` / `'hingeZ'` / `'hingeIJ'` / `'hingeZIJ' $arg...` (single hinge or both dual-$z$ hinges)<br>• `'hingeI'` / `'hingeZI'` / `'hinge1' $arg...` (end $I$)<br>• `'hingeJ'` / `'hingeZJ'` / `'hinge2' $arg...` (end $J$)<br>• `'hingeY'` / `'hingeYIJ'` / `'hingeYI'` / `'hingeY1'` / `'hingeYJ'` / `'hingeY2' $arg...` (**3D only**)<br>• `'hingeBiaxial'` / `'hingeC'` / `'hingeBiaxialIJ'` / `'hingeCIJ'` / `'hingeBiaxialI'` / `'hingeCI'` / `'hingeBiaxial1'` / `'hingeBiaxialJ'` / `'hingeCJ'` / `'hingeBiaxial2' $arg...` (**3D only**)<br>• *(no prefix)*: broadcasts `$argv` to all sections and active hinges |
| `LadrunoIMKBeam2d` | `'E'`, `'A'`, `'Iz'` (or `'I'`), `'rho'`, `'GAvy'` (or `'GAy'`, `'GAv'`), `'P0'`, `'pmAlpha'` | • `'hingeI'` / `'hinge1'` / `'hingeZi'` / `'hingeZI'` / `'hingeMzI'` / `'matZi' $arg...`<br>• `'hingeJ'` / `'hinge2'` / `'hingeZj'` / `'hingeZJ'` / `'hingeMzJ'` / `'matZj' $arg...`<br>• `'matZ'` / `'hingeZ'` / `'material'` / `'hinge' $arg...` (both ends, or `'hinge' i|j|0|1 $arg...`)<br>• `'shearHinge'` / `'shearHingeY'` / `'matShear' $arg...` |
| `LadrunoIMKBeam` (3D) | `'E'`, `'A'`, `'Iz'`, `'Iy'`, `'G'`, `'Jx'` (or `'J'`), `'rho'`, `'GAvy'` (or `'GAy'`), `'GAvz'` (or `'GAz'`), `'P0'`, `'pmAlpha'` | • `'hingeI'` / `'hinge1'` / `'hingeZi'` / `'hingeZI'` / `'hingeMzI'` / `'matZi' $arg...`<br>• `'hingeJ'` / `'hinge2'` / `'hingeZj'` / `'hingeZJ'` / `'hingeMzJ'` / `'matZj' $arg...`<br>• `'hingeYi'` / `'hingeYI'` / `'hingeMyI'` / `'matYi' $arg...`<br>• `'hingeYj'` / `'hingeYJ'` / `'hingeMyJ'` / `'matYj' $arg...`<br>• `'matZ'` / `'hingeZ' $arg...` (both $z$ hinges)<br>• `'matY'` / `'hingeY' $arg...` (both $y$ hinges)<br>• `'hingeBiaxial'` / `'hingeBiaxialI'` / `'hingeBiaxialJ' $arg...`<br>• `'shearHingeY'` / `'shearHingeZ'` / `'shearHinge' $arg...`<br>• `'material'` / `'hinge' $arg...` (all 4 hinges, or `'hinge' Zi|Zj|Yi|Yj|0|1|2|3 $arg...`) |
| `ForceBeamColumn2d`<br>`ForceBeamColumn3d` | `'rho'`, `'lch'` | • `'section' $secNum $arg...`<br>• `'sectionX' $x $arg...`<br>• `'allSections' $arg...` (or default fallback without prefix) |
| `LadrunoCohesiveHinge` | `'Mc'`, `'Gf'`, `'Kpen'`, `'P0'`, `'alphaPM'` | — |
| `LadrunoCohesiveHingeBiaxial` | `'Mcz'` (or `'Mc'`), `'Gfz'` (or `'Gf'`), `'Mcy'`, `'Gfy'`, `'P0'`, `'alphaPM'`, `'bkEta'` (or `'eta'`), `'penaltyRatio'` | — |

---

## 6. `eleResponse` & Material Response Reference

All frame and wall elements share the canonical token matcher in [`SRC/element/LadrunoResponseTokens.h`](file:///C:/Users/nmora/.gemini/antigravity/worktrees/OpenSees/review_ladrunobeam_implementation/SRC/element/LadrunoResponseTokens.h).

### 6.1 `LadrunoDispBeamColumn2d/3d` & `LadrunoForceBeamColumn2d/3d`

| Query Token(s) | 2D Size & Layout | 3D Size & Layout |
| :--- | :--- | :--- |
| `'forces'` / `'force'` / `'globalForce'` / `'globalForces'` | `6`: `[Px_i, Py_i, Mz_i, Px_j, Py_j, Mz_j]` | `12`: `[Px_i, Py_i, Pz_i, Mx_i, My_i, Mz_i, Px_j, Py_j, Pz_j, Mx_j, My_j, Mz_j]` |
| `'localForce'` / `'localForces'` | `6`: `[N_i, V_i, M_i, N_j, V_j, M_j]` | `12`: `[N_i, Vy_i, Vz_i, T_i, My_i, Mz_i, N_j, Vy_j, Vz_j, T_j, My_j, Mz_j]` |
| `'basicForce'` / `'basicForces'` / `'chordForce'` | `3`: `[N, Mz_i, Mz_j]` | `6`: `[N, Mz_i, Mz_j, My_i, My_j, T]` |
| `'basicStiffness'` | `9`: flat $3\times 3$ condensed basic stiffness | `36`: flat $6\times 6$ condensed basic stiffness |
| `'chordRotation'` / `'chordDeformation'` / `'basicDeformation'` | `3`: `[eps, thetaZ_i, thetaZ_j]` | `6`: `[eps, thetaZ_i, thetaZ_j, thetaY_i, thetaY_j, thetaX]` |
| `'plasticRotation'` / `'plasticDeformation'` / `'hingeRotation'` / `'hingeRotations'` | `3`: `[epsP, thetaZ_iP, thetaZ_jP]` | `6`: `[epsP, thetaZ_iP, thetaZ_jP, thetaY_iP, thetaY_jP, thetaXP]` |
| `'hingeSummary'` / `'hingeState'` | `4`: `[theta_i, M_i, theta_j, M_j]` | `8`: `[thetaZ_i, Mz_i, thetaZ_j, Mz_j, thetaY_i, My_i, thetaY_j, My_j]` |
| `'hingeJump'` / `'hingeJumpI'` | `1`: `[theta_i]` | `2`: `[thetaZ_i, thetaY_i]` |
| `'hingeJumpJ'` | `1`: `[theta_j]` | `2`: `[thetaZ_j, thetaY_j]` |
| `'hinge'` / `'hingeI'` / `'hingeJ'` `$matResp...` | Forwards `$matResp` to major-axis ($z$) hinge material | Forwards `$matResp` to major-axis ($z$) hinge material |
| `'hingeY'` / `'hingeYI'` / `'hingeYJ'` `$matResp...` | — | Forwards `$matResp` to minor-axis ($y$) hinge material |
| `'hingeBiaxial'` / `'hingeBiaxialI'` / `'hingeBiaxialJ'` `$matResp...` | — | Forwards `$matResp` to biaxial hinge `NDMaterial` |
| `'implexError'` / `'ImplexError'` | `1`: `[maxImplexErr]` across active hinges | `1`: `[maxImplexErr]` across active hinges |
| `'avgImplexError'` / `'AvgImplexError'` | `1`: `[avgImplexErr]` across active hinges | `1`: `[avgImplexErr]` across active hinges |
| `'integrationPoints'` | `N`: `[x_1, ..., x_N]` (physical coords $\xi_i L$) | `N`: `[x_1, ..., x_N]` (physical coords $\xi_i L$) |
| `'integrationWeights'` | `N`: `[w_1, ..., w_N]` (physical weights $w_i L$) | `N`: `[w_1, ..., w_N]` (physical weights $w_i L$) |
| `'section'` `$secNum` `$secResp...` | Forwards to 1-based section `$secNum` (`'force'`, `'deformation'`, `'fiber'`, etc.) | Forwards to 1-based section `$secNum` |
| `'sectionX'` `$x` `$secResp...` | Forwards to section nearest physical coordinate $x \in [0, L]$ ($\xi = x/L$) | Forwards to section nearest physical coordinate $x \in [0, L]$ ($\xi = x/L$) |
| `'energy'` / `'Energy'` | `1`: integrated section energy $\sum w_i L E_i$ **plus** active embedded hinge energy (`theHinge->getEnergy()`) | `1`: integrated section energy $\sum w_i L E_i$ **plus** active embedded hinge energy (`theHinge*->getEnergy()`) |
| `'localAxes'` / `'localaxes'` | `9`: `[vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z]` | `9`: `[vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z]` |
| `'dampingForce'` / `'dampingForces'` | `6` (if elemental damping active) | `12` (if elemental damping active) |

---

### 6.2 `LadrunoIMKBeam2d` & `LadrunoIMKBeam` (3D)

| Query Token(s) | 2D Size & Layout | 3D Size & Layout |
| :--- | :--- | :--- |
| `'forces'` / `'force'` / `'globalForce'` / `'globalForces'` | `6`: `[Px_i, Py_i, Mz_i, Px_j, Py_j, Mz_j]` | `12`: `[Px_i, Py_i, Pz_i, Mx_i, My_i, Mz_i, Px_j, Py_j, Pz_j, Mx_j, My_j, Mz_j]` |
| `'localForce'` / `'localForces'` | `6`: `[N_i, V_i, M_i, N_j, V_j, M_j]` | `12`: `[N_i, Vy_i, Vz_i, T_i, My_i, Mz_i, N_j, Vy_j, Vz_j, T_j, My_j, Mz_j]` |
| `'basicForce'` / `'basicForces'` / `'chordForce'` | `3`: `[N, Mz_i, Mz_j]` | `6`: `[N, Mz_i, Mz_j, My_i, My_j, T]` |
| `'basicStiffness'` | `9`: flat $3\times 3$ condensed basic stiffness | `36`: flat $6\times 6$ condensed basic stiffness |
| `'chordRotation'` / `'chordDeformation'` / `'basicDeformation'` | `3`: `[eps, thetaZ_i, thetaZ_j]` | `6`: `[eps, thetaZ_i, thetaZ_j, thetaY_i, thetaY_j, thetaX]` |
| `'hingeRotation'` / `'hingeRotations'` / `'plasticRotation'` / `'plasticDeformation'` | `2`: `[thetaH_Zi, thetaH_Zj]` | `4`: `[thetaH_Zi, thetaH_Zj, thetaH_Yi, thetaH_Yj]` |
| `'hingeMoment'` / `'hingeMoments'` | `2`: `[M_Zi, M_Zj]` | `4`: `[M_Zi, M_Zj, M_Yi, M_Yj]` |
| `'hingeSummary'` / `'hingeState'` | `4`: `[thetaH_Zi, M_Zi, thetaH_Zj, M_Zj]` | `8`: `[thetaH_Zi, M_Zi, thetaH_Zj, M_Zj, thetaH_Yi, M_Yi, thetaH_Yj, M_Yj]` |
| `'shearHinge'` *(1 arg)* | `2`: `[gamma_y, Vy]` | `4`: `[gamma_y, Vy, gamma_z, Vz]` |
| `'shearHingeY'` / `'shearHingeZ'` *(1 arg)* | `2`: `[gamma_y, Vy]` (`'shearHingeY'`) | `2`: `[gamma_y, Vy]` or `[gamma_z, Vz]` |
| `'shearHinge'` / `'shearHingeY'` / `'shearHingeZ'` `$matResp...` | Forwards `$matResp` to shear hinge `UniaxialMaterial` | Forwards `$matResp` to shear hinge `UniaxialMaterial` |
| `'hingeBiaxial'` / `'hingeBiaxialI'` / `'hingeBiaxialJ'` `$matResp...` | — | Forwards `$matResp` to biaxial hinge `NDMaterial` (or returns `[thetaZ, Mz, thetaY, My]` summary when 1 arg) |
| `'section'` `1\|2` `'force'` | `2`: `[N_end, Mz_end]` at end `1` ($I$) or `2` ($J$) | `4`: `[N_end, Mz_end, My_end, T_end]` at end `1` ($I$) or `2` ($J$) |
| `'material'` / `'hinge'` `$idx` `$matResp...` | Forwards `$matResp` to 0-based slot `$idx` (`0`=Zi, `1`=Zj, `'shear'`=shearHinge) | Forwards `$matResp` to 0-based slot `$idx` (`0`=Zi, `1`=Zj, `2`=Yi, `3`=Yj, `'biaxialI'`, `'biaxialJ'`, `'shearY'`, `'shearZ'`) |
| `'implexError'` / `'ImplexError'` | `1`: `[maxImplexErr]` across active hinges | `1`: `[maxImplexErr]` across active hinges |
| `'avgImplexError'` / `'AvgImplexError'` | `1`: `[avgImplexErr]` across active hinges | `1`: `[avgImplexErr]` across active hinges |
| `'localAxes'` / `'localaxes'` | `9`: `[vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z]` | `9`: `[vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z]` |
| `'energy'` / `'Energy'` | `1`: elastic strain energy + hinge work | `1`: elastic strain energy + hinge work |

---

### 6.3 `MVLEM`, `SFI_MVLEM`, `MVLEM_3D`, `SFI_MVLEM_3D`

| Query Token(s) | `MVLEM` (2D) | `SFI_MVLEM` (2D) | `MVLEM_3D` (3D) | `SFI_MVLEM_3D` (3D) |
| :--- | :--- | :--- | :--- | :--- |
| `'forces'` / `'force'` / `'globalForce'` / `'globalForces'` *(+ `'forceG'`/`'forcesG'` in 3D)* | `6`: `[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]` | `6`: `[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]` | `24`: 4-node global `[Fx, Fy, Fz, Mx, My, Mz] * 4` | `24`: 4-node global `[Fx, Fy, Fz, Mx, My, Mz] * 4` |
| `'localForce'` / `'localForces'` *(+ `'forceL'`/`'forcesL'` in 3D)* | `6` (2D global = local) | `6` (2D global = local) | `24`: 4-node local forces | `24`: 4-node local forces |
| `'Curvature'` / `'curvature'` | `1`: scalar wall curvature $\phi$ | `1`: scalar wall curvature $\phi$ | `1`: scalar wall curvature $\phi$ | `1`: scalar wall curvature $\phi$ |
| `'ShearDef'` / `'sheardef'` | `1`: shear deformation $D_{sh}$ | `1`: shear deformation $D_{sh}$ | `1`: shear deformation $D_{sh}$ | `1`: shear deformation $D_{sh}$ |
| `'Fiber_Strain'` / `'fiber_strain'` | `m`: vertical fiber strains $\varepsilon_y$ | *(query via panel `'strain'`)* | `m`: vertical fiber strains $\varepsilon_y$ | *(query via panel `'strain'`)* |
| `'Fiber_Stress_Concrete'` / `'fiber_stress_concrete'` | `m`: concrete stresses $\sigma_{yc}$ | — | `m`: concrete stresses $\sigma_{yc}$ | — |
| `'Fiber_Stress_Steel'` / `'fiber_stress_steel'` | `m`: steel stresses $\sigma_{ys}$ | — | `m`: steel stresses $\sigma_{ys}$ | — |
| `'Shear_Force_Deformation'` / `'shear_force_deformation'` | `2`: `[shearDef, shearFrc]` | — | `2`: `[shearDef, shearFrc]` | — |
| `'integrationWeights'` | `m`: fiber widths `[b_1, ..., b_m]` | `m`: panel widths `[b_1, ..., b_m]` | `m`: fiber widths `[b_1, ..., b_m]` | `m`: panel widths `[b_1, ..., b_m]` |
| `'localAxes'` / `'localaxes'` | `9`: identity `[1,0,0, 0,1,0, 0,0,1]` | `9`: identity `[1,0,0, 0,1,0, 0,0,1]` | `9`: rows of $3\times 3$ wall rotation $\mathbf{T}_t$ | `9`: rows of $3\times 3$ wall rotation $\mathbf{T}_t$ |
| `'material'` *(or `'RCpanel'` / `'RCPanel'` / `'RC_panel'` / `'RC_Panel'`)* `$idx` ... | `1..m` `['concrete'\|'steel'\|'shear']` `$resp` | `1..m` `$panelResp` (e.g. `'strain'`, `'stress'`) | `'material'` `1..m` `['concrete'\|'steel'\|'shear']` `$resp` | `1..m` `$panelResp` (e.g. `'strain'`, `'stress'`) |

---

### 6.4 Hinge Material Responses (`LadrunoCohesiveHinge` & `LadrunoCohesiveHingeBiaxial`)

| Material | Query Token(s) | Size & Description |
| :--- | :--- | :--- |
| `LadrunoCohesiveHinge` | `'stress'` / `'force'` | `1`: Current hinge moment $M$ |
| | `'strain'` / `'deformation'` / `'rotationJump'` / `'jump'` / `'kappa'` | `1`: Rotational jump $[\![\theta]\!]$ |
| | `'tangent'` / `'stiffness'` | `1`: Current rotational tangent $\partial M / \partial [\![\theta]\!]$ |
| | `'maxJump'` / `'kappaMax'` | `1`: `[TkappaMax]` (maximum rotational jump norm) |
| | `'dissipatedEnergy'` / `'energy'` / `'work'` | `1`: Cumulative dissipated fracture energy |
| | `'damage'` / `'D'` | `1`: Scalar cohesive damage index $D \in [0, 1]$ |
| | `'failed'` / `'open'` | `1`: `1.0` if completely fractured ($D \ge 1$), else `0.0` |
| | `'pmMultiplier'` / `'fpm'` | `1`: Current $P\text{--}M$ capacity factor $f_{pm}(N) \in [0.05, 2.5]$ |
| | `'implexError'` / `'ImplexError'` | `1`: Relative IMPLEX extrapolation error $\|\tilde{M} - M_{\text{imp}}\| / (\|M_{\text{imp}}\| + 10^{-12})$ |
| `LadrunoCohesiveHingeBiaxial` | `'stress'` / `'stresses'` | `2`: `[Mz, My]` |
| | `'strain'` / `'strains'` | `2`: `[alphaZ, alphaY]` |
| | `'momentZ'` / `'Mz'` ; `'momentY'` / `'My'` | `1`: Individual moment component $M_z$ or $M_y$ |
| | `'jumpZ'` / `'alphaZ'` ; `'jumpY'` / `'alphaY'` | `1`: Individual rotational jump $\alpha_z$ or $\alpha_y$ |
| | `'tangent'` / `'stiffness'` | `4`: Flat $2\times 2$ tangent `[Kzz, Kzy, Kyz, Kyy]` |
| | `'maxNorm'` / `'rMax'` | `1`: Maximum normalized jump norm $r_{\max}$ |
| | `'dissipatedEnergy'` / `'energy'` / `'work'` | `1`: Cumulative biaxial dissipated fracture energy |
| | `'damage'` / `'D'` | `1`: Scalar biaxial cohesive damage index $D \in [0, 1]$ |
| | `'pmMultiplier'` / `'fpm'` | `1`: Current $P\text{--}M$ capacity factor $f_{pm}(N) \in [0.05, 2.5]$ |
| | `'implexError'` / `'ImplexError'` | `1`: Relative biaxial IMPLEX error |

---

## 7. `.ladruno` HDF5 Recorder Schema & Post-Processing (`Ladruno_ElementResults.h`)

### 7.1 Recorder Command Syntax
```tcl
recorder ladruno $filePath.ladruno \
    -N displacement reaction \
    -E force localForce sectionForce sectionDeformation endForces hingeState localAxes \
    <-dt $deltaT> <-flush $flushEvery> <-compress $gzipLevel>
```

### 7.2 Element Topology & Quadrature Mapping in HDF5

Defined in [`SRC/recorder/Ladruno_ElementResults.h`](file:///C:/Users/nmora/.gemini/antigravity/worktrees/OpenSees/review_ladrunobeam_implementation/SRC/recorder/Ladruno_ElementResults.h):

| Element Class | `ElementGeometryType` | `ElementIntegrationRuleType` | Quadrature Metadata Written to `/MODEL_STAGE[i]/QUADRATURE/` |
| :--- | :--- | :--- | :--- |
| `LadrunoDispBeamColumn2d` (`33013`)<br>`LadrunoDispBeamColumn3d` (`33014`)<br>`LadrunoForceBeamColumn2d` (`33024`)<br>`LadrunoForceBeamColumn3d` (`33025`)<br>`ForceBeamColumn2d` (`73`)<br>`ForceBeamColumn3d` (`74`) | `Line_2N` (`1`) | `CustomIntegrationRule` (`100`, `custom_rule_dimension = 1`) | Queries `'integrationPoints'`, normalizes $\xi_k = 2(x_k - x_{\min}) / \Delta x - 1 \in [-1, 1]$, and writes `GP_X` (and `GP_WEIGHTS`). |
| `LadrunoIMKBeam` (`33003`)<br>`LadrunoIMKBeam2d` (`33004`) | `Line_2N` (`1`) | `Line_GaussLegendre_1` (`1`, 1 GP at midpoint) | Standard 1-point line rule. |
| `MVLEM` (`162`)<br>`SFI_MVLEM` (`163`) | `Line_2N` (`1`, 2 external nodes) | `CustomIntegrationRule` (`100`, `custom_rule_dimension = 1`) | Probes `'material'` `1..m` `setResponse()` XML metadata (`<GaussPointOutput eta="..." weight="...">`) to extract normalized panel locations and weights (`GP_X`, `GP_WEIGHTS`). |
| `MVLEM_3D` (`212`)<br>`SFI_MVLEM_3D` (`213`) | `Quadrilateral_CohesiveBand_4N` (`10`, 4 external quad nodes) | `CustomIntegrationRule` (`100`, `custom_rule_dimension = 2`) | Probes `'material'` `1..m` `setResponse()` XML metadata (`<GaussPointOutput eta="..." weight="...">`) and writes 2D custom rule (`GP_X`, `GP_Y = 0`, `GP_WEIGHTS`). |

### 7.3 `-E` Element Dataset Shapes & Semantics

All element results are stored under `/MODEL_STAGE[i]/RESULTS/ON_ELEMENTS/<ResultName>/<ElementGroup>/DATA` with shape `(N_steps, N_elements_in_group, N_components)`:

| `-E` Result Token | Internal Query | Components (`N_components`) | Notes |
| :--- | :--- | :--- | :--- |
| `force` | `'forces'` | 2D beam/wall: `6`<br>3D beam: `12`<br>3D wall (`MVLEM_3D`/`SFI_MVLEM_3D`): `24` | Global nodal resisting forces. |
| `localForce` | `'localForce'` | 2D beam/wall: `6`<br>3D beam: `12`<br>3D wall: `24` | Element-local nodal resisting forces. |
| `endForces` | `'forces'` | 2D beam: `6`<br>3D beam: `12` | Whole-element end forces (`isWholeElementItem() == true`, 1 vector per element, never multiplied by `num_gauss_points`). |
| `hingeState` | `'hingeSummary'` | 2D beam (`LadrunoDispBeamColumn2d`, `LadrunoForceBeamColumn2d`, `LadrunoIMKBeam2d`): `4` (`[th_i, M_i, th_j, M_j]`)<br>3D beam (`LadrunoDispBeamColumn3d`, `LadrunoForceBeamColumn3d`, `LadrunoIMKBeam`): `8` (`[thZ_i, Mz_i, thZ_j, Mz_j, thY_i, My_i, thY_j, My_j]`) | Whole-element hinge summary (`isWholeElementItem() == true`). Works uniformly across `LadrunoDispBeamColumn`, `LadrunoForceBeamColumn`, and `LadrunoIMKBeam`! |
| `sectionForce` | `'section' $gp 'force'` | `N_gp * nSecForce` | Per-integration-point section stress resultants (`LadrunoDispBeamColumn`, `LadrunoForceBeamColumn`, `ForceBeamColumn`; for `LadrunoIMKBeam`, queries `'section' 1 'force'`). |
| `sectionDeformation` | `'section' $gp 'deformation'` | `N_gp * nSecDef` | Per-integration-point section deformations (`LadrunoDispBeamColumn`, `LadrunoForceBeamColumn`, `ForceBeamColumn`). |
| `localAxes` | `'localAxes'` | `9`: `[vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z]` | Written once as static element orientation metadata and/or per step. |

---

## 8. OpenMP Multi-Threading (`ladrunoThreads`) & IMPLEX Guidelines

1. **Thread-Safe Element State Update & Tangent Assembly (`WP-192` / `WP-196` / `WP-200` / `WP-201`):**
   * `LadrunoDispBeamColumn2d`, `LadrunoDispBeamColumn3d`, `LadrunoForceBeamColumn2d`, `LadrunoForceBeamColumn3d`, `LadrunoIMKBeam2d`, `LadrunoIMKBeam`, `ForceBeamColumn2d`, and `ForceBeamColumn3d` override `ladrunoThreadSafeUpdate()` and `ladrunoThreadSafeForm()` (returning `true` when sections, coordinate transformations, and hinges are thread-safe).
   * Static scratch matrices/vectors inside these elements and their transformations (`LinearCrdTransf2d/3d`, `PDeltaCrdTransf2d/3d`, `CorotCrdTransf2d/3d`) and sections (`ElasticSection2d/3d`, `FiberSection2d/3d`) are `static thread_local`, guaranteeing bit-identical results under `ops.ladrunoThreads(1)`, `2`, `4`, `8+`.
   * **Uniaxial / Fiber Material Thread-Safety Requirement:** `FiberSection2d/3d` and `SectionAggregator` only report `ladrunoThreadSafe() == true` when every underlying `UniaxialMaterial` overrides `ladrunoThreadSafe() { return true; }`:
     - **Parallel-safe (`ladrunoThreadSafe() == true`):** `ElasticMaterial`, `Steel01`, `Steel02`, `Concrete01`, `Concrete02`, `ASDConcrete1D`, `ASDSteel1D`, `LadrunoUniaxialJ2`, `LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `Bilin`, `ModIMKPeakOriented`.
     - **Serial fallback (`ladrunoThreadSafe() == false`):** `Concrete04`, `SteelMPF`, `ReinforcingSteel`, `Hysteretic` (if used inside a `FiberSection`, the parent beam element automatically falls back to serial update/tangent assembly).
   * **Note on `MVLEM` / `SFI_MVLEM`:** The macro-wall elements use class-level or non-thread-marked `NDMaterial` (`FSAM`) instances and retain the base `Element::ladrunoThreadSafeUpdate() == false` serial fallback in hybrid models.
2. **IMPLEX (`algorithm Implex` / `-implex`):**
   * `LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `Bilin`, `ModIMKPeakOriented`, `LadrunoUniaxialJ2`, `ASDConcrete1D`, and `ASDSteel1D` all support `-implex` (as do `LadrunoDispBeamColumn`, `LadrunoForceBeamColumn`, and `LadrunoIMKBeam` at the element level; note that `LadrunoRebarBuckling` and `LadrunoBondSlip` do not take `-implex`).
   * When emitting an IMPLEX model from `apeGmsh`, pass `-implex` on the softening/plastic materials (and/or beam elements) and monitor element-level `'implexError'` / `'avgImplexError'` responses to drive adaptive time-stepping.

---

## 9. `apeGmsh` High-Level Python API (`ops.element.*`) Reference & Examples

In `apeGmsh` (`src/apeGmsh/opensees/element/beam_column.py` and `src/apeGmsh/opensees/element/_namespace.py`), the three Ladruno frame element classes expose typed Python keyword arguments that validate mutual-exclusion rules at model-build time and emit the exact Tcl/Live argument streams specified in Section 1.

### 9.1 `ops.element.LadrunoDispBeamColumn` (`LadrunoDispBeamColumn`)

| Python Keyword Argument | Type | Default | Emitted C++ / Tcl Flag(s) |
| :--- | :--- | :--- | :--- |
| `tag`, `i_node`, `j_node`, `transf_tag`, `integration_tag` | `int` | *Required* | Positional `tag iNode jNode transfTag integrationTag` |
| `c_mass` | `bool` | `False` | `'-cMass'` |
| `mass` | `float` | `0.0` | `'-mass', mass` (emitted when `mass > 0.0`) |
| `nl` | `bool` | `False` | `'-nl'` (mutually exclusive with all hinge flags) |
| `span_load` | `bool` | `True` | `'-noSpanLoad'` when `span_load=False` (coupled $(\beta_P, \beta_\kappa)$ bubble is active by default; automatically disabled in C++ when `-nl` is active) |
| `lch` | `float \| None` | `None` | `'-lch', lch` (must be `> 0`) |
| `hinge`, `hinge_i`, `hinge_j` | `int \| None` | `None` | `'-hinge', tag` / `'-hingeI', tag` / `'-hingeJ', tag` |
| `hinge_ij` | `tuple[int, int] \| None` | `None` | `'-hingeIJ', tag_i, tag_j` |
| `hinge_xi` | `float \| None` | `None` | `'-hingeXi', xi` ($\xi \in [0, 1]$) |
| `hinge_loc` | `float \| tuple[float, float] \| None` | `None` | `'-hingeLoc', xi_i[, xi_j]` ($\xi \in [0, 1]$) |
| `hinge_y`, `hinge_yi`, `hinge_yj` | `int \| None` | `None` | `'-hingeY', tag` / `'-hingeYI', tag` / `'-hingeYJ', tag` (3D only) |
| `hinge_yij` | `tuple[int, int] \| None` | `None` | `'-hingeYIJ', tag_i, tag_j` (3D only) |
| `hinge_biaxial`, `hinge_biaxial_i`, `hinge_biaxial_j` | `int \| None` | `None` | `'-hingeBiaxial', tag` / `'-hingeBiaxialI', tag` / `'-hingeBiaxialJ', tag` (3D only; mutually exclusive with uniaxial `hinge*`/`hinge_y*`) |
| `hinge_biaxial_ij` | `tuple[int, int] \| None` | `None` | `'-hingeBiaxialIJ', tag_i, tag_j` (3D only) |
| `pm` | `tuple[float, float] \| float \| None` | `None` | `'-pm', P0, alpha` (if `float`, emits `alpha = 1.0`) |
| `implex` | `bool` | `False` | `'-implex'` |

### 9.2 `ops.element.LadrunoForceBeamColumn` (`LadrunoForceBeamColumn`)

| Python Keyword Argument | Type | Default | Emitted C++ / Tcl Flag(s) |
| :--- | :--- | :--- | :--- |
| `tag`, `i_node`, `j_node`, `transf_tag`, `integration_tag` | `int` | *Required* | Positional `tag iNode jNode transfTag integrationTag` |
| `mass` | `float` | `0.0` | `'-mass', mass` (emitted when `mass > 0.0`) |
| `iter` | `tuple[int, float] \| None` | `None` | `'-iter', max_iters, tol` |
| `subdivide` | `int \| tuple[int, float] \| None` | `None` | `'-subdivide', max_levels[, min_sub_factor]` |
| `lch` | `float \| None` | `None` | `'-lch', lch` (must be `> 0`) |
| `total_element_lch` | `bool` | `False` | `'-totalElementLch'` |
| `shear` | `float \| tuple[float, float] \| None` | `None` | 2D: `'-shear', GAy` (`float`)<br>3D: `'-shear', GAy, GAz` (`tuple[float, float]`, or `(GA, GA)` if `float`) |
| `hinge`, `hinge_i`, `hinge_j`, `hinge_ij`, `hinge_xi`, `hinge_loc`, `hinge_y*`, `hinge_biaxial*`, `pm`, `implex` | Same as `LadrunoDispBeamColumn` | `None` / `False` | Same embedded cohesive hinge flags as `LadrunoDispBeamColumn` |

### 9.3 `ops.element.LadrunoIMKBeam` (`LadrunoIMKBeam`)

| Python Keyword Argument | Type | Default | Emitted C++ / Tcl Flag(s) |
| :--- | :--- | :--- | :--- |
| `tag`, `i_node`, `j_node`, `E`, `A`, `Iz`, `transf_tag` | `int \| float` | *Required* | 2D: `tag iNode jNode E A Iz transfTag` |
| `Iy`, `G`, `J` | `float \| None` | `None` | 3D: `tag iNode jNode E A Iz Iy G J transfTag` (all three required in 3D) |
| `mat_z`, `mat_y`, `mat_torsion` | `int \| tuple[int, int] \| None` | `None` | `'-matZ'`, `'-matY'` (3D), `'-matTorsion'` (3D) |
| `hinge_biaxial`, `hinge_biaxial_i`, `hinge_biaxial_j` | `int \| None` | `None` | `'-hingeBiaxial'`, `'-hingeBiaxialI'`, `'-hingeBiaxialJ'` (3D only; mutually exclusive with `mat_z`/`mat_y`) |
| `hinge_biaxial_ij` | `tuple[int, int] \| None` | `None` | `'-hingeBiaxialIJ', tag_i, tag_j` (3D only) |
| `shear_hinge` | `int \| None` | `None` | 2D: `'-shearHinge', tag` (`UniaxialMaterial` tag for inelastic transverse shear $V_y\text{--}\Delta u_y$) |
| `shear_hinge_y`, `shear_hinge_z` | `int \| None` | `None` | 3D: `'-shearHingeY', tag` ($V_y\text{--}\Delta u_y$), `'-shearHingeZ', tag` ($V_z\text{--}\Delta u_z$) |
| `release`, `release_z`, `release_y` | `int` | `0` | 2D: `'-release', flag` (`1`=I, `2`=J, `3`=IJ)<br>3D: `'-releaseZ', flag`, `'-releaseY', flag` |
| `shear` | `float \| tuple[float, float] \| None` | `None` | 2D: `'-shear', GAy`<br>3D: `'-shear', GAy, GAz` |
| `end_zones` | `tuple[float, float] \| None` | `None` | `'-endZones', dI, dJ` |
| `pm`, `pm_z`, `pm_y` | `tuple[float, float] \| float \| None` | `None` | `'-pm', P0, alpha` / `'-pmZ', P0, alpha` / `'-pmY', P0, alpha` |
| `c_mass`, `implex`, `mass`, `iter` | `bool \| float \| tuple` | `False` / `0.0` / `None` | `'-cMass'`, `'-implex'`, `'-mass'`, `'-iter'` |

### 9.4 Python `apeGmsh` Usage Examples

```python
from apeGmsh.opensees.model import Model

model = Model(ndm=3, ndf=6)

# 1. Force-based fiber column with Tier-1 lch, Timoshenko shear, and dual 3D biaxial hinges
model.element.LadrunoForceBeamColumn(
    tag=101,
    i_node=1,
    j_node=2,
    transf_tag=1,
    integration_tag=10,
    lch=0.45,
    shear=(8.5e7, 8.5e7),
    hinge_biaxial_ij=(501, 502),
    hinge_loc=(0.0, 1.0),
    pm=(2500.0e3, 0.35),
    subdivide=(10, 1.0e-12),
    implex=True,
)

# 2. Displacement-based fiber beam with coupled (beta_P, beta_kappa) spanLoad bubble (default)
model.element.LadrunoDispBeamColumn(
    tag=201,
    i_node=2,
    j_node=3,
    transf_tag=2,
    integration_tag=20,
    lch=0.50,
    span_load=True,  # default True; set False only to emit '-noSpanLoad'
    c_mass=True,
)

# 3. 3D concentrated-plasticity IMK beam-column with biaxial hinges, shear hinges, and rigid end zones
model.element.LadrunoIMKBeam(
    tag=301,
    i_node=4,
    j_node=5,
    E=30.0e9,
    A=0.25,
    Iz=0.0052,
    Iy=0.0052,
    G=12.5e9,
    J=0.0085,
    transf_tag=1,
    hinge_biaxial_ij=(601, 602),
    shear_hinge_y=701,
    shear_hinge_z=702,
    shear=(1.04e9, 1.04e9),
    end_zones=(0.30, 0.30),
    pm=(3000.0e3, 0.40),
    c_mass=True,
)
```

---

## 10. Cross-Implementation Accuracy & Performance Benchmarks (Opus-Verified)

All three Ladruno frame element families have been benchmarked against upstream ("vanilla") OpenSees element formulations and multi-element `zeroLength` + `equalDOF` workarounds on identical structural problems.

### 10.1 Accuracy: Single-Element Span-Load & Section Response Parity

On an **unsymmetric cracked RC T-beam** ($L = 5.0\text{ m}$, uniform gravity load $w_y = -80\text{ kN/m}$ via `eleLoad -beamUniform`, axial compression $N = -100\text{ kN}$, `Concrete01` + `Steel01` fiber section with shifting neutral axis and strong $P\text{--}M$ coupling):

| Element Formulation | Elements / Mesh | Midspan Moment $M_{\text{mid}}$ ($\text{kN}\cdot\text{m}$) | Midspan Curvature $\kappa_{\text{mid}}$ ($\text{m}^{-1}$) | End Rotation $\theta_I$ ($\text{rad}$) | Midspan Moment Error |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Refined Reference (`dispBeamColumn`)** | **16 elements** (17 nodes) | $250.65$ | $0.004798$ | $-0.007812$ | Reference |
| **Vanilla `dispBeamColumn`** | **1 element** (2 nodes) | $166.67$ | $0.003115$ | $-0.005018$ | **$-33.3\%$** (linear $\kappa(x)$ blind to $w_y$) |
| **`LadrunoDispBeamColumn` (`-noSpanLoad`)** | **1 element** (2 nodes) | $166.67$ | $0.003115$ | $-0.005018$ | **$-33.3\%$** (identical to vanilla) |
| **`LadrunoDispBeamColumn` (default `-spanLoad` bubble)** | **1 element** (2 nodes) | **$247.88$** | **$0.004744$** | **$-0.007882$** | **$-0.85\%$** ($0.00\%$ on linear sections) |
| **Vanilla `forceBeamColumn`** | **1 element** (2 nodes) | **$250.00$** | **$0.004796$** | **$-0.007856$** | **$0.00\%$** (equilibrium exact) |
| **`LadrunoForceBeamColumn`** | **1 element** (2 nodes) | **$250.00$** | **$0.004796$** | **$-0.007856$** | **$0.00\%$** (matches `forceBeamColumn` to $<10^{-13}$) |

> [!TIP]
> With a single element per span, `LadrunoDispBeamColumn`'s default coupled $(\beta_P, \beta_\kappa)$ internal bubble eliminates the $33.3\%$ midspan moment error of vanilla `dispBeamColumn` without requiring internal nodes or force-based element state iterations. When exact force equilibrium is required along the span (or with embedded hinges + Timoshenko shear), `LadrunoForceBeamColumn` gives exact equilibrium parity ($< 10^{-13}$) with vanilla `forceBeamColumn`.

### 10.2 Accuracy: Condensed Series Hinges vs. Multi-Node `zeroLength` + `equalDOF` Assemblies

| Verification Case | Ladruno Single-Element Model (2 nodes, `Plain` constraints) | Vanilla Multi-Element Assembly (`zeroLength` + `equalDOF` + `Transformation`) | Max Relative Difference |
| :--- | :--- | :--- | :---: |
| **Fiber + Shear + Cohesive Hinges** (`LadrunoForceBeamColumn2d` with `-hingeIJ` + `-shear` vs. `forceBeamColumn` + `Section Aggregator Vy` + 2 `zeroLength` + `equalDOF`) | 1 element, 2 nodes, 6 DOFs | 3 elements, 4 nodes, 12 DOFs, 4 `equalDOF` constraints | **$2.93 \times 10^{-16}$** (moment)<br>**$0.00$** (displacement) |
| **IMK Flexural + Timoshenko + Shear Hinges** (`LadrunoIMKBeam2d` with `-matZ` + `-shear` + `-shearHinge` vs. `ElasticTimoshenkoBeam` + 2 rotational `zeroLength` + 1 shear `zeroLength` + `equalDOF`) | 1 element, 2 nodes, 6 DOFs | 4 elements, 5 nodes, 15 DOFs, 6 `equalDOF` constraints | **$1.20 \times 10^{-15}$** (moment)<br>**$5.63 \times 10^{-16}$** (shear)<br>**$0.00$** (displacement) |

### 10.3 Wall-Clock Performance Benchmarks: Multi-Story Frame Pushover & Cyclic Analysis

#### Benchmark A — 2D 8-Bay × 8-Story RC Frame (136 Members, Fiber Sections + Shear + Softening End Hinges, 100 Steps)

| Model Configuration | Nodes / Elements / Constraints | Threads & Solver | Wall Time (ms) | Speedup vs. Vanilla | Limit State Reached |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Vanilla Assembly** (`forceBeamColumn` + `Section Aggregator Vy` + 2 `zeroLength` + `equalDOF`) | 353 nodes, 408 elements, `Transformation` | 1T `Newton` | `1,517.7 ms` | `1.00×` | Snap-back limit $\lambda = 5.86628$ at step 96 |
| **`LadrunoForceBeamColumn2d`** (`-hingeIJ` + `-shear`) | 81 nodes, 136 elements, `Plain` | 1T `Newton` | `377.8 ms` | **`4.02×`** | Exact same limit $\lambda = 5.86628$ at step 96 |
| **`LadrunoForceBeamColumn2d`** (`-hingeIJ` + `-shear`, OpenMP) | 81 nodes, 136 elements, `Plain` | 4T `Newton` | `273.1 ms` | **`5.56×`** | Exact same limit $\lambda = 5.86628$ at step 96 |
| **`LadrunoForceBeamColumn2d`** (`-hingeIJ` + `-shear` + `-implex`, OpenMP) | 81 nodes, 136 elements, `Plain` | 4T `Linear` (`-implex`) | **`168.8 ms`** | **`8.99×`** | Completes all 100 steps (`ok = 0`) through post-peak softening |

#### Benchmark B — 3D 10-Bay × 10-Story Frame (210 Members, Concentrated Plasticity, 150 Cyclic Steps)

| Model Configuration | Nodes / DOFs / Elements / Constraints | Threads & Solver | Wall Time (ms) | Speedup vs. Vanilla |
| :--- | :--- | :--- | :---: | :---: |
| **Vanilla Assembly** (`elasticBeamColumn` + 2 `zeroLength` + `equalDOF`) | 541 nodes (3,246 DOFs), 630 elements, `Transformation` | 1T `Newton` | `7,001.8 ms` | `1.00×` |
| **`LadrunoIMKBeam` 3D** (condensed internal hinges) | 121 nodes (726 DOFs), 210 elements, `Plain` | 1T `Newton` | `1,614.7 ms` | **`4.34×`** |
| **`LadrunoIMKBeam` 3D** (`-implex` + OpenMP) | 121 nodes (726 DOFs), 210 elements, `Plain` | 4T `Linear` (`-implex`) | **`543.5 ms`** | **`12.88×`** |

