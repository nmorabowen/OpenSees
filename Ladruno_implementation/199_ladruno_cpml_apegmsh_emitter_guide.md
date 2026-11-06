---
title: "LadrunoCPML — apeGmsh Emitter, Multi-Layer PML Shell Generator & Coupled Biot u-p Guide"
project: Ladruno
status: ready-for-apegmsh
priority: high
adr: WP 199
tags:
  - apegmsh
  - emitter
  - postprocessor
  - cpml
  - absorbing-boundaries
  - biot-up
  - wave-propagation
---

# `LadrunoCPML` (`ELE_TAG_LadrunoCPML = 33023`) — `apeGmsh` Implementation Guide

This guide specifies how **`apeGmsh`** should generate geometry/mesh layers, emit OpenSees commands, enforce build-time safety gates, and post-process results for the unified 2D/3D Complex-Frequency-Shifted Perfectly Matched Layer (**CFS-PML / C-PML**) absorbing continuum and coupled Biot $u\text{–}p$ wave-propagation element shipped in **WP-199** ([`199_ladruno_cpml_guide.md`](199_ladruno_cpml_guide.md), PRs #995 and #996).

---

## Part 1 — Architectural Comparison: `LadrunoCPML` vs. `ASDAbsorbingBoundary` vs. Legacy `PML2D/3D`

| Aspect | `ASDAbsorbingBoundary2D/3D` (ADR 0054) | Legacy `PML2D` / `PML3D` | `LadrunoCPML` (`ELE_TAG_LadrunoCPML = 33023`) |
|---|---|---|---|
| **Nodal DOFs (`ndf`)** | `ndm` (`2` in 2D, `3` in 3D) | `5` (2D) or `9` (3D) + `equalDOF` tying to `ndm` soil | **`ndm` (`2` or `3`) in single-phase; `ndm + 1` (`3` or `4`) or Taylor–Hood in `-up` mode** (zero auxiliary DOFs, zero `equalDOF` constraints) |
| **Supported Topologies** | `Q4` (2D), `H8` (3D) only | `Q4` (2D), `H8` (3D) only | **All 9 continuum topologies**: `T3, Q4, T6, BT6` (2D) and `TET4, H8, TET10, BTET10, H20` (3D) |
| **Layer Thickness** | **1 element thick** (local Lysmer–Kuhlemeyer dashpot + free-field skin) | Multi-element layer ($5\text{–}10$ elements) | **Multi-element layer** ($N_{\text{pml}} = 4\text{–}10$ elements across thickness $L_{\text{pml}}$) |
| **Boundary Orientation** | Axis-aligned (`L, R, F, K, B`) only; requires per-element `btype` string | Axis-aligned `-box` only | **Axis-aligned `-box`** (one shared `-box` definition automatically classifies faces, 2D corners, 3D edges, and 3D corners at each Gauss point!) **AND arbitrary oblique `-plane`** (`-x0, -n, -L`) |
| **Coupled Biot $u\text{–}p$ Support** | None (single-phase total stress only) | None | **Native `-up` mode** (`u-p` equal-order with McGann $\tilde{\mathbf{H}}$ stabilization or Taylor–Hood `pOrder="linear"` on `T6, BT6, TET10, BTET10, H20`; `std` and `bbar` formulations), bit-for-bit impedance-matched to `LadrunoUP` (`33017`) at $\xi = 0$ |
| **Outer-Wall Boundary Condition** | Outer ghost nodes left free (`B` base / free-field) | Outer PML nodes fixed (`u = 0`) | **Outer PML boundary nodes FIXED** (`fix ... 1 1 [1]` in single-phase; `fix ... 1 1 [1] 1` in `-up` mode) |
| **Rayleigh Damping** | Required on ASD elements (drives free-field column damping) | Optional | Optional (`alphaM, betaK, betaK0, betaKc`); in `-up` mode, stiffness-proportional Rayleigh damping is **automatically isolated strictly to the $(u, u)$ solid block** |
| **Element Assembly Cost** | $\sim 0.43\,\mu\text{s}$ (`Q4`), $\sim 7.0\,\mu\text{s}$ (`H8`) per element-step | $\sim 140\,\mu\text{s}$ (`H8` 72-DOF Newmark) | $\sim 2.42\,\mu\text{s}$ (`Q4`), $\sim 11.37\,\mu\text{s}$ (`H8`), $\sim 40.2\,\mu\text{s}$ (`H20`) per element-step (post-PR #996 hoisted GP geometry + diagonal `D_F` fast path) |

---

## Part 2 — Geometry & Mesh Generation in `apeGmsh` (`add_cpml_shell` / `add_plane_wave_box`)

### 2.1 Multi-Element Offset Shell vs. 1-Element Ghost Skin
Whereas `ASDAbsorbingBoundary` uses a 1-element-thick morphological dilation (`n_layers = 1`), `LadrunoCPML` attenuates outgoing body, evanescent, and Rayleigh/Lamb surface waves across a finite physical thickness $L_{\text{pml}} = (L_x, L_y[, L_z])$ discretized into $N_{\text{pml}} \ge 4$ elements (recommended default: `n_pml_layers = 5` or `6`).

Because `LadrunoCPML` in `-box` mode computes the local penetration depth
$$\xi_m(\mathbf{x}_g) = \max\bigl(0,\, x_{\min, m} - x_{g, m},\, x_{g, m} - x_{\max, m}\bigr)$$
directly at **every Gauss point $\mathbf{x}_g$**, **`apeGmsh` does NOT need to partition the PML shell into 5/9/17 separate `btype` physical groups (`L`, `R`, `B`, `BL`, `BLF`, ...)!**

- **Single Homogeneous Soil Layer**: All PML cells wrapping the 3 (2D: $-X, +X, -Y$) or 5 (3D: $-X, +X, -Y, +Y, -Z$) truncation faces share **one single physical group** `"CPML"` (plus one outer boundary surface/curve physical group `"CPML_Outer"` for Dirichlet fixation)!
- **Vertically Layered Soil Profile (`Axis1D` z-layers)**: Each soil layer $k$ (with material $(E_k, \nu_k, \rho_k)$ and optional Biot parameters) gets one PML physical group `"CPML_layer_{k}"` at the same elevation $z \in [z_{k,\min}, z_{k,\max}]$ (with the bottom $-Z$ PML layer inheriting the bottom-most soil layer's properties), all sharing the **same global inner-domain bounding box** `[xMin, xMax, yMin, yMax, zMin, zMax]`!

### 2.2 Free Surface (`+Y` in 2D / `+Z` in 3D) Rule
For half-space geotechnical and SSI problems where the top surface ($y = y_{\max}$ in 2D or $z = z_{\max}$ in 3D) is a free ground surface (so Rayleigh/Lamb waves propagate along the top boundary into the lateral PML walls):
- **Do NOT place a PML layer above the free surface**, AND
- **Set `yMax` (2D) or `zMax` (3D) in the `-box` argument above the top of the mesh** (e.g. `zMax_box = z_top + 10.0 * Lz`), so that $\xi_z(\mathbf{x}_g) = 0$ at the free surface $z = z_{\text{top}}$ while $\xi_z(\mathbf{x}_g) = z_{\min} - z_g > 0$ in the bottom PML layer ($z < z_{\min}$)!
- Conversely, for full-space / buried waveguide problems wrapped on all 4 (2D) or 6 (3D) sides, set `[xMin, xMax, yMin, yMax[, zMin, zMax]]` to the exact inner soil box bounds.

### 2.3 Outer Boundary Fixation (`CPML_Outer`)
The outer termination surfaces of the PML shell (`x == xMin - Lx`, `x == xMax + Lx`, `y == yMin - Ly`, `[y == yMax + Ly]`, `z == zMin - Lz`) MUST be tagged into a boundary physical group `"CPML_Outer"`:
- **Single-phase mode (`ndf = ndm`)**: fix all displacement DOFs on `"CPML_Outer"` (`fix $nodeTag 1 1` in 2D, `fix $nodeTag 1 1 1` in 3D).
- **Coupled Biot `-up` mode**:
  - Equal-order (`pOrder="equal"`, all nodes `ndf = ndm + 1`): fix both solid displacements and pore pressure (`fix $nodeTag 1 1 1` in 2D, `fix $nodeTag 1 1 1 1` in 3D) on `"CPML_Outer"`.
  - Taylor–Hood (`pOrder="linear"` on `T6, BT6, TET10, BTET10, H20`): vertex nodes on `"CPML_Outer"` have `ndf = ndm + 1` (`fix $vNode 1 1 [1] 1`), while mid-edge nodes on `"CPML_Outer"` have `ndf = ndm` (`fix $mNode 1 1 [1]`). `apeGmsh`'s per-node resolved `ndf` mask machinery (ADR 0074) handles this automatically when `fix` masks validate against resolved per-node `ndf`.

---

## Part 3 — Emitter Side (`apeGmsh.opensees.element.LadrunoCPML`)

### 3.1 Canonical Tcl & Python (`openseespy`) Syntax

#### Single-Phase Continuum Mode (`ndf = ndm`)
```tcl
# Axis-aligned -box mode (2D: 6 box args; 3D: 9 box args)
element LadrunoCPML $eleTag $n1 ... $nk \
    [-shape T6|BT6|TET10|BTET10] \
    {-mat $matTag | -E $E -nu $nu -rho $rho} \
    [-thick $t] \
    -box $xMin $xMax $yMin $yMax [$zMin $zMax] $Lx $Ly [$Lz] \
    [-m $mDeg] [-R0 $R0] [-kappaMax $kappaMax] [-alphaMax $alphaMax] \
    [-mass consistent|lumped] [-stage 0|1] [-body $bx $by [$bz]]

# Arbitrary oblique -plane mode
element LadrunoCPML $eleTag $n1 ... $nk \
    [-shape T6|BT6|TET10|BTET10] \
    {-mat $matTag | -E $E -nu $nu -rho $rho} \
    [-thick $t] \
    -plane -x0 $x0 $y0 [$z0] -n $nx $ny [$nz] -L $Lpml \
    [-m $mDeg] [-R0 $R0] [-kappaMax $kappaMax] [-alphaMax $alphaMax] \
    [-mass consistent|lumped] [-stage 0|1] [-body $bx $by [$bz]]
```

#### Coupled Biot $u\text{–}p$ Poroelastic Mode (`-up`)
```tcl
element LadrunoCPML $eleTag $n1 ... $nk \
    [-shape T6|BT6|TET10|BTET10] \
    {-mat $matTag | -E $E -nu $nu -rho $rho} \
    [-thick $t] \
    {-box $xMin $xMax $yMin $yMax [$zMin $zMax] $Lx $Ly [$Lz] | \
     -plane -x0 $x0 $y0 [$z0] -n $nx $ny [$nz] -L $Lpml} \
    -up -Kf $Kf -poro $poro -rhoF $rhoF \
    {-perm $kx [$ky [$kz]] | -permH $khx [$khy [$khz]] -gammaW $gw} \
    [-alphaB $alphaB] [-Ks $Ks] \
    [-pOrder equal|linear] [-formulation std|bbar] \
    [-stab auto [$alpha0] | -stab off | -stab manual $alpha] \
    [-fluidBody $fbx $fby [$fbz]] \
    [-m $mDeg] [-R0 $R0] [-kappaMax $kappaMax] [-alphaMax $alphaMax] \
    [-mass consistent|lumped] [-stage 0|1] [-body $bx $by [$bz]]
```

### 3.2 Recommended `apeGmsh` Typed Primitive Signature

In `src/apeGmsh/opensees/element/continuum.py` (or `absorbing.py`), expose `ops.element.LadrunoCPML`:

```python
ops.element.LadrunoCPML(
    pg: str,
    material: NDMaterial | None = None,
    *,
    E: float | None = None,
    nu: float | None = None,
    rho: float | None = None,
    thickness: float = 1.0,                 # 2D only (elided if 1.0 or 3D)
    shape: Literal["T6", "BT6", "TET10", "BTET10"] | None = None,
    # Geometry mode (mutually exclusive: supply box OR plane)
    box: tuple[float, ...] | None = None,   # 2D: (xMin, xMax, yMin, yMax, Lx, Ly)
                                            # 3D: (xMin, xMax, yMin, yMax, zMin, zMax, Lx, Ly, Lz)
    plane_x0: tuple[float, ...] | None = None,
    plane_normal: tuple[float, ...] | None = None,
    plane_length: float | None = None,
    # CFS-PML grading parameters (None -> elide; C++ defaults: m=2.0, R0=1e-6, kappaMax=5.0, alphaMax=pi)
    m_deg: float | None = None,
    R0: float | None = None,
    kappa_max: float | None = None,
    alpha_max: float | None = None,
    lumped: bool = False,
    stage: Literal[0, 1] | None = None,     # 0 for gravity/seepage stage, 1 for dynamic (default 1)
    body_force: tuple[float, ...] | None = None,
    # Coupled Biot u-p mode (-up)
    up: bool = False,
    Kf: float | None = None,
    poro: float | None = None,
    rhoF: float | None = None,
    perm: float | tuple[float, ...] | None = None,
    permH: float | tuple[float, ...] | None = None,
    gammaW: float | None = None,
    alphaB: float | None = None,
    Ks: float | None = None,
    p_order: Literal["equal", "linear"] | None = None,
    formulation: Literal["std", "bbar"] | None = None,
    stab: Literal["auto", "off"] | float | tuple[str, float] | None = None,
    fluid_body: tuple[float, ...] | None = None,
)
```

### 3.3 Topology Fan-Out & Node Ordering Across All 9 Shapes

`LadrunoCPML` accepts Gmsh connectivity **verbatim** across all 9 topologies (`node_reorder = identity`):

| `(ndm, nNodes)` | Gmsh `etype` | Topology | `-shape` token needed? | Taylor–Hood `ndf_floor_per_slot` (`up=True, p_order="linear"`) |
|---|---|---|---|---|
| `(2, 3)` | `2` (`tri3`) | `T3` | No | Equal-order only (`(3, 3, 3)`) |
| `(2, 4)` | `3` (`quad4`) | `Q4` | No | Equal-order only (`(3, 3, 3, 3)`) |
| `(2, 6)` | `9` (`tri6`) | `T6` (Lagrange) or `BT6` (Bézier) | Elide for `T6`; emit `-shape BT6` if pairing with `BezierTri6` / `LadrunoUP` | `(3, 3, 3, 2, 2, 2)` (3 vertices `ndf=3`, 3 mid-edges `ndf=2`) |
| `(3, 4)` | `4` (`tet4`) | `TET4` | No | Equal-order only (`(4, 4, 4, 4)`) |
| `(3, 8)` | `5` (`hex8`) | `H8` | No | Equal-order only (`(4,)*8`) |
| `(3, 10)` | `11` (`tet10`) | `TET10` (Lagrange) or `BTET10` (Bézier) | Elide for `TET10`; emit `-shape BTET10` if pairing with `BezierTet10` / `LadrunoUP` | `(4,)*4 + (3,)*6` (4 vertices `ndf=4`, 6 mid-edges `ndf=3`) |
| `(3, 20)` | `17` (`hex20`) | `H20` | No | `(4,)*8 + (3,)*12` (8 corners `ndf=4`, 12 mid-edges `ndf=3`) |

> **Important pairing rule (`BT6` / `BTET10` vs. `T6` / `TET10`):**
> - When the interior domain uses `LadrunoLST` (Lagrange `T6`) or standard Lagrange `TET10`, omit `-shape` (defaults to `T6` / `TET10`).
> - When the interior domain uses `BezierTri6`, `BezierTet10`, or `LadrunoUP` (`BT6` / `BTET10`), emit `-shape BT6` (2D) or `-shape BTET10` (3D) so the shared interface faces use the exact same Bernstein–Bézier shape functions and achieve bit-identical impedance matching at $\xi = 0$!

---

## Part 4 — Five Critical Emitter Rules & Build Gates

### Rule 1 — Body-Force Unit Convention Split (`single-phase` vs. `-up` mode)
To preserve 100% compatibility with both `LadrunoBrick` (single-phase) and `LadrunoUP` (coupled Biot $u\text{–}p$):
- **Single-Phase Mode (`up=False`)**: `-body $bx $by [$bz]` is a **force per unit volume** $[F/L^3]$ (i.e. $\rho \mathbf{g}$), identical to `LadrunoBrick` and `LadrunoQuad`.
- **Coupled Biot `-up` Mode (`up=True`)**: `-body $bx $by [$bz]` and `-fluidBody $fbx $fby [$fbz]` are **accelerations** $[L/T^2]$ (i.e. $\mathbf{g} = (0, 0, -9.81)$), identical to `LadrunoUP`! Inside `LadrunoCPML`, the solid rows multiply `-body` by $\rho$ and the pressure rows multiply `-fluidBody` by $\rho_f \bar{k}_i$.
- **Why this matters for `apeGmsh`**: An `apeGmsh` model pairing `LadrunoUP` in the interior with `LadrunoCPML(up=True)` in the absorbing shell passes the **exact same acceleration tuple** `body_force=(0, -9.81)` to both elements with zero factor-of-$\rho$ discrepancy at the interface!

### Rule 2 — General Unsymmetric Solver Gate in `-up` Mode
In `-up` mode, `LadrunoCPML` uses the exact same **honest-$p$ block layout** as `LadrunoUP`:
$$\mathbf{K}_{\text{eff}} = \begin{bmatrix} \mathbf{K}_{\text{eff}, uu} & -\mathbf{Q}_{\text{eff}} \\ \mathbf{0} & \mathbf{H}_{\text{eff}} \end{bmatrix}, \qquad \mathbf{C}_{\text{eff}} = \begin{bmatrix} \mathbf{C}_{\text{eff}, uu} + \mathbf{C}_{\text{Ray}, uu} & \mathbf{0} \\ +\mathbf{Q}_{\text{eff}}^T & \mathbf{S}_{\text{eff}} + \tilde{\mathbf{H}} \end{bmatrix}$$
Because $\mathbf{K}_T^{\text{dyn}} = c_1 \mathbf{K}_{\text{eff}} + c_2 \mathbf{C}_{\text{eff}} + c_3 \mathbf{M}_{\text{eff}}$ is **unsymmetric**, `apeGmsh` MUST apply the exact same **general-solver build gate** to `LadrunoCPML(up=True)` that it applies to `LadrunoUP` (ADR 0074):
- Allow-list: `UmfPack`, `SparseGEN` / `SparseGeneral`, `FullGeneral`, `BandGeneral`, `Mumps` (with unsymmetric `ICNTL(14)` / `SYM=0`), or `Pardiso` (`-matrixType 11`).
- Refuse symmetric solvers (`ProfileSPD`, `BandSPD`, `SymBandEigen`, etc.) at deck-build time.
- In **single-phase mode (`up=False`)**, $\mathbf{K}_{\text{eff}}, \mathbf{C}_{\text{eff}}, \mathbf{M}_{\text{eff}}$ are **strictly symmetric positive semi-definite**, so symmetric and SPD solvers (`ProfileSPD`, `Pardiso -matrixType 2`, `UmfPack`, `Mumps`) are all legal!

### Rule 3 — Gravity & Steady-Seepage Staged Activation (`-stage 0` $\to$ `stage 1`)
When a model has a static gravity / initial stress / steady-state Darcy seepage stage followed by a transient wave-propagation stage:
1. Emit `LadrunoCPML` with `-stage 0`.
   - At `stage 0`, all C-PML auxiliary filter states are frozen, $\mathbf{M}_{\text{eff}} = \mathbf{M}_0$, $\mathbf{K}_{\text{eff}} = \mathbf{K}_0 = \begin{bmatrix} \mathbf{K}_{0,uu} & -\mathbf{Q}_0 \\ \mathbf{0} & \mathbf{H}_0 \end{bmatrix}$, and `getResistingForce()` returns the unstretched static elastic + pore-pressure + Darcy seepage balance.
2. After `StaticAnalysis` (`analyze(1)`) converges and before starting `TransientAnalysis`, `s.activate_absorbing()` (or a new `s.activate_cpml()`) emits:
   ```tcl
   setParameter -val 1 -ele $firstCpmlTag ... $lastCpmlTag stage
   ```
   (or `parameter $pTag element $eleTag stage; updateParameter $pTag 1`).
   - On the `0 -> 1` transition, `LadrunoCPML` automatically snapshots the committed static displacement and pore-pressure state $\mathbf{u}_0 = [\mathbf{u}_{s,0}, \mathbf{p}_0]^T$ and static internal force $\mathbf{f}_0 = \mathbf{K}_0 \mathbf{u}_0$, zeroes all recursive-convolution filter states, and thereafter feeds only the **dynamic perturbation** $(\mathbf{u} - \mathbf{u}_0, \mathbf{p} - \mathbf{p}_0)$ into the C-PML filters (`f_total = f0 + f_cpml(u - u0)`). This guarantees **zero unbalance jump** at the first transient step!

### Rule 4 — Material Specification (`-mat $matTag` vs. `-E -nu -rho`)
`LadrunoCPML` accepts either:
- `-mat $matTag`: extracts $(E, \nu, \rho)$ from any linear isotropic `NDMaterial` (`ElasticIsotropic`) in the domain at creation time (and optional `-rho $rho` overrides the material density if given), OR
- `-E $E -nu $nu -rho $rho`: inline elastic parameters without needing an `nDMaterial` tag.
In `apeGmsh`, allowing `material=soil_mat` (which emits `-mat $matTag`) lets the user pass the exact same `ElasticIsotropic` material object used by `LadrunoQuad` / `LadrunoBrick` / `LadrunoUP`.

### Rule 5 — Fork Capability Gate
Add `"LadrunoCPML"` (and `"ladrunoCPML"`) to `_FORK_ONLY_ELEMENTS` in `src/apeGmsh/opensees/emitter/live.py`, and add `ELE_TAG_LadrunoCPML = 33023` to `_response_catalog.py`.

---

## Part 5 — Postprocessor & Recorder Integration

### 5.1 Supported `eleResponse` Queries on `LadrunoCPML`
| Token | Output Size | Contents |
|---|---|---|
| `"force"` / `"forces"` / `"globalForce"` / `"globalForces"` | `nDOF` | Total resisting force vector (including C-PML damping/convolution + inertia + Rayleigh forces) |
| `"stiff"` / `"stiffness"` / `"tangent"` | `nDOF * nDOF` | Row-major flattened $\mathbf{K}_{\text{eff}}(\Delta t)$ |
| `"damp"` / `"damping"` | `nDOF * nDOF` | Row-major flattened $\mathbf{C}_{\text{eff}}(\Delta t) + \mathbf{C}_{\text{Ray}}$ |
| `"mass"` | `nDOF * nDOF` | Row-major flattened $\mathbf{M}_{\text{eff}}(\Delta t)$ |
| `"dampingForce"` / `"dampingForces"` | `nDOF` | Solid-only Rayleigh damping force vector |
| `"basicForce"` / `"basicForces"` | `nDOF` | Kernel internal resisting force $\text{scale} \cdot \mathbf{R}_{\text{int}}$ (displacement + velocity convolution + kernel mass rate, without nodal $M\ddot{u}$ or Rayleigh) |
| `"energy"` / `"energies"` | `2` | `[W_ext, P_ext]` where $W_{\text{ext}} = \int \mathbf{v}^T \mathbf{R}_{\text{int}}\,dt$ (trapezoidal cumulative work absorbed by the element) and $P_{\text{ext}} = \mathbf{v}^T \mathbf{R}_{\text{int}}$ (instantaneous power) |
| `"stage"` | `1` | Current stage (`0.0` = static, `1.0` = dynamic C-PML active) |

### 5.2 `.ladruno` Recorder & `EnergyBalanceRecorder` Integration
- **`EnergyBalanceRecorder`**: `LadrunoCPML` implements the full 6-accessor contract (`getTangentStiff`, `getInitialStiff`, `getDamp`, `getMass`, `getResistingForce`, `getResistingForceIncInertia`), so `recorder EnergyBalance` automatically integrates the energy absorbed by the C-PML layer in `DW_ele + IE` with no residual leak!
- **3D Viewport / Contour Filtering**: Just as `apeGmsh` excludes `ASDAbsorbingBoundary` ghost elements from default stress/displacement contour auto-scaling so that attenuated boundary waves do not skew the color bar of the physical soil domain, tag `LadrunoCPML` physical groups as absorbing layers so the viewer can toggle their visibility on/off independently of the interior soil domain.
