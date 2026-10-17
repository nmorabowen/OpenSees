# `LadrunoMembrane` — User & Developer Guide (ADR 98)

`LadrunoMembrane` (`ELE_TAG_LadrunoMembrane = 33022`) is an arbitrary $N$-node polygonal floor/roof load-distribution and optional in-plane CST diaphragm element in the Ladruno fork of OpenSees.

It solves two coupled problems in structural frame modeling:
1. **Surface-to-Frame Gravity/Pressure Load Distribution:** Converts uniform surface pressures $q$ acting on arbitrary planar polygons (including non-convex L-shapes, oblique edges, interior `-hole` openings, partial/interior supporting `-frames`, and isolated column/wall `-supportNodes`) into moment-conserving partial trapezoidal `Beam2dPartialUniformLoad` / `Beam3dPartialUniformLoad` instances on supporting beams and direct nodal forces on isolated point supports.
2. **Optional Guyan-Condensed In-Plane Diaphragm:** Optionally provides in-plane membrane stiffness and consistent/lumped mass by triangulating the polygon with a Constrained Delaunay Triangulation (CDT) and statically condensing internal Steiner vertices onto the polygon's external boundary nodes.

---

## 1. Command Syntax

### Tcl

```tcl
element LadrunoMembrane $eleTag \
    -nodes $n1 $n2 $n3 ... \
    [-frames $f1 $f2 ...] \
    [-supportNodes $sn1 $sn2 ...] \
    [-hole $hn1 $hn2 $hn3 ...] \
    [-method <tributary|oneway|fea_poisson|fea_plate>] \
    [-spanAxis $dx $dy [$dz]] \
    [-orthoRatio $eta] \
    [-pressure $q] \
    [-thick $t -mat $ndMatTag] \
    [-stiff $E $nu $t] \
    [-rho $rho] \
    [-lumped] \
    [-meshSize $hMax]
```

### OpenSeesPy

```python
ops.element(
    "LadrunoMembrane", eleTag,
    "-nodes", *boundaryNodeTags,
    "-frames", *frameEleTags,                # optional (required unless -supportNodes given)
    "-supportNodes", *supportNodeTags,       # optional point supports (columns/walls)
    "-hole", *hole1NodeTags,                 # repeatable for multiple openings
    "-method", "tributary",                  # "tributary" | "oneway" | "fea_poisson" | "fea_plate"
    "-spanAxis", dx, dy, dz,                 # span direction for oneway or orthotropic fea_poisson
    "-orthoRatio", eta,                      # E2/E1 ratio in (0, 1] for fea_poisson (default 1.0)
    "-pressure", q,                          # baseline normal surface pressure (acts along +n)
    "-thick", t, "-mat", ndMatTag,           # optional in-plane CST diaphragm via PlaneStress NDMaterial
    "-stiff", E, nu, t,                      # optional in-plane CST diaphragm via isotropic E, nu, t
    "-rho", rho,                             # optional mass density (mass per unit volume; areal = rho*t)
    "-lumped",                               # optional row-sum lumped mass instead of Guyan consistent mass
    "-meshSize", hMax,                       # optional target CDT mesh size (default: 0.12 * L_char)
)
```

---

## 2. Arguments and Flags Reference

| Flag | Arguments | Default | Description |
|---|---|---|---|
| `-nodes` | `$n1 $n2 $n3 ...` ($\ge 3$) | *(Required)* | Ordered outer boundary node tags of the planar polygon. Must be coplanar within $10^{-3} L_{\text{char}}$ and have non-zero area. |
| `-frames` | `$f1 $f2 ...` | `[]` | Supporting 2D or 3D beam-column element tags (`ElasticBeam2d/3d`, `ModElasticBeam3d`, `DispBeamColumn2d/3d`, `LadrunoDispBeamColumn2d/3d`, `ForceBeamColumn2d/3d/Thermal/CBDI3d`, `ElasticTimoshenkoBeam2d/3d`). Frames may lie along the polygon perimeter, partially along an edge, or inside the polygon interior. Resolved lazily on first use, so `LadrunoMembrane` may be declared before or after its supporting frame elements. |
| `-supportNodes` | `$sn1 $sn2 ...` | `[]` | Isolated point-support node tags (e.g., flat-slab columns or wall piers not connected by beams). Receives direct nodal gravity reactions during `addLoad`/`getResistingForce`. May be either boundary nodes or interior nodes. |
| `-hole` | `$h1 $h2 $h3 ...` ($\ge 3$) | None | Ordered node tags defining an interior opening/cutout. May be specified multiple times for multiple openings. Hole boundaries have natural free-edge ($\mathbf{n}\cdot\mathbf{D}\nabla w = 0$) conditions and carry zero surface load. |
| `-method` | `tributary` \| `oneway` \| `fea_poisson` \| `fea_plate` | `tributary` | Surface-to-frame load distribution engine (see Section 3). |
| `-spanAxis` | `$dx $dy [$dz]` | `(1, 0, 0)` | Principal span direction vector projected onto the slab plane. Used by `oneway` (strip ray direction) and `fea_poisson` (principal axis $\hat{\mathbf{e}}_1$ when `-orthoRatio < 1`). |
| `-orthoRatio` | `$eta` $\in (0, 1]$ | `1.0` | Orthotropy ratio $\eta = D_2 / D_1$ for `fea_poisson`. As $\eta \to 10^{-4}$, `fea_poisson` smoothly approaches the one-way strip distribution along `-spanAxis`. |
| `-pressure` | `$q` | `0.0` | Reference normal surface pressure (force per unit area; acts along the polygon's Newell normal $\hat{\mathbf{n}}$). Scaled by the active `LoadPattern` factor when the membrane is included in an `eleLoad` or pattern. |
| `-thick` + `-mat` | `$t` (`> 0`), `$matTag` | `0.0`, `0` | Enables the Guyan-condensed in-plane CST diaphragm using an existing `NDMaterial` (cloned via `"PlaneStress"` or `"ThreeDimensional"` plane-stress condensation). |
| `-stiff` | `$E $nu $t` | Disabled | Direct isotropic plane-stress in-plane CST diaphragm ($E > 0$, $-1 < \nu < 0.5$, $t > 0$) without requiring a separate `nDMaterial` command. |
| `-rho` | `$rho` ($\ge 0$) | `0.0` | Volumetric mass density (so areal mass density is $\rho_{\text{area}} = \rho \cdot t$, or $\rho$ directly if $t = 0$). |
| `-lumped` | *(flag)* | `false` | Uses row-sum lumped nodal mass on the polygon boundary nodes instead of Guyan-condensed consistent CST mass. |
| `-meshSize` | `$hMax` (`> 0`) | $0.12 L_{\text{char}}$ | Maximum triangle circumradius / edge target for Chew-Ruppert CDT refinement (clamped to $[10^{-4}\sqrt{A},\; 10 L_{\text{char}}]$ with a 4000-vertex safety cap). |

---

## 3. Load-Distribution Engines (`-method`)

All four methods construct a deterministic 2D Planar Straight-Line Graph (PSLG) in the slab's local orthonormal basis $(\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2, \hat{\mathbf{n}})$, split intersecting or T-junction segments, tag acute corners ($< 60^\circ$), and refine a Constrained Delaunay Triangulation (CDT).

### 3.1 `tributary` (Default — CDT-Based Euclidean Voronoi / Yield-Line Geometric Partition)
Each point $\mathbf{x} \in \Omega$ is mapped to the nearest supporting feature $S^*$ (frame sub-segment or `-supportNodes` point) under exact 2D Euclidean distance $d(\mathbf{x}, S)$.
- CDT triangles whose 7 sample points map to the same interior frame segment are integrated in closed form via 3-point triangle quadrature.
- CDT triangles straddling a Voronoi ridge between two frame segments are clipped analytically along the bisector $d_1(\mathbf{x}) - d_2(\mathbf{x}) = 0$.
- Multi-feature apex triangles are recursively subdivided (depth 3).
- **Properties:** On rectangular two-way bays, reproduces the classical $45^\circ$ yield-line triangles ($q L_s^2 / 4$) and trapezoids ($q L_s (L_l - L_s / 2) / 2$) to $< 10^{-4}$ relative error, while remaining strictly non-negative ($w \ge 0$), leak-free around interior `-hole` openings, and well-defined on non-convex L-shapes and partial edges.

### 3.2 `oneway` (Bidirectional Ray-Cast Strip Partitioning)
Projects quadrature points $\mathbf{x}_q$ along $\pm\hat{\mathbf{s}}$ (`-spanAxis`) to find the nearest supporting intersections at distances $d_L$ and $d_R$:
- **Two-sided span:** $f_L = f_q \frac{d_R}{d_L + d_R}$, $f_R = f_q \frac{d_L}{d_L + d_R}$ (exact lever-arm equilibrium conserving $\sum F_z$, $\sum M_x$, and $\sum M_y$ pointwise).
- **Single-sided cantilever overhang:** 100% of $f_q$ is assigned to the single intersected support (contributing vertical force; cantilever root bending moment requires `fea_plate`).

### 3.3 `fea_poisson` (1-DOF CST Membrane Deflection Analogy)
Solves $-\nabla \cdot (\mathbf{D}\nabla w) = q$ on the CDT mesh with homogeneous Dirichlet boundary conditions $w = 0$ on all supporting frame sub-segments and `-supportNodes`, and natural Neumann free edges $\mathbf{n}\cdot\mathbf{D}\nabla w = 0$ on unsupported perimeter edges and `-hole` boundaries:
- Consistent reaction residual $\mathbf{R} = \mathbf{K}\mathbf{w} - \mathbf{f}_{\text{ext}}$.
- Reactions at shared vertices (corners/T-junctions) are never split 50/50; they are apportioned across incident constrained sub-segments proportionally to each adjacent triangle's integrated outward normal diffusive flux $L_{S_m}\max(0, -\hat{\mathbf{n}}_{S_m}\cdot\mathbf{D}\nabla w_e)$, and mapped to linear trapezoids via the inverse 1D consistent mass matrix $\mathbf{M}_{1D}^{-1}$.

### 3.4 `fea_plate` (3-DOF Discrete Kirchhoff Triangle Plate Bending)
Assembles the 9-DOF Batoz-Bathe-Ho (1980) DKT plate bending stiffness $\mathbf{K}_{\text{DKT}}$ for $(w, \theta_x, \theta_y)$ at every CDT vertex:
- **Non-collinear supports ($\operatorname{rank} \ge 2$, $\sigma_2 \ge 10^{-6}L_{\text{char}}$):** Supported nodes act as knife-edge simple supports ($w = 0$, rotations $\theta_x, \theta_y$ free). Signed reactions (including corner hold-down uplift forces) are preserved so $\sum F_z = qA$, $\sum M_x$, and $\sum M_y$ hold to machine precision.
- **Collinear or single-beam supports ($\operatorname{rank} < 2$):** Rotations $\theta_x = \theta_y = 0$ are clamped at supported nodes to prevent rigid-body rotation about the support line, enabling pure cantilever slabs on a single girder.

---

## 4. Multi-Pattern Loading and `eleLoad` Integration

`LadrunoMembrane` supports applying surface pressures across any number of `LoadPattern`s:

1. **Via `-pressure $q` on the element + `eleLoad -ele $eleTag -type -surfaceLoad` (or `-selfWeight $gx $gy $gz`):**
   ```python
   ops.element("LadrunoMembrane", 100, "-nodes", 1, 2, 3, 4,
               "-frames", 1, 2, 3, 4, "-pressure", -5.0)
   ops.timeSeries("Linear", 1)
   ops.pattern("Plain", 1, 1)
   ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")
   ```
2. **Multi-pattern scaling (`Dead` + `loadConst` + `Live`):**
   When multiple `LoadPattern`s act concurrently (e.g. Dead load followed by `loadConst("-time", 0.0)` and then Live load), `LadrunoMembrane::zeroLoad()` clears its internal `m_activeBeamLoads` arena and each `addLoad(theLoad, loadFactor)` call allocates fresh `Beam3dPartialUniformLoad` / `Beam2dPartialUniformLoad` instances pre-scaled by `loadFactor`. Supporting `ElasticBeam3d` and `ForceBeamColumn3d` elements never double-count across patterns.

---

## 5. Recorder and `eleResponse` Queries

Query any `LadrunoMembrane` via `ops.eleResponse(eleTag, query)` or Tcl `eleResponse $eleTag $query`:

| Query | Return Format | Description |
|---|---|---|
| `"area"` | `[A]` (length 1) | Net planar surface area $A_{\text{net}} = A_{\text{outer}} - \sum A_{\text{holes}}$. |
| `"totalLoad"` | `[Fx, Fy, Fz]` (length 3) | Total integrated global surface force vector $\sum \mathbf{F}_{\text{applied}}$ currently applied across active load patterns. |
| `"tributaryAreas"` (alias `"tributaryWidths"`) | `[A_1, ..., A_K, A_s1, ...]` (length $N_{\text{frames}} + N_{\text{supportNodes}}$) | Effective tributary area $A_k = \int_{\text{frame } k} w_k^{\text{unit}}(s)\,ds$ carried by each supporting frame in `-frames` order, followed by the tributary area of each isolated `-supportNodes` entry. |
| `"frameLoads"` | `[aL_1, bL_1, wa_1, wb_1, ...]` ($4 \times N_{\text{frames}}$) | Equivalent moment-conserving partial trapezoidal load parameters per supporting frame in `-frames` order: normalized start `aL` ($a/L$), normalized end `bL` ($b/L$), start intensity `wa`, and end intensity `wb`. |
| `"subSegmentLoads"` (alias `"beamSubLoads"`) | `[frameTag, aL, bL, wya, wyb, wza, wzb, wxa, wxb, ...]` ($9 \times N_{\text{sub}}$) | Exact multi-segment partial trapezoidal loads emitted onto supporting frames, transformed into each target frame's local coordinate axes $(\hat{x}_b, \hat{y}_b, \hat{z}_b)$ and scaled by the active global traction (or unit normal pressure before `addLoad`). Used by `apeGmsh` to reconstruct exact piecewise-quadratic shear $V(x)$ and piecewise-cubic moment $M(x)$ diagrams (see [`98_ladruno_membrane_apegmsh_emitter_guide.md`](98_ladruno_membrane_apegmsh_emitter_guide.md)). |
| `"subSegments"` | `[N_sub]` (length 1) | Number of emitted partial trapezoidal beam sub-segments. |
| `"overhangTorques"` | `[T_1, ..., T_K]` (length $N_{\text{frames}}$) | Total cantilever overhang torque about each frame's local $\hat{x}_b$ axis when `supportRank < 2` (`fea_plate` cantilever). |
| `"supportReactions"` (alias `"nodalReactions"`) | `[R_1, R_2, ..., R_S]` (length $N_{\text{supportNodes}}$) | Direct normal force reactions transferred to each isolated `-supportNodes` entry. |
| `"localAxes"` | `[e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz]` (length 9) | Flattened $3 \times 3$ orthonormal local basis $(\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2, \hat{\mathbf{n}})$ of the membrane plane. |
| `"meshStats"` | `[numVertices, numTriangles, numSubSegments, supportRank]` (length 4) | Internal 2D Chew-Ruppert CDT mesh telemetry and plate support rank. |
| `"meshVertices"` | `[X_0, Y_0, Z_0, ...]` ($3 \times N_{\text{verts}}$) | 3D global coordinates of every internal CDT mesh vertex. |
| `"meshTriangles"` | `[v0_0, v1_0, v2_0, ...]` ($3 \times N_{\text{tris}}$) | 0-based triangle vertex indices into `"meshVertices"`. |

---

## 6. Examples

### 6.1 Two-Way Rectangular Floor Bay (`tributary` vs `fea_poisson`)

```python
from _testbed import ops

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)

Lx, Ly, q = 6.0, 4.0, -10.0  # 6 m x 4 m slab, 10 kN/m^2 downward (-Z) pressure
ops.node(1, 0.0, 0.0, 0.0)
ops.node(2, Lx,  0.0, 0.0)
ops.node(3, Lx,  Ly,  0.0)
ops.node(4, 0.0, Ly,  0.0)
for n in (1, 2, 3, 4):
    ops.fix(n, 1, 1, 1, 1, 1, 1)

ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
E, G, A, J, Iy, Iz = 30e6, 12.5e6, 0.18, 0.002, 0.003, 0.003
ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, J, Iy, Iz, 1)
ops.element("elasticBeamColumn", 2, 2, 3, A, E, G, J, Iy, Iz, 1)
ops.element("elasticBeamColumn", 3, 3, 4, A, E, G, J, Iy, Iz, 1)
ops.element("elasticBeamColumn", 4, 4, 1, A, E, G, J, Iy, Iz, 1)

ops.element(
    "LadrunoMembrane", 100,
    "-nodes", 1, 2, 3, 4,
    "-frames", 1, 2, 3, 4,
    "-method", "tributary",
    "-pressure", q,
)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")

ops.system("BandGeneral")
ops.numberer("Plain")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)
ops.reactions()

# Total upward support reaction equals -q * Lx * Ly = +240.0 kN
Rz_total = sum(ops.nodeReaction(n, 3) for n in (1, 2, 3, 4))
assert abs(Rz_total - (-q * Lx * Ly)) < 1e-9

# Short beams (2, 4) carry Ly^2 / 4 = 4.0 m^2 each; long beams (1, 3) carry 8.0 m^2 each
trib_areas = ops.eleResponse(100, "tributaryAreas")
```

### 6.2 In-Plane Diaphragm with L-Shaped Polygon and Opening

```python
ops.element(
    "LadrunoMembrane", 200,
    "-nodes", 1, 2, 3, 4, 5, 6,
    "-hole", 11, 12, 13, 14,
    "-frames", 10, 11, 12, 13, 14, 15,
    "-method", "fea_poisson",
    "-pressure", -8.0,
    "-stiff", 28e6, 0.20, 0.15,   # E = 28 GPa, nu = 0.20, t = 0.15 m
    "-rho", 2.5,                  # 2.5 t/m^3 -> 0.375 t/m^2 areal mass
)
```

---

## 7. `apeGmsh` Preprocessor Contract Row

When `apeGmsh` exports a floor/roof diaphragm physical group to OpenSees Tcl or OpenSeesPy, it emits `LadrunoMembrane` according to the following normative contract:

| `apeGmsh` Field / Physical Group | Emitted `LadrunoMembrane` Syntax | Contract & Ordering Guarantees |
|---|---|---|
| `Physical Surface("Slab_<id>_<method>")` outer loop | `-nodes $n1 $n2 ... $nM` | Ordered boundary loop nodes from Gmsh curve loop `0` (CCW or CW; `LadrunoMembrane` auto-orients in local plane). No interior Gmsh surface mesh nodes need to be emitted into the OpenSees domain because `LadrunoMembrane` performs its own internal PSLG/CDT triangulation and Guyan static condensation. |
| Interior curve loops `1..H` on the same surface | `-hole $h1 $h2 ...` (one `-hole` flag per interior loop) | Ordered loop node tags of each opening/shaft cutout. |
| `Physical Line` frame elements bordering or crossing the surface | `-frames $f1 $f2 ...` | All 1D frame element tags (`elasticBeamColumn`, `forceBeamColumn`, `dispBeamColumn`, `TimoshenkoBeamColumn`) whose endpoints lie in the slab plane within $10^{-3} L_{\text{char}}$. Declaration order between `element LadrunoMembrane` and `element *Beam*` is unrestricted (`resolveAndBuildTransfer()` is lazy). |
| Isolated column/wall top nodes in the surface plane | `-supportNodes $sn1 $sn2 ...` | Point-support node tags (`ndf == 3` or `ndf == 6`) that touch the slab without an intervening beam element. |
| Distribution strategy metadata (`method`, `span_dir`, `ortho_ratio`) | `-method <mode> [-spanAxis $dx $dy $dz] [-orthoRatio $eta]` | Defaults to `-method tributary` when unspecified; emits `-spanAxis` whenever `method == "oneway"` or `ortho_ratio < 1.0`. |
| Diaphragm membrane property (`E`, `nu`, `t`, `rho`, `lumped`) | `[-stiff $E $nu $t] [-rho $rho] [-lumped]` or `[-thick $t -mat $matTag]` | Omitted for pure load-distribution-only (`t == 0`) panels; included when the floor acts as a semi-rigid in-plane CST diaphragm. |
| Surface gravity / live pressure load | `-pressure $q` + `eleLoad -ele $eleTag -type -surfaceLoad` | Emitted per `LoadPattern` so Dead, Superimposed Dead, and Live patterns combine cleanly through `m_activeBeamLoads` without double-counting. |
