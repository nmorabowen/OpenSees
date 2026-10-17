---
title: "LadrunoMembrane & Non-Uniform Frame Loads — apeGmsh Emitter & Postprocessor Guide"
project: Ladruno
status: ready-for-apegmsh
priority: medium
adr: ADR 98
tags:
  - apegmsh
  - emitter
  - postprocessor
  - membrane
  - tributary
  - beam-loads
  - etabs-import
---

# `LadrunoMembrane` & Partial Trapezoidal Frame Loads — `apeGmsh` Implementation Guide

This guide specifies how **`apeGmsh`** should emit and post-process the two capabilities shipped in **WP-098** ([`98_ladruno_membrane_load_distribution_adr.md`](98_ladruno_membrane_load_distribution_adr.md) and [`LadrunoMembrane_guide.md`](LadrunoMembrane_guide.md)):

1. **Non-Uniform (Partial / Trapezoidal) Frame Loads** (`LOAD_TAG_Beam2dPartialUniformLoad = 12`, `LOAD_TAG_Beam3dPartialUniformLoad = 121`) across all 2D and 3D OpenSees frame elements (`ElasticBeam`, `ModElasticBeam3d`, `ElasticTimoshenkoBeam`, `DispBeamColumn`, `LadrunoDispBeamColumn`, `ForceBeamColumn`, `ForceBeamColumnCBDI3d`, `ForceBeamColumn3dThermal`).
2. **ETABS-Style Load-Distributing Membrane Element** (`element LadrunoMembrane`, `ELE_TAG_LadrunoMembrane = 33022`) supporting 4 load-transfer modes (`tributary`, `oneway`, `fea_poisson`, `fea_plate`), arbitrary $N$-gon perimeter polygons (`-nodes`), interior openings (`-hole`), supporting frame elements (`-frames`), isolated column/wall point supports (`-supportNodes`), and optional Guyan-condensed in-plane CST diaphragm stiffness and mass (`-thick`, `-mat`, `-E`, `-nu`, `-rho`, `-lumped`).

---

## Part 1 — Emitter Side (`apeGmsh.opensees` & `etabs_import.py`)

### 1.1 Emitting Partial Trapezoidal Loads on Frame Elements

Any `BeamUniformLoad` primitive in `apeGmsh` now supports partial span bounds `[aOverL, bOverL]` ($0 \le a/L < b/L \le 1$) and distinct start/end intensities $(w_a, w_b)$ on **every** 2D and 3D frame element type in Ladruno OpenSees.

#### 3D Frame Elements (`ndm == 3`)
- **Full-span uniform load (3 numeric args after `-beamUniform`):**
  ```tcl
  eleLoad -ele $beamTag -type -beamUniform $Wy $Wz $Wx
  ```
  ```python
  ops.eleLoad("-ele", beam_tag, "-type", "-beamUniform", Wy, Wz, Wx)
  ```
- **Partial uniform load (5 numeric args after `-beamUniform`):**
  ```tcl
  eleLoad -ele $beamTag -type -beamUniform $Wy $Wz $Wx $aOverL $bOverL
  ```
  ```python
  ops.eleLoad("-ele", beam_tag, "-type", "-beamUniform", Wy, Wz, Wx, a_over_L, b_over_L)
  ```
- **Partial trapezoidal load (8 numeric args after `-beamUniform`):**
  ```tcl
  eleLoad -ele $beamTag -type -beamUniform $Wya $Wza $Wxa $aOverL $bOverL $Wyb $Wzb $Wxb
  ```
  ```python
  ops.eleLoad(
      "-ele", beam_tag, "-type", "-beamUniform",
      Wya, Wza, Wxa, a_over_L, b_over_L, Wyb, Wzb, Wxb,
  )
  ```
  > **Argument order note (3D):** Notice that OpenSees `Beam3dUniformLoad` places `[aOverL, bOverL]` at indices `3, 4` (immediately after `Wya, Wza, Wxa`), followed by the end intensities `Wyb, Wzb, Wxb` at indices `5, 6, 7`.

#### 2D Frame Elements (`ndm == 2`)
- **Partial trapezoidal load (6 numeric args after `-beamUniform`):**
  ```tcl
  eleLoad -ele $beamTag -type -beamUniform $Wya $Wxa $aOverL $bOverL $Wyb $Wxb
  ```
  ```python
  ops.eleLoad("-ele", beam_tag, "-type", "-beamUniform", Wya, Wxa, a_over_L, b_over_L, Wyb, Wxb)
  ```

---

### 1.2 Emitting `element LadrunoMembrane`

When an area object in `apeGmsh` (or imported from ETABS `$2k` / `.edb`) is configured as a load-distribution membrane or deck (rather than a meshed `ShellMITC4` / `ShellDKGQ` floor), `apeGmsh` does **not** need to mesh the interior of the slab panel in the global OpenSees model. It emits a single `LadrunoMembrane` element referencing the perimeter corner nodes and the supporting beam/column tags:

#### Tcl Syntax
```tcl
element LadrunoMembrane $eleTag \
    -nodes $n1 $n2 $n3 ... \
    [-frames $f1 $f2 ...] \
    [-supportNodes $s1 $s2 ...] \
    [-hole $h1 $h2 $h3 ...] \
    [-method tributary | oneway | fea_poisson | fea_plate] \
    [-spanAxis $ax $ay $az] \
    [-orthoRatio $r] \
    [-segments $nSeg] \
    [-meshSize $h] \
    [-thick $t] [-mat $matTag] [-E $E] [-nu $nu] \
    [-rho $rho] [-lumped] \
    [-pressure $p] [-b1 $b1] [-b2 $b2] [-b3 $b3]
```

#### Python (`openseespy`) Syntax
```python
ops.element(
    "LadrunoMembrane", ele_tag,
    "-nodes", *poly_node_tags,
    "-frames", *frame_ele_tags,
    "-supportNodes", *support_node_tags,   # optional
    "-hole", *hole_loop_1_node_tags,       # optional (repeatable for multiple openings)
    "-method", method,                     # "tributary" | "oneway" | "fea_poisson" | "fea_plate"
    "-spanAxis", ax, ay, az,               # optional (used by "oneway" and orthotropic FEA)
    "-orthoRatio", r,                      # optional (D22/D11 in [0.01, 100.0])
    "-segments", n_seg,                    # optional (0 = exact CDT sub-segments)
    "-meshSize", h,                        # optional (internal 2D CDT target edge size)
    "-thick", t, "-mat", mat_tag,          # optional in-plane CST diaphragm
    "-rho", rho,                           # optional mass density
    "-pressure", p,                        # optional nominal surface pressure
)
```

#### Applying Pattern Loads to `LadrunoMembrane`
Inside any `pattern Plain $patTag $tsTag { ... }`, `LadrunoMembrane` accepts standard OpenSees elemental loads and automatically distributes them as scoped `Beam3dPartialUniformLoad` / `Beam2dPartialUniformLoad` instances onto its `-frames`:
1. **Gravity / Self-Weight Multiplier (`-selfWeight`):**
   ```tcl
   eleLoad -ele $memTag -type -selfWeight $gx $gy $gz
   ```
   Applies global traction $\mathbf{q} = \text{loadFactor} \cdot w_{\text{base}} \cdot (g_x, g_y, g_z)$ where:
   - $w_{\text{base}} = \rho \cdot t$ if both $\rho > 0$ and $t > 0$,
   - else $w_{\text{base}} = \text{pressure}$ if `-pressure $p` is non-zero,
   - else $w_{\text{base}} = 1.0$ (so `eleLoad -ele $memTag -type -selfWeight 0.0 0.0 -5.0` directly applies $-5.0\text{ kPa}$ global $Z$ traction when $\rho = 0, p = 0$).
2. **Normal Surface Pressure (`-surfaceLoad`):**
   ```tcl
   eleLoad -ele $memTag -type -surfaceLoad
   ```
   Applies $\mathbf{q} = \text{loadFactor} \cdot (p_{\text{base}} \hat{n} + (b_1, b_2, b_3))$ where $p_{\text{base}} = (\text{pressure} \ne 0) ? \text{pressure} : 1.0$.
3. **Directional Surface Traction (`-beamUniform`):**
   ```tcl
   eleLoad -ele $memTag -type -beamUniform $qy_local $qz_normal $qx_local
   ```
   Applies traction components $(q_{x,\text{local}}\,\hat{e}_1 + q_{y,\text{local}}\,\hat{e}_2 + q_{z,\text{normal}}\,\hat{n}) \times \text{loadFactor}$ per load pattern (ideal when Dead and Live patterns have different uniform pressures on the same membrane without changing `-pressure`).

#### ETABS Area Object Mapping (`etabs_import.py`)
| ETABS Shell / Deck Type | Recommended `LadrunoMembrane` Configuration |
|---|---|
| `Membrane` (Two-Way Load Distribution, No Out-of-Plane Bending) | `-method tributary` (or `-method fea_poisson` for smooth variational distribution) + `-thick $t -mat $matTag` for in-plane diaphragm stiffness |
| `Membrane` / `Deck` (One-Way Load Distribution, `Local Axis 1 = alpha`) | `-method oneway -spanAxis $ax $ay $az` |
| `Plate` / `Shell-Thin` used strictly as unmeshed load-transfer panel (incl. cantilevers) | `-method fea_plate` (automatically clamps supported beam edge and transfers overhang torsion when support rank $< 2$) |
| Unmeshed pure load-distribution panel ("None" stiffness in ETABS) | Omit `-thick` and `-mat` ($K = 0, M = 0$; acts purely as a zero-stiffness load distributor) |

> **Declaration Order Note:** `LadrunoMembrane` resolves `-frames` **lazily** on first load/stiffness/response evaluation (`RP1 HIGH-1`), so `apeGmsh` can emit `element LadrunoMembrane` either before or after the frame elements.

---

## Part 2 — Postprocessor Side: Exact Internal Force Diagram Reconstruction (`V(x)`, `M(x)`, `N(x)`, `T(x)`)

### 2.1 Why End Forces Alone Are Not Enough

For any frame element $e$ of length $L$ spanning from Node $I$ ($x = 0$) to Node $J$ ($x = L$):
- Calling `ops.eleResponse(beamTag, "localForce")` returns the **12 basic+local end forces** at the element ends:
  $$\mathbf{f}_{\text{local}} = \bigl[\underbrace{N_I,\, V_{yI},\, V_{zI},\, T_I,\, M_{yI},\, M_{zI}}_{\text{Node } I\ (x = 0)},\; \underbrace{N_J,\, V_{yJ},\, V_{zJ},\, T_J,\, M_{yJ},\, M_{zJ}}_{\text{Node } J\ (x = L)}\bigr]$$
- Those 12 numbers already include the **exact fixed-end reactions** of all partial trapezoidal loads (both direct `eleLoad`s and loads injected by `LadrunoMembrane`).
- However, to plot the **internal force diagrams along the member** ($x \in [0, L]$), linearly interpolating between Node $I$ and Node $J$ is **wrong** whenever distributed loads act on the span:
  - Under a constant load $w(x) = w_0$, shear $V(x)$ is linear and moment $M(x)$ is quadratic.
  - Under a trapezoidal / triangular load $w(x) = w_a + \frac{w_b - w_a}{b - a}(x - a)$ on $[a, b]$, shear $V(x)$ is **piecewise-quadratic** and bending moment $M(x)$ is **piecewise-cubic**!

---

### 2.2 What `LadrunoMembrane` Exposes via `eleResponse`

To make post-processing completely self-contained without requiring `apeGmsh` to re-run computational geometry or local FEA, `LadrunoMembrane` exposes the following `eleResponse` queries:

| Query Token | Length | Returned `Vector` Layout & Units |
|---|---|---|
| `"subSegmentLoads"` (alias `"beamSubLoads"`) | $9 \times N_{\text{sub}}$ (min 1) | **Exact partial trapezoidal sub-loads emitted onto supporting beams**, already transformed into each target beam's local coordinate axes $(\hat{x}_b, \hat{y}_b, \hat{z}_b)$ and scaled by the active load patterns:<br>`[frameTag, aOverL, bOverL, wya, wyb, wza, wzb, wxa, wxb, ...]` |
| `"subSegments"` | $1$ | `[N_sub]` — number of emitted partial trapezoidal sub-loads. |
| `"frameLoads"` | $4 \times N_{\text{frames}}$ | Single moment-conserving equivalent trapezoid per frame in `-frames` declaration order:<br>`[aOverL_1, bOverL_1, waEq_1, wbEq_1, ...]` |
| `"overhangTorques"` | $N_{\text{frames}}$ | Total cantilever overhang torque `[T_1, ..., T_K]` about each frame's local $\hat{x}_b$ axis when `supportRank < 2` (`fea_plate` cantilever). |
| `"supportReactions"` (alias `"nodalReactions"`) | $N_{\text{supportNodes}}$ | Normal force reaction at each node in `-supportNodes`. |
| `"tributaryAreas"` (alias `"tributaryWidths"`) | $N_{\text{frames}} + N_{\text{supportNodes}}$ | Effective tributary area $A_k$ carried by each frame (`0 .. nF-1`) followed by each support node (`nF .. nF+nS-1`). Satisfies $\sum A_k = A_{\text{net}}$. |
| `"totalLoad"` | $3$ | `[Fx, Fy, Fz]` — total active global force vector applied to the membrane across all active load patterns. |
| `"area"` | $1$ | `[A_net]` — net polygon area (outer polygon minus `-hole` openings). |
| `"localAxes"` | $9$ | `[e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz]` — orthonormal in-plane basis $(\hat{e}_1, \hat{e}_2)$ and unit normal $\hat{n}$. |
| `"meshStats"` | $4$ | `[numVertices, numTriangles, numSubSegments, supportRank]` — internal CDT mesh telemetry. |
| `"meshVertices"` | $3 \times N_{\text{verts}}$ | `[X_0, Y_0, Z_0, X_1, Y_1, Z_1, ...]` — 3D global coordinates of every internal Chew-Ruppert CDT vertex (for optional viewport mesh overlay). |
| `"meshTriangles"` | $3 \times N_{\text{tris}}$ | `[v0_0, v1_0, v2_0, ...]` — 0-based triangle vertex indices into `"meshVertices"`. |
| `"stress"` / `"strain"` | $3$ | `[sxx, syy, sxy]` or `[exx, eyy, gxy]` — in-plane membrane stress/strain when `-thick` and `-mat` are active. |

---

### 2.3 Closed-Form Integration of Partial Trapezoidal Segment Loads

Suppose a 3D frame element `beamTag` of length $L$ has local end forces at Node $I$ ($x = 0$) from `f_local = ops.eleResponse(beamTag, "localForce")`:
- $N_I = f_{\text{local}}[0]$ (axial force at Node $I$, positive in $+x$ direction, i.e., compression at Node $I$)
- $V_{yI} = f_{\text{local}}[1]$ (local-$y$ force at Node $I$)
- $V_{zI} = f_{\text{local}}[2]$ (local-$z$ force at Node $I$)
- $T_I = f_{\text{local}}[3]$ (torsion at Node $I$)
- $M_{yI} = f_{\text{local}}[4]$ (moment about local $y$ at Node $I$)
- $M_{zI} = f_{\text{local}}[5]$ (moment about local $z$ at Node $I$)

And suppose the beam carries a list of partial trapezoidal distributed loads (collected from `ops.eleResponse(memTag, "subSegmentLoads")` for all `LadrunoMembrane` elements referencing `beamTag`, plus any direct `Beam3dPartialUniformLoad`s emitted on `beamTag`):
$$\mathcal{S}(\text{beamTag}) = \bigl\{\bigl(a_m,\, b_m,\, w_{ya,m},\, w_{yb,m},\, w_{za,m},\, w_{zb,m},\, w_{xa,m},\, w_{xb,m}\bigr)\bigr\}_{m=1}^{M}$$
where $a_m = (a/L)_m \cdot L$ and $b_m = (b/L)_m \cdot L$.

#### Contribution of One Partial Trapezoid $[a, b]$ at Station $x \in [0, L]$

Let $c = b - a > 0$ be the loaded segment length, and let $w(s)$ vary linearly from $w_a$ at $s = a$ to $w_b$ at $s = b$:
$$w(s) = w_a + \frac{w_b - w_a}{c}(s - a), \qquad s \in [a, b]$$

Define the **integrated shear contribution** $\Delta F(x; a, b, w_a, w_b) = \int_a^{\min(x, b)} w(s)\,ds$ and **integrated moment contribution** $\Delta M(x; a, b, w_a, w_b) = \int_a^{\min(x, b)} (x - s)\,w(s)\,ds$:

1. **Before the load ($x \le a$):**
   $$\Delta F(x) = 0, \qquad \Delta M(x) = 0$$

2. **Inside the trapezoid ($a < x < b$):**
   Let $\xi = x - a \in (0, c)$ and $k = \frac{w_b - w_a}{c}$ be the load slope:
   $$\Delta F(x) = w_a \xi + \frac{1}{2} k \xi^2 \quad \text{(quadratic in } x\text{)}$$
   $$\Delta M(x) = \frac{1}{2} w_a \xi^2 + \frac{1}{6} k \xi^3 \quad \text{(cubic in } x\text{)}$$

3. **Past the load ($x \ge b$):**
   The total segment force is $F_{\text{seg}} = \frac{1}{2}(w_a + w_b) c$ and its moment about $s = a$ is $M_{a,\text{seg}} = \frac{c^2}{6}(w_a + 2 w_b)$, so:
   $$\Delta F(x) = F_{\text{seg}}$$
   $$\Delta M(x) = F_{\text{seg}} (x - a) - M_{a,\text{seg}} = F_{\text{seg}} \left(x - a - \frac{c(w_a + 2 w_b)}{3(w_a + w_b)}\right) \quad \text{(linear in } x\text{)}$$

#### Complete Internal Force Diagrams Along $x \in [0, L]$

Summing over all partial trapezoidal segments $m \in \mathcal{S}(\text{beamTag})$ gives the **exact** internal section forces at any station $x \in [0, L]$ (matching `ForceBeamColumn3d`'s internal Gauss-Lobatto section forces to machine precision):

$$\begin{aligned}
N(x) &= -N_I - \sum_{m} \Delta F(x;\, a_m, b_m, w_{xa,m}, w_{xb,m}) \\
V_y(x) &= V_{yI} + \sum_{m} \Delta F(x;\, a_m, b_m, w_{ya,m}, w_{yb,m}) \\
V_z(x) &= V_{zI} + \sum_{m} \Delta F(x;\, a_m, b_m, w_{za,m}, w_{zb,m}) \\
T(x) &= -T_I \\
M_z(x) &= -M_{zI} + V_{yI}\,x + \sum_{m} \Delta M(x;\, a_m, b_m, w_{ya,m}, w_{yb,m}) \\
M_y(x) &= -M_{yI} - V_{zI}\,x - \sum_{m} \Delta M(x;\, a_m, b_m, w_{za,m}, w_{zb,m})
\end{aligned}$$
*(Note the minus signs on $-V_{zI}\,x - \sum \Delta M_z$ in $M_y(x)$: in a right-handed 3D frame $(\hat{x}, \hat{y}, \hat{z})$, $\hat{x} \times \hat{y} = +\hat{z}$ whereas $\hat{x} \times \hat{z} = -\hat{y}$, so OpenSees `CrdTransf3d` defines $V_{zI} = -(M_{yI} + M_{yJ})/L + p_0[3]$.)*

> **End-Point Consistency Check at $x = L$:**
> Because OpenSees computed $f_{\text{local}}$ using the exact Hermite fixed-end forces of the same partial trapezoids, evaluating the formulas above at $x = L$ automatically satisfies:
> $$N(L) = N_J, \quad V_y(L) = -V_{yJ}, \quad V_z(L) = -V_{zJ}, \quad M_z(L) = M_{zJ}, \quad M_y(L) = M_{yJ}$$
> to machine precision ($\sim 10^{-12}$) on every element type (`ElasticBeam3d`, `ModElasticBeam3d`, `ElasticTimoshenkoBeam3d`, `DispBeamColumn3d`, `LadrunoDispBeamColumn3d`, `ForceBeamColumn3d`).

---

### 2.4 Drop-In Python Reference Implementation for `apeGmsh`

The function below can be dropped directly into `apeGmsh`'s postprocessor to collect all membrane sub-loads and evaluate exact internal force diagrams at arbitrary sampling stations `x_stations`:

```python
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class PartialTrapezoidSubLoad:
    frame_tag: int
    a_over_L: float
    b_over_L: float
    wya: float
    wyb: float
    wza: float
    wzb: float
    wxa: float
    wxb: float


def collect_membrane_beam_subloads(ops, membrane_tags: list[int]) -> dict[int, list[PartialTrapezoidSubLoad]]:
    """Query all LadrunoMembrane elements for their active beam sub-segment loads."""
    by_beam: dict[int, list[PartialTrapezoidSubLoad]] = {}
    for mem_tag in membrane_tags:
        n_sub_raw = ops.eleResponse(mem_tag, "subSegments")
        if not n_sub_raw or int(round(n_sub_raw[0])) <= 0:
            continue
        n_sub = int(round(n_sub_raw[0]))
        raw = ops.eleResponse(mem_tag, "subSegmentLoads")
        if not raw or len(raw) < 9 * n_sub:
            continue
        for s in range(n_sub):
            base = 9 * s
            f_tag = int(round(raw[base + 0]))
            sub = PartialTrapezoidSubLoad(
                frame_tag=f_tag,
                a_over_L=float(raw[base + 1]),
                b_over_L=float(raw[base + 2]),
                wya=float(raw[base + 3]),
                wyb=float(raw[base + 4]),
                wza=float(raw[base + 5]),
                wzb=float(raw[base + 6]),
                wxa=float(raw[base + 7]),
                wxb=float(raw[base + 8]),
            )
            by_beam.setdefault(f_tag, []).append(sub)
    return by_beam


def _integrate_trapezoid_array(
    x: np.ndarray, a: float, b: float, wa: float, wb: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return (dF(x), dM(x)) arrays at stations x for a partial trapezoid on [a, b]."""
    c = b - a
    dF = np.zeros_like(x, dtype=float)
    dM = np.zeros_like(x, dtype=float)
    if c <= 1e-14 or (abs(wa) + abs(wb) <= 1e-15):
        return dF, dM

    k = (wb - wa) / c
    inside = (x > a) & (x < b)
    xi = x[inside] - a
    dF[inside] = wa * xi + 0.5 * k * (xi ** 2)
    dM[inside] = 0.5 * wa * (xi ** 2) + (1.0 / 6.0) * k * (xi ** 3)

    past = x >= b
    F_seg = 0.5 * (wa + wb) * c
    Ma_seg = (c ** 2) * (wa + 2.0 * wb) / 6.0
    dF[past] = F_seg
    dM[past] = F_seg * (x[past] - a) - Ma_seg
    return dF, dM


def reconstruct_beam_3d_diagrams(
    L: float,
    local_force_12: list[float] | np.ndarray,
    sub_loads: list[PartialTrapezoidSubLoad],
    num_stations: int = 51,
) -> dict[str, np.ndarray]:
    """
    Reconstruct exact 3D internal force diagrams along a frame element of length L.
    Automatically includes all segment transition breakpoints (a_m, b_m) in x!
    """
    Lf = np.asarray(local_force_12, dtype=float)
    N_I, Vy_I, Vz_I, T_I, My_I, Mz_I = Lf[0], Lf[1], Lf[2], Lf[3], Lf[4], Lf[5]

    breakpoints = [0.0, L]
    for sub in sub_loads:
        breakpoints.append(np.clip(sub.a_over_L * L, 0.0, L))
        breakpoints.append(np.clip(sub.b_over_L * L, 0.0, L))
    x = np.unique(np.concatenate([np.linspace(0.0, L, num_stations), np.array(breakpoints)]))

    N = np.full_like(x, -N_I)
    Vy = np.full_like(x, Vy_I)
    Vz = np.full_like(x, Vz_I)
    T = np.full_like(x, -T_I)
    Mz = -Mz_I + Vy_I * x
    My = -My_I - Vz_I * x

    for sub in sub_loads:
        a = sub.a_over_L * L
        b = sub.b_over_L * L
        dFx, _ = _integrate_trapezoid_array(x, a, b, sub.wxa, sub.wxb)
        dFy, dMy_z = _integrate_trapezoid_array(x, a, b, sub.wya, sub.wyb)
        dFz, dMz_y = _integrate_trapezoid_array(x, a, b, sub.wza, sub.wzb)

        N -= dFx
        Vy += dFy
        Vz += dFz
        Mz += dMy_z
        My -= dMz_y

    return {"x": x, "N": N, "Vy": Vy, "Vz": Vz, "T": T, "My": My, "Mz": Mz}
```

---

## Part 3 — Implementation Checklist for `apeGmsh`

1. **Primitive & Schema (`apeGmsh/opensees/primitives.py`):**
   - Extend `BeamUniformLoad` to accept optional `a_over_L: float = 0.0`, `b_over_L: float = 1.0`, `wyb: float | None = None`, `wzb: float | None = None`, `wxb: float | None = None`. Emit 3 args when uniform full-span, 5 args when uniform partial-span, and 8 args (3D) / 6 args (2D) when trapezoidal.
   - Add `LadrunoMembrane` element primitive accepting `nodes`, `frames`, `support_nodes`, `holes`, `method`, `span_axis`, `ortho_ratio`, `segments`, `mesh_size`, `thick`, `mat_tag`, `E`, `nu`, `rho`, `lumped`, `pressure`, `b1`, `b2`, `b3`.
2. **ETABS Import (`apeGmsh/importers/etabs_import.py`):**
   - Map ETABS Membrane / Deck / unmeshed Floor panels to `LadrunoMembrane`, populating `-frames` from the bounding/interior coplanar beams and `-method` (`tributary` for two-way membrane, `oneway` + `-spanAxis` for one-way deck).
3. **Postprocessor Diagram Engine (`apeGmsh/post/diagrams.py`):**
   - Query `ops.eleResponse(memTag, "subSegmentLoads")` after each analysis step (or persist in `results.h5` under `/Membranes/subSegmentLoads`) and pass the per-beam sub-segment list into `reconstruct_beam_3d_diagrams()` when rendering internal force diagrams ($N, V_y, V_z, T, M_y, M_z$).
4. **Optional 3D Viewport Tributary / CDT Mesh Overlay:**
   - Query `ops.eleResponse(memTag, "meshVertices")` and `ops.eleResponse(memTag, "meshTriangles")` to render the internal Chew-Ruppert CDT triangulation and tributary load arrows onto the supporting beams.
