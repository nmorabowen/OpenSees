---
title: ADR 98 — Non-Uniform Frame Loading and LadrunoMembrane Surface-to-Frame Load Distribution
project: Ladruno
status: accepted
priority: high
owner: nmora
tags:
  - implementation
  - element
  - membrane
  - beam
  - load-distribution
related:
  - "[[87_ladruno_warrant_and_architecture_adr]]"
  - "[[91_ladruno_shell_stiffness_modifiers_adr]]"
  - "[[96_ladruno_contact_passenger_dof_adr]]"
---

# ADR 98 — Non-Uniform Frame Loading and `LadrunoMembrane` Surface-to-Frame Load Distribution

## 1. What & Why

Floor and roof diaphragms in structural frame models require two complementary capabilities that upstream OpenSees lacks or implements with defects:

1. **Exact Non-Uniform (Partial Trapezoidal) Frame Element Loads (`Beam2dPartialUniformLoad` / `Beam3dPartialUniformLoad`):**
   Upstream `ElasticBeam3d`, `ModElasticBeam3d`, `DispBeamColumn2d/3d`, `ForceBeamColumnCBDI3d`, and `ElasticTimoshenkoBeam2d/3d` either drop the trapezoidal end-intensities `(wyb, wzb, wab)` in `Beam3dPartialUniformLoad`, omit `LOAD_TAG_Beam2d/3dPartialUniformLoad` from `FEM_ObjectBrokerAllClasses`, ignore moment end-releases (`release`, `releasey`, `releasez`), ignore the Timoshenko shear flexibility factor $\Phi = \frac{12EI}{GA_s L^2}$, or invert the sign of the triangular shear term in `ForceBeamColumn3d::getElementLoad` (`SECTION_RESPONSE_VZ`).
2. **Surface-to-Frame Load Distribution & Optional In-Plane Diaphragm Element (`LadrunoMembrane` + `LadrunoSurfaceLoad`):**
   Engineers need an arbitrary $N$-node polygonal membrane/floor panel (supporting non-convex L-shapes, oblique edges, interior `-hole` openings, partial/interior supporting beams `-frames`, and isolated column/wall `-supportNodes`) that deterministically distributes surface pressures $q$ onto supporting frame elements as moment-conserving partial trapezoidal `eleLoad`s under four selectable mechanics modes (`tributary`, `oneway`, `fea_poisson`, `fea_plate`), while optionally providing Guyan-condensed in-plane CST diaphragm stiffness and mass (`-stiff $E $nu $t [-rho $rho]`).

---

## 2. Normative Mathematical Specification

### 2.1 Exact Partial Trapezoidal Fixed-End Forces (Euler-Bernoulli, Timoshenko $\Phi$, and End Releases)

Let a linear trapezoidal load act over $x \in [a, b] \subset [0, L]$ with normalized bounds $\alpha = a/L$, $\beta = b/L$ ($0 \le \alpha < \beta \le 1$) and linear intensity $w(\xi) = w_a + (w_b - w_a)\frac{\xi - \alpha}{\beta - \alpha}$ for $\xi = x/L \in [\alpha, \beta]$.

Using the exact 3-point Gauss-Legendre rule on $u \in [0, 1]$ ($\xi(u) = \alpha + (\beta - \alpha)u$, exact for polynomial degree $\le 5$, eliminating subtraction cancellation when $\beta - \alpha \ll 1$):

- **Simply-supported basic reactions (`p0`):**
  $$F = \frac{1}{2}(w_a + w_b)(b - a), \qquad V_1^{\text{ss}} = L(\beta - \alpha)\sum_{g=1}^3 \omega_g w(u_g)(1 - \xi_g), \qquad V_2^{\text{ss}} = L(\beta - \alpha)\sum_{g=1}^3 \omega_g w(u_g)\xi_g$$
  $$P = \frac{1}{2}(w_{aa} + w_{ab})(b - a), \qquad P_J = L(\beta - \alpha)\sum_{g=1}^3 \omega_g w_{\text{ax}}(u_g)\xi_g, \qquad P_I = P - P_J$$
- **Fixed-end basic moments with Timoshenko shear factor $\Phi = \frac{12EI}{G A_s L^2} \ge 0$:**
  From the exact Timoshenko Hermite rotation shape functions $N_{\theta_1}^{(\Phi)}(\xi) = \frac{L}{1+\Phi}\xi(1-\xi)(1-\xi+\Phi/2)$ and $N_{\theta_2}^{(\Phi)}(\xi) = -\frac{L}{1+\Phi}\xi(1-\xi)(\xi+\Phi/2)$:
  $$M_{z1}^{\text{fixed}} = -\frac{L^2(\beta - \alpha)}{1 + \Phi_z}\sum_{g=1}^3 \omega_g w_y(u_g)\xi_g(1-\xi_g)\!\left(1 - \xi_g + \frac{\Phi_z}{2}\right)$$
  $$M_{z2}^{\text{fixed}} = +\frac{L^2(\beta - \alpha)}{1 + \Phi_z}\sum_{g=1}^3 \omega_g w_y(u_g)\xi_g(1-\xi_g)\!\left(\xi_g + \frac{\Phi_z}{2}\right)$$
  and for 3D local-$z$ loading ($x \times z = -y$, so $+w_z$ induces $+M_{y1}$ at node I and $-M_{y2}$ at node J):
  $$M_{y1}^{\text{fixed}} = +\frac{L^2(\beta - \alpha)}{1 + \Phi_y}\sum_{g=1}^3 \omega_g w_z(u_g)\xi_g(1-\xi_g)\!\left(1 - \xi_g + \frac{\Phi_y}{2}\right)$$
  $$M_{y2}^{\text{fixed}} = -\frac{L^2(\beta - \alpha)}{1 + \Phi_y}\sum_{g=1}^3 \omega_g w_z(u_g)\xi_g(1-\xi_g)\!\left(\xi_g + \frac{\Phi_y}{2}\right)$$
- **End-Release Reduction (`release`, `releasez`, `releasey` $\in \{0, 1, 2, 3\}$):**
  For any basic moment pair $(M_1^{\text{fixed}}, M_2^{\text{fixed}})$:
  - `0` (fixed-fixed): $(M_1, M_2) = (M_1^{\text{fixed}}, M_2^{\text{fixed}})$
  - `1` (pinned at I): $(M_1, M_2) = \left(0, \; M_2^{\text{fixed}} - \frac{1}{2}M_1^{\text{fixed}}\right)$
  - `2` (pinned at J): $(M_1, M_2) = \left(M_1^{\text{fixed}} - \frac{1}{2}M_2^{\text{fixed}}, \; 0\right)$
  - `3` (pinned-pinned): $(M_1, M_2) = (0, 0)$
  In the basic system, `p0 = [-P, -Vy1_ss, -Vy2_ss, -Vz1_ss, -Vz2_ss]` and `q0 = [-PJ, Mz1, Mz2, My1, My2]`; `CrdTransf` (and local-system vectors `ql0`) assemble the total end shears as $V_{y1} = -V_{y1}^{\text{ss}} + (M_{z1} + M_{z2})/L$, $V_{y2} = -V_{y2}^{\text{ss}} - (M_{z1} + M_{z2})/L$, $V_{z1} = -V_{z1}^{\text{ss}} - (M_{y1} + M_{y2})/L$, $V_{z2} = -V_{z2}^{\text{ss}} + (M_{y1} + M_{y2})/L$.
- **`ForceBeamColumn3d` Internal Section Response (`sp`):**
  For $x \in (a, b)$ with $\Delta x = x - a$ and $w_z(x) = w_{za} + \frac{w_{zb} - w_{za}}{b - a}\Delta x$:
  $$M_y(x) = V_{zI} x - \frac{1}{2}w_{za}\Delta x^2 - \frac{1}{6}(w_z(x) - w_{za})\Delta x^2, \qquad V_z(x) = \frac{dM_y}{dx} = V_{zI} - w_{za}\Delta x - \frac{1}{2}(w_z(x) - w_{za})\Delta x$$
  fixing the inverted sign on `0.5*(wzz - wz)*(x - a)` in `ForceBeamColumn3d.cpp`, `ForceBeamColumn3dThermal.cpp`, and `ForceBeamColumnCBDI3d.cpp`.

### 2.2 Four Load-Distribution Modes on `LadrunoMembrane`

All four modes operate on a deterministic 2D Planar Straight-Line Graph (PSLG) and Constrained Delaunay Triangulation (CDT) refined via Chew-Ruppert quality meshing in the slab's local plane $(\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2)$:

1. **`tributary` (CDT-Based Euclidean Voronoi / Nearest-Supported-Feature Integration):**
   Each point $\mathbf{x} \in \Omega$ is assigned to the nearest supporting feature $S^*$ (supporting frame segment or `-supportNodes` point) under Euclidean distance $d(\mathbf{x}, S)$. Triangles whose sample points map to a single interior frame segment are integrated in closed form via 3-point triangle quadrature; triangles straddling a Voronoi bisector between two frame segments are clipped along the linear bisector $d_1(\mathbf{x}) - d_2(\mathbf{x}) = 0$, and multi-feature ridge apex triangles are adaptively subdivided. On convex rectangular bays, this reproduces the exact $45^\circ$ yield-line triangles and trapezoids while remaining 100% robust on non-convex polygons, holes, free edges, and point supports.
2. **`oneway` (Bidirectional Ray-Cast Strip Partitioning):**
   Along the user-specified unit span vector $\hat{\mathbf{s}} = (\cos\theta, \sin\theta)$, rays are cast from quadrature points $\mathbf{x}_q$ in $\pm\hat{\mathbf{s}}$ to find the nearest supporting intersections at distances $(d_L, d_R)$. Two-sided strips partition $q\,dA$ by lever arm $f_L = f_q \frac{d_R}{d_L + d_R}$ and $f_R = f_q \frac{d_L}{d_L + d_R}$ (identically conserving $\sum F_z$, $\sum M_x$, and $\sum M_y$ for every point $\mathbf{x}_q$); single-sided cantilever strips assign 100% of $f_q$ to the intersected support.
3. **`fea_poisson` (1-DOF CST Membrane Deflection Analogy):**
   Solves $-\nabla \cdot (\mathbf{D}\nabla w) = q$ on the CDT mesh with homogeneous Dirichlet boundary conditions $w = 0$ on all supporting frame segments and `-supportNodes`, and natural Neumann condition $\mathbf{n}\cdot\mathbf{D}\nabla w = 0$ on unsupported/hole boundaries, where $\mathbf{D} = \mathbf{R}(\theta)\operatorname{diag}(1, \eta)\mathbf{R}(\theta)^T$ supports both isotropic ($\eta = 1$) and orthotropic ($\eta \ll 1$) load paths.
4. **`fea_plate` (3-DOF Discrete Kirchhoff Triangle Bending):**
   Assembles the 9-DOF Batoz-Bathe-Ho (1980) DKT plate bending stiffness $\mathbf{K}_{\text{DKT}}$ for $(w, \theta_x, \theta_y)$ per node. Support coordinate covariance rank is checked:
   - If $\operatorname{rank} \ge 2$ (non-collinear supports), supported nodes are knife-edge simple supports ($w = 0$, rotations free), preserving signed downward/uplift reactions so $\sum F_z = qA$, $\sum M_x$, and $\sum M_y$ are satisfied exactly.
   - If $\operatorname{rank} < 2$ (single beam or collinear supports), rotations at supported nodes are clamped ($\theta_x = \theta_y = 0$) to prevent a kinematic rigid-body mechanism, and overhanging moments are extracted.

### 2.3 Moment-Conserving Trapezoid Projection (`(F_j, M_j) -> (w_a, w_b)`)

On every constrained sub-segment $j$ of length $L_j = b_j - a_j$ belonging to parent frame $k$, the accumulated force $F_j = \int_0^{L_j} w(s)\,ds$ and first moment $M_j = \int_0^{L_j} s\,w(s)\,ds$ (or in FEM modes, $F_j = R_{j,a} + R_{j,b}$ and $M_j = L_j R_{j,b}$, which is identically the inverse of the 1D consistent mass matrix $\mathbf{M}_{1D}^{-1}$) determine the exact linear trapezoid endpoints $(w_{a,j}, w_{b,j})$:
$$w_{a,j} + w_{b,j} = \frac{2 F_j}{L_j}, \qquad w_{b,j} - w_{a,j} = \frac{12\!\left(M_j - \tfrac{1}{2} F_j L_j\right)}{L_j^2}$$
This guarantees exact conservation of total force $\sum F_z = qA$ and global overturning moments $\sum M_x, \sum M_y$ to machine precision ($< 10^{-12}$).

---

## 3. Architectural Decisions (`RP1` 3-Lens Opus Review Synthesis)

| # | Decision | Rationale |
|---|---|---|
| **D1** | **Scoped Load Instance Arena (`m_activeBeamLoads`)** | Mutating persistent `Beam3dPartialUniformLoad` instances across multiple `LoadPattern`s double-counts loads in `ElasticBeam3d` and `ForceBeamColumn3d`. `LadrunoMembrane::zeroLoad()` clears `m_activeBeamLoads`, and each `addLoad(theLoad, loadFactor)` allocates fresh `std::unique_ptr<ElementalLoad>` instances pre-scaled by `loadFactor` and passes `(rawPtr, 1.0)` to `targetFrame->addLoad()`. |
| **D2** | **Broker Registration for Tags `12` and `121`** | Register `LOAD_TAG_Beam2dPartialUniformLoad` (`12`) and `LOAD_TAG_Beam3dPartialUniformLoad` (`121`) in `FEM_ObjectBrokerAllClasses::getNewElementalLoad()` for parallel serialization (`sendSelf`/`recvSelf`). |
| **D3** | **Lazy Frame Resolution (`resolveFrames()`)** | `setDomain()` resolves only the membrane's own nodes and sets `m_framesResolved = false`. Supporting frame pointers and the unit load transfer table are resolved lazily on first `addLoad()`, `getTangentStiff()`, `getResistingForce()`, or `getResponse()`, eliminating element declaration order sensitivity. |
| **D4** | **Dynamic Mixed-NDF Offsets (`m_nodeDofOffsets`)** | `LadrunoMembrane` supports mixed 3-DOF and 6-DOF nodes plus `-supportNodes` appended to `connectedExternalNodes`, computing `m_numDOF = sum(theNodes[i]->getNumberDOF())` dynamically. |
| **D5** | **Guyan Static Condensation for In-Plane CST Diaphragm** | Internal Steiner vertices created by CDT refinement do not exist in the OpenSees `Domain`. In-plane CST stiffness $\mathbf{K}$ and consistent mass $\mathbf{M}$ are statically condensed onto the polygon's external boundary nodes via $\mathbf{G} = -\mathbf{K}_{ii}^{-1}\mathbf{K}_{ib}$, $\mathbf{K}_{\text{cond}} = \mathbf{K}_{bb} + \mathbf{K}_{bi}\mathbf{G}$, and $\mathbf{M}_{\text{cond}} = \mathbf{T}^T \mathbf{M} \mathbf{T}$. |
| **D6** | **Element-Edge Diffusive Flux Apportionment (No 50/50 Corner Split)** | Shared corner nodes are never split 50/50; reactions at shared vertices are apportioned across incident constrained sub-segments by each adjacent triangle's integrated normal diffusive flux $L_{S_m}\max(0, -\hat{\mathbf{n}}_{S_m}\cdot\mathbf{D}\nabla w_e) + \epsilon L_{S_m}(\hat{\mathbf{n}}_{S_m}^T\mathbf{D}\hat{\mathbf{n}}_{S_m})$, and projected onto each sub-segment via $\mathbf{M}_{1D}^{-1}$. |
| **D7** | **Support Rank Check in `fea_plate` (DKT)** | Support coordinate singular values $\sigma_1, \sigma_2$ determine whether supports are non-collinear ($\sigma_2 \ge 10^{-6}L_{\text{char}} \implies$ knife-edge $w=0$, rotations free) or collinear ($\operatorname{rank} < 2 \implies$ clamped rotations $\theta_x = \theta_y = 0$ with overhang torque extraction). |
| **D8** | **End-Release Partial-Load Reduction** | `ElasticBeam2d`, `ElasticBeam3d`, and `ModElasticBeam3d` apply exact propped-cantilever moment reductions (`release = 1, 2, 3`) under `Beam2d/3dPartialUniformLoad`. |
| **D9** | **Timoshenko $\Phi$-Dependent FEFs & `ForceBeamColumn3d` $V_z$ Sign Fix** | `ElasticTimoshenkoBeam2d/3d` include $\Phi = \frac{12EI}{GA_s L^2}$ in partial trapezoidal fixed-end moments; `ForceBeamColumn3d`, `Thermal`, and `CBDI3d` fix the sign of the triangular $V_z$ shear term. |
| **D10** | **4-Stage PSLG Preprocessor + Chew-Ruppert Acute Guard** | Vertices are deduplicated, segment intersections and point-on-segment incidences are split, and vertices with interior angle $\alpha < 60^\circ$ are tagged to suppress encroachment ping-pong with floor $h_{\min} = 10^{-4}\sqrt{A}$ and vertex cap. |
