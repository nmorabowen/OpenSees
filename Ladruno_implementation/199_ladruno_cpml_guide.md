# WP-199 — `LadrunoCPML`: Unified 2D/3D Convolutional Perfectly Matched Layer (CFS-PML)

## 1. Overview & Motivation

`LadrunoCPML` (`ELE_TAG_LadrunoCPML = 33023`) is a unified 2D/3D continuum wave-absorbing boundary element based on the **Non-Split Second-Order Complex Frequency-Shifted Convolutional Perfectly Matched Layer (CFS-PML / C-PML)** formulation of Matzen (2011) and Xie, Komatitsch, Martin & Matzen (2014).

Compared to local Lysmer–Kuhlemeyer dashpots (`ASDAbsorbingBoundary2D`/`3D`) and classical split-field / mixed-DOF PML elements (`PML2D`/`PML3D`):
- **Standard Displacement DOFs Only (`ndf = ndm`)**: Uses 2 DOFs/node in 2D and 3 DOFs/node in 3D, connecting directly to standard continuum elements (`FourNodeQuad`, `LadrunoQuad`, `LadrunoCST`, `LadrunoLST`, `BezierTri6`, `stdBrick`, `LadrunoBrick`, `LadrunoBrick20`, `FourNodeTetrahedron`, `TenNodeTetrahedron`, `BezierTet10`) with **no `equalDOF` or multipoint constraints** and **no auxiliary nodal DOFs**.
- **Symmetric Positive-Definite (SPD) Effective Tangents**: All C-PML memory variables live at Gauss points and are updated via exact **Piecewise-Linear Recursive Convolution (`RCFilter1D`)** cascades, yielding strictly symmetric positive-definite effective matrices $\mathbf{K}_{\text{eff}}, \mathbf{C}_{\text{eff}}, \mathbf{M}_{\text{eff}}$ compatible with standard symmetric solvers (`Pardiso`, `Mumps`, `ProfileSPD`, `BandSPD`).
- **Rayleigh Surface Waves, Evanescent Modes & Grazing Incidence**: Complex frequency shifting ($\alpha_{\max} > 0$) and real coordinate scaling ($\kappa_{\max} \ge 1$) absorb Rayleigh surface waves and near-grazing body waves (reducing Lamb's problem reflection error by $> 64\times$ compared to `ASDAbsorbingBoundary2D`) and eliminate the late-time secular drift of classical unshifted PML.
- **Unstructured & Simplicial Mesh Support (`9` Topologies)**: Because the stretched weak form operates in physical Cartesian/local-normal coordinates at each Gauss point while the isoparametric mapping $\mathbf{x}(\boldsymbol{\xi})$ is purely real and geometric, `LadrunoCPML` works out-of-the-box on both structured and unstructured meshes across all 9 topologies:
  * **2D (`ndm = 2, ndf = 2`)**: `T3` (3-node linear triangle, 1 GP), `Q4` (4-node bilinear quad, 4 GPs), `T6` (6-node quadratic Lagrange triangle, 3 GPs), `BT6` (6-node quadratic Bézier triangle, 3 GPs).
  * **3D (`ndm = 3, ndf = 3`)**: `TET4` (4-node linear tetrahedron, 1 GP), `H8` (8-node trilinear brick, 8 GPs), `TET10` (10-node quadratic Lagrange tetrahedron, 4 GPs), `BTET10` (10-node quadratic Bézier tetrahedron, 4 GPs), `H20` (20-node serendipity quadratic brick, 8-point Gauss rule).

---

## 2. Mathematical Formulation

### 2.1 CFS-PML Complex Coordinate Stretching & Polynomial Grading
In the local orthonormal frame $\mathbf{x}' = \mathbf{R}(\mathbf{x} - \mathbf{x}_{\text{ref}})$ ($\mathbf{R} \in SO(\text{ndm})$), each local coordinate $x'_m$ ($m \in \{0, \dots, \text{ndm}-1\}$) is mapped to the complex stretched coordinate $\tilde{x}'_m$ via:
$$\frac{\partial \tilde{x}'_m}{\partial x'_m} = s_m(\omega) = \kappa_m(\xi_m) + \frac{d_m(\xi_m)}{\alpha_m(\xi_m) + i\omega}$$
with inverse stretching factor:
$$s_m^{-1}(\omega) = \frac{1}{\kappa_m} - \frac{d_m / \kappa_m^2}{\hat{\alpha}_m + i\omega}, \qquad \hat{\alpha}_m = \alpha_m + \frac{d_m}{\kappa_m}$$
At each Gauss point $\mathbf{x}_g$ with normalized penetration depth $\bar{\xi}_m = \min(1, \max(0, \xi_m(\mathbf{x}_g)/L_m))$:
- For $\bar{\xi}_m = 0$: $d_m = 0$, $\kappa_m = 1$, $\alpha_m = 0$ (exact reduction to standard linear elasticity).
- For $\bar{\xi}_m > 0$:
  $$d_{0,m} = -\frac{(m_{\text{deg}} + 1) V_p \ln(R_0)}{2 L_m}, \qquad V_p = \sqrt{\frac{\lambda + 2\mu}{\rho}}$$
  $$d_m(\bar{\xi}_m) = d_{0,m} (\bar{\xi}_m)^{m_{\text{deg}}}, \qquad \kappa_m(\bar{\xi}_m) = 1 + (\kappa_{\max} - 1)(\bar{\xi}_m)^{m_{\text{deg}}}, \qquad \alpha_m(\bar{\xi}_m) = \alpha_{\max}(1 - \bar{\xi}_m)$$

### 2.2 Non-Split Second-Order Stretched Weak Form
Multiplying the frequency-domain elastodynamic momentum balance by the stretched Jacobian $J_s(\omega) = \prod_{m=0}^{\text{ndm}-1} s_m(\omega)$ and integrating by parts over $\Omega_e$ yields the symmetric non-split weak form:
$$\int_{\Omega_e} \rho\, w_i\, J_s(\omega)(i\omega)^2 \hat{u}_i \, d\Omega + \int_{\Omega_e} \Phi'_{a,j} \, C_{ijkl} \, F_{jl}(\omega) \hat{G}'_{kl} \, d\Omega = \int_{\Omega_e} w_i \hat{b}_i \, d\Omega$$
where $\Phi'_{a,j} = \sum_p R_{jp} \frac{\partial N_a}{\partial x_p}$, $G'_{kl} = \sum_b u'_{b,k} \Phi'_{b,l}$, and $F_{jl}(\omega) = \frac{J_s(\omega)}{s_j(\omega) s_l(\omega)}$ is a rational filter cascade:
- **In 2D (`ndm = 2`)**: $F_{00} = s_1 \circ s_0^{-1}$, $F_{11} = s_0 \circ s_1^{-1}$, $F_{01} = F_{10} = 1$.
- **In 3D (`ndm = 3`)**: $F_{00} = s_2 \circ s_1 \circ s_0^{-1}$, $F_{11} = s_2 \circ s_0 \circ s_1^{-1}$, $F_{22} = s_1 \circ s_0 \circ s_2^{-1}$, $F_{01} = F_{10} = s_2$, $F_{02} = F_{20} = s_1$, $F_{12} = F_{21} = s_0$.

Using $s_m(\omega) = \kappa_m + \phi_m(\omega)$ with $\phi_m(\omega) = \frac{d_m}{\alpha_m + i\omega}$ and $g_m(\omega) = (i\omega)\phi_m(\omega) = d_m - \frac{\alpha_m d_m}{\alpha_m + i\omega}$, the inertial operator $J_s(\omega)(i\omega)^2 \hat{\mathbf{u}}$ decomposes without third-order time derivatives into:
1. **Acceleration term ($\mathcal{A}_{\text{acc}}$)**: $(\kappa_0 \kappa_1 \kappa_2)\ddot{\mathbf{u}}$ (no memory variables; supports both consistent and HRZ-lumped mass $\mathbf{M}_{\text{eff}}$).
2. **Velocity term ($\mathcal{A}_{\text{vel}}$)**: $(\kappa_1 \kappa_2 g_0 + \kappa_0 \kappa_2 g_1 + \kappa_0 \kappa_1 g_2)[\dot{\mathbf{u}}]$.
3. **Displacement-inertial term ($\mathcal{A}_{\text{disp}}$)**: $(\kappa_2 (g_1 \circ g_0) + \kappa_1 (g_2 \circ g_0) + \kappa_0 (g_2 \circ g_1) + (\phi_2 \circ g_1 \circ g_0))[\mathbf{u}]$.

### 2.3 Exact Piecewise-Linear Recursive Convolution (`RCFilter1D`) & SPD Tangents
Every elementary first-order rational operator $H(\omega) = c_0 + \frac{c_1}{\beta + i\omega}$ acting on input $x(t)$ over $[t_n, t_{n+1}]$ ($\Delta t = t_{n+1} - t_n > 0$) is integrated exactly assuming piecewise-linear variation of $x(\tau)$:
$$\zeta^{n+1} = b\,\zeta^n + \phi_0\,x^n + \phi_1\,x^{n+1}, \qquad y^{n+1} = c_0\,x^{n+1} + c_1\,\zeta^{n+1}$$
where for $\theta = \beta \Delta t$:
$$b = e^{-\theta}, \qquad \phi_0 = \frac{1 - b(1 + \theta)}{\beta^2 \Delta t}, \qquad \phi_1 = \frac{\theta - 1 + b}{\beta^2 \Delta t}$$
(with 4th-order Taylor expansion for $|\theta| < 10^{-8}$). The instantaneous derivative with respect to the current trial state $x^{n+1}$ is strictly positive:
$$D_H(\Delta t) = \frac{\partial y^{n+1}}{\partial x^{n+1}} = c_0 + c_1 \phi_1(\beta, \Delta t) > 0$$
For a $K$-stage cascade $H = H^{(K)} \circ \dots \circ H^{(1)}$, the chain rule yields $D_{\text{cascade}}(\Delta t) = \prod_{k=1}^K D_{H^{(k)}}(\Delta t) > 0$. Consequently, the exact element tangent matrices:
$$\mathbf{K}_{\text{eff}}(\Delta t) = \frac{\partial \mathbf{R}_{\text{disp}}^{n+1}}{\partial \mathbf{u}_e^{n+1}}, \qquad \mathbf{C}_{\text{eff}}(\Delta t) = \frac{\partial \mathbf{R}_{\text{vel}}^{n+1}}{\partial \dot{\mathbf{u}}_e^{n+1}}, \qquad \mathbf{M}_{\text{eff}} = \frac{\partial \mathbf{f}_{\text{acc}}^{n+1}}{\partial \ddot{\mathbf{u}}_e^{n+1}}$$
are **strictly symmetric positive semi-definite (SPD)** for any $\Delta t > 0$, and the Newmark dynamic tangent $\mathbf{K}_T^{\text{dyn}} = c_1 \mathbf{K}_{\text{eff}} + c_2 \mathbf{C}_{\text{eff}} + c_3 \mathbf{M}_{\text{eff}}$ matches $\frac{d\mathbf{R}_{\text{int}}^{n+1}}{d\mathbf{u}_e^{n+1}}$ to machine precision.

---

## 3. Command Syntax (Tcl & OpenSeesPy)

```tcl
element LadrunoCPML $eleTag {*}$nodeTags \
    (-E $E -nu $nu -rho $rho | -mat $matTag <-rho $rho>) \
    (-box $xmin $xmax $ymin $ymax <$zmin $zmax> <$Lx <$Ly <$Lz>>> \
     | -plane -x0 $x0 $y0 <$z0> -n $nx $ny <$nz> -L $Lpml \
     | -plane $x0 $y0 <$z0> $nx $ny <$nz> <$Lpml>) \
    <-L $L1 <$L2 <$L3>>> <-Lx $Lx> <-Ly $Ly> <-Lz $Lz> \
    <-m $mDeg | -mDeg $mDeg> <-R0 $R0> \
    <-kappaMax $kappaMax | -kappa $kappaMax> \
    <-alphaMax $alphaMax | -alpha $alphaMax | -f0 $f0> \
    <-thick $thickness | -thickness $thickness> \
    <-lumped> <-stage $initialStage> \
    <-shape T3|Q4|T6|BT6|bezier|lagrange|TET4|H8|TET10|BTET10|H20> \
    <-body $bx $by <$bz> | -b $bx $by <$bz>>
```

```python
ops.element(
    "LadrunoCPML", eleTag, *nodeTags,
    "-E", E, "-nu", nu, "-rho", rho,          # or "-mat", matTag (ElasticIsotropic)
    "-box", xmin, xmax, ymin, ymax,           # (+ zmin, zmax in 3D), or "-plane", ...
    "-L", Lpml,                               # or "-Lx", Lx, "-Ly", Ly, ("-Lz", Lz)
    "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 5.0, "-alphaMax", 31.4159,
    "-stage", 1,                              # 0 = static elastic gravity, 1 = dynamic C-PML
    "-b", bx, by,                             # optional body force per unit volume (or -body)
)
```

### Geometry Modes
1. **`-box` Mode (Axis-Aligned Bounding Box)**:
   - Specifies the bounding box $[x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}] (\times [z_{\min}, z_{\max}])$ of the **interior physical domain**.
   - At each Gauss point $\mathbf{x}_g$, the penetration depth along axis $m$ is computed automatically as:
     $$\xi_m(\mathbf{x}_g) = \max\!\big(0,\; x_{m,\min} - x_{g,m},\; x_{g,m} - x_{m,\max}\big)$$
   - Automatically classifies face ($\xi_m > 0$ along 1 axis), edge ($\xi_m > 0$ along 2 axes), and corner ($\xi_m > 0$ along 3 axes) quadrature points without requiring separate face/edge/corner element types. To leave a free surface or symmetry plane unstretched, set its bound beyond the mesh extent (e.g. `ymax = 1e6` for a free surface at $y = 0$).
2. **`-plane` Mode (Arbitrary Planar Layer with Outward Normal $\hat{\mathbf{n}}$)**:
   - Accepts either keyword flags (`-plane -x0 $x0 $y0 <$z0> -n $nx $ny <$nz> -L $Lpml`, with `-normal` alias for `-n`) or inline positional values (`-plane $x0 $y0 <$z0> $nx $ny <$nz> <$Lpml>`).
   - Constructs a local orthonormal frame $\mathbf{R} \in SO(\text{ndm})$ with first axis along $\hat{\mathbf{n}}$ and stretches along $\xi_0(\mathbf{x}_g) = \max(0, \hat{\mathbf{n}} \cdot (\mathbf{x}_g - \mathbf{x}_0))$.

### Static Gravity to Dynamic Stage Transition (`-stage 0 -> 1`) & Body Loads
- Initialize with `-stage 0` during static gravity analysis (`LoadControl`): the element acts as a pure linear-elastic isotropic continuum element with stiffness $\mathbf{K}_0$ and zero C-PML damping.
- Body forces can be applied directly via `-b` / `-body $bx $by <$bz>` or scaled through elemental load patterns (`eleLoad -type -selfWeight` or `eleLoad -type -BrickSelfWeight`).
- After gravity converges, transition to stage `1` via:
  * `ops.parameter(pTag, "element", eleTag, "stage")` + `ops.updateParameter(pTag, 1)` (**recommended** for multi-material soil-structure models), or
  * `ops.updateMaterialStage("-material", eleTag, "-stage", 1)` (filtered by `eleTag` so soil `updateMaterialStage` calls with other material tags never hijack `LadrunoCPML` elements).
- On `0 -> 1` transition, `LadrunoCPML` latches the committed static displacement $\mathbf{u}_0$ and static internal force $\mathbf{f}_0 = \mathbf{K}_0 \mathbf{u}_0$, resets all C-PML recursive-convolution filters, and applies C-PML absorption strictly to the dynamic perturbation $(\mathbf{u} - \mathbf{u}_0, \dot{\mathbf{u}}, \ddot{\mathbf{u}})$, preserving static equilibrium to machine precision ($\le 10^{-13}$) with zero transient jump.

---

## 4. Coupled Biot $u\text{–}p$ Poroelastic Mode (`-up`)

`LadrunoCPML` natively supports saturated porous media in the $u\text{–}p$ (solid displacement $\mathbf{u}$, pore fluid pressure $p$, compression-positive) approximation of Biot poroelasticity (`-up`), pairing directly with [`LadrunoUP`](71_ladruno_up_adr.md) (`ELE_TAG_LadrunoUP = 33017`) across all 9 continuum topologies in 2D and 3D without interface constraints or auxiliary nodal DOFs.

### 4.1 Stretched Biot $u\text{–}p$ Weak Form & Recursive Convolution Cascades
In the local orthonormal frame $\mathbf{x}' = \mathbf{R}(\mathbf{x} - \mathbf{x}_{\text{ref}})$ with complex coordinate stretching $\frac{\partial}{\partial \tilde{x}'_j} = \frac{1}{s_j(\omega)} \frac{\partial}{\partial x'_j}$ and stretched Jacobian $J_s(\omega) = \prod_{m=0}^{\text{ndm}-1} s_m(\omega)$, the frequency-domain Biot $u\text{–}p$ momentum and fluid mass conservation equations are:
$$-\frac{\partial}{\partial \tilde{x}'_j}\!\left(\sigma''_{ij} - \alpha_B \delta_{ij} \hat{p}\right) + \rho (i\omega)^2 \hat{u}_i = \rho \hat{b}_i$$
$$\alpha_B \frac{\partial (i\omega \hat{u}'_j)}{\partial \tilde{x}'_j} + \frac{1}{\bar{Q}}(i\omega \hat{p}) - \frac{\partial}{\partial \tilde{x}'_j}\!\left(k'_{jl} \frac{\partial \hat{p}}{\partial \tilde{x}'_l}\right) = -\frac{\partial}{\partial x'_j}\!\left(k'_{jl} \rho_f \hat{b}'_{f,l}\right)$$
where $\rho$ is the saturated mixture mass density, $\alpha_B \in [n, 1]$ is the Biot-Willis effective stress coefficient, $\frac{1}{\bar{Q}} = \frac{n}{K_f} + \frac{\alpha_B - n}{K_s}$ is the inverse Biot storage modulus (with $\frac{1}{K_s} = 0$ when $K_s = 0$), and $k'_{jl} = \sum_m R_{jm} \bar{k}_m R_{lm}$ is the rotated mobility tensor ($\bar{k}_m = k_{h,m}/\gamma_w$ when hydraulic conductivity $k_{h,m}$ and fluid unit weight $\gamma_w$ are supplied).

Multiplying both equations by $J_s(\omega)$ and integrating by parts over $\Omega_e$ yields four stretched poroelastic operators at each Gauss point $g$:

1. **Undrained Fast P-Wave ($P_1$) Grading Speed**:
   Because the fast compressional wave ($P_1$) governs the fastest characteristic propagation speed in a saturated porous medium, the C-PML damping profile $d_{0,m}$ uses the undrained fast P-wave velocity:
   $$V_{p,\text{fast}} = \sqrt{\frac{\lambda + 2\mu + \alpha_B^2 \bar{Q}}{\rho}}, \qquad d_{0,m} = -\frac{(m_{\text{deg}} + 1)\, V_{p,\text{fast}} \ln(R_0)}{2 L_m}$$

2. **Darcy Permeability Cascade $F_{jl}^{(p)}(\omega)$**:
   The stretched Darcy flux term becomes:
   $$\int_{\Omega_e} \Phi^{(p)}_{b,j} \, k'_{jl} \, F_{jl}^{(p)}(\omega)\!\left[\Phi^{(p)}_{c,l} \hat{p}_c\right] d\Omega, \qquad F_{jl}^{(p)}(\omega) = \frac{J_s(\omega)}{s_j(\omega) s_l(\omega)}$$
   where $\Phi^{(p)}_{b,j} = \sum_m R_{jm} \frac{\partial N_b^{(p)}}{\partial x_m}$ is the local gradient of the pressure shape function $N_b^{(p)}$, and $F_{jl}^{(p)}(\omega)$ has the identical 2-stage (2D) or 3-stage (3D) rational cascade structure as the solid stiffness filter $F_{jl}(\omega)$.

3. **Hydro-Mechanical Coupling Cascades $E_j^{(p)}(\omega)$ and $E_j^{(\dot{u})}(\omega)$ (Exact Integration by Parts)**:
   Define the transverse product filter along local axis $j$:
   $$E_j(\omega) = \frac{J_s(\omega)}{s_j(\omega)} = \prod_{m \neq j} s_m(\omega) = \begin{cases} E_0 = s_1, \;\; E_1 = s_0 & (\text{2D}) \\ E_0 = s_1 \circ s_2, \;\; E_1 = s_0 \circ s_2, \;\; E_2 = s_0 \circ s_1 & (\text{3D}) \end{cases}$$
   Because $s_m(\omega)$ depends **only** on local coordinate $x'_m$, $E_j(\omega)$ is independent of $x'_j$:
   $$\frac{\partial E_j(\omega)}{\partial x'_j} \equiv 0$$
   Therefore, integrating the pore-pressure gradient term $J_s(\omega) \frac{\alpha_B}{s_j(\omega)} \frac{\partial \hat{p}}{\partial x'_j} = \alpha_B \frac{\partial}{\partial x'_j}\!\big(E_j(\omega)\hat{p}\big)$ by parts produces **zero spatial derivative of the stretching factors**, giving the exact weak coupling terms:
   - **On solid momentum rows ($\mathbf{u}$-rows)**:
     $$R_{\text{disp}, a i}^{(u p), n+1} = -\sum_g w_g \alpha_B \sum_{j=0}^{\text{ndm}-1} R_{ji} \tilde{\Phi}_{a,j} \, \tilde{p}_j^{n+1}, \qquad \tilde{p}_j(\omega) = E_j^{(p)}(\omega)[\hat{p}(\mathbf{x}_g)]$$
   - **On fluid continuity rows ($p$-rows)**:
     $$R_{\text{vel}, b}^{(p u), n+1} = \sum_g w_g \alpha_B N_b^{(p)} \sum_{j=0}^{\text{ndm}-1} \tilde{\dot{\varepsilon}}_{jj}^{n+1}, \qquad \tilde{\dot{\varepsilon}}_{jj}(\omega) = E_j^{(\dot{u})}(\omega)\!\left[\sum_a \tilde{\Phi}_{a,j} \dot{u}'_{a,j}\right]$$
   where $\tilde{\Phi}_{a,j} = \Phi'_{a,j}$ for `-formulation std` and $\tilde{\Phi}_{a,j} = \bar{\Phi}'_{a,j}$ (element-mean local gradient) for `-formulation bbar`. Since $E_j^{(p)}(\omega)$ and $E_j^{(\dot{u})}(\omega)$ share the identical instantaneous chain-rule derivative $D_{E,j}(\Delta t) = \prod_{m \neq j} D_{s_m}(\Delta t) > 0$, the discrete coupling tangent block is exact:
   $$\left(\mathbf{Q}_{\text{eff}}(\Delta t)\right)_{a i, b} = \sum_g w_g \alpha_B \left(\sum_{j=0}^{\text{ndm}-1} R_{ji} \tilde{\Phi}_{a,j} D_{E,j}(\Delta t)\right) N_b^{(p)}$$

4. **Fluid Storage Cascade $J_s^{(p)}(\omega)$**:
   The stretched fluid compressibility term $\frac{1}{\bar{Q}} J_s(\omega)(i\omega \hat{p})$ is filtered directly through the forward cascade $J_s^{(p)}(\omega) = s_{\text{ndm}-1} \circ \dots \circ s_0$ applied to the Gauss-point pressure rate $\dot{p}_g = \sum_c N_c^{(p)} \dot{p}_c$:
   $$R_{\text{vel}, b}^{(p p, S), n+1} = \sum_g w_g \frac{1}{\bar{Q}} N_b^{(p)} \, J_s^{(p)}[\dot{p}_g]^{n+1}, \qquad \left(\mathbf{S}_{\text{eff}}(\Delta t)\right)_{bc} = \sum_g w_g \frac{1}{\bar{Q}} D_{J_s}(\Delta t) N_b^{(p)} N_c^{(p)}$$

5. **Unstretched McGann Pressure Stabilization $\tilde{\mathbf{H}}$**:
   For equal-order interpolation (`-pOrder equal`), the McGann et al. (2012, 2015) pressure-Laplacian stabilization matrix:
   $$\tilde{H}_{bc} = \sum_g w_g \alpha_{\text{stab}} \nabla N_b^{(p)} \cdot \nabla N_c^{(p)}, \qquad \alpha_{\text{stab}} = \begin{cases} 0 & \text{if } \texttt{-pOrder linear} \text{ or } \texttt{-stab off} \\ \dfrac{\alpha_0 h^2}{\lambda + 2\mu} & \text{if } \texttt{-stab auto } (\text{default } \alpha_0 = 0.25) \\ \alpha_{\text{user}} & \text{if } \texttt{-stab manual } \alpha_{\text{user}} \end{cases}$$
   is added **unstretched** to the rate matrix $\mathbf{C}_{\text{eff}}$ (and to $\mathbf{R}_{\text{vel}}$ via $\tilde{\mathbf{H}}\dot{\mathbf{p}}_e$), cured of undrained $u\text{–}p$ checkerboarding without introducing artificial C-PML memory variables.

---

### 4.2 OpenSees Honest-$p$ Block Placement & Unsymmetric Solver Requirement
Matching `LadrunoUP` (`ADR-71`) and upstream `FourNodeQuadUP` / `BrickUP`, `LadrunoCPML -up` places the physical pore pressure $p$ (without time-integrating $p$ into a potential) directly into the OpenSees element matrices:
$$\mathbf{K}_{\text{eff}}(\Delta t) = \begin{bmatrix} \mathbf{K}_{uu,\text{eff}}(\Delta t) & -\mathbf{Q}_{\text{eff}}(\Delta t) \\ \mathbf{0} & \mathbf{H}_{\text{eff}}(\Delta t) \end{bmatrix}, \qquad \mathbf{C}_{\text{eff}}(\Delta t) = \begin{bmatrix} \mathbf{C}_{uu,\text{eff}}(\Delta t) + \mathbf{C}_{\text{Ray},uu} & \mathbf{0} \\ \mathbf{Q}_{\text{eff}}(\Delta t)^T & \mathbf{S}_{\text{eff}}(\Delta t) + \tilde{\mathbf{H}} \end{bmatrix}, \qquad \mathbf{M}_{\text{eff}} = \begin{bmatrix} \mathbf{M}_{uu,\text{eff}} & \mathbf{0} \\ \mathbf{0} & \mathbf{0} \end{bmatrix}$$

> [!IMPORTANT]
> **Unsymmetric Linear Solver Required in `-up` Mode**:
> Because $-\mathbf{Q}_{\text{eff}}$ sits in the upper-right $(u,p)$ block of $\mathbf{K}_{\text{eff}}$ while $+\mathbf{Q}_{\text{eff}}^T$ sits in the lower-left $(p,u)$ block of $\mathbf{C}_{\text{eff}}$, the dynamic Newmark tangent $\mathbf{K}_T^{\text{dyn}} = c_1 \mathbf{K}_{\text{eff}} + c_2 \mathbf{C}_{\text{eff}} + c_3 \mathbf{M}_{\text{eff}}$ has upper-right block $-c_1 \mathbf{Q}_{\text{eff}}$ and lower-left block $+c_2 \mathbf{Q}_{\text{eff}}^T$ ($c_1 \neq c_2$). Even in static analysis ($c_2 = c_3 = 0$), $\mathbf{K}_0$ has $-\mathbf{Q}_0$ in the upper right and $\mathbf{0}$ in the lower left. Always use a general unsymmetric solver:
> - `system UmfPack`
> - `system SparseGEN`
> - `system FullGeneral`
> - `system Mumps` (with unsymmetric `SYM=0`)

#### Solid-Only Rayleigh Damping Isolation
When element or domain `rayleigh` damping is active in `-up` mode, `LadrunoCPML` restricts $\mathbf{C}_{\text{Ray}}$ strictly to the translational solid $(u,u)$ block (`getDamp()` and `getRayleighDampingForces()`), never contaminating the fluid continuity rows ($p$-rows) with artificial $\beta_K \mathbf{H}_0 \dot{\mathbf{p}}$ or $-\beta_K \mathbf{Q}_0 \dot{\mathbf{p}}$ forces.

---

### 4.3 Equal-Order (`-pOrder equal`) vs. Taylor–Hood (`-pOrder linear`) Discretizations

| Pressure Mode | Supported Topologies | Nodal DOF Layout | Pressure Stabilization |
|---|---|---|---|
| **`-pOrder equal`** *(default on linear topologies)* | All 9 topologies: `T3, Q4, T6, BT6` (2D); `TET4, H8, TET10, BTET10, H20` (3D) | Every node $a \in \{0, \dots, n_{\text{nodes}}-1\}$ carries $\text{ndm} + 1$ DOFs: $(u_x, u_y, [u_z], p)$. | McGann pressure Laplacian active by default (`-stab auto`, $\alpha_0 = 0.25$). |
| **`-pOrder linear`** *(Taylor–Hood; auto-detected on quadratic topologies if mid-edge nodes have `ndm` DOFs)* | Quadratic topologies (`nNodes > nVertices`): `T6, BT6` (2D, $n_P = 3$); `TET10, BTET10` (3D, $n_P = 4$); `H20` (3D, $n_P = 8$) | **Vertex nodes** $0 \dots n_P-1$ carry $\text{ndm} + 1$ DOFs $(u_x, u_y, [u_z], p)$; **mid-edge nodes** $n_P \dots n_{\text{nodes}}-1$ carry $\text{ndm}$ DOFs $(u_x, u_y, [u_z])$. | Inf-sup (LBB) stable by construction; McGann stabilization is automatically disabled ($\tilde{\mathbf{H}} \equiv \mathbf{0}$). |

---

### 4.4 Command Syntax for Coupled Biot $u\text{–}p$ Mode (`-up`)

```tcl
element LadrunoCPML $eleTag {*}$nodeTags \
    (-E $E -nu $nu -rho $rho | -mat $matTag <-rho $rho>) \
    (-box $xmin $xmax $ymin $ymax <$zmin $zmax> <$Lx <$Ly <$Lz>>> \
     | -plane -x0 $x0 $y0 <$z0> -n $nx $ny <$nz> -L $Lpml \
     | -plane $x0 $y0 <$z0> $nx $ny <$nz> <$Lpml>) \
    -up -Kf $Kf -poro $n -rhoF $rhoF \
    (-perm $kx <$ky <$kz>> | -permH $khx <$khy <$khz>> -gammaW $gammaW) \
    <-alphaB $alphaB> <-Ks $Ks> \
    <-stab off|none|auto <$alpha0>|manual $alpha|$alpha> \
    <-pOrder equal|linear> <-formulation std|bbar> \
    <-body $gx $gy <$gz>> <-fluidBody $fgx $fgy <$fgz>> \
    <-L $L1 <$L2 <$L3>>> <-m $mDeg> <-R0 $R0> <-kappaMax $k> <-alphaMax $a | -f0 $f0> \
    <-thick $t> <-lumped> <-stage 0|1> <-shape ...>
```

| Flag | Description |
|---|---|
| `-up` | Enables coupled Biot $u\text{–}p$ poroelastic CFS-PML mode. |
| `-Kf $Kf` | Pore fluid bulk modulus $K_f > 0$ *(required in `-up`)*. |
| `-poro $n` | Porosity $n \in [0, 1)$ (aliases: `-n_poro`, `-porosity`) *(required in `-up`)*. |
| `-rhoF $rhoF` | Pore fluid mass density $\rho_f \ge 0$ (alias: `-rhof`) *(required in `-up`)*. |
| `-perm $kx <$ky <$kz>>` | Mobility $\bar{k}_m = \kappa_m / \mu_f \ge 0$ ($\text{m}^3\cdot\text{s}/\text{kg}$, isotropic if 1 value, orthotropic if `ndm` values; alias: `-kbar`). Mutually exclusive with `-permH` / `-gammaW`. |
| `-permH $khx <$khy <$khz>> -gammaW $gw` | Hydraulic conductivity $k_{h,m} \ge 0$ ($\text{m}/\text{s}$) and fluid unit weight $\gamma_w > 0$ ($\text{N}/\text{m}^3$), converted internally via $\bar{k}_m = k_{h,m} / \gamma_w$. Both flags must be supplied together. |
| `-alphaB $alphaB` | Biot-Willis effective stress coefficient $\alpha_B \in [n, 1]$ (alias: `-biotAlpha`; default `1.0`). |
| `-Ks $Ks` | Solid grain bulk modulus $K_s > 0$ (optional; default `0.0` $\implies$ incompressible grains $1/K_s = 0$). |
| `-stab ...` | Equal-order McGann stabilization: `off` / `none` ($\alpha_{\text{stab}} = 0$), `auto <$alpha0>` ($\alpha_{\text{stab}} = \alpha_0 h^2 / (\lambda + 2\mu)$, default $\alpha_0 = 0.25$), or `manual $alpha` ($\alpha_{\text{stab}} = \alpha$). Ignored when `-pOrder linear`. |
| `-pOrder equal\|linear` | Pressure interpolation order: `equal` (all nodes carry $\text{ndm}+1$ DOFs) or `linear` / `TH` / `taylor-hood` (vertex nodes carry $\text{ndm}+1$ DOFs, mid-edge nodes carry $\text{ndm}$ DOFs on `T6, BT6, TET10, BTET10, H20`). Auto-detected on quadratic topologies from mid-edge node DOF count if omitted. |
| `-formulation std\|bbar` | Solid volumetric strain formulation: `std` (default) or `bbar` (mean-dilatation $\bar{\mathbf{B}}$ on both solid stiffness and Biot coupling $\mathbf{Q}_{\text{eff}}$; alias: `-form`). |
| `-body $gx $gy <$gz>` | Gravitational acceleration vector $\mathbf{g}$ in `-up` mode: contributes nodal mixture weight $\rho V_a \mathbf{g}$ on solid DOFs (matching `LadrunoUP`). |
| `-fluidBody $fgx $fgy <$fgz>` | Fluid gravitational acceleration vector $\mathbf{g}_f$ for Darcy seepage right-hand side $\int_{\Omega_e} \nabla N_b^{(p)} \cdot \bar{\mathbf{k}} \rho_f \mathbf{g}_f \, d\Omega$ (defaults to `-body` when omitted). |

---

### 4.5 Complete Working Examples (`-up` Mode)

#### Example 1: Saturated Soil Domain (`LadrunoUP`) Surrounded by `LadrunoCPML -up` Absorbing Layer
```python
import opensees as ops

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)  # (ux, uy, p)

E, nu, rho_sat = 5.0e7, 0.30, 2000.0
Kf, poro, rho_f, kbar = 2.2e9, 0.35, 1000.0, 1.0e-7

# ... define nodes in interior [0, 10] x [0, 10] and PML ring of thickness L_pml = 2.0 ...
ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho_sat)

# Interior element (LadrunoUP takes positional matTag = 1 after node tags)
ops.element(
    "LadrunoUP", 1, 1, 2, 5, 4, 1,
    "-Kf", Kf, "-poro", poro, "-rhoF", rho_f, "-perm", kbar,
    "-alpha", 1.0, "-stab", "auto",
)

# Absorbing C-PML element (LadrunoCPML -up)
ops.element(
    "LadrunoCPML", 101, 2, 3, 6, 5,
    "-mat", 1, "-box", 0.0, 10.0, 0.0, 10.0, "-L", 2.0,
    "-up", "-Kf", Kf, "-poro", poro, "-rhoF", rho_f, "-perm", kbar,
    "-alphaB", 1.0, "-stab", "auto", "-stage", 1,
)

# Honest-p Biot tangent is unsymmetric -> use UmfPack / SparseGEN / Mumps(SYM=0)
ops.system("UmfPack")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("Newmark", 0.5, 0.25)
ops.algorithm("Linear")
ops.analysis("Transient")
ops.analyze(200, 0.001)
```

#### Example 2: Staged Hydrostatic Gravity + Steady Seepage Consolidation (`-stage 0 -> 1`)
```python
import opensees as ops

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# 1. Build model with -stage 0 and gravity acceleration -body 0.0 -9.81
#    In stage 0, LadrunoCPML -up acts as an unstretched static Biot u-p element
#    with exact hydrostatic pore-pressure and effective-stress equilibrium.
ops.element(
    "LadrunoCPML", 201, 1, 2, 4, 3,
    "-E", 4.0e7, "-nu", 0.25, "-rho", 2000.0,
    "-box", 0.0, 10.0, 0.0, 10.0, "-L", 2.0,
    "-up", "-Kf", 2.2e9, "-poro", 0.30, "-rhoF", 1000.0,
    "-permH", 1.0e-4, "-gammaW", 9810.0,
    "-body", 0.0, -9.81, "-stage", 0,
)
ops.parameter(9001, "element", 201, "stage")

# 2. Solve static gravity / steady seepage equilibrium
ops.system("UmfPack")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)

# 3. Transition 0 -> 1: latches committed (u0, p0) and static force f0 = K0 * [u0; p0],
#    resets all C-PML convolution filters, and absorbs only dynamic perturbations
ops.updateParameter(9001, 1)
ops.wipeAnalysis()
ops.system("UmfPack")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("Newmark", 0.5, 0.25)
ops.algorithm("Linear")
ops.analysis("Transient")
ops.analyze(100, 0.001)
```

---

## 5. Verification & Benchmark Summary (Single-Phase & Coupled Biot $u\text{–}p$ Warrant Package)

| Verification Layer / Anchor | Test File | Key Metrics & Thresholds |
|---|---|---|
| **P0 Standalone C++17 Kernel & Blind Python Oracle (Single-Phase + Biot `-up`)** | [`tests/test_ladruno_cpml_kernel_cpp.py`](../tests/test_ladruno_cpml_kernel_cpp.py) | 67 C++ self-checks + 18 single-phase + 14 coupled Biot $u\text{–}p$ multi-step cases across all 9 topologies (`T3, Q4, T6, BT6, TET4, H8, TET10, BTET10, H20`, equal-order & Taylor–Hood, `std` & `bbar`): max relative error vs. independent Python oracle = **$6.16 \times 10^{-14}$**; **10/10** internal kernel mutants (`5/5` single-phase + `5/5` Biot `-up`: `UP_NO_DARCY_STRETCH, UP_NO_COUPLING_STRETCH, UP_NO_STORAGE_STRETCH, UP_WRONG_Q_SIGN, UP_DRAINED_VP_GRADING`) killed. |
| **Element-Level Oracle Parity & FD Tangent (Single-Phase)** | [`tests/test_ladruno_cpml_waveguide_mms.py`](../tests/test_ladruno_cpml_waveguide_mms.py) | All 18 cases (`inline_mat` & `ndmat_tag`) through `opensees.pyd` match Python oracle to **$\le 10^{-11}$**; finite-difference tangent check matches $\mathbf{K}_T^{\text{dyn}}$ to **$\le 10^{-7}$**. |
| **Anchor 1: 1D Analytical Waveguide (P, S, `-box`, `-plane`, $37^\circ$ Rotated)** | [`tests/test_ladruno_cpml_waveguide_mms.py`](../tests/test_ladruno_cpml_waveguide_mms.py) | P-wave reflection $R_P = \mathbf{0.414\%} \le 1\%$; S-wave reflection $R_S = \mathbf{0.497\%} \le 1\%$; `-plane` bit-identical to `-box` ($\le 10^{-14}$); $37^\circ$ rotated waveguide matches axis-aligned to **$2.68 \times 10^{-17}$**. |
| **Anchor 2: 2D/3D MMS Patch Completeness & $O(h^2)$ Convergence** | [`tests/test_ladruno_cpml_waveguide_mms.py`](../tests/test_ladruno_cpml_waveguide_mms.py) | Multi-element unstructured patch completeness across all 9 topologies in `-box` corner and `-plane` oblique modes ($\le 10^{-11}$); $O(h^2)$ spatial convergence ($p \ge 1.85$) on unstructured `T3` and `Q4` meshes. |
| **Anchor 3a/3b: 2D Cylindrical & 3D Spherical Full-Space Radiation** | [`tests/test_ladruno_cpml_greens_lamb.py`](../tests/test_ladruno_cpml_greens_lamb.py) | **2D (`Q4` & unstructured `T3`)**: $0.60\%\text{–}1.40\%$ L2 error across $0^\circ, 45^\circ$ corner, and $68^\circ$ grazing receivers (vs. $95.8\%\text{–}173.6\%$ unabsorbed). **3D (`H8` & unstructured `TET4`)**: $0.38\%\text{–}0.84\%$ L2 error across faces, edges, and $(1,1,1)$ 3D corners ($90.6\times\text{–}199.3\times$ reduction). |
| **Anchor 4: Lamb's Problem Rayleigh Wave Head-to-Head vs. `ASDAbsorbingBoundary2D`** | [`tests/test_ladruno_cpml_greens_lamb.py`](../tests/test_ladruno_cpml_greens_lamb.py) | `LadrunoCPML` achieves **$0.60\%$ L2 / $0.62\%$ peak** reflection error vs. **$38.46\%$ L2 / $50.09\%$ peak** for `ASDAbsorbingBoundary2D` (**$64.61\times$ lower L2 error** and **$80.17\times$ lower peak error**). |
| **Anchor 5: Ecosystem, Energy Closure, 2500-Step Stability & Checkpoint** | [`tests/test_ladruno_cpml_ecosystem.py`](../tests/test_ladruno_cpml_ecosystem.py) | `EnergyBalanceRecorder` closure $\le 10^{-10}$; 2500-step late-time Lyapunov stability ($< 10^{-4}$ of peak in 2D/3D corners); `-stage 0 -> 1` gravity latching ($\le 10^{-13}$); `FileDatastore` live & fresh checkpoint round-trip ($\le 10^{-14}$); Tcl parity via `OpenSees.exe`; 26 fail-closed parser checks. |
| **Coupled Biot `-up` Suite (Gate 0 + Anchors 1–3)** | [`tests/test_ladruno_cpml_up.py`](../tests/test_ladruno_cpml_up.py) | **Gate 0**: interior ($\xi_m = 0$) equivalence vs. `LadrunoUP` across all 5 shared topologies (`T3, Q4, H8, BT6, BTET10` $\le 10^{-13}$) + 14-case oracle parity across all 9 topologies ($\le 10^{-10}$) + 4-block FD tangent ($\le 10^{-7}$). **Anchor 1**: 1D saturated soil column fast $P_1$ ($u_x$ & pore pressure $p$ $\le 1\%$) and S-wave ($\le 1\%$) absorption (`-box`, `-plane`, $37^\circ$ rotated invariance $\le 10^{-11}$). **Anchor 2**: 2D (`Q4`) & 3D (`H8`) poroelastic radiation + saturated `Q4` Lamb's half-space Rayleigh + pore-pressure wave absorption ($\le 3.5\%\text{–}4.5\%$ L2, $\ge 10\times$ reduction). **Anchor 3**: hydrostatic gravity + steady seepage `-stage 0 -> 1` latching ($\le 10^{-10}$), solid-only Rayleigh isolation ($\le 10^{-14}$), 2500-step 2D/3D corner stability, `FileDatastore` round-trip across `Q4, BT6, H8, BTET10` ($\le 10^{-14}$), and 28 fail-closed parser/`setDomain` refusals. |
| **ADR-87 D2 Mutation Gate (`CPML` Family)** | [`Ladruno_scripts/mutation_gate.py`](../Ladruno_scripts/mutation_gate.py) | `CPML=ZERO`: **$0.715$** (`103/144` killed, floor `0.600`); `CPML=SCALE`: **$0.632$** (`91/144` killed, floor `0.400`). |

---

## 6. References
1. Matzen, R. (2011). *An efficient finite element time-domain formulation for the elastic second-order wave equation: A non-split complex frequency shifted convolutional perfectly matched layer*. International Journal for Numerical Methods in Engineering, 88(11), 1157–1178.
2. Xie, Z., Komatitsch, D., Martin, R., & Matzen, R. (2014). *Improved forward wave propagation and adjoint-based sensitivity kernel calculations using a numerically stable finite-element PML*. Geophysical Journal International, 198(3), 1714–1747.
3. Roden, J. A., & Gedney, S. D. (2000). *Convolution PML (CPML): An efficient FDTD implementation of the CFS-PML for arbitrary media*. Microwave and Optical Technology Letters, 27(5), 334–339.
4. Zeng, Y. Q., He, J. Q., & Liu, Q. H. (2001). *The application of the perfectly matched layer in numerical modeling of wave propagation in poroelastic media*. Geophysics, 66(4), 1258–1266.
5. Martin, R., Komatitsch, D., & Ezziani, A. (2008). *An unsplit convolutional perfectly matched layer improved at grazing incidence for seismic wave propagation in poroelastic media*. Geophysics, 73(4), T51–T61.
6. McGann, C. R., Arduino, P., & Mackenzie-Helnwein, P. (2012). *Stabilized single-point 4-node quadrilateral element for dynamic analysis of fluid saturated porous media*. Acta Geotechnica, 7(4), 297–311.
