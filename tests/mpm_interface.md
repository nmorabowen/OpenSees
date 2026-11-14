# MPM explicit engine and incremental hypo-kernel entry — the frozen interface (WP-232b)

**Status:** phase 1 — proposed for freeze. Frozen when the orchestrator says so on the program
issue (#1075). After the freeze a change goes through the orchestrator, as an amendment logged at
the end of this file, never as a silent edit.

**Authority:** `Ladruno_implementation/232_ladruno_mpm_general_driver_adr.md` §3, §5.2, §6, §8.
Where this file is more precise than the ADR it **decides**; where it seems to contradict the ADR,
the ADR wins and this file is wrong.

**Who reads it.** The fork C++ worker (232c) codes against this file and the two case files —
`tests/mpm_cases.txt`, `tests/hypo_incremental_cases.txt` — and never reads the numpy oracle
(`tests/mpm_reference.py`, `tests/hypo_incremental_reference.py`). Everything needed to reproduce a
case is therefore here or in the case file. Written without reading anything under `SRC/`.

---

## 1. Conventions

| Item | Decision |
|---|---|
| Units | Any consistent set; the engine never converts. Every case uses SI (m, kg, s, N, Pa). |
| Sign | **Tension positive.** Compression is negative stress and negative normal strain. Gravity is a vector, given with its sign (a downward pull in 2D is `(0, -g)`). |
| Dimensions | `ndm` ∈ {1, 2, 3}. The fork engine runs `ndm` 2 (plane strain) and 3; `ndm` 1 exists in this file only to define the 1D factor of the basis (§4). Every 1D case in `mpm_cases.txt` is written as a 2D or 3D strip. |
| Matrices | Row-major. `F[i][j] = ∂x_i/∂X_j`, stored at `F[ndm*i + j]`. Same for `R`, `L`, `ΔF`. |
| Velocity gradient | `L[i][j] = ∂v_i/∂x_j`. |
| Voigt, 3D | `[xx, yy, zz, xy, yz, zx]`. Strain carries **engineering** shear (γ = 2ε); stress carries tensor components. This is the existing hypo kernel's convention, which the new entry must interoperate with: `tests/hypo_reference.py` (docstring `{00, 11, 22, 01, 12, 20}` and its Voigt helpers) and ADR-79 §3 (confirmed by the orchestrator, 2026-11-13). No other order is admitted. |
| Voigt, plane strain | `[xx, yy, xy]`, same shear convention. ε_zz = ε_yz = ε_zx = 0 by construction; σ_zz is **not carried by the engine** (it is the material's business; no comparison in this file uses it). |
| Plane-strain thickness | Masses and volumes in 2D are **full** values (area × thickness). Thickness never enters the kinematics; every 2D case uses thickness 1. |
| Indices | 0-based everywhere: nodes, cells, particles, components. |
| Time | `t_{n+1} = t_n + Δt`. The step that starts at `t_n` uses the load factor at `t_{n+1}` (ADR §5.2 step 2: `Domain::applyLoad(t + Δt)` runs before P2G). |

## 2. The grid

One uniform Cartesian background grid per engine, rebuilt (zeroed) every step (ADR §3 S1).

- **Geometry:** origin `x0` (the minimum corner, `ndm` floats), one cell size `h > 0` for every
  axis, cell counts `n = (n_x[, n_y[, n_z]])`, each ≥ 1. The grid spans `[x0_d, x0_d + n_d·h]` on
  axis d.
- **Nodes:** `(n_x+1)·(n_y+1)·(n_z+1)` of them. Node `(i, j, k)` sits at `x0 + (i·h, j·h, k·h)`,
  with `0 ≤ i ≤ n_x` etc. Its id is

  ```
  1D: id = i
  2D: id = i + (n_x+1)·j
  3D: id = i + (n_x+1)·(j + (n_y+1)·k)
  ```
  (x fastest, then y, then z).
- **Cells:** cell `(i, j, k)` spans `[x0 + i·h, x0 + (i+1)·h]` per axis, `0 ≤ i < n_x`; id
  `i + n_x·(j + n_y·k)`. A coordinate `x` lies in cell `floor((x − x0)/h)`, clamped to `n − 1` so
  the maximum face belongs to the last cell.
- **Faces** (where grid conditions live, §7): `xmin, xmax, ymin, ymax, zmin, zmax` — the node
  sets `i = 0`, `i = n_x`, `j = 0`, `j = n_y`, `k = 0`, `k = n_z`. Each has an **inward unit
  normal** `n` pointing into the grid: `xmin → +e_x`, `xmax → −e_x`, and so on.
- Interior planes are not grid conditions in v1 (open question Q3).

## 3. Particle state

Structure of arrays over N particles. `nv` = 6 (3D) or 3 (plane strain). Everything marked
*trial* has a committed copy and a trial copy; the trial copy becomes committed only at ADR §5.2
step 7 (§6.8 here).

| Array | Shape, dtype | Units | Meaning |
|---|---|---|---|
| `id` | [N] int64 | — | stable, unique, never reused |
| `body` | [N] int32 | — | body id |
| `mat` | [N] int32 | — | material tag |
| `flags` | [N] uint32 | — | reserved; 0 |
| `x` | [N, ndm] float64 | m | position (trial) |
| `v` | [N, ndm] float64 | m/s | velocity (trial) |
| `m` | [N] float64 | kg | mass, constant |
| `V0` | [N] float64 | m³ | initial volume, > 0, constant |
| `V` | [N] float64 | m³ | current volume = `J·V0` (trial) |
| `lp` | [N, ndm] float64 | m | GIMP half-lengths, constant in v1 (uGIMP, §4) |
| `F` | [N, ndm, ndm] float64 | — | deformation gradient (trial); `F = I` at t = 0 |
| `R` | [N, ndm, ndm] float64 | — | polar rotation of `F` (trial); `I` at t = 0 |
| `eps_hat` | [N, nv] float64 | — | accumulated unrotated strain, engineering shear — the material's **strain feed** (trial) |
| `sig_hat` | [N, nv] float64 | Pa | unrotated (material-frame) stress (trial) |
| `sig` | [N, nv] float64 | Pa | Cauchy stress in the global frame, `R·sig_hat·Rᵀ` (trial) |
| `sig0_hat` | [N, nv] float64 | Pa | initial stress, unrotated frame (= global at t = 0, since `R = I`) |

Derived, never stored: `ρ = m/V`, `J = det F`.

## 4. Basis functions

Two bases, chosen per engine: `gimp` (v1 production) and `linear` (the hat — kept as the
cell-crossing contrast; ADR §3 S4 says it is not production grade).

### 4.1 uGIMP, 1D factor

Uniform GIMP: the particle's characteristic function has constant half-length `lp` (it does **not**
follow `F` — that would be cpGIMP, not v1). Bardenhagen & Kober (2004), eq. for the uniform
1D weight; written out in the piecewise form of Steffen, Kirby & Berzins (2008), §2. With
`r = x_p − x_i` (particle minus node), cell size `h`, half-length `lp`, **0 < lp ≤ h/2**:

| region | `S(r)` | `dS/dr` |
|---|---|---|
| `|r| < lp` | `1 − (r² + lp²)/(2·h·lp)` | `−r/(h·lp)` |
| `lp ≤ |r| < h − lp` | `1 − |r|/h` | `−sign(r)/h` |
| `h − lp ≤ |r| < h + lp` | `(h + lp − |r|)²/(4·h·lp)` | `−sign(r)·(h + lp − |r|)/(2·h·lp)` |
| `|r| ≥ h + lp` | 0 | 0 |

`S` and `dS/dr` are continuous at every breakpoint (so region boundaries need no tie rule). The
gradient of node i's function, evaluated at the particle, is `dS/dr` (since `r = x_p − x_i`, the
derivative with respect to `x_p` is `+dS/dr`).

**Support enumeration.** On axis d, with `c = floor((x_p − x0_d)/h)`, the candidate nodes are
`c−1, c, c+1, c+2`; drop those outside `[0, n_d]`. Candidates with `S = 0` may be kept or dropped
(their weight and gradient are both 0). At most 3 per axis carry weight.

**Admissibility (as amended by Amendment 2).** A particle is admissible iff on every axis its
**centre** lies in the grid, `x0_d ≤ x_p ≤ x0_d + n_d·h`, and `0 < lp_d ≤ h/2`. Support candidates
outside `[0, n_d]` are **dropped**: their weight and gradient are lost. When the whole domain is
inside the grid (`x0_d ≤ x_p − lp_d` and `x_p + lp_d ≤ x0_d + n_d·h`), `Σ_i S_ip = 1` and
`Σ_i ∇S_ip = 0` exactly (partition of unity). A domain that overhangs a face by `δ < lp` loses the
weight `δ²/(4·h·lp)` of the missing node.

**Limitation.** Dropping the support nodes outside the grid is exact only where those nodes would
be constrained anyway: faces with a fixed normal velocity. A particle within `lp` of an
**unconstrained** grid face loses partition of unity, and mass and momentum go to the dropped
nodes. The C++ engine therefore **warns** about such particles, and a grid must extend at least
one cell plus `lp` beyond every free surface.

### 4.2 Linear hat, 1D factor

With `c` the particle's cell (clamped, §2) and `ξ = (x_p − (x0 + c·h))/h ∈ [0, 1]`: node `c` gets
`S = 1 − ξ`, `dS/dx = −1/h`; node `c+1` gets `S = ξ`, `dS/dx = +1/h`. Always exactly these two
nodes per axis, **including** when `ξ = 0` or `ξ = 1` (a zero weight with a non-zero gradient is
kept). `lp` is ignored; admissible iff `x0_d ≤ x_p ≤ x0_d + n_d·h`.

### 4.3 Tensor product (2D, 3D)

`S_ip = Π_d S_d(r_d)`; `∂S_ip/∂x_d = S_d'(r_d) · Π_{e≠d} S_e(r_e)`. Each node index tuple maps to a
node id by §2.

**Shape functions are evaluated once per step at the step-start positions `x^n`** and used for
P2G, G2P and the velocity gradient (ADR §5.2 step 5).

## 5. Seeding rule used by the cases (lattice)

The cases describe bodies, not particle lists, so that a case is reproducible from text. The rule
(ADR §7.4 item 2, "lattice first"):

- `ppc = n` particles per cell per axis. The global sub-lattice has points at
  `x0_d + (a + ½)·h/n`, `a = 0 … n·n_d − 1`, per axis.
- A sub-lattice point belongs to a body when it is inside it: `box` — `lo_d ≤ x_d ≤ hi_d` on every
  axis; `disc` (2D) — `|x − c|² < R²`; `ball` (3D) — same, strict.
- Each particle gets `V0 = (h/n)^ndm` (× thickness 1 in 2D), `m = ρ·V0`, `lp_d = h/(2n)`, `F = R = I`,
  `eps_hat = 0`, the body's initial velocity and initial stress (`sig_hat = sig = sig0_hat`).
- **Ordering and ids:** bodies in file order; inside a body, sub-lattice points in order of
  `a_x` fastest, then `a_y`, then `a_z`. Ids run 0, 1, 2, … across all bodies in that order.
- A case may also list explicit particles (`particle` lines); they come after the lattice bodies,
  with the ids they declare.

## 6. One explicit step (MUSL)

The numbering follows ADR §5.2. All sums run over the support nodes of each particle; `S_ip` and
`∇S_ip` are evaluated at `x_p^n`.

### 6.1 Step 1 — admissibility and the critical time step

1. Every particle must be admissible (§4.1, §4.2) at `x^n`; otherwise the step is **refused**
   (error naming the first offending particle id), and nothing changes.
2. Critical time step, from the committed state:

   ```
   c_p   = sqrt(M_p / ρ_p),   ρ_p = m_p / V_p
   M_p   = the largest normal diagonal entry of the particle's current material tangent
           (Voigt rows xx, yy[, zz]); for isotropic elasticity M = λ + 2μ
   Δt_c  = h / max_p ( c_p + |v_p| )            (|·| the Euclidean norm)
   ```
   The step is refused iff `Δt > s·Δt_c`, with the safety factor `s` (default 0.5, §11).
   Sources: the CFL condition (Courant, Friedrichs & Lewy 1928); the `h/(c + |v|)` form for MPM in
   de Vaucorbeil et al. (2020, §4.6); ADR §5.2 step 1 for "cell size over the largest particle wave
   speed plus particle speed". For elasticity `λ = Eν/((1+ν)(1−2ν))`, `μ = E/(2(1+ν))`, in both
   plane strain and 3D.

### 6.2 Step 2 — load factor

`λ_{n+1} = λ(t_{n+1})` from the host's time series. Case time series: `none` → λ ≡ 1;
`linear T_r` → `λ(t) = min(t/T_r, 1)`. Body force per unit mass: `b = λ_{n+1}·g`.

### 6.3 Step 3 — grid reset and P2G

```
m_i     = Σ_p S_ip m_p
p_i^n   = Σ_p S_ip m_p v_p^n                         (momentum)
fext_i  = Σ_p S_ip m_p b
fint_i  = − Σ_p V_p^n  σ_p^n · ∇S_ip                (σ Cauchy, tensor; (σ·∇S)_a = Σ_b σ_ab ∂S/∂x_b)
f_i     = fint_i + fext_i
```
Sulsky, Chen & Schreyer (1994), eqs. for the discrete momentum balance; the internal force is the
particle quadrature of `−∫ σ : ∇N dV` at the current configuration. In plane strain `σ` is the 2×2
in-plane block of `[xx, yy, xy]`.

### 6.4 Active nodes — the mass cutoff

```
m_ref = max_p m_p                  (over all particles of the engine, fixed for the run)
node i is active  ⇔  m_i > c_m · m_ref        (strict)
```
`c_m` default `1e-10` (§11). A node with `m_i = 0` is never active, even with `c_m = 0`. Inactive
nodes take `v_i^n = v_i^{n+1} = a_i = v_i* = 0` and still appear in the G2P sums (with zero
values). The rule is relative so that it does not depend on units; it is the small-mass guard of
Buzzi, Pedroso & Giacomini (2008) in a stated form.

### 6.5 Step 4 — damping, grid momentum update, grid conditions

For every active node, componentwise (d = each axis):

```
v_i^n       = p_i^n / m_i
f_id       ←  f_id − α_D · |f_id| · sign(v_id^n)          sign(0) = 0
v_i^{n+1}   = v_i^n + Δt · f_i / m_i
apply grid conditions to v_i^{n+1}                          (§7, full operator)
a_i         = ( v_i^{n+1} − v_i^n ) / Δt                    (after the conditions)
```
Local non-viscous damping in the FLAC form — the damping force opposes the velocity with a
magnitude proportional to the unbalanced force (Cundall 1987; Itasca, *FLAC Theory and
Background*, "mechanical damping"); its use in MPM for quasi-static equilibrium: Al-Kafaji (2013,
PhD thesis, Stuttgart). `α_D` ∈ [0, 1), default 0. `f` here is the total unbalanced force `fint +
fext`; damping acts before the update and before the grid conditions.

`a_i` is defined **after** the grid conditions, so a constrained component's FLIP increment is
consistent with its constrained velocity (a fixed component gives `a_id = −v_id^n/Δt`).

### 6.6 Step 5 — particle update (MUSL)

```
v_p^{n+1}   = α_F · ( v_p^n + Δt Σ_i S_ip a_i ) + (1 − α_F) · Σ_i S_ip v_i^{n+1}
x_p^{n+1}   = x_p^n + Δt · Σ_i S_ip v_i^{n+1}

p_i*        = Σ_p S_ip m_p v_p^{n+1}                      (second projection)
v_i*        = p_i* / m_i   (active nodes; 0 otherwise)
apply grid conditions to v_i*                               (§7, second-projection operator)

L_p         = Σ_i v_i* ⊗ ∇S_ip                              (L_ab = Σ_i v*_ia ∂S_ip/∂x_b)
ΔF_p        = I + Δt · L_p
F_p^{n+1}   = ΔF_p · F_p^n
(deps, R1, J1) = hypoIncrementFromF(F_p^n, F_p^{n+1})       (§8; status −1 ⇒ the step fails)
V_p^{n+1}   = J1 · V0_p
eps_hat^{n+1} = eps_hat^n + deps
sig_hat^{n+1} = material(eps_hat^{n+1})                     (§9)
sig^{n+1}     = R1 · sig_hat^{n+1} · R1ᵀ                     (Voigt push, §8.4)
R_p^{n+1}     = R1
```
`α_F` ∈ [0, 1] is the **FLIP fraction**: 1 = pure FLIP (the default — Sulsky, Zhou & Schreyer
1995's MUSL), 0 = pure PIC. The blend form is Zhu & Bridson (2005). MUSL with a second projection
of the updated particle momentum before the strain update: Sulsky, Zhou & Schreyer (1995); Buzzi
et al. (2008); de Vaucorbeil et al. (2020, §4.2, "MUSL"). The position update uses the **grid**
velocity `v^{n+1}`, not the blended particle velocity.

`ops_TheActiveElement` is 0 on every thread during this step (ADR W9; no effect on the oracle).

### 6.7 Step 6 — commit and refusal

(As amended by Amendment 1.)

- **Kernel refusal fails the step.** If any particle's kernel status was −1 in §6.6, the step
  fails **before any material commits**: every particle keeps its committed state from step n, the
  materials' trial states revert, and the engine reports the failure.
- **A material refusal at commit aborts the run.** The particle materials then commit. If any
  refuses (the WP-99 counter), the run is **aborted**, not failed and retried: a refusal cannot be
  undone once sibling points have committed (ADR W17; `Domain::commit` itself aborts the analysis
  in that case).

The oracle's elastic material never refuses.

### 6.8 Step 7 — G2P kinematics take their trial values

`x, v, F, R, V, eps_hat, sig_hat, sig` become the committed state. Then ADR step 8
(`Domain::commit`, recorders) and step 9 (profiler rows), which the oracle does not model.

## 7. Grid conditions

Conditions are attached to grid **faces** (§2) and act on the active nodes of that face.

- **`fix <face> <components>`** — the listed velocity components are set to 0 (`components` is a
  subset of `x y z`). A slip (roller) face is `fix` of the face's normal component only.
- **`wall <face> <μ>`** — a stationary rigid wall with Coulomb friction coefficient `μ ≥ 0`, after
  the grid-based contact of Bardenhagen, Brackbill & Sulsky (2000), against a rigid body. With `n`
  the face's inward normal and `v` the node's trial velocity:

  ```
  v_n = v · n
  if v_n < 0:                          (moving into the wall)
      v_t   = v − v_n n
      |v_t| > 0 ?  v ← v_t − min(μ·|v_n|, |v_t|) · v_t/|v_t|
                :  v ← v_t                                  (= 0 tangential)
  else: unchanged                       (separating or at rest: no force)
  ```
  The tangential velocity is reduced by at most μ times the removed normal velocity — Coulomb's
  law in impulse form; `μ·|v_n|·m_i/Δt` is the friction force at the node.

**Order.** Full operator (on `v^{n+1}`, §6.5): every `wall` in face order `xmin, xmax, ymin, ymax,
zmin, zmax`, then every `fix` (a fixed component wins). Second-projection operator (on `v*`,
§6.6): every `wall` applies its **no-penetration part only** (`v_n < 0 ⇒ v ← v − v_n n`, the
tangential part untouched — friction is applied once per step), then every `fix`.

Domain `fix`, `sp` and load patterns never reach the grid (ADR §5.2 step 2).

## 8. The incremental hypo-kernel entry

### 8.1 Signature (C++; the names are frozen)

```cpp
// 3D. Fn, F1 row-major 3x3; deps Voigt [xx yy zz xy yz zx], engineering shear;
// R1 row-major 3x3; returns 0, or -1 when det(F_mid) <= 0 or det(F1) <= 0
// (outputs then untouched).
int hypoIncrementFromF(const double Fn[9], const double F1[9],
                       double deps[6], double R1[9], double* J1);

// Plane strain. Fn, F1 row-major 2x2; deps [xx yy xy] engineering shear;
// R1 row-major 2x2. Same return codes.
int hypoIncrementFromFPlaneStrain(const double Fn[4], const double F1[4],
                                  double deps[3], double R1[4], double* J1);
```
In `SRC/element/solidTransformation/LadrunoHypoKernel.h`, beside the existing entry (ADR §6).

### 8.2 The mathematics (3D)

Hughes & Winget (1980) midpoint, integrated in the unrotated (Green–Naghdi) frame:

```
F_mid  = ½ (F_n + F_1)                     refuse if det F_mid ≤ 0
G      = (F_1 − F_n) · F_mid⁻¹             (the incremental displacement gradient on the midpoint configuration)
Δε     = ½ (G + Gᵀ)                        (spatial, at the midpoint)
R_mid  = polar(F_mid)                      (F = R U, R proper orthogonal, U symmetric positive definite)
Δε̂     = R_midᵀ · Δε · R_mid               (unrotated frame)
deps   = Voigt(Δε̂) with engineering shear: [Δε̂_00, Δε̂_11, Δε̂_22, 2Δε̂_01, 2Δε̂_12, 2Δε̂_20]
R1     = polar(F_1),   J1 = det F_1        refuse if J1 ≤ 0
```
This is the algebra of the fork's existing nodal-displacement entry, with `F` given directly
instead of built from nodal displacements: `tests/hypo_reference.py` (`hypo_increment`) forms
`G = (u1 − un)ᵀ·dNdX·F_mid⁻¹ = (F_1 − F_n)·F_mid⁻¹`. The polar factor is unique for `det F > 0`;
any algorithm (SVD, Higham iteration) must agree to rounding.

Objectivity check (used by the cases): if `F_1 = Q·F_n` with `Q` a rotation of angle < π, then
`G = 2(Q − I)(Q + I)⁻¹` is skew (the Cayley transform), so `deps = 0` exactly, and
`R1 = Q·polar(F_n)`.

### 8.3 Plane-strain embedding

Embed each 2×2 `F` as the 3×3 `[[F00, F01, 0], [F10, F11, 0], [0, 0, 1]]`, run §8.2, and return
`deps = [deps6[0], deps6[1], deps6[3]]`, `R1` = the upper-left 2×2 block, `J1` = the 3×3
determinant (= the 2×2 one). By the block structure `deps6[2] = deps6[4] = deps6[5] = 0` and `R1`'s
3×3 form is a rotation about z — exactly, not to tolerance.

### 8.4 Stress push

`sig = R1 · sig_hat · R1ᵀ`, as tensors. 3D Voigt: the existing `pushStress6` (row `[xx yy zz xy yz
zx]`, tensor shear). Plane strain: the 2×2 analogue on `[xx, yy, xy]`; exact for in-plane
rotations, because `sig_hat` has no xz or yz part.

## 9. The material update the oracle models

- The material receives the **total** accumulated unrotated strain `eps_hat^{n+1}` (OpenSees
  `NDMaterial` takes total trial strain), not the increment.
- **Oracle material: hypoelastic isotropic** — `sig_hat = sig0_hat + D·eps_hat`, with `D` the
  isotropic elasticity matrix in the Voigt convention of §1 (3D 6×6; plane strain 3×3, rows
  `[λ+2μ, λ, 0; λ, λ+2μ, 0; 0, 0, μ]`). The initial stress is a particle-held offset in the
  unrotated frame. For `ElasticIsotropic` the C++ engine must reproduce exactly this law; how it
  holds `sig0_hat` is its choice (open question Q2 for non-elastic materials).
- The tangent for §6.1 is `D` (so `M = λ + 2μ`).

## 10. Observables (what the cases compare)

All sums run over all particles of the engine unless a case names a body.

| Name | Definition |
|---|---|
| `com` | `Σ m x / Σ m` |
| `com_v` | `Σ m v / Σ m` |
| `momentum` | `Σ m v` |
| `KE` | `½ Σ m |v|²` |
| `SE` | `Σ ½ V · eps_hatᵀ · D · eps_hat` (elastic only; Voigt with engineering shear, current `V`) |
| `E_total` | `KE + SE` (only for cases with no damping, no gravity and no walls) |
| `energy_ratio_end` | `E_total(end) / E_total(0)` |
| `period_comx` | from `com_v_x(t_k)` recorded every step: upward zero crossings `com_v_x(t_k) < 0 ≤ com_v_x(t_{k+1})`, crossing time by linear interpolation between the two samples; the mean of the differences of consecutive crossings (needs ≥ 2 crossings) |
| `selfweight_error ρ0 g l0 y_base` | `max_p |sig_yy,p − σ_exact(Y_p)| / (ρ0·|g|·l0)` with `σ_exact(Y) = −ρ0·|g|·(l0 − (Y − y_base))`, `Y_p` the particle's **initial** y (Charlton, Coombs & Augarde 2017, the 1D column; exact for any material by equilibrium when the column is laterally confined) |
| `cell_crossing_meter ρ0 g l0 y_base` | the same residual as RMS: `sqrt(mean_p (sig_yy,p − σ_exact(Y_p))²) / (ρ0·|g|·l0)` |
| `max_speed` | `max` over every recorded step and every particle of `|v_p|` |
| `com_disp_x <body>` | `com_x(end) − com_x(0)` of the named body |
| `com_v_x <body>` | `com_v_x` of the named body after the last step |
| `translation_error <v0>` | `max_p |x_p(end) − (X_p + v0·t_end)|` (Euclidean), `X_p` the initial position |
| `F_error` | `max_p max_ab |F_ab − δ_ab|` after the last step |
| `max_abs_stress` | `max_p max_k |sig_k|` after the last step |
| `accepted` | every step of the case runs without refusal |
| `particle_fields` | per particle, by id: `x, v, F (row-major), V, sig, sig_hat, eps_hat` at the end |
| `refused` | the engine refuses the first step (§6.1) and changes nothing |
| `dt_crit` | `Δt_c` of §6.1 at the initial state |
| `active_nodes` | sorted ids of the active nodes after P2G of the first step (§6.4) |
| `shape` | for one particle: every support node id with `S`, `∂S/∂x_d` (§4) |

`shape`, `dt_crit` and `active_nodes` are unit-level: the C++ tests need a test hook exposing the
engine's basis, estimate and active set (open question Q4).

## 11. Parameters and defaults

| Parameter | Symbol | Default | Range |
|---|---|---|---|
| basis | — | `gimp` | `gimp`, `linear` |
| FLIP fraction | `α_F` | 1.0 | [0, 1] |
| local damping | `α_D` | 0.0 | [0, 1) |
| mass cutoff | `c_m` | 1e-10 | ≥ 0 |
| Δt safety factor | `s` | 0.5 | (0, 1] |
| gravity | `g` | 0 | vector |
| time series | λ(t) | `none` (≡ 1) | `none`, `linear T_r` |

The default `s = 0.5` sits below `1/√3`, the factor by which `h/c` may overestimate the
stability limit of a lumped-mass trilinear cube (its highest mode is the dilatational
checkerboard, with an effective length near `h/√ndm`; compare the characteristic-length estimates
of Flanagan & Belytschko 1981). GIMP's smoothing raises the limit; phase 2 measures it (open
question Q1).

## 12. Case files

- `tests/mpm_cases.txt` — engine cases; format in its header.
- `tests/hypo_incremental_cases.txt` — kernel cases; format in its header (the `case … endcase`
  layout of `tests/hypo_cases.txt`, so the C++ check can reuse that reader pattern).
- `pending-oracle` marks a value the phase-2 oracle fills in; inputs and tolerances are frozen now.
  If the oracle cannot meet a frozen tolerance against a closed form, phase 2 reports it to the
  orchestrator instead of loosening it.

## 13. Sources

- Bardenhagen, S.G., Brackbill, J.U., Sulsky, D. (2000). The material-point method for granular
  materials. *CMAME* 187, 529–541. — grid-based frictional contact (§7).
- Bardenhagen, S.G., Kober, E.M. (2004). The generalized interpolation material point method.
  *CMES* 5(6), 477–495. — GIMP (§4.1).
- Buzzi, O., Pedroso, D.M., Giacomini, A. (2008). Caveats on the implementation of the generalized
  material point method. *CMES* 31(2), 85–106. — MUSL, small nodal masses (§6.4, §6.6).
- Charlton, T.J., Coombs, W.M., Augarde, C.E. (2017). iGIMP: an implicit generalised interpolation
  material point method for large deformations. *Computers & Structures* 190, 108–125. — the
  self-weight column (§10).
- Courant, R., Friedrichs, K., Lewy, H. (1928). Über die partiellen Differenzengleichungen der
  mathematischen Physik. *Math. Ann.* 100, 32–74. — CFL (§6.1).
- Cundall, P.A. (1987). Distinct element models of rock and soil structure. In *Analytical and
  Computational Methods in Engineering Rock Mechanics*, Allen & Unwin, 129–163; Itasca, *FLAC —
  Theory and Background*. — local non-viscous damping (§6.5).
- Al-Kafaji, I.K.J. (2013). *Formulation of a dynamic material point method (MPM) for geomechanical
  problems*. PhD thesis, Universität Stuttgart. — local damping in MPM (paper only; Anura3D source
  is barred, ADR §2).
- de Vaucorbeil, A., Nguyen, V.P., Sinaie, S., Wu, J.Y. (2020). Material point method after 25
  years: theory, implementation, and applications. *Adv. Appl. Mech.* 53, 185–398. — the MUSL
  algorithm, the time-step estimate (§6).
- Flanagan, D.P., Belytschko, T. (1981). A uniform strain hexahedron and quadrilateral with
  orthogonal hourglass control. *IJNME* 17, 679–706. — characteristic length (§11).
- Hughes, T.J.R., Winget, J. (1980). Finite rotation effects in numerical integration of rate
  constitutive equations arising in large-deformation analysis. *IJNME* 15, 1862–1867. — the
  midpoint increment (§8.2).
- Steffen, M., Kirby, R.M., Berzins, M. (2008). Analysis and reduction of quadrature errors in the
  material point method. *IJNME* 76, 922–948. — the uGIMP piecewise formulas (§4.1).
- Sulsky, D., Chen, Z., Schreyer, H.L. (1994). A particle method for history-dependent materials.
  *CMAME* 118, 179–196. — MPM, P2G, the colliding discs (§6.3, cases).
- Sulsky, D., Zhou, S.-J., Schreyer, H.L. (1995). Application of a particle-in-cell method to
  solid mechanics. *Comput. Phys. Commun.* 87, 236–252. — MUSL with FLIP (§6.6).
- Zhu, Y., Bridson, R. (2005). Animating sand as a fluid. *ACM TOG* 24(3), 965–972. — the PIC/FLIP
  blend (§6.6).

## 14. Open questions (for the orchestrator)

- **Q1 — the 2D/3D stability limit.** §6.1 follows the ADR (`h` over wave speed). For lumped-mass
  linear cells the 2D/3D limit is nearer `h/(c√ndm)`; GIMP raises it. Phase 2 measures the
  spectral limit of the oracle on the bar and on a 2D/3D block; if `s = 0.5` is not safe, the
  default changes by amendment.
- **Q2 — initial stress for non-elastic materials.** §9's additive `sig0_hat` is exact for
  elasticity only. For `LadrunoJ2` and `DruckerPrager` the initial stress has to enter the
  material's own state (yield must see it); that is 232c's choice with 232a's Contract 1
  (initial stress field), and no case in this file depends on it.
- **Q3 — interior planes and moving walls.** v1 grid conditions live on grid faces only, and walls
  are stationary. Rigid footings (phase 1) will need interior or moving planes; a later amendment.
- **Q4 — the unit-level test hook.** `shape`, `dt_crit` and `active_nodes` need the engine to expose
  its basis, estimate and active set to a test (e.g. an `mpm debug …` sub-command or a Python-only
  hook). `dt_crit` maps onto `mpm criticalTimeStep`; the other two are 232c's design.

## Amendments

The interface was frozen at `8383758d5` (declared on #1075, 2026-11-13). Every change after that
is listed here, with a date.

### Amendment 1 — 2026-11-13 — refusals (§6.7)

Source: the fork ADR's W17, from a review of #1074. A material refusal cannot be undone once
sibling points have committed; `Domain::commit` itself aborts the analysis ("INCONSISTENT …
aborted").

- A kernel status of −1 in §6.6 fails the step **before** any material commits: every particle
  keeps its committed state, and the materials' trial states revert.
- A material refusal **at commit** (the WP-99 counter) **aborts the run**. It does not fail and
  retry the step.
- No case changes: the oracle's elastic material never refuses.

### Rulings on the open questions — 2026-11-13 (orchestrator, on #1075)

- **Q1:** `s = 0.5` stays. Phase 2 measures the 2D/3D spectral limit; any change comes by
  amendment.
- **Q2:** initial stress for non-elastic materials is 232c's design; the elastic offset of §9 stays.
- **Q3:** v1 has grid-face conditions and stationary walls only.
- **Q4:** 232c provides the test hooks.

### Amendment 2 — 2026-11-13 — admissibility (§4.1)

Found by the phase-2 oracle. The frozen rule required each particle's whole GIMP domain inside
the grid. Every body that rests on a face condition then gets refused within a few steps: a bar
fixed at `x = 0` compresses toward the face, and a self-weight column compresses into its base,
so their domains cross the face. Since conditions live on grid faces (Q3), that is the normal
case, not an edge case.

- **New rule:** a particle is admissible iff its centre lies in the grid and `0 < lp ≤ h/2`.
  Support nodes outside the grid are dropped, with their weight and gradient.
- Partition of unity still holds exactly for any particle whose whole domain is inside.
- The linear basis was already centre-based; unchanged.
- No case's frozen value changes. The basis cases keep their domains inside the grid.
- Accepted by the orchestrator (2026-11-13), with the limitation now stated in §4.1: dropping is
  exact only at faces with a fixed normal velocity. Near an unconstrained face it loses partition
  of unity, so the C++ engine warns, and grids extend one cell plus `lp` beyond every free
  surface.

### Amendment 3 — 2026-11-13 — case corrections (`tests/mpm_cases.txt`)

Found by the phase-2 oracle. Each is a fault in the case, not in the interface.

- **`incline_slide_2d`, `incline_stick_2d`:** gravity is now ramped, `series linear 0.2`.
  - Applied suddenly, gravity makes the elastic block bounce on the wall. The normal force then
    drops to zero for part of each cycle, and the "stick" block crept 0.13 m.
  - Ramped, the stick block moves 9e-5 m (tolerance 0.005 m).
  - The slide's closed form becomes `a·(T_r²/6 + T_r·(T − T_r)/2 + (T − T_r)²/2)` = 0.97678 m at
    `T = 1 s`, `T_r = 0.2 s`. The oracle gives 0.97229 (0.46%); tolerance 3%.
- **`bar_dt_limit_2d`:** the energy bound goes from `le 1.02` to `le 1.05`.
  - At 0.989·Δt_c the measured energy oscillates within ±3.3%, with the same maximum in each third
    of a 6000-step run. It is bounded, not growing.
  - 1.02 was a guess, not a closed form.
- **`particle_fields` layout:** a filled case reads `expect particle_fields eq fields tol …`,
  followed by one `field` line per particle.
- **`k0_quiet_start_2d`:** the oracle's `max_speed` is exactly 0. GIMP equilibrium is exact for a
  linear stress field, and here every number is dyadic. The `eq` line's tolerance goes from
  `0 1e-6` (relative to 0, so no tolerance at all) to `1e-12 0` (absolute), so the C++ engine's
  summation order cannot fail it. The `le 0.006` bound stays.

### Phase-2 measurements — 2026-11-13

- **Q1 (stability limit), measured on the oracle.** The largest stable Δt, as a multiple of the
  §6.1 estimate:
  - 1D GIMP bar: 1.74
  - 2D free block (ppc 2, ν = 0.3; GIMP and linear alike): 1.66
  - 3D free block: 1.31

  The default `s = 0.5` is safe by a factor of at least 2.6. No change.
- **`column_selfweight_2d_gimp` misses its frozen tolerance.** The oracle gives
  `selfweight_error` = 0.072 against the frozen 0.03. Refining does not help: 0.067 at h = 0.25, and
  0.106 at h = 0.125 over the same duration, which is not yet quiescent. The error is explicit-GIMP
  cell-crossing noise at 20% compression.
  - The linear basis gives 0.47.
  - The small-deformation column (`ρ0·g·l0/E` = 0.002) gives 0.013, which is O(lp/l0).
  - The tolerance stays frozen at 0.03. The pytest marks this check as a strict xfail until the
    orchestrator rules: change the gate, or the case. (Ruled: Amendment 4.)

### Amendment 4 — 2026-11-13 — the self-weight column is split

The orchestrator's ruling on the shortfall above. The error grows under refinement
(0.067 → 0.106), so it is a limitation of the basis, not noise. The tolerance is **not** loosened,
and the gate is **not** moved to the RMS meter.

- **New case `column_selfweight_small_2d_gimp`:** the same column with `E = 5e7`
  (`ρ0·g·l0/E` = 0.002), `dt = 0.001`, 6000 steps. Its `selfweight_error` is gated at the frozen
  0.03; the oracle gives 0.0128. **This is the slice-1 gate for the closed form.**
- **`column_selfweight_2d_gimp` (20% compression) stays a strict xfail** at 0.03. Measured
  numbers, also in the case file's comments: max error 0.072 (h = 0.5), 0.067 (h = 0.25),
  0.106 (h = 0.125); RMS meter 0.021. It is the acceptance test for phase 3's domain-tracking
  GIMP or B-spline basis.
- **The RMS `cell_crossing_meter` is an observable, not a gate.** Its `eq` lines are cross-checks
  of the C++ engine against the oracle's number. Only the GIMP < linear relation is asserted on it.
- Quirk ledger: `Ladruno_implementation/ledger/quirks/WP-232b-ugimp-fixed-domains-do-not-converge-under-large-compression.md`.
