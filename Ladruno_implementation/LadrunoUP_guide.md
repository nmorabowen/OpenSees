# LadrunoUP — user & validation guide (ADR 71, `ELE_TAG` 33017)

`LadrunoUP` is the fork's **Biot u-p saturated-porous continuum element**: one
class that carries solid displacement **u** and pore pressure **p** as coupled
nodal fields, so drained-vs-undrained behaviour *emerges* from the race between
loading rate and permeability instead of being assumed per analysis. It is the
element physics behind consolidation, staged construction on soft ground,
partially-drained seismic response, and liquefaction triggering.

One element covers the whole geometry family through internal, stateless **shape
providers**: `(ndm, nodeCount)` dispatches. v1 ships **T3, Q4, H8** (equal-order
linear p, stabilized) and Bézier **T6, Tet10** (Taylor–Hood: quadratic Bernstein
u, vertex-linear p). Solid formulation is `-formulation std|bbar`, `-geom
linear|corot` (corot = 3D large-rotation, ADR 78 — see §9b).

The load-bearing divergence from upstream OpenSees's eight UP elements is the
**honest pressure DOF**: the nodal DOF *is* p (upstream's is ∫p·dt, recorded as a
*velocity*). Ours records, constrains, and initializes p directly, matches
Abaqus/Kratos semantics, and makes **static** analysis well-posed. The price is an
**unsymmetric** effective tangent — a general solver is mandatory (§2, read it
first). This guide publishes the P1–P4 validation record and the modelling
recipes; it is the companion to the design record `71_ladruno_up_family_adr.md`.

---

## 1. Element command

```
element LadrunoUP $tag $n1 … $nk $matTag \
    <-thick $t>                 ;# 2D only (plane strain implied), default 1.0
    -Kf $Kf -poro $n -rhoF $rhof            ;# REQUIRED
    -perm $k1 $k2 <$k3>         ;# REQUIRED (ndm values), k_hydraulic / gamma_w
    <-permH $k1 $k2 <$k3> -gammaW $gw>      ;# sugar: conductivity / gamma_w internally
    <-alpha $biotAlpha> <-Ks $Ks>           ;# defaults 1.0 / <=0 (infinite grains)
    <-body $b1 $b2 <$b3>>       ;# accelerations, default 0
    <-fluidBody $f1 $f2 <$f3>>  ;# defaults to -body; drives the seepage source
    <-formulation std|bbar> <-pOrder equal|linear> <-lumped> \
    <-stab auto <$alpha0> | off | $alpha>   ;# equal-order default: auto (alpha0=0.25)
    <-dynSeepage on|off>        ;# default OFF (P4 B5 adjudication); on = research opt-in
    <-geom linear|corot>        ;# corot = 3D-only large-rotation u-p (ADR 78); finite reserved
```

`(ndm, k)` selects the provider: **(2,3) T3 · (2,4) Q4 · (2,6) Bézier T6 ·
(3,8) H8 · (3,10) Bézier Tet10**.

- **`$matTag`** — any 2D/3D `nDMaterial`. The element copies it per GP as
  `PlaneStrain` (2D) or `ThreeDimensional` (3D) and the material sees only the
  **effective stress σ′**; p never enters the constitutive update (the
  effective-stress seam). Density rides the material (`getRho` = saturated mixture
  ρ); `-rhoF` is the fluid density used *only* in the seepage body-force term.
- **`-Kf $Kf -poro $n`** — fluid bulk modulus and porosity entered **raw**; the
  element forms storage `1/Q̄ = n/K_f + (α−n)/K_s` itself. (quadUP's `bulk` is the
  pre-combined `Q̄ ≈ K_f/n`; do not pass that here.)
- **`-perm`** — permeability is `k̄ = k_hydraulic / γ_w` [L³·T/M], per axis. If
  you have a conductivity in m/s, use the sugar `-permH $k1 … -gammaW $gw` and the
  parser divides internally (they come together or not at all). This kills the most
  common quadUP user error.
- **`-alpha`/`-Ks`** — Biot–Willis α (default 1.0) and grain bulk modulus
  (default `≤0` ⇒ incompressible grains, α effectively 1). Storage requires
  `0 < n ≤ α ≤ 1` (parser-policed).
- **`-body`/`-fluidBody`** — **always-on** accelerations (upstream behaviour). A
  `SelfWeight` load pattern **replaces** them with loadFactor-scaled values for
  both the solid rows and the seepage source (`zeroLoad` restores the ctor values)
  — this is the staging knob the gravity init recipe drives (§4). All other
  elemental loads are rejected with a named error; surface tractions come via
  standard patterns.

**Legality matrix.** Quadratic shapes (2,6)/(3,10) are **Taylor–Hood only** and
**require `-pOrder linear`** (quadratic u, linear vertex p); omitting it or passing
`-pOrder equal` on a quadratic shape is **fatal** (equal-order quadratic is a
reserved axis with no theory/gate). Linear shapes accept `equal`; `linear` there is
a documented synonym (all nodes are vertices). On equal-order pairs an omitted
`-stab` defaults to `auto` (α₀=0.25) with a one-line notice, so B4's `-stab off`
leg is always an explicit opt-out. `-stab` on a TH pair is **fatal** (TH is inf-sup
stable without it). 2D is **plane strain only**.

> [!warning] Unknown flags are FATAL — a deliberate break from the family
> The rest of the fork warns-and-continues on an unrecognized flag; LadrunoUP
> **refuses**. A mistyped u-p flag (`-Kf` vs `-kf`, `-perm` vs `-permH` without
> `-gammaW`) silently changes the physics, so the parser stops rather than run the
> wrong model. Parser truth: `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp`.

**Responses** (per-GP): `stresses` (effective σ′), `stressesTotal`
(σ_total = σ′ − α·p on the normals), `porePressure` (scalar, `Np·p_nodal`), `flux`
(Darcy `−k̄·(∇p − ρ_f(b_f − ü))`, ndm comps), plus `material $gp $arg` forwarding.
Full response story in §8.

---

## 2. THE SOLVER REQUIREMENT (read before you run anything)

> [!danger] LadrunoUP models MUST name a general solver. The default is wrong, and it fails SILENTLY.
> The honest-p contract makes the effective transient tangent **unsymmetric** —
> `−Q` lives in `getTangentStiff()` (u-rows) and `+Qᵀ` in `getDamp()` (p-rows), so
> `c₁K + c₂C` never has matching off-diagonal pairs. **Symmetric-storage solvers
> (ProfileSPD, SparseSYM, MUMPS SYM=2) store only the upper triangle: one Q block
> is silently discarded at assembly, and the solve returns plausible-looking
> garbage pore pressures.**
>
> `ProfileSPD` is the **no-`system`-command default**. On the identical Terzaghi
> column, `system UmfPack` gives p ∈ [0, 5] kPa (matches the series); `system
> ProfileSPD` gives p ≈ **1e88–1e89** with every `analyze()` still returning 0.
> Nothing fails loudly. No framework hook lets an element reject an SOE.

**Always** name a general solver:

```python
ops.system("UmfPack")     # serial, first choice
# also serial-OK: SuperLU, FullGeneral, BandGeneral
# MPI targets: system("Mumps", "-ICNTL14", ...) with SYM=0 (serial Mumps is not compiled)
```

The parser prints a one-line solver notice at the first element creation (once per
process). Pinned:
`tests/test_ladruno_up_element_analytic.py::test_wrong_solver_divergence_profilespd_vs_umfpack`;
LEDGER_quirks row *"Symmetric-storage solvers silently DROP one u-p coupling
block"*. Upstream comparisons in the batteries are therefore **response-level**
equivalence gates, not matrix-level reduce-to.

---

## 3. DOF and BC semantics (the honest-p contract)

DOF order is per-node-interleaved `[ux uy (uz) p]`; the pressure DOF is slot
**ndm+1**. Consequences you model with directly:

- **Record p** with the classic `recorder Node -dof <ndm+1> disp` — p rides the
  displacement channel (not `vel`, as upstream). `ops.nodeDisp(node, ndm+1)` reads
  it live.
- **Drained boundary** = fix the p slot (`ops.fix(node, …, 1)` on slot ndm+1, or
  `sp` a prescribed head under a load pattern).
- **Impervious boundary** = leave the p slot free (the natural, no-flux BC).

> [!warning] Statics need ≥1 drained node per hydraulically-connected region — but an all-impervious static does NOT fail loudly
> With no p fixed anywhere, H is a pure-Neumann Laplacian (`H·1 = 0`) and the
> block-triangular static tangent `det(K)·det(H)` is structurally singular — the
> first thing a user building a sealed model hits. **The ADR's original "loud
> singularity" expectation is REFUTED (P1 finding):** every serial general solver
> (UmfPack / FullGeneral / BandGeneral / SparseGeneral) factorizes the singular
> system through round-off and returns **rc = 0 with an arbitrary,
> solver-dependent pressure level** — the p-RHS is consistent, so no solver reports
> the rank deficiency. The singularity is real (σ_min/σ_max < 1e-12), but the
> symptom is a *silently arbitrary p level*, not a crash. Provide a drained node.
> Pinned (strict xfail on the loud-failure claim + companion SVD test):
> `tests/test_ladruno_up_element_analytic.py::test_all_impervious_static_*`.

**Static well-posedness** also needs `k̄ > 0`. When it holds, statics are a genuine
new capability upstream cannot offer (drained solid + steady-state seepage
`H·p = f_seep`) — and statically p is *not* a Lagrange multiplier, so **equal-order
pairs cannot checkerboard on the static path** (H̃ lives only in the damp block,
never assembled by a `Static` integrator).

**Mixed-ndf regions** (dry ndf=ndm next to saturated ndf=ndm+1) are legal; tie them
with the **explicit-DOF-list** `equalDOF` on shared u-DOFs — see §6 for the worked
example, and `ndf_and_mixed_models_guide.md` for the full mixed-ndf reference (the
bare `equalDOF` sizes from the builder NDF and mis-sizes against the smaller-ndf
node).

---

## 4. Initialization recipes

There is no bulk `ic` command for p (it triggers no element update). Establishing
an initial pore-pressure field is a **recipe**, and one sequencing rule dominates.

Three routes (gated in `tests/test_ladruno_up_init_recorders.py`):

1. **Static steady-seepage stage → transient continuation.** The genuinely new
   route: solve `H·p = f_seep` under a `Static` integrator with prescribed heads
   (§3), then continue transient. Requires a general solver and ≥1 drained node.
2. **Gravity-stage transient → hydrostatic p, then shaking.** Run a few large-Δt
   Newmark steps (γ=0.6) under self-weight to settle a hydrostatic p profile
   (±1% of `ρ_f·g·z`), then continue into the excitation.
3. **Scripted `setNodeDisp $node $pDof $val -commit`** per node — direct, for
   analytic profiles.

Each route should start its follow-on transient **without a pressure shock** (the
first-step Δp is bounded — a gate in the init battery).

> [!warning] Sequencing trap: displacement-zeroing staged steps ALSO zero nodal p
> Under the honest-p contract, `InitialStateAnalysis off` → `revertToStart`
> (displacement-zeroing tricks) zero the p DOF along with the displacements.
> **Sequence any displacement-zeroing step BEFORE you establish the p field**, or
> you erase it. See `project_initdefgrad_staged` for the staged-analysis interplay.

> [!warning] `InitialStateAnalysis` requires a build with the upstream-`191c67c2d` backport
> On builds predating the backport (see the LEDGER_quirks row),
> `ops.InitialStateAnalysis("on"|"off")` heap-corrupts the process on the NEXT
> `wipe()`/model build (upstream dangling-parameter bug) **and** repeat calls
> silently never flip the `ops_InitialStateAnalysis` flag — on those builds use
> `ops.reset()` for the zeroing step (what
> `tests/test_ladruno_up_init_recorders.py` gates) and avoid ISA entirely with
> the UW soil materials. Fixed builds pass
> `tests/test_initial_state_analysis_lifetime.py`; `ops.reset()` remains the
> recommended honest-p sequencing recipe either way.

> [!warning] Staged prescribed-head `sp` uses Penalty (or Lagrange), never Transformation
> Adding a nonzero pressure `sp` **mid-analysis** under
> `constraints('Transformation')` converges cleanly (rc=0) to a **wrong** steady
> state — measured p ≈ −73 for a top head of +1 on a sealed 4-element column.
> Penalty and Lagrange are correct in both sequences (sp-from-step-1 and
> sp-added-later). Use `constraints('Penalty', 1e12, 1e12)`. This mirrors the
> family's fully-prescribed-rig Transformation trap. Pinned in
> `tests/test_ladruno_up_element_analytic.py`; LEDGER_quirks row.

---

## 5. Analysis guidance

- **Transient integration: γ ≥ 0.6.** Newmark kinematics carry a parasitic
  history root `−(1−γ)/γ` on the p-row, **damped only for γ > ½**. A γ=0.5
  consolidation march measurably diverges (~1e10 by Tv=0.5). Use the ZS84
  production set **γ=0.6, β=0.3025**; the milder **γ=0.51, β=0.2575** set (Dewoolkar
  1996) is less algorithmically damped and better against centrifuge data — B5
  pins both.
- **`-dynSeepage off` for quasi-static consolidation.** The seepage drive's `−ü`
  term is physically right for genuine dynamics but is **pure noise amplification**
  in a quasi-static run: trial accelerations of numerically-damped compressible-wave
  modes are noise, and `f_seep` integrates them. **Measured (ZS84 column, γ=0.6):**
  with `-dynSeepage off` the error *converges* 1.3e-3 → 7e-4 as Δt shrinks
  0.08 → 0.005; with `on` it **grows** 1.8e-2 → 8.7e-1 (smaller Δt is
  worse). And the B5 dynamic column showed `on` misbehaves **in genuine dynamics
  too** (wandering post-front p ≈ 1.7–2.0 vs the exact β = 0.973 plateau;
  unbounded shallow-station growth) — so since P4 the **default is `off`**;
  `-dynSeepage on` is an explicit research opt-in (SWANDYNE parity; its +G
  residual term stays FD-gated). Pinned: `tests/test_ladruno_up_element_analytic.py`
  sweep, `tests/test_ladruno_up_element_b5.py` documentation gate; LEDGER_quirks
  row; ADR §12 log 2026-07-13.
- **`-stab off` for wave propagation.** The auto-α Laplacian (h² on ṗ) targets the
  undrained/checkerboard limit; on the B5 fast-wave column it injects ~10%
  spurious deep-station p ringing (measured p = 1.08 at ξ=40 vs β = 0.973).
  Consolidation/footing runs keep `-stab auto`; **wave-propagation runs opt out
  explicitly** with `-stab off`.
- **Δt guidance.** For dynamics use ≈ CFL/2 of the fast (undrained P) wave; for
  consolidation resolve `T/100` at the target `Tv`. Under u-p, refining Δt does
  *not* refine the p field beyond the mesh — the boundary layer is spatial.
- **Convergence test: NormDispIncr, not NormUnbalance.** Under `Penalty`
  prescribed-head constraints, a `NormUnbalance` floor bites at the `H·p` penalty
  scale (huge residual entries at the constrained DOFs) and the test stalls. Use
  `ops.test("NormDispIncr", tol, maxIter)`. The two-leg upstream equivalence gates
  also run unbalance-based tests *only* against upstream (whose ∫p·dt DOF grows
  unboundedly, so disp-increment tests mis-scale there) — for our element,
  displacement-increment is the right convergence norm.
- **Validity envelope (documented, not enforced).** u-p drops the fluid
  acceleration *relative to* the solid; it covers the complete earthquake frequency
  range for `k′ < 10⁻³ m/s` (Zone II of Zienkiewicz–Chang–Bettess). Full-Biot
  effects intrude only at `k′ ≳ 10⁻³ m/s` *and* very high frequency (out of scope).

---

## 6. Shapes and pressure order

### 6.1 Equal-order linear shapes (T3, Q4, H8)

All nodes carry p; every node is `ndf = ndm+1`; stabilization is on by default
(§7). The parser enforces the builder NDF up front for these (a clear P1-style
error). This is the everyday lane for structured meshes.

**T3 honesty (REWORDED from the ADR — P2 refutation).** Crossed-diagonal T3
inherits CST's undrained locking, **but LESS than full-integration Q4**: the ADR's
"worse than Q4" framing is refuted against the std-Q4 reference. Measured
settlement ratios as ν → 0.49: crossed-T3 / std-Q4 goes 0.998 → **1.262** (the
union-jack centre bubble relieves dilatation; std-Q4 itself locks, dropping to
0.751 vs bbar). The pin's spirit survives only against the **locking-free bbar
reference** (T3x/Q4bbar 0.960 → 0.948 — T3 does lock relative to bbar). Also: the
crossed-T3 spurious pressure mode is **centre-node-localized**, not the corner
checkerboard — gate it with a cell-center metric, not the corner lattice. Pinned:
`tests/test_ladruno_up_element_t3.py`.

### 6.2 Taylor–Hood Bézier shapes (T6, Tet10) — the modeling dance

Quadratic u, linear vertex p, LBB-stable, **no stabilization**. Vertices are
pressure carriers (`ndf = ndm+1`); mid-edge nodes carry displacement only
(`ndf = ndm`). Building one is a **heterogeneous-ndf** dance (verified working):

```python
# 1. vertices under ndf = ndm+1 (pressure carriers)
ops.model("basic", "-ndm", 2, "-ndf", 3)
for t, (x, y) in vertices.items():
    ops.node(t, x, y)
# 2. RE-ISSUE model to drop the builder to ndf = ndm for the mid-edge nodes
#    (this changes the default for SUBSEQUENT nodes; existing nodes keep their ndf)
ops.model("basic", "-ndm", 2, "-ndf", 2)
for t, (x, y) in midedges.items():
    ops.node(t, x, y)
# 3. element: VERTICES first, then MID-EDGES, then matTag, then -pOrder linear
ops.element("LadrunoUP", e, v0, v1, v2, m01, m12, m20, matTag,
            *common, "-pOrder", "linear")
```
*(idiom from `tests/test_ladruno_up_element_th.py`, `_bt6_unit_square`)*

- **Node order.** BT6: `v0, v1, v2, m(v0,v1), m(v1,v2), m(v2,v0)`. BTET10:
  `v0..v3` then edges `(0,1)(1,2)(0,2)(0,3)(2,3)(1,3)`. (Pinned to
  `LadrunoUPShapes.h` / `LadrunoUP.cpp` `setDomain`.)
- **Straight sides are MANDATORY.** The Bézier providers assume an affine map. A
  mid-edge node more than `1e-6·edgeLen` off the true midpoint trips a **loud
  setDomain guard** ("…off the midpoint … STRAIGHT sides … element not
  activated") and deactivates the element (zero stiffness → the solve
  singularizes; `analyze` returns rc≠0). Gated:
  `tests/test_ladruno_up_element_th.py::test_vi_straight_side_guard_loud`.
- **BTET10 winding is accepted both ways.** Under the pinned node↔barycentric map,
  a conventionally right-handed tet evaluates `detJ < 0`; the element folds
  `dv = w·|detJ|` (the *measure* folds — the signed-J gradients stay
  orientation-correct) and prints a one-time informational note. Right-handed and
  left-handed input give identical p/u. Gated `test_v_btet10_winding_acceptance`.
- **`-stab off` is FATAL on TH** (and `-stab` in any form). TH does not need
  stabilization; the parser refuses it. Practical caveat: because of this, a mixed
  mesh of TH and equal-order elements **cannot share one argument list** — the
  equal-order elements want `-stab`, the TH ones reject it. Build the two families
  with separate arg lists.

### 6.3 Mixed-ndf interface (the guide's canonical example)

A dry region (standard `quad`, ndf=2) stacked on a saturated BT6-TH region
(ndf=3), interface nodes **duplicated** and tied on the u-DOFs with the
**explicit-DOF-list** `equalDOF` (⟨FW-F10⟩ — the bare form mis-sizes):

```python
# wet vertices retain; dry corners constrained; tie u1,u2 only (p left uncoupled)
ops.equalDOF(4, 20, 1, 2)      # wet vertex 4  <- dry node 20 (dup of node 4's position)
ops.equalDOF(3, 21, 1, 2)      # wet vertex 3  <- dry node 21
```
*(full topology in `tests/test_ladruno_up_element_th.py::test_iv_equaldof_mixed_ndf_interface`)*

Interface u-continuity is exact; the saturated-region p stays free and finite;
static gravity and a transient smoke both converge. This is the pattern for
frame-on-soil / dry-cover-on-saturated-core models.

---

## 7. Stabilization (equal-order lanes)

Equal-order pairs must satisfy inf-sup or pressure **checkerboards** and volumetric
locking appear at the undrained/incompressible limit. The cure is the pressure-
Laplacian stabilization `S* = S + H̃`, `H̃ = ∫(∇N_p)ᵀ α ∇N_p dV`, in the damp p-p
slot. **On by default** for equal-order shapes (McGann lineage); a one-line notice
prints when `-stab` is omitted so `-stab off` is always an explicit opt-out.

**Auto formula (papers-verified, `-stab auto <α0>`):**

```
alpha = alpha0 * h^2 / (K_s + 4*G_s/3)
```

- `α₀` default **0.25**, supported range **[0.1, 0.5]** (McGann 2012/2015; a
  once-per-process warning prints outside it). `α₀ ≤ 0` de-stabilizes → fatal.
- `K_s, G_s` are the **drained skeleton** moduli, extracted from the material's
  **initial elastic tangent** per element per material — so the value is correct
  after `updateMaterialStage`/plasticity, not frozen at construction (a fix over
  SSP's frozen `Kstab`).
- **`h` = largest EDGE length** (not the diagonal). This definition back-calculates
  McGann's own footing example (α₀ ≈ 0.254 with h = 3 m); pinned across the family
  (Q4/H8 side, simplex largest vertex-pair).

**Calibration lineage (B4).** The McGann checkerboard footing pins **α = 6.8×10⁻⁵**
(α₀ ≈ 0.25, h = 3 m) as the clean value; the element's `auto` computes it to
within measurement (6.686×10⁻⁵, ±2%). A **twin-run identity** — `-stab auto`
vs `-stab <that manual value>` giving bit-identical fields — is the proof the auto
extraction fires. Pinned: `tests/test_ladruno_up_element_equiv.py` (2D),
`tests/test_ladruno_up_element_h8.py` (3D cube-footing analog).

**Checkerboard diagnosis.** The CB index `CB = |⟨p, χ⟩| / (‖p‖·‖χ‖)` with χ the
alternating lattice mode; the Laplacian-roughness variant `CB_lap` is the primary
metric (captures the full oscillation, monotone across the α₀ sweep). Measured
suppression (B4 2D): `CB_lap` off=0.603 → α₀=0.25: 0.126 → α₀=0.5: 0.0946.
`-stab off` checkerboards; `auto` is clean. Equal-order needs this; **TH does not**
(it is inf-sup stable by construction — the same instability appears on
*fully-integrated* Q1-P1 comparators, so it is inf-sup-driven, not a one-point-solid
artifact).

> [!warning] Over-stabilization on slivers
> `h²` scaling makes H̃ vanish under refinement (no over-stabilization on
> well-shaped meshes), **but** a sliver element with one long edge inflates `h` and
> can over-smear its pressure. There are no over-stabilization diagnostics in the
> papers; keep aspect ratios sane, and if a graded/distorted mesh looks
> over-diffused in p, sweep α₀ downward.

---

## 8. Responses and recorders

| Response | Content | Length |
|---|---|---|
| `stresses` | **effective** σ′ (material stress, tension-positive) | per-GP |
| `stressesTotal` | σ_total = **σ′ − α·p** on the normal components (compression-positive p) | per-GP |
| `porePressure` | scalar `Np·p_nodal` at each GP | **per-GP** (nGP values) |
| `flux` | Darcy `−k̄·(∇p − ρ_f(b_f − ü))` | ndm per GP |
| `material $gp $arg` | forwarded to the GP material (family idiom) | material-defined |

`porePressure` is per-GP (nGP values), *not* per-carrier — matters on TH shapes
where nGP ≠ nP (e.g. BT6 reports 3 GP values). Response IDs 1–4 avoid the plane
family's 21 (`stressZZ`); anything else returns a null response.

**Recorder story:**

- **Classic node recorder works unmodified:** `recorder Node -dof <ndm+1> disp`
  reads p off the displacement channel. Gated
  `tests/test_ladruno_up_element_th.py::test_ii_recorder_pressure_and_eleresponse_lengths`.
- **`Ladruno_NodeResults` PressureSource is contract-aware:** for a LadrunoUP node
  it reads the **disp slot** (honest-p); for an upstream UP element it keeps the
  legacy **vel slot**. Flag-free, automatic, mixed models supported. Gated in
  `tests/test_ladruno_up_init_recorders.py` (⟨FW-F4⟩).
- **MPCO stays upstream-only.** The frozen MPCORecorder PRESSURE result keeps the
  vel-assumption (change-budget: our-hooks-only) — do not use MPCO's pressure
  channel for LadrunoUP; use the node recorder or PressureSource.

---

## 9. Validation record (benchmarks B1–B5)

| Gate | What it proves | Test | Tolerance / result |
|---|---|---|---|
| **B1** | ZCB80 periodic-load layer at the undrained limit (Q̄=10⁴ + 10⁹) — TH no-checkerboard, quadratic-u | `test_ladruno_up_element_th_b1.py` | Q̄=10⁴ leg \|p̂\| L2 **1.47%**, \|û\| 0.06%; Q̄=10⁹ leg no-checkerboard + monotone L2 convergence |
| **B2** | ZS84 Example-1 consolidation column + γ=0.5 vs 0.6 oscillation | `test_ladruno_up_element_analytic.py` | ~1e-3 rel vs the 1D series (`-dynSeepage off`) |
| **B3** | Boone–Ingraffea poroelastic column (the only α≠1 gate) | `test_ladruno_up_element_analytic.py` | hard numbers u(0⁺)=0.254 mm, u(∞)=0.079 mm, p(0⁺)=0.410 MPa |
| **B4** | McGann checkerboard footing — stabilization gate | `test_ladruno_up_element_equiv.py` (2D), `_h8.py` (3D) | `CB_lap` off/auto suppression + **α=6.8e-5** twin-run identity (6.686e-5, ±2%) |
| **B5** | Simon–Zienkiewicz–Paul dynamic column (full-Biot closed form) + step-load p-oscillation under both Newmark sets | `test_ladruno_up_element_b5.py` | σ̂ front + π̂ plateau hard gates at stations ξ > 40 |

Plus: two-leg upstream equivalence vs **quadUP** (P1) and **brickUP** (P2) — leg 1
tight (γ=½/β=¼, `-dynSeepage off`, machine-identical, 3.6e-15/1.1e-15) and leg 2
production (γ=0.6, Δt-halving mutual convergence); 1D/3D Terzaghi vs series;
SSPquadUP pressure-block cross-check; PDMY liquefaction column vs quadUP
(`test_ladruno_up_element_pdmy.py`, exercising the `updateMaterialStage` transport +
K₀/H̃ cache-dirty path); MP smoke + serial DB round-trip
(`test_ladruno_up_mp_smoke.py`, `_th.py`).

**Benchmark caveats:**

- **B1 needs graded meshes.** A uniform mesh **cannot** pass B1 — the undrained
  layer has a ~5 cm boundary layer that a graded BT6-TH column resolves and a
  uniform one does not. The graded mesh is load-bearing, not cosmetic.
- **B1 quadratic-u is proven by exactness, not by a rate.** The ADR's "quadratic-u
  convergence rate > 2" pin is REFUTED: u inherits the pressure boundary-layer rate
  (~1.0–1.2) *through coupling*. The replacement (stronger) proof is
  **parabola-exactness** — the TH u-field reproduces an exact P2 field to 1.5e-14.
- **B1 density note.** B1 lists grain + water densities but the gated quantities
  are ρ-insensitive; treating ρ=2000 as a single mixture value is self-consistent.
- **B5 errata.** The B5 oracle deliberately KEEPS four paper errata (eq 43 front
  sign, eq 44 σ′=σ−π, Table I spike branch, Table II vs Fig-3 Q swap). **Do not
  "fix them back"** — see ADR §7.1 for the exact list and the algebraic
  justification. B5's `û(0,τ)` is **measured-first** (⟨UP-5⟩): a 1D u-p oracle
  pre-quantifies the full-Biot-vs-u-p discrepancy; if > 0.5% it demotes û to a
  documented comparison (the slow-P-wave boundary layer ξ ≲ 40 is dropped by
  construction, so hard gates live only at ξ > 40).

---

## 9b. Geometric nonlinearity: `-geom corot` (ADR 78)

`-geom corot` opens large **rotation** + geometric (load) stiffness on the 3D
lanes (H8, Bézier Tet10) while reusing the small-strain material library
unchanged (PDMY01 et al. — the whole point; `-geom finite` stays reserved
because the LogStrain lifting adaptor assumes a linear-elastic inner law).
What it buys and what it does not:

- **Buys:** the wedge rotation a bearing mechanism needs; total-stress
  prestress in the geometric stiffness (pore pressure feeds K_geo); frame-
  consistent u-p coupling with an **incremental storage-coupling residual**
  `Qᵀ·Δu_d/Δt` (ADR 78 §3.3 — a velocity-form coupling would accrue a
  systematic, Q̄-amplified compressive pore pressure from the finite-rotation
  chord; the incremental form makes a rigidly rotating saturated block
  generate *exactly zero* excess pore pressure, gated at solver tolerance);
  seepage gravity drive that follows the rotated skeleton.
- **Does not buy:** large deformational strain (a developed shear band is
  outside corot), porosity/permeability evolution with volume change (deferred
  with `finite`), or 2D (`-geom corot` in 2D is fatal — a planar node cloud
  has no polar rotation).
- Composes with `-formulation std|bbar`, `-pOrder`, `-stab`, `-lumped`.
  `stresses`/`stressesTotal` report in the **corotated material frame**;
  `darcyFlux` is pushed forward to the **global frame**.
- **The `-b` trap when twinning against a dry model:** LadrunoUP `-body` is an
  ACCELERATION (× the material's saturated ρ); `LadrunoBrick -b` is a FORCE
  PER UNIT VOLUME. Convert (`b_brick = ρ_sat·b_UP`, or buoyant γ′ for a
  drained twin) and always check the summed base reaction against the expected
  weight — it catches a mix-up instantly and nothing else does.

Battery: `tests/test_ladruno_up_element_corot.py` (rigid-rotation gates,
rotated-frame hydrostatic seepage, corot→linear consolidation, and the
drained-equivalence gate against a dry `LadrunoBrick -geom corot` twin).

## 10. Limitations and roadmap

- **v1 solid part is `-formulation std|bbar`; `-geom linear|corot` (corot 3D
  lanes, ADR 78).** `-geom finite` u-p is a reserved research phase (material-
  side blocker, ADR 78 §1.1); EAS/URI/hourglass are not u-p axes.
- **Equal-order quadratic is a reserved axis** (no theory/gate) — quadratic shapes
  are Taylor–Hood only.
- **bbar-on-TH is legal but ungated.** The formulation axis was validated at P1/P2
  (linear shapes); bbar mean-dilatation on a quadratic Bézier shape parses and runs
  but has no battery behind it.
- **Hybrid u/p (incompressibility) is reserved mode 2 (P6)** — same Q, pressure
  providers, and DOF plumbing, differing in the p-row law (algebraic compliance, a
  genuine static path) and a deviatoric material contract. Ships under the same ELE
  tag 33017 behind a mode flag.
- **Explicit u-p is research (P7).** The pressure equation is parabolic and the
  p-rows carry zero mass, so the existing explicit stack cannot advance them. The
  near-term explicit answer for undrained problems stays material-level
  (`FluidSolidPorousMaterial`). Own ADR when it opens.
- **Growth axis (P5, by demand):** Q9 and H20 providers add via one provider +
  tests (the ADR-72 `LadrunoBrick20` H20 basis is the natural donor) — no new
  class, tag, or broker rows.
- **PyLiq1/TzLiq1/QzLiq1 interop gap:** those materials read mean effective stress
  via *friend access to FourNodeQuad(UP) specifically* and will not see LadrunoUP;
  keep upstream quads under py/tz/qz interfaces until a public `meanEffectiveStress`
  response lands (small follow-up if demand appears).

## 11. Files

`SRC/element/ladrunoUP/` (`LadrunoUP.{h,cpp}`, `OPS_LadrunoUP.cpp`,
`LadrunoUPKernel.h`, `LadrunoUPShapes.h`) ·
`tests/test_ladruno_up_element_{analytic,equiv,h8,t3,th,th_b1,b5,pdmy}.py`,
`tests/test_ladruno_up_init_recorders.py`, `tests/test_ladruno_up_mp_smoke.py` ·
`Ladruno_implementation/71_ladruno_up_family_adr.md` (design record) ·
`ndf_and_mixed_models_guide.md` (mixed-ndf reference).
