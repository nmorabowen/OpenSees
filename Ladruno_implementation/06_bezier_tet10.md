---
title: Bézier tetrahedral element (BezierTet10)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - element
  - adr
---

# Bézier tetrahedral element (BezierTet10)

> **ADR** (Architecture Decision Record). The 3D sibling of [[04_bezier_elements]]
> (`BezierTri6`), deferred there under **D10** ("Tet10 gets its own ADR, same
> recipe"). Captures the v1 scope for porting the Kadapa (2018) quadratic Bézier
> **tetrahedron** (§5, Eq. 57) into the `ladruno` fork. The Tri6 ADR is the
> template — this doc records only what **changes** for 3D and inherits the rest.
> Companions: [[mpco_ladruno_element_contract]] (recorder side — already specifies
> the tet case), [[Ladruno_explicit_roadmap]] (the explicit integrator this element
> feeds), [[project_bezier_element]] (ordering decision + repo layout).

## What

`BezierTet10` — a 10-node quadratic **Bézier** tetrahedron (Bernstein basis) for
3D linear elastostatics and elastodynamics, after Kadapa (2018, IJNME 117:543 §5).
It reuses standard quadratic-Lagrange (T10) tet meshes via the same per-edge
control-point map (Eq. 11–13 applied edge-wise), delegates the constitutive law to
any 3D `NDMaterial`, and offers the optional consistent **B-bar** formulation
(Eq. 45, full 6-component 3D form). Its defining property is the **all-positive
lumped mass** `ρVe/10·diag[1₁₀,1₁₀,1₁₀]` (Eq. 57), which makes quadratic tets
usable in explicit dynamics.

**In scope (v1):**
- Pure-displacement and B-bar formulations, full 3D (3 DOF/node, 6 stress comps).
- Lumped (default) and consistent (`-cMass`) mass.
- **Straight-sided elements only** (mid-edge nodes at edge midpoints ⇒ control
  points = nodes, P = X).
- Self-weight / body-force loads via the standard load-pattern path.
- Recorder self-declaration via the element contract (`basisInfo`, `family=bernstein`,
  `topology=tet`, `paramDomain=bary`, `rational=0`, `numCtrl=10`).

**Out of scope (deferred — inherits Tri6's deferrals):**
- **Curved** Bézier faces/edges and the Eq. 14 Dirichlet control-point↔DOF mapping.
- Rational/NURBS tets, higher order.
- A proper surface traction load (the `-pressure` volume hack carries over only if ported).
- Selectable quadrature (fixed 4-point is exact for straight-sided).
- The explicit time integrator (Chung–Lee β=13/12) — element-side is ready;
  the integrator belongs to [[Ladruno_explicit_roadmap]].

## Why

- **Explicit dynamics in 3D need a non-negative lumped mass.** This is the paper's
  *primary* 3D motivation: standard quadratic Lagrange tets have negative
  mass-matrix entries → unusable lumped diagonals → no explicit integration.
  Bernstein non-negativity gives each node exactly `ρVe/10` (Eq. 57).
- **Near-incompressible locking in 3D.** The 3D B-bar (Eq. 45) gives optimal
  convergence as ν→0.5 — the headline 3D example is the thick-walled sphere under
  internal pressure (§6.2), the natural validation case here. This is where the
  element earns its keep: linear tets (TET4) lock badly; quadratic-Bézier B-bar
  (TETB10B) stays optimal independent of ν.
- **Uses existing mesh generators.** apeGmsh/Gmsh already emit quadratic-Lagrange
  tets; the edge-wise map turns them into Bézier elements at negligible cost.
- **The Tri6 recipe is proven.** v1 Tri6 is merged (`ELE_TAG 33000`, PR #6), validated
  against the paper and on real apeGmsh meshes. Tet10 is the same machine in 3D.

## Where

New code (all under **`SRC/element/bezierTetrahedron/`** in the `ladruno` fork —
a new sibling dir, paralleling `bezierTriangle/`):

| File | Role |
|------|------|
| `BezierTet10.{h,cpp}` | the element |
| `OPS_BezierTet10.cpp` | Tcl/Python factory (`element BezierTet10 …`) |
| `CMakeLists.txt` | subdir build |

Modify (fork-local):
- `SRC/classTags.h` — add `ELE_TAG_BezierTet10` (shipped as `33001` in the ladruno
  private band ≥33000, next after Tri6's `33000`; the early sub-300 candidate `273` was dropped).
- `SRC/interpreter/OpenSeesElementCommands.cpp` — dispatch `BezierTet10`/`bezierTet10`.
- `SRC/element/CMakeLists.txt` — `add_subdirectory(bezierTetrahedron)`.

Reference (copy idioms):
- `SRC/element/bezierTriangle/BezierTri6.{h,cpp}` — **the primary template** (static
  return Matrix/Vector pattern, control-point map, B-bar, lumped mass, `basisInfo`
  response trio, SelfWeight path, straight-sided guard).
- `SRC/element/tetrahedron/TenNodeTetrahedron.cpp` — **the node-ordering + Voigt
  oracle** (see D2): its `shp3d` defines the canonical OpenSees T10 edge order; its
  `computeB` defines the 3D Voigt convention `{xx,yy,zz,xy,yz,zx}` — *identical to
  Kadapa Eq. 26–27*, so the B-bar Eq. 45 maps in directly.

No new external dependency; no compilation-journal entry needed.

## Decisions

D1, D3, D4, D6, D7, D8, D11 **inherit verbatim from [[04_bezier_elements]]**
(Kadapa linear map not IGA; constitutive via `NDMaterial`; lumped mass default +
`-cMass`; SelfWeight load path; recorder via the contract; home = the fork; B-bar
is a `-bbar` flag not a separate element). The decisions that **change for 3D**:

| # | Decision | Rationale | Consequence |
|---|----------|-----------|-------------|
| **D2′** | **Node order = OpenSees `TenNodeTetrahedron` order**, byte-identical: 4 vertices, then mid-edges N5=(1‑2) N6=(2‑3) N7=(1‑3) N8=(1‑4) **N9=(3‑4) N10=(2‑4)** (the jaabell/Larenas `shp3d` order, deliberately N9↔N10-swapped to match Gmsh). | Must be interchangeable with existing T10 meshes + Gmsh output, exactly as D2 demands for T6. **Verified from `shp3d` derivatives**; still to confirm Gmsh native tet10 == this order (O11). | Recorder reuses the Lagrange T10 node order, swaps only the basis to Bernstein. Bernstein N: vertices `Nᵢ=ξᵢ²`, edges `N=2ξₐξ_b` in the same edge order. |
| **D5′** | **B-bar is the 3D headline; no plane-stress guard.** Full 6-row Eq. 45 split across all three normal strains. | Tet10 is pure 3D — there is no plane-stress degeneracy to guard (the D5 guard was 2D-only). B-bar is the *reason* for the 3D element (§6.2). | `-bbar` always valid. Plane-strain/axisymmetric variants of Eq. 45 are N/A here. |
| **D9′** | **v1 = straight-sided only** (P = X), same as D9 — but the guard checks **all 6 mid-edge** nodes against their edge midpoints. | 4-pt tet quadrature is exact, BCs interpolatory, map is trivial when P=X. | `setDomain` warns/rejects if any of the 6 mid-edge nodes deviates from its edge midpoint. Curved faces + Eq. 14 = follow-up. |
| **D10′** | **Quadrature = 4-point tet rule** (degree-2 exact) for stiffness, internal force, **and** the B-bar volume average (Eq. 46). Consistent mass uses a higher-order rule. | Kadapa §6: "four-point rule for tetrahedral elements." Remark 3: a 1-point rule for `ε̄_vol` drops convergence one order — use the full 4-pt. | Same `GP_PARAM` (4×3 barycentric) / `GP_WEIGHT` (4× `1/24`, summing to `Ve_ref=1/6`) emitted to the recorder. Consistent mass needs ~degree-4 (e.g. an 11- or 15-pt rule) for exactness. |

**Element constants (vs Tri6):** `NEN=10`, `NDOF=3`, `NELD=30`, `NSTRESS=6`,
`NGAUSS=4` (stiffness/force/B-bar avg), `NGAUSS_MASS=` higher-order (consistent
mass). Static returns: `K_return(30,30)`, `M_return(30,30)`, `P_return(30)`.

## How (design sketch)

Mechanically identical to `BezierTri6`, lifted to 3D. Per-GP loop:

1. **Shape funcs / derivatives** — quadratic Bernstein on barycentric (ξ₁,ξ₂,ξ₃,ξ₄),
   ξ₄=1−ξ₁−ξ₂−ξ₃. Vertices `Nᵢ=ξᵢ²` (i=1..4), edges `N=2ξₐξ_b` in the D2′ edge order.
   Derivatives w.r.t. the 3 free coords (chain through ξ₄). Cross-check against the
   `TenNodeTetrahedron::shp3d` analytic derivatives (Lagrange→Bernstein differs but
   the node→edge map must match).
2. **Jacobian** 3×3 from control points (P=X for straight-sided) → `detJ`, `J⁻¹`,
   physical derivatives `dN/dx,dN/dy,dN/dz`.
3. **B matrix** 6×30, submatrix `Bₐ` = Kadapa Eq. 38 (rows `{xx,yy,zz,xy,yz,zx}`).
4. **B-bar** (if `-bbar`): volume-averaged `B̄₁,B̄₂,B̄₃` (Eq. 46, over the 4 GPs) →
   `B̄ₐ` = Eq. 45; the three normal-strain rows get the `(/3)` deviatoric split, the
   three shear rows are unchanged.
5. **K / Fint / mass** — `K=Σ wᵢ|Jᵢ| Bᵀ D B`, `Fint=Σ wᵢ|Jᵢ| Bᵀσ`, lumped
   `m=ρVe/10` per DOF (Eq. 57) or consistent `∫ρNᵀN`.
6. **Control-point map** — `computeControlPoints()` applies Eq. 11/12 per edge:
   `P_edge = 2[X_edge − 0.25 Xₐ − 0.25 X_b]` for the two endpoints (a,b) of each of
   the 6 edges; vertices `P=X`. Straight-sided ⇒ P=X (the guard enforces it).

### Recorder side (the contract — **zero recorder code**)

Per [[mpco_ladruno_element_contract]], `BezierTet10` is recorded by the
`MPCO_Ladruno` recorder with **no recorder-side edits** — it self-declares via three
`setResponse` hooks, copied from Tri6 and re-dimensioned. The contract's *Ordering
convention* section already specifies the tet case explicitly:

- `"basisInfo"` → emits an `ElementBasis` tag, attrs: `topology="tet"`,
  `family="bernstein"`, `paramDomain="bary"`, `rational=0`, `numCtrl=10`,
  `numGP=4`, `orderU=orderV=orderW=2` (simplex order). Returns a non-null sentinel.
- `"integrationPoints"` → `Matrix(4,3)` — the 4 GP barycentric coords (store 3 of
  the 4 area coords, `PARAM_DOMAIN="bary"`).
- `"integrationWeights"` → `Vector(4)` — the 4 Gauss weights.
- `rational=0` ⇒ **`controlPointWeights` / `CTRL_WEIGHT` omitted entirely** (same as Tri6).

Result tree (`stress`/`strain` at the 4 GPs → 6 comps each = 24 values; `gaussCoord`
→ 4×3 = 12; `force` → 30; material delegation per GP) mirrors Tri6's `getResponse`.
**The recorder reconstructs `x(ξ)=Σ Nᵢ(ξ)Xᵢ` from BASIS+QUADRATURE universally** — no
per-class table touched. This is the inverted-dispatch payoff: the element ships, the
recorder is untouched.

### Test bed (mirrors Tri6's three tracks)

**Track C — Python oracle (build-free, write first).** There is **no Tet10 Python
reference yet** (Tri6 had `bezier_tri6_element.py`; Tet10 must be authored). Plan:
`bezierFEM/python examples/bezier_tet10_element.py` with the same self-checks Tri6
passed to machine precision — partition of unity, Bernstein non-negativity, exact
volume `∫NₐdΩ=Ve/10`, consistent mass = analytical, equal lumped `ρVe/10`,
constant-stress patch test, 6 rigid-body null-space modes. Then a
`validate_thick_sphere.py` (the paper's §6.2 3D example) reproducing the radial
displacement / `σ_rr`,`σ_θθ` and the **TET4-vs-TETB10-vs-TETB10B locking study at
ν→0.5** — the 3D analogue of `validate_cooks_membrane.py`. Validates the
*formulation* before any C++ is built.

**L0 — C++ smoke (`test_bezier_tet10_smoke.py`).** After fork registration + build:
patch test exact (σ=D·ε at all 4 GPs under a linear displacement field), `gaussPoint`
response, **lumped-mass positivity** (all 30 diagonals > 0 — the whole point), `-bbar`
path runs and is symmetric/PD. Single-element, no mesh generator. Models Tri6's
`test_bezier_smoke.py` (which was 4/4).

**Integration — apeGmsh-driven (deferred until it compiles).** Gmsh emits a
straight-sided T10 tet mesh → extract nodes + connectivity → drive `BezierTet10`
directly (apeGmsh py3.11, our `opensees.pyd` py3.12 — the two-env direct-drive path
proven for Tri6). **First gate: confirm Gmsh native tet10 node order == D2′ order**
(verify max mid-node deviation ~1e-15, as Tri6 did). Benchmark: thick-walled sphere
(§6.2) or a 3D cantilever vs beam theory. The `apeSees` typed registry won't include
`BezierTet10` (same O9 gap as Tri6) — direct-drive is the path.

**Regression.** A committed `test_bezier_tet10.py` under
`Ladruno_scripts/bezier_tests/` (where Tri6's regression lives), kept green.

## Risks / open questions

> [!question]
> **O11 — node-order reconciliation (the #1 correctness risk).** D2′ is read off the
> `TenNodeTetrahedron::shp3d` derivatives (vertices then edges 1‑2/2‑3/1‑3/1‑4/3‑4/2‑4,
> post the Larenas N9↔N10 swap). **Must be triple-confirmed against (a) Gmsh's native
> tet10 ordering and (b) what apeGmsh emits**, before trusting any multi-element result —
> a wrong edge order silently produces a valid-looking but wrong stiffness. The Tri6
> integration test caught this for T6 (deviation 1.78e-15); replicate exactly for T10.
> If Gmsh disagrees with `TenNodeTetrahedron`, the element must match **the mesh source
> apeGmsh uses**, and the discrepancy gets a permutation map + a loud comment.

> [!question]
> **O12 — consistent-mass quadrature order.** Lumped `ρVe/10` (the default, D4) needs
> only the 4-pt rule. **Consistent** mass `∫ρNᵀN` integrates degree-4 polynomials and
> needs a higher-order tet rule (≥11-pt) for exactness on straight-sided elements —
> pick and document it (Tri6 used a 6-pt Dunavant rule for its degree-4 consistent
> mass). Low risk (consistent mass is opt-in) but get it exact, not approximate.

> [!question]
> **O5/O6/O7, O8, O9 — inherited deferrals.** Curved faces + Eq. 14 Dirichlet map,
> real surface traction, selectable quadrature (O5–O7); the Chung–Lee explicit
> integrator (O8, → [[Ladruno_explicit_roadmap]]); apeGmsh first-class typed primitive
> (O9). All deferred *exactly as in Tri6* — do not implement piecemeal.

- **Class-tag selection** — shipped as `ELE_TAG_BezierTet10 = 33001` in the fork's
  `classTags.h` (ladruno private band ≥33000; the early candidate `273` was dropped).
- **Validation bar** — "done" gated by Kadapa §6.2: thick-sphere convergence, the
  TET4/TETB10/TETB10B locking study at ν=0.49999, lumped-mass positivity, plus the
  constant-stress patch test.
- **Backwards compat** — lumped-mass default (D4) means `eigen`/dynamics scripts that
  expect consistent mass must pass `-cMass`; flag in the element docs.

## Deferred follow-ups (parked — not active)

| Item | Ref | Trigger to revive |
|------|-----|-------------------|
| **Curved tet faces** + Eq-14 Dirichlet control-point↔DOF mapping | D9′ / O5 | a model needs genuinely curved Bézier faces (the element warns on curved meshes) |
| **Real surface traction** (replace any `-pressure` volume hack) | O6 | a real surface-pressure load case |
| **Selectable quadrature** (`-nip`) | O7 | curved elements (4-pt is exact for straight-sided) |
| **Rational / NURBS tet, higher order** | — | an IGA or exact-geometry need; own ADR |
| **Explicit integrator** (Chung–Lee β=13/12) — global `TransientIntegrator`, not element-side | O8 | Kadapa's exact damped scheme. Element is explicit-ready today with stock `CentralDifference` → [[Ladruno_explicit_roadmap]] |
| **apeGmsh first-class typed primitive** (`ops.element.BezierTet10`) | O9 | want bridge emission/run instead of direct-drive |

## Implementation log

- 2026-05-30 — **ADR drafted.** Sibling of [[04_bezier_elements]] under its D10.
  Tet10 math sourced from Kadapa (2018) §5: lumped mass Eq. 57 (`ρVe/10`), 4-pt
  quadrature, 3D B-bar Eq. 45, control-point map Eq. 11–13 (edge-wise). **Node order
  (D2′) and Voigt convention read off `TenNodeTetrahedron`** (`shp3d` edge order
  1‑2/2‑3/1‑3/1‑4/3‑4/2‑4 post Larenas swap; `computeB` rows `{xx,yy,zz,xy,yz,zx}` ==
  Kadapa Eq. 26–27). Recorder side is zero-code via the element contract (tet case
  already specified). Open: O11 (Gmsh node-order reconciliation), O12 (consistent-mass
  rule). **Next: write the Python oracle (`bezier_tet10_element.py`) — Track C before
  any C++**, exactly as Tri6 sequenced it.
- 2026-05-30 — **Track C Python oracle DONE + validated (23/23 green).** Authored
  `bezierFEM/python examples/bezier_tet10_element.py` (3D sibling of `bezier_tri6_element.py`).
  Verifies to machine precision: Bernstein partition-of-unity + nonnegativity + FD-checked
  derivatives; **D2′ node placement** (vertices N=Lᵢ², edges N=2LₐL_b in the `TenNodeTetrahedron`
  edge order 1‑2/2‑3/1‑3/1‑4/3‑4/2‑4); exact volume; **consistent mass == analytical simplex
  integral** (1e-12) and **all-positive**; **lumped mass = ρVe/10** exactly equal per DOF (Eq. 57);
  stiffness symmetric with **all 6 rigid-body modes** in the null space (1e-13); **constant-stress
  patch test exact to 2e-16** for *both* displacement and B-bar; B-bar SPD + patch-exact.
  **Formulation validated.** Two findings to carry to the C++ port: (1) **integrate with `|detJ|`**
  (not signed) — `BᵀDB`/`NᵀN`/`J⁻¹` are orientation-independent, so K/M/F are correct for *either*
  vertex handedness; pragmatic hedge for **O11** (Gmsh-tet orientation vs D2′ still unverified — a
  consistently-negative detJ is reversed handedness, not a tangle). (2) **consistent-mass quadrature
  needs the 4×4×4 collapsed-Duffy rule** (64-pt, degree-7), *not* 27-pt: the (1-a)² Duffy Jacobian
  pushes the degree-4 `NₐNᵦ` integrand to degree 6 in the first collapsed variable — **resolves O12**
  (production 4-pt is for K/F/B-bar-avg; mass needs the higher rule). **Next: `validate_thick_sphere.py`
  (§6.2 convergence + TET4/TETB10/TETB10B locking study), then fork registration (new dir
  `SRC/element/bezierTetrahedron/`, `ELE_TAG 273`, dispatch, CMake) + L0 C++ smoke.**
- 2026-05-30 — **Track C validation DONE — `validate_thick_sphere.py` PASSES (both benchmarks).**
  Self-contained (no gmsh; controlled D2′ ordering; sparse assemble/solve). **Benchmark 1 — thick
  sphere (§6.2):** octant shell of straight-sided T10 tets (barycentric-refined spherical triangle ×
  radial layers → conforming prism-split tets), consistent internal-pressure load on inner faces,
  symmetry BCs. Both formulations converge **monotonically to the Lamé solution** at **rate ~1.8**
  (n=2→5). **KEY FINDING:** with **straight-sided** tets the faceted-boundary **geometry error
  (O(h²)) DOMINATES and masks Poisson locking** — displacement and B-bar are nearly identical on the
  sphere (rate ~2, not the optimal ~3). This is *expected* and is why Kadapa uses **curved** elements
  for the sphere (our deferred **O5**); the curved-geometry convergence study is genuinely blocked on
  O5. **Benchmark 2 — near-incompressible bending beam (FLAT geometry ⇒ straight tets EXACT, the 3D
  analogue of Tri6's Cook's membrane):** isolates locking cleanly. Result is textbook: at ν=0.3 disp
  & B-bar agree (δ 3.98 vs 4.23, peak |p_hyd| ~0.02-0.03); at ν=0.4999 the **displacement form
  develops a spurious hydrostatic-pressure spike ~6× physical** (peak |p_hyd| 0.186 vs B-bar 0.036,
  **ratio 5.2×**) while B-bar deflection is **ν-insensitive** (2.8%). **Confirms the paper's actual
  claim: for QUADRATIC elements displacement locking is mild in *displacement* (Kadapa: adequate to
  ν≈0.48) — B-bar's payoff is the *stress/pressure field*** (so a tip-deflection-ratio pass criterion
  is wrong for quadratic; the spurious-pressure metric is the right one, mirroring the Tri6 finding).
  **Formulation fully validated. Track C complete. Next: fork registration (new dir
  `SRC/element/bezierTetrahedron/`, `ELE_TAG 273`, dispatch, CMake) + L0 C++ smoke.**
- 2026-05-30 — **C++ element WRITTEN + compile-clean.** Ported the validated oracle to
  `SRC/element/bezierTetrahedron/{BezierTet10.{h,cpp},OPS_BezierTet10.cpp,CMakeLists.txt}`,
  mirroring `BezierTri6` lifted to 3D (NEN=10, NDOF=3, NELD=30, NSTRESS=6, NGAUSS=4; static
  return matrices, control-point map, B/B-bar, lumped/consistent mass, `basisInfo` response trio,
  SelfWeight, straight-sided guard on all 6 mid-edges). Both oracle findings baked in: integrate
  with **`fabs(detJ)`** (O11 hedge — `computeJacobian` returns signed detJ, computes J⁻¹ for either
  sign, callers use `|detJ|`); consistent mass uses the **4×4×4 collapsed-Duffy rule** built inline
  from a 1D 4-pt Gauss-Legendre table (O12). Voigt {xx,yy,zz,xy,yz,zx} == `TenNodeTetrahedron`.
  Wiring (4 points): `classTags.h` `ELE_TAG_BezierTet10 = 33001`; `add_subdirectory(bezierTetrahedron)`
  in `SRC/element/CMakeLists.txt`; fwd-decl + `functionMap` "BezierTet10"/"bezierTet10" in
  `OpenSeesElementCommands.cpp`. `OPS_BezierTet10` is plain C++ linkage (not `extern "C"`).
  **Standalone `cl /c` compile-check of both new TUs = exit 0, zero warnings** (flags from the main
  tree's `build.ninja`, `-I<worktree>\SRC` prepended so classTags.h tag 273 resolves; recipe from
  [[project_mpco_ladruno]]). L0 smoke test written: `test_bezier_tet10_smoke.py` (gitignored root
  helper) — 3D patch test (disp + B-bar), gaussPoint response, lumped-mass eigen positivity.
  **Next: full OpenSeesPy build in the worktree + run the smoke test (the L0 gate).**
- 2026-05-30 — **BUILT + L0 SMOKE PASSES (4/4).** Full `build.bat OpenSeesPy` in the worktree
  (exit 0; fresh `dist/bin/opensees.pyd`, both `bezierTetrahedron/*.cpp.obj` compiled and linked
  into the DLL). Copied the main tree's `mumps-install/` into the worktree to skip the MUMPS step
  (FetchContent still rebuilt MUMPS — known ~20-min cost, harmless). `test_bezier_tet10_smoke.py`
  run with the BUILD python (`pythoncore-3.12-64`, no venv boot `.pth`) + `os.add_dll_directory`
  on `dist/bin` for MKL → **4/4 PASS:** patch test σ_xx=1.346154 **exact** at all 4 GPs
  (σ_yy/zz=0.576923 exact); gaussPoint 12 vals; lumped-mass default → 3 positive eigenfreqs;
  3D B-bar patch test exact. (One transient hiccup: the smoke print used a `≈`/`~` char that
  cp1252 chokes on — run with `PYTHONIOENCODING=utf-8`; element math was never affected.)
  **C++ element fully validated end-to-end. v1 element COMPLETE. Next: apeGmsh integration
  (verify Gmsh tet10 order == D2′ — resolves O11) + PR to ladruno (mirror Tri6 PR #6).**
- 2026-05-30 — **RECORDER SIDE WIRED (MPCO_Ladruno inverted dispatch) — round-trip 15/15.**
  Finding: the recorder side was a *documented stub* — Tri6 only ever implemented the **element**
  half of the contract (`basisInfo`/`integrationPoints`/`integrationWeights`); the recorder's
  consumption was reserved as `TODO(Step 3)` in `MPCOL_ElementResults.h` and **never implemented**
  (neither Bézier element was referenced anywhere in the recorder). Implemented Step 3, serving
  Tri6 + Tet10 + every future contract element: **(1)** a `BasisInfoCaptureStream : DummyStream`
  (captures the `ElementBasis` attrs) + `captureBasisInfo`/`captureResponseMatrix`/`...Vector`
  helpers; **(2)** `writeModelElements` now probes `basisInfo` first → overrides `FAMILY=bernstein`,
  `ORDER=[2]`, writes the element's own `GP_PARAM`/`GP_WEIGHT`, legacy class-tag table as fallback;
  **(3)** added `BezierTri6→Triangle_6N` and `BezierTet10→Tetrahedron_10N` (+ real 3-pt/4-pt Gauss
  rules) to `getGeometryAndIntRuleByClassTag` so the engine gets correct node/GP counts. **Two
  latent bugs in the ported recorder fixed:** (a) `getCustomGaussPointLocations` expects
  `integrationPoints`→Vector but the contract returns a Matrix → sidestepped via the real Gauss
  rule (non-custom); (b) **HDF5 cannot create a group under a dataset** — the ported `QUADRATURE`
  child-group code crashes (`group::create(dset_id,…)`); rewrote to store `GP_PARAM`/`GP_WEIGHT`
  as **array attributes** on the element dataset (like the existing `GP_X`). Verified: model runner
  `recorder_tet10_model.py` (build py) + `recorder_tet10_check.py` (venv h5py) → **15/15** —
  `FAMILY=bernstein`, `TOPOLOGY=tet`, `NUM_CTRL=10`, `NUM_GP=4`, `ORDER=[2]`, `PARAM_DOMAIN=bary`,
  `CONNECTIVITY (1,11)`, `GP_PARAM` = exact 4-pt tet rule, `GP_WEIGHT`=[1/24]×4. **Class tag moved
  to the ladruno private band: `ELE_TAG_BezierTet10 = 33001` (≥33000, avoids upstream collision;
  was 273). Macro-based code picks it up on rebuild — current binary still shows tag 273.**
- 2026-05-30 — **REBASED onto current `origin/ladruno` → recorder wiring SUPERSEDED by the canonical
  merged implementation; Tet10 wiring is now a 1-block add.** The branch was 40+ commits behind; the
  rebase brought in the *real* recorder work that my pre-rebase session had reinvented in parallel
  (worse): **PR #16** (the exact HDF5 group-under-dataset fix), **PR #18** ("Wire BezierTri6 into
  MPCO_Ladruno + consume basisInfo self-declaration" — the inverted dispatch already merged, with a
  generic `getElementBasisInfo` capture, a multi-dim `getCustomGaussPointLocations` that reads
  `integrationPoints` as a `[nGP×dim]` Matrix verbatim, and **`MODEL/ELEMENTS/<name>` as a GROUP**
  holding `CONNECTIVITY` + a `QUADRATURE` child group with `GP_PARAM [nGP×ndir]` + `GP_WEIGHT [nGP]`
  datasets — the canonical bucket-as-group schema), **PR #21** (BezierTri6 272→33000 band), **PR #23**
  (`mpco-quadrature-global-gp`). **Action:** discarded all my redundant recorder edits, rebased
  (classTags conflict resolved: keep BezierTri6=33000 + add BezierTet10=33001). **Tet10 recorder
  wiring is now exactly ONE block** in `getGeometryAndIntRuleByClassTag`: `BezierTet10 →
  Tetrahedron_10N + CustomIntegrationRule + custom_rule_dimension=3` (mirrors BezierTri6's dim=2).
  Everything else is generic and reused — the element already emits the contract (basisInfo
  topology=tet/numCtrl=10/numGP=4/orderU=2, integrationPoints `Matrix(4,3)`, integrationWeights
  `Vector(4)`); ladruno's recorder reads `orderU` into a single `order` and replicates it per
  GP_PARAM column → `ORDER=[2,2,2]`. **How results are stored (canonical, adopted):** bucket-as-group
  — `MODEL/ELEMENTS/<33001-BezierTet10[...]>/{CONNECTIVITY, QUADRATURE/{GP_PARAM[4×3],GP_WEIGHT[4]}}`,
  FAMILY=bernstein on the group. Updated `recorder_tet10_check.py` to assert that group schema.
  Ledgers + banner updated. Rebuild (classTags band ⇒ near-full) in progress; re-verify L0 smoke +
  recorder round-trip after. **NEXT: rebuild verify, then PR to ladruno + apeGmsh O11 check.**
- 2026-05-31 — **REBASED AGAIN onto current `origin/ladruno` (PR #58) — the recorder was RENAMED
  MPCO_Ladruno → Ladruno mid-session.** During the previous rebuild `origin/ladruno` advanced 32
  commits (PR #26→#58): `61eb2be9e` renamed the whole recorder (`MPCORecorderLadruno→LadrunoRecorder`,
  `MPCOL_*→Ladruno_*`, cmd `mpcoLadruno→ladruno`, `GENERATOR→"Ladruno"`) + added a profiler subsystem
  (P7/P8/viewer). **The Tet10 wiring survived cleanly:** git rename-detection auto-applied my
  `BezierTet10` registration to the renamed `Ladruno_ElementResults.h` (still
  `Tetrahedron_10N + CustomIntegrationRule + dim=3`); `classTags.h` auto-merged. Resolved 6 doc/banner
  conflicts (keep origin's comprehensive ledger/README rows + add the BezierTet10 rows; regenerated
  the banner). Updated the recorder test scripts for the rename: `recorder('ladruno', …)`,
  `GENERATOR=='Ladruno'`. **The element + recorder mechanism are rename-agnostic** — only the cmd name,
  GENERATOR string, and file names changed. Full rebuild + re-verify in progress.**
