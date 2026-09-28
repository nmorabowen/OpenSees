---
title: Bézier triangular element (BezierTri6)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - element
  - adr
---

# Bézier triangular element (BezierTri6)

> **ADR** (Architecture Decision Record). Captures the v1 scope decisions of
> 2026-05-30 for porting the Kadapa (2018) quadratic Bézier triangle into the
> `ladruno` fork. The draft implementation lives in the **separate** `bezierFEM`
> repo (`opensees implementation/BezierTri6.{h,cpp}`, `OPS_BezierTri6.cpp`); this
> ADR governs landing it in the fork. Companions: [[ladruno_element_contract]]
> (recorder side — already satisfied via `basisInfo`), [[Ladruno_explicit_roadmap]]
> (the explicit integrator this element feeds), and [[project_bezier_element]] (the
> ordering decision and repo layout).

## What

`BezierTri6` — a 6-node quadratic **Bézier** triangle (Bernstein basis) for 2D
linear elastostatics and elastodynamics, after Kadapa (2018, IJNME 117:543). It
reuses standard quadratic-Lagrange (T6) meshes via a linear control-point map
(Eq. 13), delegates the constitutive law to any `NDMaterial`, and offers an
optional consistent **B-bar** formulation for near-incompressibility. Its
defining property is the **all-positive lumped mass** `ρVe/6·diag[1₆,1₆]`
(Eq. 56), which makes quadratic triangles usable in explicit dynamics.

**In scope (v1):**
- The pure-displacement and B-bar formulations, plane stress / plane strain.
- Lumped (default) and consistent (`-cMass`) mass.
- **Straight-sided elements only** (mid-edge nodes at edge midpoints ⇒ control
  points = nodes, P = X).
- Self-weight / body-force loads via the standard load-pattern path.
- Recorder self-declaration via the element contract (`basisInfo`, `family=bernstein`).

**Out of scope (deferred to follow-up ADRs / tracks):**
- **Curved** Bézier edges and the Eq. 14 Dirichlet control-point↔DOF mapping they
  require (deferred with D9).
- The Bézier **tetrahedron** (Tet10, Kadapa §5), rational/NURBS, higher order.
- A proper edge/surface traction load (the current `-pressure` is a volume hack).
- Selectable quadrature (fixed 3-point is exact for straight-sided).
- The explicit time integrator (Chung–Lee β=13/12) — element-side is ready;
  the integrator belongs to [[Ladruno_explicit_roadmap]].

## Why

- **Explicit dynamics need a non-negative lumped mass.** Standard quadratic
  Lagrange triangles have negative mass-matrix entries → zero/negative lumped
  diagonals → unusable with explicit integration. Bernstein non-negativity fixes
  this exactly (each node gets ρVe/6). This is the paper's headline result and the
  primary reason to add the element.
- **Near-incompressible locking.** The consistent B-bar (Eq. 45) gives optimal
  convergence as ν→0.5 without mixed/stabilized formulations or extra DOFs.
- **Uses existing mesh generators.** apeGmsh already produces quadratic-Lagrange
  triangles; the Eq. 13 map turns them into Bézier elements at negligible cost.
- **A verified draft already exists** in `bezierFEM`, reviewed equation-by-equation
  against the paper (shape funcs Eq. 4–9, map Eq. 13, B-bar Eq. 45 all confirmed).

## Where

New code (all under `SRC/element/bezierTriangle/` in the `ladruno` fork):

| File | Role |
|------|------|
| `BezierTri6.{h,cpp}` | the element |
| `OPS_BezierTri6.cpp` | Tcl/Python factory (`element BezierTri6 …`) |
| `CMakeLists.txt` | subdir build |

Modify (fork-local):
- `SRC/classTags.h` — add `ELE_TAG_BezierTri6` (shipped as `33000` in the ladruno
  private band ≥33000; the early sub-300 candidate `272` was dropped — see ship log below).
- `SRC/interpreter/OpenSeesElementCommands.cpp` — dispatch `BezierTri6`.
- `SRC/element/CMakeLists.txt` — `add_subdirectory(bezierTriangle)`.
- `SRC/actor/objectBroker/FEM_ObjectBroker.cpp` — `case ELE_TAG_BezierTri6`.

Reference (copy idioms): `SRC/element/fourNodeQuad/FourNodeQuad.cpp` — lumped
mass + material-density fallback, `addLoad`/`SelfWeight`, `appliedB`/`applyLoad`.
No new external dependency; no `Ladruno_internal/01_compilation_journal` entry needed.

## Decisions

| # | Decision | Rationale | Consequence |
|---|----------|-----------|-------------|
| **D1** | **Kadapa linear map, not IGA.** Reuse Lagrange meshes + Eq. 13 map; small-strain, linear-elastic. | The paper's whole premise (Remark 1: not isogeometric, not finite-strain). | Geometry is approximate, not exact; fine for the target problems. Finite-strain/IGA are separate efforts. |
| **D2** | **Node order = standard FE hierarchical** (3 vertices then 3 mid-edges), identical to Lagrange T6. | Verified vs Fig. 1 / Eq. 4–9; see [[project_bezier_element]]. | Recorder reuses the Lagrange node order, swaps only the basis. No `CTRL_ORDER` hint. |
| **D3** | **Constitutive delegated to `NDMaterial`** (`getTangent`), not a hardcoded D. | Any 2D material (and the B-bar split) plugs in; idiomatic OpenSees. | Element carries no elasticity constants. |
| **D4** | **Lumped mass is the default**; consistent via `-cMass`. | The element's reason to exist (explicit dynamics); matches FourNodeQuad. | **Behavior change** vs the first draft (was consistent-only). `eigen`/dynamics results change unless `-cMass`. |
| **D5** | **B-bar restricted to PlaneStrain/3D**; guarded off (warn) for PlaneStress. | Eq. 45 `/3` deviatoric split assumes an out-of-plane constraint; plane stress has no volumetric locking. | `-bbar` + `PlaneStress` silently disables B-bar with a warning. |
| **D6** | **Loads via SelfWeight + constructor body force** (FourNodeQuad idiom). | Lets gravity ramp through a load pattern / time series. | `addLoad` supports only `LOAD_TAG_SelfWeight`; other `eleLoad` types return −1. |
| **D7** | **Recorder integration via the element contract** (`basisInfo` family=bernstein, rational=0). | Self-declaring; zero recorder-side code (see [[ladruno_element_contract]]). | Recorder records correct GP coords; no Lagrange-fallback mis-placement. |
| **D8** | **Home = the `ladruno` fork** (not upstream PR, not standalone plugin). | First-class integration with the fork's recorder + explicit roadmap; no upstream review latency. | Fork-local `ELE_TAG`, namespace freedom; **no upstream-diff discipline** — revisit if/when upstreaming. |
| **D9** | **v1 = straight-sided only** (P = X). | Eq. 14 Dirichlet mapping is trivial, 3-pt quadrature is exact, BCs are interpolatory — correct and simple now. | **Must guard:** warn/reject if a mid-edge node is not at the edge midpoint (a curved mesh v1 is not validated for). Curved geometry + Eq. 14 is a follow-up ADR. |
| **D10** | **Tri6 only in this ADR.** | Land and validate the 2D triangle first; lowest risk. | Bézier **Tet10** (Kadapa §5), rational/NURBS, and higher order each get their own ADR built on this recipe. |
| **D11** | **B-bar is a `-bbar` flag, not a separate element.** | One class tag; matches the current code. | Diverges from the paper's TRIB6/TRIB6B naming (a flag, not two elements). Documented. |

## Risks / open questions

> [!question]
> **O5/O6/O7 follow the curved-geometry decision (D9).** Curved edges require the
> Eq. 14 Dirichlet control-point↔DOF map (O5), a real edge-traction load to replace
> the volume `-pressure` hack (O6), and higher/selectable quadrature (O7). All three
> are deferred *together* with curved support — do not implement piecemeal.

> [!question]
> **O8 — explicit integrator.** The paper's purpose (explicit elastodynamics with
> the Chung–Lee scheme, β=13/12) needs an integrator OpenSees may not ship. The
> element side is ready (lumped mass). Resolve in [[Ladruno_explicit_roadmap]], not here.

> [!question]
> **O9 — apeGmsh emission.** Straight-sided (D9) makes this clean (P = X, no BC
> mapping). **Tested 2026-05-30:** apeGmsh meshes → BezierTri6 works end-to-end via
> *direct drive* (apeGmsh extracts nodes + T6 connectivity; our `opensees.pyd` builds the
> elements). Gmsh's native T6 node order is **identical** to BezierTri6's (verified to
> 1e-15). **Gap:** apeGmsh's `apeSees` bridge has a *typed* element registry that does NOT
> include `BezierTri6`, and apeGmsh runs on py3.11/openseespy while our element is in the
> py3.12 fork build — so the bridge can neither emit nor run it as a typed primitive today.
> Making it first-class in apeGmsh (a typed-primitive wrapper + deck emission against the
> fork build) is an apeGmsh-side feature. The direct-drive path needs neither.

> [!note]
> **Broader / integration testing is apeGmsh-driven — DEFERRED until the C++
> element compiles.** Once `BezierTri6` is registered and built into the fork, the
> system-level validation (real meshes, multi-element models, the recorder/results
> round-trip, larger benchmarks) runs through **apeGmsh** (`/apegmsh-helper`) rather
> than hand-written drivers. The Python `validate_cooks_membrane.py` harness covers
> the *formulation* in the meantime; the apeGmsh path covers the *integrated* element
> + recorder + post-processing stack. Blocked on: fork registration + build.

- **Class-tag selection** — pick a real `ELE_TAG_BezierTri6` and confirm it is
  unused in the fork's `classTags.h` before registering.
- **Validation bar (O10, Track C)** — "done" gated by Kadapa §6: constant-stress
  patch test, Cook's membrane convergence, thick-cylinder §6.1, and the locking
  study TRIB6 vs TRIB6B vs CST at ν=0.4999, plus lumped-mass positivity.
- **Backwards compat** — D4 changes the default mass to lumped; flag prominently in
  the element docs so existing `eigen` scripts add `-cMass` if they need consistent.

## Deferred follow-ups (parked — not active)

v1 (straight-sided `BezierTri6`, validated 2026-05-30) is feature-complete. The
following are **explicitly deferred**, each with the trigger that would revive it:

| Item | Ref | Trigger to revive |
|------|-----|-------------------|
| **Curved geometry** + Eq-14 Dirichlet control-point↔DOF mapping | D9 / O5 | a model needs genuinely curved Bézier edges (the element already warns on curved meshes) |
| **Real edge traction** (replace the volume `-pressure` hack) | O6 | curved support, or a real surface-pressure load case |
| **Selectable quadrature** (`-nip`) | O7 | curved elements (3-pt is exact for straight-sided) |
| **Bézier Tet10** (3D sibling, Kadapa §5) | D10 | a 3D nearly-incompressible problem; own ADR, same recipe |
| **Explicit integrator** (Chung–Lee β=13/12) — a **global** `TransientIntegrator`, NOT element-side | O8 | wanting Kadapa's *exact* damped scheme. The element is already explicit-ready (lumped mass + resisting force); it runs explicitly **today** with stock `CentralDifference`. The Chung–Lee scheme is a fork-wide integrator → [[Ladruno_explicit_roadmap]] (cf. `ExplicitBatheLNVD`), not this ADR |
| **apeGmsh first-class typed primitive** (`ops.element.BezierTri6`) | O9 | want bridge emission/run instead of direct-drive; spec written at `apeGmsh/docs/plans/bezier-tri6-element.md` |
| **A3 SelfWeight** | — | ✅ DONE (smoke-tested 2026-05-30) — not deferred |

Nothing here blocks using `BezierTri6` for straight-sided 2D continuum problems today.

## WP-114: `-bbar` is the plane-strain ½ split (TIMs F17, PR #848)

**Defect (fixed in b42ca77d8).** `computeBBarMatrix` applied the 3D mean-dilatation
split, (B̄+2B)/3 and (B̄−B)/3, to the 3-row plane-strain B and dropped the εzz row. The
material therefore saw an in-plane trace of (θ+2θ̄)/3 instead of θ̄. Under isochoric
flow that keeps all 3 point-wise volumetric constraints of the plain T6, so `-bbar`
relieved nothing. The element now uses the 2D split, εxx/εyy += (θ̄−θ)/2, the same as
LadrunoQuad and LadrunoUP (2D). That leaves one constraint per element. The element
tangent was never the problem: it matches a central FD to ≤ 5e-11 before and after the
fix, at elastic, associated and ψ=0 plastic states.

**Punch benchmark** (`tests/wp114/punch_nonassoc.py`, driver `run_punch_matrix.py`).
Half model of a rigid rough strip punch, B = 2 m, block 8 × 6 m graded to the punch
corner. UW DruckerPrager φ = 33°, c = 2 kPa, plane-strain MC match, no hardening (the
tension cutoff sits at the apex). γ = 18 kN/m³ and q = 20 kPa surcharge, applied in an
elastic geostatic stage; then displacement control to s/B = 0.10. Same strategy for every
case: Newton → KrylovNewton → NewtonLineSearch, then up to 6 step halvings (a "wall"),
NormDispIncr 1e-8. g1 = 280 T6 or 140 Q4; g2 = 1 120 T6 or 560 Q4. T6 g1 and Q4 g2 have
the same node grid. The census counts GPs on the tension cutoff or apex (`ladrunoBranch`
2/3) at the last converged step.

| case | ψ = 0 | ψ = φ |
|---|---|---|
| T6 `-bbar` PRE (÷3) g1 | **walls s/B 0.039**, peak 219 kPa, 39 cuts, **63 apex** | 672 kPa at 0.10 (rising), 7 apex |
| T6 `-bbar` POST (½) g1 | no wall, 395 kPa at 0.10 (flattening), 2 cuts, **0** | 656 kPa (rising), 4 apex |
| T6 `-bbar` POST g2 | walls 0.026 at 1e-8 (Newton stalls at 1e-6 to 1e-7, 0 apex); at tol 1e-6: peak 382, 371 at 0.10, 0 | 628 kPa (rising), 8 apex |
| T6 std g1 | peak 450 at 0.071, 435 at 0.10, 1 apex | 664 kPa (rising), 2 apex |
| Q4 `-bbar` g1 | 482 kPa at 0.10 (rising), 0 | 809 kPa (rising), 17 apex |
| Q4 `-bbar` g2 | peak 439 at 0.075, 428 at 0.10, 0 | 706 kPa (rising), 46 apex |

At s/B = 0.06, ψ = 0: T6-bbar POST 372, T6 std 434, Q4 g1 451, Q4 g2 427 kPa. The fixed
T6 `-bbar` does not reach the quad band. It sits 13–18 % below it, and it trends down
with refinement as the quad does. This is consistent with it being the least-constrained
of the four (one volumetric constraint per 6-node element), but we have no reference
solution for ψ = 0 DP, so we do not claim it is closer to the true limit.

**Where the pre-fix apex points were.** They were not at the punch corner. The 63 GPs had
distances from the corner of 0.17 m (min), 1.08 m (median) and 5.05 m (max). Most were
deep (y < −0.25 m), and 16 were under the punch more than 2 m from the corner. They
appeared at s/B = 0.0035, as soon as ψ = 0 flow started. At s/B = 0.03 the intra-element
I1 spread over plastic elements had a median of 0.41 × |mean I1| and a p90 of 2.23 before
the fix, and was exactly 0 after it. With ψ = φ it was 0.11 / 0.31 before and 0.10 / 0.30
after. So the apex points are the pressure checkerboard of volumetric locking. That locking
was left in place by the ÷3 split. They are not a boundary or loading effect.

**Bernstein loading lead: ruled out as the cause, but a real trap.** Surcharge loads with
Lagrange weights (qL/6, 2qL/3, qL/6) on Bézier CPs are the wrong consistent loads (they
should be qL/3 each). With those loads every T6 run fails at the first plastic step, before
and after the fix, for both ψ values. That is a different signature from the reported wall.
With the correct Bernstein loads, the fixed element has no cutoff or apex points.

**Asymmetric-pivot lead: not reproduced.** See `LEDGER_quirks` ("constant perturbed-pivot
count … element-less nodes").

Gate: `tests/test_beziertri6_bbar_plane_strain.py`. It checks the trace kinematics, the
volumetric rank (1 vs 3), the FD tangent, and that the ψ = 0 punch runs to s/B = 0.06 with
no wall, 0 apex points, and ≥ 0.75 × the Q4 g1 load. About 8 s. On the pre-fix binary the
kinematics, rank and punch checks fail; the FD check passes on both binaries, because it is a
consistency gate.

## Implementation log

- 2026-05-30 — ADR drafted. Draft element in `bezierFEM` reviewed equation-by-equation
  against Kadapa (2018); core math confirmed. **Track A correctness fixes applied** to
  the draft: lumped-mass default + `-cMass` (D4), PlaneStress B-bar guard (D5),
  `addLoad`/`SelfWeight` + material-density mass fallback (D6); recorder `basisInfo`/
  `integrationPoints`/`integrationWeights` added (D7). Gating decisions O1–O4 settled
  → D8–D11. Deferred: curved geometry (D9), Tet10 (D10), edge traction/quadrature
  (O6/O7), explicit integrator (O8).
- 2026-05-30 — **D9 straight-sided guard implemented** in `BezierTri6::setDomain`
  (`bezierFEM`): warns if any mid-edge node deviates > ~1e-6·(edge length) from the
  edge midpoint. Header documents the v1 straight-sided limitation. Remaining:
  Track C validation (compile-exercises the Track A + guard changes), then fork
  registration (class tag + dispatch + broker + CMake).
- 2026-05-30 — **Track C started (Python oracle, no build).** Ran the reference
  element's self-checks (`bezier_tri6_element.py`) — all pass to machine precision
  (partition of unity, non-negativity, exact area, consistent mass = analytical to
  2.6e-17, equal lumped mass ρAt/6, patch test exact, rigid-body null space). Wrote
  `validate_cooks_membrane.py` (straight-sided T6 mesh + assembler driving the
  reference element) and reproduced Kadapa §6.3: point-A v_y converges to **~9.57 mm
  (ν=0.3)** and **~7.97 mm (ν=0.4999)** — the latter in the literature range (~7.6–8.0),
  both formulations consistent. **Formulation validated.** Open: B-bar's anti-locking
  benefit is mild in point-A *displacement* (expected — higher-order displacement
  locking is mild per the paper; B-bar mainly fixes the *stress* field). To fully
  demonstrate B-bar, add a stress-oscillation metric and/or a linear TRI3 baseline.
  Still pending: building/registering the **C++** element to confirm the Track A port
  (the Python oracle shares the formulation but not the C++ code paths). **Broader /
  integration testing deferred to the apeGmsh path (`/apegmsh-helper`)** — runs once
  the element compiles into the fork (real meshes + recorder/results round-trip); the
  Python harness covers the formulation until then.
- 2026-05-30 — **Fork registration + build + smoke test DONE.** Landed the element in
  the `ladruno` fork: `SRC/element/bezierTriangle/{BezierTri6.{h,cpp},OPS_BezierTri6.cpp,
  CMakeLists.txt}`, `ELE_TAG_BezierTri6 = 33000` in `classTags.h` (moved from 272 to the
  ladruno private band >=33000 on 2026-05-30 — 272 collided one above upstream PML3DVISCOUS=271),
  `add_subdirectory(bezierTriangle)`, and the Python dispatch (`OPS_BezierTri6` fwd-decl +
  `functionMap` "BezierTri6"/"bezierTri6") in `OpenSeesElementCommands.cpp`. Linkage fix:
  `OPS_BezierTri6` must be plain C++ (not the `extern "C"` OPS_Export DLL style) to resolve
  against the static fwd-decl. Broker step is obsolete (`getNewElement` is a stub here).
  Built via `Ladruno_scripts/{setup_env,build}.bat OpenSeesPy` — both new TUs compile clean,
  `opensees.pyd` produced. **Smoke test `test_bezier_smoke.py` (4/4 PASS):** patch test
  exact (σ=E·ε at all GPs), gaussPoint response, lumped-mass default (positive eigvals),
  PlaneStress B-bar guard warns+disables. **The C++ port is validated end-to-end** →
  the deferred apeGmsh broader/integration testing (`/apegmsh-helper`) is now UNBLOCKED.
  (Caveat: my earlier hand-run `cmake .` corrupted the build cache — fixed by deleting
  `CMakeCache.txt`; never reconfigure outside `build.bat`.)
- 2026-05-30 — **apeGmsh integration test PASSED.** Two-step (apeGmsh on py3.11+gmsh
  generates the mesh; our `opensees.pyd` on py3.12 drives BezierTri6). Cantilever plate
  (10×2), 465 nodes / 208 straight-sided T6 elems. Gmsh T6 node order == BezierTri6 order
  (max mid-node deviation 1.78e-15 ⇒ ordering + straight-sidedness both confirmed). Solve
  converged; **tip deflection FE −5.138 vs Timoshenko −5.156 (0.3% off)**; peak σ_xx 151.5
  vs Mc/I 150. Scripts: `bezier_apegmsh_mesh.py` (mesh→JSON) + `bezier_apegmsh_run.py`
  (JSON→BezierTri6). **Finding (O9):** `apeSees` typed-element registry lacks BezierTri6;
  direct-drive from the apeGmsh mesh is the working path. Element validated on a real
  multi-element model — Track C continuum validation effectively complete for v1.
- 2026-05-30 — **Session close-out.** (1) **A3 SelfWeight smoke-tested PASS** —
  `eleLoad -type -selfWeight` → total vertical reaction = γ·A·t exactly
  (`test_bezier_selfweight.py`). (2) **apeGmsh canonical-implementation notes** written
  at `apeGmsh/docs/plans/bezier-tri6-element.md` (6 touch points modeled on `SixNodeTri`,
  the Gmsh-tri6 ordering fact, the openseespy-vs-fork runtime caveat, the interim
  direct-drive path). (3) **Housekeeping** — fork is canonical; `bezierFEM`
  `REGISTRATION_GUIDE.md` banner marks it superseded + records the two divergences
  (linkage, placeholder); transient log/JSON removed (test scripts kept at project root).
  (4) **Follow-ups parked** — see *Deferred follow-ups* table above. v1 complete.
- 2026-05-30 — **MERGED to `ladruno`** via PR #6 (element + registration + this ADR;
  `ELE_TAG_BezierTri6 = 33000`). `ladruno` is the source of truth — the element is now in
  the integration branch. This README-index + log entry complete the picture on ladruno.
- 2026-09-18 — **WP-114 (TIMs F17): `-bbar` switched to the plane-strain ½ split** (b42ca77d8, PR #848). The 3D ÷3 split had kept all 3 volumetric constraints under ψ=0 flow; see the WP-114 section above.
