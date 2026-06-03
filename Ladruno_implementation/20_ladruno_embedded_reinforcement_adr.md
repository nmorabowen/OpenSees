# ADR 20 — Embedded Reinforcement + Bond-Slip for 3D RC

**Status:** proposed, **revised after adversarial review** (workflow `wf_32a898cd`,
64 agents: 26 confirmed + 23 partially-confirmed findings, 6 refuted, 2 blockers).
Overall verdict: **proceed-with-changes**. This revision folds in the must-fix list.
**Companions:** [[19_rc3d_modeling_recipe]] (gap map), [[09_ladruno_brick]] /
[[06_bezier_tet10]] (concrete hosts), [[14_ladruno_rebar_buckling_adr]] /
[[13_ladruno_uniaxial_j2_adr]] (rebar materials), [[solid_transformation_wrapper]]
(corot/finite host kinematics). Memory: `project_rc3d_embedment_scoping`.

---

## 1. Context & problem

The fork can mesh 3D concrete (`LadrunoBrick` hex, `BezierTet10` tet) with a triaxial
plastic-damage concrete (`ASDConcrete3D`) and discrete-bar steel (`LadrunoUniaxialJ2` +
`LadrunoRebarBuckling`). The missing layer is **coupling discrete reinforcement to the
concrete continuum** under **both implicit and explicit** analysis.

Today's options are all deficient: `ASDEmbeddedNodeElement` / `EmbeddedBeamInterfaceP`
are **penalty** (implicit-only in practice); `EmbeddedBeamInterfaceL` is Lagrange; the
only auto-embed helper (`generateInterfacePoints`) is **Tcl-only** and carries confirmed
silent-failure bugs (§ D3); none provide a **1D rebar bond-slip** law.

## 2. Goals / non-goals

**Scope note (post-gate, §8):** these goals describe the full RC-3D *capability*. After the
gates, the **v1 deliverables are the conformal path** (recipe + bar-layout helper +
adaptive-stepping wrapper + confinement note — *no new C++*); the dual-mode embedded
element + bond-slip material below are the **v2** deliverables. "v1" in the goal line means
the capability as a whole, not the first shippable increment.

**Goals — BOTH implicit and explicit are in scope** (drives the dual-mode D2; the
implicit cyclic pushover on confined columns is a primary target, so R1/R7 must be
*measured*, not assumed):
- Discrete rebar (truss) embedded in a non-matching solid mesh.
- **Perfect-bond** mode AND **bond-slip** mode.
- Implicit (static + transient) AND explicit (CentralDifferenceLadruno / ExplicitBathe).
- **Host-agnostic** (hex `LadrunoBrick` + tet `BezierTet10`); **v1 = straight-sided hosts**.
- openseespy-native; set-to-set user API via apeGmsh.

**Non-goals (v1, now explicit with encoded scope-limit tests):** *quantitative* dowel
action / transverse bar shear — the default `truss` is axial-only, and the opt-in `beam`
rebar (D6) gives only an *approximate, mesh-limited* transverse contribution (true dowel =
local bearing at ~few `d_b`, needs a fine near-crack mesh or a dowel interface law → v2);
distributed Gauss-point bond (node-lumped first, O(h²));
**cyclic bond-slip degradation (→ v2)**; **cover-controlled splitting bond failure (→ v1.1,
xfail test)**; finite-strain tributary-length stretch correction (→ v2); multi-axial rebar;
**Mode T under corot/finite hosts (→ use Mode P)**; smeared reinforcement.

## 3. Decisions

### D1 — Discrete embedded reinforcement; conformal-mesh proven-out FIRST as a gate
Rebar = independent 1D `corotTruss` elements (carrying `LadrunoUniaxialJ2` /
`LadrunoRebarBuckling`) located inside host solids. Only this lets the fork's discrete-bar
models (buckling, fatigue, fracture) express themselves, and it decouples the meshes.
Scope-kill objections (smeared covers v1; conformal covers everything) were **refuted** —
the smeared stack (FSAM/MVLEM) is 2D-only, and conformal meshing is genuinely infeasible
for dense cages/joints. **But the embed-vs-conformal head-to-head is now a PRE-BUILD GATE
(Gate 1, §6), not a trailing validation item:** prove conformal (apeGmsh bar-layout +
`corotTruss` + `zeroLength` perfect-bond, no new C++) handles regular members and *measure
the empirical infeasibility boundary* (bar density / cage geometry) before building the
element. Element justified only past that boundary.

### D2 — One element, two modes (preserves the O(N) explicit solve)
`LadrunoEmbeddedRebar` (classTag **33005**, Element band) exposes:
- **Mode T (transformation / perfect bond, implicit-only):** rebar node translational DOFs
  slaved to host interpolation `u_r = Σ_i N_i(ξ_r) u_i^host` via an `MP_Constraint` resolved
  by the **Transformation** handler → rebar DOFs eliminated. No penalty.
- **Mode P (penalty / bond-slip, explicit-capable + implicit):** rebar keeps own DOFs + mass;
  anisotropic coupling assembles `K = Bᵀ diag(k_t, k_t, k_axial) B`,
  `B = [−I₃ | N₁I₃ … N_nI₃]`, with transverse `k_t` = volume-scaled penalty and axial
  `k_axial = dτ/ds · (πd) · L_trib` (perfect-bond axial ⇒ large penalty). NOT isotropic
  `iK·BᵀB` (would violate virtual-work symmetry).

**Why two modes (corrected wording):** the split **preserves the O(N) explicit
(diagonal-mass) solve**, it is not a free preference. `CentralDifferenceLadruno` does a
trivial **diagonal** M⁻¹ solve (`system Diagonal`, [CentralDifferenceLadruno.cpp:186]). A
transformation MPC condenses rebar DOFs into host DOFs (`TransformationDOF_Group::getTangent`
applies `Tᵀ M T` via `addMatrixTripleProduct`) → **densifies the mass** → `DiagonalSOE` would
silently assemble wrong inertia. Penalty keeps M diagonal (adds stiffness only) — exactly
why LS-DYNA CBIS uses penalty for explicit. **Lagrange multipliers with a localized bordered
solve also preserve diagonal mass at ~O(8N) and are a deferred alternative (R1).**
- **GUARD (must-fix):** `setDomain()` hard-errors when an explicit integrator is paired with
  Mode T / Transformation handler (prevents the silent dense-`TᵀMT`-into-`DiagonalSOE` bug).
- **dt_cr cost:** penalty `K_p` enters the coupled rebar-host eigensolve
  ([CriticalTimeStep.cpp:226]); stiff `K_p` on light rebar can cut `dt_cr` 2–10×. Heuristic
  `K_p ≲ 10·K_host`; verify via `-cfl`; recover the step with **rebar mass scaling** + the
  R1 audit below.
- Do NOT collapse to a single penalty mode: penalty-implicit conditioning at `K~1e18` is
  unproven (pre-build patch-test sweep in §6 decides if single-mode is viable).

### D3 — Host-agnostic isoparametric inverse-map (guarded; do NOT reuse legacy as-is)
Element stores `(hostEleTag, ξ)` and evaluates `N_i(ξ)` each step; inverse-map
`X = Σ N_i(ξ) X_i` solved at **build time** (reference config). **The legacy
`Tcl_generateInterfacePoints.cpp::invIsoMapping` carries CONFIRMED BLOCKERS and must NOT be
reused unguarded:** on non-convergence it returns −1 but leaves `inBounds=true`, its `opserr`
is commented out, and both call sites ignore the return code → silently embeds garbage ξ;
`Matrix::Solve` has no singularity/condition guard (det(R)→0 ⇒ inf/nan silently); Newton
tolerance is hardcoded **absolute** 1e-10 (won't converge for curved geometry). The new
inverse-map MUST: use **relative** Newton tolerance; check `det(R)`/condition before solve
(or SVD least-squares); check the return code at the call site; set `inBounds=false` + warn
on non-convergence; define an explicit **out-of-bounds policy** (reject-with-error vs snap).
**v1 supports straight-sided hosts only; the Bernstein (BezierTet10 curved) inverse-map is
deferred** and must be added with its own convergence study.

### D4 — 1D bond-slip τ–s uniaxial material (`LadrunoBondSlip`, MAT_TAG 33002)
**v1 = MC2010 monotonic backbone only** (cyclic → v2). τ vs slip s, exposing
`τ_max, s1, s2, s3, τ_f, α` as inputs (do not hardcode the regime), **plus `G_F^bond`** for
D4.3 regularization — a **7th explicit user input** (NOT derived from τ_max/s1, so the
constructor contract is unambiguous; the post-peak slope is set to dissipate `G_F^bond/lch`).
Drives Mode P's axial slot.
- **D4.1 Initial-slip regularization (must-fix):** the power-law `τ = τ_max(s/s1)^α` has
  `dτ/ds → ∞` at s→0 → singular first-iteration tangent. Use a **linear segment for
  `s < s0` (≈0.1·s1)** or a tangent cap `K_max`.
- **D4.2 Solution control (must-fix):** the descending branch has **negative tangent** →
  load control diverges. Past `τ_max` require **DisplacementControl / ArcLength / IMPLEX**
  (load control on the softening branch is unsupported; IMPLEX may reuse the ASDConcrete3D
  infrastructure).
- **D4.3 Objectivity (must-fix; retracts the earlier false claim):** node-lumped `L_trib`
  bond softening is **NOT** mesh-objective. Regularize the post-peak slope by a **bond
  fracture energy `G_F^bond`** scaled by characteristic length (`lch_ref/lch`, mirroring
  ASDConcrete3D), and add a **rebar-refinement convergence leg** to §6. v1 is first-order
  O(h²) (≳6–8 segments per development length); distributed Gauss bond → v2.
- **Scope:** v1 = well-confined pull-out; **cover-controlled splitting → v1.1** (xfail a
  small-cover splice test; do not silently over-predict bond).
- **Unit contract:** `F = τ · (πd) · L_trib` (perimeter × tributary length).

### D5 — Two-layer API; explicit `-mode` (no element-time integrator auto-detection)
- **Primitive:** `element('LadrunoEmbeddedRebar', tag, rebarNode, hostEleTag, ξ, η, ζ, '-mode', 'T'|'P', '-bond', matTag, '-perimeter', p, '-ltrib', L, '-kt', kt)`.
- **Generator (apeGmsh):** `g.reinforce(host=<set>, bars=<set>, bond='cebfip'|'perfect', mode='T'|'P')` — runs the guarded inverse-map at build time and emits primitives with an **explicit `-mode`** (apeGmsh knows the intended integrator). OpenSees elements cannot query the integrator at `setDomain`, so **no auto-detection**; instead `setDomain` adds a **consistency guard** that errors on integrator/mode mismatch (catches integrator swap + implicit↔explicit chaining). classTags 33005 / 33002 are collision-free (per-registry namespaces — objection refuted).

### D6 — Rebar element is USER-SELECTABLE: `truss` (default) or `beam` (opt-in)
The embedding couples **translations only** for both element types (`u_r = Σ N_i u_i^host`),
so element choice is orthogonal to the embedding machinery — only the extra beam DOFs and a
twist guard differ. This hands the axial-vs-flexural trade-off to the user (the Abaqus
`*EMBEDDED ELEMENT` model: constrain translations, leave rotations to the element).

- **`truss` (default) — `corotTruss`, axial-only.** Cheapest, explicit-clean, no rotational
  DOF/inertia. Covers the bulk; the bar's transverse stiffness comes from the Mode P
  penalty tie, not the element.
- **`beam` (opt-in, advanced) — `dispBeamColumn`/`forceBeamColumn`, translations-coupled /
  rotations-free.** Adds the bar's **distributed flexural/transverse stiffness** to the host
  (an *approximate* dowel/confinement contribution). Two caveats, both recorded as
  implementation requirements:
  - **Quantitative dowel is mesh-limited (not free).** Real dowel is a local bearing
    mechanism at ~3–6 `d_b`, far below typical solid mesh size; a beam slaved to a coarse
    host smears it ⇒ you get *a* transverse stiffness, **not** the right dowel capacity
    without a fine near-crack mesh or a dedicated dowel/bearing interface law (→ v2). Bending
    rotations are well-posed (the beam self-stiffens given imposed end translations).
  - **Bar-axis twist is a spurious zero-energy mode — MUST be stabilized.** The torsional
    sub-block `[[GJ/L,−GJ/L],[−GJ/L,GJ/L]]` is singular and the host has no rotation field to
    resist rigid twist about the bar axis ⇒ singular tangent. Fix: constrain the bar-axis
    rotation (SP) or add a small drilling stiffness. (Bending/translation DOFs are fine.)
  - **No double-count with bond-slip:** for `beam`, the embedding ties translations
    (perfect-bond transverse) and bond-slip acts **axially only** — the beam already supplies
    transverse stiffness; do not also apply a transverse bond spring.

Mode P transverse handling (applies to `truss`; for `beam` the element supplies it):
- Mode P transverse hold uses **volume-scaled** penalty (`k_t ≈ α·E_host·∛V`, mirror
  `ASDEmbeddedNodeElement`'s `iK·√V`) + the anisotropic `K` of D2 — a bare `corotTruss` is
  rank-1 (`EA/L·n⊗n`), so a bar parallel to a host face leaves a near-floating transverse
  DOF without this. Add a Zone-A **spurious-transverse-mode / tangent-rank** test.
- **Mode T restricted to small-rotation / linear hosts** (R5 promoted to a decision): the
  frozen reference-config `N_i(ξ)` drifts under large host rotation; corot/finite hosts
  **default to Mode P**. `setDomain` errors if Mode T meets a corot/finite host. (A `beam`
  rebar's own frame also goes stale under large host rotation — same restriction applies.)
- finite-strain `L_trib *= λ` (stretch) deferred; v1 documents the small-strain (host axial
  strain ≲ 2%) bond assumption.

## 4. Reference-code alignment (confinement claim CORRECTED)
- **Abaqus** `*EMBEDDED ELEMENT`: DOF-elim perfect bond, translations-only constraint (rotations of an embedded beam left free) → our Mode T + the D6 `beam` option; **no native bond-slip** → we add it.
- **LS-DYNA** `*CONSTRAINED_BEAM_IN_SOLID`: penalty + bond + mass scaling for explicit → our Mode P.
- **DIANA** bond-slip reinforcement: own-DOF rebar + 1D τ–s (CEB-FIP) → our Mode P + D4.
- **Confinement (corrected):** ASDConcrete3D *has* a Lubliner plastic-damage model with a Kc
  triaxial parameter and plastic strain ([ASDConcrete3DMaterial.h:27,90,342,370]) — so the
  triaxial surface advances faster under confining pressure. **BUT it has no tunable dilation
  angle, and the backbone peak stress is set by the user's compression curve, so whether
  confined `fcc` actually *rises* is uncertain and backbone-dependent — it is NOT a settled
  "emergent, no fc input" property.** The earlier claim (and the recipe-doc "dilation angle is
  the knob") was wrong: there is no dilation input. **This is resolved empirically by Gate 2
  (§6): a standalone triaxial-vs-Mander unit test, no rebar, run BEFORE any embedded-column
  test.** If `fcc` error >~10%, either calibrate the compression backbone per confinement
  level (Mander-wrapper workflow) or adopt a concrete model with a proper flow rule.

## 5. Implementation sketch (registration is a Definition-of-Done checklist)
- `SRC/element/ladrunoEmbeddedRebar/{LadrunoEmbeddedRebar.{cpp,h},OPS_LadrunoEmbeddedRebar.cpp,CMakeLists.txt}`.
- `SRC/material/uniaxial/LadrunoBondSlip.{cpp,h},OPS_LadrunoBondSlip.cpp`.
- `SRC/classTags.h`: add `ELE_TAG_LadrunoEmbeddedRebar 33005`, `MAT_TAG_LadrunoBondSlip 33002`.
- Register: forward-decl + `functionMap` entry in `OpenSeesElementCommands.cpp` (element) and
  `OpenSeesUniaxialMaterialCommands.cpp` (material); **broker** `case 33005`/`case 33002` in
  `FEM_ObjectBrokerAllClasses.cpp` (after LadrunoIMKBeam2d / LadrunoRebarBuckling groupings) +
  a **serialization round-trip test** (BezierTet10 precedent).
- Ship-time: `LEDGER_implementations.md` rows; `stamp_headers.py` on the new files;
  `banner_features.txt` line + `patch_banner.py`.

## 6. Validation plan — GATES first, then build
**Gate 1 (no element code) — embed-vs-conformal:** conformal column (apeGmsh layout +
`corotTruss` + `zeroLength` bond) at low/dense/mixed bar counts (implicit + explicit). Pass ⇒
conformal is viable for regular members (embed becomes a complex-cage tool); fail ⇒ measure
the infeasibility boundary, then build the element.
**Gate 2 (no rebar) — concrete confinement:** standalone triaxial compression (3D solid,
lateral p = 0 → 0.3·fc) through ASDConcrete3D vs Mander/Richart; quantify emergent `fcc` and
ductility. >~10% error ⇒ block the "emergent confinement" claim; calibrate or change model.
(LEDGER_quirks entry either way.)

Then, on the element:
1. **Inverse-map unit** — recover ξ for random points in straight-sided hex + tet; **distorted-hex / near-edge / outside-all-hosts** cases prove non-convergence is *caught*, not silently embedded.
2. **Single-bar pull-out (bond-slip)** — DisplacementControl; vs CEB-FIP backbone; first-step convergence at s≈0 (D4.1); **rebar-refinement objectivity leg** (D4.3).
3. **Lap-splice** — 3a well-confined (pass); **3b small-cover (xfail, splitting out of scope)**.
4. **Tension-stiffening RC prism** — Mode T vs Mode P head-to-head; document Mode T over-smears cracking (perfect bond into bulk GPs).
5. **Perfect-bond patch** — constant-strain on a **distorted** (non-parallelepiped) host (Mode T is exact only for affine; quantify error).
5.5 **Mode P penalty calibration** — sweep `K_p` 0.1–100·K_host; dt_cr decay; spurious-transverse-mode / tangent-rank test (bar ∥ host face); transverse-ringing under mass scaling.
6. **Implicit↔explicit consistency** — Mode T implicit vs Mode P explicit (mass-scaled) same quasi-static result; energy closure with the R1 mass-scaling audit.
7. **Confined column, cyclic** — Mode T then Mode P; confinement already certified by Gate 2.
8. **Large-host-rotation Mode-T drift** — quantify ξ-staleness; confirms the D6 Mode-T-linear-only restriction.

## 7. Risks (R5 promoted to a decision; R1/R2/R4 updated)
- **R1 (explicit/mass) — open, must measure:** does Mode P + rebar mass scaling preserve dynamics? Add a **per-element mass-scaling audit column** (γ) to the EnergyBalanceRecorder so a falsely-closed balance can't hide scaled-mass pollution; §6 item 6 checks it. Lagrange-bordered solve is the deferred diagonal-mass-preserving alternative.
- **R2 (constraint handler):** Transformation assumes a single MP per node — add pre-assembly validation rejecting (a) a retained host node that is itself constrained and (b) a constrained node referenced by >1 MP (rebar + rigidDiaphragm); both currently silent-wrong. Document the single-MP limitation.
- **R3 (inverse-map):** mitigated by the D3 guards + §6 item 1.
- **R4 (bond objectivity):** resolved by D4.3 fracture-energy regularization + §6 item 2 convergence leg.
- **R5 (corot consistency) → DECISION:** Mode T is small-rotation/linear-host only; corot/finite hosts use Mode P (D6).
- **R6 (scope):** converted to Gate 1.
- **R7 (dual-mode trap):** mitigated — explicit `-mode` (D5) + a loud API note when a mode is forced; the pre-build patch-test sweep decides single-vs-dual.
- **R8 (confinement calibration) — RESOLVED by Gate 2 ([[21_rc3d_validation_gates]], 2026-06-03):** ASDConcrete3D *does* emerge confined `fcc` within **5% of Mander** for p/fc∈[0,0.20] (and the ductility gain), via the Lubliner Kc surface — the "can't confine" blocker is refuted. Residual caveat: no tunable dilation-angle input, so the amount is set by Kc + the compression backbone (validate, don't pre-inflate fc). v1 still has no bar-restraint DOFs / interface-damage coupling.
- **Gates passed (see [[21_rc3d_validation_gates]]):** Gate 1 — conformal meshing (shared-node rebar, no new C++) is viable for **regular grid-aligned cages**; *measured* to converge at 4, 8 **and 12** longitudinal bars (pure count is not the wall). bars=4·n locks bars to mesh node-lines, so the embedded element's justification is **non-grid bar positions, ties/hoops, and the O(N²·nz) mesh-cost** — *inferred*, still wants a hoop/non-grid demonstration. Build order: ship the conformal recipe + bar-layout helper first; build `LadrunoEmbeddedRebar` for what conformal can't mesh. Both gates needed KrylovNewton + step-halving ⇒ the adaptive-stepping wrapper is a confirmed prerequisite.

---

## 8. Post-gate decision & v1 build order (this ADR's conclusion)

Both pre-build gates ran ([[21_rc3d_validation_gates]], 2026-06-03) and reshape the plan:
the concrete host is trustworthy for confinement (Gate 2: emergent `fcc` within 5% of
Mander), and conformal meshing already covers regular ≤8-bar members with **zero new C++**
(Gate 1). The embedded element is therefore **not the first thing to build**. Final ordering:

1. **v1a — Conformal RC recipe + parametric bar-layout helper (no new C++).** Generalize
   `rc3d_gates/gate1_conformal_column.py` into a documented recipe + an apeGmsh/openseespy
   helper that lays out longitudinal bars on mesh node-lines (`corotTruss` on shared nodes)
   for regular columns/beams/walls (≤8 bars). Covers the common 80%.
2. **v1b — Blessed adaptive-stepping wrapper (Python).** Both gates needed KrylovNewton +
   recursive step-halving; ship the tested wrapper proc (try-step → halve-on-fail →
   grow-on-success) so RC pushovers run unattended. Not a C++ integrator.
3. **v1c — Concrete confinement note + LEDGER_quirks.** Record Gate 2: emergent confinement
   is real via the Lubliner Kc surface (no fc inflation), but there is **no dilation-angle
   input** — calibrate Kc + the compression backbone; validate per the Gate-2 unit test.
4. **v2 — `LadrunoEmbeddedRebar` (33005) + `LadrunoBondSlip` (33002)** for the cases
   conformal can't mesh (non-grid bar positions, hoops/ties, beam-column joints — pure
   longitudinal count is not the trigger; 12 grid-bars converge conformally), built to the
   D1–D6 decisions with every must-fix guard (D2 setDomain explicit/Transformation guard;
   D3 guarded inverse-map; D4 initial-slip regularization + disp-control/IMPLEX +
   fracture-energy regularization, cyclic→later, splitting xfail; D5 explicit `-mode`;
   D6 volume-scaled anisotropic penalty + Mode T linear-host-only).

**Scope locked:** v1 = items 1–3 (no new C++ classes); v2 = item 4 (the two new classes).
The embedded element ships only after the conformal path is in users' hands and a dense-cage
case has demonstrably out-grown it.

## 9. SHIPPED so far (2026-06-03) + the next refactor (own session)

**Shipped & CI-green on `ladruno`:** `LadrunoBondSlip` (MAT 33002, PR #168) and
`LadrunoEmbeddedRebar` (ELE 33005, PR #169 + test-fix #171). The element is **Mode P**
(penalty coupling: perfect-bond-via-penalty + bond-slip); **Mode T deferred** (needs a
multi-retained constraint the stock `MP_Constraint` can't express). Current parse:
`element LadrunoEmbeddedRebar tag rebarNode nHost h1..hN -shape N1..NN -dir dx dy dz (-bond matTag [-bondScale bs] | -perfect kAxial) [-kt kt]`.

**Shape-fn-from-host refactor — IMPLEMENTED (2026-06-03, PR [#175](https://github.com/nmorabowen/OpenSees/pull/175), CI-only build):**
the user no longer has to re-supply host nodes or `-shape` weights — both can come from
the host *object*. OpenSees had **no shape-fn API on `Element`** (verified), so we added one.
Both steps landed entirely in the parser + base-class virtual; **`LadrunoEmbeddedRebar.{h,cpp}`
itself is UNCHANGED** (it still just stores `Nshape` — the host query happens in the parser,
keeping the element host-agnostic):

- **Step 1 — `-host eleTag` form (parser-only).** The parser resolves the host's nodes via
  `OPS_GetDomain()->getElement(eleTag)->getExternalNodes()` (the external-node list is complete
  at CONSTRUCTION, before `setDomain`, so this is a parse-time copy). Host must be defined before
  the rebar element. The explicit `nHost h1..hN` form is **kept** — Zone-A unit tests need a fake
  host (bare nodes) and arbitrary hosts need the escape hatch. Disambiguated by peeking the token
  after `rebarNode`: `-host` ⇒ host-element form, else `atoi` ⇒ explicit `nHost`.
- **Step 2 — `getInterpolationWeights` + `-xi` form.** Added
  `virtual int getInterpolationWeights(const Vector& xi, Vector& N)` to **`Element.{h,cpp}`**
  (default `return -1`; ledgered in `LEDGER_vanilla_files` + `// Ladruno` markers — vanilla
  base-class edit, vtable change, recompile-all but additive). Overridden on `LadrunoBrick`
  (trilinear `0.125·∏(1+ξ_I ξ)`, N size 8) and `BezierTet10` (quadratic Bernstein via the
  private `shapeFunctions`, N size 10). The parser's new `-xi x1..x_ndm` form (requires `-host`)
  calls `hostEle->getInterpolationWeights(xi, Nshape)` — single source of truth; errors with a
  "supply -shape" hint if the host returns −1. `-shape` still works for any host.
- **Grammar now:** `element LadrunoEmbeddedRebar tag rebarNode {nHost h1..hN | -host eleTag}
  {-shape N1..NN | -xi x1..x_ndm} -dir dx dy [dz] (-bond matTag [-bondScale bs] | -perfect kAxial) [-kt kt]`.
- Zone-A tests added (`test_ladrunoEmbeddedRebar_element.py`): real unit-cube `LadrunoBrick`
  host with `-host` (node auto-resolve), `-xi 0 0 0` (centroid → all 0.125), and an off-centroid
  `-xi` whose corner force split matches the trilinear formula (proves N really comes from the host).
- The **inverse-map** (global bar point → ξ) stays in the apeGmsh generator regardless — that
  is the irreducible point-location step; the only question is whether it emits `N` or `ξ`. With
  `-xi` it can now emit `ξ` and let the host element own the shape evaluation.

**CI/test discipline learned (apply next session):** ladruno auto-merge gates on the
classTag+manifest check but **NOT** on Zone-A pytest — a red test CAN land (it did, #169),
so watch the build and fix-forward. Run `ci/check_classtags.py` + `ci/check_manifest.py`
locally (pure-Python) before pushing — a new classTag needs a `testbed/manifest.yaml` row or
the fast gate fails. Element coupling tests must leave **≥1 free DOF** (don't fix-host AND
sp-all-rebar — zero free DOFs aborts the linear solver, exit 255). New C++ verifies via the
CI Zone-A build (the worktree isn't locally build-configured).
