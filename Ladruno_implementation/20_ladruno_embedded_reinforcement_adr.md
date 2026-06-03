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
> **STATUS UPDATE (2026-06-03): Mode T is DEFERRED INDEFINITELY.** The two-mode framing
> below is kept for the record, but the active design has moved on: Mode T (DOF-elimination
> via a multi-retained `MP_Constraint`) is parked indefinitely — it is implicit-only,
> perfect-bond-only, and needs a constraint object the stock `MP_Constraint` can't express
> (see §9). All forward development is on the **Mode P family**, now generalized to an
> **`-enforce {penalty|al|nitsche}`** flag (see **§10**). `transformation` survives only as a
> future, unscheduled `-enforce` value. Where this ADR says "Mode T / Mode P", read
> "`-enforce transformation` (deferred) / `-enforce penalty` (active)".

`LadrunoEmbeddedRebar` (classTag **33005**, Element band) exposes:
- **Mode T (transformation / perfect bond, implicit-only) — DEFERRED INDEFINITELY (§10):**
  rebar node translational DOFs slaved to host interpolation `u_r = Σ_i N_i(ξ_r) u_i^host`
  via an `MP_Constraint` resolved by the **Transformation** handler → rebar DOFs eliminated.
  No penalty. Parked: implicit-only + perfect-bond-only + needs a multi-retained MPC.
- **Mode P (penalty / bond-slip, explicit-capable + implicit) — the active path:** rebar keeps own DOFs + mass;
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
  R1 audit below. **CORRECTED in §10.6:** this line is inaccurate as-built — the element
  returns zero mass ([LadrunoEmbeddedRebar.cpp:262]) so the *per-element* `CriticalTimeStep`
  scan skips it entirely (it never sees `K_p`); the real limit lives only in the assembled
  coupled system and is currently a **false-safe** (collapse can exceed 30×). §10.6 specs the
  bipenalty fix that makes it visible + tunable.
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
(penalty coupling: perfect-bond-via-penalty + bond-slip); **Mode T deferred INDEFINITELY**
(see §10 — implicit-only, perfect-bond-only, needs a multi-retained constraint the stock
`MP_Constraint` can't express; all forward work is the Mode P `-enforce` family). Current parse:
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


---

## 10. Mode-P improvement roadmap & the `-enforce` generalization

> **Status (2026-06-03):** `LadrunoEmbeddedRebar` (ELE 33005) ships **Mode P** only —
> penalty coupling, perfect-bond-via-penalty + axial bond-slip (`LadrunoBondSlip`, MAT 33002).
> This section is the forward plan: it replaces the legacy `-mode {T|P}` switch with an
> `-enforce {penalty|al|nitsche|transformation}` flag and lays out the kernels, the shared
> infrastructure they inherit, the explicit-`dt_cr` story, and the large-rotation fix, in
> sequenced build order.

### 10.1 The `-enforce` flag and the two-bucket architecture

The constraint is one and the same kinematic tie at every enforcement level:
`c(u) = g = u_r − Σ_i N_i(ξ) u_i^host = 0`, with the shipped discrete operator
`B = [ I | −N_1 I … −N_M I ]`, gap split `s = g·dir`, `g_t = g − s·dir`, traction
`t = F_axial(s)·dir + k_t·g_t`, tangent `K = BᵀD B`, `D = k_axial·dir⊗dir + k_t·(I−dir⊗dir)`
(`LadrunoEmbeddedRebar.cpp:181-235`). Penalty, augmented-Lagrangian (AL), and Nitsche are the
three classical Lagrangian treatments of *that one* tie (Wriggers 2006, *Computational Contact
Mechanics* Ch.6; Belytschko et al. 2014, *Nonlinear FE* §6.8); they differ only in the per-step
`(t, D)` kernel inside `formBandTraction` plus a small amount of extra strategy state.

This yields a clean **two-bucket split**:

- **Bucket A — own-DOF additive family `{penalty, al, nitsche}`.** ONE element, one strategy
  enum chosen in the constructor. All three share: external nodes, `getNumDOF`, zero mass
  (`getMass`→`M0≡0`, explicit-capable), the B-matrix kinematics, the volume-scaled penalty
  auto-scale, co-rotation of `dir`, serialization, responses, and — critically — they inherit
  `dt_cr` (`CriticalTimeStep.cpp:205/226` calls only `getMass`/`getTangentStiff`) and energy
  accounting (`EnergyBalanceKernel.h` consumes only `getResistingForce`/`getMass`) **for free**.
- **Bucket B — DOF-elimination `{transformation}`.** A genuinely different object family: a
  multi-retained `MP_Constraint` resolved by the Transformation handler, eliminating the rebar
  DOFs, with **no element K**. Implicit-only (it densifies `TᵀMT` → wrong inertia in
  `DiagonalSOE`, D2/ADR:90), perfect-bond-only, and requires a multi-retained constraint object
  the stock single-retained `MP_Constraint` cannot express (§9).

**Decision.** `-enforce` supersedes `-mode`: `penalty ← P`, `transformation ← T`
(**deferred indefinitely** — it survives only as a future, unscheduled flag value), plus new
bucket-A kernels `al` and `nitsche`. `-mode P`/`-mode T` become deprecated aliases.
`-enforce transformation` is a **hard parse-time error** ("deferred indefinitely, ADR 20 §10;
use penalty|al") until the bucket-B builder exists. The setDomain explicit-integrator guard
generalizes from "`mode==T && explicit`" to "`strategy ∈ bucket B && explicit`"; bucket A needs
no guard (all stiffness-only, mass diagonal). `-perfect`/`-bond` remain orthogonal — they select
the **axial law**, while `-enforce` selects the **constraint treatment**.

### 10.2 Shared Mode-P infrastructure (highest leverage — build first)

All of penalty/al/nitsche inherit these. This is the highest-value work because it is written
once and pays off three times.

**(a) Auto-scaled penalty.** Replace the raw `kt` default with a stiffness-matched
`k_t = α · scale · lch`, the dimensionally-correct generalization of ASD's geometric
`iK = m_K·√V` (`ASDEmbeddedNodeElement.cpp:849/1296`). A continuum host contributes a nodal
stiffness `~E·lch`, so the penalty should be a dimensionless multiple of it (Wriggers §6.3;
de Souza Neto et al. 2008 penalty regularization): constraint residual `‖g‖ ~ t/k_t = O(σ/(αE))`
shrinks while `κ(K) ~ α·κ(K_host)` stays bounded. **Recommended host-agnostic scale:**
`lch = host->getCharacteristicLength()` (LadrunoBrick `∛V` at `LadrunoBrick.cpp:3190`,
BezierTet10 `∛6V`; base `Element::getCharacteristicLength` min-edge fallback) and
`scale = ‖host->getInitialStiff()‖_∞ / lch` — self-calibrates to the host's actual tangent,
needs no user `E`, evaluated once at `setDomain`. Default `α = 1e3` (‖g‖≈1e-3, κ≤1e10);
floor `α≳10` (below it the bar floats transversely — D6 rank-1 note); ceiling `α≲1e4` for
implicit Krylov solvers. **Only the transverse `k_t` is auto-scaled** — bond-slip axial comes
from `bondMat->getInitialTangent()·bondScale` (a physical stiffness, untouched). *Wiring (the
§9 host path makes this feasible now):* store `hostEleTag` so the host can be re-fetched at
`setDomain` (`OPS_LadrunoEmbeddedRebar.cpp:101` already resolves it but discards it); add a
`-kt auto [α]` sentinel; explicit `nHost` form with no host falls back to own-node min-edge lch
and should pass a numeric `kt`. **Loud documentation note:** ASD's `m_K=1e18` is *not* E-scaled —
porting users must NOT pass `1e18` as `α` (κ blow-up). *Effort: moderate.*

**(b) Energy accounting.** The recorder already absorbs coupling work *blindly*:
`IE += Fᵀv = tᵀġ = F_axial·ṡ + k_t·g_tᵀġ_t` (orthogonal split, `dir ⊥ (I−dir⊗dir)`). The defect
is classification — the artificial penalty work `½k_t‖g_t‖²` (and, in perfect-bond, the entire
`½k_axial s²`) is booked as physical internal energy, silently inflating IE and producing a
*falsely-closed* balance (R1's concern has a penalty twin). **Fix (moderate, openseespy-native,
no vanilla edit):** the element reports `eleResponse` channels — `penaltyEnergy` (artificial),
`bondEnergy`/`bondDissipation` (physical, the τ–s hysteresis `∮τds`, single-sourced from a new
`LadrunoBondSlip` `dissipatedEnergy` response), `constraintViolation` (`‖g_t‖`) — accumulated in
`commitState()` with the recorder's trapezoidal rule; users bin embedded elements in their own
`MeshRegion` and subtract `Σ penaltyEnergy` for the net physical balance. *Caveat:* the
element's committed-increment integral matches the recorder's `Fᵀv` only as `Δt→0`; decide
whether `penaltyEnergy` is a *diagnostic* or a *subtractand* and match conventions accordingly.
**Future-proofing (vanilla, defer until a 2nd kernel lands):** promote to
`virtual double Element::getArtificialEnergyRate(){return 0;}` so AL/Nitsche (each with its own
artificial split — Nitsche's stabilization term has one too) net uniformly in the kernel. References:
Wriggers Ch.6 (penalty energy = constraint-violation energy → 0 only as `k→∞`); Belytschko
§6 (artificial-energy fraction <5–10% as a quality metric); LS-DYNA SLEOUT/GLSTAT interface-energy
precedent. *Effort: moderate.*

**(c) Co-rotated frame.** *Required before AL/Nitsche*, else they inherit the same staleness.
The frozen `dir` (set once from reference config, `LadrunoEmbeddedRebar.cpp:61-65`) is the only
true large-rotation defect — see §10.5. Shared upgrade living in `formBandTraction`. *Effort: moderate.*

### 10.3 Enforcement kernel — penalty hardening (`-enforce penalty`)

The shipped kernel. "Hardening" here means finishing the auto-scale (10.2a) so the default
`k_t` is physically conditioned rather than a raw `1e12`, and exposing a `penalty`/`kt`
`setResponse` reporting the resolved value + cap-active flag (for the §6 item 5.5 sweep). No new
mechanics. Hook points: `formBandTraction` (`:145`), `getTangentStiff`/`getResistingForce`
(`:206/:181`), `setResponse` (`:387`). References: Wriggers §6.3; de Souza Neto §6. *Effort: quick-win
(rides on 10.2a).*

### 10.4 Enforcement kernel — augmented Lagrangian (`-enforce al`)

**Mechanics.** Add a per-element multiplier `λ ∈ R^ndm` (committed state): traction becomes
`t = λ + D·g` (Simo & Laursen 1992; Wriggers §6.4). The element tangent is **identical** to
penalty (`K = BᵀDB`) — `λ` is constant within an inner solve, entering only the residual as
`Bᵀλ`. The Uzawa update `λ ← λ + D·g` fires **once per converged step inside `commitState()`**
(the natural per-step variant — `Element::commitState` is driver-called exactly once
post-convergence, `:127`), needing **no analysis-core change**. The converged gap contracts
geometrically `g_{k+1}=ρ g_k`, `ρ = K_struct/(K_struct+K_p)<1` (Bertsekas 1982), so **`g→0` at
*moderate* `K_p`** — no `K_p→∞`. This is the headline win: AL reaches near-exact constraint
satisfaction at `K_p ≈ K_host` (exactly the D2 `K_p ≲ 10·K_host` sweet spot), so it (i) dodges
the unproven `K~1e18` conditioning risk, (ii) *enlarges* explicit `dt_cr` versus a stiff
perfect-bond penalty, and (iii) drives the artificial penalty energy →0, fixing the 10.2b IE
pollution in the converged state.

**Bond-slip interaction (clean split):** augment **only the transverse projection** `(I−dir⊗dir)g`;
leave the axial slot to `bondMat` so the intended τ–s dissipation is untouched.

**Hook points.** `commitState` (`:127`, the Uzawa point), `getResistingForce` (`:181`,
`t(k)+=λ(k)`), tangent **unchanged**, `revertToStart/LastCommit` (`:139/:134`, reset `λ`),
`sendSelf/recvSelf` (`:271-373`, grow payload by ndm + the BezierTet10 round-trip test),
`setResponse` ("augLambda"/"gap"). Per-step Uzawa is invisible to the inner `CTest` (λ frozen
within a step); the small first-iteration residual jump at λ-update is absorbed by KrylovNewton
(the ADR's blessed solver). References: Simo & Laursen (1992) *Comput. Struct.* 42:97; Wriggers
§6.4.2; Bertsekas (1982) Ch.2; Laursen (2002) §3.3 (explicit-compatible, mass stays diagonal).
**This is the strategic next kernel** — high value, low risk, no missing infrastructure.
*Effort: moderate.*

### 10.5 Large-rotation resolution (crisp)

The "frozen `N_i`" worry is a **red herring**: `N_i(ξ)` at a fixed *material* `ξ` is
deformation-independent (de Souza Neto §4, isoparametric map of a fixed parametric point), and
`g = u_r − Σ N_i u_host` is a current-config relative displacement — under rigid host rotation
`Q`, both terms transform identically so `g → Qg` (objective). The **coupling/tie is therefore
already exactly objective**, and Mode P inherits **no MPC-linearization error** (the geometric
nonlinearity is carried exactly to penalty tolerance every step) — this is precisely why Mode T
is the one with the small-rotation restriction, *not* Mode P.

**What actually drifts:** the *constitutive split* against the frozen `dir`. With `dir` fixed,
`s = (Qg)·dir_0 ≠ g·dir_0`, so a pure rigid rotation registers spurious axial slip and a bogus,
non-objective traction `t` (it does not rotate as `Qt`; violates frame indifference, Crisfield
1997 Vol.2 Ch.17). **Cure:** co-rotate `dir` (and only `dir`) per step from current host
geometry, then `s = g·dir(t)`, `g_t = g − s·dir(t)`, `t → Qt` is restored. **Recommended v1
path (Option i, host-agnostic, no vanilla edit):** store two natural coords `ξ_a, ξ_b`
straddling the bar, recompute `dir(t) = normalize( Σ N_i(ξ_b)x_i − Σ N_i(ξ_a)x_i )` from current
node positions via the existing `getInterpolationWeights` called twice — captures host rigid
rotation exactly (secant smears curvature, negligible for a bar within one host element). The
principled alternative (Option ii, `dir = R_host·dir_0` via a polar decomposition of `F_host`)
needs a NEW `getInterpolationGradients`/`getDeformationGradient(ξ)` Element virtual — research-grade,
re-opens the §9 base-class surface.

**Verdict:** with the co-rotated `dir`, **the small-ROTATION restriction can be dropped** for
Mode P (default the corot flag OFF for bit-identical back-compat, ON for finite-rotation use). The
small-STRETCH `L_trib *= λ` correction (D6) stays deferred — co-rotation must not be oversold as
full finite-strain support. Mass stays zero, so the fix is explicit-safe. Implicit consistent
tangent gains a `∂dir/∂u` geometric (generally non-symmetric, follower-type) term; **drop it in
v1** (standard EICR practice, exact for explicit, converges for moderate per-step rotation under
the v1b step-halving wrapper). References: Crisfield Vol.2 Ch.17; Belytschko §4 (`F = Σ x_i⊗∇_X N_i`);
de Souza Neto §3 (objectivity). *Effort: moderate.* New LEDGER_quirk: anisotropic embedded coupling
needs a co-rotated axis; ASD's isotropic `√V` penalty does not — which is why ASD recomputes from
*reference* coords and is still objective.

### 10.6 Explicit `dt_cr` (correct the framing; bipenalty vs localized mass scaling)

**The ADR D2 claim that "penalty `K_p` enters the coupled rebar-host eigensolve at
`CriticalTimeStep.cpp:226`" is FALSE as built.** `computeCriticalTimeStep` is a strictly
**per-element** scan (`:200`): it reads `getMass()` first (`:205`), and because the element
returns `M0≡0` (`:262`), `mPositive=false` (`:91`) → DSYGVX skipped, DGGEV runs with `B=0` →
every `β≈0` → all eigenpairs rejected → `λ_max=−1` → the element is **silently non-contributing**
(`:246`). The host element's per-element tangent does not contain the coupling K (it lives on the
rebar element, assembled globally), so the host scan misses it too. **Net: the penalty mode is
invisible to the reported `dt_cr` — a dangerous false-safe.** The real limit lives only in the
assembled coupled system: `ω²_penalty ≈ k_eff/m_min`, `dt_penalty = 2√(m_min/k_eff)`, with
`m_min = min(m_rebar, min_i N_i·m_host,i)`. For `k_eff=10·k_host` and a light steel bar
(`m_rebar/m_host~0.01`) the collapse can **exceed 30×**, well past the ADR's "2–10×". This must
be stated plainly: `-cfl` users can pick an unstable step today.

**Two fixes:**

- **Bipenalty (Askes & Hetherington 2010, *Proc. R. Soc. A* 467:1205; Caleyron et al. 2013) —
  recommended.** Pair the stiffness penalty `k_p` with a *mass* penalty `m_p = k_p/(β·ω_host)²`
  added to the **rebar block of `getMass()` only** (lump the full `m_p` on the slave node — the
  faithful `m_p·BᵀB` has host off-diagonal `N_iN_j` blocks that **would** densify M and corrupt
  `DiagonalSOE`; the lumped form does NOT — this is a hard implementation constraint, not a
  gloss). This bounds the spurious constraint frequency to `β·ω_host` (`β~2–5`), so `dt_cr` is
  set by a chosen physical ratio, *decoupled from* `k_p`. Bonus: with a nonzero rebar-block
  diagonal, the DGGEV fallback (`CriticalTimeStep.cpp:137-179`) now returns
  `λ_max = k_eff/m_p = (β·ω_host)²` → the element **self-reports** `dt = 2/(β·ω_host)`, turning
  the invisible instability into a visible, user-tunable bound. Gate on `-enforce penalty`
  (set `m_p=0` for `al`; reuse for `nitsche`).
- **Localized mass scaling + γ audit (R1, the fallback).** Scale the *rebar element* mass to
  recover the step. Inferior to bipenalty: it pollutes *all* rebar modes (real axial-wave
  dynamics), needs the EnergyBalanceRecorder γ-column to certify added/true mass stays small,
  and the γ has no natural home on a zero-mass coupling element (it lives on the corotTruss).

**Open `dt_cr` audit items:** (i) where does `ω_host` come from? — user input, or
`~c_dilatational/∛V` approximation; (ii) confirm the DGGEV `λ_max` selection picks the
finite-β rebar mode and is not masked by host-DOF zero-β pairs (needs a 2-node toy-pencil unit
test); (iii) multiple embeds on one rebar node would double-count `m_p` — needs per-node
de-duplication; (iv) Rayleigh `β_K` on the bounded penalty mode could still shrink the damped
`dt` (`:249-255`) — decide whether the element zeroes its own Rayleigh factors.
References: Belytschko §6.7; Hughes (1987) Ch.9; Hetherington & Askes (2009). *Effort: moderate.*

### 10.7 Enforcement kernel — Nitsche (`-enforce nitsche`)

**Mechanics.** Replace the arbitrary penalty spring with a variationally-consistent flux:
`t = −⟨t_h⟩ + γ g`, where the consistency term `t_h·dir = (σ_h·m̂)·dir` is the host stress at
ξ resolved on the bar axis (Nitsche 1971; Hansbo & Hansbo 2002; Annavarapu, Hautefeuille & Dolbow
2012). Three terms: consistency `−∫(t_h·dir)δs` (removes the `1/γ` accuracy error), adjoint/symmetry
`−∫(t_h(δu)·dir)s` (symmetric variant → optimal L2, symmetric tangent), stabilization `+∫γ sδs`
(the demoted penalty). Accuracy becomes **γ-independent for any `γ ≥ γ_min ~ C·E_host/h`** (an
inverse estimate — the same `√V/∛V` h-scaling the fork already ships, just without the consistency
terms); `γ` only buys coercivity.

**Feasibility wall (the reason for the rating).** Two pieces are **missing** in OpenSees:
(1) no `getInterpolationGradients`/∂N analogue to the §9 `getInterpolationWeights` — needed for
`ε_h(ξ) = Σ ∂N_i⊗u_i`; (2) no API to evaluate host material stress/tangent at an *arbitrary* ξ
(the cached-Response idiom — LadrunoBrick `damageResponse`, `:287-294` — only reaches Gauss-point
slots). Phase-0 cheap route: sample the nearest GP (O(h) consistency error). Rigorous route: a new
`virtual int Element::getInterfaceFlux(ξ, dir, traction, dTr_dUhost)` overridden per host —
**re-coupling the element to host types**, breaking the host-agnostic stance the fork deliberately
adopted. Additional risks: codim-2 mismatch (a bar is codim-2, classical Nitsche is codim-1 — the
"flux" and the inverse constants are non-standard, **little published precedent**, must be
*validated not cited*); softening host breaks constant `γ_min` (ASDConcrete3D's degrading `D_h`
→ tangent-tracking `γ` à la Annavarapu weighted flux); symmetric-vs-nonsymmetric trade (some
OpenSees solver paths assume symmetry). The transverse block has essentially no precedent — likely
keep it pure penalty/AL and apply Nitsche to the axial bond only. References as above + Burman &
Hansbo (2012) ghost penalty for bar-near-face degenerate cuts. **Benchmark AL first — it achieves
the same γ-independence with *no* host-stress/∂N query.** *Effort: research-grade.*

### 10.8 Sequenced build order

1. **Quick-win — penalty hardening + auto-scale (§10.2a, §10.3). ✅ SHIPPED (PR #177).**
   `-kt auto [-ktAlpha α]` resolves `kt = α·max|K_host(i,i)|` (≈ α·E·lch) **lazily on first
   assembly** (not `setDomain` — dodges element setDomain-ordering) from the §9 host pointer;
   `eleResponse "kt"` exposes the resolved value. Default `α = 1e3`. Opt-in (numeric `-kt`
   unchanged), so existing models are bit-identical.
2. **Quick-win — energy channels (§10.2b). ✅ SHIPPED (PR #177).** `LadrunoBondSlip` `energy`
   /`dissipatedEnergy` response (cumulative ∫τ ds, trapezoidal, committed + serialized);
   element `penaltyEnergy` (artificial ½kt|gt|²[+½k s² perfect-bond]), `constraintViolation`
   (|gt|), `bondEnergy` (= bondScale·material work, single-sourced + bond-law-agnostic via a
   cached material sub-response). Region-net post-processor left to the user (openseespy-side).
3. **Strategic — co-rotated `dir` (§10.2c/§10.5). ✅ SHIPPED (PR #180).** `-corot` recomputes
   the bar axis each step as the secant from the embed point to a point B along the bar
   (`-xiB`/`-shapeB`), using current host node positions, so the axial/transverse split stays
   frame-objective under large host rotation. Default OFF (frozen `-dir`, bit-identical).
   `eleResponse "dir"` reports the working axis; v1 omits the `∂dir/∂u` tangent term (EICR).
   Headline test green: rigid host z-rotation → axis follows `Q·dir0` with `-corot`, frozen
   without. **This lifts the small-ROTATION restriction for Mode P** (small-STRETCH `L_trib·=λ`
   still deferred). New LEDGER_quirk: anisotropic embedded coupling needs a co-rotated axis;
   ASD's isotropic penalty does not.
4. **Strategic — `-enforce` flag + AL kernel (§10.1, §10.4). ✅ SHIPPED (PR #181).**
   `-enforce {penalty|al}` (default penalty); `nitsche` rejected (not built), `transformation`
   rejected (deferred indefinitely). AL adds a per-element multiplier `lambda` to the traction
   (`t = penalty/bond traction + lambda`) with the **tangent unchanged** (`K = BᵀDB`; lambda
   constant within an inner solve) and a per-step Uzawa update `lambda += Δλ` in `commitState`
   — `Δλ` = the PERFECT-BOND penalty traction only (transverse `kt·g_t` always; axial
   `kAxialPerfect·s` only when there is no bond law), so bond-slip's physical τ–s axial is
   never driven to zero. Inherits the co-rotated `dirCur` from item 3. Responses `augLambda`
   /`gap`; `lambda` serialized. Tests green: AL drives the perfect-bond gap → 0 at moderate
   `kt` (the multiplier carries the load) where penalty leaves `P/kt`; AL leaves bond-slip
   axial untouched; `-enforce penalty` bit-identical to default. No analysis-core change
   (commitState is driver-called once per converged step). **This is the near-exact,
   well-conditioned, explicit-safe perfect-bond path** — the strategic payoff of the roadmap.
   *Review fix:* for **bond-slip + `-corot` + `al`**, the transverse multiplier is re-projected
   onto the current transverse plane each step so a rotating `dirCur` can't leak it into a
   spurious axial force on the τ–s slot; perfect-bond keeps the full 3D multiplier. Open (→§10.9):
   AL+corot is exact only to O(per-step rotation) like the dropped ∂dir/∂u term, and per-step
   Uzawa under **cyclic** load is path-approximate (v1 is monotonic) — both want a v2 test leg.
5. **Strategic — bipenalty `dt_cr` (§10.6).** Rebar-block lumped `m_p`, `-wcap β`, self-reported
   `dt_cr` response. Gated on `-enforce penalty`. *Independent of 3/4 but pairs naturally with
   explicit AL runs.* Correct the D2 framing in this ADR as part of the PR.
6. **Research-grade — Nitsche (§10.7), separate PR.** Phase-0 `getInterpolationGradients` vanilla
   virtual (ledgered, mirrors §9) → Phase-1 host-stress query (cheap GP-sample first) → Phase-2
   kernel behind the flag, `xfail`/`skip` until the host-flux API lands. *Justify against AL on
   conditioning + `dt_cr`, not against pure penalty.*

### 10.9 Open questions / validation legs

- **AL convergence rate:** is per-step (single) Uzawa's `g→0` rate adequate for an RC pull-out
  *end-state*, or is a thin Python in-step re-augmentation wrapper needed? (Per-step first; add
  only if a gate shows residual end-of-step gap.)
- **Auto-scale on a force-based/condensed host:** confirm `getInitialStiff` is valid pre-analysis
  and non-zero-norm (LadrunoBrick `:486` and BezierTet10 are fine); guard a zero-norm return.
- **Nodal-mass `dt_cr` cap ordering:** the rebar corotTruss mass may not be on the node at the
  embedded element's `setDomain` — prefer the explicit `-wcap`/`-dtcr` budget for v1; lazy
  first-`update()` gathering is a v2 refinement.
- **Bipenalty:** lumping-induced constraint-mode-shape error must not leak into real low-frequency
  bar dynamics (Zone-A check); DGGEV `λ_max` selection toy-pencil test; multi-embed de-dup.
- **Co-rotation:** confirm `getTrialDisp()` is populated at the `CentralDifferenceLadruno`
  predictor/corrector call order; degenerate-axis (`|x_b−x_a|→0`) fallback-to-`dir_0` + warn.
- **Nitsche:** γ-independence sweep `γ∈[γ_min, 100γ_min]` (the headline claim, penalty fails by
  construction); γ_min coercivity vs `E_host/h`; softening-host robustness; head-to-head conditioning
  + `dt_cr` vs penalty and AL; ghost-penalty bar-near-face case.
- **Energy:** decide whether `penaltyEnergy` is a recorder *subtractand* (must match `Fᵀv`
  convention) or an independent diagnostic; finite-`k` caveat (perfect-bond residual `s` is partly
  real elastic bond compliance, only fully "artificial" as `k→∞`).
- **Whether Nitsche is worth re-coupling to host types at all** — if AL covers the consistency /
  conditioning need, Nitsche may join `transformation` in permanent deferral, leaving
  `-enforce {penalty|al}` as the real bucket A.
