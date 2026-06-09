---
title: "ADR 30 — LadrunoProjectionHandler: momentum-conserving explicit constraint projection (DYNA-style)"
project: Ladruno
type: ADR / implementation spec
status: SCOPED — formulation locked on paper, no code; Phase 0 falsification tests defined
related:
  - "[[05_robust_central_difference]]"                # the host integrator (CentralDifferenceLadruno, 33003)
  - "[[24_ladruno_coupling_constraints_adr]]"         # coupling family; D3b mass-condensation idea subsumed here
  - "[[28_ladruno_distributing_coupling_rbe3_adr]]"   # RBE3 element (penalty/AL) — the element-route sibling
  - "[[29_ladruno_kinematic_coupling_rbe2_adr]]"      # RBE2 element (penalty/AL/bipenalty) — ditto
  - "[[23_ladruno_embedded_node_adr]]"                # penalty embed whose dt_cr pain motivated this
  - "[[LEDGER_implementations]]"
tags: [adr, constraints, explicit-dynamics, projection, momentum-conserving, central-difference, ls-dyna, abaqus, constraint-handler]
updated: 2026-06-09
---

# ADR 30 — `LadrunoProjectionHandler` (explicit constraint projection)

**Status:** scoped, no code. This ADR locks the formulation, the architecture seams, the
v1 scope, and the validation battery. It is the handler-level realization of the
"LS-DYNA-style rigid-body / tied-constraint integrator" idea parked as **D3b** in
[[24_ladruno_coupling_constraints_adr|ADR 24]].

> [!info] What this is, in one line
> A new constraint handler + a small hook in `CentralDifferenceLadruno` that enforces
> `MP_Constraint`s in explicit dynamics by **mass-orthogonal projection of the solved
> accelerations** — exact, momentum-conserving, Δt-neutral — instead of penalty springs
> (which eat `dt_cr`) or equation-level condensation (which breaks the diagonal-mass
> recipe). = how LS-DYNA enforces `*CONSTRAINED_LINEAR` (Theory §28.2, `ü_c = [LᵀML]⁻¹LᵀF`),
> tied contact (Theory §26.9, mass/force redistribution), and nodal constraint sets
> (Theory §28.1, `a = ΣMᵢaᵢ/ΣMᵢ`). = the kinematic-contact half of Abaqus/Explicit.

---

## 1. Driver & goal

The fork's two flagship tracks **fight each other** today:

1. **Explicit dynamics** (`CentralDifferenceLadruno` 33003, ExplicitBathe/LNVD): the whole
   recipe is *lumped mass → `system Diagonal` → `algorithm Linear` → dt < dt_cr*. One
   trivial diagonal solve per step.
2. **Constraints/coupling** (rigidDiaphragm, equalDOF, rigidLink, RBE2/RBE3 elements,
   embedded nodes): every available enforcement path damages recipe #1:

   - **Penalty** (handler or element): adds `α·CᵀC` stiffness on the *full* space →
     raises `ω_max` → collapses `dt_cr`. This is the documented pain of
     [[23_ladruno_embedded_node_adr]] ("penalty embed kills explicit dt") and the reason
     RBE2/RBE3 grew the bipenalty escape hatch. A penalty soft enough to leave `dt_cr`
     alone is too soft to enforce a *rigid* tie — LS-DYNA itself only uses Δt-bounded
     penalty for *joints* (compliant by nature, Theory Eq. 25.30), never for ties.
   - **Lagrange**: zero diagonal block → impossible with a diagonal SOE. Dead on arrival
     for explicit.
   - **Transformation**: kinematically exact, but `TᵀMT` puts **off-diagonal mass
     coupling** between the retained DOFs of a constraint (the transport/rotary terms of
     a rigid link or diaphragm). A diagonal SOE silently *drops* those off-diagonals at
     assembly → **wrong inertia, wrong physics, no warning** (to be proven in Phase 0,
     test T6). With a banded solver it is correct but forfeits the O(n) diagonal solve.
     Plus the standing upstream defects: one-MP-per-constrained-node silently corrupts the
     DOF map, segfault at zero free DOFs (LEDGER_quirks), no MP chains.

**Goal:** constrained explicit dynamics that is simultaneously (i) **exact** (rigid ties
truly rigid), (ii) **Δt-neutral** (constraints can only *lower* `ω_max`, never raise it),
(iii) **momentum-conserving**, (iv) compatible with the untouched `Diagonal` SOE recipe,
and (v) loudly diagnostic where upstream is silently wrong.

**Non-driver:** implicit analyses. Penalty / Lagrange / Transformation / Auto remain the
right tools there; this handler refuses non-explicit integrators (§5.3).

---

## 2. Formulation (locked)

### 2.1 The projection

Explicit CD solves `a_raw = M⁻¹ r` with diagonal `M`. A set of linear constraints
restricts admissible accelerations to a subspace: stack each constraint group's relation
`u_c = C u_r` (twice differentiated: `a_c = C a_r`, valid because `C` is constant — §6
non-goal for time-varying `C`) and write the group's full DOF vector through the
retained ones:

```
a_group = L a_r ,      L = [ I ]    (retained rows)
                           [ C ]    (constrained rows)
```

The enforced acceleration is the **M-orthogonal projection** of the raw one onto that
subspace — i.e. minimize `‖a − a_raw‖²_M` subject to `a = L a_r`:

```
a_r = (Lᵀ M_g L)⁻¹ Lᵀ M_g a_raw,g          (dense solve, size = n_retained of the group)
a_c = C a_r
```

with `M_g = diag` of the group's lumped masses. Untouched DOFs keep `a_raw`. Special
cases recover the LS-DYNA forms exactly:

- `equalDOF` (C = I, one retained DOF): `a = (m_r a_r + m_c a_c)/(m_r + m_c)` — the
  mass-weighted common acceleration of `*CONSTRAINED_NODE_SET` (Theory Eq. 28.1).
- general linear MP: `[LᵀML]⁻¹Lᵀ(M a_raw) = [LᵀML]⁻¹LᵀF` — `*CONSTRAINED_LINEAR`
  (Theory §28.2), since `M a_raw = r` is the assembled force.

**Constraint force recovery** (free byproduct): `f_tie = M (a_raw − a_proj)` per DOF —
exposed as a query in Phase 3.

### 2.2 Properties (the load-bearing claims)

- **Momentum:** `Lᵀ M (a_proj,g − a_raw,g) = 0` by construction — the generalized
  momentum conjugate to every retained direction is conserved; for translational ties
  this is literal linear-momentum conservation.
- **Exactness & drift:** the leap-frog updates are linear (`v_{n+1/2} = v_{n−1/2} + Δt·aₙ`,
  `u_{n+1} = uₙ + Δt·v_{n+1/2}`), so if `u₀, v₀` satisfy the constraints and every
  `a` (including the starter `a₀`) is projected, then **u and v remain on the constraint
  manifold to machine precision for all time** — no drift term exists for constant `C`.
  Initial-condition compliance is *checked* at `domainChanged()` (warn + optional
  projection of `u₀, v₀`, `-projectICs`).
- **Stability (Δt-neutrality):** projection restricts the dynamics to a subspace; by
  eigenvalue interlacing (Rayleigh quotient on a subspace), `ω_max(constrained) ≤
  ω_max(unconstrained)` — so the existing per-element `dt_cr` estimate from
  `CriticalTimeStep` stays **conservative**. Contrast penalty, which adds stiffness on
  the full space and *raises* `ω_max`. Verified empirically by T4 (dt sweep at
  0.99/1.01·dt_cr with and without ties).
- **Energy:** for a tie satisfied from t=0 the projection removes nothing (the violating
  component is zero) — no spurious dissipation; closure checked with
  `EnergyBalanceRecorder` (T9).

### 2.3 Constraint groups = connected components

MPs sharing nodes must be enforced **jointly**, not sequentially (LS-DYNA's
last-evaluated-wins is the failure mode we refuse to inherit — Vol I p. 923). At
`handle()` time, build the undirected constraint graph (nodes = DOFs, edges = MP rows)
and take connected components; each component becomes one group with one stacked `L`.
A rigid diaphragm (N slaves, one master) is then a *single* group whose dense solve is
3×3 (the master's in-plane DOFs) — small, local, cheap.

**Chains** (a constrained node retained elsewhere) are *composable* in principle
(substitute `C` matrices, as Abaqus does); **v1 refuses them with a named error**
identifying both MPs — already infinitely better than upstream Transformation's silent
DOF-map corruption. Composition is Phase 4.

### 2.4 SP interaction

- Homogeneous SPs (`fix`): excluded from the equation set, PlainHandler-style. A retained
  DOF that is SP-fixed is removed from the group's free set (its `a ≡ 0` enters the
  projection as a known zero — equivalently, delete that column of `L`).
- A **slave DOF that is also SP-fixed** = overconstraint → named error at `handle()`.
- Non-homogeneous SP / `imposedMotion`: **not in v1** (§6). The DYNA-style answer
  (overwrite `a` on those DOFs with the prescribed acceleration before projecting the
  rest) is Phase 4.

### 2.5 Massless DOFs

`LᵀML` singular ⟺ a whole retained direction is massless → **named error** at first
projection build ("massless retained set", echoing the RBE2 massless-DOF-scan lesson).
A massless *slave* DOF is mechanically fine (its raw `a` is irrelevant — the projection
overwrites it), **but** the diagonal SOE will refuse the zero-diagonal solve before we
ever get there. v1: detect at `handle()` (mass scan over slave DOFs), error with the
recipe ("add nodal mass or use the penalty path"). A SOE-cooperative skip is Phase 4
(needs a vanilla touch; not worth it for v1). → Open question O1.

---

## 3. Architecture & seams

### 3.1 Two fork-owned pieces, zero vanilla seam changes

1. **`LadrunoProjectionHandler`** (`SRC/analysis/handler/`, HANDLER_TAG **33001** —
   handler band is 1–6 upstream; fork private band ≥33000, first handler entry).
   `handle()`:
   - creates plain `DOF_Group`s / `FE_Element`s for everything (PlainHandler-style;
     homogeneous SPs excluded from equations; **MPs deliberately NOT enforced at
     assembly** — slave DOFs keep their own equations and their own diagonal mass);
   - builds the constraint graph → groups → per-group `(dof lists, C)` (§2.3) and runs
     the full diagnostic battery (§5.2);
   - at `doneNumberingDOF()` (equation numbers now known) freezes each group's equation
     indices into a `LadrunoConstraintProjector` object.
2. **A small consumer interface, not a concrete downcast**:
   `LadrunoProjectionConsumer` (abstract, one method:
   `setConstraintProjector(LadrunoConstraintProjector*)`), implemented by
   `CentralDifferenceLadruno` in v1. The handler pushes it: the base class's
   `setLinks()` already hands the handler the `Integrator&`
   (`ConstraintHandler.h:56-58`), so the handler `dynamic_cast`s to the *interface* —
   never to the concrete integrator class — and pushes. This makes P4 adoption by
   ExplicitBathe/LNVD a pure integrator-side change (implement the interface, add the
   per-sub-step insertion points); the handler is untouched. The projection is thus
   **explicit-family-portable but deliberately not integrator-agnostic**: implicit
   schemes are excluded on principle (D6 — projecting Newton iterates without also
   projecting residual+tangent breaks constrained equilibrium; doing it consistently
   *is* the Transformation method, which already exists). Coupled-explicit schemes
   with damping on the LHS (e.g. `NewmarkExplicit`'s `(M+γΔtC)`) would change the
   projection weight from `M` to the LHS matrix — out of scope; the interface contract
   assumes an M-only LHS. **If the cast fails → hard, early, named error**
   ("LadrunoProjectionHandler requires an explicit projection-aware integrator; use
   constraints Transformation/Penalty for implicit analyses"). Both classes are
   fork-owned → no `AnalysisModel`/vanilla edit for the seam. (classTags.h +
   commands-registration are the only vanilla touches; standard LEDGER_vanilla rows.)

### 3.2 The projector object

`LadrunoConstraintProjector` (header-mostly, fork-authored, LADRUNO header stamp):
- per group: equation-index ID, `C` matrix, factored `LᵀML` (dense, tiny, LAPACK
  `dgetrf/dgetrs` or hand-rolled ≤6×6);
- `buildMass()`: lazy, on first use after each `domainChanged()` — reads diag(M) from
  the assembled DiagonalSOE *after* the integrator's M-only `formTangent()` (the
  integrator calls `projector->buildMass(theLinSOE)` right after its first factor; mass
  doesn't change between domain changes);
- `project(Vector& a)`: in-place, per group, the §2.1 update;
- `tieForces(const Vector& a_raw, const Vector& a_proj)`: Phase 3 query.

### 3.3 Insertion points in `CentralDifferenceLadruno` (3 lines of behavior)

The integrator's flow (CentralDifferenceLadruno.cpp) gains exactly:
1. **starter** (`newStep()` first-step block, after `theLinSOE->solve()`,
   ~line 432): project `a₀` before `Aprev = X` and before the
   `v_{−1/2} = v₀ − Δt/2·a₀` seed — keeps the back-half-step velocity on the manifold;
2. **main path** (`update(U)`, before the `Vfull` computation ~line 493): copy `U`,
   project, use the projected vector for `Vfull`, `setResponse`, and `Aprev`;
3. **`domainChanged()`**: IC-compliance check (`g(u₀)`, `C v₀ − v_c`) → warn /
   `-projectICs`; invalidate the cached mass.

No change to `formEleTangent` (M-only stays), no change to the SOE, no change to the
leap-frog algebra. ExplicitBathe/LNVD adoption is Phase 4 (same hook, more insertion
points per sub-step).

### 3.4 Command

```
constraints LadrunoProjection <-verbose> <-projectICs> <-icTol $tol>
```

`Ladruno` prefix = collision-proof and greppable, per house convention. `-verbose`
prints the group report (per-group size, retained set, mass totals, `dt_cr` note) in
the spirit of AutoConstraintHandler's report.

---

## 4. What v1 enforces (scope, locked)

| Constraint | v1 | Mechanism |
|---|---|---|
| homogeneous SP (`fix`, `fixX/Y/Z`) | ✅ | equation exclusion (Plain-style) |
| `equalDOF` / `equalDOF_Mixed` | ✅ | projection, mass-weighted average |
| `rigidLink -bar` | ✅ | projection (translation block) |
| `rigidLink -beam`, `rigidDiaphragm` | ✅ | projection with transport columns of `C` (small-rotation `C`, frozen — same kinematics as upstream MP, enforced exactly) |
| `equationConstraint` (EQ_Constraint, upstream 5/2025) | Phase 3 | same `L` machinery, multi-retained group |
| non-homogeneous SP / `imposedMotion` | Phase 4 | prescribed-`a` overwrite, then project the rest |
| MP chains (constrained node retained elsewhere) | refused, named error | composition = Phase 4 |
| time-varying `C`, finite-rotation links | ❌ non-goal | needs per-step `L` rebuild + drift control; element route (RBE2) covers it meanwhile |
| `Pressure_Constraint` | ❌ non-goal | PFEM-only, separate pipeline |
| RBE2/RBE3 eliminable block | Phase 4 | route the elements' rigid part through the projector, retire bipenalty where eliminable |
| OpenSeesSP/MP partition-crossing MPs | ❌ non-goal v1 | groups must be partition-local; refuse + message |

---

## 5. Diagnostics — the shared sub-component

Building §2.3's constraint graph *is* 90% of a `verifyConstraints` facility, which the
review of 2026-06-09 identified as the #1 robustness gap vs LS-DYNA (error-stop on tied
conflicts, Theory p. 543) and Abaqus (automatic overconstraint resolution). So:

### 5.1 Ship the graph as its own small library
`LadrunoConstraintGraph` (used by the handler; reusable later by a standalone
`verifyConstraints` command and by the other handlers if upstreamed).

### 5.2 Named errors/warnings at `handle()` (all tested, T5)
- slave DOF under **two MPs** (upstream Transformation: silent corruption) → error;
- **chain** (slave is retained elsewhere) → error naming both MP tags;
- **SP ∩ slave DOF** → error;
- **redundant cycle** (u1=u2, u2=u3, u3=u1) → error (rank check per group);
- **zero free DOFs** after exclusion (upstream: segfault) → clean error;
- **massless retained set / massless slave DOF** (§2.5) → error with recipe;
- duplicate SPs on one DOF → warning (Auto-handler parity).

---

## 6. Validation battery (Zone-A pytest, locked tolerances)

| # | Test | Pass criterion |
|---|---|---|
| T1 | two-mass `equalDOF` tie, impulse on one mass | combined-mass analytic trajectory; linear momentum error < 1e-12 |
| T2 | same model: CD+Projection vs CD+**Transformation+BandGen** (correct reference) | trajectories match < 1e-10 |
| T3 | rigid-diaphragm slab + 4 columns, base excitation | vs implicit Newmark+Transformation reference (T2a-tier tolerance); diaphragm gap ≡ 0 to machine eps |
| T4 | dt sweep 0.99·dt_cr (stable) / 1.01·dt_cr (diverges) **with and without ties** | tie does not change the stability boundary (interlacing claim §2.2) |
| T5 | conflict battery (§5.2, 7 cases) | each refused with its named error; no segfault, no silent pass |
| T6 | **Phase-0 falsification**: rigidDiaphragm + CD + Transformation + `system Diagonal` vs BandGen | expected: Diagonal drops `TᵀMT` off-diagonals → wrong response. Outcome → LEDGER_quirks row either way |
| T7 | `rigidLink -beam` with rotary inertia, free vibration | frequencies vs assembled-reference < 1e-8 rel |
| T8 | tie-force recovery `M(a_raw−a_proj)` | matches high-α penalty reference force within 0.1% |
| T9 | energy closure on T3 via `EnergyBalanceRecorder` | drift < 1e-3 of peak KE over the record |

---

## 7. Phases

- **P0 — falsify & baseline (no SRC change):** run T6; baseline Transformation+BandGen
  explicit cost vs Diagonal on a ~50k-DOF model (quantifies the win); confirm
  DiagonalSOE zero-diagonal behavior for §2.5.
- **P1 — core:** handler (Plain-style assembly + graph + diagnostics) + projector +
  CDL hook; `equalDOF` only. T1, T2, T4, T5 green.
- **P2 — general C:** rigidLink/rigidDiaphragm transport columns, SP-fixed-retained
  column deletion, group rank checks. T3, T6, T7 green.
- **P3 — queries & EQ:** tie-force recorder query; `EQ_Constraint` groups; `-verbose`
  report polish. T8, T9 green.
- **P4 — deferred:** prescribed-motion overwrite; MP-chain composition; ExplicitBathe/
  LNVD adoption; RBE2/RBE3 eliminable-block routing (retire bipenalty where possible);
  SOE-cooperative massless-slave skip.

Each phase = one PR, **fresh branch off `ladruno`** (stranded-commit lesson: never pile
onto a merged PR branch).

---

## 8. Bookkeeping obligations (same-PR, per house rules)

- `LEDGER_implementations.md`: row for `LadrunoProjectionHandler` (HANDLER 33001,
  + projector/graph files) at P1; status flips per phase.
- `LEDGER_vanilla_files.md`: `classTags.h` (+1 define), handler registration in the
  commands file (`// Ladruno` comment at the edit).
- `LEDGER_quirks.md`: T6 outcome (Transformation+Diagonal mass-coupling drop) — this is
  a quirk worth recording even if we never ship the handler.
- `Ladruno_scripts/stamp_headers.py`: add the new files to GLOBS, re-run.
- `banner_features.txt` + `patch_banner.py`: one line at P2 (first user-visible
  shipped state).
- Zone-A manifest row + classTag gate.

---

## 9. Open questions

- **O1 (massless slaves, §2.5):** is handle()-time refusal acceptable for v1, or do we
  need the SOE-cooperative skip immediately? Proposed: refuse in v1; revisit if the
  embedded-node Phase-4 routing needs massless interface nodes.
- **O2 (`-projectICs` default):** warn-only vs auto-project when `g(u₀) ≠ 0`. Proposed:
  warn + abort by default (silent state mutation is the house anti-pattern), flag to
  opt into projection.
- **O3 (handler name):** `LadrunoProjection` vs plain `Projection`. Proposed: Ladruno
  prefix (collision-proof); revisit only if upstreaming.

---

## 10. Decision record

| # | Decision | Choice | Why |
|---|---|---|---|
| D1 | enforcement mechanism | **acceleration projection**, not Δt-bounded penalty | bounded penalty (DYNA Eq. 25.30) is for compliant joints; ties need exactness AND Δt-neutrality — only projection gives both |
| D2 | vs Transformation in explicit | new handler | `TᵀMT` off-diagonals are incompatible with the Diagonal-SOE recipe (silently dropped, T6); plus upstream's silent one-MP limit / 0-free-DOF segfault |
| D3 | where the projection lives | integrator post-solve hook behind a `LadrunoProjectionConsumer` interface; handler builds the operator | matches the math (projection acts on `a`); handler already receives the Integrator via `setLinks` → zero vanilla seam; interface (not concrete downcast) keeps the handler untouched when ExplicitBathe/LNVD adopt it (P4) |
| D4 | simultaneous vs sequential group enforcement | connected-component groups, joint solve | refuses DYNA's last-evaluated-wins fragility; groups are tiny (≤ retained-DOF count) |
| D5 | chains | refuse with named error (v1) | diagnostic-over-silent; Abaqus-style composition deferred (P4) |
| D6 | scope | explicit-only, CDL-first | implicit already well served; CDL is fork-owned so the seam is ours to cut |
