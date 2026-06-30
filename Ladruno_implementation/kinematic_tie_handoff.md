---
title: "LadrunoTie (ADR-62) — handoff: concept VALIDATED, build the auto-generator next"
project: Ladruno
type: handoff
status: P0 (numpy) + P1 (real-solver) GREEN; the LadrunoTie auto-generator (C++ + build) is the remaining work
related:
  - "[[62_ladruno_kinematic_mesh_tie_adr]]"      # the spec
  - "[[30_ladruno_explicit_constraint_projection_adr]]"   # the enforcement (SHIPPED handler)
  - "[[projection_handler_handoff]]"             # the projection handler's API + limits
  - "[[39_ladruno_contact_domain_adr]]"          # bucket-sort + closest-point projection to reuse
  - "[[61_ladruno_contact_bipenalty_adr]]"       # the SHELVED penalty route this replaces
tags: [adr, constraints, explicit-dynamics, mesh-tie, projection, kinematic, handoff]
---

# LadrunoTie — handoff (next session: build the generator)

## TL;DR

A non-conforming **explicit mesh-tie** is enforced KINEMATICALLY (not by penalty) by emitting
`u_s = Σ N_i u_{m,i}` as constraints and letting the **shipped** `LadrunoProjectionHandler`
(ADR-30) enforce them — exact, `dt_cr`-neutral, momentum-clean, no fictitious mass. ADR-61's
penalty routes (SOFT penetration, bipenalty ~100× mass) are replaced.

**The concept is fully validated (this session, two oracles green):**
- `Ladruno_implementation/kinematic_tie_validation/proto_p0_kinematic_tie.py` — the projection
  math (numpy): exact, idempotent, energy-clean, + the penalty-`dt_cr`-collapse contrast.
- `Ladruno_implementation/kinematic_tie_validation/proto_p1_kinematic_tie_opensees.py` — the
  REAL solver (shipped binary, no build): a weighted multi-master non-conforming tie
  `u_4 = 0.7 u_2 + 0.3 u_3`, tie error **6.7e-18** / 400 steps, `dt_cr` tied==untied (ratio
  **1.000000**), load-carrying.

**It already works TODAY with no new code** (this is the key point): per slave node, hand-emit
```python
ops.equationConstraint(s, dof, 1.0, m1, dof, -N1, m2, dof, -N2, m3, dof, -N3)
ops.constraints("LadrunoProjection"); ops.system("Diagonal")
ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear")
```
**The remaining work = a setup-time C++ generator** that automates the GEOMETRY (pair slave
nodes → master facet, compute `N_i(ξ̄)`, emit the constraints). Pure ergonomics; correctness
is settled.

---

## The task for next session (P1 productization)

Build **`LadrunoTie`** — a standalone constraint generator (decided: API option **(a)**,
ADR-62 OQ-1). Scope = **solid–solid, node-to-segment collocation** (ADR-62 OQ-4 = P1 only;
defer shells + integral-mortar). At model-build it:

1. **Pairs** each slave-surface node to the nearest master facet (closest-point projection).
2. **Evaluates** the master facet shape functions `N_i(ξ̄)` at the projection point.
3. **Emits**, per slave node per translational DOF, an `EQ_Constraint`
   `1·u_s + Σ_i (−N_i)·u_{m,i} = 0` (exactly what `equationConstraint` builds).

User then runs `constraints LadrunoProjection` (shipped). Done.

### Where to reuse code (do NOT re-derive)
- **Closest-point projection + facet shape weights** — the contact engine already does this:
  `SRC/analysis/handler/LadrunoContactFE.cpp::segmentActive(...)` fills the oriented normal,
  the gap, AND `N[nps]` (the master shape weights at the projection) + the B-operator. That is
  precisely the per-slave-node geometry the generator needs. The oracle-validated projection
  is `LadrunoContactProjection` (ADR-39). Lift the projection+`N[]` computation; you do NOT
  need the penalty force/tangent.
- **Broad phase** — `SRC/domain/contact/LadrunoContactBucketSort.*` (ADR-39) to find the
  candidate master facet per slave node (reference config; a tie is static so one pass).
- **Constraint emission** — `SRC/domain/constraints/EQ_Constraint.{h,cpp}`. Construct one per
  slave node (coef vector = `[1, −N_1, …, −N_nps]` over DOFs `[s, m_1, …, m_nps]`) and
  `theDomain->addEQ_Constraint(...)` (mirror what the `equationConstraint` command does —
  grep `Py_ops_equationConstraint` in `SRC/interpreter/PythonWrapper.cpp:941` and the Tcl
  twin in `TclWrapper.cpp` for the exact construction + domain-add call).
- **Command registration** — add `LadrunoTie` to `SRC/interpreter/PythonWrapper.cpp` +
  `TclWrapper.cpp` (mirror an existing geometry/command pair). Suggested signature:
  `LadrunoTie -slaveNodes {…} -masterFacets {…}` (or `-masterSurface` by element set);
  internally resolve facets → nodes.

### `equationConstraint` signature (verified this session)
`equationConstraint(cNode, cDOF, cCoef, rNode, rDOF, rCoef, …)` — variadic; enforces
`cCoef·u_c + Σ rCoef·u_r = 0` ⇒ `u_c = −(Σ rCoef·u_r)/cCoef`. For a tie use `cCoef=1`,
`rCoef_i = −N_i`. Confirmed enforced to `<1e-12` by the projection handler (it accepts the
multi-master EQ row — ADR-30 P3, `tests/test_adr30_projection_p3.py::test_EQ_*`).

---

## GOTCHAS hit this session (save the next session the debugging)

1. **`dt_cr` reads ELEMENT mass, not nodal `mass`.** `criticalTimeStep()` returned −1 ("no
   element produced a finite estimate") with `ops.mass(...)` nodal mass — the CriticalTimeStep
   kernel uses `element->getMass()`. Give elements density (`-rho`) so they carry mass (the
   assembled diagonal then also feeds the projection handler's "tied DOFs need lumped mass"
   requirement). This is the same nodal-vs-element-mass split that motivates the SOFT scheme's
   `ladrunoBuildNodalMass`.
2. **`criticalTimeStep()` needs `-cfl` AND a first step.** Pass `integrator
   CentralDifferenceLadruno -cfl`, then `analyze(1, dt0)` once with a SAFE `dt0` (≪ dt_cr) to
   trigger `domainChanged` (which computes dt_cr), THEN query. Querying before any step → −1.
3. **2D Truss defeats the dt_cr eigensolve** (transverse zero-stiffness direction → no finite
   estimate). Use a pure **1D (ndm=1,ndf=1)** axial chain for the bar oracle, or a real
   solid/shell for the patch test.
4. **Projection handler refuses MP-chains + double-constrained DOFs** (`projection_handler_
   handoff.md`). The generator MUST guarantee node-disjoint slave/master surfaces and ONE
   facet per slave node (collocation) → no chains. Detect-and-refuse otherwise (hand off to
   SOFT). This is BLOCKER-1 in ADR-62.
5. **Tied DOFs need lumped mass** (handler keeps slave DOFs in the equation set). Real solid
   nodes have it; check + named-refuse a massless tied DOF (BLOCKER-2).
6. **IC on the manifold** (BLOCKER-3, ADR-62 OQ-3 STILL OPEN): a non-conforming as-built
   interface starts with the slave off the master facet → the handler aborts unless ICs satisfy
   `u_s = Σ N_i u_{m,i}`. Decide: snap the slave onto the facet at emission (`-projectICs`
   semantics) vs require conforming-at-interface geometry. **Needs sign-off before coding the
   IC path.**

---

## Build + run recipe (this machine)

- **Run the validated oracles** (no build): Python 3.12 at
  `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe` (has numpy). The real-
  solver oracle bootstraps the shipped binary itself (`os.add_dll_directory` +
  `sys.path.insert` on `…\OpenSees\dist\bin`; `LADRUNO_OPENSEES_QUIET=1`;
  `sys.stdout.reconfigure(encoding="utf-8")` for the Windows console). No oneAPI setvars
  needed to RUN.
- **Build** (for the new command): `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`
  — invoke via the **PowerShell tool** (`cmd /c …`), NOT the Bash tool (build-gotcha,
  [[project_bezier_charlen]]). If configure fails with ZLIB-not-found, it's the system
  CMake-4.3 shadow — use the conan cmake by full path ([[project_build_env_cmake43_conan_zlib]]).
  ~30 min. The built `opensees.pyd` lands in `dist/bin`.
- **Test**: add `tests/test_ladrunoTie_*.py` (Zone-A, `pytest.mark.zone_a`, `from _testbed
  import ops`). Mirror `tests/test_adr30_projection_p3.py`.

---

## P1 validation to add (after the generator builds)

- **Solid–solid patch test**: two non-conforming solid blocks tied across a flat interface;
  a uniform stress field must transmit EXACTLY across the tie (constant-stress patch). Reuse
  `_testbed/fem_checks.check_constant_stress`.
- **Split-bar equivalence**: a bar split by a non-conforming tie matches the monolithic mesh
  (disp + dt_cr).
- **Penetration vs SOFT**: the kinematic tie gap is **0**; the SOFT tie leaves
  `δ/h = ε_iface·CFL²/(4·α_m·SOFSCL)` (the ADR-61 P0 closed form). Tabulate the contrast.
- Momentum + energy are already covered by the shipped projection-handler tests (T8/T9).

---

## Bookkeeping (when the generator ships)

- `LEDGER_implementations.md` — row: *LadrunoTie kinematic mesh-tie generator* (emits
  `EQ_Constraint`s for `LadrunoProjectionHandler`), **no new class tag**, status, PR.
- `LEDGER_vanilla_files.md` — none expected (public Domain `addEQ_Constraint` API).
- `LEDGER_quirks.md` — the gotchas above (dt_cr-element-mass; the projection-handler
  requirements the generator must refuse-and-hand-off on).
- `classTags.h` — no change (free contact ELE slot is 33016 if a tagged object is ever needed).
- Banner — a `banner_features.txt` line for LadrunoTie.

## Decided / open
- **DECIDED**: API = standalone `LadrunoTie` (OQ-1a); collocation node-to-segment (OQ-2);
  `-dtcr` n/a (no penalty); scope = solid–solid P1 (OQ-4).
- **OPEN (needs sign-off)**: OQ-3 IC handling (snap vs require-conforming) — gotcha #6.

## Status of the PRs
- **PR #449** (`guppi/adr62-kinematic-mesh-tie` → `ladruno`): ADR-62 + ADR-61-finalized +
  the P0/P1 oracles + this handoff. Merge it to land the record. (PR #445 already merged the
  original ADR-61.)
