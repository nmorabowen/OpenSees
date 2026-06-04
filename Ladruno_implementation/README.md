---
title: Ladruno implementation plans
project: Ladruno
tags:
  - index
  - planning
  - implementation
---

# Ladruno implementation plans

This folder holds two kinds of doc:

1. **Ledgers** (`LEDGER_*.md`) — the always-current build-control record of what
   the fork has changed. Update these as part of every feature PR.
2. **Plans** — forward-looking design docs for functionality we want to add.

## Build-control ledgers

On a fast development track we keep three running ledgers so we never lose track
of what diverged from upstream:

- [[LEDGER_vanilla_files]] — every **upstream file we touched**, why, and the PR.
  The rebase-onto-upstream checklist.
- [[LEDGER_implementations]] — every **new feature/file we authored**, its class
  tag, and the PR. Mirrors the splash-banner feature list.
- [[LEDGER_quirks]] — **OpenSees gotchas** we learned the hard way.

> [!important] Banner ↔ ledger sync
> The splash banner prints the active-feature list. Its source of truth is
> `Ladruno_scripts/banner_features.txt`; every `shipped` row in
> [[LEDGER_implementations]] should have a matching line there. After editing,
> run `python Ladruno_scripts/patch_banner.py` and rebuild — the script
> regenerates the `FEATURES-START/END` blocks in `tclMain.cpp` (Tcl) and
> `PythonModule.cpp` (openseespy/mp).

## Element selection & usage

- [[ladruno_continuum_elements_guide]] — **Continuum Elements — Modeling &
  FE-Selection Guide**: the single decision desk for picking between
  `BezierTri6` (2D), `BezierTet10` (3D tet) and `LadrunoBrick` (3D hex) —
  selection axes, a decision procedure, per-element intended-use profiles, and
  cross-cutting modeling guidance. Links down to the per-element references.
  (reference)
- [[LadrunoBrick_reference]] — the brick's living theory/implementation/usage
  reference (the deep doc the selection guide points to). (reference)

## Plans

Forward-looking planning docs for new functionality we want to add to this OpenSees fork. Each plan lives in its own file and walks through:

- **What** — feature description, scope, non-goals
- **Why** — motivation, what user problem it solves
- **Where** — files in `OpenSees/SRC/` that need to change, similar existing implementations to reference
- **How** — design sketch, API, integration points, testing strategy
- **Risks** — what could go wrong, dependencies, open questions

## Conventions

- One file per implementation, numbered loosely by priority. Rename files freely as priorities shift — don't worry about preserving filenames as stable links until a plan is actually being executed.
- Use Obsidian-flavored Markdown: frontmatter, wikilinks, fenced code blocks. Plain Markdown also renders elsewhere.
- Cross-link to [[../Ladruno_internal/01_compilation_journal|the compilation journal]] when a plan touches the build (e.g. needs a new dependency, new compile flag, new CMake target).
- A plan starts as a sketch — incomplete sections are fine. Mark unresolved questions with `> [!question]` callouts so they stand out in Obsidian.
- When a plan is implemented, add a final `## Implementation log` section at the bottom with commits / dates / surprises, then move the file to [[../Ladruno_internal/]] under a `implemented_<name>.md` name. The plan stops being forward-looking and becomes history.

## Plan template

When starting a new plan, copy [[_template]] and rename.

## Plans

- [[03_ladruno_recorder]] — **Ladruno**: modular recorder fork (sibling of the frozen MPCORecorder), apeGmsh-native `.ladruno` schema, global + envelope results. (ADR, draft)
  - [[ladruno_schema_v1]] — the on-disk HDF5 schema spec (self-describing BASIS/QUADRATURE for Bézier + Belytschko). (draft)
  - [[ladruno_element_contract]] — the element-side `setResponse` contract elements implement to be recorded. (draft)
- [[04_bezier_elements]] — **BezierTri6**: 6-node quadratic Bézier triangle (Kadapa 2018) — non-negative lumped mass for explicit dynamics + consistent B-bar. v1 = straight-sided Tri6, **merged** (`ELE_TAG 272`, PR #6); implements the [[ladruno_element_contract]]. (ADR)
  - [[bezier_apegmsh_integration]] — how apeGmsh meshes drive BezierTri6 (direct-drive today; typed-primitive deferred). Regression test: `Ladruno_scripts/bezier_tests/test_bezier_tri6.py`.
- [[06_bezier_tet10]] — **BezierTet10**: 10-node quadratic Bézier tetrahedron (Kadapa 2018 §5) — the 3D sibling of [[04_bezier_elements]] (deferred there under D10). Non-negative lumped mass `ρVe/10` (Eq. 57) for explicit 3D dynamics + 3D B-bar for near-incompressibility. v1 = straight-sided. (ADR, draft)
- **Ladruno brick element(s)** — our own higher-order hexahedral element(s), the solid-side sibling of BezierTri6 (planned). Will implement the [[ladruno_element_contract]] for zero-edit Ladruno recording, with non-negative lumped mass for explicit dynamics and B-bar/assumed-strain against volumetric locking. Scope/order TBD — plan doc to be written. (draft, no plan file yet)
- [[ladruno_apegmsh_contract]] — **apeGmsh feature reference**: the fork-only features apeGmsh emits/reads, with the canonical command and apeGmsh touch-points for each. The companion to [[LEDGER_implementations]] (which is authoritative for tags + PRs). (reference)
- [[19_ladruno_rc_shell_adr]] — **Ladruno RC shell stack**: a header-only `LadrunoRCKernel.h` (cloning [[10_ladruno_j2_plasticity|LadrunoJ2Kernel]]'s "one core, many views" pattern) that adds MCFT compression softening + degrading aggregate-interlock shear + tension stiffening to `ASDConcrete3D`'s plastic-damage spine, delivered as an order-5 `PlateFiber` `nDMaterial` view that drops into the **unmodified** `ASDShellQ4` + `LayeredShellFiberSection` seam. 5-phase path; Phase 1 closes the squat-wall in-plane-shear gap with zero element/section edit. Designed via a 6-dimension design-panel + adversarial workflow (β-on-strength-axis is the blocking Phase-1 gate). (ADR, draft)

## Reference guides (shipped features)

User-facing, living reference docs for features already on `ladruno` (theory →
architecture → OpenSees implementation → usage), distinct from the forward-looking
plans above:

- [[Ladruno_materials_guide]] — **the material catalog**: every fork-authored
  constitutive material (the J2 plasticity core, the finite-strain & staged
  wrappers, the steel/rebar overlays), organized by family with theory, OpenSees
  command, and use case. The single entry point for materials; links the per-material
  guides below.
- [[finite_strain_trifecta_guide]] — **the large-deformation stack**: how the
  element geometry layer (`-geom corot|finite`), the Hencky material wrapper
  (`nDMaterial LogStrain`), and the constitutive law (`LadrunoJ2`) compose into
  finite-strain elastoplasticity. The single entry point; links the three per-leg
  guides below.
- [[LadrunoBrick_reference]] — the unified hex element (formulations + geometry seams).
- [[LadrunoJ2_guide]] — combined-hardening von Mises (J2) `nDMaterial`.
- [[LadrunoUniaxialJ2_guide]] — the uniaxial J2 twin (fibers/truss/zeroLength).
- [[LadrunoJ2Finite_guide]] — finite-strain-native combined J2 (co-rotating backstress).
- [[LogStrain_guide]] — the Hencky log-strain finite-strain material adaptor.
- [[LadrunoStaged_guide]] — the `Staged*` family (`InitDefGrad` finite / `StagedStrain` small).
- [[LadrunoLemaitreDamage_guide]] — the Lemaitre ductile-damage mode on the J2 family.
- [[LadrunoRebarBuckling_guide]] — reinforcing-bar buckling overlay `uniaxialMaterial`.
- [[LadrunoBondSlip_guide]] — 1D bond-slip τ–s `uniaxialMaterial` for embedded rebar.
- [[solid_transformation_wrapper]] — the solid geometry-method layer (linear/corot/finite).
- [[09_finite_strain_material_wrapper]] — the log-strain (`LogStrain`) adaptor.
- [[18_finite_strain_validation_report]] — the finite-strain V&V execution record.

## Companion folder

- [[../Ladruno_internal/README|Ladruno_internal]] — internal docs about the existing build and patches.
