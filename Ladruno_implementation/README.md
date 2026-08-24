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

> [!warning] Known stale content, audited 2026-08-23 (at `b17e8bd82`)
> Two classes of rot are flagged in place rather than silently fixed, so nobody
> re-discovers them:
>
> - **12 shipped ADRs still carry a pre-implementation `status:`** — 09, 19, 20
>   (stabilized arc-length), 30 (explicit constraint projection), 31, 32, 35, 36, 45,
>   52, 75, 78 (parallel contact). Each now opens with a `> [!warning]` banner.
>   **Rule: trust [[LEDGER_implementations]] for *does it work*, the ADR for *why*.**
>   A CI check cross-validating frontmatter `status:` against the ledger's Status cell
>   would stop this recurring; none exists yet.
> - **12 unresolved `[[NN_topic]]` links in [[Ladruno_explicit_roadmap]]** are forward
>   references to section docs that were never written, under a local numbering that is
>   *not* the ADR numbering. Explained in a banner at the top of that file.
>
> Everything else that this audit found was fixed outright: the `ELE_TAG_LadrunoUP`
> 33017 "RESERVED" cells, the superseded contact ELE-33016 reservation (in the ledger
> and in ADR-61/62), a `Fix` bullet stranded in the wrong [[LEDGER_quirks]] entry by a
> squash merge, and renumbering-damaged links in ADR 42, 40 and 25.

## Element selection & usage

- [[ladruno_continuum_elements_guide]] — **Continuum Elements — Modeling &
  FE-Selection Guide**: the single decision desk for picking between
  `BezierTri6` (2D), `BezierTet10` (3D tet) and `LadrunoBrick` (3D hex) —
  selection axes, a decision procedure, per-element intended-use profiles, and
  cross-cutting modeling guidance. Links down to the per-element references.
  (reference)
- [[LadrunoBrick_reference]] — the brick's living theory/implementation/usage
  reference (the deep doc the selection guide points to). (reference)

## What a plan contains

Forward-looking planning docs for new functionality we want to add to this OpenSees fork. Each plan lives in its own file and walks through:

- **What** — feature description, scope, non-goals
- **Why** — motivation, what user problem it solves
- **Where** — files in `OpenSees/SRC/` that need to change, similar existing implementations to reference
- **How** — design sketch, API, integration points, testing strategy
- **Risks** — what could go wrong, dependencies, open questions

## Conventions

- One file per implementation, numbered loosely by priority.
- **A number is a permanent identifier once the doc is cited from `SRC/` or `tests/`.** Never
  reuse a number for a new subject and never renumber an existing file: there are ~1,100 numbered
  wikilinks in this folder and ~2,700 bare `ADR NN` mentions across ~300 files in `SRC/` (plus
  ~1,000 in `tests/`), and the bare form is exactly the form that cannot be mechanically rewritten.
  Allocate the next free integer instead. (An earlier renumber already left
  `42_ladruno_buckling_adr.md` pointing at a *different, plausible-looking* ADR — fixed, but that
  is the failure mode.) Superseding a decision is done with `supersedes:` / `amends:` frontmatter,
  not by taking over its number.
- **Bare `ADR NN` is only legal where `NN` is unique.** These numbers currently carry more than
  one doc — cite them with a qualifier or a full wikilink:

  | # | resolves to |
  |---|---|
  | **78** | parallel MPI contact · quadratic mortar tie · `LadrunoUP -geom corot` (+ an apeGmsh companion) |
  | **20** | stabilized arc-length · brick EAS stabilization · embedded reinforcement |
  | **30** | explicit constraint projection · parallel numberer |
  | **31** | LadrunoConcrete3D · robust solve driver |
  | **19** | brick EAS Simo–Rifai · RC shell |
  | 06 · 14 | see the filenames — three files each |
  | 04 · 08 · 09 · 10 · 21 · 22 · 52 · 80 | an ADR plus a companion (`_handoff`, `_findings`, `_validation_gates`, …) whose suffix disambiguates it |

  Pairing a handoff or gates doc with its ADR under one number is deliberate and fine. Two
  unrelated ADRs sharing a number is not — that is what the qualifier rule is for.
- Use Obsidian-flavored Markdown: frontmatter, wikilinks, fenced code blocks. Plain Markdown also renders elsewhere.
- Cross-link to [[../Ladruno_internal/01_compilation_journal|the compilation journal]] when a plan touches the build (e.g. needs a new dependency, new compile flag, new CMake target).
- A plan starts as a sketch — incomplete sections are fine. Mark unresolved questions with `> [!question]` callouts so they stand out in Obsidian.
- When a plan is implemented, add a final `## Implementation log` section at the bottom with commits / dates / surprises, **and flip its `status:` frontmatter to match its row in [[LEDGER_implementations]]**. The file stays here — a shipped ADR is the permanent record of *why*, and `SRC/` cites it by name.
  > **Note.** This rule previously said to move the file to `../Ladruno_internal/` as
  > `implemented_<name>.md`. That was never once done (0 of 74 ADRs; `Ladruno_internal/` holds no
  > `implemented_*.md` at all), and moving files would break the ~2,700 `ADR NN` citations in
  > `SRC/`. The rule now describes what the fork actually does.
- **Keep frontmatter `status:` current.** It drifts badly: at `b17e8bd82` at least 21 ADRs whose
  ledger row says `shipped` still had frontmatter reading `draft` / `proposed` / `NO code`. Trust
  the ledger for *does it work*, the ADR for *why*. A CI check cross-validating the two would stop
  the drift.

## Plan template

When starting a new plan, copy [[_template]] and rename.

## Currently-active plans

> **This is a curated selection, not the folder index.** The folder holds ~119 numbered docs
> (~74 ADRs plus execution artifacts); this list is the handful worth reading next. For the full
> picture use `ls *_adr.md`, and use [[LEDGER_implementations]] — which carries class tags, PR
> numbers and status — as the authoritative catalogue of what exists.

- [[03_ladruno_recorder]] — **Ladruno**: modular recorder fork (sibling of the frozen MPCORecorder), apeGmsh-native `.ladruno` schema, global + envelope results. (ADR, draft)
  - [[ladruno_schema_v1]] — the on-disk HDF5 schema spec (self-describing BASIS/QUADRATURE for Bézier + Belytschko). (draft)
  - [[ladruno_element_contract]] — the element-side `setResponse` contract elements implement to be recorded. (draft)
- [[04_bezier_elements]] — **BezierTri6**: 6-node quadratic Bézier triangle (Kadapa 2018) — non-negative lumped mass for explicit dynamics + consistent B-bar. v1 = straight-sided Tri6, **merged** (`ELE_TAG 33000`, PR #6); implements the [[ladruno_element_contract]]. (ADR)
  - [[bezier_apegmsh_integration]] — how apeGmsh meshes drive BezierTri6 (direct-drive today; typed-primitive deferred). Regression test: `Ladruno_scripts/bezier_tests/test_bezier_tri6.py`.
- [[06_bezier_tet10]] — **BezierTet10**: 10-node quadratic Bézier tetrahedron (Kadapa 2018 §5) — the 3D sibling of [[04_bezier_elements]] (deferred there under D10). Non-negative lumped mass `ρVe/10` (Eq. 57) for explicit 3D dynamics + 3D B-bar for near-incompressibility. v1 = straight-sided. (ADR, draft)
- [[72_ladruno_second_order_brick_adr]] — **ADR 72 — LadrunoBrick20**: 20-node serendipity quadratic hex (`ELE_TAG 33018` planned; 33017 went to ADR-71 LadrunoUP), `-formulation {std|uri}` (27-pt full with a reduce-to-`Twenty_Node_Brick` anchor / 2×2×2 Barlow-point reduced, no hourglass control needed), HRZ-only `-lumped` (row-sum has −M/8 corner masses), implicit-first with explicit permitted-but-discouraged; contact faces + `-geom finite` + H27 sibling deferred. Fills the long-standing "higher-order hex" placeholder row. (ADR, draft)
- [[ladruno_apegmsh_contract]] — **apeGmsh feature reference**: the fork-only features apeGmsh emits/reads, with the canonical command and apeGmsh touch-points for each. The companion to [[LEDGER_implementations]] (which is authoritative for tags + PRs). (reference)
- [[58_ladruno_rigid_body_adr]] — **ADR 58 — RigidBody DomainComponent + SO(3) integrator**: promotes explicit-roadmap §5.5 to a numbered decision-capture. A new top-level `DomainComponent` kind (parallel to `Element`/`MP_Constraint`) owning 6 DOFs + inertia tensor + a Lie-group rotational integrator, with LS-DYNA §25 CoM mass-condensation so the explicit/diagonal-mass path stays clean. Gates the joint family (roadmap §5.6) and is the missing "free rigid body" home for AEM contact debris ([[51_ladruno_element_removal_adr|ADR 51]]). (ADR / scoping, no code, draft)
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
