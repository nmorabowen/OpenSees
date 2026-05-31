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
- **Ladruno brick element(s)** — our own higher-order hexahedral element(s), the solid-side sibling of BezierTri6 (planned). Will implement the [[ladruno_element_contract]] for zero-edit Ladruno recording, with non-negative lumped mass for explicit dynamics and B-bar/assumed-strain against volumetric locking. Scope/order TBD — plan doc to be written. (draft, no plan file yet)
- [[ladruno_apegmsh_contract]] — **apeGmsh feature reference**: the fork-only features apeGmsh emits/reads, with the canonical command and apeGmsh touch-points for each. The companion to [[LEDGER_implementations]] (which is authoritative for tags + PRs). (reference)

## Companion folder

- [[../Ladruno_internal/README|Ladruno_internal]] — internal docs about the existing build and patches.
