---
title: Ladruno implementation plans
project: Ladruno
tags:
  - index
  - planning
  - implementation
---

# Ladruno implementation plans

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

- [[03_mpco_ladruno]] — **MPCO_Ladruno**: modular recorder fork (sibling of the frozen MPCORecorder), apeGmsh-native `.ladruno` schema, global + envelope results. (ADR, draft)
  - [[mpco_ladruno_schema_v1]] — the on-disk HDF5 schema spec (self-describing BASIS/QUADRATURE for Bézier + Belytschko). (draft)
  - [[mpco_ladruno_element_contract]] — the element-side `setResponse` contract elements implement to be recorded. (draft)
- [[04_bezier_elements]] — **BezierTri6**: 6-node quadratic Bézier triangle (Kadapa 2018) — non-negative lumped mass for explicit dynamics + consistent B-bar. v1 = straight-sided Tri6, **merged** (`ELE_TAG 272`, PR #6); implements the [[mpco_ladruno_element_contract]]. (ADR)
  - [[bezier_apegmsh_integration]] — how apeGmsh meshes drive BezierTri6 (direct-drive today; typed-primitive deferred). Regression test: `Ladruno_scripts/bezier_tests/test_bezier_tri6.py`.

## Companion folder

- [[../Ladruno_internal/README|Ladruno_internal]] — internal docs about the existing build and patches.
