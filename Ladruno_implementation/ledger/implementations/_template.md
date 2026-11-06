---
title: Ledger — our own implementations
project: Ladruno
tags:
  - ledger
  - features
  - implementation
---

# Ledger — Ladruno implementations (new code we authored)

Every brand-new file / feature the Ladruno fork adds on top of vanilla
OpenSees: what it is, its class tag, where it lives, and the PR(s) that
shipped it. This is the counterpart to [[LEDGER_vanilla_files]] (which tracks
edits to *pre-existing* upstream files).

## Conventions

- **One row per feature**, not per file. List the owning files in the *Files* cell.
- *Status*: `shipped` (merged to `ladruno`), `draft` (plan/WIP), or `frozen`.
- **The banner feature list mirrors this ledger.** Every `shipped` feature
  should have a matching line in `Ladruno_scripts/banner_features.txt`. After
  editing that file run `python Ladruno_scripts/patch_banner.py` and rebuild.
- Class tags live in `SRC/classTags.h`; keep them recorded here so we never
  collide on a tag.
- **Reserved tags:** a tag pre-allocated by an ADR but not yet implemented is
  marked `— RESERVED, not yet built`; it is recorded here to prevent collisions
  but does **not** appear in `SRC/classTags.h` until the implementation merges.
  Class-tag bands are per-registry (Element / nDMaterial / uniaxial / Integrator /
  Recorder each have their own 33000-space), so the same number in two registries
  is not a collision.
- When a forward-looking plan in this folder ships, move it to
  `Ladruno_internal/implemented_<name>.md` and add/flip its row here.

## Ledger

| Feature | Kind | Class tag | Files | Status | PR(s) |
|---|---|---|---|---|---|
<!-- ledger:rows table -->

## Documentation / ADR PRs + shipped-feature build history

The first bullets shaped the design but shipped docs only — recorded for
traceability. The later bullets (DispBeamColumn / CohesiveHinge /
CohesiveHingeBiaxial / RigidBody) are the **detailed build logs** of features
that DID ship source; each now also has a summary row in the feature table
above, so the table stays the authoritative shipped-feature index.

<!-- ledger:history -->
