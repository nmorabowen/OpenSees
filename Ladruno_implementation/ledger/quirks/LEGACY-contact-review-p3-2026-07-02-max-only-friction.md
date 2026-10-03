---
wp: LEGACY
title: "Contact-review P3 (2026-07-02): τmax-only friction was an UNBOUNDED bond; validation choke points; mortar isomap scale fix"
date: 2026-07-02
legacy_seq: 128
---
### Contact-review P3 (2026-07-02): τmax-only friction was an UNBOUNDED bond; validation choke points; mortar isomap scale fix
- **τmax-only (HIGH-2):** `-tauMax` without `-mu`/`-cohesion` is the unified cone `min(μN+c, τmax) =
  min(0, τmax) = 0` — a ZERO cone radius (free slip). The kernel's `cap≤0` branch returned the RAW
  ELASTIC traction ("byte-identity is the guard's job") with a consistent stick tangent, so the config
  silently became an UNBOUNDED elastic bond — the FE guards treat `tauMax>0` as a friction REQUEST and
  the handler auto-defaults `kt`, so it fell exactly between the tested branches (Tresca is tested as
  μ=0 WITH cohesion). FIX at two layers: the kernel `cap≤0` branch now FREE-SLIPS (zero traction, slip
  absorbs the motion, zero tangent block — `frictionReturnMap` + `frictionTangentBlock` +
  the numpy mirror in `proto_a1_friction.py`), and the command surface REFUSES `-tauMax`/`-edgeTauMax`
  without a `-mu`/`-cohesion` (the sanctioned shear-capped BOND is `-cohesion`, optionally + `-tauMax`);
  on the NTS lane `-tauMax` is refused outright as mortar-only (it was silently INERT there — the
  positional-μ NTS cone has no shear cap; adversarial-gate MINOR-1, fail-loud like `-geomtan`).
  cap>0 branches byte-unchanged (166k-case gate fuzz bitwise-identical + oracle byte-guards pass on
  BOTH pre/post kernels). Oracle: `contact_prototypes/proto_friction_validation.cpp` (12 checks,
  6 FAIL pre-fix). GATE-SURFACED FOLLOW-UP (pre-existing, not fixed here): mortar AUTO-orientation
  silently degenerates for COINCIDENT conforming facets — `orientDir = scen − mcen = 0` and the
  facet-normal flip test never fires, leaving the sign to winding luck; pass `-outward` for coincident
  interfaces (the shipped c2_1 tests always do). Candidate warn/refuse in the review PR-5 batch.
- **Validation choke points (previously silent):** duplicate CONTACT tags now refused across
  contact/mortar/contactPlane (`contactTagInUse` — the tag is the leading key of every path-state
  store: a duplicate ALIASED friction slots last-writer-wins and ping-ponged the re-emit fingerprint
  every handle); `kt<0` refused (mirrors the kn guard); faceted-surface connectivity must be a whole
  number of segments (`size % nps == 0` — the trailing partial segment was silently DROPPED);
  a missing node or a 2D (`-ndm 2`) node in a contact surface now SKIPS the whole contact LOUDLY at
  handle() (`ladrunoSurfaceNodesOk` — previously silent dead pairs, and a 2D node read `getCrds()(2)`
  OUT OF BOUNDS, unchecked in release).
- **Mortar `inverseIsomap2D` (the PR-1 gate follow-up):** absolute `tolR=1e-13` on a LENGTH-unit
  residual (aux-plane UV ~ facet size h, noise floor ~eps·h) ⇒ GPs silently SKIPped for h ≳ 2000
  (378/400 dead at h=5e4) ⇒ `-mortar` contact and `LadrunoTie -mortar` quietly lost their integration
  at mm-unit-building scales. FIX = the PR-1 parametric-step escape (`|dξ|+|dη| < 1e-10` after the
  update). In-solver gate: `tests/test_contact_review_p3_validation.py` mortar press at h=1 AND h=5e4
  settle to the SAME penalty prediction (pre-fix: free fall at h=5e4).
