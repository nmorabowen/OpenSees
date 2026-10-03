---
wp: LEGACY
title: "Contact-review P5 (2026-07-02): per-contact NormalField sign re-key + the hygiene batch"
date: 2026-07-02
legacy_seq: 130
---
### Contact-review P5 (2026-07-02): per-contact NormalField sign re-key + the hygiene batch
- **Shared-master `-smoothNormal` sign (MED, the headline):** the ADR-63 nodal-normal field was keyed
  by MASTER SURFACE tag, but its frozen global outward SIGN is a per-CONTACT datum (parity with the
  faceted per-contact `orientDir`). A SECOND contact sharing the master (a two-sided baffle: slaves on
  BOTH faces of one declared plate) inherited the FIRST contact's frozen sign — `signCaptured` skipped
  the whole vote, so even its own explicit `-outward` was IGNORED — and its slaves were read as
  penetrating from the wrong side and EJECTED through the plate (silent pass-through). FIX =
  `theNormalFields` keyed by CONTACT tag (each contact votes + freezes its own sign; the topological
  σ/sharedEdge cache duplicates per contact on a shared master — cheap, fingerprint-cached). Gate:
  `tests/test_contact_review_p5_percontact.py::test_p5_two_sided_baffle_per_contact_sign`.
- **One-time warnings were PROCESS-lifetime:** every contact warning latch was a function-local
  `static bool warned` ⇒ after the first warning anywhere, every later model in the same interpreter
  (wipe/new model) stayed SILENT, and a second contact with the same defect was silenced by the first.
  FIX = per-(contactTag, topic) latches on the engine (`LadrunoContactDomain::warnOnce`, reset
  naturally on wipe since `Domain::clearAll` deletes the engine; kept across `revertToStart` — a reset
  re-run repeats the same configuration). The rigid-plane adapter now carries its contactPlane tag for
  the latch key.
- **`-visc` was NOT statics-inert:** the "v≡0 in statics ⇒ byte-identical" claim ignored that trial
  velocities are STATE — a static stage AFTER a transient one (or `setNodeVel`) keeps the last
  committed velocities, so the dashpot injected a constant unphysical force with NO matching tangent
  (`StaticIntegrator::formEleTangent` never calls `addCtoTang`) and silently shifted the static
  equilibrium. FIX = the dashpot is DISABLED under a `StaticIntegrator` (dynamic_cast gate in
  `LadrunoContactFE::viscousActive`, warn once per contact); genuinely-zero velocities made the old
  force exactly 0.0, so the gate is byte-identical for a pure static run (the shipped D2 static test).
- **Edge-edge first-capture sign could coin-flip:** the E2 body-fixed sign is captured ONCE from
  `sign(n·orientDir)` — with the reference nearly ⟂ n the capture was a numerical coin flip that then
  BOUND the pair for its whole life. FIX = the H2-style conditioning guard (mirrors
  `normalOriented`'s `|proj| < 1e-12·|ref|` refusal) DEFERS the capture (the eval is treated as
  no-contact; retried at the next, better-conditioned config) — and the defer is LOUD (gate lens-A):
  on a symmetric rig the reference can stay ⟂ n forever (e.g. an axis-aligned `-outward` exactly in
  the crossing plane), so a one-time per-contact warning names the pair inert and points at
  `-outward`. Also closes the degenerate-normal hole where `edgeGap`'s failure path left `n`
  uninitialized but the caller kept assembling.
- **`LadrunoEdgeKernel` TAU_LEN was an absolute floor (the PR-1/PR-3 unit-trap class):** `a ≤ 1e-9` on
  a SQUARED length ⇒ every edge shorter than ~3e-5 length units was silently DEGENERATE (edge-edge
  contact dead at small length units). FIX = RELATIVE floor `TAU_LEN_REL = 1e-12` (squared-length
  ratio vs the LONGER edge of the pair; both-zero ⇒ degenerate). Oracle mirrors
  (`proto_e0_closest_point.py` / `proto_e2_penalty.py`) updated in lockstep, all green.
- **Tri/quad shared-edge guard 2× factor:** `onSharedInteriorEdge`'s parametric `edgeTol` (0.1) was
  calibrated on the QUAD span of 2 (5% of span) but applied verbatim to the TRI span of 1 (10%) — the
  P2.1 ownership guard fired twice as aggressively on tri meshes. FIX = tri branch uses `edgeTol/2`
  (same fraction-of-span; quad byte-unchanged — the calibrated P2.1 ridge gates are quad meshes).
- **Coincident/staggered zero-gap mortar AUTO orientation (the PR-3 gate follow-up):** the kernel
  signs the facet normal by flipping it toward `orientDir = scen − mcen`, so the hazard is the
  NORMAL COMPONENT of that vector vanishing — a COINCIDENT conforming pair (`orientDir ≈ 0`) AND a
  zero-gap STAGGERED non-conforming pair (`orientDir` tangential, O(h) — the canonical mortar
  interface; gate lens-B) both leave the flip test to roundoff ⇒ the mortar sign is WINDING LUCK
  (an unlucky winding = attractive contact, no diagnostic). Now WARNS once per contact when
  `|n̂_m·orientDir| ≤ 1e-6·h` at a proximate pair (`|orientDir| ≤ 2h`; far coplanar pairs stay
  silent); `-tie` exempt — its full 3-vec bond is sign-independent and coincident is its common
  case. Pass `-outward` for zero-gap interfaces.
