---
title: "ADR 72 P2 — signed test specs (task 2.2 gate-semantics review)"
status: SIGNED 2026-07-13 (Fable orchestrator) — implementation (2.3) may start; any
  deviation from an assertion below must be reported back for re-adjudication,
  not silently adapted
amendments: re-signed 2026-07-13 during implementation (same signer) —
  (1) S0 mass equality relaxed 1e-14→1e-12 — M is extracted as (K+c3M)−K from
  printA pairs, so the floor is |K|·eps cancellation noise, not eps on M (the
  in-element M is bitwise formulation-independent by code path);
  (2) S1 subspace metric restated as sin(theta_max) < 1e-7 via SVD — arccos
  cannot resolve 1e-8 rad in double precision;
  (3) S8 advisory-under-uri limited to the never-fires branch — the fires-once
  branch is pinned by the P1 battery and the PROCESS-once static makes it
  unrepeatable in a shared pytest process (code path formulation-independent);
  (4) S9 gate set at ratio < 0.8 — the eleResponse round-trip adds equal
  absolute overhead to both sides, diluting the pure-assembly 0.3–0.4×;
  (5) S5 re-signed from measured physics (4×1×1 probe, mesh-refinement sweep
  n=4/8/16): the ν=0.3 ≥0.98 control and the ν=0.4999 uri∈[0.95,1.05] pin were
  over-tight — at n=4 BOTH formulations sit ≈0.95–0.965 at ν=0.3 (shared
  discretization/reference error, gap 0.012, both → ~0.99 on refinement), and
  at ν=0.4999 uri=0.891/std=0.749 (uri halves the error but is "not a mixed
  element", ADR §3.4 verbatim). Re-signed gates pin the CONTRAST: ν=0.3 both
  ≥0.95 AND |uri−std| ≤ 0.03 (no formulation gap when ν is moderate);
  ν=0.4999 std ≤ 0.90 (the lock — escalation pin UNCHANGED), uri ≥ 0.85, and
  uri − std ≥ 0.10 (the relief). Measured margins: 0.965/0.953 (gap 0.012);
  0.749 ✓, 0.891 ✓, gap 0.142 ✓. NOT an ADR §3.4 errata — §3.4 already
  states partial relief; the spec pin was stricter than the ADR's own claim;
  (6) 2026-07-19 (same signer): S9 gate relaxed to < 0.95 — a loaded-box
  full-suite replay measured 0.834 on ~20 ms wall-clock samples
  (noise-dominated under contention; the 0.8 gate was CI-safe but not
  load-safe). The assert now only pins "uri is cheaper than std"; the
  measured ratio stays REPORT (0.63 idle-box, guide records it).
---

# ADR 72 P2 `-formulation uri` — signed test specifications

Companion to [[_adr72_implementation_plan]] task 2.2 and the ADR §6 P2 gate row.
These specs ARE the §3.2 spurious-mode contract; a wrong-but-passing test here
is worse than none. Every assertion below is a *raising* pytest assertion, not
a logged number, unless explicitly marked REPORT.

**The P2 anchor.** Upstream has no reduced-integration 20-node element, so the
P1 reduce-to anchor does not exist for `uri`. The correctness anchor is the P0
sympy oracle (`tests/hex20_reference.py`): numpy-assembled K/f/σ at the GP8
table (z-fastest lexicographic, the Brick order) vs the element via printA /
recorders, on BOTH the unit cube and the P1 distorted-hex fixture. Everything
else (census, communicability, pathology, bending, locking) is *contract
semantics* layered on top of that anchor.

**Shared fixtures.** Reuse the P1 harness (`_hex20_nodes`, `_static_rig`,
printA K-extraction, distorted-hex corner set). Debt (c) applies: import node/GP
tables from `tests/hex20_reference.py` — no re-transcription. Elastic material
E=1000, ν=0.25 unless stated. Eigen census = `scipy/numpy eigh` on the printA
matrix (symmetric), zero-mode count = #{λ_i : λ_i < 1e-12·λ_max}.

---

## S0 — oracle anchor (headline, before any contract test)

- **Mesh:** single element, (a) unit cube, (b) the P1 distorted hex.
- **Assert:** element `uri` K (printA) == oracle K assembled at GP8, rel-max
  ≤ 1e-12, both meshes. Resisting force under an imposed nodal displacement
  vector (same harness as P1 reduce-to force test) == oracle B·σ integration,
  ≤ 1e-12. Per-GP stress (element recorder / response, 8 stations) == oracle
  σ(GP_L), ≤ 1e-12, station-by-station in GP8 order — this pins the
  materialPointers[L] ↔ GP8[L] pairing.
- **Mass:** `uri` consistent M (printA under Newmark, P1 pattern) == `std` M
  == upstream `Twenty_Node_Brick` M, ≤ 1e-14 rel. Mass is 27-pt and
  formulation-independent BY DESIGN — this is the assertion that proves it.

## S1 — mode census (ADR gate: "single element: rank 48, exactly 6 zero-energy modes")

- **Mesh:** single unit-cube element, unrestrained, `uri`.
- **Assert:** eigen census of K = **exactly 12** zero modes (6 RBM + 6
  spurious); separation λ12/λ13 ≥ 1e10 (P1 saw 6e13 for the std 6/7 split).
- **Catalogue:** project the 12-dim null basis onto the 6-dim rigid-body space
  (built from nodal coords); the 6-dim orthogonal complement = the spurious
  modes. Assert the complement's principal angles vs the P0 oracle's
  8-pt null-space complement (recomputed via hex20_reference, NOT a checked-in
  dump) are all < 1e-8 rad. This is the "catalogued" clause — same modes as P0,
  through the element path.
- **Control:** same census on `std` = exactly 6 (already a P1 test; re-assert
  here side-by-side so the contrast lives in one test).
- **Distorted hex:** census repeats on the distorted fixture — still exactly 12
  (rank deficiency is topological, not geometric). REPORT the 12 eigenvalues.

## S2 — non-communicability pin (ADR gate: "restrained ≥2×2×2 block: zero spurious global modes")

- **Mesh:** 2×2×2 block (8 `uri` elements, unit cubes), shared serendipity
  nodes deduped by coordinate (extend the P1 two-element mesher).
- **BCs:** *statically determinate* 6-DOF restraint (3-2-1 scheme on three
  corner nodes of one face) — deliberately minimal so nothing but element rank
  can remove mechanisms. NOT a fixed face (a fixed face could mask a mechanism).
- **Assert:** eigen census of the restrained global K = **zero** modes below
  1e-10·λ_max; AND λ_min(uri) ≥ 0.05·λ_min(std on the identical mesh+BCs) —
  the softest uri mode is a physical mode, not a near-mechanism.
- **Single-element contrast (the Bathe Fig. 5.46 half):** one `uri` element
  under the SAME 3-2-1 determinate restraint → census shows ≥ 1 (expect 6)
  zero modes = outright unstable; `std` under the same restraint → zero.
  The pair of assertions IS the non-communicability statement: deadly alone,
  cured by any 3D assembly.

## S3 — single-stack pathology (ADR gate: "demonstrated & documented — the honest counter-example")

- **Mesh:** 1×1×4 stack of `uri` unit cubes (one element thick both transverse
  directions — the documented propagation topology), base face z=0 FULLY fixed
  (realistic cantilever usage, not a contrived support).
- **Primary assert (eigen):** λ_min(uri-stack)/λ_max(uri-stack) < 1e-3 ·
  [λ_min(std-stack)/λ_max(std-stack)] — i.e. the uri stack carries a
  near-mechanism at least three orders softer than the softest physical mode
  of the identical std stack. Expect machine-zero-ish; the 1e-3 floor keeps
  the assertion honest if the fixed base leaks slight stiffness into the mode.
- **Mode-shape documentation:** REPORT (test log, not assert) the offending
  eigenvector's tip-face pattern (corner vs mid-edge alternation) — the living
  documentation clause.
- **std control:** λ_min(std-stack) corresponds to the physical cantilever
  bending mode: REPORT its Rayleigh-quotient frequency vs beam theory (sanity,
  ±20% — coarse mesh).
- **Escalation trigger (task 2.4):** if the uri stack shows NO soft mode
  (ratio ≥ 1e-3), ADR §3.2's propagation claim is contradicted → STOP, full
  adversarial gate + ADR errata. The test failing IS the alarm; do not weaken
  the assertion to make it pass.

## S4 — coarse bending (ADR gate: Oñate Fig. 8.23 replication, tip ≥ 0.98·analytic)

- **Mesh:** cantilever L=10, b=h=1, **ν=0** (kills anticlastic curvature so
  beam theory is the honest reference), E=1000. 5×1×1 `LadrunoBrick20 uri`
  (one element through thickness — the coarse-bending claim).
- **Load:** end shear P as the *consistent* nodal distribution of a uniform
  traction on the tip face (use eleLoad/nodal equivalents from N; document).
- **Reference:** Timoshenko tip deflection δ = PL³/3EI + PL/(κGA), κ=5/6.
- **Assert:** mean tip-face deflection(uri) ≥ 0.98·δ; also run `std` same mesh
  — expect it passes too (quadratic bending is the element's raison d'être;
  the contrast is elsewhere) — assert std ≥ 0.98 as well.
- **Comparative table (REPORT + two coarse asserts):** same cantilever meshed
  with 5×1×1 `LadrunoBrick std` (H8) → assert < 0.90·δ (parasitic-shear lock
  demonstrated); `LadrunoBrick eas` 5×1×1 → REPORT (expect between; "needs
  refinement" is the story, don't over-pin a fork-internal number).

## S5 — near-incompressible relief (ADR gate: "ν=0.4999: uri relieves, std locks")

- **Mesh:** same cantilever family as S4 (bending is where the volumetric
  constraint bites the serendipity hex; a fully confined block has a linear
  exact solution and shows nothing), 4×1×1, E=1000, **ν=0.4999**.
- **Reference:** same Timoshenko δ(ν=0.4999) (G from E,ν).
- **Assert (amendment 5 — re-signed from the measured 4×1×1 physics):**
  ν=0.3 control: both formulations ≥ 0.95·δ(0.3) AND |uri−std| ≤ 0.03 (no
  formulation gap at moderate ν — the shared ~3.5% shortfall is coarse-mesh
  discretization + beam-theory reference error, converging to ~0.99 at n=16);
  ν=0.4999: std/δ ≤ 0.90 (the lock, per AUG §28.1.1 + the r≈0.44 constraint
  count — the ESCALATION pin, unchanged), uri/δ ≥ 0.85, uri−std ≥ 0.10 (the
  relief contrast — uri halves the error but is not a mixed element, §3.4).
- **Escalation trigger:** if std does NOT lock (>0.90) at this mesh, that
  contradicts ADR §3.4 — report to orchestrator (2.4 decision: try the
  distorted variant / more slender before declaring errata; do NOT silently
  tighten ν or slenderness until it passes — the spec change must be
  re-signed).
- **REPORT:** checkerboard diagnostic — GP-pressure alternation across
  neighboring GPs for std at ν=0.4999 (ADR §3.4's diagnostic), logged only.

## S6 — Barlow superconvergence (ADR gate: "uri GP stresses beat std GP stresses vs analytic bending")

- **Mesh:** the S4 ν=0 cantilever (5×1×1), both formulations.
- **Field:** σ_xx(x,y) = M(x)·y/I from beam theory at each GP's GLOBAL
  coordinates (element GP order via the recorder round-trip of S7).
- **Assert:** RMS relative error of σ_xx over interior elements' GPs
  (exclude the support element — St-Venant): uri(8 GP) < std(27 GP). Each set
  scored against the analytic field at its OWN stations.
- **REPORT:** both RMS values (expect uri ~O(h³)-class vs std generic).

## S7 — recorder seam, debt (a) (basisInfo/custom rule; formulation-aware GP metadata)

- **uri element → LadrunoRecorder round-trip** (P1 test pattern): exactly **8**
  GP stations; GLOBAL_GP_COORDS == isoparametric map of the GP8 table (z-fastest
  lexicographic) through the mesh geometry, ≤ 1e-12; per-GP σ matches S0's
  closed form. Mechanism: the element answers the `basisInfo` probe (BezierTet10
  pattern — topology=hex, numCtrl=20, numGP formulation-dependent) and/or
  `integrationPoints`/`integrationWeights`; the recorder must NOT fall through
  to the hardwired 33018→GaussLegendre_3 arm when the element says 8.
- **std regression:** the P1 27-GP round-trip must still pass unchanged
  (whatever seam is built, std metadata stays GL3-identical byte-for-byte).
- **Battery hygiene:** the P1 unconditional `NUM_GP==27` pin becomes
  formulation-aware (27 std / 8 uri) — grep-level acceptance: no bare 27 GP
  count left that would silently pass for a uri element.

## S8 — surface & serialization regression

- Parser: `-formulation uri` (and alias `reduced`) now ACCEPTED — replace
  `test_formulation_uri_rejected_until_p2` with an acceptance test asserting
  Print/JSON reports `uri`. `-hourglass` stays a hard error (unchanged test).
  `-lumped` still parses + getMass falls back (P3 unchanged).
- sendSelf/recvSelf round-trip of a `uri` element preserves ordinal 1 — the F8
  wire-coercion defense must be RETIRED for ordinal 1 (it exists to protect
  against claiming-uri-while-computing-std; after P2 uri computes uri). Assert
  the round-tripped element's K == original K (printA), and Print says uri.
- U1 damage advisory: fires identically under uri (probe is
  formulation-independent) — one parametrized case.

## S9 — cost (ADR gate: "assembly time ratio ≈ 0.3–0.4× std") — REPORT + loose gate

- **Bench:** one moderately sized block (e.g. 4×4×4 elements), time N repeated
  tangent-forming analyze steps, uri vs std, same mesh, wall clock.
- **Assert (loose, CI-safe):** ratio < 0.7. REPORT the measured ratio with the
  0.3–0.4 expectation noted. Timing asserts tighter than this are CI flake
  bait; the honest number lives in the log and the guide.

---

## Debts folded into 2.1/2.3 (from ADR §6 P2 row)

- (a) recorder basisInfo seam → S7 (element + recorder + battery pin).
- (b) blocked/sparse BᵀDB via per-node 3×3 blocks — do it while the assembly
  is touched for uri; re-anchored by S0 (≤1e-12 vs oracle) AND the P1 std
  reduce-to suite re-run green (the refactor must not move std by 1 ulp beyond
  the existing gates).
- (c) battery imports `tests/hex20_reference.py` tables (sympy is in Zone-A) —
  applies to the whole P2 file and retrofits the P1 file's private
  `_shape_h20`/`_gp27_brcshl` copies.
