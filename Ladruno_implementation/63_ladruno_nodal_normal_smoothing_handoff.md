# ADR-63 nodal-normal smoothing — handoff (P2 onward)

Companion to [[63_ladruno_nodal_normal_smoothing_adr]] (the design + the folded pre-code gate + the
4-lens P1 review). This is the resume point after **P0+P1 shipped to `ladruno` via PR #457**
(merge `fd82af5b0`, 2026-07-01).

---

## Status snapshot

- **On `ladruno` now:** the `-smoothNormal` opt-in averaged nodal-normal field for the **NTS** contact
  lane. Resolves ADR-60 **R3** (curved-master ridge flip) + the NTS slice of ADR-41 **Q-NORMAL** (C0
  normal across junctions). **OFF by default ⇒ byte-identical; no new classTag.**
- **Q-TANGENT is RESOLVED** (symmetric-first): the MVP ships the symmetric **frozen-field `kn·BᵀB`**
  tangent. The full `∂n_smooth/∂u` is the conditionally-required, gated **P3**.
- **Validation shipped:** oracle `proto_nodal_normal_selfcheck.cpp` 21/21; in-solver
  `tests/test_adr63_smoothnormal_p1.py` 5/5; ADR-39/41/57/60 battery 147 OFF byte-identical.
- **Ledgers/ADRs updated in the PR:** LEDGER_implementations row (status `draft`), LEDGER_vanilla row
  (`-smoothNormal` parser — the only vanilla touch), 2 LEDGER_quirks, ADR-47 #4a marked pulled-forward,
  ADR-60 R3 flipped to resolved. **Banner line NOT added** (held until the feature graduates draft→
  shipped after P2; add a `banner_features.txt` line + `patch_banner.py` then).

---

## Where the code lives (entry points the next session needs)

- **`SRC/domain/contact/LadrunoContactNormalField.h`** (new, header-only, OpenSees-free) — the geometry:
  - `propagateOrientation(mTags,nSeg,nps,sigma)` — flood-fill manifold winding `σ_s` (topological). Returns
    `Status` (OK / NON_MANIFOLD / DISCONNECTED / NON_ORIENTABLE / CLOSED). **The edge-adjacency map built
    here is the key asset P2 reuses for facet ownership** (it already knows which edges are shared vs boundary).
  - `voteSign(nSeg,nps,segCoords,sigma,seed,&conf)` — decides the GLOBAL outward sign ONCE + its confidence.
  - `nodalNormals(mTags,nSeg,nps,segCoords,sigma,G,segNodalNorm)` — area-weighted (Newell) nodal normals with
    the FROZEN sign `G`; relative-degeneracy leaves a cancelling node zero.
  - `smoothNormal(nps,N,nodalNorm,n)` — the shape-fn blend + renormalize; false on a degenerate blend.
- **`SRC/domain/contact/LadrunoContactProjection.h`** — `evalSegmentSmooth(...)` (the smoothed sibling of
  `evalSegment`; facet projection + in-bounds/penetration gate UNCHANGED; falls back to `normalOriented(refDir)`
  on a degenerate blend). **This is where the P2 facet-ownership guard goes.**
- **`SRC/domain/contact/LadrunoContactDomain.{h,cpp}`** — the `NormalField` store per master-surface tag
  (cached `σ` + fingerprint + `segNodalNorm` + `signCaptured`/`globalSign`/`signConf`); `setNormalField` /
  `getSegNodalNorm` / `getNormalFieldSignConf` / `clearNormalFields`; `Contact.smoothNormal` flag.
- **`SRC/analysis/handler/LadrunoContactHandler.cpp`** — the `if (ct.smoothNormal)` field-build block (deformed
  coords, once-captured global `smoothSeed`, refuse + ill-conditioned-sign + downgrade warnings, missing-node
  refuse) + the per-segment `setSmoothNormals` install; the candidate-loop `orientDir = smoothSeed` override.
- **`SRC/analysis/handler/LadrunoContactFE.{h,cpp}`** — `setSmoothNormals` + frozen `nodalNorm[4][3]` +
  `useSmoothNormal`; `segmentActive` branches to `evalSegmentSmooth`; the faceted B3 `∂n/∂u` is suppressed
  under smoothing (`consistentNormal && !useSmoothNormal`).
- **`SRC/interpreter/OpenSeesOutputCommands.cpp`** — `-smoothNormal` parse (NTS-only; refused with `-mortar`).

---

## P2 backlog (priority order)

### P2.1 — Facet ownership at a sharp ridge (HIGHEST — the one real limitation)
The deferred bring-up limitation. A slave sitting AT a sharp convex ridge has its closest-point land on the
**shared edge** of the *adjacent* facet (barycentric ≈ 0, marginally in-bounds) → that neighbor reads a large
spurious penetration → big ejecting force → instability. The faceted path only *incidentally* prunes this (at a
90° ridge the neighbor normal is ~⟂ `orientDir` → the `normalOriented` fail-safe kills it). The blunt
interior-margin gate was **tried and REVERTED** (it also killed free-edge contacts and removed a
spuriously-*helpful* neighbor force → destabilized other geometries).

**Designed fix (not yet built):** *interior ownership by SHARED edge only.* `propagateOrientation` already
builds the edge-adjacency map — it knows which of a segment's edges are **shared (interior, 2 segments)** vs
**boundary (free, 1 segment)**. Thread a per-segment "shared-edge mask" through the field store → into
`evalSegmentSmooth`. Reject the contact **only** when the projection's nearest edge is a *shared* edge and the
barycentric coord is within a small band of it (the neighbor's interior owns that region). A projection near a
*free/boundary* edge is NOT rejected (that facet legitimately owns its free edge — this is why the blunt margin
failed). At the exact ridge both facets briefly reject (a negligible dead-band; the smooth field keeps the
normal continuous either side). Oracle: the tri-3 and quad convex-ridge "press into the facet" cases (the ones
the current tests deliberately avoid) must hold. Cross-check against ADR-57 #4b edge-handoff (this is adjacent).

### P2.2 — Friction under the smoothed normal
`segmentActive` already threads the smoothed `n` into the friction slip `gTvec`
(`LadrunoFrictionKernel::tangentPart(drel,n,gTvec)`) and the friction tangent — so friction *should* compose,
but it is untested with `-smoothNormal + -mu`. Add a dragged-block-over-a-ridge-with-friction gate: traction
continuous across the junction, no spike, energy balance closed.

### P2.3 — Compose with re-emit (lift the R3 `-outward` caveat for the full combo)
The field is already built from the same deformed feed and invalidated by the same R1 membership fingerprint as
re-emit. The remaining work is a **live gate**: `-smoothNormal + -reemit + -mu + curved master` sustaining
contact across ≥3 crossings under sliding (the ADR-60 "exposed combo", now closeable). Watch the frozen global
sign holds across re-handles (it's captured once, keyed by surface tag — verify a re-mesh re-captures correctly).

### P2.4 — Curved-implicit convergence tripwire (from the tangent-posture review)
BLOCKER-TANGENT-POSTURE's P2 gate: run a **genuinely curved** (not piecewise-flat) implicit master under Newton
with `-smoothNormal -consistentNormal` and record iteration count vs the faceted `-consistentNormal` baseline.
If the frozen-field tangent **converges comparably** → P3 stays deferred *with evidence*. If it **stalls** (not
just slows) → promote P3 to required + document curved-implicit as explicit-first until P3.

---

## P3 (gated) — full `∂n_smooth/∂u` consistent tangent
Only if P2.4 trips. `-consistentNormalSmooth`: differentiate the nodal-normal average + the shape-fn
interpolation through the mesh motion → a non-symmetric, 1-ring-wide (~10–17 node) stencil. **Implicit-only —
MUST NOT be assembled under an explicit integrator** (the P2a getTangent-returning-K trap). Mirror the existing
`normalGeomTangent` (B3) oracle discipline: analytic == FD ~1e-8 on a warped patch. Converged answer unchanged
(tangent-independent). See the ADR §D4 / BLOCKER-TANGENT-POSTURE.

---

## Locked decisions — DO NOT re-litigate
- **Symmetric-first tangent** (Q-TANGENT resolved via the 2-iteration adversarial loop; steelman-B conceded).
- **Global sign is FROZEN once** at first OK handle (F1 fix) — never re-decided from a deformed config.
- **Refuse** non-manifold / disconnected / non-orientable / **closed** / missing-node surfaces → faceted `-outward`.
- **NTS lane only.** Mortar GP-normal smoothing + edge-edge (own `signN`) are fenced (own future passes).
- **Facet projection + active set unchanged** by smoothing — only the normal *direction* changes (so ADR-60 R5
  slide-off is untouched). P2.1's ownership guard is the ONE sanctioned active-set change, and it's smoothed-path-only.

---

## Gotchas / lessons (so you don't rediscover them)
- **Local build recipe (worktree has no build tree):** junction `mumps-install/archive/src/build` from the main
  checkout (`New-Item -ItemType Junction`); `set LADRUNO_OPENSEES_QUIET=1` (MANDATORY — the boot `.pth` banner
  corrupts CMake's Python probe); `Ladruno_scripts\build.bat OpenSeesPy`. Incremental rebuild ~minutes.
- **Run tests against THIS worktree's pyd:** `python -S Ladruno_scripts\_run_adr63_tests.py` (the `-S` skips the
  boot `.pth` that would pre-import a *different* worktree's stale pyd; the runner wires `dist/bin` + site-packages
  by hand). The runner is UNTRACKED (machine-local, ADR-60 convention) — recreate it if missing (recipe in the
  ADR-60 memory note).
- **Debugging confounder (cost me two rebuilds):** a process-global `static int` debug counter in `segmentActive`
  spanned TWO `ops` models in a probe that ran outward-then-smooth → the `use=0` lines were the OUTWARD run;
  worse, the freelist reused the SAME `this` addresses across `ops.wipe()` → false pointer-match. **Debug the
  smooth case IN ISOLATION** (one model per process).
- **The "holds the ridge" test can't distinguish** smooth-field-active from the faceted+global-seed fallback
  (both hold). Use the `use=` flag / a direct normal probe to confirm the smooth path is actually engaged.
- The master centroid for the auto seed must be over **UNIQUE** nodes (the flat `mTags` double-counts shared
  ridge nodes → biased centroid → sign flip). Already fixed; don't regress.

---

## Oracle + tests
- Build-free: `Ladruno_implementation/contact_prototypes/proto_nodal_normal_selfcheck.cpp`
  (`g++ -std=c++11 -I ../../SRC/domain/contact ... && ./a.out`, 21/21).
- In-solver: `tests/test_adr63_smoothnormal_p1.py` (5/5): quad convex-ridge R3 fix + faceted pass-through
  companion; flat-master byte-consistency; tri-3 chain + consistency (flat, to dodge the P2.1 pathology).
