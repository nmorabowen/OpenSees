---
title: "ADR 72 — implementation plan (P0–P3): task breakdown, Opus-4.8 delegation map, review gates"
status: plan — companion to [[72_ladruno_second_order_brick_adr]]; update the checkboxes as PRs land
---

# ADR 72 implementation plan — LadrunoBrick20

Execution companion to [[72_ladruno_second_order_brick_adr]] (which owns the
*decisions*; this file owns the *work*). Four PRs off `ladruno`, one per phase,
strictly sequenced P0 → P1 → P2 → P3.

**Delegation legend** (per [[feedback_prefer_opus48_agents]], updated for the
Fable tier):

| Mark | Meaning |
|---|---|
| **[OPUS]** | Opus 4.8 subagent implements solo — a hard correctness anchor (upstream reduce-to, numpy oracle, exact pinned fractions) catches its mistakes, so no human/Fable review of the code itself is needed beyond the anchor passing |
| **[OPUS-V]** | independent Opus 4.8 *verifier* agent — a SEPARATE agent from the implementer, run where a defect would silently poison the anchors themselves |
| **[FABLE]** | main-session work: orchestration, sign-offs, gate-semantics review, PR/CI shepherding, conditional adversarial escalation |
| **[USER]** | a decision only the user can make |

The principle: **implementation is delegable wherever an anchor is stronger
than a reviewer; review effort goes to the few places that anchor everything
else** (node/GP ordering, oracle independence, test semantics) — a defect there
passes every downstream test while making them meaningless.

---

## 0. Standing constraints (every phase — bake into every subagent prompt)

- **Build**: from PowerShell, `cmd /c Ladruno_scripts\build.bat OpenSees`
  (never via the Bash tool — [[project_bezier_charlen]]); all-5-targets only
  with NO args; close VS Code if DLLs lock ([[project_installer_dll_lock]]).
- **Test env**: worktree runners MUST use `python -S` + manual paths +
  `assert opensees.__file__` (the boot-`.pth` stale-pyd trap,
  [[project_opensees_test_env]]). Mirror `tests/test_ladrunoBrick_element.py`
  harness patterns (K/M extraction, distorted-hex fixtures).
- **PR hygiene**: one PR per branch; verify `state==OPEN` before every push
  ([[feedback_stranded_commits_after_automerge]]); auto-merge gates are NOT the
  Zone-A battery — **WATCH `gh pr checks`** until green
  ([[feedback_ladruno_ci_gating_and_openseespy_peek]]); PR body via
  `--body-file -` ([[feedback_gh_pr_body_stdin]]).
- **New files**: stamp headers (`stamp_headers.py`, add to GLOBS —
  [[feedback_always_stamp_header]]); ledgers in the SAME PR
  ([[LEDGER_implementations]] + [[LEDGER_vanilla_files]] for every vanilla
  touch, `// Ladruno` comments); banner line at P1 ship
  (`banner_features.txt` → `patch_banner.py`).
- **Code conventions to inherit** (from the LadrunoBrick adversarial history /
  [[LEDGER_quirks]]): number-aware option parsing =
  `OPS_GetDoubleInput` + `OPS_ResetCurrentInputArg(-1)`, never
  `OPS_GetString`+`strtod` (openseespy PyFloat trap); materials cloned via
  `getCopy("ThreeDimensional")` (generic `getCopy()` may `exit(-1)` —
  [[project_initdefgrad_staged]]); `Formulation` enum ordinals {STD=0, URI=1}
  are **serialized — append-only forever**; loop bounds from one `NGP`
  constant (the upstream `setParameter i<4` bug class).
- **Tag frontier**: 33018 assumed; **re-check at P0 land time** against ADR-71
  `LadrunoUP` (33017, drafted same day on branch `guppi/focused-gould-1fa488`)
  — whichever P0 lands second takes the next free slot and fixes its ADR.

---

## P0 — pure kernel + oracle + tag reservation (PR-1, size S/M)

Goal: `LadrunoHex20Shape.h` proven against an independent oracle before any
element code exists. ADR §6 P0 row = the gate list; this section = the how.

| # | Task | Owner |
|---|---|---|
| 0.1 | **Extract the canonical orderings**: read the 20-node natural-coordinate node order and the 27-pt GP order + weights out of upstream (`Twenty_Node_Brick.cpp` / `shp3dv.cpp` `brcshl`) and the fork recorder tables (`Ladruno_ElementResults.h:1266` rule + `:1418` node block). Produce a one-page ordering memo (node k ↔ (ξ,η,ζ); GP L ↔ (ξ,η,ζ,w)) committed as a comment block in the kernel header. | **[OPUS]** |
| 0.2 | **Ordering cross-verification** — independent agent re-derives the memo from the sources WITHOUT seeing 0.1's memo, then diff. Any mismatch blocks the phase. *This is the poison-risk item: a wrong ordering passes every kernel self-test and silently breaks the P1 reduce-to + recorder pairing.* | **[OPUS-V]** (blocking) |
| 0.3 | **`SRC/element/ladrunoBrick/LadrunoHex20Shape.h`** — pure doubles, no OpenSees deps (the `LadrunoMassLumping.h` idiom): serendipity N/∂N∂ξ at arbitrary ξ (ADR §3.1 formulas), node table (0.1 order), 27-pt (brcshl order) + 8-pt (`Hexahedron_GaussLegendre_2` lexicographic, z-fastest — the Brick order) tables, Jacobian/inverse/det, 6×60 B fill (Voigt {xx,yy,zz,xy,yz,zx}), consistent-mass block, row-sum + HRZ lump helpers (HRZ via `Ladruno::hrzLumpRaw` reuse). | **[OPUS]** |
| 0.4 | **Numpy oracle** — layout per ADR-70 P0 (`git show 6347ab76a --stat` for the pattern): shape functions derived **symbolically (sympy) from the serendipity definition, NOT transcribed from the C++**; runs the full ADR §6 P0 gate list: partition-of-unity / Kronecker / Σ∂N=0 / degree-2 completeness; distorted-hex J vs numpy; K rank 54@27pt, 48@8pt + the 6 spurious eigenvectors dumped & catalogued; mass diagonals >0; **HRZ cube fractions == 7/248 & 2/31 exactly; row-sum corner == −M/8 demonstrated**; GP tables == 0.1 memo; C++-kernel-vs-oracle to ~1e-12 via the standalone g++ driver. | **[OPUS]** |
| 0.5 | **Oracle-independence audit** — verifier reads the oracle and confirms zero formula transcription from the kernel header (sympy derivation only), and that every ADR §6 P0 gate is actually asserted (not just computed). | **[OPUS-V]** (blocking) |
| 0.6 | 14-pt Irons rule: numerical rank check in the oracle (confirm rank-sufficient → keep the reserved `i14` note; else edit ADR §7). | **[OPUS]** (fold into 0.4) |
| 0.7 | Tag reservation `#define ELE_TAG_LadrunoBrick20 33018` + comment; adjust the "33017-33019 free" comment; LEDGER_vanilla_files row. **Frontier re-check vs ADR-71 first.** | **[FABLE]** |
| 0.8 | PR-1 assembly, CI watch, plan-file checkbox update. | **[FABLE]** |

- [x] P0 landed via [PR #548](https://github.com/nmorabowen/OpenSees/pull/548) (R1 ✔ ordering blind-verify; R2 ✔ oracle audit incl. Patch-1 hardening; oracle 7/7; C++ cross-check 8/8; 14-pt verdict: rank-sufficient, `i14` note stands)

---

## P1 — element `-formulation std` + registration + recorder (PR-2, size L)

Goal: the 27-pt element, correctness anchored by the **reduce-to
`Twenty_Node_Brick` ~1e-12** gate. ADR §6 P1 row = gates.

| # | Task | Owner |
|---|---|---|
| 1.1 | **`LadrunoBrick20.{h,cpp}`** — copy the `LadrunoBrick` skeleton, resize to 20 nodes / `NGP=27`, **strip** formulations to `std` (URI enum present, parser-rejected until P2), strip geom seam / hourglass / ssp / eas machinery entirely. Keep & resize: state cycle over `materialPointers[27]`, `formResidAndTangent` = plain BᵀσB Gauss loop via the P0 kernel, consistent mass (27-pt, formulation-independent), `zeroLoad`/`addLoad(SelfWeight)`/`appliedB`, `setResponse` tree (`material`/`stresses`/`strains` per-GP + `stress3D6`/`strain3D6`), `getCharacteristicLength`=∛V, `Print` JSON, `sendSelf`/`recvSelf` (mirror LadrunoBrick layout: nodes, formulation ordinal, massType, b[], per-GP materials + Damping), `setDamping` with `Damping* theDamping[27]`, `displaySelf` (corner-brick render, copy upstream), fixed `setParameter` loop bounds. | **[OPUS]** |
| 1.2 | **`OPS_LadrunoBrick20.cpp`** factory — 20 nodes + matTag + `-formulation {std}` + `-lumped` (accepted, routes to massType; HRZ path activates P3 — until then error "not yet implemented" to keep API honest) + `-b` + `-damp`; number-aware parsing per §0; `-hourglass` → hard error with the ADR §2.2 one-line rationale. | **[OPUS]** |
| 1.3 | **Registrations** (vanilla, strictly additive, `// Ladruno`): `FEM_ObjectBrokerAllClasses.cpp`; interpreter + Tcl element dispatch; `ladrunoBrick/CMakeLists.txt`; `vtktypes[33018]=*_QUADRATIC_HEXAHEDRON` in VTK/VTKHDF/PVD/Gmsh recorders; fork-owned `Ladruno_ElementResults.h` tag lists (rule `Hexahedron_GaussLegendre_3`). | **[OPUS]** |
| 1.4 | **Zone-A battery `tests/test_ladrunoBrick20_element.py`** — the ADR §6 P1 gates: reduce-to (K, resisting force, consistent M ~1e-12, distorted hex, mixed loads, per-GP stress pairing), 2-element distorted patch test (all 6 strain components), rank-54/6-RBM, SelfWeight/`-b` vs closed form **incl. the corner-share sign assertion**, lch=∛V, LadrunoRecorder round-trip (GP coords vs `GLOBAL_GP_COORDS`), sendSelf/recvSelf round-trip, parser-rejection cases (`-hourglass`, `uri`, `-lumped`). | **[OPUS]** |
| 1.5 | **Diff review before PR**: `/code-review` (medium) + the vanilla-footprint checklist (every vanilla file additive-only + ledgered + `// Ladruno`; serialization ordinal freeze noted in the header comment). Per [[feedback_adversarial_gate_when]] NO full adversarial gate — the reduce-to anchor is stronger than a reviewer here. | **[FABLE]** (runs the review; findings back to the Opus implementer) |
| 1.6 | Banner line + LEDGER_implementations row + guide stub (`LadrunoBrick20` section in [[LadrunoBrick_reference]] or its own guide — decide by size); apeGmsh contract row TBD→command+ordering note. | **[OPUS]** |
| 1.7 | **[USER] decision — RESOLVED 2026-07-10: BOTH.** Guide documents the softening caveat AND the element emits a one-time `opserr` advisory at `setDomain` when the attached material exposes a "damage" response (probe via the cached-Response pattern from `LadrunoBrick::damageResponse`; advisory only, run proceeds). Folded into tasks 1.1 (element) + 1.4 (a test asserting the advisory fires once for ASDConcrete3D-class materials and never for elastic/J2) + 1.6 (guide). | ~~[USER]~~ done |
| 1.8 | PR-2 assembly, build + battery green locally, CI watch. | **[FABLE]** |

- [ ] P1 landed (PR #___)

---

## P2 — `-formulation uri` (PR-3, size M)

Goal: the 2×2×2 production configuration, with the spurious-mode contract made
legible. ADR §6 P2 row = gates.

| # | Task | Owner |
|---|---|---|
| 2.1 | uri path: 8-pt loop over the P0 kernel table, `materialPointers[8]`, factory accepts `uri` (ordinal 1), mass STAYS 27-pt, per-GP response tree sized 8, recorder rule `Hexahedron_GaussLegendre_2` keyed by formulation. | **[OPUS]** |
| 2.2 | **Gate-semantics review of the P2 test designs BEFORE implementation** — the mode census (exactly 6 zero-energy eigenvectors, catalogued), the restrained-block non-communicability pin, the **single-stack pathology demonstration**, the Barlow superconvergence comparison, and the ν=0.4999 relief test. *A wrong-but-passing test here is worse than none: these tests ARE the ADR §3.2 contract.* Deliverable: a half-page test spec each (mesh, BCs, load, assertion, tolerance) signed off before 2.3 starts. | **[FABLE]** (blocking) |
| 2.3 | Implement the battery per the signed specs + the coarse-bending Oñate Fig. 8.23 replication (1-layer cantilever ≥0.98·analytic; `LadrunoBrick std` locks, `eas` needs refinement — comparative table in the test log) + assembly-cost timing (~0.3–0.4× std). | **[OPUS]** |
| 2.4 | **Conditional escalation**: if the mode census, communicability, or pathology results contradict ADR §3.2 predictions → STOP, full Opus adversarial gate on the uri formulation + ADR errata. Otherwise skip (per ADR gate policy). | **[FABLE]** (trigger decision) |
| 2.5 | Guide: the std-vs-uri selection table (eigen/single-stack/point-loads/soft-support → std; smooth production ≥2 elements thick → uri). `/code-review` (low) on the diff. | **[OPUS]** impl, **[FABLE]** review |
| 2.6 | PR-3 assembly + CI watch. | **[FABLE]** |

- [ ] P2 landed (PR #___)

---

## P3 — dynamics: HRZ `-lumped` + explicit proof-of-life (PR-4, size M)

Goal: positive lumped mass through the element path; one honest explicit run.
ADR §6 P3 row = gates.

| # | Task | Owner |
|---|---|---|
| 3.1 | `-lumped` = HRZ: `getMass()` massType-1 path builds the 27-pt consistent block then `Ladruno::hrzLump` with the all-translational `dofDir` ID; serialize massType (already in the P1 layout). Un-error the factory `-lumped`. | **[OPUS]** |
| 3.2 | Dynamics battery: eigen bracketing (consistent above / HRZ below analytic bar+beam frequencies), element-path HRZ vector == P0 cube fractions, `criticalTimeStep()` returns the 60-DOF pencil value, explicit wave bar under `CentralDifferenceLadruno` at pencil Δt (stable; **measured Δt ratio vs equal-node H8 mesh reported** ≈0.8 ballpark), energy-balance closure with the hourglass channel asserted zero. Consult the explicit-dynamics skill + [[project_ladruno_mass_scaling]]/[[project_energy_balance_feature]] quirks (velocity-IC gates, Δt reseed) when writing the wave test. | **[OPUS]** |
| 3.3 | Guide finalize: the "explicit permitted-but-discouraged" wording (ADR §3.6), HRZ accuracy caveat (Cook p. 373), pointer to Bézier/H8-uri for real explicit work. `/code-review` (low). | **[OPUS]** impl, **[FABLE]** review |
| 3.4 | PR-4 assembly + CI watch; move ADR §9 log entries; consider the ADR → `Ladruno_internal/implemented_*` move once P4 items are formally parked. | **[FABLE]** |

- [ ] P3 landed (PR #___)

---

## Review map (the one-glance answer)

**Opus 4.8 builds solo** (anchors carry correctness): 0.1, 0.3, 0.4, 0.6 —
1.1–1.4, 1.6 — 2.1, 2.3, 2.5 — 3.1–3.3. That is ~90% of the code and tests.

**Requires review — blocking:**

| Gate | What | Who | Why it can't be skipped |
|---|---|---|---|
| R1 (P0) | node/GP ordering cross-verification (0.2) | independent Opus verifier | a silent mismatch passes every kernel test and poisons reduce-to + recorder pairing |
| R2 (P0) | oracle-independence audit (0.5) | independent Opus verifier | an oracle transcribed from the implementation proves nothing |
| R3 (P1) | `/code-review` medium + vanilla-footprint checklist (1.5) | Fable | vanilla touches + serialization layout are the only P1 surfaces the reduce-to anchor doesn't cover |
| R4 (P2) | gate-semantics sign-off of the five P2 test specs (2.2) | Fable, before test code | these tests ARE the spurious-mode contract; wrong-but-green is the worst outcome |

**Conditional review:** R5 (P2, 2.4) — full Opus adversarial gate ONLY if the
mode census / pathology results contradict ADR §3.2.

**Light review:** `/code-review` low on the P2/P3 diffs (2.5, 3.3).

**User decisions:** U1 — softening-material advisory: warn vs docs-only (1.7,
due at P1). U2 — PR merges are auto once CI is green; no other user gate.

---

## Orchestration log (loop state — newest first)

- **2026-07-13, iteration 11: #564 MERGED (8e8b623b3) — P0+P1+Gmsh runway
  COMPLETE, loop CLOSED.** Shipped this loop: #548 (P0 kernel+oracle+33018),
  #561 (P1 std element + battery + R3 hardening), #564 (GmshRecorder hex20
  type-17 + permutation, incl. upstream Twenty_Node_Brick). Small debt: the
  two #564 ledger rows still read "_PR pending (GmshRecorder hex20 fix)_" —
  fold into the next ledger-backfill batch (the #560 convention). Remote
  branches guppi/xenodochial-lamarr-adr72-p1 + guppi/gmsh-hex20-fix left in
  place (squash-merged). **P2 (`uri`) is NOT started** — fresh phase, gates +
  R3-deferred debts pinned in ADR §6 P2 row; P3 dynamics gates incl. the
  betaK-Rayleigh clobber assert; P4 demand-gated.

- **2026-07-13, iteration 10: P1 MERGED (#561, squash 4fa7b8ff2) → Gmsh PR-3 = #564 in CI.**
  Three CI rounds on #561: (r1) stale-branch conflict vs 8 landed PRs →
  merged ladruno, keep-both broker/Tcl rows, banner regen — NB the 4
  "whole-file conflicts" were line-ending artifacts, resolved by taking
  --theirs + re-applying our additive lines; (r2) Zone-A FAIL =
  **ASDConcrete3D HardeningLawStorage poisoning** (static store-if-absent
  registry keyed by material TAG, survives wipe(); our advisory test's tag-1
  toy law poisoned test_ladrunoBrick_asdconcrete's mesh-objectivity ratio to
  7.03 — repro'd locally by test-order sequencing, fixed with file-unique tag
  337218 + LEDGER_quirks entry; rule: file-unique ASD tags in every test
  file); (r3) merged #562/#563 cleanly (pinned a betaK-Rayleigh clobber gate
  into the ADR P3 row — Brick20 should be structurally immune post-F12/F13),
  all green → squash-merged (repo does NOT allow gh auto-merge). Then
  guppi/gmsh-hex20-fix merged onto new ladruno (clean), pushed, **PR #564**
  opened (delta: GmshRecorder.{h,cpp} + ledger + round-trip test), CI watch
  armed. NEXT: #564 green → squash-merge → phase pause; P2 uri is a fresh
  phase (R3-deferred debts pinned in ADR §6 P2 row).

- **2026-07-10, iteration 9b: Gmsh hex20 side-fix COMPLETE (user-fired chip).**
  Isolated-worktree agent, commit `ac63a879a` on LOCAL branch
  `guppi/gmsh-hex20-fix` (based on origin/ladruno b530f0854, NOT pushed):
  GmshRecorder 20-node hex → MSH type 17 (`GMSH_HEXAHEDRON_20`) + write-order
  permutation [8,11,16,9,17,10,18,19,12,15,13,14] (independently re-derived
  from shp3dv L/M/N + the Gmsh reference manual; corners pass-through), for
  BOTH `Twenty_Node_Brick` and 33018; round-trip pytest 2/2 green + real-gmsh
  oracle (detJ ≡ 1/8, volume 1.0); ledger rows incl. upstreamable table.
  Siblings TwentyNodeBrick / _u_p_U / TotalLagrangianFD20NodeBrick still carry
  the type-12 bug (documented, out of scope). SEQUENCING: after PR-2 merges →
  rebase (ledger-row conflict expected, keep-both), push, open PR-3.

- **2026-07-10, iteration 9: R3 review DONE → fix round in flight.**
  R3 = 8-angle review (line-scan / copy-strip audit / cross-file tracer /
  reuse / simplification / efficiency / altitude / conventions), all reported.
  Implementer deviations adjudicated: dual-scalar U1 probe ACCEPTED (ADR §7
  noted incl. lazy-allocation limitation), loud getMass fallback ACCEPTED,
  printA-based gates ACCEPTED, guide WRITTEN by orchestrator. Verified
  findings → fix-round agent (Opus, running) applying F1–F17: headline **F1
  per-GP detJ≤0 guard** (dropped vs the upstream anchor — silent dv sign-flip
  on curved second-order meshes), F2 factory 3D-capability probe (null-clone
  segfault), F3 advisories per-element→process-static (O(N) opserr flood),
  F4 parse-time `-lumped` notice, F5 argv guard, F6 trailing `-damp` error,
  F7 updateParameter aggregation, F8 URI wire-coercion guard, F9 setDomain
  ordering, **F10 GmshRecorder 33018 line REMOVED** (MSH type 12 = 27-node
  hex per GmshRecorder.h's own comment + gmsh hex20 node-order mismatch →
  H20 gmsh output deferred; task chip spawned for the identical upstream
  `Twenty_Node_Brick` bug), F11 per-instance geometry cache, F12 M0 mass
  cache + momentum-pass separation, F13 residual-path gating, F14 massless
  early-out, F15 getInitialStiff dedup, F16 dead-code sweep, F17 marker.
  Deferred-with-teeth: ADR §6 P2 row now carries the basisInfo recorder seam,
  blocked BᵀDB, and battery-oracle-import debts; guide rows added (gmsh,
  embedment). NEXT: fixer report → verify gates (reduce-to must stay ≤1e-12)
  → PR-2 assembly (commit all incl. U1 docs, push, `gh pr create --base
  ladruno`, CI watch, stranded-commit discipline).

- **2026-07-10, iteration 8: P1 IMPLEMENTED — all gates green, awaiting R3.**
  Implementer delivered tasks 1.1–1.4 + 1.6 on branch
  `guppi/xenodochial-lamarr-adr72-p1` (uncommitted — R3 review first). New:
  `LadrunoBrick20.{h,cpp}` (std-only, 27-pt, P0-kernel-consuming; URI ordinal 1
  reserved+parser-rejected; U1 advisory at setDomain), `OPS_LadrunoBrick20.cpp`
  (-hourglass hard error; -lumped parsed, getMass errors "lands in P3" + falls
  back consistent), registrations (broker + functionMap + classic-Tcl table +
  CMakeLists + 4 viz vtktypes + Ladruno_ElementResults 33018 in the 20N/GL3
  dispatch), `tests/test_ladrunoBrick20_element.py` (18/18 PASS). **Reduce-to
  residuals (rel-max): disp 6.8e-15, force 5.5e-15, per-GP stress 3.2e-15,
  K 8.1e-16, Newmark K+c3M 8.0e-16** — three orders below the ~1e-12 gate.
  Rank: exactly 6 zero modes, λ6/λ7 separation 6e13. Corner-share sign
  pathology asserted (−V/8 corners / +V/6 mids). Classic-Tcl smoke PASS
  (OpenSees.exe). P0 kernel harness still 8/8 (local box needed
  `pip install sympy`); sibling LadrunoBrick battery 20/20;
  check_manifest/check_classtags OK. Bookkeeping: banner line + patch_banner,
  LEDGER_implementations row, 8 LEDGER_vanilla rows (_PR pending_), manifest
  row planned→active + dispatch{tcl,functionMap}=true + LadrunoQuad-style tcl
  WAIVED wording; headers stamped. **DEVIATION for R3:** the U1 probe
  discriminates on the damage-channel SHAPE (advise only for dual-scalar
  size≥2 = ASDConcrete3D {d+,d−} / LadrunoConcrete3D {ω_t,ω_c}); plain
  LadrunoJ2 exposes a size-1 Lemaitre "damage" response even with the law OFF,
  so a bare present/absent probe would violate the "never for LadrunoJ2" gate —
  LadrunoJ2 itself untouched. NEXT: R3 (/code-review medium + vanilla-footprint
  checklist), then commit + PR off ladruno.

- **2026-07-10, iteration 7: P0 MERGED (#548, squash 564eda1a6) → P1 launched.**
  CI round 2 all green (Zone-A 14m incl. the hex20 g++ cross-check on Linux;
  manifest gate 22s). Squash-merged; branch `guppi/xenodochial-lamarr-adr72-p1`
  cut from updated ladruno (carries the U1 doc edits, uncommitted). **P1
  implementer agent RUNNING** (tasks 1.1-1.4 + 1.6: element + factory +
  registrations + battery + bookkeeping incl. the manifest-row update and the
  U1 one-time advisory; build+test loop on the Windows worktree). NEXT on its
  report: R3 = /code-review (medium) + vanilla-footprint checklist, then PR-2.
  NB PR #549 (ADR-70 docs) also landed on ladruno meanwhile — no overlap.

- **2026-07-10, iteration 6: PR #548 opened; first CI round FAILED → fixed.**
  Two failures, both diagnosed + fixed in 42da0bcf0: (1) **G9 manifest gate**
  — every Ladruno classTag needs a `testbed/manifest.yaml` row (NEW standing
  constraint: add the manifest row in the SAME commit as any classTags.h tag;
  §0 should be read as including this) → added LadrunoBrick20 row
  (status=planned, pytest=test_hex20_kernel_cpp.py, tcl WAIVED until P1);
  (2) **Zone-A collection ImportError** — CI env had no `sympy` → added to
  the Zone-A pip line (ladruno.yml) + `pytest.importorskip("sympy")` guard in
  the harness. Local gates re-verified: check_manifest OK, check_classtags
  OK, harness 8/8. U1 RESOLVED by user: **BOTH** (guide caveat + one-time
  opserr advisory; folded into P1 tasks 1.1/1.4/1.6; ADR §7 updated,
  uncommitted — rides with the P1 branch). CI re-watch armed. NEXT: on green
  + merge, cut the P1 branch from updated ladruno and spawn the P1
  implementer.

- **2026-07-10, iteration 5: R2 PASS → PR-1 assembled.** R2 audit verdict:
  PASS on all five items — formula independence (mutation-tested: an injected
  corner-derivative transcription bug misses by ~10¹⁰× the tolerance), table
  provenance vs primary sources (20/20 nodes vs L/M/N; GP27/GP8 vs
  brcshl/recorder; Irons-14 verified by moment conditions), full gate coverage
  (every ADR §6 P0 gate = a raising assertion; tolerances judged sound: K rtol
  1e-11 is the right floor for -O2 vs LAPACK, the P1 reduce-to is the strict
  1e-12 gate), the G3 fix correct-and-non-vacuous for the committed seed, and
  all duplicated driver constants have live drift-failure paths. Applied R2
  **Patch 1** (deterministic top-degree even distinguisher in G3 — kills the
  lucky-seed vacuity risk); skipped optional Patch 2 (low value, noted).
  Re-ran after the patch: oracle 7/7 + cross-check 8/8 PASS. **P0 review
  gates R1 + R2 both closed. PR-1 committed/pushed from this branch → CI
  watch.**

- **2026-07-10, iteration 4: P0 GATES ALL GREEN; R2 in flight; PR-1 staged.**
  P0 finisher reported: driver `tests/hex20_kernel_check.cpp` + harness
  `tests/test_hex20_kernel_cpp.py` written; oracle G1–G7 PASS; C++ cross-check
  8/8 PASS (~1e-11/1e-12); NO header algebra bugs; header stamped (dir glob
  already covered it); ONE oracle fix (G3 probe: per-variable degree is wrong
  for the non-product Irons rule → total-degree probe; rule moments
  hand-verified). **14-pt verdict: rank 54 CONFIRMED — ADR §7 `i14` note
  stands.** Key pinned numbers: ranks 54/48/54, 6 spurious modes catalogued,
  row-sum corner −M/8, HRZ 7/248 & 2/31 positive+conserving, VOL27==VOL8.
  Orchestrator: LEDGER_vanilla_files row added for the classTags.h edit
  (_PR pending_); regenerable `hex20_spurious_modes.txt` removed from the
  tree. R2 (oracle-independence + gate-coverage audit) RUNNING. NEXT on R2
  PASS: commit changeset (ADR + plan + kernel + oracle + driver + harness +
  classTags + README + ledger) on this branch, push, `gh pr create` →
  `ladruno` (body via --body-file -), fix the ledger PR link, tick the P0
  checkbox, WATCH `gh pr checks`. On R2 FAIL: apply its patches first.

- **2026-07-10, iteration 3: R1 SIGNED OFF + tag reserved.** R1 blind verifier
  reported **FULLY CONSISTENT** — its independent derivation (from
  Twenty_Node_Brick.cpp L1847-1873, shp3dv.cpp shap3dv L19-215 + brcshl
  L219-393, Brick.cpp L124-125/L536-542) matches the kernel's ordering memo
  row-for-row: all 20 nodes, all 27 GPs + weights (b=√(3/5), {5/9,8/9}³), 8-pt
  ζ-fastest lexicographic, exact Abaqus C3D20 convention, recorder tables
  agree at every point. Orchestrator diffed R1's tables against
  LadrunoHex20Shape.h NODE_XI/GP27/GP8: **zero mismatches → R1 gate PASS.**
  Task 0.7 DONE: `ELE_TAG_LadrunoBrick20 = 33018` reserved in classTags.h
  (+ a 33017 LadrunoUP claim note, + 33019 pencilled for H27).
  Still running: the P0 finisher agent (driver + pytest harness + gates to
  green + stamp + 14-pt verdict). NEXT on its completion: spawn R2
  (oracle-independence audit), then PR-1 assembly (ledger row for classTags.h,
  checkbox, gh pr create off a fresh branch → ladruno, CI watch).

- **2026-07-10, iteration 2 (P0 IN FLIGHT).** Classifier recovered. During the
  outage the orchestrator STAGED (written, unrun): the kernel
  `SRC/element/ladrunoBrick/LadrunoHex20Shape.h` (ordering memo + tables +
  serendipity N/dN + J/B/M/rowsum/dofDirs/elastic-K) and the oracle
  `tests/hex20_reference.py` (sympy-derived, gates G1–G7 incl. exact −1/8 &
  7/248 / 2/31 fractions, Irons 14-pt self-validating rule). Now RUNNING in
  background: (a) **P0 finisher agent** — writes `tests/hex20_kernel_check.cpp`
  + `tests/test_hex20_kernel_cpp.py` (mirror fs2d harness: g++ -O2 -std=c++17,
  -I ladrunoBrick + analysis/integrator, -DLADRUNO_MASSLUMPING_STANDALONE),
  runs oracle + cross-check to green, stamps header, reports 14-pt verdict;
  (b) **R1 blind verifier** — independent ordering derivation from
  shp3dv/brcshl + recorder, no access to the implementation. NEXT on their
  completion: orchestrator diffs R1 memo vs kernel memo (R1 sign-off), spawns
  **R2 oracle-independence audit**, then task 0.7 tag reservation edit
  (33018) + PR-1 assembly. NOTE for R2: kernel AND oracle were both authored
  by the orchestrator (outage fallback) — the audit checks the oracle's
  sympy-derivation independence and gate coverage as specified.

- **2026-07-10, iteration 1 (INCOMPLETE — infra outage).** /loop started. Done:
  task 0.7 frontier re-check **PASSED** (ELE 33018/33019 free; the classTags
  33018/33019 hits are ND/LADRUNO registries — per-registry bands, no
  collision). Orchestrator adjudication notes for the R1 diff captured:
  `shp3dv.cpp` header documents the node pattern (corners 1–4 lower CCW from
  (−1,−1,−1) / 5–8 upper; mid-edges 9–12 bottom ring 1-2,2-3,3-4,4-1; 13–16 top
  ring; 17–20 vertical 1-5,2-6,3-7,4-8 — the Abaqus C3D20 convention), L/M/N
  tables at `shp3dv.cpp:75-79`, variable-node correction at `:129-211`;
  `brcshl` GP arrays RA/SA/TA at `:251-277` (corners → edge-mids → faces
  +r,+s,+t,−r,−s,−t → centroid; weights {5/9,8/9}³) — MATCHES the recorder
  `Hexahedron_GaussLegendre_3` table. BLOCKED: the P0 implementer + R1
  verifier agent spawns (and ScheduleWakeup) failed repeatedly on a
  `claude-opus-4-8[1m]` safety-classifier outage — no state-changing tool
  launches available. **Next action on resume: spawn the P0 implementer
  (tasks 0.1/0.3/0.4/0.6) and the R1 blind verifier (prompts as specified in
  P0 rows above + briefing packet), then R2 after the implementer completes.**

## Subagent briefing packet (attach to every Opus prompt)

1. The ADR: `Ladruno_implementation/72_ladruno_second_order_brick_adr.md`
   (decisions + §6 gates are non-negotiable).
2. This plan file (the task being executed + §0 standing constraints).
3. Templates/anchors: `SRC/element/ladrunoBrick/LadrunoBrick.{h,cpp}` +
   `OPS_LadrunoBrick.cpp` (the skeleton), `SRC/element/brick/Twenty_Node_Brick.cpp`
   (the reduce-to anchor + brcshl source), `SRC/recorder/Ladruno_ElementResults.h`
   (rules/orderings), `SRC/analysis/integrator/LadrunoMassLumping.h` (HRZ),
   `tests/test_ladrunoBrick_element.py` (test harness patterns),
   ADR-70 P0 commit `6347ab76a` (oracle layout).
4. The quirks that bite: build via `cmd /c build.bat`; `python -S` test
   runner; number-aware parsing; `getCopy("ThreeDimensional")`; append-only
   serialized ordinals; stamp headers; ledgers in-PR.
