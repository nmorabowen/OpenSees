# Dispatch briefs — one per work package

The orchestrator copies a brief, fills in the `<…>` fields, and dispatches it with the roster's
`model` and `effort` and `isolation: "worktree"`. Briefs for apeGmsh work (232d, 232e, R's pin
bump) are handed to apeGmsh's T coordinator instead (SKILL.md, "apeGmsh — only through its T
coordinator"). Until WP-232 merges, the skill, the ADR and the
plan exist only on `origin/wp/232-mpm-general-driver`: workers read them with
`git show <ref>:<path>`, where `<ref>` is that branch or, after the merge, `origin/ladruno`.

Every brief starts with the same header:

```text
You are a worker on the Ladruno MPM program (WP-232). Read, in order, from <ref>:
  .claude/skills/ladruno-mpm/SKILL.md            (the "Worker rules" section — binding)
  .claude/skills/ladruno-mpm/references/architecture.md
  .claude/skills/ladruno-mpm/references/fork-checklist.md   (fork WPs)
  Ladruno_implementation/232_ladruno_mpm_general_driver_adr.md  (§ named below)
Program issue: #<program>. Your PR: #<pr or "open it on your first push">.
Last handoff's Next: <paste>.
```

---

## 232a — the contracts and the apeGmsh ADR (opus · high)

```text
<header>  ADR sections: §7 (all), §4 rows "LadrunoRecorder and the .ladruno format", W13.
Goal: the two contracts as field-by-field schemas, fixture files each way, and the apeGmsh ADR
that admits them — so 232c (fork writer) and 232d (apeGmsh readers) can build in parallel.

Fork half — branch wp/232a-mpm-contracts from origin/ladruno:
  1. Ladruno_implementation/232a_mpm_contracts.md — Contract 1 (particle set) and Contract 2
     (the .ladruno particle group — an additive change: FORMAT_MINOR_VERSION 1 → 2,
     FORMAT_VERSION stays 1, because apeGmsh refuses any other major): every field with name, dtype, shape,
     units, required/optional, default, validation rule; versioning; id rules (stable, never
     reused, removal and addition); output cadence and field selection; precision (float64
     positions); chunking and compression; grid metadata. Read the format docs first:
     Ladruno_implementation/03_ladruno_recorder.md, 164_ladruno_recorder_write_path.md,
     165_ladruno_recorder_stage_semantics.md.
  2. tests/mpm_contract_fixtures.py (python3.12, h5py + numpy) writing two tiny fixtures:
     tests/data/mpm/particle_set_v1.h5 and tests/data/mpm/particles_v1m2.ladruno (a minimal
     .ladruno shell with INFO/FORMAT_VERSION = 1, FORMAT_MINOR_VERSION = 2 and the particle
     group; a few particles, a few outputs, one removal). Cover the fork's own format tooling:
     Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py and make_synthetic.py.
  3. tests/test_mpm_contracts.py — a schema validator that accepts both fixtures and refuses
     each malformed variant (missing field, negative volume, reused id, wrong dtype).

apeGmsh half — a DRAFT ONLY, in the fork branch. You cannot write in the apeGmsh checkout; the
orchestrator hands the draft to apeGmsh's T coordinator, which lands it:
  4. Ladruno_implementation/232a_apegmsh_adr_draft.md — an apeGmsh ADR in apeGmsh's own format
     (read, read-only, from apeGmsh origin/main: .claude/skills/apegmsh-adr-docs/SKILL.md and two
     recent ADRs under architecture/decisions/), Status: Proposed, its number left as <next>.
     It decides: the particle-set placement (a FEMData stream in model.h5 with the
     compose-streams obligations, or a sidecar as an explicit exception to ADR 0112 D1); that
     results are the .ladruno particle group (ADRs 0122, 0131); the particle-cloud result
     ontology; the three.js points layer (ADRs 0112 D8, 0116); emit-time tag resolution
     (charter P4); deck-only staged runs (ADR 0121). Only apeGmsh's maintainer may accept it.

Done when: both schemas complete; fixtures generated and validated (pytest green); the apeGmsh ADR
draft complete; ledger gate clean. Report per the Worker rules.
```

## 232b — the numpy oracle (opus · high)

```text
<header>  ADR sections: §3, §5.2, §6, §8. Branch wp/232b-mpm-oracle from origin/ladruno.
Goal: an independent reference implementation of the explicit MPM step and of the incremental
hypo-kernel entry, and a frozen interface the C++ worker codes against.
Independence: read NOTHING under SRC/. Use the ADR, the literature, and numpy (python3.12).
Conventions: the hypo kernel's Voigt order {00, 11, 22, 01, 12, 20}, engineering shear in strain —
from tests/hypo_reference.py (test-side, allowed) and ADR-79 §3.

Push FIRST (the interface freeze — 232c waits on it):
  1. tests/mpm_interface.md — particle state layout; grid indexing; GIMP weights and gradients
     (1D/2D/3D); P2G and G2P formulas; the MUSL order (ADR §5.2 steps 3–7); the PIC/FLIP blend;
     the mass cutoff; local non-viscous damping; the critical-Δt estimate; the incremental kernel
     entry signature (Fn, F1) → (Δε in the unrotated frame, R1, J1); sign and Voigt conventions.
  2. tests/mpm_cases.txt — machine-readable cases with expected numbers (inputs, step count,
     tolerances), and tests/hypo_incremental_cases.txt for the kernel entry.
Then:
  3. tests/mpm_reference.py and tests/hypo_incremental_reference.py — the oracles.
  4. tests/test_mpm_reference.py — the oracle checked against closed forms: 1D vibrating bar
     (frequency; PIC vs FLIP dissipation), elastic self-weight column with local damping
     (σ = ρ0·g·(l0 − Y), Charlton et al. 2017), rigid translation and rotation of a stressed body
     (zero spurious stress), two-body elastic impact (momentum, energy), block on a frictional
     incline, critical Δt (stable just below the estimate), the cell-crossing meter (a number).

Done when: the interface and cases are pushed and the orchestrator has declared them frozen;
`python3.12 -m pytest tests/test_mpm_reference.py` green; every case reproducible from the
cases files alone. Report per the Worker rules.
```

## 232c — fork slice 1 (opus · high; never the agent that wrote 232b)

```text
<header>  ADR sections: §4, §5.1, §5.2, §5.4, §6, §8 (S1 rows), §9.
Branch wp/232c-mpm-explicit-slice from origin/ladruno (after 232b's freeze; rebase after 232b
merges). Read tests/mpm_interface.md and tests/mpm_cases.txt from 232b's branch; never
tests/mpm_reference.py.
Goal: slice 1 of the fork — the explicit engine, hosted as a LoadPattern, end to end from a
Contract-1 file to a .ladruno particle group, in all three interpreter engines.

Build (see fork-checklist.md for every registration and gate):
  - SRC/domain/pattern/ladrunoMPM/ — the host (LadrunoMPMPattern), the engine (grid, active
    nodes, GIMP, P2G/G2P, MUSL, PIC/FLIP, mass cutoff, local damping, critical Δt, grid planes,
    frictional wall), the particle store and its Contract-1 reader; wired like
    ladrunoPorousOverlay (its CMakeLists.txt and the add_subdirectory line).
  - The incremental entry in SRC/element/solidTransformation/LadrunoHypoKernel.h, plus
    tests/hypo_incremental_kernel_check.cpp against tests/hypo_incremental_cases.txt (the
    g++ header-check pattern of tests/hypo_kernel_check.cpp).
  - The `mpm` command: one LadrunoCommandTable.h row; no analysis bridge in slice 1 (mode F has no
    analysis); refused in OpenSeesSP / MP.
  - The .ladruno particle group writer (FORMAT_MINOR_VERSION 2; FORMAT_VERSION stays 1), matching
    232a's fixture layout, with ladruno_format.py and make_synthetic.py updated.
  - The stage path: `parameter <t> loadPattern <mpmTag> …` (classic Tcl's `pattern` keyword is
    only `pattern <tag> lambda`).
  - Admitted materials and loud refusals (architecture.md, "Materials"); a commit refusal aborts
    the run (W17).
  - Tests: tests/test_ladruno_mpm_*.py for every S1 row of ADR §8, marked zone_a, run in classic Tcl
    and Python; h5py added to Zone-A's pip line if a test needs it; code that builds with GCC.
  - The mutation gate: the MPM family through the existing mechanism (fork-checklist.md, "Tests").
  - Warrant package: manifest row (feature_class: loadpattern), implementations fragment,
    vanilla fragment for every upstream file touched, quirks learned, banner line, header
    stamps, the guide Ladruno_implementation/232_ladruno_mpm_guide.md, and the ADR-87 D4 handoff
    pack for the owner's out-of-family verdict.
Batch edits; ask the orchestrator for each build (all five targets) and each mutant build.

Done when: every S1 test green against 232b's cases in classic Tcl and Python; the MPM family's
ZERO and SCALE mutants red; all ci/ gates clean; Zone-A green on a dispatched run; ladrunoBuild
matches your head on both the exe and the pyd. Report per the Worker rules.
```

## 232d — apeGmsh seeding, readers and the points layer (opus · high; sonnet · medium for seeding plumbing)

```text
<header, apeGmsh repo>  ADR sections: §7. Starts only after the owner accepted the apeGmsh ADR.
Work in nmorabowen/apeGmsh from origin/main under its own rules; load its skills first:
`apegmsh` (usage), `apegmsh-viewer-results`, `apegmsh-viewers-change`,
`apegmsh-viewers-visual-check`, and `apegmsh-adr-docs` if you touch docs.
Build against 232a's fixtures (tests/data/mpm/ in the fork), not against a fork build.
  1. Lattice seeding: n particles per cell per axis inside each body (Gmsh inside test or the
     body's mesh); positive volumes; GIMP half-lengths cell/(2n); a particles-per-cell
     uniformity report; plane-strain thickness. Write Contract 1 at emit, with material tags
     resolved from labels by the emitter. No Gauss-point seeding in this WP.
  2. The .ladruno particle reader in src/apeGmsh/results/readers/_ladruno.py and a TypeScript
     reader beside apeGmshViewer/src/reader/ladruno*.ts.
  3. A THREE.Points layer with per-output position buffers, coloured by stress or state, after
     apeGmshViewer/src/renderer/nodes.ts (makeNodePoints); the visual check with evidence.

Done when: the fixture .ladruno loads and animates in apeGmshViewer; seeding writes a valid
Contract-1 file and its uniformity report; apeGmsh CI green; the visual check posted; the card on
board #1203 updated. Report per the Worker rules.
```

## R — the fork release and the pin (the owner, with the orchestrator)

Not a worker brief. After 232c merges: follow the release recipe in `Ladruno_scripts/README.md`
(build and verify, then tag; the workflow makes the release; publish the installer). Then an
apeGmsh PR bumps `FORK_PIN` — `tag`, `sha`, `repo` and `asset`, all four together — and apeGmsh's
`live-fork` lane must go green on the new pin.

## 232e — the apeGmsh bridge (sonnet · medium)

```text
<header, apeGmsh repo>  ADR sections: §7.1, §7.4(5). Starts only after R.
Load `apegmsh-bridge-feature` and follow its checklist: typed primitives for the `mpm`
commands; the `ladruno_fork` marker on every test; a row in
docs/concepts/backend-capabilities.md; capability probes that return None when inconclusive;
staged MPM runs deck-only (live runs are flat only); the particle file path brace-quoted in
emitted Tcl; the emitter writes Contract 1's material column (charter P4).

Done when: apeGmsh CI and the live-fork lane green on the pinned release; a deck emitted by
apeGmsh runs a slice-1 case end to end and its .ladruno opens in apeGmshViewer. Report per the
Worker rules.
```

## 232f and later

Each phase of ADR §9 is started with `/ladruno-mpm 232f` (and so on). The orchestrator writes its
brief at kickoff, appends it here in the same PR, and re-decides model and effort; phase 1
(implicit through mode A) also follows `.claude/skills/ladruno-new-element/SKILL.md` for the
particle element.
