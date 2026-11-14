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

## 232f — phase 1: implicit dynamic and quasi-static through mode A

Kickoff 2026-11-14 (orchestrator, Opus 5.5). Branch `wp/232f-mpm-implicit-mode-a`, one draft PR
for the WP. Two workers, run as the 232b / 232c pair was: the verification worker never reads the
mode-A C++, and the C++ worker never reads the verification worker's reference scripts — only its
cases file.

**Model and effort, re-decided at kickoff.** C++ worker: opus · high (the roster's fork C++ row;
the WP is new element and driver numerics, not a re-plumb). Verification worker: opus · high, a
different agent (the roster's oracle row: its references are the WP's warrant). Ready-flip critic:
fable · high (the authors are Opus). `/code-review high`. No out-of-family verdict is pre-booked:
232c's was waived by the owner (#1075, 2026-11-14); ask the owner at ready-flip.

**Pre-registered gates** (fixed here, before any mode-A number exists; a band is never widened
after a run — a miss is reported as a miss, with its diagnosis):

| Row (ADR §8, P1) | Gate |
|---|---|
| large-deflection cantilever, quasi-static (gravity ramp) | tip deflection within 5 % of the reference at the finest of three resolutions; error non-increasing across them |
| large-deflection cantilever, transient (Newmark) | first-peak tip deflection and first period within 5 % of the transient reference |
| non-associated mild-cone footing (φ_ps = 27.47°, plane strain, exact q₀·N_q = 138.9 kPa) | the peak / plateau at the finest of three resolutions (cells across B × particles per cell) within [0.90, 1.10] × 138.9 kPa; all three reported. FE's own monotone sequence 1.0020 / 0.9517 / 0.9260 (ADR-79) is context, not the gate |
| elastic column under self-weight, implicit quasi-static | the slice-1 small-compression case's closed form σ = ρ₀g(l₀ − Y) at its frozen tolerance (`tests/mpm_cases.txt`, Amendment 4; the 20 % case stays an XFAIL — quirk `WP-232b-ugimp-fixed-domains-…`) |
| plate with a hole, stated d/W | peak σ_θθ / σ_net within 5 % of Howland's factor at that d/W |
| 1D vibrating bar, implicit (Newmark average acceleration, small Δt) | frequency within 1 % of closed form — the transient re-seed (W1) leg |
| failed-step retry | a step that fails before commit is cut and retried; G2P count = accepted steps; final kinematics equal (rel 1e-12) to a run that used the cut Δt from the start |
| no-cell-crossing re-seed | a step in which no particle changes cell equals the re-seeded reference (rel 1e-12); a debug switch that skips the unconditional `domainChange()` on such a step gives a measurably different (wrong) answer |
| recorder synchrony | particle-group positions and stresses of every output come from the same accepted step |
| restart (database) | save at step n, restore, continue to N: particle state equal (rel 1e-10) to the uninterrupted run, in both shipped engines |
| refusals | non-whitelisted integrator, `LadrunoHHT`, `LadrunoGeneralizedAlpha`, `VariableTransient`, arc-length, displacement control on a grid node — each loud, none a wrong answer |

### 232f-V — verification references (opus · high; never reads the mode-A C++)

```text
<header>  ADR sections: §5.3, §5.4, §8 (the P1 rows), §9 phase 1. Your PR: #<pr>.
Branch: work detached on origin/wp/232f-mpm-implicit-mode-a; commit; `git pull --rebase origin
wp/232f-mpm-implicit-mode-a`; `git push origin HEAD:wp/232f-mpm-implicit-mode-a`. Another worker
pushes C++ to the same branch; you touch only the files below.
Independence: read nothing under SRC/domain/pattern/ladrunoMPM/ and no new mode-A element source.
You MAY read and run the fork's FE elements and the shipped explicit `mpm` engine (its guide:
Ladruno_implementation/232_ladruno_mpm_guide.md) — they are different code paths.
Goal: the expected numbers that gate 232f, each from a reference that is not self-referential
(closed form, an independent numpy solve, or cross-code / a different fork code path), with the
model-form gap between reference and MPM stated.
Deliver:
  1. tests/mpm_implicit_cases.txt — machine-readable, in the style of tests/mpm_cases.txt: per
     case the geometry, material, loading, grid (h, cells, basis), particles per cell, step count
     and Δt, the observable and the expected value, and the pre-registered tolerance copied from
     the gate table in .claude/skills/ladruno-mpm/references/wp-briefs.md ("232f"). Never widen
     a band; if one is physically unreachable, say why in your report — the orchestrator decides.
     Cases: cantilever QS and transient (three resolutions), the mild-cone footing (three
     resolutions; rough rigid footing as a prescribed-velocity grid condition — not a
     DisplacementControl integrator on a grid node, which ADR §5.3 forbids), the self-weight
     column (implicit QS leg of the slice-1 small-compression case), plate with a hole (state
     d/W; the particle set written as a Contract-1 file via tests/mpm_contracts.py), and the 1D
     vibrating bar under Newmark average acceleration.
  2. tests/mpm_implicit_reference.py — the scripts that produced every expected number: e.g. an
     elastica / beam ODE solve in numpy for the cantilever plus a fork FE run (2D plane strain or
     a 3D slab; say which element and why its kinematics are comparable to the hypo kernel at
     these strains), Prandtl–Reissner N_q for the footing with the ADR-79 FE sequence as context,
     Howland's factor at your d/W (cite the table you used), the bar's closed-form frequency.
  3. tests/test_mpm_implicit_reference.py — checks the references themselves (closed forms
     reproduced, ODE converged, FE run converged in its own mesh) — never against MPM.
Every FE number names the build hash (`ladrunoBuild`), the deck, the output path and the date.
python3.12; ruff F811 clean; the fork-checklist gates before every push.
Done when: the cases file covers every row above with a reference and a tolerance; your reference
test is green; reproducible from the cases file alone. Report per the Worker rules.
```

### 232f-C — mode A in C++ (opus · high; a different agent from 232f-V)

```text
<header>  ADR sections: §4 (assets, W1–W17), §5.1, §5.3, §5.4, §6, §8 (P1 rows), §9 phase 1.
Also follow .claude/skills/ladruno-new-element/SKILL.md for the particle element. Your PR: #<pr>.
Branch: as 232f-V (detached on origin/wp/232f-mpm-implicit-mode-a, rebase, push HEAD to it).
Read tests/mpm_implicit_cases.txt when it lands; never tests/mpm_implicit_reference.py.

ROUND 1 — stop and report after it (the interface freeze; the orchestrator reviews it):
  tests/mpm_implicit_interface.md, at most ~3 pages, deciding with code anchors quoted:
  - the particle element (provisional LadrunoMPMPoint): state, nodes = the support of its cell
    (GIMP support spans neighbours), in-place re-pointing and how it refetches Node*;
    getResistingForce / getTangentStiff with the hypo-kernel entry (frozen positions, F_trial
    from grid increments, the geometric term), mass model (one, per the element guide), Rayleigh;
  - grid nodes as Domain nodes: lifetime, tag scheme (never reuse a retired tag, W5), and the
    empty-support treatment for phase 1 (excluded or ghost-regularized: pick one, keep it
    switchable for M0b, and say what it changes in the answer);
  - the per-step sequence: P2G of committed disp (zero), vel, accel onto Domain nodes; the
    unconditional C++ Domain::domainChange() (never the user command, W4); analyzeStep(dt)
    (W10) / StaticAnalysis::analyze(1,false) + flushRecorders (W16); G2P in commitState (W11,
    must not fail); the driver's cut-and-retry with nothing to undo; W17 abort;
  - integrator whitelist by getClassTag (Newmark, HHT, GeneralizedAlpha); VariableTransient by
    dynamic_cast; LadrunoHHT / LadrunoGeneralizedAlpha refused; two bridges for classic Tcl and
    Python as LadrunoStaggeredAnalyze has;
  - the command surface added to `mpm` (mode selection, prescribed-velocity grid conditions with
    a reaction observable for the footing, `mpm analyze` in mode A, any `info` counters the tests
    need — e.g. a G2P count — and the debug switch for the re-seed test);
  - restart: what the element and the host serialize, broker entries, and lifting the slice-1
    database refusal in mode A only;
  - mutation: the MPM family's ZERO / SCALE hooks reaching the mode-A internal force and G2P.

ROUND 2 — on the orchestrator's "frozen":
  Implement, batching edits; ask the orchestrator for each build (all five targets) and each
  mutant build — never start build.bat yourself. Tests tests/test_ladruno_mpm_implicit_*.py
  (zone_a), classic Tcl and Python, reading tests/mpm_implicit_cases.txt; GCC-clean code.
  Registrations per fork-checklist.md and the element guide (classTag, broker, both element
  dispatch sites if the element is user-creatable, manifest row feature_class element, ledger
  fragments WP-232-mpm-implicit-mode-a.md, vanilla rows + `// Ladruno` markers, banner, stamps).
  Update the guide Ladruno_implementation/232_ladruno_mpm_guide.md with mode A.
Done when: every row of the 232f gate table green (or a reported miss with diagnosis, never a
widened band); MPM ZERO and SCALE mutants red on a mode-A oracle; ci/ gates clean; Zone-A green on
a dispatched run; ladrunoBuild matches your head on the exe and the pyd. Report per the Worker
rules.
```
