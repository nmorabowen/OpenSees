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
| non-associated mild-cone footing (φ_ps = 27.47°, plane strain, exact q₀·N_q = 138.9 kPa) | the peak / plateau at the finest of three resolutions (cells across B × particles per cell) within [0.90, 1.10] × 138.9 kPa; all three reported. FE's own monotone sequence 1.0020 / 0.9517 / 0.9260 (ADR-79) is context, not the gate **Outcome (2026-11-15):** unreachable at r3 even by the fork's FE (B-bar 0.882x, locked 1.114x); recorded as a strict xfail, and the owner chose to redesign the gate before phase 4 (ADR 232, amendment after the §8 table). **Redesigned at 232j kickoff (2026-11-17): row L4 of "232j".** |
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

## 232g — phase 2: the M0b measurements

Kickoff 2026-11-15 (orchestrator, Opus 5.5 — the program's orchestrator session). The roster's
"max when declaring M0b" is met by a Fable critic pass at max that checks the application of the
table below before the ADR amendment is written. Branch `wp/232g-mpm-m0b-measurements`, one draft
PR for the WP, one worker. Done when the four measurements are recorded and the mode-C decision is
written into the ADR as an amendment.

**Model and effort, re-decided at kickoff.** Worker: opus · medium. The roster's harness row is
sonnet · medium, but this harness must design 3D mode-A models that converge at 10⁵ particles and
add profiler scopes in `SRC/` (one vanilla file); the rules it serves are fixed here, so medium.
Ready-flip critic: fable · high (the author is Opus), briefed with the data and the pre-registered
table, not with the orchestrator's reading. `/code-review medium`. No out-of-family verdict is
pre-booked (measurement and instrumentation, no new numerics); ask the owner at ready-flip.

**Pinned conditions** (every number carries them): one build, its `ladrunoBuild` printed by the
harness from the binary it runs; machine i7-13700KF (8 P + 8 E cores), 64 GB; Python 3.12;
`MKL_NUM_THREADS` and `OMP_NUM_THREADS` stated per run (implicit: 8; explicit: 1, plus one row at
the default); no build or other heavy job on the machine during a timed run (the orchestrator
launches them); median of 3 repeats with the spread; warm-up steps excluded and stated.

**Definitions.**
- *Rebuild* t_rb: everything a mode-A step spends because the Domain changed — the host's
  preparation (P2G, grid-node churn, re-pointing, SP rebuild, the `domainChange()` raise), the
  analysis's `domainChanged()` (handler, numberer, `setSize` and graph), and every `soe.symbolic`
  that the structure change triggers. Symbolic analysis counts as rebuild: a fixed-structure scheme
  (mode C) would reuse it.
- *Newton iteration* t_it: formTangent + linearSolve numeric (factor + solve) + formUnbalance +
  update + convTest, per iteration, without the symbolic phase.
- R = t_rb / t_it; k = Newton iterations per accepted step; f = t_rb / (wall time of the step).

**Pre-registered decisions** (fixed here, before any M0b number exists; a threshold is never moved
after a run — an outcome that fits no row is reported as such and goes to the owner):

| # | Measurement | Decides | Rule |
|---|---|---|---|
| M1 | t_rb, t_it, R, k, f at N = 10³, 10⁴, 10⁵ particle-elements, 3D; quasi-static (LoadControl) and transient (Newmark); elastic and Drucker–Prager; PARDISO and UmfPack; each empty-node treatment | **whether mode C is ever needed** | Read on the default (PARDISO, `renormalize` 5e-2) at every N, taking the worst of QS / transient and elastic / DP. **C-closed:** f ≤ 0.25 and R ≤ 1 at every N, and R does not grow with N (log-log slope of R against N ≤ 0.1) → C moves from deferred to rejected for implicit, reopened only by a new measurement in phase 11 (parallel). **C-cheaper-first:** f > 0.25, but the largest rebuild sub-phase is avoidable without C (numbering, `setSize` / graph or symbolic analysis on a step whose connectivity did not change — `mpm info connectivityChanges` — or host bookkeeping) → C is not needed; a targeted rebuild optimisation is filed as its own issue. **C-needed:** f > 0.5 at N = 10⁵ and symbolic analysis of a changed graph is ≥ 50 % of t_rb → C is scheduled and the owner places it in §9. **C-deferred:** anything else → C stays deferred, the numbers are recorded, revisit at phase 11. UmfPack decides nothing alone: if it falls in a worse row than PARDISO, the guide recommends PARDISO for mode A. |
| M1-t | the three empty-node treatments (`renormalize` 5e-2 default, `exclude`, `ghost`), same models | whether the default treatment changes | `exclude` anchors free surfaces (232f study): cost only, never a candidate. `ghost` replaces the default only if its answer is within 0.5 % of `renormalize`'s (the case's tip or settlement observable) **and** it is ≥ 25 % cheaper per step at every N → proposed as an amendment for the owner. Otherwise the default stands. |
| M2 | explicit (mode F) throughput in particle-steps per second at N = 10³, 10⁴, 10⁵, 3D, `ElasticIsotropic` → `DruckerPrager` → `LadrunoJ2`, threads 1 (and the default for elastic); the material update's share of the step | baselines; **cannot reopen C** (F owns explicit, ADR §5.5) | SANISAND is not admitted at particles in v1 (`LadrunoMPMPattern.cpp:375-382`), so its row is *projected*: its per-call cost from a non-MPM single-point driver, composed with F's measured non-material cost per particle. The composition is validated on DP, measured both ways: if the DP projection misses the DP measurement by > 25 %, the SANISAND row is reported as unreliable. If F's non-material share exceeds 70 % for elastic at 10⁵ on one thread, a performance issue on F's transfer loops is filed. The SANISAND / DP cost ratio is recorded for phase 5. |
| M3 | memory per particle-element: the slope of private bytes against N (10⁴ → 10⁵) after the first step, mode A and mode F, elastic and DP, with what holds it (element, material copy, `FE_Element` and its matrices, grid share) | the memory side of C | Mode A ≤ 8 kB per particle-element → no memory case for C. > 32 kB → a memory case for C at 10⁶ particles (> 32 GB) — unless one storage item that is ≥ 50 % of it is removable without C (per-element dense matrices that a shared scratch would replace), which is then filed as its own issue instead. Between → recorded, revisit at phase 11. |
| M4 | the 232f data owed: Newton iterations per step on the cantilever QS case (`tests/mpm_implicit_cases.txt`, r1–r3); the tip deflection at `massCutoff` 1e-6 and 1e-10 against 5e-2, under `renormalize` (a divergence and its step are data); formTangent's share of t_it | whether 5e-2 biases the answer; the declined 232f review item 15 | Data, not a gate. If `renormalize` converges at 1e-6 or 1e-10 at all three resolutions and its tip differs from 5e-2's by more than the 232f gate tolerance (5 %), the default is biased → an issue for the owner; otherwise the default stands. If formTangent is > 30 % of t_it at 10⁵, item 15 (the tangent B-build cleanup) is filed as an issue. |

**Then** the orchestrator applies the table, a Fable critic at max checks the application, and the
ADR gets "Amendment (232g, <date>) — M0b result" at §5.5 and an implementation-log line.

```text
<header>  ADR sections: §5.3 (and its 232f amendments), §5.5, §8 (mutation note), §9 phase 2.
Your PR: #<pr>. Branch: work detached on origin/wp/232g-mpm-m0b-measurements; commit; `git pull
--rebase origin wp/232g-mpm-m0b-measurements`; `git push origin HEAD:wp/232g-mpm-m0b-measurements`.
Read the "232g" section of .claude/skills/ladruno-mpm/references/wp-briefs.md first: its
definitions and pre-registered table bind you. You never apply the table: you produce the numbers
and the orchestrator decides.
Goal: a harness that produces every number in that table, reproducibly, from one build.

ROUND 1 (stop and report after it; the orchestrator builds, then reviews a 10^3 smoke run):
  1. Instrumentation, profiler scopes only. No behaviour change: OPS_PROFILE_SCOPE is a no-op
     when the profiler is off; quote SRC/utility/profiler/ProfilerMacros.h:77 in the PR.
     - SRC/analysis/analysis/StaticAnalysis.cpp: a "domainChanged" scope around the call in
       analyze() and dc.* scopes inside domainChanged(), mirroring DirectIntegrationAnalysis.cpp
       (:212, :433-461). A vanilla file: `// Ladruno WP-232g` markers and a row in
       Ladruno_implementation/ledger/vanilla/WP-232g-mpm-m0b-measurements.md.
     - SRC/domain/pattern/ladrunoMPM/LadrunoMPMImplicit.cpp: "mpm.a.prepare" with children for
       P2G, grid-node churn, re-pointing and the SP rebuild; "mpm.a.g2p" in the batched G2P.
     - Mode F (LadrunoMPMPattern.cpp:653-714) has "mpm.advance": add "mpm.f.material" around
       the particle material update so M2's share is measured, not inferred (if that loop is
       threaded, say how you measured the share instead). Check that the profiler's per-thread
       trees neither double-count nor drop threaded scopes.
  2. Ladruno_scripts/mpm_m0b/: the harness (python3.12, opensees.pyd, following
     tests/ladruno_mpm_implicit_harness.py's conventions): m1_rebuild.py, m2_explicit.py,
     m3_memory.py, m4_cantilever.py and run_all.py, each writing JSON under
     Ladruno_scripts/mpm_m0b/results/ with ladrunoBuild, the pinned conditions, the model
     parameters, raw per-step profiler totals, and the derived quantities exactly as the brief
     defines them. Memory through ctypes GetProcessMemoryInfo (psutil is not installed; h5py
     is). Every run in a fresh process.
     Models: 3D cases that converge at 10^3-10^5 particle-elements and cross cells during the
     timed steps (report connectivityChanges per step), e.g. a soft block under a gravity ramp
     or a 3D extrusion of the 232f cantilever, plus a DP variant that yields. M2: a 3D granular
     column collapse (DP yields) and its elastic and J2 twins. The SANISAND per-call driver is
     the same driver for ElasticIsotropic, DP and SANISAND (the composition is validated on DP).
     State every choice and its reason in a README.
  3. `--smoke` runs N = 10^3 with one repeat in under 2 minutes per script. You cannot run it
     before the orchestrator's build.
  Report: worktree path and head sha; files; the cases and why; anything in the table you think
  cannot be measured as defined (say why, and propose the closest measurable quantity: the
  orchestrator decides).

ROUND 2 (after the build and the smoke review): fix what the smoke showed. The orchestrator
launches the full runs (WMI, logged) and gives you the results directory. You write
Ladruno_implementation/232g_mpm_m0b_measurements.md: method, pinned conditions, one table per
measurement with every number, the sub-phase breakdown of t_rb, the scaling fits, and an
"Applying the pre-registered table" section stating which row each measurement falls in,
labelled as your proposal (the orchestrator's application is the decision).
Never start build.bat or a full run yourself. Ledger: an implementations fragment for the
harness (status shipped, no class tags), the vanilla fragment, a quirk entry for anything
learned. python3.12; ruff F811 clean; the fork-checklist gates before every push. Report per the
Worker rules.
```

## 232i — phase 3: the B-spline basis and the manufactured-solution orders

Kickoff 2026-11-16 (orchestrator, Opus 5.5 — the program's orchestrator session, per #1075).
Branch `wp/232i-mpm-bspline-mms`, one draft PR for the WP, one C++ worker. Base `ladruno` at
`6ee968060` (release `ladruno-v20261116`). ADR §8 row "method of manufactured solutions — spatial
convergence order per basis" (P3) and §9 phase 3. Done when the gate table below is green (or a
reported miss with its diagnosis, never a widened band), the warrant package is in the PR, and the
ADR carries "Amendment (232i, <date>)" recording the basis, its boundary treatment and the
contract decision.

**Model and effort, re-decided at kickoff.** C++ worker: opus · high (the roster's fork C++ row: a
new basis in both modes, a boundary treatment, a per-particle body-force path, and its MMS
harness). No separate verification worker: a manufactured solution is its own oracle (a wrong
body force or a wrong reading of the kernel shows up as an order collapse, not as agreement), and
the MMS derivation is checked by an independent finite-difference residual test (below). Ready-flip
critic: fable · high (the author is Opus), briefed with the diff, this section and the ADR rows —
not with the orchestrator's reading. `/code-review high`. No out-of-family verdict is pre-booked
(232c's and 232f's were waived); ask the owner at ready-flip.

**Decisions fixed at kickoff** (the worker may object in its Round-1 report; the orchestrator
decides):

- **K1 — degrees.** Quadratic (`bspline2`) is required; cubic (`bspline3`) is in scope only if it
  falls out of the same code path (a degree parameter; `AxisSupport` already holds 4 nodes per
  axis). Cubic never delays the quadratic gate.
- **K2 — Contract 1 is not changed in 232i.** Its §1.7 makes any new allowed `BASIS` value a
  `FORMAT_VERSION` bump, which is apeGmsh work through T. A B-spline grid is chosen on the deck
  only (`mpm grid … -basis bspline2`); `mpm particles -file` into a B-spline grid is refused loudly
  (the file's `GRID/BASIS` differs from the deck's). Contract 1 v2 is filed as a follow-up for T.
- **K3 — Contract 2** (the `.ladruno` particle group writes `GRID/BASIS`): the orchestrator asks T
  at kickoff whether apeGmsh's merged particle reader refuses an unknown `BASIS` value. The worker
  implements whichever the orchestrator relays before Round 2 (write the new value and amend
  `232a_mpm_contracts.md` §2, or hold the value behind a contract amendment T schedules).
  **Answered by T, 2026-11-16 (apeGmsh main `df1869079`):** write the new value and amend
  Contract 2 §2 in this PR, listing the exact allowed values. Until apeGmsh releases a reader
  that takes the merged contract text, the TypeScript viewer refuses `"bspline2"` by name
  (`ladrunoParticles.ts` `readGrid`, C2-GRID). The Python reader passes the value through and
  nothing reads it. Neither shows wrong values. Holding the value back would make the file lie
  about its grid. Contract 1 v2, the reader change and the bridge's `-basis` value are
  nmorabowen/apeGmsh#3013.
- **K4 — the boundary treatment is the Round-1 decision.** Today a face's `fix`/`wall` zeroes the
  face nodes and the ghost layer alike (interface A5). That keeps linear reproduction for GIMP
  (its support near a face never reaches the ghost layer), but not for a uniform quadratic
  B-spline, whose support does: zeroed ghosts leave an O(h) boundary layer. Candidates: mirrored
  ghost values (odd for a fixed component, even for a free one), a boundary-modified
  (open-knot / clamped) basis at the grid-box faces (Steffen, Kirby & Berzins 2008; Gan et al.
  2018), or the A5 zeroing kept and its order loss reported. The choice must work in mode F and in
  mode A (where a mirrored ghost is an MP constraint with factor −1 under the user's handler — say
  whether that is admissible, or why a modified basis avoids it). Literature only: no GPL / LGPL /
  unlicensed source (Worker rule 5).

**Pre-registered gates** (fixed here, before any B-spline number exists; a band is never widened
after a run — a miss is reported as a miss, with its diagnosis).

The manufactured solution, after Wallstedt & Guilkey (2008): unit square (unit cube for the 3D
row), `u_i(X, t) = A sin(2π X_i) sin(c π t)` for each axis i, rollers on every face (the normal
component fixed; F stays diagonal, so there is no rotation, no shear stress and no shear
traction), `ElasticIsotropic` through the hypo kernel. The body force is
`b = a − (1/ρ₀) ∇_X · P` with the stress the kernel actually integrates (quote the kernel's strain
measure and the stress it returns from the interface documents before deriving; for a coaxial
stretch the Hughes–Winget increments sum to the logarithmic strain up to O(Δt²)). Error
`e = RMS_p |u_p − u(X_p, t_end)| / A` over all particles at `t_end` = one period. Order `p` =
least-squares slope of `log e` against `log h` over the three finest rungs. Ladder
`1/h = 8, 16, 32, 64` (add 128 if a run costs under 10 min); two particles per cell per axis (a
3×3 leg reported). Explicit Δt fixed across the ladder at ≤ 0.1 × the finest rung's critical Δt;
**validity check:** halving Δt at the finest rung changes `e` by < 5 % — if not, Δt is reduced
until it holds (a validity condition on the evidence, not a band). FLIP with α = 1 for the gated
legs (PIC reported).

| Row | Gate |
|---|---|
| G1 basis properties, each new basis, 2D and 3D | at ≥ 10⁴ random admissible positions, including every position within 2h of each face, after the K4 treatment: `S ≥ 0`; partition of unity `\|Σ S − 1\| ≤ 1e-13`; `\|Σ ∇S\| ≤ 1e-11 / h`; linear reproduction `\|Σ S xᵢ − x_p\| ≤ 1e-12 h` and `\|Σ xᵢ ⊗ ∇S − I\| ≤ 1e-11`; the 1D cardinal values in the interior (quadratic: `S(0) = 3/4`, `S(±h/2) = 1/2` and `1/8`; cubic: `S(0) = 2/3`, `S(±h) = 1/6`) to 1e-14; C¹ continuity across knots by a finite-difference sweep of `dS` (max jump ≤ 1e-8 / h). In the g++ core check and through `mpm debug shape`. |
| G2 MMS, mode F, 2D plane strain, small amplitude (max \|F_ii − 1\| = 1e-3) | `bspline2`: p ≥ 1.8 and `e` non-increasing along the whole ladder |
| G3 MMS, mode F, 2D plane strain, large amplitude (max \|F_ii − 1\| = 0.1; particles cross cells at the finest rungs — report crossings per rung) | `bspline2`: p ≥ 1.5 and `e` non-increasing |
| G4 MMS, mode A, quasi-static, 2D, small amplitude (`u(X, λ) = λ A sin(2π X_i)`, λ ramped in ≥ 4 load steps, body force updated per step) | `bspline2`: p ≥ 1.8 and `e` non-increasing |
| G5 cubic, if delivered (K1) | G2–G4 at the same thresholds for `bspline3` |
| G6 mutation | the MPM family's ZERO and SCALE mutants each turn G2 (or G3) red |
| G7 contracts and refusals | `-basis bspline2` round-trips through `mpm info summary`, the `.ladruno` output (per K3) and the database restart in mode A; `mpm particles -file` into a B-spline grid refused (K2); an inadmissible particle refused as today |
| G1-K4 (added at the freeze, before any number) | folded `S′` within 2h of each `fix` face at the G1 tolerances: odd component `Σ S′ (x_ι − x_face)_n = (x_p − x_face)_n`; even component `Σ S′ = 1` and `Σ S′ x_t = x_{p,t}`; `S′ ≥ 0` |
| G8 images (added at the freeze) | mode F, `bspline2`: a body symmetric about an interior grid plane, run on the full grid, against its half on a grid whose face is that plane with a roller (`fix`, images): particle positions and stresses equal to rel 1e-10 after ≥ 100 steps, at small and large amplitude |
| R1 reported, not gated | G2 and G3 for `gimp` and `linear` (literature expectation: GIMP near 2 at small amplitude; linear erratic once particles cross cells) — a GIMP small-amplitude order below 1.5 is filed as an issue against shipped code; PIC legs; the 3×3 leg; the 3D small-amplitude MMS at `1/h = 4, 8, 16`; mode A transient (Newmark average acceleration) MMS; the slice-1 self-weight column at both compressions (the 20 % case is the uGIMP XFAIL) per basis at three resolutions; the cell-crossing meter (S1 row) per basis, as a number |

**Freeze (2026-11-16, orchestrator, on `tests/mpm_bspline_interface.md` @ `b25586c8b`).**
Frozen as B1–B12 with the following conditions:
- K4 = mirrored ghosts by images on `fix` faces of B-spline grids (B5). `wall` stays A5 and is
  first order there, documented.
- Amendments 1 (`mpm info basis`), 2 (G1-K4 above) and 4 (cubic in both modes) are accepted.
- Amendment 3: the gate keeps `e` at one period. `e(T/4)` and `max_t e(t)` are reported, not
  gated.
- **F1 — K3 per T:** write `bspline2` and `bspline3`. Amend C2-GRID in `232a_mpm_contracts.md` §2
  to list exactly `gimp | linear | bspline2 | bspline3`, with no minor bump, and say that a reader
  older than the amendment refuses the value by name (the safe outcome). Split `check_grid`'s value
  set per contract. G7's `.ladruno` leg is a round trip.
- **F2 — memory:** raising `MAX_DOF` must not grow per-element storage for `gimp`, `linear` or
  `bspline2`. Size by the actual support, or keep cubic on its own pool. Report bytes per
  particle-element at 10⁴ (the M0b M3 method) before and after: growth over 5 % for GIMP is a
  defect.
- **F3 — unchanged answers:** B8's per-component path must reproduce GIMP and linear bit for bit.
  The whole existing MPM suite (slice 1, the 232f gates, contracts, lifecycle, the g++ core check)
  passes unchanged, with the counts in the report.
- **F4 — G8**, the images-against-full-domain equality, is required.

**Critic round (2026-11-16, Fable · high, verdict AMEND at `da23a1fed`; each finding verified
in-tree by the orchestrator).** These are added before the final warrant run. No band changes.
- **A1 (major, a defect).** A face-layer node of a `fix` face is its own image, so its mass is
  one-sided. The active-node cutoff (`LadrunoMPMImplicit.cpp:554-561`; `LadrunoMPMCore.h:683-685`
  and `:749-752`) then drops in the half model a node the full model keeps. The critic's probe:
  3-14 % stress difference in 2D and 9.6 % in a 3D octant at mode A's default cutoff; up to 53 %
  in mode F at `massCutoff 0.05`. Fix: on B-spline grids the cutoff tests the face-layer node's
  mass times 2^k, where k is the number of `fix` faces through the node. **G8 gains two legs:**
  G8-A (mode A, Newmark, default `renormalize` 5e-2) and G8-F-cut (mode F at `massCutoff 0.05`),
  each with a body reaching a free surface near the face. The equality bar is unchanged:
  rel 1e-10.
- **A2 (major, the mutation margin).** On the Zone-A rungs SCALE survives the one-period MMS
  tests, so the family's SCALE score (floor 0.38, `mutation_gate.py:134`) has one test of margin.
  Added a Zone-A mutation oracle, which is not an accuracy gate: the quarter-period peak error
  `e(T/4) ≤ 0.05` (5 % of the amplitude) for `bspline2` G2 at 1/h = 32. Re-measure the MPM
  family's ZERO and SCALE scores with `mutation_gate.py` against their floors.
- **A3 (minor).** `verdict()` fails when a gate row lacks its criterion or has `p = None`
  (evidence void, never a pass). The results test recomputes `p` and monotonicity from `e` and
  `h`, fails when G2-G4 are missing, and checks that the rows' `ladrunoBuild` is one build. Local
  orders are reported per row. The guide describes a face that fixes only tangential components as
  an antisymmetry plane (odd tangential, even normal image).

```text
<header>  ADR sections: §4 (walls W1–W17), §5.2, §5.3 (and its 232f amendments), §5.4, §6, §8
(the P3 row, the mutation note), §9 phase 3. Your PR: #<pr>.
Branch: work detached on origin/wp/232i-mpm-bspline-mms; commit; `git pull --rebase origin
wp/232i-mpm-bspline-mms`; `git push origin HEAD:wp/232i-mpm-bspline-mms`.
Read the "232i" section of .claude/skills/ladruno-mpm/references/wp-briefs.md first: its
decisions K1–K4 and its gate table bind you. Never widen a band.
Independence: as in 232c, read tests/mpm_interface.md, tests/mpm_implicit_interface.md and the
cases files; never tests/mpm_reference.py or tests/mpm_implicit_reference.py.

ROUND 1 — stop and report after it (the orchestrator freezes it):
  tests/mpm_bspline_interface.md, at most ~3 pages, deciding with code anchors quoted:
  - the basis: the 1D functions and derivatives per degree, the support rule (which nodes, from
    which cell), admissibility, how it plugs into Basis / axisSupport / particleSupport
    (LadrunoMPMCore.h:74-238) and every site that switches on the basis today (grep it: the
    command parser, the LadrunoMPMIO.cpp BASIS read and write, LadrunoMPMImplicit.cpp support and
    serialization, LadrunoMPMPattern.cpp summary and output);
  - K4: the boundary treatment, with the reason, what it does to `fix` and `wall` faces and to
    free ghost nodes, in mode F and in mode A, and what it changes for GIMP and linear (nothing,
    or say why);
  - mode A: the particle element's node set under the new support (3 or 4 nodes per axis),
    re-pointing, and any size assumption (27 / 64 nodes) in LadrunoMPMImplicit.cpp;
  - the MMS body-force path: a generic per-particle body force (no MMS formula in C++), settable
    from Python each step as one flat list, honoured by P2G in mode F and by the particle element
    in mode A; its lifecycle under `reset` / `wipe` / restart; Tcl may lack it (say so in the
    guide). Also the per-particle initial velocity the MMS needs;
  - the kernel facts the MMS rests on, quoted from the interface documents: the strain measure
    the hypo kernel accumulates, the stress it returns (Cauchy or Kirchhoff), and the
    plane-strain embedding;
  - K2 and K3 as the orchestrator relays them; the `mpm grid` syntax.

ROUND 2 — on the orchestrator's "frozen":
  Implement, batching edits; ask the orchestrator for each build (all five targets) and each
  mutant build — never start build.bat yourself.
  - tests/mpm_mms.py: the manufactured fields, accelerations and body force (numpy; sympy is not
    installed — derive by hand, show the derivation in the docstring).
    tests/test_mpm_mms_reference.py checks the derivation independently: a centred
    finite-difference evaluation of ρ₀ a − ∇_X · P − ρ₀ b at random (X, t) is ≤ 1e-6 relative,
    never against MPM.
  - tests/test_ladruno_mpm_bspline_*.py (zone_a), Python, reading a cases file
    tests/mpm_bspline_cases.txt that holds each row's parameters and its gate copied from the
    table; a classic Tcl test for G1 and G7 (both shipped engines). Gated legs run in Zone-A only
    if they fit its time budget (say what each costs); the full ladder is the warrant run, which
    the orchestrator launches, with results JSON committed under tests/data/mpm/bspline/ carrying
    ladrunoBuild, Δt, crossings per rung, e per rung and p.
  - The g++ core check (tests/mpm_engine_check.cpp, tests/test_mpm_engine_core_cpp.py) extended
    for G1.
  Registrations per fork-checklist.md: the ledger fragments WP-232i-mpm-bspline-mms.md
  (implementations; vanilla rows + `// Ladruno` markers for any upstream file), a banner line,
  the guide (Ladruno_implementation/232_ladruno_mpm_guide.md: the basis, its boundary treatment,
  the body-force command, what Tcl lacks), stamps. Draft the ADR amendment text in your report;
  the orchestrator writes it.
Done when: every gated row green (or a reported miss with diagnosis), R1 reported in
Ladruno_implementation/232i_mpm_bspline_mms.md (method, tables per row, the orders, the
crossings), ci/ gates clean, ladrunoBuild matching your head on the exe and the pyd. Report per
the Worker rules.
```

## 232j — phase 4: anti-locking in rate form

Kickoff 2026-11-17 (orchestrator, Opus 5.5, started by the program orchestrator per #1075).
Branch `wp/232j-mpm-antilocking`, one draft PR for the WP. Base `ladruno` at `1d2a692ba` (232i
#1099 merged as `0b26fdb12`; release `ladruno-v20261116` predates 232i). ADR §9 phase 4. The
footing gate owed by the 232f amendment is redesigned first. The ADR carries it as
"Amendment (232j kickoff, 2026-11-17)", written with this section and before any phase-4 MPM
number. Done when the gate table below is green (or each miss is reported with its diagnosis, never
a widened band), the warrant package is in the PR, and the ADR carries "Amendment (232j, <date>)"
recording the scheme, its modes and bases, and its default.

Read before dispatch: #1103 (232i review follow-ups), #1095 (mode-A crash in a step-cut retry, a
hazard for the footing ladder), #1096 (the tangent rebuilds B nen² times), #1097, #1098 (mode-A
memory, 27-node supports), #1090. None is in scope unless the anti-locking rewrites the same code,
and none may regress.

**Model and effort, re-decided at kickoff.**
- C++ worker: opus · high (the roster's fork C++ row). The work is new kinematics in both modes, a
  consistent mode-A tangent and possibly a new element coupling.
- Verification worker 232j-V: opus · high, a different agent (the roster's oracle row). It
  produces L3's FE references and never reads the anti-locking C++.
- Ready-flip critic: fable · high, from the other family. Its review plus `/code-review high` is
  the review of record (ADR-87 D4, as for 232c, 232f, 232g and 232i). No owner verdict is
  pre-booked.

**Decisions fixed at kickoff** (the C++ worker may object in its Round-1 report; the orchestrator
decides):

- **K1 — the scheme is the Round-1 decision.** Two candidates, plus any other the worker finds:
  - the cell-averaged volumetric increment, an F-bar in rate form (de Souza Neto et al. 1996;
    for implicit MPM, Coombs, Charlton, Cortis & Augarde 2018);
  - a projection of the volumetric increment to the grid and back, of reduced degree for
    B-splines (B-bar / F-bar projection; Hughes 1980; Moutsanidis et al. 2020).

  The criteria:
  - L1 exactness on homogeneous fields;
  - all four bases (`gimp`, `linear`, `bspline2`, `bspline3`);
  - the same scheme in mode F and mode A;
  - a consistent mode-A tangent;
  - the locality and memory cost in mode A (K4);
  - behaviour when particles change cell between steps.

  Use the literature only: no GPL, LGPL or unlicensed source (Worker rule 5). Verify every
  citation before it enters the repo.
- **K2 — rate form.** The modification acts on the step's increment: ΔF in mode A, the velocity
  gradient or ΔF in mode F. It never acts on the total F. The hypo kernel's incremental entry
  (Fₙ, ΔF → strain increment and rotation, ADR §6) is kept. In plane strain the volumetric
  correction acts in-plane (d = 2), and the out-of-plane ΔF_zz stays 1.
- **K3 — opt-in, default off.** It is a deck-level setting; the syntax is decided in Round 1. With
  it off, every answer is bit for bit today's. Flipping the default is the owner's later call.
- **K4 — mode-A coupling.** If the scheme makes a particle's internal force depend on other
  particles' supports, Round 1 decides the element structure: either an enlarged particle-element
  node set or one element per cell. It quotes the node and DOF counts, the `FE_Element` pool limit
  (`MAX_NUM_DOF` = 64, #1098) and the bytes per particle-element it implies.
- **K5 — contracts are unchanged** unless Round 1 shows the output must record the setting. That
  change goes to T, never into apeGmsh from here.
- **K6 — the footing ladder (L4, L3).** MPM and FE get the same solution controls, per accepted
  step: KrylovNewton (NormUnbalance 1e-6, 40), NewtonLineSearch (1e-6, 60), then KrylovNewton
  (1e-5, 80), which counts as "relaxed"; on failure, halve the step. This is the ladder
  `tests/mpm_implicit_reference.py::fe_footing` used. Failed attempts and relaxed steps are
  reported per run. The C++ worker maps the ladder onto `mpm analyze` in mode A (its cut-and-retry,
  or per-step algorithm switching in the harness) and says how.

**Pre-registered gates** (fixed here, before any phase-4 MPM number exists; a band is never
widened after a run — a miss is reported as a miss, with its diagnosis). "On" means anti-locking
on, "off" means off.

| Row | Gate |
|---|---|
| L1 kinematics (both modes, every basis, on) | (a) An affine grid velocity field, deviatoric and volumetric, gives every particle a modified increment equal to the unmodified one, rel 1e-12. (b) The S1 rigid-body rotation of a stressed body shows zero spurious stress at its frozen tolerance. (c) In plane strain, ΔF_zz = 1 exactly. Checked in the g++ core check where the scheme is header-level, and through the engines. |
| L2 mode-A tangent (on) | Element tangent against a centred finite difference (FD converged over ε 1e-8..1e-10, as 232f measured). Elastic undeformed: ‖K − K_FD‖/‖K‖ ≤ 1e-5. Elastic at 0.2 % strain: at most 2× the off value at the same state. |
| L3 von Mises strip footing, closed form (mode A, QS, `gimp`, ppc 2 per axis) | The domain, grid, staging and rough `gridvel` footing of `footing_mildcone_r1..r3` (4 / 8 / 16 cells across B). Weightless soil, no surcharge. Material: `LadrunoJ2`, perfectly plastic, E = 1.0e5 kPa, ν = 0.3, σ_y = 100√3 kPa (k = 100 kPa). Push 2 mm per step to s = 0.1 m. The observable is q_lim, the mean footing pressure over s ∈ [0.03, 0.05] m. **Gate:** at r3, q_lim / ((2 + π) k) ∈ [0.97, 1.08], with (2 + π) k = 514.16 kPa. **Admission:** the band holds only if 232j-V's FE (LadrunoBrick B-bar, same domain, staging and ladder) lands inside it at 16 across B, measured before any MPM number. If it does not, the row's gate becomes L4's form: MPM within ±10 % of that FE value at each rung. The FE std (locked) legs are reported. If they sit inside the band at r3, the row is recorded as non-discriminating at r3. |
| L4 non-associated mild-cone footing, redesigned (mode A, QS, `gimp`, ppc 2) | `footing_mildcone_r1..r3` unchanged except for anti-locking on and K6's ladder. **Gate:** at each of r1, r2 and r3, q_max is within ±10 % of the FE B-bar q_max on the same domain at that resolution (`fe_footing.json`, build `f37e2ff19`): 141.43, 130.66, 122.53 kPa. The bands are [127.29, 155.58], [117.60, 143.73], [110.28, 134.78]. Each run must reach s = 0.2 m. A stop is a miss. Reported: q/exact (the retired [0.90, 1.10] × 138.907 band's verdict), the plateau, s at q_max, failed attempts and relaxed steps. **Outcome (2026-11-17, build `4df5764b7`):** r1 129.23 and r2 120.52 pass; r3 is a miss, stopping at s = 0.094 m with q_max 115.64 inside its band. It is Newton stagnation of non-associated DP, which the FE shows too. The owner accepted the miss; the cause is tracked in #1107. |
| L5 orders kept (on) | MMS G2 (mode F, small amplitude) and G4 (mode A, QS), `bspline2`, the 232i ladder and definitions: p ≥ 1.8 and e non-increasing. |
| L6 unchanged answers (off) | The whole existing MPM suite passes unchanged (report the counts). Off reproduces 232i's warrant figures bit for bit on the rows re-run. Per-particle-element memory with off grows ≤ 5 % for `gimp` (10⁴, the M0b M3 method). With on, it is reported. |
| L7 mutation | The MPM family's ZERO and SCALE scores (`mutation_gate.py`) at or above their floors, 0.46 and 0.38. Plus a Zone-A **switch oracle**: a small locking-sensitive case, such as L3 at r1 if it fits Zone-A's time, where on is within L3's admitted gate and off is outside it. The oracle goes red if the averaging is the identity. |
| L8 surfaces | The setting round-trips through `mpm info summary` and the mode-A database restart (the pack version bumps; a v2 database restores with off). Classic Tcl and Python give the same answer on one L1 case. Unsupported combinations (by basis, mode or face) are refused loudly. |
| R reported, not gated | Off against on for L3, L4 and L5. The L4 sequence trend. `bspline3` footing legs where the cost allows. L3 at r1 in mode F with local damping, on against off. The 232i self-weight column at 20 % compression per basis, on. The granular-column runout, on against off. Wall time per step, on against off. |

`test_footing_gate` (the 232f strict xfail on the old band) is replaced by the L4 test. The old
band's verdict stays in the results record and is never re-gated.

**Freeze (2026-11-17, orchestrator, on `tests/mpm_antilocking_interface.md` @ `ba201819e`).**
Frozen as J1–J13, before any phase-4 MPM number, with these conditions:
- **K1 = J1:** a cell-centre F-bar on the step increment, one rule for all four bases and both
  modes. The quadrature weight is the carried volume, J9. The critic is asked to check that
  choice against de Souza Neto's actual-J weight.
- **Amendments 1–5 are accepted.**
  - L1(a) covers particles whose support is not truncated; the truncated result is reported.
  - L8 round-trips through `mpm info antiLocking`.
  - **L7's switch oracle at r1** uses L3's fallback form. On must lie within ±10 % of the FE
    B-bar q_lim at r1 (570.10 kPa; band [513.09, 627.11]). Off must lie outside it. The reason:
    L3's closed-form band is an r3 gate, and the FE B-bar itself sits at 1.109× at r1.
  - The citation is corrected.
  - Mode A refuses `fbar` with symmetric-storage SOEs, PARDISO's symmetric mode included.
- **Q1, the halving floor:** the FE stops when a sub-step falls below 1e-6 m of footing
  displacement (`fe_footing`'s `dsc < 1e-6`). That is h < 5e-4 in the cases' pseudo-time
  (v = 0.002 per unit). Put it in the `ladder` key, and count halvings.
- **Q2:** cases default to on unless they say `antilocking off`. Accepted. Every off leg is
  explicit.
- **Q3:** the setting is fixed per epoch, like `mpm mode`. Accepted.
- **Q4:** measure the switch oracle's Zone-A cost in Round 2. If both legs together exceed
  ~3 min, propose a smaller case before any number.
- **F1 — the fallback fraction is reported.** J5's fallback is the identity wherever a cell has
  an inactive corner. Every L3 and L4 run reports the number of averaged and fallback
  particle-steps, and the fallback share inside the footing's plastic zone (x ≤ 3B, y ≥ −2B).
  A miss with a large fallback share there is diagnosed as that, never re-tuned.
- **F2 — L3's MPM legs at r2 and r3 wait for its admission.** The admission is the FE B-bar
  at 16 across B, which the orchestrator's FE run records. Until then, no L3 run at r2 or r3.


**Critic round (2026-11-17, Fable · high, verdict AMEND at `d09ee3ee6`; each finding verified
in-tree by the orchestrator).** The warrant run on `d09ee3ee6` was stopped after L4 r2. That build
gave:
- L5: G2 p = 2.040 and G4 p ≈ 2.0, both non-increasing.
- L3 at r1, reported: q_lim 550.0 = 1.070×.
- L4 r1: q_max 127.39, inside [127.29, 155.58] by 0.08 %.
- L4 r2: **stopped** at t = 74.5 on the halving floor, with q_max 119.06 before the stop. That is
  a miss under "a stop is a miss".
- L4 r3 was not run.

These numbers are context for the corrected code, not its warrant. No band changes.

- **C1 (major, a defect, reproduced by the orchestrator).** In mode F, `fbar` does not conserve
  energy.
  - The probe: a soft free block (E = 2e5, ν = 0.45, gravity, FLIP α = 1). KE_max is 503 with
    `fbar` off and 1.33e4 with it on, still growing. The same happens at 0.05 dt_crit and with
    ppc 3. With v0 = 1 m/s, particles leave the grid.
  - It stays stable when the block is stiff or confined.
  - The cause: J1–J9 drive the stress with ℓ̄ = ℓ + βI_d, β = (1/d)(tr ℓ_c − tr ℓ). But the force
    w̄σ̄·g_a is not work-conjugate to ℓ̄, so the scheme has no energy balance.
  - **Decision: the force becomes the conjugate one, in both modes** (one scheme, K1):
    f_a = w̄ [σ̄ g_a + p̄_d (ĝ_a − g_a)], with p̄_d = (1/d) Σ_{i<d} σ̄_ii, and ĝ_a = ∇N_a(x_c)·ΔF_c⁻¹
    on the cell's corners and 0 elsewhere. Then Σ_a f_a·v_a = w̄ σ̄ : ℓ̄ exactly. This is the
    B-bar force (Hughes 1980), carried in rate form beside the F-bar kinematics.
  - J9 and J10 are re-derived for it, and the L2 FD check re-verifies the tangent.
  - If the conjugate form cannot pass L9 in mode F, `fbar` is refused in mode F. That refusal is
    recorded as a miss, with the ADR amendment saying so; mode A keeps `fbar` if its L9 leg passes.
- **C2 (major).** Mode F's L1 did not run the production path: `affineCheck` re-assembled ΔF and
  ΔF_c itself. The fix: one shared helper that both `Engine::advance` and `affineCheck` call, or a
  mode-F affine run through `advance` itself.
- **C3 (major).** The symmetric-storage refusal has gaps. It misses SProfileSPD and the
  distributed SPD SOEs, and it misses `system Mumps -matrixType 1|2`, which needs a matType
  getter. The refusal message must not recommend Mumps without its unsymmetric type.
- **C4 (minor).** J5's fallback almost never fires: 0 in every warrant record. The guide's
  "free surface's outer cells" text is corrected to what the code does. F1 is still reported.
- **C5 (minor, L2 strengthened, thresholds unchanged).** New rows:
  - ν = 0.3 and 0.45;
  - 3D;
  - `bspline2` with a `fix`-face image (the B8 path);
  - strained states at 1 % and 5 %, at the same bars: ≤ 1e-5 undeformed, ≤ 2× off when
    strained.
- **C6 (minor).**
  - `judge` "out" must fail on NaN.
  - `debug packVersion` must not be sticky.
  - L1b rotates a pre-stressed body, as the gate says.
  - A mode-A debug hook must not block a later `mpm antiLocking` before the first step.
- **New gate L9 — energy (pre-registered here, before any number on the conjugate form).**
  The critic's probe: 2D free block 1 × 1.5 m on `fix ymin xy`, h = 0.25, ppc 2, ρ = 1800,
  ElasticIsotropic E = 2e5, gravity. The bases are `gimp`, `linear` and `bspline2`, with
  ν = 0.3 and 0.45.
  - **Mode F, FLIP α = 1, dt = 0.15 dt_crit, 4000 steps:** KE_max(on) ≤ 1.5 × KE_max(off), and
    no particle leaves the grid. The variant with v0 = 1 m/s runs at 0.05 dt_crit for 9000 steps
    under the same bar.
  - **Mode A, Newmark average acceleration, 400 steps of the off run's stable Δt:**
    KE_max(on) ≤ 1.5 × KE_max(off).
  - **The S1 explicit energy-ledger case, with `fbar` on:** its frozen tolerance holds.
- **The L4 r2 stop:** diagnose it on the conjugate build (which particles, which material state,
  why the ladder fails where the FE ladder survived), never by changing a band or the ladder.

### 232j-V — verification references for L3 (opus · high; never reads the anti-locking C++)

```text
<header>  ADR sections: §8 (the footing row and its 232f and 232j-kickoff amendments), §9
phase 4. Your PR: #<pr>.
Branch: work detached on origin/wp/232j-mpm-antilocking; commit; `git pull --rebase origin
wp/232j-mpm-antilocking`; `git push origin HEAD:wp/232j-mpm-antilocking`. Another worker pushes
C++ to the same branch. You touch only tests/mpm_implicit_reference.py (a new leg),
tests/test_mpm_implicit_reference.py, tests/data/mpm/implicit/fe_footing_j2.json and a new cases
block in tests/mpm_antilocking_cases.txt (create it; the C++ worker reads it).
Independence: read nothing the C++ worker adds under SRC/domain/pattern/ladrunoMPM/. You may read
and run the fork's FE elements.
Binary: <pyd path>, built by the orchestrator; print ops.ladrunoBuild() into every record.
Goal: L3's FE references, before any MPM number exists.
  1. A leg `leg_footing_j2` beside `leg_footing`, reusing `fe_footing`'s slab, staging, rough
     footing and ladder. Material: LadrunoJ2, perfectly plastic. Quote its factory line
     (SRC/material/nD/LadrunoJ2.cpp:54) and derive K and G from E = 1.0e5, nu = 0.3; s0 = 100
     sqrt(3). No surcharge, weightless. Push to s = 0.1 m at 2 mm per step. Legs: B-bar and std at
     4 / 8 / 16 across B.
     Output per leg: q(s), q_max, q_lim = the mean of q over s in [0.03, 0.05] m, failed attempts,
     relaxed steps and wall time.
  2. The closed form (2 + pi) k = 514.16 kPa. Cite it, and state why the result does not depend
     on roughness for a weightless Tresca/von Mises strip (Prandtl vs Hill), and the plane-strain
     relation k = sigma_y / sqrt(3) for von Mises.
  3. Report the admission verdict for L3: is the B-bar q_lim at 16 across B inside
     [0.97, 1.08] x 514.16? Is std inside? Never adjust the band.
  4. Write the cases block for footing_j2_r1..r3 in the style of tests/mpm_implicit_cases.txt,
     with the gate copied from the "232j" table and the FE values as `ref` lines.
  5. tests/test_mpm_implicit_reference.py checks the closed form and the FE record's shape. It
     never checks against MPM.
Long FE runs: tell the orchestrator the command and it launches them. Do not run anything over
10 min yourself.
python3.12; ruff F811 clean; the fork-checklist gates before every push. Report per the Worker
rules.
```

### 232j-C — anti-locking in C++ (opus · high; a different agent from 232j-V)

```text
<header>  ADR sections: §4 (walls W1–W17), §5.2, §5.3 and its 232f amendments, §5.4, §6, §8 (the
footing row and its amendments, including 232j kickoff; the mutation note), §9 phase 4. Your PR:
#<pr>.
Branch: work detached on origin/wp/232j-mpm-antilocking; commit; `git pull --rebase origin
wp/232j-mpm-antilocking`; `git push origin HEAD:wp/232j-mpm-antilocking`.
Read the "232j" section of .claude/skills/ladruno-mpm/references/wp-briefs.md first: K1–K6 and
the gate table bind you. Never widen a band.
Independence: as before, read tests/mpm_interface.md, tests/mpm_implicit_interface.md,
tests/mpm_bspline_interface.md and the cases files. Never read tests/mpm_reference.py or
tests/mpm_implicit_reference.py.
Context: issues #1095, #1096, #1098 and #1103 touch code you will change. Do not fix them unless
your change rewrites that code; never regress them.

ROUND 1 — stop and report after it (the orchestrator freezes it):
  tests/mpm_antilocking_interface.md, at most ~3 pages, deciding with code anchors quoted:
  - K1: the scheme. Give its equations in rate form for each basis, and the averaging domain
    (cell, patch or grid projection) and its weights. Show why L1(a) holds exactly. Say what
    happens when a particle changes cell or a cell holds one particle, and how an empty or
    one-particle cell is handled.
  - Mode F: where in the step it enters (after the grid update, before the material call), and
    its cost.
  - Mode A (K4): the dependence of a particle's internal force on other DOFs, and the element
    structure that carries it (enlarged node set or per-cell element). Give the node and DOF
    counts per basis in 2D and 3D against MAX_NUM_DOF, the bytes per particle-element (the M0b M3
    method), the consistent tangent including the geometric term, re-pointing, the empty-node
    treatment, cut-and-retry with nothing to undo, and the restart pack.
  - K3: the syntax, the default and the refusals.
  - K6: how the footing ladder maps onto mode A.
  - K5: whether the output must record the setting.
  - The harness: the cases file tests/mpm_antilocking_cases.txt (232j-V adds the L3 block; you
    add L1, L2, L4, L5 and L7), and the test files.

ROUND 2 — on the orchestrator's "frozen":
  Implement, batching edits. Ask the orchestrator for each build (all five targets) and each
  mutant build; never start build.bat yourself.
  Tests tests/test_ladruno_mpm_antilocking_*.py (zone_a, within its time budget; say what each
  costs). The full L3 and L4 ladders and the L5 ladders are the warrant run, which the
  orchestrator launches. Results JSON go under tests/data/mpm/antilocking/ and carry
  ladrunoBuild, the ladder counts and per-rung values. Replace test_footing_gate with the L4
  test, and keep the old band as a reported value.
  Classic Tcl for L1 and L8. The g++ core check extended for L1 where it applies.
  Registrations per fork-checklist.md: the ledger fragments WP-232j-mpm-antilocking.md
  (implementations; vanilla rows and `// Ladruno` markers for any upstream file), a banner line,
  the guide (Ladruno_implementation/232_ladruno_mpm_guide.md: the setting, its scheme, its
  limits, when to use it), stamps.
  Draft the ADR amendment text in your report; the orchestrator writes it.
Done when: every gated row is green (or a reported miss with diagnosis); the reported rows are in
Ladruno_implementation/232j_mpm_antilocking.md (method, tables, off against on); the ci/ gates
are clean; ladrunoBuild matches your head on the exe and the pyd. Report per the Worker rules.
```

## 232k — phase 5: soil models at particles, the SANISAND cell-crossing meter first

Kickoff 2026-11-17 (orchestrator, Opus 5.5, started by the program orchestrator per #1075).
Branch `wp/232k-mpm-soil-models`, one draft PR for the WP. Base `ladruno` at `3e738308c` (232j
#1105 merged; 232i and 232j are not in a release, the last being `ladruno-v20261116`). ADR §9
phase 5. The ADR carries "Amendment (232k kickoff, 2026-11-17)", written with this section and
before any soil-model MPM number. Done when the gate table below is green (or each miss is
reported with its diagnosis, never a widened band), the warrant package is in the PR, and the ADR
carries "Amendment (232k, <date>)" with the admission-table rows filled.

**Facts checked at kickoff** (source-verified on `3e738308c`; the workers re-quote, never recall):
- **NorSand is not on `ladruno`** (`grep -ril norsand SRC/` is empty). `LadrunoNORSAND` is the
  draft #896 (WP-144). 232k does not pull it in. NorSand, its matched-material Berkeley track and
  the §8 tailings flow-slide row are **deferred until #896 lands** (amendment).
- **Scope:** `LadrunoSANISAND` (33019; `getCopy` gives 33020 3D / 33021 PlaneStrain, order 6 / 3;
  integration schemes 1 ModifiedEuler (deck default), 2 BackwardEuler_CPPM, 129 SAS-ME, and
  ADR-92 IMPL-EX opt-in) and `LadrunoHS` (33029; one class, `dim_` modes; `-smallStrain` /
  `LadrunoHSSmall`; `-implex`; `-sigmat`). PM4Sand and PDMY stay out.
- **Threads (W14):** `ManzariDafalias::ladrunoThreadSafeUpdate` returns false (`:1235-1269`), the
  SANISAND subclass delegates (`LadrunoSANISAND.cpp:2741-2758`), and `LadrunoHS` inherits
  `Material::ladrunoThreadSafeUpdate()` = false. Soil particles are serial.
- **Refusals (W17):** both declare commit refusals (SANISAND: `cpp:4623-4695`, `:4884-4926`; HS:
  `cpp:690-745`, `:832-834`). A refusal at particles aborts the run. The §8 P5 row (refusal at
  commit aborts cleanly, with the first refusing material) moves here, per review row 6.
- **Initial stress:** SANISAND has none (`initialize()` zeroes σ, `ManzariDafalias.cpp:1044-1047`);
  HS has `-sig0` / `-sig03` (uniform per tag). MPM interface Q2 refuses an initial stress on a
  non-elastic particle material (`LadrunoMPMPattern.cpp:442-447`, `:505-507`, `:531-534`,
  `:968-976`): the offset is added outside the material and would bypass its yield check.
- **Stages:** SANISAND's stage is a class-wide static `mElastFlag` (`ManzariDafalias.h:557`); id 1
  ("updateMaterialStage", `theInt`) and id 5 ("materialState", `theDouble`). The flip's
  `Elastic2Plastic()` (the p < 3·Pmin reset, `:5789-5794`) runs only on instances that
  `updateParameter` reaches; the others flip lazily through `ladrunoRunStageFlipOnce`, which does
  not call it (`LadrunoSANISAND.cpp:3390-3405`, `:5054`). HS keeps a mutex-guarded tag→stage map
  and seeds p_p, γ_p from the committed stress at 0→1 (`LadrunoHS.cpp:54-73`, `:517-577`).
  Mode F's `parameter … loadPattern` reaches every particle material
  (`LadrunoMPMPattern.cpp:347-364`); mode A's `updateMaterialStage` stops at the first element that
  answers (`MaterialStageParameter.cpp:79-81`). **The two modes therefore flip SANISAND
  differently today.**
- **Low p:** SANISAND `-Pmin` (default 1e-3 P_atm), `-Presidual` (0), `-pRe` (0, off), and the
  SAS-ME tension cutoff (WP-152, off). SAS-ME refuses a committed start with p + p_r ≤ 0
  (`LadrunoSANISANDSasME.cpp:72-75`). HS: p_min = max(1e-6, 1e-3 p_ref), `-sigmat`; it refuses
  |Δε|² > 25. ADR-93's floors are "available, not adopted; the owner decides".
- **Cost:** M0b (232g) SANISAND/DP plastic cost ratio **341 (222–477) at a strain increment of
  2.56e-4, 647 (462–920) at 9.44e-4**: 0.61 ms and 1.16 ms per plastic call. The increment a
  particle actually sees was not measurable while SANISAND was not admitted; 232k measures it.

**Model and effort, re-decided at kickoff.**
- C++ worker 232k-C: opus · high (the roster's fork C++ row). Admission, staging in both modes,
  the refusal path, the harness and the meter's MPM leg.
- Verification worker 232k-V: opus · high, a different agent. The FE references on the same
  domains and the one-point drivers. It never reads the MPM soil C++.
- Ready-flip critic: fable · high, from the other family, plus `/code-review high`: the review of
  record (ADR-87 D4, as for 232c, 232f, 232g, 232i and 232j). No owner verdict is pre-booked.
- Runs: esmeralda (the owner's instruction). A fresh `git clone --depth 1 --branch
  wp/232k-mpm-soil-models` into `~/ladruno_build_test/OpenSees_232k`; never touch
  `~/ladruno_build_test/OpenSees`. Soil runs are serial processes, run side by side (at most 20 at
  once on its 24 cores). The orchestrator launches every build and every long run.

**Decisions fixed at kickoff** (a worker may object in its Round-1 report; the orchestrator
decides):

- **K1 — admission.** Both materials join `admitMaterial`'s table for `ThreeDimensional` and
  `PlaneStrain`, in both modes, with every option that reads no element length (none does).
  Options with a mode-specific hazard are decided in Round 1 and either admitted with a test or
  refused loudly: SANISAND schemes 1, 2, 129; IMPL-EX on both (its time source `ops_Dt` at
  particles); HS `-implex`, `-sig0`/`-sig03` (K3). The §6 table's columns (type strings, hypo use,
  stage path, process statics, threads, length consumer) are filled for both rows.
- **K2 — serial.** Soil particles run serial (W14). A requested thread count above 1 with a soil
  particle is honoured as 1 with one warning, or refused: Round 1 picks, and says which the
  existing non-allowlisted materials already get.
- **K3 — geostatic state by staging only.** Q2 stays: no initial-stress offset on soil particles.
  The geostatic path is stage 0 (elastic) under the gravity ramp, then the flip to stage 1. HS's
  own `-sig0`/`-sig03` (inside the material, so inside its yield check) is admitted only if Round 1
  shows it is uniform-state-correct at particles; it is never the geostatic route.
- **K4 — one stage flip, identical in both modes.** Every soil particle in mode A and mode F
  reaches the same post-flip state as an FE Gauss point with the same committed stress, including
  SANISAND's `Elastic2Plastic` reset. Round 1 decides the mechanism (for example, mode A forwarding
  the stage to every particle element, or one `mpm`-level stage walk that both modes share) and
  quotes why `MaterialStageParameter`'s first-answer stop and the lazy flip cannot leave a particle
  un-reset. SANISAND's stage flag is process-static: an FE SANISAND and an MPM SANISAND in one
  process share it; the table's "process statics" column says so and the manual warns.
- **K5 — no floor decisions.** No default changes to `-Pmin`, `-Presidual`, `-pRe`, the SAS-ME
  cutoff or HS's p_min: those are the owner's (ADR-93). Gated cases keep their confinement (a
  surcharge or a confined domain) so that no gate depends on a floor choice. Free-surface cases are
  reported, with counts of particle-steps at the floor, clamps and refusals.
- **K6 — refusals abort.** A refusal in a gated case is a miss. The refusal test uses a real
  refusing material with a deterministic trigger (Round 1 picks: HS's |Δε|² > 25, or a SAS-ME start
  at p + p_r ≤ 0), in both modes: the run aborts, the particles are unmoved, the message names the
  material and its tag, and a later `mpm analyze` is refused as after `injectRefusal`.
- **K7 — mode A tangent.** SANISAND's and HS's tangents are unsymmetric unless IMPL-EX makes them
  symmetric. Mode A refuses symmetric-storage systems with a soil particle unless Round 1 shows
  the tangent is symmetric for the chosen options (reuse 232j's refusal list, its C3).
- **K8 — contracts unchanged** unless Round 1 shows the particle output must carry material
  internal variables (SANISAND's α, e, fabric; HS's hardening). Such a change goes to T, never
  into apeGmsh from here; the gates do not need it (the harness reads them through `mpm` queries).
- **K9 — the footing ladder and `fbar`.** Footing rows use 232j's K6 ladder (KrylovNewton
  1e-6/40 → NewtonLineSearch 1e-6/60 → KrylovNewton 1e-5/80 relaxed → halving, floor 5e-4 in
  pseudo-time) for MPM and FE alike, `mpm antiLocking fbar` on (FE reference B-bar), off reported.

**Pre-registered gates** (fixed here, before any soil-model MPM number; a band is never widened
after a run, a miss is reported as a miss with its diagnosis). References are the fork's own FE on
the same domain, a one-point FE driver, or closed form — never MPM, and never ADR-79's
1108–1152 kPa anchor. SANISAND deck: the fork's `_build_confined` parameters
(`tests/test_ladruno_sanisand.py:590`, kPa), scheme 1 with `-maxSubsteps` as the deck uses it.
HS deck: the oedometer/footing deck of `tests/test_ladruno_hs_material.py`. Both are quoted, with
their lines, in the cases file.

| Row | Gate |
|---|---|
| **K0 SANISAND cell-crossing meter — runs FIRST** (both modes; `gimp` and `bspline2`) | **Case:** the laterally confined column (2D plane strain, one cell wide, `fix xmin x`, `fix xmax x`), l0 = 10 m, h = 0.5, ppc 2, ρ = 2000, gravity. Stage 0 gravity ramp → flip → stage-1 hold → **translation**: the floor is a `gridvel` box under the column with v = (free, −v) and `track y` (`tests/mpm_implicit_interface.md:157`), so the stressed column slides down its frictionless walls through at least 4 cells while its exact stress stays static. Mode A QS (`LoadControl`, h/8 per step) and mode F (local damping as the S1 column; v small enough that inertia is at most 1 % of ρ g l0 — Round 1 states v and checks it). **Reference:** σ_yy = −ρ g (l0 − Y) exactly (equilibrium); σ_xx and σ_zz from 232k-V's FE column on the same domain, staging and material (no translation), at the particle's initial depth Y. **Meter:** m = RMS over particles of ‖σ_p − σ_ref(Y_p)‖ / (ρ g l0), over the in-plane components and σ_zz, after the full translation. Each material (ElasticIsotropic, DruckerPrager, LadrunoSANISAND, LadrunoHS) runs the identical case, plus a no-translation twin (v = 0). **Gate (admission of SANISAND for everything after K0):** in each mode and basis, A = m(SANISAND) / m(DP) ≤ 2, or m(SANISAND) ≤ 0.01. **Reported:** the twins, the crossing counts per particle, A for elastic and HS, and SANISAND's spurious state drift (Δe and ‖Δα‖ across the translation, against the twin). **A miss stops SANISAND** at this row: the orchestrator raises it to the owner before any further SANISAND MPM run. HS continues. |
| S1 one-point identity (mode A QS, and mode F where homogeneous) | A homogeneous block under a prescribed affine grid motion (every particle sees the same ΔF): drained triaxial compression, oedometer, and a 3-cycle isochoric shear, each from the stage flip. Reference: 232k-V's one-point FE driver on the same strain path. Where a fork FE element drives the same hypo kernel entry (`-geom hypo`, ADR-79), particle σ and internal state match to rel 1e-8; otherwise small strain at ‖ε‖ ≤ 5e-4, rel 2e-3 (Round 1 says which, quoting the kernel entry). Both materials, 3D and plane strain. |
| S2 K0 column against FE (mode A QS; `gimp`, `bspline2`) | K0's column without translation, then a surcharge equal to ρ g l0 on the top row (the surcharge mechanism of `footing_mildcone`). Gate: σ_yy against closed form within 0.03 of ρ g l0 (the S1 self-weight bar), and σ_xx RMS against 232k-V's FE column on the same domain within 0.02 of ρ g l0, at the end of stage 0, after the flip, and after the surcharge. Both materials. |
| S3 strip footing against FE (mode A QS, `gimp`, ppc 2, K9) | The 232j L4 domain and rough `gridvel` footing at r1 / r2 (4 / 8 cells across B) with a 20 kPa surcharge (K5's confinement), stage 0 gravity → flip → push to s = 0.05 B. Observables q(s) at s/B = 0.01, 0.02, 0.05. Gate: each within ±10 % of 232k-V's FE B-bar on the same domain, staging and ladder at that rung; each run reaches s = 0.05 B (a stop is a miss). HS at r1 and r2; SANISAND at r1 and r2 if K0 passes. r3 is reported only, and only if Round 1's cost estimate puts it under 12 h on esmeralda; the orchestrator decides at the freeze, before any number. |
| S4 stage parity (both modes) | After the flip on the S2 column (stage-0 state identical by construction), every particle's committed σ and internal state are bit for bit equal between mode A and mode F built on the same state, and equal to the FE Gauss point with the same committed stress, including particles with p < 3·Pmin (SANISAND's reset). |
| S5 refusal (both modes; the §8 P5 row) | K6's deterministic refusal, per material: the run aborts, the particles are unmoved (positions, F and σ bit for bit at the last accepted step), the message names material and tag, and a later `mpm analyze` is refused. Classic Tcl and Python both. |
| S6 restart (mode A) | The database round trip with SANISAND and HS particles mid-stage-1: the continued run equals the uninterrupted one bit for bit. |
| S7 refusals and options | Loud refusal for: symmetric-storage SOEs with an unsymmetric soil tangent (K7); each option Round 1 refuses; a soil particle with an initial-stress offset (Q2). K2's thread rule tested. `mpm info` reports the material, scheme and stage. |
| S8 unchanged answers | The whole existing MPM suite passes unchanged (counts reported); the 232j warrant figures re-run on L1 and L4 r1 bit for bit. |
| S9 mutation | The MPM family's ZERO and SCALE scores at or above 0.46 and 0.38. Plus a Zone-A soil oracle: S1's triaxial at small size, which goes red if the stage flip or the material update is skipped at particles. |
| S10 engines | Classic Tcl and Python give the same S1 triaxial answer. |
| R reported, not gated | The cost per particle-step of SANISAND and HS at particles, with the strain increment particles actually see (closes M0b's 222–920 range). The granular-column collapse (Lube) with SANISAND and HS in mode F: runout against DP, floor and clamp counts, refusals. The S3 legs with `fbar` off. HSSmall's S1 cyclic shear. |

**Order.** Round 1 (both workers, in parallel) → the orchestrator's freeze → 232k-C implements
K1–K4 and the meter harness → **the K0 meter runs before any other soil MPM case** → K0's verdict
(a miss goes to the owner) → S1–S10 → the warrant run on one build.

**Freeze (2026-11-17, orchestrator, on `tests/mpm_soil_interface.md` @ `766cc1450` and 232k-V's
Round-1 report).** Frozen as D1–D23, before any soil-model MPM number. The orchestrator verified the
objections in-tree: `track` moves only the lower bound (`LadrunoMPMImplicit.cpp:661`); mode F never
reads `gridVels` (its only uses are the B-spline ghost check and a clear, `LadrunoMPMPattern.cpp:669`,
`:873`); a mode-F trial refusal reverts and fails the step (`LadrunoMPMCore.h:990-995`). Rulings:

- **Units and case (D19).** kN, t, m, s, kPa: ρ = 2.0 t/m³ (the kickoff's "2000" is kg/m³),
  g = 9.81, ρ g l0 = 196.2 kPa. The K0 column is **two cells wide** (B-spline grids need 2 cells
  per axis); the problem stays 1D. No `fix ymin`; the floor is the D17 `profile … translate`
  gridvel; translation D = 2.5 m (5 cells). Accepted as amendments to the K0 row; no band changes.
- **D17 accepted:** `mpm gridvel … profile <t0> <T_r> <t1> translate <axis>`; `start`/`track`
  unchanged.
- **D18 accepted, and it is a defect fix:** mode F applies `gridvel` (both projections) and refuses
  `pforce` loudly. Filed as an issue for the released binaries; the fix lands here, with D18's test.
- **D20/D21 accepted** as the K0 timelines. The mode-F inertia ratio (|dP_y/dt| / W over 100 steps)
  is gated at ≤ 1 % at the measurement; its maximum is reported. The meter is read after the
  ramp-down and settle (D21: local damping is not Galilean invariant).
- **D23 accepted** as the pre-registered K5 contingency: a gated K0 run with any floor or clamp
  particle-step is not admitted, and the row re-runs with a bonded elastic cap (one cell, 20 kPa),
  with its FE reference.
- **K0 references:** 232k-V's FE column on the same domain for all four materials (its Round 1
  shows it equals K0 = ν/(1−ν) for SANISAND's stage 0 and K0nc for HS; DP is staged with
  `materialState` 0 first). One source for every leg.
- **D7–D11 (K4) accepted:** one host walk, three routes (`mpm stage`, `parameter … loadPattern`
  registering the host, mode A's `updateMaterialStage` handing the host the registration), and the
  D10 guard refusing `mpm analyze` on a stage mismatch.
- **S4 corrected (a fact, not a band):** SANISAND's reset fires at tr σ < 3·Pmin, i.e. **p < Pmin**
  (`ManzariDafalias.cpp:5789`), not p < 3·Pmin as the kickoff text said. The S2 column never reaches
  it (top p ≈ 1.6 kPa). S4 adds a **low-g twin** of the S2 column, g chosen so that at least one
  particle row sits at p < Pmin at the flip; the test asserts ≥ 1 reset particle in both modes. Its
  reference is the one-point FE driver flipped by a route that reaches that point (direct
  `updateParameter`), so `Elastic2Plastic` runs there too.
- **S5 amended (D12–D14):** the W17 abort leg is SANISAND `-implex -maxSubsteps 1` on one stage-1
  shear step, both modes, after the one-point FE confirms the trigger fires. HS has no reachable
  commit refusal at particles: its leg is the trial refusal (γ = 6 simple shear), which must fail
  the step loudly, name material and tag, keep the last accepted step, and **not** abort.
- **S1 tolerance under substepping (pre-registered now):** triaxial and oedometer at rel 1e-8
  against `LadrunoBrick -geom hypo -formulation std`. The 3-cycle shear is at rel 1e-8 up to the first step
  at which the MPM and FE substep counts differ (reported); from that step on it is at rel 1e-3
  (10 × the deck's ModifiedEuler TolE 1e-4).
- **S3 SANISAND: the FE deck is unreachable as kicked off.** 232k-V: with the deck's defaults
  (scheme 1, `-maxSubsteps 0`, TolE 1e-4) the FE B-bar r1 cannot take its first push step on the
  K9 ladder (15 failed attempts to the floor). Settled before any MPM number by a bounded FE probe
  (each leg ≤ 2 h on esmeralda), in this fixed order; **the first variant whose FE B-bar reaches
  s = 0.05 B at r1 and r2 on the K9 ladder becomes the S3 SANISAND deck**, for MPM and FE alike:
  1. scheme 1 with `-maxSubsteps 1000` (#792's cap);
  2. scheme 129 (SAS-ME) with its defaults;
  3. scheme 2 with `-cppmOnFail refuse` (WP-130);
  4. scheme 1 with `-implex -maxSubsteps 1000` (ADR-92).
  If none reaches, S3 SANISAND becomes reported only, with 232k-V's diagnosis, and the orchestrator
  raises it to the program orchestrator before any S3 SANISAND MPM run. K0, S1 and S2 keep the
  deck defaults.
- **S3 rungs:** HS r1 and r2 gated, r3 run and reported. SANISAND r1 gated; r2 gated if its FE
  reaches (MPM cost 8–15 h, accepted, run side by side); r3 not run (30–60 h).
- **Costs (§7) accepted** as upper bounds; the R row measures the real increments.
- **232k-V's FE flips:** `updateMaterialStage` reaches one element's points; the rest flip lazily
  without `Elastic2Plastic`. That is inert on the gated FE domains (η = 0.86 < M_c, p ≫ Pmin), and
  232k-V records it in each FE JSON.

### 232k-V — FE references (opus · high; never reads the MPM soil C++)

```text
<header>  ADR sections: §6 (the materials table), §8 (the cell-crossing meter row, the footing row
and its 232f / 232j amendments), §9 phase 5, §11 (the value proposition), and the amendment
"232k kickoff". Your PR: #<pr>.
Branch: work detached on origin/wp/232k-mpm-soil-models; commit; `git pull --rebase origin
wp/232k-mpm-soil-models`; `git push origin HEAD:wp/232k-mpm-soil-models`. Another worker pushes
C++ to the same branch. You touch only tests/mpm_soil_reference.py, tests/test_mpm_soil_reference.py,
tests/data/mpm/soil/*.json and the reference blocks of tests/mpm_soil_cases.txt (create it; the
C++ worker reads it and adds the MPM blocks).
Independence: read nothing under SRC/domain/pattern/ladrunoMPM/. You may read and run the fork's
FE elements and materials, and tests/mpm_implicit_reference.py's fe_footing (you extend it).
Read the "232k" section of .claude/skills/ladruno-mpm/references/wp-briefs.md first: the facts,
K1–K9 and the gate table bind you. Never widen a band.
Binary: <pyd path>, built by the orchestrator; print ops.ladrunoBuild() into every record.

ROUND 1 — report before running anything long (at most 300 words + a plan):
  - the two decks, quoted with file:line (SANISAND `_build_confined`, HS oedometer/footing), in
    the units the MPM cases use; the stage sequence each deck uses (command, id, theInt/theDouble);
  - for S1: whether a fork FE element drives the same hypo kernel entry as the particles
    (`-geom hypo`, ADR-79): quote it, or say small strain applies;
  - the FE column (K0, S2) and footing (S3) models: element, B-bar, mesh = the MPM grid, staging,
    ladder; and a cost estimate per run (use the M0b ratio, 341 / 647 against DP, and your own
    timing of one stage-1 FE step);
  - anything in the gate table you consider unreachable by the FE itself (as 232f's footing was),
    with the number that shows it — before any MPM number exists.
ROUND 2 — on the orchestrator's "frozen":
  1. K0 references: the FE column σ_xx(Y), σ_zz(Y) for each material, at the end of stage 0, after
     the flip and after the hold; JSON with the build hash.
  2. S1 one-point drivers: triaxial, oedometer, 3-cycle isochoric shear, both materials, 3D and
     plane strain; the strain paths written so the MPM harness can impose the same ΔF.
  3. S2 FE column with surcharge; S3 FE B-bar footing r1 / r2 (r3 if the freeze says so) with the
     K9 ladder, q(s) at s/B 0.01 / 0.02 / 0.05, failed attempts, relaxed steps, wall time.
  4. tests/test_mpm_soil_reference.py: shape and closed-form checks; never against MPM.
Long FE runs: give the orchestrator the command; it launches them on esmeralda. Nothing over
10 min yourself. python3.12; ruff F811 clean; the fork-checklist gates before every push. Report
per the Worker rules.
```

### 232k-C — soil particles in C++ (opus · high; a different agent from 232k-V)

```text
<header>  ADR sections: §4 (W1–W17, especially W9, W12, W14, W17), §5.2 step 6, §5.3 and its
amendments, §6, §8 (the meter row, the P5 refusal row), §9 phase 5, §11, and the amendments
"232j" and "232k kickoff". Your PR: #<pr>.
Branch: work detached on origin/wp/232k-mpm-soil-models; commit; `git pull --rebase origin
wp/232k-mpm-soil-models`; `git push origin HEAD:wp/232k-mpm-soil-models`.
Read the "232k" section of .claude/skills/ladruno-mpm/references/wp-briefs.md first: the facts,
K1–K9 and the gate table bind you. Never widen a band.
Independence: read tests/mpm_interface.md, tests/mpm_implicit_interface.md,
tests/mpm_bspline_interface.md, tests/mpm_antilocking_interface.md and the cases files. Never read
tests/mpm_reference.py, tests/mpm_implicit_reference.py or tests/mpm_soil_reference.py.
Context: #1095, #1096, #1098, #1103, #1106 and #1107 touch code near yours. Do not fix them unless
your change rewrites that code; never regress them.

ROUND 1 — stop and report after it (the orchestrator freezes it):
  tests/mpm_soil_interface.md, at most ~3 pages, deciding with code anchors quoted:
  - K1: the admission rows for both materials; each option admitted (with its test) or refused
    (with its message); IMPL-EX's `ops_Dt` at particles in each mode.
  - K2: the thread rule, and what non-allowlisted materials get today.
  - K3: whether HS `-sig0`/`-sig03` is admitted, and why.
  - K4: the stage-flip mechanism for both modes; why no particle can miss SANISAND's
    Elastic2Plastic reset or HS's 0→1 seeding; the process-static warning.
  - K6: the deterministic refusal trigger per material, and the abort path you exercise.
  - K7: the SOE refusal list.
  - K8: whether the output must change.
  - The K0 meter: v and the step for each mode, the inertia check (at most 1 % of ρ g l0 in mode
    F), the domain height that allows at least 4 cells of translation, the run cost per material
    from the M0b ratio, and how the harness queries particle σ and SANISAND's e and α.
  - The cases file tests/mpm_soil_cases.txt (232k-V adds the reference blocks; you add the MPM
    blocks) and the test files, with each one's Zone-A cost.
ROUND 2 — on the orchestrator's "frozen":
  Implement K1–K4 and the meter harness first, batching edits; ask the orchestrator for each
  build (all five targets) and each mutant build; never start build.bat yourself. The K0 meter is
  the first soil MPM number: report it and stop for the orchestrator's verdict. Then S1–S10.
  Tests tests/test_ladruno_mpm_soil_*.py (zone_a within its time budget). The full K0, S2 and S3
  runs are the warrant run, which the orchestrator launches on esmeralda. Results JSON under
  tests/data/mpm/soil/ carry ladrunoBuild, the counts and per-row values.
  Classic Tcl for S5 and S10. Registrations per fork-checklist.md: the ledger fragments
  WP-232k-mpm-soil-models.md (implementations; vanilla rows and `// Ladruno` markers for any
  upstream file), a banner line, the guide (Ladruno_implementation/232_ladruno_mpm_guide.md: the
  soil rows, staging, the floors as the owner's, cost, serial), stamps.
  Draft the ADR amendment text in your report; the orchestrator writes it.
Done when: every gated row is green (or a reported miss with diagnosis); the reported rows are in
Ladruno_implementation/232k_mpm_soil_models.md; the ci/ gates are clean; ladrunoBuild matches your
head on the exe and the pyd. Report per the Worker rules.
```
