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
| non-associated mild-cone footing (φ_ps = 27.47°, plane strain, exact q₀·N_q = 138.9 kPa) | the peak / plateau at the finest of three resolutions (cells across B × particles per cell) within [0.90, 1.10] × 138.9 kPa; all three reported. FE's own monotone sequence 1.0020 / 0.9517 / 0.9260 (ADR-79) is context, not the gate **Outcome (2026-11-15):** unreachable at r3 even by the fork's FE (B-bar 0.882x, locked 1.114x); recorded as a strict xfail, and the owner chose to redesign the gate before phase 4 (ADR 232, amendment after the §8 table). |
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
