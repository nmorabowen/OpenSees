---
title: "ADR 90 orchestration plan — agents, effort, PR sequence, and the orchestrator's protocol"
project: Ladruno
type: execution plan (companion to the planning brief)
status: "KICKED OFF 2026-09-04 — owner approved §1 table, draft PR-1/PR-2, CP1 ownership; WP-A and WP-B launched. Number 90 allocated (88 = PR #778, 89 = proposed Track T)"
owner: nmora
orchestrator: "Claude Fable 5.1 (this session's model) — coordinates; does not merge"
related:
  - "[[_adr90_regularization_planning_brief]]"   # the WHAT — read first
  - "[[87_ladruno_depth_with_width_adr]]"         # warrant package + branch/PR rules in force
  - "[[testbed/00_canonical_testbed]]"
tags: [planning, orchestration, adr, regularization, duvaut-lions, tims]
updated: 2026-09-04
---

# ADR 90 orchestration plan

> [!abstract] One line
> Five work packages, five PRs, one conditional. The critical path is **A → C → D**
> (oracle decides the architecture → wrapper → acceptance case). Compute, not code, is the
> long pole: two full builds in WP-C and hours of slow-tier runs in WP-D. The orchestrator
> owns builds, reviews, PR state and checkpoints; agents own code and never end a turn with a
> build running.

## 0. What is deliberately different from the usual way

| Usual | This plan | Why |
|---|---|---|
| Code first, gate later | **A0 oracle decides D2 before any C++** (two-track vs true DL) | Brief R1 (BLOCKING); ADR-59 died on an unverified seam claim |
| ADR written once, reviewed by merge | **ADR text gets a 3-lens adversarial pass before WP-C opens** | Same convention that saved #709 (3 lenses) and sank ADR-59 (5 critics) — cheaper before code |
| Agent launches build, session moves on | **Agents never end a turn with a build running**; builds are `Monitor`-ed to the artifact | Memory: an agent is *never* re-woken by its own build finishing |
| Micro-PRs, auto-merge | ADR-87: `wp/<n>-<slug>`, **draft day one, one merge event, merge commit, owner merges** | Protection is live; drafts skip the 40-min Zone-A build until ready |
| Whole sweep in-session | **Slow legs run detached (nightly Zone-B / background) with wakeups**; sessions stop at checkpoints | Memory: sessions ran too long; "continue" = one phase |
| Fork decides alone | **Two formal TIMs handoffs** — spec (§5 of the brief) at kickoff, *verified* after WP-D | Vault 65 D5: the acceptance case is co-owned |

## 1. Agent and effort selection

Model roles (the Agent tool's `model` field; effort = `/code-review` level):

| Role | Model | Rationale |
|---|---|---|
| **Orchestrator** — decisions, ADR final synthesis, D2 call, checkpoints, PR state | **Fable** (this session) | Holds the only context that spans the vault, the sweeps and the brief |
| **Constitutive oracle (A0)** and **C++ wrapper (WP-C)** | **Opus** | Correctness-critical numerics (blend, tangent, `sendSelf` protocol, static-buffer traps); precedent: PR-2/PR-4 of the TIMs triage were Opus |
| **Adversarial critics** on the ADR text (3 lenses) and the out-of-family verdict on WP-C | **Opus** ×3, background, independent prompts | Critics must not share the author's context |
| **Acceptance decks (WP-D)** design + width metric | **Opus** | Physics interpretation (note-82's checkerboard lesson) |
| **Prerequisite fixes (WP-B)**, ledger/manifest/guide wiring, stale-doc fixes, TIMs handoff doc (WP-E) | **Sonnet** | Well-specified, pattern-following work; precedent: PR-3 of the triage |
| **Read-only sweeps** during any WP | **Explore** (Sonnet), "very thorough" | Keeps the orchestrator's context lean |
| Haiku | not used | Nothing here is mechanical enough to be safe |

Review effort per PR: WP-A `medium` (python + docs) → WP-B `medium` → **WP-C `high` + out-of-family verdict agent** → **WP-D `high` + a physics-lens Opus critic on the interpretation** → WP-E `low`.

## 2. Work packages and PRs

Dependency graph: `A ∥ B` → `C` (needs A's D2 decision; B merged or rebased in) → `D` → `E`.
`F` only if A0 says two-track is inadequate.

### WP-A — ADR 90 + numpy oracle · `wp/90a-adr-oracle` · PR-1

| | |
|---|---|
| Content | `90_ladruno_viscoplastic_regularization_adr.md` (from the brief + A0 result, alternatives table §3, claims §4, D-list §7); `tests/_testbed/duvaut_lions_ref.py` (two-track DL over J2 + DP point models; **A0 1-D bar** N = 20/40/80/160, imperfect centre element, two-track vs true DL); `tests/test_duvaut_lions_oracle.py` (PV1–PV6 clones + A0 width gates, seconds, Zone-A); the two `_adr90_*` docs; `LEDGER_implementations` row; `README` index line |
| Agent | **Opus** general-purpose for the oracle; ADR draft by Opus from the brief, **final synthesis by the orchestrator**; 3 Opus critics (constitutive / FEM-numerics / strategy) on the ADR draft |
| Build | **none** (no C++) |
| Gates | oracle green; A0: width ∝ h at τ = 0, converges at τ > 0; two-track-vs-true-DL difference quantified on a hardening law; critics' BLOCKING findings folded |
| Wall time | oracle ~1 day; ADR + review ~1 day; Zone-A on ready ~40 min (unavoidable, runs once) |
| **Checkpoint 1** | **D2 decision** (generic two-track ✓ / needs true DL → WP-F) — the owner decides on the A0 numbers; WP-C does not open before it |
| Merge | merge commit, owner |

### WP-B — SANISAND plane-strain prerequisites · `wp/90b-sanisand-ps-prereqs` · PR-2

| | |
|---|---|
| Content | cover `LadrunoSANISANDPlaneStrain::getCopy(void)`; fix `ManzariDafaliasPlaneStrain` null-ctor classTag (`LEDGER_quirks.md:4826`) with a broker round-trip test; **F8 stale docs** (`LadrunoConcrete3D_guide.md`, ADR-31 line 91); vanilla-file ledger rows for anything under `UWmaterials/` |
| Agent | **Sonnet** general-purpose |
| Build | one (`OpenSees OpenSeesPy`), in its own worktree |
| Gates | Zone-A; `test_ladruno_sanisand.py` + `test_fspm_over_manzari_family.py` green; `check_manifest` |
| Wall time | code ½ day; build ~40 min; Zone-A ~40 min |
| Merge | small fix → squash allowed |

Runs **in parallel with WP-A** from kickoff; its build is the background task while A0 is written.

### WP-C — the wrapper · `wp/90c-duvaut-lions-wrapper` · PR-3 (family PR, full warrant)

| | |
|---|---|
| Content | `SRC/material/nD/LadrunoDuvautLions.{h,cpp}` (ND tag **33022**, modelled on `StagedStrainNDMaterial` 33014): adopting ctor, forwarding, `getCopy(type)` only, copy-not-alias, `setParameter` forward on tag miss, `ops_Dt` with `dt ≤ 0 ⇒ inviscid`, tangent blend, τ as `Parameter`, responses `overstress/beta/dt`, non-uniform-Δt diagnostic; registration ×3 + CMake; Tcl + Python; `classTags.h` with the Ladruno comment; **`-DLADRUNO_MUTATE_DUVAUTLIONS`**; `ledger/ladruno_duvaut_lions.yml`; guide stub; banner line; g++ byte check vs the oracle fixture; Zone-A tests: byte gate (material point + `LadrunoBrick`), PV3 overstress, non-tautology, stage-forwarding over `LadrunoSANISAND`, database round-trip, MP wire (`test_fspm…` twin), mutation gate |
| Agent | **Opus** general-purpose — **one agent for the whole WP**, continued via `SendMessage` across review rounds so its context survives |
| Build | **two**: normal + mutation flag (the gate must turn the suite red). Plan ~2 h of build time; batch *all* C++ before the first build |
| Gates | ADR-87 warrant complete; `/code-review high`; out-of-family verdict (`reviews/handoff_adr90.md` → `reviews/adr90_verdict.md`); Zone-A green on the ready head |
| Wall time | code 1–2 days; builds 2 h; review rounds 1–2 days |
| **Checkpoint 2** | warrant complete → flip to ready → Zone-A → **owner merges** |

### WP-D — acceptance case · `wp/90d-acceptance-case` · PR-4

| | |
|---|---|
| Content | Legs **A1** (DP `quad PlaneStrain`), **A2** (`LadrunoSANISANDPlaneStrain`), **A3** (`LadrunoSANISAND` on `LadrunoBrick -formulation bbar`, one element thick), **A4** (R3 Prandtl regression at the declared De); G5-style positive + negative pairing; De × {½, 1, 2} sweep with matched-pair collapse; width metric (OQ5) settled on A0's analytic bar first; capacity three-clause rule; `ladrunoBuild()` stamped; **one leg in Zone-A slow tier, the sweep in Zone-B nightly** |
| Agent | **Opus** general-purpose for deck + metric; runs are launched by the **orchestrator** (background, wakeups), not held open by the agent |
| Build | none (tests only) unless WP-F fired — runs on the **merged** PR-3 binary |
| Gates | C4, C5, C6 of the brief; A4 unchanged within the De family; physics-lens critic on the interpretation |
| Wall time | decks 2 days; **runs: hours per leg** (R3 alone is 61 min; A3 × 3 meshes × {τ=0, τ>0} × 3 De ≈ a day of compute) |
| **Checkpoint 3** | the numbers, interpreted, before anything is written to the vault |

### WP-E — TIMs handoff · `wp/90e-tims-handoff` · PR-5

`Ladruno_implementation/90_ladruno_duvaut_lions_tims_report.md` in the shape of the ADR-86
TIMs report (findings as *what it does / measured evidence / what it means for you*), plus the
apeGmsh emitter note (`-tau` flag, De reporting fields, parser ordering rule). **Sonnet**. Docs
only. The vault note that records the acceptance case as *verified* is TIMs' to write; the fork
hands over the draft.

### WP-F — conditional · `wp/90f-sanisand-tau` · PR-3b

Only if Checkpoint 1 says two-track is inadequate: a material-specific `-tau` inside
`LadrunoSANISAND` (true DL on $\sigma$, $\alpha$, $z$, using the `protected` base members),
cross-validated against the wrapper on the A0 law. **Opus**. One build. Same warrant as WP-C.

## 3. Orchestrator protocol (what the orchestrator does, in order, every time)

**Kickoff (once)**
1. `git worktree list` and `gh pr list` — no duplicate lane (memory: two agents worked #753/#754 at once).
2. Confirm with the owner: opening draft PRs is a publishing action — get an explicit yes per PR.
3. Send TIMs the brief's §5 as the OQ5 draft; ask for OQ2's numbers. Do not wait on them for WP-A/B/C.

**Per WP**
4. One worktree per live WP (Agent `isolation: worktree`, or `git worktree add` on `wp/<n>-<slug>` from fresh `ladruno`). **Never build in the shared checkout.**
5. Spawn the WP agent with: the brief, the WP row above, and the standing rules (§4). Background by default; the orchestrator continues.
6. **Builds are the orchestrator's** or are `Monitor`-ed by the agent to the artifact — never "launched and forgotten". Verify `dist/bin/opensees.pyd` mtime and `ops.ladrunoBuild()` before any test claim; a stale `tests/opensees.pyd` shadows a fresh build (BUILD_GOTCHAS §4b).
7. Open the **draft PR day one** with the agent's first commit; the PR description is the acceptance record and grows with the WP.
8. Review loop: `/code-review <effort>` → findings → `SendMessage` to the *same* agent (context intact) → re-review. Independent critics as separate background agents.
9. Warrant check against ADR-87 D9 before flipping to ready; then wait for Zone-A on the ready head. **The owner merges. Agents do not.**

**Checkpoints (the orchestrator stops and reports; "continue" means the next phase only)**
- **CP1** after A0: D2 decision. — **CP2** WP-C warrant complete: merge decision. — **CP3** WP-D numbers: interpretation before the vault hears anything.

**Long waits** — builds (40 min), Zone-A (40 min), slow legs (hours): use `Monitor` on the log/artifact or a scheduled wakeup, not polling and not a blocked session. Windows traps on record: `Start-Process` builds die on session restart (launch via WMI + Monitor the log); Git Bash `cmd /c build.bat` is a silent no-op; MSVC C1001 is transient — rerun.

## 4. Standing rules handed to every agent (verbatim in the prompt)

1. Read `_adr90_regularization_planning_brief.md` §1–§2 before touching anything; do not re-open §1.
2. Batch **all** C++ before the first build; edit with Write/Edit, never heredocs.
3. If you launch a build, `Monitor` it to `dist/bin/opensees.pyd` and **do not end your turn while it runs**; report the pyd mtime and `ops.ladrunoBuild()` with every test result.
4. Build `OpenSees OpenSeesPy` (the exe alone leaves the pyd stale).
5. Every number you report names build hash, driver, output path, date.
6. Do not open, flip, or merge PRs; do not write to the vault; do not touch `C:\Program Files\Ladruno\`. Report back.
7. `stdBrick` swallows material return codes — use `LadrunoBrick` for anything that must fail loudly.
8. `ManzariDafalias` family: copy returned vectors (static buffers); `getCopy(type)` only; forward `setParameter` on tag miss; `revertToLastCommit` is a stub.
9. Ledgers are part of the PR: `LEDGER_implementations` row, `LEDGER_vanilla_files` for any upstream file (with a `// Ladruno` comment), `LEDGER_quirks` for anything learned.
10. Commit trailer `Co-Authored-By: Claude <model> <noreply@anthropic.com>` per the harness; no `--no-verify`.

## 5. Calendar (order-of-magnitude, serial compute)

| Week | Milestone |
|---|---|
| 1 | WP-A oracle + A0 (days 1–2) ∥ WP-B (days 1–2); ADR draft + 3-lens review (days 3–4); **CP1** (day 4); PR-1/PR-2 ready → Zone-A → merge (day 5) |
| 2 | WP-C code (days 1–2), two builds + Zone-A, review rounds, verdict; **CP2** (day 5) |
| 3 | WP-D decks (days 1–2); runs detached (days 2–5); **CP3**; WP-E draft |
| 4 | buffer: WP-F if fired; TIMs verification; PR-4/PR-5 merge |

~3 weeks nominal, 4 with WP-F. The single largest schedule risk is OQ2 arriving late (band-width
targets); mitigation: WP-D runs with provisional bands and re-fits, curves stored whole (the
vault 65 D6 pattern).

## 6. Go/no-go items for the owner before kickoff

1. Approve the model/effort table (§1) — or reassign.
2. Approve opening **draft PR-1 and PR-2** on `nmorabowen/OpenSees` at kickoff (publishing action).
3. Confirm the orchestrator may send the §5 acceptance-case draft to TIMs (external message) — or you send it.
4. Confirm CP1 is a decision **you** make, not the orchestrator.
