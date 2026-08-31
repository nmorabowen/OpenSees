---
title: Depth with width — verification, workflow, and warrant program
project: Ladruno
status: accepted
priority: high
owner: nmora
tags:
  - adr
  - governance
  - verification
  - process
aliases:
  - ADR-87
  - depth with width
  - warrant program
updated: 2026-08-31
---

# ADR-87 — Depth with width

> Number **87** on purpose: `ladruno` already carries ADR-86 (SANISAND,
> PRs 768–772). Do not reuse 86. An earlier, unrelated draft briefly held
> this number ("warrant the kernel", never committed, withdrawn
> 2026-08-30); this document replaces it and inherits nothing from it
> except the evidence it verified.
>
> This is a **governance + verification** ADR. It adds no class tag. It
> defines a program: how the fork keeps growing (width) while making
> every claim it ships independently checkable (depth).
>
> **Accepted 2026-08-30** by PR #773, which did the WP-1 work rather than
> merely adding this file (branch protection live, banner cut to
> shipped-only ≤70 cols, LadrunoUP collapsed to one truthful row). The
> acceptance bar it set for itself is in *How* below.

## What

A standing operating rule plus four work tracks:

**The rule — depth is the admission ticket for width.** New families
keep coming. Every new family PR arrives with a verification manifest,
a mutation gate, and a user guide, or it does not merge. The existing
surface is retrofitted through the same bar, family by family. There is
**no freeze**.

**Not in scope:** deleting `Ladruno*` objects, halting apeGmsh
integration, or rewriting history.

## Why

The 2026 agentic campaign produced a real kernel: continuum, explicit
integration, contact/ties, porous media, solvers, and a set of
upstream silent-failure repairs — much of it anchored to external
oracles (Simon 1984 closed form incl. 4 documented paper errata, ZS84,
ZCB80, Terzaghi, two-leg equivalences vs stock) and with a record that
retracts its own claims when measurement refutes them (CMS §32,
ADR-79 vs ADR-78 §1, ADR-80 G3).

But the *warrant* infrastructure has not kept pace with the surface
(verified against the tree 2026-08-30):

- `LEDGER_implementations.md` is 484,515 characters; 3 of 4 LadrunoUP
  rows carry "Element still not built" on the same line as "SHIPPED
  (P0–P4 complete)". `merge=union` duplicates in-place edits
  (`[[../Ladruno_internal/WORKFLOW_GOTCHAS]]` §2a).
- The banner is 67 lines; 44 exceed 70 columns; the contact line is
  1,362 characters.
- The merge path cannot refuse: auto-merge lands in minutes, a
  non-compiling PR has merged (#175), `gh pr merge --auto` has
  squashed before Zone-A registered, and a Tcl suite once reported
  green as "OK (0 checks passed)" (since fixed in #653 —
  `check_tcl_results.py` fails on `total == 0`).
- Verification citations are dominated by same-model-family panels.
  They find real bugs; they are not an independent laboratory.

The gap between what we claim (banner, "shipped") and what a third
party could check is the project's main liability. Closing it is work;
this ADR schedules that work without stopping the feature program.

## Decisions

### D0 — No freeze. Warrant-first.

A new element / material / handler / integrator family is welcome iff
its opening PR carries: (a) a verification manifest (D1), (b) a
mutation gate (D2), (c) a user guide stub, and (d) both-interpreter
wiring or an explicit deferral note. "While we're here" modes ride the
same rule: they extend the manifest or they wait.

### D1 — Verification manifest per feature.

One machine-readable YAML per family under
`Ladruno_implementation/ledger/<feature>.yml`:

```yaml
feature: LadrunoUP
kind: element
class_tag: 33017
status: shipped          # shipped | draft | frozen | failed
owner: nmora             # a human name, not an agent
oracles:                 # at least one non-self-referential entry
  - type: closed_form    # closed_form | mms | cross_code | two_leg
    ref: "Simon 1984, B5 dossier"
    test: tests/test_ladruno_up_element.py
mutation_gate: tests/test_mutation_up.py   # D2
guide: Ladruno_implementation/LadrunoUP_guide.md
interpreters: [tcl, python]
parallel: serial-only    # or: mp-verified, sp-verified
prs: [547, 557, 677]
```

The human-readable ledger table is **generated** from these files
(script in `Ladruno_scripts/`). Session narrative lives in the ADR
implementation logs, not the ledger. This dissolves both the 484k
monolith and the `merge=union` duplicate-row footgun: one file per
feature, one status field, no in-place table edits.

### D2 — Mutation gates: green must be falsifiable.

For every family called `shipped`, CI proves the tests test the
physics: a mutation build (`-DLADRUNO_MUTATE_<FAMILY>` — stress kernel
returns zero, or tangent returns identity) must turn that family's
suite **red**. A suite that stays green under physics deletion is a
defect of the same severity as a wrong stress.

Mechanics: one CMake option per family guarded to touch only fork
files; a CI lane builds N mutants and asserts `pytest -m <family>`
fails for each. Rolled out top-down: continuum, UP, contact, explicit
integrators, SANISAND first.

### D3 — CI that can refuse is the merge criterion.

- Zone-A (build + `pytest -m zone_a`) becomes a **required**
  branch-protection check on `ladruno`. A green classTag/manifest
  check is not a merge.
- Draft PRs do not merge. `gh pr merge --auto` is not used on
  family-level PRs (the stale-PR "no checks" trap,
  `WORKFLOW_GOTCHAS` §6).
- Banner lint in CI: every line ≤70 columns, maps to a ledger entry
  with `status: shipped` and a passing mutation gate. The splash is an
  enforced claim, not a changelog.
- CI jobs that cannot run (self-hosted lanes with zero registered
  runners: Zone-B, cross-tier leak) are either moved to hosted runners
  or deleted along with the claims that cite them.

### D4 — Out-of-family review is the citable review.

Same-family agent panels remain bug-finding tools; they are logged
when they find something and cited as nothing else. For family-level
PRs and ADR acceptances, the workflow produces a review handoff pack
(`reviews/handoff_<adr>.md`: the ADR, key diffs, specific questions)
that the owner runs through an out-of-family model (SuperGrok today);
the verdict is committed as `reviews/<adr>_verdict.md` beside the ADR.
Precedent: apeGmsh `internal_docs/handoff_adr0094_grok.md`.

### D5 — Depth artifacts: MMS, cross-code, derivations.

Retrofit program over the existing kernel, in priority order:

1. **MMS + convergence-order tables** for the continuum family and
   LadrunoUP: manufactured solutions, observed vs theoretical order,
   one log-log figure per element lane. The canonical depth artifact.
2. **Cross-code dossier**: 3–5 canonical problems (Terzaghi column,
   Prandtl–Reissner bearing, contact patch, wall pushover) run as
   stock-OpenSees / fork / Abaqus-theory reference values; one delta
   table. Instruments already exist (abaqus-theory suite, two-leg
   habit).
3. **Derivation notes** per family: weak form, return map, first
   integrator step — short, committed beside the guide. The durable
   form of "a named human can derive it"; seed chapters of a theory
   manual.

### D6 — Show both axes: the matrix, the gallery, the release.

- The fork README's centerpiece becomes a **feature × verification
  matrix**: rows = families; columns = oracle, MMS, cross-code,
  mutation gate, guide, both interpreters, parallel. Generated from
  the D1 manifests. Width is the row count; depth is the fill.
- A worked-example gallery: one apeGmsh → fork → validated-result
  example per family (the LadrunoUP / PorousOverlay guides prove the
  pattern).
- Tagged releases: `ladruno-v1.x` = installer + generated matrix +
  verification report. A rolling HEAD cannot be warranted; a tagged
  snapshot can.
- One upstream slice: 2–3 vanilla silent-failure fixes PRed to
  `OpenSees/OpenSees` (no 33000 band required). External acceptance is
  validation no internal process can manufacture.

### D7 — Nightly observatory.

3–5 canonical decks run nightly on hosted CI: accuracy vs pinned
reference plus wall-clock, tracked over time. Catches physics
regressions and silent performance decay; produces the longitudinal
evidence that the kernel is stable.

### D8 — Close the named residuals.

Standing debts, each raising warrant when closed: the LadrunoUP ledger
contradiction (D1 collapses it); ADR-80 `getTangForce`; ADR-84's owed
Cerro Lindo M3 re-run; SANISAND through the D0 bar; CMS finishes P4 or
takes an honest `frozen` status; finite-sliding `-reemit`
finish-or-freeze. Generalize ADR-84's `strict_convergence` pattern
(exhaustion-accept, f-decreasing exits, stress-scaled tolerances) as
the fork-wide convention for constitutive returns.

### D9 — One merge per work package; continuous work inside it.

The unit of merge to `ladruno` is the **work package / ADR phase** —
larger than the 2026-campaign micro-PR (758 merged PRs no human read),
smaller than a mega-branch that drifts for months:

- One branch per WP from fresh `ladruno`, named `wp/<n>-<slug>`.
- A **draft PR opens on day one**: the fast static gates (classTag +
  manifest + recorder oracle) run on every push; the full Zone-A build
  deliberately skips drafts (a ~40-min build per push would be waste)
  and fires when the PR is flipped to ready (`ready_for_review` is
  wired in the workflow trigger) — or on demand mid-campaign via
  `workflow_dispatch`. Drafts cannot merge, so the one deliberate
  merge event always carries a real Zone-A run on its head commit.
  Agents commit continuously to the WP branch in a dedicated worktree;
  no intermediate PRs, no `--auto`.
- Ready = the warrant package is complete: manifest + mutation gate +
  guide (D0), Zone-A green, and (family-level) the committed
  out-of-family verdict (D4). The PR description is the acceptance
  record.
- **Merge commits for WP PRs; squash only for small fixes.** Squash
  breaks `git branch --merged` ancestry detection — the mechanism that
  produced the current zombie-branch pile — so WP history merges
  intact.
- Hotspots: classTag reservations are claimed in the manifest at WP
  kickoff; per-feature ledger files (WP-5) remove the ledger/banner
  conflict class.
- Small bugfix/docs PRs remain legal but wait for Zone-A like
  everything else.

### D10 — Branch hygiene: sweep, triage, graveyard.

Census 2026-08-30: 192 local branches (141 ancestry-merged), 685
remote (326 ancestry-merged), 20 worktrees — and squash merges mean
ancestry *undercounts* the dead.

- Enable GitHub "automatically delete head branches" (fixes the
  future), then delete all ancestry-merged branches, local and remote.
- For the remainder, script the squash-zombie check: a branch whose
  head was consumed by a merged PR (`gh pr list --head <branch>
  --state merged`) is dead — delete.
- True stragglers (unique commits, no merged PR) are the dangerous
  class (ADR-79 §9 stranded five days on a merged branch): each gets a
  row in `Ladruno_internal/branch_graveyard.md` — branch, SHA,
  verdict (*rescued as PR* / *recorded and dropped*) — before
  deletion. The owner reviews the triage report before any unmerged
  branch dies.
- Worktrees: prune those on merged/dead branches; standing rule is one
  worktree per live WP.

## Work packages (planning skeleton)

| WP | Track | Scope | Deliverable |
|---|---|---|---|
| WP-0 | D10 | Branch sweep: auto-delete setting, ancestry-merged deletions, squash-zombie script, straggler triage report | clean branch list + graveyard file |
| WP-1 | D3 | Branch protection: Zone-A required, no draft merges; banner cut to shipped-only ≤70 cols; fix LadrunoUP ledger rows | settings + one docs/banner PR |
| WP-2 | D2 | Mutation-gate framework + first gate (continuum) | CMake option pattern, CI lane, 1 red-proof |
| WP-3 | D2 | Mutation gates: UP, contact, explicit, SANISAND | 4 more gates |
| WP-4 | D1 | Manifest schema + generator + top-10 family manifests | `ledger/*.yml`, generated table |
| WP-5 | D1 | Ledger split: retire the monolith, per-feature files | restructured `Ladruno_implementation/` |
| WP-6 | D5.1 | MMS: continuum family | convergence report + figures |
| WP-7 | D5.1 | MMS: LadrunoUP | convergence report + figures |
| WP-8 | D5.2 | Cross-code dossier (3–5 problems) | delta table + decks |
| WP-9 | D4 | Review handoff pack generator + first Grok cycle | `reviews/` convention |
| WP-10 | D7 | Nightly observatory (decks + hosted CI + trend plots) | workflow + dashboard |
| WP-11 | D6 | Matrix generator + README + example gallery skeleton | generated matrix |
| WP-12 | D6 | `ladruno-v1.0` tagged release | installer + report |
| WP-13 | D6 | Upstream slice (2–3 vanilla fixes) | PRs to OpenSees/OpenSees |
| WP-14 | D8 | Residual burn-down (UP rows via WP-1; ADR-80/84/CMS/`-reemit`) | per-item PRs |
| WP-15 | D5.3 | Derivation notes, top-5 families | theory chapters |

Sequencing: WP-1 first (a day; flips the repo to refuse-by-default).
Then WP-2..5 (the warrant infrastructure), then WP-6..9 in parallel
with normal feature work under D0. WP-10..15 as capacity allows.
Agent assignment and effort level per WP are decided at kickoff and
recorded in the implementation log below.

## Where

| Change | File |
|---|---|
| This decision | `Ladruno_implementation/87_ladruno_depth_with_width_adr.md` |
| Manifests + generated table | `Ladruno_implementation/ledger/` (new) |
| Traps promoted to policy | `[[../Ladruno_internal/WORKFLOW_GOTCHAS]]` |
| Banner source | `Ladruno_scripts/banner_features.txt` + new lint |
| CI | `.github/workflows/ladruno.yml` + repo branch-protection settings |
| Reviews | `reviews/` (new, repo root or `Ladruno_implementation/reviews/`) |

## How (acceptance of this ADR)

A PR that only adds this file is not acceptance. Acceptance is WP-1
done in the same PR (or an immediately following one referenced here):

1. Frontmatter `status: accepted`.
2. Zone-A required check live on `ladruno` (or documented as blocked
   with the owner's screenshot of the settings attempted).
3. Banner: shipped-only, ≤70 columns, contradiction-free.
4. LadrunoUP collapsed to one truthful status.
5. Ledger row for ADR-87 + README index line + WORKFLOW_GOTCHAS
   policy note (the standard build-control wiring).

## Risks

- **The manifest becomes ceremony.** Mitigation: D2 — a manifest
  without a red-proving mutation gate does not count as `shipped`.
- **Mutation builds bloat CI.** Mitigation: mutants build only the
  fork target and run only the family's marker; one lane, matrix
  strategy, cached toolchain.
- **Out-of-family review has no API** (SuperGrok is manual).
  Mitigation: the handoff pack makes the manual step 5 minutes; it is
  required only at family level, not per PR.
- **Width outruns the retrofit.** Mitigation: D0 binds new width to
  the same bar, so the unwarranted backlog can only shrink.
- **The owner does not accept.** Then nothing changes, and the gap
  keeps growing at the measured rate (banner +10 lines in the week of
  2026-08-23 → 30).

## Non-goals

- A freeze on new families (the withdrawn draft's D2; explicitly
  rejected here).
- A constitutive audit of every 33000-band object in one campaign
  (D5/D8 retrofit it incrementally).
- Renaming `-consistanttan` (ABI; alias later).

## Implementation log

- 2026-08-30 — drafted, replacing the withdrawn "warrant the kernel"
  draft after owner review (freeze rejected; warrant-first adopted).
  Evidence base re-verified against the tree the same day: ledger
  484,515 chars; 3/4 LadrunoUP rows self-contradictory; banner 67
  lines / 44 over 70 cols / max 1,362 chars; Zone-A workflow present
  but not required. Status `proposed`; WP table awaits agent/effort
  assignment.
- 2026-08-30 — D9 (WP-sized merges, day-one draft PR, merge commits)
  and D10 (branch sweep + graveyard) added at the owner's direction,
  replacing the micro-PR/auto-merge flow. Branch census the same day:
  192 local / 685 remote branches, 141 / 326 ancestry-merged, 20
  worktrees; squash merges undercount the dead. WP-0 added.
- 2026-08-30 — **WP-0 launched + WP-1 executed** (owner directive:
  "launch WP-0 and WP-1"). Settings applied live via API:
  `delete_branch_on_merge=true`; branch protection on `ladruno` with
  required checks `Zone-A (Ubuntu)` + `classTag + manifest gates`,
  `enforce_admins=true` (self-hosted lanes deliberately NOT required —
  zero registered runners would deadlock every merge; D3's
  fix-or-delete still owed). `ladruno` previously had **no protection
  at all** (API 404). Banner rebuilt: 57 features, 0 lines over 70
  cols (was 44 over; contact line 1,362 → 60 chars; duplicate
  Quad/CST line deduped; the unfinished finite-sliding `-reemit`
  claim removed from the splash). LadrunoUP collapsed 4 rows → 1
  (stale P0/P2/P3 `merge=union` snapshots deleted; "element not
  built" tail replaced by a dated dedup note). Protocol recorded in
  project `CLAUDE.md` (PRs and branches) + `WORKFLOW_GOTCHAS` policy
  header. Status flipped to `accepted`; the WP-1 PR is the accepting
  PR, run under D9 itself (`wp/1-accept-adr87`, day-one draft).
  **WP-1 MERGED as [#773](https://github.com/nmorabowen/OpenSees/pull/773)**
  — Zone-A green (16m18s), landed as a true merge commit
  (`68ab8baab`, two parents), head branch auto-deleted. D9's first
  live exercise behaved as designed.
- 2026-08-30 — **WP-0 COMPLETE** (D10 sweep). 769 refs deleted:
  132 local + 337 remote ancestry-merged, then 293 remote + 37 local
  **squash-zombies** (proof = a merged PR used the branch as head —
  the necessary phase, since squash merges make `git branch --merged`
  under-report). `origin/*` 656 → 26, local 193 → 24; `jaabell` (29)
  untouched; every exclusion verified intact afterwards. 13 stragglers
  triaged in `[[../Ladruno_internal/branch_graveyard]]`, **none
  deleted** (D10 reserves that for owner sign-off): 5 `up/*` upstream-
  campaign branches keep, 8 drop. Both subagent "rescue"
  recommendations were **overturned on verification** — the SOE branch
  is the losing lane of the #753/#754 duplicate race (the landed test
  is +437 lines richer), and RebarBuckling shipped and has since
  gained v2; the 46-commit stiff-soil branch is safe to drop only
  because those commits live on `jaabell` AND the model already ships
  (`StiffSoil*.h`). Verify subagent verdicts before acting on them.
  Also fixed: the main checkout's local `ladruno` ref was 3 commits
  stale, which silently under-detects merged branches.
- 2026-08-31 — **WP-0 close-out**: the 8 `drop`-verdict stragglers were
  deleted after owner sign-off (SHAs recorded in
  `[[../Ladruno_internal/branch_graveyard]]` first, and re-verified
  against the remote immediately before deletion). `origin/*` is now 18
  real branches: `ladruno`, the 5 `up/*` campaign branches, and 12 held
  by live worktrees.
- 2026-08-31 — **WP-2: mutation-gate framework SHIPPED + first
  (CONTINUUM) gate live.** Four design decisions, each made against a
  specific failure mode:
  (a) **Call sites in the Element API accessors**, not the internals —
  `getTangentStiff` / `getInitialStiff` / `getResistingForce` /
  `getResistingForceIncInertia` are the whole of what an analysis can
  observe about an element, so 11 call sites cover LadrunoBrick's SIX
  assembly paths (std/bbar, URI, physical-hourglass, SSP, EAS, finite,
  hypo) plus LadrunoBrick20, and a future formulation is covered the
  day it is written.
  (b) **Always compiled, never `#ifdef`-ed out** — mode NONE folds every
  branch to `if (0)`, so the default build is unaffected while the
  mutation code stays type-checked. Scaffolding that can rot into a
  no-op is worse than none: it reports a safety it is not delivering.
  (c) **A fail-loud identity verb.** `ladrunoMutation` (Python + Tcl +
  classic Tcl) reports "none" or e.g. `CONTINUUM=ZERO`, and the driver
  REFUSES to score a binary that disagrees with `--expect`. Without it
  a silently-failed mutant build would re-run the previous binary, every
  test would pass, and the gate would report the exact OPPOSITE of the
  truth — "this suite cannot detect deleted physics" — with no error
  anywhere. Same idiom as `ladrunoBuild` for stale-`.pyd`.
  (d) **The survivor list is the deliverable**, not the score: the tests
  that pass with the physics deleted are the actionable output. Modes:
  `ZERO` (the canonical return-0), `SCALE` x1.5 (the diagnostic one — a
  plausible-but-wrong answer, so a suite green under SCALE is a smoke
  test), `IDENT` (tangent only; most tests correctly survive, so read a
  low score as "no tangent coverage", not as a bug).
  Verified offline (no build needed): header compiles and its logic is
  right (ZERO nulls, IDENT builds a true identity, NONE is inert); the
  mutation-string builder emits exactly `none` / `CONTINUUM=ZERO` /
  `CONTINUUM=SCALE,UP=IDENT`; scoring excludes baseline-failures and
  both gate directions exit correctly; all four CMake paths behave,
  including **a typo'd family being FATAL rather than a silently
  unmutated build**; both patched elements pass `g++ -fsyntax-only`;
  static CI gates green.
  **Caught during WP-2 and worth keeping:** a family's test glob must
  select exactly the code that carries its call sites — the first
  CONTINUUM glob pulled in LadrunoBrick20 tests before that element was
  gated, which would have reported unmutated-element tests as "tests
  that ignore the physics". Fixed by gating Brick20 too; the driver now
  also refuses to run an ungated family so a meaningless 0.0 can never
  be read as a damning result.
  **NOT yet measured:** the score itself needs the three-build CI lane,
  and `workflow_dispatch` only offers workflows that already live on the
  default branch — so the first real number arrives after this merges.
  The SCALE floor is deliberately left ungated until it can be set from
  a measured baseline instead of a guess.
- 2026-08-31 — **WP-2b: the first mutation-lane run FAILED, and it earned
  its keep.** Run 33365505674 died in the baseline step and exposed two
  defects, the second serious:
  (i) the driver probed and ran pytest from the repo ROOT while CI copies
  `opensees.so` into `tests/` — `ModuleNotFoundError` on a perfectly good
  build. Zone-A works only because it does `cd tests` first. Fixed with a
  `--module-dir` (default `tests/`) used as the cwd for BOTH the probe and
  pytest; test files now resolve as names relative to it.
  (ii) **the gate could not fail.** GitHub's default shell is `bash -e` —
  errexit but NOT pipefail — so `python ... | tee` returned tee's exit
  code and a FAILED gate would have reported success. The tell was inside
  the failed run itself: the SCALE score step reported `success` while its
  input file did not exist. Fixed with `set -o pipefail` and explicit
  `${PIPESTATUS[0]}` propagation; the SCALE step now separates "the mutant
  never ran" (skip) from "scoring itself broke" (fail).
  Defect (ii) is this ADR's own thesis reproduced INSIDE the gate meant to
  enforce it: a check that runs, reports green, and cannot say no. It
  survived review because the step *looked* like it gated. In Actions, a
  pipe into `tee` is a gate-killer unless pipefail is set — worth carrying
  into every future CI lane.
  Verified locally before re-running: the probe now finds a module in
  `--module-dir` and the identity-mismatch guard fires (synthetic module);
  `false | tee` measured to exit 0 under `bash -e` while `${PIPESTATUS[0]}`
  correctly captures 3.
