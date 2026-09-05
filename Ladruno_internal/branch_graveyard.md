---
title: Branch graveyard — what was swept and why
project: Ladruno
tags:
  - internal
  - workflow
  - git
---

# Branch graveyard

Per [[../Ladruno_implementation/87_ladruno_depth_with_width_adr|ADR-87]] D10:
no branch with unique commits is deleted silently. Every one gets a row here
first — branch, head SHA, and the evidence behind the verdict — so a deleted
ref can still be found and reasoned about later.

Ancestry-merged and squash-merged branches are **not** listed individually:
they carry no unique content by definition (a merged PR is the proof), and
listing 769 of them would defeat the purpose of a readable record. The counts
are in the sweep log below.

## Sweep 1 — 2026-08-30 (WP-0, ADR-87 D10)

Before: 193 local, 656 `origin/*`, 20 worktrees. After: 24 local, 26
`origin/*`. The `jaabell` remote (29 branches, upstream lineage) was out of
scope and untouched.

| Phase | Basis | Refs deleted |
|---|---|---|
| 1 | local ancestry-merged (`git branch -d`, safe mode) | 132 |
| 2 | remote ancestry-merged (`git push --delete`) | 337 |
| 3 | squash-zombies — a **merged PR** used the branch as head | 293 remote + 37 local |
| 4 | stragglers — no merged PR (triage only, no deletion) | 0 |

Why phase 3 exists: the fork squash-merged for most of 2026, so
`git branch --merged` under-reports the dead. A branch whose head was consumed
by a merged PR is dead regardless of what ancestry says.

Protection verified intact after the sweep: `ladruno`, the WP branch, and all
16 worktree-checked-out branches survived; zero `gh` lookup errors.

### Stragglers — verdicts (owner sign-off pending)

Thirteen branches had no merged PR. Each verdict below was **re-verified
against the tree** (subagent recommendations were not taken on trust; two were
overturned).

| Branch | Head | Verdict | Evidence |
|---|---|---|---|
| `claude/charming-joliot-b678fc` | `2017c48ff` | drop | Losing lane of the #753/#754 duplicate race ([[../Ladruno_internal/WORKFLOW_GOTCHAS]]). The landed `tests/test_printa_unsized_soe.py` is **richer** (+437 lines) and already skips absent-`system` rows with a more careful crash-vs-skip rule. *(Overturned from "rescue".)* |
| `guppi/ladruno-rebar-buckling-impl` | `363e74893` | drop | `LadrunoRebarBuckling.{h,cpp}` are in the tree and have since gained v2 cyclic re-straightening (`6db43fc09`). The branch holds only the superseded original. *(Overturned from "rescue".)* |
| `merged-Explicit-ASDP-master` | `0abe2d5cd` | drop | 46 unique commits incl. Hardening/Stiffening Soil work — but those commits are reachable from the **`jaabell` remote** (preserved upstream) and the model already ships here (`SRC/material/nD/ASDPlasticMaterial3D/**/StiffSoil*.h`). 1,736 commits behind. |
| `guppi/explicitbathe-sms` | `9c1c9b675` | drop | Its `ExplicitBatheSMS` 33009/33010 design was deliberately superseded by the ADR-52 W1-E2 6→1 collapse (#419); those tags are recorded as deprecated aliases in `classTags.h`. |
| `docs/apegmsh-fork-online` | `15a5f0c5c` | drop | Docs-only scoping notes, 683 commits behind; superseded by later apeGmsh integration docs. |
| `feature/explicit-dtcr-hardening` | `f67007755` | drop | Single commit, 1,468 behind; the explicit `dt_cr` machinery has been rebuilt since (ADR-36/38, overlay-aware advisory). |
| `fix/ladruno-ledger-conflict` | `351b1b614` | drop | Fixed conflict markers in a ledger that no longer has them; 1,404 behind. |
| `feature/nmb_compile` | `58f6b3915` | drop | Conan-2 build tweaks from 2025-04; superseded by `Ladruno_scripts/build.bat`. 2,153 behind. |
| `up/00-tennodetet-shp3d` | `b84972dda` | **keep** | Live upstream PR-campaign branch (fresh-branch port to `jaabell`). Not squash-merge candidates by design — see [[../Ladruno_implementation/upstream_pr_campaign]]. |
| `up/01-portability-crash-fixes` | `bdc1b2f52` | **keep** | Same campaign. |
| `up/02-quad-tri-rho-serialization` | `4e4780e82` | **keep** | Same campaign (the element-rho serialization fix). |
| `up/03-gmsh-hex20` | `b4939b3d0` | **keep** | Same campaign. |
| `up/04-domain-clearall-eq-leak` | `b5f2b398a` | **keep** | Same campaign. |

The eight `drop` rows are recorded but **not yet deleted** — D10 reserves that
for owner sign-off. Deleting them later needs no further investigation; this
table is the record.

### Operational findings

- The main checkout's local `ladruno` ref was 3 commits behind `origin/ladruno`
  and has been refreshed. A stale local ref silently under-detects merged
  branches — always sweep against `origin/ladruno`.
- All 20 worktrees point at commits already merged into `origin/ladruno`, so
  each is a `git worktree remove` candidate **once its owning session is done**
  (merged ≠ session finished — this session's own worktree is in that list).
  Not acted on.
- Going forward `delete_branch_on_merge` is on, so merged heads self-delete and
  a sweep this large should never be needed again.

## Sweep 2 — 2026-09-04 (worktree prune, ADR-87 D10)

Before: 20 registered worktrees, of which **17 were dead** (zero commits ahead of
`origin/ladruno`, 78–885 behind, twelve still carrying a 1–2 GB `build/`+`dist/`),
plus one unregistered orphan folder and the main checkout parked on a merged
branch with a stranded commit. After: 4 worktrees — the main checkout and three
live WP worktrees.

The census that found this (run per worktree): branch, `git rev-list
--left-right --count origin/ladruno...HEAD`, `git status --porcelain | wc -l`,
last commit date, `build/` + `dist/bin/opensees.pyd` present. "0 ahead" is the
safe-to-remove test; a dirty tree on a 0-ahead branch is the **stranded-work**
class and is rescued first.

### Rescued

| What | Where it sat | Rescued to | Notes |
|---|---|---|---|
| `a71cdb197` build(cmake): threaded-layer opt-in for the Linux PARDISO | pushed to `claude/pardiso-linux` **after** #781 had merged it (WORKFLOW_GOTCHAS §1) | `wp/rescue-pardiso-linux-threaded` = `e5593cdfc`, **draft PR #782** | `cherry-pick -x` onto fresh `ladruno`, clean |
| ADR-80 P3 tangent-predictor WIP: `TransformationFE`/`DOF_Group` getTangForce, `LadrunoLoadControl`, gates, `80d` verdict doc, test (20 files, +504) | uncommitted in `.claude/worktrees/sp-strengthening-9c0ddc` on a merged branch, 229 behind | `origin/rescue/adr80-p3-tangent-predictor-wip` = `5b2debaac` | unbuilt snapshot; rebase before use |
| LadrunoBrick/LadrunoQuad k-stab damage-scaling plasticity + `tests/test_ladrunoBrick_kstab_plasticity.py` (+304) | uncommitted in `pardisio-profiling-0a03b1`, 268 behind | `origin/rescue/brick-kstab-plasticity-wip` = `76594384d` | unbuilt snapshot |
| ADR-79 tet10 bearing backbone outputs + `build_mesh_tet10.py` | uncommitted in `cool-haibt-781905`, 265 behind | `origin/rescue/adr79-bearing-tet10-backbone` = `530adc838` | data + one script |
| ADR-75b lane-B thread-count profile results + build logs | uncommitted in `mumps-opensees-study-f833bf`, 444 behind | `origin/rescue/adr75b-laneB-thread-results` = `d2d7d24d7` | data |

**`rescue/*` branch rule:** a rescue branch is a snapshot, never a PR head. It
therefore never auto-deletes — list them with
`git branch -r --list 'origin/rescue/*'`. The next agent on that lane
cherry-picks or rebases the snapshot onto fresh `ladruno`, then deletes the
rescue branch in the same PR.

### Dropped without rescue (verified)

| Item | Evidence |
|---|---|
| `_wt_sp_findings` staged `80_ladruno_sp_imposition_strengthening_adr.md` | **older** than the tree's copy (17+/74− against `ladruno`) |
| build logs in `compiled-version-install-db4e20`, `OpenSees-eas`, `OpenSees-cmsp4` | logs only |
| orphan folder `.claude/worktrees/adr78-p4-battery` (8,807 files, unregistered, Aug 12–14) | 4,534 byte-identical to `ladruno`; 2,740 differ by CRLF only; 1,419 untracked = a MUMPS build tree + `contact_p4/_out` logs; of the 114 real diffs, 113 versions exist in `ladruno` history and the last (`PARDISOSymLinSOE.cpp`) is the **pre-#754** `exit(-1)` version. Nothing unique. Deleted. |

### Worktrees removed (17) and local branches deleted

`.claude/worktrees/`: `compiled-version-install-db4e20`, `cool-haibt-781905`,
`ecstatic-cannon-a822ae`, `mumps-opensees-study-f833bf`,
`pardisio-profiling-0a03b1`, `sp-strengthening-9c0ddc`.
Siblings: `_wt_sp_findings`, `OpenSees-adr77review`, `OpenSees-adr78`,
`OpenSees-cmsp4`, `OpenSees-eas`, `OpenSees-hypo`, `OpenSees-integrators`,
`OpenSees-p1`, `OpenSees-rigidbody`, `OpenSees-robust`, `OpenSees-ssi-toolchain`.

Local branches deleted (all 0 ahead of `origin/ladruno`, or rescued above):
`claude/sp-prescribed-displacement-findings`, `claude/compiled-version-install-db4e20`,
`claude/latest-pr-ci-issue-0ef542`, `claude/cms-p4-verdict`, `claude/adr49a-c06-caveat`,
`claude/adr78-contact-tcl-registration`, `claude/adr78-handoff`, `claude/up-geom-comment-accuracy`,
`integ-work`, `feat/adr58-rigidbody-pytest-banner`, `guppi/robust-creduction`,
`claude/geometric-locking-coupled-816302`, `claude/adr75-p1-pardiso`,
`claude/tier-a-kstab-damage-scaling-46d9af`, `claude/sp-strengthening-9c0ddc`.

Left as found: the **main checkout** on `claude/pardiso-linux` (another agent
was drafting ADR-88 there, untracked, during the sweep — do not switch its
branch under them); `origin/claude/pardiso-linux` (its stranded commit is now
#782, so the remote branch can go once #782 lands); the empty directory
`.claude/worktrees/ecstatic-cannon-a822ae` (registration removed, directory
held open by another process — `rmdir` it when free).
