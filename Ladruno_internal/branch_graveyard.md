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
