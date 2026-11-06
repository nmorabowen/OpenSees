---
title: Ledger — vanilla files we touched
project: Ladruno
tags:
  - ledger
  - provenance
  - upstream
---

# Ledger — vanilla OpenSees files we touched

Every upstream ("vanilla") OpenSees file the Ladruno fork modifies, **why** we
touched it, and the **PR** it landed in. This is the provenance record: if we
ever rebase onto a newer upstream, this table is the list of edits to re-apply
and re-verify.

## Conventions

- **Vanilla = pre-existing upstream file.** Brand-new files we author live in
  [[LEDGER_implementations]] instead — do not list them here.
- Keep edits minimal and marked. Every Ladruno edit in a vanilla file carries a
  `// Ladruno ...` comment so `grep -rn "Ladruno" SRC/` reconstructs this table.
- One row per (file, PR). If the same file is touched by several PRs, add a row
  per PR so the history stays linear.
- PR numbers are on the fork: `github.com/nmorabowen/OpenSees` (branch `ladruno`).
- When you touch a new vanilla file, **add the row in the same PR**.

## Deleted vanilla files — how to resolve the next upstream sync

Four upstream files are **deleted**, not modified (the `DELETED` rows below).
When someone next merges `OpenSees/OpenSees` into `ladruno` and upstream has
touched one of them, git raises a **modify/delete conflict**. It is not a
warning that we did something wrong:

> **Resolution is always `git rm` — keep them deleted.** Restoring one puts
> `makeWIN.bat` back in the repo root, which is exactly the failure this fork
> spent a PR removing (agents build through it and then test a stale or
> differently-linked binary). Deliberate divergence, recorded here.

Low blast radius, measured 2026-09-08: `makeWIN.bat`, `makeMac.sh` and
`conanfile2.py` have **2 upstream commits each, ever**;
`OpenSeesAWS-Ubuntu22.04.sh` about the same. Last sync from `OpenSees/OpenSees`
was 2026-04-26. So this is one conflict, once, if ever.

**After every upstream sync, refresh the upstream manifest** (WP-122). The CI step
"header stamp covers GLOBS" treats any tracked `SRC` source missing from
`Ladruno_scripts/upstream_src_manifest.txt` as fork-authored, so upstream's new files
turn it red until you run
`git fetch upstream master && python Ladruno_scripts/stamp_headers.py --refresh-upstream-manifest`
and commit the result. Do not add upstream files to GLOBS to silence it.

This costs nothing in the other direction: the upstream PR campaign builds every
package as a *fresh branch off `jaabell/ladruño`* with files copied in
(`upstream_pr_campaign.md` — "our git history is not portable. No
cherry-picking"), so a deletion on `ladruno` can never reach a port branch.

**Deliberately NOT deleted**, though they are also unused build systems:
`Win32/` (204 upstream commits), `Win64/` (418, last touched 2026-02-19),
`MAKES/` (83), `Makefile`, `Dockerfile`, `docker/`. The first three are actively
maintained upstream, so deleting them would mean a real conflict on *every*
future sync — and unlike `makeWIN.bat` in the root, nobody reaches for a Visual
Studio solution file by accident.

## Ledger

| Vanilla file | Why touched | PR |
|---|---|---|
<!-- ledger:rows main -->

> [!note] Upstreamable bugfixes
> Some PRs fix genuine upstream bugs (not fork-only features) and are candidates
> to send back to OpenSeesFramework. Track those in the table below so we know
> what to PR upstream.

| Vanilla file | Upstreamable fix | Fork PR |
|---|---|---|
<!-- ledger:rows upstreamable -->
