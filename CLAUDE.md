# OpenSees — Ladruno fork: working rules for agents

This is a fork of OpenSees (`github.com/nmorabowen/OpenSees`, default branch
`ladruno`). On top of upstream we add fork-only features. We move fast, so we
keep a **build-control** record of everything that diverges from upstream.

## Build-control ledgers — keep them current (REQUIRED)

Three ledgers live in `Ladruno_implementation/`. Updating them is part of the
work, not an afterthought — do it **in the same PR** as the change:

- **`LEDGER_vanilla_files.md`** — touched an *upstream* file? Add a row: file,
  why, PR. Mark the edit in source with a `// Ladruno ...` comment so the table
  is reconstructable via `grep -rn "Ladruno" SRC/`.
- **`LEDGER_implementations.md`** — added a *new* feature/file we authored? Add a
  row: feature, kind, class tag, files, status, PR. Record class tags here so we
  never collide (`SRC/classTags.h`).
- **`LEDGER_quirks.md`** — learned an OpenSees gotcha? Write it down so the next
  agent doesn't rediscover it.

## Splash-banner feature list — keep it in sync

The splash banner prints an active-feature list under the LADRUNO ASCII art.
**Do not hand-edit the C strings.** Source of truth is
`Ladruno_scripts/banner_features.txt` (one feature per line). Workflow:

1. Add/edit a line in `Ladruno_scripts/banner_features.txt`.
2. `python Ladruno_scripts/patch_banner.py` — regenerates the `FEATURES-START/END`
   blocks in `SRC/tcl/tclMain.cpp` and `SRC/interpreter/PythonModule.cpp`.
3. Rebuild: `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`.

Every `shipped` row in `LEDGER_implementations.md` should have a matching banner
line. (Banner art itself comes from `banner_ASCII.txt` → `BANNER-START/END`.)

## Doc folders

- `Ladruno_implementation/` — ledgers (above) + forward-looking feature plans.
- `Ladruno_internal/` — build history / compilation journal; deep build and
  toolchain detail. See `BUILD_GOTCHAS.md` (env/runtime workarounds: Python 3.12
  ABI, batch traps, CMake-4.3 shadow, MUMPS, test bootstrap, installer DLL-lock)
  and `WORKFLOW_GOTCHAS.md` (PR/CI/git traps on this fast-auto-merging fork:
  stranded commits, stale-PR "no checks", header stamp, vanilla footprint).
- `Ladruno_scripts/` — build, installer, banner, and test tooling.

## PRs and branches (ADR-87 D9/D10 — in force)

Base fork PRs on `ladruno` (the default branch), not `main`. See the global
`~/.claude/CLAUDE.md` for the stacked-PR `--base` pitfall.

The micro-PR/auto-merge era is over. The working shape:

- **One branch per work package**, named `wp/<n>-<slug>`, cut from fresh
  `ladruno`. One worktree per live WP.
- **Open a draft PR on day one.** Zone-A runs on every push; a draft cannot
  merge. Commit continuously to the WP branch — no intermediate PRs.
- **One deliberate merge event**: flip to ready only when the warrant package
  is complete (verification manifest + mutation gate + guide for family work,
  Zone-A green). The owner merges; agents do not.
- **Merge commits for WP PRs; squash only for small fixes.** Never
  `gh pr merge --auto`.
- Branch protection is live on `ladruno` (Zone-A + classTag/manifest checks
  required, `enforce_admins` on) and merged head branches auto-delete. Do not
  work around a red check — fix it.
- Don't accumulate long-lived personal branches: rescue or record-and-drop
  per `Ladruno_internal/branch_graveyard.md` (ADR-87 D10).
