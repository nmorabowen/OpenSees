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
  toolchain detail.
- `Ladruno_scripts/` — build, installer, banner, and test tooling.

## PRs

Base fork PRs on `ladruno` (the default branch), not `main`. See the global
`~/.claude/CLAUDE.md` for the stacked-PR `--base` pitfall.
