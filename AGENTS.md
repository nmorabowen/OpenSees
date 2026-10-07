# OpenSees — Ladruno fork: working rules for agents

This is a fork of OpenSees (`github.com/nmorabowen/OpenSees`, default branch
`ladruno`). On top of upstream we add fork-only features. We move fast, so we
keep a **build-control** record of everything that diverges from upstream.

## Task guides — read the matching one before starting

Short checklists for specific kinds of work. Each item links to the ledger entry
that explains it; the guides point into the ledgers, they do not replace them.

| Doing this | Read first |
|---|---|
| Adding or changing an element (`SRC/element/`) | [`.claude/skills/ladruno-new-element/SKILL.md`](.claude/skills/ladruno-new-element/SKILL.md) |
| Adding or changing a material (`SRC/material/`) | [`.claude/skills/ladruno-new-material/SKILL.md`](.claude/skills/ladruno-new-material/SKILL.md) |
| Changing the profiler or monitor viewers (`Ladruno_tools/`) | [`.claude/skills/ladruno-viewer-tools/SKILL.md`](.claude/skills/ladruno-viewer-tools/SKILL.md) |

The mechanical items are enforced by `python ci/check_quirk_patterns.py` (a
static gate in `ladruno.yml`). When a quirk names a greppable pattern, add it
to that script instead of adding another paragraph to a guide.

Lint rules, tests and registrations are named by slug or WP, never by the next
free number (a rule is a `RULES` slug with its own `ci/test_quirk_<slug>.py`);
ruff F811 and `ci/check_duplicate_defs.py` fail a duplicate definition. Never add
`merge=union` to code or tests.

## Building — there is exactly ONE way (REQUIRED)

From a fresh `cmd.exe` at the fork root:

```cmd
call Ladruno_scripts\setup_env.bat
Ladruno_scripts\build.bat               :: all 5 targets -> dist\
Ladruno_scripts\build.bat installer     :: ...and wrap dist\ into the setup.exe
Ladruno_scripts\build.bat clean installer   :: the full release recipe
```

`Ladruno_scripts/` is the **only** build system in this fork. The upstream
recipes (`makeWIN.bat`, `makeMac.sh`, `OpenSeesAWS-Ubuntu22.04.sh`,
`conanfile2.py`) and the legacy zip packager (`make_installer.ps1`) were
**deleted** — they had drifted (wrong MUMPS path, lp64-vs-ilp64 MKL, no
`dist/`, no installer) and agents kept reaching for them because they sat in
the repo root. If a doc still cites one, that doc is stale: fix it.

Four rules that cost real time when broken:

1. **The build tree is `build\build\Release`**, not `build\Release` — Conan's
   `cmake_layout` nests it. A hand-rolled `cmake --build build/Release …`
   configures a *second, different* cache; you then "build" successfully and
   test a stale binary. Always go through `build.bat`.
   (`Ladruno_internal/02_esmeralda_linux_build_guide.md` legitimately uses
   `build/Release` — different machine, different layout. Not a contradiction.)
2. **Naming targets skips the rest.** `build.bat OpenSees` leaves
   `dist\bin\opensees.pyd` stale, so pytest tests an old binary while Tcl tests
   the new one. Check mtimes, or just build everything.
3. **Packaging requires a full build.** A 4-of-5 build drops `dist\openseesmp`
   entirely and the installer then ships a MIXED install, silently.
   `build.bat installer` refuses explicit targets for this reason.
4. **Build in a worktree, never the shared main checkout** — its branch moves
   under you, and commits/builds leak onto another agent's branch.

Artifacts: `dist\bin\` (`OpenSees.exe`, `OpenSeesSP.exe`, `OpenSeesMP.exe`,
`opensees.pyd`) + `dist\openseesmp\` (`openseesmp.pyd` with its own Intel MPI
runtime). The installer lands in
`Ladruno_files\Ladruno_OpenSees_<version>_setup.exe`.

Details: [BUILDING.md](BUILDING.md) (what/how),
[Ladruno_scripts/README.md](Ladruno_scripts/README.md) (which script does
what), [Ladruno_internal/BUILD_GOTCHAS.md](Ladruno_internal/BUILD_GOTCHAS.md)
(why it broke last time).

## Build-control ledgers — keep them current (REQUIRED)

Three ledgers, kept as **per-WP fragments** (WP-161) so no two PRs write the
same file. Updating them is part of the work, not an afterthought — do it **in
the same PR** as the change:

- **Ledgers are per-WP fragments:** write
  `Ladruno_implementation/ledger/<kind>/WP-<nnn>-<slug>.md`, and never edit
  `LEDGER_*.md`, which is generated (the committed files are stubs; the gate
  fails on an edited one). Format: `Ladruno_implementation/ledger/README.md`.
- **Read** with `rg <pattern> Ladruno_implementation/ledger/` or
  `python ci/ledger.py build` (writes `Ladruno_implementation/ledger/_build/`).
- **A follow-up to an older entry edits that WP's fragment**; never duplicate it.
- **`vanilla/`** — touched an *upstream* file? One fragment per PR, one row per
  file: file, why, PR. Mark the edit in source with a `// Ladruno ...` comment
  so the rows are reconstructable via `grep -rn "Ladruno" SRC/`.
- **`implementations/`** — added a *new* feature/file we authored? A row:
  feature, kind, class tag, files, status, PR; put the tags in `class_tags` so
  the gate checks them against `SRC/classTags.h`.
- **`quirks/`** — learned an OpenSees gotcha? One `### ` entry, so the next
  agent doesn't rediscover it.

## GitHub issues are the shared memory (REQUIRED)

Live state lives on GitHub, never in a session's private memory: a WP's state
and handoffs on its **draft PR**, a multi-WP program on its **`program`
issue**, everything open on the pinned **board** issue (#937), each out-of-scope bug
or gap in **its own issue**. Boot from the board and the last handoff comment;
before the session ends, post a handoff and refresh the board row. Private
memory keeps preferences and pointers ("program X is #n"), not open-PR lists,
job ids or "resume here" notes. Something learned that must not be
rediscovered graduates from an issue into a ledger quirk or a lint rule.
Protocol, labels and templates: [Ladruno_internal/ISSUES.md](Ladruno_internal/ISSUES.md).

## Splash-banner feature list — keep it in sync

The splash banner prints an active-feature list under the LADRUNO ASCII art.
**Do not hand-edit the C strings.** Source of truth is
`Ladruno_scripts/banner_features.txt` (one feature per line). Workflow:

1. Add/edit a line in `Ladruno_scripts/banner_features.txt`.
2. `python Ladruno_scripts/patch_banner.py` — regenerates the `FEATURES-START/END`
   blocks in `SRC/tcl/tclMain.cpp` and `SRC/interpreter/PythonModule.cpp`.
3. Rebuild: `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`.

Every `shipped` implementations fragment should have a matching banner line
(name it in the fragment's `banner:` and the gate checks it exists). (Banner art itself comes from `banner_ASCII.txt` → `BANNER-START/END`.)

## Doc folders

- `Ladruno_implementation/` — ledger fragments in `ledger/` (above) + forward-looking feature plans.
- `Ladruno_internal/` — build history / compilation journal; deep build and
  toolchain detail. See `BUILD_GOTCHAS.md` (env/runtime workarounds: Python 3.12
  ABI, batch traps, CMake-4.3 shadow, MUMPS, test bootstrap, installer DLL-lock)
  and `WORKFLOW_GOTCHAS.md` (PR/CI/git traps on this fast-auto-merging fork:
  stranded commits, stale-PR "no checks", header stamp, vanilla footprint)
  and `ISSUES.md` (issues as shared memory: board, handoffs, labels).
- `Ladruno_scripts/` — build, installer, banner, and test tooling.
- `.claude/skills/` — the task guides above (tracked; the rest of `.claude/` is
  session scratch and stays ignored).

## PRs and branches (ADR-87 D9/D10 — in force)

Base fork PRs on `ladruno` (the default branch), not `main`. See the global
`~/.claude/CLAUDE.md` for the stacked-PR `--base` pitfall.

The micro-PR/auto-merge era is over. The working shape:

- **One branch per work package**, named `wp/<n>-<slug>`, cut from fresh
  `ladruno`. One worktree per live WP.
- **Open a draft PR on day one.** Drafts run the fast static gates on every
  push (the full Zone-A build intentionally skips drafts; trigger it on
  demand with `workflow_dispatch`); a draft cannot merge. Commit
  continuously to the WP branch — no intermediate PRs.
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
