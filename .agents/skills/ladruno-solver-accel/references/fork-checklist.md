# WP-260 Fork Checklist — Registrations, Gates, Builds & Claude Review

## 1. Static CI Gates (Run Before Every Push)

From the worktree root:

```bash
python ci/check_classtags.py
python ci/check_manifest.py
python ci/check_ladruno_commands.py
python ci/check_ledger_fragments.py
python ci/check_quirk_patterns.py
python ci/check_duplicate_defs.py
python -m ruff check --no-cache --select F811 ci tests Ladruno_scripts
python Ladruno_scripts/stamp_headers.py --check
python -m pytest -q ci/test_quirk_*.py ci/test_duplicate_defs.py
```

## 2. Building on Windows (Exactly One Way)

Run in the worktree root (never simultaneously across two worktrees):

```cmd
cmd.exe /c "cd /d <worktree> && call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat"
```

- Always build all 5 targets (never partial targets when testing Python + Tcl).
- Build tree is `build\build\Release`; artifacts land in `dist\bin\` and `dist\openseesmp\`.
- Verify `ops.ladrunoBuild()` (Python) and `ladrunoBuild` (Tcl) match your worktree `HEAD` SHA before trusting any benchmark or test result.

## 3. Ledger Fragments (WP-161)

Never edit `LEDGER_implementations.md`, `LEDGER_quirks.md`, or `LEDGER_vanilla_files.md`. Add per-WP fragments:
- `Ladruno_implementation/ledger/implementations/WP-260<x>-<slug>.md`
- `Ladruno_implementation/ledger/vanilla/WP-260<x>-<slug>.md` (mark every upstream edit with `// Ladruno ...`)
- `Ladruno_implementation/ledger/quirks/WP-260<x>-<slug>.md`

## 4. Registrations (When Adding New Classes or Commands)

- `SRC/classTags.h` (fork band `>= 33000`, right registry)
- `Ladruno_implementation/testbed/manifest.yaml`
- `SRC/interpreter/LadrunoCommandTable.h` (if adding Tcl/Python commands)
- `CMakeLists.txt` + `FEM_ObjectBrokerAllClasses.cpp` (where required for SP/MP/restart)
- `Ladruno_scripts/stamp_headers.py` (`GLOBS` + `--stamp`)
- `Ladruno_scripts/banner_features.txt` + `python Ladruno_scripts/patch_banner.py` (for shipped features)

## 5. Warrant Package for Intensive Claude Review

Before requesting final review / flipping Draft PR to Ready:
1. All static CI gates green.
2. Local 5-target build green; `ladrunoBuild` SHA matches PR `HEAD`.
3. `pytest -m zone_a` green (both `OpenSees.exe` and `opensees.pyd` exercised where applicable).
4. Explicit ADR-68 gate proof (`G-BYTE` MD5/signbit identity, `G-EQUIV` oracle, or `G-OPT` + fail-closed checks; TSan for threaded changes).
5. Before/after `-perStep` + `-deep` benchmark table pasted in the PR description.
6. Zone-A Ubuntu CI dispatched (`gh workflow run ladruno.yml --ref <branch>`) and green (`gh pr checks <pr>`).
7. `Ladruno_internal/ISSUES.md` §3 `### Handoff` comment posted on the PR and board `#937` refreshed.
