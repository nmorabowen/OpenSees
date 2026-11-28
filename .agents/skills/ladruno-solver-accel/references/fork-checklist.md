# WP-260 Fork Checklist — Registrations, Gates, Builds & Review

## 1. Static CI Gates (Run Before Every Push)

The same commands `.github/workflows/ladruno.yml` runs on every PR, drafts included. From the worktree root:

```bash
pip install ruff==0.16.10 pytest h5py numpy        # versions CI pins
python ci/check_duplicate_defs.py
ruff check --no-cache --select F811 ci tests Ladruno_scripts
python ci/check_classtags.py
python ci/check_ledger_fragments.py
pytest -q ci/test_ledger.py
python ci/ledger.py build
python ci/check_manifest.py
(cd Ladruno_scripts/ladruno_recorder_tests && pytest -q test_ladruno_harness.py)
pytest -q -p no:cacheprovider tests/test_mpm_contracts.py
python Ladruno_scripts/stamp_headers.py --check
pytest -q ci/test_quirk_*.py ci/test_duplicate_defs.py
python ci/check_quirk_patterns.py
pytest -q ci/test_check_viewer_ledger.py
python ci/check_viewer_ledger.py --base origin/ladruno --head HEAD
```

If `ladruno.yml` changes, this list follows it, not the other way round.

### 1b. MKL edits — oneMKL 2024.2.2 floor

Required for any edit under `SRC/system_of_eqn/linearSOE/pardiso/`, a `_PARDISO` or `_LADRUNO_MKL_FEAST` block, or a new MKL routine/macro/`iparm` value (`pardiso_64`, `mkl_rci.h` `dfgmres`, `iparm[1|10|12|23|24|27]`):

- Every symbol must exist in **oneMKL 2024.2.2** (esmeralda, the oldest fork MKL). Zone-A configures without MKL, so only `ci/check_mkl_compat.py` (WP-148) compiles these sources against it. It runs as its own job in `ladruno.yml` on every push, drafts included, on Linux against the 2024.2.2 headers and libraries (`--build-dir build/mkl --self-test`); it cannot run on the Windows box. A green local Windows build proves nothing here: confirm that job is green on the PR before claiming the MKL edit works.
- Changes inside the WP-132 `-deterministic` (CNR) block of `PARDISOGenLinSolver.cpp` re-prove CNR bit identity; new parallel `iparm` options are refused under `-deterministic` unless proved.

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

A follow-up to a WP-260 quirk edits that quirk's fragment (`#### Follow-up (WP-260<x>)`); if a quirk heading changes, update the matching `Quirks:` pointer in `SKILL.md` (lint rule `pointers` checks it).

## 4. Registrations (When Adding New Classes or Commands)

- `SRC/classTags.h` (fork band `>= 33000`, right registry)
- `Ladruno_implementation/testbed/manifest.yaml`
- `SRC/interpreter/LadrunoCommandTable.h` (if adding Tcl/Python commands)
- `CMakeLists.txt` + `FEM_ObjectBrokerAllClasses.cpp` (where required for SP/MP/restart)
- `Ladruno_scripts/stamp_headers.py` (`GLOBS` + `--stamp`)
- `Ladruno_scripts/banner_features.txt` + `python Ladruno_scripts/patch_banner.py` (for shipped features)

## 5. Warrant Package Before Ready-Flip

1. All static CI gates of §1 green (and §1b for MKL edits).
2. Local 5-target build green; `ladrunoBuild` SHA matches PR `HEAD`.
3. `pytest -m zone_a` green (both `OpenSees.exe` and `opensees.pyd` exercised where applicable).
4. Explicit ADR-68 gate proof (`G-BYTE` MD5/signbit identity, `G-EQUIV` oracle, or `G-OPT` + fail-closed checks; TSan for threaded changes).
5. Before/after profile table pasted in the PR description, measured on a quiet machine with a load probe per run (`review-protocol.md` §0 rule 5): the D1 decks the WP was picked on, method of the `260-w0` brief §2 (`1t_deep` for phase shares, `1t_coarse` and `nt_coarse` for wall time and scaling, default thread environment, step 1 dropped, interleaved reps, median and spread), phase speedup `s` and whole-run gain.
6. Zone-A Ubuntu CI dispatched (`gh workflow run ladruno.yml --ref <branch>`) and green (`gh pr checks <pr>`).
7. `### Review request` posted with its evidence table (`references/review-protocol.md` §0–§1: every claim has an artifact, every oracle names what it compares against, CI status copied from `gh pr checks`, approvals linked); Claude's verdict is `ready`, every finding answered by number.
8. The PR states the plan §0.5 numbers it was picked on (p, s, projected G) next to the measured whole-run gain.
9. Family-level PR (per D7): `Ladruno_implementation/reviews/handoff_260<x>.md` written for the owner's ADR-87 D4 run; `reviews/adr260<x>_verdict.md` committed.
10. `Ladruno_internal/ISSUES.md` §3 `### Handoff` comment posted on the PR and board `#937` refreshed.
