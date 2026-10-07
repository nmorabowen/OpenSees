# The Windows self-hosted runner (WP-179)

Zone-A is Ubuntu without MKL. Every test with a Windows-only leg (Pardiso, FEAST,
the MSVC byte-identity baselines, the boot-`.pth` wiring) runs only on a Windows
build. Until a runner exists they have **no CI** (issue #934), and a drift can sit
unseen for weeks (WP-175).

The runner is a developer desktop, so the Windows gates are **on demand only**:
`.github/workflows/ladruno_windows.yml` has no schedule and no PR trigger, and it
is not a required check. Nothing runs unless you dispatch it with the runner up.

Registering the runner is an **owner** step. It needs a GitHub registration
token and runs a process under your account, so agents cannot do it.

## 0. Security first: this repository is PUBLIC

A pull request from a fork can edit a workflow file and point a job at
`runs-on: self-hosted`, which would run outside code on your desktop. GitHub's
own guidance is that self-hosted runners and public repositories need care.
Before registering:

1. **Settings → Actions → General → "Approval for running fork pull request
   workflows from contributors"**: choose **"Require approval for all external
   contributors"**. A fork PR's workflows then never start without your click.
2. Never add a `pull_request` / `pull_request_target` trigger to a job that runs
   on `self-hosted`. The fork's jobs only use `workflow_dispatch` and `schedule`,
   which only people with write access can cause.
3. Run the runner **interactively** (`run.cmd`) only while you want CI, not as an
   always-on service (§2). The attack surface is then only the windows you open.

## 1. Prerequisites (already true on the build desktop)

- The documented build works from a fresh `cmd.exe`:
  `call Ladruno_scripts\setup_env.bat` then `Ladruno_scripts\build.bat`
  (oneAPI, MSVC, Conan, the CMake pin, MUMPS).
- `py -3.12` works and its `site-packages` has `pytest`, `numpy`, `h5py`, `gmsh`
  ([[BUILD_GOTCHAS]] §0, §4).
- For the exact ADR-97 gate-4 leg, the host is registered in
  `Ladruno_implementation/adr97_oracle/baselines/hosts.json` (WP-177). The current
  desktop (`win32 | AMD Ryzen AI 7 PRO 350 w/ Radeon 860M`) is. Any other host
  runs the 1e-8 relative leg, and `test_gate4_this_host_has_its_own_baseline`
  SKIPs with the registration recipe.

## 2. Register the runner (once)

1. GitHub → **Settings → Actions → Runners → New self-hosted runner →
   Windows, x64**. Keep that page open; it shows a short-lived token.
2. In PowerShell, in a folder **outside** any repo, OneDrive or SeaDrive (for
   example `C:\actions-runner`), run the page's *Download* commands.
3. Configure with ONLY the `ladruno-win` label:

   ```powershell
   .\config.cmd --url https://github.com/nmorabowen/OpenSees --token <TOKEN-FROM-THE-PAGE> --name $env:COMPUTERNAME --labels ladruno-win --work _work
   ```

   Answer **N** to "run as service". A service runs as `NETWORK SERVICE`,
   which has no access to your oneAPI activation, your Conan cache
   (`%USERPROFILE%\.conan2`) or your `py -3.12` install, so every build would fail.

   **Why `ladruno-win` and not `ladruno-perf`:** the nightly jobs in
   `ladruno.yml` (`zone-b-nightly`, `cross-tier-nightly`) target `ladruno-perf`
   and fire at 06:00 UTC. Without that label they stay dormant. Add it later
   (`.\config.cmd remove`, then reconfigure with `--labels ladruno-win,ladruno-perf`)
   only if you want them on this desktop.

## 3. Use it

```powershell
cd C:\actions-runner; .\run.cmd          # leave this window open while CI runs
```

From any shell:

```bash
gh workflow run ladruno_windows.yml --ref <branch>                       # platform tests (default)
gh workflow run ladruno_windows.yml --ref <branch> -f scope=zone_a -f slow=true
gh run watch
```

- `scope=platform` runs every `test_*.py` with a non-portable platform branch,
  listed by `python ci/check_quirk_patterns.py --list-platform-tests`. Nobody
  maintains that list: a new win32-only test joins it automatically.
- The run's check appears on that commit, and therefore on its PR.
- Close `run.cmd` (Ctrl+C) when done. Queued dispatches wait for the runner.

The first run builds from scratch (about 40 min with MUMPS). Later runs reuse
`_work\OpenSees\OpenSees\build` (`clean: false`) and take minutes, unless the
C++ changed.

## 4. What a green run proves, and the exit codes

Three checks, each closing a failure met in WP-175..178:

1. `CMakeLists.txt` is touched, so CMake reconfigures and re-stamps
   `ladrunoBuild()` with this commit ([[BUILD_GOTCHAS]] "ladrunoBuild lags").
2. `build.bat` byte-verifies `dist\` against the build and fails on a locked or
   stale copy (WP-178, [[BUILD_GOTCHAS]] §5c).
3. `Ladruno_scripts/ci_run_pytest.py` runs under `python -S` (no boot `.pth`),
   pins `dist\bin`, and requires `ladrunoBuild() == github.sha`.

| Failure | Meaning |
|---|---|
| job stays *Queued* | the runner is not running (`run.cmd`) or lacks the `ladruno-win` label |
| build step: `could not copy … open in another process` | something on the desktop has `_work\…\dist\bin\opensees.pyd` loaded; close it and re-dispatch |
| exit **90** | the engine is the wrong file, or `ladrunoBuild()` ≠ the commit (stale `dist\` or a stamp that did not refresh) |
| exit **91** | the launcher was started without `-S` |
| `--list-platform-tests returned nothing` | the selector broke; the job refuses to report a green run of zero tests |

## 5. Open items

- Flip the `# ci-coverage: local-only` notes on the platform tests to a CI kind
  (and add an `on-demand-windows` kind to lint L8) only **after** the first green
  dispatch on this runner. Claiming coverage before that would be false.
- `zone-b-nightly`'s perf step still runs a bare `python` (not pinned). It only
  matters if you opt into `ladruno-perf`.
