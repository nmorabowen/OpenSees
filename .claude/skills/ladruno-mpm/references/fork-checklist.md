# Fork checklist — registrations, gates, builds, warrant

For fork WPs (232b, 232c, 232f+). The general house rules are `AGENTS.md`; the element and
material task guides (`.claude/skills/ladruno-new-element/`, `ladruno-new-material/`) also apply
when a WP adds an element (mode A's particle element) or changes a material.

## Registrations for a new fork class or command

- [ ] **classTag** in `SRC/classTags.h`, fork band ≥ 33000, in the right registry (the host is a
      PATTERN tag; `PATTERN_TAG_LadrunoPorousOverlay` = 33022). Numbers are per registry: an
      equal number in another registry is not a collision, and the comment says so.
- [ ] **Manifest row** in `Ladruno_implementation/testbed/manifest.yaml`: `feature_class:
      loadpattern` for the host (the overlay's row is the model), `element` for the mode-A
      particle; `pytest:` a real test; `tcl:` a Tcl test or `WAIVED(<reason>)`;
      `core_files:` every upstream file touched; `dispatch:` both registries.
- [ ] **Command:** one `LADRUNO_COMMAND` row in `SRC/interpreter/LadrunoCommandTable.h`; nothing
      registered by hand in `commands.cpp`, `TclWrapper.cpp` or `PythonWrapper.cpp`.
      If the host is also exposed as a `pattern` type, wire it like the overlay:
      `OpenSeesPatternCommands.cpp` (Python) and the Tcl pattern command.
- [ ] **CMake:** a `CMakeLists.txt` in the new directory and its `add_subdirectory` line, the
      way `SRC/domain/pattern/ladrunoPorousOverlay` is wired.
- [ ] **Object broker** (`FEM_ObjectBrokerAllClasses.cpp`): needed for database restart and for
      OpenSeesSP / MP. v1 refuses both, so the entry waits for phase 1's restart test; say so in
      the PR.
- [ ] **Upstream files** you edit carry a `// Ladruno …` marker at the edit, and each gets a row
      in your WP's `ledger/vanilla/` fragment.
- [ ] **Header stamps:** add the new files to `GLOBS` in `Ladruno_scripts/stamp_headers.py`, run it.
- [ ] **Banner:** a line in `Ladruno_scripts/banner_features.txt`, then
      `python Ladruno_scripts/patch_banner.py`; name it in the implementations fragment's `banner:`.

## Ledger fragments (per WP; never edit `LEDGER_*.md`)

- `Ladruno_implementation/ledger/implementations/WP-<nnn>-<slug>.md` — front matter `wp`, `title`,
  `pr`, `date`, `status`, `class_tags: ["SYMBOL=value"]`, `banner`, `section: "table"`, then the row.
- `ledger/vanilla/WP-<nnn>-<slug>.md` — one row per upstream file: file, why, PR.
- `ledger/quirks/WP-<nnn>-<slug>.md` — one `###` entry per gotcha learned. A follow-up to an
  older entry edits THAT fragment (`#### Follow-up (WP-<nnn>)`).

## Gates — run before every push (the forms Zone-A uses)

```bash
python3.12 ci/check_classtags.py
python3.12 ci/check_manifest.py
python3.12 ci/check_ladruno_commands.py
python3.12 ci/check_ledger_fragments.py
python3.12 ci/check_quirk_patterns.py
python3.12 ci/check_duplicate_defs.py
python3.12 -m ruff check --no-cache --select F811 ci tests Ladruno_scripts   # CI pins ruff 0.16.10
python3.12 Ladruno_scripts/stamp_headers.py --check
python3.12 -m pytest -q ci/test_quirk_*.py ci/test_duplicate_defs.py
```

## Build — the orchestrator launches it, in the worker's worktree

There is one build system: `Ladruno_scripts\`. Build all five targets (naming targets leaves
`opensees.pyd` stale); the tree is `build\build\Release`. Launch through WMI so the build survives
the tool's time limit and a session restart, and log to a file:

```powershell
$root = '<worker worktree>'; $log = "$root\build_wp232.log"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{
  CommandLine = "cmd.exe /c `"cd /d $root && call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat > $log 2>&1`"" }
```

Then watch the log with a bounded wait (`=== Done`, `Build of … failed`, or a stall). Check that
the log file exists before reading a build as done. Git Bash `cmd /c build.bat` is a silent no-op.
An MSVC C1001 internal compiler error is transient: rebuild. Afterwards, `ladrunoBuild` must
print the worker's head on both `dist\bin\OpenSees.exe` and `dist\bin\opensees.pyd`.

## Tests

- pytest files `tests/test_ladruno_mpm_*.py` carry `pytestmark = pytest.mark.zone_a`. Zone-A runs
  only `pytest -m zone_a`, on Ubuntu with a GCC build (`.github/workflows/ladruno.yml:287`), so an
  unmarked test never gates, and code must build with GCC as well as MSVC.
- Zone-A installs no h5py. A test that reads a `.ladruno` or particle-set file needs `h5py` added
  to that workflow's pip line in the same PR.
- Engine parity runs in the two shipped engines: classic Tcl (`OpenSees.exe`) and Python
  (`opensees.pyd`). The DL Tcl engine is registered but never instantiated.
- Oracles: `tests/*_reference.py` with `tests/*_cases.txt`; kernel header checks as
  `tests/*_kernel_check.cpp` (the pattern of `tests/hypo_kernel_check.cpp`).
- Mutation gate (ADR-87 D2), through the existing family mechanism: add `MPM` to the
  `LADRUNO_MUTATE_FAMILY` list in `CMakeLists.txt` (a vanilla-ledgered line), to the options of
  `.github/workflows/ladruno_mutation.yml` and to `Ladruno_scripts/mutation_gate.py`; hook the
  engine through `SRC/Ladruno_mutation.h`. `ZERO` sabotages the P2G internal force and `SCALE` the
  G2P position increment (ADR §8). The orchestrator builds the mutant with
  `Ladruno_scripts/mutation_build.bat`, or dispatches the CI lane; score with `mutation_gate.py`.
  Ordering hazards (G2P after `Domain::commit`; mode A's conditional `domainChange()`) are guarded
  by regression tests, not by mutants.

## Warrant package (ADR-87) — the PR flips to ready only with all of it

- [ ] oracles that are not self-referential (closed form, MMS, cross-code or two-leg), each with
      its test;
- [ ] the mutation gate: the MPM family's ZERO and SCALE mutants each turn their oracle red;
- [ ] the guide, `Ladruno_implementation/232_ladruno_mpm_guide.md`, updated for what the WP adds;
- [ ] both shipped engines exercised (classic Tcl and Python);
- [ ] for a family PR (232c and later): the out-of-family verdict and its handoff pack, committed
      per ADR-87 D4 (`reviews/handoff_<adr>.md`, `reviews/<adr>_verdict.md`);
- [ ] the manifest row, the ledger fragments and the banner line;
- [ ] Zone-A green on a dispatched run (`gh workflow run ladruno.yml --ref <branch>`; a ready-flip
      starts none).

## Release (R) — the owner's

`Ladruno_scripts/README.md`, "Release recipe": build and verify, then the maintainer tags
`ladruno-v<version>`; the release workflow builds and publishes. Agents never create the tag or
the release.
