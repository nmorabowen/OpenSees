# `Ladruno_scripts/` — the fork's only build system

There used to be two ways to compile this fork. There is now one. The upstream
recipes (`makeWIN.bat`, `makeMac.sh`, `OpenSeesAWS-Ubuntu22.04.sh`,
`conanfile2.py`) and the legacy zip packager (`make_installer.ps1`) were deleted
— they had drifted out of date and, because `makeWIN.bat` sat in the repo root,
they were what agents found first.

## The build chain

From a fresh `cmd.exe` at the fork root:

```cmd
call Ladruno_scripts\setup_env.bat        :: vcvars64 + 3x oneAPI vars.bat
Ladruno_scripts\build.bat                 :: Conan -> CMake -> Ninja -> dist\
Ladruno_scripts\build.bat installer       :: ...then wrap dist\ into a setup.exe
```

| Invocation | Effect |
|---|---|
| `build.bat` | Build all 5 targets, refresh `dist\` |
| `build.bat clean` | Wipe `build/`, `install/`, `dist/` first |
| `build.bat rebuild` | Wipe `build/` only |
| `build.bat OpenSees OpenSeesPy` | Build just those targets *(leaves the rest stale)* |
| `build.bat installer` | Full build, then `build_inno_installer.ps1` |
| `build.bat clean installer` | The full release recipe |

`installer` is **refused** alongside explicit targets: `build.bat` refreshes
`dist\` only for the targets it built, so packaging a partial build produces a
plausible-looking but mixed `setup.exe`.

### Outputs

| Path | Contents |
|---|---|
| `dist\bin\` | `OpenSees.exe`, `OpenSeesSP.exe`, `OpenSeesMP.exe`, `opensees.pyd`, Tcl + MKL runtime |
| `dist\openseesmp\` | `openseesmp.pyd` + its own Intel MPI runtime (deliberately separate from `bin\`) |
| `Ladruno_files\Ladruno_OpenSees_<version>_setup.exe` | The wizard installer |

### The build tree is `build\build\Release`

Conan's `cmake_layout` nests it. `build\Release` is **not** the build tree on
Windows. A hand-rolled `cmake --build build/Release …` configures a second,
independent cache — it will appear to succeed while you go on to test a stale
binary. Go through `build.bat`.

The one legitimate exception is Linux: `Ladruno_internal/02_esmeralda_linux_build_guide.md`
configures `-B build/Release` directly on the Esmeralda cluster. Different
machine, different layout — not a contradiction, and not a Windows recipe.

## Script inventory

### Build & package
| Script | Purpose |
|---|---|
| `install_prereqs.ps1` | One-time per machine: VS 2022, oneAPI, CMake, Ninja, Python 3.12, Conan, Inno Setup |
| `setup_env.bat` | Loads the toolchain into the current shell; refuses with the exact missing tool |
| `build.bat` | The build driver (Conan → CMake → Ninja → `dist\`) |
| `opensees-msvc-static.profile` | Conan profile used by `build.bat` |
| `build_inno_installer.ps1` | Wraps `dist\` into `Ladruno_files\*_setup.exe` via `iscc.exe` |
| `installer.iss` | Inno Setup script (venv-picker wizard) driving the above |

### Banner
| Script | Purpose |
|---|---|
| `banner_features.txt` | **Source of truth** for the splash feature list — edit this, not the C strings |
| `banner_ASCII.txt` | The LADRUNO ASCII art |
| `patch_banner.py` | Regenerates the `FEATURES-START/END` + `BANNER-START/END` blocks in `tclMain.cpp` and `PythonModule.cpp` |

### Gates & tests
| Script | Purpose |
|---|---|
| `run_zone_a.ps1` | The Zone-A no-regression suite (the CI gate) |
| `mutation_build.bat` / `mutation_gate.py` | Mutation-testing gate for the warrant package (ADR-87) |
| `stamp_headers.py` | Header provenance stamps |
| `bezier_tests/`, `rigidbody_tests/`, `ladruno_recorder_tests/`, `robust_solve_tests/`, `zfp_benchmark/` | Per-feature test batteries |
| `verify_*.tcl` | Tcl parity checks (classic / modal / FEAST / ARPACK-MP-MUMPS) |

### Environment wiring
| Script | Purpose |
|---|---|
| `wire_pyenv.ps1`, `wire_venv_pth.py` | Point a venv at a built/installed `opensees.pyd` |

Everything else (`ladruno_solve.py`, `robust_drive.py`, `analyze_augmented.py`,
`*_study.py`, `*_figures.py`, `make_profiler_dump.py`, …) is analysis tooling,
not part of the build.

## Editing the `.bat` files

They are **CRLF**. `cmd.exe` parses batch files by byte offset, so a file
silently converted to LF fails with garbage like `'M' is not recognized as an
internal or external command` (that is `REM` losing its first two bytes). Git
Bash's `cat -A` and `sed` hide the CRs, so they will not warn you. After editing
a `.bat` with any Unix-side tool, check:

```bash
file Ladruno_scripts/build.bat   # must say "with CRLF line terminators"
```

## See also

- [`BUILDING.md`](../BUILDING.md) — the end-user procedure (what and how)
- [`Ladruno_internal/BUILD_GOTCHAS.md`](../Ladruno_internal/BUILD_GOTCHAS.md) — env/runtime workarounds
- [`Ladruno_internal/01_compilation_journal.md`](../Ladruno_internal/01_compilation_journal.md) — why the toolchain looks like this
