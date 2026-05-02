# Building Ladruno OpenSees

End-user procedure for compiling this fork from source on Windows. All commands are run from the fork's root directory (the directory containing this `BUILDING.md`).

For the *why* behind toolchain choices, patches, and gotchas, see [Ladruno_internal/01_compilation_journal.md](Ladruno_internal/01_compilation_journal.md). This file only documents the *what* and *how*.

## TL;DR

```cmd
:: from a fresh cmd.exe in the fork root
call Ladruno_scripts\setup_env.bat
Ladruno_scripts\build.bat clean
```

About 10–15 min later, four binaries land in `dist\bin\`:

| Binary | What it is |
|---|---|
| `OpenSees.exe` | Sequential Tcl interpreter |
| `OpenSeesSP.exe` | Partitioned-domain MPI Tcl interpreter |
| `OpenSeesMP.exe` | Multi-interpreter MPI |
| `opensees.pyd` | Python module (OpenSeesPy, Python 3.12 ABI) |

## Prerequisites (one-time per machine)

| Tool | How to install | Used for |
|---|---|---|
| Visual Studio 2022 Community | [Microsoft installer](https://visualstudio.microsoft.com/), select the C++ workload + Windows SDK | `cl.exe`, `link.exe` |
| Intel oneAPI Base + HPC | [Intel installer](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html) | `ifx.exe` (Fortran), MKL, Intel MPI |
| CMake 3.23+ | `winget install Kitware.CMake` | Build generator |
| Ninja | `winget install Ninja-build.Ninja` | Fast build executor |
| Python 3.12 | `winget install Python.Python.3.12` (or via the py launcher) | OpenSeesPy ABI target |
| Conan 2.x | `py -3.12 -m pip install --user conan` | Fetches HDF5, Tcl, Zlib, Eigen3 |
| Inno Setup 6 (optional) | `winget install JRSoftware.InnoSetup` | Building the wizard `setup.exe` |
| MUMPS 5.5.1 source | Download `MUMPS_5.5.1.tar.gz` from [mumps-solver.org](https://mumps-solver.org/) and place at `mumps-archive/mumps_src.tar.gz` | Sparse direct solver for SP/MP |

If any tool is missing, `setup_env.bat` prints exactly which one and refuses to continue.

## The build chain

```
setup_env.bat   →   build.bat [mode]
   │                    │
   │                    ├── builds MUMPS 5.5.1 (first run only)
   │                    ├── conan install (fetches HDF5/Tcl/Zlib/Eigen)
   │                    ├── cmake configure
   │                    ├── ninja build (≈1860 compile units)
   │                    ├── cmake --install
   │                    └── copy curated binaries to dist/bin
   │
   └── activates: VS 2022 Developer env, Intel oneAPI (compiler + MKL + MPI),
       Conan on PATH, NoDefaultCurrentDirectoryInExePath cleared
```

`setup_env.bat` is idempotent — re-running it in the same shell is safe but unnecessary. `build.bat` *requires* it to have run first; if the toolchain isn't activated, the build dies immediately with a clear error.

## Build modes

`Ladruno_scripts\build.bat <mode>` accepts:

| Mode | What it does | When to use |
|---|---|---|
| *(no arg)* | Incremental — Ninja rebuilds only changed targets | Default for daily edits |
| `clean` | Wipes `build/`, `install/`, `dist/`, then full rebuild | After toolchain change, Python ABI bump, or when something inexplicably breaks |
| `rebuild` | Wipes `build/` only (keeps Conan/MUMPS caches), full rebuild | Faster reset than `clean`; reuses MUMPS and Conan downloads |
| `<target>` | Builds just one target. Valid: `OpenSees`, `OpenSeesSP`, `OpenSeesMP`, `OpenSeesPy` | Quick iteration on a single front end |

## Output

```
dist/                         ← gitignored; regenerable
└── bin/
    ├── OpenSees.exe
    ├── OpenSeesSP.exe
    ├── OpenSeesMP.exe
    ├── opensees.pyd
    └── *.dll                 ← Intel MKL/iomp runtimes (co-located so opensees.pyd resolves them)
```

`dist/bin/` is the canonical run location. Binaries elsewhere (`build/`, `install/`) are intermediate.

## Wiring OpenSeesPy into a Python venv

For `import opensees` to work, the venv needs a `.pth` file pointing at `dist/bin/`. Two ways:

### Quick (developer use)

```powershell
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\wire_pyenv.ps1
```

Creates a venv at `<fork-root>\opensees_venv` (Python 3.12), drops the `.pth`, runs a smoke test. Override target:

```powershell
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\wire_pyenv.ps1 -VenvPath D:\envs\my-venv -Force
```

### Distributable (end-user use)

Build a packaged installer that any Windows user can run — see *Packaging for distribution* below.

## Packaging for distribution

Two parallel installer formats, both drawing from the same `dist/` and writing to `Ladruno_files/`.

### Wizard installer (`setup.exe`) — recommended for non-developers

```powershell
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\build_inno_installer.ps1
```

Produces `Ladruno_files\Ladruno_OpenSees_<version>_setup.exe`: a single self-contained Inno Setup wizard with a venv-picker page, optional PATH entry, and a `.pth`-writing helper. Default install location is `Program Files\Ladruno\OpenSees` (UAC-elevated). Requires Inno Setup 6 installed locally.

### PowerShell installer (`install.ps1` + `.zip`) — fallback for CI / power users

```powershell
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\make_installer.ps1
```

Produces `Ladruno_files\install.ps1` and `Ladruno_OpenSees_<version>.zip`. Distribute both files together. Receivers run:

```powershell
powershell -ExecutionPolicy Bypass -File install.ps1 `
    -InstallPath C:\Tools\Ladruno -VenvPath C:\Tools\Ladruno\venv -AddToPath -Yes
```

The `-Yes` flag suppresses prompts — useful for unattended/silent installs.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `setup_env.bat` says `[MISSING] cl.exe` | VS 2022 not installed or not on the standard path | Reinstall VS, or pass `-VsPath <path>` if your install is non-standard |
| `setup_env.bat` says `[MISSING] ifx.exe` | Intel oneAPI HPC kit not installed | Install Intel oneAPI HPC toolkit (the Base toolkit alone is not enough — Fortran lives in HPC) |
| `setup_env.bat` says `[MISSING] conan.exe` | Conan not installed under Python 3.12 | `py -3.12 -m pip install --user conan` |
| Build fails with `MUMPS_5.5.1.tar.gz: No such file` | Source archive not placed at expected path | Download from [mumps-solver.org](https://mumps-solver.org/) and put at `mumps-archive\mumps_src.tar.gz` |
| `import opensees` → `ImportError: DLL load failed` | Venv Python ABI doesn't match the `.pyd` | Use Python 3.12 for any venv that wires this build. Check with `python -c "import sys; print(sys.version)"` |
| `OpenSees` not found after wizard install | Existing shell has stale `PATH` from before the install | Open a fresh `cmd.exe` (or `$env:Path = [Environment]::GetEnvironmentVariable('Path','User')` in the current PowerShell) |
| `OpenSeesSP.exe` / `OpenSeesMP.exe` says `mpiexec not found` | Intel MPI not on PATH at runtime | `setup_env.bat` adds it; or run from any shell where `setvars.bat` has been called |

For deeper diagnostic context (the eight CMakeLists patches, MUMPS LP64-vs-ILP64 decision, oneAPI activation pitfalls), see the [compilation journal](Ladruno_internal/01_compilation_journal.md).
