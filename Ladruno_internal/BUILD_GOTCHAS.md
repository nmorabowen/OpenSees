---
title: Ladruno build / environment gotchas
project: Ladruno
tags:
  - internal
  - build
  - environment
---

# Build & environment gotchas

The "why" behind the workarounds baked into `setup_env.bat` / `build.bat`, plus
the runtime/test-env facts you need *after* a successful build. The fixes are in
the scripts; this file explains them so a fresh clone (or an upstream re-clone)
isn't blind when something breaks.

> Ported from machine-local agent memory so it survives a clone on another
> machine. The narrative "why a CMake patch exists" lives in
> [[01_compilation_journal]]; this file is the env/runtime layer around it.

---

## 0. Python version — the build targets CPython 3.12 (NOT 3.11)

**Authoritative current state:** the built `dist/bin/opensees.pyd` is linked
against **`python312.dll`** (verify with `dumpbin /dependents`). It loads ONLY
under CPython 3.12.

- Build Python: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`
  (py-launcher / Microsoft Store install; has numpy + pytest). Must be passed
  explicitly as `-DPython_EXECUTABLE=` so CMake doesn't pick up some other Python
  on PATH.
- **Conan** runs under 3.11 (`...\Python311\Scripts\conan.exe`) and still works
  fine as a *build-tool driver* even though the engine targets 3.12. If 3.11 is
  removed, reinstall via `py -3.12 -m pip install --user conan` and update
  `setup_env.bat`.

> NOTE: parts of [[01_compilation_journal]] still say "Python 3.11" — those lines
> describe an earlier era. The pyd ABI target is **3.12** as of the current build.
> See §4 for the test bootstrap.

---

## 1. Four Windows batch traps worked around in `setup_env.bat`

Environment-specific quirks with **misleading symptoms** (they won't appear on a
clean Windows install or on Mac/Linux). All four are already mitigated in
`setup_env.bat`; this is why those mitigations exist:

1. **No `if (...)` blocks when paths contain `(x86)`.** The oneAPI path
   `C:\Program Files (x86)\Intel\oneAPI\` has parens; inside a parenthesized `if`
   body in cmd.exe those parens prematurely close the body. Use goto-style error
   labels instead of `if not exist X (echo err)`.
2. **Intel `setvars.bat` iteration is broken on this machine.**
   `call "...\oneAPI\setvars.bat" intel64 vs2022` emits
   `'vars.bat' is not recognized` for every component. Sidestep by calling each
   component's `env\vars.bat` directly:
   ```bat
   call "%ONEAPI_ROOT%\compiler\latest\env\vars.bat" intel64
   call "%ONEAPI_ROOT%\mkl\latest\env\vars.bat" intel64
   call "%ONEAPI_ROOT%\mpi\latest\env\vars.bat"
   ```
   (The `vars.bat is not recognized` symptom *looks* like a broken Intel install
   but is unrelated.)
3. **`vswhere.exe` must be on PATH for Conan/TCL nmake builds.** It lives at
   `C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe`, not on
   PATH by default. Without it, Conan's `conanvcvars.bat` and TCL's `makefile.vc`
   can't locate VS → TCL fails with `nmakehlp not recognized` then
   `rules.vc(1070): syntax error`. `setup_env.bat` prepends the Installer dir.
4. **`NoDefaultCurrentDirectoryInExePath=1` is set in the user environment.**
   This blocks cmd from auto-searching `.` for executables. TCL 8.6.11's
   `makefile.vc` builds `nmakehlp.exe` in cwd then invokes it as bare `nmakehlp`
   (no `.\`), which fails. `setup_env.bat` clears it **for the build session only**
   (`set "NoDefaultCurrentDirectoryInExePath="`) — the user/machine setting is
   untouched. Don't change it globally.

---

## 2. CMake 4.3 (system) shadows the conan-pinned 3.31.12 → ZLIB-not-found

**Symptom:** `cmake` configure dies at `find_package(ZLIB REQUIRED)` → "Could NOT
find ZLIB", preceded by `-- NOT USING CONAN (System Search)`. A failed configure
also deletes `build/build/Release/CMakeFiles/rules.ninja`, so `ninja` then errors
"loading 'CMakeFiles\rules.ninja'".

**Root cause:** the project pins **CMake 3.31.12** via Conan, but historically
`setup_env.bat` took whatever `cmake.exe` was first on PATH. A system CMake
upgrade to **4.3.0** (`C:\Program Files\CMake\bin\`) shadows it, and CMake 4.3
resolves builtin module-mode `FindZLIB` instead of Conan 2's generated config →
ZLIB not found. Purely the cmake version on PATH, **not** a CMakeLists bug.

**Permanent fix — shipped in `build.bat` (PR #225):** two steps run automatically
after `conan install`:
- **Step 2b (pin):** greps the conan cache for the pinned `cmake.exe`
  (`dir /s /b "%USERPROFILE%\.conan2\cmake.exe" | findstr "\p\bin\cmake.exe"` — no
  hash hardcoded) and prepends its dir to PATH so every `cmake` call resolves to
  3.31.x. Echoes the version; falls back to PATH cmake with a loud warning.
- **Step 2c (heal):** if `CMakeCache.txt`'s `CMAKE_COMMAND` isn't the conan copy
  (a shadow configure poisoned it), deletes JUST the cache file — clean
  reconfigure, compiled `.obj` survive.

So you no longer hand-recover. **If you ever do** need a manual recovery, invoke
the conan cmake by FULL PATH (e.g.
`C:\Users\nmora\.conan2\p\cmakefb7f5258f2df7\p\bin\cmake.exe`) so `cl`/`ifx`/`ninja`
stay on PATH from setup_env while cmake itself is 3.31.x.

> **Do NOT** prepend the conan cmake dir via `set "PATH=%CM%;%PATH%"` in a cmd
> one-liner — `%PATH%` expands at PARSE time (before setup_env runs) and wipes the
> vcvarsall additions (`cl` vanishes). Use full-path invocation.
>
> **Do NOT** try to auto-pin cmake inside `setup_env.bat` (attempted, abandoned):
> (1) editing the .bat with LF-only line endings makes cmd silently SKIP the block
> (cmd needs CRLF); (2) after the oneAPI `vars.bat` load, `for /f` and pipes fail
> to spawn child cmd. The `build.bat` approach above is the supported one. Simplest
> alternative permanent fix: deprioritise/uninstall system CMake 4.3.

---

## 3. MUMPS 5.5.1 — always pass `-Dintsize64=OFF` explicitly

`build.bat` builds MUMPS once into `mumps-install/` (from
`mumps-archive/mumps_src.tar.gz`, auto-downloaded from `mumps-solver.org` — the
old `mumps.enseeiht.fr` URL is dead — consumed by a `scivision/mumps` v5.5.1.5
clone). Configure: Ninja + ifx/cl, `-Darith=d`, `BUILD_SHARED_LIBS=OFF`,
`CMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded` (matches OpenSees `/MT`).

**The trap:** `build.bat` passes `-Dintsize64=OFF` **even though OFF is the
default**, because scivision/mumps' `option(intsize64 ...)` is **sticky in CMake's
cache**. If a previous configure ran with `ON`, then you delete `mumps-build/` and
reconfigure *without* the flag, `FetchContent` retains the old subbuild cache and
still applies `-i8`. Symptom: an **ABI split inside one struct** — C writes 4
bytes for `id.job`, Fortran reads 8 for `INSTANCE_NUMBER` — surfacing at runtime
as **"Instance Error 1 in DMUMPS_F77"** with `INSTANCE_NUMBER` a huge garbage
integer.

Always pass `-Dintsize64=OFF` explicitly, or wipe `mumps-build/` **and** its
`_deps/` subtree. With OFF: Fortran INTEGER = 4 bytes, C `MUMPS_INT` = `int`,
Intel MPI LP64 (32-bit handles), and OpenSees links the `lp64` ScaLAPACK variants
(`mkl_scalapack_lp64` / `mkl_intel_lp64` / `mkl_blacs_intelmpi_lp64`). Only flip
to ILP64 if matrices exceed ~2³¹ nonzeros — non-trivial (MUMPS + OpenSees
ScaLAPACK linkage + the ILP64 impi wrapper all have to flip together). If a user
reports "Instance Error 1" or garbage MUMPS solutions, check this is still OFF and
ScaLAPACK is linking `lp64`.

---

## 4. Test bootstrap — run pytest under CPython 3.12

There is no reliable persistent `.pth` on this box (the `opensees_venv` one is
3.11 and stale — it points at a non-existent `OpenSees_Compile\dist\bin` path).
Run tests with the 3.12 interpreter and bootstrap the DLL dir explicitly:

```python
import os, sys
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
import opensees as ops      # set LADRUNO_OPENSEES_QUIET=1 to mute the banner
```

- The MKL DLLs already sit in `dist/bin`; `os.add_dll_directory(DIST)` is enough
  — you do **not** need to source oneAPI setvars just to *run* (only to *build*).
- `opensees_venv` (Python 3.11) **cannot** import the 3.12-linked pyd
  ("DLL load failed"); that error is the 3.11-vs-3.12 ABI mismatch, not a missing
  MKL/MPI dll. Use `pythoncore-3.12-64\python.exe`; find others via `py -0p`.
- Run a battery:
  `& <py3.12> -m pytest tests\test_X.py -v -p no:cacheprovider`.
- **gmsh coexists**: `pip install gmsh` (→ 4.15.2, SDK wheel with bundled
  binaries) imports fine in the same 3.12 env and runs in ONE process with
  opensees (mesh, `gmsh.finalize()`, then build the ops model). This is the RAW
  gmsh SDK, not the `apeGmsh` wrapper.

---

## 5. Installer DLL-lock on upgrade ("DeleteFile failed; code 5")

Re-running the Inno installer to UPGRADE `C:\Program Files\Ladruno\OpenSees`
throws `DeleteFile failed; code 5. Access is denied` on a `bin\` DLL (seen on
`mkl_intel_thread.2.dll`). It is a **file lock, not a permissions/elevation
problem.**

**Root cause:** a previous installer wired `opensees_venv` with a *boot* `.pth`
that **runs `import opensees` on every interpreter startup**. So every Python in
that venv (VS Code's Black/isort/pylint `lsp_server.py`, linters, Jupyter
kernels) loads `opensees.pyd` + the install's MKL DLLs and pins them. VS Code's
Python extension respawns its formatter server within ~1s if killed — whack-a-mole.

**Diagnose:** `Get-Process | ? { $_.Modules.ModuleName -match 'mkl_intel_thread' }`
(no admin needed for your own session). Culprit is typically
`...\opensees_venv\Scripts\python.exe` running an `ms-python.*\lsp_server.py`.

**Fixes (most robust first):**
1. Close VS Code → installer **Retry** → reopen. Foolproof.
2. Keep VS Code open: a ~25s loop that `Stop-Process` any `python.exe` whose
   CommandLine matches `lsp_server.py` (stateless formatter — safe; **never** match
   ipykernel/Jupyter kernels, which hold user state), and click Retry in the window.
3. Windows restart (heavier).

**Permanent fix (do once):** convert the `opensees_venv` wiring from the boot
`.pth` (active `import`) to a *passive* `.pth` (adds the dist path to `sys.path`
but does NOT import on startup). Then VS Code's Python tools no longer auto-load
the install DLLs and upgrades stop getting locked. See `wire_venv_pth.py` /
`wire_pyenv.ps1` in `Ladruno_scripts`.

---

## 6. Building the all-features installer

The all-features installer is built off `ladruno` HEAD in a throwaway detached
worktree with a **short path** (e.g.
`C:\Users\nmora\Github\OpenSees_Compile\ladruno-build`) to dodge MAX_PATH, with
MUMPS seeded from an existing build to skip the ~20 min step. **`build.bat` builds
all 5 targets only when invoked with NO args** — passing target names restricts
to just those.
