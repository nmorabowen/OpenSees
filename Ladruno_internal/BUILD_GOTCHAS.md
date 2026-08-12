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

**Permanent fix — DONE in the tool (#735), no longer a manual step.**
`wire_venv_pth.py` now generates a *passive* boot module. It still registers
search locations at startup (`sys.path`, `add_dll_directory`, process-local
`PATH`) — none of which LOAD anything — but the `openseespy` alias that used to
do an eager `import opensees` right there is now resolved lazily by a
`sys.meta_path` finder. The first thing to load `opensees.pyd` is the user's own
import, so VS Code's Python tools no longer pin the install DLLs.
**Re-run the wirer once per venv to pick it up:**

```
<venv>\Scripts\python.exe Ladruno_scripts\wire_venv_pth.py "C:\Program Files\Ladruno\OpenSees\bin" "C:\Program Files\Ladruno\OpenSees\openseesmp"
```

**Measure it with `tasklist /m opensees.pyd`, as a set difference against a live
baseline** — a bare venv Python must add NO new holder, while one whose script
imports opensees adds some (that contrast is the only thing that proves the probe
works at all). Two probes that LOOK right and silently lie:

* `Get-Process($pid).Modules` — enumerated **7** modules for a Python that had
  demonstrably imported the .pyd (its own stdout printed the Program-Files path),
  and reports a clean `0` held whether or not the file is pinned. It will tell you
  the lock is gone when it is not.
* `tasklist /m X /fi "PID eq N"` — misses, because the venv launcher re-spawns the
  real interpreter under a **different PID** than `Start-Process` hands back. Run
  it unfiltered and diff the PID set.

The same change retires the companion `LEDGER_quirks` entry ("An INSTALLED
Ladruno hijacks `import opensees`"): with nothing imported at startup, a
worktree's own `sys.path.insert(0, dist/bin)` wins again, and
`LADRUNO_OPENSEES_BIN` is demoted from crutch to override.

**It does NOT save you from a process that legitimately imported the engine** —
an analysis run, a pytest session, its `multiprocessing` workers. The passive
`.pth` stops *bare* interpreters from holding the DLLs; a real user of OpenSees
still holds them, correctly.

### 5b. `/SUPPRESSMSGBOXES` turns "could not close applications" into a SILENT abort

A scripted upgrade (`setup.exe /VERYSILENT /SUPPRESSMSGBOXES`) returned **exit 5**
and changed nothing, twice, while `install.log` showed only the PREVIOUS run's
lines — easy to misread as "the installer is broken". It was not. Inno's
RestartManager found six Python processes using files it was about to replace,
failed to shut them down, and `/SUPPRESSMSGBOXES` answered the resulting
Abort/Retry/Ignore box with its default, **Abort**:

```
RestartManager found an application using one of our files: Python   (x6)
Shutting down applications using our files.
Some applications could not be shut down.
Defaulting to Abort for suppressed message box (Abort/Retry/Ignore)
User canceled the installation process.
```

**Read `/LOG=<file>`, not `install.log`,** to see this: `install.log` is written by
our own `[Code]` at `ssPostInstall`, so an install that aborts BEFORE that never
writes a line, and its stale contents describe an earlier run. Exit 5 = user
cancelled during install; with `/SUPPRESSMSGBOXES` that "user" is the default
button.

Nothing is damaged when this happens — the abort is before any file is replaced,
and the existing install keeps working (verify with `ladrunoBuild`). Options: wait
for the runs to finish; `/NOCLOSEAPPLICATIONS` (then locked files fail
individually instead); or install elsewhere with `/DIR=`.

**If you do test-install to a scratch `/DIR=`, repair the uninstall key
afterwards.** `AppId` is FIXED in `installer.iss`, so every install — scratch ones
included — overwrites the machine's single Add/Remove Programs entry under
`HKLM\...\Uninstall\{8C8E2E87-...}_is1`. Test-installing to a temp folder silently
repoints `InstallLocation` / `UninstallString` / `QuietUninstallString` /
`Inno Setup: App Path` at that folder, and deleting it then leaves a dangling
entry and a real install with no entry at all. Same key also makes the wizard's
default destination "sticky", so the next interactive run offers the temp path.

---

## 6. Building the all-features installer

The all-features installer is built off `ladruno` HEAD in a throwaway detached
worktree with a **short path** (e.g.
`C:\Users\nmora\Github\OpenSees_Compile\ladruno-build`) to dodge MAX_PATH, with
MUMPS seeded from an existing build to skip the ~20 min step. **`build.bat` builds
all 5 targets only when invoked with NO args** — passing target names restricts
to just those.

---

## 7. Where a new subsystem's sources live — two sanctioned patterns

Two build-organization patterns exist in the fork; **pick by dependency weight
and maturity, and know the CI consequence of each.** This is a real fork
convention question that came up during the ADR-1000 (LadrunoCMS) review.

**Pattern A — integrate-and-always-build (the default; use this unless B's
triggers apply).** Drop the new `.cpp`/`.h` into the nearest existing `SRC/`
subdir and add them to that subdir's `target_sources(OPS_<Existing> ...)`. The
classes compile unconditionally into a stock `OPS_*` object library; only
platform/accelerator *linkage* is guarded (a CMake option flips a macro, the
heavy path degrades to a stub). Precedents:
- **FEAST** (ADR-43) → `SRC/system_of_eqn/eigenSOE/`, always compiled into
  `OPS_SysOfEqn`; only the MKL FEAST kernel is behind `LADRUNO_MKL_FEAST_LINUX`
  (OFF) and stubs out where MKL is absent. `eigen -feast` always exists.
- **Contact engine** → `SRC/domain/contact/` + `SRC/analysis/handler/`, always
  compiled into `OPS_Domain` / `OPS_Analysis`, no CMake option.

**Pattern B — isolated opt-in library (`OFF` by default).** A new top-level
`add_library(OPS_<Name> STATIC)` created only `if(<OPTION>)`, its own
`CMakeLists.txt` pulled in by a gated `add_subdirectory`, linked into only the
targets that can use it. **Sole precedent: `OPS_LadrunoCMS`** (ADR-1000):
`option(LADRUNO_CMS ... OFF)`, its own `SRC/system_of_eqn/ladrunoCMS/CMakeLists`,
linked into MP/PyMP only, with `_LADRUNO_CMS` gating both the sources and the
`eigen -ladrunoCMS` interpreter wiring.

**Use Pattern B only when ALL of these hold** (CMS is the model): the subsystem
(i) hard-requires heavy deps that the default build does NOT already carry (CMS:
MPI + MUMPS + ScaLAPACK + METIS as `PUBLIC` link reqs, `FATAL_ERROR` if absent —
so it cannot degrade to a stub the way FEAST's MKL path does); (ii) is usable in
only a subset of targets (CMS is MP-only); and (iii) is large and pre-shipment
("remediation" maturity). Otherwise use Pattern A — a stubbable accelerator or a
small always-usable feature does **not** justify a separate library.

**The cost you take on with Pattern B — and must pay down:** an `OFF`-by-default
library is **never compiled by the default CI build**, so its entire core rots
without a compile signal. When ADR-1000 shipped, only `LadrunoCMSOptions.cpp`
(g++-compiled by the `zone_a` core test) had any CI coverage; the other ~6k
lines (Hierarchy, Lanczos, MUMPS wrapper, SubspaceRefiner, EigenSOE/Solver) had
none. **Any Pattern-B subsystem MUST ship with a dedicated CI lane that
configures `-D<OPTION>=ON` (+ its `_BUILD_TESTS`) and at least compiles the
library + runs its standalone checks.** For CMS this is the nightly self-hosted
Zone-B job (it has the oneAPI MPI+MUMPS toolchain; the Ubuntu Zone-A runner does
not). Build it locally the same way: `build.bat` configures with
`-DMUMPS_DIR`/`-DMUMPS_INCLUDE_DIR` already; add `-DLADRUNO_CMS=ON
-DLADRUNO_CMS_BUILD_TESTS=ON` to the configure and the ~7 CMS TUs compile
incrementally (the rest of the tree is cached).

## 8. oneAPI 2026 ships NO Fortran → three cascading build failures (ADR-75)

**Symptom (any one of three, depending on how far you get):**

1. `CMAKE_Fortran_COMPILER: ifx ... is not a full path and was not found in the PATH`
   — at the **MUMPS configure** (build.bat Step 1).
2. `LINK : fatal error LNK1104: cannot open file 'ifconsol.lib'`
   — when **linking `OpenSeesMP` / `OpenSeesSP` / `OpenSeesPyMP`**.
3. `CMAKE_Fortran_COMPILER: ...\compiler\latest\bin\ifx.exe is not a full path to an
   existing compiler tool` — at the **main OpenSees configure**.

**Root cause — one event, three faces.** oneAPI **2026.x ships no Fortran at all**: Intel split it
into a separate package, so `compiler\2026.1` has neither `ifx.exe` nor `ifconsol.lib`. An update
re-points the **`compiler\latest` junction**, and everything downstream breaks:

- no `ifx` ⇒ failure (1);
- no `ifconsol.lib` ⇒ failure (2), because **MKL's ScaLAPACK/BLACS import libraries still carry a
  `/DEFAULTLIB:ifconsol` directive**. **Serial targets never link ScaLAPACK, so they keep building
  fine** — which is precisely why this stays invisible until somebody builds an MP target;
- CMake caches the **resolved absolute** compiler path, so a cache written while `latest` was 2025.x
  now names a file that no longer exists ⇒ failure (3). Note **`-DCMAKE_Fortran_COMPILER=ifx` does
  NOT fix (3)** — the cached absolute path wins.

This was observed **live**: an update flipped `latest` 2025.1 → 2026.1 *between two builds in one
session* (the earlier log shows `Fortran compiler identification is IntelLLVM 2025.1.1` from
`latest\bin`, the later one finds nothing there).

**Fixes, all shipped (`Ladruno_scripts/`):**
- `setup_env.bat` §2a — prefer `compiler\latest`, but if it has **no `ifx.exe`** fall back to the
  newest installed compiler that does and source **that** `vars.bat`. Hard error (with a pointer to
  the separate Fortran package) if none has it.
- `setup_env.bat` §2b — independently, if the chosen compiler has no `ifconsol.lib`, prepend the lib
  dir of one that does. Kept separate because a compiler *can* have `ifx` but not `ifconsol`.
- `build.bat` §2d — delete `CMakeCache.txt` when its `CMAKE_Fortran_COMPILER` no longer exists on
  disk. Compiled objects survive (ninja re-runs only what its hashes changed).

**Diagnosis one-liner** (run after `setup_env.bat`, in `cmd`):
```
for %P in (ifx.exe) do @echo %~$PATH:P
for %P in (ifconsol.lib) do @echo %~$LIB:P
```
Both must print a real path. If `ifx` resolves under `2025.x` while `latest` is `2026.x`, the
fallback is doing its job.

**Note for these scripts:** `setup_env.bat` must **not** `setlocal` (it exports the environment to
its caller), so **delayed expansion `!VAR!` is unavailable** — write these blocks flat with labels,
never as a parenthesized read-after-write.

---

## 9. MKL bumps its DLL SONAME → a hardcoded copy list ships nothing (ADR-75)

Companion to §8. That one is about the oneAPI **compiler** package moving under
you; this one is about the **MKL runtime** doing the same thing, and it is worse
because it produces a *successful-looking build* that fails on the user's machine
rather than an error on yours.

Intel version-stamps the MKL runtime DLLs and **bumps the stamp across
releases**. Going 2025.1 → 2026.1:

| base | oneMKL 2025.1 | oneMKL 2026.1 |
|---|---|---|
| `mkl_core`, `mkl_intel_thread`, `mkl_def`, `mkl_avx2`, `mkl_avx512`, `mkl_mc3` | `.2.dll` | **`.3.dll`** |
| `mkl_scalapack_lp64`, `mkl_blacs_intelmpi_lp64` | `.2.dll` | `.2.dll` (unchanged) |

Note the two groups move **independently** — you cannot assume one suffix for
all of MKL, which is exactly why globbing per base name is the only safe form.

`build.bat` used to mirror these into `dist\bin\` and `dist\openseesmp\` by
**full filename**, each behind an `if exist` guard. After the upgrade every
guard missed, so the copy loop **silently copied nothing** and left the
*previous* version's `.2.dll` files sitting in `dist\`. The shipped package then
either:

- **fails to load** — `mkl_core.3.dll not found`, because the freshly linked
  import libraries reference the `.3` SONAME; or
- **loads the wrong MKL** — the stale `.2` DLLs satisfy an older binary, so you
  ship something you never tested.

No error, no warning: the `if exist` guard swallowed the whole failure.

**Mitigation** (`build.bat`): mirror by **base name**, globbing `<base>.*.dll`
via `MKL_RUNTIME_DLLS`, with `del /q dist\...\mkl_*.dll` before copying so a
bump cannot leave two generations behind. `libiomp5md.dll` gets the same
treatment — prefer `compiler\latest\bin`, else the newest compiler dir that
actually has it.

**Rules of thumb:**
1. Never name an Intel runtime DLL by its stamped filename in a script.
2. Never guard a **required** copy with a bare `if exist` — that converts a
   missing dependency into a silent packaging bug. Guard optional things only.
3. After any MKL upgrade, verify the shipped artifact, not just the build:
   `dist\bin` and the installer must contain the CPU-dispatch kernels
   (`mkl_def` / `mkl_avx2` / `mkl_avx512` / `mkl_mc3`). Missing those does not
   fail the build — MKL just silently falls back to a generic kernel.

### Rollback

Both MKL generations stay on disk. To revert, repoint the symlink at
`...\Intel\oneAPI\mkl\latest` back to `...\mkl\2025.1` (admin shell) and
rebuild — the scripts pick up whatever `latest` resolves to.

### A note on measured performance

The 2025.1 → 2026.1 upgrade is **not** justified by a measured speedup. Observed
cross-session variance on the dev box (~33% at identical settings — UmfPack
n=20 @ 4 threads measured 8.29 s in one session and 11.00 s in another) exceeds
every 2025.1-vs-2026.1 delta recorded, so the two versions are indistinguishable
on the evidence gathered. Treat any single-session A-vs-B MKL comparison here as
uninformative unless both are benched back-to-back in one session.

One mechanism *was* settled while chasing this: UmfPack speeds up 1.35x from
1→4 MKL threads, confirming it runs inside threaded MKL BLAS, whereas MKL
PARDISO uses its own supernodal kernels and is untouched by BLAS-layer work.

## 10. `LADRUNO_CMS=ON` never linked — static + dynamic MKL in one target

**Symptom.** Any `LADRUNO_CMS=ON` build dies at the *link* step (the whole
library compiles first, which makes it look like a late/unrelated failure):

```
mkl_core.lib(mkl_verbose.obj) : error LNK2005: mkl_serv_verbose_mode already
    defined in mkl_intel_thread_dll.lib(mkl_intel_thread.3.dll)
mkl_core.lib(mkl_xerbla.obj)  : error LNK2005: mkl_serv_xerbla already defined ...
ladruno_cms_mumps_check.exe : fatal error LNK1169: one or more multiply defined
```

**Cause.** `SRC/system_of_eqn/ladrunoCMS/CMakeLists.txt` linked BOTH
`${SCALAPACK_LIBRARIES}` — the **static sequential** MKL stack the root
`CMakeLists` builds for the parallel targets (`mkl_scalapack_lp64` +
`mkl_intel_lp64` + `mkl_sequential` + `mkl_core` + `mkl_blacs_intelmpi_lp64`) —
and `${LAPACK_LIBRARIES}`, which on any box where `find_package(LAPACK)`
resolves to oneAPI MKL expands to the **dynamic** interface
(`mkl_intel_lp64_dll` + `mkl_intel_thread_dll` + `mkl_core_dll` + `libiomp5md`).
MKL does not support mixing the static and dynamic layers in one link. The two
sets collide on `mkl_serv_*`.

`OpenSeesSP`/`OpenSeesMP` link `${SCALAPACK_LIBRARIES}` and **not**
`${LAPACK_LIBRARIES}` — that is why only the CMS library tripped this, and why
it propagated (the offending list was `PUBLIC`) to any MP target built with
`LADRUNO_CMS=ON`. Whether it bites is environment-dependent: it needs a box
where `find_package(LAPACK)` finds MKL, so the same tree can link on one machine
and fail on another. This is the concrete reason the ADR-1000 P3/P4 execution
plans were authored "without a runtime build".

**Fix.** Drop `${LAPACK_LIBRARIES}` from `OPS_LadrunoCMS` and re-add it only
`if(NOT MKL_FOUND)` (the non-MKL box, where `SCALAPACK_LIBRARIES` is
user-supplied and carries no BLAS). Same rule for any future Pattern-B library
(#7): **take the BLAS/LAPACK layer from exactly one source, and match whichever
one the targets you link into already use.**

**A second, independent defect in the same lane.** `build.bat`'s CMS run loop
echoed `--- running ..._check (mpiexec -n 4) ---` *inside* a parenthesized
`for`-body. `cmd` closed the block on that inner `)` and failed with `--- was
unexpected at this time` — **after** a fully successful build, so the lane
reported failure having never run a check. Escape inner parens as `^( ^)` (this
is trap #5 in the family of batch traps in #1). Both defects had to be fixed
before a single CMS numerical check had ever executed in this build tree.

## 11. Never `sed -i` a `.bat` file — it strips CRLF and cmd silently miscompiles

`sed -i` (Git Bash) rewrites the whole file with **LF-only** line endings. `cmd.exe`
requires **CRLF** in batch files: with bare LFs, single-line commands still work,
but multi-line constructs silently fall apart. A one-character edit to a `for`
list produced:

```
'do' is not recognized as an internal or external command,
'not' is not recognized as an internal or external command,
```

and a stray interactive Python prompt — because `for %%C in (...) do (` and the
`if not defined ...` lines were no longer parsed as blocks. The `git diff` looked
perfect (two changed characters); nothing in it hints at the real damage. `git`
also warns `LF will be replaced by CRLF the next time Git touches it`, which is
easy to dismiss as routine noise — on a `.bat` it is the actual error message.

**Rule:** edit `.bat` files with the Edit tool (it preserves existing endings) or
a CRLF-aware editor, never `sed -i`/`tr`. To repair one:

```powershell
$t = [IO.File]::ReadAllText($p) -replace "`r`n","`n" -replace "`n","`r`n"
[IO.File]::WriteAllText($p, $t, (New-Object Text.UTF8Encoding($false)))
```

Verify with a bare-LF count (must be 0), not by eye — and re-run the script,
because a `git diff` cannot show you this class of breakage.

## 12. Running an MP binary from the raw build tree — three DLL traps (ADR-1000 §32)

The curated installer layout resolves its own runtime. A freshly built
`build\build\Release\OpenSeesMP.exe` does not, and fails in three different ways
that all look like "the build is broken". Each one cost time in the ADR-1000 §32
session; none is a build problem.

**a. `c0000135` on every rank = a missing DLL, and it is usually the MKL
generation.** `dumpbin /DEPENDENTS OpenSeesMP.exe` is the fast answer. This fork
currently imports `mkl_core.3.dll` / `mkl_intel_thread.3.dll` — the **2026.1**
generation. `mkl\2025.1\bin` ships `mkl_core.2.dll`, so sourcing the 2025.1 MKL
vars gets you a clean environment with the wrong SONAME (the §9 failure mode,
seen from the running end instead of the packaging end). Source
`mkl\2026.1\env\vars.bat`.

**b. `C:\Program Files\Ladruno\OpenSees\bin` is on the machine PATH and shadows
the build tree.** The *installed* distribution ships its own `libiomp5md.dll`
plus MKL/MUMPS DLLs. A raw build-tree exe picks them up and aborts with:

```
OMP: Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already
initialized.
```

Do **not** reach for `KMP_DUPLICATE_LIB_OK=TRUE` — that is the documented unsafe
workaround and it leaves you measuring an unknown mixture of two runtimes. Strip
the installed directory from `PATH` for the run instead. The abort is the benign
outcome here; the dangerous one is a *silent* stale-DLL shadow when the versions
happen to be compatible, which would quietly invalidate any performance number.

**c. `PMPI_Init` fails with `MPIDI_OFI_mpi_init_hook`** when `mpi\latest\bin` is
on `PATH` but the rest of the Intel MPI environment is not — libfabric and
`I_MPI_ROOT` come from `mpi\latest\env\vars.bat`, not from the `bin` directory.
Call the `vars.bat`, do not hand-assemble `PATH`.

Plus the non-DLL trap: **`TCL_LIBRARY` must point at the conan Tcl**
(`%USERPROFILE%\.conan2\p\b\tcl*\p\lib\tcl8.6`) or the interpreter starts with no
commands at all and the first line of the deck dies with
`invalid command name "wipe"`.

A working recipe (PowerShell), which imports the real component environments
rather than guessing at `PATH`:

```powershell
$O = "C:\Program Files (x86)\Intel\oneAPI"
cmd /c "call `"$O\compiler\2026.1\env\vars.bat`" >nul & call `"$O\mkl\2026.1\env\vars.bat`" >nul & call `"$O\mpi\latest\env\vars.bat`" >nul & set" |
  ForEach-Object { if ($_ -match '^([^=]+)=(.*)$') { Set-Item "env:$($matches[1])" $matches[2] -EA SilentlyContinue } }
$env:PATH = (($env:PATH -split ';') | Where-Object { $_ -notmatch 'Ladruno\\OpenSees\\bin' }) -join ';'
$env:TCL_LIBRARY = "C:\Users\nmb\.conan2\p\b\tcl1fa6686758830\p\lib\tcl8.6"
```
---

## 13. Building on a non-VS2022 box (VS 2026 / v18 verified 2026-07)

Everything above assumes the original VS **2022** machine. The fork now also
builds on VS **2026 (v18, MSVC 19.5x)**. Three deltas were needed; two are baked
into the scripts, one is a per-shell env var:

1. **`setup_env.bat` finds vcvars64 via `vswhere`** (not a hardcoded 2022 path):
   `vswhere -latest -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64
   -property installationPath`, falling back to the old 2022 Community path. Works
   across editions/years unchanged.
2. **The conan profile's MSVC version is env-overridable.**
   `opensees-msvc-static.profile` renders
   `compiler.version={{ os.getenv("LADRUNO_MSVC_VERSION", "194") }}`. Default `194`
   preserves the VS 2022 box. On a VS 2026 box **`set LADRUNO_MSVC_VERSION=195`**
   before `build.bat`, else conan dies in the zlib `generate()` step with
   `VS non-existing installation: Visual Studio 17` (it maps 194→VS17, which isn't
   installed). Confirm the target maps: `conan profile detect` prints
   `Found msvc 18` / `compiler.version=195` on a v18 box.
3. **Toolchain-on-PATH quirks on a fresh winget bootstrap:**
   - System CMake from `winget Kitware.CMake` is **4.x** — the §2 hazard. Install
     the pinned line via `py -3.12 -m pip install --user "cmake~=3.31.0"`; it lands
     in `%APPDATA%\Python\Python312\Scripts` (the dir `setup_env` prepends for
     conan), so it wins for the pre-conan MUMPS step too.
   - `winget Ninja-build.Ninja` installs to
     `%LOCALAPPDATA%\Microsoft\WinGet\Packages\Ninja-build.Ninja_*\` with **no**
     Links symlink → not on PATH. Prepend that dir before `build.bat`.

With `LongPathsEnabled=1` (default on Win11) the deep `.claude\worktrees\...`
build tree is fine — the `CMAKE_NINJA_FORCE_RESPONSE_FILE=ON` flag in `build.bat`
already covers the command-line-length side.
