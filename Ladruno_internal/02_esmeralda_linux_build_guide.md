---
title: esmeralda Linux build guide
project: Ladruno
tags:
  - internal
  - build
  - linux
---

# Building & running ladruno OpenSees on esmeralda (Linux)

**Status:** proven 2026-06-11 — full build from `ladruno` HEAD (`5865f22`) compiled
clean and passed a physics smoke test (elastic truss, `nodeDisp` = PL/EA exact).

`esmeralda` is an Ubuntu 22.04 Linux box reachable as `ssh esmeralda` (key-based,
`BatchMode` works, so scripts can drive it non-interactively). This guide records
what's on the machine, where the built binaries live, how to run them, and how to
rebuild when `ladruno` moves.

---

## 1. The machine

| Item | Value |
|---|---|
| OS | Ubuntu 22.04.1 LTS (Jammy), kernel 6.8 |
| CPU / RAM | 24 cores / 31 GB |
| Home | `/mnt/deadmanschest/nmorabowen` (NVMe, ~710 GB free) |
| C/C++/Fortran | gcc, g++, gfortran 11.4 (gcc-12 also installed) |
| CMake | 3.22.1 (system) — above our 3.16 floor, **below 3.23** (see §5) |
| MPI | OpenMPI 4.1.2 + `libopenmpi-dev` (`mpicc`/`mpicxx`/`mpirun` on PATH) |
| BLAS/LAPACK | apt MKL 2020 packages, OpenBLAS, reference LAPACK — all present |
| Tcl | tcl8.6 + `tcl8.6-dev` (system), plus the conan tcl 8.6.11 used by the build |
| Python | 3.10.12 + `python3-dev` (configure found Development/Embed for OpenSeesPy) |
| Conan | **not system-installed** — lives in the venv below |

## 2. Where everything lives on esmeralda

Everything is under `~/ladruno_build_test/`:

| Path | What it is |
|---|---|
| `~/ladruno_build_test/OpenSees/` | clone of `github.com/nmorabowen/OpenSees`, branch `ladruno` |
| `~/ladruno_build_test/OpenSees/build/Release/OpenSees` | **the built binary** (~40 MB, Release) |
| `~/ladruno_build_test/opensees.sh` | **run wrapper** — sets `TCL_LIBRARY` and execs the binary (use this) |
| `~/ladruno_build_test/conan_venv/` | Python venv holding Conan 2.29 (`./conan_venv/bin/conan`) |
| `~/ladruno_build_test/conan_install.log` | dependency-install log from the proven build |
| `~/ladruno_build_test/cmake_configure.log` | configure log |
| `~/ladruno_build_test/build_opensees.log` | full compile log (ends `BUILD_OK`) |
| `~/.conan2/p/tcl*/p/lib/tcl8.6/` | conan Tcl runtime scripts (`init.tcl`) the binary needs |

## 3. How to run it

### The wrapper (recommended)

```bash
ssh esmeralda
~/ladruno_build_test/opensees.sh model.tcl        # batch
~/ladruno_build_test/opensees.sh                  # interactive OpenSees > prompt
```

The wrapper exists because of one gotcha: the binary links conan's Tcl 8.6.11,
which can't find its `init.tcl` startup scripts on its own. Run bare, it dies with
*"Can't find a usable init.tcl"*. The wrapper locates the scripts in the conan
cache and exports `TCL_LIBRARY` before exec'ing the real binary.

### Running the binary directly

If you'd rather not use the wrapper (e.g. inside another script):

```bash
export TCL_LIBRARY=$(dirname $(find ~/.conan2/p -path '*tcl*' -name init.tcl | head -1))
~/ladruno_build_test/OpenSees/build/Release/OpenSees model.tcl
```

### One-shot from Windows (no interactive login)

```powershell
ssh esmeralda '~/ladruno_build_test/opensees.sh /path/to/model.tcl'
```

Pair with `scp` to push the model and pull results:

```powershell
scp model.tcl esmeralda:/tmp/
ssh esmeralda '~/ladruno_build_test/opensees.sh /tmp/model.tcl'
scp esmeralda:/tmp/results.out .
```

### Sanity check

```bash
ssh esmeralda '~/ladruno_build_test/opensees.sh /tmp/smoke.tcl'
```

`/tmp/smoke.tcl` (left on the machine) is a 1-element elastic truss; expected
output is `DISP=0.001` (= PL/EA with P=1, L=1, E=1000, A=1). The splash also
prints the git SHA the binary was built from — use it to confirm which commit
you're running.

## 4. How it was built (the proven recipe)

This mirrors the Zone-A Linux CI job (`.github/workflows/ladruno.yml`):
Conan supplies HDF5 1.14.0 / Tcl 8.6.11 / zlib 1.3.1 / Eigen 3.4.0, CMake
consumes the conan toolchain.

```bash
# one-time setup
mkdir -p ~/ladruno_build_test && cd ~/ladruno_build_test
python3 -m venv conan_venv
./conan_venv/bin/pip install conan                      # Conan 2.x
export PATH=$HOME/ladruno_build_test/conan_venv/bin:$PATH
git clone --depth 1 --branch ladruno https://github.com/nmorabowen/OpenSees.git
cd OpenSees
conan profile detect --force                            # gcc 11, gnu17, Release

# dependencies (~1 min — all four come as prebuilt conan-center binaries)
conan install . --build=missing

# configure + build (~30 min wall at -j20)
cmake -S . -B build/Release \
  -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/Release --target OpenSees -j20
```

Configure highlights worth knowing: profiler HDF5 writer enabled (HDF5 1.14.0),
Python 3.10 Development/Embed found (so the `OpenSeesPy` target is configured),
no ZLIB drama (the Windows CMake-4.3-shadowing quirk doesn't exist here).

## 5. Gotchas

- **CMake 3.22 < 3.23 → no presets.** Conan prints `cmake --preset conan-release`
  as the suggested command; that fails on esmeralda. Use the explicit
  `-DCMAKE_TOOLCHAIN_FILE=...` form shown above.
- **`TCL_LIBRARY` is mandatory at runtime** (§3). Symptom: *"Can't find a usable
  init.tcl"* then `invalid command name "model"`.
- **The build looks hung at 98%.** It isn't. `OPS_Material` ends with a few huge
  translation units compiling single-threaded for ~10+ minutes after the parallel
  phase drains. Check `pgrep -c cc1plus` — 1 running compiler means it's working.
- **Conan is venv-local.** A bare `conan` on a fresh shell won't resolve; either
  `export PATH=$HOME/ladruno_build_test/conan_venv/bin:$PATH` or call
  `~/ladruno_build_test/conan_venv/bin/conan` directly.

## 6. Rebuilding when `ladruno` moves

```bash
ssh esmeralda
cd ~/ladruno_build_test/OpenSees
git fetch origin ladruno && git reset --hard origin/ladruno
export PATH=$HOME/ladruno_build_test/conan_venv/bin:$PATH
conan install . --build=missing      # no-op unless deps changed
cmake --build build/Release --target OpenSees -j20
```

Incremental rebuilds only recompile what changed; the 30-minute figure is for a
cold build. (The clone is `--depth 1`, so `git fetch` + `reset --hard` is the
update path — don't expect full history locally.)

## 7. Not yet built (but expected to work)

- **`OpenSeesPy`** — the target is configured (Python 3.10 dev found); build with
  `cmake --build build/Release --target OpenSeesPy -j20`. Note the resulting
  module imports under esmeralda's Python 3.10, not the Windows 3.12 test env.
- **`OpenSeesSP` / `OpenSeesMP`** — OpenMPI 4.1.2 + dev headers are present, so
  the parallel interpreters are plausible on those 24 cores. Untested; expect to
  pass the MPI toolchain through conan/CMake flags and validate separately.
