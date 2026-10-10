---
wp: LEGACY
title: "wp/build-scripts-unify -- 4 vanilla row(s)"
files: ["`makeWIN.bat`", "`makeMac.sh`", "`OpenSeesAWS-Ubuntu22.04.sh`", "`conanfile2.py`"]
table: "main"
legacy_seq: [1, 2, 3, 4]
---
| `makeWIN.bat` | **DELETED.** Upstream's Windows build recipe (`conan install` + a hand-written `cmake` line into `build/Release`). Kept drifting from the fork's real chain — stale MUMPS path (`..\..\mumps\build`), lp64/static MKL against the fork's ilp64/dynamic, no `dist/` staging, no installer — while sitting in the repo root where it was the first thing anyone found. Agents repeatedly built through it (or copied its `build/Release` path) and then tested a stale or differently-linked binary. Superseded entirely by `Ladruno_scripts/build.bat`. | wp/build-scripts-unify |
| `makeMac.sh` | **DELETED.** Upstream macOS build recipe. The fork builds on Windows (`Ladruno_scripts/build.bat`) and Linux/Esmeralda (`Ladruno_internal/02_esmeralda_linux_build_guide.md`); this was an unmaintained third entry point with no fork coverage. | wp/build-scripts-unify |
| `OpenSeesAWS-Ubuntu22.04.sh` | **DELETED.** Upstream AWS/Ubuntu provisioning + build script. Same reason: an unmaintained parallel build recipe. The fork's Linux path is the Esmeralda guide. | wp/build-scripts-unify |
| `conanfile2.py` | **DELETED.** Upstream's second Conan 2 recipe, sitting next to the `conanfile.py` that `build.bat` actually uses. Two peer-looking conanfiles with no marker saying which is live. | wp/build-scripts-unify |
