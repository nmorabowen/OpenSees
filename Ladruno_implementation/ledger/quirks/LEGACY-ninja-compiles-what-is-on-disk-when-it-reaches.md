---
wp: LEGACY
title: "Ninja compiles what is on disk WHEN IT REACHES THE TU, so editing source during a build silently mixes trees"
legacy_seq: 343
---
## Ninja compiles what is on disk WHEN IT REACHES THE TU, so editing source during a build silently mixes trees

**Cost 2026-08-27, ADR-86 PR-3, ~20 minutes.** A full build was launched from committed HEAD to
establish a baseline, and source edits were made while it ran, on the assumption that the
material TUs had already compiled (the log tail showed `1964/1968`). They had not — ninja
interleaves targets, and `LadrunoSANISAND.cpp` compiled *after* the edit. The resulting binary
was HEAD **plus some** of the edits, and the "baseline" battery run against it reported a
failure that looked like a real defect on `ladruno` HEAD.

This is the inverse of `86_ladruno_sanisand_handoff` §1's stale-binary trap and it reads exactly
the same from the outside: a green or red result about a tree that never existed.

- **Never edit `SRC/` while a build is running.** A progress counter is not a per-file cursor.
- The tell is a warning or behaviour in the run that is **newer than the build you launched**.
- For iterating on one material, skip `build.bat` entirely:
  `cmake --build build\build\Release --target OpenSeesPy -j 8` then copy `OpenSeesPy.dll` to
  `dist\bin\opensees.pyd` — **~30 s** against ~8 minutes, because it skips the CMake configure
  and the MUMPS step (which re-runs on this box even for a one-line C++ change).
