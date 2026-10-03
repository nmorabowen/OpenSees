---
wp: LEGACY
title: "A running analysis HOLDS dist/bin/opensees.pyd, so a rebuild relinks the DLL, SILENTLY fails to copy it, and exits 0 -- you then test the OLD binary"
legacy_seq: 363
---
## A running analysis HOLDS `dist/bin/opensees.pyd`, so a rebuild relinks the DLL, SILENTLY fails to copy it, and exits 0 -- you then test the OLD binary

**Measured 2026-09-06, ADR-92 P1.**

`build.bat` ends with `copying OpenSeesPy.dll -> opensees.pyd`. If ANY process has
that `.pyd` loaded -- a background leg run with `LADRUNO_DIST_BIN` pointing at that
worktree's `dist/bin`, a stray `python -c "import opensees"` -- Windows refuses the
overwrite. The build **prints the copying line, exits 0, and leaves the old
`opensees.pyd` in place**. Measured: `build/build/Release/OpenSeesPy.dll` at
`01:53:37`, `dist/bin/opensees.pyd` still at `01:13:43`, forty minutes stale, with
three `python3.12` processes holding it.

- **Why it is dangerous rather than annoying:** every downstream gate then runs the
  PREVIOUS binary and attributes its behaviour to the change you just made. In this
  case the next step would have been a BVP gate on a freshly fixed guard, and the
  answer would have been the unfixed one, reported with a build hash that looked
  plausible because it *was* a real build of this branch -- just the wrong one.
- **Detection, and do it every time:** compare `ops.ladrunoBuild()` against
  `git rev-parse HEAD`, AND compare the mtimes of `build/build/Release/OpenSeesPy.dll`
  and `dist/bin/opensees.pyd`. The hash alone is not enough -- a stale `.pyd` from an
  earlier commit on the SAME branch reports a hash that passes a casual glance.
- **Fix:** finish or kill everything holding the module, then copy the DLL over the
  `.pyd` by hand (no need to rebuild -- the link already succeeded), or rebuild once
  the lock is clear.
- Related: the installer's DLL-lock entry, and `ladruno-build-targets-exe-vs-pyd`
  (building only `OpenSees` leaves the `.pyd` stale, which fails the same way and is
  ALSO invisible unless the mtime is checked).
