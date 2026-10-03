---
wp: LEGACY
title: "Copy-Item PRESERVES the source's LastWriteTime, so restoring a file from a backup copy leaves ninja convinced the object is up to date — the \"rebuild\" silently…"
legacy_seq: 318
---
### `Copy-Item` PRESERVES the source's LastWriteTime, so restoring a file from a backup copy leaves ninja convinced the object is up to date — the "rebuild" silently keeps testing the OLD code
- **Bites:** the standard A/B pattern for a negative control -- save a fixed source aside,
  `git checkout --` it to get the broken version, rebuild, measure, then `Copy-Item` the fixed
  version back and rebuild again. The **restore** build is a no-op: `Copy-Item` stamps the
  destination with the SOURCE file's timestamp (that of the aside copy, taken *before* the
  negative-control build), so the restored `.cpp` is OLDER than the `.obj` ninja produced from
  the broken one and ninja skips it. The binary you then test is still the broken build, and
  it reads as "my fix does not work". Hit exactly this on 2026-08-18: the pytest gate reported
  an already-fixed file as crashing again.
- **Why it is nastier than an ordinary stale build:** every honest signal says the restore
  worked. `git diff` shows the fix present, the build script exits 0, and the log even shows a
  relink (other targets moved). Only that file's own compile line is missing from the log,
  which is not a thing anyone looks for. It also inverts the usual reading of a red gate --
  the test is right and the binary is lying, so the instinct to go debug the test is wrong.
- **Workaround/status (2026-08-18):** after restoring a file from a copy, **touch it** --
  `(Get-Item path).LastWriteTime = Get-Date` -- then rebuild; or restore with `git checkout --`
  / `git stash pop`, which write fresh mtimes. Cheap verification: `grep` the build log for
  that file's own compile line, or check that the artifact's mtime is newer than the source's.
  `cp -p` and `robocopy` preserve timestamps by default too. Sibling of the stale-`.pyd` trap
  already recorded for `build.bat`.
