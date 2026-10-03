---
wp: LEGACY
title: "Copy-Item preserves the source's timestamp, so restoring a mutation backup can leave ninja compiling the MUTANT"
legacy_seq: 346
---
## `Copy-Item` preserves the source's timestamp, so restoring a mutation backup can leave ninja compiling the MUTANT

**Cost 2026-08-27, ADR-86 follow-up, ~10 minutes of chasing a "failing" test that was correct.**

The mutation-testing loop is: back the file up, apply a mutation, rebuild, run the tests, restore,
rebuild, re-run. Restoring with PowerShell's `Copy-Item` **preserves the source file's
`LastWriteTime`**, so the restored source carried the timestamp of when the backup was TAKEN —
which is *earlier* than the object file built from the mutant. Ninja compares mtimes, saw the
source as older than its `.obj`, and **skipped the recompile**. The suite then reported the
post-restore build as still failing, which reads exactly like "the fix does not work".

Measured: restored source at `19:10`, `ManzariDafalias.cpp.obj` at `19:10:44` from the mutant, and
`fastbuild` printed `FASTBUILD: OK` having compiled nothing. The link step still ran and refreshed
`OpenSeesPy.dll`'s mtime, so even the artifact timestamp looked current.

- **After restoring any file, `touch` it before rebuilding** — or restore with something that
  stamps the current time (`cat backup > file`, or the Write/Edit tools).
- **Confirm the recompile, do not infer it.** Grep the build log for the specific TU:
  `... | Select-String "ManzariDafalias"` must show a `Building CXX object` line. `FASTBUILD: OK`
  on its own means the *link* succeeded, not that your file was compiled.
- This is the same family as `86_ladruno_sanisand_handoff` §1's stale-binary trap and the
  edit-during-a-build trap above, and it presents identically: a green or red result about a tree
  that is not the one on disk.
