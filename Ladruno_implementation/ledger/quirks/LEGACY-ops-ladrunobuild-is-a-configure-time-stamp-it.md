---
wp: LEGACY
title: "ops.ladrunoBuild() is a CONFIGURE-time stamp — it LAGS after an incremental rebuild"
legacy_seq: 448
---
### `ops.ladrunoBuild()` is a CONFIGURE-time stamp — it LAGS after an incremental rebuild
- **Bites:** you edit C++, run `Ladruno_scripts\build.bat <targets>`, and the new
  binary reports the hash of an *older* commit. Every evidence run in WP-99 did
  this: the round-0 binary reported `bab19cfae` while `HEAD` was `c0c31f977`, and
  the round-1 binary reported `c0c31f977` while `HEAD` was `fa042bf51`. If you
  paste that into a PR as "the binary this was measured on", you have understated
  what you tested by one or more commits — and if you were checking *for* a stale
  binary, you would have concluded the opposite of the truth.
- **Why:** `CMakeLists.txt:200-207` captures the hash in an `execute_process`
  running `git log -1 --format=%H` — **at configure time**, into a cached
  `GIT_VERSION` that becomes a compile definition. An incremental `build.bat` run
  does not re-run CMake configure, so the cached value is reused no matter how
  many commits have landed since. The `.pyd`/`.exe` mtimes *are* fresh; only the
  stamp is stale.
- **Workaround/status (2026-09-14):** before an evidence run, force a
  reconfigure — `touch CMakeLists.txt` (or `Ladruno_scripts\build.bat clean`,
  which is the guaranteed way) — or state the lag explicitly and prove the
  binary behaviourally instead: assert on something the new code emits and the
  old code cannot (WP-99 used the new `Domain::commit() - N integration point(s)
  REFUSED this commit` line, which exists in neither of the two candidate older
  commits, plus the widened 6-slot `implexRefusals` response). The stamp is still
  the right first check for a *grossly* stale build (see the memory entry
  "ladrunoBuild provenance command"); it just cannot resolve one commit.
  See `Ladruno_internal/BUILD_GOTCHAS.md`.
