---
wp: PR-877
title: "LadrunoProbeCommitScope: throwaway shadow material copies must not trip the real commit latch (WP concrete3d-hang-diagnosis, #877)"
pr: "#877"
date: 2026-10-10
---
### `LadrunoProbeCommitScope`: throwaway shadow material copies must not trip the real commit latch (WP concrete3d-hang-diagnosis, [#877](https://github.com/nmorabowen/OpenSees/pull/877))
- **Bites:** `LadrunoBrick` / `LadrunoQuad` commit hourglass-SHADOW material copies through the same `commitState()` as their real points; a refusal in a probe (K&R coarse, LadrunoQuad ssp, 0.6 mm) failed the real commit while every real Gauss point was healthy. An element wraps every probe commit in a `LadrunoProbeCommitScope` (RAII; `ladrunoNoteCommitRefusal()` is a no-op while one is open); `tests/test_ladruno_probe_commit_guard.py` parses the sources (brace-matched, comments/strings stripped) and fails if a `commitHgShadows()` call has no scope in scope. NDTest's direct `commitState` is scoped the same way; `wipe` and `revertToLastCommit` clear the counter / the material's `lastStatus`. The counters are plain non-atomic statics (serial commit).
