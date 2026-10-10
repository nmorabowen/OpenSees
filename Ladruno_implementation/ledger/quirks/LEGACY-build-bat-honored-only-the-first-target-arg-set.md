---
wp: LEGACY
title: "build.bat honored only the FIRST target arg (set MODE=%1) — build.bat OpenSeesPy OpenSeesPyMP silently built ONLY serial. And a parallel A/B compare can falsel…"
legacy_seq: 113
---
### `build.bat` honored only the FIRST target arg (`set MODE=%1`) — `build.bat OpenSeesPy OpenSeesPyMP` silently built ONLY serial. And a parallel A/B compare can falsely PASS on TWO diverged runs.
- **The build trap (FIXED):** `Ladruno_scripts/build.bat` took `set "MODE=%1"`, so extra target args were dropped — the MP `.pyd` was never rebuilt and a STALE pre-edit binary ran (old warnings + the un-scaled fallback → the run diverged). Fixed to route the whole non-`clean`/`rebuild` arg list (`set MODE=%*`). Symptom to recognize: the run prints OLD `opserr` warning text you already changed, or the built `.pyd` mtime is older than your edits — always confirm the binary is fresh and the build log shows `Step 4: Building targets: <your target>` + your changed TUs recompiling.
- **The validation trap (FIXED):** a stale-binary divergence made BOTH `np=1` and `np=2` overflow to the SAME garbage (`-2.1e+179`), so a tip-disp A/B comparator saw `diff=0` and FALSELY PASSED. A parallel correctness compare MUST reject non-finite / unphysically large output BEFORE comparing (`compare*.py` now guard `|disp|>1.0` and `isfinite`). Cross-check against an INDEPENDENT serial reference (serial `DiagonalSOE`+`consistentPCG`), not only `np=1`-vs-`np=2` of the same new code. Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]] V5.
