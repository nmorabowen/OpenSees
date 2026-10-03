---
wp: LEGACY
title: "LadrunoStaggeredAnalyze overrides -subcycle while driving — the driver syncs the fluid every step by construction"
legacy_seq: 186
---
### `LadrunoStaggeredAnalyze` overrides `-subcycle` while driving — the driver syncs the fluid every step by construction
- **Bites:** an overlay built with `-subcycle N>1` (or `-subcycle auto`) accumulates the Δu window across N commits under plain `analyze`, but under the ADR-73 P2 driver `LadrunoStaggeredAnalyze` the fluid is advanced (and committed) at EVERY driver step — the configured N is ignored. If you expect the same subcycled fluid cadence you configured, you don't get it while driving.
- **Why:** the driver's whole point is the iterated per-step fixed-point solve; a multi-commit accumulation window is incompatible with re-solving the solid against the current-step fluid iterate. The latched `onDomainCommit` (SM_MARCH branch) does `commitFluid` + counter reset only — no window accumulation, no extra `advanceTrial`.
- **Workaround/status (2026-07-17, ADR-73 P2):** intended, not a bug. A one-time advisory prints ("LadrunoStaggeredAnalyze overrides -subcycle while driving") when the driven overlay has `subcycleN>1` or `-subcycle auto`. A pending window at driver entry is caught up with an early fs1 sync (`catchUpPendingWindow`) so the first driver advance doesn't pair a multi-commit Δu with the single driver dt. Counters are zeroed at each latched commit, so a post-driver plain `analyze` restarts its window cleanly.
