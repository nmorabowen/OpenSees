---
wp: WP-232
title: "DirectIntegrationAnalysis::analyze sub-steps internally: one analyze(dt) call can commit several times"
date: 2026-11-13
---
### `DirectIntegrationAnalysis::analyze` sub-steps internally — one `analyze(dt)` call can commit several times

- **Bites:** a driver that wraps `analyze(1, dt)` and assumes one call is one committed step — for example one that does per-step work before and after the call (an MPM grid reset and particle update, WP-232). When a step fails and `numSubLevels != 0`, `analyze` retries through `analyzeSubLevel` with smaller steps and commits each, so the wrapper's per-step work runs once around several commits.
- **Why:** `DirectIntegrationAnalysis::analyze` calls `analyzeSubLevel(1, dT)` only when `numSubLevels != 0` (`DirectIntegrationAnalysis.cpp:176-177`); it recurses (`:280`, `:293`) and commits each successful sub-step.
- **Workaround / status (2026-11-13, WP-232 — disclosed):** a driver that must own the step boundary calls `analyzeStep(dt)` directly and does its own step-size cutting; per-step work that must not repeat on a retry belongs in the commit slot. This applies to transient analyses only: `StaticAnalysis` has no `analyzeStep` and no sub-levels (`StaticAnalysis.h:66`), so a static driver calls `analyze(1, false)` and flushes recorders itself.
