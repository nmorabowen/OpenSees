---
wp: LEGACY
title: "wipe / Domain::clearAll() did NOT clear EQ_Constraints (upstream bug) — leaks across models"
legacy_seq: 105
---
### `wipe` / `Domain::clearAll()` did NOT clear EQ_Constraints (upstream bug) — leaks across models
- **Bites:** `EQ_Constraint` (the `equationConstraint` command, a later upstream addition) was never wired into `Domain::clearAll()` — it clears `theSPs/thePCs/theMPs` but omitted `theEQs`. So `ops.wipe()` (and any `clearAll`) LEAVES equation constraints in the domain; the next model silently inherits them. This is invisible with the stock handlers (none of them iterate `getEQs()` in the common path), so it sat latent. It surfaces the moment a handler reads `getEQs()` — `LadrunoProjectionHandler` (ADR-30 P3) does, and a stale EQ from a prior model then mis-assembles a constraint group (wrong groups / partition-guard refusal / wrong projection) in the NEXT analysis. Found 2026-06-20 by the ADR-30 P3 full-suite regression: the EQ test poisoned every subsequent projection test ONLY in combined runs (passed in isolation) — the classic test-ordering signature of leaked global domain state.
- **Fix:** add `theEQs->clearAll();` to `Domain::clearAll()` (one line, mirrors `theMPs->clearAll()`). Upstreamable. Ledgered in LEDGER_vanilla. **General lesson:** a "passes alone, fails in a combined pytest run" ordering failure = leaked global OpenSees state; suspect a container `wipe`/`clearAll` doesn't clear (here EQ_Constraints). Run new constraint/handler features in a COMBINED suite, not just in isolation, to catch it.
