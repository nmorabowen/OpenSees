# ADR 39 — P1b design: LadrunoContactDomain + surface + lifecycle hooks

> Pre-code design for P1b (on the proven P1a injection plumbing). Still ZERO
> force — P1b proves the **lifecycle + Domain-attachment + parser**, narrow phase
> is P2. Grounded in source read 2026-06-21. Parent: `39_..._adr.md`; P1a:
> `_adr39_p1_design.md`; loop: `_adr39_loop_state.md`.

## P1b goal

On the green P1a injection: add the **Domain-owned `LadrunoContactDomain`** (holds
surfaces + pair state, survives `AnalysisModel::clearAll`), the **commit/revert
lifecycle hooks** the design gate mandated, the **`contactSurface`/`contact`
parser**, and make the handler **rebind** its adapters to ContactDomain pairs by
key in `handle()`. Adapters stay zero-force (P2 adds the narrow phase).

## Components

1. **`LadrunoContactDomain`** — new `SRC/domain/contact/LadrunoContactDomain.{h,cpp}`.
   - Owned by `Domain` (raw ptr, nullable). Holds `LadrunoContactSurface`s, the
     contact definitions, and the `LadrunoContactPair` state map (keyed by
     (slaveNode, masterSeg)). Path-dependent state (gap0, friction) lives HERE.
   - `commit()` / `revertToLastCommit()` — commit/restore all pair state.
   - `buildAdapters()` — called from the handler's `handle()`; (P1b) returns the
     adapter set (still empty connectivity / zero force); (P2) per-segment.
   - Survives `domainChanged` because that calls `AnalysisModel::clearAll`, NOT
     `Domain::clearAll` (verified: `DirectIntegrationAnalysis.cpp:412`).
   - `sendSelf/recvSelf` (P1b stub → 0; carries the contact serialization later).

2. **`LadrunoContactSurface`** — `SRC/domain/contact/LadrunoContactSurface.{h,cpp}`.
   - From a **node-set** (slave) or an **element-face-set** (master segments).
   - P1b: stores node tags + segment connectivity; resolves Node*/coords at
     setDomain. No projection yet.

3. **`Domain` vanilla edits** (`// Ladruno` ADR-39, ledger):
   - `Domain.h`: `LadrunoContactDomain *theContactDomain;` member (init 0 in every
     ctor), + `setLadrunoContactDomain(LadrunoContactDomain*)` /
     `getLadrunoContactDomain()` accessors; delete in dtor.
   - `Domain::commit()` (`Domain.cpp:2185`, before `return 0;`): insert
     `if (theContactDomain != 0) theContactDomain->commit();`.
   - `Domain::revertToLastCommit()` (`Domain.cpp:2211`, after node/ele revert,
     before `applyLoad`): insert
     `if (theContactDomain != 0) theContactDomain->revertToLastCommit();`.
   - `Domain::clearAll()`: do NOT delete theContactDomain here on a normal
     domainChanged path — only on a true model wipe (mirror how other Domain-owned
     objects are handled; verify the wipe path). Decide: clear pairs but keep
     surfaces, or full reset. (Open Q-P1b-1.)

4. **Handler rebind** — `LadrunoContactHandler::handle()` asks
   `theDomain->getLadrunoContactDomain()->buildAdapters()` for the adapter set
   (instead of the P1a hardcoded single adapter) and injects each. If no
   ContactDomain (null) → no contact adapters (pure Plain behaviour, byte-identical).

5. **Parser** — `contactSurface`/`contact` commands (openseespy
   `OpenSeesCommands.cpp` + classic-Tcl `commands.cpp`, mirror the dual wiring the
   code gate required for `constraints`). They create the `LadrunoContactDomain`
   on the active `Domain` (lazily) + add surfaces / contact pairs. (≈the ~6
   wrapper edits the P1 gate flagged — they land HERE, not P1a.)

## Open questions for the (light) P1b review

> [!question] Q-P1b-1 (ContactDomain lifetime) — **RESOLVED**
> Delete/reset `theContactDomain` in `Domain::clearAll()` (`Domain.cpp:1041`, the
> `ops.wipe()` path) — `if (theContactDomain) { delete theContactDomain;
> theContactDomain = 0; }` next to the container `clearAll()` calls. This MIRRORS
> the ADR-30 `theEQs->clearAll()` leak-fix at `Domain.cpp:1054` (upstream-omitted
> cleanup that leaked into the next model). `domainChanged` correctly does NOT
> touch it (it runs `AnalysisModel::clearAll`, a different method). Also null-init
> in every `Domain` ctor + delete in `~Domain`.

> [!question] Q-P1b-2 (adapter↔pair rebind key stability)
> Adapters are destroyed each `handle()` and rebuilt; pair state on ContactDomain
> is keyed by (slaveNodeTag, masterSegId). Node tags are stable; the surface must
> assign master-segment ids stably (by face-set insertion order, fixed at
> setDomain). Confirm at P1b code time. (Lower risk — still zero-force in P1b.)

> [!question] Q-P1b-3 (cross-target lib wiring) — **RESOLVED**
> Put `LadrunoContactDomain` + `LadrunoContactSurface` in `SRC/domain/contact/`
> (target `OPS_Domain`, new `CMakeLists.txt` + `add_subdirectory(contact)` in
> `SRC/domain/CMakeLists.txt`). No new cycle: `Domain` is in `OPS_Domain`, and the
> handler (`OPS_Analysis`) ALREADY references `Domain` (LadrunoProjectionHandler
> includes Domain.h) ⇒ `OPS_Analysis→OPS_Domain` linkage exists, so the handler
> can reference `LadrunoContactDomain` the same way. The `LadrunoContactFE` adapter
> STAYS in `SRC/analysis/handler/` (OPS_Analysis) for P1b — it's an FE_Element,
> assembly-side; only the Domain-resident state + surfaces move to OPS_Domain.

## P1b acceptance gates

- Null ContactDomain (no `contactSurface`/`contact`): byte-identical to stock
  (the nullable-ptr guard) — under CDL + Newmark.
- Define a contact surface pair (zero-force adapters): model still runs;
  `commit`/`revertToLastCommit` hooks fire (instrument a counter) and are no-op
  on state at zero force; bitwise-identical to P1a.
- `ops.wipe()` then rebuild a fresh model: no leak, no stale ContactDomain (Q-P1b-1).
- Forced `domainChanged` (remove unrelated element): ContactDomain + pair state
  survive; adapters rebind by key (Q-P1b-2).
