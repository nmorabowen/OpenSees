---
wp: LEGACY
title: "The bare this->getMass(); side-effect idiom still lives in the cache-LESS LadrunoCST/LadrunoCSTPair — four pre-planted traps for the next mass-cache extension"
legacy_seq: 240
---
### The bare `this->getMass();` side-effect idiom still lives in the cache-LESS `LadrunoCST`/`LadrunoCSTPair` — four pre-planted traps for the next mass-cache extension
- **Bites:** `LadrunoCST.cpp` (`addInertiaLoadToUnbalance`, `getResistingForceIncInertia`) and `LadrunoCSTPair.cpp` (same two) call `this->getMass();` purely to refill class-static `K`, then read `K(i,i)` directly — byte-for-byte the idiom that silently corrupted Quad/LST inertia when the G2 cache landed there (a hit skips the formation; `K` holds the last *tangent*). Correct today only because these two elements have no cache; the ledger row for `LadrunoMassCache` lists them as "NOT applied (trivial formation)", which invites exactly the future extension that would fire all four sites at once.
- **Workaround/status:** in-source `// Ladruno (ADR-77 review wave): DO NOT add the G2 cache here without first rewriting this to consume getMass()'s return` breadcrumbs at all four sites. If the cache is ever extended: fix the callers FIRST (the Quad/LST pattern, PR #650), and remember the Zone-A dynamics/Rayleigh battery — not A/B equality — is the oracle that catches it. *2026-07-27 (ADR-77 review wave).*
