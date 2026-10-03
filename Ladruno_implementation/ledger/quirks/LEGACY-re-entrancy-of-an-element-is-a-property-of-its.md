---
wp: LEGACY
title: "Re-entrancy of an element is a property of its whole call graph, and of its CONFIGURATION — not of its class"
legacy_seq: 452
---
### Re-entrancy of an element is a property of its whole call graph, and of its CONFIGURATION — not of its class

Two traps that the obvious "audit the element class" framing misses:

1. **The material decides.** `ManzariDafalias` is re-entrant under `IntScheme 1`
   (ModifiedEuler: no static work arrays, no `Matrix::Invert`) and is NOT under
   `IntScheme 2` or `4`. Same class, same object, opposite answers — which is why
   WP-107 made the allowlist a **runtime virtual** (`ladrunoThreadSafeUpdate()`)
   rather than a class-tag table. The same holds one level up: `LadrunoQuad` is
   re-entrant for `std`/`bbar`/`ssp` under `-geom linear` and is not for `eas` or
   `-geom finite`.
2. **A diagnostic ledger can block threading even when the physics is clean.**
   `LadrunoSANISAND` under `-implex` keeps a process-wide accumulator
   (`LadrunoImplexGlobals`) whose `sumError`/`maxError` are **floating point**.
   Adding an atomic does not fix it: a threaded sum changes the reported average's
   last bits with the thread count, which is exactly the determinism the threaded
   loop exists to preserve. The right answer was to refuse, not to "make it
   thread-safe".
