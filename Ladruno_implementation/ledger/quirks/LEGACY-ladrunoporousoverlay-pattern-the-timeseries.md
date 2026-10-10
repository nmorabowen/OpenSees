---
wp: LEGACY
title: "LadrunoPorousOverlay pattern: the TimeSeries/load factor is IGNORED by design — the overlay owns its force amplitudes"
legacy_seq: 182
---
### LadrunoPorousOverlay pattern: the TimeSeries/load factor is IGNORED by design — the overlay owns its force amplitudes
- **Bites:** attaching a `timeSeries` to `pattern LadrunoPorousOverlay ...` (or expecting `loadConst`-style factor scaling) does nothing: the injected nodal forces are always `+Q·p_committed` at full amplitude — the pore-pressure field, not a factored load, is the amplitude. Silent expectation mismatch if you try to "ramp" the overlay.
- **Why:** the overlay is a domain ENGINE riding the LoadPattern plumbing (H5DRM precedent); its `applyLoad(time)` ignores `time` and any series. A one-time informational notice prints if a series was assigned. Note the python/Tcl surface cannot even attach a series to it structurally (verified 1.E-ii, 2026-07-14).
- **Workaround/status:** ramp the SOLID loads (they live in ordinary patterns); stage the fluid via `-pInit` / `-staticMode`. ADR-73 §4.1; P1 battery gates the bit-exactness of "factor changes nothing".
