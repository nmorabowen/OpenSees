---
wp: LEGACY
title: "Overlay -layer moduli overrides use a > 0 sentinel — layerNu 0 is unreachable via the parameter route (rejected loudly)"
legacy_seq: 188
---
### Overlay `-layer` moduli overrides use a `> 0` sentinel — `layerNu 0` is unreachable via the parameter route (rejected loudly)
- **Bites:** trying to set a per-layer Poisson ratio of exactly 0.0 through `parameter $p loadPattern $tag layerNu $i` + `updateParameter` fails with a warning, even though nu = 0 is a physically legal value (and the global `nu` accepts it).
- **Why:** the `Layer` struct encodes "inherit the overlay-global value" as `nu <= 0` (P1 sentinel, serialized that way); a stored layer nu of 0.0 would be silently re-interpreted as "unset" by `resolveCellModuli`, turning the update into a no-op. The P2 panel (robustness-7) flagged the silent path; the fix rejects it loudly instead.
- **Workaround/status (2026-07-17, ADR-73 P2):** set the GLOBAL `nu` to 0 (parameter id `nu`) and leave the layer inheriting, or use a tiny positive value. Changing the sentinel to an explicit per-field override flag would touch the serialized layer payload — deferred until a real user needs layered nu = 0.
