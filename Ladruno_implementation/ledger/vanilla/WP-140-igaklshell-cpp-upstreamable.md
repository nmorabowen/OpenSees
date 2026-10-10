---
wp: WP-140
title: "880 -- upstreamable-table row(s)"
pr: "#880"
files: ["`SRC/element/IGA/IGAKLShell.cpp`"]
table: "upstreamable"
legacy_seq: [674]
---
| `SRC/element/IGA/IGAKLShell.cpp` | `// ladruno-lint: sequence-ok` (WP-140) — COMMENT ONLY, no code change: the lint L7 waiver on `res = this->getResistingForce() + this->getMass() * NodalAccelerations;` in `getResistingForceIncInertia`. Audited harmless: `getMass()` writes only `*mass`, never the `*resid` `getResistingForce()` returns, so either call order gives the same value. Owner chose the waiver over a rewrite or an IGA exemption (2026-09-27). Stale-waiver detection only scans fork-stamped files, so this vanilla waiver is not re-checked automatically. | [#880](https://github.com/nmorabowen/OpenSees/pull/880) |
