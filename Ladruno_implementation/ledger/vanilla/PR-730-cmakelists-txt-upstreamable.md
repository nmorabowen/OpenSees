---
wp: PR-730
title: "730 -- upstreamable-table row(s)"
pr: "#730"
files: ["`CMakeLists.txt`"]
table: "upstreamable"
legacy_seq: [425]
---
| `CMakeLists.txt` | `// Ladruno` (ADR-78 P1, backfilled — #730/#731 shipped without ledger rows): introduce `OPS_CONTACT_PER_TARGET_SOURCES` holding `LadrunoContactAbort.cpp`, and add it to all five targets, so the contact fatal exit is compiled with each target's parallel defines instead of once with none. | [#730](https://github.com/nmorabowen/OpenSees/pull/730) |
