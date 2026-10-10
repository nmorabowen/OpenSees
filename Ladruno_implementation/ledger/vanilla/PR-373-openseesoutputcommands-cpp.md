---
wp: PR-373
title: "373 -- 1 vanilla row(s)"
pr: "#373"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [56]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-41 C2.0: `contactSurface` gains the `-slave-segments nps tag…` kind (faceted slave for the mortar `D`); the `contact` command gains a `-mortar` selector + `-epsN auto\|val / -augTol / -maxAug / -ngp` options → `addMortarContact` (stored separately from NTS, inert until C2.1). The kn pre-parse now treats ANY `-`-prefixed token (not just `-outward`) as "kn omitted" so `contact … -mortar …` with no numeric kn parses. `ladrunoContactInfo` appends a 4th element `numMortarContacts` (callers indexing [0..2] unaffected; `test_adr39_contact_p1` length-asserts updated to `[0,0,0,0]`). | [#373](https://github.com/nmorabowen/OpenSees/pull/373) |
