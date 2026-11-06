---
wp: PR-757
title: "// Ladruno ADR-85 T3 (2D mortar -thickness): FOUR edit blocks in ladrunoContactImpl -- (1) the option-table comment, (2) the -thickness <h> parser branch (one…"
pr: "#757"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "upstreamable"
legacy_seq: [476]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno ADR-85 T3` (2D mortar `-thickness`): FOUR edit blocks in `ladrunoContactImpl` -- (1) the option-table comment, (2) the `-thickness <h>` parser branch (one double, `h > 0` validated at parse), (3) the post-loop `-mortar`-only validation (mirrors the `-tauMax` rule), (4) the `addMortarContact` call-site 34th argument. Dimension routing stays at handle() (a parse-time surface peek cannot know the pair dimension); a 3D mortar pair with `h != 1` draws the handler's named FATAL. 3D decks are byte-identical: the flag defaults to 1.0 and every consumer multiplies by it only inside the 2D injection branch. Measured at ship: `contact_dump` bit-identical x2 (`B0F8F770...81E4`), 3D battery 142 passed N unchanged. PR: [#757](https://github.com/nmorabowen/OpenSees/pull/757) -- ADR-85 T3. |
