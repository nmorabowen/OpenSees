---
wp: PR-706
title: "706 -- upstreamable-table row(s)"
pr: "#706"
files: ["`.gitattributes`"]
table: "upstreamable"
legacy_seq: [406]
---
| `.gitattributes` | `# --- Ladruno: re-enable EOL normalization on the shared registration surface ---`: an override block appended at the END of the file (last matching line wins) turning `text` back ON for the eight files every fork feature touches to register a class — `SRC/interpreter/OpenSeesCommands.{cpp,h}`, `SRC/interpreter/{Python,Tcl}Wrapper.cpp`, `SRC/interpreter/OpenSees{Element,Pattern}Commands.cpp`, `SRC/tcl/commands.cpp`, `SRC/classTags.h`. WHY: the inherited upstream `.gitattributes` is 5,345 lines of which **5,323 are per-file `-text`** (a CVS/SVN→git conversion artifact); `-text` means git stores bytes VERBATIM and normalizes nothing, so those files have no canonical line ending and any contributor whose editor differs rewrites the whole file — how #700 recorded 10,543 changed lines for 63 real ones. All eight targets were already `i/lf` in the blob, so this is a **ZERO-DIFF** change (`git add --renormalize` touches no source byte); verified A/B by CRLF-ifying a covered file (staged diff **empty**) vs the still-`-text` `OpenSeesCommandsTcl.cpp` (**188/188** whole-file churn). Deliberately NOT a repo-wide renormalize: ~447 blobs are still CRLF and rewriting them would be permanent whole-file divergence from upstream. | [#706](https://github.com/nmorabowen/OpenSees/pull/706) |
