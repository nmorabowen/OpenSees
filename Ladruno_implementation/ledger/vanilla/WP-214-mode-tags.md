---
wp: WP-214
title: "modalProperties report: append the Ladruno mode-kind tags"
date: 2026-10-26
pr: "#1026"
table: "main"
files: ["`SRC/domain/domain/DomainModalProperties.cpp`"]
---
| `SRC/domain/domain/DomainModalProperties.cpp` | `// Ladruno WP-214`: `print_internal` and `printDict` ask `LadrunoModeKinds::report(OPS_GetDomain(), eigenvalues())` for the WP-213 mode kinds and append its strings: a `  # R` / `  # Z` trailing comment on each mode row of sections 2 and 6–10, `# NOTE:` lines after the section-2 table, and a `modeKind` key in `-return`. The hook returns empty strings (and no key) unless a valid registry record matches the report's eigenvalues bit for bit, so without `residualModes` / Ritz modes the report is byte-identical. Decoration from outside was not possible: the rows are written straight to `opserr` / the file stream. | #1026 |
