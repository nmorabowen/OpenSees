---
wp: WP-203
title: "eigen -feast: trim node eigenvectors with the eigenvalues; -traceStop"
pr: "#1010"
date: 2026-10-22
files: ["`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/tcl/commands.cpp`"]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-203`: the FEAST reconcile trimmed only the domain eigenvalue list to the found count, leaving every node's eigenvector matrix at the analysis-loop cap (>=128 columns), so `modalProperties` aborted with "inconsistent eigenvector matrix and ID" (DomainModalProperties.cpp:589) under every constraint handler. Now calls `FeastEigenSOE::trimDomainModes` (eigenvalues AND node eigenvectors). Also parses the new `eigen -feast ... -traceStop` (FEAST now stops on the residual by default; see the WP-203 quirk). | [#1010](https://github.com/nmorabowen/OpenSees/pull/1010) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-203`: same reconcile fix and `-traceStop` for classic Tcl `eigen -feast`. | [#1010](https://github.com/nmorabowen/OpenSees/pull/1010) |
