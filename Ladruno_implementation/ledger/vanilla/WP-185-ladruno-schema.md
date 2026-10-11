---
wp: WP-185
pr: "#972"
title: "WP-185 vanilla rows"
date: 2026-10-11
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/tcl/commands.cpp`"]
table: "main"
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno WP-185`: `#include "LadrunoSchema.h"` and `OPS_LadrunoSchema()` (the DL `ladrunoSchema` entry point, next to `OPS_LadrunoBuild`). | [#972](https://github.com/nmorabowen/OpenSees/pull/972) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-185`: declares `OPS_LadrunoSchema()`. | [#972](https://github.com/nmorabowen/OpenSees/pull/972) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-185`: `#include "../interpreter/LadrunoSchema.h"`, the forward declaration and the classic-Tcl bridge `ladrunoSchema` (registered through the WP-168 hook by its table row). | [#972](https://github.com/nmorabowen/OpenSees/pull/972) |
