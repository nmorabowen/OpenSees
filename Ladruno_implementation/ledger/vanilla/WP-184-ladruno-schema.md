---
wp: WP-184
title: "WP-184 vanilla rows"
date: 2026-10-11
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/tcl/commands.cpp`"]
table: "main"
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno WP-184`: `#include "LadrunoSchema.h"` and `OPS_LadrunoSchema()` (the DL `ladrunoSchema` entry point, next to `OPS_LadrunoBuild`). | — |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-184`: declares `OPS_LadrunoSchema()`. | — |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-184`: `#include "../interpreter/LadrunoSchema.h"`, the forward declaration and the classic-Tcl bridge `ladrunoSchema` (registered through the WP-168 hook by its table row). | — |
