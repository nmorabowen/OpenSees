---
wp: LEGACY
title: "(pre-campaign; PR untracked) -- upstreamable-table row(s)"
files: ["`SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp`"]
table: "upstreamable"
legacy_seq: [342]
---
| `SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp` | **UNMARKED, previously UNLEDGERED** (same audit): file-scope global `SuperLUStat_t stat` collides with MSVC's POSIX `stat()` (via `<sys/stat.h>`/`<io.h>`) → renamed `superlu_stat`. Windows build portability fix; upstreamable. Add the `// Ladruno` marker. | (pre-campaign; PR untracked) |
