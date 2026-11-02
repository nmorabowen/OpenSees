---
wp: WP-221
title: "Asymmetry guard query + -strictSymmetry: upstream files touched"
date: 2026-11-02
pr: "TBD"
files: ["SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.h", "SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.cpp", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.h", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.cpp", "SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp", "SRC/interpreter/OpenSeesCommands.cpp", "SRC/interpreter/OpenSeesCommands.h", "SRC/tcl/commands.cpp", "SRC/tcl/commands.h"]
---
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.h` | WP-221: `asymWarned` becomes a `LadrunoSymmetryGuard symGuard` member; `getSymmetryGuard()`, `setStrictSymmetry()`, `solve()` override | TBD |
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.cpp` | WP-221: addA guard records into `symGuard` and ignores the sampling budget under strict; `solve()` refuses a latched strict system | TBD |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.h` | WP-221: same members as the PARDISO SOE | TBD |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp` | WP-221: same guard/strict/solve changes as the PARDISO SOE | TBD |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.cpp` | WP-221: `OPS_ProfileSPDLinDirectSolver` parses `-strictSymmetry` (it read no options before; unknown ones now warn) | TBD |
| `SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp` | WP-221: stop with -1 when the LinearSOE's `solve()` fails; the return was stored and never read (ido -1) or not stored (ido 1), so eigen iterated on an unsolved X. Upstream defect | TBD |
| `SRC/interpreter/OpenSeesCommands.cpp` | WP-221: `-strictSymmetry` in `OPS_PARDISOGenLinSolver`; `OPS_LadrunoSymmetryGuard` | TBD |
| `SRC/interpreter/OpenSeesCommands.h` | WP-221: `OPS_LadrunoSymmetryGuard` declaration | TBD |
| `SRC/tcl/commands.cpp` | WP-221: `-strictSymmetry` in classic `system Pardiso` / `system ProfileSPD`; `ladrunoSymmetryGuard` proc on `theSOE` | TBD |
| `SRC/tcl/commands.h` | WP-221: `ladrunoSymmetryGuard` declaration | TBD |
