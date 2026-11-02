---
wp: WP-221
title: "Asymmetry guard query + -strictSymmetry: upstream files touched"
date: 2026-11-02
pr: "#1059"
files: ["SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.h", "SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.cpp", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.h", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp", "SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.cpp", "SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp", "SRC/system_of_eqn/eigenSOE/ArpackSOE.h", "SRC/system_of_eqn/eigenSOE/ArpackSOE.cpp", "SRC/system_of_eqn/linearSOE/LinearSOE.h", "SRC/interpreter/OpenSeesCommands.cpp", "SRC/interpreter/OpenSeesCommands.h", "SRC/tcl/commands.cpp", "SRC/tcl/commands.h"]
---
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.h` | WP-221: `asymWarned` becomes a `LadrunoSymmetryGuard symGuard` member; `getSymmetryGuard()`, `setStrictSymmetry()`, `solve()` override | #1059 |
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.cpp` | WP-221: addA guard records into `symGuard` and ignores the sampling budget under strict; `solve()` refuses a latched strict system | #1059 |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.h` | WP-221: same members as the PARDISO SOE | #1059 |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp` | WP-221: same guard/strict/solve changes as the PARDISO SOE | #1059 |
| `SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.cpp` | WP-221: `OPS_ProfileSPDLinDirectSolver` parses `-strictSymmetry` (it read no options before; unknown ones now warn) | #1059 |
| `SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp` | WP-221: stop with -1 when the LinearSOE's `solve()` fails; the return was stored and never read (ido -1) or not stored (ido 1), so eigen iterated on an unsolved X. Upstream defect | #1059 |
| `SRC/system_of_eqn/eigenSOE/ArpackSOE.h` | WP-221: `LadrunoSymmetryGuard` member, `solve()` override, `getSymmetryGuard()` | #1059 |
| `SRC/system_of_eqn/eigenSOE/ArpackSOE.cpp` | WP-221: every element K of every eigen is mirror-checked (not the `-shift*M` addM routes through addA); warn once; `solve()` refuses when the linear system is `-strictSymmetry` | #1059 |
| `SRC/system_of_eqn/linearSOE/LinearSOE.h` | WP-221: `virtual int getStrictSymmetry() const` (default 0), read by ArpackSOE | #1059 |
| `SRC/interpreter/OpenSeesCommands.cpp` | WP-221: `-strictSymmetry` in `OPS_PARDISOGenLinSolver`; `OPS_LadrunoSymmetryGuard` | #1059 |
| `SRC/interpreter/OpenSeesCommands.h` | WP-221: `OPS_LadrunoSymmetryGuard` declaration | #1059 |
| `SRC/tcl/commands.cpp` | WP-221: `-strictSymmetry` in classic `system Pardiso` / `system ProfileSPD`; `ladrunoSymmetryGuard` proc on `theSOE` | #1059 |
| `SRC/tcl/commands.h` | WP-221: `ladrunoSymmetryGuard` declaration | #1059 |
