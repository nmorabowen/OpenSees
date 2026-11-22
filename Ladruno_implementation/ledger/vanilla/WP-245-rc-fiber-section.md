---
wp: WP-245
title: "LadrunoNDFiberSection3d and LadrunoUniaxialFiber: registration and a protected NDFiberSection3d constructor"
pr: "TBD"
date: 2026-11-22
files: ["SRC/material/section/NDFiberSection3d.h", "SRC/material/section/NDFiberSection3d.cpp", "SRC/classTags.h", "SRC/material/section/CMakeLists.txt", "SRC/material/section/Makefile", "SRC/Makefile", "SRC/material/nD/CMakeLists.txt", "SRC/interpreter/OpenSeesSectionCommands.cpp", "SRC/interpreter/OpenSeesNDMaterialCommands.cpp", "SRC/material/section/TclModelBuilderSectionCommand.cpp", "SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp", "SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp"]
---
| `SRC/material/section/NDFiberSection3d.h` | WP-245: a protected `NDFiberSection3d(tag, classTag, numFibers, alpha, compCentroid)` so the subclass `LadrunoNDFiberSection3d` registers its own classTag. Additive; stock constructors and behaviour unchanged | TBD |
| `SRC/material/section/NDFiberSection3d.cpp` | WP-245: the body of that constructor (a copy of the empty-section constructor with the classTag passed through) | TBD |
| `SRC/classTags.h` | WP-245: `SEC_TAG_LadrunoNDFiberSection3d 33001`, `ND_TAG_LadrunoUniaxialFiber 33030` | TBD |
| `SRC/material/section/CMakeLists.txt` | WP-245: `LadrunoNDFiberSection3d.cpp/.h`, `LadrunoFiberShearKernel.h` | TBD |
| `SRC/material/section/Makefile` | WP-245: `LadrunoNDFiberSection3d.o` | TBD |
| `SRC/Makefile` | WP-245: `LadrunoNDFiberSection3d.o` | TBD |
| `SRC/material/nD/CMakeLists.txt` | WP-245: `LadrunoUniaxialFiber.cpp/.h` | TBD |
| `SRC/interpreter/OpenSeesSectionCommands.cpp` | WP-245: `section LadrunoNDFiber` (3D only); sets the active `NDFiberSection3d` pointer so the stock `fiber`/`patch`/`layer` commands fill it | TBD |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | WP-245: `nDMaterial LadrunoUniaxialFiber` | TBD |
| `SRC/material/section/TclModelBuilderSectionCommand.cpp` | WP-245: classic-Tcl `section LadrunoNDFiber $tag opts { fibres }` (`TclCommand_addLadrunoNDFiberSection`: the shared fail-closed option parser, the stock fibre block, no `-GJ`), a `currentSectionIsLadruno` flag reset by the stock fibre-section entry, and the 3D build branch that constructs the Ladruno section | TBD |
| `SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp` | WP-245: classic-Tcl `nDMaterial LadrunoUniaxialFiber` | TBD |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | WP-245: broker cases for SEC 33001 and nD 33030 | TBD |
