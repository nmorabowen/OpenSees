---
wp: WP-248
title: "model testNDMaterial: registration in the model command and the build lists"
pr: "#1142"
date: 2026-11-23
files: ["SRC/modelbuilder/tcl/myCommands.cpp", "SRC/tcl/CMakeLists.txt", "SRC/Makefile", "SRC/modelbuilder/tcl/Makefile", "Win64/proj/tcl/tcl.vcxproj"]
---
| `SRC/modelbuilder/tcl/myCommands.cpp` | WP-248: `model testNDMaterial` (no count; it commits only on `ndTest commit`) builds a `TclNDMaterialTester`; the include. Additive; the stock testers are unchanged | #1142 |
| `SRC/tcl/CMakeLists.txt` | WP-248: `TclNDMaterialTester.cpp` | #1142 |
| `SRC/Makefile` | WP-248: `TclNDMaterialTester.o` | #1142 |
| `SRC/modelbuilder/tcl/Makefile` | WP-248: `TclNDMaterialTester.o` | #1142 |
| `Win64/proj/tcl/tcl.vcxproj` | WP-248: `TclNDMaterialTester.cpp/.h` | #1142 |
