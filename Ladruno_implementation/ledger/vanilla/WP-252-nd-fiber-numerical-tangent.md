---
wp: WP-252
title: "LadrunoNDFiber -tangent numerical: the classic Tcl builder passes the option through"
pr: "#1148"
date: 2026-11-23
files: ["SRC/material/section/TclModelBuilderSectionCommand.cpp"]
---
| `SRC/material/section/TclModelBuilderSectionCommand.cpp` | WP-252: the classic-Tcl `section LadrunoNDFiber` reads `-tangent`/`-tangentStep` through the shared parser and passes them to the Ladruno section's constructor (two statics beside `currentLadrunoSpec`). Stock sections unchanged | #1148 |
