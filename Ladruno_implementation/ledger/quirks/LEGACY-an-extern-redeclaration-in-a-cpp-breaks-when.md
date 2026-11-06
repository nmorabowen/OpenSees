---
wp: LEGACY
title: "An extern redeclaration in a .cpp breaks when the definition becomes thread_local"
legacy_seq: 454
---
### An `extern` redeclaration in a `.cpp` breaks when the definition becomes `thread_local`

Making `ops_TheActiveElement` `thread_local` compiled everywhere that included
`OPS_Globals.h` / `G3Globals.h`, and failed with MSVC **C2370 "redefinition;
different storage class"** in the two files that had declared it themselves:
`LadrunoDispBeamColumn2d.cpp:60` and `3d.cpp:56`. If you change the storage class
of a global in this tree, grep for local `extern` copies of it — the headers are
not the whole story.
