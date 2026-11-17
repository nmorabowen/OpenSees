---
wp: WP-199
date: 2026-11-16
title: "A g++-built harness exe segfaults on Windows: it loads Git for Windows' libstdc++-6.dll (WP-199)"
---
### A g++-built harness exe segfaults on Windows: it loads Git for Windows' libstdc++-6.dll (WP-199)

**Symptom.** `tests/test_ladruno_cpml_kernel_cpp.py` went red on a box with
WinLibs MinGW-W64 g++ 15.2 on PATH. `--self-check` printed most lines and then
the oracle run stopped, rc 139 / 0xC0000005, right after
`FD K_T_dyn BTET10-RH-plane-lump`. MSVC builds passed, `-O0` g++ builds passed,
and a 64 MB stack changed nothing. It looked like UB in the CPML kernel. It was not.

**Cause.** A default g++ build imports `libstdc++-6.dll`, and Windows resolves
that DLL through PATH. On this box Git for Windows' `C:\Program Files\Git\mingw64\bin`
comes before the WinLibs `bin`, so a binary compiled against GCC 15.2 headers runs
against Git's older libstdc++. gdb places the crash inside that DLL, in
`std::basic_ifstream` construction (`TokenStream::loadFile`, the first ifstream
the program opens). The `FD K_T_dyn ...` line is only where buffered stdout
stopped. The `-O0` build passed only because it inlined different code, not
because the kernel had a bug.

**Fix (this PR).** On Windows the test links the GCC runtime statically
(`-static`), so the exe carries its own libstdc++ and PATH order no longer matters.
Check: `objdump -p <exe> | grep "DLL Name"` lists no `libstdc++-6.dll`.

**Still exposed.** About 20 other `tests/test_*.py` harnesses compile with
`shutil.which("g++")` and no `-static`, so they hit the same trap on any box
with this PATH order. A crash inside `libstdc++-6.dll` from a `...\Git\mingw64\bin`
path is this quirk, not your kernel. (2026-11-16)
