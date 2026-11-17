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

#### Follow-up: every g++ harness links statically, and the gate enforces it (2026-11-16)

All 23 harness tests now append `*cxx_runtime_flags()` from `tests/_cxx.py`
(`["-static"]` on Windows, `[]` elsewhere) to their g++ compile command, the
CPML one included. The grep in the paragraph above found 21. The new rule found
two more that look up g++ with single quotes and a `path=` argument:
`test_ladruno_membrane_kernel_cpp.py` and `test_mutation_membrane.py`.

The quirk gate rule `gxx-static` (`ci/check_quirk_patterns.py`, self-test
`ci/test_quirk_gxx_static.py`) fails any `tests/**/*.py` that calls
`which("g++" | "c++" | "clang++")` unless its code uses `cxx_runtime_flags` or
passes `"-static"`. It reads the AST, so comments, docstrings and a bare import
of the helper do not count. It has no waiver. A new harness copied from an old
one without the flag goes red in the static gate, not on someone's Windows box.

Static linking also fixes the sibling failure described in `tests/conftest.py`
(`_native_path_resync`): an exe that dies at load with 0xC0000135 because the
inherited PATH lost the compiler's `bin`. A static exe needs no MinGW DLL on PATH.
