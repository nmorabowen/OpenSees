"""Link flags for the standalone g++ harnesses (WP-199 follow-up).

About twenty `tests/test_*.py` files compile a small C++ driver with
`shutil.which("g++")` (or c++ / clang++) and run it. On Windows a default
MinGW build imports `libstdc++-6.dll`, `libgcc_s_seh-1.dll` and
`libwinpthread-1.dll`, and Windows resolves them through PATH. When Git for
Windows' `C:\\Program Files\\Git\\mingw64\\bin` precedes the compiler's own
`bin`, an exe built against GCC 15 headers runs against Git's older libstdc++
and segfaults (0xC0000005), usually while constructing a `std::ifstream`.
Linking the GCC runtime statically makes PATH order irrelevant.

Every g++ harness compile appends ``*cxx_runtime_flags()`` to its command; the
quirk gate (`check_quirk_patterns.py --only gxx-static`) fails a harness that
does not. LEDGER quirk: "loads Git for Windows' libstdc++-6.dll".

Top-level module on purpose (like `_engine.py`): it imports only `os`, so the
harnesses stay runnable standalone.
"""
import os


def cxx_runtime_flags():
    """Extra link flags for a g++/clang++ harness build: ``["-static"]`` on Windows, else none.

    Append them after the sources and ``-o``: they are link flags.
    """
    return ["-static"] if os.name == "nt" else []
