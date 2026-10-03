---
wp: PR-712
title: "712 -- upstreamable-table row(s)"
pr: "#712"
files: ["`SRC/interpreter/PythonModule.cpp`"]
table: "upstreamable"
legacy_seq: [410]
---
| `SRC/interpreter/PythonModule.cpp` | `// Ladruno`: **harden module init/shutdown against same-DLL re-import** (the LEDGER_quirks "re-import ... crashes Python AT SHUTDOWN" entry). `moduledef.m_size > 0` makes CPython RE-RUN `PyInit_opensees` when the module is re-imported after a `sys.modules.pop` of the same pyd; the init unconditionally called `Py_AtExit(cleanupFunc)`, so cleanup registered twice, and `cleanupFunc` deref'd the static `PythonModule* module` BEFORE its null check and never nulled it after `delete` — the second atexit call at `Py_FinalizeEx` wiped a dangling pointer (0xC0000005/0xC0000409, heap-state dependent). Two-part fix: `Py_AtExit` behind a `static bool` (register exactly once per TU), and a null-safe idempotent `cleanupFunc` (`if (module) { wipe; delete; module = 0; }`). One TU per Python target (`PythonMPIModule.cpp` re-includes this file), so `openseesmp` inherits the fix. Gated by the 3-line repro: `python3.12 -c "import sys, opensees; sys.modules.pop('opensees'); import opensees"` — exit −1073741819 pre-fix, 0 post-fix. | [#712](https://github.com/nmorabowen/OpenSees/pull/712) |
