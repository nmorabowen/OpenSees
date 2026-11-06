---
wp: ADR-66
title: "The ladruno_opensees.pth boot module pins ONE worktree's pyd — a fresh build in ANOTHER worktree is silently ignored (ADR-66 P5.1)"
legacy_seq: 146
---
## The `ladruno_opensees.pth` boot module pins ONE worktree's pyd — a fresh build in ANOTHER worktree is silently ignored (ADR-66 P5.1)

`Ladruno_scripts/wire_venv_pth.py` writes `_ladruno_opensees_boot.py` into the py-3.12
site-packages with the generating checkout's `dist\bin` HARD-CODED, `sys.path.insert(0)`-ed at
interpreter startup, and — the sharp edge — an EAGER `import opensees` (for the
`openseespy` aliasing), so `opensees` is already in `sys.modules` before any test bootstrap or
`PYTHONPATH` entry can win. **Symptom:** you build a NEW element in worktree B, the build exits 0,
the pyd timestamp is fresh — and pytest says `element type X is unknown`, because the import came
from worktree A (check `opensees.__file__` FIRST when a freshly-built symbol is "unknown").
**Bypass without touching the other session's wiring:** set `PMI_RANK=1` in the child env (the boot's
MPI guard skips the eager import + aliasing) and `sys.path.insert(0, <your dist\bin>)` +
`os.add_dll_directory` + PATH-prepend in a small driver BEFORE importing pytest
(the P5.1 `run_gates.py` pattern). Re-running `wire_venv_pth.py` re-pins instead, but stomps the
sibling session. **2026-08-10: superseded for the common case** — `set LADRUNO_OPENSEES_BIN=<your
dist\bin>` before importing wins WITHOUT stomping the sibling session or needing the PMI_RANK trick
(the boot module itself now checks the env var first); see the "An INSTALLED Ladruno hijacks" entry
below for the fix detail.
