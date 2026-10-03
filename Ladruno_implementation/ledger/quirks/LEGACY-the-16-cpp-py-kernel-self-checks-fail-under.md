---
wp: LEGACY
title: "The 16 *_cpp.py kernel self-checks FAIL under plain pytest and PASS under pytest -s — it's stdout capture breaking subprocess, not a regression"
legacy_seq: 268
---
### The 16 `*_cpp.py` kernel self-checks FAIL under plain `pytest` and PASS under `pytest -s` — it's stdout capture breaking `subprocess`, not a regression
- **Bites:** a full Zone-A sweep reports `16 failed, 1742 passed`, and every failure is a C++ kernel self-check (`test_hypo_kernel_cpp`, `test_ladrunoJ2_*_cpp`, `test_logstrain*_cpp`, `test_ladrunoConcrete3D_material`, `test_ladrunoCMS_core_cpp`, `test_ladruno_up_kernel_cpp`, `test_ladrunoRCConcrete_{reg,tensstiff}_cpp`, …). It reads as if a source change just broke every material kernel at once. It did not.
- **Why:** those tests shell out to `g++` via `subprocess.run(..., capture_output=True)`. Under pytest's default capture the parent's stdio handles are not real console handles, and Windows fails the spawn with **`OSError: [WinError 50] The request is not supported`** at `subprocess.py:1416` — *before* any C++ is compiled or any OpenSees code runs. `g++` itself is fine (`/c/msys64/mingw64/bin/g++`, MSYS2 rev8 15.2.0) and the same `subprocess.run` succeeds from a plain `python3.12 -c`.
- **Tell:** the failure is an `OSError` at spawn time, never an assertion or a numeric mismatch. A real kernel regression fails on *values*.
- **Workaround/status (2026-08-04):** run those tests with **`pytest -s`** (capture disabled) — `16 passed` on the identical binary that "failed" 16. So the honest Zone-A total is **1758 passed, 0 real failures**. Verify a suspected kernel regression with `-s` **before** believing it. Two further modules (`test_ladruno_overlay_{driver,physics}.py`) cannot even be *collected* here because `matplotlib` is absent — also unrelated; `--ignore` them or install matplotlib.
