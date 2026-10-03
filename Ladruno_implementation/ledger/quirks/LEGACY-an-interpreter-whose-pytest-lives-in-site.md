---
wp: LEGACY
title: "An interpreter whose pytest lives in site-packages loses it under python -S: the byte-identity child process fails before running any deck (harness)"
legacy_seq: 539
---
### An interpreter whose `pytest` lives in site-packages loses it under `python -S`: the byte-identity child process fails before running any deck (harness)
- **Bites:** `tests/test_ladruno_sanisand_sasme.py::test_existing_schemes_byte_identical` with an interpreter that keeps pytest only in its site-packages: Esmeralda's `~/ladruno_build_test/conan_venv/bin/python`, and the nmora desk's `pythoncore-3.12-64`.
  - The test spawns `sys.executable -S` (the Windows `-S` trap); `-S` drops site-packages, where that venv keeps pytest; `wp129_sanisand_byteid` imports `test_ladruno_sanisand`, which imports pytest → `ModuleNotFoundError`.
- **Rule:** Read that failure as an environment artifact unless its message is a row mismatch. Put the interpreter's site-packages on `PYTHONPATH` (the child builds `sys.path` from it) and the test runs for real.
- **Workaround/status:** documented (WP-152 review, 2026-09-29): with site-packages on `PYTHONPATH` it PASSES on the nmora desk at 7f1562c81 (111 s). The fix belongs to the harness, not to the material.
