---
wp: LEGACY
title: "A missing OPTIONAL dependency does not skip two tests — it aborts the ENTIRE suite at collection"
legacy_seq: 293
---
### A missing OPTIONAL dependency does not skip two tests — it aborts the ENTIRE suite at collection
- **Bites:** `tests/test_ladruno_overlay_{driver,physics}.py` import the ADR-71 frozen toy for its SOLVERS; that toy does `import matplotlib` at module scope. On a box without matplotlib the import raises during COLLECTION, and pytest treats collection errors as fatal: `!!!! Interrupted: 2 errors during collection !!!!`, `9 skipped, 2 errors`. All ~2000 other tests never ran. `pytest tests/` was simply unusable, and had been for as long as the box lacked matplotlib.
- **Why it evades the usual guards:** a module-level `pytest.skip(allow_module_level=True)` (which this repo uses correctly for the build gate and the gmsh gate) contains the damage; a bare `import` does not. The failure also names only the two modules, so it reads as "two broken files" rather than "the suite cannot run".
- **What catches it:** guard transitive optional deps with `pytest.importorskip("<dep>")` BEFORE the import that needs them, so the blast radius is the module rather than the run. And periodically run the whole suite on a machine that does NOT have the research-only extras — the gap only shows up there.
- **NB the fix belongs in the TEST, not the frozen module:** `meshless_p_toy.py` is ADR-cited and marked DO NOT MODIFY, so the import is not lazified at source. Install matplotlib to actually run those tests; the guard only stops them taking everything else down with them.
- **Related:** the whole `zone_b` tier (129 tests) self-skips without gmsh via `tests/conftest.py`. That one is BY DESIGN and correct — the contrast is the point: a declared, named skip is fine; an unhandled import is not.
