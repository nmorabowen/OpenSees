---
wp: LEGACY
title: "LadrunoRecorder shell per-layer stress: the verb is a DOTTED token, and material.fiber.* needed a fix"
legacy_seq: 47
---
### LadrunoRecorder shell per-layer stress: the verb is a DOTTED token, and `material.fiber.*` needed a fix
- **Bites (1) — syntax:** the `-E` element verb is a single **dot-joined** token that the
  recorder splits on `'.'` (`recorder('ladruno', f, '-E', 'section.fiber.stress')`), NOT
  space-separated args. Passing `'-E','section','fiber','stress'` stores only `["section"]`
  (the rest are parsed as later options), which then **segfaults** (see bite 3).
- **Bites (2) — the real bug:** for a layered shell, the *obvious* per-layer verb
  `material.fiber.stress` silently emitted **no element bucket**. Root cause: the request
  builder only set `do_all_fibers` for `fiber` under `do_all_sections`, and the
  `do_all_materials` path had no fiber-index expansion — so the fiber id was never
  substituted and `setResponse` returned null for every element. `section.fiber.stress`
  worked because the recorder swaps section->material for shells and runs the section
  fiber-expansion. The section-level read itself was always fine
  (`LayeredShellFiberSection`/`MembranePlateFiberSection` answer `fiber <i> <resp>` and tag
  `FiberOutput`; `eleResponse(tag,'material','1','fiber','1','stress')` works directly).
- **Fix:** extend the `fiber` trigger to `do_all_sections || do_all_materials` and give the
  `do_all_materials` branch the same per-(gp,layer) expansion as `do_all_sections` (driven
  by the shared `elem_ngauss_nfiber_info` discovery table). Now `material.fiber.stress`
  emits the per-layer bucket **byte-identical** to `section.fiber.stress` (verified: 4 GP ×
  3 layers × 5 comp = 60 cols, maxdiff 0.0). Regression gate `SHELL LAYER STRESS`
  (`shell_layer_model.py`/`shell_layer_check.py`) in `run_regression.bat`.
- **Bites (3) — latent crash (FIXED):** a bare `-E section` / `-E material` (the keyword
  with NO sub-verb) segfaulted — `request_mod` was `["section",""]` (argc=2), so the element
  was queried as `["section"/"material", <id>]` and forwarded a **zero-length** arg list
  (`&argv[2]`, argc-2==0) to its section's `setResponse`, which derefs `argv[0]`. Affected
  any shell or beam. **Fix:** guard both non-fiber `setResponse` call sites in
  `initElementSources` — only call when `argc > (int)<keyword>_id_placeholder_index + 1`
  (i.e. at least one sub-verb token follows the id). A bare verb now emits no bucket instead
  of crashing. Regression gate `SHELL BARE VERB` (`shell_bare_verb_model.py` /
  `shell_bare_verb_check.py`). Learned + fixed 2026-06-03.
- **Known cosmetic limitation:** the shell per-layer COLUMN_MAP records `fiber_id=-1`,
  `section_tag=-1`, and `UnknownStress` component names (layer identity is flattened into a
  running `gauss_id`). This is IDENTICAL on the already-working `section.fiber.stress` path
  (not introduced here) — a metadata refinement for later, the stress *values* are correct.
  (Adversarial-review note: a reviewer claimed this *collapses* layers into one
  `normalize_element` key and loses data — REFUTED: the 3 layers map to distinct `gauss_id`
  0..11, so the parity dict has 4·12·5=240 distinct entries, no loss.)
