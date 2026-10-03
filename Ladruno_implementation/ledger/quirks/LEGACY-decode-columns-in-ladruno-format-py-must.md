---
wp: LEGACY
title: "decode_columns in ladruno_format.py must flatten COLUMN_MAP arrays — the recorder writes them 2-D [k×1]"
legacy_seq: 39
---
### `decode_columns` in `ladruno_format.py` must flatten COLUMN_MAP arrays — the recorder writes them 2-D `[k×1]`
- **Bites:** `int(mult[i])` (and the other per-block scalars) in `decode_columns`
  throws `TypeError: only 0-dimensional arrays can be converted to Python scalars`
  under numpy 2.x when reading a **real recorder** `.ladruno`.
- **Why:** the C++ writer stores each per-block COLUMN_MAP array via
  `createAndWrite(vec, k, 1)` → **2-D `[k×1]`**, so `arr[i]` is a `(1,)`-array, not
  a scalar. `make_synthetic.py` writes them 1-D `[k]`, which masked it; and the
  element-**parity** gate keys results by flat column index and never calls
  `decode_columns`, so no test exercised it on real 2-D output until PR #45's
  element-envelope checker.
- **Status:** fixed (PR #45) — `decode_columns` `.reshape(-1)`s GAUSS_ID/SECTION_TAG/
  FIBER_ID/NUM_COMP/MULTIPLICITY (LEVELS is consumed via `np.atleast_1d`, left as-is).
  Learned 2026-05-31.
