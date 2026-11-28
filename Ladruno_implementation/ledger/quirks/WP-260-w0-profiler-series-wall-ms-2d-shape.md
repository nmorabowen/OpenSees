---
wp: WP-260
title: "Profiler HDF5 series/wall_ms dataset has 2D shape (N, 1), not 1D (N,)"
date: 2026-11-26
pr: "#1158"
---
### Profiler HDF5 `series/wall_ms` dataset has 2D shape `(N, 1)`, not 1D `(N,)`

- **Bites:** Python scripts reading `runs/<run_id>/series/wall_ms` from `ops.profiler("report", ...)` HDF5 output and iterating `[float(x) for x in ds[:]]` directly — raises `TypeError: only 0-dimensional arrays can be converted to Python scalars` on NumPy ≥ 1.25 / Python 3.12.
- **Why:** The C++ profiler writes per-step series datasets as 2D `(N, 1)` HDF5 datasets rather than 1D `(N,)` vectors.
- **Fix/status:** Always flatten with `.reshape(-1)` before scalar conversion (`[float(x) for x in grp["series"]["wall_ms"][:].reshape(-1)]`). Learned 2026-11-26 (#1158).
