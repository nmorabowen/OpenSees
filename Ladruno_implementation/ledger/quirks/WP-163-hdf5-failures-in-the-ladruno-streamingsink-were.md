---
wp: WP-163
title: "HDF5 failures in the Ladruno StreamingSink were silent: a result could vanish, and a duplicate request doubled the rows (WP-163)"
date: 2026-10-03
---
### HDF5 failures in the Ladruno `StreamingSink` were silent: a result could vanish, and a duplicate request doubled the rows (WP-163)
- **Bites:** no return code was checked in `createTimeSeries3d` / `begin()` / `appendSlab3d`. A failed `H5Dcreate`
  (disk full, a > 4 GiB chunk, a name clash) left a DATA-less group and every later `accept()` returned quietly —
  only HDF5-DIAG noise (or nothing). `-N displacement displacement` made two sinks on one group: the second's
  create failed, it still marked itself initialized, then appended into the first's group → DATA 2T rows,
  TIME/STEP every step twice.
- **Workaround/status:** fixed (WP-163 R5): checked creates/appends, one error by result name, channel stopped;
  a pre-existing result group is refused; a short buffer skips the whole step (DATA/TIME/STEP stay aligned);
  the id axis is tiled above a 1 GiB slab. HDF5 chunks must stay < 4 GiB and every chunk dim ≤ its max dim
  (a `[T×0×C]` dataset with chunk 1 cannot be created — an empty channel now writes nothing).
