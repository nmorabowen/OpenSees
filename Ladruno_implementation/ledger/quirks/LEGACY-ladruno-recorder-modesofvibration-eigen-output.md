---
wp: LEGACY
title: "Ladruno recorder modesOfVibration (eigen) output writes no data — modal DATA group/dataset collision"
legacy_seq: 44
---
### Ladruno recorder `modesOfVibration` (eigen) output writes no data — modal `DATA` group/dataset collision
- **Bites:** `recorder ladruno -N modesOfVibration` after `ops.eigen(n)` creates the
  `MODEL_STAGE[*]/RESULTS/ON_NODES/MODES_OF_VIBRATION(U)` group with a valid schema
  (ID, COMPONENTS) but `DATA` stays empty `(0, nNodes, nComp)` and no `MODE_k`
  datasets appear; HDF5-DIAG "can't synchronously write data / Write failed" errors
  fire. Happens for ANY model (reproduced with a bare elasticBeamColumn portal —
  NOT the known fiber-section `writeSections` noise), under both `ops.record()` and
  an `analyze()` step.
- **Why:** `LadrunoRecorder::recordModeChannel` (LadrunoRecorder.cpp ~1692) calls
  `ch.sink->begin()`, which (StreamingSink) creates `.../MODES_OF_VIBRATION(U)/DATA`
  as a **chunked dataset** (normal time-series layout), then tries to
  `h5::group::create` a **group** at `DATA/STEP_<step>` with `MODE_k` datasets
  under it (the MPCO modal layout). A group cannot be created beneath an existing
  dataset → the HDF5 calls fail, no modal data is written.
- **Status:** **FIXED 2026-05-31.** `recordModeChannel` no longer calls the StreamingSink
  `begin()`. It now owns the modal init, mirroring frozen `ResultRecorderModesOfVibration::record`:
  once per stage (idempotent via `H5Lexists`) it creates the result group
  (`h5::group::createResultGroup`) + `ID` dataset + `DATA` **group**, then per step writes
  `DATA/STEP_<step>/MODE_<k>` datasets with MODE/LAMBDA/OMEGA/FREQUENCY/PERIOD attrs. The
  validator `ladruno_format.py::_check_data_shape` was taught the modal layout (DATA = group
  of STEP_<step> groups of MODE_<k> datasets). **Modal eigenvectors now match frozen mpco to
  1e-12**; the EIGEN gate is promoted into the counted regression battery and `eigen_check.py`
  does a real modal value-parity diff vs the ref `.mpco`. Files: `SRC/recorder/LadrunoRecorder.cpp`,
  `ladruno_format.py`, `eigen_check.py`, `run_regression.bat`. Found + fixed 2026-05-31 by the
  new eigen coverage gate (the test scheme catching, then confirming the fix of, a real bug).
