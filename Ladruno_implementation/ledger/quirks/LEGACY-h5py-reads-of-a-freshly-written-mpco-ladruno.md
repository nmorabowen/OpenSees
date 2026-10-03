---
wp: LEGACY
title: "h5py reads of a freshly-written .mpco/.ladruno HANG on HDF5 file locking"
legacy_seq: 40
---
### h5py reads of a freshly-written `.mpco`/`.ladruno` HANG on HDF5 file locking
- **Bites:** a venv-python checker calling `h5py.File(path, "r")` on a file the
  build-python just wrote (in the same gate run) **hangs indefinitely** — no error,
  just blocks (e.g. `parity_check.py` stuck forever).
- **Why:** HDF5's default file locking; the writer's lock/superblock state isn't
  cleared promptly on a synced/Temp FS, so the reader blocks acquiring the lock.
- **Workaround:** set `HDF5_USE_FILE_LOCKING=FALSE` in the checker's environment
  (`$env:HDF5_USE_FILE_LOCKING="FALSE"`) → opens instantly (parity 80/80). Apply to
  ALL `.ladruno`/`.mpco` read steps after a recorder run. Learned 2026-05-31.
