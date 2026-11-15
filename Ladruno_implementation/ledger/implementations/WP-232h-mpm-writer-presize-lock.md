---
wp: WP-232h
title: "WP-232h — MPM particle writer: presized root, guard lock (#1085)"
date: 2026-11-15
status: "draft"
section: "table"
---
| **WP-232h — MPM particle writer: presized root header and the guard lock** (#1085; ports WP-233's recorder fixes to `mpm output`). The particle file's root is presized and its links stay compact (`plist::presizeLinks`), so a kill inside an `mpm reset` layout can tear only the stage being laid out, never the root and every earlier stage. HDF5 file locking is off for the particle file and the writer holds the recorder's guard lock (Windows `LockFileEx` at 2^62, POSIX `flock(LOCK_EX)`) for the file's life, so a reader with locking on is refused cleanly and never breaks the SWMR writer. | MPM IO robustness | **none — no new classes** | `SRC/domain/pattern/ladrunoMPM/LadrunoMPMIO.cpp`, `SRC/recorder/Ladruno_Hdf5.h`, `Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py`, `tests/`, `Ladruno_implementation/232a_mpm_contracts.md §2.5` | draft | — |
