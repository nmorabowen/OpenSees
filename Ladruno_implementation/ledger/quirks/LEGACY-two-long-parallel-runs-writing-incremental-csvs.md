---
wp: LEGACY
title: "Two long parallel runs writing incremental CSVs: re-running one leg's name TRUNCATES the live file under it, and the survivor keeps writing at its old offset"
legacy_seq: 255
---
### Two long parallel runs writing incremental CSVs: re-running one leg's name TRUNCATES the live file under it, and the survivor keeps writing at its old offset
- **Bites:** any campaign whose legs stream results to `f = open(path, "w")` +
  `flush()` per step and run for hours in parallel — the `hypo_bearing` legs, and
  the same idiom in the perf testbeds. Starting a short smoke of a leg whose
  full run is ALREADY in flight reopens the same path with `"w"`. The smoke
  truncates the file to zero; the live process still holds its own descriptor at
  (say) byte 12000, so its next write lands there and the OS zero-fills the gap.
  Measured result: a 15 753-byte CSV that was 14 589 NUL bytes, holding the
  smoke's 11 rows, then padding, then the real run's tail — the live leg's
  entire early backbone (through s/B = 3.37%, ~80 min of compute, including the
  1% and 2.5% checkpoints) simply gone.
- **Tell:** NUL bytes in a CSV; a first data row whose settlement is *larger*
  than rows further down; a file far bigger than its line count justifies. It
  does NOT crash, and the live process's end-of-run summary still prints correct
  numbers from its in-memory rows — so the loss is silent unless you read the file.
- **Rule:** two guards, both cheap. (1) A capped/smoke run must write a
  DIFFERENT filename (`backbone_<leg>__smoke.csv`), so a smoke can never address
  a real leg's output. (2) Refuse to open an existing output file modified within
  ~180 s — that means another process is actively appending — with an explicit
  env override for the case where you know it is dead. Note Windows does not
  block the second open, so nothing but your own guard prevents this.
  *2026-07-30 (ADR-79 locking leg).*
