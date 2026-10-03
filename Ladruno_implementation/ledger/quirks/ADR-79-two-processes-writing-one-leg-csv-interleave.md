---
wp: ADR-79
title: "Two processes writing one leg CSV interleave rows and tear a line — the ADR-79 runner's 180 s guard exists for this, and a new testbed driver that omits it wil…"
legacy_seq: 438
---
### Two processes writing one leg CSV interleave rows and tear a line — the ADR-79 runner's 180 s guard exists for this, and a new testbed driver that omits it will lose a leg
- **Bites:** a leg's CSV has more rows than the run reported steps, a line with
  the wrong field count in the middle, and a `wall_s` column that jumps
  backwards. The reducer then quietly reports a curve neither process wrote.
- **Why:** two invocations of the same leg (an accidental double launch of a
  batch script) open the same `out/f10_<leg>.csv` in `"w"` mode and both keep
  writing at their own offsets. `hypo_bearing/README.md` already records this for
  the ADR-79 runner, whose fix is to refuse a CSV another process touched in the
  last 180 s (`ADR79_FORCE=1` overrides).
- **Measured (ADR-92 F10, legs L and M):** the physics was unaffected — re-run
  single-process, leg M reproduced `s/B = 0.011286` to the digit — but both wall
  times were wrong (M 115 s contended vs 169 s alone) and leg L's CSV carried a
  7-field line at row 2767.
- **Workaround/status (2026-09-14):** the F10 driver now carries the same guard
  (`F10_FORCE=1` overrides) and its reducer drops torn lines. **Copy the guard
  into any new testbed runner that writes one file per named leg** — a batch
  script that can be launched twice is not a hypothetical.
