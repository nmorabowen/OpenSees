---
wp: LEGACY
title: "ops.logFile(path, '-noEcho') is how a Python driver captures opserr warnings — but you must redirect it away again before reading the file"
legacy_seq: 359
---
## `ops.logFile(path, '-noEcho')` is how a Python driver captures `opserr` warnings — but you must redirect it away again before reading the file

**Found 2026-09-05, ADR-90 WP-A2.**

Material-level warnings that matter for a soil deck — `ManzariDafalias`'s
`"stage-switch stress ratio ... exceeded the bounding surface M_c ... (Outside Bounding!)"`
M_c-inflation notice, and the low-p `"mean stress p = ... is below the floor ... CLAMPING"` notice
— go to `opserr`, not to Python. A driver that wants to ASSERT on them, rather than hope a human
reads the console, can do:

```python
ops.logFile(log_path, '-noEcho')      # opserr -> file, console silent
...  run the analysis ...
ops.logFile(other_path, '-noEcho')    # release the handle
n_outside = open(log_path, errors='ignore').read().count('Outside Bounding')
```

- The positive control that the capture is live: the file also contains the material's own
  construction echo, so an empty file means the redirect did not take, not that nothing warned.
- `-noEcho` is what silences the console; without it the file is written *and* echoed.
