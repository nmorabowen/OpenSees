---
wp: WP-138
title: "A failed analyze() reverts the domain — the failing step's trial iterates are gone; to see the first one, commit it on purpose with test FixedNumIter 1 at a wa…"
legacy_seq: 510
---
### A failed `analyze()` reverts the domain — the failing step's trial iterates are gone; to see the first one, commit it on purpose with `test FixedNumIter 1` at a wall (WP-138)
- **Bites:** a post-mortem of a step that walled wants the strain increment the material could not integrate, but `StaticAnalysis::analyze()` calls `revertToLastCommit()` on every failure path, so after the return every material's trial state equals the committed one (the TIMs F20a observation that `substeps` read 0 after a failed step is the same effect).
- **Workaround:** at the wall, where the run is over anyway, read the committed state, then `integrator LoadControl -ds_first; test FixedNumIter 1; algorithm Newton; analyze 1`. The test always "converges", so the first Newton iterate is COMMITTED and its strains can be read; the difference to the committed state is that iterate's real Δε. If the material refuses at iterate 1 the step still fails; halve ds and retry. Destroys the equilibrium state, so only at the end of a run (`footing_ab.py` post-mortem). Cumulative per-instance censuses (`substepStats`, `sasStats`) survive the revert and do count the failed attempts. *2026-09-27.*
