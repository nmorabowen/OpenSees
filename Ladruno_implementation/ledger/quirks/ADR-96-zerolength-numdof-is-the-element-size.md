---
wp: ADR-96
title: "ZeroLength::numDOF is the element size everywhere, not a \"count check\" (ADR-96)"
legacy_seq: 405
---
## `ZeroLength::numDOF` is the element size everywhere, not a "count check" (ADR-96)

`setDomain()` dispatches `numDOF`/`elemType` on `(dimension, ndf)` pairs and
every accessor, the `t1d` transformation, `d0`/`v0`, `commitSensitivity` and the
responses loop to `numDOF` or `numDOF/2`. Three sites subtract whole nodal
vectors (`disp2 - disp1` in `setDomain`, `update`, `getResponse`), which throws
on a (3,4) pair before any count check is reached. "Relax the count check" is
therefore not a one-line change; the passenger scatter (ADR-96 D4) is the
minimal one that keeps the vanilla path byte-identical.
