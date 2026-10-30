---
wp: WP-216
title: "TransformationFE overflowed a 64-int static buffer on elements with more than 64 nodes"
date: 2026-10-27
pr: "#1032"
---
### An element with more than 64 nodes under `constraints Transformation` corrupted the heap

- **Bites:** a LadrunoMembrane floor with its wall and column nodes as
  `-diaphragmNodes` easily has > 64 nodes. Under Transformation every
  assembly (K, M, residual; `analyze`, `eigen`) wrote past a static buffer;
  the process died later with 0xc0000374 (heap corruption), typically at the
  next `wipe`, far from the cause. Penalty, Lagrange and Plain were fine.
  Results computed before the crash may also be affected.
- **Why:** `TransformationFE::dofData` holds one DOF count per NODE
  (`numDOFs.setData(dofData, numGroups)`) but was sized `MAX_NUM_DOF` = 64,
  the per-DOF pool limit. The pooled tangent and residual objects already had
  a size check; this buffer did not. Upstream has the same code.
- **Fix/status:** fixed in WP-216 (grown to the largest node count; the thread
  scratch too). Regression: `test_transformation_fe_element_with_more_than_64_nodes`,
  which is NOT deterministic: heap corruption need not crash, so it can pass on
  an unfixed build. The gate tower crashed on every run before the fix.
  When a crash shows up only at `wipe`, suspect a static buffer sized by a
  different count than the one that indexes it. Learned 2026-10-27 (WP-216).
