---
wp: WP-222
title: "LadrunoMembrane's eager transfer build dominated a tower's model build"
pr: "#1061"
date: 2026-11-02
---
### `o.element('LadrunoMembrane', ...)` takes ~30 ms each and a 46-storey tower spends ~80 s in `apeSees.run` before any analysis

- **Where the time goes.** `Domain::addElement` → `update()` → `resolveAndBuildTransfer()` runs the
  whole transfer build eagerly as soon as the supporting frames exist: the validation CDT mesh
  (twice: once in the parser's `validateNodesAndBasis`, once in `setDomain`), the gravity
  `distribute_surface_load`, the diaphragm CDT mesh and `condense_lst`. The cost is O(1) in model
  size per membrane (YOO RevB T2, 2795 membranes: 80.9 s; padding the domain with 400 k dummy
  nodes changed nothing). In `distribute_surface_load` the `tributary` path called the
  nearest-feature lambda ~29 M times, each call doing a `std::hypot` per candidate and two heap
  allocations.
- **WP-222 fix, bit-identical:** (a) a memo of the three pure stages keyed on the exact bit image
  of their inputs (`LadrunoMembraneMemo.h`; repeated storeys are bit-identical in the panel's own
  2D frame, so T2's 2795 membranes need 1038 entries over the three stages); (b) the
  nearest-feature query screens candidates by squared distance and evaluates `hypot` only on the
  tie shortlist, allocation-free. T2: 80.9 s → 5.8 s; (b) alone 44.4 s.
- **Traps for the next edit of the kernel:**
  - The fast `nearest_features` returns a **reference to a reused buffer**. Never hold its `tied`
    across a second call (no call site does today).
  - A change to any input of a memoised stage must go into that stage's key, or the memo returns
    another panel's result. The keys are built next to the calls in `resolveAndBuildTransfer` and
    `buildPlaneBasis`. The stages must stay pure: no globals, no randomness, no dependence on the
    element tag. Tags are mapped after the memo.
  - `-DLADRUNO_MEMBRANE_KERNEL_REFERENCE` compiles the pre-WP-222 query and containers.
    `tests/test_wp222_membrane_build_memo.py` checks both builds, and `-noMemo`, bit for bit. Run
    it after touching either.
