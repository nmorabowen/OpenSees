---
wp: LEGACY
title: "...and the obvious repair — propagate any < 0 — is WRONG. Measured."
legacy_seq: 374
---
### ...and the obvious repair — propagate any `< 0` — is WRONG. Measured.

The first ADR-86b cut did exactly that and **broke two long-green gates**:
`test_ladrunoBrick_asdconcrete_bend.py::test_notched_bend_mesh_objectivity` and
`::test_ssp_hourglass_energy_fraction_in_cracked_band` died at
`Domain::update - domain failed in update`, load factor 605, after months of passing.

`ASDConcrete3DMaterial::setTrialStrain` **does** return negative — out of `compute()`
(`ASDConcrete3DMaterial.cpp:1607+`) — but it means *"an inner iteration missed, here is my best
state"*, not *"I did not integrate this increment"*. That is precisely the case the fork already
settled the other way in ADR-33/34 (the `LadrunoDispBeamColumn3d` `solveHingeJump` entry above):
**return 0 with the last iterate plus a loud warning**, because a failure code there makes softening
analyses fragile — a fixed-increment run dies at a kink the global Newton would have walked through.

**OpenSees has no convention separating the two meanings**, so ADR-86b introduced one for the fork:
`LADRUNO_MATERIAL_REFUSED` (`SRC/material/LadrunoMaterialStatus.h`, value `-33086`) means *the
increment was NOT integrated and the committed state is unchanged*. `LadrunoBrick::update()`
propagates **only that exact value**; every other non-zero code keeps whatever treatment it already
had. Adopting it cost **zero** behaviour change anywhere else, and the two ASD gates went green
again.

- **Rule for any future fork material with a real failure mode:** returning a bare `-1` buys you
  nothing (elements swallow it), and returning `-1` *and* teaching an element to honour `< 0` buys
  you the ASD regression. Return the sentinel.
- **`updateHypo()` / `updateFinite()` keep their pre-existing `< 0` tests.** They predate ADR-86b
  and were not touched; only the std/b-bar, SSP and the two URI paths use the sentinel.
