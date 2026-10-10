---
wp: PR-877
title: "Re-evaluating a COMMITTED LadrunoRCConcrete / ASDConcrete3D state whose last step was actively damaging is NOT idempotent — a database round trip (or any resta…"
pr: "#877"
date: 2026-10-10
---
### Re-evaluating a COMMITTED LadrunoRCConcrete / ASDConcrete3D state whose last step was actively damaging is NOT idempotent — a database round trip (or any restart) changes the stress
- **Bites:** `database_roundtrip` with a stress probe fails by ~2 % on σ_yy (0.3803 → 0.3738) after a path whose final step still grew the tensile damage (e.g. transverse compression added to a cracked state: the Lubliner tensile measure grows through √J2). The same path ending on an elastic step round-trips bit-for-bit.
- **Why:** the committed effective stress is the post-reduction one, `σ̃_c = (1−d_t^pl)·ST + (1−d_c^pl)·SC`, not the elastic predictor that produced the committed thresholds. Re-running the return map at the committed strain decomposes `σ̃_c` (not the predictor), so the nominal stress `(1−d̄_t)·ST(σ̃_c)` differs by a factor `(1−d_t^pl)` on the tensile cone. Inherited from `ASDConcrete3D::compute`; not introduced by C2 (it showed up there first because `-crackedNu 0` removes the Poisson unloading that made the default path's last step elastic).
- **Workaround:** end a round-trip test path on an elastic (unloading) step; do not expect restart-from-commit to reproduce the trial stress of an actively-loading step. `tests/test_ladrunoRCConcrete_crackednu_betac.py::test_serialization_carries_the_c2_params` documents it. 2026-09-26.
