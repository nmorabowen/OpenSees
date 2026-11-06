---
wp: WP-144
title: "Domain::recvSelf calls update() on every element after its recvSelf, so a material's TRIAL stress/tangent on the wire is dead data (WP-144 G2)"
date: 2026-11-06
---
### `Domain::recvSelf` calls `update()` on every element after its recvSelf, so a material's TRIAL stress/tangent on the wire is dead data (WP-144 G2)
- **Bites:** `Domain::recvSelf` (Domain.cpp ~4273; `Domain::addElement` does the same at ~497) calls `element->update()` for every restored element, which re-integrates each Gauss point from its restored COMMITTED state to the strain of the restored nodal displacements. For a restored material that is a ZERO increment, so the committed state, the trial STATE and every counter survive, but the trial stress and the trial TANGENT are RECOMPUTED (stress within a last bit; the tangent is the elastic one, not the saved consistent plastic one).
  - LadrunoNorSand first sent them anyway (42 doubles, `WIRE_LEN` 128). No openseespy call can tell a wire that carries them from one that does not: the round-trip tests could only hold the committed state (bit-exact) and a continuation from it.
- **Rule:** put on the wire only what the host's post-restore `update()` does NOT recompute: the committed State, the strains, the latch and the refusal record. Rebuild the trial stress/tangent in `recvSelf` from the restored trial State so a restored point is self-consistent before that `update()`. A test that "proves" a round trip through the tangent proves nothing here.
- **Workaround/status:** FIXED (WP-144): `LadrunoNorSand::sendSelf` sends 86 doubles (`WIRE_LEN` 86, one Vector under the dbTag; FE_Datastore size-uniqueness needs nothing else, see "FE_Datastore keys a sent Vector by its SIZE"), `recvSelf` recomputes `sigT`/`CT` from `sT`. The latch and the refusal cause DO cross the wire (`test_latch_survives_a_database_round_trip`).
