---
wp: WP-255
title: "LadrunoForceBeamColumn: \"evalBulkState() -- could not invert flexibility\" after a shear failure, section strains ~1e14"
pr: "#1151"
date: 2026-11-25
---
### LadrunoForceBeamColumn: "evalBulkState() -- could not invert flexibility" after a shear failure, section strains ~1e14 (WP-255 adds `-monolithic`)
- **Bites:** a `LadrunoForceBeamColumn3d` with `LadrunoNDFiber` sections of `LadrunoMCFTStrip` fibres and `LadrunoMCFTSteel` bars (WP-254), cycled past the shear failure (Sezen & Moehle 2002 2CLD12, 2CLH18; Lynn et al. 1996 3CMD12, 3CMH18). The analysis stops with `could not invert flexibility`; sub-stepping, other algorithms and smaller tolerances do not help.
- **Why:** the nested force-based state determination needs the section flexibility `f_s = K_s⁻¹` at every inner iteration. Once every strip of a section is on a flat branch and every bar is past fracture, `K_s` is exactly singular: `f_s` does not exist, the inner iterate runs away (section strains ~1e14) and the element flexibility cannot be inverted. Making the singular branch fall through does not complete the runs; the iteration itself is the limit.
- **Fix/status (WP-255):** `element LadrunoForceBeamColumn … -monolithic` solves the section strains and end forces in one Newton on `K_s` and never inverts a section; only the whole element matrix must be regular. With `-monoAccept 10e3 10e6` all 8 lab tests complete. `evalBulkState` returning -1 on a singular trial flexibility, instead of falling to its initial-flexibility branches, is a separate defect of the nested path.
- See [[200_ladruno_frame_elements_apegmsh_emitter_guide]] §2.3 and ADR `255_lfbc_monolithic_solve.md`.
