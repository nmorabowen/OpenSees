---
wp: WP-142
title: "ASDSteel1D -auto_regularization in a shell reads the SHELL's size — getCharacteristicLength()/2, a quarter of the shortest edge under ASDShellQ4 EAS — not the…"
legacy_seq: 522
---
### ASDSteel1D `-auto_regularization` in a shell reads the SHELL's size — `getCharacteristicLength()/2`, a quarter of the shortest edge under ASDShellQ4 EAS — not the bar direction or the tie spacing (WP-142)
`ASDSteel1DMaterial::setTrialStrain` sets `params.lch_element = ops_TheActiveElement->getCharacteristicLength()/2` on every call (`ASDSteel1DMaterial.cpp:2150`). Wrapped in `PlateRebar` inside a LayeredShell, the active element is the shell: `ASDShellQ4::getCharacteristicLength` returns the shortest node-to-node distance, halved again when EAS is on (the default), so the bar sees a quarter of the shortest element edge — isotropic, blind to the bar direction and to the tie/stirrup spacing. With `-auto_regularization` that length scales the fracture softening (`epl_max = eupl + eupl·16r/(2·lch_element)`) and switches the buckling RVE's elastic correction (`lch_element > length`); without the flag it is unused.
- **Rule:** for bars in shells use `-buckling $lch` with `$lch` = the tie/stirrup spacing (the RVE takes `lch/2`) and do NOT pass `-auto_regularization`; otherwise mesh refinement changes the bar's softening and buckling.
- **Workaround/status:** ⚠️ documented from source (WP-142); no gate. *2026-09-27.*
