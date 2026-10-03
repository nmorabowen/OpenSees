---
wp: LEGACY
title: "AllASDInternalVariableTypes.h and AllASDHardeningFunctions.h have NO include guard"
legacy_seq: 414
---
### `AllASDInternalVariableTypes.h` and `AllASDHardeningFunctions.h` have NO include guard

Including either directly in a translation unit that also includes
`ASDPlasticMaterial3D.h` (which pulls both in) is a redefinition storm. Relevant to
any standalone syntax-check / pre-flight translation unit.
