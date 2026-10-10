---
wp: LEGACY
title: "CrdTransf::getLocalAxes base default zeros the axes — wiring localAxes on an element whose transform doesn't override it emits a *degenerate* frame"
legacy_seq: 46
---
### `CrdTransf::getLocalAxes` base default zeros the axes — wiring `localAxes` on an element whose transform doesn't override it emits a *degenerate* frame
- **Bites:** when extending the Ladruno `"localAxes"` (response id 30) coverage to
  more beam elements, `DispBeamColumn2dInt` *looks* trivially wireable — it owns a
  `crdTransf` and the copy-paste pattern compiles. But its transform is
  `LinearCrdTransf2dInt`, which does **not** implement `getLocalAxes`, so it falls
  through to the base `CrdTransf::getLocalAxes` (`SRC/coordTransformation/CrdTransf.cpp:126`)
  that simply `Zero()`s xAxis/yAxis/zAxis and returns 0.
- **Why it matters:** the recorder's `writeModelLocalAxes` records *any* element that
  answers `"localAxes"`. A non-null response carrying an all-zeros frame would be
  written/quaternion-converted as a **degenerate** orientation — strictly worse than
  the current behaviour, where a silent (no-response) element falls back to a clean
  identity quaternion.
- **Status/workaround:** `DispBeamColumn2dInt` is **deliberately left unwired**. To
  wire it correctly, first give `LinearCrdTransf2dInt` a real `getLocalAxes` (mirror
  `LinearCrdTransf2d::getLocalAxes`), then add the id-30 response. The standard-transform
  beams (Elastic/Force/Disp/Mixed/GradientInelastic Beam(Column)2d/3d) are all safe —
  their LinearCrdTransf2d/3d, PDelta*, and Corot* transforms all override `getLocalAxes`.
  Learned 2026-06-03.
