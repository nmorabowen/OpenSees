---
wp: LEGACY
title: "Ladruno recorder nested result-group names need the <display> parent pre-created"
legacy_seq: 38
---
### Ladruno recorder nested result-group names need the `<display>` parent pre-created
- **Bites:** writing an element result whose `schema.name` is nested
  (`"<display>/<bucket>"`, e.g. `stress/204-FourNodeQuad[201:0:0]`) under a parent
  that does NOT already contain the intermediate `<display>` group → the
  `H5Gcreate` in `createResultGroup` returns an invalid handle and the datasets
  underneath silently fail (HDF5-DIAG noise / missing data).
- **Why:** the recorder's group-creation property list (`h_group_proplist`) is a
  GCPL with `CRT_ORDER_TRACKED|INDEXED` but carries **no create-intermediate-group
  LINK property**, and `createResultGroup` passes `lcpl=H5P_DEFAULT`. `H5Gcreate`
  only creates the *final* path component; any intermediate must already exist.
- **Status/workaround:** the time-series `StreamingSink` path is safe because the
  recorder pre-creates `ON_ELEMENTS/<display>` at init. The `EnvelopeSink` is
  self-contained, so it pre-creates each prefix segment of `m_name` inside
  `writeEnvelope` (the per-flush `H5Ldelete` only deletes the leaf link, so the
  intermediate persists). This is why element envelopes were latently broken until
  PR #45 — only flat node/domain names had ever been written. Learned 2026-05-31.
