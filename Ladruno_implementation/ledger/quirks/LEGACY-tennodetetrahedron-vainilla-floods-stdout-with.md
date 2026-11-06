---
wp: LEGACY
title: "TenNodeTetrahedron (vainilla) floods stdout with one line per shape-function entry -- silently kills large-mesh MPI jobs on a shared cluster long before any re…"
legacy_seq: 382
---
### `TenNodeTetrahedron` (vainilla) floods stdout with one line per shape-function entry -- silently kills large-mesh MPI jobs on a shared cluster long before any real memory limit is hit
- **Bites:** any MPI run of `TenNodeTetrahedron` past a few hundred thousand elements. Discovered
  during ADR-88's T3 cross-element H5DRM validation (`drm_load_pattern/86_drm_free_field_all_elements.ipynb`):
  the `TenNodeTetrahedron_h2.5` run (565,248 elements, 4 Gauss points, 4x10 shape-function table)
  was assumed to be dying from Mumps OOM (the same SIGKILL/exit-137 signature as three sibling
  runs in the same batch) -- root-caused instead by inspecting the raw SLURM `.out` file: **86.6
  million lines, 761 MB**, all bare floating-point numbers, one per shape-function entry per Gauss
  point per element (565,248 x 4 x 4 x 10 ~= 90.4M, matching the observed count almost exactly).
- **Why:** `TenNodeTetrahedron.cpp` (`computeBasis`/shape-function precompute loop, upstream
  commit `887ea413ef`, Jose A. Abell, 2024-03-25) has an unconditional `std::cout << shp[p][q] <<
  std::endl;` right next to the intended `Shape[p][q][count] = shp[p][q];` assignment -- almost
  certainly a leftover interactive debug line that was never gated behind a verbosity flag or
  removed. Harmless on the small decks this element was previously exercised with (a few thousand
  elements); catastrophic at the mesh sizes a real DRM free-field model needs -- the sheer I/O
  volume (and/or whatever buffering/flush behavior it triggers under MPI on a shared filesystem)
  is enough by itself to get the job SIGKILLed well before Mumps ever gets a chance to run out of
  memory, masquerading as the OOM failures already documented for tet10-scale meshes elsewhere in
  this fork's history (see `06_drm_free_field_bezier_tet.ipynb` in `soil_model_01_ATLAS/history/`).
- **Workaround/status (2026-08-30, ADR-88 -- FIXED):** the stray `std::cout` line removed and
  marked `// Ladruno` (vanilla file, per the fork's edit-marking discipline) — see
  `SRC/element/tetrahedron/TenNodeTetrahedron.cpp`. `Shape[p][q][count] = shp[p][q]` (the real
  computation the loop exists for) is untouched. Rebuilt `OpenSeesMP` incrementally (single
  `.cpp` recompile + relink, confirmed via the build log) on the ADR-88 isolated build
  (`/mnt/deadmanschest/pxpalacios/opensees_tmp/`, branch `adr88-h5drm-higher-order-elements`)
  and redeployed to `bin/OpenSeesMP`; `TenNodeTetrahedron_h2.5` resubmitted (job 145758) to
  confirm the fix under real MPI load. **General lesson**: when several sibling jobs in the same
  batch die with an identical SIGKILL signature, do not assume they share one root cause — check
  each one's own raw output before applying the same fix (more nodes, more memory) to all of them;
  here 3 of 4 were genuine OOM candidates and 1 was an unrelated, much cheaper bug hiding behind
  the same exit code.
