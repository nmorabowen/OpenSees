---
wp: WP-232a
title: "The Ladruno HDF5 attribute writer skips empty strings and empty vectors, so an attribute whose value is \"\" is ABSENT in a .ladruno (WP-232a)"
date: 2026-11-13
---
### The Ladruno HDF5 attribute writer skips empty strings and empty vectors, so an attribute whose value is `""` is ABSENT in a `.ladruno` (WP-232a)
- **Bites:** `h5::attribute::write(obj, name, const std::string&)` in `SRC/recorder/Ladruno_Hdf5.h`
  returns 0 without writing when the string is empty, and the `std::vector<int>` / `std::vector<double>`
  overloads do the same for an empty vector. A result group whose `DIMENSION` is `""` (the modal
  groups, a dimensionless MPM field such as `DEFORMATION_GRADIENT`) has no `DIMENSION` attribute at
  all. A reader that does `attrs["DIMENSION"]` raises `KeyError`; a fixture written by h5py with
  `attrs["DIMENSION"] = ""` hides the bug, because h5py does write the empty string.
- **Also:** the same wrapper writes "scalar" ints and doubles as 1-element 1-D arrays (`writei` /
  `writed` with `data_size = 1`) and strings as fixed-length `NULLTERM` in a 1-element array
  (`writes`), not as scalar dataspaces or variable-length strings. `ladruno_format._attr` already
  squeezes the arrays; apeGmsh's readers must accept both encodings.
- **Workaround/status:** treat an absent string attribute as `""` wherever `""` is a legal value
  (232a contract §0), and write test fixtures through `make_synthetic._ai/_ad/_as`, which mimic the
  C++ wrapper (including the skip), so readers are built against what the fork really writes.
