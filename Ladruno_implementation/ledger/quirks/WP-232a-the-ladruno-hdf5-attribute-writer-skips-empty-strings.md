---
wp: WP-232a
title: "The Ladruno HDF5 attribute writer skips empty strings and empty vectors, so an attribute whose value is \"\" is ABSENT in a .ladruno (WP-232a)"
pr: "#1077"
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
- **Datasets too:** the `h5::dataset::createAndWrite` overloads (`Ladruno_Hdf5.h:241-282`) return
  `HID_INVALID` for an empty `std::vector` or an empty string and create no dataset, so a zero-length
  dataset (an empty output, an empty registry) needs a direct `H5Dcreate`, as the WP-165 empty-stage
  branch of `writeModelNodes` does. The wrapper also has **no int64 dataset writer** (only
  `createAndWritei1/i2`, `H5T_STD_I32LE`) and **no string-array dataset writer** (`createAndWrites`
  writes one string). The MPM particle group needs int64 `ID` / `OFFSETS` / registry `ID` and a
  string `BODIES/NAME` dataset (232a contract §2.4-§2.5), so 232c must add those writers.
- **h5py trap when mimicking the C++ strings:** writing a fixed-length `NULLTERM` type of exactly
  `len(s)` bytes through h5py's normal path (`attrs.create(..., dtype=<NULLTERM type>)` with a numpy
  `S` array) converts NULLPAD → NULLTERM and **drops the last byte** (`"NONE"` reads back `"NON"`),
  because the size leaves no room for a terminator. The C++ `H5Awrite(attr, atype, …)` passes the
  file type as the memory type, so nothing is converted. Mimic it with the low-level
  `h5a.create(...).write(data, mtype=tid)`, as `make_synthetic._as` does.
- **Workaround/status:** treat an absent string attribute as `""` wherever `""` is a legal value
  (232a contract §0), and write test fixtures through `make_synthetic._ai/_ad/_as`, which mimic the
  C++ wrapper (including the skip), so readers are built against what the fork really writes.
