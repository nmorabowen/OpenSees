---
title: "WP-232a — MPM contracts: the particle set (Contract 1) and the .ladruno particle group (Contract 2)"
project: Ladruno
status: "draft — schema for review; no C++"
owner: nmora
wp: "WP-232a (program WP-232, issue #1075)"
related:
  - "[[232_ladruno_mpm_general_driver_adr]] (§7.2, §7.3, §4 'LadrunoRecorder and the .ladruno format', W13)"
  - "[[ladruno_schema_v1]] (the .ladruno layout this extends)"
  - "[[03_ladruno_recorder]] · [[164_ladruno_recorder_write_path]] · [[165_ladruno_recorder_stage_semantics]]"
  - "[[232a_apegmsh_adr_draft]] (the apeGmsh side; placement of Contract 1)"
tags: [mpm, contract, hdf5, ladruno-format, schema, wp-232a]
updated: 2026-11-13
---

# WP-232a — MPM contracts

Two HDF5 contracts between apeGmsh and the fork, written so that the fork's C++ writer and reader
(232c) and apeGmsh's Python and TypeScript readers (232d) can each be built from this page alone.

| | Contract 1 — particle set | Contract 2 — particle results |
|---|---|---|
| Direction | apeGmsh → fork | fork → apeGmsh |
| File | a standalone HDF5 file, `<stem>.particles.h5` (placement: the apeGmsh ADR draft, D1) | a group inside a `.ladruno` file: `FORMAT_VERSION = 1`, `FORMAT_MINOR_VERSION = 2` |
| Read by | `mpm particles -file {path}` (232c) | `_ladruno.py`, `apeGmshViewer/src/reader/ladruno*.ts` (232d) |
| Owner of the layout | the fork (this page and its validator) | the fork (as all of `.ladruno`, schema §7b.1, apeGmsh ADR 0122 D4) |
| Executable spec | `tests/mpm_contracts.py::validate_particle_set` | `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py::validate` (its `_validate_particles`); `tests/mpm_contracts.py::validate_particle_results` calls it |
| Fixture | `tests/data/mpm/particle_set_v1.h5` | `tests/data/mpm/particles_v1m2.ladruno` |

The fixtures are written by `tests/mpm_contract_fixtures.py` and checked by
`tests/test_mpm_contracts.py`, which also refuses one malformed variant per rule. Where this page and
the validator disagree, the validator is a bug: fix both in one PR.

**Ownership.** The fork owns both layouts. The fork's reader is what admits or refuses a particle
set, and other producers than apeGmsh (a user's own script) must be able to write one. apeGmsh is
the canonical producer of Contract 1 and the canonical consumer of Contract 2; it changes either
only through a fork PR, as for the rest of `.ladruno` (ADR 0122 D4). This is a proposal: the owner
confirms it (open question Q1).

---

## 0. Conventions common to both contracts

- **Exact dtypes.** Every dataset has exactly the dtype listed, little-endian (`<i4`, `<i8`, `<u4`,
  `<f8`, `<f4`). A reader refuses another dtype; it never casts silently. (h5py writes `<f8` for a
  `numpy.float64` array; C++ writes `H5T_IEEE_F64LE` / `H5T_STD_I64LE` / `H5T_STD_I32LE`.)
- **Attribute encodings.** A scalar attribute may be stored with a scalar dataspace *or* a 1-D
  dataspace of length 1; a string attribute may be fixed-length (`NULLTERM`) or variable-length, ASCII
  or UTF-8. Readers accept all of these. The reason: Python writers (h5py) produce scalar
  dataspaces and variable-length strings, while the fork's C++ wrapper
  (`SRC/recorder/Ladruno_Hdf5.h`, `h5::attribute::writei/writed/writes`) produces 1-element arrays
  and fixed-length strings. The Contract-1 fixture uses the first encoding, the Contract-2 fixture
  the second, so each reader is built against what its producer really writes.
- **Empty strings are absent.** The fork's `h5::attribute::write(obj, name, const std::string&)`
  returns without writing when the string is empty (and the `std::vector` overloads skip an empty
  vector). So an attribute whose value is `""` is **absent** in a fork-written file. Every string
  attribute below whose value can be empty says "absent ≡ `""`" (quirk
  `WP-232a-the-ladruno-hdf5-attribute-writer-skips-empty-strings`).
- **Integers that a browser reads.** `int64` values are kept within ±(2⁵³ − 1), so h5wasm's
  `BigInt64Array` converts to `Number` exactly; the viewer already refuses an unsafe integer
  (`ladruno.ts`, "is not a safe integer").
- **Storage.** Writers chunk every per-particle dataset along axis 0 (target ≈ 1 MiB per chunk,
  the WP-164 rule) and compress it with shuffle + deflate (level 1 by default). Readers must not
  depend on the chunk shape and must also accept contiguous datasets. Only the `shuffle`, `deflate`
  (gzip) and `fletcher32` filters are allowed: h5wasm and the fork's runtime-loaded HDF5 both read
  them, and nothing else is guaranteed (no szip, lzf, blosc, zfp, scale-offset) (C1-FILTER,
  C2-FILTER).
- **Units.** Neither file converts units. Both carry the producer's declared consistent system
  (`UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME`); mass is in `F·T²/L`, density in `F·T²/L⁴`, stress in
  `F/L²`. The fork copies them from Contract 1 into Contract 2.
- **Stress.** Cauchy stress, **tension positive**, in the global axes, Voigt order
  `11,22,33,12,23,13` in 3D and `11,22,33,12` in 2D plane strain (σ₃₃ is the out-of-plane normal
  stress, needed for a K0 start; the shear is the tensor shear σᵢⱼ — for stress tensor and
  "engineering" shear are the same number). This is OpenSees' `nDMaterial` 3D order
  (`ElasticIsotropic3D`: xx, yy, zz, xy, yz, zx); the plane-strain row adds σ₃₃ to OpenSees'
  3-component plane-strain vector.
- **Plane strain.** For `ndm = 2`, every volume and mass is per the declared thickness `t`
  (`PLANE_STRAIN_THICKNESS`): `VOLUME` is area × t, in `L³`, so `mass = density × VOLUME` holds in
  both dimensions.

---

## 1. Contract 1 — the particle set (apeGmsh → fork)

One HDF5 file, `FORMAT = "ladruno-mpm-particles"`, `FORMAT_VERSION = 1`. All paths below are
from the file root `/`. A file with any rule violated is refused by `mpm particles` with a message
that names the rule (the codes `C1-…` below; 232c's refusal text should carry them).

### 1.1 Root attributes

| Name | Type | Req. | Default | Rule |
|---|---|---|---|---|
| `FORMAT` | string | yes | — | `== "ladruno-mpm-particles"` (C1-FMT) |
| `FORMAT_VERSION` | int32 | yes | — | `== 1`; any other value is refused (C1-VER) |
| `GENERATOR` | string | yes | — | non-empty, e.g. `"apeGmsh"` (C1-GEN) |
| `SET_ID` | string | yes | — | a uuid4 in canonical 36-character lowercase form, minted when the file is written; the fork copies it into Contract 2 (C1-SETID) |
| `CREATED_UTC` | string | no | absent | ISO 8601, `YYYY-MM-DDThh:mm:ssZ` (informative only) |
| `NDM` | int32 | yes | — | `2` or `3` (C1-NDM) |
| `NUM_PARTICLES` | int64 | yes | — | `N ≥ 1`, equal to the length of every `PARTICLES/*` dataset (C1-N) |
| `PLANE_STRAIN_THICKNESS` | float64 | iff `NDM = 2` | — | finite, `> 0`, in `L`; must be absent when `NDM = 3` (C1-THK) |
| `UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME` | string | yes | — | non-empty unit symbols of one consistent system, e.g. `m`, `kN`, `s` (C1-UNIT) |
| `UNIT_SYSTEM` | string | no | absent | apeGmsh's unit-system id when it has one (`kN_m_s`, `tf_m_s`, `kgf_m_s`, `N_mm_s`, `N_m_s`, `kip_in_s`, `dyne_cm_s`; `apeGmshViewer/src/state/units.ts`); informative |
| `MASS_SOURCE` | string | yes | — | `"mass"` or `"density"`: which of `PARTICLES/MASS` and `PARTICLES/DENSITY` is present (C1-MASS) |
| `STRESS_SIGN` | string | yes | — | `== "tension-positive"` (C1-SIGN) |
| `VOIGT_ORDER` | string | yes | — | `== "11,22,33,12,23,13"` when `NDM = 3`, `== "11,22,33,12"` when `NDM = 2` (C1-VOIGT) |

`STRESS_SIGN` and `VOIGT_ORDER` admit one value each. They are in the file so that a producer
states the convention it used and a reader refuses a file written to another one, rather than
reading it wrongly.

### 1.2 `GRID` — the background grid (attributes on the group `/GRID`)

| Name | Type | Shape | Req. | Rule |
|---|---|---|---|---|
| `ORIGIN` | float64 | `[NDM]` | yes | finite; the grid's minimum corner, in `L` (C1-GRID) |
| `CELL_SIZE` | float64 | `[NDM]` | yes | finite, `> 0` per axis, in `L` (C1-GRID) |
| `NUM_CELLS` | int32 | `[NDM]` | yes | `≥ 1` per axis (C1-GRID) |
| `BASIS` | string | scalar | yes | `"gimp"` or `"linear"` (C1-BASIS) |

The grid box is `[ORIGIN, ORIGIN + NUM_CELLS·CELL_SIZE]` per axis, closed. Nodes sit at
`ORIGIN + i·CELL_SIZE`, `i = 0 … NUM_CELLS` per axis. The engine runs on this grid and the viewer
draws it; the seeding was done against it. A grid that differs from the one the deck's `mpm grid`
declares is refused by the fork (232c), not by the file validator.

### 1.3 `BODIES` — names for the body ids (datasets in `/BODIES`)

| Name | Type | Shape | Req. | Rule |
|---|---|---|---|---|
| `ID` | int32 | `[B]` | yes | unique, `≥ 0` (C1-BODY) |
| `NAME` | string (vlen UTF-8) | `[B]` | yes | unique, non-empty: the body's label in apeGmsh (C1-BODY) |

`B ≥ 1`. Every `PARTICLES/BODY_ID` value is one of `BODIES/ID` (C1-BODYREF). The names reach
Contract 2 unchanged, so the viewer names bodies without apeGmsh (apeGmsh ADR 0112 D2c).

### 1.4 `PARTICLES` — one row per particle (datasets in `/PARTICLES`)

| Name | Type | Shape | Units | Req. | Default | Rule |
|---|---|---|---|---|---|---|
| `ID` | int64 | `[N]` | — | yes | — | strictly increasing (so unique), `0 ≤ id ≤ 2⁵³ − 1` (C1-ID) |
| `POSITION` | float64 | `[N, NDM]` | `L` | yes | — | finite; the particle and its GIMP domain `x ± l` lie in the grid box (C1-POS) |
| `VOLUME` | float64 | `[N]` | `L³` | yes | — | finite, `> 0`; zero or negative is refused (C1-VOL) |
| `MASS` | float64 | `[N]` | `F·T²/L` | iff `MASS_SOURCE = "mass"` | — | finite, `> 0` (C1-MASS) |
| `DENSITY` | float64 | `[N]` | `F·T²/L⁴` | iff `MASS_SOURCE = "density"` | — | finite, `> 0`; the fork computes `mass = DENSITY × VOLUME` once, at load (C1-MASS) |
| `GIMP_HALF_LENGTH` | float64 | `[N, NDM]` | `L` | yes | — | finite, `0 < l ≤ CELL_SIZE/2` per axis; ignored by the engine when `BASIS = "linear"` (C1-GIMP) |
| `VELOCITY` | float64 | `[N, NDM]` | `L/T` | no | `0` | finite (C1-FINITE) |
| `STRESS` | float64 | `[N, 6]` (3D), `[N, 4]` (2D) | `F/L²` | no | `0` | finite; §0 convention (C1-FINITE) |
| `BODY_ID` | int32 | `[N]` | — | yes | — | in `BODIES/ID` (C1-BODYREF) |
| `MATERIAL_TAG` | int32 | `[N]` | — | yes | — | `> 0`; a tag **resolved at emit** from the body's material label (charter P4). The fork refuses a tag that names no `nDMaterial` defined before `mpm particles`, or one outside the admitted-material table (232c) (C1-MAT) |
| `FLAGS` | uint32 | `[N]` | — | no | `0` | v1 defines no bit: every value must be 0, and the fork refuses a set bit rather than ignore it (C1-FLAGS) |

- Every dataset has length `N = NUM_PARTICLES` on axis 0 (C1-N), and exactly the dtype and the
  second dimension listed (C1-DTYPE, C1-SHAPE). Values are finite (C1-FINITE). A missing required
  attribute, group or dataset anywhere in Contract 1 is C1-REQ.
- **No other dataset** may sit in `PARTICLES` (C1-UNKNOWN): an unknown per-particle field could
  carry meaning, so a reader refuses it instead of dropping it. A new field is a
  `FORMAT_VERSION` bump.
- **For lattice seeding** the GIMP domain fills the particle's volume:
  `VOLUME = ∏(2·l)` in 3D and `t·∏(2·l)` in 2D. This is not a rule (Gauss-point seeding, later,
  breaks it); the validator reports a mismatch as a warning only.

### 1.5 `PROVENANCE` (group, attributes only)

Producer-written, lowercase names. Readers ignore what they do not know.

| Name | Type | Req. | Meaning |
|---|---|---|---|
| `apegmsh_version` | string | when `GENERATOR = "apeGmsh"` | the apeGmsh version that wrote the file |
| `apegmsh_commit` | string | when `GENERATOR = "apeGmsh"` | 40 lowercase hex characters, or `"unknown"` |
| `session_id` | string | when `GENERATOR = "apeGmsh"` | the session's uuid4 (`/meta/session_id` of the paired `<stem>.h5`) |
| `fem_hash`, `model_hash` | string | no | copied byte-equal from the paired `<stem>.h5` `/meta/lineage`, when it has them |

(C1-PROV checks the three required ones.)

### 1.6 `SEEDING` (optional group, attributes only)

| Name | Type | Shape | Meaning |
|---|---|---|---|
| `METHOD` | string | scalar | `"lattice"`, `"gauss"` or `"custom"` |
| `PARTICLES_PER_CELL` | int32 | `[NDM]` | `n` per axis, for `"lattice"` |

Informative; the uniformity report is the producer's (232d).

### 1.7 Versioning

`FORMAT_VERSION = 1`. A reader refuses every other version. Any new field, any changed meaning, and
any new allowed value of an existing attribute bumps the version. There is no minor version for
Contract 1: the fork's reader refuses unknown fields, so nothing additive is silently safe.

---

## 2. Contract 2 — the particle group in `.ladruno` (fork → apeGmsh)

The particle group is an **additive** change to schema v1: `INFO/FORMAT_VERSION` stays **1** and
`INFO/FORMAT_MINOR_VERSION` goes from 1 to **2**. Everything in `ladruno_schema_v1.md` holds
unchanged; minor 2 adds `MODEL/PARTICLES` and `RESULTS/ON_PARTICLES` to a `MODEL_STAGE`.

### 2.1 The version rule

- **No major bump.** Schema §9: "Additive changes (new optional groups) keep the version"; a
  semantic or additive revision moves `FORMAT_MINOR_VERSION` (the fork writes both attributes,
  `LadrunoRecorder.cpp:989-994`). apeGmsh treats a `FORMAT_VERSION` bump as a breaking layout change
  and opens only `(1,)` (`results/schema/_versions.py:73-90`, `_ladruno.py:633-637`,
  `apeGmshViewer/src/reader/ladruno.ts:56`), while it tolerates additive changes within a version and
  already reads the minor (`ladruno.ts:172`, `:342`). A major bump would make every new `.ladruno`
  unreadable by today's apeGmsh. (An earlier draft of this page bumped `FORMAT_VERSION` to 2; the
  orchestrator's review withdrew that on 2026-11-13.)
- **`FORMAT_MINOR_VERSION = 2`** means "minor 1 (schema §2a: per-stage `KIND`, stage barriers), plus
  the particle group of this section". Minors are cumulative.
- **A file that holds a particle group writes minor `≥ 2`** (C2-MINOR). Whether the recorder also
  writes 2 on FE-only files is 232c's choice; the recommendation is no, so FE-only output stays
  byte-identical and the two fork checks that assert minor `== 1`
  (`Ladruno_scripts/ladruno_recorder_tests/stage_kind_check.py:99`,
  `tests/test_ladruno_recorder_stage_semantics.py:183`) keep holding. If 232c writes 2 everywhere, it
  relaxes both to `>= 1` in the same PR.
- **What an older reader does.** A reader written for minor 1 opens a minor-2 file and skips
  `ON_PARTICLES` and `MODEL/PARTICLES`, which it does not know. On a particle-only file it shows an
  empty model (zero nodes). That is the accepted price of an additive change; a reader that wants
  particles checks `FORMAT_MINOR_VERSION >= 2` and then looks for the groups.
- **apeGmsh** keeps `LADRUNO_SUPPORTED_FORMAT_VERSIONS = (1,)`; its readers read the particle group
  when the minor is `≥ 2` and the group is present (232d). apeGmsh's own `.ladruno` writers
  (`_ladruno_writer.py`, the MPCO converter) keep writing minor 1; they never write particles.
- **The fork's executable spec.** `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py` stays
  `FORMAT_VERSION = 1`, gains `FORMAT_MINOR_VERSION_PARTICLES = 2`, and `validate()` runs the particle
  checks of this section on every stage that has a particle group (`_validate_particles`). Its
  reader gains `particle_stages()` and `particle_output()`. `make_synthetic.build(path,
  particles=True)` adds a particle group to the synthetic FE file (the shared-file case), and the
  harness `test_ladruno_harness.py` checks both. 232a lands these; nothing else in those files changes.

### 2.2 `INFO` in a file with particles

As schema §2, with: `FORMAT_VERSION = 1`; `FORMAT_MINOR_VERSION = 2`; `SPATIAL_DIM = 2 | 3` (C2-NDM) equal to
the particle set's `NDM` (and to the nodes' dimension when the file also holds FE output);
`PARTITIONED = 0` — the v1 driver is serial, and a file with particles that is one part of a
partitioned set is refused (C2-PART). `STORED_PRECISION` keeps its schema-§2 meaning for
`ON_NODES` / `ON_ELEMENTS`; a particle field's precision is its own dataset's dtype (§2.6).

### 2.3 The `MODEL_STAGE` that holds particles

- `MODEL_STAGE[<id>]` follows schema §2a for its id and attributes. `KIND` is `"transient"` for the
  explicit engine (mode F) and for mode-A transient steps, `"static"` for mode-A quasi-static steps;
  never `"eigen"` (C2-KIND). A particle stage has no running `Analysis` object in mode F, so a writer
  that shares `LadrunoRecorder`'s stage logic must set `KIND` itself (232c).
- `MODEL/NODES` is present with `ID` and `COORDINATES`, zero-length when the run has no FE nodes;
  `MODEL/ELEMENTS` is present, possibly empty (schema §1 shape, as WP-165's empty stage).
- **A new stage opens** at a user `domainChange` barrier (schema §2a, rule 2) and whenever the grid
  changes (§2.4 is per stage). Particle removal or addition **never** opens a stage: it is recorded
  in the per-output id lists (§2.5).
- `MODEL/PARTICLES` and `RESULTS/ON_PARTICLES` are both present or both absent in a stage (C2-PAIR).
  One particle group per stage in minor 2.

### 2.4 `MODEL_STAGE[<id>]/MODEL/PARTICLES` — per-stage particle metadata

**Attributes**

| Name | Type | Req. | Rule |
|---|---|---|---|
| `SOURCE_SET_ID` | string | yes | `SET_ID` of the Contract-1 file the particles were loaded from; with several files, the first one's, and `SOURCE_SET_IDS` lists all (C2-SRC) |
| `SOURCE_SET_IDS` | string | no | a JSON array of every loaded `SET_ID`, when there is more than one |
| `PATTERN_TAG` | int32 | yes | the host pattern's tag (`LadrunoMPMPattern`); `-1` in mode A, which has none |
| `MODE` | string | yes | `"explicit"`, `"implicit"` or `"quasi-static"` (C2-MODE) |
| `PLANE_STRAIN_THICKNESS` | float64 | iff `SPATIAL_DIM = 2` | copied from Contract 1 (C2-THK) |
| `UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME` | string | yes | copied from Contract 1 (C2-UNIT) |
| `UNIT_SYSTEM` | string | no | copied when Contract 1 has it |

**`GRID`** — subgroup with the same four attributes as Contract 1 §1.2 (`ORIGIN`, `CELL_SIZE`,
`NUM_CELLS`, `BASIS`), constant for the stage (C2-GRID).

**`BODIES`** — subgroup with `ID` (int32 `[B]`) and `NAME` (string `[B]`), copied from Contract 1
(C2-BODY). In a fork-written file `NAME` is a fixed-length string dataset; readers accept fixed or
variable length.

**The registry** — datasets directly in `MODEL/PARTICLES`, one row per particle that appears in any
output of the stage:

| Name | Type | Shape | Units | Rule |
|---|---|---|---|---|
| `ID` | int64 | `[M]` | — | strictly increasing (C2-REG) |
| `BODY_ID` | int32 | `[M]` | — | in `BODIES/ID` (C2-REG) |
| `MATERIAL_TAG` | int32 | `[M]` | — | `> 0` (C2-REG) |
| `MASS` | float64 | `[M]` | `F·T²/L` | finite, `> 0`; constant per particle (MPM conserves it), so it is stored once here, not per output (C2-REG) |

The registry is extensible on axis 0: a particle added during the stage is appended. Because added
particles take ids above every id used before (§2.7), appending keeps `ID` sorted.

### 2.5 `MODEL_STAGE[<id>]/RESULTS/ON_PARTICLES` — the outputs

One output is the state of the live particles after one step. Outputs are **ragged**: their particle
count can change (removal, addition), so they are concatenated along one axis with an offset table —
the CSR form schema §3.4 already uses for connectivity — instead of schema §7.0's fixed
`[T × nIds × nComp]`.

**Attribute:** `PARTITION_REDUCTION = "NONE"` (schema §7.1; serial only in minor 2) (C2-RED).

**Index datasets** (directly in `ON_PARTICLES`; `T` outputs, `S = OFFSETS[T]` rows in all):

| Name | Type | Shape | Rule |
|---|---|---|---|
| `STEP` | int32 | `[T]` | the Domain commit tag at the output; strictly increasing (C2-AXES) |
| `TIME` | float64 | `[T]` | the Domain time at the output; finite, non-decreasing (C2-AXES) |
| `OFFSETS` | int64 | `[T + 1]` | `OFFSETS[0] = 0`, non-decreasing, `OFFSETS[T] = S` (C2-OFF) |
| `ID` | int64 | `[S]` | output `k`'s particle ids are `ID[OFFSETS[k] : OFFSETS[k+1]]`, strictly increasing within an output, each one in the registry (C2-ID) |

`T ≥ 1` (C2-AXES). All four are extensible (`maxshape = None` on axis 0).

**Field groups** — one subgroup per recorded field, `ON_PARTICLES/<FIELD>/`:

| Item | Type | Rule |
|---|---|---|
| attr `DISPLAY_NAME` | string | non-empty |
| attr `COMPONENTS` | string | comma-separated labels; their count is `nComp` (C2-COMP) |
| attr `DIMENSION` | string | absent ≡ `""` (dimensionless) |
| attr `DESCRIPTION` | string | absent ≡ `""` |
| attr `TYPE` | int32 | `0` (Generic, `ResultType` in `Ladruno_Types.h`) |
| attr `DATA_TYPE` | int32 | `0` Scalar, `1` Vectorial, `2` Tensorial (`ResultDataType`) |
| attr `PARTITION_REDUCTION` | string | `"NONE"` |
| dataset `DATA` | float64 or float32 | `[S, nComp]`, row `r` belongs to particle `ID[r]` (C2-DATA) |

A missing `DISPLAY_NAME`, `COMPONENTS`, `TYPE`, `DATA_TYPE` or `PARTITION_REDUCTION` is C2-ATTR.
The field set is fixed for a stage: every field's `DATA` has exactly `S` rows. `POSITION` is
required (C2-POS); every other field is optional. Field names match `^[A-Z][A-Z0-9_]*$`, or
`^STATE_[A-Za-z0-9_]+$` for material state (C2-NAME). Every subgroup of `ON_PARTICLES` is a field;
no dataset other than the four index datasets sits there (C2-UNKNOWN).

**Reading output `k` of field `F`** (Python; the TypeScript reader is the same with hyperslabs):

```python
a, b = OFFSETS[k], OFFSETS[k + 1]
ids  = ON_PARTICLES["ID"][a:b]                   # sorted
vals = ON_PARTICLES[F]["DATA"][a:b]              # [b - a, nComp]
row  = np.searchsorted(PARTICLES["ID"][...], ids)  # registry rows: body, material, mass
```

A field's full history for one particle is a gather over outputs; a dense `[T, M, nComp]` view with
NaN for absent particles is a reader convenience, not the stored layout.

### 2.6 The field catalogue

Reserved names. When one of these is written, its `COMPONENTS` must be exactly as listed (C2-CAT).
All values are in the global axes. POSITION is float64 always; every other field is float64 unless
the user asks for float32 output (a lossy opt-in, as `-precision f32` on the recorder) (C2-DTYPE).

| Field | `COMPONENTS` (ndm = 2) | `COMPONENTS` (ndm = 3) | `DIMENSION` | `DATA_TYPE` | Meaning |
|---|---|---|---|---|---|
| `POSITION` | `x,y` | `x,y,z` | `L` | 1 | current position (required) |
| `VELOCITY` | `Vx,Vy` | `Vx,Vy,Vz` | `L/t` | 1 | particle velocity after G2P |
| `ACCELERATION` | `Ax,Ay` | `Ax,Ay,Az` | `L/t^2` | 1 | mode A only |
| `STRESS` | `sigma11,sigma22,sigma33,sigma12` | `sigma11,sigma22,sigma33,sigma12,sigma23,sigma13` | `F/L^2` | 2 | Cauchy, tension positive (§0) |
| `VOLUME` | `V` | `V` | `L^3` | 0 | current volume, `J·V₀` (2D: area × t) |
| `GIMP_HALF_LENGTH` | `lx,ly` | `lx,ly,lz` | `L` | 1 | current half-lengths, when they evolve |
| `DEFORMATION_GRADIENT` | `F11,F12,F21,F22` | `F11,F12,F13,F21,F22,F23,F31,F32,F33` | absent (`""`) | 2 | row-major; in 2D `F33 = 1` is not stored |
| `STATE_<token>` | the material's response labels | same | as the material says | 0 | a material state variable (below) |

- **Values are finite** in every field except `STATE_*` (C2-FINITE).
- **`STATE_<token>`**: `<token>` is the material response token the user named in `mpm output`
  (232c's guide lists the admitted ones, e.g. an equivalent plastic strain). A particle whose material
  does not answer the token gets `NaN` in every component — the only place `NaN` is allowed. New
  tokens need no version bump.
- A new reserved field (e.g. a pore pressure in phase 6) is additive: it is a new optional group with
  its own name, which a minor-2 reader that does not know it must skip, not refuse. Readers skip
  unknown field groups; the validator refuses only malformed ones.

### 2.7 Particle ids — stable, never reused

- A particle keeps the id it had in Contract 1 for the whole run, in every stage and output.
- **Removal** (a particle leaves the grid, or is deleted): its id is absent from every later output
  of the run, in this stage and in later stages. It never reappears (C2-REUSE). Its registry row
  stays in the stage where it lived.
- **Addition** (inflow, a second `mpm particles` load): the engine gives each new particle an id
  above every id the run has used, and appends it to the registry. A second Contract-1 file whose
  ids meet any id the run has used is refused (232c).
- `mpm reset` returns the particles to their loaded state; a writer opens a new stage after it, and
  the ids are the loaded ones again (the reset is a new run for this rule).

### 2.8 Output cadence and field selection (what the writer does; the command is 232c's)

- **Never every step by default.** 10⁶ particles × about 20 doubles is about 160 MB per output.
  `mpm output` selects every `N` steps or every `Δt_out` of time, and the fields; `POSITION` is
  always written.
- **First and last.** At the first `mpm analyze` of a stage the writer outputs the state before the
  first step (so the viewer has the initial configuration), and at the end of every `mpm analyze`
  it outputs the last step if the cadence skipped it, so the final state is never lost.
- **Synchrony.** An output is written from `Domain::commit`'s recorder pass or after it, so
  positions and stresses come from the same step (ADR §5.2 step 8, W11).
- **Flush.** On the recorder's wall-clock cadence (`-flush`, WP-164) and at stage end and close.
- **Chunks.** `ID` and each `DATA`: rows per chunk `≈ 2²⁰ / (nComp · itemsize)`, shuffle + deflate
  at the `-compress` level (default 1). `STEP`, `TIME`, `OFFSETS`: 1024 per chunk.

### 2.9 Executable spec and where this goes when 232c lands

- `ladruno_format.validate(path)` checks every `C2-…` rule above on each stage with a particle
  group (`_validate_particles`), together with every schema-v1 invariant of the rest of the file.
  `tests/mpm_contracts.py::validate_particle_results(path)` is a thin wrapper that also insists the
  file has a particle group at all. apeGmsh SHOULD import `ladruno_format`, as for v1 (schema §9).
- When 232c lands the writer, its format-document PR moves §2 into `ladruno_schema_v1.md` (a new
  §7.5 `ON_PARTICLES`, a §2 line for minor 2 and a §9 versioning line), with this page linked from
  there.

---

## 3. Round trip between the contracts

The fixtures are one model seen from both sides: `particles_v1m2.ladruno` is what the fork would write
after loading `particle_set_v1.h5`. `tests/test_mpm_contracts.py` checks the bridge between them:

| Contract 1 | becomes in Contract 2 |
|---|---|
| `SET_ID` | `MODEL/PARTICLES@SOURCE_SET_ID` |
| `NDM`, `PLANE_STRAIN_THICKNESS`, `UNIT_*` | `INFO/SPATIAL_DIM`, `MODEL/PARTICLES@PLANE_STRAIN_THICKNESS`, `@UNIT_*` |
| `GRID@*` | `MODEL/PARTICLES/GRID@*` |
| `BODIES/{ID, NAME}` | `MODEL/PARTICLES/BODIES/{ID, NAME}` |
| `PARTICLES/{ID, BODY_ID, MATERIAL_TAG}` | the registry, same rows |
| `PARTICLES/MASS`, or `DENSITY × VOLUME` | registry `MASS` |
| `PARTICLES/POSITION`, `VOLUME`, `VELOCITY`, `STRESS` | output 0's `POSITION`, `VOLUME`, `VELOCITY`, `STRESS` |

## 4. Fixtures

`python3.12 tests/mpm_contract_fixtures.py` rewrites both files. The model: 2D plane strain,
`t = 1 m`, `kN_m_s`; a 4 × 4 grid of 0.5 m cells, GIMP; two bodies — `soil` (8 particles, material
tag 2, K0 initial stress, `ρ = 2.0`) and `block` (4 particles, tag 1, `ρ = 7.85`) — seeded on a 2 × 2
lattice per cell. Contract 2 holds one transient stage with 3 outputs (steps 0, 10, 20); particle 11
is removed after output 1; the fields are `POSITION`, `VELOCITY`, `STRESS`, `VOLUME`,
`DEFORMATION_GRADIENT` (no `DIMENSION`, as the fork writes it) and `STATE_eqPlasticStrain` (`NaN`
for the block's elastic particles). The token `eqPlasticStrain` is illustrative; 232c's guide names
the real ones.

## 5. Open questions (for the owner and the apeGmsh maintainer)

- **Q1 — ownership of Contract 1.** Proposed: the fork owns it (its reader admits or refuses), and
  apeGmsh changes it through a fork PR. The alternative, apeGmsh owning it as one of its zones, would
  bring apeGmsh's zone rules (semver key in `/meta`, its int32 policy) to a file the fork parses.
- **Q2 — the ADR's `FORMAT_VERSION` 1 → 2.** ADR §7.3 and architecture.md say the version goes
  from 1 to 2. This page keeps 1 and moves the minor to 2 (§2.1, the orchestrator's ruling of
  2026-11-13); the ADR text needs the matching amendment. Also open: whether the recorder writes
  minor 2 on FE-only files (recommended: no).
- **Q3 — int64 ids.** The ADR fixes int64; the viewer converts with a safe-integer check, and the
  rule `id ≤ 2⁵³ − 1` keeps that exact. apeGmsh's int32 policy for its ADR-0112 zones does not apply
  to either file (neither is a `model.h5` zone) — confirmed only if the apeGmsh ADR draft's D1 is
  accepted.
