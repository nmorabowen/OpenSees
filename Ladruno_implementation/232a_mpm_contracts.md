---
title: "WP-232a — MPM contracts: the particle set (Contract 1) and the .ladruno particle group (Contract 2)"
project: Ladruno
status: "draft — schema for review; no C++"
owner: nmora
wp: "WP-232a (program WP-232, issue #1075, PR #1077)"
related:
  - "[[232_ladruno_mpm_general_driver_adr]] (§7.2, §7.3, §4 'LadrunoRecorder and the .ladruno format', W13)"
  - "[[ladruno_schema_v1]] (the .ladruno layout this extends)"
  - "[[03_ladruno_recorder]] · [[164_ladruno_recorder_write_path]] · [[165_ladruno_recorder_stage_semantics]]"
  - "[[232a_apegmsh_adr_draft]] (the apeGmsh side; placement of Contract 1)"
  - "tests/mpm_interface.md on wp/232b-mpm-oracle (the engine's state layout and Voigt order)"
tags: [mpm, contract, hdf5, ladruno-format, schema, wp-232a]
updated: 2026-11-13
---

# WP-232a — MPM contracts

Two HDF5 contracts between apeGmsh and the fork, written so that the fork's C++ writer and reader
(232c) and apeGmsh's Python and TypeScript readers (232d) can each be built from this page alone.

| | Contract 1 — particle set | Contract 2 — particle results |
|---|---|---|
| Direction | apeGmsh → fork | fork → apeGmsh |
| File | a standalone HDF5 file, `<deck stem>.particles.h5`, written beside the deck that names it (apeGmsh ADR 0134 draft, D1) | a group inside a `.ladruno` file: `FORMAT_VERSION = 1`, `FORMAT_MINOR_VERSION = 2` |
| Read by | `mpm particles -file {path}` (232c) | `_ladruno.py`, `apeGmshViewer/src/reader/ladruno*.ts` (232d) |
| Owner of the layout | the fork (this page and its validator) | the fork (as all of `.ladruno`, schema §7b.1, apeGmsh ADR 0122 D4) |
| Executable spec | `tests/mpm_contracts.py::validate_particle_set` | `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py::validate` (its `_validate_particles`); `tests/mpm_contracts.py::validate_particle_results` calls it |
| Fixtures | `tests/data/mpm/particle_set_v1.h5` | `tests/data/mpm/particles_v1m2.ladruno` (all particles), `tests/data/mpm/particles_mixed_v1m2.ladruno` (an FE stage, then a particle stage) |

`tests/mpm_contract_fixtures.py` writes the fixtures and `tests/test_mpm_contracts.py` checks them.
That test also refuses at least one malformed variant for **every rule code on this page**: it
reads the codes from this file and fails on a code with no case. Where this page and the validator
disagree, the validator is a bug; fix both in one PR.

**Ownership.** The fork owns both layouts. The fork's reader admits or refuses a particle set, and a
producer other than apeGmsh (a user's own script) must be able to write one. apeGmsh is the
canonical producer of Contract 1 and the canonical consumer of Contract 2. It changes either one only
through a fork PR, as for the rest of `.ladruno` (ADR 0122 D4). The owner has not confirmed this yet
(Q1).

---

## 0. Conventions common to both contracts

- **Exact dtypes.** Every dataset has exactly the dtype listed, little-endian (`<i4`, `<i8`, `<u4`,
  `<f8`, `<f4`). A reader refuses any other dtype; it never casts silently.
- **Integer attributes are int32** unless a table says int64 (only Contract 1's `NUM_PARTICLES`).
- **Attribute encodings.** A reader accepts a scalar attribute stored either with a scalar dataspace
  or with a 1-D dataspace of length 1. It accepts a string attribute either fixed-length or
  variable-length, ASCII or UTF-8. Python writers (h5py) produce scalar dataspaces and
  variable-length strings. The fork's C++ wrapper (`SRC/recorder/Ladruno_Hdf5.h`,
  `h5::attribute::writei/writed/writes`) produces 1-element arrays and fixed-length `NULLTERM`
  strings of exactly the string's length (`:82`). The Contract-1 fixture uses the first encoding and
  the Contract-2 fixture the second, so each reader is built against what its producer really
  writes.
- **Empty strings are absent.** The fork's string writer returns without writing when the string is
  empty, and the vector writers skip an empty vector (`Ladruno_Hdf5.h:258-279`). So in a fork-written
  file, an attribute whose value is `""` is **absent**. Every string attribute below whose value can
  be empty says "absent ≡ `""`". See the quirk
  `WP-232a-the-ladruno-hdf5-attribute-writer-skips-empty-strings`, which also lists the dataset
  writers 232c must add.
- **Integers that a browser reads.** `int64` values stay within ±(2⁵³ − 1). h5wasm's
  `BigInt64Array` then converts to `Number` exactly, and the viewer already refuses an unsafe
  integer (`ladruno.ts`, "is not a safe integer").
- **Storage.** Writers chunk every per-particle dataset along axis 0, aiming for about 1 MiB per
  chunk (the WP-164 rule). They compress it with shuffle + deflate, level 1 by default. Readers must
  not depend on the chunk shape and must also accept contiguous datasets. Only the `shuffle`,
  `deflate` (gzip) and `fletcher32` filters are allowed, because h5wasm and the fork's runtime-loaded
  HDF5 both read them. Nothing else is guaranteed: no szip, lzf, blosc, zfp or scale-offset
  (C1-FILTER, C2-FILTER).
- **Units.** Neither file converts units. Both carry the producer's declared consistent system
  (`UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME`). Mass is in `F·T²/L`, density in `F·T²/L⁴`, stress in
  `F/L²`. The fork copies the units from Contract 1 into Contract 2.
- **Plane strain.** For `ndm = 2`, every volume and mass is per the declared thickness `t`
  (`PLANE_STRAIN_THICKNESS`): `VOLUME` is area × t, in `L³`, so `mass = density × VOLUME` holds in
  both dimensions.
- **The grid box** is `[ORIGIN, ORIGIN + NUM_CELLS·CELL_SIZE]` per axis, closed, with a tolerance of
  10⁻⁹ times the largest extent. Contract 1 requires each particle and its GIMP domain `x ± l` to lie
  in it (C1-POS). Contract 2 requires every output position to lie in it (C2-POS).

### 0.1 Stress and strain: the index map

Cauchy stress is **tension positive** and in the global axes. The Voigt order is the hypo kernel's,
and 232b's engine carries the same order (`tests/mpm_interface.md` §1, `nv = 6` in 3D and `nv = 3` in
plane strain). Stress carries the tensor shear; strain carries the **engineering** shear (γ = 2ε),
as the engine's `eps_hat` does.

| Index | 3D stress | 3D strain | Plane-strain stress | Plane-strain strain |
|---|---|---|---|---|
| 0 | σ₁₁ = σxx | ε₁₁ | σ₁₁ = σxx | ε₁₁ |
| 1 | σ₂₂ = σyy | ε₂₂ | σ₂₂ = σyy | ε₂₂ |
| 2 | σ₃₃ = σzz | ε₃₃ | σ₁₂ = σxy | γ₁₂ = 2ε₁₂ |
| 3 | σ₁₂ = σxy | γ₁₂ = 2ε₁₂ | — | — |
| 4 | σ₂₃ = σyz | γ₂₃ | — | — |
| 5 | σ₁₃ = σzx | γ₁₃ | — | — |

The labels write `13` for the zx component, as the fork's element stress labels do. It is the same
component as 232b's `zx` (σ₁₃ = σ₃₁).

In plane strain, **σ₃₃ is a separate, optional dataset** (`STRESS33`), not a fourth Voigt
component. The engine does not carry it (232b §1). It is the material's business:

- **Contract 1** may carry it, so that a K0 initial state is complete. The v1 engine refuses a
  nonzero σ₃₃ on a particle whose material is not elastic, until 232c's initial-stress mechanism
  supports it (232b Q2).
- **Contract 2** writes it only for a material that exposes it. No v1 material does, so v1 files do
  not have it.
- In 3D, σ₃₃ is in `STRESS` and `STRESS33` must not appear (C1-STRESS33; in Contract 2 it is a C2-CAT
  error).

---

## 1. Contract 1 — the particle set (apeGmsh → fork)

One HDF5 file, with `FORMAT = "ladruno-mpm-particles"` and `FORMAT_VERSION = 1`. All paths below
are from the file root `/`. `mpm particles` refuses a file that breaks any rule, and its message
names the rule's code; 232c's refusal text should carry the code. A missing required attribute,
group or dataset anywhere in Contract 1 is C1-REQ.

**Placement and path resolution.**

- **Who writes it.** apeGmsh's deck emitter writes the file **next to the deck it emits**, as
  `<deck stem>.particles.h5`, in the same call that writes the deck, whenever the model declares MPM
  bodies. Emitting `runs/case1.tcl` writes `runs/case1.particles.h5`.
- **How the deck names it.** Each deck emitter generates the `mpm particles -file …` line for its
  own target. The line anchors the path to the deck's own folder:
  - Tcl: `mpm particles -file [file join [file dirname [info script]] {case1.particles.h5}]`
  - Python: `os.path.join(os.path.dirname(os.path.abspath(__file__)), "case1.particles.h5")`
- **How the fork opens it.** Exactly as given. A relative path resolves from the process's working
  directory, as for any other file argument; the deck's anchoring is what makes it work from
  anywhere. The engine records the absolute path it opened (232c).
- **The `-file` line is not a `command()` row.** It stays out of apeGmsh's `/opensees/commands` and
  `model_hash`. A replay reproduces the particle content under a **new `SET_ID`**, not a
  byte-identical file.
- **Its stem is the deck's.** The model archive (`<stem>.h5`) need not share it; the two files pair
  through `PROVENANCE` (§1.5).
- **Overwrites.** The emitter replaces an existing `<deck stem>.particles.h5` only when that file is
  a Contract-1 file written by an apeGmsh emit: `FORMAT = "ladruno-mpm-particles"` and
  `GENERATOR = "apeGmsh"`. Any other file refuses the emit, and a refused set means no deck is
  written either.
- **Serial only.** The fork refuses `mpm` in OpenSeesSP and OpenSeesMP, so apeGmsh refuses MPM on a
  partitioned or per-rank emit. There is one writer.
- **No live path.** MPM is deck-only.

All of this is apeGmsh ADR 0134 (draft) D1.

### 1.1 Root attributes

| Name | Type | Req. | Default | Rule |
|---|---|---|---|---|
| `FORMAT` | string | yes | — | `== "ladruno-mpm-particles"` (C1-FMT) |
| `FORMAT_VERSION` | int32 | yes | — | `== 1`; every other value is refused (C1-VER) |
| `GENERATOR` | string | yes | — | non-empty, e.g. `"apeGmsh"` (C1-GEN) |
| `SET_ID` | string | yes | — | a uuid4 in canonical 36-character lowercase form, minted when the file is written; the fork copies it into Contract 2 (C1-SETID) |
| `CREATED_UTC` | string | no | absent | ISO 8601, `YYYY-MM-DDThh:mm:ssZ`; informative only |
| `NDM` | int32 | yes | — | `2` or `3` (C1-NDM) |
| `NUM_PARTICLES` | **int64** | yes | — | `N ≥ 1`, equal to the length of every `PARTICLES/*` dataset (C1-N) |
| `PLANE_STRAIN_THICKNESS` | float64 | iff `NDM = 2` | — | finite and `> 0`, in `L`; must be absent when `NDM = 3` (C1-THK) |
| `UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME` | string | yes | — | non-empty unit symbols of one consistent system, e.g. `m`, `kN`, `s` (C1-UNIT) |
| `UNIT_SYSTEM` | string | no | absent | apeGmsh's unit-system id when it has one (`kN_m_s`, `tf_m_s`, `kgf_m_s`, `N_mm_s`, `N_m_s`, `kip_in_s`, `dyne_cm_s`; `apeGmshViewer/src/state/units.ts`); informative |
| `MASS_SOURCE` | string | yes | — | `"mass"` or `"density"`: which one of `PARTICLES/MASS` and `PARTICLES/DENSITY` is present (C1-MASS) |
| `STRESS_SIGN` | string | yes | — | `== "tension-positive"` (C1-SIGN) |
| `VOIGT_ORDER` | string | yes | — | `== "11,22,33,12,23,13"` when `NDM = 3`, `== "11,22,12"` when `NDM = 2` (§0.1) (C1-VOIGT) |

The integer attributes must have the widths listed (C1-DTYPE). `STRESS_SIGN` and `VOIGT_ORDER`
each admit one value. They are in the file so that a producer states the convention it used, and a
reader refuses a file written to another convention instead of reading it wrongly.

### 1.2 `GRID` — the background grid (attributes on the group `/GRID`)

| Name | Type | Shape | Req. | Rule |
|---|---|---|---|---|
| `ORIGIN` | float64 | `[NDM]` | yes | finite; the grid's minimum corner, in `L` (C1-GRID) |
| `CELL_SIZE` | float64 | `[NDM]` | yes | finite and `> 0` on each axis, in `L` (C1-GRID) |
| `NUM_CELLS` | int32 | `[NDM]` | yes | `≥ 1` on each axis (C1-GRID) |
| `BASIS` | string | scalar | yes | `"gimp"` or `"linear"` (C1-BASIS) |

Nodes sit at `ORIGIN + i·CELL_SIZE`, for `i = 0 … NUM_CELLS` on each axis. The file allows a
different cell size per axis. 232b's v1 engine has one cell size `h` for all axes (its §2), so 232c
refuses a set whose `CELL_SIZE` components differ, and a grid that differs from the deck's
`mpm grid`. Both are engine refusals, not file rules.

### 1.3 `BODIES` — names for the body ids (datasets in `/BODIES`)

| Name | Type | Shape | Req. | Rule |
|---|---|---|---|---|
| `ID` | int32 | `[B]` | yes | unique and `≥ 0` (C1-BODY) |
| `NAME` | string | `[B]` | yes | unique and non-empty: the body's label in apeGmsh (C1-BODY) |

`B ≥ 1`. Every `PARTICLES/BODY_ID` value is one of `BODIES/ID` (C1-BODYREF). The names reach
Contract 2 unchanged, so the viewer can name bodies without apeGmsh (apeGmsh ADR 0112 D2c).

### 1.4 `PARTICLES` — one row per particle (datasets in `/PARTICLES`)

| Name | Type | Shape | Units | Req. | Default | Rule |
|---|---|---|---|---|---|---|
| `ID` | int64 | `[N]` | — | yes | — | strictly increasing (so unique), `0 ≤ id ≤ 2⁵³ − 1` (C1-ID) |
| `POSITION` | float64 | `[N, NDM]` | `L` | yes | — | inside the grid box together with its GIMP domain (§0) (C1-POS) |
| `VOLUME` | float64 | `[N]` | `L³` | yes | — | `> 0`; zero or negative is refused (C1-VOL) |
| `MASS` | float64 | `[N]` | `F·T²/L` | iff `MASS_SOURCE = "mass"` | — | `> 0` (C1-MASS) |
| `DENSITY` | float64 | `[N]` | `F·T²/L⁴` | iff `MASS_SOURCE = "density"` | — | `> 0`; the fork computes `mass = DENSITY × VOLUME` once, at load (C1-MASS) |
| `GIMP_HALF_LENGTH` | float64 | `[N, NDM]` | `L` | yes | — | `0 < l ≤ CELL_SIZE/2` on each axis; the engine ignores it when `BASIS = "linear"` (C1-GIMP) |
| `VELOCITY` | float64 | `[N, NDM]` | `L/T` | no | `0` | — |
| `STRESS` | float64 | `[N, 6]` (3D), `[N, 3]` (2D) | `F/L²` | no | `0` | the §0.1 index map (C1-SHAPE) |
| `STRESS33` | float64 | `[N]` | `F/L²` | no; 2D only | `0` | σ₃₃ in plane strain (§0.1); absent in 3D (C1-STRESS33) |
| `BODY_ID` | int32 | `[N]` | — | yes | — | in `BODIES/ID` (C1-BODYREF) |
| `MATERIAL_TAG` | int32 | `[N]` | — | yes | — | `> 0`; a tag **resolved at emit** from the body's material label (charter P4). The fork refuses a tag that names no `nDMaterial` defined before `mpm particles`, or one outside the admitted-material table (232c) (C1-MAT) |
| `FLAGS` | uint32 | `[N]` | — | no | `0` | v1 defines no bit, so every value must be 0. The fork refuses a set bit rather than ignore it (C1-FLAGS) |

- Every dataset has length `N = NUM_PARTICLES` on axis 0 (C1-N), and exactly the dtype and second
  dimension listed (C1-DTYPE, C1-SHAPE). Every float value is finite (C1-FINITE).
- **No other dataset** may sit in `PARTICLES` (C1-UNKNOWN). An unknown per-particle field could
  carry meaning, so a reader refuses it instead of dropping it. A new field needs a
  `FORMAT_VERSION` bump.
- **For lattice seeding**, the GIMP domain fills the particle's volume: `VOLUME = ∏(2·l)` in 3D and
  `t·∏(2·l)` in 2D. This is not a rule, because Gauss-point seeding (later) breaks it. The
  validator reports a mismatch as a warning only.

### 1.5 `PROVENANCE` (group, attributes only)

Producer-written attributes with lowercase names. Readers ignore the ones they do not know.

| Name | Type | Req. | Meaning |
|---|---|---|---|
| `apegmsh_version` | string | when `GENERATOR = "apeGmsh"` | the apeGmsh version that wrote the file |
| `apegmsh_commit` | string | when `GENERATOR = "apeGmsh"` | 40 lowercase hex characters, or `"unknown"` |
| `session_id` | string | when `GENERATOR = "apeGmsh"` | the session's uuid4 (`/meta/session_id` of the paired `<stem>.h5`) |
| `fem_hash`, `model_hash` | string | no | copied byte-equal from the paired `<stem>.h5` `/meta/lineage`, when it has them |

C1-PROV checks the three required attributes.

### 1.6 `SEEDING` (optional group, attributes only)

| Name | Type | Shape | Meaning |
|---|---|---|---|
| `METHOD` | string | scalar | `"lattice"`, `"gauss"` or `"custom"` |
| `PARTICLES_PER_CELL` | int32 | `[NDM]` | `n` per axis, for `"lattice"` |

This group is informative. The uniformity report is the producer's (232d).

### 1.7 Versioning

`FORMAT_VERSION = 1`, an integer the fork owns. It is not one of apeGmsh's semver zone keys (ADRs
0023 and 0113 do not govern it), and the apeGmsh ADR draft names that exception. A reader refuses
every other version. Any new field, any changed meaning, and any new allowed value of an existing
attribute bumps the version. Contract 1 has no minor version, because the fork's reader refuses
unknown fields, so no addition is silently safe.

---

## 2. Contract 2 — the particle group in `.ladruno` (fork → apeGmsh)

The particle group is an **additive** change to schema v1. `INFO/FORMAT_VERSION` stays **1** and
`INFO/FORMAT_MINOR_VERSION` goes from 1 to **2**. Everything in `ladruno_schema_v1.md` holds
unchanged. Minor 2 adds `MODEL/PARTICLES` and `RESULTS/ON_PARTICLES` to a `MODEL_STAGE`.

### 2.1 The version rule: the minor is a feature flag

- **No major bump.** Schema §9 keeps `FORMAT_VERSION` for additive changes, and §3.4 is the
  precedent: CSR connectivity was added without a version bump, because old readers skip groups they
  do not know. apeGmsh treats a `FORMAT_VERSION` bump as a breaking change and opens only `(1,)`
  (`results/schema/_versions.py:73-90`, `_ladruno.py:633-637`,
  `apeGmshViewer/src/reader/ladruno.ts:56`). A major bump would therefore make every new `.ladruno`
  unreadable by today's apeGmsh.
- **The minor is used here as a feature flag.** Old readers skip the particle groups anyway, as
  §3.4's CSR groups are skipped, so safety does not need the minor. It is set so that a reader can
  tell from one attribute read whether to look for particles. Schema §9 reserves the minor for
  revisions a v1 reader still parses. `FORMAT_MINOR_VERSION = 2` means "minor 1 (schema §2a:
  per-stage `KIND`, stage barriers), plus this section". Minors are cumulative. The fork writes
  both attributes (`LadrunoRecorder.cpp:989-994`), and the viewer already reads the minor
  (`ladruno.ts:172`, `:342`).
- **A file that holds a particle group writes minor `≥ 2`** (C2-MINOR). Whether the recorder also
  writes 2 on FE-only files is 232c's choice. The recommendation is no: FE-only output then stays
  byte-identical, and the two fork checks that assert minor `== 1`
  (`Ladruno_scripts/ladruno_recorder_tests/stage_kind_check.py:99`,
  `tests/test_ladruno_recorder_stage_semantics.py:183`) keep holding. If 232c writes 2 everywhere, it
  relaxes both to `>= 1` in the same PR.
- **What an older reader does**, by scope. A particle-only stage is written as the recorder's empty
  stage (§2.3).
  - **A file with at least one FE stage** opens in today's apeGmsh Python reader and ignores the
    particle groups. So does any file opened with `model_h5=`. The fixture
    `particles_mixed_v1m2.ladruno` is that case.
  - **An all-particle file** is **refused** by today's Python reader. `LadrunoReader.fem()` skips
    every `EMPTY_PARTITION` stage, and `Results.from_ladruno` then raises "carries no MODEL group …
    Pass model_h5=". The fixture `particles_v1m2.ladruno` is that case.
  - **Today's TypeScript reader** does not read `EMPTY_PARTITION`. It throws a whole-file layout
    error on any stage `KIND` outside its `KIND_MAP`, so a particle stage's `KIND` is pinned (§2.3).
    With the pin, it opens both kinds of file and shows no particles.

  This is acceptable until apeGmsh's reader slice (M2) reads particles. A reader that wants particles
  checks `FORMAT_MINOR_VERSION >= 2` and then looks for the groups.
- **apeGmsh** keeps `LADRUNO_SUPPORTED_FORMAT_VERSIONS = (1,)` (232d). Its own `.ladruno` writers
  (`_ladruno_writer.py`, the MPCO converter) keep minor 1, since they never write particles.
- **The fork's executable spec**, `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py`:
  - stays `FORMAT_VERSION = 1` and gains `FORMAT_MINOR_VERSION_PARTICLES = 2`;
  - its `validate()` runs this section's checks on every stage that has a particle group
    (`_validate_particles`);
  - its reader gains `particle_stages()` and `particle_output()`.

  `make_synthetic.write_particle_group` is the reference writer, with the C++ encodings.
  `make_synthetic.build(path, particles=True)` covers the case where FE results and particles share
  one file (the base synthetic file writes minor 1). `test_ladruno_harness.py` checks both.

### 2.2 `INFO` in a file with particles

As schema §2, with these values:

- `FORMAT_VERSION = 1` and `FORMAT_MINOR_VERSION = 2`.
- `SPATIAL_DIM = 2 | 3` (C2-NDM), equal to the particle set's `NDM`, and to the nodes' dimension when
  the file also holds FE output.
- `PARTITIONED = 0`. The v1 driver is serial, so a file with particles that is one part of a
  partitioned set is refused (C2-PART).

`STORED_PRECISION` keeps its schema-§2 meaning for `ON_NODES` / `ON_ELEMENTS`. A particle field's
precision is the dtype of its own dataset (§2.6).

### 2.3 The `MODEL_STAGE` that holds particles

- `MODEL_STAGE[<id>]` follows schema §2a for its id and its `STEP` / `TIME` / `KIND` attributes.
- **`KIND` is pinned.** It is `"transient"` for explicit runs (mode F) and implicit-dynamic runs
  (mode A), and `"static"` for quasi-static runs (mode A). It is never empty and never any other
  value: absent, `""`, `"eigen"` or anything else is C2-KIND. Today's apeGmsh TypeScript reader
  throws a whole-file layout error on any other `KIND` (§2.1). **`KIND` and `MODE` agree**
  (C2-MODEKIND):

  | `MODE` | `KIND` |
  |---|---|
  | `explicit` (mode F) | `transient` |
  | `implicit` (mode A) | `transient` |
  | `quasi-static` (mode A) | `static` |

  Mode F has no running `Analysis` object. A writer that shares `LadrunoRecorder`'s stage logic must
  therefore set `KIND` itself (232c), and must not leave the recorder's provisional `"static"` on an
  explicit stage.
- **`PARTICLE_EPOCH`** (int32, required on a particle stage): 0 at the start of a run, incremented
  by each `mpm reset`, and never decreasing in stage order (C2-EPOCH). The id rules of §2.7 apply
  within an epoch.
- **A particle-only stage is LadrunoRecorder's empty stage** (C2-EMPTY). When the run has no FE
  nodes, the writer produces exactly what the recorder writes for an empty stage
  (`LadrunoRecorder.cpp`, `writeModel` at about 1099-1118 and the WP-165 MP-8 branch of
  `writeModelNodes` at about 1189-1227), and then adds the particle groups:
  - the stage attribute `EMPTY_PARTITION = 1`;
  - `MODEL/NODES/ID` as int32 `[0 × 1]` and `MODEL/NODES/COORDINATES` as float64 `[0 × ndm]`;
  - an empty `MODEL/ELEMENTS`;
  - empty `RESULTS/ON_NODES`, `ON_ELEMENTS`, `ON_DOMAIN` and `ON_REGIONS`.

  Today's apeGmsh readers skip such a stage as empty. A stage that has FE nodes must not set
  `EMPTY_PARTITION`.
- **When a new stage opens:** at a user `domainChange` barrier (schema §2a, rule 2), at a grid change
  (§2.4 is per stage), and after `mpm reset` (a new epoch). Removing or adding a particle **never**
  opens a stage; the per-output id lists record it (§2.5).
- `MODEL/PARTICLES` and `RESULTS/ON_PARTICLES` are both present or both absent in a stage (C2-PAIR).
  A stage has at most one particle group in minor 2.

### 2.4 `MODEL_STAGE[<id>]/MODEL/PARTICLES` — per-stage particle metadata

**Attributes**

| Name | Type | Req. | Rule |
|---|---|---|---|
| `SOURCE_SET_ID` | string | when a Contract-1 file was loaded | `SET_ID` of the Contract-1 file the particles came from; with several files, the first one's. It is absent when no file was loaded (C2-SRC) |
| `SOURCE_SET_IDS` | string | no | a JSON array of every loaded `SET_ID`, when there is more than one (C2-SRC) |
| `PATTERN_TAG` | int32 | yes | the host pattern's tag (`LadrunoMPMPattern`); `-1` in mode A, which has no host pattern (C2-PATTERN) |
| `MODE` | string | yes | `"explicit"`, `"implicit"` or `"quasi-static"` (C2-MODE) |
| `PLANE_STRAIN_THICKNESS` | float64 | iff `SPATIAL_DIM = 2` | copied from Contract 1 (C2-THK) |
| `UNIT_LENGTH`, `UNIT_FORCE`, `UNIT_TIME` | string | yes | copied from Contract 1 (C2-UNIT) |
| `UNIT_SYSTEM` | string | no | copied when Contract 1 has it |

**`GRID`** is a subgroup with the same four attributes as Contract 1 §1.2 (`ORIGIN`, `CELL_SIZE`,
`NUM_CELLS`, `BASIS`), constant for the stage (C2-GRID).

**`BODIES`** is a subgroup with `ID` (int32 `[B]`) and `NAME` (string `[B]`), copied from Contract 1
(C2-BODY). In a fork-written file `NAME` is a fixed-length string dataset. Readers accept fixed or
variable length.

**The registry** is a set of datasets directly in `MODEL/PARTICLES`, with one row for every particle
that appears in any output of the stage:

| Name | Type | Shape | Units | Rule |
|---|---|---|---|---|
| `ID` | int64 | `[M]` | — | strictly increasing (C2-REG) |
| `BODY_ID` | int32 | `[M]` | — | in `BODIES/ID` (C2-REG) |
| `MATERIAL_TAG` | int32 | `[M]` | — | `> 0` (C2-REG) |
| `MASS` | float64 | `[M]` | `F·T²/L` | finite and `> 0`; constant per particle (MPM conserves mass), so it is stored once here and not per output (C2-REG) |

The registry is extensible on axis 0, and a particle added during the stage is appended to it. Added
particles take ids above every id used in the epoch (§2.7), so appending keeps `ID` sorted.

### 2.5 `MODEL_STAGE[<id>]/RESULTS/ON_PARTICLES` — the outputs

One output is the state of the live particles after one step. Outputs are **ragged**: their particle
count can change through removal or addition. They are therefore concatenated along one axis with an
offset table, the CSR form schema §3.4 already uses for connectivity, instead of schema §7.0's fixed
`[T × nIds × nComp]`.

**Attribute:** `PARTITION_REDUCTION = "NONE"` (schema §7.1; serial only in minor 2) (C2-RED).

**Index datasets** (directly in `ON_PARTICLES`; `T` outputs, `S = OFFSETS[T]` rows in all):

| Name | Type | Shape | Rule |
|---|---|---|---|
| `STEP` | int32 | `[T]` | the Domain commit tag at the output; strictly increasing (C2-AXES) |
| `TIME` | float64 | `[T]` | the Domain time at the output; finite (C2-AXES). It need **not** increase: `loadConst -time 0` and `setTime` may restart time inside a stage, as on `ON_NODES`. Readers key outputs by index or `STEP`, never by `TIME` |
| `OFFSETS` | int64 | `[T + 1]` | `OFFSETS[0] = 0`, never decreasing, `OFFSETS[T] = S` (C2-OFF) |
| `ID` | int64 | `[S]` | output `k`'s particle ids are `ID[OFFSETS[k] : OFFSETS[k+1]]`, strictly increasing within an output, each one in the registry (C2-ID) |

`T ≥ 1` (C2-AXES). All four index datasets are extensible (`maxshape = None` on axis 0).

**Field groups.** Each recorded field is a subgroup, `ON_PARTICLES/<FIELD>/`:

| Item | Type | Rule |
|---|---|---|
| attr `DISPLAY_NAME` | string | non-empty |
| attr `COMPONENTS` | string | comma-separated labels, none empty or repeated; their count is `nComp` and equals the `DATA` column count (C2-COMP) |
| attr `DIMENSION` | string | absent ≡ `""` (dimensionless) |
| attr `DESCRIPTION` | string | absent ≡ `""` |
| attr `TYPE` | int32 | `0` (Generic, `ResultType` in `Ladruno_Types.h`) |
| attr `DATA_TYPE` | int32 | `0` Scalar, `1` Vectorial, `2` Tensorial (`ResultDataType`) |
| attr `PARTITION_REDUCTION` | string | `"NONE"` |
| dataset `DATA` | float64 or float32 | 2-D `[S, nComp]`; row `r` belongs to particle `ID[r]` (C2-DATA) |

- A missing `DISPLAY_NAME`, `COMPONENTS`, `TYPE`, `DATA_TYPE` or `PARTITION_REDUCTION` is C2-ATTR.
- The field set is fixed for a stage, so every field's `DATA` has exactly `S` rows.
- `POSITION` is required (C2-POS), and every other field is optional.
- Field names match `^[A-Z][A-Z0-9_]*$`, or `^STATE_[A-Za-z0-9_]+$` for material state (C2-NAME).
- Every subgroup of `ON_PARTICLES` is a field. No dataset other than the four index datasets sits
  there (C2-UNKNOWN).

**Append order, and the reader's rule after a crash.** For each output the writer appends, in this
order:

1. every field's `DATA` rows;
2. the `ID` rows;
3. `STEP` and `TIME`;
4. last, `OFFSETS[T+1]`, which commits the output.

The registry rows for added particles are appended before step 1. A closed file satisfies every rule
above. A file left by a crash can hold a partial output beyond the last committed one. A reader then
uses `T' = min(len(STEP), len(TIME), len(OFFSETS) − 1)` outputs, keeps only the outputs whose
`OFFSETS[k+1]` does not exceed `len(ID)` or any field's `len(DATA)`, and ignores rows beyond
`OFFSETS[T']`. The validator checks closed files and requires `OFFSETS[T] = len(ID)`.

**Reading output `k` of field `F`** (Python; the TypeScript reader does the same with hyperslabs):

```python
a, b = OFFSETS[k], OFFSETS[k + 1]
ids  = ON_PARTICLES["ID"][a:b]                     # sorted
vals = ON_PARTICLES[F]["DATA"][a:b]                # [b - a, nComp]
row  = np.searchsorted(PARTICLES["ID"][...], ids)  # registry rows: body, material, mass
```

A field's full history for one particle is a gather over outputs. A dense `[T, M, nComp]` view, with
NaN where a particle is absent, is a reader convenience, not the stored layout.

### 2.6 The field catalogue

The names below are reserved. When one of these fields is written, its `COMPONENTS` must be exactly
as listed (C2-CAT). All values are in the global axes, except `STRAIN` (below). `POSITION` is always
float64 (C2-DTYPE). Every other field is float64 unless the user asks for float32 output, a lossy
opt-in like `-precision f32` on the recorder (C2-DTYPE).

| Field | `COMPONENTS` (ndm = 2) | `COMPONENTS` (ndm = 3) | `DIMENSION` | `DATA_TYPE` | Meaning |
|---|---|---|---|---|---|
| `POSITION` | `x,y` | `x,y,z` | `L` | 1 | current position (required); inside the grid box (§0) (C2-POS) |
| `VELOCITY` | `Vx,Vy` | `Vx,Vy,Vz` | `L/t` | 1 | particle velocity after G2P |
| `ACCELERATION` | `Ax,Ay` | `Ax,Ay,Az` | `L/t^2` | 1 | mode A only |
| `STRESS` | `sigma11,sigma22,sigma12` | `sigma11,sigma22,sigma33,sigma12,sigma23,sigma13` | `F/L^2` | 2 | Cauchy stress, tension positive, in the §0.1 order |
| `STRESS33` | `sigma33` | — (refused in 3D, C2-CAT) | `F/L^2` | 0 | σ₃₃ in plane strain, only from a material that exposes it (none in v1) |
| `STRAIN` | `eps11,eps22,gamma12` | `eps11,eps22,eps33,gamma12,gamma23,gamma13` | absent (`""`) | 2 | reserved: the accumulated unrotated strain feed, 232b's `eps_hat` (material frame, engineering shear) |
| `VOLUME` | `V` | `V` | `L^3` | 0 | current volume, `J·V₀` (2D: area × t) |
| `GIMP_HALF_LENGTH` | `lx,ly` | `lx,ly,lz` | `L` | 1 | current half-lengths, when they evolve |
| `DEFORMATION_GRADIENT` | `F11,F12,F21,F22` | `F11,F12,F13,F21,F22,F23,F31,F32,F33` | absent (`""`) | 2 | row-major; in 2D `F33 = 1` is not stored |
| `STATE_<token>` | the material's response labels | same | as the material says | 0 | a material state variable (below) |

- **Values are finite** in every field except `STATE_*` (C2-FINITE).
- **`STATE_<token>`**: `<token>` is the material response token the user names in `mpm output`
  (232c's guide lists the admitted ones, e.g. an equivalent plastic strain). A particle whose
  material does not answer the token gets `NaN` in every component; this is the only place `NaN` is
  allowed. New tokens need no version bump.
- **A new reserved field** (e.g. a pore pressure in phase 6) is additive: a new optional group with
  its own name. A minor-2 reader that does not know it skips it rather than refusing it. The
  validator refuses only malformed groups.

### 2.7 Particle ids — stable, never reused within an epoch

- A particle keeps the id it had in Contract 1 for the whole epoch, in every stage and every output.
- **Removal** means its id is absent from every later output of the epoch, in this stage and in later
  stages of the same epoch. It never reappears (C2-REUSE), and its registry row stays in the stage
  where it lived. **Removal is forward design:** the v1 engine never removes a particle. It refuses
  a step in which a particle leaves the grid. The rule is in place for the phases that delete
  particles (erosion, outflow).
- **Addition** (inflow, or a second `mpm particles` load) gives each new particle an id above every id
  the epoch has used, appended to the registry. 232c refuses a second Contract-1 file whose ids meet
  any id the epoch has used.
- **`mpm reset`** returns the particles to their loaded state and starts a new epoch: the writer
  opens a new stage with `PARTICLE_EPOCH` one higher. The loaded ids are legal again, and the
  validator restarts its id sets when the epoch changes.

### 2.8 Output cadence and field selection (what the writer does; the command is 232c's)

- **Never every step by default.** 10⁶ particles × about 20 doubles is about 160 MB per output.
  `mpm output` selects outputs every `N` steps or every `Δt_out` of time, and selects the fields.
  `POSITION` is always written.
- **First and last outputs.** At the first `mpm analyze` of a stage, the writer outputs the state
  before the first step, so the viewer has the initial configuration. At the end of every
  `mpm analyze`, it outputs the last step if the cadence skipped it, so the final state is never
  lost.
- **Synchrony.** An output is written from `Domain::commit`'s recorder pass or after it, so positions
  and stresses come from the same step (ADR §5.2 step 8, W11).
- **Flush.** On the recorder's wall-clock cadence (`-flush`, WP-164), and at stage end and close.
- **Chunks.** For `ID` and each `DATA`: about `2²⁰ / (nComp · itemsize)` rows per chunk, with shuffle +
  deflate at the `-compress` level (default 1). For `STEP`, `TIME` and `OFFSETS`: 1024 entries per
  chunk.

### 2.9 Executable spec and where this goes when 232c lands

- `ladruno_format.validate(path)` checks every `C2-…` rule above on each stage with a particle group
  (`_validate_particles`), together with every schema-v1 invariant of the rest of the file.
  `tests/mpm_contracts.py::validate_particle_results(path)` is a thin wrapper that also requires the
  file to have a particle group. apeGmsh SHOULD import `ladruno_format`, as for v1 (schema §9).
- When 232c lands the writer, its format-document PR moves §2 into `ladruno_schema_v1.md` (a new §7.5
  `ON_PARTICLES`, a §2 line for minor 2 and a §9 versioning line) and links this page from there.

---

## 3. Round trip between the contracts

The fixtures are one model seen from both sides. `particles_v1m2.ladruno` is what the fork would write
after loading `particle_set_v1.h5`; `tests/test_mpm_contracts.py` checks the mapping below.
`particles_mixed_v1m2.ladruno` carries the same particle stage after an FE stage.

| Contract 1 | becomes in Contract 2 |
|---|---|
| `SET_ID` | `MODEL/PARTICLES@SOURCE_SET_ID` |
| `NDM`, `PLANE_STRAIN_THICKNESS`, `UNIT_*` | `INFO/SPATIAL_DIM`, `MODEL/PARTICLES@PLANE_STRAIN_THICKNESS`, `@UNIT_*` |
| `GRID@*` | `MODEL/PARTICLES/GRID@*` |
| `BODIES/{ID, NAME}` | `MODEL/PARTICLES/BODIES/{ID, NAME}` |
| `PARTICLES/{ID, BODY_ID, MATERIAL_TAG}` | the registry, same rows |
| `PARTICLES/MASS`, or `DENSITY × VOLUME` | registry `MASS` |
| `PARTICLES/POSITION`, `VOLUME`, `VELOCITY`, `STRESS` | output 0's `POSITION`, `VOLUME`, `VELOCITY`, `STRESS` |
| `PARTICLES/STRESS33` | the material's initial σ₃₃; in Contract 2 only if the material exposes it |

## 4. Fixtures

`python3.12 tests/mpm_contract_fixtures.py` rewrites both files.

- **The model:** 2D plane strain, `t = 1 m`, `kN_m_s`. The grid is 4 × 4 cells of 0.5 m, with GIMP.
  There are two bodies, seeded on a 2 × 2 lattice per cell:
  - `soil`: 8 particles, material tag 2, `ρ = 2.0`, a K0 initial stress. Its σ₃₃ is in `STRESS33`.
  - `block`: 4 particles, material tag 1, `ρ = 7.85`.
- **The fixture tests the file, not the engine.** The v1 engine would refuse the soil's nonzero σ₃₃
  unless tag 2 is elastic (§0.1).
- **Contract 2** is one transient stage: the recorder's empty stage plus the particle groups, with
  `PARTICLE_EPOCH 0` and outputs at steps 0, 10 and 20. Particle 11 is removed after output 1; this
  is valid in the file even though the v1 engine never removes a particle.
- **The mixed fixture**, `particles_mixed_v1m2.ladruno`, holds two stages:
  - `MODEL_STAGE[0]`: `make_synthetic`'s static FE stage, with its nodes, elements and results;
  - `MODEL_STAGE[1]`: the particle stage above.

  This is the case where today's readers open the file and ignore the particles (§2.1).
- **Fields:** `POSITION`, `VELOCITY`, `STRESS`, `STRAIN`, `VOLUME`, `DEFORMATION_GRADIENT` (no
  `DIMENSION`, as the fork writes it) and `STATE_eqPlasticStrain`. The last is `NaN` for the block's
  elastic particles. The token `eqPlasticStrain` is illustrative; 232c's guide names the real ones.

## 5. Open questions

- **Q1 — ownership of Contract 1.** Proposed: the fork owns it, because its reader admits or refuses,
  and apeGmsh changes it through a fork PR. The alternative is apeGmsh owning it as one of its zones,
  which would bring apeGmsh's zone rules (a semver key in `/meta`, its int32 policy) to a file the
  fork parses.
- **Q2 — minor 2 on FE-only files.** Recommended: no (§2.1). This is 232c's call.
- **Q3 — int64 ids.** The ADR fixes int64. The viewer converts them with a safe-integer check, and
  the rule `id ≤ 2⁵³ − 1` keeps that conversion exact. apeGmsh's int32 policy for its ADR-0112 zones
  does not apply to either file, since neither is a `model.h5` zone. That holds only if the apeGmsh
  ADR draft's D1 is accepted.

## 6. Change log

- **2026-11-13, first draft (PR #1077).** Both contracts, the fixtures, the validators and the apeGmsh
  ADR draft. During the draft, the orchestrator ruled that Contract 2 does not bump `FORMAT_VERSION`:
  it is minor 2.
- **2026-11-13, review rulings applied:**
  - plane-strain `STRESS` is `[11, 22, 12]`, with σ₃₃ as the optional `STRESS33`, and the index map is
    stated (§0.1);
  - `STRAIN` is reserved;
  - a particle-only stage is the recorder's empty stage (C2-EMPTY);
  - `PARTICLE_EPOCH` is added, and `mpm reset` starts a new epoch (C2-EPOCH);
  - removal is marked as forward design;
  - `TIME` is only required to be finite;
  - `PATTERN_TAG` has its own code, and `MODE` must agree with `KIND`;
  - C2-COMP is emitted;
  - the grid-box rule, the append order and the crash rule are stated;
  - integer attributes are int32 unless stated;
  - `SOURCE_SET_ID` is optional;
  - fixture strings are `NULLTERM`, and the base synthetic file writes minor 1;
  - every rule code has a refusal case.
- **2026-11-13, apeGmsh review of ADR 0134 (PR #2976, `changes`) applied:**
  - the older-reader claim is scoped (§2.1);
  - a particle stage's `KIND` is pinned, and an absent or empty `KIND` is C2-KIND;
  - a second fixture, `particles_mixed_v1m2.ladruno`, is added;
  - the placement is one rule: beside the emitted deck as `<deck stem>.particles.h5`, deck-only
    (§1);
  - seeding runs in apeGmsh's emitter, from FEMData (ADR 0134 draft D2).
- **2026-11-13, second apeGmsh review of ADR 0134 (`changes`) applied:**
  - the particle file has an overwrite rule of its own (Contract-1 provenance), and a refused set
    means no deck;
  - MPM emits are serial only;
  - the `-file` line is per deck target and not a `command()` row, so a replay gets a new `SET_ID`;
  - both emitters anchor the path to the deck's folder, and the fork opens the path as given (§1).
- **Gate before code.** The owner must ratify apeGmsh ADR 0134 (`232a_apegmsh_adr_draft.md`, through
  apeGmsh's T coordinator and maintainer) before 232c's particle-set reader and `.ladruno` writer,
  and before 232d, start (ADR §9, "Gate before code (232a)").
