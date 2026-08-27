---
title: LadrunoSANISAND — Manzari-Dafalias with settable, wired and echoed low-stress constants
project: Ladruno
status: implemented (PR-1); PR-2 and PR-3 outstanding
priority: high
owner: nmora
tags:
  - implementation
  - adr
  - material
  - nd-material
  - soil
  - sanisand
  - manzari-dafalias
  - low-confinement
---

# LadrunoSANISAND (`nDMaterial LadrunoSANISAND`)

> **ADR. PR-1 implemented 2026-08-26 — see the Log; §5 test 5 carries a correction box.** A thin **subclass of `ManzariDafalias`** —
> ND classTags **33019 / 33020 / 33021** (3D and PlaneStrain wrappers) — whose only v1
> difference is that the two low-stress constants `m_Presidual` and `m_Pmin`, currently
> hardcoded in `initialize()`, become **optional constructor arguments, carried on the
> wire, and echoed at construction**. Defaults `0.0` and `1.0e-3 * P_atm`.
> **Vanilla `ManzariDafalias` is not edited.** Motivated by a source audit against
> Dafalias & Manzari (2004) which found three low-stress devices in the code that are
> not in the paper, one of which is currently making **serial and parallel runs of the
> same deck use different materials**. §7 carries two adjacent defects found in the same
> reading that are *not* part of this ADR.

## What

`ManzariDafalias::initialize()` (`SRC/material/nD/UWmaterials/ManzariDafalias.cpp:835-836`)
assigns, unconditionally and with no way to override:

```cpp
m_Pmin      = 1.0e-4 * m_P_atm;   //  0.0101 kPa at P_atm = 101
m_Presidual = 1.0e-2 * m_P_atm;   //  1.01   kPa at P_atm = 101
```

`m_Presidual` is declared in the header as *"small residual pressure (due to cohesion)"*
and is threaded through **thirty** mean-stress evaluations as
`p = one3*GetTrace(stress) + m_Presidual`. It is not a constructor argument, not
settable from the deck, not echoed in the input, and **not on the wire**.

This ADR adds `LadrunoSANISAND : public ManzariDafalias` in which both become inputs.
It changes no existing code path: the only upstream files PR-1 touches are registries.

## Why

### 1. It is an apparent cohesion, and it is large at low confinement

The residual pressure acts as an apparent cohesion of `c = p_r·tanφ ≈ 0.95 kPa`. Its
relative effect on the returned stress ratio is `p_r/p'`, so it is a sub-2 % curiosity
at 50–200 kPa and dominant below about 10 kPa. Measured on a single-element drained
triaxial (build `25a0647`, `stdBrick`/`LadrunoBrick`, IntScheme 1, 400 steps), against
the model's **own** internal identity `η_end = M^b(ψ_end)` — a condition the formulation
satisfies by construction, so any departure is the instrument and no experiment is
needed to see it:

| p₀ (kPa) | 100 | 10 | 3 | **1** | 0.5 | 0.2 |
|---|---|---|---|---|---|---|
| departure from its own bounding surface | +1.0 % | +3.7 % | ~+8 % | **+20.1 %** | ~+30 % | ~+52 % |

The invariant `err × p'_end ≈ p_r` reproduces across the whole sweep, which identifies
the constant without needing any of the reporting project's context: at p₀ = 1 kPa,
`err × p'_end = 1.0136` against `p_r = 1.0100`, i.e. 0.36 %.

Because η saturates at a geometric ceiling in a triaxial test, a "+52 % in stress ratio"
understates it badly in physical terms — the model returns mobilised friction angles up
to **83°** where its own `M^b` says 47°, and the deviator returned barely changes while
the confinement that should produce it falls fiftyfold. A strength independent of
confinement is a cohesion.

**Every leg converges and reports success.** This is a silent-wrong-answer defect, not a
robustness one.

### 2. It reaches the state parameter, not just the stress return

`GetPSI()` is called with the same `p + m_Presidual`. So the residual pressure perturbs
**ψ**, the variable the whole critical-state formulation is organised around, and
through ψ it displaces `M^b`, `M^d` and the dilatancy together. Dafalias & Manzari
(2004) Eq. 9–10 define `M^b = M·exp(-n_b ψ)` and `M^d = M·exp(n_d ψ)` with
`ψ = e - e_c(p)` on the **true** mean stress. Nothing in the paper puts a residual
pressure there.

### 3. Serial and parallel are already different materials — this is the sharp one

`sendSelf` serialises into `Vector data(97)`; `m_Pmin` occupies the last slot, `data(96)`
(`:636`, `:718`). **`m_Presidual` has no slot.**

On the broker path, `FEM_ObjectBrokerAllClasses` constructs via the class-tag/null
constructor, which sets `m_P_atm = 0.0` and *then* calls `initialize()` — so
`m_Presidual = 1.0e-2 * 0.0 = 0.0`. `recvSelf` restores `m_P_atm` from `data(9)` but
**never re-runs `initialize()`**, and there is no slot to restore the residual from.

Consequently **an `OpenSeesMP` worker or a database-restored material runs at
`p_r = 0` while the serial process beside it runs at 1.01 kPa.** Worse, it is
path-dependent: `revertToStart()` (`:465-476`) calls `initialize()` unconditionally
outside an initial-state analysis, by which time `m_P_atm` is real, so the constant
reappears mid-life.

Nothing warns. Any MP-vs-serial comparison on this material at low confinement is
currently invalid.

### 4. Two independent adaptations of this same code already do it the paper's way

Both descend from the UW implementation, so where they differ the difference is a choice
somebody made rather than an ambiguity in the paper.

`SRC/material/nD/UANDESmaterials/SAniSandMS.cpp` (Abell, Liu, Pisanò) — **zero
occurrences of `Presidual`**. It computes `p = one3*GetTrace(stress)` throughout and
clamps at `m_Pmin` with the deviator-preserving rebuild.

`NTUASand02.cpp` (D. N. Gorini, supplied to the reporting project) —

```cpp
// set minimum allowable p and residual p. Why set at this value?
m_Pmin = 1.0e-3 * m_P_atm;
// Presidual=0 because coehsion is not considered
m_Presidual = 0.0;
```

The zero is documented as deliberate and for the correct reason: a cohesionless sand has
no cohesion. This ADR's defaults are exactly these two values.

### 5. Why a new class rather than changing the default in place

Changing `ManzariDafalias`'s own default would alter every existing deck's results
silently, and exposing the constant on the vanilla class would require a wire-format
change (`data(97)` → `data(99)`) that breaks mixed-build parallel and database restore.
A new class leaves every archived result reproducible under the name it was measured
with, and makes "which soil did this number describe" a question about a class name
rather than about a build hash.

## Decisions

| | Decision | |
|---|---|---|
| **D1** | `LadrunoSANISAND : public ManzariDafalias`, ND tags 33019/33020/33021 | Accepted |
| **D2** | v1 difference is **only** `m_Presidual`/`m_Pmin` as optional args, wired and echoed. Defaults `0.0` and `1.0e-3*P_atm` | Accepted |
| **D3** | **Vanilla `ManzariDafalias` is not edited** in PR-1 | Accepted |
| **D4** | **Subclass**, not a full copy, and **not** `virtual` promotions — see §4.1, the RO hazard | Accepted |
| **D5a** | The `D_factor` sigmoid's **shape** — whether it should exist once `m_Presidual = 0`, and whether its constants are right — is **deferred**. Family-wide; §7.2 | Accepted |
| **D5b** | The sigmoid's **unit dependence is fixed in vanilla**, at all four sites, in **PR-2** — *not* in the subclass, and not deferred. It is a **no-op** in the calibrated units. §7.2.1 | Accepted, amended 2026-08-26 |
| **D6** | The `m_e_init` elastic modulus is **not** corrected in v1 — it moves a calibrated quantity (§7.3) | Accepted |
| **D7** | Behaviour-changing items (`mTolR`, the interpolant) go in PR-2, each with a ledger row and its own commit | Accepted |

## How

### 4.1 Route: subclass + flag seams. Why not the other two

`m_Presidual` and `m_Pmin` are **protected data members read at run time** by every base
integrator (header protected block starts at `ManzariDafalias.h:108`). A subclass
therefore never has to override an integrator to change them — it only has to win the
last write. Every method it must control for that is **already virtual** through
`NDMaterial`/`MovableObject`: `revertToStart`, `sendSelf`, `recvSelf`, both `getCopy`
forms, `Print`.

**A full copy** (the `SAniSandMS` route) would duplicate 5229 lines, inherit every
remaining defect invisibly, double the cost of each future fix, and forfeit the one
property the reporting project's experiment needs — being able to state that two runs
differ in exactly one constant and are otherwise the same code.

**Promoting base methods to `virtual` must not be done, and this is the trap.**
`ManzariDafaliasRO` (`SRC/material/nD/UWmaterials/ManzariDafaliasRO.h:93-97`) *shadows*
`initialize()` and two `GetElasticModuli` overloads with identical signatures. Adding
`virtual` to those would convert every one of those shadows into an **override**, so the
base integrators would begin calling Ramberg-Osgood elasticity in every existing RO
deck. A change intended to be inert for other users would silently change their results.

Where a later stage genuinely needs different behaviour inside a **non-virtual** base
method (`integrate()`, `explicit_integrator()`, `ModifiedEuler()` are all non-virtual and
bind statically inside the base TU), the instrument is a **flag seam**: a protected
member on the base, defaulting to the vanilla value, read at the defect site, set only
from the derived constructor. One ledgered line, vanilla bit-identical. PR-2 uses this
for `mTolR`:

```cpp
// ManzariDafalias.cpp:1320
double T = 0.0, dT = 1.0, dT_min = 1e-6, TolE = mHonorTolRInME ? mTolR : 1e-4; // Ladruno
```

### 4.2 New files

| File | Contents |
|---|---|
| `SRC/material/nD/LadrunoSANISAND.{h,cpp}` | class + `OPS_LadrunoSANISAND(void)` parser |
| `SRC/material/nD/LadrunoSANISAND3D.{h,cpp}` | wrapper from `ManzariDafalias3D` (~140 lines) |
| `SRC/material/nD/LadrunoSANISANDPlaneStrain.{h,cpp}` | wrapper from `ManzariDafaliasPlaneStrain` |
| `tests/test_ladruno_sanisand.py` | §5 |

Fork-authored ND materials live flat in `SRC/material/nD/` (the `LadrunoJ2` convention),
not in `UWmaterials/` or `UANDESmaterials/`.

**The wrappers must preserve the sign flip.** `ManzariDafalias3D.cpp:80` is
`mEpsilon = -1.0 * strain_from_element; // -1.0 is for geotechnical sign convention`,
with the matching flip on stress out. The material is compression-positive internally
and tension-positive at the element face. A wrapper that loses this is catastrophic and
silent.

### 4.3 Class tags

```c
#define ND_TAG_LadrunoSANISAND             33019
#define ND_TAG_LadrunoSANISAND3D           33020
#define ND_TAG_LadrunoSANISANDPlaneStrain  33021
```

Highest ND tag currently in use is `ND_TAG_LadrunoRCFiniteStrain = 33018`, so these are
free **in the ND registry**. `LADRUNO_TAG_ComplexEigen = 33019`, `ELE_TAG_LadrunoSolidShell`
and `ELE_TAG_LadrunoCSTPair` are numerically equal but in other registries — deliberately
not collisions, per the convention already documented at `classTags.h:592`, `:643`,
`:945-946`. Add the same note. **Re-grep at merge**: `ND_TAG_LadrunoRCConcrete` carries
`(moved 33014->33015: StagedStrain took 33014 on ladruno)`, so tags have moved mid-flight
before.

### 4.4 The six overrides, and what each prevents

| Override | What it prevents |
|---|---|
| `initialize()` shadow — calls base, then re-applies | The base recomputing `p_r` from `P_atm` |
| `revertToStart()` — replicate the `ops_InitialStateAnalysis` guard at `:468` | `:472` restoring 1.01 kPa **mid-analysis**. Not optional |
| `sendSelf`/`recvSelf` — base `data(97)` unchanged, plus one extra `Vector(4)` | §3 above: worker and serial running different soils |
| `getCopy()` and `getCopy(const char*)` | Base `:392-410` **hardcodes the vanilla wrappers** and would strip the settings at every Gauss point |
| `Print` | A record that cannot state what it ran |

`m_Pmin` takes a sentinel (`< 0` → resolve to `1.0e-3 * P_atm`) because `P_atm` is not
known until the base constructor has run.

**The echo is per construction, not once per process.** A latched message is observable
only by whichever test runs first — see the isolation note in §5.

> **Refined during PR-1, by measurement.** "Per construction" taken literally means *per Gauss
> point*, because `getCopy(const char*)` runs a full constructor for every integration point:
> measured at 513 echo lines for a 64-element mesh, i.e. ~400k lines (~83 MB of stderr) on a
> 50k-element model before the first analysis step. The requirement this clause actually encodes
> is that the message must not be **latched**; echoing once per `nDMaterial` command satisfies it
> in full. PR-1 therefore guards the echo on `getClassTag() == ND_TAG_LadrunoSANISAND` — the
> deck-level object the parser builds announces itself every time, in every process, while its
> Gauss-point copies (which carry the wrapper tags) stay quiet. `Print()` is deliberately
> **not** guarded, so a per-Gauss-point record is still available on demand.

### 4.5 Deck syntax

First 18 positional args and the 5 optional ones are **identical to
`ManzariDafalias`**, so a deck migrates by renaming the command:

```tcl
nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m \
    $h0 $ch $nb $A0 $nd $z_max $cz $Rho <$IntScheme $TanType $JacoType $TolF $TolR> \
    <-Presidual $pr> <-Pmin $pmin> <-honorTolR 0>
```

> `-honorTolR` accepts **only `0`** in PR-1. The seam it would control is one line in
> `ManzariDafalias.cpp:1320`, and PR-1 does not edit vanilla (D3/D7) — so `-honorTolR 1`
> is a **hard parse error**, not a silent no-op. It becomes `0|1` when PR-2 lands the seam.

### 4.6 Registration — five sites, not one

Single-path registration is this fork's most common silent failure; the symptom is a
feature that works in `openseespy` and fails in `OpenSees.exe` with an error that reads
like a broken install. `SAniSandMS` is the working template (correctly wired in both).

1. `SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp` — `extern` near `:138`,
   dispatch block near `:709` (classic `OpenSees.exe`)
2. `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` — decl near `:51`,
   `nDMaterialsMap.insert(...)` near `:173`
3. `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — includes + 3 `case`s beside
   the Manzari block (`~:2539-2554`)
4. `SRC/runtime/runtime/TclPackageClassBroker.cpp` — mirror the Manzari rows
5. `SRC/material/nD/CMakeLists.txt` — sources and headers

### 4.7 Contract items, same PR

`LEDGER_implementations.md` row (tags, files, status, PR); `LEDGER_vanilla_files.md` rows
for each registry file touched, with `// Ladruno` markers in source;
`LEDGER_quirks.md` rows for the two quirks this reading found (the parallel `p_r = 0`
path-dependence of §3; the `revertToStart` reset trap); one line in
`Ladruno_scripts/banner_features.txt` then `python Ladruno_scripts/patch_banner.py`
(never hand-edit the C strings); PR based on `ladruno`.

**Two items this section originally omitted, both found during PR-1. The first is a
blocking CI gate.**

- **`Ladruno_implementation/testbed/manifest.yaml` — one row per class tag, so three.**
  `ci/check_manifest.py` (gate G9) ERRORS (exit 1) on any Ladruno tag in `SRC/classTags.h`
  with no manifest row, keyed on `tag_symbol`. Omitting these leaves `ladruno` HEAD red and
  every later PR inherits the failure (`Ladruno_internal/WORKFLOW_GOTCHAS.md` §3). Per that
  same section the `#define`, the manifest row and the ledger row go in **one commit**.
- **`Ladruno_scripts/stamp_headers.py`** — add the new files' glob to `GLOBS` and run it;
  the stamped LADRUNO header on every new authored file is non-optional
  (`WORKFLOW_GOTCHAS.md` §5). `--check` must exit 0.
- Note also that the banner regen writes **two vanilla files**
  (`SRC/tcl/tclMain.cpp`, `SRC/interpreter/PythonModule.cpp`), so it needs its own
  `LEDGER_vanilla_files.md` row. Easy to miss because the files are generated, not edited.

## 5. Tests — `tests/test_ladruno_sanisand.py`

Model the file on `tests/test_manzari_safety_pack.py`. **The trap that file documents
governs the design:** `mElastFlag` is `static` (`ManzariDafalias.cpp:58`) and every
Manzari-family constructor resets it, so every test must `wipe` and call
`updateMaterialStage` explicitly and no test may rely on stage state surviving a
construction.

Gorini's calibrated parameter set, for a self-contained deck:

```
G0 264.32   nu 0.3129   e_init 0.6944   Mc 1.33090   c 0.71    lambda_c 0.027
e0 0.83     ksi 0.45    P_atm 101.0     m 0.005      h0 1.3    ch 0.968
nb 3.5      A0 0.05     nd 5.75         z_max 12.5   cz 1100.0 Rho 2.0
```

1. **`test_vanilla_equivalence` — the load-bearing gate.** Same strain-driven stack twice:
   `ManzariDafalias` vs `LadrunoSANISAND -Presidual [expr 1.0e-2*$Patm] -Pmin [expr 1.0e-4*$Patm] -honorTolR 0`.
   Committed stresses equal to ~1e-12 relative. This is what licenses the sentence *"the
   two differ in exactly one constant"*. If it does not pass, PR-1 is not finished.
2. **`test_echo_and_defaults(capfd)`** — construct with defaults; assert the echo names
   `p_residual = 0` and `p_min = 1e-3*P_atm`. Plus a **subprocess smoke** running a
   3-line Tcl deck through `OpenSees.exe`, or the classic-Tcl wiring is untested.
3. **`test_revert_to_start_preserves_settings`** — analyse, `reset()`, re-analyse, assert
   identical. A subclass missing the `revertToStart` override fails **only** here.
4. **`test_db_roundtrip_carries_presidual`** — FileDatastore save/restore; assert the
   restored material continues with `p_r = 0`. MP 2-rank smoke where the zone allows.
5. **`test_presidual_is_the_low_p_defect`** — drained shear at `p0 = 0.01*P_atm`; assert
   vanilla departs from its own `M^b` identity by ~+20 % and the `p_r = 0` run by
   materially less. **State as a band, not a number** — the band absorbs step-count and
   scheme sensitivity, which is real (see the measured numbers below).

   > **CORRECTED 2026-08-26, during PR-1, by measurement.** This item originally read
   > *"the `D_factor` sigmoid still suppresses 59 % of dilatancy at that pressure (§7.2),
   > so `p_r = 0` is not expected to give zero."* **That reasoning is wrong and the
   > prediction is wrong.** `D_factor` scales the dilatancy `D`, which governs how the
   > state *travels* toward the bounding surface; it does not perturb the identity
   > `η = M^b(ψ)`, which holds once the state is *there* (`⟨b,n⟩ → 0`, i.e. `K_p → 0`).
   > The instrument `err = η/M^b − 1` measures the `p_r` wedge and nothing else.
   > Measured on the PR-1 battery at `p0 = 1.01 kPa`: vanilla **+18.13 %**
   > (`err·p_end = 0.9903` against `p_r = 1.0100`), `p_r = 0` **+0.075 %**
   > (`err·p_end = 0.0025`). So `p_r = 0` *does* give zero to measurement.
   > `D_factor` shows up instead in the very different `p_end` the two legs reach
   > (3.29 vs 5.46 kPa) — the path, not the identity. Independently adjudicated from
   > Eq. 9–10 and the code's own `GetStateDependent` by a second model.
   >
   > **Also learned here, and sharper than §8 risk 6:** `p_r` is doubling as a numerical
   > regulariser at these pressures. The `p_r = 0` leg *fails* at 400 steps where vanilla
   > survives, and needs 1200; under a strain-controlled isochoric variant at the same
   > `p0` it collapses onto the `m_Pmin` floor and reports `η = 0`. Adopting `p_r = 0`
   > will therefore change **which decks converge**, not only what answers they give.
   > This test is `@pytest.mark.slow` (measured 21.4 s) for that reason, not for cost.

## 6. Backward compatibility

| Change | Existing decks | Disclosure |
|---|---|---|
| PR-1 (new class) | **None** — vanilla untouched | Ledger, banner |
| PR-2 `mTolR` seam | None (flag defaults false; vanilla bit-identical). Migrated decks asking for tight `TolR` on scheme 1 start **getting** it — slower, more accurate | Echo prints the honoured tolerance |
| PR-2 clamp warning | Diagnostics only | PR text |
| PR-2 interpolant fix | **Yes**, but smaller than this row originally claimed and confined to ONE of three sites. MEASURED in PR-2: committed stress moves **2.8e-6 relative** on the G1 deck, and on this ADR's own instrument at `p0 = 1.01 kPa` the `M^b` departure moves **+18.1297 % -> +18.1258 %** (0.004 percentage points). The row's original "~0.012 % on `M^b`" was never reproduced. **`RungeKutta4`/`RungeKutta45` do not move at all**: their loops assign `NextVoidRatio` but read a separate `CurVoidRatio`, so the fix corrects the *returned* void ratio there without touching stress. Own commit for bisection | `LEDGER_vanilla_files` row |
| `-Pmin` default `1e-3` vs vanilla `1e-4` | New class only, but it moves the clamp trigger `p < m_Pmin + m_Presidual` | Echo; A/B protocols must pin `-Pmin` |

## 7. Adjacent defects found in the same reading — NOT part of this ADR

### 7.1 `OPS_SAniSandMSMaterial` drops optional arguments

`SRC/material/nD/UANDESmaterials/SAniSandMS.cpp`, lines ~134 and ~143. `numArgs` is read
before the tag is consumed, so it is `20 + k` for `k` optionals; the code does
`numData = numArgs - 19` (should be `- 20`), consumes `min(numData,3)` ints, then
`numData -= 5` (should be `-= 3`, matching the ints just read).

| k | ints/doubles wanted | as-is | |
|---|---|---|---|
| 0–3 | correct | correct | harmless — the over-read fails and `TclInterpreter::getInt` returns `-1` **before writing**, leaving the default |
| **4** | 3 / 1 | 3 / **0** | **TolF silently dropped** |
| **5** | 3 / 2 | 3 / **1** | **TolR silently dropped** |

`ManzariDafalias` is correct only because its count lines up (tag + 18 doubles = 19) and
it reads all five optionals as doubles in one call, casting the first three. The MS
version looks adapted from it: a 19th parameter was added and the `- 19` never updated.

Two-line fix, but a **behaviour change**: `SAniSandMS`'s integrator *does* honour
`mTolR` (`:1484`, `TolE = mTolR` — it does **not** have the `ModifiedEuler` hardcode),
so decks passing 4–5 optionals are currently running error control at the `1e-7` default
while their input says otherwise. Needs its own regression test.

### 7.2 The `D_factor` dilatancy sigmoid is dimensional — family-wide

Present identically in **all three** implementations (`ManzariDafalias`, `SAniSandMS`,
`NTUASand02`), and in `ManzariDafalias` at **four sites**: `GetStateDependent` (`:4776`)
plus three replicas inside the analytical Jacobian (`:3021`, `:3836`, `:4301`, each with
its derivative `temp1`).

```cpp
if (p < 0.05 * m_P_atm)
    D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
```

The trigger `0.05 * m_P_atm` scales with the unit system. **The sigmoid does not** —
`7.2713` multiplies a raw stress and is a bare literal, so it must carry units of
1/stress. At a true confinement of 1 kPa the factor is **0.410 in kPa, 1.000 in Pa
(never fires), 0.0005 in MPa (total suppression)**. Three different materials from one
input file, and OpenSees has no unit system to catch it.

### 7.2.1 The units fix is a no-op in kPa, and therefore belongs in vanilla (D5b)

Rewrite as `D_factor = 1/(1 + exp(a - b*p/m_P_atm))` with **`b = 7.2713 * 101.0`  (write it as the **product**, not as the literal `734.401`: `734.401/101` is `7.27129703`, a 4.1e-7 relative shift, so the rounded literal silently breaks the no-op that is this fix's whole justification. `(7.2713*101.0)/101.0 == 7.2713` bit-for-bit.)**
(`a` unchanged at 7.6349). Wherever `P_atm = 101.0` kPa this reproduces the shipped
function **exactly** — `b*p/P_atm` *is* `7.2713*p` there — measured at 0.00 % change for
every pressure in 0.2 … 5 kPa. It changes behaviour only for a deck that declared its
atmosphere in some other unit, and for that deck the present behaviour is the defect.

Two consequences follow, and the second is the one that decides scope.

**Disclosure owed.** A deck declaring `P_atm = 101.3` rather than `101.0` shifts by
**1.3 %** in `D_factor` at 1 kPa and about 1 % at 0.5 kPa, because its declared
atmosphere differs from the one the bare constant implicitly assumed. Small, real, and it
goes in the PR text rather than being absorbed.

**It must go in vanilla, not in `LadrunoSANISAND`.** Because the fix is a no-op in
kilopascals, putting it in the subclass would change nothing for any kPa user — which is
every user we know of, including both sibling implementations. The people it helps are
the pascal and megapascal users, and they are reachable only through
`ManzariDafalias` itself. So this is a **vanilla repair in PR-2**, carrying a
`LEDGER_vanilla_files.md` row and a `// Ladruno` marker, and it is *not* a
behaviour-changing item to be gated behind evidence the way §7.3 is.

The same reasoning applies to the sibling copies: `SAniSandMS.cpp:2627` carries the
identical sigmoid with the identical bare constants (one site there, not four), and so
does `NTUASand02`. Fixing all three is one coordinated change or three independent ones;
either way the argument is the same.

**The shape question stays open (D5a).** With `m_Presidual = 0`, is this sigmoid still
needed at all? The two look like overlapping patches for the same problem applied at
different times, and removing one may make the other unnecessary — or may reveal it as
the better device of the two. That is a modelling question for the model's authors, and
it is not settled by non-dimensionalising the constants.

### 7.3 Elastic `G` uses `m_e_init`, not the current void ratio

Dafalias & Manzari (2004) p. 623 give `G = G0·p_at·(2.97-e)²/(1+e)·(p/p_at)^½` with `e`
the **current** void ratio. All three `GetElasticModuli` overloads
(`ManzariDafalias.cpp:4539-4600`, six occurrences) use `m_e_init`. The function is handed
the current void ratio as parameter `en` and does not use it; the **correct line sits
commented out** three lines above (`:4553`); and `b0` in the same file uses the current
`e`. `SAniSandMS` (`:2295`) and `NTUASand02` both use `en` correctly.

Not fixed here because it moves a **calibrated** quantity — the reporting project's
`G0 = 264.32` was fitted against the frozen form, so correcting `G` without revisiting
`G0` would preserve one error by means of another.

## 8. Risks, ranked

1. **An un-audited `initialize()` call site** restoring 1.01 kPa post-construction. Known
   sites are the four constructors (`:205`, `:279`, `:332`, `:384`) and `revertToStart`
   (`:472`). Grep for new ones on the PR checklist; test 3 is the tripwire.
2. **Single-path registration** — five sites (§4.6), two gates.
3. **Static `mElastFlag`, inherited and not fixable by a subclass.** Constructing any
   Manzari-family material resets the stage flag for every instance in the process.
   Existing quirks-ledger rule (explicit `updateMaterialStage` per tag) still governs.
4. **Two variables at low p** if `-Pmin` is not pinned in an A/B comparison.
5. **Tag renumbering by a sibling PR.** Re-grep at merge.
6. **Stall-pattern shift** — `p_r` enters the clamp guard `p < m_Pmin + m_Presidual`, so
   with it at zero the (already non-monotonic) low-p stalls may move. No gate should
   depend on an extreme-low-p leg completing.

## 9. Suggested agent staffing

Offered because the failure modes here are almost all **silent** — the wrong thing
compiles, runs, converges and reports success. The rule applied: *model tier tracks the
cost of a silent error, not the difficulty of the code*, with the release valve that once
a stage is gated by a fail-before/pass-after test, a cheaper tier can write it.

| Stage | Suggested | Effort |
|---|---|---|
| Class + parser (§4.4) | Opus | high |
| Wrappers (§4.2) | Sonnet | medium |
| Registration ×5 (§4.6) | Sonnet | high |
| Tests (§5) | Opus | high |
| Ledgers + banner (§4.7) | Sonnet | low |
| **Adversarial pre-PR review** | **A different model from the authors** | xhigh |

The last row is the one worth keeping. During the audit that produced this ADR, a
second-model pass over source the first model had already read carefully returned three
corrections, one of which **reversed a stated conclusion** — §3 above, which had
initially been recorded as harmless. That was not a difficulty failure; it was a
same-mind failure.

## Log

- 2026-08-26 — ADR written from a source audit of `ManzariDafalias` against Dafalias &
  Manzari (2004) and against `SAniSandMS` and `NTUASand02`. Not yet implemented. The
  reporting project's decision record is its own ADR 66; this document is self-contained
  and does not depend on it.
- 2026-08-26 — **Amended (D5 -> D5a/D5b, new §7.2.1):** the `D_factor` sigmoid's *unit*
  dependence is repaired in **vanilla, in PR-2** (`b = 7.2713*101.0` on `p/m_P_atm`, all
  four sites; an exact no-op at `P_atm = 101`), because the fix is a no-op in kPa and so
  reaches only the Pa/MPa users, who are reachable only through `ManzariDafalias` itself.
  The sigmoid's *shape* question stays open (D5a).
- 2026-08-26 — **PR-1 implemented.** All six overrides ship; five registration sites wired;
  three class tags landed. Battery `tests/test_ladruno_sanisand.py` green (6 passed /
  2 skipped default; 7 / 1 with `--runslow`), `test_manzari_safety_pack.py` still green
  alongside in both file orders. **G1 passes bit-identically**, with a non-vacuity
  companion at 5.58e-2. Corrections applied to this document from the implementation and
  from an adversarial second-model review: the `initialize()` assignment lines
  (`:834-835` -> `:835-836`), `revertToStart` (`:468-474` -> `:465-476`), the RO shadow
  site (`ManzariDafaliasRO.h:49` -> `:93-97`), §4.6's "four sites" (it lists five), §4.7's
  omission of `manifest.yaml` and the header stamp, and — the substantive one — §5 test 5's
  `D_factor` rationale, which was refuted by measurement (see the correction box there).
  §7.1's `SAniSandMS.cpp` line cites `~134`/`~143` were challenged during implementation
  and are **correct as written**; the proposed `:132`/`:141` correction was itself wrong.
  Known gap shipped knowingly: the **PlaneStrain lane (33021) is wired and compiling but
  has no test**; the manifest row says so rather than claiming coverage.
- 2026-08-26 — **PR-1, third pass: two adversarial reviews acted on.** One pass charged to refute
  every claim of fact in the docs, one running a 17-mutation campaign against the battery. **No
  defect was found in the shipped C++ behaviour.** What they found was (a) false claims —
  "MP smoke self-arming" (it skips unconditionally), an incomplete "four sites" fix, an
  irreproducible convergence anecdote — all corrected; and (b) that **CI ran 5 of the 8 tests, and
  the 3 that did not run carried the unique coverage**: `dist/bin` never exists on the runner, so
  the classic-Tcl smoke skipped on every push. Arming it turned Zone-A red on
  `Can't find a usable init.tcl` — `Ladruno_internal/WORKFLOW_GOTCHAS.md` §7 exactly — and the fix
  is the same `TCL_LIBRARY` export the two neighbouring Tcl steps already carried, plus a
  self-discovering retry in the test so the gate is portable. **That gate had never once executed.**
  Three coverage holes then closed, each proven by mutation rather than argued: the PlaneStrain
  lane, `Print`, and the per-Gauss-point echo (§4.4's refinement box). Still open and recorded, not
  hidden: `LadrunoSANISAND3D::getCopy(void)` is uncovered, and `-Pmin` has no *behavioural* gate —
  a real one needs a deck at `p` near `m_Pmin`, i.e. the extreme-low-p regime §8 risk 6 says no
  gate should depend on, so it belongs with **PR-2's clamp diagnostic**, which supplies an
  observable to assert on instead of inferring the clamp from stress output.
