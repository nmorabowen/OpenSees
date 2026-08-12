# asdconcrete_tag_cache — ASDConcrete3D silently reuses the FIRST hardening law registered under a material tag

**Status: REPORTED, NOT FIXED (decision 2026-08-12).** `ASDConcrete3DMaterial.cpp`
is upstream and the fix was deliberately deferred rather than taken here. This
directory is the evidence, so nobody has to re-derive it.

**Live consequence:** `tests/test_ladrunoBrick_asdconcrete_bend.py::test_notched_bend_mesh_objectivity`
FAILS in a full-suite run and PASSES in isolation. That is expected until the
material is fixed — **do not "repair" that test**; it is correctly detecting this
bug. See the §Why the test fails section below.

## The defect

`ASDConcrete3DMaterial::HardeningLawStorage` is a process-global singleton keyed
**only by material tag**, and `store()` writes only into an empty slot
(`ASDConcrete3DMaterial.cpp:1263-1276`):

```cpp
auto& item = m_tension[hl.tag()];
if (item == nullptr)
    item = std::make_shared<HardeningLaw>(hl);   // first one wins, forever
```

It is a function-local `static`, so **`ops.wipe()` cannot clear it** and it lives
for the whole process.

`HardeningLaw::regularize()` reads that cache on every auto-regularized material
(`:947-957`):

```cpp
double lch_scale = lch > 0.0 ? lch_ref / lch : 0.0;
if (!m_fracture_energy_is_bounded || lch_scale <= 0.0 || lch_scale == 1.0)
    return;                 // <-- guard that hides the bug when lch == lch_ref
deRegularize();             // <-- recover(m_tag, m_type): overwrites THIS law
                            //     with whatever was stored first for this tag
```

So the **first** `ASDConcrete3D` defined with tag *N* owns that tag's tension and
compression backbones for the life of the process. Every later material with the
same tag silently inherits them.

### Preconditions (all three required)

1. `-autoRegularization` is on (otherwise `regularize()` is never called), **and**
2. `lch != lch_ref` (otherwise the `lch_scale == 1.0` guard returns before
   `deRegularize()`), **and**
3. some earlier `ASDConcrete3D` in the same process used the same tag with a
   different backbone.

Miss any one and the bug is invisible — which is exactly why a first attempt at a
reproducer, using a unit cube with `lch_ref = 1.0`, came back clean.

### Why this matters outside the test suite

A softening-curve parameter sweep in one Python session — the natural way to write
one, reusing tag 1 and calling `ops.wipe()` between cases — returns the **first**
curve for every case. Nothing warns. The results look entirely plausible.

## The measurement

`python repro_tag_cache.py` (needs `LADRUNO_OPENSEES_BIN` pinned to a build; no
gmsh, no pytest, one 8-node cube):

```
EXPERIMENT  same tag 1, ductile first then brittle
  tag1 ductile W = 0.090456
  tag1 brittle W = 0.090456
  identical      = True     <-- BUG
CONTROL     the two backbones really are different (fresh tags)
  tag21 ductile W = 0.090456
  tag22 brittle W = 0.030812
  differ          = True
```

The control is load-bearing: it proves the two backbones genuinely differ by
**2.9x**, so the identical pair in the experiment can only be the cache.

Two instrument bugs were caught before the numbers were believed, both by
checking a control rather than a verdict:

* **Peak load cannot see this.** Peak reaction is `FT * area = 3.0` for BOTH
  backbones, so the first version of this script read `2.999999` four times and
  its "differ" check passed on 1e-9 float noise. The softening *branch* is what
  differs, so the script integrates dissipated work instead.
* **`side=2.0` is required, not cosmetic.** With a unit cube and `lch_ref=1.0`
  the `lch_scale == 1.0` guard returns before `deRegularize()`, the cache is
  never read, and the bug does not reproduce.

## Why the notched-bend test fails

`tests/test_ladrunoBrick_asdconcrete.py` (single-element, runs FIRST alphabetically)
registers tag **1** with a brittle backbone — `Te=[ET0, 1.0e-3, 8.0e-3]`,
`Ts=[FT, 1.0, 0.05]`.

`tests/test_ladrunoBrick_asdconcrete_bend.py` then asks for tag **1** with a
deliberately **ductile** backbone — `Te=[ET0, 2.0e-3, 1.5e-2]`,
`Ts=[FT, 1.8, 0.15]` — whose own comment says it is shallow-softening "so
load-point displacement control + step-cutting traverses the full crack". It
silently gets the brittle one, the softening is far sharper than the deck was
designed for, and the coarse solve stalls:

```
AssertionError: coarse: only reached d=0.034 (<0.05)
```

Note the failure never reaches the mesh-objectivity assertions at all — it dies
on the reachability precondition, which reads like a solver problem rather than a
material problem. That is what made it expensive to attribute.

Confirmed both ways:

| bend deck material tag | result |
|---|---|
| 1 (collides with the single-element suite) | 1 failed, `only reached d=0.034` |
| 91 (unique) | 4 passed |

Which of the three tests in the single-element file leaks:

| test | registers ASDConcrete3D tag 1 | leaks |
|---|---|---|
| `test_char_length_is_cbrt_volume` | no — `ElasticIsotropic` only | clean |
| `test_lch_handshake_mesh_objectivity` | yes | **leaks** |
| `test_tier_a_damage_scaled_kstab` | yes | **leaks** |

## If someone fixes it later

The two candidate fixes, both in `HardeningLawStorage`:

1. **Replace on mismatch** — `store()` overwrites when the incoming law differs
   from the cached one. Smallest change, keeps the sharing optimisation the cache
   exists for (many Gauss-point copies of one material tag share one law object).
2. **Key by tag + law content** — different backbones under one tag coexist.
   Cleaner, slightly larger.

Either is a `LEDGER_vanilla_files` row (`ASDConcrete3DMaterial.cpp` is upstream),
and `repro_tag_cache.py` is the gate: `identical` must flip to `False` while the
control keeps reading `True`.
