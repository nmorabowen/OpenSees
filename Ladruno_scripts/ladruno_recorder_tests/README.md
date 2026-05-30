# MPCO_Ladruno test harness (test-first)

The **executable form** of the `.ladruno` schema spec, written *before* the C++ recorder
so the spec is validated and the implementers have a concrete target. Pure Python
(h5py + numpy) — no OpenSees build needed.

Specs this implements: [`../../Ladruno_implementation/mpco_ladruno_schema_v1.md`](../../Ladruno_implementation/mpco_ladruno_schema_v1.md)
(layout) and the test plan in [`../../Ladruno_implementation/03_mpco_ladruno.md`](../../Ladruno_implementation/03_mpco_ladruno.md)
(§ Testing).

## Files

| File | Role | Test-plan layer |
|---|---|---|
| `ladruno_format.py` | `.ladruno` reader, `validate()` (the executable spec), and canonical normalizers | L2 |
| `ladruno_basis.py` | basis functions (lagrange/bernstein) + global-coord reconstruction `x(ξ)=ΣRᵢ(ξ)Xᵢ` | L3 |
| `make_synthetic.py` | builds a tiny spec-conformant `.ladruno` — self-test fixture **and** reference artifact for the C++ implementers | — |
| `test_ladruno_harness.py` | pytest: validator accepts/rejects, normalizers extract, reconstruction matches | L0/L2/L3 |
| `test_parity.py` | **(Phase 3)** diffs a real `.ladruno` vs legacy `.mpco` to 1e-12 | L1 |

## Run

```powershell
& C:/Users/nmora/venv/opensees_venv/Scripts/python.exe -m pytest -q
# write the reference fixture to inspect by hand:
& C:/Users/nmora/venv/opensees_venv/Scripts/python.exe make_synthetic.py reference.ladruno
```

## How it plugs into the build (the diamond)

- **Now (test-first):** `validate_ladruno` + reader + reconstruction self-test against
  the synthetic fixture. ✅ 10/10 passing.
- **Phase 1 gate:** the C++ skeleton's first (empty) file must pass `validate()`.
- **Phase 3 gate (L1 parity):** `test_parity.py` runs a canonical model with **both**
  `recorder mpco` and `recorder mpcoLadruno`, normalizes each
  (`normalize_nodal`/`normalize_element` here; `STKO_to_python` for the legacy side), and
  asserts value equality to 1e-12. A component-name alias map handles the renames
  (`localForce`→`axial_force`, …).

## Reference fixture model

`make_synthetic.py` builds a 2-D model exercising the load-bearing schema features in one
small file: a unit-square `FourNodeQuad` (stress, standard 2×2 rule) + a 2-node
`DispBeamColumn2d` with a 3-fiber section (`section.fiber.stress`, multiplicity-compressed
COLUMN_MAP) and a `LOCAL_AXES` frame. If the real recorder's output disagrees with this
file's layout, one of them is wrong — that's the point.

## Note on COLUMN_MAP/COMP_NAMES

`COMP_NAMES` is an **attribute** on the `COLUMN_MAP` group (one newline-separated string,
one line per block), not a dataset. The other COLUMN_MAP fields
(`LEVELS/GAUSS_ID/SECTION_TAG/FIBER_ID/NUM_COMP/MULTIPLICITY`) are datasets.
