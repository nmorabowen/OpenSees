# Ladruno test harness (test-first)

The **executable form** of the `.ladruno` schema spec, written *before* the C++ recorder
so the spec is validated and the implementers have a concrete target. Pure Python
(h5py + numpy) — no OpenSees build needed.

Specs this implements: [`../../Ladruno_implementation/ladruno_schema_v1.md`](../../Ladruno_implementation/ladruno_schema_v1.md)
(layout) and the test plan in [`../../Ladruno_implementation/03_ladruno_recorder.md`](../../Ladruno_implementation/03_ladruno_recorder.md)
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
  `recorder mpco` and `recorder ladruno`, normalizes each
  (`normalize_nodal`/`normalize_element` here; `STKO_to_python` for the legacy side), and
  asserts value equality to 1e-12. A component-name alias map handles the renames
  (`localForce`→`axial_force`, …).

## Reference fixture model

`make_synthetic.py` builds a 2-D model exercising the load-bearing schema features in one
small file: a unit-square `FourNodeQuad` (stress, standard 2×2 rule) + a 2-node
`DispBeamColumn2d` with a 3-fiber section (`section.fiber.stress`, multiplicity-compressed
COLUMN_MAP) and a `LOCAL_AXES` frame. If the real recorder's output disagrees with this
file's layout, one of them is wrong — that's the point.

## Live regression battery (`run_regression.bat`)

Once OpenSeesPy is built (`Ladruno_scripts\build.bat OpenSeesPy`), the battery
drives the **real** recorder on many element types/outputs and value-diffs each
`.ladruno` against a frozen `.mpco` (or, for envelope/local-axes/energy/bezier,
against the schema + an analytic/kernel oracle). It points `TMP` at the freshly
built `dist\bin` so the `opensees.pyd` and its MKL DLLs are co-located — that is
why **all** gates run, not just the few that hardcode a separate MKL dir.

```powershell
Ladruno_scripts\ladruno_recorder_tests\run_regression.bat
```

Models run under the **build python** (matches the pyd ABI); checks run under the
**venv python** (h5py + numpy). Each model emits both recorders' files into
`_devrun\`. Gates (all green as of this writing):

| Gate | Element(s) / output | Check kind |
|---|---|---|
| NODAL PARITY | truss frame — displacement, reactionForce | vs mpco 1e-12 |
| ELEMENT PARITY quad/beam | FourNodeQuad stress; dispBeamColumn `section.fiber.stress` | vs mpco 1e-12 |
| MULTI-STAGE | model rebuild across stages | vs mpco 1e-12 |
| STANDARD QUAD | Quad/Tri31/stdBrick/Quad9/Tri10/Hex20 GP geometry | analytic oracle |
| NODAL / ELEMENT ENVELOPE | `-envelope` min/max/absmax | self-consistency |
| LOCAL AXES 2D/3D | beam `LOCAL_AXES` frame quaternions | reconstruction |
| ENERGY BALANCE | `-G energy` ON_DOMAIN + ON_REGIONS | schema + sidecar kernel |
| BEZIER TRI6 | fork BezierTri6 element | schema + mpco |
| **TRUSS/ZEROLENGTH FORCE** | Truss `axialForce`/`basicForces` + zeroLength `force` | vs mpco 1e-12 |
| **FRAME3D localForce** | elasticBeamColumn3d + forceBeamColumn3d (fiber) | vs mpco 1e-12 |
| **SHELL forces** | ShellMITC4 stress resultants | vs mpco 1e-12 |
| **EIGEN modes** | `-N modesOfVibration` eigenvectors (KIND=eigen, `DATA/STEP_k/MODE_k`) | vs mpco 1e-12 |
| **TCL FLAG ORDER** | classic Tcl exe, `-G energy <tag>` in a NON-final position (`flag_order_model.tcl`) | OK marker + all channels present |

(The bold rows are the new element-type coverage. `mp_parallel` is openseesmp-only
and is run separately.)

> **Note — the EIGEN gate caught (and now guards) a real recorder bug.** The
> modal/eigen write path originally wrote **no** eigenvector data: `recordModeChannel`
> called the StreamingSink `begin()` (which makes `DATA` a chunked *dataset*) then
> tried to create a *group* at `DATA/STEP_k/MODE_k` — a name collision. Fixed
> 2026-05-31 (recordModeChannel now owns the modal init: result group + `ID` + `DATA`
> *group*, mirroring the frozen MPCO layout). The gate now diffs modal eigenvectors
> vs mpco to 1e-12. See `LEDGER_quirks.md`.

## Note on COLUMN_MAP/COMP_NAMES

`COMP_NAMES` is an **attribute** on the `COLUMN_MAP` group (one newline-separated string,
one line per block), not a dataset. The other COLUMN_MAP fields
(`LEVELS/GAUSS_ID/SECTION_TAG/FIBER_ID/NUM_COMP/MULTIPLICITY`) are datasets.
