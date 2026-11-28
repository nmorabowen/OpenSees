# WP-260 Wave 0 (`260-w0`) — Phase-0 Profiler Scripts

Reproducible Phase-0 measurement harness for the WP-260 Solver & Wave Acceleration program (`#1156`).
Zero `SRC/` modifications; scripts follow `.agents/skills/ladruno-solver-accel/references/wp-briefs.md` §1–§6.

## Dependencies
- Python 3.12 matching the ABI of `dist/bin/opensees.pyd` (`py -3.12`)
- `h5py` and `numpy` installed in the Python 3.12 environment (used in the **parent** process only;
  child processes import `opensees` and write HDF5 via C++ so `opensees.pyd` never co-loads a second HDF5 DLL)
- 5-target release build at `dist/bin/` (`ladrunoBuild()` identical between `OpenSees.exe` and `opensees.pyd`)

## Decks (`decks.py`)
All 5 D1 families (`S` and `L` sizes):
- `d1a_brick_sanisand` (`S`, `L`): `LadrunoBrick -formulation bbar` + `LadrunoSANISAND`, `Pardiso` + `KrylovNewton -maxDim 20 -reuseFactor 3 -reuseIter 8`, `Newmark`
- `d1b_brick_concrete3d` (`S`, `L`): `LadrunoBrick -formulation bbar` + `LadrunoConcrete3D`, `Pardiso` + `KrylovNewton -maxDim 20 -reuseFactor 3 -reuseIter 8`, `DisplacementControl`
- `d2_cpml_wave` (`S`, `L >= 1e5 DOF`): `LadrunoBrick -formulation std` + `LadrunoCPML` + `H5DRM`, `Pardiso` + `Linear -factorOnce`, `Newmark`
- `d3a_explicit_cdl_cpml` (`S`, `L`): `CentralDifferenceLadruno`, `LadrunoBrick -formulation uri -lumped` + `LadrunoCPML -lumped`, `Diagonal`
- `d3b_explicit_bathe_cpml` (`S`, `L`): `ExplicitBathe`, `LadrunoBrick -formulation uri -lumped` + `LadrunoCPML -lumped`, `Diagonal`
- `d3c_explicit_sms_cons` (`S`): `CentralDifferenceSMSConsistent`, `Pardiso`, `LadrunoBrick -formulation bbar`
- `d4_mpm_modeA` (`S`): Mode-A implicit MPM (`mpm grid 232`, `mpm mode implicit`), `Pardiso`, `Newmark`
- `d5_mcft_2cld12` (`S`): `LadrunoForceBeamColumn -monolithic` + `LadrunoNDFiber` + `LadrunoMCFTStrip -cyclicDegradation`, `Pardiso`

## Usage

### 1. Smoke test (small decks + `Linear -factorOnce == Newton` equivalence check)
```bat
py -3.12 Ladruno_scripts\wp260_profiles\run_profiles.py --smoke --verify-d2-linear --export-dir %TEMP%\wp260w0_smoke_json
py -3.12 Ladruno_scripts\wp260_profiles\summarize.py --from-json %TEMP%\wp260w0_smoke_json --out-results %TEMP%\wp260w0_smoke_RESULTS.md
```

### 2. Full ranking sweep (run AFTER `#1153` & `#1154` merge, on a quiet machine)
```bat
py -3.12 Ladruno_scripts\wp260_profiles\run_profiles.py --full --export-dir Ladruno_scripts\wp260_profiles\results\<sha9>
py -3.12 Ladruno_scripts\wp260_profiles\summarize.py --from-json Ladruno_scripts\wp260_profiles\results\<sha9> --out-results Ladruno_scripts\wp260_profiles\RESULTS.md
```

## Outputs
- **Raw `.h5` files (local only, never committed):** `%TEMP%\wp260w0_prof\*.h5`
- **Committed reproducible snapshot:** `Ladruno_scripts/wp260_profiles/results/<sha9>/{manifest.json,trees.json,dump_trees.txt}` (post-#1153/#1154 ranking run)
- **Generated ranking report:** `Ladruno_scripts/wp260_profiles/RESULTS.md` (`summarize.py --from-json` rebuilds it cleanly without `.h5` files or a binary build)
