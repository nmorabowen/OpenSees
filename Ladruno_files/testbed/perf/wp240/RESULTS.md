# WP-240: LadrunoDynamicRelaxation, measured 2026-11-19

Windows desktop, MSVC build of `Ladruno_scripts\build.bat` (all five targets).
Before: `ladruno` f9864cfec. After: `wp/240-dr-single-update-mpi-relax` 22b0573b9.

## 1. One element update per step: same answer, time per step

`dr_single_update_bench.py run <dist/bin> <out.json>` with each build, then
`compare`. Wall time per step is the best of three marches of 399 steps, serial
engine (`opensees.pyd`), default thread settings.

| Model | Values compared | Differ | Before (ms/step) | After (ms/step) | Ratio |
|---|---|---|---|---|---|
| `brick_j2`: 16x16x16 LadrunoBrick, LadrunoJ2, top face driven by `sp`, Transformation handler | 14741 | 0 | 68.621 | 49.146 | 1.40 |
| `stdbrick_j2`: 16x16x8 stdBrick, J2Plasticity, same loading | 7805 | 0 | 52.069 | 35.957 | 1.45 |
| `quad_j2`: 64x64 LadrunoQuad plane strain, LadrunoJ2, edge load | 8452 | 0 | 23.464 | 17.425 | 1.35 |

The values are every nodal displacement after 400 steps plus the residual norm
and the kinetic energy, compared as hex floats. The 18 observables of
`tests/test_wp240_dr_single_update.py` (kinetic and viscous damping) are also
bit-identical between the two builds.

José Abell measured 76 -> 13 ms per step on a plastic 8192-Gauss-point stdBrick
model with a Mohr-Coulomb material (jaabell/OpenSees c8b0f3bc4). The gain here
is smaller. These models use J2 materials, and they run with the default
KE-peak refresh of M* and the stability-margin probe, which the change does not
touch. The share of each part of the step was not measured.

## 2. OpenSeesMP: 1, 2 and 4 processes against serial

Model of `tests/tcl/wp240_dr_mp.tcl`: 8 LadrunoBrick elements of LadrunoJ2 along
x, split by x, loaded to 1.1 x yield, `-noAutoRefresh`, 400 steps. Residual norm
and kinetic energy as each process reports them.

After (22b0573b9), `OpenSees.exe` / `OpenSeesMP.exe`; `opensees.pyd` /
`openseesmp.pyd` print the same digits:

| Run | Step 50 residual | Step 50 KE | Step 400 residual | Step 400 KE |
|---|---|---|---|---|
| serial | 15.023531144663 | 1.96631851444749 | 0.191240200058889 | 7.2740895203218e-07 |
| 1 process | 15.023531144663 | 1.96631851444749 | 0.191240200058889 | 7.2740895203218e-07 |
| 2 processes (both) | 15.0235311446629 | 1.96631851444749 | 0.191240200056459 | 7.2740895202312e-07 |
| 4 processes (all four) | 15.0235311446629 | 1.96631851444749 | 0.191240200056459 | 7.2740895202312e-07 |

Over the 400 steps every value is within 1e-9 of the serial one, relative to
the largest value of its column (the tolerance of `tests/test_wp240_dr_mp.py`).
At step 400 the residuals differ by 1.3e-11 of their value.

Before (f9864cfec), `OpenSeesMP.exe`, same model: each process reports its own
numbers and the run does not relax.

| Run | Process | Step 50 residual | Step 50 KE | Step 400 residual | Step 400 KE |
|---|---|---|---|---|---|
| 2 processes | 0 | 3.72739306177107 | 0.000820203573867786 | 1.69447406929974 | 0.000109341022747634 |
| 2 processes | 1 | 11.1221227552692 | 1.79099459840461 | 2.86024246932345 | 10.8505835945052 |
| 4 processes | 0 | 0.874310584699984 | 5.63505915752928e-05 | 0.252894980510327 | 3.00640511786832e-05 |
| 4 processes | 1 | 1.06119120537298 | 0.000155982201158822 | 0.516345483139077 | 0.237134894348007 |
| 4 processes | 2 | 1.71537352830878 | 0.000285858250479147 | 1.11058534824041 | 0.00167046673705741 |
| 4 processes | 3 | 11.9065151956387 | 1.3833184809323 | 1.32011714430995 | 0.000174312803317877 |

Not run: more than 4 processes, a cluster, Linux, OpenSeesSP.
