# ADR-75 P1 gate — Lane-B desktop solver baseline (MEASURED)

Run 2026-07-24 on `nmb` desktop, `dist/bin/opensees.pyd` (2026-06-25), python3.12 `-S`,
**MKL/OMP threads pinned to 1**, median of 7 interleaved rounds, machine verified free of
compilers/builders before starting. Model: Lane-B `15³` LadrunoBrick + LadrunoJ2, 3375 elements,
~11.5k DOF, implicit static, 15 LoadControl steps, Newton, `numberer RCM`, `constraints Plain` —
the same model as the ADR-40b dominance report.

## Results

| Solver | median (s) | vs UmfPack | min–max (s) | spread | tip `ux` | rel err |
|---|---|---|---|---|---|---|
| **UmfPack** | **22.711** | **1.00×** | 22.543 – 28.192 | 24.9% | 2.271760382e+00 | — |
| SparseSYM | 47.648 | 2.10× | 46.525 – 94.646 | 101.0% | 2.271760382e+00 | **0.0** |
| SuperLU | 78.521 | 3.46× | 77.260 – 80.264 | 3.8% | 2.271760382e+00 | **0.0** |
| Mumps | — | — | — | — | **UNAVAILABLE** | `WARNING unknown system type Mumps` |
| Pardiso | — | — | — | — | **UNAVAILABLE** | `WARNING unknown system type Pardiso` |

## Findings

1. **UmfPack is the serial champion — 22.71 s is the bar `system Pardiso` must beat (P1 gate).**
   It is 2.1× faster than SparseSYM and 3.5× faster than SuperLU on this model.

2. **ADR-75's central claim is now empirically confirmed, not just statically inferred.** Both
   `Mumps` and `Pardiso` fail at *runtime* with `WARNING unknown system type` in the serial build —
   exactly as the CMake/`#ifdef _MUMPS` analysis predicted. The desktop regime genuinely has **no
   threaded sparse-direct solver today**.

3. **Correctness cross-check passed exactly.** All three available solvers returned a
   **bit-identical** tip displacement (`rel err = 0.0`), so the timings compare like for like and no
   solver is buying speed with accuracy.

4. **⚠ "Symmetric ⇒ ~2× faster" is contradicted by the only symmetric solver we can measure.**
   `SparseSYM` — a *symmetric* solver on a symmetric tangent — is **2.1× SLOWER** than unsymmetric
   UmfPack. Implementation quality dominates the storage-format advantage. This is a concrete
   datapoint **reinforcing the ADR-75 §12 hedge** ("sparse-symmetric time savings often <2×;
   measure"): PARDISO's `mtype ±2` must be *measured*, never assumed to deliver 2×. (Caveat:
   SparseSYM is a dated implementation, so this bounds the *claim*, not symmetric factorization in
   principle.)

## Measurement quality note

The box's known ±30% wall swing (PHASE0_RECIPE.md) showed up: UmfPack 24.9% and SparseSYM 101%
max-vs-min spread, from one or two contaminated rounds each. **The medians are nonetheless robust** —
timing contamination is strictly *one-sided* (background load can only slow a run, never speed it),
and every median sits within 0.5–2.5% of its own minimum, the signature of a clean central value with
a few high outliers. SuperLU's tight 3.8% spread confirms the harness itself is not the noise source.
Median-of-7 did exactly the job it is in the policy for.

## Reproduce

```bash
cd Ladruno_files/testbed/perf/phase1 && python3.12 -S laneB_solver_bench.py
```

Raw data: `laneB_solver_baseline.json` (all 7 samples per solver retained).
