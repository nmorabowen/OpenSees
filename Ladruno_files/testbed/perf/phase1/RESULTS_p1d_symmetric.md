# ADR-75 P1d — symmetric PARDISO (`-matrixType`): MEASURED

**Verdict: symmetric wins on both axes, and it is the memory lever BLR was supposed to be.**
1.35× faster than unsymmetric PARDISO and **−42% peak memory**, exactly — not approximately.

Date: 2026-07-25 · worktree `adr-75-session-handoff-29effa` · same binary for every row.

## What was measured

`system Pardiso -matrixType 0|1|2` (aliases `-symmetric` = 2, `-spd` = 1), where the SOE stores the
**upper triangle only** and PARDISO factorizes with `mtype` 11 / 2 / −2 respectively.

Model: Lane B — 15³ `LadrunoBrick` + `LadrunoJ2` (associated flow ⇒ **genuinely symmetric
tangent**), 11,520 DOF, implicit static, 15 `LoadControl` steps, Newton, RCM.
Harness: `laneB_solver_bench.py` (median of 5 **interleaved** rounds) driven by `sweep_p1d.sh`;
memory by `p1d_memory.py`; Tcl path by `p1d_tcl_smoke.tcl`.

## 1. Time

| MKL threads | solver | median s | vs UmfPack | **vs unsym PARDISO** | tip u_x rel err |
|---|---|---|---|---|---|
| 1 | UmfPack | 18.6426 | 1.00× | — | 0.0 |
| 1 | Pardiso (unsym) | 14.4015 | 1.29× | 1.00× | 0.0 |
| 1 | Pardiso `-matrixType 2` | 11.4641 | 1.63× | **1.26×** | 0.0 |
| 1 | Pardiso `-matrixType 1` | 10.3004 | 1.81× | **1.40×** | 0.0 |
| 4 | UmfPack | 17.6229 | 1.00× | — | 0.0 |
| 4 | Pardiso (unsym) | 12.1545 | 1.45× | 1.00× | 2.0e-16 |
| 4 | Pardiso `-matrixType 2` | **8.9878** | **1.96×** | **1.35×** | 0.0 |
| 4 | Pardiso `-matrixType 1` | 9.0499 | 1.95× | 1.34× | 0.0 |

**Read the "vs unsym PARDISO" column, not "vs UmfPack".** The UmfPack anchor here (17.6 s @4T) is
*not* the locked P1 baseline (22.711 s) — that number came from a different session on a differently
loaded machine, and P1's own PARDISO-vs-UmfPack ratio was 1.71× where this run shows 1.45×. Only the
**interleaved within-run** comparison (same process, same session, A/B/A/B rounds) is trustworthy at
this precision. The cross-session UmfPack numbers should be treated as ±20%.

## 2. Memory — the point of the exercise

`system Pardiso -stats` (new in P1d) dumps PARDISO's own counters: `iparm[14]` peak symbolic,
`iparm[15]` permanent, `iparm[16]` peak numeric, total peak = `max(iparm[14], iparm[15]+iparm[16])`.

| | nnz(A) stored | nnz in factors | **TOTAL PEAK** |
|---|---|---|---|
| unsymmetric (mtype 11) | 818,892 | 8,464,734 | **105.43 MB** |
| symmetric (mtype −2) | 415,206 (**−49.3%**) | 4,459,964 (**−47.3%**) | **61.31 MB (−41.8%)** |
| SPD (mtype 2) | 415,206 | 4,459,964 | **56.42 MB (−46.5%)** |

The stored-nnz figure is exact by construction: `(818892 − 11520)/2 + 11520 = 415206`.

**Contrast with P2b's BLR result — this is the important part.** BLR shrank the *stored factors*
21.8% but moved *peak* memory only 4.6%, because peak is dominated by the active frontal/working
space, not by what is stored. Symmetric half-storage shrinks **both**: the fronts themselves are
half-size, so the fit/no-fit number drops ~42%. And it does so **exactly** — no accuracy tolerance,
no Newton degradation, legal on the byte-identical/oracle paths where BLR is forbidden.

Against the P1c capability wall (UmfPack OOM at 86,490 DOF; PARDISO solved 136,080 DOF), a 42% cut
in peak is the single largest memory lever ADR-75 has produced.

## 3. Correctness

- Tip displacement **bit-identical to UmfPack** (rel err 0.0) for every symmetric configuration at
  every thread count — one 2.0e-16 on the *unsymmetric* row, i.e. within a single ULP. Half-storage
  is exact on a symmetric tangent, as it must be.
- Tcl parity (`p1d_tcl_smoke.tcl`, elastic `stdBrick` cantilever, `OpenSees.exe`): `UmfPack`,
  `Pardiso`, `-matrixType 2`, `-matrixType 1`, `-symmetric`, `-spd` all return
  `ux = 1.035149231914e+00` — identical to 12 digits.
- `matrixType 0` is byte-identical to the pre-P1d code path (the half-storage filter is a
  `!halfStore ||` short-circuit; nothing else in the unsymmetric path changed).

## 4. Why the default stays UNSYMMETRIC

Despite winning on both axes here, `-matrixType 0` remains the default:

1. **Half-storage on an unsymmetric tangent silently solves the wrong system.** Only the `col >= row`
   half of each element matrix is read — no averaging, no detection. The fork has genuinely
   unsymmetric tangents: contact (`LadrunoContact`), non-associated flow, follower loads,
   `LadrunoUP`. Those must not be silently symmetrized by a default.
2. ADR-75 §5 quotes the Kratos precedent as *"default symmetric but expose the override"*. That
   holds for a code whose element library is symmetric by construction; this fork's is not.
3. The win is one explicit token away, and `-stats` now makes it self-evidently worth taking.

**Recommendation to model authors:** if your tangent is symmetric (elastic, J2/associated plasticity,
no contact), use `-matrixType 2`. Prefer **2 over 1**: at 4 threads they are within 0.7% on time,
`2` is only 8% worse on peak memory, and `1` (SPD, Cholesky, no pivoting) **fails outright on an
indefinite tangent** — which any softening or buckling model eventually is. The `-4` zero-pivot error
now says so explicitly. `1`'s single-thread edge (1.40× vs 1.26×) is real but not worth the fragility.

## 5. Reproduce

```bash
cd Ladruno_files/testbed/perf/phase1
./sweep_p1d.sh                                        # time, 1 + 4 threads
OPS_PYD=<worktree>/dist/bin python3.12 -S p1d_memory.py   # peak memory, -stats
<worktree>/dist/bin/OpenSees.exe p1d_tcl_smoke.tcl        # Tcl parity
```

## 6. Two traps banked while measuring this

- **`-matrixType` as a STRING silently costs you the solver.** openseespy
  `system('Pardiso','-matrixType','2')` fails `OPS_GetIntInput`; the factory returned 0, `theSOE`
  stayed null, and OpenSees fell back to the **default `ProfileSPDLinSOE`** with only a warning. On a
  3D solid mesh that is catastrophically slow — a Lane-B bench ran 25 minutes before being killed,
  looking merely "slow" rather than wrong. Fixed twice over: the parse failure now degrades to
  unsymmetric instead of returning 0 (which is what the warning always claimed), and the sweep script
  keeps the full log instead of grepping four lines out of it.
- **MKL 2026.1 renamed its runtime DLLs `.2.dll` → `.3.dll`**, so `build.bat`'s hardcoded staging list
  copied nothing and `import opensees` died with "DLL load failed". Written up as
  `BUILD_GOTCHAS.md §9` (fixed concurrently by #627 — base-name globbing + a stale-DLL purge; P1d only adds the new `mkl_avx10` kernel).
