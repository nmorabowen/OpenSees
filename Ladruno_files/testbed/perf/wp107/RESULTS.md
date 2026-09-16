# WP-107 results — threaded `Domain::update()` (ADR-75b loop A), desktop-scoped

Binary: `wp/107-omp-state-determination`, built by `Ladruno_scripts\build.bat`
(`LADRUNO_OPENMP=ON`). Box: 24 logical CPUs, Windows 11. **`MKL_NUM_THREADS=1`
pinned on every run** (ADR-75b §3 P-6: otherwise the solver's own ~1 ULP jitter
masquerades as an assembly race).

Deck: `wp107_strip_bench.py` — plane-strain `LadrunoQuad -formulation bbar`,
`constraints Transformation`, `numberer RCM`, `system Pardiso`, `algorithm Newton`,
gravity then prescribed settlement. Default load case is the **laterally confined**
one (`--load oedometer`); see the note in the bench about why the footing case is
not an instrument.

> **Caveat on every wall number below:** the box was concurrently running two
> sibling wp/106 SANISAND jobs. Speed-ups are therefore lower bounds, and the
> bit-identity results (which are what actually gate the WP) are unaffected.

---

## 1. The gate — loop A's fraction of step (serial, deep profiler, 1 thread)

ADR-75b's gate is **a single loop ≥ 40 % of step**, evaluated per loop, never on
the aggregate kernel fraction (§12.4's own correction).

| deck | elements / Gauss pts | nDOF | **loop A** | loop B | loop C | solve | gate |
|---|---|---|---|---|---|---|---|
| `LadrunoQuad -bbar` + **LadrunoSANISAND** (IntScheme 1, TanType 0) | 6 400 / 25 600 | 12 798 | **51.23 %** | 15.80 % | 8.02 % | 17.06 % | **PASS** |
| `LadrunoQuad -bbar` + **ElasticIsotropic** | 14 400 / 57 600 | 28 798 | **7.80 %** | 26.07 % | 12.38 % | 47.46 % | **FAIL** |

Compare G-L3 (ADR-75b §13, `stdBrick` + `J2Plasticity` under MUMPS at np=16):
loop A was **0.26 %** at 540 675 DOF, a ~42x gate failure.

**So §1 of the WP note is confirmed: the desktop + expensive-kernel regime is a
genuinely different regime, and loop A passes the gate there.** Kernel vs scatter
inside loop A on the SANISAND deck: 50.23 % kernel / 1.96 % scatter — i.e. almost
all of it is threadable work, unlike loop B (38 % of it is scatter).

Per-element cost, for the fork/join regression check: `LadrunoQuad` classTag 33007
at **7.70 µs/ele** in loop A on the SANISAND deck, against a barrier cost ADR-75b
§7 measured at ~19.5 µs for a whole region. No regression risk at these sizes.

**And this is the WP's central tension:** the only material that makes loop A
51 % of step is the one §3 shows is not thread-safe.

---

## 2. Bit-identity and speed-up — `ElasticIsotropicPlaneStrain2D` (allowlisted)

14 400 elements / 57 600 Gauss points, 12 steps, 3 repeats per thread count.
Curves compared as an exact string compare of `repr()`'d doubles.

| threads | per-step wall (s, min of 3) | speed-up | curve bit-identical | max abs diff | refused? |
|---|---|---|---|---|---|
| 1 | 0.16557 | 1.00x | YES | 0.000e+00 | no |
| 2 | 0.16518 | 1.00x | YES | 0.000e+00 | no |
| 4 | 0.15146 | 1.09x | YES | 0.000e+00 | no |
| 8 | 0.14872 | 1.11x | YES | 0.000e+00 | no |

**12/12 runs bit-identical, `maxdiff` exactly 0.000e+00.** The threaded loop
demonstrably ran (the one-time `ladrunoThreads: element update loop THREADED on N
threads` line is present, and absent at 1 thread).

The speed-up is ~1.1x and that is the **correct** answer: loop A is 7.80 % of this
deck, so Amdahl caps it at **1.08x (8T)**. Threading buys nothing here because
there is nothing to buy. Reported to make the point that the mechanism is sound
and the *deck* decides whether it pays.

---

## 3. `LadrunoSANISAND` — measured, then WITHDRAWN

### 3.1 What was measured before the defect was found

6 400 elements / 25 600 Gauss points, IntScheme 1, TanType 0, `ds = 0.002`,
15 steps, 3 repeats:

| threads | per-step wall (s, min of 3) | speed-up | curve bit-identical | max abs diff |
|---|---|---|---|---|
| 1 | 0.27401 | 1.00x | YES | 0.000e+00 |
| 2 | 0.24163 | 1.13x | YES | 0.000e+00 |
| 4 | 0.20485 | 1.34x | YES | 0.000e+00 |
| 8 | 0.20280 | 1.35x | YES | 0.000e+00 |

Against the Amdahl ceiling from loop A = 51.23 % (2T 1.34x, 4T 1.62x, 8T 1.81x),
that is 84 % / 83 % / 75 % of the ceiling. Backing the realized loop-A speed-up out
of the totals gives **1.30x at 2T, 1.97x at 4T, 2.03x at 8T** — the element loop
itself saturates near **2x** and stops, so beyond ~4 threads the binding constraint
is loop A's own parallel efficiency, not Amdahl. (Two sibling jobs were on the box;
treat 2x as a floor.)

**These numbers are withdrawn as a shippable result.** They were produced on a
configuration the WP now refuses, and they were produced on a load path that had
barely entered plasticity — see §3.2.

### 3.2 Why it is refused: the audit was clean and it segfaults anyway

`ManzariDafalias` under IntScheme 1 (ModifiedEuler) has **no function-scope static
anywhere on its update call graph**, and no `Matrix::Solve`/`Invert` on it. It was
allowlisted on that basis and returned the 12/12 bit-identical table above.

Then a harder load path (`ds = 0.02`, or TanType 2 at `ds = 0.002` — both simply
put more Gauss points into the plastic branch) **segfaults**: 4/4 at 4 threads.

| hypothesis | test | result |
|---|---|---|
| the fork's subclass | vanilla `nDMaterial ManzariDafalias` (`--mat manzari`) | crashes identically, 4/4 |
| solver interaction | `BandGeneral` vs `Pardiso` | both crash |
| worker-thread stack overflow | `KMP_STACKSIZE=64M` (libiomp5md **is** our runtime) | no change, 3/4 still crash |
| `opserr` from inside the parallel region | `-Pmin 1e-8` so the low-`p` clamp never warns | crashes with **zero** warnings emitted |
| the element | same element + `ElasticIsotropicPlaneStrain2D`, 10 000 elements, 8 threads | **6/6 clean, bit-identical** |
| the elastic branch of the same material | gravity stage (`mElastFlag == 0`) | clean |

Root cause **not located** ⇒ `ManzariDafalias::ladrunoThreadSafeUpdate()` returns
`false` unconditionally, and a threaded run on such a deck falls back to the serial
loop **loudly**, naming the element tag:

```
WARNING ladrunoThreads: element 1 (classTag 33007) is not on the WP-107
thread-safe allowlist -- running the element loop SERIAL.
```

Verified after the change: the previously-crashing configuration is 4/4 clean at
`--threads 4`, with identical curves and no `THREADED` line.

**The transferable lesson**, banked in `LEDGER_quirks.md`: a `static` grep is an
audit, not a re-entrancy proof; and a clean threaded run on a mild load path proves
nothing. ADR-75b §7's ThreadSanitizer protocol item is the next tool, and it needs
clang/gcc — Esmeralda, not this desktop.

---

## 4. Serial-path regression

| check | result |
|---|---|
| `LADRUNO_OPENMP=OFF` build (`set LADRUNO_NO_OPENMP=1`) compiles | yes |
| OFF build, `ladrunoThreads 4` | warns "built WITHOUT LADRUNO_OPENMP … stays SERIAL", never threads |
| OFF vs ON at 1 thread, elastic deck (14 400 ele) | **byte-identical** |
| OFF vs ON at 1 thread, SANISAND deck (6 400 ele) | **byte-identical** |
| `ladrunoThreads` query / set / clamp-at-1 | correct; `0` clamps to 1 with a warning |
| pytest `test_ladruno_sanisand{,_integrator,_responses,_responsetype,_implex}`, `test_ladrunoQuad_eas`, `test_ladrunoquad_finite`, `test_fourNodeQuad_T0`, `test_quad_tri_rho_db_restart`, `test_ladrunoQuad_sanisand_implex_commit_refusal`, `test_ladrunoTie_mortar_quad8` | **128 passed, 2 skipped, 2 xfailed** |

## 5. Not measured / not run

- **No "before" baseline binary was built.** These suites were run on the WP
  binary only; the serial-path evidence is the OFF-vs-ON byte-identity above plus
  the suites being green, not an A/B against a separately built `ladruno` tip.
- **Full Zone-A was not run** (deliberately — it is a 27-minute sweep).
- **ThreadSanitizer was not run.** It is the tool that would find §3.2's defect and
  it needs clang/gcc.
- **No cluster / MPI measurement.** `PartitionedDomain` and `Subdomain` refuse by
  construction; ADR-75b §11 q6 (hybrid MPI+threads) is untouched.
- **Loops B/C were not threaded**, so there is no "update+formTangent" table and no
  measured serialised-assembly cap. The blocker is structural (`FE_Element::theTangent`
  is a class-wide pool) and is argued in the WP note §4, not measured here.
- **Wall times are contended.** Two sibling wp/106 jobs ran throughout.
