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

> **The wall numbers in §2 were RE-MEASURED on an idle box** after the red-team
> review (finding S4) and the first table was withdrawn. The original run was
> taken while two sibling wp/106 SANISAND jobs were on the box and was captioned
> "speed-ups are therefore lower bounds" — which has the **wrong sign**: on an
> idle box the numbers get *worse*, not better, because contention was suppressing
> the 8-thread oversubscription penalty as much as it was inflating the serial
> baseline. §2 below is the idle-box re-measure. Bit-identity (which is what
> actually gates the WP) was unaffected either way.

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

14 400 elements / 57 600 Gauss points, 12 steps, 3 repeats per thread count,
**idle box**, `MKL_NUM_THREADS=1`, `system Pardiso`, `--h 0.05 --steps 12
--mat elastic`. Re-measured after red-team S4; see the caveat at the top of this
file for why the original table was withdrawn.

| threads | per-step wall (s, min of 3) | mean of 3 | speed-up | full field bit-identical | announced THREADED |
|---|---|---|---|---|---|
| 1 | 0.16014 | 0.17301 | 1.00x | YES | no (serial path, by design) |
| 2 | 0.15550 | 0.16865 | **1.03x** | YES | yes |
| 4 | 0.15474 | 0.17034 | **1.03x** | YES | yes |
| 8 | 0.17050 | 0.19819 | **0.94x — a REGRESSION** | YES | yes |

**12/12 runs bit-identical**, and the oracle is now the FULL FIELD (red-team S6):
one md5 per step over every node's `nodeDisp` and `nodeReaction` at `repr()`
precision plus one element's stress vector. Across all twelve runs there is
exactly **one distinct digest**, `2ddbd45b00ffcc207b374600bbe43195`.

### 2.1 Say the 8-thread result plainly: on this deck, 8 threads is slower than serial

It is **0.94x**, i.e. a 6 % regression on the min-of-3 and 15 % on the mean, and
it is reproducible (the red team measured 0.98x independently on the same deck).
Two effects, both structural:

* **Loop A is only 7.80 % of this deck's step** (§1). Amdahl therefore caps the
  whole-step win at **1.08x** no matter how many threads are thrown at it — so
  the entire available prize is ~8 %, and anything that costs more than that in
  overhead turns the exercise negative.
* **Per-region OpenMP overhead grows with the thread count, and the region is
  entered on every Newton iteration.** Fork/join plus the `schedule(dynamic,8)`
  work-queue contention across 8 threads is a fixed per-iteration tax that is not
  amortised by 14 400 cheap elastic elements; at 24 logical CPUs the last threads
  land on SMT siblings and on cores already carrying the (sequential-MKL) solve.
  The serial pre-work in front of the loop — the element snapshot, the O(nEle)
  virtual allowlist re-audit, and the 3·O(nNode) `Node` trial-state pre-pass — is
  paid once per iteration regardless of thread count (red-team N5) and eats into
  the 8 % as well.

The original table reported **1.11x at 8 threads, which exceeds this deck's own
Amdahl ceiling of 1.08x** and should have been flagged as noise when it was
written. It was not, and that is the lesson worth banking: a speed-up above the
ceiling you computed yourself is a measurement defect, not a result.

The honest summary of §2 is therefore: **the mechanism is correct and costs
nothing at 1 thread, and on the only deck it is allowed to run today it buys
about 3 % at 2–4 threads and loses money at 8.** The deck that would pay
(51.2 % loop A, §1) is refused by §3.

### 2.2 What is NOT measured here, and cannot be — the threaded failure path

Red-team S5. `Domain::ladrunoThreadedUpdate()` carries a deterministic
failure-reporting path: a `critical` that tracks the lowest serial index among
failing elements, an extra diagnostic line, and step-cut parity with the serial
loop. **None of it has ever executed at any thread count, and no deck can make it
execute today.** The only allowlisted material is `ElasticIsotropicPlaneStrain2D`,
whose four `setTrialStrain*` overloads `return 0` unconditionally, and
`LadrunoQuad::update()` returns nothing but the sum of those codes. The path is
therefore **dead code until a second material is allowlisted**, and review item 3
("is failure reporting deterministic under 4 threads?") is answered by
construction and by reading, not by experiment. Recorded rather than papered
over; the cheapest way to close it later is a test-only always-failing
allowlisted `NDMaterial`.

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
| pytest `tests/test_wp107_threaded_update.py` (NEW — the WP-107 warrant, red-team B1) | **18 passed** |
| mutation gate: `LadrunoQuad::shp`/`shpBar` `thread_local` → `static`, rebuilt | pre-existing suite **128 passed** (blind to it); new file **7 failed** |

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
- ~~**Wall times are contended.**~~ WITHDRAWN and re-measured on an idle box
  after red-team S4 — see §2 and the caveat at the top. The "lower bound" framing
  was wrong in sign.
- **The threaded failure path is unreachable** with the current allowlist, so it
  is untested and untestable — §2.2 (red-team S5).
- **No MPI measurement, and none is possible:** WP-107 now REFUSES to thread on
  any `_PARALLEL_INTERPRETERS` / `_PARALLEL_PROCESSING` binary (red-team B2), so
  OpenSeesSP / OpenSeesMP / OpenSeesPyMP run the serial loop unconditionally.
