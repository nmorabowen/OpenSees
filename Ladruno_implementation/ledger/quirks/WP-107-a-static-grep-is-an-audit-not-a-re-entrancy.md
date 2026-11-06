---
wp: WP-107
title: "A static grep is an audit, not a re-entrancy proof — measured on ManzariDafalias (WP-107)"
legacy_seq: 459
---
### A `static` grep is an audit, not a re-entrancy proof — measured on ManzariDafalias (WP-107)

This is the most useful thing WP-107 found, and it cost the WP its headline
payoff, so it is worth stating bluntly.

ADR-75b §5.1 sizes the threading hazard as "~5,600 function-/file-scope
`static Matrix|Vector|ID` declarations across 587 files", which frames
re-entrancy as a *grep problem*. It is not. `ManzariDafalias` under `IntScheme 1`
(ModifiedEuler) has **no function-scope static anywhere on its update call
graph** —

    integrate -> explicit_integrator -> ModifiedEuler
              -> {GetElastoPlasticTangent, Stress_Correction,
                  IntersectionFactor, GetStateDependent, GetStiffness}

every `static Vector`/`Matrix` in that file is in `RungeKutta45`, `NewtonIter`,
`getPStrain` or `sendSelf`/`recvSelf`, all off that path — and no
`Matrix::Solve`/`Invert` on it either. It was allowlisted for WP-107's threaded
`Domain::update()` on exactly that basis.

**It segfaults.** On a 6400-element `LadrunoQuad -bbar` plane-strain deck at 4
threads, reproducibly (4/4), as soon as the PLASTIC branch is exercised in
volume. The elastic branch (gravity stage, `mElastFlag == 0`) is clean and
bit-identical, which is why a mild load path ran 12/12 clean with a byte-exact
curve before a harder one was tried — the most dangerous possible result.

What the experiments rule out, none of which changed the outcome:

| hypothesis | test | result |
|---|---|---|
| the fork's subclass | vanilla `nDMaterial ManzariDafalias` | crashes identically, 4/4 |
| solver interaction | `BandGeneral` vs `Pardiso` | both crash |
| worker-thread stack overflow | `KMP_STACKSIZE=64M` (libiomp5 is the runtime) | no change |
| `opserr` from inside the parallel region | `-Pmin 1e-8` so the clamp never warns | crashes with **zero** warnings emitted |
| the element | same element + `ElasticIsotropicPlaneStrain2D`, 10 000 elements, 8 threads | **6/6 clean, bit-identical** |
| a benign FP-order difference | — | it is a segfault, not a last-bits difference |

Root cause **not located**. The family is therefore refused by
`ManzariDafalias::ladrunoThreadSafeUpdate()` returning false unconditionally.

**A second round excluded three more hypotheses** (full table in
[[75b_ladruno_threaded_assembly_adr]] §14.1), and one of them reframes the
problem: serializing the **entire `theEle->update()`** with `#pragma omp critical`
— so the threads exist and enter the region but never run an update concurrently
— **still faults 0/4**. Serializing the whole of `ManzariDafalias::integrate()`
likewise changes nothing, and `OMP_STACKSIZE`/`KMP_STACKSIZE` at 256 MB changes
nothing. **So this is not a data race between element updates**, which is what
every hypothesis up to that point had assumed. What is left is something about
running this particular update path on an OpenMP worker thread at all — the
elastic path on the same worker thread, same element, same counts, is clean 6/6.

The practical blocker for going further on this box: no `cdb`/`WinDbg`/`procdump`
is installed and the Release build emits **no PDBs**, so a faulting frame could
not be obtained. That is the first thing to fix next time, not another hypothesis.

Three things to carry forward:

1. **A located-but-unfixed hazard is worse than an un-audited one**, because the
   audit creates confidence. WP-107's allowlist defaults to `false` precisely so
   that being wrong is expensive to do rather than free.
2. **A clean threaded run on a mild load path proves nothing.** The bit-identity
   gate has to be run on a deck that exercises the branch you care about; ours
   passed 12/12 with `maxdiff = 0.000e+00` on a configuration that had barely
   entered plasticity.
3. **The next tool is ThreadSanitizer, not another grep.** ADR-75b §7's
   correctness protocol already says so ("the only tool that finds the misses
   `grep` cannot"); this is the measurement that proves the protocol item is
   load-bearing rather than belt-and-braces. TSan needs clang/gcc, so it means
   the Esmeralda/Linux path, not the MSVC desktop one.
