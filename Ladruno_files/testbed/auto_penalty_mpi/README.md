# auto_penalty_mpi — `constraints Auto` sizes its penalty from ONE rank

Regression harness for the ADR-78 follow-up fix: `AutoConstraintHandler`'s
cross-rank `MPI_Allreduce` had been compiled out of every binary since the
handler was contributed.

## The bug

`constraints Auto` picks its penalty from the order of magnitude of the mean
element-diagonal stiffness:

```
PVAL = 10^( round(log10(KAVG)) + oom )        oom defaults to 3
```

`KAVG` is meant to be a **global** quantity — upstream wrote an `MPI_Allreduce`
to make it one. But `AutoConstraintHandler.cpp` compiles into the `OPS_Analysis`
OBJECT library, which is built **once with no parallel define** and folded into
every target. Its `#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)`
block therefore compiled to nothing in every binary, `OpenSeesMP` included.

So under MPI each rank averaged only its **own** elements. A rank holding the
soft part of a model picks a penalty orders of magnitude below the one a rank
holding the stiff part picks — and `PVAL` is precisely the value used for any
constrained node with no local elements, which is the partition-straddling
interface case. Nothing is printed when this happens.

Same trap ADR-78 P1 hit with its first `MPI_Abort`, and ADR-02 calls "the core
blocker" for `OpenSeesCommands.cpp`. Fixed the same way: the reduction moved to
`SRC/analysis/handler/LadrunoAutoPenaltyReduce.cpp`, listed in
`OPS_MPI_PER_TARGET_SOURCES` and compiled per target with that target's defines.
`AutoConstraintHandler.cpp` is now left with **no** `#ifdef _PARALLEL_*` at all,
which is what keeps the trap from silently reopening.

## The model

Deliberately trivial, and built so the three possible answers are four orders of
magnitude apart — no tolerance argument, no rounding ambiguity.

| population | elements | KAVG | PVAL |
|---|---|---|---|
| soft | 100 trusses, `k = 1e2` | `1e2` | **`1.000000e+05`** |
| stiff | 1 truss, `k = 1e6` | `1e6` | **`1.000000e+09`** |
| whole | both | `1e4` | **`1.000000e+07`** |

The whole-model `KAVG` lands on exactly `1e4`: `(100*2*1e2 + 1*2*1e6)/202 = 1e4`.
Each truss contributes exactly 2 diagonal entries above the handler's
`1e-12*kmax` tolerance, hence the counts 200 and 2.

Each rank also declares an **elementless** `equalDOF` pair. That is the probe:
the handler accumulates a per-node penalty from the elements attached to a node,
and falls back to `m_global_penalty` when a constrained node has none — exactly
what a replicated interface node whose elements live on the other rank looks
like. So `PVAL` is not merely printed here, it is the value actually used.

## Reference values

| run | expected |
|---|---|
| `serial.tcl soft` | one `PVAL = 1.000000e+05` |
| `serial.tcl stiff` | one `PVAL = 1.000000e+09` |
| `serial.tcl whole` | one `PVAL = 1.000000e+07` |
| `mp.tcl` on 2 ranks | **two** `PVAL = 1.000000e+07` |

**Pre-fix signature — measured, not predicted:** the 2-rank run prints
`1.000000e+09` and `1.000000e+05`, one per rank, each sized from its own
elements. Obtained by mutating `ladrunoAutoPenaltyReduce()` to `return 1`
immediately (which is exactly what the dead `#ifdef` amounted to), rebuilding
`OpenSeesMP` alone, and re-running: the three serial rows stayed green and only
the 2-rank row flipped. `run.py` names that case explicitly when it sees it, so
a regression reads as itself rather than as a generic number mismatch.

Re-run that mutation if you touch this path. A green harness proves nothing
until the broken version fails — the same lesson ADR-78 P1 learned when its
first `MPI_Abort` compiled clean and did nothing.

The checker asserts that **both** `PVAL` lines equal the whole-model value
rather than attributing a line to a rank: the two ranks' verbose reports
interleave on stdout in an arbitrary order, and "both are right" is both
stronger and order-independent.

The serial rows matter as much as the parallel one. They are what make the
2-rank result a *fix* rather than a coincidence: `1e5` and `1e9` are shown to be
the values the two partitions produce alone, so a 2-rank run printing `1e7`
twice can only have reduced across ranks.

## Running

```
python run.py <tree>\build\build\Release
```

It runs all four cases and prints a PASS/FAIL table. It names the binaries by
full path on purpose — a bare `OpenSees` on PATH has resolved to another
session's build on this machine more than once.

## Trap

A hung MPI job leaves orphaned ranks alive that hold the binary open, and the
*next* link fails with `LNK1104: cannot open file 'OpenSeesMP.exe'`. If a link
fails that way, look for stray `OpenSeesMP.exe` / `hydra_pmi_proxy` processes
before touching the build system.
