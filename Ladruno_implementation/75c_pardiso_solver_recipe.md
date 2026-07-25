# `system Pardiso` — implementation & usage guide (for the apeGmsh team)

**Status:** usage guide for a **shipped** feature. Everything here is on `ladruno`
today and measured, not planned.

**Audience:** whoever wires solver selection into apeGmsh's OpenSees export /
live backend. You do not need to read the ADR to use this — but if you want the
decision record it is [[75_ladruno_sparse_direct_strategy_adr]], with the
measurements in `Ladruno_files/testbed/perf/phase1/RESULTS_p1*.md`.

**One-line summary:** on a desktop, `system Pardiso` replaces `system UmfPack`
and is **1.7×–3.4× faster** (the gap grows with model size) and solves models
UmfPack cannot fit in memory at all.

---

## TL;DR — what to emit

```python
ops.system('Pardiso')                     # <- correct default for almost everything
```

Threading comes from the **environment**, not from a solver option:

```bash
set MKL_NUM_THREADS=4
```

That is the whole recommendation for a first integration. The options below are
opt-in tuning; the default is the safe one on purpose.

---

## 1. When it applies

| Target | Verb |
|---|---|
| `OpenSees.exe`, `openseespy` (serial/desktop) | **`system Pardiso`** |
| `OpenSeesSP` / `OpenSeesMP` (MPI) | `system Mumps` — Pardiso refuses with an explicit error |

PARDISO is shared-memory threaded (OpenMP). It is wired into the **serial**
targets only; the parallel targets link the sequential MKL layer, so there is no
shared-memory win to be had there. Asking for `Pardiso` under `OpenSeesMP`
now fails loudly rather than silently continuing on whatever solver was set
before — but don't rely on that: branch on the target when you generate the deck.

A build without MKL simply does not register the verb; fall back to `UmfPack`.

---

## 2. Options

| Option | Value | Default | What it does |
|---|---|---|---|
| `-matrixType` | `0` \| `1` \| `2` | `0` | `0` unsymmetric (full CSR) · `1` SPD (`LL^T`) · `2` symmetric indefinite (`LDL^T`) |
| `-symmetric` | *(flag)* | — | alias for `-matrixType 2` |
| `-spd` | *(flag)* | — | alias for `-matrixType 1` |
| `-krylov` | int `L` | off | reuse the previous factorization as a **preconditioner** (CGS); stops at `10^-L` |
| `-stats` | *(flag)* | off | dump PARDISO's own peak-memory counters + CGS decisions, once per pattern |

Same spelling in Tcl and Python. Numbering deliberately matches
`system Mumps -matrixType`.

---

## 3. Which configuration for which model

```
Is the tangent symmetric?
├─ NO  (contact, non-associated flow, follower loads,
│       corotational transforms, LadrunoUP)   ->  system Pardiso
│                                                 (+ -krylov 6 if ≳25k DOF)
└─ YES (associated-flow plasticity, elasticity,
        standard solid/shell models)
   ├─ and it stays positive definite         ->  -matrixType 1  (+ -krylov 6)
   └─ softening / buckling / limit points    ->  -matrixType 2   (NO -krylov)
```

**If you are unsure whether the tangent is symmetric, use the default.** Being
wrong in the symmetric direction is a *silent wrong answer*; being wrong in the
unsymmetric direction only costs speed. See trap 1.

---

## 4. What it buys (measured, Lane B = 3D solid + J2, 4 threads)

| Change | Effect |
|---|---|
| `UmfPack` → `Pardiso` | **1.71×** at 11.5k DOF, **3.40×** at 51k — the gap *compounds with size* |
| … capability | UmfPack **ran out of memory at 86k DOF**; PARDISO solved it, and 136k DOF in 69 s |
| `+ -matrixType 1\|2` | a further **~1.25×**, and **−42% peak memory** — exact, not approximate |
| `+ -krylov 6` | a further **1.51×** at 51k DOF under full Newton |
| `-matrixType 1 -krylov 6` | **1.57×** vs plain `Pardiso` (the two levers compose sublinearly) |

Two things worth internalising:

- **Never generalise a solver verdict from a small model.** Every ratio above
  grows with N. A 5k-DOF trial will make PARDISO look barely worth it, and a
  50k-DOF production model will not.
- **`-krylov` only pays under full Newton** (and other tangent-refreshing
  algorithms). Under `ModifiedNewton`/`Initial` the factorization is already
  being reused by a different mechanism and `-krylov` has nothing to add.

---

## 5. Traps — read this section before shipping a deck generator

### Trap 1 — `-matrixType 1|2` on an unsymmetric tangent silently solves the WRONG system

Half-storage keeps only the `col >= row` half of each element matrix. If the
assembled tangent is genuinely unsymmetric, the lower half is **discarded** and
the run solves the *reflected* system. It converges. It produces a plausible
answer. It is wrong.

Unsymmetric in this fork means, at least: **contact**, **non-associated flow**
(e.g. `LadrunoConcrete3D`), **follower loads**, **corotational transforms**, and
**`LadrunoUP`**.

There is a detector — it compares each element matrix against its mirror and
warns once — but **it samples only the first 3 tangent assemblies of a pattern.**
A tangent that *becomes* unsymmetric later (contact closing at step 50,
non-associated flow after first yield) will **not** be caught. Treat
`-matrixType 1|2` as an assertion you are making, not a check the code performs.

### Trap 2 — `-matrixType 1` (SPD) fails on any indefinite tangent

Softening, buckling, or limit-point models have negative eigenvalues. `1` is a
Cholesky factorization and will stop with a zero-pivot error (`-4`). Use `2`.
The error message says this explicitly.

### Trap 3 — pass numeric option values as **ints**, not strings

```python
ops.system('Pardiso', '-matrixType', 2)      # correct
ops.system('Pardiso', '-matrixType', '2')    # parse fails -> falls back, WARNS
```

The Pardiso factory now degrades to unsymmetric with a warning rather than
returning null. That specific catastrophe is fixed here — but the general
pattern is not fixed everywhere in OpenSees: a factory that returns null leaves
the SOE unset and the analysis silently continues on **`ProfileSPDLinSOE`**,
which on a 3D solid mesh is catastrophically slow. The run then looks *slow*,
never *wrong*, so it gets blamed on the model. Emit ints.

### Trap 4 — `-krylov` is not available with `-matrixType 2`

Intel documents the CG preconditioner for symmetric **positive definite** only.
`-matrixType 2` (symmetric indefinite) warns and falls back to a direct solve.
So the softening/limit-point class — the one that most wants solver speed —
cannot use this lever. Not a bug; a documented scope limit.

### Trap 5 — `-krylov` can change the POST-PEAK equilibrium branch

`-krylov` is bit-identical to a direct solve throughout the physically
meaningful range. Once a softening structure passes its **limit point** under
`LoadControl`, equilibrium is non-unique and the run can land on a different
branch — or give up at a different step — than the same deck without `-krylov`.
Both answers are non-physical there, so this is a reproducibility hazard rather
than a correctness defect. **Keep `-krylov` off when the deliverable is a
post-peak path** (progressive-collapse / AEM work).

### Trap 6 — "perturbed pivots" is not a cosmetic warning

```
WARNING PARDISOGenLinSolver: PARDISO perturbed N pivot(s) ...
```

PARDISO *replaces* a pivot below a threshold rather than failing, so a
near-singular tangent returns `error == 0` and a solution to a matrix that is
**not** the one assembled. If you see this, treat a stalling Newton as a
near-singular tangent, not as a solver hiccup. The symmetric path uses a
threshold five orders looser than the unsymmetric one, so `-matrixType 1|2`
widens the window.

### Trap 7 — threaded PARDISO is not byte-reproducible run-to-run

**This one will bite apeGmsh's regression tests directly.** With
`MKL_NUM_THREADS > 1`, running the *same deck with the same binary twice* can
give results differing in the last bit. Measured: a 14³ Lane-B model at 4
threads returned **two distinct tip displacements across 10 runs, a 5/5 split**
(~1 ULP). At `MKL_NUM_THREADS=1` it is 10/10 identical. It is size-dependent —
small models are stable even threaded, which is exactly why it went unnoticed.

This is MKL's threaded factorization, not a fork bug, and the drift is last-bit,
so it is a **reproducibility** problem, not an accuracy one.

**If you compare outputs byte-for-byte in CI, pin `MKL_NUM_THREADS=1` for those
runs.** Use a tolerance (1e-12 relative is ample) for threaded runs. Do not
"fix" a flaky golden-file test by widening it blindly — check the thread count
first.

---

## 6. Verifying it actually engaged

```python
ops.system('Pardiso', '-matrixType', 1, '-stats')
```

prints once per sparsity pattern:

```
PARDISO stats: n=... nnz(A)=... mtype=2 (symmetric, upper-triangle CSR)
  peak symbolic  = ... MB
  permanent      = ... MB
  peak numeric   = ... MB   (factor + solve)
  TOTAL PEAK     = ... MB   <- the fit/no-fit number
  nnz in factors = ...
  perturbed pivots = 0   refinement steps = 2
```

Check `mtype` is what you asked for. **`TOTAL PEAK` is the number that decides
whether a model fits** — use it, not `nnz`, when sizing. With `-krylov` you also
get one line per solve showing whether CGS answered it or PARDISO refactored.

Threading is *not* reported. Verify it externally (CPU utilisation, or time a
run at `MKL_NUM_THREADS=1` vs `4`). `iparm[2]` does **not** set the thread count
— only the environment does.

---

## 7. Suggested apeGmsh defaults

```python
def solver_directives(model, target):
    if target in ("OpenSeesSP", "OpenSeesMP"):
        return ["Mumps", "-ICNTL14", 200]

    # 0 unsymmetric (always safe) / 1 SPD / 2 symmetric indefinite
    mat_type = 0
    if model.tangent_is_symmetric:
        mat_type = 2 if model.has_softening else 1   # traps 1 and 2

    args = ["Pardiso"]
    if mat_type:
        args += ["-matrixType", mat_type]            # ints, not strings (trap 3)

    # -krylov pays only under a tangent-refreshing algorithm, only at scale,
    # and is not legal for mat_type 2 (trap 4). Off for post-peak work (trap 5).
    if (mat_type != 2 and model.ndof >= 25_000
            and model.algorithm == "Newton"
            and not model.wants_post_peak_path):
        args += ["-krylov", 6]

    return args
```

`model.tangent_is_symmetric` must be **False** whenever the model contains
contact, `LadrunoUP`, non-associated flow, follower loads, or corotational
transforms (trap 1). If apeGmsh cannot determine that reliably, hard-code the
default `["Pardiso"]` — it is correct for every model, just not the fastest for
some.

---

## 8. Environment note

For live apeGmsh-on-fork runs use the `opensees_env` (Python 3.12) environment
and point `APEGMSH_OPENSEES_BIN` at the fork build; the module imports as
`opensees`. The MKL runtime DLLs must sit beside `opensees.pyd` — the fork's
`build.bat` stages them into `dist/bin` automatically, so ship that whole
directory, not just the `.pyd`.
