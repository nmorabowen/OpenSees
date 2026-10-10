---
wp: LEGACY
title: "Matrix::Solve and Matrix::Invert run on a PROCESS-WIDE scratch buffer that they FREE AND REALLOCATE"
legacy_seq: 451
---
### `Matrix::Solve` and `Matrix::Invert` run on a PROCESS-WIDE scratch buffer that they FREE AND REALLOCATE

`Matrix::matrixWork` / `Matrix::intWork` (`SRC/matrix/Matrix.cpp:51-52`) are class
statics. `Matrix::Solve(Vector&,Vector&)` (`:373`), `Solve(Matrix&,Matrix&)`
(`:461`) and `Invert` (`:567`) each contain the same block: if the matrix is
bigger than the fixed work area, `delete [] matrixWork; matrixWork = new
double[dataSize];` — then they copy `data` into it and factor in place.

So any code path reachable from a threaded loop that calls either method is not
merely a data race on a shared buffer: it is a **use-after-free**, because one
thread can free the buffer another thread is mid-factorization on. A grep for
`static` inside the element or material file finds nothing — the static is three
directories away in `SRC/matrix/`.

This is the concrete reason WP-107's allowlist refuses `ManzariDafalias`
`IntScheme 2` / `4` (their `NewtonSol*`/`NewtonIter*` call both) and
`LadrunoQuad -formulation eas` (its static condensation inverts `Kaa`), even
though nothing about those paths *looks* shared. The `Matrix` **constructors**
also lazily allocate that buffer, which is benign in practice only because
thousands of `Matrix` objects are built during model construction on the master
thread — do not rely on that if a threaded phase is ever added before model
build completes.
