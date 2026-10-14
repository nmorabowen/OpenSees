---
wp: WP-192
title: "Making an element/material thread-safe: static scratch -> thread_local, and the four traps that survive the conversion"
date: 2026-10-14
---
### Making an element or material thread-safe: `static` scratch -> `thread_local`, and the traps that survive the conversion

- **What bites.** OpenSees elements and materials return references to `static Matrix` /
  `static Vector` scratch (function scope and class scope) and share class-wide pools
  (`FE_Element::theMatrices`, `TransformationFE::modMatrices`/`dataBuffer`,
  `Element::theMatrices` for Rayleigh, `Matrix::matrixWork`). Two elements on two threads
  then write the same buffer: silent wrong answers or NaN, not a crash.
- **The conversion that works (WP-192).** `static X` -> `static thread_local X` (and the
  class-static definition `X C::x(...)` -> `thread_local X C::x(...)`): zero bytes per
  element, identical arithmetic on one thread, so the serial path stays byte-identical.
  The shared POOLS are not converted; they get a per-thread twin used only inside a threaded
  region (`ladruno_inThreadedRegion()`), because their pointers are cached per object at
  construction.
- **Trap 1 — a static sized by its FIRST caller.** `static Matrix BB(6, m_num_dofs)` in one
  ASDAbsorbingBoundary3D kernel was sized once and never resized; per-thread copies are each
  sized by a different first caller. Resize per call (WP-192 did; the serial code was
  already reading past the end for any later element with another dof count).
- **Trap 2 — a size DECISION that reads the shared size.** `Matrix::addMatrixTripleProduct`
  takes a different code path when its work area is too small, and `Invert` passes the work
  size to LAPACK as `lwork`. A per-thread work area must answer with the SHARED logical size
  (never grown inside a region), or a worker takes another path and the bits change.
- **Trap 3 — statements captured in the wrong branch.** Wrapping a growth block in
  `if (in region) {...} else {...}` swallowed the `x = b;` that followed it in
  `Matrix::Solve(Matrix&, Matrix&)`: serial unchanged, threaded ForceBeamColumn3d inverted
  garbage. The full-field threaded-vs-serial hash caught it on the second step;
  `LADRUNO_THREADS_LOOPS=1|2|4|8` bisected it to the update loop in one run.
- **Trap 4 — a process-wide accumulator in commit.** ASDConcrete3D/1D add every point's
  IMPL-EX error to `GlobalParameters` (max, sum, count) in `commitState`, and
  ASDAbsorbingBoundary3D adds to the `ABSORB_LEAK` energy channel. A threaded commit
  races on them and reorders an FP sum. Serialise (a mutex inside the region) or keep that
  element's commit serial (`ladrunoThreadSafeCommit() == false`).
- **How to test.** Byte-identity at N threads is only meaningful on a deterministic solver:
  OpenSeesMP + MUMPS is not run-to-run reproducible (WP-188 quirk), so compare threaded vs
  serial on a serial Pardiso deck with `MKL_NUM_THREADS=1` (or UmfPack), every node every
  step at 17 digits.
