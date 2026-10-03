---
wp: LEGACY
title: "A singular tangent can pass on Windows and fail on Linux — a green local run is not proof the model is well-posed"
legacy_seq: 248
---
### A singular tangent can pass on Windows and fail on Linux — a green local run is not proof the model is well-posed
- **Bites:** `tests/test_ladruno_response_tokens.py::test_distributing_coupling_tokens` tied a `LadrunoDistributingCoupling` reference node to two COLLINEAR hex corners, leaving the reference rotation about that line unconstrained. Same source, same `system BandGeneral` + `numberer RCM`: the **Windows/MKL build returned `analyze == 0`**, Zone-A on **Ubuntu failed the step** with `BandGenLinLapackSolver::solve() -factorization failed, matrix singular U(i,i) = 0`. The LAPACK behind `BandGeneral` differs between the two (MKL vs the CI's reference/OpenBLAS) and so does the pivot at which it calls a matrix singular.
- **Consequence:** a model that is *actually* degenerate can develop and pass locally for as long as you only run it on Windows, then fail the moment CI builds it. The green local run proved nothing; the Linux failure was the correct answer. Same family as the ADR-76 gate ("a singular matrix must not report SUCCESS") — that one is about the reported status of a known-singular solve, this one is about the two platforms disagreeing on whether the solve is singular at all.
- **Rule for new test models:** when an element PRINTS a well-posedness warning, treat it as an error in a test. This element said `reference rotation about axis (1, 0, 0) is unconstrained (degenerate independent set)` on the very run that "passed" — the diagnostic was right and the assertion was wrong. For `LadrunoDistributingCoupling` / `LadrunoKinematicCoupling` specifically, the independent node set must SPAN (non-collinear for a rotation-carrying reference node); check `nKept` == the number of rotation axes you expect.
- *2026-07-28 (recorder-token consistency sweep).*
