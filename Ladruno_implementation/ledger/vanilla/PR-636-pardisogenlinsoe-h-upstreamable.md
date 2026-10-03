---
wp: PR-636
title: "636 -- upstreamable-table row(s)"
pr: "#636"
files: ["`SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.{h,cpp}`"]
table: "upstreamable"
legacy_seq: [373]
---
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.{h,cpp}` | `// Ladruno ADR-75 P1f`: **`addA` scatter — linear scan → binary search.** The stock inner loop rescanned the whole CSR row for every one of `idSize²` element entries, i.e. `O(idSize² × rowlen)`; for a 3D brick that is 24×24 lookups each walking a ~81-entry row. ADR-75b's L3-0 profile measured it at **1699 ms of an ~11.9 s step, 1.28× slower than UmfPack's** equivalent scatter, and it only grew in relative terms as P1a/P1d/P1e made the *solve* cheaper. New `ops_pardiso_findCol()` binary search; legal because the ascending-CSR invariant is **enforced**, not assumed — `setSize` already rejects the matrix (`size = 0`, `return -1`) if any row's columns are not strictly ascending. **Measured (Lane B, 4 threads, interleaved A/B against a snapshot of the pre-fix binary): 1.098× at 26.5k DOF, `-matrixType 1` 1.029×.** Exactness verified at `MKL_NUM_THREADS=1`, where old and new agree **bit-for-bit over 10 runs each** (the 4-thread comparison is useless for this — see the determinism quirk). Also adds a **free miss detector**: the `k >= 0` test already existed, so `else missing = 1` costs a never-taken branch and turns "an element entry has no CSR slot ⇒ its stiffness is silently discarded" from undetectable into a once-per-SOE warning; the old linear scan simply fell off the end of the row. | [#636](https://github.com/nmorabowen/OpenSees/pull/636) |
