---
wp: PR-782
title: "782 -- upstreamable-table row(s)"
pr: "#782"
files: ["`CMakeLists.txt` (root)"]
table: "upstreamable"
legacy_seq: [498]
---
| `CMakeLists.txt` (root) | `# Ladruno (TIMs, 2026-09-03, second opt-in)`: **`option(LADRUNO_MKL_PARDISO_LINUX_THREADED OFF)`** — swaps the Linux PARDISO opt-in's MKL layer from `mkl_sequential` to `mkl_gnu_thread` + `libgomp`; refuses to configure if `LAPACK_LIBRARIES`/`BLAS_LIBRARIES` still name `mkl_sequential`, and (#782 verification) refuses `LADRUNO_MKL_FEAST_LINUX=ON` beside it, because the FEAST block resolves `mkl_sequential` into the same cached `LADRUNO_MKL_SEQ` and the threaded lookup would silently no-op. Default OFF ⇒ every existing build byte-identical. Verified on esmeralda (spack oneAPI MKL 2024.2): links `libmkl_gnu_thread.so.2` + `libgomp.so.1`, no sequential layer; 16³/24³ stdBrick block matches UmfPack to all printed digits; 46 875 DOF solve 3.39 s → 1.44 s at 1 → 8 threads. | [#782](https://github.com/nmorabowen/OpenSees/pull/782) |
