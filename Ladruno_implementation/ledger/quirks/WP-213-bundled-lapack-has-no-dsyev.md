---
wp: WP-213
title: "Linux link: undefined reference to dsyev_ (bundled LAPACK subset)"
date: 2026-10-26
pr: "#1023"
---
### Linux Zone-A link fails with `undefined reference to 'dsyev_'` although Windows builds

- **Bites:** a new fork file calls a LAPACK routine (here `dsyev_` in
  `LadrunoRitzKernel.cpp`); the Windows build links (MKL provides all of
  LAPACK), the Ubuntu Zone-A build fails at the final `OpenSees` link.
- **Why:** the Linux build links the LAPACK subset compiled from
  `OTHER/LAPACK/` (`OTHER/LAPACK/CMakeLists.txt` lists the files). It has
  `dsyevx`, `dsygvx`, `dggev`, `dgeev`, `dgesvd`, ... but not `dsyev`.
  `LadrunoCMSLocalLanczos.cpp` also calls `dsyev_`, so it is not a precedent.
- **Fix/status:** check `OTHER/LAPACK/CMakeLists.txt` before calling a LAPACK
  routine from fork code; `LadrunoRitzKernel` now uses `dsyevx` with
  `RANGE = 'A'` and `ABSTOL = 2·DBL_MIN`. Learned 2026-10-26 (WP-213).
