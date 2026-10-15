---
wp: WP-194
title: "OpenSeesMP + system Mumps on 2 ranks is not byte-reproducible run to run"
date: 2026-10-14
pr: "#982"
---
### Two runs of the same OpenSeesMP deck with `system Mumps` give different last bits (same binary, same input)
- **Bites:** any byte-identity gate or A/B that compares OpenSeesMP + MUMPS outputs from ONE run per arm. Seen on the WP-188 frame deck (`tests/wp188_factor_reuse/frame_gm.tcl`, 300 steps, 2 ranks, Windows, Intel MPI, MUMPS 5.9.0, MKL sequential, 1 thread/rank), on the pre-WP-194 binary as much as on the WP-194 one: elastic, 1 run in 10 gives a different roof-history md5; fiber, 6 to 7 runs in 10 differ from the modal one (last digits only, relative ~1e-16 to 1e-15; the iteration count never changed). `MKL_CBWR=AUTO,STRICT` does not remove it.
- **Why (inferred, not traced):** with distributed entry (`ICNTL(18)=3`) MUMPS gathers each rank's triplets and sums duplicates in message-arrival order, so the floating-point summation order of shared entries can change between runs. Not a fork bug and not WP-194's: the stock double path does it with the WP-194 option absent.
- **Do:** compare OpenSeesMP + MUMPS results as distributions (several runs per arm) or with a tolerance (1e-12 relative is ample); never conclude "byte-identical" or "changed" from n = 1, the same lesson as the threaded PARDISO entry (`LEGACY-threaded-mkl-pardiso-is-not-byte-reproducible`). WP-190's `test_stats_off_prints_nothing_new_and_same_answer` compares two single runs exactly and can flake for this reason.
