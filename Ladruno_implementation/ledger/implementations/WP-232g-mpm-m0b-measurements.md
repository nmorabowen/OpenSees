---
wp: WP-232g
title: "WP-232g: the M0b measurement harness (ADR 232 phase 2)"
pr: "#1089"
date: 2026-11-15
status: "shipped"
class_tags: []
---
| **WP-232g — M0b measurement harness** ([[232_ladruno_mpm_general_driver_adr]] §9 phase 2, §5.5; the pre-registered table in `.claude/skills/ladruno-mpm/references/wp-briefs.md` "232g") — `Ladruno_scripts/mpm_m0b/`: M1 (mode-A rebuild t_rb against one Newton iteration t_it, R, k, f, the rebuild's sub-phases, per (N, quasi-static / transient, elastic / DruckerPrager, PARDISO / UmfPack, renormalize / exclude / ghost) on a laterally confined 3D block), M2 (mode-F throughput and the material share on a 3D granular column collapse, ElasticIsotropic / DruckerPrager / LadrunoJ2; the SANISAND row projected from a single-point stdBrick driver and validated on DruckerPrager), M3 (private bytes per particle-element against N, modes A and F, with a profiler `-memory` attribution), M4 (the 232f cantilever: Newton iterations, tip at three mass cutoffs, formTangent's share). Every run a fresh process with the thread environment pinned; `--pyd` selects the binary; `--smoke` = N 10³, one repeat. Instrumentation is profiler scopes only: `StaticAnalysis` / `DirectIntegrationAnalysis` `dc.*`, mode A's `mpm.a.*`, mode F's deep `mpm.f.material` (per-call `ElemScope` in the material adapter). | Tooling + profiler scopes | — | `Ladruno_scripts/mpm_m0b/{README.md,common.py,models.py,m1_rebuild.py,m2_explicit.py,m3_memory.py,m4_cantilever.py,run_all.py}`, `SRC/domain/pattern/ladrunoMPM/{LadrunoMPMImplicit.cpp,LadrunoMPMPattern.cpp}`, `SRC/analysis/analysis/{StaticAnalysis.cpp,DirectIntegrationAnalysis.cpp}` | shipped | [#1089](https://github.com/nmorabowen/OpenSees/pull/1089) |
