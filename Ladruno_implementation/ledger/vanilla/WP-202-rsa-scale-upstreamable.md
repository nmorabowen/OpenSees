---
wp: WP-202
title: "responseSpectrumAnalysis -scale was parsed and silently ignored"
date: 2026-10-22
table: "upstreamable"
files: ["`SRC/analysis/analysis/ResponseSpectrumAnalysis.cpp`", "`SRC/analysis/analysis/ResponseSpectrumAnalysis.h`"]
---
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.cpp` | `// Ladruno WP-202`: `-scale` was stored in `m_scale` and never used, so every spectral ordinate ignored it (no warning). `getSa()` now returns `m_scale * getSaUnscaled(T)`, so the per-mode (`analyze`), single-mode (`-mode`) and `-combine` paths all honour it, for both the `-Tn/-Sa` list and the timeSeries forms. Upstream bug, candidate for upstream. Test `tests/test_wp202_rsa_scale.py`. | WP-202 |
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.h` | `// Ladruno WP-202`: declares `getSaUnscaled(double)`. | WP-202 |
