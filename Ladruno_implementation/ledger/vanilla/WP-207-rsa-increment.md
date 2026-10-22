---
wp: WP-207
title: "responseSpectrumAnalysis -increment: per-mode RSA about the committed (gravity) state"
date: 2026-10-22
pr: "TBD"
table: "main"
files: ["`SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}`"]
---
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}` | `// Ladruno WP-207`: opt-in `-increment [-sign ±1] [-revert]`. `solveMode()` sets trial = committed + sign·u_k (instead of replacing it); `endMode()` updates the domain and calls `Domain::record()` but does NOT commit; `beginMode()` reverts to the last commit first, and `-revert` reverts at the end. Standalone `responseSpectrumAnalysis -revert` reverts the domain. `setIncrement()` + in-class members (`m_increment=false`, `m_sign`, `m_revert_at_end`) keep the ctor untouched. `-increment` absent ⇒ byte-identical to stock. | TBD |
