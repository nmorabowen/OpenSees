---
wp: WP-212
title: "responseSpectrumAnalysis: mass-shortfall report and guard (-requireMass), -combine displacements-only warning"
date: 2026-10-26
pr: "#1021"
table: "main"
files: ["`SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}`"]
---
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}` | `// Ladruno WP-212`: new private `checkMass(mp)`, called after `check()` on every path (`analyze()`, `analyze(mode)`, `analyzeCombined()`, so `-increment` too). It prints the cumulative effective-mass ratio of the excited direction once per call (hidden by `LADRUNO_OPENSEES_QUIET`), warns below 90% (ASCE 7-22 §12.9.1.1, NEC-SE-DS §6.2.2e), and returns −1 when the direction has no free mass or every Γ of it is zero (the WP-211 diaphragm-mass signature). The new option `-requireMass $ratio` in (0, 1] turns the warning into an error (`setRequiredMassRatio()` + in-class `m_require_mass=0`; ctor untouched). `analyzeCombined()` prints one warning per call: `-combine` combines nodal displacements only, so forces read afterwards are not modally combined. Only stderr changes: results are byte-identical when the mass is complete (ADR 45 G-F). | #1021 |
