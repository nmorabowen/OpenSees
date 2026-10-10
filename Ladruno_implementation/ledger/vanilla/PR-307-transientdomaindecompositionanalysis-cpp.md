---
wp: PR-307
title: "307 -- 1 vanilla row(s)"
pr: "#307"
files: ["`SRC/analysis/analysis/TransientDomainDecompositionAnalysis.cpp`"]
table: "main"
legacy_seq: [70]
---
| `SRC/analysis/analysis/TransientDomainDecompositionAnalysis.cpp` | `// Ladruno` ADR-30 (P2): `domainChanged()` already checked `handle()` + `Integrator::domainChanged()` but DROPPED the `doneNumberingDOF()` return (overwritten by `setSize`). Added the `<0` ⇒ return −1 check so a handler diagnostic is not silently swallowed under distributed transient analysis. Defense-in-depth (v1 projection handler is partition-interior only). | [#307](https://github.com/nmorabowen/OpenSees/pull/307) |
