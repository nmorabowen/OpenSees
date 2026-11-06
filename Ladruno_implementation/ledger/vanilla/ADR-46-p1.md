---
wp: ADR-46
title: "ADR46 P1 -- 1 vanilla row(s)"
files: ["`SRC/domain/domain/Domain.{h,cpp}`"]
table: "main"
legacy_seq: [260]
---
| `SRC/domain/domain/Domain.{h,cpp}` | `// Ladruno` ADR46 P1: (1) domain-level copy of the Rayleigh factors + `getRayleighDampingFactors()` getter — `setRayleighDampingFactors` was cascade-and-discard (elements/nodes get the values, the Domain kept nothing), so no analysis-side consumer could read back what damping the model carries; the complex-modal Route-A closed form (and later ADR 44) needs exactly that. (2) `getNumEigenvalues()` non-exiting presence probe — `getEigenvalues()` **exit(-1)s** when never set (kernel-killer). (3) `clearAll()` now resets the Rayleigh copy AND the upstream `theEigenvalues`/`theEigenvalueSetTime` — the spectrum surviving `wipe()` was a latent upstream leak the domain-coupled `complexEigen` turned observable (P1 Opus gate CRITICAL-1/2/3; see [[LEDGER_quirks]]). All strictly additive except the clearAll resets. | ADR46 P1 |
