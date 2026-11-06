---
wp: ADR-46
title: "ADR46 P2, P3 -- 1 vanilla row(s)"
files: ["`SRC/domain/node/Node.{h,cpp}`"]
table: "main"
legacy_seq: [257]
---
| `SRC/domain/node/Node.{h,cpp}` | `// Ladruno` ADR46 P2: `getNumEigenvectors()` non-exiting presence probe — `getEigenvectors()` **exit(0)s** when unset (fully-fixed nodes never enter the eigen analysis), so the reduced-operator projection must be able to probe before gathering (zero rows = the fixed node's physical mode shape). Additive (decl + 4-line body; Matrix is forward-declared in Node.h so the body lives in the .cpp). **P3:** complex mode-shape storage — `setNumComplexEigenvectors`/`setComplexEigenvector(mode, re, im)`/`getComplexEigenvectors{Re,Im}`/`getNumComplexEigenvectors` (Re/Im Matrix pair mirroring the real `theEigenvectors`; in-class-initialized members so the 4 constructors stay untouched; destructor deletes; NOT serialized — complexEigen is serial-only). | ADR46 P2, P3 |
