---
wp: PR-739
title: "739 -- 2 vanilla row(s)"
pr: "#739"
files: ["`SRC/domain/domain/Domain.cpp`", "`SRC/domain/domain/Domain.h`"]
table: "main"
legacy_seq: [37, 38]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno` ADR-78 P2: `sendSelf`/`recvSelf` carry the contact engine's DEFINITIONS — `domainData` ID 17→19 (slot 17 = packed size, 0 when no engine; slot 18 = the Vector's dbTag), the definitions Vector sent after the Parameters, recv rebuilds an empty engine / VERIFIES a populated one. Plus `dbContact` member init at all 7 db-tag sites. **Stream-format change:** a pre-P2 (or upstream) database cannot be restored by a P2+ build — see [[LEDGER_quirks]]. Response byte-identity unaffected (analysis path untouched). | [#739](https://github.com/nmorabowen/OpenSees/pull/739) |
| `SRC/domain/domain/Domain.h` | `// Ladruno` ADR-78 P2: `int dbContact` member (dbTag for the contact-engine definitions Vector, sibling of dbNod/dbEle/…). | [#739](https://github.com/nmorabowen/OpenSees/pull/739) |
