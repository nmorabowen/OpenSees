---
wp: ADR-78
title: "A fork database written pre-ADR-78-P2 cannot be restored by a P2+ build (domainData 17→19) — and the failure names nothing useful (ADR 78 P2)"
legacy_seq: 296
---
### A fork database written pre-ADR-78-P2 cannot be restored by a P2+ build (domainData 17→19) — and the failure names nothing useful (ADR 78 P2)
- **Bites:** `Domain::sendSelf/recvSelf` grew the leading `domainData` ID from 17 to 19 slots (slot 17 = the contact-engine definitions Vector's packed size, slot 18 = its dbTag) so `database File` save/restore carries contact definitions. `FileDatastore` keys its record files BY OBJECT SIZE (`.IDs.17` vs `.IDs.19`), so a P2+ build restoring a pre-P2 database (or an upstream one) looks for a 19-slot record that does not exist and fails with the generic `Domain::recv - channel failed to recv the initial ID` — nothing says "format changed".
- **Why it evades the usual guards:** both builds are green on their own round-trips; only the CROSS-build restore breaks, and saved databases usually outlive the build that wrote them by exactly long enough to forget this.
- **Rule:** a saved OpenSees database is build-lineage-scoped. Re-save after upgrading across P2 (#this PR); do not archive `database File` outputs as long-term state.
