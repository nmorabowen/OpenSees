---
wp: WP-232
title: "PartitionedDomain::analysisStep never tells its Subdomains the domain changed: the notification loop's condition is impossible (upstream)"
date: 2026-11-13
---
### `PartitionedDomain::analysisStep` never tells its Subdomains the domain changed — the notification loop's condition is impossible (upstream bug)

- **Bites:** OpenSeesSP / MP runs in which a subdomain's topology changes mid-analysis and something relies on the partitioned domain propagating `domainChange()` back to every subdomain.
- **Why:** `PartitionedDomain.cpp:1147-1156` — inside `if (domainChangedAnySubdomain == true)`, the propagation loop is `while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false))`, which is false on entry, so `theSub->domainChange()` never runs; only the partitioned domain marks itself changed (`:1148`). Related: `PartitionedDomain::analysisStep` calls the base `Domain::analysisStep` (`:1159`), but `Subdomain::analysisStep` (`Subdomain.cpp:808`) and `ShadowSubdomain::analysisStep` override it and delegate to their decomposition analysis, so a hook placed in `Domain::analysisStep` never fires inside a subdomain.
- **Workaround / status (2026-11-13, WP-232 — found by review, confirmed in-tree, not fixed):** upstream code with no `// Ladruno` marker. The one-token fix (`== true` in the loop condition) has an unmeasured effect on SP / MP runs; fix it in a WP with an SP / MP deck that changes topology mid-run, not as a drive-by.
