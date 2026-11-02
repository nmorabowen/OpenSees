---
wp: WP-219
title: "Concrete02 returns the committed response at the committed strain"
date: 2026-11-01
table: "upstreamable"
files: ["SRC/material/uniaxial/Concrete02.cpp"]
---
| `SRC/material/uniaxial/Concrete02.cpp` | WP-219: at `\|deps\| < DBL_EPSILON` set `sig = sigP; e = eP` instead of keeping the previous trial's values; the stress jump trapped force-based elements in an endless cycle. Upstream defect | TBD |
