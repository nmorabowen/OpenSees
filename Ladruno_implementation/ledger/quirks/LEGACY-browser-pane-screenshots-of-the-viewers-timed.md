---
wp: LEGACY
title: "Browser-pane screenshots of the viewers timed out in two sessions — assert through the DOM (JS eval) instead"
legacy_seq: 480
---
### Browser-pane screenshots of the viewers timed out in two sessions — assert through the DOM (JS eval) instead
- **Bites:** 2026-05-31 (P8, #51): the preview tool's screenshot call kept timing out on the profiler viewer, while DOM queries through JS eval and the console log worked; the session noted "assert via DOM queries, not screenshots". 2026-07-04 (#484): the next viewer session tried screenshots again and hit the same timeout ("still flaky"). The second loss was already written down.
- **Workaround/status:** make the DOM the evidence — count the icicle's SVG `rect` frames, read legend and table text, check the console for errors — and treat a screenshot as a bonus when it works. Cause unknown; not re-tested in WP-121. *Moved from session memory, WP-121.*
