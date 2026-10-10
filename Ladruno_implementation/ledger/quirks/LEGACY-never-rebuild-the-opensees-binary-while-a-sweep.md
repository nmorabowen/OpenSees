---
wp: LEGACY
title: "Never rebuild the OpenSees binary while a sweep is running"
legacy_seq: 228
---
### Never rebuild the OpenSees binary while a sweep is running
- **Bites:** a multi-mode sweep re-execs the binary once per mode (`openseesmp.sh` → `exec`). Rebuilding mid-job silently swaps the executable between modes, so the A/B comparison spans two different binaries and nothing in the output says so.
- **Workaround/status:** hold rebuilds until `squeue` is clear, or build to a distinct path and point the wrapper at it explicitly. *2026-07-26 (ADR-75 P2h).*
