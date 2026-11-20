---
wp: WP-242
title: "ASDShellQ4/T3 answer localAxes (material frame for MODEL/LOCAL_AXES)"
date: 2026-11-19
files: ["SRC/element/shell/ASDShellQ4.cpp", "SRC/element/shell/ASDShellT3.cpp"]
---
| `SRC/element/shell/ASDShellQ4.cpp` | WP-242: new `localAxes` element response (id 5, 9 packed cosines e1\|e2\|e3) = the section/material frame: reference local system rotated by `m_angle` about its normal, so `-local` (projected on the plane) or the default mid-edge axis. Lets `recorder ladruno` write `MODEL/LOCAL_AXES` for the shells, which the viewer needs to place crack vectors. Additive; ids 1-4 untouched | — |
| `SRC/element/shell/ASDShellT3.cpp` | WP-242: same `localAxes` response (id 5); default e1 = P2 - P1 | — |
