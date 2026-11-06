---
wp: WP-217
title: "IMPL-EX clock in ModIMKPeakOriented and Bilin"
date: 2026-10-29
files: ["SRC/material/uniaxial/ModIMKPeakOriented.cpp", "SRC/material/uniaxial/ModIMKPeakOriented.h", "SRC/material/uniaxial/Bilin.cpp", "SRC/material/uniaxial/Bilin.h"]
---
| `SRC/material/uniaxial/ModIMKPeakOriented.cpp` | WP-217: `-implex` factor from `LadrunoImplexClock` (frozen per step) instead of the live `ops_Dt`; `-implexDt step\|pseudo`; clock bit packed in the implex slot | TBD |
| `SRC/material/uniaxial/ModIMKPeakOriented.h` | WP-217: `implexClock_` member, `setImplexClock`/`getImplexClock` | TBD |
| `SRC/material/uniaxial/Bilin.cpp` | WP-217: as ModIMKPeakOriented | TBD |
| `SRC/material/uniaxial/Bilin.h` | WP-217: as ModIMKPeakOriented | TBD |
