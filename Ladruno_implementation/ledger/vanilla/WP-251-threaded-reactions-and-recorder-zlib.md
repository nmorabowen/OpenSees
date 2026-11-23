---
wp: WP-251
title: "WP-251: threaded reaction loop, its per-element opt-out, recorder zlib link"
date: 2026-11-23
files: ["SRC/domain/domain/Domain.cpp", "SRC/domain/domain/Domain.h", "SRC/element/Element.h", "SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.h", "SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.h", "SRC/recorder/CMakeLists.txt"]
---
| `SRC/domain/domain/Domain.cpp` | WP-251: `calculateNodalReactions` tries `ladrunoThreadedReactions(flag)` first (parallel element forces, serial-order scatter, bit-identical); the serial loop is unchanged and runs whenever it declines | — |
| `SRC/domain/domain/Domain.h` | WP-251: declares `ladrunoThreadedReactions(int)` | — |
| `SRC/element/Element.h` | WP-251: `virtual bool ladrunoThreadSafeReaction() const` (default = `ladrunoThreadSafeForm()`) | — |
| `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.h` | WP-251: `ladrunoThreadSafeReaction()` false — it overrides `addResistingForceToNodalReaction` | — |
| `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.h` | WP-251: `ladrunoThreadSafeReaction()` false — it overrides `addResistingForceToNodalReaction` | — |
| `SRC/recorder/CMakeLists.txt` | WP-251: `Ladruno_ChunkCompress.{h,cpp}`; `OPS_Recorder` links `ZLIB::ZLIB` (it calls `compress2`) | — |
