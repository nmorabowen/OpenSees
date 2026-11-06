---
wp: LEGACY
title: "Graph::getVertexPtr cost is set by the STORAGE the graph was built on — ArrayOfTaggedObjects is O(1) for dense tags, MapOfTaggedObjects is O(log n)"
legacy_seq: 208
---
### `Graph::getVertexPtr` cost is set by the STORAGE the graph was built on — `ArrayOfTaggedObjects` is O(1) for dense tags, `MapOfTaggedObjects` is O(log n)
- **Bites:** an RCM BFS (or any per-vertex loop calling `getVertexPtr`) on a `MapOfTaggedObjects`-backed graph pays a `std::map::find` per read — ~1.6×10⁹ `find`s ≈ 25-45 min at 19 M nodes, the whole post-T0 numberer residual. It reads like "RCM is just slow"; it's actually the storage choice underneath `getVertexPtr`. `AnalysisModel::getDOFGroupGraph` builds on a `MapOfTaggedObjects`, so every parallel-numberer merge inherited it.
- **Why:** `Graph::getVertexPtr` → `TaggedObjectStorage::getComponentPtr`. `ArrayOfTaggedObjects` returns `theComponents[tag]` directly (O(1)) when the tag sits at its own index — true for the dense 0..N-1 tags graphs use; `MapOfTaggedObjects` is a red-black tree (O(log n)) always.
- **Workaround/status:** the T1 lever — `LadrunoParallelNumberer` builds its merged graph on **owned `ArrayOfTaggedObjects` storage** + a `tag → Vertex*` mirror, so RCM's BFS reads are O(1) and edge inserts skip `Graph::addEdge`'s two lookups entirely (`Vertex::addEdge` direct). numberDOF 15.5 → 6.4 s at 2.0 M ([#594](https://github.com/nmorabowen/OpenSees/pull/594)). Bit-identity holds because the adjacency `ID::insert` is a sorted set (insertion-order-canonical) and dense tags iterate ascending in BOTH storages. Rule: if you build a Graph you will read per-vertex, build it on array storage and confirm tags are dense. *2026-07-22.*
