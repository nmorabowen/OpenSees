---
wp: WP-126
title: "PARTITION_REDUCTION on Ladruno recorder output (WP-126)"
pr: "#861"
section: "table"
legacy_seq: 12
---
| **`PARTITION_REDUCTION` on Ladruno recorder output (WP-126)** ([[126_partition_reaction_reduction]]) — under OpenSeesMP each rank's reaction (and unbalanced load) at a node on a partition interface is that rank's partial; apeGmsh's stitch kept the first copy ((0, 10) where the serial reaction is (20, 30)) and `-envelope` of a partial was unrecoverable. The recorder now writes `PARTITION_REDUCTION` = `NONE` | `SUM` | `UNSUPPORTED` on every result group (from `ResultSource::partitionReduction()`; energyBalance = UNSUPPORTED), refuses `-envelope` of a non-NONE source in partitioned runs (warning) and warns that partitioned energy is not mergeable. Schema §7.1 + validator + apeGmsh contract updated; zone_a `tests/test_ladruno_partition_reduction.py` (4) and the openseesmp two-rank gate `mp_reaction_*.py` (stitch per contract = serial). apeGmsh stitch fix = separate apeGmsh PR. From WP-120 OQ1. | recorder fix + schema addition | — | `SRC/recorder/{Ladruno_ResultIO.h, Ladruno_DomainResults.h, Ladruno_Sinks.{h,cpp}, LadrunoRecorder.cpp}`, `Ladruno_scripts/ladruno_recorder_tests/{ladruno_format.py, test_ladruno_harness.py, mp_reaction_model.py, mp_reaction_check.py, README.md}`, `tests/test_ladruno_partition_reduction.py`, schema v1 §7.1, `ladruno_apegmsh_contract.md`, `Ladruno_implementation/wp126_partition_reactions/*` | **ready (owner merges)** | #861 |
