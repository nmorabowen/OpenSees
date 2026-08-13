# ADR-78 P4 — parallel-contact validation battery

The fork-side validation battery for ADR-78 (parallel MPI contact), phase P4:

    python run.py <build-dir>

`<build-dir>` holds `OpenSees.exe` / `OpenSeesMP.exe` (normally `dist\bin` of
the tree you built). Outputs land in `./_out` (wiped per run). See the
docstring in `run.py` for the case list and `pound.tcl` / `tie.tcl` headers
for the deck shapes; measured results are recorded in the ADR's §P4 log
(`Ladruno_implementation/78_ladruno_parallel_contact_adr.md`).

Companion harnesses: `../contact_parallel/` (P0), `../contact_p2/` (P2 —
the full soft-refusal matrix + DB round-trip; this battery re-pins only the
partitioned-soft refusal on its own pounding deck).
