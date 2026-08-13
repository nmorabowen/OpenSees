# ADR-78 P2 gate harness

Gates both P2 deliverables ([[78_ladruno_parallel_contact_adr]] §P2 LOG):
the `-soft`-under-partitioning refusal and the contact-definitions
serialization through `Domain::sendSelf/recvSelf`.

    python run.py <build-dir>

`<build-dir>` holds `OpenSees.exe` / `OpenSeesMP.exe` / `opensees.pyd`
(`dist\bin` of the tree you built, or `build\build\Release`). Named
explicitly on purpose — a bare `OpenSees` on PATH has resolved to another
session's build on this machine more than once.

| case | deck | expects |
|---|---|---|
| serial-auto | `../contact_parallel/serial.tcl` | runs; reference values |
| mp-auto | `../contact_parallel/mp.tcl`, np2 | response parity vs serial ≤1e-12 (measured 1.6e-14) — the "-kn auto needs no change" gate |
| serial-soft | `serial_soft.tcl` | RUNS — refusal must not fire at np=1 |
| mp1-soft | `serial_soft.tcl` under `OpenSeesMP -n 1` | RUNS — refusal is np-keyed |
| mp2-soft | `mp_soft.tcl`, np2 | named FATAL refusal; tail line never prints |
| mp2-soft-mortar | `mp_soft_mortar.tcl`, np2 | same refusal in the lane P1 never guarded (pre-P2 this RUNS silently — the honest mutation) |
| db-roundtrip | `db_roundtrip.py` | contact defs survive save→wipe→restore exactly; live-restore verify path silent; vanilla model round-trips too |
| db-all-lanes | `db_roundtrip_all_lanes.py` | the review follow-up: ALL FOUR definition lanes (NTS+friction/geomtan, mortar tie, rigid plane+visc, mortar friction+edge-edge block) round-trip exactly in one model; the packed Vector is byte-deterministic; six per-field mutations (kt, mu, edgeKn, epsTie, plane kn, plane normal) each provably change the packed bytes; the live-restore verify instrument stays silent on a match and warns on a mutated definition set |

Negative control: run the same driver against a pre-P2 build
(`C:\Program Files\Ladruno\OpenSees\bin` at the time of writing) —
`mp2-soft` / `mp2-soft-mortar` / `db-roundtrip` must FAIL there
(mp2-soft dies on the P1 slave-mass abort instead of the P2 refusal;
mp2-soft-mortar runs to completion; db-roundtrip restores a contact-free
model, rt=0.0).

Deck-design note that cost a round: both MP soft decks carry `rho > 0` ON
PURPOSE. With a massless model the P1 zero-slave-mass abort fires first and
the refusal test can "pass" via the wrong abort. With mass, the P2 refusal
is the only abort available.
