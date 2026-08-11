# contact_parallel — the ADR-78 verification harness

Every number quoted in ADR-78 P0 / P0.5 / P1 came from these decks. They are here
so P2 and P4 can re-run them instead of re-deriving a reference, and so a
regression has something to fail against.

## The model

Two stacked unit `stdBrick`s (E = 2e7, ν = 0), NTS penalty contact, `-kn auto`,
`-outward 0 0 1`, lateral DOFs fixed so it reduces to a 1-D column. The bottom
block is the **master**, the top block the **slave**, and the interface starts
with a 1e-3 **interpenetration** — with an exactly zero gap the penalty is
inactive at iteration 1, the top block has a free rigid-body mode, and Newton
oscillates (measured: no convergence in 50 iterations, residual stuck at 1.7e4).
That is a property of the test model, not of contact.

The serial answer is analytically exact, which is what makes it a reference:
each block compresses `PL/EA` = 5e-4, and the base reaction equals the applied
10 kN.

## Reference values

| | serial | 2-rank |
|---|---|---|
| `w15` (top) | −5.624999999999997e−3 | −5.624999999999919e−3 |
| `w11` (interface) | −5.124999999999998e−3 | −5.124999999999914e−3 |
| `w5` (master face) | −5.000000000000000e−4 | −5.000000000000034e−4 |
| `ΣR` (base) | 10000.0 | 10000.0 |

Explicit lane (`CentralDifferenceLadruno`, 200 × 5e-5 s): serial and 2-rank are
**bit-identical** — `w15` 0.0037458245314017234, `w11` 0.003692142684617696,
`w5` −4.5278662184168605e−05.

Rank 0's **ghost** `w11` must equal rank 1's native `w11` bit-for-bit. That
single equality is the whole ADR-78 premise: a contact force assembled into a
ghost DOF reaches the owning rank through the distributed SOE.

## Files

* `serial.tcl` / `mp.tcl` — the reference pair. `mp.tcl` is the shape apeGmsh
  ADR 0092 emits: contact on the master-owning rank only, the slave interface
  ghosted as `node` + replayed `fix`, **no ghost mass** (INV-6/INV-7).
* `mp_noghost.tcl` — mutation: the ghost declarations removed. Must FATAL and
  tear the job down in ~1 s. Before ADR-78 P1 this warned, limped to
  `analyze=-3`, hung for minutes and left orphaned ranks alive.
* `probe.tcl` — is the contact command surface registered at all? This is what
  caught P0.a (`OpenSeesMP` linked `commands.cpp`, which had no contact verbs).
* `query.tcl` + `serial_model.tcl` — smoke the eight query/augment commands the
  P0 deck never calls, against a live contact.
* `model.py` / `run_serial.py` / `run_mp.py` / `run_exp.py` / `explicit.py` —
  the Python-lane twins (`opensees` / `openseesmp` modules), including the
  explicit lane and the `P0_MUTATE` / `P0_GHOSTMASS` mutations.
* `pin.py` — **read this before trusting any Python result.** A bare
  `import opensees` in the dev venv resolves through `ladruno_opensees.pth`,
  which on this machine pointed at a build in an unrelated session's scratchpad.
  Serial and MP must come from the same build or the comparison is meaningless.
  Pin with `LADRUNO_OPENSEES_BIN` / `LADRUNO_OPENSEESMP_BIN` **before** Python
  starts and assert `mod.__file__`.

## Running

Tcl (adjust the build dir):

```
set TCL_LIBRARY=<tree>\dist\lib\tcl8.6
<build>\OpenSees.exe serial.tcl
mpiexec -n 2 <build>\OpenSeesMP.exe mp.tcl
mpiexec -n 2 <build>\OpenSeesMP.exe mp_noghost.tcl      # must abort, not hang
```

## Two traps that cost real time

* **A hung MPI job leaves orphaned ranks alive.** They hold the binary open and
  the *next* link fails with `LNK1104: cannot open file 'OpenSeesMP.exe'`. If a
  link fails that way, look for stray `OpenSeesMP.exe` / `hydra_pmi_proxy`
  processes before touching the build system. P1's teardown exists so this stops
  happening, but a pre-P1 binary can still do it.
* **A green run proves nothing until a broken one fails.** The first
  `MPI_Abort` fix compiled cleanly, left both happy paths bit-identical, and was
  entirely inert — the `#ifdef` was in a TU compiled once with no parallel
  define. Only `mp_noghost.tcl` caught it, by hanging exactly as before. Run the
  mutation every time.
