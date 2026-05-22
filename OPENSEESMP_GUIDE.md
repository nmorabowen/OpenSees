# Using `openseesmp` — MPI-parallel OpenSees in Python

> The pip world calls this "openseespymp". **In this Ladruno build the import
> name is `openseesmp`** (it coexists with the sequential `import opensees`
> and `import openseespy.opensees` in the same venv).

---

## 1. What it is — and what it is *not*

`openseesmp` is the **`_PARALLEL_INTERPRETERS`** model:

- You launch **N independent OpenSees interpreters** with `mpiexec -n N`.
- **Every rank runs your entire Python script**, each with its **own separate
  model/domain**.
- Ranks tell themselves apart with **`ops.getPID()`** (0 … N-1) and
  **`ops.getNP()`** (= N).

This is ideal for **embarrassingly-parallel structural workloads**: parametric
sweeps, IDA, Monte-Carlo / sensitivity, running many ground motions or many
designs at once — *one independent analysis per rank*.

> **It is NOT domain decomposition.** `openseesmp` does **not** split one big
> model across ranks. A single huge model partitioned across processes is the
> *SP* paradigm, which is **not available in Python** in this build — use the
> Tcl `OpenSeesSP.exe` for that. (`ops.partition()` is intentionally a no-op /
> error here.)

---

## 2. Requirements

| Requirement | Why |
|---|---|
| **Python 3.12** | `openseesmp.pyd` is built against the CPython 3.12 ABI. Any other Python version fails to import. |
| **The bundled `mpiexec`** | Ships in `…\openseesmp\mpiexec.exe`. It is the exact Intel MPI version OpenSees + MUMPS were built against. Using a *different* `mpiexec`/`impi.dll` causes MUMPS *“Instance Error 1”* or an `MPI_Init` abort. |
| **No oneAPI install needed** | The `…\openseesmp\` folder is self-contained: `impi.dll`, `libfabric.dll`, the Hydra launchers and the MKL runtime are all bundled. |

---

## 3. Quick start

### A. If you installed via the Ladruno installer (recommended)

Pick **“create a new venv”** (or point at an existing **Python 3.12** venv) in
the wizard. The installer wires that venv so `import openseesmp` just works —
no `PATH`/`sys.path` fiddling.

```bat
REM <install> = e.g. C:\Program Files\Ladruno\OpenSees  (or your per-user dir)
REM <venv>    = the venv the installer created/wired, e.g. <install>\opensees_venv

"<install>\openseesmp\mpiexec.exe" -n 4 "<venv>\Scripts\python.exe" driver.py
```

That's it. `driver.py` does `import openseesmp as ops`.

### B. Straight from a dev build (no installer)

```powershell
$mp  = 'C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\openseesmp'
$py  = 'C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe'
$env:PATH = "$mp;$env:PATH"          # so openseesmp.pyd finds impi.dll / libfabric / MKL
& "$mp\mpiexec.exe" -n 4 $py "$mp\test_openseesmp.py"
```

`dist\openseesmp\test_openseesmp.py` is a ready-made smoke test — run it first
to confirm everything works (it prints a distinct `[rank k/4]` line per rank).

> Tip: `set LADRUNO_OPENSEES_QUIET=1` suppresses the per-rank ASCII banner so
> parallel output stays clean.

---

## 4. The programming model

```python
import openseesmp as ops

pid = ops.getPID()      # this rank's id: 0 .. nprocs-1
nprocs = ops.getNP()    # total number of ranks (the -n value)
```

Every rank executes the whole script. You make ranks do *different* work by
slicing your job list with `pid`/`nprocs`:

```python
jobs = [...]                                  # your full list of cases
for i in range(pid, len(jobs), nprocs):       # rank-strided split
    case = jobs[i]
    ...build & analyze a model for `case`...
```

**Golden rule:** anything written to disk (recorders, result files, logs) must
include `pid` in its name so ranks never overwrite each other.

---

## 5. Worked example — parametric sweep, one case per rank

```python
# driver.py   ->   mpiexec -n 8 python driver.py
import openseesmp as ops

pid, nprocs = ops.getPID(), ops.getNP()
stiffnesses = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

for i in range(pid, len(stiffnesses), nprocs):
    k = stiffnesses[i]
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1); ops.load(2, 10.0)
    ops.system('ProfileSPD')          # a *sequential* solver: each rank
    ops.numberer('Plain')             # solves its own independent model
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

    with open(f'out_case{i}_rank{pid}.txt', 'w') as f:   # pid in the name!
        f.write(f'k={k}\tux={ops.nodeDisp(2, 1):.6f}\n')

print(f'[rank {pid}/{nprocs}] finished', flush=True)
```

`mpiexec -n 8 python driver.py` → 8 cases solved concurrently, 8 result files,
then post-process them in plain (serial) Python afterwards.

---

## 6. A more realistic pattern — a ground-motion / IDA suite

```python
import os, glob, openseesmp as ops

pid, nprocs = ops.getPID(), ops.getNP()
motions = sorted(glob.glob(r'ground_motions\*.AT2'))   # same list on every rank
scales  = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
cases   = [(m, s) for m in motions for s in scales]    # full (gm, scale) grid

os.makedirs('results', exist_ok=True)
for idx in range(pid, len(cases), nprocs):
    gm, sf = cases[idx]
    tag = f'{os.path.splitext(os.path.basename(gm))[0]}_sf{sf:g}'
    ops.wipe()
    # ... build your structure ...
    # ... define the GM with timeSeries/pattern using `gm`, scaled by `sf` ...
    ops.recorder('Node', '-file', f'results/roof_{tag}_r{pid}.out',
                 '-time', '-node', 999, '-dof', 1, 'disp')
    # ... transient analysis loop ...
    print(f'[rank {pid}/{nprocs}] done {tag}', flush=True)
```

Run: `…\openseesmp\mpiexec.exe -n 16 <python> ida_driver.py`. Collapse-check /
fragility post-processing reads the `results/` files serially afterwards.

---

## 7. Collecting results across ranks

Two approaches, in order of recommendation:

1. **File-per-rank, then post-process (recommended).** Simple, robust,
   restart-friendly. Every rank writes `…_r{pid}…` files; a separate serial
   script aggregates them when the run is done. This is the workflow shown
   above and the one to start with.

2. **In-run message passing (advanced).** The same MPI commands as Tcl
   `OpenSeesMP` are available:
   `ops.getPID()`, `ops.getNP()`, `ops.barrier()`,
   `ops.send(...)`, `ops.recv(...)`, `ops.Bcast(...)`.
   `getPID`/`getNP` are verified here; if you use `send`/`recv`/`Bcast`,
   sanity-check the exact argument form with a 2-rank toy first (they mirror
   the Tcl OpenSeesMP commands). Use this only when ranks must exchange data
   *during* the run.

---

## 8. Launch command reference

```text
Installed (wired venv):
  "<install>\openseesmp\mpiexec.exe" -n N "<venv>\Scripts\python.exe" driver.py

Dev build:
  set PATH=...\dist\openseesmp;%PATH%
  ...\dist\openseesmp\mpiexec.exe -n N <python3.12>.exe driver.py

  N            number of ranks (processes). Start with N = physical cores.
  driver.py    your script; does `import openseesmp as ops`
```

`mpiexec` extras you may want: `-genv LADRUNO_OPENSEES_QUIET 1` (quiet banner on
all ranks), `-outfile-pattern out.%r.log` / `-errfile-pattern err.%r.log`
(per-rank stdout/stderr capture).

---

## 9. Gotchas & troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ImportError: DLL load failed` on `import openseesmp` | Wrong Python (must be **3.12**), or the `openseesmp` folder isn't reachable. Use the installer-wired venv, or prepend `…\openseesmp` to `PATH`. |
| `Fatal error in PMPI_Init … MPIDI_OFI_mpi_init_hook` (abort 1090959) | `libfabric.dll` not found. Launch with the **bundled** `…\openseesmp\mpiexec.exe` and make sure `…\openseesmp` is on the process `PATH` (the installer-wired venv does this automatically). |
| MUMPS *“Instance Error 1”* under `-n > 1` | A mismatched `impi.dll`/`mpiexec` (different Intel MPI than OpenSees was built with). Always use the bundled `mpiexec`. |
| Every rank prints `PID 0 / NP 1` | You ran `python driver.py` **without** `mpiexec` (that's a single rank). Launch via the bundled `mpiexec -n N`. |
| Ranks overwrite each other's output | Put `ops.getPID()` in every filename / recorder path. |
| Banner spam (N copies) | `set LADRUNO_OPENSEES_QUIET=1` (or `mpiexec -genv LADRUNO_OPENSEES_QUIET 1`). |
| Want to split one big model across ranks | Not supported in Python (that's SP). Use the Tcl `OpenSeesSP.exe`, or restructure as independent per-rank models. |
| `import opensees` **and** `import openseesmp` in one script | Don't. Use `opensees` for serial scripts and `openseesmp` for `mpiexec` scripts — one OpenSees runtime per process. |

---

## 10. One-page reference

```python
import openseesmp as ops
r, n = ops.getPID(), ops.getNP()         # rank, world size
for i in range(r, len(jobs), n):         # rank-strided work split
    ops.wipe()
    ...                                   # build + analyze job i
    write(f'result_{i}_r{r}.dat', ...)    # pid in EVERY output name
ops.barrier()                             # optional sync point
```

```bat
"<…>\openseesmp\mpiexec.exe" -n 8 "<venv>\Scripts\python.exe" driver.py
```

- Independent interpreter per rank — **not** domain decomposition.
- Python **3.12** only · use the **bundled `mpiexec`** · `pid` in all filenames.
- Serial / non-MPI work → plain `import opensees`.
- Smoke test: `dist\openseesmp\test_openseesmp.py`.

---
*Ladruno OpenSees — `openseesmp` is the `_PARALLEL_INTERPRETERS` Python module
(Patch 9). For the build/architecture rationale see
`Ladruno_implementation/02_openseespymp.md`.*
