# Running parallel OpenSees in Python — User Guide

*For people who installed **Ladruno OpenSees** with the setup wizard. No
building, no developer tools, no Intel software required — everything you
need was installed for you.*

---

## 1. What you have and where it is

The installer placed everything under one folder (your **install folder**).
By default that is:

- **`C:\Program Files\Ladruno\OpenSees`** — if you installed for all users, or
- **`%LOCALAPPDATA%\Programs\Ladruno\OpenSees`** — if you installed just for you.

Not sure which? Open **Settings ▸ Apps**, find **“Ladruno OpenSees”**, or
re-open the installer — it shows the path. In the examples below, replace
`<install>` with your actual folder.

Inside `<install>` the parts you care about:

| Folder | What it is |
|---|---|
| `<install>\opensees_venv\` | The ready-to-use Python environment the wizard created and wired for you (only if you chose “create a venv”). |
| `<install>\openseesmp\` | The **parallel** Python module **and its own `mpiexec.exe`** — always launch parallel jobs with this `mpiexec`. |
| `<install>\bin\` | Sequential `opensees` for Python + the command-line `OpenSees.exe`. |

If during setup you pointed at **your own** virtual environment instead of
letting the wizard create one, use that venv’s folder wherever this guide
says `<install>\opensees_venv`.

---

## 2. 30-second check that it works

Open a terminal (Command Prompt or PowerShell) and run **one line** (adjust
`<install>`):

```bat
"<install>\openseesmp\mpiexec.exe" -n 4 "<install>\opensees_venv\Scripts\python.exe" -c "import openseesmp as o; print('rank', o.getPID(), 'of', o.getNP())"
```

You should see **four lines**, one per rank, like:

```
rank 0 of 4
rank 1 of 4
rank 2 of 4
rank 3 of 4
```

Four *different* rank numbers = parallel OpenSees is working. (If every line
says `rank 0 of 1`, see Troubleshooting.)

---

## 3. How it works (the one idea to understand)

When you launch with `mpiexec -n N`, you get **N independent copies** of
OpenSees running at once. **Each copy runs your whole script**, and each has
its **own separate model**.

Your script asks *“which copy am I?”* and does a different piece of the work:

```python
import openseesmp as ops
me   = ops.getPID()   # 0, 1, 2, … (this copy's number)
total = ops.getNP()   # how many copies in total (the N you launched with)
```

So `openseesmp` is perfect when you have **many independent runs** — e.g. the
same structure under 50 ground motions, or a design swept over 200 parameter
values. Each rank just handles a share of the list.

> It does **not** make one giant model solve faster by splitting it across
> cores. It runs **many models at the same time**. (For one model, the normal
> `import opensees` is what you want.)

### Defining nodes & elements: which rank “owns” what?

This is the question people coming from parallel FEM ask first, so to be
crystal clear:

**You do *not* assign a structure’s nodes/elements to ranks. There is no
shared model.** Every `node`, `element`, `fix`, `load`, … command you call
goes into **this rank’s own private model only**. Tags are *local*: rank 0’s
`node 1` and rank 1’s `node 1` are two unrelated nodes in two unrelated
models. Nothing is connected, summed, or solved jointly across ranks. There
is no global stiffness matrix spanning ranks.

So the rank number chooses **which complete model this copy builds** — never
**which piece of one structure** it builds.

❌ **Wrong mental model** (this does *not* work here):

```python
# DON'T: trying to split ONE frame across ranks
me = ops.getPID()
if me == 0:
    # nodes/elements 1..100  (left half of the building)
elif me == 1:
    # nodes/elements 101..200 (right half)
# ...expecting OpenSees to couple the halves and solve one frame -> NO.
# Each rank just gets a disconnected half-model. Nothing couples them.
```

✅ **Right mental model** — full model on every rank, rank picks the *case*:

```python
me, total = ops.getPID(), ops.getNP()
motions = ['gm01', 'gm02', 'gm03', 'gm04', 'gm05', 'gm06', 'gm07', 'gm08']
for i in range(me, len(motions), total):
    gm = motions[i]
    ops.wipe()
    # ---- build the ENTIRE structure here (all nodes, all elements) ----
    # ---- apply ground motion `gm`, run, write result_<gm>_rank<me>.txt
```

Rule of thumb: **build the whole structure inside every rank**, exactly as you
would for a single serial run; the only thing that changes between ranks is
the *input* (which load case / ground motion / parameter set), selected with
`ops.getPID()`.

> **What if I really need one huge model split across machines/cores?**
> That is *domain decomposition* (each rank owns a subdomain of one coupled
> model). It is **not available from Python** in this build — the `partition`
> command is intentionally inert here. Use the Tcl **`OpenSeesSP.exe`**
> (in `<install>\bin\`), which auto-partitions one model across ranks, or
> re-cast the problem as many independent runs as above.

---

## 4. Your first parallel script

Save this as `hello_parallel.py`:

```python
import openseesmp as ops

me, total = ops.getPID(), ops.getNP()

# A 1-DOF spring whose stiffness depends on which rank we are.
k = float(me + 1)
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
ops.node(1, 0.0); ops.node(2, 1.0)
ops.fix(1, 1)
ops.uniaxialMaterial('Elastic', 1, k)
ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1); ops.load(2, 10.0)
ops.system('ProfileSPD'); ops.numberer('Plain'); ops.constraints('Plain')
ops.integrator('LoadControl', 1.0); ops.algorithm('Linear')
ops.analysis('Static'); ops.analyze(1)

print(f"rank {me}/{total}:  k={k}  ux={ops.nodeDisp(2, 1):.4f}", flush=True)
```

Run it (replace `<install>`):

```bat
"<install>\openseesmp\mpiexec.exe" -n 4 "<install>\opensees_venv\Scripts\python.exe" hello_parallel.py
```

Each rank solves its *own* spring and prints a different answer
(`ux = 10/k`).

---

## 5. A real example — run many cases at once

This is the pattern you’ll actually use: a list of cases, split across ranks.

```python
# sweep.py
import openseesmp as ops

me, total = ops.getPID(), ops.getNP()

cases = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]   # e.g. scale factors

# This rank handles every Nth case starting at its own number:
for i in range(me, len(cases), total):
    factor = cases[i]
    ops.wipe()
    # ----- build your model here, using `factor` -----
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 1.0); ops.fix(1, 1)
    ops.uniaxialMaterial('Elastic', 1, 1000.0)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1); ops.load(2, 10.0 * factor)
    ops.system('ProfileSPD'); ops.numberer('Plain'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear')
    ops.analysis('Static'); ops.analyze(1)

    # IMPORTANT: put the rank number in every output filename so the
    # parallel copies never overwrite each other.
    with open(f"result_case{i}_rank{me}.txt", "w") as fh:
        fh.write(f"factor={factor}  ux={ops.nodeDisp(2, 1):.6f}\n")

print(f"rank {me}/{total}: done", flush=True)
```

```bat
"<install>\openseesmp\mpiexec.exe" -n 8 "<install>\opensees_venv\Scripts\python.exe" sweep.py
```

Eight cases run simultaneously and write eight result files.

---

## 6. Getting your results back

The simple, reliable way (recommended): **each rank writes its own files**
(always with the rank number in the name, e.g. `..._rank0.txt`), and you
read them all back with a normal, non-parallel script afterwards:

```python
# collect.py  — run normally:  <venv>\Scripts\python.exe collect.py
import glob
for f in sorted(glob.glob("result_case*_rank*.txt")):
    print(f, "->", open(f).read().strip())
```

(Advanced users can exchange data between ranks live with
`ops.barrier()`, `ops.send(...)`, `ops.recv(...)`, `ops.Bcast(...)` —
but the write-files-then-collect approach above is simplest and is the
recommended starting point.)

---

## 7. Command cheat-sheet

```bat
REM Always: <bundled mpiexec>  -n <copies>  <wired venv python>  <your script>

"<install>\openseesmp\mpiexec.exe" -n 8 "<install>\opensees_venv\Scripts\python.exe" myjob.py
```

- **`-n 8`** = run 8 parallel copies. A good starting value is the number of
  physical CPU cores on your machine.
- Quieter output (hide the banner on every rank): add
  `-genv LADRUNO_OPENSEES_QUIET 1` right after `mpiexec.exe`.
- For **non-parallel** work, just activate the venv and use `import opensees`
  the normal way — no `mpiexec` needed.

Activate the venv for normal interactive use:

```bat
"<install>\opensees_venv\Scripts\activate.bat"      REM Command Prompt
```
```powershell
& "<install>\opensees_venv\Scripts\Activate.ps1"    # PowerShell
```

---

## 8. Troubleshooting

| What you see | What to do |
|---|---|
| Every line says **`rank 0 of 1`** | You ran plain `python myjob.py` (one copy). You must launch through the bundled **`mpiexec -n N`** as shown above. |
| `ImportError` / `DLL load failed` on `import openseesmp` | Use the venv the wizard wired (`<install>\opensees_venv\Scripts\python.exe`), not some other Python. The module needs **Python 3.12**. |
| `mpiexec.exe` “is not recognized” | Use the **full path** to the bundled one: `"<install>\openseesmp\mpiexec.exe"`. Don’t use an `mpiexec` from any other MPI install. |
| A fatal `MPI_Init` / `libfabric` error | You launched with the wrong `mpiexec`. Always use `<install>\openseesmp\mpiexec.exe` together with the wired venv’s `python.exe`. |
| Results files look mixed up / overwritten | Put `ops.getPID()` (the rank number) into **every** output filename and recorder path. |
| “I want to split one big model across cores” | That’s a different mode and isn’t available from Python here. Either run many independent models (what this guide shows), or ask for the Tcl `OpenSeesSP.exe` workflow. |
| It worked once, then I moved/renamed the install folder | The venv was wired to the original path. Re-run the installer (or re-wire the venv) after moving things. |

---

## 9. Quick recap

1. Write a normal OpenSees-Python script. Add at the top:
   `import openseesmp as ops` and read `ops.getPID()` / `ops.getNP()`.
2. Use the rank number to pick which slice of your work this copy does.
3. Put the rank number in every output filename.
4. Launch with the **bundled** `mpiexec` and the **wired venv** python:
   `"<install>\openseesmp\mpiexec.exe" -n N "<install>\opensees_venv\Scripts\python.exe" myjob.py`
5. Collect the per-rank result files with a normal script afterward.

That’s the whole workflow. For ordinary single runs, keep using
`import opensees` as usual — both live in the same environment.
