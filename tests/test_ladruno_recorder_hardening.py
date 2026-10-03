"""Ladruno recorder hardening (WP-163): R1 eigen gate, R2 failed rebuild, R5 silent loss.

Found by the 2026-10-03 red/blue review
(Ladruno_implementation/163_ladruno_recorder_hardening_roadmap.md):

  * R1 -- `eigen` -> `wipe` -> new model with `-N modesOfVibration` -> `analyze`
    killed the process: `wipe` never resets the interpreter's numEigen, the recorder
    gated its modal path on it, and Domain::getEigenvalues() exit(-1)'d on the fresh
    domain. The gate now asks the domain (getNumEigenvalues, ADR46).
  * R2 -- a MODEL_STAGE rebuild that failed part-way (here: every node of the `-R`
    region removed -> "no nodes to write") returned before the old sources were
    released, and the stamp was already committed, so the next step called
    getResponse() on Response objects of DELETED elements. Sources are now released
    first and recording is suspended (one error) until the next domain change.
  * R5 -- HDF5 create/append failures were ignored: a result could vanish from the
    file with only HDF5-DIAG noise. A duplicate request (`-N displacement
    displacement`) is the reproducible case: the second sink's create failed, it
    still marked itself initialized, then appended into the FIRST sink's group,
    doubling its rows. It is now refused once, by name, and DATA keeps T rows.

R1 and R2 kill or corrupt the interpreter on the unfixed binary, so they run in a
SUBPROCESS and assert a marker printed after the last command.
"""
import os
import subprocess
import sys
import textwrap

import pytest

from _testbed import ops

h5py = pytest.importorskip("h5py")
pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

DONE = "WP163_DONE"


def _run_child(body, tmp_path):
    """Run `body` (uses `ops`) in a fresh interpreter on THIS build of opensees."""
    moddir = os.path.dirname(os.path.abspath(ops.__file__))
    script = "\n".join([
        "import os, sys",
        f"sys.path.insert(0, {moddir!r})",
        f"os.add_dll_directory({moddir!r})" if hasattr(os, "add_dll_directory") else "",
        "import opensees as ops",
        textwrap.dedent(body),
        f"print({DONE!r}, flush=True)",
    ])
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1", HDF5_USE_FILE_LOCKING="FALSE")
    # -S: a boot .pth could pre-import a DIFFERENT opensees build
    # (LEDGER_quirks "An INSTALLED Ladruno hijacks `import opensees`").
    return subprocess.run([sys.executable, "-S", "-c", script], env=env,
                          capture_output=True, text=True, timeout=180,
                          cwd=str(tmp_path))


# Two x-direction trusses with lumped x-masses: 2 free DOFs, a well-posed
# generalized eigenproblem for -fullGenLapack.
TRUSS_MODEL = """
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1); ops.mass(2, 1.0, 0.0)
ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1); ops.mass(3, 1.0, 0.0)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)
"""

STATIC = """
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0); ops.analysis("Static")
"""


def test_r1_modes_after_wipe_do_not_kill_the_process(tmp_path):
    out = str(tmp_path / "r1.ladruno")
    body = TRUSS_MODEL + """
ops.eigen("-fullGenLapack", 1)       # leaves the interpreter's numEigen = 1
""" + TRUSS_MODEL + f"""
ops.recorder("ladruno", {out!r}, "-N", "displacement", "modesOfVibration")
""" + STATIC + """
ops.analyze(1)                       # fresh domain: no spectrum -> no modal write
ops.wipeAnalysis()
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.eigen("-fullGenLapack", 1)       # a real spectrum now exists
ops.record()                         # -> modes are written
ops.wipe()
"""
    r = _run_child(body, tmp_path)
    log = r.stdout + r.stderr
    assert DONE in r.stdout, f"process died (rc={r.returncode}):\n{log}"
    assert "Eigenvalues were never set" not in log, log
    with h5py.File(out, "r") as f:
        modes = []
        f.visit(lambda name: modes.append(name)
                if name.rsplit("/", 1)[-1].startswith("MODE_") else None)
    assert modes, "modesOfVibration wrote no MODE_<k> dataset after the real eigen"


def test_r2_failed_stage_rebuild_suspends_instead_of_dangling(tmp_path):
    out = str(tmp_path / "r2.ladruno")
    body = f"""
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
# region 1: two trusses (nodes 1-2-3); a separate truss 10-11 keeps the model alive
ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1)
ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1)
ops.node(10, 0.0, 5.0); ops.fix(10, 1, 1)
ops.node(11, 1.0, 5.0); ops.fix(11, 0, 1)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)
ops.element("Truss", 3, 10, 11, 1.0, 1)
ops.region(1, "-ele", 1, 2)
ops.recorder("ladruno", {out!r}, "-R", 1, "-N", "displacement", "-E", "force")
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(11, 1.0, 0.0)
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
ops.integrator("LoadControl", 0.1); ops.analysis("Static")
ops.analyze(2)
# remove the whole region: its elements, its nodes' SPs, its nodes
ops.remove("element", 1); ops.remove("element", 2)
ops.remove("sp", 1, 1); ops.remove("sp", 1, 2)
ops.remove("sp", 2, 2); ops.remove("sp", 3, 2)
ops.remove("node", 1); ops.remove("node", 2); ops.remove("node", 3)
ops.analyze(3)                       # rebuild fails once, then 2 more commits
ops.wipe()
"""
    r = _run_child(body, tmp_path)
    log = r.stdout + r.stderr
    assert DONE in r.stdout, f"process died (rc={r.returncode}):\n{log}"
    assert log.count("could not be written; recording is suspended") == 1, log
    # the first stage is intact: 2 rows of displacement and of element force
    with h5py.File(out, "r") as f:
        stages = sorted(k for k in f if k.startswith("MODEL_STAGE"))
        first = f[stages[0]]
        assert first["RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0] == 2


def test_r5_duplicate_request_is_refused_not_doubled(tmp_path, capfd):
    out = str(tmp_path / "r5.ladruno")
    exec(TRUSS_MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "displacement")
    exec(STATIC, {"ops": ops})
    ops.analyze(3)
    ops.wipe()
    log = "".join(capfd.readouterr())
    assert "already exists in this MODEL_STAGE" in log, log
    assert log.count("already exists in this MODEL_STAGE") == 1, log   # once, not per step
    with h5py.File(out, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
        g = f[f"{stage}/RESULTS/ON_NODES/DISPLACEMENT"]
        assert g["DATA"].shape[0] == 3, g["DATA"].shape          # was 6 (2 sinks)
        assert g["TIME"].shape[0] == 3 and g["STEP"].shape[0] == 3
