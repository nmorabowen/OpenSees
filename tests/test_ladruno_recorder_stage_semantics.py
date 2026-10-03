"""Ladruno recorder stage semantics (WP-165).

Found by the 2026-10-03 red/blue review (Ladruno_implementation/163_...roadmap.md):

  * R6 -- every move of the domain-change stamp started a new MODEL_STAGE (full
    model copy, envelope + energy reset), but the stamp moves for much more than
    topology: a pattern holding SPs, `eleLoad`, contact `-reemit` re-sorts. Now a
    stamp move rebuilds only when the node/element/pressure-constraint SET changed.
  * R3 -- the energy integrals restarted at every MODEL_STAGE with a rate x t
    jump, and were integrated only over the `-T` samples. They now live in the
    recorder and are integrated on every commit.
  * MP-9 -- INFO carries RUN_ID (+ RUN_ID_SCOPE) so a reader can reject stale
    part files from an earlier run.
"""
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from _testbed import ops

h5py = pytest.importorskip("h5py")
pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def _truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 1.0e3)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1); ops.mass(2, 1.0, 0.0)
    ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1); ops.mass(3, 1.0, 0.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 2, 3, 1.0, 1)


def _stages(path):
    with h5py.File(path, "r") as f:
        return sorted(k for k in f if k.startswith("MODEL_STAGE"))


def _static(dlam):
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("LoadControl", dlam); ops.analysis("Static")


def test_r6_sp_pattern_does_not_start_a_new_stage(tmp_path):
    """A new pattern with an imposed displacement bumps the stamp (Domain::addSP_
    Constraint) but changes no node/element: one MODEL_STAGE, rows continue."""
    out = str(tmp_path / "r6.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.pattern("Plain", 2, 1)
    ops.sp(3, 1, 0.01)
    assert ops.analyze(2) == 0
    ops.wipe()
    stages = _stages(out)
    assert len(stages) == 1, stages
    with h5py.File(out, "r") as f:
        assert f[f"{stages[0]}/RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0] == 4


def test_r6_topology_change_still_starts_a_new_stage(tmp_path):
    out = str(tmp_path / "r6b.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.node(99, 5.0, 5.0); ops.fix(99, 1, 1)          # a new node: topology changed
    assert ops.analyze(2) == 0
    ops.wipe()
    assert len(_stages(out)) == 2, _stages(out)


def _dynamic_energy(out, extra, steps_a=40, steps_b=40, add_node=False):
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement", "-G", "energy", *extra)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    assert ops.analyze(steps_a, 0.01) == 0
    if add_node:
        ops.node(99, 5.0, 5.0); ops.fix(99, 1, 1)
    assert ops.analyze(steps_b, 0.01) == 0
    ops.wipe()
    rows = []
    with h5py.File(out, "r") as f:
        for st in sorted((k for k in f if k.startswith("MODEL_STAGE")),
                         key=lambda k: int(k[len("MODEL_STAGE["):-1])):
            g = f[f"{st}/RESULTS/ON_DOMAIN/energyBalance"]
            rows.append((g["TIME"][...], g["DATA"][:, 0, :]))   # [T x 6] KE,IE,DW,ULW,RES,ERR
    return rows


def test_r3_energy_continues_across_a_stage_change(tmp_path):
    rows = _dynamic_energy(str(tmp_path / "r3.ladruno"), [], add_node=True)
    assert len(rows) == 2, len(rows)
    (t1, e1), (t2, e2) = rows
    ie_last, ie_next = e1[-1, 1], e2[0, 1]
    ulw_last, ulw_next = e1[-1, 3], e2[0, 3]
    # the next commit adds ONE step of work; it used to restart near rate x t_total
    step_ulw = abs(e1[-1, 3] - e1[-2, 3])
    assert abs(ulw_next - ulw_last) < 5.0 * step_ulw + 1e-12, (ulw_last, ulw_next, step_ulw)
    assert ie_next > 0.5 * ie_last, (ie_last, ie_next)
    # and the closure still holds: |RES| small next to the external work
    assert abs(e2[-1, 4]) < 0.05 * abs(e2[-1, 3]) + 1e-12, e2[-1]


def test_r3_energy_does_not_depend_on_T_sampling(tmp_path):
    every = _dynamic_energy(str(tmp_path / "e1.ladruno"), ["-T", "nsteps", 1])
    sparse = _dynamic_energy(str(tmp_path / "e5.ladruno"), ["-T", "nsteps", 5])
    (t_a, e_a), (t_b, e_b) = every[0], sparse[0]
    # compare at the common sample times
    common = np.intersect1d(np.round(t_a, 9), np.round(t_b, 9))
    assert len(common) >= 10, common
    ia = np.isin(np.round(t_a, 9), common)
    ib = np.isin(np.round(t_b, 9), common)
    np.testing.assert_allclose(e_b[ib][:, 1:4], e_a[ia][:, 1:4], rtol=1e-9, atol=1e-12)


def test_mp9_info_carries_a_run_id(tmp_path):
    out = str(tmp_path / "mp9.ladruno")
    moddir = os.path.dirname(os.path.abspath(ops.__file__))
    script = "\n".join([
        "import os, sys",
        f"sys.path.insert(0, {moddir!r})",
        f"os.add_dll_directory({moddir!r})" if hasattr(os, "add_dll_directory") else "",
        "import opensees as ops",
        textwrap.dedent(f"""
            ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
            ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0); ops.fix(1, 1, 1); ops.fix(2, 0, 1)
            ops.uniaxialMaterial("Elastic", 1, 1.0e3); ops.element("Truss", 1, 1, 2, 1.0, 1)
            ops.recorder("ladruno", {out!r}, "-N", "displacement")
            ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 1.0, 0.0)
            ops.system("FullGeneral"); ops.integrator("LoadControl", 1.0)
            ops.algorithm("Linear"); ops.analysis("Static"); ops.analyze(1); ops.wipe()
        """),
    ])
    launcher = ("SLURM_JOB_ID", "SLURM_STEP_ID", "OMPI_MCA_ess_base_jobid",
                "PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK",
                "SLURM_NTASKS", "SLURM_PROCID")
    env = {k: v for k, v in os.environ.items() if k not in launcher}
    env.update(LADRUNO_RUN_ID="run-abc", LADRUNO_OPENSEES_QUIET="1",
               HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([sys.executable, "-S", "-c", script], env=env,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr

    def s(v):
        v = v.flat[0] if hasattr(v, "flat") else v
        return v.decode() if isinstance(v, bytes) else str(v)

    with h5py.File(out, "r") as f:
        assert s(f["INFO"].attrs["RUN_ID"]) == "run-abc"
        assert s(f["INFO"].attrs["RUN_ID_SCOPE"]) == "user"
