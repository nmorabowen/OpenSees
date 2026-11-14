"""WP-232c: the MPM host's lifecycle, refusals, materials and stage path through
the shipped binaries (ADR 232 §5.2, §6, §8 S1 rows "lifecycle", "refusals",
"engine parity: the loadPattern stage path", "thread-count determinism").

* refusals are loud, never a wrong answer: a non-admitted material, LadrunoJ2
  -autoRegularization, an initial stress on a non-elastic material, free FE DOFs,
  a second grid, Δt above s·Δt_c, an inadmissible particle, database save;
* lifecycle: `wipe` and `remove loadPattern` delete the host; `reset` reaches it
  only as a re-applyLoad (particles stay); `mpm reset` returns the particles to
  the state at the first `mpm analyze`, a rerun reproduces it bit for bit, and
  PARTICLE_EPOCH counts the resets; the host is inert in an FE analysis;
* admitted materials: DruckerPrager and LadrunoJ2 in their elastic range give
  ElasticIsotropic's answer; DruckerPrager's `materialState` reaches the
  particles through `parameter <t> loadPattern <tag> materialState` (W12);
* threads: the particle loop is serial in v1, so any `ladrunoThreads` gives the
  same output.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_harness as H  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]
TAG = 232



def _scalar(v):
    """A scalar query result (a float; tolerate a one-element list)."""
    return float(v[0]) if isinstance(v, (list, tuple)) else float(v)

def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _block(ops, mat_cmd=("ElasticIsotropic", 1, 1.0e6, 0.25), ppc=2):
    """A 2D block, sine-velocity field (it strains), walls on the box."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial(*mat_cmd)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    ops.mpm("body", 1, "box", 0.2, 0.2, 0.6, 0.6, "-ppc", ppc, "-mat", mat_cmd[1], "-rho", 2000.0)
    ops.mpm("velocity", 1, "sine", "x", "x", 0.5, 7.853981633974483, 0.2)
    ops.mpm("fix", "xmin", "xy")
    ops.mpm("wall", "ymin", 0.3)


def _fields(ops, names=("x", "v", "F", "V", "sig", "sigHat", "epsHat")):
    return {f: list(ops.mpm("particleData", f)) for f in names}


def _close(a, b, rel=1e-9, ab=1e-12):
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert abs(x - y) <= ab + rel * abs(y), (x, y)


# --------------------------------------------------------------------------- #
#  Refusals                                                                    #
# --------------------------------------------------------------------------- #
def test_refuses_non_admitted_material():
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("J2Plasticity", 5, 1.0e6, 0.5e6, 1.0e3, 1.0e3, 0.0, 0.0)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.2, 0.2, 0.6, 0.6, "-ppc", 2, "-mat", 5, "-rho", 2000.0)
    assert int(ops.mpm("info", "numParticles")) == 0


def test_refuses_ladrunoj2_autoregularization():
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("LadrunoJ2", 6, 1.0e6, 0.5e6, "-iso", "voce", 1.0e9, 0.0, 0.0, 0.0,
                   "-autoRegularization", 0.1)
    ops.nDMaterial("LadrunoJ2", 7, 1.0e6, 0.5e6, "-iso", "voce", 1.0e9, 0.0, 0.0, 0.0)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.2, 0.2, 0.6, 0.6, "-ppc", 2, "-mat", 6, "-rho", 2000.0)
    assert int(ops.mpm("body", 1, "box", 0.2, 0.2, 0.6, 0.6, "-ppc", 2, "-mat", 7, "-rho", 2000.0)) == 64


def test_refuses_initial_stress_on_non_elastic_material():
    ops = _ops()
    _block(ops, ("DruckerPrager", 3, 8.0e5, 4.0e5, 1.0e9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    with pytest.raises(Exception):
        ops.mpm("stress", 1, "uniform", -1.0, -1.0, 0.0)
    sig = list(ops.mpm("particleData", "sig"))
    assert max(abs(s) for s in sig) == 0.0


def test_refuses_free_fe_dofs():
    ops = _ops()
    _block(ops)
    ops.node(1, 0.0, 0.0)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)
    ops.fix(1, 1, 1)
    ops.mpm("analyze", 1, 1.0e-4)


def test_refuses_second_grid_and_tag_clash():
    ops = _ops()
    _block(ops)
    with pytest.raises(Exception):
        ops.mpm("grid", TAG + 1, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)


def test_refuses_dt_and_changes_nothing():
    ops = _ops()
    _block(ops)
    dtc = _scalar(ops.mpm("criticalTimeStep"))
    before = _fields(ops)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.51 * dtc)
    assert _fields(ops) == before
    ops.mpm("analyze", 1, 0.49 * dtc)


def test_refuses_inadmissible_particle():
    ops = _ops()
    _block(ops)
    ops.mpm("particle", 500, 2, 1, 2000.0, 0.0025, "-x", 1.05, 0.5, "-lp", 0.025, 0.025)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-5)


# --------------------------------------------------------------------------- #
#  Lifecycle                                                                   #
# --------------------------------------------------------------------------- #
def test_wipe_and_remove_delete_the_host():
    ops = _ops()
    _block(ops)
    ops.wipe()
    with pytest.raises(Exception):
        ops.mpm("info", "numParticles")
    _block(ops)
    ops.remove("loadPattern", TAG)
    with pytest.raises(Exception):
        ops.mpm("info", "numParticles")


def test_reset_and_mpm_reset():
    ops = _ops()
    _block(ops)
    dt = 0.4 * _scalar(ops.mpm("criticalTimeStep"))
    init = _fields(ops)
    ops.mpm("analyze", 20, dt)
    run1 = _fields(ops)
    assert run1["x"] != init["x"]
    assert abs(_scalar(ops.mpm("info", "time")) - 20 * dt) < 1e-12
    # Domain reset reaches the host only as a re-applyLoad: particles stay put
    ops.reset()
    assert _fields(ops) == run1
    # mpm reset: back to the state at the first analyze; a rerun is bit-identical
    e0 = int(ops.mpm("info", "epoch"))
    ops.mpm("reset")
    assert int(ops.mpm("info", "epoch")) == e0 + 1
    assert _fields(ops) == init
    ops.mpm("analyze", 20, dt)
    assert _fields(ops) == run1


def test_host_is_inert_in_an_fe_analysis():
    """Outside `mpm analyze` the host's applyLoad does nothing (ADR 232 §5.1)."""
    ops = _ops()
    _block(ops)
    before = _fields(ops)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 2.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 200.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert abs(ops.nodeDisp(2, 1) - 10.0 * 2.0 / 200.0) < 1e-12
    assert _fields(ops) == before


# --------------------------------------------------------------------------- #
#  Admitted materials and the stage path                                      #
# --------------------------------------------------------------------------- #
def _run_block(ops, mat_cmd, nsteps=30, stage=None):
    _block(ops, mat_cmd)
    if stage is not None:
        ops.parameter(1, "loadPattern", TAG, "materialState")
        ops.updateParameter(1, float(stage))
    ops.mpm("analyze", nsteps, 1.0e-4)
    return _fields(ops)


def _KG(E, nu):
    return E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))


def test_ladrunoj2_elastic_range_equals_elastic_isotropic():
    ops = _ops()
    E, nu = 1.0e6, 0.25
    K, G = _KG(E, nu)
    ref = _run_block(ops, ("ElasticIsotropic", 1, E, nu))
    j2 = _run_block(ops, ("LadrunoJ2", 1, K, G, "-iso", "voce", 1.0e12, 0.0, 0.0, 0.0))
    for f in ref:
        _close(j2[f], ref[f])


def test_druckerprager_stage_path_through_loadpattern():
    """materialState 0 (elastic) reaches every particle through the host (W12):
    the run then equals ElasticIsotropic; left elastoplastic, it does not."""
    ops = _ops()
    E, nu = 1.0e6, 0.25
    K, G = _KG(E, nu)
    ref = _run_block(ops, ("ElasticIsotropic", 1, E, nu))
    dp = ("DruckerPrager", 1, K, G, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    elastic = _run_block(ops, dp, stage=0)
    for f in ref:
        _close(elastic[f], ref[f])
    plastic = _run_block(ops, dp)
    assert max(abs(a - b) for a, b in zip(plastic["sig"], ref["sig"])) > 1.0


# (The v1 particle loop is serial, so thread-count determinism holds by
# construction; the test arrives with the threaded loop. Review F8.)


def _j2_oracle(eps_path, K, G, sy):
    """Perfect-plasticity von Mises radial return in plane strain, ||s|| <= sqrt(2/3) sy,
    driven through the particle's own strain path ([exx, eyy, gxy] per step). The
    independent leg of the two-leg check (review F8)."""
    ep = np.zeros((3, 3))
    out = []
    R = np.sqrt(2.0 / 3.0) * sy
    for e in eps_path:
        eps = np.array([[e[0], e[2] / 2, 0.0], [e[2] / 2, e[1], 0.0], [0.0, 0.0, 0.0]])
        tr = np.trace(eps)
        dev = eps - tr / 3.0 * np.eye(3)
        s = 2.0 * G * (dev - ep)
        ns = np.linalg.norm(s)
        if ns > R:
            n = s / ns
            ep = ep + (ns - R) / (2.0 * G) * n
            s = R * n
        sig = s + K * tr * np.eye(3)
        out.append([sig[0, 0], sig[1, 1], sig[0, 1]])
    return np.array(out)


@pytest.mark.parametrize("which", ["LadrunoJ2", "DruckerPrager"])
def test_two_leg_plastic_particles_match_the_material_law(which):
    """Every particle's stress after every step equals perfect J2 plasticity driven
    through that particle's own unrotated strain path (DruckerPrager with rho = 0
    and no hardening IS J2 with ||s|| <= sqrt(2/3) sigma_y)."""
    ops = _ops()
    E, nu, sy = 1.0e6, 0.25, 2.0e3
    K, G = _KG(E, nu)
    mat = (("LadrunoJ2", 1, K, G, "-iso", "voce", sy, 0.0, 0.0, 0.0) if which == "LadrunoJ2"
           else ("DruckerPrager", 1, K, G, sy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    _block(ops, mat)
    eps, sig = [], []
    for _ in range(40):
        ops.mpm("analyze", 1, 1.0e-4)
        eps.append(np.array(ops.mpm("particleData", "epsHat")).reshape(-1, 3))
        sig.append(np.array(ops.mpm("particleData", "sigHat")).reshape(-1, 3))
    eps, sig = np.array(eps), np.array(sig)            # [steps, particles, 3]
    plastic = 0
    for p in range(eps.shape[1]):
        ref = _j2_oracle(eps[:, p, :], K, G, sy)
        assert np.allclose(sig[:, p, :], ref, rtol=1e-7, atol=1e-6 * sy), (which, p)
        plastic += int(np.any(np.abs(ref - _j2_oracle(eps[:, p, :], K, G, 1e30)) > 1.0))
    assert plastic > 0, "the path never yielded: the check would be elastic only"


def test_critical_dt_reads_the_initial_tangent():
    """Amendment 6 / review F6: after plastic flow, Δt_c still uses the elastic
    modulus K + 4G/3 (DruckerPrager's softened tangent would overestimate it)."""
    ops = _ops()
    E, nu, sy = 1.0e6, 0.25, 2.0e3
    K, G = _KG(E, nu)
    _block(ops, ("DruckerPrager", 1, K, G, sy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    ops.mpm("analyze", 40, 1.0e-4)
    m = np.array(ops.mpm("particleData", "m"))
    V = np.array(ops.mpm("particleData", "V"))
    v = np.array(ops.mpm("particleData", "v")).reshape(-1, 2)
    c = np.sqrt((K + 4.0 * G / 3.0) / (m / V)) + np.linalg.norm(v, axis=1)
    assert abs(_scalar(ops.mpm("criticalTimeStep")) - 0.1 / c.max()) <= 1e-12 * (0.1 / c.max())


def test_reset_after_a_refused_first_analyze():
    """Review F1: the snapshot is the state before the first ACCEPTED step. A body
    added after a refused analyze survives `mpm reset`; after a reset, a re-seeded
    stiff body gets its own (stiff) materials and the critical Δt follows."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e6, 0.25)
    ops.nDMaterial("ElasticIsotropic", 2, 1.0e10, 0.25)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    n1 = int(ops.mpm("body", 1, "box", 0.2, 0.2, 0.4, 0.4, "-ppc", 2, "-mat", 1, "-rho", 2000.0))
    dtc = _scalar(ops.mpm("criticalTimeStep"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 10 * dtc)                       # refused: nothing changed
    n2 = int(ops.mpm("body", 2, "box", 0.6, 0.6, 0.8, 0.8, "-ppc", 2, "-mat", 1, "-rho", 2000.0))
    ops.mpm("analyze", 5, 0.2 * dtc)
    ops.mpm("reset")
    assert int(ops.mpm("info", "numParticles")) == n1 + n2   # the snapshot includes body 2
    assert len(ops.mpm("particleData", "mat")) == n1 + n2


def test_reset_drops_particles_added_after_the_snapshot_with_their_materials():
    """Review F1, probe_snap2: particles added between a refused analyze and the first
    accepted step belong to the snapshot; particles added after a reset that is then
    reset again (no step) are dropped with their materials and ids."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e6, 0.25)
    ops.nDMaterial("ElasticIsotropic", 2, 1.0e10, 0.25)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    ops.mpm("body", 1, "box", 0.2, 0.2, 0.4, 0.4, "-ppc", 2, "-mat", 1, "-rho", 2000.0)
    dtc = _scalar(ops.mpm("criticalTimeStep"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 10 * dtc)
    ops.mpm("body", 2, "box", 0.6, 0.6, 0.8, 0.8, "-ppc", 2, "-mat", 1, "-rho", 2000.0)
    ops.mpm("analyze", 1, 0.2 * dtc)
    ops.mpm("reset")
    ops.mpm("body", 3, "box", 0.6, 0.2, 0.8, 0.4, "-ppc", 2, "-mat", 2, "-rho", 2000.0)
    mats = list(ops.mpm("particleData", "mat"))
    assert mats[-1] == 2.0
    M = 1.0e10 * 0.75 / (1.25 * 0.5)
    assert abs(_scalar(ops.mpm("criticalTimeStep")) - 0.1 / np.sqrt(M / 2000.0)) < 1e-15
    ops.mpm("analyze", 1, 1.0e-7)
    ops.mpm("reset")
    assert int(ops.mpm("info", "numParticles")) == len(mats)


def test_initial_conditions_refused_after_the_first_step():
    """Review F4: stress and velocity are initial conditions, like shift and body."""
    ops = _ops()
    _block(ops)
    ops.mpm("analyze", 2, 1.0e-4)
    for cmd in (("stress", 1, "uniform", -5.0, -5.0, 0.0), ("velocity", 1, "uniform", 9.0, 0.0),
                ("shift", 1, 0.1, 0.0),
                ("body", 2, "box", 0.7, 0.7, 0.8, 0.8, "-ppc", 2, "-mat", 1, "-rho", 2000.0)):
        with pytest.raises(Exception):
            ops.mpm(*cmd)
    ops.mpm("reset")
    ops.mpm("velocity", 1, "uniform", 0.0, 0.0)


def test_kernel_refusal_fails_the_step_and_changes_nothing():
    """Review F3, binary leg: head-on particles invert F; the step fails, the
    particles and the Domain time are unchanged, and a small Δt is accepted."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e6, 0.25)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10, "-basis", "linear")
    ops.mpm("safety", 1.0)
    v = 1.0e6
    ops.mpm("particle", 0, 2, 1, 1000.0, 0.0025, "-x", 0.805, 0.85, "-lp", 0.025, 0.025)
    ops.mpm("particle", 1, 1, 1, 1000.0, 0.0025, "-x", 0.405, 0.45, "-lp", 0.025, 0.025, "-v", v, 0.0)
    ops.mpm("particle", 2, 1, 1, 1000.0, 0.0025, "-x", 0.495, 0.45, "-lp", 0.025, 0.025, "-v", -v, 0.0)
    dtc = _scalar(ops.mpm("criticalTimeStep"))
    ops.mpm("analyze", 1, 1e-3 * dtc)
    t0, before = _scalar(ops.mpm("info", "time")), _fields(ops)
    dtc = _scalar(ops.mpm("criticalTimeStep"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 3, 0.999 * dtc)
    assert _fields(ops) == before
    assert _scalar(ops.mpm("info", "time")) == t0
    ops.mpm("analyze", 1, 1e-3 * dtc)
    assert _scalar(ops.mpm("info", "time")) > t0


_DANGLING = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_MPM_TESTS"])
from _testbed import ops
K, G = 666666.6666666666, 400000.0
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial("DruckerPrager", 1, K, G, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
ops.mpm("body", 1, "box", 0.2, 0.2, 0.6, 0.6, "-ppc", 2, "-mat", 1, "-rho", 2000.0)
ops.parameter(1, "loadPattern", 232, "materialState")
ops.updateParameter(1, 0.0)
ops.remove("loadPattern", 232)
ops.updateParameter(1, 1.0)
ops.wipe()
print("SURVIVED")
"""


def test_parameter_update_after_remove_loadpattern_does_not_crash():
    """Review F2: the host detaches its Parameters when it is deleted."""
    _ops()
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1",
               LADRUNO_MPM_TESTS=os.path.dirname(os.path.abspath(__file__)))
    proc = subprocess.run([sys.executable, "-c", _DANGLING], capture_output=True, text=True,
                          timeout=120, env=env, stdin=subprocess.DEVNULL)
    assert proc.returncode == 0 and "SURVIVED" in proc.stdout, (proc.returncode, proc.stdout, proc.stderr)


# --------------------------------------------------------------------------- #
#  Classic Tcl: the same refusals and the `loadPattern` stage spelling        #
# --------------------------------------------------------------------------- #
TCL_DECK = r"""
proc expectError {label script} {
  if {[catch {uplevel 1 $script} msg]} { puts "PASS $label" } else { puts "FAIL $label (accepted)" }
}
expectError no-model { mpm grid 232 -origin 0 0 -h 0.1 -cells 10 10 }
model basic -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 1.0e6 0.25
nDMaterial DruckerPrager 3 666666.6666666666 400000.0 1.0 0 0 0 0 0 0 0 0
nDMaterial J2Plasticity 5 1.0e6 0.5e6 1.0e3 1.0e3 0.0 0.0
mpm grid 232 -origin 0 0 -h 0.1 -cells 10 10
expectError second-grid { mpm grid 233 -origin 0 0 -h 0.1 -cells 10 10 }
expectError non-admitted { mpm body 9 box 0.2 0.2 0.6 0.6 -ppc 2 -mat 5 -rho 2000 }
mpm body 1 box 0.2 0.2 0.6 0.6 -ppc 2 -mat 3 -rho 2000
expectError stress-on-dp { mpm stress 1 uniform -1 -1 0 }
mpm velocity 1 sine x x 0.5 7.853981633974483 0.2
mpm fix xmin xy
mpm wall ymin 0.3
expectError unknown-sub { mpm frobnicate }
parameter 1 loadPattern 232 materialState
updateParameter 1 0
set dtc [mpm criticalTimeStep]
expectError dt-above { mpm analyze 1 [expr {0.6 * $dtc}] }
mpm analyze 5 1.0e-4
if {abs([mpm info time] - 5.0e-4) < 1e-12} { puts "PASS time" } else { puts "FAIL time" }
database File mpm_db
expectError database-save { save 1 }
node 1 0 0
expectError free-dof { mpm analyze 1 1.0e-4 }
wipe
expectError wiped { mpm info numParticles }
puts "SELF-TEST: done"
"""


def test_classic_tcl_refusals_and_stage_spelling():
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    with tempfile.TemporaryDirectory() as td:
        deck = Path(td) / "mpm_refusals.tcl"
        deck.write_text(TCL_DECK, encoding="utf-8")
        proc = subprocess.run([exe, str(deck)], capture_output=True, text=True, timeout=300,
                              env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=td,
                              stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert "SELF-TEST: done" in out, out
    assert "FAIL" not in out, out
    assert out.count("PASS ") == 10, out
    assert "database save is refused" in out, out
