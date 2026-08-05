"""PDMY01 (``PressureDependMultiYield``) plastic-strain recorder channels.

PDMY01 is hypoelastic multi-surface plasticity: it stores nested yield surfaces,
a PPZ (dilation) state and the current stress, but it has never stored a plastic
strain.  The fork adds three recordable channels driven by one accumulator in
``accumPlasticStrain()``, which splits each *sub-step's* strain increment with
the SAME pressure-dependent moduli the elastic predictor in ``setTrialStress()``
just used::

    d_eps^p = d_eps - dev(d_sigma)/(2G) - (d_p/(3B)) I

  * ``plasticStrain``            -- the tensor (engineering shear; 3 comps in 2D)
  * ``equivalentPlasticStrain``  -- int sqrt(2/3 de^p_dev : de^p_dev)  (von Mises
                                    convention, matching LadrunoJ2's ``ebarP``)
  * ``volumetricPlasticStrain``  -- trace, signed (the dilatancy signal)

Driver: a 1x1x2 stack of ``stdBrick`` under a fully prescribed *uniform* strain
path (three roller faces + prescribed normal displacement on the three positive
faces).  Every Gauss point sees the same state, so the material response can be
differenced step by step and reconstructed independently in numpy.

The gates:

  G1  stage 0 (elastic) accrues exactly zero plastic strain.
  G2  ebarP is monotone non-decreasing (it is a path length).
  G3  an elastic unload window after yielding accrues ~0 -- this is the sharp
      one: it proves the hypoelastic inversion is exact, moduli and all, and it
      is the check that would fail on any G/B mix-up.
  G4  a full numpy reconstruction of the plastic strain from the recorded
      (stress, strain) history matches the engine.  When the model takes one
      sub-step per load step, the reconstruction is EXACT, so this is asserted
      tightly (1e-9 relative).
  G5  volumetric plastic strain is contractive on this low-stress-ratio path,
      as a contractive sand (contract=0.05, eta well below the PT line) must be.

G4 is exact only because the oracle reproduces the model's own moduli, which
needs PDMY01's *actual* defaults, not the documented ones:
``OPS_PressureDependMultiYield`` sets ``param[21] = .3`` for cohesion while its
own usage string says ``cohesi (=.5)`` (upstream typo), and ``setUpSurfaces``
stores ``residualPress`` POSITIVE while ``refPressure`` is negative.  Getting
either wrong perturbs G by ~0.15%, which is invisible in a stress check but
throws the plastic strain -- a small difference of two large numbers -- off by
a factor of 5.  See LEDGER_quirks.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST)

pytestmark = [pytest.mark.zone_a]   # hand-built mesh, no gmsh

# medium-loose sand, UCSD 16-arg signature, nd=3, liquefaction OFF
P = dict(rho=2.0, G=5.5e4, B=1.5e5, phi=33.0, gammaPeak=0.1, refPress=80.0,
         d=0.5, PT=27.0, contract=0.05, dil1=0.0, dil2=0.0,
         liq1=0.0, liq2=0.0, liq3=0.0, NYS=20)

_E0 = 1.5e-4      # lateral strain at t=1 (keeps p' near refPress)
_DELTA = 3.0      # axial/lateral ratio - 1; sets the (sub-failure) stress ratio
_NSTEP = 240
_T_FLIP = 0.25    # pseudo-time of updateMaterialStage 0 -> 1

_XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]


def _build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, 0.5 * k)
    ops.nDMaterial("PressureDependMultiYield", 1, 3, P["rho"], P["G"], P["B"],
                   P["phi"], P["gammaPeak"], P["refPress"], P["d"], P["PT"],
                   P["contract"], P["dil1"], P["dil2"], P["liq1"], P["liq2"],
                   P["liq3"], P["NYS"])
    for e in (1, 2):
        b = 4 * (e - 1)
        ops.element("stdBrick", e, b+1, b+2, b+3, b+4, b+5, b+6, b+7, b+8, 1)

    # rollers on the three negative faces
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)

    # prescribed uniform strain on the three positive faces
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -_E0)
            if y == 1.:
                ops.sp(n, 2, -_E0)
            if k == 2:
                ops.sp(n, 3, -_E0 * (1. + _DELTA))

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / _NSTEP)
    ops.analysis("Static")


def _probe():
    return (np.array(ops.eleResponse(1, "material", "1", "stress")[:6]),
            np.array(ops.eleResponse(1, "material", "1", "strain")[:6]),
            np.array(ops.eleResponse(1, "material", "1", "plasticStrain")),
            ops.eleResponse(1, "material", "1", "equivalentPlasticStrain")[0],
            ops.eleResponse(1, "material", "1", "volumetricPlasticStrain")[0])


@pytest.fixture(scope="module")
def run():
    _build()
    hist = []
    n_flip = int(round(_T_FLIP * _NSTEP))
    for s in range(_NSTEP):
        assert ops.analyze(1) == 0, f"load step {s} failed"
        hist.append(_probe())
        if s + 1 == n_flip:
            ops.updateMaterialStage("-material", 1, "-stage", 1)

    n_load = len(hist)
    ops.integrator("LoadControl", -1.0 / _NSTEP)          # elastic unload leg
    for s in range(_NSTEP // 8):
        assert ops.analyze(1) == 0, f"unload step {s} failed"
        hist.append(_probe())

    return dict(
        sig=np.array([h[0] for h in hist]),
        eps=np.array([h[1] for h in hist]),
        epl=np.array([h[2] for h in hist]),
        ebar=np.array([h[3] for h in hist]),
        evol=np.array([h[4] for h in hist]),
        n_flip=n_flip, n_load=n_load)


def test_sendself_recvself_round_trip(tmp_path):
    """`sendSelf`/`recvSelf` carry the accumulators, and the yield surfaces
    still land on their upstream offsets.

    The accumulators are APPENDED past the surfaces (slots 70+n*8 .. 76+n*8),
    so the vanilla layout stays a strict prefix of ours -- the minimal-and-
    additive rule for vanilla files.  An earlier revision inserted them at 70
    and pushed every surface to 77; that was self-consistent (send and recv
    agreed) and therefore invisible to a pure round-trip.

    Mutation-verified: shifting ONLY the recvSelf surface offset by one slot
    (70 -> 71) is caught by gate (b) below -- the mobilised ratio goes to `inf`,
    because the surface size gets read out of the middle of a stress vector.
    Gate (c) then drives the model forward as a second, independent net.

    Honest scope: this gates the format against future edits.  It would not
    have caught the offset choice itself, which was a layout-policy defect, not
    a correctness one.

    Uses plain `Newton`, not the module default `KrylovNewton`: restoring a
    database resets the algorithm, and Krylov's accumulator restarts with it.
    PDMY commits `currentStress = trialStress` from the LAST getStress(), which
    is the second-to-last Newton iterate, so the committed state is
    iteration-path dependent -- a restarted Krylov run diverges by ~6e-4
    relative even though every byte of material state round-tripped. Memoryless
    Newton takes the identical path and reproduces to ~4e-14, which is what lets
    this assert a real tolerance instead of a decorative one.
    """
    def _to_plastic_state():
        _build()
        ops.algorithm("Newton")           # memoryless -- see the docstring
        n_flip = 40
        for s in range(80):
            assert ops.analyze(1) == 0
            if s + 1 == n_flip:
                ops.updateMaterialStage("-material", 1, "-stage", 1)

    def _snapshot():
        return dict(
            stress=np.array(ops.eleResponse(1, "material", "1", "stress")),
            strain=np.array(ops.eleResponse(1, "material", "1", "strain")),
            epl=np.array(ops.eleResponse(1, "material", "1", "plasticStrain")),
            ebar=ops.eleResponse(1, "material", "1", "equivalentPlasticStrain")[0],
            evol=ops.eleResponse(1, "material", "1", "volumetricPlasticStrain")[0])

    # reference: uninterrupted run, 10 steps past the save point
    _to_plastic_state()
    for _ in range(10):
        assert ops.analyze(1) == 0
    reference = _snapshot()

    # same run, but save/restore at the 80-step mark before continuing
    _to_plastic_state()
    before = _snapshot()
    db = str(tmp_path / "pdmy")
    ops.database("File", db)
    ops.save(1)
    ops.restore(1)
    after = _snapshot()

    # (a) every channel survives the round-trip bit-identically
    for k in ("epl", "ebar", "evol"):
        assert np.array_equal(np.asarray(before[k]), np.asarray(after[k])), \
            f"{k} changed across save/restore"
    assert np.array_equal(before["stress"][:6], after["stress"][:6])
    assert np.array_equal(before["strain"], after["strain"])

    # (b) the mobilised stress ratio is derived from committedSurfaces[NYS].size(),
    #     so it is a direct probe of the surface block landing on the right offset
    assert abs(before["stress"][6] - after["stress"][6]) < 1e-12, \
        "mobilised ratio moved -- the surface block did not round-trip"

    # (c) the real gate: keep driving.  Wrong surface offsets survive (a) and (b)
    #     but diverge as soon as the model has to use the surfaces again.
    for _ in range(10):
        assert ops.analyze(1) == 0
    continued = _snapshot()
    for k in ("ebar", "evol"):
        assert abs(continued[k] - reference[k]) <= 1e-10 * max(abs(reference[k]), 1e-30), \
            f"{k} diverged after restore: {continued[k]:.12e} vs {reference[k]:.12e}"
    # atol scaled to each tensor's own magnitude: the zero components here carry
    # ~1e-24 of noise against a ~1e-7 signal, and atol=0 would gate on that noise
    assert np.allclose(continued["epl"], reference["epl"], rtol=1e-10,
                       atol=1e-10 * np.abs(reference["epl"]).max())
    assert np.allclose(continued["stress"][:6], reference["stress"][:6], rtol=1e-10,
                       atol=1e-10 * np.abs(reference["stress"][:6]).max())


def test_channels_exist_and_are_sized(run):
    assert run["epl"].shape[1] == 6          # 3D -> 6 components
    assert np.isfinite(run["ebar"]).all()
    assert np.isfinite(run["evol"]).all()


def test_g1_elastic_stage_accrues_nothing(run):
    """Stage 0 is linear elastic: no plastic strain, exactly."""
    n = run["n_flip"]
    assert np.abs(run["ebar"][:n]).max() == 0.0
    assert np.abs(run["epl"][:n]).max() == 0.0
    assert np.abs(run["evol"][:n]).max() == 0.0


def test_g2_equivalent_plastic_strain_is_monotone(run):
    """ebarP is a path length -- it can never decrease, including on unload."""
    d = np.diff(run["ebar"])
    assert d.min() >= -1e-18, f"ebarP decreased by {d.min():.3e}"
    assert run["ebar"][run["n_load"] - 1] > 1e-9, "no plasticity was mobilised"


def test_g3_elastic_unload_window_accrues_nothing(run):
    """The sharp one: after yielding, reverse the path.  PDMY drops to
    activeSurfaceNum = 0 and the sub-steps are purely elastic, so the
    hypoelastic inversion must return exactly zero -- G, B and all."""
    eb = run["ebar"]
    n = run["n_load"]
    d = np.diff(eb[n - 1:])
    scale = eb[n - 1]
    assert scale > 1e-9
    assert d.max() / scale < 1e-10, (
        f"elastic unload accrued {d.max():.3e} against ebarP {scale:.3e}")


def test_g4_matches_independent_numpy_reconstruction(run):
    """Rebuild d_eps^p from the recorded (stress, strain) history using the
    model's own pressure-dependent moduli.  One sub-step per load step here,
    so the reconstruction is exact, not approximate."""
    sig, eps = run["sig"], run["eps"]
    n0, n1 = run["n_flip"], run["n_load"]

    sin_phi = np.sin(np.radians(P["phi"]))
    m_nys = 6. * sin_phi / (3. - sin_phi)
    # OPS_PressureDependMultiYield param[21]; the usage string's ".5" is wrong
    residual_press = max(2. * 0.3 / m_nys, 0.0001 * 101.)   # stored POSITIVE
    ref_press = -P["refPress"]                              # stored NEGATIVE

    def factor(p):
        return max(((p - residual_press) / (ref_press - residual_press))
                   ** P["d"], 1e-10)

    def tensorial(v):
        w = v.copy()
        w[3:] /= 2.0
        return w

    ep = np.zeros(6)
    ebar = 0.0
    for s in range(n0, n1):
        dsig = sig[s] - sig[s-1]
        deps = tensorial(eps[s]) - tensorial(eps[s-1])
        f = factor(sig[s-1][:3].sum() / 3.)
        two_g, three_b = 2. * P["G"] * f, 3. * P["B"] * f
        dp = dsig[:3].sum() / 3.
        dep = deps.copy()
        dep[:3] -= (dsig[:3] - dp) / two_g + dp / three_b
        dep[3:] -= dsig[3:] / two_g
        dev = dep.copy()
        dev[:3] -= dep[:3].sum() / 3.
        ebar += np.sqrt(2. / 3. * ((dev[:3] ** 2).sum() + 2. * (dev[3:] ** 2).sum()))
        ep += dep

    got = run["ebar"][n1 - 1]
    assert abs(ebar - got) <= 1e-9 * abs(got), \
        f"ebarP oracle {ebar:.9e} vs engine {got:.9e}"

    ep_eng = np.concatenate([ep[:3], 2. * ep[3:]])
    assert np.allclose(ep_eng, run["epl"][n1 - 1], rtol=1e-9,
                       atol=1e-9 * np.abs(ep_eng).max())
    assert abs(ep[:3].sum() - run["evol"][n1 - 1]) <= 1e-9 * abs(ep[:3].sum())


def test_2d_reduced_view():
    """ndm=2 returns the 3-component [xx, yy, gxy] view, matching what the
    ``strain`` channel already does -- a separate branch in getResponse()."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in enumerate([(0., 0.), (1., 0.), (1., 1.), (0., 1.)], start=1):
        ops.node(t, x, y)
    ops.nDMaterial("PressureDependMultiYield", 1, 2, P["rho"], P["G"], P["B"],
                   P["phi"], P["gammaPeak"], P["refPress"], P["d"], P["PT"],
                   P["contract"], P["dil1"], P["dil2"], P["liq1"], P["liq2"],
                   P["liq3"], P["NYS"])
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, -_E0)
    ops.sp(3, 1, -_E0)
    ops.sp(3, 2, -_E0 * (1. + _DELTA))
    ops.sp(4, 2, -_E0 * (1. + _DELTA))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / 40)
    ops.analysis("Static")

    for s in range(40):
        assert ops.analyze(1) == 0, f"2D step {s} failed"
        if s == 9:
            ops.updateMaterialStage("-material", 1, "-stage", 1)

    epl = ops.eleResponse(1, "material", "1", "plasticStrain")
    strain = ops.eleResponse(1, "material", "1", "strain")
    assert len(epl) == 3, f"expected the 3-component 2D view, got {len(epl)}"
    assert len(epl) == len(strain), "plasticStrain must match the strain view"
    ebar = ops.eleResponse(1, "material", "1", "equivalentPlasticStrain")
    evol = ops.eleResponse(1, "material", "1", "volumetricPlasticStrain")
    assert len(ebar) == 1 and len(evol) == 1
    assert ebar[0] > 0.0, "2D run mobilised no plasticity"
    assert np.isfinite(epl).all()


def test_g5_contractive_path_gives_contractive_plastic_volume(run):
    """contract=0.05 with eta well below the phase-transformation line: the
    plastic volumetric strain must be compactive (negative in this
    compression-negative convention), never dilative."""
    ev = run["evol"][run["n_load"] - 1]
    assert ev < 0.0, f"expected plastic compaction, got {ev:.3e}"
    assert abs(ev) < abs(run["ebar"][run["n_load"] - 1]), \
        "volumetric plastic strain should be the smaller part here"
