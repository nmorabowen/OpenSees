"""LadrunoRCConcrete (classTag 33015) -- Ladruno WP-241: the spatial add-on plane, the damage
crack gate, and the condensed initial tangent. Zone-A, the material as integrated in OpenSees.

* ``-addonPlane spatial``: beta, tension stiffening and the cracked-nu latch read eps1 and p1
  from the FULL 3x3 strain tensor instead of the x-y block. Needed for beam fibres, whose
  shear sits in gamma_xy AND gamma_xz (the BeamFiber view, a later WP of program #1131). The fixed-crack
  interlock is x-y only and is refused with it.
* ``-crackGate damage`` (the spatial default): the cracked-nu latch and the TS floor open only
  once tensile damage has started. With the strain gate (eps1 >= ft/E, the membrane default)
  Poisson expansion opens them in plain compression -- the WP-241 probe: a plane-stress point
  in uniaxial compression latched at ~0.42 fc and the step diverged.
* ``getInitialTangent`` of the condensed views (PlaneStress, PlateFiber) is the sigma_33 = 0
  Schur complement, so ``Newton -initial`` iterates against E/(1 - nu^2).

The kernel tangents of the new code are gated by the g++ driver
tests/_testbed/rc_spatial_gpp.cpp (test_ladrunoRCConcrete_spatial_cpp.py). With the new flags
at their defaults the material is unchanged (the rest of test_ladrunoRCConcrete_*.py).
"""
import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]


def _rc(tag, *flags, cls="LadrunoRCConcrete"):
    ops.nDMaterial(cls, tag, E, NU,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC, *flags)


def _cube(flags, disp, nsteps=60):
    """Homogeneous strain on one stdBrick: every nodal displacement prescribed, u = disp(x, y, z)
    scaled by the load factor. Returns the stress history (6 comps) and the nuCracked history."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    _rc(1, *flags)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        u = disp(x, y, z)
        for d in range(3):
            ops.sp(t, d + 1, u[d])
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    sig, lat = [], []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
        sig.append(list(ops.eleResponse(1, "stresses"))[0:6])
        lat.append(ops.eleResponse(1, "material", 1, "nuCracked")[0])
    return sig, lat


_MCFT = ["-beta", "-tensStiff", "vc", "-crackedNu", 0.0]


@pytest.mark.t1
def test_spatial_reduces_to_membrane_on_an_xy_path():
    """An x-y state whose largest principal strain lies in x-y (eps_zz = 0 < eps1, no transverse
    shear): the spatial plane must give the membrane answer (strain gate on both)."""
    def disp(x, y, z):
        return (-2.0e-3 * x + 1.0e-3 * y, 3.0e-3 * y, 0.0)
    m, lm = _cube(_MCFT, disp)
    s, ls = _cube(_MCFT + ["-addonPlane", "spatial", "-crackGate", "strain"], disp)
    for a, b in zip(m, s):
        assert b == pytest.approx(a, rel=1e-10, abs=1e-9)
    assert lm == ls and lm[-1] == 1.0


@pytest.mark.t1
def test_spatial_is_invariant_under_xy_xz_shear_and_membrane_is_not():
    """Pure shear gamma_xy vs gamma_xz, same history. Spatial: tau_xy == tau_xz and
    sigma_yy == sigma_zz (y <-> z swap). Membrane reads only x-y, so the x-z shear gets no
    beta, no TS floor, no latch: the responses differ (break-once evidence for the gate)."""
    g = 4.0e-3
    for plane, invariant in (("spatial", True), ("membrane", False)):
        flags = _MCFT + ["-addonPlane", plane, "-crackGate", "strain"]
        a, la = _cube(flags, lambda x, y, z: (g * y, 0.0, 0.0))
        b, lb = _cube(flags, lambda x, y, z: (g * z, 0.0, 0.0))
        worst = max(max(abs(p[3] - q[5]), abs(p[1] - q[2]), abs(p[0] - q[0])) for p, q in zip(a, b))
        if invariant:
            assert worst < 1e-9, worst
            assert la == lb
        else:
            assert worst > 0.05, worst


def _plane_stress_compression(flags, eyy_end=-1.8e-3, nsteps=60):
    """One PlaneStress quad, uniaxial compression in y, laterally free (the WP-241 probe).
    Returns (steps converged, nuCracked history, sigma_yy history)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _rc(1, *flags)
    for n, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(n, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(3, 2, eyy_end)
    ops.sp(4, 2, eyy_end)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    lat, syy = [], []
    for k in range(nsteps):
        if ops.analyze(1) != 0:
            return k, lat, syy
        lat.append(ops.eleResponse(1, "material", 1, "nuCracked")[0])
        syy.append(ops.eleResponse(1, "material", 1, "stress")[1])
    return nsteps, lat, syy


@pytest.mark.t1
def test_damage_gate_keeps_plane_stress_compression_uncracked():
    """The WP-241 probe as a regression. Strain gate (default): Poisson expansion reaches
    eps_cr = ft/E at |eyy| = eps_cr/nu = 5e-4 and the cracked-nu latch fires in plain
    compression (no crack). Damage gate: the Lubliner tensile measure is zero without a
    tensile principal stress, so nothing latches and the run reaches ~0.9 fc' with the
    elastic-Poisson response."""
    n_s, lat_s, _ = _plane_stress_compression(_MCFT)
    assert (n_s < 60) or (max(lat_s) == 1.0), "strain gate no longer reproduces the probe defect"
    n_d, lat_d, syy_d = _plane_stress_compression(_MCFT + ["-crackGate", "damage"])
    assert n_d == 60
    assert max(lat_d) == 0.0
    assert min(syy_d) < -0.8 * max(CS)


@pytest.mark.t1
def test_damage_gate_still_latches_in_tension():
    """Uniaxial tension past the backbone peak (eps = 1e-4): the damage gate latches (one step
    after the threshold passes the peak; the predictor needs nu before the step's damage)."""
    _, lat = _cube(_MCFT + ["-addonPlane", "spatial"],
                     lambda x, y, z: (6.0e-4 * x, 0.0, 0.0), nsteps=60)
    first = lat.index(1.0)
    assert 0 < first < 60
    assert lat[first - 1] == 0.0


@pytest.mark.t1
def test_plane_stress_initial_tangent_is_condensed():
    """Elastic plane-stress quad under a small nodal load, Newton -initial with a 2-iteration
    budget. With the sigma_33 = 0 condensed initial tangent the linear problem converges in one
    correction (the 2nd iterate is ~0). With the raw C0 block (lambda + 2 mu, 6 % too stiff
    in plane) the iterates contract by ~0.06 per step and miss 1e-12 in two."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _rc(1)
    for n, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(n, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, 0.5)
    ops.load(4, 0.0, 0.5)              # sigma_yy = 1 MPa, below ft = 3
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 2, 0)
    ops.algorithm("Newton", "-initial")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 2) == pytest.approx(1.0 / E, rel=1e-9)


def test_serialization_carries_the_wp241_params():
    """Schema v7: addonPlane and crackGate survive a database round trip. The path ends on
    an elastic unloading step (a re-evaluated damaging step is not idempotent in this kernel,
    see test_serialization_carries_the_c2_params)."""
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _CUBE.items():
            ops.node(t, *c)
        _rc(1, *_MCFT, "-addonPlane", "spatial", "-crackGate", "damage")
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t, (x, y, z) in _CUBE.items():
            ops.sp(t, 1, 1.0e-3 * z)       # gamma_xz: only the spatial plane sees this shear
            ops.sp(t, 2, 0.0)
            ops.sp(t, 3, 0.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Penalty", 1.0e14, 1.0e14)
        ops.test("NormDispIncr", 1.0e-10, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 0.1)
        ops.analysis("Static")
        for _ in range(10):
            assert ops.analyze(1) == 0
        ops.integrator("LoadControl", -0.1)   # elastic unloading step
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
    database_roundtrip(build, probe_nodes=[7], ndf=3, dbname="rc_wp241",
                       probe_fn=lambda: list(ops.eleResponse(1, "stresses"))[0:6]
                       + list(ops.eleResponse(1, "material", 1, "nuCracked")))


@pytest.mark.parametrize("cls", ["LadrunoRCConcrete", "LadrunoRCFiniteStrain"])
@pytest.mark.parametrize("flags", [
    ["-addonPlane", "spatial", "-interlock"],
    ["-addonPlane", "spatial", "-cyclic"],
    ["-addonPlane", "spatial", "-xcrack"],
    ["-addonPlane", "spatial", "-shearRetention", "rots"],
    ["-addonPlane", "plane"],
    ["-crackGate", "stress"],
    ["-addonPlane"],
])
def test_parser_refuses(cls, flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    with pytest.raises(Exception):
        _rc(1, *flags, cls=cls)


@pytest.mark.parametrize("cls", ["LadrunoRCConcrete", "LadrunoRCFiniteStrain"])
def test_parser_accepts(cls):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _rc(1, "-addonPlane", "spatial", "-crackGate", "damage", "-beta", "-tensStiff", "vc", cls=cls)
    _rc(2, "-addonPlane", "membrane", "-crackGate", "damage", "-interlock", cls=cls)
