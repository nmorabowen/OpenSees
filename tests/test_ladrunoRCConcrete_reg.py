"""LadrunoRCConcrete Phase-3b crack-band (Bazant-Oh) regularization (`-autoRegularization`,
classTag 33015) — Zone-A battery, the material AS INTEGRATED IN OPENSEES.

`-autoRegularization $lch_ref` rescales the post-peak softening so the regularized specific
fracture energy is g_reg = G_f0 * (lch_ref/lch) — i.e. the dissipated energy DENSITY * lch (the
physical fracture energy) is mesh-objective. lch is the element's getCharacteristicLength()
unless `-lch` is supplied. Default off => baseline-identical. The regularize() math is proven
first in the numpy oracle (rc_shell_ref.py run_R1) and a standalone g++ build (rc_reg_gpp.cpp);
here we verify the C++ material end-to-end on a uniaxial-stress stdBrick, using `-lch` to set the
"element size" deterministically.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
# tension backbone with a clear softening branch all the way to zero stress
TE = [0.0, 0.0001, 0.0010, 0.0040]
TS = [0.0, 3.0,    0.5,    0.0]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0, 0.999]
LCH_REF = 50.0


def _rc(tag, auto_reg=False, lch=0.0, lch_ref=LCH_REF):
    args = ["LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if auto_reg:
        args += ["-autoRegularization", lch_ref]
        if lch > 0.0:
            args += ["-lch", lch]
    ops.nDMaterial(*args)


def _run(mat_fn, eps_target, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.25, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 1, eps_target / nsteps)
    ops.analysis("Static")
    eps, sig = [], []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"reg analyze failed (target {eps_target})"
        ops.eleResponse(1, "forces")
        eps.append(ops.nodeDisp(2, 1))
        sig.append(list(ops.eleResponse(1, "stresses"))[0])
    return np.array(eps), np.array(sig)


_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # np.trapz removed in numpy 2.0


def _postpeak_energy(eps, sig):
    ipk = int(np.argmax(sig))
    return float(_trapz(sig[ipk:], eps[ipk:]))       # dissipated energy density past the peak


def test_reg_reduce_noop():
    """-autoRegularization with lch == lch_ref is a no-op (lch_scale==1 quick return) => the
    trajectory matches the un-regularized material bit-for-bit."""
    base = _run(lambda t: _rc(t), 4.0e-3, 120)
    reg = _run(lambda t: _rc(t, auto_reg=True, lch=LCH_REF), 4.0e-3, 120)
    for (e0, s0), (e1, s1) in zip(zip(*base), zip(*reg)):
        assert abs(s1 - s0) <= 1.0e-7 * (abs(s0) + 1.0), (e0, s0, s1)


def test_reg_energy_mesh_objectivity():
    """The post-peak dissipated energy DENSITY * lch must be mesh-objective: g_reg(lch)*lch ==
    G_f0*lch_ref, constant across element sizes. Drive identical uniaxial tension at lch in
    {50,25,12.5} (scale 1,2,4 — all stretch, no gmin clamp) and compare g_density*lch."""
    prods = []
    for lch in (50.0, 25.0, 12.5):
        eps, sig = _run(lambda t, L=lch: _rc(t, auto_reg=True, lch=L), 0.02, 400)
        g = _postpeak_energy(eps, sig)
        prods.append(g * lch)
    ref = prods[0]
    for p in prods:
        assert p == pytest.approx(ref, rel=0.05), prods
    # smaller lch must give a more ductile (higher post-peak energy density) response
    e25, s25 = _run(lambda t: _rc(t, auto_reg=True, lch=25.0), 0.02, 400)
    e50, s50 = _run(lambda t: _rc(t, auto_reg=True, lch=50.0), 0.02, 400)
    assert _postpeak_energy(e25, s25) > 1.5 * _postpeak_energy(e50, s50)


def test_reg_element_lch_path():
    """Without -lch the material must latch the element's getCharacteristicLength() and still
    run + soften (the production path inside an element)."""
    eps, sig = _run(lambda t: _rc(t, auto_reg=True, lch_ref=2.0), 4.0e-3, 150)
    assert sig.max() == pytest.approx(3.0, rel=0.15)        # tensile peak ~ ft
    assert sig[-1] < 0.5 * sig.max()                        # softened
    assert sig[-1] >= -1.0e-6                               # not pathological


def test_reg_parser_guard():
    """-autoRegularization requires lch_ref > 0."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoRCConcrete", 1, E, NU, "-Ce", *CE, "-Cs", *CS,
                       "-Te", *TE, "-Ts", *TS, "-autoRegularization", 0.0)
