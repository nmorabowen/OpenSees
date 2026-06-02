"""LadrunoBrick — single-point output mirror (the redundant-eval fix).

The single-integration-point formulations (`eas`, `uri+stiffness`, `uri+viscous`)
now evaluate the constitutive model ONCE at the centroid (material slot 0) instead
of pushing the same centroid strain to all 8 material copies. Slots 1-7 are no
longer strained; the per-Gauss-point output is *mirrored* from slot 0 in
getResponse / the material response channel.

These tests pin the observable contract of that change:
  * single-point formulations report IDENTICAL stress/strain at all 8 GPs, equal
    to the element-average (centroid) value, and physically correct (≠ 0);
  * the per-GP `material` response channel returns the centroid value for ANY gp;
  * full-integration formulations (std/bbar) still report genuinely DISTINCT
    per-GP stresses under a non-uniform (bending) load — i.e. the mirror is
    applied ONLY to the single-point formulations, not universally.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


def _gp_stress_blocks(tag=1):
    """eleResponse 'stresses' -> 48 values = 8 GP blocks of 6. Returns list of 8."""
    s = list(ops.eleResponse(tag, "stresses"))
    return [s[6 * g:6 * g + 6] for g in range(8)]


def _uniaxial_cube(form_args, sig0=2.0, E=1000.0, nu=0.3):
    """Single hex, 1/8-symmetry restraints, uniform uniaxial-x traction → a
    constant (uniform) stress state in every formulation."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, sig0 / 4.0, 0.0, 0.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0


@pytest.mark.parametrize("form_args", [
    ["-formulation", "eas"],
    ["-formulation", "uri", "-hourglass", "stiffness"],
    ["-formulation", "uri", "-hourglass", "viscous"],
])
def test_singlepoint_mirrors_all_gauss_points(form_args):
    """Single-point formulations report the centroid stress at every GP: all 8
    blocks identical, nonzero (σ11≈sig0), and equal to stress3D6 (the average)."""
    sig0 = 2.0
    _uniaxial_cube(form_args, sig0=sig0)
    blocks = _gp_stress_blocks()
    # σ11 is recovered and nonzero
    assert abs(blocks[0][0] - sig0) < 1e-6, f"{form_args}: σ11 {blocks[0][0]} != {sig0}"
    # every GP block equals slot 0 (the mirror)
    for g in range(1, 8):
        for j in range(6):
            assert abs(blocks[g][j] - blocks[0][j]) < 1e-12, (
                f"{form_args}: GP{g+1} comp{j} {blocks[g][j]} != centroid {blocks[0][j]}"
            )
    # and equals the element-average channel
    avg = list(ops.eleResponse(1, "stress3D6"))
    for j in range(6):
        assert abs(avg[j] - blocks[0][j]) < 1e-9, f"{form_args}: stress3D6 != centroid"


@pytest.mark.parametrize("form_args", [
    ["-formulation", "eas"],
    ["-formulation", "uri", "-hourglass", "stiffness"],
])
def test_material_channel_routes_to_centroid(form_args):
    """The per-GP `material` response channel returns the live centroid material
    (slot 0) for ANY requested Gauss point in a single-point formulation."""
    _uniaxial_cube(form_args, sig0=2.0)
    s1 = list(ops.eleResponse(1, "material", 1, "stress"))
    for gp in (2, 5, 8):
        sg = list(ops.eleResponse(1, "material", gp, "stress"))
        assert sg, f"{form_args}: material channel gp{gp} returned nothing"
        for a, b in zip(s1, sg):
            assert abs(a - b) < 1e-12, f"{form_args}: gp{gp} stress != gp1 (not mirrored)"


def test_full_integration_keeps_distinct_gauss_points():
    """Guard that the mirror is NOT applied to full integration: under a bending
    load, std reports genuinely DISTINCT per-GP stresses (the tension/compression
    sides differ). If this fails, the single-point mirror leaked into std/bbar."""
    L, E, nu, P = 4.0, 1000.0, 0.0, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # one slender hex, clamp x=0 face, transverse z tip load -> bending gradient
    coords = {1: (0, 0, 0), 2: (L, 0, 0), 3: (L, 1, 0), 4: (0, 1, 0),
              5: (0, 0, 1), 6: (L, 0, 1), 7: (L, 1, 1), 8: (0, 1, 1)}
    for t, c in coords.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "std")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.0, 0.0, P / 4.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    blocks = _gp_stress_blocks()
    # σ11 (bending fiber stress) must differ across GPs (top vs bottom fibers)
    s11 = [b[0] for b in blocks]
    assert max(s11) - min(s11) > 1e-6, (
        f"std reports uniform σ11 across GPs ({s11}) — the single-point mirror "
        f"must NOT apply to full integration"
    )
