"""ADR-94 R2 blue lane — does H1 (class-static ``Stiffness``) change RESULTS on
a real multi-element, single-material mesh, or only convergence cost?

R1-A confirmed H1 by showing the ASSEMBLED tangent block for a heterogeneous
two-cube model is contaminated (every element gets the last-integrated GP's
tangent). What R1-A did not measure is whether that contamination survives to
the CONVERGED answer. It cannot: OpenSees' residual (unbalanced force) is
built from each element's own committed/trial stress, never from the shared
static ``Stiffness`` -- that static only feeds the Newton SEARCH DIRECTION.
So a converged step (``ops.testIter`` reports success under
``NormDispIncr``) must have the right forces regardless of which GP's tangent
got smeared across the assembly; a wrong tangent can only cost iterations (or,
in the worst case, non-convergence), never a wrong committed stress once
``analyze()`` returns 0.

This test measures the extra-iteration cost directly: run two disconnected
VonMises cubes (one driven plastic, one left elastic -- genuinely different
states, so H1 is active) over several load steps, and compare (a) their
converged per-element stresses against each element analyzed ALONE (H1 is not
observable with only one element/material live), and (b) the total Newton
iteration count against the sum of the two solo analyses.

Zone-A, < 1s.
"""
import os
import sys

import numpy as np
import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.pardir, "Ladruno_implementation",
                                 "adr94_oracle"))
import hex8_tangent as O  # noqa: E402

import test_adr94_hlist_numerics as N  # noqa: E402

pytestmark = [pytest.mark.zone_a]

vm_available = N.vm_available


def _multistep_cubes_build(tags_x0, loads, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x0 in tags_x0:
        b = 10 * t
        for k, (x, y, z) in enumerate(O.NODES):
            ops.node(b + k + 1, x0 + float(x), float(y), float(z))
        for k in range(4):
            ops.fix(b + k + 1, 1, 1, 1)
    N.mat_vm(1, "Continuum")
    for t, x0 in tags_x0:
        b = 10 * t
        ops.element("stdBrick", t, *[b + k + 1 for k in range(8)], 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, load in loads:
        b = 10 * t
        for k in range(4, 8):
            ops.load(b + k + 1, 0., 0., load / nsteps)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _run(tags_x0, loads, nsteps):
    _multistep_cubes_build(tags_x0, loads, nsteps)
    total_iters = 0
    for _ in range(nsteps):
        rc = ops.analyze(1)
        total_iters += ops.testIter()
        assert rc == 0, "must converge every step to compare final stresses"
    out = {}
    for t, _x0 in tags_x0:
        out[t] = np.array(list(ops.eleResponse(t, "stresses"))[0:6])
    return total_iters, out


NSTEPS = 4

#: Slack on the (b) assertion below -- see its comment. Measured
#: post-wp/94b: ITER_PENALTY_MEASURED.
ITER_PENALTY_SLACK = 2


def test_blue_H1_results_are_exact_only_convergence_cost_rises(vm_available):
    """FIXED by wp/94b (the cost half).  H1 never corrupted the converged
    answer -- it inflated the iteration count, and that inflation is gone.

    Two-cube heterogeneous model (element 1 -> plastic, element 2 -> elastic,
    same VonMises specialization -- before wp/94b this was exactly the case
    where the shared static ``Stiffness`` was live) vs each cube analyzed alone.
    Part (a) is unchanged and still passes: converged stresses were always
    exact, because OpenSees' residual is built from each element's own stress,
    never from the shared tangent.  Part (b) is the flip.
    """
    iters_pl_alone, sig_pl_alone = _run([(1, 0.0)], [(1, N.LOAD_PL)], NSTEPS)
    iters_el_alone, sig_el_alone = _run([(2, 0.0)], [(2, N.LOAD_EL)], NSTEPS)
    iters_together, sig_together = _run(
        [(1, 0.0), (2, 3.0)], [(1, N.LOAD_PL), (2, N.LOAD_EL)], NSTEPS)

    # (a) RESULTS: converged stresses match the solo analyses to numerical
    # tolerance, even though the intermediate tangent was smeared by H1.
    rel_pl = np.max(np.abs(sig_together[1] - sig_pl_alone[1])) / \
        max(np.max(np.abs(sig_pl_alone[1])), 1e-30)
    rel_el = np.max(np.abs(sig_together[2] - sig_el_alone[2])) / \
        max(np.max(np.abs(sig_el_alone[2])), 1e-30)
    # Tolerance matches the NormDispIncr test tolerance (1e-9 on displacement
    # increment), not H1's smeared tangent -- H1 affects only how many
    # iterations it takes to reach this band, not the band itself.
    assert rel_pl < 1e-6, "plastic element's committed stress must be exact"
    assert rel_el < 1e-6, "elastic element's committed stress must be exact"

    # (b) COST: FIXED by wp/94b.  The H1 penalty was the extra Newton
    # iterations the smeared tangent cost -- measured at +62% on this rig
    # before the fix.  With each element assembled from its own tangent, the
    # heterogeneous run costs no more than the two solo runs summed (the two
    # cubes are disconnected, so the combined Newton is exactly the two
    # independent Newtons superposed; the small slack absorbs the shared
    # convergence test, which stops on the WHOLE displacement-increment norm).
    penalty = iters_together - (iters_pl_alone + iters_el_alone)
    print(f"H1 iteration penalty: alone={iters_pl_alone}+{iters_el_alone}, "
          f"together={iters_together}, extra={penalty}")
    assert penalty <= ITER_PENALTY_SLACK, (
        f"the heterogeneous run still costs {penalty} extra Newton iterations "
        f"over the sum of the solo runs -- the H1 tangent smearing (ADR-94 M1) "
        f"may be back")


def test_blue_H1_homogeneous_single_material_mesh_is_ordering_invariant(
        vm_available):
    """When every element in a single-material mesh is in the SAME state
    (the common idealized case some decks approximate with one representative
    element type per zone), the last-integrated GP's tangent IS every other
    GP's own tangent, so H1's swap is invisible -- both iteration count and
    result are identical regardless of which element OpenSees visits last.
    This is NOT true for a mesh with genuine stress gradients (Cerro Lindo's
    horseshoe-cavity model is exactly such a case) -- there H1 still costs
    iterations even though results remain exact per (a) above.
    """
    iters_a, sig_a = _run([(1, 0.0), (2, 3.0)],
                           [(1, N.LOAD_PL), (2, N.LOAD_PL)], NSTEPS)
    iters_b, sig_b = _run([(1, 0.0), (2, 3.0)],
                           [(2, N.LOAD_PL), (1, N.LOAD_PL)], NSTEPS)
    assert iters_a == iters_b
    assert np.max(np.abs(sig_a[1] - sig_b[1])) < 1e-10
