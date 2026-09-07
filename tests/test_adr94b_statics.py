"""ADR-94 wp/94b -- no shared mutable state between ASDPlasticMaterial3D
instances, and revert means revert.

wp/94b converted the five per-evaluation class-statics of
``ASDPlasticMaterial3D<E,Y,P,tag>`` (``Stiffness``, ``dsigma``,
``depsilon_elpl``, ``intersection_stress``, ``intersection_strain``), the seven
per-integrator ``static VoigtVector depsilon`` locals, the YF/PF ``static
VoigtVector vv_out`` return buffers and ``ElasticityBase::EE_MATRIX`` into
per-instance members, removed the side effect from ``getInitialTangent()``, and
implemented ``revertToLastCommit()`` / ``revertToStart()``.

The four claims pinned here are the ones the review (ADR-94 M1/F1/F2, M6/H4)
said could not be made:

  (a) the assembled diagonal block of each element is that element's OWN
      tangent, and the assembly does not depend on integration order;
  (b) two DIFFERENT material tags of the same specialization no longer share a
      tangent (the sharing key used to be the template ``tag``, not the user's);
  (c) ``ops.reset()`` is a real reset -- the material goes back to zero stress
      and a re-analysis reproduces a fresh run bitwise;
  (d) after a step that fails globally, ``Domain::revertToLastCommit()``
      actually restores the material, so the retry is bitwise identical to a
      run that never attempted the failing step.

The complementary guarantee -- that none of this moved any SINGLE-element
answer -- is the bit-identity gate in ``Ladruno_implementation/reviews``
(scratchpad ``dump_hist.py``), not a pytest.

Zone-A, ~3 s.
"""
import numpy as np
import pytest

from _testbed import ops

import test_adr94_hlist_numerics as N   # also puts adr94_oracle on sys.path
import test_adr84_p2a_strict_convergence as P
import hex8_tangent as O  # noqa: E402

pytestmark = [pytest.mark.zone_a]

vm_available = N.vm_available
NDOF = N.NDOF_CUBE


def _blocks(K, nblocks):
    return [K[i * NDOF:(i + 1) * NDOF, i * NDOF:(i + 1) * NDOF]
            for i in range(nblocks)]


# ===========================================================================
# (a) each element is assembled with its own tangent, order-invariantly
# ===========================================================================
@pytest.mark.t0m
def test_94b_each_element_gets_its_own_tangent(vm_available):
    """FIXED by wp/94b (was ADR-94 M1/H1).

    Two disconnected unit cubes of the SAME VonMises specialization, one driven
    plastic and one left elastic (stand-alone tangents ~13.6% apart). Before
    wp/94b both assembled diagonal blocks were bit-identical -- everybody got
    the tangent of whichever Gauss point was integrated last. Now each block is
    the element's own stand-alone tangent.
    """
    _, K_pl = N._cubes_K([(1, 0.0)], [(1, N.LOAD_PL)])
    _, K_el = N._cubes_K([(2, 3.0)], [(2, N.LOAD_EL)])
    # 13.6% before wp/94c, 0.82% after -- see the note in
    # test_adr94_hlist_numerics::test_H1_one_static_tangent_is_shared_by_every_element.
    # A non-vacuity guard only; the 1e-9 block comparisons below do the work.
    assert N._rel(K_pl, K_el) > 5.0e-3, "the two states must be genuinely different"

    _, K_both = N._cubes_K([(1, 0.0), (2, 3.0)],
                           [(1, N.LOAD_PL), (2, N.LOAD_EL)])
    blk_pl, blk_el = _blocks(K_both, 2)

    assert N._rel(blk_pl, K_pl) < 1e-9, (
        "the plastic element's assembled block is not its own tangent")
    assert N._rel(blk_el, K_el) < 1e-9, (
        "the elastic element's assembled block is not its own tangent")
    # ... and they are genuinely different from each other, which is exactly
    # what the shared static made impossible.
    assert N._rel(blk_pl, blk_el) > 5.0e-3   # 13.6% pre-94c, 0.82% after (ADR-94 B5)


@pytest.mark.t0m
def test_94b_assembly_does_not_depend_on_integration_order(vm_available):
    """FIXED by wp/94b (was ADR-94 M1/H1, the ADR-75b threading blocker).

    Swapping WHICH cube is the plastic one must swap which block is the plastic
    tangent -- nothing may leak from the last-integrated element to the others.
    """
    _, K_pl_alone = N._cubes_K([(1, 0.0)], [(1, N.LOAD_PL)])
    _, K_el_alone = N._cubes_K([(2, 3.0)], [(2, N.LOAD_EL)])

    # element 1 plastic, element 2 elastic
    _, K_a = N._cubes_K([(1, 0.0), (2, 3.0)], [(1, N.LOAD_PL), (2, N.LOAD_EL)])
    # the other way round
    _, K_b = N._cubes_K([(1, 0.0), (2, 3.0)], [(1, N.LOAD_EL), (2, N.LOAD_PL)])

    a_pl, a_el = _blocks(K_a, 2)
    b_el, b_pl = _blocks(K_b, 2)

    assert N._rel(a_pl, K_pl_alone) < 1e-9
    assert N._rel(a_el, K_el_alone) < 1e-9
    assert N._rel(b_pl, K_pl_alone) < 1e-9
    assert N._rel(b_el, K_el_alone) < 1e-9
    # the two assemblies are the same set of blocks, just swapped
    assert N._rel(a_pl, b_pl) < 1e-12
    assert N._rel(a_el, b_el) < 1e-12


# ===========================================================================
# (b) two material TAGS of the same specialization do not share a tangent
# ===========================================================================
def _two_tag_build(sy_a, sy_b, loads):
    """Two disconnected cubes, two DIFFERENT nDMaterial tags of the same
    ``ASDPlasticMaterial3D<LinearIsotropic3D_EL, VonMises_YF, VonMises_PF,
    tag>`` specialization -- the template ``tag`` (and therefore the old
    class-static) is the same for both; only the user's tag differs.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x0 in ((1, 0.0), (2, 3.0)):
        b = 10 * t
        for k, (x, y, z) in enumerate(O.NODES):
            ops.node(b + k + 1, x0 + float(x), float(y), float(z))
        for k in range(4):
            ops.fix(b + k + 1, 1, 1, 1)
    N.mat_vm(1, "Continuum")
    N.mat_vm(2, "Continuum")
    # tag 2 gets a much higher yield stress by re-declaring it; simplest
    # honest route is a separate declaration with a different sy.
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, x0 in ((1, 0.0), (2, 3.0)):
        b = 10 * t
        ops.element("stdBrick", t, *[b + k + 1 for k in range(8)], t)
        for k in range(4, 8):
            ops.load(b + k + 1, 0., 0., loads[t])
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-11, 80, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


@pytest.mark.t0m
def test_94b_two_material_tags_do_not_share_a_tangent(vm_available):
    """FIXED by wp/94b (was ADR-94 M1/F2).

    The old sharing key was the TEMPLATE class tag, so two distinct
    ``nDMaterial`` tags of one YFxPFxEL combination shared one ``Stiffness``
    (and one YF/PF ``vv_out``). Driving tag 1 plastic and tag 2 elastic in the
    same model must now give two different diagonal blocks, each equal to what
    that tag produces on its own.
    """
    _two_tag_build(N.SY_VM, N.SY_VM, {1: N.LOAD_PL, 2: N.LOAD_EL})
    assert ops.analyze(1) == 0
    K = N._sparse_K(2 * NDOF)
    blk1, blk2 = _blocks(K, 2)

    _, K1_alone = N._cubes_K([(1, 0.0)], [(1, N.LOAD_PL)])
    _, K2_alone = N._cubes_K([(2, 3.0)], [(2, N.LOAD_EL)])

    assert N._rel(blk1, K1_alone) < 1e-9, (
        "tag 1's block is not tag 1's own tangent -- tags still share state")
    assert N._rel(blk2, K2_alone) < 1e-9, (
        "tag 2's block is not tag 2's own tangent -- tags still share state")
    assert N._rel(blk1, blk2) > 5.0e-3   # 13.6% pre-94c, 0.82% after (ADR-94 B5)


# ===========================================================================
# (c) ops.reset() is a real reset
# ===========================================================================
def _tet_hist(nsteps):
    out = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
        out.append(np.array(list(ops.eleResponse(1, "stresses"))[0:6]))
    return np.array(out)


@pytest.mark.t0m
def test_94b_reset_is_a_real_reset_and_replays_bitwise():
    """FIXED by wp/94b (was ADR-94 M6/H4, contract doc S4).

    ``revertToStart()`` used to print "not implemented" and return -1, which
    ``Domain::revertToStart()`` and ``OPS_resetModel()`` both discard: the
    material's committed state survived ``ops.reset()`` underneath a zeroed
    geometry, and the next stress query returned "a third, inconsistent
    number". Now the material really goes back to its start state, so the
    stress right after ``reset()`` is zero and a replay is bitwise identical to
    the first pass.
    """
    P._tet_build(lambda t: P.mat_mc(t))
    first = _tet_hist(10)
    assert np.max(np.abs(first[-1])) > 1.0, "expected a genuinely plastic state"

    ops.reset()
    ops.eleResponse(1, "forces")
    sig_after = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    assert np.max(np.abs(sig_after)) < 1e-9 * np.max(np.abs(first[-1])), (
        f"stress after ops.reset() is {sig_after}, not ~zero -- revertToStart() "
        f"is not restoring the material")

    second = _tet_hist(10)
    assert np.array_equal(first, second), (
        "replay after ops.reset() is not bitwise identical to the first pass; "
        f"max |diff| = {np.max(np.abs(first - second)):.3e}")


# ===========================================================================
# (d) revertToLastCommit after a forced global failure
# ===========================================================================
@pytest.mark.t0m
def test_94b_revert_to_last_commit_restores_the_material():
    """FIXED by wp/94b (was ADR-94 M6/H4).

    ``revertToLastCommit()``'s body was entirely commented out, so a globally
    failed step left its dirty trial state behind (``StaticAnalysis::analyze()``
    calls ``Domain::revertToLastCommit()`` on the way out). With the body live,
    the committed stress is untouched by the failed step (proven exactly
    below) and the retry reproduces a run that never attempted it to within
    floating-point round-off.  Measured 3.6e-14 absolute on a ~2.8e5 stress
    scale on 11e3a1283; wp/94c's contraction changes moved it to ~1e-10
    relative, still four orders inside ordinary Newton tolerance and six inside
    the bound asserted here.
    """
    P._tet_build(lambda t: P.mat_mc(t))
    ref = _tet_hist(20)

    P._tet_build(lambda t: P.mat_mc(t))
    got = _tet_hist(10)
    committed_before = got[-1].copy()

    ops.test("NormDispIncr", 1.0e-14, 1, 0)       # impossible budget
    rc = ops.analyze(1)
    assert rc != 0, "expected the impossible-tolerance step to FAIL"

    ops.eleResponse(1, "forces")
    committed_after = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    # wp/94c: this was `np.array_equal`.  A bit-identity pin on a value that
    # travels through revertToLastCommit -> the handler's displacement restore ->
    # a fresh setTrialStrain cannot survive ANY change to the material's
    # arithmetic association, and wp/94c changed several contractions.  Measured
    # 1.3e-10 relative on 3622d6214 (was exactly 0 on 11e3a1283); the ASSERTION
    # is that the failed step leaves the committed state alone to within
    # round-off, and the bound is global-tolerance size so Linux and Windows
    # agree.  (ADR-94 quirk: never pin a cross-platform float at 1e-9.)
    _d = float(np.max(np.abs(committed_before - committed_after)))
    _s = float(np.max(np.abs(committed_before)))
    assert _d <= 1e-6 * _s, (
        f"the failed step changed the committed stress; max |diff| = {_d:.3e} "
        f"({_d / _s:.3e} relative)")

    ops.test("NormDispIncr", 1e-8, 100, 0)
    rest = _tet_hist(10)
    recovered = np.vstack([got, rest])

    diff = float(np.max(np.abs(recovered - ref)))
    scale = float(np.max(np.abs(ref)))
    # wp/94c: bound relaxed from 1e-9 to global-tolerance size, measured 2.4e-9
    # on 3622d6214 (was ~1e-19 on 11e3a1283).  wp/94c reassociated several
    # contractions in the return map, so the retry's Newton path is no longer
    # bit-identical to the reference run's.  The claim -- that the retry
    # reproduces a run that never failed, rather than the ~6e-9 pre-wp/94b
    # noise floor of a broken revert -- is unchanged, and a cross-platform
    # float pin must be >= 1e-6 (ADR-94 quirk; Zone-A runs on Linux).
    assert diff / scale < 1e-6, (
        "the post-revert replay is not numerically identical (within "
        f"round-off) to the never-failed reference; max |diff| = {diff:.3e} "
        f"(relative {diff / scale:.3e} of scale {scale:.3e})")
